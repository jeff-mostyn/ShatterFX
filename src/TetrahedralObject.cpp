#include <iostream>

#include "TetrahedralObject.h"

// ------------------------------------------------------------
// 
// -------------------- TETRAHEDRAL OBJECT --------------------
// 
// ------------------------------------------------------------

TetrahedralObject::TetrahedralObject(
	std::unique_ptr<MaterialData> a_matData, 
	float a_energyConsumptionPercent,
	float a_energySpreadFactor
) : m_min(FLT_MAX, FLT_MAX, FLT_MAX), m_max(FLT_MIN, FLT_MIN, FLT_MIN)
{
	m_tets = std::vector<Tetrahedron *>();
	// m_min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	// m_max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);

	m_matData = std::move(a_matData);
	m_energyConsumptionPercent = a_energyConsumptionPercent;
	m_energySpreadFactor = a_energySpreadFactor;
	ComputeMaterialMatrix();
}

TetrahedralObject::~TetrahedralObject()
{
	for (auto tet : m_tets)
	{
		delete tet;
	}
}

void checkTetMinMax(Tetrahedron *tet, vec3 &min, vec3 &max)
{
	min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);
	for (vec3 point : tet->m_points)
	{
		min[0] = point[0] < min[0] ? point[0] : min[0];
		min[1] = point[1] < min[1] ? point[1] : min[1];
		min[2] = point[2] < min[2] ? point[2] : min[2];
		max[0] = point[0] > max[0] ? point[0] : max[0];
		max[1] = point[1] > max[1] ? point[1] : max[1];
		max[2] = point[2] > max[2] ? point[2] : max[2];
	}
}

void TetrahedralObject::AddTet(std::vector<vec3> a_points)
{
	Tetrahedron *tet = new Tetrahedron(a_points, this);

	m_tets.push_back(tet);

	// check whether each point in the tetrahedron has been accounted for yet
	// if not, add it and its index to unordered_map for tracking, and add it to the vector of points
	for (const vec3 &point : tet->m_points) {
		if (!m_pointIndices.count(point)) {
			int nextIndex = m_points.size();
			m_pointIndices[point] = nextIndex;
			m_points.push_back(point);
		}
	}

	vec3 min, max;
	checkTetMinMax(tet, min, max);
	m_min[0] = min[0] < m_min[0] ? min[0] : m_min[0];
	m_min[1] = min[1] < m_min[1] ? min[1] : m_min[1];
	m_min[2] = min[2] < m_min[2] ? min[2] : m_min[2];
	m_max[0] = max[0] > m_max[0] ? max[0] : m_max[0];
	m_max[1] = max[1] > m_max[1] ? max[1] : m_max[1];
	m_max[2] = max[2] > m_max[2] ? max[2] : m_max[2];
}

void TetrahedralObject::DumpPoints() {
	m_points.clear();
	m_pointIndices.clear();
	m_min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	m_max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);
}

void TetrahedralObject::Draw(GU_Detail* gdp)
{
	// name handle is so we can assign each tet a name attribute for its fragment
	// this should allow each fragment to function as an independent group for physics
	GA_RWHandleS nameHandle(gdp->addStringTuple(GA_ATTRIB_PRIMITIVE, "name", 1));

	int i = 0;
	for (TetFragment* frag : m_frags) {
		for (Tetrahedron* tet : frag->m_tets) {
			tet->Draw(gdp, i, nameHandle);
		}

		i++;
	}
}

const std::vector<vec3> TetrahedralObject::GetPointsSingleton() {
	return m_points;
}

const std::vector<Tetrahedron *> TetrahedralObject::GetTets() {
	return m_tets;
}

const Eigen::MatrixXf* TetrahedralObject::GetMaterialMatrix() {
	return &m_materialMatrix;
}

const Eigen::VectorXf* TetrahedralObject::GetDisplacementVector() {
	return &m_pointDisplacementVector;
}

const int TetrahedralObject::GetVertexIndex(vec3 a_vertex) {
	return m_pointIndices[a_vertex];
}

vec3 TetrahedralObject::GetMin()
{
	return m_min;
}

vec3 TetrahedralObject::GetMax()
{
	return m_max;
}

float TetrahedralObject::GetTotalEnergy() {
	float sum = 0;
	for (auto& tet : m_tets) {
		sum += tet->m_W;
	}

	return sum;
}

/// <summary>
/// Entry point to the calculation. This will start function flow for running fracture sims.
/// </summary>
/// <param name="a_dir">Force direction</param>
/// <param name="a_mag">Impact force magnitude</param>
/// <param name="a_location">Desired impact location in 3D space. Closest vertex will be determined</param>
void TetrahedralObject::RegisterImpact(vec3 a_dir, float a_mag, vec3 a_location) {
	// input sanitization
	a_dir = a_dir.Normalize();

	// computation
	Eigen::VectorXf f_global = Eigen::VectorXf::Zero(3 * m_points.size());

	// get closest point to impact
	int closestIndex = -1;
	float minDist = FLT_MAX;

	for (int i = 0; i < m_points.size(); ++i) {
		float dist = (m_points[i] - a_location).SqrLength();

		if (dist < minDist) {
			minDist = dist;
			closestIndex = i;
		}
	}

	// Apply force at that point
	if (closestIndex >= 0) {
		vec3 force = a_dir * a_mag;

		f_global(3 * closestIndex + 0) = force[0];
		f_global(3 * closestIndex + 1) = force[1];
		f_global(3 * closestIndex + 2) = force[2];
	}

	// Compute vector of vertex displacements
	m_pointDisplacementVector = SolveFEM(f_global);

	//std::cout << "Displacements (first few): " << m_pointDisplacementVector.head(10).transpose() << std::endl;

	for (Tetrahedron* tet : m_tets) {
		tet->ComputeStrainTensor();
		tet->ComputeStrainEnergy();
	}
}

void TetrahedralObject::GenerateFragments(std::vector<vec3> sites) {
	m_frags.clear();

	// Map to group tetrahedra by their voxel grid cell
	std::unordered_map<int, TetFragment*> cellToFragment;

	for (Tetrahedron* tet : m_tets)
	{
		int closest = -1;
		float closestDist = FLT_MAX;
		for (int i = 0; i < sites.size(); i++) {
			float dist = (tet->GetCenterOfMass() - sites[i]).SqrLength();
			if (dist < closestDist) { 
				closest = i; 
				closestDist = dist;
			}
		}

		if (cellToFragment.find(closest) == cellToFragment.end())
		{
			cellToFragment[closest] = new TetFragment();
			m_frags.push_back(cellToFragment[closest]);
		}

		cellToFragment[closest]->m_tets.push_back(tet);
	}

}

void TetrahedralObject::ComputeMaterialInformation() {
	ComputeMaterialMatrix();
	ComputeGlobalStiffnessMatrix();
}

std::vector<vec3> TetrahedralObject::GenerateFractureSites(vec3 a_impactPoint) {
	const int maxSites = 50; // clamp to avoid infinite loop
	std::vector<vec3> sites;
	float totalEnergy = GetTotalEnergy() * m_energyConsumptionPercent;
	float fractureEpislon = 0.005;

	std::cout << "Energy to be spent: " << totalEnergy << std::endl;

	// Step 1: Sample initial sites using weighted energy distribution
	std::list<Fracture> candidates;
	std::map<vec3, float> fracSiteEnergyContrib;
		
	// Generate initial candidates
	for (auto& tet : m_tets) {
		vec3 side1, side2;
		side1 = tet->m_points[1] - tet->m_points[0];
		side2 = tet->m_points[2] - tet->m_points[0];

		float sideArea = 0.5 * side1.Cross(side2).Length();

		// if deformation energy in tetrahedron is enough to fracture, log the fracture point
		// and compute distance weighted energy
		if (tet->m_W > sideArea * m_matData->fractureToughness) {
			candidates.push_back(Fracture(tet->GetCenterOfMass(), tet->m_W));
		}
	}


	// Step 2: Calculate energy contributed by each tet to each fracture
	// and use to limit number of fractures to those such that the energy spent does not
	// exceed the total energy of the system, but that the fractures that happen are the highest
	// energy ones
	float boundingRadius = max(m_max.Length(), m_min.Length());
	float sigma_global = boundingRadius * m_energySpreadFactor; // affects how far energy spreads from impact

	float radiiTotal = 0.f;
	for (const auto& tet : m_tets) {
		radiiTotal += (tet->m_points[0] - tet->GetCenterOfMass()).Length();
	}
	float avgRadius = radiiTotal / m_tets.size();
	float sigma_local = 2.f * avgRadius; // affects how many sites a tet can influence

	// initialize map
	for (const auto &candidate : candidates) {
		fracSiteEnergyContrib[candidate.location] = 0;
	}

	// calculate combined energy applied to each fracture site from nearby tets
	for (const auto& tet : m_tets) {
		float distFromImpact = (tet->GetCenterOfMass() - a_impactPoint).Length();

		// Global falloff based on impact distance (e.g. Gaussian)
		float globalWeight = std::exp(-(distFromImpact * distFromImpact) / (2.0f * sigma_global * sigma_global));
	
		float weightedEnergy = tet->m_W * globalWeight;

		// for each fracture, assign energy contributed by this tet
		for (std::list<Fracture>::iterator it = candidates.begin(); it != candidates.end(); ++it) {
			float d_local = (tet->GetCenterOfMass() - (*it).location).Length();

			// Local kernel (e.g. Gaussian again, or could be inverse-square)
			float localWeight = std::exp(-(d_local * d_local) / (2.0f * sigma_local * sigma_local));

			fracSiteEnergyContrib[(*it).location] += weightedEnergy * localWeight;
		}
	}

	// sort candidates by energy
	float currentHighestEnergy = FLT_MAX;
	float nextHighestEnergy = 0;
	vec3 nextHighestCandidate = vec3Zero;
	std::vector<vec3> sortedCandidates = std::vector<vec3>();
	std::vector<float> sortedEnergies = std::vector<float>();
	
	for (int i = 0; i < candidates.size(); i++) {
		for (auto& candidate : candidates) {
			float candEnergy = fracSiteEnergyContrib[candidate.location];
			if (candEnergy < currentHighestEnergy && candEnergy >= nextHighestEnergy) {
				nextHighestEnergy = candEnergy;
				nextHighestCandidate = candidate.location;
			}
		}

		currentHighestEnergy = nextHighestEnergy;
		sortedCandidates.push_back(nextHighestCandidate);
		sortedEnergies.push_back(nextHighestEnergy);

		nextHighestEnergy = 0;
	}

	// pare "extra energy" fractures if any exist
	float contributedEnergy = 0;
	for (int i = 0; i < sortedCandidates.size(); i++) {
		std::cout << "Energy expenditure [" << i << "]: " << sortedEnergies[i] << std::endl;
		// if total contributed energy plus energy of next fracture is less than the total
		// system energy, add the current fracture to the sites
		if (contributedEnergy + sortedEnergies[i] <= totalEnergy && sortedEnergies[i] > fractureEpislon) {
			std::cout << "Energy spent!" << std::endl;
			sites.push_back(sortedCandidates[i]);
			contributedEnergy += sortedEnergies[i];
		}
		else {
			break;
		}
	}

	// Step 3: Run Lloyd’s algorithm to compute CVD in 3D
	//ComputeCVD(sites, 3); // 5 iterations of Lloyd

	// Step 4: Assign each tet to nearest site
	for (int t = 0; t < m_tets.size(); ++t) {
		vec3 com = m_tets[t]->GetCenterOfMass();
		float minDist2 = FLT_MAX;
		int closest = 0;

		int i = 0;
		for (auto& sites : sites) {
			float d2 = (com - sites).SqrLength();
			if (d2 < minDist2) {
				minDist2 = d2;
				closest = i;
			}
			i++;
		}
	}

	return sites;
}

void TetrahedralObject::ComputeMaterialMatrix() {
	float E = m_matData->stiffness;
	float v = m_matData->strainRatio;

	// this is the 3D version of the material matrix D for Hooke's Law
	Eigen::MatrixXf D{
		{ 1.f - v, v, v, 0.f, 0.f, 0.f },
		{ v, 1.f - v, v, 0.f, 0.f, 0.f },
		{ v, v, 1.f - v, 0.f, 0.f, 0.f },
		{ 0.f, 0.f, 0.f, (1.f - (2.f * v) / 2.f), 0.f, 0.f },
		{ 0.f, 0.f, 0.f, 0.f, (1.f - (2.f * v) / 2.f), 0.f },
		{ 0.f, 0.f, 0.f, 0.f, 0.f, (1.f - (2.f * v) / 2.f) }
	};
	D *= (E / ((1.f + v) * (1.f - (2.f * v))));

	m_materialMatrix = D;
}

/// <summary>
/// Computes global stiffness matrix for the tetrahedral object. It is important that this is only called after all tetrahedrons are added and processed.
/// Additionally, ONLY CALL WHEN NECESSARY. This is an expensive function
/// </summary>
void TetrahedralObject::ComputeGlobalStiffnessMatrix() {
	// the size of the global stiffness matrix is 3n x 3n, where n is the number of points in the construct
	// this represents the degrees of freedom (x, y, z) for each point
	int degreeOfFreedomCount = m_points.size() * 3;

	Eigen::SparseMatrix<float> K_global(degreeOfFreedomCount, degreeOfFreedomCount);

	// triplets are a data structure that represents row, column, and value for a place in a matrix
	// they are good for inserting data into a sparse matrix all at once, so we'll collect everything in a vector of them and then insert
	std::vector<Eigen::Triplet<float>> triplets;

	for (size_t t = 0; t < m_tets.size(); t++) {
		const Eigen::MatrixXf& K_local = m_tets[t]->K_e;
		const std::vector<vec3>& tetVertices = m_tets[t]->m_points;

		// Loop through the local stiffness matrix and map to global
		// i and j are up to 4, because they are multiplied by 3 to get the "start" of each point's data in the local stiffness matrix, which is 12x12
		// those are 12x12 because the tetrahedrons have 4 points * 3 degrees of freedom.
		// di and dj are "offsets" that will specify the x, y, or z of each point
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				for (int di = 0; di < 3; ++di) {
					for (int dj = 0; dj < 3; ++dj) {
						// based on the index of the tetrahedrons' vertices in the vertex list for the whole object
						// assign a row and column in the global matrix for this data to be stored in
						int global_row = m_pointIndices[tetVertices[i]] * 3 + di;
						int global_col = m_pointIndices[tetVertices[j]] * 3 + dj;

						float value = K_local(i * 3 + di, j * 3 + dj);
						if (value != 0.0) {
							triplets.emplace_back(global_row, global_col, value);
						}
					}
				}
			}
		}
	}

	// Set the values in the sparse matrix
	// setFromTriplets will apparently automatically handle accumulation, if that is needed
	K_global.setFromTriplets(triplets.begin(), triplets.end());

	m_globalStiffness = K_global;
}

void TetrahedralObject::ComputeCVD(std::vector<vec3>& sites, int iterations) {
	for (int iter = 0; iter < iterations; ++iter) {
		
		// Assign tets to nearest site
		std::vector<std::vector<Tetrahedron*>> clusters(sites.size());

		for (Tetrahedron* tet : m_tets) {
			vec3 center = tet->GetCenterOfMass();

			int closestIndex = 0;
			float minDist = FLT_MAX;

			for (int i = 0; i < sites.size(); ++i) {
				float dist = (center - sites[i]).SqrLength();

				if (dist < minDist) {
					minDist = dist;
					closestIndex = i;
				}
			}

			clusters[closestIndex].push_back(tet);
		}

		// Update each site to the average center of mass of its assigned tets
		for (int i = 0; i < sites.size(); ++i) {
			if (clusters[i].empty()) {
				continue;
			}

			vec3 avg = vec3Zero;

			for (Tetrahedron* tet : clusters[i]) {
				avg += tet->GetCenterOfMass();
			}
			sites[i] = avg / static_cast<float>(clusters[i].size());
		}
	}
}

Eigen::VectorXf TetrahedralObject::SolveFEM(const Eigen::VectorXf& a_Force) {
	Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;
	solver.analyzePattern(m_globalStiffness);
	solver.factorize(m_globalStiffness);
	return solver.solve(a_Force);
};