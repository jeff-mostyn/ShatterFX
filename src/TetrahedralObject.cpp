#include "TetrahedralObject.h"

// ------------------------------------------------------------
// 
// -------------------- TETRAHEDRAL OBJECT --------------------
// 
// ------------------------------------------------------------

TetrahedralObject::TetrahedralObject(std::unique_ptr<MaterialData> a_matData) : m_min(FLT_MAX, FLT_MAX, FLT_MAX), m_max(FLT_MIN, FLT_MIN, FLT_MIN)
{
	m_tets = std::vector<Tetrahedron *>();
	// m_min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	// m_max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);

	m_matData = std::move(a_matData);
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
	//for (Tetrahedron *tet : m_tets)
	//{
	//	tet->Draw(gdp);
	//}

	for (TetFragment* frag : m_frags)
	{
		for (Tetrahedron* tet : frag->m_tets)
		{
			tet->Draw(gdp);
		}
	}
	//for (Tetrahedron* tet : m_frags[0]->m_tets)
	//{
	//	tet->Draw(gdp);
	//}
}

const std::vector<vec3> TetrahedralObject::GetPointsSingleton() {
	return m_points;
}

const std::vector<Tetrahedron *> TetrahedralObject::GetTets() {
	return m_tets;
}

const Eigen::MatrixXf TetrahedralObject::GetMaterialMatrix() {
	return m_materialMatrix;
}

vec3 TetrahedralObject::GetMin()
{
	return m_min;
}

vec3 TetrahedralObject::GetMax()
{
	return m_max;
}

void TetrahedralObject::GenerateFragments(float cellSize)
{
	m_frags.clear();

	// Map to group tetrahedra by their voxel grid cell
	std::unordered_map<std::string, TetFragment*> cellToFragment;

	for (Tetrahedron* tet : m_tets)
	{
		vec3 center = tet->GetCenterOfMass();

		// Snap to grid
		int cx = static_cast<int>((center[0] - m_min[0]) / cellSize);
		int cy = static_cast<int>((center[1] - m_min[1]) / cellSize);
		int cz = static_cast<int>((center[2] - m_min[2]) / cellSize);

		// Create a unique key for this grid cell
		std::string key = std::to_string(cx) + "_" + std::to_string(cy) + "_" + std::to_string(cz);

		if (cellToFragment.find(key) == cellToFragment.end())
		{
			cellToFragment[key] = new TetFragment();
			m_frags.push_back(cellToFragment[key]);
		}

		cellToFragment[key]->m_tets.push_back(tet);
	}

}

void TetrahedralObject::MoveFragments(float distanceFromCenter)
{
	// Step 1: Compute center of the full mesh
	vec3 objectCenter = 0.5f * (m_min + m_max);

	// Step 2: Loop through fragments
	for (TetFragment* frag : m_frags)
	{
		// Compute average center of all tets in the fragment
		vec3 fragCenter = vec3Zero;
		int count = 0;
		for (Tetrahedron* tet : frag->m_tets)
		{
			fragCenter += tet->GetCenterOfMass();
			count++;
		}
		if (count == 0) continue;
		fragCenter /= static_cast<float>(count);

		// Step 3: Compute direction from object center to fragment center
		vec3 direction = fragCenter - objectCenter;
		float length = direction.Length();
		if (length == 0.0f) continue; // Prevent divide by zero
		direction /= length; // Normalize

		// Step 4: Move each point in the fragment
		for (Tetrahedron* tet : frag->m_tets)
		{
			for (vec3& pt : tet->m_points)
			{
				pt += direction * distanceFromCenter;
			}
		}
	}
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

void TetrahedralObject::ComputeMaterialMatrix() {
	float E = m_matData->stiffness;
	float v = m_matData->strainRatio;

	// this is the 3D version of the material matrix D for Hooke's Law
	Eigen::MatrixXf D {
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