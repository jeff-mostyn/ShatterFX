#include "TetrahedralObject.h"

// -----------------------------------------------------
// -------------------- TETRAHEDRON --------------------
// -----------------------------------------------------
vec3 Tetrahedron::GetCenterOfMass()
{
	if (m_points.size() == 4)
	{
		return 0.25 * (m_points[0] + m_points[1] + m_points[2] + m_points[3]);
	}
	else
	{
		return vec3Zero;
	}
}

float Tetrahedron::TetrahedralVolume()
{
	Eigen::Matrix3f A;
	// Calculate the Jacobian
	// the columns of this matrix are the tetrahedrons 1st point coords, subtracted from its 2nd, 3rd, and 4th points coords.
	// the rows correspond to x, y, z
	vec3 diff10 = m_points[1] - m_points[0];
	vec3 diff20 = m_points[2] - m_points[0];
	vec3 diff30 = m_points[3] - m_points[0];

	A << diff10[0], diff20[0], diff30[0],
		diff10[1], diff20[1], diff30[1],
		diff10[2], diff20[2], diff30[2];

	// assign calculated Jacobian to struct member variable J
	J = A;

	// volume is the determinant of the above matrix (whatever the name for that matrix is) divided by 6
	V = J.determinant() / 6.0;
	return V;
}

void Tetrahedron::Draw(GU_Detail *gdp)
{
	for (int i = 0; i < 4; i++)
	{
		vec3 p1 = m_points[i];
		vec3 p2 = m_points[(i + 1) % 4];
		vec3 p3 = m_points[(i + 2) % 4];

		GU_PrimPoly *poly = GU_PrimPoly::build(gdp, 3, GU_POLY_CLOSED);
		GA_Offset ptoff1 = poly->getPointOffset(0);
		GA_Offset ptoff2 = poly->getPointOffset(1);
		GA_Offset ptoff3 = poly->getPointOffset(2);

		gdp->setPos3(ptoff1, UT_Vector3(p1[0], p1[1], p1[2]));
		gdp->setPos3(ptoff2, UT_Vector3(p2[0], p2[1], p2[2]));
		gdp->setPos3(ptoff3, UT_Vector3(p3[0], p3[1], p3[2]));
	}
}

void Tetrahedron::ComputeCoefficients() {
	a = vector<float>();
	b = vector<float>();
	c = vector<float>();
	d = vector<float>();

	// Use Cramer's rule to solve system of 4 equations of form AX = B
	// equations are the shape function N = a + bx + cy + dz
	// coefficients will be calculated by determinant(Amod) / determinant(A) where A =
	//
	// | 1 x0 y0 z0 |
	// | 1 x1 y1 z1 |
	// | 1 x2 y2 z2 |
	// | 1 x3 y3 z3 |
	//
	// and Amod is A, except one of its columns (which one depends on if we're calculating for coeff a, b, c, or, d)
	// is replaced by...
	//
	// B = [1 0 0 0], [0 1 0 0], [0 0 1 0], or [0 0 0 1] (which one depends on if we're getting coeffs for 1st, 2nd, 3rd, or 4th point)
	//
	Eigen::Matrix4f A;
	A << 1.f, m_points[0][0], m_points[0][1], m_points[0][2],
		1.f, m_points[1][0], m_points[1][1], m_points[1][2],
		1.f, m_points[2][0], m_points[2][1], m_points[2][2],
		1.f, m_points[3][0], m_points[3][1], m_points[3][2];

	// iterate to calculate a, b, c, d, for each point in m_points. i determines which point, inner loop j determines a, b, c, d
	// the location of 1 in vector B corresponds to the first, second, third, fourth points.
	for (int i = 0; i < 4; i++)
	{
		Eigen::Vector4f B = {0, 0, 0, 0};
		B[i] = 1;

		// iterate across matrix A, substituting each column with B. first column gets us a, 2nd b, 3rd c, etc.
		for (int j = 0; j < 4; j++)
		{
			Eigen::Matrix4f A_copy = A;
			A_copy.col(j) = B;

			switch (j)
			{
			case 0:
				a.push_back(A_copy.determinant() / A.determinant());
				break;
			case 1:
				b.push_back(A_copy.determinant() / A.determinant());
				break;
			case 2:
				c.push_back(A_copy.determinant() / A.determinant());
				break;
			case 3:
				d.push_back(A_copy.determinant() / A.determinant());
				break;
			}
		}
	}

	// with the coefficients, we can compute the strain-displacement matrix of the tetrahedron
	// it is 12 x 6 and generally of the form
	// | b1 00 00 b2 00 00 b3 00 00 b4 00 00 |
	// | 00 c1 00 00 c2 00 00 c3 00 00 c4 00 |
	// | 00 00 d1 00 00 d2 00 00 d3 00 00 d4 |
	// | c1 b1 00 c2 b2 00 c3 b3 00 c4 b4 00 |
	// | 00 d1 c1 00 d2 c2 00 d3 c3 00 d4 c4 |
	// | d1 00 b1 d2 00 b2 d3 00 b3 d4 00 b4 |
	Eigen::MatrixXf Strain_Disp{
		{ b[0], 0.00, 0.00, b[1], 0.00, 0.00, b[2], 0.00, 0.00, b[3], 0.00, 0.00 },
		{ 0.00, c[0], 0.00, 0.00, c[1], 0.00, 0.00, c[2], 0.00, 0.00, c[3], 0.00 },
		{ 0.00, 0.00, d[0], 0.00, 0.00, d[1], 0.00, 0.00, d[2], 0.00, 0.00, d[3] },
		{ c[0], b[0], 0.00, c[1], b[1], 0.00, c[2], b[2], 0.00, c[3], b[3], 0.00 },
		{ 0.00, d[0], c[0], 0.00, d[1], c[1], 0.00, d[2], c[2], 0.00, d[3], c[3] },
		{ d[0], 0.00, b[0], d[1], 0.00, b[1], d[2], 0.00, b[2], d[3], 0.00, b[3] }
	};

	B = Strain_Disp;

	// Calculate the stiffness matrix for this tetrahedron
	// K_e = V * B^T * D * B
	// V - Volume
	// B - Strain-displacement matrix
	// D - Material stiffness matrix
	K_e = V * B.transpose() * m_myObj->GetMaterialMatrix() * B;
}

// ------------------------------------------------------------
// -------------------- TETRAHEDRAL OBJECT --------------------
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

	//for (TetFragment* frag : m_frags)
	//{
	//	for (Tetrahedron* tet : m_frags->m_tets)
	//	{
	//		tet->Draw(gdp);
	//	}
	//}
	for (Tetrahedron* tet : m_frags[0]->m_tets)
	{
		tet->Draw(gdp);
	}
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