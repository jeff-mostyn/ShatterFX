#include "TetrahedralObject.h"

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

double Tetrahedron::TetrahedralVolume()
{
	// the columns of this matrix are the tetrahedrons 1st point coords, subtracted from its 2nd, 3rd, and 4th points coords.
	// the rows correspond to x, y, z
	Eigen::Matrix3f A;

	vec3 diff10 = m_points[1] - m_points[0];
	vec3 diff20 = m_points[2] - m_points[0];
	vec3 diff30 = m_points[3] - m_points[0];

	A << diff10[0], diff20[0], diff30[0],
		diff10[1], diff20[1], diff30[1],
		diff10[2], diff20[2], diff30[2];

	// volume is the determinant of the above matrix (whatever the name for that matrix is) divided by 6
	return A.determinant() / 6.0;
}

void Tetrahedron::ComputeCoefficients() {
	a = vector<double>();
	b = vector<double>();
	c = vector<double>();
	d = vector<double>();

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
	for (int i = 0; i < 4; i++) {
		Eigen::Vector4f B = { 0, 0, 0, 0 };
		B[i] = 1;

		// iterate across matrix A, substituting each column with B. first column gets us a, 2nd b, 3rd c, etc.
		for (int j = 0; j < 4; j++) {
			Eigen::Matrix4f A_copy = A;
			A_copy.col(j) = B;

			switch (j) {
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
}


TetrahedralObject::TetrahedralObject() : m_min(FLT_MAX, FLT_MAX, FLT_MAX), m_max(FLT_MIN, FLT_MIN, FLT_MIN)
{
	m_tets = std::vector<Tetrahedron *>();
	// m_min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	// m_max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);
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

void TetrahedralObject::AddTet(std::vector<vec3> a_points) {
	Tetrahedron* tet = new Tetrahedron(a_points);

	m_tets.push_back(tet);

	for (const vec3 &point : tet->m_points)
	{
		if (m_pointSet.count(point) == 0)
		{
			m_pointSet.insert(point);
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

void TetrahedralObject::DumpPoints()
{
	m_pointSet.clear();
	m_min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	m_max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);
}

std::set<vec3> TetrahedralObject::GetPointsSingleton()
{
	return m_pointSet;
}

std::vector<Tetrahedron *> TetrahedralObject::GetTets()
{
	return m_tets;
}

vec3 TetrahedralObject::GetMin()
{
	return m_min;
}

vec3 TetrahedralObject::GetMax()
{
	return m_max;
}