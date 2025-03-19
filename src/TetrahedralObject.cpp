#include "TetrahedralObject.h"


double Tetrahedron::TetrahedralVolume() {
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


TetrahedralObject::TetrahedralObject() {
	m_tets = std::vector<Tetrahedron*>();
}

TetrahedralObject::~TetrahedralObject() {
	for (auto tet : m_tets) {
		delete tet;
	}
}

void TetrahedralObject::AddTet(std::vector<vec3> a_points) {
	Tetrahedron* tet = new Tetrahedron();
	tet->m_points = std::vector<vec3>(a_points);

	m_tets.push_back(tet);

	for (const vec3 &point : tet->m_points) {
		if (m_pointSet.count(point) == 0) {
			m_pointSet.insert(point);
		}
	}
}

void TetrahedralObject::DumpPoints() {
	m_pointSet.clear();
}

std::set<vec3> TetrahedralObject::GetPointsSingleton() {
	return m_pointSet;
}

std::vector<Tetrahedron*> TetrahedralObject::GetTets() {
	return m_tets;
}
