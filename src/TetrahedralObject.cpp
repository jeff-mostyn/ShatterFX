#include "TetrahedralObject.h"

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
	tet->points = std::vector<vec3>(a_points);

	m_tets.push_back(tet);

	for (const vec3 &point : tet->points) {
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
