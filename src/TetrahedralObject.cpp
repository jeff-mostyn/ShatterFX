#include "TetrahedralObject.h"

TetrahedralObject::TetrahedralObject() {
	m_points = std::vector<vec3*>();
	m_tets = std::vector<Tetrahedron*>();
}

TetrahedralObject::~TetrahedralObject() {
	for (auto p : m_points) {
		delete p;
	}

	/*for (auto p : m_pointsSingleton) {
		delete p;
	}*/

	for (auto tet : m_tets) {
		delete tet;
	}
}

void TetrahedralObject::AddPoint(vec3* a_point) {
	m_points.push_back(a_point);
	if (m_pointSet.count(*a_point) == 0) {
		//m_pointsSingleton.push_back(a_point);
		m_pointSet.insert(*a_point);
	}
}

void TetrahedralObject::AddTet(std::vector<exint> a_indices) {
	Tetrahedron* tet = new Tetrahedron();
	tet->pointIndices = std::vector<exint>(a_indices);
}

void TetrahedralObject::DumpPoints() {
	for (vec3* vec : m_points) {
		delete vec;
	}
	m_points.clear();

	/*for (vec3* vec : m_pointsSingleton) {
		delete vec;
	}
	m_pointsSingleton.clear();*/
}

std::vector<vec3*> TetrahedralObject::GetPoints() {
	return m_points;
}

std::set<vec3> TetrahedralObject::GetPointsSingleton() {
	return m_pointSet;
}

std::vector<Tetrahedron*> TetrahedralObject::GetTets() {
	return m_tets;
}
