#include "TetrahedralObject.h"

TetrahedralObject::TetrahedralObject() : m_min(FLT_MAX, FLT_MAX, FLT_MAX), m_max(FLT_MIN, FLT_MIN, FLT_MIN)
{
	m_tets = std::vector<Tetrahedron*>();
	//m_min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	//m_max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);
}

TetrahedralObject::~TetrahedralObject() {
	for (auto tet : m_tets) {
		delete tet;
	}
}

void checkTetMinMax(Tetrahedron* tet, vec3 &min, vec3 &max) {
	min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);
	for (vec3 point : tet->points)
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
	Tetrahedron* tet = new Tetrahedron();
	tet->points = std::vector<vec3>(a_points);

	m_tets.push_back(tet);

	for (const vec3 &point : tet->points) {
		if (m_pointSet.count(point) == 0) {
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

void TetrahedralObject::DumpPoints() {
	m_pointSet.clear();
	m_min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	m_max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);
}

std::set<vec3> TetrahedralObject::GetPointsSingleton() {
	return m_pointSet;
}

std::vector<Tetrahedron*> TetrahedralObject::GetTets() {
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