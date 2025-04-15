#include <iostream>

#include "TetrahedralObject_Simple.h"

// ------------------------------------------------------------
// 
// -------------------- TETRAHEDRAL OBJECT --------------------
// 
// ------------------------------------------------------------

TetrahedralObject_Simple::TetrahedralObject_Simple() : m_min(FLT_MAX, FLT_MAX, FLT_MAX), m_max(FLT_MIN, FLT_MIN, FLT_MIN)
{
	m_tets = std::vector<Tetrahedron_Simple*>();
}

TetrahedralObject_Simple::~TetrahedralObject_Simple()
{
	for (auto tet : m_tets)
	{
		delete tet;
	}
}

void checkTetMinMax(Tetrahedron_Simple* tet, vec3& min, vec3& max)
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

void TetrahedralObject_Simple::AddTet(std::vector<vec3> a_points)
{
	Tetrahedron_Simple* tet = new Tetrahedron_Simple(a_points, this);

	m_tets.push_back(tet);

	// check whether each point in the tetrahedron has been accounted for yet
	// if not, add it and its index to unordered_map for tracking, and add it to the vector of points
	for (const vec3& point : tet->m_points) {
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

void TetrahedralObject_Simple::DumpPoints() {
	m_points.clear();
	m_pointIndices.clear();
	m_min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	m_max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);
}

void TetrahedralObject_Simple::Draw(GU_Detail* gdp)
{
	for (Tetrahedron_Simple*tet : m_tets)
	{
		tet->Draw(gdp);
	}
	//for (Tetrahedron* tet : m_frags[0]->m_tets)
	//{
	//	tet->Draw(gdp);
	//}
}

const std::vector<vec3> TetrahedralObject_Simple::GetPointsSingleton() {
	return m_points;
}

const std::vector<Tetrahedron_Simple*> TetrahedralObject_Simple::GetTets() {
	return m_tets;
}

const int TetrahedralObject_Simple::GetVertexIndex(vec3 a_vertex) {
	return m_pointIndices[a_vertex];
}

vec3 TetrahedralObject_Simple::GetMin()
{
	return m_min;
}

vec3 TetrahedralObject_Simple::GetMax()
{
	return m_max;
}