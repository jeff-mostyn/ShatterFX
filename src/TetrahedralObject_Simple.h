#pragma once

#include <set>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <SYS/SYS_Types.h>
#include <SOP/SOP_Node.h>
#include <GU/GU_PrimPoly.h>

#include "vec.h"

class TetrahedralObject_Simple;

// ------------------------------------------
// ---------- FORWARD DECLARATIONS ----------
// ------------------------------------------
struct Tetrahedron_Simple {
	TetrahedralObject_Simple* m_myObj;
	std::vector<vec3> m_points;

	Tetrahedron_Simple(std::vector<vec3> a_points, TetrahedralObject_Simple* a_obj)
	{
		m_myObj = a_obj;
		m_points = std::vector<vec3>(a_points);
	}

	void Draw(GU_Detail* gdp) {
		for (int i = 0; i < 4; i++)
		{
			vec3 p1 = m_points[i];
			vec3 p2 = m_points[(i + 1) % 4];
			vec3 p3 = m_points[(i + 2) % 4];

			GU_PrimPoly* poly = GU_PrimPoly::build(gdp, 3, GU_POLY_CLOSED);
			GA_Offset ptoff1 = poly->getPointOffset(0);
			GA_Offset ptoff2 = poly->getPointOffset(1);
			GA_Offset ptoff3 = poly->getPointOffset(2);

			gdp->setPos3(ptoff1, UT_Vector3(p1[0], p1[1], p1[2]));
			gdp->setPos3(ptoff2, UT_Vector3(p2[0], p2[1], p2[2]));
			gdp->setPos3(ptoff3, UT_Vector3(p3[0], p3[1], p3[2]));
		}
	}
};

// ----------------------------------
// ---------- DECLARATIONS ----------
// ----------------------------------

class TetrahedralObject_Simple
{
public:
	TetrahedralObject_Simple();
	~TetrahedralObject_Simple();

	// manipulator functions
	void AddTet(std::vector<vec3> a_points);

	// output functions
	void DumpPoints();
	void Draw(GU_Detail* gdp);

	// accessor functions
	const std::vector<vec3> GetPointsSingleton();
	const std::vector<Tetrahedron_Simple*> GetTets();
	const int GetVertexIndex(vec3 a_vertex);
	vec3 GetMin();
	vec3 GetMax();

private:
	std::vector<vec3> m_points;
	std::unordered_map<vec3, int> m_pointIndices;
	std::vector<Tetrahedron_Simple*> m_tets;
	vec3 m_min, m_max; // maintain minimum and maximum ranges in the set
};