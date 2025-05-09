#pragma once

#include <Eigen/Dense>
#include <vector>

#include "vec.h"
#include "TetrahedralObject.h"

// ------------------------------------------
// ---------- FORWARD DECLARATIONS ----------
// ------------------------------------------
class TetrahedralObject;

// ----------------------------------
// ---------- DECLARATIONS ----------
// ----------------------------------
struct Tetrahedron
{
	TetrahedralObject* m_myObj;
	std::vector<vec3> m_points;
	std::vector<float> a, b, c, d; 	// coefficients for shape function
	Eigen::Matrix3f J;				// Jacobian
	Eigen::MatrixXf B;				// strain-displacement matrix
	Eigen::MatrixXf K_e;			// Element Stiffness Matrix
	Eigen::Matrix3f m_StrainTensor;	// Tetrahedron strain tensor
	Eigen::Matrix3f m_StressTensor;	// Tetrahedron stress tensor
	float m_W = 0;
	float m_V = 0;					// Volume

	Tetrahedron(std::vector<vec3> a_points, TetrahedralObject* a_obj)
	{
		m_myObj = a_obj;
		m_points = std::vector<vec3>(a_points);
		TetrahedralVolume();
		ComputeCoefficients();
	}
	vec3 GetCenterOfMass();
	Eigen::Vector3f GetVertexDisplacements(int index);
	float TetrahedralVolume();
	void Draw(GU_Detail* gdp, int fragmentId, GA_RWHandleS nameHandle);

	// computation functions
	void ComputeCoefficients();
	void ComputeStrainTensor();
	void ComputeStrainEnergy();
};