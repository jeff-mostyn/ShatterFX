#pragma once

#include <set>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <SYS/SYS_Types.h>
#include <SOP/SOP_Node.h>
#include <GU/GU_PrimPoly.h>

#include "Tetrahedron.h"
#include "vec.h"

// ------------------------------------------
// ---------- FORWARD DECLARATIONS ----------
// ------------------------------------------
class Tetrahedron;

// ----------------------------------
// ---------- DECLARATIONS ----------
// ----------------------------------
struct MaterialData {
	float stiffness;	// Young's Modulus
	float strainRatio;	// Poisson Ratio
};

struct TetFragment
{
	//vec3 m_min, m_max; // maintain minimum and maximum range of fragment
	std::vector<Tetrahedron*> m_tets;
	//vec3 GetMin();
	//vec3 GetMax();
};

class TetrahedralObject
{
public:
	TetrahedralObject(std::unique_ptr<MaterialData> a_matData);
	~TetrahedralObject();

	// manipulator functions
	void AddTet(std::vector<vec3> a_points);

	// output functions
	void DumpPoints();
	void Draw(GU_Detail* gdp);

	// accessor functions
	const std::vector<vec3> GetPointsSingleton();
	const std::vector<Tetrahedron *> GetTets();
	const Eigen::MatrixXf* GetMaterialMatrix();
	const Eigen::VectorXf* GetDisplacementVector();
	const int GetVertexIndex(vec3 a_vertex);
	vec3 GetMin();
	vec3 GetMax();
	void RegisterImpact(vec3 a_dir, float a_mag, vec3 a_location);
	
	// Generation/Math functions
	void GenerateFragments(float cellSize);
	void MoveFragments(float distanceFromCenter);
	void ComputeMaterialInformation();

private:
	std::vector<vec3> m_points;
	std::unordered_map<vec3, int> m_pointIndices;
	std::vector<Tetrahedron *> m_tets;
	std::vector<TetFragment*> m_frags;
	std::unique_ptr<MaterialData> m_matData;
	Eigen::MatrixXf m_materialMatrix;
	Eigen::SparseMatrix<float> m_globalStiffness;
	Eigen::VectorXf m_pointDisplacementVector;
	vec3 m_min, m_max; // maintain minimum and maximum ranges in the set

	void ComputeMaterialMatrix();
	void ComputeGlobalStiffnessMatrix();
	Eigen::VectorXf SolveFEM(const Eigen::VectorXf& a_force);
};