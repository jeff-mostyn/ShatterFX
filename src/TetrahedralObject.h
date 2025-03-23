#include <set>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <SYS/SYS_Types.h>
#include <SOP/SOP_Node.h>
#include <GU/GU_PrimPoly.h>

#include "vec.h"

// ------------------------------------------
// ---------- FORWARD DECLARATIONS ----------
// ------------------------------------------
class TetrahedralObject;

// ----------------------------------
// ---------- DECLARATIONS ----------
// ----------------------------------
struct MaterialData {
	float stiffness;	// Young's Modulus
	float strainRatio;	// Poisson Ratio
};

struct Tetrahedron
{
	TetrahedralObject* m_myObj;
	std::vector<vec3> m_points;
	std::vector<float> a, b, c, d; 	// coefficients for shape function
	Eigen::Matrix3f J;				// Jacobian
	Eigen::MatrixXf B;				// strain-displacement matrix
	Eigen::MatrixXf K_e;			// Element Stiffness Matrix
	double V = 0;					// Volume

	Tetrahedron(std::vector<vec3> a_points, TetrahedralObject* a_obj)
	{
		m_myObj = a_obj;
		m_points = std::vector<vec3>(a_points);
		TetrahedralVolume();
		ComputeCoefficients();
	}
	vec3 GetCenterOfMass();
	double TetrahedralVolume();
	void Draw(GU_Detail *gdp);
	void ComputeCoefficients();
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
	const Eigen::MatrixXf GetMaterialMatrix();
	vec3 GetMin();
	vec3 GetMax();
	
	// Generation/Math functions
	void GenerateFragments(float cellSize);
	void ComputeGlobalStiffnessMatrix();

private:
	std::vector<vec3> m_points;
	std::unordered_map<vec3, int> m_pointIndices;
	std::vector<Tetrahedron *> m_tets;
	std::vector<TetFragment*> m_frags;
	std::unique_ptr<MaterialData> m_matData;
	Eigen::VectorXf m_materialMatrix;
	Eigen::SparseMatrix<float> K_global;
	vec3 m_min, m_max; // maintain minimum and maximum ranges in the set

	void ComputeMaterialMatrix();
};