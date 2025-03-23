#include <set>
#include <vector>

#include <Eigen/Dense>
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
	TetrahedralObject* myObj;
	std::vector<vec3> m_points;
	std::vector<float> a, b, c, d; 	// coefficients for shape function
	Eigen::Matrix3f J;				// Jacobian
	Eigen::MatrixXf B;				// strain-displacement matrix
	double V = 0;					// Volume

	Tetrahedron(std::vector<vec3> a_points, TetrahedralObject* a_obj)
	{
		myObj = a_obj;
		m_points = std::vector<vec3>(a_points);
		//ComputeCoefficients();
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
	void AddTet(std::vector<vec3> a_points);
	void DumpPoints();
	std::set<vec3> GetPointsSingleton();
	std::vector<Tetrahedron *> GetTets();
	vec3 GetMin();
	vec3 GetMax();
	void GenerateFragments(float cellSize);
	void Draw(GU_Detail *gdp);

private:
	std::set<vec3> m_pointSet; // we maintain a set so we make sure there are no duplicate points
							   // do we need this? would it be better to have duplicate points so things can be broken more easily?
	std::vector<Tetrahedron *> m_tets;
	std::vector<TetFragment*> m_frags;
	std::unique_ptr<MaterialData> m_matData;
	vec3 m_min, m_max; // maintain minimum and maximum ranges in the set
};