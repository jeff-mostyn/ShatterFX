#include <set>
#include <vector>

#include <Eigen/Dense>
#include <SYS/SYS_Types.h>

#include "vec.h"

struct Tetrahedron {
	std::vector<vec3> m_points;
	std::vector<double> a, b, c, d;

	vec3 GetCenterOfMass();
	double TetrahedralVolume();
};

class TetrahedralObject {
public:
	TetrahedralObject();
	~TetrahedralObject();
	void AddTet(std::vector<vec3> a_points);
	void DumpPoints();
	std::set<vec3> GetPointsSingleton();
	std::vector<Tetrahedron*> GetTets();
	vec3 GetMin();
	vec3 GetMax();

private:
	std::set<vec3> m_pointSet;	// we maintain a set so we make sure there are no duplicate points
								// do we need this? would it be better to have duplicate points so things can be broken more easily?
	std::vector<Tetrahedron*> m_tets;
	vec3 m_min, m_max;			// maintain minimum and maximum ranges in the set
};