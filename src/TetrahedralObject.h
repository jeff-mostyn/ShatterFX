#include <set>
#include <vector>
#include <SYS/SYS_Types.h>

#include "vec.h"

struct Tetrahedron {
	std::vector<exint> pointIndices;
};

class TetrahedralObject {
public:
	TetrahedralObject();
	~TetrahedralObject();
	void AddPoint(vec3* a_point);
	void AddTet(std::vector<exint> a_indices);
	void DumpPoints();
	std::vector<vec3*> GetPoints();
	std::set<vec3> GetPointsSingleton();
	std::vector<Tetrahedron*> GetTets();

private:
	std::vector<vec3*> m_points;
	//std::vector<vec3*> m_pointsSingleton;
	std::set<vec3> m_pointSet;	// we maintain a set so we make sure there are no duplicate points
								// do we need this? would it be better to have duplicate points so things can be broken more easily?
	std::vector<Tetrahedron*> m_tets;
};