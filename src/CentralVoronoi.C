


#include <UT/UT_DSOVersion.h>
//#include <RE/RE_EGLServer.h>


#include <UT/UT_Math.h>
#include <UT/UT_Interrupt.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <CH/CH_LocalVariable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>


#include <limits.h>
#include "CentralVoronoi.h"
using namespace HDK_Sample;

// Help is stored in a "wiki" style text file. 
// See the sample_install.sh file for an example.

///
/// newSopOperator is the hook that Houdini grabs from this dll
/// and invokes to register the SOP.  In this case we add ourselves
/// to the specified operator table.
///
void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(
	    new OP_Operator("CusCentroidalVoronoi",			// Internal name
			    "CVD",						// UI name
			     SOP_CVD::myConstructor,	// How to build the SOP
			     SOP_CVD::myTemplateList,	// My parameters
			     1,				// Min # of sources
			     1,				// Max # of sources
			     SOP_CVD::myVariables,	// Local variables
			     OP_FLAG_GENERATOR)		// Flag it as generator
	    );
}


// DECLARE PARAMETERS (arg1: internal, arg2: descriptive)

static PRM_Name boundsName("bounds", "Voronoi Region Bounds");
static PRM_Name cellSizeName("cellSize", "Voronoi Cell Size");


// SET PARAMETER DEFAULTS

static PRM_Default boundsDefault[]{
	PRM_Default(1.0),
	PRM_Default(1.0),
	PRM_Default(1.0)
};
static PRM_Default cellSizeDefault(1.0);

// USE PARAM NAMES AND PARAMS TO INITIALIZE

PRM_Template SOP_CVD::myTemplateList[] = {
	PRM_Template(
		PRM_FLT,		// Declare a 3-float parameter
		3,              // Number of components
		&boundsName,    // Parameter name
		boundsDefault),
	PRM_Template(PRM_FLT, 1, &cellSizeName, &cellSizeDefault),

    PRM_Template()
};


// Here's how we define local variables for the SOP.
enum {
	VAR_PT,		// Point number of the star
	VAR_NPT		// Number of points in the star
};

CH_LocalVariable
SOP_CVD::myVariables[] = {
    { "PT",	VAR_PT, 0 },		// The table provides a mapping
    { "NPT",	VAR_NPT, 0 },		// from text string to integer token
    { 0, 0, 0 },
};

bool
SOP_CVD::evalVariableValue(fpreal &val, int index, int thread)
{
    // myCurrPoint will be negative when we're not cooking so only try to
    // handle the local variables when we have a valid myCurrPoint index.
    /*
	if (myCurrPoint >= 0) {
	// Note that "gdp" may be null here, so we do the safe thing
	// and cache values we are interested in.
		switch (index)
		{
			case VAR_PT:
				val = (fpreal) myCurrPoint;
				return true;
			case VAR_NPT:
				val = (fpreal) myTotalPoints;
				return true;
			default:
			// do nothing ;
		}
    }
	*/
    // Not one of our variables, must delegate to the base class.
    return SOP_Node::evalVariableValue(val, index, thread);
}

OP_Node *
SOP_CVD::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {
    return new SOP_CVD(net, name, op);
}

SOP_CVD::SOP_CVD(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op) {
	// declare local var defaults?
}

SOP_CVD::~SOP_CVD() {}

unsigned
SOP_CVD::disableParms()
{
    return 0;
}

OP_ERROR
SOP_CVD::cookMySop(OP_Context &context)
{
	// ------------------------------------------------------------------
	// ------------------------ INPUT GUARDS ----------------------------
	// ------------------------------------------------------------------
	
	// make sure to lock input prior to work
	if (lockInput(0, context) >= UT_ERROR_ABORT) {
		addError(SOP_MESSAGE, "Failed to lock input geometry.");
		return error();
	}

	// collect input geometry as a const, because it should not be modified directly
	const GU_Detail* inputGeo = SOP_CVD::inputGeo(0, context);
	if (!inputGeo) {
		addError(SOP_MESSAGE, "No input geometry found.");
		return error();
	}

	// Ensure the input contains polygonal data
	bool hasPolygons = false;
	for (GA_Iterator it(inputGeo->getPrimitiveRange()); !it.atEnd(); ++it) {
		if (inputGeo->getPrimitive(*it)->getTypeId() == GEO_PRIMPOLY) {
			hasPolygons = true;
			break;
		}
	}
	if (!hasPolygons) {
		addError(SOP_MESSAGE, "Input must be a polygonal mesh.");
		return error();
	}

	fpreal now = context.getTime();

	// PUT YOUR CODE HERE
	// DECLARE FIELDS FOR PARAMETERS      
    //    NOTE : [ALL-CAPS] is a function that you need to use and it is declared in the header file to update your values instantly while cooking 
	vec3 bounds;
	bounds = BOUNDS(now);

	float cellSize;
	cellSize = CELL_SIZE(now);

	std::cout << "---------------------------------------------" << std::endl;

	// PROCESS DATA HERE

	PopulateVoronoiPoints(bounds, cellSize);
	PrintVoronoiPoints();


	std::cout << "---------------------------------------------" << std::endl;
	
	
	// DRAW GEOMETRY
    /* These were all examples, not sure exactly what they doooo*/
		/*float		 rad, tx, ty, tz;
		int			 divisions, plane;
		int			 xcoord =0, ycoord = 1, zcoord =2;
		float		 tmp;
		UT_Vector4		 pos;
		GU_PrimPoly		*poly;
		int			 i;*/
	UT_Interrupt	*boss;

    // Check to see that there hasn't been a critical error in cooking the SOP.
    if (error() < UT_ERROR_ABORT) {
		boss = UTgetInterrupt();

		// Start the interrupt server
		if (boss->opStart("Building GEO")) {
			// ----------------------------------------------------------------
			// ------------------------ PROCESSING ----------------------------
			// ----------------------------------------------------------------

			// make a handle to an editable copy of the input geometry
			gdp->clearAndDestroy();
			duplicateSource(0, context, gdp);
		}

		//draw points
		/*
		for (float i = 0; i < voronoiPoints.size(); i++) {
			for (float j = 0; j < voronoiPoints[i].size(); j++) {
				for (float k = 0; k < voronoiPoints[i][j].size(); k++) {

					//poly = GU_PrimPoly::build(gdp, 2, GU_POLY_OPEN);
					GA_Offset ptoff = gdp->appendPoint();
					vec3 point = voronoiPoints[i][j][k];
					//DrawVoronoiEdge(gdp, vec3(0, 0, 0), point);
					gdp->setPos3(ptoff, UT_Vector3(point[0], point[1], point[2]));
					//GA_Offset p0 = poly->getPointOffset(0);
					//GA_Offset p1 = poly->getPointOffset(1);
					
					//UT_Vector3 start(0, 0, 0);
					//UT_Vector3 end(point[0], point[1], point[2]);
					//gdp->setPos3(p0, start);
					//gdp->setPos3(p1, end);
				}
			}
		}*/

		//DrawVoronoiCells(gdp);

		// Tell the interrupt server that we've completed. Must do this
		// regardless of what opStart() returns.
		boss->opEnd();
    }

	unlockInput(0);

    return error();
}

void SOP_CVD::PopulateVoronoiPoints(vec3 a_bounds, float a_cellSize) {
	voronoiPoints.clear();

	int xMax = a_bounds[0] / a_cellSize;
	int yMax = a_bounds[1] / a_cellSize;
	int zMax = a_bounds[2] / a_cellSize;

	// first we compute points in a range 2 larger in each axis. We will use these to average our approximate centroids.
	std::vector<std::vector < std::vector<vec3>>> tmpVoronoiPoints;
	for (float i = 0; i <= xMax+1; i++) {
		tmpVoronoiPoints.push_back(std::vector<std::vector<vec3>>());
		for (float j = 0; j <= yMax+1; j++) {
			tmpVoronoiPoints[i].push_back(std::vector<vec3>());
			for (float k = 0; k <= zMax+1; k++) {
				// generate a [0, 1] noise value in the given cell
				// then scale by cell size to get offset from corner
				vec3 noiseVal = Noise3D(vec3(i, j, k)) * a_cellSize;

				vec3 node = vec3(a_cellSize * (i - 1.0), a_cellSize * (j - 1.0), a_cellSize * (k - 1.0)) + noiseVal;

				tmpVoronoiPoints[i][j].push_back(node);
			}
		}
	}

	// from the above temp matrix, we will use the average of the current cell and the six in cardinal directions to compute
	// the approximate centroid. This is in lieu of more expensive methods like Delaunay triangulation or Qhull
	for (float i = 1; i <= xMax; i++) {
		voronoiPoints.push_back(std::vector<std::vector<vec3>>());
		// using i less 1 because we're starting at one higher index to account for larger temp matrix
		for (float j = 1; j <= yMax; j++) {
			voronoiPoints[i - 1].push_back(std::vector<vec3>());
			for (float k = 1; k <= zMax; k++) {
				vec3 average = (
					tmpVoronoiPoints[i][j][k] +
					tmpVoronoiPoints[i+1][j][k] +
					tmpVoronoiPoints[i-1][j][k] +
					tmpVoronoiPoints[i][j+1][k] +
					tmpVoronoiPoints[i][j-1][k] +
					tmpVoronoiPoints[i][j][k+1]	+ 
					tmpVoronoiPoints[i][j][k-1] 
				) / 7.0;

				// using i and j less 1 because we're starting at one higher index to account for larger temp matrix
				voronoiPoints[i - 1][j - 1].push_back(average);
			}
		}
	}
}

// This is not range [0-1], gotta work on that
vec3 SOP_CVD::GeneratePerlinNoise(vec3 a_cell) {
	UT_Noise noise;
	noise.setType(UT_Noise::ALLIGATOR); // Set noise type to Perlin
	noise.setSeed(1234);

	UT_Vector3 noiseVec;
	UT_Vector3 pos(a_cell[0], a_cell[1], a_cell[2]);

	noise.turbulence(pos, 1.0, noiseVec, 0.7, 1.0);

	std::cout << "Turbulence at position " << pos << ": " << noiseVec << std::endl;

	return vec3(noiseVec.x(), noiseVec.y(), noiseVec.z());
}

// This is modeled after some glsl shader noise functions
vec3 SOP_CVD::Noise3D(vec3 a_cell) {
	vec3 nonFractional = vec3(
		abs(sin(Dot(a_cell, vec3(343.7, 151.1, 934.2)))),
		abs(sin(Dot(a_cell, vec3(678.9, 432.1, 342.8)))),
		abs(sin(Dot(a_cell, vec3(845.3, 473.2, 432.7))))
	) * 200419.35;

	vec3 fractional = vec3(
		nonFractional[0] - (long)nonFractional[0],
		nonFractional[1] - (long)nonFractional[1],
		nonFractional[2] - (long)nonFractional[2]
	);

	return fractional;
}

// prints voronoi points in "slices"
const void SOP_CVD::PrintVoronoiPoints() {
	for (float i = 0; i < voronoiPoints.size(); i++) {
		for (float j = 0; j < voronoiPoints[i].size(); j++) {
			for (float k = 0; k < voronoiPoints[i][j].size(); k++) {

				std::cout << voronoiPoints[i][j][k] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "-----------" << std::endl;
	}
}

void SOP_CVD::DrawVoronoiEdge(GU_Detail* gdp, vec3 p1, vec3 p2)
{
	GU_PrimPoly* poly = GU_PrimPoly::build(gdp, 2, GU_POLY_OPEN);
	GA_Offset ptoff1 = poly->getPointOffset(0);
	GA_Offset ptoff2 = poly->getPointOffset(1);

	gdp->setPos3(ptoff1, UT_Vector3(p1[0], p1[1], p1[2]));
	gdp->setPos3(ptoff2, UT_Vector3(p2[0], p2[1], p2[2]));

}

void SOP_CVD::DrawVoronoiCells(GU_Detail* gdp) {
	// For a 3D Voronoi diagram, we need to compute the Delaunay triangulation first
	// Then construct cells from the dual of that triangulation

	// This is a simplified approach for demonstration
	for (int i = 0; i < voronoiPoints.size(); i++) {
		for (int j = 0; j < voronoiPoints[i].size(); j++) {
			for (int k = 0; k < voronoiPoints[i][j].size(); k++) {
				vec3 p1 = voronoiPoints[i][j][k];

				// Check all potential neighboring cells in a reasonable radius
				for (int ni = std::max(0, i - 1); ni <= std::min((int)voronoiPoints.size() - 1, i + 1); ni++) {
					for (int nj = std::max(0, j - 1); nj <= std::min((int)voronoiPoints[i].size() - 1, j + 1); nj++) {
						for (int nk = std::max(0, k - 1); nk <= std::min((int)voronoiPoints[i][j].size() - 1, k + 1); nk++) {
							// Skip self
							if (i == ni && j == nj && k == nk)
								continue;

							vec3 p2 = voronoiPoints[ni][nj][nk];

							// Find midpoint between the two points
							vec3 midpoint = (p1 + p2) * 0.5;

							// Direction from p1 to p2
							vec3 dir = p2 - p1;
							float dist = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
							dir = dir / dist;  // Normalize

							// For simplicity in this example, we'll just draw short line segments 
							// at the midpoint perpendicular to the line connecting the two points
							// In a full implementation, you'd compute intersections of these planes

							// Perpendicular directions (simplified for visualization)
							vec3 perp1, perp2;

							// Find perpendicular vectors
							if (fabs(dir[0]) > fabs(dir[1])) {
								perp1 = vec3(-dir[2], 0, dir[0]);
							}
							else {
								perp1 = vec3(0, -dir[2], dir[1]);
							}

							float len = sqrt(perp1[0] * perp1[0] + perp1[1] * perp1[1] + perp1[2] * perp1[2]);
							perp1 = perp1 / len;  // Normalize

							// Cross product for third perpendicular vector
							perp2[0] = dir[1] * perp1[2] - dir[2] * perp1[1];
							perp2[1] = dir[2] * perp1[0] - dir[0] * perp1[2];
							perp2[2] = dir[0] * perp1[1] - dir[1] * perp1[0];

							// Draw a small "cross" at the midpoint to visualize the cell boundary
							float edgeLength = dist * 0.25;  // Adjust size as needed

							vec3 edge1Start = midpoint - perp1 * edgeLength;
							vec3 edge1End = midpoint + perp1 * edgeLength;
							DrawVoronoiEdge(gdp, edge1Start, edge1End);

							vec3 edge2Start = midpoint - perp2 * edgeLength;
							vec3 edge2End = midpoint + perp2 * edgeLength;
							DrawVoronoiEdge(gdp, edge2Start, edge2End);
						}
					}
				}
			}
		}
	}
}