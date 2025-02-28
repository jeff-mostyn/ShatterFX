


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
			     0,				// Min # of sources
			     0,				// Max # of sources
			     SOP_CVD::myVariables,	// Local variables
			     OP_FLAG_GENERATOR)		// Flag it as generator
	    );
}


// DECLARE PARAMETERS (arg1: internal, arg2: descriptive)

static PRM_Name boundsName("bounds", "Voronoi Region Bounds");
static PRM_Name cellSizeName("cellSize", "Voronoi Cell Size");


// SET PARAMETER DEFAULTS

static PRM_Default boundsDefault[]{
	PRM_Default(0.0),
	PRM_Default(0.0),
	PRM_Default(0.0)
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
	fpreal now = context.getTime();

	// PUT YOUR CODE HERE
	// DECLARE FIELDS FOR PARAMETERS      
    //    NOTE : [ALL-CAPS] is a function that you need to use and it is declared in the header file to update your values instantly while cooking 
/*
	float angle;
	angle = ANGLE(now);

	float stepSize;
	stepSize = STEPSIZE(now);

	float iterations;
	iterations = ITERATIONS(now);

	UT_String grammar;
	GRAMMAR(grammar, now);*/
	float cellSize;
	cellSize = CELL_SIZE(now);

	vec3 bounds;
	bounds = BOUNDS(now);
	std::cout << "---------------------------------------------" << std::endl;

	std::cout << "bounds" << std::endl;
	std::cout << bounds[0] << ", " << bounds[1] << ", " << bounds[2] << std::endl;

	// PROCESS DATA HERE

	vec3 noiseVal = GeneratePerlinNoise(bounds);
	std::cout << "noise" << std::endl;
	std::cout << noiseVal[0] << ", " << noiseVal[1] << ", " << noiseVal[2] << std::endl;

	std::cout << "---------------------------------------------" << std::endl;
	// DRAW GEOMETRY
    float		 rad, tx, ty, tz;
    int			 divisions, plane;
    int			 xcoord =0, ycoord = 1, zcoord =2;
    float		 tmp;
    UT_Vector4		 pos;
    GU_PrimPoly		*poly;
    int			 i;
    UT_Interrupt	*boss;

    // Check to see that there hasn't been a critical error in cooking the SOP.
    if (error() < UT_ERROR_ABORT) {
		boss = UTgetInterrupt();
		gdp->clearAndDestroy();

		// Start the interrupt server
		if (boss->opStart("Building GEO")) {
      
		}

		// Tell the interrupt server that we've completed. Must do this
		// regardless of what opStart() returns.
		boss->opEnd();
    }

    return error();
}

void SOP_CVD::PopulateVoronoiPoints(vec3 a_bounds, float a_cellSize) {
	voronoiPoints.clear();

	int xMax = a_cellSize / a_bounds[0];
	int yMax = a_cellSize / a_bounds[1];
	int zMax = a_cellSize / a_bounds[2];

	for (float i = 0; i < xMax; i++) {
		for (float j = 0; j < yMax; j++) {
			for (float k = 0; k < zMax; k++) {
				// generate a [0, 1] noise value in the given cell
				// unfortunately not working and Idk why yet (just returns 0?)
				vec3 noiseVal = GeneratePerlinNoise(vec3(i, j, k));
			}
		}
	}
}

vec3 SOP_CVD::GeneratePerlinNoise(vec3 a_cell) {
	UT_Noise noise;
	noise.setType(UT_Noise::ALLIGATOR); // Set noise type to Perlin

	UT_Vector3 noiseVec;
	UT_Vector3 pos;
	pos.x() = a_cell[0];
	pos.y() = a_cell[1];
	pos.z() = a_cell[2];

	std::cout << "input" << std::endl;
	std::cout << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;

	noise.turbulence(pos, 1.0, noiseVec, 1.0, 1.0);

	std::cout << "output" << std::endl;
	std::cout << noiseVec.x() << ", " << noiseVec.y() << ", " << noiseVec.z() << std::endl;

	return vec3(noiseVec.x(), noiseVec.y(), noiseVec.z());
}

