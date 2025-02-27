


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

//
// Help is stored in a "wiki" style text file. 
//
// See the sample_install.sh file for an example.
//
// NOTE : Follow this tutorial if you have any problems setting up your visual studio 2008 for Houdini 
//  http://www.apileofgrains.nl/setting-up-the-hdk-for-houdini-12-with-visual-studio-2008/


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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//PUT YOUR CODE HERE
//You need to declare your parameters here
//Example to declare a variable for angle you can do like this :
//static PRM_Name		angleName("angle", "Angle");

static PRM_Name angleName("angle", "Angle");
static PRM_Name stepSizeName("stepSize", "Step Size");
static PRM_Name iterationsName("iterations", "Iterations");
static PRM_Name grammarName("grammar", "Grammar");







//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				     ^^^^^^^^    ^^^^^^^^^^^^^^^
//				     internal    descriptive version


// PUT YOUR CODE HERE
// You need to setup the initial/default values for your parameters here
// For example : If you are declaring the inital value for the angle parameter
// static PRM_Default angleDefault(30.0);	

static PRM_Default angleDefault(30.0);
static PRM_Default stepSizeDefault(1.0);
static PRM_Default iterationsDefault(2.0);
static PRM_Default grammarDefault(0, "");









////////////////////////////////////////////////////////////////////////////////////////

PRM_Template
SOP_CVD::myTemplateList[] = {
// PUT YOUR CODE HERE
// You now need to fill this template with your parameter name and their default value
// EXAMPLE : For the angle parameter this is how you should add into the template
// PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &angleName, &angleDefault, 0),
// Similarly add all the other parameters in the template format here
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &angleName, &angleDefault, 0),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &stepSizeName, &stepSizeDefault, 0),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &iterationsName, &iterationsDefault, 0),
	PRM_Template(PRM_STRING, PRM_Template::PRM_EXPORT_MIN, 1, &grammarName, &grammarDefault, 0),





/////////////////////////////////////////////////////////////////////////////////////////////

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
			/* do nothing */;
		}
    }
    // Not one of our variables, must delegate to the base class.
    return SOP_Node::evalVariableValue(val, index, thread);
}

OP_Node *
SOP_CVD::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {
    return new SOP_CVD(net, name, op);
}

SOP_CVD::SOP_CVD(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op) {
    myCurrPoint = -1;	// To prevent garbage values from being returned
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
	// Decare the necessary variables and get always keep getting the current value in the node
	// For example to always get the current angle thats set in the node ,you need to :
	//    float angle;
	//    angle = ANGLE(now)       
    //    NOTE : ANGLE is a function that you need to use and it is declared in the header file to update your values instantly while cooking 
	LSystem myplant;

	float angle;
	angle = ANGLE(now);

	float stepSize;
	stepSize = STEPSIZE(now);

	float iterations;
	iterations = ITERATIONS(now);

	UT_String grammar;
	GRAMMAR(grammar, now);


	///////////////////////////////////////////////////////////////////////////

	//PUT YOUR CODE HERE
	// Next you need to call your Lystem cpp functions 
	// Below is an example , you need to call the same functions based on the variables you declare
    // myplant.loadProgramFromString("F\nF->F[+F]F[-F]";  
    // myplant.setDefaultAngle(30.0f);
    // myplant.setDefaultStep(1.0f);

	myplant.loadProgramFromString((std::string)grammar);
	myplant.setDefaultAngle(angle);
	myplant.setDefaultStep(stepSize);

	///////////////////////////////////////////////////////////////////////////////

	// PUT YOUR CODE HERE
	// You the need call the below function for all the genrations ,so that the end points points will be
	// stored in the branches vector , you need to declare them first

	std::vector<LSystem::Branch> branches;

	//for (int i = 0; i < generations ; i++)
	for (int i = 0; i < iterations ; i++) {
		  myplant.process(i, branches);
	}

	///////////////////////////////////////////////////////////////////////////////////


	// Now that you have all the branches ,which is the start and end point of each point ,its time to render 
	// these branches into Houdini 
    

	// PUT YOUR CODE HERE
	// Declare all the necessary variables for drawing cylinders for each branch 
    float		 rad, tx, ty, tz;
    int			 divisions, plane;
    int			 xcoord =0, ycoord = 1, zcoord =2;
    float		 tmp;
    UT_Vector4		 pos;
    GU_PrimPoly		*poly;
    int			 i;
    UT_Interrupt	*boss;

    // Since we don't have inputs, we don't need to lock them.

    divisions  = 5;	// We need twice our divisions of points
    myTotalPoints = divisions;		// Set the NPT local variable value
    myCurrPoint   = 0;			// Initialize the PT local variable



    // Check to see that there hasn't been a critical error in cooking the SOP.
    if (error() < UT_ERROR_ABORT)
    {
	boss = UTgetInterrupt();
	if (divisions < 4)
	{
	    // With the range restriction we have on the divisions, this
	    //	is actually impossible, but it shows how to add an error
	    //	message or warning to the SOP.
	    addWarning(SOP_MESSAGE, "Invalid divisions");
	    divisions = 4;
	}
	gdp->clearAndDestroy();

	// Start the interrupt server
	if (boss->opStart("Building LSYSTEM"))
	{
        // PUT YOUR CODE HERE
	    // Build a polygon
	    // You need to build your cylinders inside Houdini from here
		// TIPS:
		// Use GU_PrimPoly poly = GU_PrimPoly::build(see what values it can take)
		// Also use GA_Offset ptoff = poly->getPointOffset()
		// and gdp->setPos3(ptoff,YOUR_POSITION_VECTOR) to build geometry.

		float radius = 0.1;
		int sides = 6;
		float tinc = M_PI * 2 / (float)sides;

		// these will hold points of start and end of cylinders
		std::vector<UT_Vector3> startpts, endpts;

		for (int i = 0; i < branches.size(); i++) {
			// clear old cylinder points
			startpts.clear();
			endpts.clear();

			// calculate tangent and bitangent so branches are oriented the right way / maintain thickness
			UT_Vector3 dir; 
			dir(0) = (branches[i].second - branches[i].first)[0];
			dir(1) = (branches[i].second - branches[i].first)[1];
			dir(2) = (branches[i].second - branches[i].first)[2];

			UT_Vector3 tan = (fabs(dir[0]) > 0.99) ? UT_Vector3(0, 1, 0) : UT_Vector3(1, 0, 0);
			tan -= dir * (dir.dot(tan));
			tan.normalize();

			UT_Vector3 bit = dir;
			bit.cross(tan);
			bit.normalize();

			// create circle for start of branch
			GU_PrimPoly* start = GU_PrimPoly::build(gdp, sides, GU_POLY_CLOSED);
			for (int j = 0; j < sides; j++) {
				// Since we expect the local variables to be used in specifying
				// the radii, we have to evaluate the channels INSIDE the loop
				// through the points...
				float tmp = (float)j * tinc;

				UT_Vector3 center = UT_Vector3(branches[i].first[0], branches[i].first[1], branches[i].first[2]);
				UT_Vector3 posStart = center 
					+ (tan * radius * SYScos(tmp)) 
					+ (bit * radius * SYSsin(tmp));
				startpts.push_back(posStart);

				GA_Offset ptoff = gdp->appendPointOffset(); // Create a new point in the geometry detail and get its offset
				start->setVertexPoint(j, ptoff); // Associate the vertex with this point
				gdp->setPos3(ptoff, posStart);   // Set the position of the new point
			}

			// create circle for end of branch
			GU_PrimPoly* end = GU_PrimPoly::build(gdp, sides, GU_POLY_CLOSED);
			for (int j = 0; j < sides; j++) {
				// Since we expect the local variables to be used in specifying
				// the radii, we have to evaluate the channels INSIDE the loop
				// through the points...
				float tmp = (float)j * tinc;

				UT_Vector3 center = UT_Vector3(branches[i].second[0], branches[i].second[1], branches[i].second[2]);
				UT_Vector3 posEnd = center
					+ (tan * radius * SYScos(tmp))
					+ (bit * radius * SYSsin(tmp));	
				endpts.push_back(posEnd);

				GA_Offset ptoff = gdp->appendPointOffset(); // Create a new point in the geometry detail and get its offset
				start->setVertexPoint(j, ptoff); // Associate the vertex with this point
				gdp->setPos3(ptoff, posEnd);   // Set the position of the new point
			}

			// create sides of branch
			for (int k = 0; k < sides; k++) {
				int next = k + 1 < sides ? k + 1 : 0;
				GU_PrimPoly* side = GU_PrimPoly::build(gdp, 4, GU_POLY_CLOSED);

				GA_Offset ptoff0 = gdp->appendPointOffset();
				side->setVertexPoint(0, ptoff0);
				gdp->setPos3(ptoff0, startpts[k]);

				GA_Offset ptoff1 = gdp->appendPointOffset();
				side->setVertexPoint(1, ptoff1);
				gdp->setPos3(ptoff1, startpts[next]);
				
				GA_Offset ptoff2 = gdp->appendPointOffset();
				side->setVertexPoint(2, ptoff2);
				gdp->setPos3(ptoff2, endpts[next]);

				GA_Offset ptoff3 = gdp->appendPointOffset();
				side->setVertexPoint(3, ptoff3);
				gdp->setPos3(ptoff3, endpts[k]);
			}
		}


		////////////////////////////////////////////////////////////////////////////////////////////

	    // Highlight the star which we have just generated.  This routine
	    // call clears any currently highlighted geometry, and then it
	    // highlights every primitive for this SOP. 
	    select(GU_SPrimitive);
	}

	// Tell the interrupt server that we've completed. Must do this
	// regardless of what opStart() returns.
	boss->opEnd();
    }

    myCurrPoint = -1;
    return error();
}

