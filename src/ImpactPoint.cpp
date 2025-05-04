


//#include <UT/UT_DSOVersion.h>
//#include <RE/RE_EGLServer.h>

#include <CH/CH_LocalVariable.h>
#include <GA/GA_Primitive.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <GU/GU_PrimTetrahedron.h>
#include <OBJ/OBJ_Node.h>
#include <OP/OP_Director.h>
#include <OP/OP_Network.h>
#include <OP/OP_Node.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Math.h>

#include <limits.h>
#include "ImpactPoint.h"
#include "TetrahedralObject.h"
using namespace HDK_Sample;

// Help is stored in a "wiki" style text file. 
// See the sample_install.sh file for an example.

// DECLARE PARAMETERS (arg1: internal, arg2: descriptive)
static PRM_Name forceMagName("forceMag", "Impact Force Magnitude (N)");
static PRM_Name forceDirName("forceDir", "Impact Force Direction");
static PRM_Name forceLocName("forceLoc", "Impact Force Location");
static PRM_Name pickButtonName("getpos", "Calculate Fracture");
static PRM_Name filenameName("outputFile", "Output File (no suffix)");
static PRM_Name exportButtonName("export", "Export File");

// SET PARAMETER DEFAULTS

static PRM_Default forceMagDefault(100.0);
static PRM_Default forceDirDefault[]{
	PRM_Default(-1.0),
	PRM_Default(0.0),
	PRM_Default(0.0)
};
static PRM_Default forceLocDefault[]{
	PRM_Default(0.5),
	PRM_Default(0.0),
	PRM_Default(0.0)
};
static PRM_Default exportFileDefault(0, "out");

// DECLARE PARAMETER RANGES
static PRM_Range forceMagRange(PRM_RANGE_RESTRICTED, 1.0, PRM_RANGE_RESTRICTED, 1500.0);

// USE PARAM NAMES AND PARAMS TO INITIALIZE

PRM_Template SOP_Impact::myTemplateList[] = {
	PRM_Template(PRM_FLT, 1, &forceMagName, &forceMagDefault, nullptr, &forceMagRange),
	PRM_Template(
		PRM_FLT,		// Declare a 3-float parameter
		3,              // Number of components
		&forceDirName,    // Parameter name
		forceDirDefault),
	PRM_Template(
		PRM_FLT,		// Declare a 3-float parameter
		3,              // Number of components
		&forceLocName,    // Parameter name
		forceLocDefault),
	PRM_Template(PRM_CALLBACK, 1, &pickButtonName, 0, 0, 0, SOP_Impact::PickCallback),
	PRM_Template(
		PRM_STRING,       
		1,                
		&filenameName,         
		&exportFileDefault
	),
	PRM_Template(PRM_CALLBACK, 1, &exportButtonName, 0, 0, 0, SOP_Impact::ExportCallback),

	PRM_Template()
};


// Here's how we define local variables for the SOP.
enum {
	VAR_PT,		// Point number of the star
	VAR_NPT		// Number of points in the star
};

CH_LocalVariable
SOP_Impact::myVariables[] = {
	{ "PT",	VAR_PT, 0 },		// The table provides a mapping
	{ "NPT",	VAR_NPT, 0 },		// from text string to integer token
	{ 0, 0, 0 },
};

bool
SOP_Impact::evalVariableValue(fpreal& val, int index, int thread)
{
	// Not one of our variables, must delegate to the base class.
	return SOP_Node::evalVariableValue(val, index, thread);
}

OP_Node*
SOP_Impact::myConstructor(OP_Network* net, const char* name, OP_Operator* op) {
	return new SOP_Impact(net, name, op);
}

SOP_Impact::SOP_Impact(OP_Network* net, const char* name, OP_Operator* op)
	: SOP_Node(net, name, op) {
	// declare local var defaults?
}

SOP_Impact::~SOP_Impact() {}

unsigned
SOP_Impact::disableParms()
{
	return 0;
}

OP_ERROR
SOP_Impact::cookMySop(OP_Context& context)
{
	// ------------------------------------------------------------------
	// ------------------------ INPUT GUARDS ----------------------------
	// ------------------------------------------------------------------
	context = context;

	// make sure to lock input prior to work
	if (lockInput(0, context) >= UT_ERROR_ABORT) {
		addError(SOP_MESSAGE, "Failed to lock input geometry.");
		return error();
	}

	// collect input geometry as a const, because it should not be modified directly
	const GU_Detail* inputGeo = SOP_Impact::inputGeo(0, context);
	if (!inputGeo) {
		addError(SOP_MESSAGE, "No input geometry found.");
		return error();
	}

	// Ensure the input contains tetrahedral data
	bool hasPolygons = false;
	for (GA_Iterator it(inputGeo->getPrimitiveRange()); !it.atEnd(); ++it) {
		if (inputGeo->getPrimitive(*it)->getTypeId() == GEO_PRIMTETRAHEDRON) {
			hasPolygons = true;
			break;
		}
	}
	if (!hasPolygons) {
		addError(SOP_MESSAGE, "Input must be a tetrahedral mesh.");
		return error();
	}

	fpreal now = context.getTime();

	// PUT YOUR CODE HERE
	// DECLARE FIELDS FOR PARAMETERS      
	//    NOTE : [ALL-CAPS] is a function that you need to use and it is declared in the header file to update your values instantly while cooking 
	// NOW THAT WE HAVE HANDLES ON ALL VARIABLES, WE ALSO WILL DO CHECKS AND MAKE SURE THEY STAY IN SAFE RANGES
	
	float forceMag;
	forceMag = FORCE_MAG(now);

	vec3 forceDir;
	forceDir = FORCE_DIR(now);
	if (forceDir.Length() == 0) {
		forceDir = { 1.0, 1.0, 1.0 };
	}

	vec3 forceLoc;
	forceLoc = FORCE_LOC(now);


	// PROCESS DATA HERE
	// Normalize the direction vector
	UT_Vector3 normal(forceDir[0], forceDir[1], forceDir[2]);
	normal.normalize();

	UT_Vector3 center(forceLoc[0], forceLoc[1], forceLoc[2]);

	// DRAW GEOMETRY
	UT_Interrupt* boss;
	GU_PrimPoly* poly;
	GA_Offset p0, p1;

	// Check to see that there hasn't been a critical error in cooking the SOP.
	if (error() < UT_ERROR_ABORT) {
		boss = UTgetInterrupt();

		// Start the interrupt server
		if (boss->opStart("Building GEO")) {
			// ------------------------------------------------------------
			// ------------------------ OUTPUT ----------------------------
			// ------------------------------------------------------------

			// make a handle to an editable copy of the input geometry
			gdp->clearAndDestroy();
			duplicateSource(0, context, gdp);
		}

		// Compute two perpendicular vectors using Gram-Schmidt
		UT_Vector3 xdir, ydir;

		// Choose an arbitrary vector that's not colinear
		UT_Vector3 arbitrary(1, 0, 0);
		if (fabs(normal.dot(arbitrary)) > 0.99)  // Too aligned, choose another
			arbitrary = UT_Vector3(0, 1, 0);

		// First perpendicular vector
		xdir = arbitrary - (arbitrary.dot(normal)) * normal;
		xdir.normalize();

		// Second perpendicular vector
		ydir = normal;
		ydir.cross(xdir);
		ydir.normalize();
		
		// create cone group
		GA_PointGroup* coneGroup = gdp->newPointGroup("cone");
		
		//define cone parameters
		float height = 0.5f;
		float radius = 0.25f;
		int numSides = 16;

		// Compute base center
		UT_Vector3 baseCenter = center + -1 * normal * height;

		// draw cone
		for (int i = 0; i < numSides; ++i) {
			int next = (i + 1) % numSides;

			// Create quad face
			poly = GU_PrimPoly::build(gdp, 3, GU_POLY_CLOSED);

			// Add points for the quad
			GA_Offset p0 = poly->getPointOffset(0);
			coneGroup->addOffset(p0);
			GA_Offset p1 = poly->getPointOffset(1);
			coneGroup->addOffset(p1);
			GA_Offset p2 = poly->getPointOffset(2);
			coneGroup->addOffset(p2);

			float theta1 = 2 * M_PI * i / numSides;
			float theta2 = 2 * M_PI * next / numSides;
			UT_Vector3 dir1 = cos(theta1) * xdir + sin(theta1) * ydir;
			UT_Vector3 pos1 = baseCenter + radius * dir1;
			UT_Vector3 dir2 = cos(theta2) * xdir + sin(theta2) * ydir;
			UT_Vector3 pos2 = baseCenter + radius * dir2;

			// Set positions
			gdp->setPos3(p0, center);
			gdp->setPos3(p1, pos1);
			gdp->setPos3(p2, pos2);
		}

		// Tell the interrupt server that we've completed. Must do this
		// regardless of what opStart() returns.
		boss->opEnd();
	}

	unlockInput(0);

	return error();
}

int SOP_Impact::PickCallback(void* data, int index,
	float time, const PRM_Template*)
{
	// get SOP and context
	SOP_Impact* sop = static_cast<SOP_Impact*>(data);
	OP_Context myContext(time);

	// get inputs
	float forceMag;
	forceMag = sop->FORCE_MAG(time);

	vec3 forceDir;
	forceDir = sop->FORCE_DIR(time);
	if (forceDir.Length() == 0) {
		forceDir = { 1.0, 1.0, 1.0 };
	}
	forceDir.Normalize();

	vec3 forceLoc;
	forceLoc = sop->FORCE_LOC(time);

	// use this link for help setting up next node https://www.sidefx.com/docs/hdk/_h_d_k__node_intro__working_with_nodes.html
	
	// Variables
	OP_Network* parent = (OP_Network*)(sop->getParent()); // the network in which this node exists
	OP_Node* fileNode = parent->findNode("shatterExportFile1");
	OP_Node* fractureNode;
	OP_Node* fractureNodeExisting = parent->findNode("CVD1");
	OP_Node* input;

	// remove the CVD1 node if it already exists
	if (fractureNodeExisting) {
		parent->destroyNode(fractureNodeExisting);
	}

	// ---------------------------------------------------------
	//					Create & Attach CVD Node
	// ---------------------------------------------------------
	// create node
	fractureNode = parent->createNode("CusCentroidalVoronoi", "CVD1");
	
	// set parameters
	fractureNode->setFloat("forceMag", 0, time, forceMag);

	fractureNode->setFloat("forceDir", 0, time, forceDir[0]);
	fractureNode->setFloat("forceDir", 1, time, forceDir[1]);
	fractureNode->setFloat("forceDir", 2, time, forceDir[2]);

	fractureNode->setFloat("forceLoc", 0, time, forceLoc[0]);
	fractureNode->setFloat("forceLoc", 1, time, forceLoc[1]);
	fractureNode->setFloat("forceLoc", 2, time, forceLoc[2]);

	// connect the node
	input = parent->findNode("tetconform1");  // find /obj/geo1/tetConform1 as relative path
	if (input) {
		fractureNode->setInput(0, input);       // set first input to /obj/tetConform
	}

	// ---------------------------------------------------------
	//						Arrange Nodes
	// ---------------------------------------------------------
	// now that done we're done connecting it, position it relative to its inputs
	fractureNode->moveToGoodPosition();
	fractureNode->setRender(1);
	fractureNode->setDisplay(1);
	sop->forceRecook(); // Trigger cookMySop again

	if (fileNode) {
		fileNode->setInput(0, fractureNode);
		fileNode->moveToGoodPosition();
	}

	sop->setCurrent(0);
	fractureNode->setCurrent(1);

	return 1;
}

int SOP_Impact::ExportCallback(void* data, int index,
	float time, const PRM_Template*) {
	
	// get SOP and context
	SOP_Impact* sop = static_cast<SOP_Impact*>(data);
	OP_Context myContext(time);

	UT_String outputFile;
	sop->EXPORT_FILE(outputFile, time);

	OP_Network* parent = (OP_Network*)(sop->getParent());
	OP_Node* fractureNode = parent->findNode("CVD1");
	OP_Node* fileNode;
	OP_Node* fileNodeExisting = parent->findNode("shatterExportFile1");
	

	// ---------------------------------------------------------
	//			Create File Node or Assign Existing One
	// ---------------------------------------------------------
	if (!fileNodeExisting) {
		fileNode = parent->createNode("file", "shatterExportFile1");
	}
	else {
		fileNode = fileNodeExisting;
	}

	// if success, set node up and cook it
	if (fileNode) {
		// set parameters
		string fullPath = "$HIP/" + (std::string)outputFile + ".obj";

		PRM_Parm* filemodeParm = &fileNode->getParm("filemode");
		if (filemodeParm != nullptr) {
			fileNode->setInt("filemode", 0, time, 2);
		}

		fileNode->setString(fullPath, CH_STRING_LITERAL, "file", 0, time);

		if (fractureNode) {
			fileNode->setInput(0, fractureNode);
			fileNode->moveToGoodPosition();
			fileNode->forceRecook();
		}
	}

	// set as current node, but doesn't seem to run it
	fileNode->setCurrent(1);

	return 0;
}