


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
#include <GU/GU_PrimTetrahedron.h>
#include <GA/GA_Primitive.h>


#include <limits.h>
#include "CentralVoronoi.h"
#include "ImpactPoint.h"
#include "TetrahedralObject.h"
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
	table->addOperator(
		new OP_Operator("ImpactPoint",			// Internal name
			"Impact",						// UI name
			SOP_Impact::myConstructor,	// How to build the SOP
			SOP_Impact::myTemplateList,	// My parameters
			1,				// Min # of sources
			1,				// Max # of sources
			SOP_Impact::myVariables,	// Local variables
			OP_FLAG_GENERATOR)		// Flag it as generator
	);
}


// DECLARE PARAMETERS (arg1: internal, arg2: descriptive)

static PRM_Name stiffnessName("stiffness", "Young's Modulus (GPa)");
static PRM_Name strainRatioName("strainRatio", "Poisson Ratio");
static PRM_Name energySpreadName("energySpread", "Energy Spread Factor");
static PRM_Name energyPercentName("energyConsumptionPercent", "Energy Consumption Pct");
static PRM_Name fractureToughnessName("fractureToughness", "Fracture Toughness (J/m^2)");
static PRM_Name forceMagName("forceMag", "Impact Force Magnitude (N)");
static PRM_Name forceDirName("forceDir", "Impact Force Direction");
static PRM_Name forceLocName("forceLoc", "Impact Force Location");
static PRM_Name filenameName("outputFile", "Output File (no suffix)");
static PRM_Name exportButtonName("export", "Export File");
static PRM_Name physicsSimButtonName("physics", "Simulate Physics");
static PRM_Name explodedViewButtonName("explode", "Explode View");

// SET PARAMETER DEFAULTS

static PRM_Default stiffnessDefault(40.);
static PRM_Default strainRatioDefault(0.18);
static PRM_Default fractureToughnessDefault(15.0);
static PRM_Default energySpreadDefault(0.3);
static PRM_Default energyPercentDefault(0.15);
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
static PRM_Range youngsModulusRange(PRM_RANGE_RESTRICTED, 0.01, PRM_RANGE_RESTRICTED, 1200.0);
static PRM_Range poissonRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_RESTRICTED, 0.4999);
static PRM_Range fractureToughnessRange(PRM_RANGE_RESTRICTED, 10.0, PRM_RANGE_RESTRICTED, 500);
static PRM_Range energySpreadRange(PRM_RANGE_RESTRICTED, 0.01, PRM_RANGE_RESTRICTED, 2.0);
static PRM_Range energyPercentRange(PRM_RANGE_RESTRICTED, 0.1, PRM_RANGE_RESTRICTED, 1.0);
static PRM_Range forceMagRange(PRM_RANGE_RESTRICTED, 1.0, PRM_RANGE_RESTRICTED, 1500.0);

// USE PARAM NAMES AND PARAMS TO INITIALIZE

static PRM_Name materialPropsHeaderName("header_matProps", "Material Properties");
static PRM_Name forceParamsHeaderName("header_forceParams", "Force Paramters");

PRM_Template SOP_CVD::myTemplateList[] = {
	PRM_Template(PRM_LABEL, 1, &materialPropsHeaderName),
	PRM_Template(
		PRM_FLT, 
		1, 
		&stiffnessName, 
		&stiffnessDefault, 
		nullptr, 
		&youngsModulusRange,
		nullptr, nullptr, 0,
		"Ability of the material to withstand deformation under stress. In other words, it is the relationship between applied pressure and distortion.\nHigher values indicate stiffness, and lower values pliability."
	),
	PRM_Template(
		PRM_FLT, 
		1, 
		&strainRatioName, 
		&strainRatioDefault, 
		nullptr, 
		&poissonRange,
		nullptr, nullptr, 0,
		"Describes the ability of the material to expand in directions perpendicular to direction of pressure.\nHigher values resist distortion, while lower lower values deform easily."
	),
	PRM_Template(
		PRM_FLT, 
		1, 
		&fractureToughnessName, 
		&fractureToughnessDefault, 
		nullptr, 
		&fractureToughnessRange,
		nullptr, nullptr, 0,
		"Ability of a material to resist crack propagation. When exceeded, crack propagation becomes rapid and unlimited."
	),
	PRM_Template(
		PRM_FLT, 
		1, 
		&energySpreadName, 
		&energySpreadDefault, 
		nullptr, 
		&energySpreadRange,
		nullptr, nullptr, 0,
		"Controls how far from a facture seed point energy is 'collected' to determine how much energy that fracture will release.\nHigher values generally result in fewer fractures, as energy is consumed faster."
	),
	PRM_Template(
		PRM_FLT, 
		1, 
		&energyPercentName, 
		&energyPercentDefault, 
		nullptr, 
		&energyPercentRange,
		nullptr, nullptr, 0,
		"Controls cap of how much of total internal energy is released via fractures. Higher values correspond to more fractures."
	),
	PRM_Template(PRM_FLT, 1, &forceMagName, &forceMagDefault, nullptr, &forceMagRange),
		
	PRM_Template(PRM_LABEL, 1, &forceParamsHeaderName),
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
	PRM_Template(
		PRM_STRING,
		1,
		&filenameName,
		&exportFileDefault
	),
	PRM_Template(PRM_CALLBACK, 1, &exportButtonName, 0, 0, 0, SOP_CVD::ExportCallback),
	PRM_Template(PRM_CALLBACK, 1, &physicsSimButtonName, 0, 0, 0, SOP_CVD::PhysicsSimCallback),
	PRM_Template(PRM_CALLBACK, 1, &explodedViewButtonName, 0, 0, 0, SOP_CVD::ExplodedViewCallback),

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
	float stiffness;
	stiffness = STIFFNESS(now) * pow(10, 9); // we are inputting Young's Modulus in GPa, need to convert to Pa

	float strainRatio;
	strainRatio = STRAIN_RATIO(now);

	float fractureToughness;
	fractureToughness = FRACTURE_TOUGHNESS(now);

	float energySpread;
	energySpread = ENERGY_SPREAD(now);

	float energyPercent;
	energyPercent = ENERGY_PERCENT(now);

	float forceMag;
	forceMag = FORCE_MAG(now);

	vec3 forceDir;
	forceDir = FORCE_DIR(now);
	if (forceDir.Length() == 0) {
		forceDir = { 1.0, 1.0, 1.0 };
	}

	vec3 forceLoc;
	forceLoc = FORCE_LOC(now);
	
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
			// ------------------------------------------------------------
			// ------------------------ OUTPUT ----------------------------
			// ------------------------------------------------------------

			// make a handle to an editable copy of the input geometry
			gdp->clearAndDestroy();
			duplicateSource(0, context, gdp);
		}

		// ----------------------------------------------------------------------------------
		// ------------------------ PROCESS OBJECT IN GDP BUFFER ----------------------------
		// ----------------------------------------------------------------------------------

		// DYNAMIC MEMORY - DO NOT LEAVE FUNCTION WITHOUT DELETING
		std::unique_ptr<MaterialData> matData = std::make_unique<MaterialData>();
		matData->stiffness = stiffness;
		matData->strainRatio = strainRatio;
		matData->fractureToughness = fractureToughness;

		TetrahedralObject* obj = new TetrahedralObject(std::move(matData), energyPercent, energySpread);
		obj->DumpPoints();

		GA_Iterator primIter(gdp->getPrimitiveRange());
		for (; !primIter.atEnd(); ++primIter)
		{
			GA_Offset primOffset = primIter.getOffset();
			const GA_Primitive* prim = gdp->getPrimitive(primOffset);

			// Check if the primitive is a tetrahedron
			if (prim->getTypeId() == GEO_PRIMTETRAHEDRON)
			{
				const GEO_PrimTetrahedron* tetra = static_cast<const GEO_PrimTetrahedron*>(prim);
				std::vector<vec3> vertices = std::vector<vec3>();

				// Access the four vertices of the tetrahedron
				for (int i = 0; i < 4; ++i)
				{
					GA_Offset vertOffset = tetra->getVertexOffset(i);

					GA_Offset pointOffset = tetra->getPointOffset(i);
					UT_Vector3 pos = gdp->getPos3(pointOffset);

					vec3 point = vec3(pos.x(), pos.y(), pos.z());
					vertices.push_back(point);
				}

				// add the vertices for this tetrahedron to the object
				obj->AddTet(vertices);

			}
		}

		// flush gdp to render our stuff
		gdp->clearAndDestroy();

		std::cout << "-----------" << std::endl;

		std::cout << "set size: " << obj->GetPointsSingleton().size() << std::endl;

		std::cout << "-----------" << std::endl;

		// Compute all information to generate voronoi points and whatnot
		obj->ComputeMaterialInformation();
		obj->RegisterImpact(forceDir, forceMag, forceLoc);

		vector<vec3> fractureSites = obj->GenerateFractureSites(forceLoc);

		if (fractureSites.size() > 0) {
			// try to draw generated points
			for (float i = 0; i < fractureSites.size(); i++) {
				GA_Offset ptoff = gdp->appendPoint();
				vec3 point = fractureSites[i];
				gdp->setPos3(ptoff, UT_Vector3(point[0], point[1], point[2]));
			}

			obj->GenerateFragments(fractureSites);
			obj->Draw(gdp);
			gdp->consolidatePoints(0.001f);
		}

		delete obj;


		// Tell the interrupt server that we've completed. Must do this
		// regardless of what opStart() returns.
		boss->opEnd();
    }

	unlockInput(0);

    return error();
}

int SOP_CVD::ExportCallback(void* data, int index,
	float time, const PRM_Template*) {
	// Get SOP and context
	SOP_CVD* sop = static_cast<SOP_CVD*>(data);
	OP_Context myContext(time);

	// Set up Input for Output File Name
	UT_String outputFile;
	sop->EXPORT_FILE(outputFile, time);

	// Get Parent Network
	OP_Network* parent = (OP_Network*)(sop->getParent());
	

	// ---------------------------------------------------------
	//			Create File Node or Assign Existing One
	// ---------------------------------------------------------
	OP_Node* fileNode;
	OP_Node* fileNodeExisting = parent->findNode("shatterExportFile1");

	// Create or Assign Existing Node
	if (!fileNodeExisting) {
		fileNode = parent->createNode("file", "shatterExportFile1");
	}
	else {
		fileNode = fileNodeExisting;
	}

	// if success, set node up and cook it
	if (fileNode) {
		// Set up parameters

		// Full path to output file
		string fullPath = "$HIP/" + (std::string)outputFile + ".obj";

		// Set File Mode to Write Files
		PRM_Parm* filemodeParm = &fileNode->getParm("filemode");
		if (filemodeParm != nullptr) {
			fileNode->setInt("filemode", 0, time, 2);
		}

		// Set full path for output file location + name
		fileNode->setString(fullPath, CH_STRING_LITERAL, "file", 0, time);

		// Set input as CVD Fracture Geometry
		if (sop) {
			fileNode->setInput(0, sop);
			fileNode->moveToGoodPosition();
			fileNode->forceRecook();
		}
	}

	// Set as current node and run it
	sop->setCurrent(0);
	fileNode->setCurrent(1);
	fileNode->setRender(1);
	fileNode->setDisplay(1);

	return 0;
}

int SOP_CVD::PhysicsSimCallback(void* data, int index,
	float time, const PRM_Template*)
{
	// Get SOP and context
	SOP_Impact* sop = static_cast<SOP_Impact*>(data);
	OP_Context myContext(time);

	// Get Parent network and CVD fracture node
	OP_Network* parent = (OP_Network*)(sop->getParent());
	OP_Node* fractureNode = parent->findNode("CVD1");
	

	// ---------------------------------------------------------
	//			Create Fracture Group Node or Assign Existing One
	// ---------------------------------------------------------
	// Create or Assign Existing Node
	OP_Node* groupNode;
	OP_Node* groupNodeExisting = parent->findNode("fractureGroup");

	if (!groupNodeExisting) {
		groupNode = parent->createNode("groupcreate", "fractureGroup");
	}
	else {
		groupNode = groupNodeExisting;
	}

	// if success, set node up and cook it
	if (groupNode) {
		// Set up parameters

		// Set group type to points
		PRM_Parm* groupParm = &groupNode->getParm("grouptype");
		if (groupParm != nullptr) {
			groupNode->setInt("grouptype", 0, time, 1);
		}

		// Set input to CVD fracture node
		if (fractureNode) {
			groupNode->setInput(0, fractureNode);
			groupNode->moveToGoodPosition();
			groupNode->forceRecook();
		}
	}


	// ---------------------------------------------------------
	//			Create Rest Position Node or Assign Existing One
	// ---------------------------------------------------------
	OP_Node* restNode;
	OP_Node* restNodeExisting = parent->findNode("rest1");

	// Create or Assign Existing Node
	if (!restNodeExisting) {
		restNode = parent->createNode("rest", "rest1");
	}
	else {
		restNode = restNodeExisting;
	}

	// if success, set node up and cook it
	if (restNode) {
		// Set input to group node
		if (groupNode) {
			restNode->setInput(0, groupNode);
			restNode->moveToGoodPosition();
			restNode->forceRecook();
		}
	}

	
	// ---------------------------------------------------------
	//			Create Assemble Node or Assign Existing One
	// ---------------------------------------------------------
	OP_Node* assembleNode;
	OP_Node* assembleNodeExisting = parent->findNode("assemble1");

	// Create or Assign Existing Node
	if (!assembleNodeExisting) {
		assembleNode = parent->createNode("assemble", "assemble1");
	}
	else {
		assembleNode = assembleNodeExisting;
	}

	// if success, set node up and cook it
	if (assembleNode) {
		// Set up parameters

		// Set Create Name Attribute to False
		PRM_Parm* assembleParm1 = &assembleNode->getParm("newname");
		if (assembleParm1 != nullptr) {
			assembleNode->setInt("newname", 0, time, 0);
		}

		// Set Create Packed Primitives to True
		PRM_Parm* assembleParm2 = &assembleNode->getParm("pack_geo");
		if (assembleParm2 != nullptr) {
			assembleNode->setInt("pack_geo", 0, time, 1);
		}

		// Set input to Rest Node
		if (restNode) {
			assembleNode->setInput(0, restNode);
			assembleNode->moveToGoodPosition();
			assembleNode->forceRecook();
		}
	}


	// ---------------------------------------------------------
	//			Get Cone Collider from Impact Node
	// ---------------------------------------------------------
	OP_Node* coneGroupNode;
	OP_Node* coneGroupNodeExisting = parent->findNode("coneGroup");
	OP_Node* impactNode = parent->findNode("ImpactPoint1"); // is there a better way? since this is user generated

	// Create or Assign Existing Node
	if (!coneGroupNodeExisting) {
		coneGroupNode = parent->createNode("groupcreate", "coneGroup");
	}
	else {
		coneGroupNode = coneGroupNodeExisting;
	}

	// if success, set node up and cook it
	if (coneGroupNode) {
		// Set up parameters

		// Set group type to points
		PRM_Parm* groupParm1 = &coneGroupNode->getParm("grouptype");
		if (groupParm1 != nullptr) {
			coneGroupNode->setInt("grouptype", 0, time, 1);
		}

		// Set group to cone
		string cone = "cone";
		coneGroupNode->setString(cone, CH_STRING_LITERAL, "basegroup", 0, time);

		// Set input to impact Node if it exists
		if (impactNode) {
			coneGroupNode->setInput(0, impactNode);
			coneGroupNode->moveToGoodPosition();
			coneGroupNode->forceRecook();
		}
	}


	// ---------------------------------------------------------
	//			Create Blast Node or Assign Existing One
	// ---------------------------------------------------------
	OP_Node* blastNode;
	OP_Node* blastNodeExisting = parent->findNode("coneBlast");

	// Create or Assign Existing Node
	if (!blastNodeExisting) {
		blastNode = parent->createNode("blast", "coneBlast");
	}
	else {
		blastNode = blastNodeExisting;
	}

	// if success, set node up and cook it
	if (blastNode) {
		// Set up parameters

		// set group to cone group
		string cone = "cone";
		blastNode->setString(cone, CH_STRING_LITERAL, "group", 0, time);

		// Negate Blast node to keep only cone group points
		PRM_Parm* blastParm = &blastNode->getParm("negate");
		if (blastParm != nullptr) {
			blastNode->setInt("negate", 0, time, 1);
		}

		// set blast node input to group node
		if (coneGroupNode) {
			blastNode->setInput(0, coneGroupNode);
			blastNode->moveToGoodPosition();
			blastNode->forceRecook();
		}
	}


	// ---------------------------------------------------------
	//			Create Transform Node or Assign Existing One
	// ---------------------------------------------------------
	OP_Node* transformNode;
	OP_Node* transformNodeExisting = parent->findNode("coneTransform");

	// Create or Assign Existing Node
	if (!transformNodeExisting) {
		transformNode = parent->createNode("xform", "coneTransform");
	}
	else {
		transformNode = transformNodeExisting;
	}

	// if success, set node up and cook it
	if (transformNode) {
		// Set up parameters
		
		// Set transform node input to blast node
		if (blastNode) {
			transformNode->setInput(0, blastNode);
			transformNode->moveToGoodPosition();
			transformNode->forceRecook();
		}
	}

	// ---------------------------------------------------------
	//			Create RBD Solver Node or Assign Existing One
	// ---------------------------------------------------------
	OP_Node* rbdNode;
	OP_Node* rbdNodeExisting = parent->findNode("rbdbulletsolver1");

	// Create or Assign Existing Node
	if (!rbdNodeExisting) {
		rbdNode = parent->createNode("rbdbulletsolver", "rbdbulletsolver1");
	}
	else {
		rbdNode = rbdNodeExisting;
	}

	// if success, set node up and cook it
	if (rbdNode) {
		// Set up parameters

		// Set Collision Type to Deforming to use keyframe animation on transform node
		PRM_Parm* rbdParam1 = &rbdNode->getParm("collision_initialstate");
		if (rbdParam1 != nullptr) {
			rbdNode->setInt("collision_initialstate", 0, time, 2);
		}

		// Add ground plane
		PRM_Parm* rbdParam2 = &rbdNode->getParm("useground");
		if (rbdParam2 != nullptr) {
			rbdNode->setInt("useground", 0, time, 1);
		}

		// Move ground plane down - right now hard coded might want to calculate a good value based on size of object
		PRM_Parm* rbdParam3 = &rbdNode->getParm("ground_posy");
		if (rbdParam3 != nullptr) {
			rbdNode->setInt("ground_posy", 0, time, -2);
		}

		// Set input 0 to assemble node to get fractured geometry
		if (assembleNode) {
			rbdNode->setInput(0, assembleNode);
		}

		// Set input 3 to transform node to get cone collision geometry and any animation key framed
		if (transformNode) {
			rbdNode->setInput(3, transformNode);
		}

		rbdNode->moveToGoodPosition();
		rbdNode->forceRecook();
	}

	// set as current node and run it
	sop->setCurrent(0);
	rbdNode->setCurrent(1);
	rbdNode->setRender(1);
	rbdNode->setDisplay(1);

	return 0;
}

int SOP_CVD::ExplodedViewCallback(void* data, int index,
	float time, const PRM_Template*) {
	// Get SOP and context
	SOP_CVD* sop = static_cast<SOP_CVD*>(data);
	OP_Context myContext(time);

	// Get Parent Network
	OP_Network* parent = (OP_Network*)(sop->getParent());


	// ---------------------------------------------------------
	//			Create Exploded View Node or Assign Existing One
	// ---------------------------------------------------------
	OP_Node* explodedViewNode;
	OP_Node* explodedViewNodeExisting = parent->findNode("explodedview1");

	// Create or Assign Existing Node
	if (!explodedViewNodeExisting) {
		explodedViewNode = parent->createNode("explodedview::2.0", "explodedview1");
	}
	else {
		explodedViewNode = explodedViewNodeExisting;
	}

	// if success, set node up and cook it
	if (explodedViewNode) {
		// Set input as CVD Fracture Geometry
		if (sop) {
			explodedViewNode->setInput(0, sop);
			explodedViewNode->moveToGoodPosition();
			explodedViewNode->forceRecook();
		}
	}

	// Set as current node and run it
	sop->setCurrent(0);
	explodedViewNode->setCurrent(1);
	explodedViewNode->setRender(1);
	explodedViewNode->setDisplay(1);

	return 0;
}