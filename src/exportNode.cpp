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

#include "exportNode.h"
using namespace HDK_Sample;

// -----------------------------------------------------
// 
// -------------------- TETRAHEDRON --------------------
// 
// -----------------------------------------------------
vec3 Tet::GetCenterOfMass()
{
	if (m_points.size() == 4)
	{
		return 0.25 * (m_points[0] + m_points[1] + m_points[2] + m_points[3]);
	}
	else
	{
		return vec3Zero;
	}
}

// Help is stored in a "wiki" style text file. 
// See the sample_install.sh file for an example.

// DECLARE PARAMETERS (arg1: internal, arg2: descriptive)

// SET PARAMETER DEFAULTS

// DECLARE PARAMETER RANGES

// USE PARAM NAMES AND PARAMS TO INITIALIZE

PRM_Template SOP_Export::myTemplateList[] = {
	PRM_Template()
};


// Here's how we define local variables for the SOP.
enum {
	VAR_PT,		// Point number of the star
	VAR_NPT		// Number of points in the star
};

CH_LocalVariable
SOP_Export::myVariables[] = {
	{ "PT",	VAR_PT, 0 },		// The table provides a mapping
	{ "NPT",	VAR_NPT, 0 },		// from text string to integer token
	{ 0, 0, 0 },
};

bool
SOP_Export::evalVariableValue(fpreal& val, int index, int thread)
{
	// Not one of our variables, must delegate to the base class.
	return SOP_Node::evalVariableValue(val, index, thread);
}

OP_Node*
SOP_Export::myConstructor(OP_Network* net, const char* name, OP_Operator* op) {
	return new SOP_Export(net, name, op);
}

SOP_Export::SOP_Export(OP_Network* net, const char* name, OP_Operator* op)
	: SOP_Node(net, name, op) {
	// declare local var defaults?
}

SOP_Export::~SOP_Export() {}

unsigned
SOP_Export::disableParms()
{
	return 0;
}

OP_ERROR
SOP_Export::cookMySop(OP_Context& context)
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
	const GU_Detail* inputGeo = SOP_Export::inputGeo(0, context);
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

		setUpTetFragsByName(gdp);

		moveFragments();

		gdp->consolidatePoints(0.001f);

		clearAll();

		setUpTetFragsByNameSecondPass(gdp);

		moveFragmentsSecondPass();

		forceRecook();
		// Tell the interrupt server that we've completed. Must do this
		// regardless of what opStart() returns.
		boss->opEnd();
	}

	unlockInput(0);

	return error();
}

void SOP_Export::setUpTetFragsByName(const GU_Detail* gdp)
{
	vec3 min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	vec3 max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);

	GA_ROHandleS nameHandle(gdp->findPrimitiveAttribute("name"));
	if (!nameHandle.isValid())
	{
		std::cout << "No 'name' attribute found on tetrahedra." << std::endl;
		return;
	}

	for (GA_Iterator it(gdp->getPrimitiveRange()); !it.atEnd(); ++it)
	{
		// Filter for only tetrahedra (GA_PRIMPOLY with 4 points, or better: specific tetrahedral primitive type)
		const GA_Primitive* prim = gdp->getPrimitive(*it);
		if (!prim || prim->getTypeId() != GEO_PRIMTETRAHEDRON)
			continue;
		
		// Gets 4 vertices of Tetrehedron
		std::vector<vec3> vertices = std::vector<vec3>();
		std::vector<GA_Offset> offsets = std::vector<GA_Offset>();

		for (int i = 0; i < 4; ++i)
		{
			//GA_Offset vertOffset = prim->getVertexOffset(i);

			GA_Offset pointOffset = prim->getPointOffset(i);
			UT_Vector3 pos = gdp->getPos3(pointOffset);

			vec3 point = vec3(pos.x(), pos.y(), pos.z());

			// update min max
			min[0] = point[0] < min[0] ? point[0] : min[0];
			min[1] = point[1] < min[1] ? point[1] : min[1];
			min[2] = point[2] < min[2] ? point[2] : min[2];
			max[0] = point[0] > max[0] ? point[0] : max[0];
			max[1] = point[1] > max[1] ? point[1] : max[1];
			max[2] = point[2] > max[2] ? point[2] : max[2];

			vertices.push_back(point);
			offsets.push_back(pointOffset);
		}

		// Create Tetrahedron Object
		Tet* tet = new Tet(vertices, offsets);

		// Get Tet Fragment name
		UT_StringHolder name = UT_StringHolder(nameHandle.get(*it));

		// Check if Fragment Exists
		if (!m_frags.count(name))
		{
			// If not create Fragment
			TetFrag* frag = new TetFrag(name);

			// Add frag to m_frags map
			m_frags[name] = frag;
		}

		// Add Tet to fragment
		m_frags[name]->m_tets.push_back(tet);
	}

	m_center = 0.5f * (min + max);
	return;
}

void SOP_Export::moveFragments()
{
	// Step 1: Loop through fragments
	for (const auto& frag : m_frags)
	{
		// Compute average center of all tets in the fragment
		vec3 fragCenter = vec3Zero;
		int count = 0;
		for (Tet* tet : frag.second->m_tets)
		{
			fragCenter += tet->GetCenterOfMass();
			count++;
		}
		if (count == 0) continue;
		fragCenter /= static_cast<float>(count);

		// Step 3: Compute direction from object center to fragment center
		vec3 direction = fragCenter - m_center;
		float length = direction.Length();
		if (length == 0.0f) continue; // Prevent divide by zero
		direction /= length; // Normalize

		// Step 4: Store direction in map
		m_mapFrag2Direction[frag.first] = direction;

		// Step 4: Move each point in the fragment
		for (Tet* tet : frag.second->m_tets)
		{
			for (int i = 0; i < 4; i++)
			{
				tet->m_points[i] += direction * 1.0f;
				gdp->setPos3(tet->m_pointOffsets[i], 
					UT_Vector3(tet->m_points[i][0], tet->m_points[i][1], tet->m_points[i][2]));
			}
		}
	}
}

UT_Set<UT_StringHolder> SOP_Export::getAllTetNames(const GU_Detail* gdp)
{
	UT_Set<UT_StringHolder> nameSet;

	GA_ROHandleS nameHandle(gdp->findPrimitiveAttribute("name"));
	if (!nameHandle.isValid())
	{
		std::cout << "No 'name' attribute found on tetrahedra." << std::endl;
		return nameSet;
	}

	for (GA_Iterator it(gdp->getPrimitiveRange()); !it.atEnd(); ++it)
	{
		// Filter for only tetrahedra (GA_PRIMPOLY with 4 points, or better: specific tetrahedral primitive type)
		const GA_Primitive* prim = gdp->getPrimitive(*it);
		if (!prim || prim->getTypeId() != GEO_PRIMTETRAHEDRON)
			continue;

		nameSet.insert(UT_StringHolder(nameHandle.get(*it)));
	}

	return nameSet;
}

void SOP_Export::clearAll()
{
	for (auto frag : m_frags)
	{
		delete frag.second;
	}

	m_frags.clear();

	for (auto vert : m_mapPoint2VertInfo)
	{
		delete vert.second;
	}

	m_mapPoint2VertInfo.clear();
}

void SOP_Export::setUpTetFragsByNameSecondPass(const GU_Detail* gdp)
{
	vec3 min = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	vec3 max = vec3(FLT_MIN, FLT_MIN, FLT_MIN);

	GA_ROHandleS nameHandle(gdp->findPrimitiveAttribute("name"));
	if (!nameHandle.isValid())
	{
		std::cout << "No 'name' attribute found on tetrahedra." << std::endl;
		return;
	}

	for (GA_Iterator it(gdp->getPrimitiveRange()); !it.atEnd(); ++it)
	{
		// Filter for only tetrahedra (GA_PRIMPOLY with 4 points, or better: specific tetrahedral primitive type)
		const GA_Primitive* prim = gdp->getPrimitive(*it);
		if (!prim || prim->getTypeId() != GEO_PRIMTETRAHEDRON)
			continue;

		// Gets 4 vertices of Tetrehedron
		std::vector<vec3> vertices = std::vector<vec3>();

		// Get Tet Fragment name
		UT_StringHolder name = UT_StringHolder(nameHandle.get(*it));

		for (int i = 0; i < 4; ++i)
		{
			//GA_Offset vertOffset = prim->getVertexOffset(i);

			GA_Offset pointOffset = prim->getPointOffset(i);
			UT_Vector3 pos = gdp->getPos3(pointOffset);

			vec3 point = vec3(pos.x(), pos.y(), pos.z());

			// update min max
			min[0] = point[0] < min[0] ? point[0] : min[0];
			min[1] = point[1] < min[1] ? point[1] : min[1];
			min[2] = point[2] < min[2] ? point[2] : min[2];
			max[0] = point[0] > max[0] ? point[0] : max[0];
			max[1] = point[1] > max[1] ? point[1] : max[1];
			max[2] = point[2] > max[2] ? point[2] : max[2];

			vertices.push_back(point);

			if (!m_mapPoint2VertInfo.count(point)) {
				VertInfo* vert = new VertInfo(point, pointOffset, name);
				m_mapPoint2VertInfo[point] = vert;
			}
		}

		// Create Tetrahedron Object
		Tet* tet = new Tet(vertices);

		// Check if Fragment Exists
		if (!m_frags.count(name))
		{
			// If not create Fragment
			TetFrag* frag = new TetFrag(name);

			// Add frag to m_frags map
			m_frags[name] = frag;
		}

		// Add Tet to fragment
		m_frags[name]->m_tets.push_back(tet);
	}

	m_center = 0.5f * (min + max);
	return;
}

void SOP_Export::moveFragmentsSecondPass()
{
	// Loop through unique points
	for (const auto& point : m_mapPoint2VertInfo)
	{
		vec3 newPos = point.second->position + (m_mapFrag2Direction[point.second->fragment] * -1.0f);
		gdp->setPos3(point.second->offset,
			UT_Vector3(newPos[0], newPos[1], newPos[2]));
	}
}