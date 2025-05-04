
#ifndef __SHATTERFX_EXPORT_PLUGIN_h__
#define __SHATTERFX_EXPORT_PLUGIN_h__

#include <SOP/SOP_Node.h>

#include "vec.h"
struct Tet
{
    std::vector<vec3> m_points;
    std::vector<GA_Offset> m_pointOffsets;

    Tet(std::vector<vec3> a_points)
    {
        m_points = std::vector<vec3>(a_points);
    }

    Tet(std::vector<vec3> a_points, std::vector<GA_Offset> a_pointOffsets)
    {
        m_points = std::vector<vec3>(a_points);
        m_pointOffsets = std::vector<GA_Offset>(a_pointOffsets);
    }

    vec3 GetCenterOfMass();
    //void Draw(GU_Detail* gdp, int fragmentId, GA_RWHandleS nameHandle);
};

struct TetFrag
{
    UT_StringHolder name;
    std::vector<Tet*> m_tets;

    TetFrag(UT_StringHolder n)
    {
        name = n;
        m_tets = std::vector<Tet*>();
    }

    ~TetFrag()
    {
        for (auto tet : m_tets)
        {
            delete tet;
        }
    }
};

struct VertInfo
{
    vec3 position;
    GA_Offset offset;
    UT_StringHolder fragment;

    VertInfo(vec3 p, GA_Offset o, UT_StringHolder f)
    {
        position = p;
        offset = o;
        fragment = f;
    }
};

namespace HDK_Sample {
    class SOP_Export : public SOP_Node {
    public:
        static OP_Node* myConstructor(OP_Network*, const char*,
            OP_Operator*);

        /// Stores the description of the interface of the SOP in Houdini.
        /// Each parm template refers to a parameter.
        static PRM_Template		 myTemplateList[];

        /// This optional data stores the list of local variables.
        static CH_LocalVariable	 myVariables[];

    protected:

        SOP_Export(OP_Network* net, const char* name, OP_Operator* op);
        virtual ~SOP_Export();

        /// Disable parameters according to other parameters.
        virtual unsigned		 disableParms();


        /// cookMySop does the actual work of the SOP computing, in this
        /// case, a LSYSTEM
        virtual OP_ERROR		 cookMySop(OP_Context& context);

        /// This function is used to lookup local variables that you have
        /// defined specific to your SOP.
        virtual bool		 evalVariableValue(
            fpreal& val,
            int index,
            int thread);
        // Add virtual overload that delegates to the super class to avoid
        // shadow warnings.
        virtual bool		 evalVariableValue(
            UT_String& v,
            int i,
            int thread)
        {
            return evalVariableValue(v, i, thread);
        }

    private:
        std::unordered_map<UT_StringHolder, TetFrag*> m_frags;
        vec3 m_center;
        std::unordered_map<vec3, VertInfo*> m_mapPoint2VertInfo;
        std::unordered_map<UT_StringHolder, vec3> m_mapFrag2Direction;

        void setUpTetFragsByName(const GU_Detail* gdp);
        void setUpTetFragsByNameSecondPass(const GU_Detail* gdp);
        void moveFragments();
        void moveFragmentsSecondPass();
        UT_Set<UT_StringHolder> getAllTetNames(const GU_Detail* gdp);
        void clearAll();
    };
} // End HDK_Sample namespace

#endif