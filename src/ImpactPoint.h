

#ifndef __SHATTERFX_IMPACT_PLUGIN_h__
#define __SHATTERFX_IMPACT_PLUGIN_h__

//#include <GEO/GEO_Point.h>
//

#include <SOP/SOP_Node.h>

#include "vec.h"

namespace HDK_Sample {
    class SOP_Impact : public SOP_Node {
    public:
        static OP_Node* myConstructor(OP_Network*, const char*,
            OP_Operator*);

        /// Stores the description of the interface of the SOP in Houdini.
        /// Each parm template refers to a parameter.
        static PRM_Template		 myTemplateList[];

        /// This optional data stores the list of local variables.
        static CH_LocalVariable	 myVariables[];

    protected:

        SOP_Impact(OP_Network* net, const char* name, OP_Operator* op);
        virtual ~SOP_Impact();

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
        /// The following list of accessors simplify evaluating the parameters
        /// of the SOP.

        // PUT YOUR CODE HERE
        // Here you need to declare functions which need to be called from the .C file to 
        // constantly update the cook function, these functions help you get the current value that the node has
        // Example : To declare a function to fetch angle you need to do it this way 
        // fpreal  ANGLE(fpreal t)     { return evalFloat("angle", 0, t); }

        fpreal FORCE_MAG(fpreal t) {
            return evalFloat("forceMag", 0, t);
        }

        vec3 FORCE_DIR(fpreal t) {
            float x = evalFloat("forceDir", 0, t);
            float y = evalFloat("forceDir", 1, t);
            float z = evalFloat("forceDir", 2, t);
            return vec3(x, y, z);
        }

        vec3 FORCE_LOC(fpreal t) {
            float x = evalFloat("forceLoc", 0, t);
            float y = evalFloat("forceLoc", 1, t);
            float z = evalFloat("forceLoc", 2, t);
            return vec3(x, y, z);
        }

        void EXPORT_FILE(UT_String& _str, fpreal t) {
            return evalString(_str, "outputFile", 0, t);
        }

        static int PickCallback(void* data, int index,
            float time, const PRM_Template*);

        static int ExportCallback(void* data, int index,
            float time, const PRM_Template*);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        /// Member variables are stored in the actual SOP, not with the geometry
        /// In this case these are just used to transfer data to the local 
        /// variable callback.
        /// Another use for local data is a cache to store expensive calculations.

        // NOTE : You can declare local variables here like this 
        bool hasSelectedPoint = false;
    };
} // End HDK_Sample namespace

#endif
