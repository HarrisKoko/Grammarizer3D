#include "generateGrammar.h"

MStatus GenerateGrammarCmd::doIt(const MArgList& args) {
    if (!LoadMeshCmd::hasLoadedGraph()) {
        return MS::kFailure;
    }

    Graph& graph = LoadMeshCmd::getLoadedGraph(); // assumes copy constructor works. if issues, try making without.

    




    return MS::kSuccess;
}

void* GenerateGrammarCmd::creator() {
    return new GenerateGrammarCmd();
}