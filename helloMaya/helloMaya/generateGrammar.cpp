#include "generateGrammar.h"

MStatus GenerateGrammarCmd::doIt(const MArgList& args) {
    if (!LoadMeshCmd::hasLoadedGraph()) {
        return MS::kFailure;
    }

    Graph& graph = LoadMeshCmd::getLoadedGraph(); // assumes copy constructor works. if issues, try making without.
    
    // Insert code here for grammar generation (saves to graph static variable)
    auto rules = graph.generateRules(8);
    Graph blank = Graph();
    auto opt = blank.applyRandomReplacementRule(rules,false);
    graph = opt.value_or(graph);

    return MS::kSuccess;
}

void* GenerateGrammarCmd::creator() {
    return new GenerateGrammarCmd();
}