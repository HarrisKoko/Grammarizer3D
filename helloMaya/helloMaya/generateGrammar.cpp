#include "generateGrammar.h"

std::vector<std::pair<Graph, Graph>> GenerateGrammarCmd::rules = std::vector<std::pair<Graph, Graph>>();

MStatus GenerateGrammarCmd::doIt(const MArgList& args) {
    if (!LoadMeshCmd::hasLoadedGraph()) {
        return MS::kFailure;
    }

    Graph& graph = LoadMeshCmd::getLoadedGraph(); // assumes copy constructor works. if issues, try making without.
    
    // Insert code here for grammar generation (saves to graph static variable)
    this->rules = graph.generateRules(8);
    Graph blank = Graph();
    auto opt = blank.applyRandomReplacementRule(rules,false);
    graph = opt.value_or(graph);
    LoadMeshCmd::setLoadedGraph(graph);
    return MS::kSuccess;
}

void* GenerateGrammarCmd::creator() {
    return new GenerateGrammarCmd();
}