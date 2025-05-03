#include "generateGrammar.h"

std::vector<std::pair<Graph, Graph>> GenerateGrammarCmd::rules = std::vector<std::pair<Graph, Graph>>();

MSyntax GenerateGrammarCmd::newSyntax() {
    MSyntax syntax;
    syntax.addArg(MSyntax::kLong); 
    syntax.addArg(MSyntax::kLong); 
    return syntax;
}


MStatus GenerateGrammarCmd::doIt(const MArgList& args) {
    if (!LoadMeshCmd::hasLoadedGraph()) {
        return MS::kFailure;
    }
    int num_rules = 10;
    int num_iters = 1;
    if (args.length() >= 1) {
        MGlobal::displayInfo("Arg 0 (rules): " + MString() + args.asInt(0));
        num_rules = args.asInt(0);
    }
    if (args.length() >= 2) {
        MGlobal::displayInfo("Arg 1 (iters): " + MString() + args.asInt(1));
    }

    Graph& graph = LoadMeshCmd::getLoadedGraph(); 
    
    this->rules = graph.generateRules(num_rules);
    Graph blank = Graph();
    //for loop here using num_iters
    auto opt = blank.applyRandomReplacementRule(rules,false);
    graph = opt.value_or(graph);
    LoadMeshCmd::setLoadedGraph(graph);

    return MS::kSuccess;
}

void* GenerateGrammarCmd::creator() {
    return new GenerateGrammarCmd();
}