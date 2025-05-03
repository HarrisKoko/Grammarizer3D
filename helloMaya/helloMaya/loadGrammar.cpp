#include "loadGrammar.h"
#include "LoadMesh.h"


void* LoadGrammarCmd::creator() {

	return new LoadGrammarCmd();
}


MStatus LoadGrammarCmd::doIt(const MArgList& args)  {
	// args 0 should have file path from mel plugin dialog.
	if (args.length() < 1) {
		return MS::kFailure;
	}
	Graph& graph = LoadMeshCmd::getLoadedGraph();
	MString filepath = args.asString(0);
	// Go from Mstring to Qstring for jsonreader loadgrammar from file function.
	QString qfilepath = QString::fromUtf8(filepath.asChar());
	// Use grammar from return of function and apply to graph. then overwrite graph just like generateGrammar.
	std::vector<std::pair<Graph, Graph>> rules = JSONReader::LoadGrammarFromFile(qfilepath);

	Graph blank = Graph();
	auto opt = blank.applyRandomReplacementRule(rules, false);
	graph = opt.value_or(graph);
	LoadMeshCmd::setLoadedGraph(graph);


	return MS::kSuccess;
}