#include "loadGrammar.h"



void* LoadGrammarCmd::creator() {

	return new LoadGrammarCmd();
}


MStatus LoadGrammarCmd::doIt(const MArgList& args)  {
	return MS::kSuccess;
}