#ifndef LOAD_GRAM_CMD_H
#define LOAD_GRAM_CMD_H
#include "jsonreader.h"
#include <maya/MPxCommand.h>
#include <maya/MArgList.h>
class LoadGrammarCmd : public MPxCommand
{
public:
	static void* creator();
	MStatus doIt(const MArgList& args) override;
private:

};

#endif 