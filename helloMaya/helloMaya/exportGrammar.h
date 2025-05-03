#ifndef EXPORT_GRAM_CMD_H
#define EXPORT_GRAM_CMD_H
#include "jsonreader.h"
#include <maya/MPxCommand.h>
#include <maya/MArgList.h>
class ExportGrammarCmd : public MPxCommand
{
public:
	static void* creator();
	MStatus doIt(const MArgList& args) override;
private:

};

#endif 