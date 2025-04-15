#ifndef GEN_GRAM_CMD_H
#define GEN_GRAM_CMD_H

#include <maya/MPxCommand.h>
#include <maya/MArgList.h>
#include "LoadMesh.h"

class GenerateGrammarCmd : public MPxCommand {
public:
    MStatus doIt(const MArgList& args) override;
    static void* creator();
};

#endif 
