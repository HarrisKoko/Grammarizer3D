#ifndef GEN_MODELS_CMD_H
#define GEN_MODELS_CMD_H

#include <maya/MPxCommand.h>
#include <maya/MArgList.h>
#include "LoadMesh.h"
#include <maya/MFnMesh.h>
#include <maya/MPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MGlobal.h>
#include <maya/MDagModifier.h>
#include <maya/MSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MStatus.h>
#include <vector>


class GenerateModelsCmd : public MPxCommand {
public:
    MStatus doIt(const MArgList& args) override;
    static void* creator();
    MStatus createMesh(const std::vector<MPoint>& vertices,
        const std::vector<int>& faceCounts,
        const std::vector<int>& faceConnects);

};

#endif 
