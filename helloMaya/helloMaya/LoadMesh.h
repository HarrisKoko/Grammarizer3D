#ifndef LOAD_MESH_CMD_H
#define LOAD_MESH_CMD_H

#include <maya/MPxCommand.h>
#include <maya/MArgList.h>

class LoadMeshCmd : public MPxCommand {
public:
    MStatus doIt(const MArgList& args) override;
    static void* creator();
};

#endif // LOAD_MESH_CMD_H
