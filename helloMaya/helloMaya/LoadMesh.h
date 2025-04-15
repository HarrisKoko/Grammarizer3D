#ifndef LOAD_MESH_CMD_H
#define LOAD_MESH_CMD_H

#include <maya/MPxCommand.h>
#include <maya/MArgList.h>
#include "graph.h"

class LoadMeshCmd : public MPxCommand {
public:
    static void* creator();
    MStatus doIt(const MArgList& args) override;

    static Graph& getLoadedGraph();
    static bool hasLoadedGraph(); 

private:
    static Graph s_loadedGraph;
    static bool s_hasLoadedGraph;
};

#endif // LOAD_MESH_CMD_H
