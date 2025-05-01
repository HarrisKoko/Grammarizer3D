#ifndef LOAD_MESH_CMD_H
#define LOAD_MESH_CMD_H

#include <maya/MPxCommand.h>
#include <maya/MArgList.h>
#include "graph.h"

#include <maya/MFnMesh.h>
#include <maya/MItSelectionList.h>
#include <maya/MSelectionList.h>
#include <maya/MGlobal.h>
#include <maya/MPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MDagPath.h>
#include <iostream>
#include <maya/MObjectArray.h>          
#include <maya/MPlug.h>                 
#include <maya/MPlugArray.h>            
#include <maya/MFnDependencyNode.h> 
#include <maya/MObjectHandle.h>


class LoadMeshCmd : public MPxCommand {
public:
    static void* creator();
    MStatus doIt(const MArgList& args) override;

    static Graph& getLoadedGraph();
    static bool hasLoadedGraph(); 
    MObjectArray shaders;
    MIntArray shaderIndices;

private:
    static Graph s_loadedGraph;
    static bool s_hasLoadedGraph;
};

#endif // LOAD_MESH_CMD_H
