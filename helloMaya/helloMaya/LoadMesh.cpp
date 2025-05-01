#include "LoadMesh.h"
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


Graph LoadMeshCmd::s_loadedGraph  = Graph();
bool LoadMeshCmd::s_hasLoadedGraph = false;

Graph& LoadMeshCmd::getLoadedGraph() {
    return s_loadedGraph;
}

bool LoadMeshCmd::hasLoadedGraph() {
    return s_hasLoadedGraph;
}

MStatus LoadMeshCmd::doIt(const MArgList& args) {
    // Get the selected mesh in Maya
    MSelectionList selection;
    MGlobal::getActiveSelectionList(selection);

    // Makes iterable list excluding anything not a kMesh.
    MItSelectionList iter(selection, MFn::kMesh);

    // If no object selected, show error.
    if (iter.isDone()) {
        MGlobal::displayError("No mesh selected.");
        return MS::kFailure;
    }

    // loop through meshes
    for (; !iter.isDone(); iter.next()) {
        // Get dagPath (required for making MFnMesh) for current mesh
        MDagPath dagPath;
        iter.getDagPath(dagPath);

        // Get MFnMesh from dagPath
        MFnMesh meshFn(dagPath);
        MGlobal::displayInfo("Mesh loaded: " + meshFn.name());

        // Get Shaders
        // Get connected shaders
        MObjectArray shaders;
        MIntArray shaderIndices; // Maps each face to a shader index
        meshFn.getConnectedShaders(dagPath.instanceNumber(), shaders, shaderIndices);

        // Print material info per shading group
        for (unsigned int i = 0; i < shaders.length(); ++i) {
            MObject shadingGroup = shaders[i];

            // Find the surface shader connected to the shading group
            MPlug surfaceShaderPlug = MFnDependencyNode(shadingGroup).findPlug("surfaceShader", true);
            MPlugArray connections;
            surfaceShaderPlug.connectedTo(connections, true, false);

            if (connections.length() > 0) {
                MObject shaderNode = connections[0].node();
                MFnDependencyNode shaderFn(shaderNode);

                MGlobal::displayInfo("Material " + shaderFn.name() + " is connected to shading group " + MFnDependencyNode(shadingGroup).name());

                // Example: try to get the color plug
                MPlug colorPlug = shaderFn.findPlug("color", true);
                if (!colorPlug.isNull()) {
                    float r, g, b;
                    colorPlug.child(0).getValue(r);
                    colorPlug.child(1).getValue(g);
                    colorPlug.child(2).getValue(b);
                    char buf[128];
                    snprintf(buf, sizeof(buf), "  Color: (%.2f, %.2f, %.2f)", r, g, b);
                    MGlobal::displayInfo(buf);
                }
            }
        }



        // Get vertices 
        MPointArray verts;
        meshFn.getPoints(verts, MSpace::kWorld);

        // Get face information
        MIntArray faceCounts, faceConnects;
        meshFn.getVertices(faceCounts, faceConnects);

        // Print vertex data
        for (unsigned int i = 0; i < verts.length(); ++i) {
            MString msg = MString("") + i + ": (" +
                verts[i].x + ", " +
                verts[i].y + ", " +
                verts[i].z + ")";
            MGlobal::displayInfo(msg);
        }

        // Print face data
        MGlobal::displayInfo("Faces:");
        int idx = 0;
        for (unsigned int i = 0; i < faceCounts.length(); ++i) {
            char buffer[128];
            snprintf(buffer, sizeof(buffer), "Face %u:", i);  
            MString faceMsg(buffer);
            for (int j = 0; j < faceCounts[i]; ++j) {
                faceMsg += " ";
                faceMsg += faceConnects[idx++];
            }
            MGlobal::displayInfo(faceMsg);
        }


        std::cout << "Loaded mesh with " << verts.length() << " vertices and "
            << faceCounts.length() << " faces.\n";

        // Convert vertex data to glm::vec3
        std::vector<glm::vec3> vertexPositions;
        for (unsigned int i = 0; i < verts.length(); ++i) {
            vertexPositions.emplace_back(verts[i].x, verts[i].y, verts[i].z);
        }

        // Convert face data to vector of {vertex indices, face index} pairs
        std::vector<std::pair<std::vector<int>, int>> faceData;
        idx = 0;
        for (unsigned int i = 0; i < faceCounts.length(); ++i) {
            std::vector<int> faceVerts;
            for (int j = 0; j < faceCounts[i]; ++j) {
                faceVerts.push_back(faceConnects[idx++]);
            }
            faceData.emplace_back(faceVerts, static_cast<int>(i));
        }

        // Construct the graph
        s_loadedGraph = Graph(vertexPositions, faceData);
        s_hasLoadedGraph = true;


    }

    return MS::kSuccess;
}

void* LoadMeshCmd::creator() {
    return new LoadMeshCmd();
}
