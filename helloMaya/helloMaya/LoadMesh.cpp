#include "LoadMesh.h"
#include <maya/MFnMesh.h>
#include <maya/MItSelectionList.h>
#include <maya/MSelectionList.h>
#include <maya/MGlobal.h>
#include <maya/MPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MDagPath.h>
#include <iostream>
#include "graph.h"

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
        Graph graph(vertexPositions, faceData);

    }

    return MS::kSuccess;
}

void* LoadMeshCmd::creator() {
    return new LoadMeshCmd();
}
