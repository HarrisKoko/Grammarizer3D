#include "LoadMesh.h"


// Static graph object that stores the loaded graph
Graph LoadMeshCmd::s_loadedGraph = Graph();
// Boolean flag to track if the graph has been loaded
bool LoadMeshCmd::s_hasLoadedGraph = false;

MObjectArray LoadMeshCmd::shaders = MObjectArray();

// Returns the loaded graph
Graph& LoadMeshCmd::getLoadedGraph() {
    return s_loadedGraph;
}

// Returns whether a graph has been loaded
bool LoadMeshCmd::hasLoadedGraph() {
    return s_hasLoadedGraph;
}

MObjectArray& LoadMeshCmd::getShaders() {
    return shaders;
}


MStatus LoadMeshCmd::doIt(const MArgList& args) {
    // Get the active selection list in Maya
    MSelectionList selection;
    MGlobal::getActiveSelectionList(selection);

    // Set up the iterator to go through selected mesh objects
    MItSelectionList iter(selection, MFn::kMesh);

    // If no object is selected, display an error
    if (iter.isDone()) {
        MGlobal::displayError("No mesh selected.");
        return MS::kFailure;
    }

    // Loop through each selected mesh
    for (; !iter.isDone(); iter.next()) {
        // Get the dagPath (required for creating MFnMesh) for the current mesh
        MDagPath dagPath;
        iter.getDagPath(dagPath);

        // Create MFnMesh from the dagPath
        MFnMesh meshFn(dagPath);
        MGlobal::displayInfo("Mesh loaded: " + meshFn.name());

        // Get the connected shaders for the current mesh
        //Shaders: list of unique shaders in the mesh
        //shaderIndicies: length of faces in order and stores index to shader
        MIntArray shaderIndices;
        meshFn.getConnectedShaders(dagPath.instanceNumber(), shaders, shaderIndices); //this sets shader indices to the index in shaders for that face

        for (unsigned int i = 0; i < shaders.length(); ++i) {
            MObject shadingGroup = shaders[i];

            // Find the surface shader connected to the shading group
            MPlug surfaceShaderPlug = MFnDependencyNode(shadingGroup).findPlug("surfaceShader", true);
            MPlugArray connections;
            surfaceShaderPlug.connectedTo(connections, true, false);

            if (connections.length() > 0) {
                MObject shaderNode = connections[0].node();
                MFnDependencyNode shaderFn(shaderNode);

                // Display shader info directly without storing unique shaders
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


        // Retrieve vertex positions of the mesh
        MPointArray verts;
        meshFn.getPoints(verts, MSpace::kWorld);

        // Get face information (counts and vertex indices for each face)
        MIntArray faceCounts, faceConnects;
        meshFn.getVertices(faceCounts, faceConnects);

        // Print the vertex data for debugging
        for (unsigned int i = 0; i < verts.length(); ++i) {
            MString msg = MString("") + i + ": (" +
                verts[i].x + ", " +
                verts[i].y + ", " +
                verts[i].z + ")";
            MGlobal::displayInfo(msg);
        }

        // Print the face data (vertex indices for each face)
        MGlobal::displayInfo("Faces:");
        int idx = 0;
        for (unsigned int i = 0; i < faceCounts.length(); ++i) {
            char buffer[128];
            snprintf(buffer, sizeof(buffer), "Face %u:", i);
            MString faceMsg(buffer);
            // Loop through the vertices of the face and append them to the message
            for (int j = 0; j < faceCounts[i]; ++j) {
                faceMsg += " ";
                faceMsg += faceConnects[idx++];
            }
            MGlobal::displayInfo(faceMsg);
        }

        std::cout << "Loaded mesh with " << verts.length() << " vertices and "
            << faceCounts.length() << " faces.\n";

        // Convert vertex data to glm::vec3 for further processing
        std::vector<glm::vec3> vertexPositions;
        for (unsigned int i = 0; i < verts.length(); ++i) {
            vertexPositions.emplace_back(verts[i].x, verts[i].y, verts[i].z);
        }

        // Convert face data into pairs of vertex indices and corresponding shader index
        std::vector<std::pair<std::vector<int>, int>> faceData;
        idx = 0;
        for (unsigned int i = 0; i < faceCounts.length(); ++i) {
            std::vector<int> faceVerts;
            for (int j = 0; j < faceCounts[i]; ++j) {
                faceVerts.push_back(faceConnects[idx++]);
            }

            // Get the shader index associated with the current face
            int shaderIndex = shaderIndices[i];  // Use shaderIndices to get the shader index for this face
            faceData.emplace_back(faceVerts, shaderIndex);  // Store the vertex indices and the shader index as a pair
        }

        // Create a graph object using the vertex positions and face data
        s_loadedGraph = Graph(vertexPositions, faceData);
        s_hasLoadedGraph = true;
    }

    return MS::kSuccess;
}

// Creator method for the LoadMeshCmd class
void* LoadMeshCmd::creator() {
    return new LoadMeshCmd();
}
