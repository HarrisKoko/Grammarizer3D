#include "generateModels.h"

MStatus GenerateModelsCmd::createMesh(const std::vector<MPoint>& vertices,
    const std::vector<int>& faceCounts,
    const std::vector<int>& faceConnects)
{
    MStatus status;

    // Convert to MPointArray and MIntArray
    MPointArray mPoints;
    for (const auto& pt : vertices) {
        mPoints.append(pt);
    }

    MIntArray mFaceCounts;
    for (int count : faceCounts) {
        mFaceCounts.append(count);
    }

    MIntArray mFaceConnects;
    for (int idx : faceConnects) {
        mFaceConnects.append(idx);
    }

    // Create mesh
    MFnMesh meshFn;
    MObject meshObj = meshFn.create(
        mPoints.length(),
        mFaceCounts.length(),
        mPoints,
        mFaceCounts,
        mFaceConnects,
        MObject::kNullObj,
        &status
    );

    if (!status) {
        MGlobal::displayError("Failed to create mesh");
        return status;
    }

    // Rename the mesh
    MFnDagNode dagNode(meshObj);
    dagNode.setName("pGeneratedMesh");

    return MS::kSuccess;
}




MStatus GenerateModelsCmd::doIt(const MArgList& args) {
    if (!LoadMeshCmd::hasLoadedGraph()) {
        return MS::kFailure;
    }

    Graph& graph = LoadMeshCmd::getLoadedGraph(); // Assumes copy constructor works. if issues, try making without.

    // Get vertices and faces from graph
    std::vector<MPoint> vertices;  
    std::vector<int> faceCounts;  
    std::vector<int> faceConnects; 

    // Clears Maya window
    MGlobal::executeCommand("delete all");

    return createMesh(vertices, faceCounts, faceConnects);

}

void* GenerateModelsCmd::creator() {
    return new GenerateModelsCmd();
}