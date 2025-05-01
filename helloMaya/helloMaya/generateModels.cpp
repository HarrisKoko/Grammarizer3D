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

    Graph& graph = LoadMeshCmd::getLoadedGraph();

    //std::vector<glm::vec3> verts = graph.getCachedPositions();

    
    std::map<unsigned int, glm::vec3> s;
    double d;
    std::vector<glm::vec3> verts = graph.samplePositions(s, d);
    std::vector<std::pair<std::vector<int>, int>> faces;
    std::vector<glm::vec3> faceColors({ {0.f,0.f,0.f},{1.f,1.f,1.f} });

    std::map<Primitive*, int> primIndex;

    // NOTE only works in case with edges all being parts of faces obviously
    for (unsigned int i = 0; i < graph.primitives.size(); ++i) {
        // for (const uPtr<Primitive>& p : m_currGraph.primitives) {
        Primitive* p = graph.primitives.at(i).get();
        //verts.push_back(p->cachedPos.value_or(glm::vec3(0)));
        primIndex.insert({ p, i });
    }


    std::set<std::pair<Primitive*, HalfEdgeGraph>> explored;
    for (unsigned int i = 0; i < graph.primitives.size(); ++i) {
        Primitive* startP = graph.primitives.at(i).get();
        for (unsigned int j = 0; j < startP->halfEdges.size(); ++j) {
            HalfEdgeGraph curHE = startP->halfEdges.at(j);
            if (curHE.angle < 0.f) {
                if (explored.contains({ startP, curHE })) {
                    continue;
                }
                std::pair<std::vector<int>, int> face;
                face.second = curHE.faces.at(0);
                //negative -> counterclockwise direction -> outwards
                Primitive* curP = startP;
                Primitive* nextP = curP->connections.at(j);
                do {
                    explored.insert({ curP, curHE });
                    // explored.insert({{curP, curHE},true});
                    // facePrims.push_back(curP);
                    face.first.push_back(primIndex.at(curP));
                    // std::vector<std::pair<float, int>> connToSort;
                    if (nextP == nullptr) {
                        break;
                    }
                    std::set<std::pair<float, int>> potentialNext;
                    for (unsigned int k = 0; k < nextP->halfEdges.size(); ++k) {
                        if (curHE.faces.at(0) == nextP->halfEdges.at(k).faces.at(0) &&
                            curHE.faceNormals.at(0) == nextP->halfEdges.at(k).faceNormals.at(0)) {
                            // on same plane in same rotational direction
                            // but might not necessarily be next step since could have multiple HEs on this plane; hence sort in set by angle
                            potentialNext.insert({ nextP->halfEdges.at(k).angle, k });
                        }
                    }
                    std::pair<float, int> nextSel = *potentialNext.begin();
                    curP = nextP;
                    nextP = curP->connections.at(nextSel.second);
                    curHE = curP->halfEdges.at(nextSel.second);
                } while (curP != startP);

                if (nextP == nullptr) {
                    continue;
                }


                faces.push_back(std::move(face));
            }
        }
    }






    // Get vertices and faces from graph
    std::vector<MPoint> vertices;  
    std::vector<int> faceCounts;  
    std::vector<int> faceConnects; 


    for (auto &face : faces) {
        faceCounts.push_back(face.first.size());
        faceConnects.insert(faceConnects.end(), face.first.begin(), face.first.end());
    }
    for (const glm::vec3& v : verts) {
    vertices.emplace_back(v.x, v.y, v.z);
}


    // Clears Maya window
    MGlobal::executeCommand("delete all");

    return createMesh(vertices, faceCounts, faceConnects);

}

void* GenerateModelsCmd::creator() {
    return new GenerateModelsCmd();
}