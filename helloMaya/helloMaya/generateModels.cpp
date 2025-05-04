#include "generateModels.h"
#include <maya/MFnSet.h>

#include <maya/MFnMesh.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnSet.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MObjectArray.h>
#include <maya/MGlobal.h>
#include <maya/MDagPath.h>
#include <maya/MFnDagNode.h>
#include <maya/MItDag.h>
#include <maya/MFnComponentListData.h>




MStatus GenerateModelsCmd::createMesh(const std::vector<MPoint>& vertices,
    const std::vector<int>& faceCounts,
    const std::vector<int>& faceConnects,
    const std::vector<int>& shaderIndices)
{
    MStatus status;
    MGlobal::displayInfo("Before mesh creation");

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

    // Create the mesh
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
        MGlobal::displayError("Failed to create mesh.");
        return status;
    }

    MFnDagNode dagNode(meshObj);
    dagNode.setName("pGeneratedMesh");

    // Get shader list
    MObjectArray shaderList = LoadMeshCmd::getShaders();
    if (shaderList.length() == 0) {
        MGlobal::displayWarning("No shaders provided.");
        return MS::kSuccess;
    }

    // Log shader list for debugging
    MGlobal::displayInfo("Shader List Length: " + MString() + shaderList.length());
    for (unsigned int i = 0; i < shaderList.length(); ++i) {
        MGlobal::displayInfo("Shader " + MString() + i + ": " + shaderList[i].apiTypeStr());
    }

    // Get the DAG path to the mesh
    MDagPath meshPath;
    status = MDagPath::getAPathTo(meshObj, meshPath);
    if (!status) {
        MGlobal::displayError("Failed to get DAG path to mesh.");
        return status;
    }

    // Assign each face to the correct shader set
    for (unsigned int faceIdx = 0; faceIdx < shaderIndices.size(); ++faceIdx) {
        int shaderIdx = shaderIndices[faceIdx];

        if (shaderIdx < 0 || shaderIdx >= static_cast<int>(shaderList.length())) {
            MGlobal::displayWarning("Invalid shader index for face " + MString() + faceIdx);
            continue;
        }

        MObject shadingGroup = shaderList[shaderIdx];
        MFnSet shadingSet(shadingGroup, &status);
        if (!status) {
            MGlobal::displayWarning("Invalid shading group for face " + MString() + faceIdx);
            continue;
        }

        // Create a face component
        MFnSingleIndexedComponent compFn;
        MObject compObj = compFn.create(MFn::kMeshPolygonComponent, &status);
        if (!status) {
            MGlobal::displayWarning("Failed to create component for face " + MString() + faceIdx);
            continue;
        }

        compFn.addElement(faceIdx);

        // Add member using MDagPath + component
        status = shadingSet.addMember(meshPath, compObj);
        if (!status) {
            MGlobal::displayWarning("Failed to add face " + MString() + faceIdx + " to shading group.");
        }
    }

    return MS::kSuccess;
}






MStatus GenerateModelsCmd::doIt(const MArgList& args) {
    MGlobal::displayInfo("start do it");
    if (!LoadMeshCmd::hasLoadedGraph()) {
        return MS::kFailure;
    }
    Graph& graph = LoadMeshCmd::getLoadedGraph();

    //std::vector<glm::vec3> verts = graph.getCachedPositions();

    
    std::map<unsigned int, glm::vec3> s;
    double d;
    MGlobal::displayInfo("Before sample positions");
    std::vector<glm::vec3> verts = graph.samplePositions(s, d);
    std::vector<std::pair<std::vector<int>, int>> faces;

    std::map<Primitive*, int> primIndex;

    // NOTE only works in case with edges all being parts of faces obviously
    for (unsigned int i = 0; i < graph.primitives.size(); ++i) {
        // for (const uPtr<Primitive>& p : m_currGraph.primitives) {
        Primitive* p = graph.primitives.at(i).get();
        //verts.push_back(p->cachedPos.value_or(glm::vec3(0)));
        primIndex.insert({ p, i });
    }

    MGlobal::displayInfo("Before for loop");
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
                MGlobal::displayInfo("Before Face searching");
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


                            //awkward and probably some way to simplify but just doing:
                            // convert angle back to direction
                            // convert to angle within plane with face normal (0) (rather than normal with lowest valued xyz)
                            glm::vec3 direction(cos(nextP->halfEdges.at(k).angle), 0.0f, sin(nextP->halfEdges.at(k).angle));
                            glm::vec3 forward = glm::vec3(0, 1, 0);
                            glm::vec3 normal;
                            auto& norms = nextP->halfEdges.at(k).faceNormals;
                            if (norms.at(0).x == norms.at(1).x) {
                                if (norms.at(0).y == norms.at(1).y) {
                                    if (norms.at(0).z == norms.at(1).z) {
                                        normal = norms.at(0);

                                    }
                                    else {
                                        normal = (norms.at(0).z < norms.at(1).z) ? norms.at(0) : norms.at(1);
                                    }
                                }
                                else {
                                    normal = (norms.at(0).y < norms.at(1).y) ? norms.at(0) : norms.at(1);
                                }
                            }
                            else {
                                normal = (norms.at(0).x < norms.at(1).x) ? norms.at(0) : norms.at(1);
                            }

                            glm::vec3 turnDirection = glm::cross(normal, forward);
                            if (turnDirection != glm::vec3(0)) {
                                float turnAngle = glm::acos(glm::dot(normal, forward));

                                glm::mat4 rotation = glm::rotate(-turnAngle, turnDirection);

                                direction = glm::vec3(rotation * glm::vec4(direction, 1));

                            }


                            normal = nextP->halfEdges.at(k).faceNormals.at(0);
                            turnDirection = glm::cross(normal, forward);
                            if (turnDirection != glm::vec3(0)) {
                                float turnAngle = glm::acos(glm::dot(normal, forward));

                                glm::mat4 rotation = glm::rotate(turnAngle, turnDirection);

                                direction = glm::vec3(rotation * glm::vec4(direction, 1));
                            }
                            float theta = (direction.z > 0) ? glm::acos(direction.x) : -glm::acos(direction.x);


                            potentialNext.insert({ theta, k });
                        }
                    }
                    if (potentialNext.size() == 0) {
                        break; // shouldn't happen but just in case?
                    }

                    float invAngle;
                    {
                        glm::vec3 direction(cos(curHE.angle), 0.0f, sin(curHE.angle));
                        glm::vec3 forward = glm::vec3(0, 1, 0);
                        glm::vec3 normal;
                        auto& norms = curHE.faceNormals;
                        if (norms.at(0).x == norms.at(1).x) {
                            if (norms.at(0).y == norms.at(1).y) {
                                if (norms.at(0).z == norms.at(1).z) {
                                    normal = norms.at(0);

                                }
                                else {
                                    normal = (norms.at(0).z < norms.at(1).z) ? norms.at(0) : norms.at(1);
                                }
                            }
                            else {
                                normal = (norms.at(0).y < norms.at(1).y) ? norms.at(0) : norms.at(1);
                            }
                        }
                        else {
                            normal = (norms.at(0).x < norms.at(1).x) ? norms.at(0) : norms.at(1);
                        }

                        glm::vec3 turnDirection = glm::cross(normal, forward);
                        if (turnDirection != glm::vec3(0)) {
                            float turnAngle = glm::acos(glm::dot(normal, forward));

                            glm::mat4 rotation = glm::rotate(-turnAngle, turnDirection);

                            direction = glm::vec3(rotation * glm::vec4(direction, 1));

                        }


                        normal = curHE.faceNormals.at(0);
                        turnDirection = glm::cross(normal, forward);
                        if (turnDirection != glm::vec3(0)) {
                            float turnAngle = glm::acos(glm::dot(normal, forward));

                            glm::mat4 rotation = glm::rotate(turnAngle, turnDirection);

                            direction = glm::vec3(rotation * glm::vec4(direction, 1));
                        }
                        float theta = (direction.z > 0) ? glm::acos(direction.x) : -glm::acos(direction.x);

                        if (theta < 0) {
                            invAngle = theta + M_PI;
                        }
                        else {
                            invAngle = theta - M_PI;
                        }
                        
                    }
                    

                    /*potentialNext.rend();
                    for (auto it = potentialNext.rbegin(); it != potentialNext.rend(); ++it) {
                        std::pair<float, int> p = *it;
                    }*/
                        
                    //std::pair<float, int> nextSel = *(potentialNext.begin());
                    std::pair<float, int> nextSel = *std::prev(potentialNext.end());
                    // I BELIEVE this is what to do:
                    // pick greatest angle which is less than invAngle
                    // if none less than, pick the greatest overall
                    // so should pick the next one around clockwise
                    // I THINK the rotation might be reversed from what I'm expecting w/ y-axis being inverted from how I'm imagining?
                    // so do I want to go a different order then? inverse of that 

                    // TODO might be I need to go off of face [0] instead of [1]? not sure if makes a difference but if it is reversed fromw hat I think maybe do that
                    for (std::pair<float, int> p : potentialNext) {
                        /*if (invAngle + ANGLE_EPSILON < p.first) {
                            nextSel = p;
                            break;
                        }*/

                        if (p.first + ANGLE_EPSILON < invAngle) {
                            //if (p.first + ANGLE_EPSILON < invAngle) {
                            nextSel = p;
                        } else {
                            break;
                        }

                        //if ((p.first > M_PI - ANGLE_EPSILON * 2 && invAngle < 0) ? glm::epsilonEqual(p.first - float(M_PI) * 2, invAngle, ANGLE_EPSILON) :
                        //    ((invAngle > M_PI - ANGLE_EPSILON * 2 && p.first < 0) ? glm::epsilonEqual(p.first, invAngle - float(M_PI) * 2, ANGLE_EPSILON) :
                        //        glm::epsilonEqual(p.first, invAngle, ANGLE_EPSILON))) {
                        ////if (glm::epsilonEqual(p.first, invAngle, ANGLE_EPSILON)) {
                        //    break;
                        //}
                        //else {
                        //    nextSel = p;
                        //}
                    }
                    curP = nextP;
                    nextP = curP->connections.at(nextSel.second);
                    curHE = curP->halfEdges.at(nextSel.second);
                } while (curP != startP);

                if (nextP == nullptr) {
                    continue;
                }


                faces.push_back(std::move(face)); 
                // i think this already sets face label. 
                // face is a pair of the vector of positions and the int face label which is now the index in LoadMeshCmd::getShaders() list?
            }
        }
    }


    // Get vertices and faces from graph
    std::vector<MPoint> vertices;  
    std::vector<int> faceCounts;  
    std::vector<int> faceConnects;
    std::vector<int> faceLabels;


    for (auto &face : faces) {
        faceCounts.push_back(face.first.size());
        faceConnects.insert(faceConnects.end(), face.first.begin(), face.first.end());
        faceLabels.push_back(face.second);
    }
    for (const glm::vec3& v : verts) {
    vertices.emplace_back(v.x, v.y, v.z);
}


    // Clears Maya window
    MGlobal::executeCommand("delete all");

    return createMesh(vertices, faceCounts, faceConnects, faceLabels);

}

void* GenerateModelsCmd::creator() {
    MGlobal::displayInfo("Created CMD");
    return new GenerateModelsCmd();
}