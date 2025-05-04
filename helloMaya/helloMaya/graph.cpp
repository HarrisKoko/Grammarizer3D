#include "graph.h"
#include "jsonreader.h"
#include <iostream>
#include <random>

// #include <eigen-3.4.0/>
// #include <Eigen/Dense>
// #include <Eigen/SparseQR>
// #include <Eigen/Dense>
// #include <Eigen/Sparse>
// #include <Eigen/SparseQR>
#include <QFileDialog>
// #include <coin/ClpSimplex.hpp>
// #include <ortools/linear_solver/linear_solver.h>
// #include "ortools/linear_solver/linear_solver.h"
// #include <soplex.hpp>
// #include <lp_lib.h>
// #include <ortools_export.h>


unsigned int Graph::lastID{0};

Graph::Graph() : ID(lastID++) {
#if !THREE_DIMENSIONAL
    this->boundaryString.push_back(positive);
#endif
    // TODO boundary graph I think doesn't need anything for this?
    setText(QString::number(ID));
}

void Graph::resetIDs(){
    Graph::lastID = 0;
}

// trying a constructor from a general form which should be derivable from Maya's MFnMesh
// TODO note getConnectedShaders will be used to figure out labels I think (materials/shaders used correspond to 1 label each)
//  also: getPoints, getVertices (?) getPolygonVertices, getEdgeVertices
//  each face's data taken in as pair of vector of vertex indices and an integer label
//  TODO note could always change what the label is stored as
// TODO maybe use int[2] instead of pair<int,int>?
Graph::Graph(std::vector<glm::vec3> vertexPositions, std::vector<std::pair<std::vector<int>, int>> faceData)
: ID(lastID++)
{
    setText(QString::number(ID));
#if THREE_DIMENSIONAL

    std::vector<std::vector<std::tuple<unsigned int, int, int, glm::vec3, glm::vec3>>> adjacentVerts(vertexPositions.size());

#else
    std::vector<std::vector<std::array<int,3>>> adjacentVerts(vertexPositions.size());
#endif
    // std::vector<std::vector<int>> adjacentVerts(vertexPositions.size());
    // for (const std::pair<int,int> &edge : edgeVertexIndices) {
    //     adjacentVerts.at(edge.first).push_back(edge.second);
    //     adjacentVerts.at(edge.second).push_back(edge.first);
    // }
    // TODO figure out how want to store face labels? I guess (if we assume not allowing face-less edges)
    // can iterate over faceData and loop through vertices in order to construct adjacentVerts and have adjacentVerts store trios of ints instead of just 1 int
    //  make sure to skip adding to adjacent verts if already have BUT update one that's already in there to have the other face label.
    //  I think that'd work?

    // IGNORING edgeVertexIndices for now; TODO might just delete
    //  yeah deleted, was "std::vector<std::pair<int,int>> edgeVertexIndices"
    // TODO for now using label = -1 to mean no face there (edge of open mesh) -> should not draw faces with label -1 (just background color)
    for (const std::pair<std::vector<int>,int> &face : faceData) {
        glm::vec3 curFaceNormal;
        for (unsigned int i = 0; i < face.first.size(); ++i) {
            int start = face.first.at(i);
            int mid = face.first.at((i + 1) % face.first.size());
            int end = face.first.at((i + 2) % face.first.size());

            glm::vec3 startPos = vertexPositions.at(start);
            glm::vec3 midPos = vertexPositions.at(mid);
            glm::vec3 endPos = vertexPositions.at(end);

            curFaceNormal = glm::normalize(glm::cross(glm::normalize(endPos - midPos), glm::normalize(startPos - midPos)));
            if (curFaceNormal != glm::vec3(0)) {
                break;
            }
        }
        for (unsigned int i = 0; i < face.first.size(); ++i) {
            unsigned int start = face.first.at(i);
            unsigned int end = face.first.at((i + 1) % face.first.size());
            // TODO I think technically might be more efficient to just have this loop go to < size-1 then do final case after rather than modulo-ing? IDK how it compiles

            // TODO could use another loop for these below but I think since just 2 probably fine to just have two very similar bits of code
            bool yetToBeAdded = true;
#if THREE_DIMENSIONAL
            auto &startVec = adjacentVerts.at(start);
            for (unsigned int j = 0; j < startVec.size(); ++j) {
                if (std::get<0>(startVec.at(j)) == end) {
                    std::get<1>(startVec.at(j)) = face.second;
                    std::get<3>(startVec.at(j)) = curFaceNormal;
                    yetToBeAdded = false;
                    break;
                }
            }
            if (yetToBeAdded) {
                startVec.push_back({end, face.second, -1, curFaceNormal, glm::vec3(3)});
            }

            yetToBeAdded = true;
            auto &endVec = adjacentVerts.at(end);
            for (unsigned int j = 0; j < endVec.size(); ++j) {
                if (std::get<0>(endVec.at(j)) == start) {
                    std::get<2>(endVec.at(j)) = face.second;
                    std::get<4>(endVec.at(j)) = curFaceNormal;
                    yetToBeAdded = false;
                    break;
                }
            }
            if (yetToBeAdded) {
                endVec.push_back({start, -1, face.second, glm::vec3(3), curFaceNormal});
            }
#else

            std::vector<std::array<int,3>> &startVec = adjacentVerts.at(start);
            // TODO could try using find_if for this I think but I know how to construct it with for better than I know how to use that
            for (unsigned int j = 0; j < startVec.size(); ++j) {
                if (startVec.at(j).at(0) == end) {
                    startVec.at(j).at(1) = face.second;
                    yetToBeAdded = false;
                    break;
                }
            }
            if (yetToBeAdded) {
                startVec.push_back({end, face.second, -1});
            }

            yetToBeAdded = true;
            std::vector<std::array<int,3>> &endVec = adjacentVerts.at(end);
            for (unsigned int j = 0; j < endVec.size(); ++j) {
                if (endVec.at(j).at(0) == start) {
                    endVec.at(j).at(2) = face.second;
                    yetToBeAdded = false;
                    break;
                }
            }
            if (yetToBeAdded) {
                endVec.push_back({start, -1, face.second});
            }
#endif



            // NOTE requires faceData to have vertices in consistently either clockwise or counterclockwise order
            //  since need one of the adjacent faces to have one vert first and the other the other vert first when cycling in order
            // TODO if Maya doesn't give us output that way will have to change
            // If this is using counterclockwise, then I believe [1] should be left label and [2] right label
            //  shouldn't really matter which as long as it's consistent throughout


        }
    }

    for (unsigned int i = 0; i < vertexPositions.size(); ++i) {
        this->primitives.push_back(mkU<Primitive>());
    }

    for (unsigned int i = 0; i < vertexPositions.size(); ++i) {
        uPtr<Primitive> &p = this->primitives.at(i);

        glm::vec3 iPos = vertexPositions.at(i);

        // p->cachedPos = iPos; // just for displaying base graph well; not used in construction from grammar

        // add halfEdge for each adjacent vertex
        for (const auto &adjData : adjacentVerts.at(i)) {
            // for (const std::array<int,3> &adjData : adjacentVerts.at(i)) {
#if THREE_DIMENSIONAL
            glm::vec3 dir = glm::normalize(vertexPositions.at(std::get<0>(adjData)) - iPos);

            glm::vec3 forward = glm::vec3(0,1,0);
            glm::vec3 normal;
            if (std::get<3>(adjData).x == std::get<4>(adjData).x) {
                if (std::get<3>(adjData).y == std::get<4>(adjData).y) {
                    if (std::get<3>(adjData).z == std::get<4>(adjData).z) {
                        normal = std::get<3>(adjData);

                    } else {
                        normal = (std::get<3>(adjData).z < std::get<4>(adjData).z) ? std::get<3>(adjData) : std::get<4>(adjData);
                    }
                } else {
                    normal = (std::get<3>(adjData).y < std::get<4>(adjData).y) ? std::get<3>(adjData) : std::get<4>(adjData);
                }
            } else {
                normal = (std::get<3>(adjData).x < std::get<4>(adjData).x) ? std::get<3>(adjData) : std::get<4>(adjData);
            }

            glm::vec3 turnDirection = glm::cross(normal, forward);
            if (turnDirection != glm::vec3(0)) {
                float turnAngle = glm::acos(glm::dot(normal, forward));

                glm::mat4 rotation = glm::rotate(turnAngle, turnDirection);

                dir = glm::vec3(rotation * glm::vec4(dir, 1));
            //} else if (glm::dot(normal, forward) < 0) {
                // dir.x = -dir.x; // TODO not sure if that's totally right
            }
            float theta = (dir.z > 0) ? glm::acos(dir.x) : -glm::acos(dir.x);

            // rounding before so can use map method more reliably
            if (theta > 3.141) {
                theta -= 2 * M_PI;
            }
            theta = std::round(theta * 180.f) / 180.f;
            // ^TODO probably make dir.y and have this take in 2d vec? or generaliaze to something like vector for angle and leave in 3d
            // TODO maybe flip around which way is negative to match paper diagrams better (since -z is up on our display rn) but fine as is as long as consistent
            p->halfEdges.push_back(HalfEdgeGraph(theta, {{std::get<1>(adjData), std::get<2>(adjData)}}, {{std::get<3>(adjData), std::get<4>(adjData)}}));
            p->connections.push_back(this->primitives.at(std::get<0>(adjData)).get());
#else

            glm::vec3 dir = glm::normalize(vertexPositions.at(adjData.at(0)) - iPos);
            float theta = (dir.z > 0) ? glm::acos(dir.x) : -glm::acos(dir.x);
            // ^TODO probably make dir.y and have this take in 2d vec? or generaliaze to something like vector for angle and leave in 3d
            // TODO maybe flip around which way is negative to match paper diagrams better (since -z is up on our display rn) but fine as is as long as consistent
            p->halfEdges.push_back(HalfEdgeGraph(theta, {{adjData.at(1), adjData.at(2)}}));
            p->connections.push_back(this->primitives.at(adjData.at(0)).get());
#endif

        }
    }

    // this constructor should result in complete graphs (no incomplete edges in input geometry), so boundary string is just:
#if !THREE_DIMENSIONAL
    this->boundaryString.push_back(positive);
#endif

}

Graph::Graph(const Graph &other)
: ID(lastID++)
{
    setText(QString::number(ID));
    this->isSingleEdge = other.isSingleEdge;
    this->face_colors = other.face_colors;
    this->boundaryString = other.boundaryString;
    this->primitives = std::vector<uPtr<Primitive>>();
    std::map<Primitive*, Primitive*> correspondences;
    for (const uPtr<Primitive>& p : other.primitives) {
        uPtr<Primitive> newPrim = mkU<Primitive>();
        // newPrim->halfEdges = p->halfEdges; // TODO can't recall if this copies right
        newPrim->halfEdges = std::vector<HalfEdgeGraph>(p->halfEdges);
        newPrim->cachedPos = p->cachedPos;
        correspondences.insert({p.get(), newPrim.get()});
        this->primitives.push_back(std::move(newPrim));
    }

    for (const auto& [p, newPrim] : correspondences) {

    // }

    // for (const uPtr<Primitive>& p : other.primitives) {
        try {
            // Primitive* newPrim = correspondences.at(p.get());
            for (unsigned int j = 0; j < p->connections.size(); ++j) {
                Primitive* connect = p->connections.at(j);
                if (connect == nullptr) {
                    newPrim->connections.push_back(nullptr);
                } else {
                    newPrim->connections.push_back(correspondences.at(connect));
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "graph copy constructor fail" << std::endl;
            std::cerr << e.what() << std::endl;
        }
    }
#if THREE_DIMENSIONAL

    this->boundaryGraphElements = std::vector<uPtr<BoundaryElem>>();
    std::map<BoundaryElem*, BoundaryElem*> elemCorrespondences;
    for (const uPtr<BoundaryElem>& elem : other.boundaryGraphElements) {
        // uPtr<BoundaryElem> newElem = mkU<BoundaryElem>(elem->he);
        // newElem->nextTurns = elem->nextTurns;

        uPtr<BoundaryElem> newElem = mkU<BoundaryElem>(*elem);

        elemCorrespondences.insert({elem.get(), newElem.get()});
        this->boundaryGraphElements.push_back(std::move(newElem));
    }

    for (const uPtr<BoundaryElem>& newElem : this->boundaryGraphElements) {
        newElem->next = elemCorrespondences.at(newElem->next);
        newElem->prev = elemCorrespondences.at(newElem->prev);
    }

    // this->sortBoundaryGraphElements();

    this->boundaryGeneric = other.boundaryGeneric;
#endif


}

Graph &Graph::operator=(const Graph & other)
{
    setText(QString::number(other.ID));
    // Not sure why it can't get this automatically from the copy constructor being defined? But seems to work fine this way.
    if (this == &other) {
        return *this;
    }
    this->isSingleEdge = other.isSingleEdge;
    this->face_colors = other.face_colors;
    this->boundaryString = other.boundaryString;
    this->primitives = std::vector<uPtr<Primitive>>();
    std::map<Primitive*, Primitive*> correspondences;
    for (const uPtr<Primitive>& p : other.primitives) {
        uPtr<Primitive> newPrim = mkU<Primitive>();
        newPrim->halfEdges = std::vector<HalfEdgeGraph>(p->halfEdges); // TODO can't recall if this copies right
        newPrim->cachedPos = p->cachedPos;
        correspondences.insert({p.get(), newPrim.get()});
        this->primitives.push_back(std::move(newPrim));
    }

    for (const uPtr<Primitive>& p : other.primitives) {
        Primitive* newPrim = correspondences.at(p.get());
        // for (Primitive* connect : p->connections) {
        for (unsigned int j = 0; j < p->connections.size(); ++j) {
            Primitive* connect = p->connections.at(j);
            if (connect == nullptr) {
                newPrim->connections.push_back(nullptr);
            } else {
                newPrim->connections.push_back(correspondences.at(connect));
            }
        }
    }

#if THREE_DIMENSIONAL

    this->boundaryGraphElements = std::vector<uPtr<BoundaryElem>>();
    std::map<BoundaryElem*, BoundaryElem*> elemCorrespondences;
    for (const uPtr<BoundaryElem>& elem : other.boundaryGraphElements) {
        uPtr<BoundaryElem> newElem = mkU<BoundaryElem>(*elem);

        elemCorrespondences.insert({elem.get(), newElem.get()});
        this->boundaryGraphElements.push_back(std::move(newElem));
    }

    for (const uPtr<BoundaryElem>& newElem : this->boundaryGraphElements) {
        newElem->next = elemCorrespondences.at(newElem->next);
        newElem->prev = elemCorrespondences.at(newElem->prev);
    }

    // this->sortBoundaryGraphElements();

    this->boundaryGeneric = other.boundaryGeneric;

#endif


    return *this;
}

// Graph::Graph(std::unique_ptr<Primitive> prim)
// {


//     this->primitives.push_back(std::move(prim));
//     // this->boundaryString.push_back()
//     // TODO I think want to make sure halfEdges sorted in order of increasing angle value (-180 degrees -> 180 degrees), then boundary string = half edges all in order then positive turn
//     //  but since using same indices for halfEdges and connections I guess should do that when constructing those elsewhere?
//     //  perhaps easier just to make primitive instead store pairs of HEs and connections
//     //  then can just apply std::sort with order based on each .first.angle

//     // TODO actually I think don't quite want this as constructor perhaps, maybe instead just pass in the vector of HEs? since don't want connections in the single-primitive graph which this is meant to construct
// }

Graph::Graph(const std::vector<HalfEdgeGraph>& primHEs)
    : ID(lastID++)
{
    setText(QString::number(ID));
    std::unique_ptr<Primitive> prim = mkU<Primitive>();

    prim->halfEdges = std::vector<HalfEdgeGraph>(primHEs); // TODO am I correct in thinking that already works by copying?
    prim->connections = std::vector<Primitive *>(primHEs.size());

    std::sort(prim->halfEdges.begin(), prim->halfEdges.end());
    // TODO not sure if I want to sort there or just use a sorted copy type of thing? I think fine to do here since already know no connections to mess up indices of.
    // note should sort based on the operation< defined in HalfEdgeGraph, so should be ascending angle (then ascending face labels just to make order well defined but not relevant since should never have two of same angle in same primitive)
#if THREE_DIMENSIONAL

    std::vector<uPtr<BoundaryElem>> yetToAdd;
    for (const HalfEdgeGraph& he : prim->halfEdges) {
        BoundaryHE bh;
        bh.h = he;
        bh.primIndex = 0;
        // elem->he = bh;
        uPtr<BoundaryElem> elem = mkU<BoundaryElem>(bh);

        // elem->next = nullptr;
        // elem->prev = nullptr;

        bool placed = false;
        std::pair<int, glm::vec3> leftFace = {he.faces[0], he.faceNormals[0]};
        std::pair<int, glm::vec3> rightFace = {he.faces[1], he.faceNormals[1]};
        for (unsigned int i = this->boundaryGraphElements.size(); i > 0; --i) {
            // for (unsigned int i = 0; i < this->boundaryGraphElements.size(); ++i) {
            const uPtr<BoundaryElem>& otherElem = this->boundaryGraphElements.at(i - 1);
            std::pair<int, glm::vec3> otherLeftFace = {otherElem->he.h.faces[0], otherElem->he.h.faceNormals[0]};
            std::pair<int, glm::vec3> otherRightFace = {otherElem->he.h.faces[1], otherElem->he.h.faceNormals[1]};
            if (otherRightFace == leftFace) {
                if (otherLeftFace == rightFace) {
                    // same face
                    //  place in order of sorted HEs, so I guess place after?
                    this->boundaryGraphElements.insert(this->boundaryGraphElements.begin() + i, std::move(elem));
                } else {
                    // place before
                    this->boundaryGraphElements.insert(this->boundaryGraphElements.begin() + i - 1, std::move(elem));
                }
                placed = true;
                break;
            } else {
                if (otherLeftFace == rightFace) {
                    // place after
                    this->boundaryGraphElements.insert(this->boundaryGraphElements.begin() + i, std::move(elem));
                    placed = true;
                    break;
                } else {
                    // not next to this edge
                }
            }
        }
        if (!placed) {
            if (this->boundaryGraphElements.size() == 0) {
                this->boundaryGraphElements.push_back(std::move(elem));
            } else {
                // place later, once neighbor is in vector (might be a better way to do this but eh)
                yetToAdd.push_back(std::move(elem));
            }
        }
    }

    while (yetToAdd.size() > 0) {
        unsigned int curSize = yetToAdd.size();

        // for (uPtr<BoundaryElem>& elem : yetToAdd) {
        for (unsigned int j = 0; j < yetToAdd.size(); ++j) {
            uPtr<BoundaryElem>& elem = yetToAdd.at(j);
            // for (uPtr<BoundaryElem>& elem : yetToAdd) {
            std::pair<int, glm::vec3> leftFace = {elem->he.h.faces[0], elem->he.h.faceNormals[0]};
            std::pair<int, glm::vec3> rightFace = {elem->he.h.faces[1], elem->he.h.faceNormals[1]};
            for (unsigned int i = this->boundaryGraphElements.size(); i > 0; --i) {
                // for (unsigned int i = 0; i < this->boundaryGraphElements.size(); ++i) {
                const uPtr<BoundaryElem>& otherElem = this->boundaryGraphElements.at(i - 1);
                std::pair<int, glm::vec3> otherLeftFace = {otherElem->he.h.faces[0], otherElem->he.h.faceNormals[0]};
                std::pair<int, glm::vec3> otherRightFace = {otherElem->he.h.faces[1], otherElem->he.h.faceNormals[1]};
                if (otherRightFace == leftFace) {
                    if (otherLeftFace == rightFace) {
                        // same face
                        //  place in order of sorted HEs, so I guess place after?
                        this->boundaryGraphElements.insert(this->boundaryGraphElements.begin() + i, std::move(elem));
                        yetToAdd.erase(yetToAdd.begin() + j);
                        --j;
                    } else {
                        // place before
                        this->boundaryGraphElements.insert(this->boundaryGraphElements.begin() + i - 1, std::move(elem));
                        yetToAdd.erase(yetToAdd.begin() + j);
                        --j;
                    }
                    break;
                } else {
                    if (otherLeftFace == rightFace) {
                        // place after
                        this->boundaryGraphElements.insert(this->boundaryGraphElements.begin() + i, std::move(elem));
                        yetToAdd.erase(yetToAdd.begin() + j);
                        --j;
                        break;
                    } else {
                        // not next to this edge
                    }
                }
            }
        }

        if (yetToAdd.size() == curSize) {
            // fail, don't think should be possible
            break;
        }
    }


    for (unsigned int i = 0; i < this->boundaryGraphElements.size(); ++i) {
        this->boundaryGraphElements.at(i)->next = this->boundaryGraphElements.at((i + 1) % this->boundaryGraphElements.size()).get();
        this->boundaryGraphElements.at((i + 1) % this->boundaryGraphElements.size())->prev = this->boundaryGraphElements.at(i).get();
    }

    this->sortBoundaryGraphElements();

#else
    for (const HalfEdgeGraph& he : prim->halfEdges) {
        BoundaryHE bh;
        bh.h = he;
        bh.primIndex = 0;
        this->boundaryString.push_back(bh);
    }
    this->boundaryString.push_back(positive);
#endif

    this->primitives.push_back(std::move(prim));
}

bool checkMatch(std::set<std::pair<Primitive*,Primitive*>>* encountered, Primitive* fromPrim, Primitive* thisPrim, std::vector<std::tuple<Primitive*, HalfEdgeGraph, Primitive*, unsigned int>>* oCorrespondences) {
// bool checkMatch(std::map<Primitive*,Primitive*>* encountered, Primitive* fromPrim, Primitive* thisPrim) {
    if (thisPrim == nullptr || fromPrim == nullptr) {
        // TODO maybe should make return true if both nullptr?
        return false;
    }
    if (fromPrim->sameHEs(*thisPrim)) {
        encountered->insert({fromPrim, thisPrim});
        for (unsigned int j = 0; j < fromPrim->connections.size(); ++j) {
            Primitive* from2 = fromPrim->connections.at(j);
            for (unsigned int k = 0; k < thisPrim->connections.size(); ++k) {
                if (fromPrim->halfEdges.at(j) == thisPrim->halfEdges.at(k)) {
                    Primitive* this2 = thisPrim->connections.at(k);
                    if (from2 == nullptr) {
                        // TODO add to correspondence list
                        // STORE:
                            // fromPrim;
                            // fromPrim->halfEdges.at(j);
                            // this2;
                        // then afterwards, when doing substitution, iterate through boundary strings;
                        //  when reaching element with index leading to stored fromPrim,
                        //   go to prim in to which has index in corresponding place in string
                        // then for its newly copied node w/ HE which matches stored HE labels, set connection to this2
                        //  and set this2 connection to that node

                        // TODO not totally sure if should just skip this whole thing in both nullptr case or not
                        if (this2 != nullptr) {
                            unsigned int connectionIndex = 0;
                            for (unsigned int i = 0; i < this2->connections.size(); ++i) {
                                if (this2->connections.at(i) == thisPrim) {
                                    connectionIndex = i;
                                    break;
                                }
                            }

                            oCorrespondences->push_back({fromPrim, fromPrim->halfEdges.at(j), this2, connectionIndex});

                            // TODO maybe need to add check that no this2 is encountered both with a nullptr and non-nullptr pair? not sure if possible
                            encountered->insert({nullptr, this2}); // TODO not sure if this is a case that ever matters so not sure if this accounts for it right
                        }
                    } else {
                        // TODO do I also need to check in other direction?
                        //  feels like maybe need to in which case map maybe not best approach (or need 2 maps?)
                        //  but not sure that's actually the case since IDK if a graph structure where just one direction is wrong
                        //      I guess if fromPrim should go to a prim we haven't encountered yet but thisPrim goes to one we have encountered this is wrong but can that happen? I feel like maybe it can but already would cause a false elsewhere in structure?
                        //  yeah I think I do need to check both directions.
                        // encountered->find()
                        // encountered->find()
                        auto foundFrom = std::find_if(encountered->begin(), encountered->end(), [&from2](std::pair<Primitive*,Primitive*> a) {
                            return a.first == from2;
                        });
                        // auto found = encountered->find({from2,this2});
                        if (foundFrom == encountered->end()) {
                            auto foundTo = std::find_if(encountered->begin(), encountered->end(), [&this2](std::pair<Primitive*,Primitive*> a) {
                                return a.second == this2;
                            });
                            if (foundTo == encountered->end()) {
                                // not already encountered

                                if (!checkMatch(encountered, from2, this2, oCorrespondences)) {
                                    return false;
                                }
                            } else {
                                // mismatch
                                return false;
                            }
                        } else {
                            // already encountered, check if same pair
                            // TODO do I need to also do the other find here? I think not an issue
                            if (this2 != foundFrom->second) {
                                // mismatch
                                return false;
                            }
                        }
                    }
                    break;
                }
            }
        }
        return true;
    } else {
        return false;
    }
}

// TODO should make sure if this works with a rule that lacks an edge between two points that are both in from and have an edge in this
//  I think it won't at the moment; I think that's what's causing connection issues
//  after that should try to add the checks for irreducible descendants and should generally optimize the checks in hierarchy construction since can get pretty slow
std::optional<Graph> Graph::applyReplacementRule(const Graph& from, const Graph& to)
{

    Graph result = Graph(*this);

    // NOTE: I think don't need to worry about making this function be able to be applied to an incomplete graph
    //  from and to can and usually are not complete (have nullptr connections--half-edges able to be glued to)
    //  but 'this' shouldn't ever be assuming proper use of algorithm (since we always start with starter rules which produce complete graphs, all replacements apply to complete graphs)
    //  so not worrying about potential updating of boundary string of 'this' since should just be {positive} always when at this step of algorithm

    // std::cout << from.isSingleEdge << " -> " << to.isSingleEdge << std::endl;
    // TODO attempt to apply replacement
    //  search through for subgraph matching "from"
    //      matches if primitives in from have same HE as prims in this w/ appropriate connections
    //      nullptr connections in from can connect to anything in this; update substituted primitives to have same connection as this did
    //      construct new subgraph w/ connections from ^
    //          but TODO how do I figure out which connections in from match those in to? I guess could pass in separately since can find that in construction of grammar I think
    //  return true if finds an appropriate substitution
    try {
        if (from.primitives.size() == 0) {
            // starter rule case
            //  constructs a separate graph I guess? but stored in the same Graph object


            std::map<Primitive*, Primitive*> toToThisPrimitives; // TODO make a name that is less horrible maybe
            for (const uPtr<Primitive>& p : to.primitives) {
                uPtr<Primitive> newPrim = mkU<Primitive>();
                newPrim->halfEdges = p->halfEdges; // TODO can't recall if this copies right
                toToThisPrimitives.insert({p.get(), newPrim.get()});
                result.primitives.push_back(std::move(newPrim));
            }

            for (const uPtr<Primitive>& p : to.primitives) {
                Primitive* newPrim = toToThisPrimitives.at(p.get());
                // for (Primitive* connect : p->connections) {
                for (unsigned int j = 0; j < p->connections.size(); ++j) {
                    Primitive* connect = p->connections.at(j);
                    if (connect == nullptr) {
                        newPrim->connections.push_back(nullptr);
                        // std::cout << "BAD " << p->connections.size() << " " << newPrim->connections.size() << std::endl;

                        // return false; // should never have nullptr in starter rule case
                        return {};
                    } else {
                        newPrim->connections.push_back(toToThisPrimitives.at(connect));
                    }
                }
                // std::cout << "STA " << p->connections.size() << " " << newPrim->connections.size() << std::endl;

            }

            // TODO need to make tests for this being right

            // return true;
            return result;


            // TODO another case to consider but not sure how to express: single edge (i.e. graph contains just two half edges in opposite directions) I believe is considered valid in paper examples
        } else if (from.isSingleEdge) {
            // SPECIAL CASE: DOES NOT actually have any primitives in it
            //  but representing as having a single primitive to store the data (for better display purposes)
            //      should have exactly the two opposing edges

            // TODO test this case

            if (from.primitives.size() != 1 || from.primitives.at(0)->halfEdges.size() != 2) {
                // requirements for structure of singleEdge case Graph object
                // return false;
                return {};
            }
            if (to.isSingleEdge || to.primitives.size() == 0) {
                // should never be able to construct rules like this I believe but just in case
                // note if both from and to are single edge would have to be the same single edge, so pointless to apply
                // return false;
                return {};
            }

            std::vector<std::pair<Primitive*,Primitive*>> potentialPairs;
            std::vector<std::pair<int,int>> potentialIndices;
            for (const uPtr<Primitive>& curPrim : result.primitives) {
                for (unsigned int i = 0; i < curPrim->halfEdges.size(); ++i) {
                    if (curPrim->halfEdges.at(i) == from.primitives.at(0)->halfEdges.at(0) && curPrim->connections.at(i) != nullptr) {
                        Primitive* adjPrim = curPrim->connections.at(i);
                        for (unsigned int j = 0; j < adjPrim->halfEdges.size(); ++j) {
                            if (adjPrim->halfEdges.at(j) == from.primitives.at(0)->halfEdges.at(1) && adjPrim->connections.at(j) == curPrim.get()) {
                                // found a pair (curPrim, adjPrim) which has an edge between them matching "from"
                                potentialPairs.push_back({curPrim.get(),adjPrim});
                                potentialIndices.push_back({i,j});

                            }
                        }
                    }
                }
            }


            if (potentialPairs.size() > 0) {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<> distrib(0, potentialPairs.size() - 1);
                int selectedStart = distrib(gen);

                Primitive* curPrim = potentialPairs.at(selectedStart).first;
                Primitive* adjPrim = potentialPairs.at(selectedStart).second;
                int i = potentialIndices.at(selectedStart).first;
                int j = potentialIndices.at(selectedStart).second;

                //  construct copy of set of vertices from "to", connect to this pair

                std::map<Primitive*, Primitive*> toToThisPrimitives; // TODO make a name that is less horrible maybe
                for (const uPtr<Primitive>& p : to.primitives) {
                    uPtr<Primitive> newPrim = mkU<Primitive>();
                    newPrim->halfEdges = p->halfEdges; // TODO can't recall if this copies right
                    toToThisPrimitives.insert({p.get(), newPrim.get()});
                    result.primitives.push_back(std::move(newPrim));
                }


                for (const uPtr<Primitive>& p : to.primitives) {
                    Primitive* newPrim = toToThisPrimitives.at(p.get());
                    for (unsigned int k = 0; k < p->connections.size(); ++k) {
                        Primitive* connect = p->connections.at(k);
                        if (connect == nullptr) {
                            HalfEdgeGraph he = p->halfEdges.at(k);
                            // boundary of "to" graph
                            if (he == curPrim->halfEdges.at(i)) {
                                adjPrim->connections.at(j) = newPrim;
                                newPrim->connections.push_back(adjPrim);
                            } else {
                                curPrim->connections.at(i) = newPrim;
                                newPrim->connections.push_back(curPrim);

                            }
                        } else {
                            newPrim->connections.push_back(toToThisPrimitives.at(connect));
                        }
                    }
                    // std::cout << "EDG " << p->connections.size() << " " << newPrim->connections.size() << std::endl;
                }


                return result;//true;
            }


            // return false;
            return {};
        } else {
            if (result.primitives.size() == 0) {
                // return false;
                return {};
            }
            // Primitive* searchStart = from.primitives.at(0).get();

            // find matching subgraphs of "this" w/ "from"
            std::vector<Primitive*> potentialStarts;
            std::vector<std::vector<Primitive*>> potentialMatches;
            std::vector<std::vector<std::tuple<Primitive*, HalfEdgeGraph, Primitive*, unsigned int>>> potentialBoundaryCorrespondences;

            std::vector<Primitive*> primsToMatch;
            for (const uPtr<Primitive>& p : from.primitives) {
                primsToMatch.push_back(p.get());
            }
            for (unsigned int i = 0; i < result.primitives.size(); ++i) {
                Primitive* curPrim = result.primitives.at(i).get();
                // std::vector<std::pair<Primitive*, Primitive*>> fromThisCorrespondence;
                // IDEA: map with key = primitive in from, value = matching primitive in this
                // std::map<Primitive*,Primitive*> encountered;
                std::vector<Primitive*> unmatchedPrims = primsToMatch;

                std::set<std::pair<Primitive*,Primitive*>> encountered;

                // trying recursively since might be easier to write

                std::vector<std::tuple<Primitive*, HalfEdgeGraph, Primitive*, unsigned int>> boundaryCorrespondences;

                bool isMatched = true;
                // TODO need to test that this works on disconnected "from" graphs, don't have a good way to at the moment; get back to once hierarchy setup working
                while (isMatched && unmatchedPrims.size() > 0) {
                    isMatched = checkMatch(&encountered, unmatchedPrims.at(0), curPrim, &boundaryCorrespondences);
                    for (const std::pair<Primitive*,Primitive*>& pair : encountered) {
                        // unmatchedPrims.clear()
                        std::erase_if(unmatchedPrims, [&pair](Primitive* curPrim) {
                            return curPrim == pair.first;
                        });
                    }
                }

                if(isMatched) {
                    std::vector<Primitive*> matches;
                    for (auto [keyFrom,valThis]: encountered) {
                        if (keyFrom != nullptr) {
                            matches.push_back(valThis);
                        }
                    }
                    bool hitSelf = false;
                    // TODO not actually sure we always want to skip these cases
                    // TODO hitting some map::at out of range now, not sure where from, TODO will look at tomorrow
                    //  Doesn't always hit it so hard to reproduce
                    for (const auto& [fromPrim, fromHE, thisPrim, thisConnectIndex] : boundaryCorrespondences) {
                        for (Primitive* p : matches) {
                            if (p == thisPrim) {
                                hitSelf = true;
                                break;
                            }
                        }
                        if (hitSelf) {
                            break;
                        }
                    }
                    if (!hitSelf) {
                        potentialStarts.push_back(curPrim);
                        potentialMatches.push_back(matches);
                        potentialBoundaryCorrespondences.push_back(boundaryCorrespondences);
                    }
                }
                // if (curPrim->sameHEs(*searchStart)) {
                //     // perform search out to check if rest of from matches elements in this
                //     // TODO how do I track where the nullptr substitutions should be though?
                //     //  I think perhaps I need to have boundary strings working first?
                //     //  make boundary strings contain pointers to the primative+HE index the HE label originates from
                //     //  then can use from and to's boundary strings to compare which primative's connection corresponds to which
                //     //      since all the nullptr connections should correspond to the boundary string

                //     for (unsigned int j = 0; j < curPrim->connections.size(); ++j) {

                //     }
                //     // thinking:
                //     //  if all same HEs
                //     //      store pair
                //     //      for each non-null connection:
                //     //          if already encountered
                //     //              check if same pair
                //     //          else
                //     //              repeat
                //     //      for each null connection:
                //     //          store correspondence of which he and primitive at in from w/ where the connection in this goes
                //     // I guess then need to store some sort of way to see if primitives already explored
                //     //  IDEA: store pairs of corresponding primitives encountered in from and this
                // }
                // perhaps I ought to just use some existing library for this but not sure structure-wise how to incorporate
                //  subgraph isomorphism problem
                //  double pushout graph rewriting
                //  though efficiency isn't really a concern here so don't need a super-optimized approach
                //      more a question of what's easier to implement. easier to fit to exact needs when writing by self but harder to make algorithm correct obviously

                // since want replacements to be able to apply to any of the matching subgraphs in graph, should fo matching checks for all of them then randomly decide between the successful results, I guess
            }

            // TODO track which primitives were not explored by the search and repeat for those?
            //  the rules can generate non-connected graph structures I believe so this is necessary
            //  first implementing with just connected structures for now, so still TODO this
            //      I think: rather than just starting at 0, have another loop
            //      perform same search for matches but ignore all already matched primitives in both from and this


            if (potentialMatches.size() == 0) {
                // return false;
                return {};
            }


            //  line up boundary strings for use in comparison
            // TODO ^ then can use same index to check corresponding primitives
#if THREE_DIMENSIONAL
            // from.sortBoundaryGraphElements();
            // actually pre-sorting in glue operations now I think makes sense, so should already be aligned assuming sort correctly makes consistent order
#else
            unsigned int fromStringStart = 0;
            while (fromStringStart < from.boundaryString.size()) {
                bool stringMatches = true;
                for (unsigned int i = 0; i < from.boundaryString.size(); ++i) {
                    const Turn* tTo = std::get_if<Turn>(&to.boundaryString.at(i));
                    const Turn* tFrom = std::get_if<Turn>(&from.boundaryString.at((fromStringStart + i) % from.boundaryString.size()));
                    if (tTo == nullptr && tFrom == nullptr) {
                        // half edge symbol
                        const BoundaryHE* bTo = std::get_if<BoundaryHE>(&to.boundaryString.at(i));
                        const BoundaryHE* bFrom = std::get_if<BoundaryHE>(&from.boundaryString.at((fromStringStart + i) % from.boundaryString.size()));
                        if (bTo->h != bFrom->h) {
                            stringMatches = false;
                            break;
                        }
                    } else {
                        // turn symbol
                        if (tTo == nullptr || tFrom == nullptr || *tTo != *tFrom) {
                            stringMatches = false;
                            break;
                        }
                    }

                }
                if (stringMatches) {
                    break;
                } else {
                    ++fromStringStart;
                }
            }
            if (fromStringStart >= from.boundaryString.size()) {
                // nonmatching boundary strings -> failure
                // return false;
                return {};
            }
#endif

            // select which match to apply replacement to
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> distrib(0, potentialStarts.size() - 1);
            int selectedStart = distrib(gen);


            // replace selected subgraph with "to"
            //  remove matched primitives from "this" (all those in potentialMatches[selectedStart])


            if (to.isSingleEdge && potentialBoundaryCorrespondences.at(selectedStart).size() != 2) {
                // shouldn't happen (single edge graph should always have 2 boundary HEs) but just in case
                // return false;
                return {};
            }



            // TODO test this case
            if (to.isSingleEdge) {
                for (Primitive* p : potentialMatches.at(selectedStart)) {
                    // TODO maybe check if returned value == 1 here: // <- I think I meant if it did actually delete it exactly once
                    std::erase_if(result.primitives, [&p](uPtr<Primitive>& curPrim) {
                        return curPrim.get() == p;
                    });
                }

                auto& [fromPrim, fromHE, thisPrim, thisConnectIndex] = potentialBoundaryCorrespondences.at(selectedStart).at(0);
                auto& [fromPrim2, fromHE2, thisPrim2, thisConnectIndex2] = potentialBoundaryCorrespondences.at(selectedStart).at(1);
                thisPrim->connections.at(thisConnectIndex) = thisPrim2;
                thisPrim2->connections.at(thisConnectIndex2) = thisPrim;
                // return true;
                return result;
            }

            //  create new primitives corresponding to all elements of "to"
            //      TODO maybe should have a function that does this separately? since also doing same thing in the from.primitives.size() == 0 case
            //      I think: iterate over to.primitives; make new uPtr<Primitive> w/ same halfEdges
            //      store map of pointers to elements of to.primitives and the new primitives
            //      iterate over to.primitives again and set connections appropriately
            //          then maybe here also do the updating of nullptr ones based on potentialBoundaryCorrespondences
            //          -> update connections in/to new primitives using potentialBoundaryCorrespondences and boundary strings
            std::map<Primitive*, Primitive*> toToThisPrimitives; // TODO make a name that is less horrible maybe
            std::vector<uPtr<Primitive>> toAdd;
            for (const uPtr<Primitive>& p : to.primitives) {
                uPtr<Primitive> newPrim = mkU<Primitive>();
                newPrim->halfEdges = p->halfEdges; // TODO can't recall if this copies right
                toToThisPrimitives.insert({p.get(), newPrim.get()});
                toAdd.push_back(std::move(newPrim));

            }

            for (Primitive* p : potentialMatches.at(selectedStart)) {
                // check for case of connection from newly created primtive being to another newly created primitive which it isn't connected to in "to"
                // WAIT actually is this just a case we should skip? seems like the rule being applied will not actually change graph structure.
                //  IS THERE SOME CASE WHERE WE SHOULD DO THIS? TODO
                //  if not can just discount any case where thisPrim is in matches, right?
                // for (auto& [fromPrim, fromHE, thisPrim, thisConnectIndex] : potentialBoundaryCorrespondences.at(selectedStart)) {
                //     if (thisPrim == p) {
                //         for (unsigned int i = 0; i < from.boundaryString.size(); ++i) {
                //             const BoundaryHE* bFrom = std::get_if<BoundaryHE>(&from.boundaryString.at((fromStringStart + i) % from.boundaryString.size()));
                //             if (bFrom != nullptr && from.primitives.at(bFrom->primIndex).get() == fromPrim) {
                //                 const BoundaryHE* bTo = std::get_if<BoundaryHE>(&to.boundaryString.at(i));
                //                 Primitive* bToPrim = to.primitives.at(bTo->primIndex).get();
                //                 // for (unsigned int j = 0; j < bToPrim->halfEdges.size(); ++j) {
                //                 //     if (bToPrim->halfEdges.at(j) == fromHE) {
                //                 //         // thisPrim = toT
                //                 //     }
                //                 // }
                //                 // return false;
                //                 // thisPrim = toToThisPrimitives.at(bToPrim);
                //                 // break;
                //             }
                //         }

                //     }
                // }


                std::erase_if(result.primitives, [&p](uPtr<Primitive>& curPrim) {
                    return curPrim.get() == p;
                });
            }

            // TODO I think there's a one-line way to do this?
            for (uPtr<Primitive>& p : toAdd) {
                result.primitives.push_back(std::move(p));
            }



            for (const uPtr<Primitive>& p : to.primitives) {
                Primitive* newPrim = toToThisPrimitives.at(p.get());
                // for (Primitive* connect : p->connections) {
                for (unsigned int j = 0; j < p->connections.size(); ++j) {
                    Primitive* connect = p->connections.at(j);
                    if (connect == nullptr) {
                        HalfEdgeGraph he = p->halfEdges.at(j);
                        // boundary of "to" graph
                        // using boundary string, find matching primitive in "from" graph
                        // then potentialBoundaryCorrespondences tells which primitive in "this" to connect with
                        // TODO I'm not sure actually if potentialBoundaryCorrespondences[selectedStart] can include primitives which we're erasing through this algorithm
                        //  I don't THINK it should be able to but not 100% sure
                        //  If it is an issue then can add check of nullptr in the pair comparisons above? but not sure if it is a case that can happen
#if THREE_DIMENSIONAL
                        for (unsigned int i = 0; i < from.boundaryGraphElements.size(); ++i) {
                            BoundaryElem* bTo = to.boundaryGraphElements.at(i).get();
                            if (to.primitives.at(bTo->he.primIndex) == p) {
                                BoundaryElem* bFrom = from.boundaryGraphElements.at(i).get();
                                Primitive* bFromPrim = from.primitives.at(bFrom->he.primIndex).get();
                                bool matched = false;
                                for (const auto& [fromPrim, fromHE, thisPrim, thisConnectIndex] : potentialBoundaryCorrespondences.at(selectedStart)) {
                                    if (bFromPrim == fromPrim && he == fromHE) {

                                        newPrim->connections.push_back(thisPrim);
                                        thisPrim->connections.at(thisConnectIndex) = newPrim;

                                        // TODO fix bugs: sometimes thisPrim->connections is empty. this is where it crashes but not sure yet where the issue starts
                                        //  TODO uncertain if that still happens with this approach
                                        matched = true;
                                        break;
                                    }
                                }
                                if (matched) {
                                    break;
                                }

                            }
                        }
#else

                        for (unsigned int i = 0; i < from.boundaryString.size(); ++i) {
                            const BoundaryHE* bTo = std::get_if<BoundaryHE>(&to.boundaryString.at(i));
                            if (bTo != nullptr && to.primitives.at(bTo->primIndex) == p) {
                                const BoundaryHE* bFrom = std::get_if<BoundaryHE>(&from.boundaryString.at((fromStringStart + i) % from.boundaryString.size()));
                                // to.primitives.at(bTo->primIndex);
                                Primitive* bFromPrim = from.primitives.at(bFrom->primIndex).get();
                                bool matched = false;
                                for (const auto& [fromPrim, fromHE, thisPrim, thisConnectIndex] : potentialBoundaryCorrespondences.at(selectedStart)) {
                                    if (bFromPrim == fromPrim && he == fromHE) {
                                        // if (toToThisPrimitives.contains(thisPrim)) {

                                        // }

                                        newPrim->connections.push_back(thisPrim);
                                        thisPrim->connections.at(thisConnectIndex) = newPrim;

                                        // TODO fix bugs: sometimes thisPrim->connections is empty. this is where it crashes but not sure yet where the issue starts
                                        matched = true;
                                        break;
                                    }
                                }
                                if (matched) {
                                    break;
                                }
                            }
                        }
#endif

                    } else {
                        newPrim->connections.push_back(toToThisPrimitives.at(connect));
                    }
                }
                // std::cout << "NOR " << p->connections.size() << " " << newPrim->connections.size() << std::endl;
            }

            // TODO this might have some issue that doesn't appear on my computer? Will test later

            // return true;
            return result;
        }
    } catch (const std::exception& e) {
        std::cerr << "rule application fail" << std::endl;
        std::cerr << e.what() << std::endl;
        // TODO try to fix the bug causing failures but for now can just skip those fail cases
    }
    return {};
    // return false;
}

std::optional<Graph> Graph::applyRandomReplacementRule(const std::vector<std::pair<Graph, Graph>>& grammar, bool bidirectional, bool skipStarters)
{
    // TODO check this works



    std::vector<std::pair< Graph, Graph >> grammarShuffled;
    for (const auto& rule : grammar) {
        grammarShuffled.push_back(rule);
        // TODO I think can just do grammarShuffled = gramamr now
    }
    std::random_device rd; // TODO might need to have declared outside but idk
    std::mt19937 g(rd());
    std::shuffle(grammarShuffled.begin(), grammarShuffled.end(), g);

    std::vector<bool> order = {true};
    if (bidirectional) {
        order.push_back(false);
        std::uniform_real_distribution<float> flipChance(0.f, 1.f);

        //dunno if right or a good idea but make other direction less likely?
        if (flipChance(g) > 0.5) {
            std::shuffle(order.begin(), order.end(), g);
        }
        // TODO does g need to be different?
    }

    // TODO I think this shouldn't have bias in terms of which end up picked but not 100% sure
    // ones that can apply more will be picked more frequently if others can't always apply but I think that's fine

    // shuffles rules and tries in order -> picks random rules without repetition until finding a valid rule to apply
    std::optional<Graph> result;
    for (const std::pair< Graph, Graph>& rule : grammarShuffled) {
        if (skipStarters && rule.first.primitives.size() == 0) {
            continue;
        }
        for (bool LtoR : order) {
            if (LtoR) {
                // left to right
                result = applyReplacementRule(rule.first, rule.second);
                // if (result.has_value()) {
                    // std::cout << "success" << std::endl;
                    // return result;
                // }
                // std::cout << "fail" << std::endl;
            } else {
                // right to left
                result = applyReplacementRule(rule.second, rule.first);
                // if (result.has_value()) {
                    // std::cout << "success" << std::endl;
                    // return true;
                    // return result;
                // }
                // std::cout << "fail" << std::endl;
            }
            if (result.has_value()
                && !(result->isIsomorphicTo(*this))) {
                //attempt to find positions
                // skipping rule applications which don't do anything
                std::map<unsigned int,glm::vec3> cachedPositionMap;
                std::vector<unsigned int> unfreedIndices;
                for (unsigned int i = 0; i < result->primitives.size(); ++i) {
                    Primitive* p = result->primitives.at(i).get();
                    if (p->cachedPos.has_value()) {
                        cachedPositionMap.insert({i,p->cachedPos.value()});
                        unfreedIndices.push_back(i);
                    }
                }

                double posError = -1.0;
                const double minError = 1.0e-2;

                std::random_device rd;
                std::mt19937 g(rd());
                std::shuffle(unfreedIndices.begin(), unfreedIndices.end(), g);
                std::vector<glm::vec3> potentialPositions;
                double maxLen = 10;
                // while (maxLen <= 40) {
                    while (posError < 0 || posError > minError) {


                        // TODO for dimensional limits: probably should try to make random distribution bias towards 0? so can hit whole range but tends to spread out less
                        // TODO maybe maxlen should be separate from the max that it samples from? so can be longer but doesn't deliberately make longer sorta thing
                        potentialPositions = result.value().samplePositions(cachedPositionMap, posError, 1, maxLen, -30, 30);

                        if (posError >= 0 && posError <= minError) {
                            // successful
                            break;
                        }
                        if (unfreedIndices.size() > 0) {
                            // free a vertex position

                            unsigned int toFree = unfreedIndices.at(0);
                            cachedPositionMap.erase(toFree);
                            unfreedIndices.erase(unfreedIndices.begin());
                            // maxLen *= 2;
                        } else {
                            // give up on this rule
                            break;
                        }
                    }
                //     maxLen *= 2.0;
                // }
                    // posError = 0;
                if (posError >= 0.0 && posError <= minError) {
                    // update cached positions

                    for (unsigned int i = 0; i < potentialPositions.size(); ++i) {
                        result.value().primitives.at(i)->cachedPos = potentialPositions.at(i);
                    }

                    return result;

                }

                result = {};
            }
        }
    }
    return {};
}


std::vector<Graph> Graph::branchGlue(const Graph &other) const
{
    // TODO could define splice and use that plus loop glue it sounds like but might be easier doing separate


    std::vector<Graph> results;

#if THREE_DIMENSIONAL
    // if (other.primitives.size() == 1) {

    // )
    // std::vector<BoundaryElem*>
    for (const uPtr<BoundaryElem>& otherElem : other.boundaryGraphElements) {
        for (unsigned int k = 0; k < this->boundaryGraphElements.size(); ++k) {
            const uPtr<BoundaryElem>& thisElem = this->boundaryGraphElements.at(k);
        // for (const uPtr<BoundaryElem>& thisElem : this->boundaryGraphElements) {
            if (thisElem->he.h == -otherElem->he.h) {
                // mirrored pair of half edges
                // try to glue here

                Graph gluedGraph(*this);
                // for now assuming branch gluing with other having a single primitive since other cases irrelevant anyway
                uPtr<Primitive> newPrim = mkU<Primitive>(*other.primitives.at(otherElem->he.primIndex));
                // uPtr<Primitive> newPrim = mkU<Primitive>(*other.primitives.at(0));
                // TODO check if default copy constructor works? might not update ID but that's fine
                Primitive* thisPrimCopy = gluedGraph.primitives.at(thisElem->he.primIndex).get();
                for (unsigned int i = 0; i < newPrim->halfEdges.size(); ++i) {
                    if (newPrim->halfEdges.at(i) == otherElem->he.h) {
                        // newPrim->connections.at(otherElem->he.primIndex) = gluedGraph.primitives.at()
                        for (unsigned int j = 0; j < thisPrimCopy->halfEdges.size(); ++j) {
                            if (thisPrimCopy->halfEdges.at(j) == thisElem->he.h) {
                                newPrim->connections.at(i) = thisPrimCopy;
                                thisPrimCopy->connections.at(j) = newPrim.get();
                                break;
                            }
                        }
                        break;
                    }
                }

                // update boundary graph
                BoundaryElem* tE = gluedGraph.boundaryGraphElements.at(k).get();
                BoundaryElem* oE = nullptr;
                unsigned int oE_index;
                // copy boundary from other
                // std::map<BoundaryElem*,BoundaryElem*> corr; // I think can just use next index since assuming other is a single primitive, should(?) be more efficient
                unsigned int thisBoundarySize = gluedGraph.boundaryGraphElements.size();
                unsigned int otherBoundarySize = other.boundaryGraphElements.size();
                for (unsigned int i = 0; i < otherBoundarySize; ++i) {
                    uPtr<BoundaryElem> newElem = mkU<BoundaryElem>(*other.boundaryGraphElements.at(i));

                    if (otherElem == other.boundaryGraphElements.at(i)) {
                        oE = newElem.get();
                        oE_index = i + thisBoundarySize;
                    }
                    newElem->he.primIndex += this->primitives.size();
                    gluedGraph.boundaryGraphElements.push_back(std::move(newElem));
                }
                for (unsigned int i = 0; i < otherBoundarySize; ++i) {
                    BoundaryElem* target = gluedGraph.boundaryGraphElements.at(thisBoundarySize + i).get();
                    BoundaryElem* nextElem = gluedGraph.boundaryGraphElements.at(thisBoundarySize + (i + 1) % otherBoundarySize).get();
                    target->next = nextElem;
                    nextElem->prev = target;
                }
                // rotate two pairs of edges, remove island

                if (oE == nullptr) {
                    continue;
                }
                tE->prev->next = oE->next;
                oE->next->prev = tE->prev;
                tE->next->prev = oE->prev;
                oE->prev->next = tE->next;

                gluedGraph.boundaryGraphElements.erase(gluedGraph.boundaryGraphElements.begin() + oE_index);
                gluedGraph.boundaryGraphElements.erase(gluedGraph.boundaryGraphElements.begin() + k);


                gluedGraph.primitives.push_back(std::move(newPrim));
                gluedGraph.sortBoundaryGraphElements();

                results.push_back(std::move(gluedGraph));
            }
        }
    }

#else

    for (unsigned int j = 0; j < other.boundaryString.size(); ++j) {
        const BoundaryHE* otherStart = std::get_if<BoundaryHE>(&other.boundaryString.at(j));
        if (otherStart != nullptr) {
            std::vector<std::variant<BoundaryHE,Turn>> substring;
            for (unsigned int i = 1; i < other.boundaryString.size(); ++i) {
                std::variant<BoundaryHE,Turn> val = other.boundaryString.at((j + i) % other.boundaryString.size());
                BoundaryHE* bVal = std::get_if<BoundaryHE>(&val);
                if (bVal != nullptr) {
                    bVal->primIndex += this->primitives.size();
                    // TODO not sure bVal actually mutates val
                }
                substring.push_back(val);

            }
            for (unsigned int i = 0; i < boundaryString.size(); ++i) {
                const BoundaryHE* subStart = std::get_if<BoundaryHE>(&boundaryString.at(i));
                if (subStart != nullptr && subStart->h == -(otherStart->h)) {

                    Graph gluedGraph = *this;
                    std::map<Primitive*, Primitive*> otherToNewPrimitive;
                    for (const uPtr<Primitive>& p : other.primitives) {
                        uPtr<Primitive> newPrim = mkU<Primitive>();
                        newPrim->halfEdges = p->halfEdges; // TODO can't recall if this copies right
                        otherToNewPrimitive.insert({p.get(), newPrim.get()});
                        gluedGraph.primitives.push_back(std::move(newPrim));
                    }


                    for (const uPtr<Primitive>& p : other.primitives) {
                        Primitive* newPrim = otherToNewPrimitive.at(p.get());
                        for (unsigned int j = 0; j < p->connections.size(); ++j) {
                            Primitive* connect = p->connections.at(j);
                            if (connect == nullptr) {
                                newPrim->connections.push_back(nullptr);
                            } else {
                                newPrim->connections.push_back(otherToNewPrimitive.at(connect));
                            }
                        }
                    }
                    Primitive* p1 = gluedGraph.primitives.at(subStart->primIndex).get();
                    Primitive* p2 = gluedGraph.primitives.at(otherStart->primIndex + this->primitives.size()).get();

                    for (unsigned int j = 0; j < p1->connections.size(); ++j) {
                        if (p1->halfEdges.at(j) == subStart->h) {
                            p1->connections.at(j) = p2;
                            break;
                        }
                    }
                    for (unsigned int j = 0; j < p2->connections.size(); ++j) {
                        if (p2->halfEdges.at(j) == otherStart->h) {
                            p2->connections.at(j) = p1;
                            break;
                        }
                    }
                    gluedGraph.boundaryString.at(i) = Turn::negative;

                    unsigned int startIndex = i;
                    if (subStart->h.angle < 0) {
                        // subStart == a-bar
                        // a-bar -> negative,B
                        // insert after i
                        ++startIndex;
                    }
                    // else
                    // subStart == a
                    // a -> B,negative
                    //insert before i
                    gluedGraph.boundaryString.insert(gluedGraph.boundaryString.begin() + startIndex, substring.begin(), substring.end());

                    gluedGraph.cancelTurns();

                    results.push_back(gluedGraph);
                }
            }
        }
    }
#endif
    return results;
}


std::vector<Graph> Graph::loopGlue() const
{
    std::vector<Graph> results;
#if THREE_DIMENSIONAL
    // TODO 4/16
    // I THINK what I need to do:
    // loop gluing logic
    //  think done. note obviously has impossible geometry cases allowed rn--but I think fine? just makes algorithm slower to find a valid solution?
    // same boundary check
    //  I think done
    // rule application boundary correspondence
    //  actually should just need to swap from std::variant to .he mostly and otherwise similar?
    //  I think done!
    // set normals in graph constructor
    //  I think done


    // position sampling
    // figure out what to do about not getting starting rule. it's because always before the last loop glue for the starting rule it produces a rule of the two facing edges to a single edge b/c we don't check the turns



    // Note for now letting connect to any other mirrored HE regardless of location in graph (no turn checking)

    // TODO probably better unordered_map
    std::map<HalfEdgeGraph, std::vector<unsigned int>> elemsByHE;
    for (unsigned int i = 0; i < this->boundaryGraphElements.size(); ++i) {
        BoundaryElem* elem = this->boundaryGraphElements.at(i).get();


        HalfEdgeGraph invGraph = -elem->he.h;

        if (elemsByHE.contains(invGraph)) {
            std::vector<unsigned int>& matches = elemsByHE.at(invGraph);
            // for (unsigned int j = 0; j < matches.size(); ++j) {
            for (unsigned int j : matches) {
                // matching HEs

                // skipping gluing HE on same primitive together since that'll always fail to make valid geometry

                // TODO maybe some degree of turn checking is possible and would help even if not going fully precise?
                //  like at least track within plane if a path has turned % 360
                if (elem->he.primIndex != this->boundaryGraphElements.at(j)->he.primIndex) {
                    Graph gluedGraph(*this);
                    BoundaryElem* tE = gluedGraph.boundaryGraphElements.at(i).get();
                    BoundaryElem* oE = gluedGraph.boundaryGraphElements.at(j).get();

                    tE->prev->next = oE->next;
                    oE->next->prev = tE->prev;
                    tE->next->prev = oE->prev;
                    oE->prev->next = tE->next;

                    for (unsigned int k = 0; k < gluedGraph.primitives.at(tE->he.primIndex)->halfEdges.size(); ++k) {
                        if (gluedGraph.primitives.at(tE->he.primIndex)->halfEdges.at(k) == tE->he.h) {
                            gluedGraph.primitives.at(tE->he.primIndex)->connections.at(k) = gluedGraph.primitives.at(oE->he.primIndex).get();
                            break;
                        }
                    }

                    for (unsigned int k = 0; k < gluedGraph.primitives.at(oE->he.primIndex)->halfEdges.size(); ++k) {
                        if (gluedGraph.primitives.at(oE->he.primIndex)->halfEdges.at(k) == oE->he.h) {
                            gluedGraph.primitives.at(oE->he.primIndex)->connections.at(k) = gluedGraph.primitives.at(tE->he.primIndex).get();
                            break;
                        }
                    }


                    gluedGraph.boundaryGraphElements.erase(gluedGraph.boundaryGraphElements.begin() + i);
                    gluedGraph.boundaryGraphElements.erase(gluedGraph.boundaryGraphElements.begin() + j);

                    gluedGraph.sortBoundaryGraphElements();

                    results.push_back(std::move(gluedGraph));

                }

            }
        }

        elemsByHE[elem->he.h].push_back(i);
    }




#else
    for (unsigned int i = 0; i < boundaryString.size(); ++i) {
        // std::vector<std::variant<BoundaryHE, Turn>> boundaryString;
        const BoundaryHE* subStart = std::get_if<BoundaryHE>(&boundaryString.at(i));
        if (subStart != nullptr) {
            // loop glue cases:
            //  subStart = a, next term is a-bar
            //  subStart = a-bar, next 3 terms are {negative, a, positive}
            if (subStart->h.angle < 0.f) {
                // subStart == a-bar

                if (boundaryString.size() >= 4) {
                    std::vector<std::variant<BoundaryHE,Turn>> substring = {
                        boundaryString.at(i),
                        boundaryString.at((i + 1) % boundaryString.size()),
                        boundaryString.at((i + 2) % boundaryString.size()),
                        boundaryString.at((i + 3) % boundaryString.size()),
                    };

                    BoundaryHE* potentialA = std::get_if<BoundaryHE>(&substring.at(2));
                    Turn* t1 = std::get_if<Turn>(&substring.at(1));
                    Turn* t2 = std::get_if<Turn>(&substring.at(3));
                    if (potentialA != nullptr && t1 != nullptr && t2 != nullptr &&
                        *t1 == Turn::negative && *t2 == Turn::positive && potentialA->h == -(subStart->h)) {
                        Graph gluedGraph = *this;

                        Primitive* p1 = gluedGraph.primitives.at(subStart->primIndex).get();
                        Primitive* p2 = gluedGraph.primitives.at(potentialA->primIndex).get();

                        for (unsigned int j = 0; j < p1->connections.size(); ++j) {
                            if (p1->halfEdges.at(j) == subStart->h) {
                                p1->connections.at(j) = p2;
                                break;
                            }
                        }
                        for (unsigned int j = 0; j < p2->connections.size(); ++j) {
                            if (p2->halfEdges.at(j) == potentialA->h) {
                                p2->connections.at(j) = p1;
                                break;
                            }
                        }

                        if (i < boundaryString.size() - 3) {
                            gluedGraph.boundaryString.erase(gluedGraph.boundaryString.begin() + i, gluedGraph.boundaryString.begin() + i + 4);
                        } else {
                            unsigned int endIndex = (i + 3) % boundaryString.size() + 1;
                            gluedGraph.boundaryString.erase(gluedGraph.boundaryString.begin() + i, gluedGraph.boundaryString.end());
                            gluedGraph.boundaryString.erase(gluedGraph.boundaryString.begin(), gluedGraph.boundaryString.begin() + endIndex);
                        }

                        gluedGraph.cancelTurns();

                        results.push_back(gluedGraph);
                    }
                }
            } else {
                // subStart == a
                unsigned int nextIndex = (i + 1) % boundaryString.size();
                const BoundaryHE* aBar = std::get_if<BoundaryHE>(&boundaryString.at(nextIndex));
                if (aBar != nullptr && aBar->h == -(subStart->h) ) {

                    Graph gluedGraph = *this;

                    Primitive* p1 = gluedGraph.primitives.at(subStart->primIndex).get();
                    Primitive* p2 = gluedGraph.primitives.at(aBar->primIndex).get();


                    for (unsigned int j = 0; j < p1->connections.size(); ++j) {
                        if (p1->halfEdges.at(j) == subStart->h) {
                            p1->connections.at(j) = p2;
                            break;
                        }
                    }
                    for (unsigned int j = 0; j < p2->connections.size(); ++j) {
                        if (p2->halfEdges.at(j) == aBar->h) {
                            p2->connections.at(j) = p1;
                            break;
                        }
                    }

                    // TODO check I'm doing erasing right
                    if (nextIndex > i) {
                        gluedGraph.boundaryString.erase(gluedGraph.boundaryString.begin() + nextIndex);
                        gluedGraph.boundaryString.erase(gluedGraph.boundaryString.begin() + i);
                    } else {
                        gluedGraph.boundaryString.erase(gluedGraph.boundaryString.begin() + i);
                        gluedGraph.boundaryString.erase(gluedGraph.boundaryString.begin() + nextIndex);
                    }

                    gluedGraph.cancelTurns();

                    results.push_back(std::move(gluedGraph)); // std::move shouldn't be needed but I think more efficient assuming it works? since don't need it anymore anyway

                }
            }


        }
    }
#endif
    return results;
}

void Graph::cancelTurns()
{
    if (boundaryString.size() < 2) {
        return;
    }
    for (unsigned int i = 0; i < boundaryString.size(); ++i) {
        const Turn* first = std::get_if<Turn>(&boundaryString.at(i));
        const Turn* second = std::get_if<Turn>(&boundaryString.at((i + 1) % boundaryString.size()));
        if (first != nullptr && second != nullptr &&
            ((*first == Turn::positive && *second == Turn::negative) ||
             (*first == Turn::negative && *second == Turn::positive))) {
            if (i < boundaryString.size() - 1) {
                this->boundaryString.erase(this->boundaryString.begin() + i + 1);
                this->boundaryString.erase(this->boundaryString.begin() + i);
            } else {
                this->boundaryString.erase(this->boundaryString.begin() + i);
                this->boundaryString.erase(this->boundaryString.begin());
            }
            --i;
        }
    }
}

bool Graph::hasSubgraph(const Graph &other) const
{
    std::vector<Primitive*> potentialStarts;
    std::vector<std::vector<Primitive*>> potentialMatches;

    std::vector<Primitive*> primsToMatch;
    for (const uPtr<Primitive>& p : other.primitives) {
        primsToMatch.push_back(p.get());
    }
    for (unsigned int i = 0; i < this->primitives.size(); ++i) {
        Primitive* curPrim = this->primitives.at(i).get();
        std::vector<Primitive*> unmatchedPrims = primsToMatch;

        std::set<std::pair<Primitive*,Primitive*>> encountered;


        std::vector<std::tuple<Primitive*, HalfEdgeGraph, Primitive*, unsigned int>> boundaryCorrespondences;

        bool isMatched = true;
        while (isMatched && unmatchedPrims.size() > 0) {
            isMatched = checkMatch(&encountered, unmatchedPrims.at(0), curPrim, &boundaryCorrespondences);
            // TODO maybe make version that doesn't have the extra stuff tracked, just to reduce overhead
            for (const std::pair<Primitive*,Primitive*>& pair : encountered) {
                std::erase_if(unmatchedPrims, [&pair](Primitive* curPrim) {
                    return curPrim == pair.first;
                });
            }
        }

        if(isMatched) {

            return true;
        }

    }
    return false;
}


bool checkMatchExact(std::set<std::pair<Primitive*,Primitive*>>* encountered, Primitive* fromPrim, Primitive* thisPrim) {
    if (thisPrim == nullptr || fromPrim == nullptr) {
        // TODO maybe should make return true if both nullptr?
        return false;
    }
    if (fromPrim->sameHEs(*thisPrim)) {
        encountered->insert({fromPrim, thisPrim});
        for (unsigned int j = 0; j < fromPrim->connections.size(); ++j) {
            Primitive* from2 = fromPrim->connections.at(j);
            for (unsigned int k = 0; k < thisPrim->connections.size(); ++k) {
                if (fromPrim->halfEdges.at(j) == thisPrim->halfEdges.at(k)) {
                    Primitive* this2 = thisPrim->connections.at(k);
                    if (from2 == nullptr || this2 == nullptr) {
                        if (from2 != this2) {
                            // if either null both must be to be isomorphic
                            return false;
                        }
                    } else {
                        auto foundFrom = std::find_if(encountered->begin(), encountered->end(), [&from2](std::pair<Primitive*,Primitive*> a) {
                            return a.first == from2;
                        });
                        // TODO can simplify these checks some I think
                        if (foundFrom == encountered->end()) {
                            auto foundTo = std::find_if(encountered->begin(), encountered->end(), [&this2](std::pair<Primitive*,Primitive*> a) {
                                return a.second == this2;
                            });
                            if (foundTo == encountered->end()) {
                                // not already encountered

                                if (!checkMatchExact(encountered, from2, this2)) {
                                    return false;
                                }
                            } else {
                                // mismatch
                                return false;
                            }
                        } else {
                            // already encountered, check if same pair
                            if (this2 != foundFrom->second) {
                                // mismatch
                                return false;
                            }
                        }
                    }
                    break;
                }
            }
        }
        return true;
    } else {
        return false;
    }
}


bool Graph::isIsomorphicTo(const Graph &other) const
{

    // TODO I should probably rewrite this to reduce redundancy since performance is an issue

    if (other.primitives.size() != this->primitives.size()) {
        return false;
    }
    std::vector<Primitive*> primsToMatch;
    for (const uPtr<Primitive>& p : other.primitives) {
        primsToMatch.push_back(p.get());
    }
    for (unsigned int i = 0; i < this->primitives.size(); ++i) {
        Primitive* curPrim = this->primitives.at(i).get();
        std::vector<Primitive*> unmatchedPrims = primsToMatch;

        std::set<std::pair<Primitive*,Primitive*>> encountered;


        // std::vector<std::tuple<Primitive*, HalfEdgeGraph, Primitive*, unsigned int>> boundaryCorrespondences;

        bool isMatched = true;
        while (isMatched && unmatchedPrims.size() > 0) {
            isMatched = checkMatchExact(&encountered, unmatchedPrims.at(0), curPrim);
            for (const std::pair<Primitive*,Primitive*>& pair : encountered) {
                std::erase_if(unmatchedPrims, [&pair](Primitive* curPrim) {
                    return curPrim == pair.first;
                });
            }
        }

        if(isMatched) {

            return true;
        }

    }
    return false;

}

bool Graph::sameBoundaryString(const Graph &other) const
{

#if THREE_DIMENSIONAL
    // TODO 4/16
    if (this->boundaryGraphElements.size() != other.boundaryGraphElements.size()) {
        return false;
    }

    if (this->boundaryGraphElements.size() == 0) {
        return true;
    }

    return this->getBoundaryGeneric() == other.getBoundaryGeneric();

    // could do like the other graph isomoprhism check
    //  but I think another potential way to do:
    //  sort boundaryGraphElements in some consistent way
    //   probably should do outside this function and reuse to make faster
    // NOTE rn since I have a special exception for same-primitive loop gluing (which should be part of a wider exception but the rest of it not implemented yet) can match a boundary string with different gluing possibilities rn but I think that shouldn't matter?
    // std::vector<BoundaryElem*> boundaryElemsToMatch;
    // for (const uPtr<BoundaryElem>& ePtr : this->boundaryGraphElements) {
    //     boundaryElemsToMatch.push_back(ePtr.get());
    // }


    auto cmp = [](const BoundaryElem* a, const BoundaryElem* b) {
        return a->he.h < b->he.h;
    };
    std::set<BoundaryElem*, decltype(cmp)> boundaryElemsToAdd(cmp);
    std::set<BoundaryElem*, decltype(cmp)> boundaryElemsToAddOther(cmp);

    auto loopCmp = [](const std::vector<BoundaryElem*>& a, const std::vector<BoundaryElem*>& b) {
        if (a.size() == b.size()) {
            for (unsigned int i = 0; i < a.size(); ++i) {
                if (a.at(i)->he.h != b.at(i)->he.h) {
                    return a.at(i)->he.h < b.at(i)->he.h;
                }
            }
        }
        return a.size() < b.size();
    };
    std::set<std::vector<BoundaryElem*>, decltype(loopCmp)> loops;
    std::set<std::vector<BoundaryElem*>, decltype(loopCmp)> loopsOther;

    for (const uPtr<BoundaryElem>& ePtr : this->boundaryGraphElements) {
        boundaryElemsToAdd.insert(ePtr.get());
    }

    while (boundaryElemsToAdd.size() > 0) {
        BoundaryElem* loopStart = *boundaryElemsToAdd.begin();
        std::vector<BoundaryElem*> newLoop;
        // newLoop.push_back(loopStart);
        // boundaryElemsToAdd.erase(loopStart);
        BoundaryElem* curElem = loopStart;
        do {
            newLoop.push_back(curElem);
            boundaryElemsToAdd.erase(curElem);
            curElem = curElem->next;
        } while (curElem != loopStart);
        loops.insert(std::move(newLoop));
    }

    for (const uPtr<BoundaryElem>& ePtr : other.boundaryGraphElements) {
        boundaryElemsToAddOther.insert(ePtr.get());
    }

    while (boundaryElemsToAddOther.size() > 0) {
        BoundaryElem* loopStart = *boundaryElemsToAddOther.begin();
        std::vector<BoundaryElem*> newLoop;
        BoundaryElem* curElem = loopStart;
        do {
            newLoop.push_back(curElem);
            boundaryElemsToAddOther.erase(curElem);
            curElem = curElem->next;
        } while (curElem != loopStart);
        loopsOther.insert(std::move(newLoop));
    }

    if (loops.size() != loopsOther.size()) {
        return false;
    }

    auto iter1 = loops.begin();
    const auto iter1End = loops.end();
    auto iter2 = loopsOther.begin();
    while (iter1 != iter1End) {
        const std::vector<BoundaryElem*>& v1 = *iter1;
        const std::vector<BoundaryElem*>& v2 = *iter2;
        if (v1.size() != v2.size()) {
            return false;
        }
        for (unsigned int i = 0; i < v1.size(); ++i) {
            if (v1.at(i)->he.h != v2.at(i)->he.h) {
                return false;
            }
        }

        std::advance(iter1, 1);
        std::advance(iter2, 1);
    }


    return true;
#else
    if (this->boundaryString.size() != other.boundaryString.size()) {
        return false;
    }
    unsigned int length = this->boundaryString.size();
    for (unsigned int otherStart = 0; otherStart < length; ++otherStart) {
        bool isMatched = true;
        for (unsigned int i = 0; i < length; ++i) {
            // TODO I think can just compare std::variant with ==/!=? make sure works
            if (this->boundaryString.at(i) != other.boundaryString.at((otherStart + i) % length)) {
                isMatched = false;
                break;
            }
        }
        if (isMatched) {
            return true;
        }
    }
    return false;
#endif
}

void Graph::sortBoundaryGraphElements() {
    auto cmp = [](const BoundaryElem* a, const BoundaryElem* b) {
        return a->he.h < b->he.h;
    };
    std::set<BoundaryElem*, decltype(cmp)> boundaryElemsToAdd(cmp);

    auto loopCmp = [](const std::vector<BoundaryElem*>& a, const std::vector<BoundaryElem*>& b) {
        if (a.size() == b.size()) {
            for (unsigned int i = 0; i < a.size(); ++i) {
                if (a.at(i)->he.h != b.at(i)->he.h) {
                    return a.at(i)->he.h < b.at(i)->he.h;
                }
            }
        }
        return a.size() < b.size();
    };
    std::set<std::vector<BoundaryElem*>, decltype(loopCmp)> loops;

    std::map<BoundaryElem*, unsigned int> elemIndices;
    std::vector<unsigned int> indexOrder;

    for (unsigned int i = 0; i < this->boundaryGraphElements.size(); ++i) {
        const uPtr<BoundaryElem>& ePtr = this->boundaryGraphElements.at(i);
        // for (const uPtr<BoundaryElem>& ePtr : this->boundaryGraphElements) {
        boundaryElemsToAdd.insert(ePtr.get());
        elemIndices.insert({ePtr.get(), i});
    }

    while (boundaryElemsToAdd.size() > 0) {
        BoundaryElem* loopStart = *boundaryElemsToAdd.begin();
        std::vector<BoundaryElem*> newLoop;
        // newLoop.push_back(loopStart);
        // boundaryElemsToAdd.erase(loopStart);
        BoundaryElem* curElem = loopStart;
        do {
            newLoop.push_back(curElem);
            boundaryElemsToAdd.erase(curElem);
            curElem = curElem->next;
        } while (curElem != loopStart);
        loops.insert(std::move(newLoop));
    }

    this->boundaryGeneric.clear();
    for (const std::vector<BoundaryElem*>& loop : loops) {
        std::vector<HalfEdgeGraph> heLoop;
        for (BoundaryElem* ptr : loop) {
            indexOrder.push_back(elemIndices.at(ptr));
            heLoop.push_back(ptr->he.h);
        }
        this->boundaryGeneric.push_back(heLoop);
    }

    std::vector<uPtr<BoundaryElem>> originalOrder = std::move(this->boundaryGraphElements);

    this->boundaryGraphElements.clear();

    for (unsigned int i : indexOrder) {
        this->boundaryGraphElements.push_back(std::move(originalOrder.at(i)));
    }


    // for (unsigned int i = 0; i < this->boundaryGraphElements.size(); ) {
    //     std::vector<HalfEdgeGraph> curLoop;
    //     BoundaryElem* startElem = this->boundaryGraphElements.at(i).get();
    //     BoundaryElem* curElem = startElem;
    //     do {
    //         curLoop.push_back(curElem->he.h);
    //         curElem = curElem->next;
    //         ++i;
    //     } while (curElem != startElem);
    //     this->boundaryGeneric.push_back(std::move(curLoop));
    // }



}

const std::vector<std::vector<HalfEdgeGraph> > &Graph::getBoundaryGeneric() const
{
    return this->boundaryGeneric;
    // std::vector<std::vector<HalfEdgeGraph>> loops;
    // for (unsigned int i = 0; i < this->boundaryGraphElements.size(); ) {
    //     std::vector<HalfEdgeGraph> curLoop;
    //     BoundaryElem* startElem = this->boundaryGraphElements.at(i).get();
    //     BoundaryElem* curElem = startElem;
    //     do {
    //         curLoop.push_back(curElem->he.h);
    //         curElem = curElem->next;
    //         ++i;
    //     } while (curElem != startElem);
    //     loops.push_back(std::move(curLoop));
    // }
    // return loops;
}

// #include <unordered_map>
#define outputProgress 0
#define INCLUDE_SINGLE_EDGE 0
// TODO I'm not sure if all of the boundary strings are maintained totally correctly or not. looking at graphs used in rules there's some that intuitively seem like they don't match to me but I'm not sure if they actually do or not. hard to tell what counts as what turn visually personally
//  nevermind I think the case I was worried about is right; did it out on paper and seems correct
std::vector<std::pair<Graph, Graph>> Graph::generateRules(const std::vector<Primitive *> &primitives,  std::vector<glm::vec3> faceColors, unsigned int maxSteps)
{
    std::vector<std::pair<Graph,Graph>> result;

    // std::vector<std::vector<Graph>> hierarchy;
    // std::vector<Graph> hierarchy;

    std::map<std::vector<std::vector<HalfEdgeGraph>>, Graph> graphsByBoundary{};
    // std::unordered_map<std::vector<std::vector<HalfEdgeGraph>>, Graph> graphsByBoundary{};

    // TODO maybe just represent as single std::vector<Graph>? only previous tier matters I think since only used to construct next tier; otherwise just have all of them in order probably fine
    std::vector<Graph> tier1;

    std::map<HalfEdgeGraph, std::vector<std::tuple<Graph, unsigned int, unsigned int>>> heToPrims;
    //std::map<HalfEdgeGraph, std::vector<const Graph&>> heToPrims;

    std::vector<Graph> prevTier;
    std::set<HalfEdgeGraph> encounteredHalfEdges;
    for (const Primitive* p : primitives) {
        for (const HalfEdgeGraph& h : p->halfEdges) {
            encounteredHalfEdges.insert(h);
        }
    }
#if INCLUDE_SINGLE_EDGE
    for (const HalfEdgeGraph& h : encounteredHalfEdges) {
        Graph g = Graph(std::vector<HalfEdgeGraph>({h, -h}));
        g.isSingleEdge = true;
        g.face_colors = faceColors;
        // bool notYetAdded = true;
        // for (const Graph& oldG : hierarchy) {
        //     if (g.sameBoundaryString(oldG)) {
        //         // or could use isomorphism, dunno which more efficient (just avoiding hitting both HE and its opposite)
        //         notYetAdded = false;
        //         break;
        //     }
        // }
        // if (notYetAdded) {
        //     hierarchy.push_back(g);
        // }
        if (!graphsByBoundary.contains(g.getBoundaryGeneric())) {
            graphsByBoundary.insert({g.getBoundaryGeneric(), g});
        }
    }
#endif
    for (const Primitive* p : primitives) {
        // Graph g = Graph();
        Graph g = Graph(p->halfEdges);
        g.face_colors = faceColors;
        // g.face_colors.push_back(glm::vec3(0,0,0));
        // g.face_colors.push_back(glm::vec3(1,1,1));
        // TODO I GUESS should pass face colors in too? just placeholder rn

        // ignore redundant primitives
        bool noMatches = true;
        for (const Graph& existingGraph : tier1) {
            if (g.isIsomorphicTo(existingGraph)) {
                noMatches = false;
                break;
            }
        }
        if (noMatches) {
            tier1.push_back(g);
            graphsByBoundary.insert({g.getBoundaryGeneric(), g});
            prevTier.push_back(g);

            for (unsigned int j = 0; j < g.boundaryGraphElements.size(); ++j) {
            //for (const uPtr<BoundaryElem>& e : g.boundaryGraphElements) {
                const uPtr<BoundaryElem>& e = g.boundaryGraphElements.at(j);
                //heToPrims[-(e->he.h)].push_back(g);
                unsigned int edgeIndex;

                for (unsigned int i = 0; i < g.primitives.at(0)->halfEdges.size(); ++i) {
                    if (e->he.h == g.primitives.at(0)->halfEdges.at(i)) {
                        edgeIndex = i;
                        break;
                // TODO I THINK THIS IS WRONGLY SET
                        // TODO actually is j just always the same as edgeIndex is supposed to be?
                        // I don't actually see why it's wrongly set but yeah maybe can just use j
                        // nevermind j is wrong. not sure what my original issue with this was
                    }
                }
                HalfEdgeGraph heInv2 = -(e->he.h); 
                // //this^ should work I think? doesn't seem to but seems an issue elsewhere
                // // it is some issue with indexing though not finding in map so some value not lining up or similar
                // // alternatively try below approach
                // my brain can't work well enough rn finish later TODO 2:06am
                /*HalfEdgeGraph heInv;
                Primitive* p2 = p->connections.at(edgeIndex);
                for (unsigned int i = 0; i < p2->halfEdges.size(); ++i) {
                    if (p2->connections.at(i) == p) {
                        heInv = p2->halfEdges.at(i);
                        break;
                    }

                }

                if (heInv != heInv2) {
                    std::cout << "nonmatching inv" << std::endl;
                }*/
                heToPrims[heInv2].push_back({ g, edgeIndex, j });
            }
            
        }
        // TODO realize also might need to make rules for primitive to single edge case once single-edge-no-vertex graphs added
        //  NOTE IDK if primitives that are single edges should be occurring in input really but no reason not to make possible I suppose (collinear edge of face)
        //  Actually IDK. maybe just leave those out entirely since not really benefical to graph structure.
        //      not handling for now but TODO may add later; I don't think it's really beneficial though
        // fo (const Graph& )

        // tier1.push_back(g);
        // tier1.push_back(std::move(g));
        // tier1.insert(tier1.end(), g);
    }
    // hierarchy.push_back(tier1);
    // hierarchy.insert(hierarchy.end(), tier1.begin(), tier1.end());

    Graph emptyGraph{};
    emptyGraph.face_colors = faceColors; // shouldn't matter but just for consistency doing this

    std::vector<Graph> usedGraphs;

#if outputProgress
    QString fileName = QFileDialog::getSaveFileName(nullptr, "Save JSON File", "~/../../../../jsons", "JSON Files (*.json)");
#endif

    // TODO figure out end conditions
    // std::vector<Graph> prevTier = tier1;
    for (unsigned int iter = 0; iter < maxSteps; ++iter) {
        std::vector<Graph> newTier;
        std::map<std::vector<std::vector<HalfEdgeGraph>>, Graph> newGraphsByBoundary{};

        // TODO should add checks to avoid redundancy? general graph isomorphism test I guess should do
        for (const Graph& g : prevTier) {


            //for (const Graph& gPrim : tier1) {
            //    std::vector<Graph> branchGlued = g.branchGlue(gPrim);

            //    for (Graph& newGraph : branchGlued) {
            //        // bool noMatches = true;
            //        // only need to compare current tier since earlier ones all have fewer primitives and hence can't be isomorphic to this
            //        // for (const Graph& existingGraph : newTier) {
            //        //     if (newGraph.sameBoundaryString(existingGraph) && newGraph.isIsomorphicTo(existingGraph)) {
            //        //         noMatches = false;
            //        //         break;
            //        //     }
            //        // }
            //        // if (noMatches) {
            //        //     newTier.push_back(newGraph);
            //        // }
            //        const auto& bound = newGraph.getBoundaryGeneric();
            //        // if (!newGraphsByBoundary.contains(bound) || !newGraph.isIsomorphicTo(newGraphsByBoundary.at(bound))) {
            //        // auto bound = newGraph.getBoundaryGeneric();
            //        if (!newGraphsByBoundary.contains(bound) || !newGraph.isIsomorphicTo(newGraphsByBoundary.at(bound))) {
            //            newTier.push_back(newGraph);
            //            // if (bound != newGraph.getBoundaryGeneric()) {
            //            //     std::cout << "AAA" << std::endl;
            //            // }
            //            newGraphsByBoundary.insert({newGraph.getBoundaryGeneric(), (newTier.back())});
            //        }
            //    }
            //    // newTier.insert(newTier.end(), branchGlued.begin(), branchGlued.end());
            //}
            ////std::vector<Graph> loopGlued = g.loopGlue();
            std::vector<Graph> glued = g.loopGlue();

            for (unsigned int i = 0; i < g.boundaryGraphElements.size(); ++i) {
            /*for (const uPtr<BoundaryElem>& e : g.boundaryGraphElements) {*/
                try {
                    const uPtr<BoundaryElem>& e = g.boundaryGraphElements.at(i);
                    if (heToPrims.contains(e->he.h)) {
                        for (const auto& [primGraph, heIndex, beIndex] : heToPrims.at(e->he.h)) {
                            //for (const Graph& prim : heToPrims.at(e->he.h)) {
                            Graph gluedGraph(g);
                            Primitive* gPrim = gluedGraph.primitives.at(e->he.primIndex).get();
                            uPtr<Primitive> newPrim = mkU<Primitive>(*primGraph.primitives.at(0));
                            newPrim->connections.at(heIndex) = gPrim;
                            // TODO not great having a loop for this but not sure how I'd store; also not a very long loop
                            for (unsigned int j = 0; j < gPrim->halfEdges.size(); ++j) {
                                if (gPrim->halfEdges.at(j) == e->he.h) {
                                    gPrim->connections.at(j) = newPrim.get();
                                    break;
                                }
                            }

                            // update boundary graph
                            // TODO
                            //BoundaryElem* oE = nullptr;
                            //unsigned int oE_index;
                            // TODO

                            unsigned int gBoundarySize = gluedGraph.boundaryGraphElements.size();
                            unsigned int primBoundarySize = primGraph.boundaryGraphElements.size();

                            for (unsigned int j = 0; j < primBoundarySize; ++j) {
                                uPtr<BoundaryElem> newElem = mkU<BoundaryElem>(*primGraph.boundaryGraphElements.at(j));

                                newElem->he.primIndex += g.primitives.size();
                                gluedGraph.boundaryGraphElements.push_back(std::move(newElem));

                            }

                            for (unsigned int j = 0; j < primBoundarySize; ++j) {
                                BoundaryElem* target = gluedGraph.boundaryGraphElements.at(gBoundarySize + j).get();
                                BoundaryElem* nextElem = gluedGraph.boundaryGraphElements.at(gBoundarySize + (j + 1) % primBoundarySize).get();
                                target->next = nextElem;
                                nextElem->prev = target;
                            }
                            unsigned int oE_index = beIndex + gBoundarySize;
                            BoundaryElem* tE = gluedGraph.boundaryGraphElements.at(i).get();
                            BoundaryElem* oE = gluedGraph.boundaryGraphElements.at(oE_index).get();

                            tE->prev->next = oE->next;
                            oE->next->prev = tE->prev;
                            tE->next->prev = oE->prev;
                            oE->prev->next = tE->next;


                            gluedGraph.boundaryGraphElements.erase(gluedGraph.boundaryGraphElements.begin() + oE_index);
                            gluedGraph.boundaryGraphElements.erase(gluedGraph.boundaryGraphElements.begin() + i);


                            gluedGraph.primitives.push_back(std::move(newPrim));
                            gluedGraph.sortBoundaryGraphElements();

                            glued.push_back(std::move(gluedGraph));
                        }
                    }
                    else {
                        std::cerr << g.boundaryGraphElements.at(i)->he.h.angle << std::endl;
                    }
                }
                catch (const std::exception& e) {
                    std::cerr << "Glue exception" << std::endl;
                    std::cerr << e.what() << std::endl;
                }
            }

            // TODO maybe should use some sort of std::move thing here? dunno best way to
            // newTier.insert(newTier.end(), loopGlued.begin(), loopGlued.end());
            //for (Graph& newGraph : loopGlued) {
            for (Graph& newGraph : glued) {
                // bool noMatches = true;
                // for (const Graph& existingGraph : newTier) {
                //     if (newGraph.sameBoundaryString(existingGraph) && newGraph.isIsomorphicTo(existingGraph)) {
                //         noMatches = false;
                //         break;
                //     }
                // }
                // if (noMatches) {
                //     newTier.push_back(newGraph);
                // }
                const auto& bound = newGraph.getBoundaryGeneric();
                // if (!newGraphsByBoundary.contains(bound) || !newGraph.isIsomorphicTo(newGraphsByBoundary.at(bound))) {
                // auto bound = newGraph.getBoundaryGeneric();
                if (!newGraphsByBoundary.contains(bound) || !newGraph.isIsomorphicTo(newGraphsByBoundary.at(bound))) {
                    newTier.push_back(newGraph);

                    newGraphsByBoundary.insert({newGraph.getBoundaryGeneric(), (newTier.back())});

                }
            }

        }

        // for ()
        // I think going to make the whole tier first before checking for boundary equivalences in hierarchy? since then can make sure not to repeat multiple of same graph
        prevTier = std::vector<Graph>();


        std::vector<Graph> toAdd;
        for (Graph& newGraph : newTier) {

#if !THREE_DIMENSIONAL
            float curTurn = 0;
            bool skip = false;
            // TODO this isn't really doing it right I think but curious if it works at all
            for (unsigned int j = 0; j < newGraph.boundaryString.size() * 2; ++j) {
                const Turn* t = std::get_if<Turn>(&newGraph.boundaryString.at((j) % newGraph.boundaryString.size()));
                if (t == nullptr) {
                    curTurn = 0;
                } else {
                    if (*t == Turn::positive) {
                        curTurn += 180;
                    } else {
                        curTurn -= 180;
                    }
                    if (abs(curTurn) > 360) {
                        skip = true;
                        break;
                    }
                }
            }
            if (skip) {
                usedGraphs.push_back(newGraph);
                continue;
            }
#endif

            // TODO move check earlier so don't need to check ones that are descendants in the noMatches checks?
            bool isUsedDescendant = false;
            for (const Graph& otherGraph : usedGraphs) {
                if (newGraph.hasSubgraph(otherGraph)) {
                    isUsedDescendant = true;
                    break;
                }
            }

            if (!isUsedDescendant) {
                if (newGraph.sameBoundaryString(emptyGraph)) {
                // if (newGraph.boundaryString.size() == 1) {
                    // TODO I think don't need to actually check that the one symbol is positive turn? since should never be able to make a string without a positive turn in it I believe
                    //   changed to a way that's clearer in use although will be less efficient, can change back if desired
                    // complete graph
                    result.push_back({emptyGraph, newGraph});
                    usedGraphs.push_back(newGraph);
                    std::cout << "Added starter rule" << std::endl;
                } else {
                    bool ruleAdded = false;
                    // for (const Graph& otherGraph : hierarchy) {
                    //     // TODO should make spliced version of this. I think should actually include support for both so can control whether produces connected graph
                    //     if (newGraph.sameBoundaryString(otherGraph)) {
                    //         result.push_back({otherGraph, newGraph});
                    //         ruleAdded = true;
                    //         break;
                    //         // TODO uncertain if this is meant to always break here. I think should theoretically never hit more than one since all earlier ones after the first shouldn't be added to hierarchy
                    //     }
                    // }
                    auto bound = newGraph.getBoundaryGeneric();
                    if (graphsByBoundary.contains(bound)) {

                        //  trying having a check here for if newGraph can generate positions w/ very loose constraints.
                        //  if not, remove but don't make rule
                        //  theoretically could check at every time a graph generated and would make hierarchy smaller but time consuming (though also time consuming having big hierarchy) and riskier since will probably over-remove cases (which this also will but fewer chances to fail)
                        if (newGraph.isIsomorphicTo(graphsByBoundary.at(bound))) {

                            // TODO find why this isn't being caught earlier
                            // !!!!!!!!!!!


                            // Graph g = graphsByBoundary.at(bound);
                            // std::cout << newGraph.isIsomorphicTo(graphsByBoundary.at(bound)) << std::endl;
                            std::cout << "Repeated element in tier " << iter << std::endl;
                        } else {
                            double err;
                            std::map<unsigned int,glm::vec3> emptyPositionMap;
                            newGraph.samplePositions(emptyPositionMap, err, 1, 20, -10, 10, 10000, 0.5);
                            if (err != -1) {
                                result.push_back({graphsByBoundary.at(bound),newGraph});
                            }
                        }
                        ruleAdded = true;
                    }
                    // TODO figure out how to handle the special case of single-edged graph
                    //  if (newGraph.boundaryString.size() == 3) {
                    //      ...
                    //  }
                    if (ruleAdded) {

                        // TODO should this also get rid of descendants of the otherGraph? I don't think so since can have multiple that make valid rules
                        usedGraphs.push_back(std::move(newGraph));
                    } else {
                        // TODO these shouldn't copy probably, just reference, dunno how to avoid
                        //  TODO maybe prevtier as vector of pointers? or all as vectors of pointers?
                        prevTier.push_back(std::move(newGraph));
                        graphsByBoundary.insert({newGraph.getBoundaryGeneric(), (prevTier.back())});

                        // hierarchy.push_back(std::move(newGraph));
                        // TODO is it better to not add to hierarchy until after? with hierarchy push here can have rules that go between two graphs at same level but I'm unclear if that's intended or not in the algorithm as laid out in the paper. seems to make smaller grammars at least
                        // toAdd.push_back(std::move(newGraph));
                    }
                    // TODO check for guarantee of incomplete descendants and remove if so
                    //  not technically required at this moment (since we have arbitrary cutoff rn) but long term good to have

                    // TODO do we want to skip descendants of graphs already used as L of rules?
                    //  should be easy: store list of graphs that weren't added to hierarchy separately, then when adding one check if it has any of those as a subgraph
                    // TODO yeah want to do that. additional note:
                    //      STEP that I BELIEVE is optional ("algorithm can finish without it" but can make grammar have fewer rules)
                    //          when graph added to hierarchy, check if it can be used to simplify another graph
                    //              only case where this should matter I THINK is "reducing a graph with its descendants"
                    //          in that case, we want to remove all descendants besides the ones used to reduce it from the hierarchy
                    //      hence in the case where we're not doing that I THINK the idea is we remove all of its decesndants from the hierarchy
                    //      but not super explicitly specified in paper as far as I can tell
                    //          but will try it
                }
            }
        }
        // for (Graph &g : toAdd) {
        //     hierarchy.push_back(std::move(g));
        // }

        // hierarchy.push_back(newTier);
#if outputProgress
//progress saving for testing
        std::vector<std::array<Graph*,2>> gram;
        for (auto& [first, second] : result) {
            std::array<Graph*,2> gArr{&first,&second};
            gram.push_back(gArr);
        }
        JSONReader::WriteGrammarToFile(fileName.chopped(5) + QString::number(i) + ".json", gram);
#endif

        // prevTier = newTier;

    }

    // TODO handle generation 0 special case

    // TODO actual use of hierarchy to produce rules

    // TODO really should try to add some culling of incomplete-descendant cases since hierarchy does grow fast when more than a few prmitive types

    return result;
}

// TODO seems in 3d case examples I've tried to generate a lot of rules which don't do anything useful
//  I think are technically valid rules though in how defined right now, just e.g. things that swap around connections but don't add anything ? dunno if a way to avoid

std::vector<std::pair<Graph, Graph> > Graph::generateRules(unsigned int maxSteps, bool includeWholeStarter) const
{
    // generateRules()
    std::vector<Primitive*> primPtrs;
    for (const uPtr<Primitive>& p : this->primitives) {
        primPtrs.push_back(p.get());
    }
    std::vector<std::pair<Graph, Graph>> rules = generateRules(primPtrs, this->face_colors, maxSteps);
    if (includeWholeStarter) {
        rules.push_back({Graph(),*this});
    }
    return rules;
}





std::vector<Graph> Graph::generateHierarchy(const std::vector<Primitive *> &primitives, int steps, std::vector<glm::vec3> faceColors)
{
    // TODO should maybe update to have consistent behavior with actual generateGrammar function? but was just for testing purposes anyway


    std::vector<Graph> result;
    std::vector<std::vector<Graph>> hierarchy;
    std::vector<Graph> tier1;
    for (const Primitive* p : primitives) {
        Graph g = Graph(p->halfEdges);
        g.face_colors = faceColors;

        // ignore redundant primitives
        bool noMatches = true;
        for (const Graph& existingGraph : tier1) {
            if (g.isIsomorphicTo(existingGraph)) {
                noMatches = false;
                break;
            }
        }
        if (noMatches) {
            tier1.push_back(g);
        }

        result.push_back(g);
    }
    hierarchy.push_back(tier1);

    std::vector<Graph> prevTier = tier1;
    for (int i = 0; i < steps; ++i) {
        std::vector<Graph> newTier;
        for (const Graph& g : prevTier) {
            for (const Graph& gPrim : tier1) {
                std::vector<Graph> branchGlued = g.branchGlue(gPrim);

                for (Graph& newGraph : branchGlued) {
                    bool noMatches = true;
                    // only need to compare current tier since earlier ones all have fewer primitives and hence can't be isomorphic to this
                    for (const Graph& existingGraph : newTier) {
                        if (newGraph.isIsomorphicTo(existingGraph)) {
                            noMatches = false;
                            break;
                        }
                    }
                    if (noMatches) {
                        newTier.push_back(newGraph);
                        result.push_back(newGraph);
                    }
                }
                // newTier.insert(newTier.end(), branchGlued.begin(), branchGlued.end());
            }
            std::vector<Graph> loopGlued = g.loopGlue();
            // TODO maybe should use some sort of std::move thing here? dunno best way to
            // newTier.insert(newTier.end(), loopGlued.begin(), loopGlued.end());
            for (Graph& newGraph : loopGlued) {
                bool noMatches = true;
                for (const Graph& existingGraph : newTier) {
                    if (newGraph.isIsomorphicTo(existingGraph)) {
                        noMatches = false;
                        break;
                    }
                }
                if (noMatches) {
                    newTier.push_back(newGraph);
                    result.push_back(newGraph);
                }
            }

        }
        hierarchy.push_back(newTier);
        prevTier = newTier;

    }
    return result;
}

std::vector<Graph> Graph::generateHierarchy(int steps) const
{
    std::vector<Primitive*> primPtrs;
    for (const uPtr<Primitive>& p : this->primitives) {
        primPtrs.push_back(p.get());
    }
    return generateHierarchy(primPtrs, steps, this->face_colors);
}



// int main2() {
//     ClpSimplex model;
//     int numRows = 2, numCols = 2;

//     std::vector<double> objective = {1.0, 2.0};
//     std::vector<double> constraintMatrix = {1.0, 1.0, 2.0, 1.0};
//     std::vector<double> rhs = {4.0, 5.0};
//     std::vector<double> lowerBounds = {0.0, 0.0};
//     std::vector<double> upperBounds = {INFINITY, INFINITY};
//     std::vector<double> rowLower(numRows, -INFINITY);
//     std::vector<double> rowUpper = rhs;
//     std::vector<char> sense(numRows, 'L');

//     // ClpObjective
//     // model.setObjective(objective.data(), true);
//     // model.setObjective
//     model.setObjectiveCoefficient(0, 1.0);
//     model.setObjectiveCoefficient(1, 2.0);

//     for (int i = 0; i < numCols; i++) {
//         model.setColumnLower(i, lowerBounds[i]);
//         model.setColumnUpper(i, upperBounds[i]);
//     }

//     for (int i = 0; i < numRows; i++) {
//         model.setRowLower(i, rowLower[i]);
//         model.setRowUpper(i, rowUpper[i]);
//     }

//     for (int i = 0; i < numRows; i++) {
//         std::vector<int> indices(numCols);
//         std::vector<double> values(numCols);
//         for (int j = 0; j < numCols; j++) {
//             indices[j] = j;
//             values[j] = constraintMatrix[i * numCols + j];
//         }
//         model.addRow(indices.size(), indices.data(), values.data(), rhs[i], sense[i]);
//     }

//     model.primal();
//     double *solution = model.primalColumnSolution();
//     for (int i = 0; i < numCols; i++) {
//         printf("x%d = %f\n", i + 1, solution[i]);
//     }

//     return 0;
// }




// via https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
// not sure what best (quickest) test for such is especially if can be made faster with things we can assume about structure
// TODO this one only works in 2D
bool onSegment(const glm::vec3& p, const glm::vec3& q, const glm::vec3& r)
{
    if (q.x <= glm::max(p.x, r.x) && q.x >= glm::min(p.x, r.x) &&
        q.z <= glm::max(p.y, r.z) && q.z >= glm::min(p.z, r.z))
        return true;

    return false;
}
int orientation(const glm::vec3& p, const glm::vec3& q, const glm::vec3& r)
{
    int val = (q.z - p.z) * (r.x - q.x) -
              (q.x - p.x) * (r.z - q.z);
    if (val == 0) return 0;
    return (val > 0)? 1: 2;
}

bool doIntersect(const glm::vec3& p1, const glm::vec3& q1, const glm::vec3& p2, const glm::vec3& q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}




#if 1

std::vector<glm::vec3> Graph::samplePositions(const std::map<unsigned int, glm::vec3>& setValues, double& oError, float minEdgeLength, float maxEdgeLength, float minPosition, float maxPosition, unsigned int maxTries, float cosMin) const
{
    oError = -1;
    if (this->primitives.size() == 0) {
        return std::vector<glm::vec3>(0);
    }
    // std::vector<std::pair<int,Primitive*>> unmatchedPrims;
    std::map<Primitive*, int> allPrimIndices;

    // TODO maybe use some map
    // std::vector<Primitive*> unmatchedPrims;
    for (unsigned int i = 0; i < this->primitives.size(); ++i) {
        // unmatchedPrims.push_back({i, this->primitives.at(i).get()});
        allPrimIndices.insert({this->primitives.at(i).get(), i});
    }

    std::map<Primitive*, int> unmatchedPrimIndices = allPrimIndices;

    std::vector<glm::vec3> result;
    // std::vector<glm::vec3> result(this->primitives.size());


    // Eigen::MappedSparseMatrix<float> system();

    // system.insert()
    // matrix:
    // idea: columns -> (x,y,z) coordinates of vertices in order then pairs of edge lengths?
    // e.g.
    //      [ p0x ]
    //      [ p0y ]
    //      [ p0z ]
    //      [ p1x ]
    //      [ p1y ]
    //      [ p1z ]
    //      [ p2x ]
    //      [ p2y ]
    //      [ p2z ]
    //      [ l01 ]
    //      [ l02 ]
    //      [ l12 ]

    // X height = A width = primitives.size() * 3 + #edges
    // B = all 0s

    // unsigned int edgeIndex = 0;
    // unsigned int edgeIndex = this->primitives.size() * 3;
    // each edge length used just once, so just increment position each time one used
    // unsigned int rowsAdded = 0;


    // std::vector<Eigen::Triplet<float>> tripletList;
    // tripletList.reserve(this->...)

    // lprec* lp;

    // lp = make_lp(0,4);


    // delete_lp(lp);




    // return(ret);

    // operations_research::MPVariable* v;
    // std::vector<operations_research::MPVariable*> positionVariables;
    // std::vector<operations_research::MPVariable*> edgeLengthVariables;

    // operations_research::MPSolver* m;

    // operations_research::MPSolver::CreateSolver("SCIP");
    // uPtr<operations_research::MPSolver> solver(operations_research::MPSolver::CreateSolver("GLOP"));
    // operations_research::MPSolver solver("MP_Solver", operations_research::MPSolver::GLOP_LINEAR_PROGRAMMING);
    // uPtr<operations_research::MPSolver> solver(operations_research::MPSolver::CreateSolver("GLOP"));

    // for (unsigned int i = 0; i < this->primitives.size(); ++i) {
    //     positionVariables.push_back(solver->MakeNumVar(minPosition, maxPosition, "Prim " + std::to_string(i) + " x"));
    //     positionVariables.push_back(solver->MakeNumVar(minPosition, maxPosition, "Prim " + std::to_string(i) + " y"));
    //     positionVariables.push_back(solver->MakeNumVar(minPosition, maxPosition, "Prim " + std::to_string(i) + " z"));
    // }

    // operations_research::MPObjective* const objective = solver->MutableObjective();

    // std::vector<std::pair<int,int>> indexPairs;
    std::vector<std::tuple<int,int,glm::vec3>> constraints;
    std::map<int,std::pair<int,glm::vec3>> constraintMap;

    // std::vector<glm::vec3>

    // TODO I think might have some issue with constraintMap that makes it less reliable than it should be? IDK how I had it set up totally so something to look into

    // TODO probably skip p0 and make always 0,0,0
    for (const auto& [prim, index] : unmatchedPrimIndices) {
        // TODO try: add just edges that have non-null connections and angle in upper half of values (>= 0)
        //  others are redundant
        for (unsigned int i = 0; i < prim->connections.size(); ++i) {
            if (prim->connections.at(i) != nullptr && prim->halfEdges.at(i).angle >= 0) {
                unsigned int otherIndex = allPrimIndices.at(prim->connections.at(i));

                glm::vec3 direction(cos(prim->halfEdges.at(i).angle), 0.0f, sin(prim->halfEdges.at(i).angle));

#if THREE_DIMENSIONAL
                glm::vec3 forward = glm::vec3(0,1,0);
                glm::vec3 normal;
                auto& norms = prim->halfEdges.at(i).faceNormals;
                if (norms.at(0).x == norms.at(1).x) {
                    if (norms.at(0).y == norms.at(1).y) {
                        if (norms.at(0).z == norms.at(1).z) {
                            normal = norms.at(0);

                        } else {
                            normal = (norms.at(0).z < norms.at(1).z) ? norms.at(0) : norms.at(1);
                        }
                    } else {
                        normal = (norms.at(0).y < norms.at(1).y) ? norms.at(0) : norms.at(1);
                    }
                } else {
                    normal = (norms.at(0).x < norms.at(1).x) ? norms.at(0) : norms.at(1);
                }

                glm::vec3 turnDirection = glm::cross(normal, forward);
                if (turnDirection != glm::vec3(0)) {
                    float turnAngle = glm::acos(glm::dot(normal, forward));

                    glm::mat4 rotation = glm::rotate(-turnAngle, turnDirection);

                    direction = glm::vec3(rotation * glm::vec4(direction, 1));
                // } else if (glm::dot(normal, forward) < 0) {
                    // direction.x = -direction.x; // TODO not sure if that's totally right
                        // unnecessary it seems?
                }
#endif

                constraints.push_back({index, otherIndex, direction});

                if (!constraintMap.contains(index) || constraintMap.at(index).first > otherIndex) {
                    constraintMap.insert_or_assign(index,std::pair<int,glm::vec3>({otherIndex, direction}));
                }
                if (!constraintMap.contains(otherIndex) || constraintMap.at(otherIndex).first > index) {
                    constraintMap.insert_or_assign(otherIndex,std::pair<int,glm::vec3>({index, -direction}));
                }

                // operations_research::MPVariable* newEdgeVar = solver->MakeNumVar(minPosition, maxPosition, "Edge " + std::to_string(edgeIndex));

                // // TODO maybe make within some epsilon?
                // operations_research::MPConstraint* const c0 = solver->MakeRowConstraint(0.0,0.0);
                // operations_research::MPConstraint* const c1 = solver->MakeRowConstraint(0.0,0.0);
                // operations_research::MPConstraint* const c2 = solver->MakeRowConstraint(0.0,0.0);
                // c0->SetCoefficient(positionVariables.at(index * 3), -1);
                // c0->SetCoefficient(positionVariables.at(otherIndex * 3), 1);
                // c0->SetCoefficient(newEdgeVar, direction.x);

                // c1->SetCoefficient(positionVariables.at(index * 3 + 1), -1);
                // c1->SetCoefficient(positionVariables.at(otherIndex * 3 + 1), 1);
                // c1->SetCoefficient(newEdgeVar, direction.y);

                // c2->SetCoefficient(positionVariables.at(index * 3 + 2), -1);
                // c2->SetCoefficient(positionVariables.at(otherIndex * 3 + 2), 1);
                // c2->SetCoefficient(newEdgeVar, direction.z);


                // operations_research::MPConstraint* const c3 = solver->MakeRowConstraint(minEdgeLength, maxEdgeLength);
                // c3->SetCoefficient(newEdgeVar, 1);


                // objective->SetCoefficient(newEdgeVar, 1);

                // edgeLengthVariables.push_back(newEdgeVar);



                // add to matrix
                // A.insert(rowsAdded, index * 3) = -1;
                // A.insert(rowsAdded+1, index * 3+1) = -1;
                // A.insert(rowsAdded+2, index * 3+2) = -1;

                // A.insert(rowsAdded, otherIndex * 3) = 1;
                // A.insert(rowsAdded+1, otherIndex * 3+1) = 1;
                // A.insert(rowsAdded+2, otherIndex * 3+2) = 1;

                // A.insert(rowsAdded, edgeIndex) = direction.x;
                // A.insert(rowsAdded+1, edgeIndex) = direction.y;
                // A.insert(rowsAdded+2, edgeIndex) = direction.z;

                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded,     index * 3,     -1));
                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded + 1, index * 3 + 1, -1));
                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded + 2, index * 3 + 2, -1));

                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded,     otherIndex * 3,     1));
                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded + 1, otherIndex * 3 + 1, 1));
                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded + 2, otherIndex * 3 + 2, 1));

                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded,     edgeIndex, direction.x));
                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded + 1, edgeIndex, direction.y));
                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded + 2, edgeIndex, direction.z));


                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded + 3,     edgeIndex, 1));

                // rowsAdded += 3;
                // ++edgeIndex;

            }
        }

    }


    // objective->SetMinimization();
    // const operations_research::MPSolver::ResultStatus result_status = solver->Solve();
    // // result_status

    // std::cout << objective->Value() << std::endl;
    // for (const auto& v : positionVariables) {
    //     std::cout << v->name() << " = " << v->solution_value() << std::endl;
    // }
    // for (const auto& v : edgeLengthVariables) {
    //     std::cout << v->name() << " = " << v->solution_value() << std::endl;
    // }



    // Eigen::SparseMatrix<float> A(rowsAdded, edgeIndex);
    // A.setFromTriplets(tripletList.begin(), tripletList.end());

    // // Eigen::SparseVector<float> B(rowsAdded + 3);
    // Eigen::VectorXf B = Eigen::VectorXf::Zero(rowsAdded);
    // // for (unsigned int i = 3; i < rowsAdded; i += 4) {
    // //     B[i] = 1;
    // // }

    // // Eigen::SparseQR<Eigen::SparseMatrix<float>> solver;

    // // TODO no clue best ordering

    // std::cout << A << std::endl;
    // std::cout << B << std::endl;

    // Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<Eigen::SparseMatrix<float>::StorageIndex>> solver;

    // solver.compute(A);
    // Eigen::VectorXf X = solver.solve(B);

    // std::cout << X << std::endl;



    // for (unsigned int i = 0; i < this->primitives.size(); ++i) {
    //     glm::vec3 position(X[i * 3], X[i * 3 + 1], X[i * 3 + 2]);
    //     result.push_back(position);
    // }





    // while (unmatchedPrimIndices.size() > 0) {

    //     auto startData = unmatchedPrimIndices.begin();
    //     Primitive* startP = startData->first;

    //     std::vector<Primitive*> encounteredPrims;
    //     encounteredPrims.push_back(startP);

    //     result.at(startData->second) = glm::vec3(0);
    //     // TODO set start positions of each graph section; not sure how wanna do (just random I guess?)

    //     // TODO this basically is a system of linear equations, maybe try using a library built for those?

    //     for (Primitive* p : encounteredPrims) {
    //         std::erase_if(unmatchedPrimIndices, [&p](std::pair<Primitive*, int> curPrim) {
    //             return curPrim.first == p;
    //         });
    //     }

    // }



    // TODO WILL TRY as an optimization problem instead actually, since want to restrict range
    // could also just try randomly setting values and checking against constraints which I guess is what the paper does but not sure if they have like a smarter way of picking values or if just totally random
    // ClpModel....
    // ClpMatrixBase
    // CoinPackedMatrix m;





    // double objective[] = { /* objective function coefficients */ };
    // double lower[] = { /* variable lower bounds */ };
    // double upper[] = { /* variable upper bounds */ };
    // double rhsLower[] = { /* constraint lower bounds */ };
    // double rhsUpper[] = { /* constraint upper bounds */ };
    // int numCols = /* number of columns */;
    // int numRows = /* number of rows */;

    // ClpSimplex  model;
    // int status;
    // // // if (argc<2)
    // // model.loadProblem(m, lower, upper, objective, rhsLower, rhsUpper,
    // //                   nullptr, nullptr, numCols, numRows);
    // // model.loadProblem(matrix, lower, upper, objective, rhsLower, rhsUpper,
    // //                   nullptr, nullptr, numCols, numRows);
    //     // status=model.readMps("../../Mps/Sample/p0033.mps");
    // // status = model.loadProblem(m)
    // // else
    //     // status=model.readMps(argv[1]);
    // if (!status) {
    //     model.primal();
    // }

    // const double infinity = solver->infinity();


    // main2();
    // return result;

    // RANDOM SAMPLE TRY
    // TODO MAKE BETTER
    std::random_device rd;
    std::mt19937 gen(rd());
    // std::uniform_int_distribution<> distrib(0, potentialPairs.size() - 1);
    std::uniform_real_distribution<float> dist(minPosition, maxPosition);
    std::uniform_real_distribution<float> edgeDist(minEdgeLength + 0.5f, maxEdgeLength - 1.f);
    // std::uniform_real_distribution<float> edgeDist(minEdgeLength, maxEdgeLength);

    // bool notDone = true;
    // const float sampleEpsilon = 0.0001f;
    // const float cosMin = 0.995f;
    // const float cosMin = 0.98f;

    const glm::vec3 zeroVector(0);
    const unsigned int primSize = this->primitives.size();
    // int triesLeft = 30000;
    std::vector<glm::vec3> startTry(primSize);
    if (setValues.size() == 0) {
        startTry.at(0) = glm::vec3(0);
    } else {
        for (const auto& [index, pos] : setValues) {
            startTry.at(index) = pos;
        }
    }
    unsigned int tries;
    for (tries = 0; tries < maxTries; ++tries) {
    // while (--triesLeft >= 0) {
        // result = std::vector<glm::vec3>(primSize);

        result = startTry;
        // result.at(0) = glm::vec3(0);
        for (unsigned int i = 1; i < primSize; ++i) {

            const std::pair<int,glm::vec3>& mapped = constraintMap.at(i);
            if (mapped.first < i) {
                // other already mapped
                float len = edgeDist(gen);
                result.at(i) = result.at(mapped.first) - len * mapped.second;
                // TODO should probably have this change up order some rather than just same order as list
                // if
            } else {
            // if ()
#if THREE_DIMENSIONAL
                result.at(i) = glm::vec3(dist(gen), dist(gen), dist(gen));
#else
                result.at(i) = glm::vec3(dist(gen), 0.f, dist(gen));
#endif
            }

        }

        bool fail = false;
        for (const auto& [id1, id2, dir] : constraints) {
            glm::vec3 pos1 = result.at(id1);
            glm::vec3 pos2 = result.at(id2);
            glm::vec3 diff = pos2 - pos1;
            // length check
            float d = glm::length(diff);
#if 1
            if (d < minEdgeLength) {
#else
            if (d < minEdgeLength || d > maxEdgeLength) {
#endif
                fail = true;
                break;
            }
            // direction check
            glm::vec3 nDiff = glm::normalize(diff);

            float dot = glm::dot(nDiff, dir);
            if (dot < cosMin) {
                fail = true;
                break;
            }

            // if (nDiff.x >= dir.)

            //TODO face intersection check instead
            // intersection check
            // TODO maybe run all intersection checks after all position checks done? idk if really matters
#if 0
            for (const auto& [idB1, idB2, dirB] : constraints) {
                if (id1 != idB1 && id1 != idB2 && id2 != idB1 && id2 != idB2) {
                    glm::vec3 posB1 = result.at(idB1);
                    glm::vec3 posB2 = result.at(idB2);
                    // glm::vec3 diffB = posB2 - posB1;
                    // glm::vec3 p = glm::cross(diff, diffB);
                    // if (p != zeroVector) {
                    //     // not parallel, lines intersect; check if in segment range
                    //     //  TODO will have to test different when in 3D (also need to test face intersection too then)
                    //     glm::vec3 startDiff = posB1 - pos1;
                    //     float testVal = glm::dot(glm::cross(startDiff, posB1), p) / glm::dot(p,p);
                    //     if (testVal > 0 && testVal < 1) {
                    //         fail = true;
                    //         break;
                    //     }
                    // }
                    // TODO make sure this test is right
                    if (doIntersect(pos1, pos2, posB1, posB2)) {
                        fail = true;
                        break;
                    }
                }
            }
            if (fail) {
                break;
            }
#endif
        }
        if (!fail) {


            break;
        }

    }


    if (tries >= maxTries) {
        // if (triesLeft < 0) {
        oError = -1;
        // result = std::vector<glm::vec3>(this->primitives.size());

        // for (unsigned int i = 0; i < result.size(); ++i) {
        //     result.at(i) = glm::vec3(i,0,i);
        // }
    } else {
        oError = 0.0; // TODO some better measure of error so can control precision of whole system
    }




    // TODO idea: try eigen approach again but each time it generates set of values if any edge length constraints broken try adding a row that sets one edge to a random length, repeat
    // or see fi can get basis and use that? https://stackoverflow.com/questions/54402199/what-is-the-best-way-to-solve-a-system-of-linear-equations-with-infinite-solutio

    return result;
}




#else

// Trying out an alternative approach
//  solve linear equation system and repeatedly set edge lengths randomly each time a constraint fails
// runs much faster and seems to give OK results but some with misalignment I'm not sure how it's occurring
//  (not sure if limitation of this technique or just a bug on my part)
//      the other one has some alignment issues but that's a controllable parameter on my part (the epsilon for checking angle closeness), not sure why this one has it
//  figured out how to add check for precision of solution, and looks pretty good now. Might try sticking with this approach or adding an option for which to use?
//   though still plan to add intersection check, saving between steps of algorithm, and some sort of way to separate multiple subgraphs when the option for multiple starters is on
//    this is much faster but because it's less random it seems possible that this would be more likely to get "stuck" in regards especially to the intersection constraint which is something I don't think I can check within the solver itself
//      'least I've no idea if one could set up an linear equation to check for the intersection. I don't think so
//    the latter aspect there is something handled by the other approach inherently (all placed randomly anyway besides the first position of the whole graph object, while this solver gives 0,0,0 for one per each subgraph) but uncertain how to deal with in this one
// for both approaches still need to add some intersection checks for constraints
std::vector<glm::vec3> Graph::samplePositions(const std::map<unsigned int, glm::vec3>& setValues, double& oError, float minEdgeLength, float maxEdgeLength, float minPosition, float maxPosition) const
{
    oError = -1.0;
    if (this->primitives.size() == 0) {
        return std::vector<glm::vec3>(0);
    }
    std::map<Primitive*, int> allPrimIndices;

    // TODO maybe use some map
    // std::vector<Primitive*> unmatchedPrims;
    for (unsigned int i = 0; i < this->primitives.size(); ++i) {
        // unmatchedPrims.push_back({i, this->primitives.at(i).get()});
        allPrimIndices.insert({this->primitives.at(i).get(), i});
    }

    std::map<Primitive*, int> unmatchedPrimIndices = allPrimIndices;

    std::vector<glm::vec3> result;
    // std::vector<glm::vec3> result(this->primitives.size());


    // Eigen::MappedSparseMatrix<float> system();

    // system.insert()
    // matrix:
    // idea: columns -> (x,y,z) coordinates of vertices in order then pairs of edge lengths?
    // e.g.
    //      [ p0x ]
    //      [ p0y ]
    //      [ p0z ]
    //      [ p1x ]
    //      [ p1y ]
    //      [ p1z ]
    //      [ p2x ]
    //      [ p2y ]
    //      [ p2z ]
    //      [ l01 ]
    //      [ l02 ]
    //      [ l12 ]

    // X height = A width = primitives.size() * 3 + #edges
    // B = all 0s

    // unsigned int edgeIndex = 0;
    unsigned int edgeIndex = this->primitives.size() * 3;
    unsigned int edgeStartIndex = edgeIndex;
    // each edge length used just once, so just increment position each time one used
    //  now used again later but doesn't care which one it's setting per se
    unsigned int rowsAdded = 0;


    std::vector<Eigen::Triplet<double>> tripletList;
    // tripletList.reserve(this->...)

    // std::vector<std::pair<int,int>> indexPairs;
    std::vector<std::tuple<int,int,glm::vec3>> constraints;
    std::map<int,std::pair<int,glm::vec3>> constraintMap;

    // std::vector<glm::vec3>

    // TODO probably skip p0 and make always 0,0,0
    for (const auto& [prim, index] : unmatchedPrimIndices) {
        // TODO try: add just edges that have non-null connections and angle in upper half of values (>= 0)
        //  others are redundant
        for (unsigned int i = 0; i < prim->connections.size(); ++i) {
            if (prim->connections.at(i) != nullptr && prim->halfEdges.at(i).angle >= 0) {
                unsigned int otherIndex = allPrimIndices.at(prim->connections.at(i));
                glm::vec3 direction(cos(prim->halfEdges.at(i).angle), 0.0f, sin(prim->halfEdges.at(i).angle));

                // constraints.push_back({index, otherIndex, direction});

                // if (!constraintMap.contains(index) || constraintMap.at(index).first > otherIndex) {
                //     constraintMap.insert_or_assign(index,std::pair<int,glm::vec3>({otherIndex, direction}));
                // }
                // if (!constraintMap.contains(otherIndex) || constraintMap.at(otherIndex).first > index) {
                //     constraintMap.insert_or_assign(otherIndex,std::pair<int,glm::vec3>({index, -direction}));
                // }

                // add to matrix

                tripletList.push_back(Eigen::Triplet<double>(rowsAdded,     index * 3,     -1));
                tripletList.push_back(Eigen::Triplet<double>(rowsAdded + 1, index * 3 + 1, -1));
                tripletList.push_back(Eigen::Triplet<double>(rowsAdded + 2, index * 3 + 2, -1));

                tripletList.push_back(Eigen::Triplet<double>(rowsAdded,     otherIndex * 3,     1));
                tripletList.push_back(Eigen::Triplet<double>(rowsAdded + 1, otherIndex * 3 + 1, 1));
                tripletList.push_back(Eigen::Triplet<double>(rowsAdded + 2, otherIndex * 3 + 2, 1));

                tripletList.push_back(Eigen::Triplet<double>(rowsAdded,     edgeIndex, -direction.x));
                tripletList.push_back(Eigen::Triplet<double>(rowsAdded + 1, edgeIndex, -direction.y));
                tripletList.push_back(Eigen::Triplet<double>(rowsAdded + 2, edgeIndex, -direction.z));


                // tripletList.push_back(Eigen::Triplet<float>(rowsAdded + 3,     edgeIndex, 1));

                rowsAdded += 3;
                ++edgeIndex;

            }
        }

    }

    if (tripletList.size() == 0) {
        return std::vector<glm::vec3>(0);
    }

    unsigned int setRow = rowsAdded;
    for (const auto& [idx, pos] : setValues) {
        tripletList.push_back(Eigen::Triplet<double>(rowsAdded, idx * 3, 1));
        tripletList.push_back(Eigen::Triplet<double>(rowsAdded + 1, idx * 3 + 1, 1));
        tripletList.push_back(Eigen::Triplet<double>(rowsAdded + 2, idx * 3 + 1, 1));
        rowsAdded += 3;
    }


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(minPosition, maxPosition);
    std::uniform_real_distribution<double> edgeDist(minEdgeLength, maxEdgeLength);



    // // Eigen::SparseVector<float> B(rowsAdded + 3);
    // // for (unsigned int i = 3; i < rowsAdded; i += 4) {
    // //     B[i] = 1;
    // // }

    // Eigen::SparseQR<Eigen::SparseMatrix<float>> solver;

    // // TODO no clue best ordering

    // std::cout << A << std::endl;
    // std::cout << B << std::endl;

    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::StorageIndex>> solver;


    // std::cout << X << std::endl;
    unsigned int maxRows = rowsAdded + edgeIndex - edgeStartIndex;
    Eigen::SparseMatrix<double> A(maxRows, edgeIndex);
    A.conservativeResize(rowsAdded, edgeIndex);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    Eigen::VectorXd B = Eigen::VectorXd::Zero(maxRows);
    B.conservativeResize(rowsAdded);

    for (const auto& [idx, pos] : setValues) {
        B[setRow] = pos.x;
        B[setRow + 1] = pos.y;
        B[setRow + 2] = pos.z;
        setRow += 3;
    }

    unsigned int startRowsAdded = rowsAdded;

    Eigen::VectorXd X;
    bool unsatisfied = true;
    bool tryFail = false;
    const unsigned int maxTries = 120;
    unsigned int increaseAt = 15;
    const unsigned int increaseAtStep = 5; //disabled
    const double increaseBy = 10.0;
    double errorLimit = 1.0e-10;
    double relativeError = -1.0;
    unsigned int tries;
    for (tries = 0; tries < maxTries; ++tries) {
        do {
            A.makeCompressed();
            solver.compute(A);
            if (solver.info() != Eigen::Success) {
                tryFail = true;
                break;
            }
            X = solver.solve(B);
            if (solver.info() != Eigen::Success) {
                tryFail = true;
                break;
            }
            unsatisfied = false;
            for (unsigned int i = edgeStartIndex; i < edgeIndex; ++i) {

                if (X[i] < minEdgeLength || X[i] > maxEdgeLength) {
                    unsatisfied = true;
                    A.conservativeResize(rowsAdded+1, edgeIndex);
                    B.conservativeResize(rowsAdded+1);
                    A.insert(rowsAdded, i) = 1;
                    // B[rowsAdded++] = 1.2;
                    B[rowsAdded++] = edgeDist(gen);
                    // ++rowsAdded
                    break;
                }
            }
            if (rowsAdded > maxRows) {
                tryFail = true;
                break;
            }
        } while (unsatisfied);
        relativeError = (A*X - B).norm() / B.norm();
        if (relativeError > errorLimit || tryFail) {
            rowsAdded = startRowsAdded;
            A.conservativeResize(rowsAdded, edgeIndex);
            B.conservativeResize(rowsAdded);
            if (tries < maxTries - 1) {
                tryFail = false;
                if ((tries + 1) % increaseAt == 0) {
                    errorLimit *= increaseBy;
                    // increaseAt += increaseAtStep;
                }
            }
            // std::cout << tries << std::endl;
        } else {
            break;
        }
    }
    // if (relativeError < errorLimit) {
    std::cout << tries << " ";
    std::cout << relativeError << std::endl;
    // }


    // std::cout << X << std::endl;

    if (!tryFail) {
        oError = relativeError;
        for (unsigned int i = 0; i < this->primitives.size(); ++i) {
            glm::vec3 position(X[i * 3], X[i * 3 + 1], X[i * 3 + 2]);
            result.push_back(position);
        }
    }




    return result;



    // TODO idea: try eigen approach again but each time it generates set of values if any edge length constraints broken try adding a row that sets one edge to a random length, repeat
    // or see if can get basis and use that? https://stackoverflow.com/questions/54402199/what-is-the-best-way-to-solve-a-system-of-linear-equations-with-infinite-solutio

}

#endif


std::vector<glm::vec3> Graph::getCachedPositions()
{
    std::vector<glm::vec3> result;
    for (const uPtr<Primitive>& p : this->primitives) {
        if (p->cachedPos.has_value()) {
            result.push_back(p->cachedPos.value());
        } else {
            result.push_back(glm::vec3(0.f)); // TODO I guess should invoke samplepositions and return based on that in case where no value set
        }
    }
    return result;
}



// int main() {
//     // Create a CLP model object
//     ClpSimplex model;

//     // Define the problem dimensions
//     int numRows = 2;  // number of constraints
//     int numCols = 2;  // number of variables

//     // Set up the objective coefficients (e.g., minimize: 1x1 + 2x2)
//     double objective[] = {1.0, 2.0};  // coefficients for x1 and x2

//     // Set up the constraint matrix for the linear constraints:
//     // Ax <= b, where A is the constraint matrix, x is the variable vector, and b is the RHS vector
//     // Constraint 1: x1 + 2*x2 <= 4
//     // Constraint 2: x1 + x2 <= 5
//     double constraintMatrix[] = {1.0, 2.0,    // first row of coefficients (for constraint 1)
//                                  1.0, 1.0};   // second row of coefficients (for constraint 2)

//     // Right-hand side (RHS) values for the constraints
//     double rhs[] = {4.0, 5.0};  // RHS for constraint 1 and 2

//     // Set up the lower and upper bounds for the variables (x1, x2 >= 0)
//     double lowerBounds[] = {0.0, 0.0};
//     double upperBounds[] = {INFINITY, INFINITY};  // No upper bound for the variables

//     // Set up the model (objective function, constraints, and bounds)
//     model.loadProblem(numCols, numRows,
//     // model.loadProblem(numCols, numRows, constraintMatrix, nullptr, nullptr, objective, nullptr, nullptr,
//                       // lowerBounds, upperBounds, rhs, nullptr);

//     // Solve the problem using the primal simplex method
//     model.primal();

//     // Get and print the results (optimal values for x1 and x2)
//     double *solution = model.primalColumnSolution();
//     for (int i = 0; i < numCols; i++) {
//         printf("x%d = %f\n", i + 1, solution[i]);
//     }

//     return 0;
// }


