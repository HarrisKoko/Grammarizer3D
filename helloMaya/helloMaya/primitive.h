#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <set>
#include <vector>
#include <QListWidgetItem>
#include <glm/gtc/epsilon.hpp>
#include "utils.h"

#define THREE_DIMENSIONAL 1

#define ANGLE_EPSILON 0.005f
class HalfEdgeGraph {
public:
    // int id;
    float angle;
    std::array<int,2> faces;
    // int faces[2];

#if THREE_DIMENSIONAL
    // int volumes[2]; // in front of normal, behind normal // TODO actually I think other order might be more intuitive? inside then outside volume effectively
    // ^ optional
    // IGNORING FOR NOW, just going to do inside/outside volumes since no way to label volumes w/in maya mesh I think anyway
    // glm::vec3 faceNormal;
    std::array<glm::vec3,2> faceNormals = {{glm::vec3(0,1,0), glm::vec3(0,1,0)}};
    // glm::vec3 faceNormals[2];

    // HalfEdgeGraph& operator=(const HalfEdgeGraph& other) {
    //     if (this == &other)
    //         return *this;

    //     this->angle = other.angle;
    //     this->faces[0] = other.faces[0];
    //     this->faces[1] = other.faces[1];
    //     this->faceNormals[0] = other.faceNormals[0];
    //     this->faceNormals[1] = other.faceNormals[1];
    //     return *this;
    // }
#endif


    // TODO maybe ought to be in a cpp file instead
    HalfEdgeGraph operator-() const {
        HalfEdgeGraph result;
        if (this->angle < 0) {
            result.angle = this->angle + M_PI;
        } else {
            result.angle = this->angle - M_PI;
        }
        result.faces[0] = this->faces[1];
        result.faces[1] = this->faces[0];
#if THREE_DIMENSIONAL
        // result.volumes[0] = this->volumes[0]; // TODO think this is right?
        // result.volumes[1] = this->volumes[1];
        // result.faceNormal = this->faceNormal;
        // result.faceNormals[0] = this->faceNormals[0];
        // result.faceNormals[1] = this->faceNormals[1];
        result.faceNormals[0] = this->faceNormals[1];
        result.faceNormals[1] = this->faceNormals[0];
// TODO does this need two face normals stored? or do we make two HalfEdgeGraphs per outgoing edge?
//  easier with just always having two I think; can have edges with just one but can just represent with two of the same
// NOTE not going to account for case with more than two faces at a single edge (non-manifold) since not really useful for our purposes
// TODO unclear whether can just use a string for manifold meshes or if need to exapand to graph
// how do we represent turns? turns have facenormal? or do they have two?
// if in graph form, maybe:
// primitive contains vector of FaceGraphs which contain vector of HalfEdgeGraphs
// primitive's connections also then vector of vector of Primitive*
// then like a separate string per face? but I think they can have connectivity between faces matter so that doesn't quite work probably but can't think of exactly an example as to why
// MAYBE: all one string, but the substring comparison can skip over edges not w/ the same face normals? does that work?
// I think fine to have two face normals but need to always choose same one as basis for angle? but should carry same info from both knowing angle + two normals?
// need to decide on a consistent way to order normals I guess but then should work fine to always just base off of [0] for example?
//  oh wait but should probably have normals have same ordering as face labels in the arrays in these. so I guess not off of [0] but have some sorting of them which we use when finding angle
//  so whichever of [0] or [1] is considered "earlier" in that sorting (e.g. something like sort by increasing x then y then z would give -1,0,0 as basis for angle compared to sqrt(2)/2,sqrt(2)/2,0)
// if within one string, I guess can order HEs per primitive as such:
//  within same 2 faces & normals (which should only be when both totally on one face): sort as before, increasing angle
//  then the sets of HEs from that are sorted such that each shares 1 face (+normal) with adjacent set (in consistent order which I guess should be face to left (face 0?) matching face to right (face 1) on next HE)
// TODO in a manifold mesh is it possible to have edges go around in such a way that forces alternating between -bar edges (negative directions) and positive?
// might make sense convenience-wise to store faces as pair of integer and vec3 (label + normal)
// self: see paper notes. the branch gluing seems fine so far in just string but not yet sure if breaks down when loop gluing
#endif
        return result;
    }



    bool operator==(const HalfEdgeGraph& t) const {
        /*if (this->angle + ANGLE_EPSILON > M_PI) {

        }
        else {

        }*/
        const float PI2 = M_PI * 2;
        return (

            ((this->angle > M_PI - ANGLE_EPSILON * 2 && t.angle < 0) ? glm::epsilonEqual(this->angle - PI2, t.angle, ANGLE_EPSILON) :
                ((t.angle > M_PI - ANGLE_EPSILON * 2 && this->angle < 0) ? glm::epsilonEqual(this->angle, t.angle - PI2, ANGLE_EPSILON) :
                    glm::epsilonEqual(this->angle, t.angle, ANGLE_EPSILON))) &&
            //abs(t.angle - this->angle) >= PI2 - ANGLE_EPSILON ? 
            //    (t.angle > this->angle ? 
            //        glm::epsilonEqual(this->angle, t.angle - PI2, ANGLE_EPSILON) :
            //        glm::epsilonEqual(this->angle - PI2, t.angle, ANGLE_EPSILON))
            //        : glm::epsilonEqual(this->angle, t.angle, ANGLE_EPSILON) &&
        //return (glm::epsilonEqual((this->angle + ANGLE_EPSILON > M_PI ? this->angle - 2* M_PI : this->angle, t.angle + ANGLE_EPSILON > M_PI ? t.angle - 2 * M_PI : t.angle, ANGLE_EPSILON) &&
        //return glm::epsilonEqual(this->angle, t.angle, ANGLE_EPSILON) &&
#if THREE_DIMENSIONAL
                // this->volumes[0] == t.volumes[0] &&
                // this->volumes[1] == t.volumes[1] &&
            glm::all(glm::epsilonEqual(this->faceNormals[0], t.faceNormals[0], ANGLE_EPSILON)) &&
            glm::all(glm::epsilonEqual(this->faceNormals[1], t.faceNormals[1], ANGLE_EPSILON)) &&
#endif
                this->faces[0] == t.faces[0] &&
                this->faces[1] == t.faces[1]);
    }

    bool operator<(const HalfEdgeGraph& t) const {
        // TODO figure out how to handle volume/normal? this is mostly just for sorting which is plane-dependent
#if THREE_DIMENSIONAL
        //if (this->faceNormals[0] == t.faceNormals[0]) {
            //if (this->faceNormals[1] == t.faceNormals[1]) {
                //if (this->angle <= t.angle + ANGLE_EPSILON &&
                    //this->angle >= t.angle - ANGLE_EPSILON) {
        const float PI2 = M_PI * 2;

        // TODO I think this might just not work in maps for end cases. maybe pre-round so that the values being used always cap out at pi - epsilon? MORNING TODO
        if (glm::all(glm::epsilonEqual(this->faceNormals[0], t.faceNormals[0], ANGLE_EPSILON))) {
            if (glm::all(glm::epsilonEqual(this->faceNormals[1], t.faceNormals[1], ANGLE_EPSILON))) {
                //if (glm::epsilonEqual(this->angle, t.angle, ANGLE_EPSILON)) {
                if ((this->angle > M_PI - ANGLE_EPSILON * 2 && t.angle < 0) ? glm::epsilonEqual(this->angle - PI2, t.angle, ANGLE_EPSILON) :
                    ((t.angle > M_PI - ANGLE_EPSILON * 2 && this->angle < 0) ? glm::epsilonEqual(this->angle, t.angle - PI2, ANGLE_EPSILON) :
                        glm::epsilonEqual(this->angle, t.angle, ANGLE_EPSILON))) {
                //if (abs(t.angle - this->angle) >= PI2 - ANGLE_EPSILON ?
                //    (t.angle > this->angle ?
                //        glm::epsilonEqual(this->angle, t.angle - PI2, ANGLE_EPSILON) :
                //        glm::epsilonEqual(this->angle - PI2, t.angle, ANGLE_EPSILON))
                //    : glm::epsilonEqual(this->angle, t.angle, ANGLE_EPSILON)) {
                    if (this->faces[0] == t.faces[0]) {
                        return this->faces[1] < t.faces[1];
                    } else {
                        return this->faces[0] < t.faces[0];
                    }
                } else {
                    return ((this->angle > M_PI - ANGLE_EPSILON * 2 && t.angle < 0) ?
                        this->angle - PI2 + ANGLE_EPSILON < t.angle :
                        ((t.angle > M_PI - ANGLE_EPSILON * 2 && this->angle < 0) ?
                            this->angle + PI2 + ANGLE_EPSILON < t.angle :
                            this->angle + ANGLE_EPSILON < t.angle));
                    //return (abs(t.angle - this->angle) >= PI2 - ANGLE_EPSILON ?
                    //    (t.angle > this->angle ?
                    //        (this->angle + ANGLE_EPSILON < t.angle - PI2) :
                    //        (this->angle - PI2 + ANGLE_EPSILON < t.angle))
                    //    : (this->angle + ANGLE_EPSILON < t.angle));
                    //return (this->angle + ANGLE_EPSILON < t.angle);
                }
            } else {
                if (glm::epsilonEqual(this->faceNormals[1].x, t.faceNormals[1].x, ANGLE_EPSILON)) {
                    if (glm::epsilonEqual(this->faceNormals[1].y, t.faceNormals[1].y, ANGLE_EPSILON)) {
                        return (this->faceNormals[1].z + ANGLE_EPSILON < t.faceNormals[1].z);
                    }
                    return (this->faceNormals[1].y + ANGLE_EPSILON < t.faceNormals[1].y);
                }
                return (this->faceNormals[1].x + ANGLE_EPSILON < t.faceNormals[1].x);

            }
        }
        else {
            if (glm::epsilonEqual(this->faceNormals[0].x, t.faceNormals[0].x, ANGLE_EPSILON)) {
                if (glm::epsilonEqual(this->faceNormals[0].y, t.faceNormals[0].y, ANGLE_EPSILON)) {
                    return (this->faceNormals[0].z + ANGLE_EPSILON < t.faceNormals[0].z);
                }
                return (this->faceNormals[0].y + ANGLE_EPSILON < t.faceNormals[0].y);
            }
            return (this->faceNormals[0].x + ANGLE_EPSILON < t.faceNormals[0].x);
        }
#else
        if (this->angle == t.angle) {
            if (this->faces[0] == t.faces[0]) {
                return this->faces[1] < t.faces[1];
            } else {
                return this->faces[0] < t.faces[0];
            }
        } else {
            return (this->angle < t.angle);
        }
#endif
    }
    // TODO maybe just make this into a whole class now that I'm defining more stuff in it? though it's just comparison operators
};

class Primitive : public QListWidgetItem
{
private:
    static unsigned int lastID;
public:
    Primitive();
    unsigned int ID;
    std::vector<HalfEdgeGraph> halfEdges;
    std::vector<Primitive*> connections;

    static void resetIDs();

    // std::vector<std::pair<HalfEdgeGraph, Primitive*>> halfEdges;
    //  TODO try using this way instead
    //  could even be a std::set then maybe?
    // std::set<std::pair<HalfEdgeGraph, Primitive*>> halfEdgeData;
    //  want ordered but also want to refer to by index I think and not sure if can do that with sets
    //      so perhaps instead use a vector but sort it?
    //      or wait maybe don't need to use index. just use pointers to the stored pairs
    //          then boundary string elements: pointer to primitive, pointer to pair
    //          when finding a match in replacement attempt, can find corresponding primitive and its connection/lack thereof from other subgraph by following pointer from same spot in its 'string'

    bool sameHEs(const Primitive& p) const;

    // TODO need still to replace use of halfEdges and connections with new set<pair> structure
    // TODO maybe need to make a constructor that takes in the old form and makes the new one? since often in other code does work better making separately.
    //   or just need it to be a vector rather than a set but again other parts would work better sorted.
    //   or since the constructor in Graph works without needing that change, maybe just change how others structured? eh but the JSON stuff kinda needs all constructed before setting connections anyway
    // eh thinking on it more maybe leave as is and only do sorting in the generation of boundary string

    std::optional<glm::vec3> cachedPos = {};

};

#endif // PRIMITIVE_H
