#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <set>
#include <vector>
#include <QListWidgetItem>
#include "utils.h"

#define ANGLE_EPSILON 0.00001f
struct HalfEdgeGraph {
    // int id;
    float angle;
    int faces[2];

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
        return result;
    }

    bool operator==(const HalfEdgeGraph& t) const {
        return (this->angle <= t.angle + ANGLE_EPSILON &&
                this->angle >= t.angle - ANGLE_EPSILON &&
                this->faces[0] == t.faces[0] &&
                this->faces[1] == t.faces[1]);
    }

    bool operator<(const HalfEdgeGraph& t) const {
        if (this->angle == t.angle) {
            if (this->faces[0] == t.faces[0]) {
                return this->faces[1] < t.faces[1];
            } else {
                return this->faces[0] < t.faces[0];
            }
        } else {
            return (this->angle < t.angle);
        }
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

    bool sameHEs(const Primitive& p);

    // TODO need still to replace use of halfEdges and connections with new set<pair> structure
    // TODO maybe need to make a constructor that takes in the old form and makes the new one? since often in other code does work better making separately.
    //   or just need it to be a vector rather than a set but again other parts would work better sorted.
    //   or since the constructor in Graph works without needing that change, maybe just change how others structured? eh but the JSON stuff kinda needs all constructed before setting connections anyway
    // eh thinking on it more maybe leave as is and only do sorting in the generation of boundary string

};

#endif // PRIMITIVE_H
