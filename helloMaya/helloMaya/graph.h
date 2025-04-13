#ifndef GRAPH_H
#define GRAPH_H

#include "primitive.h"
#include "utils.h"

struct BoundaryHE {
    // IDEA: store pointer to primitive and index in its halfEdges vector rather than pointer directly to HE
    //  better allows manipulation of connections based on boundary string
    //  equality tests of boundary strings still done by looking at just HE itself
    //  TODO maybe instead store Primitive* p and pointer to HalfEdgeGraph then find in p when need connection index? eh
    // ah, but how does one update pointer when combining strings then? hmm
    //  opposite sort of approach: store half edge itself and an index of the primitive in the primitives vector
    //  then when combining graphs can simply add on the length of the other graph's primitive vector for one of them, then concatenate vectors.
    // Primitive* p;
    // unsigned int index;
    // std::pair<HalfEdgeGraph, Primitive*> he;
    HalfEdgeGraph h;
    unsigned int primIndex;

    bool operator==(const BoundaryHE& other) const {
        return (h == other.h);
    }
};

enum Turn {positive, negative};


class Graph : public QListWidgetItem
{
private:
    static unsigned int lastID;
public:
    Graph();
    unsigned int ID;
    static void resetIDs();
    Graph(std::vector<glm::vec3>, std::vector<std::pair<std::vector<int>,int>>);

    Graph(const Graph& other);

    Graph& operator=(const Graph& other);
    // Graph(uPtr<Primitive>);
    Graph(const std::vector<HalfEdgeGraph> &primHEs);

    std::vector<uPtr<Primitive>> primitives;
    // std::vector<> boundary_strings // TODO
    //  I think perhaps boundary strings make most sense to do as a vector of 'HalfEdgeGraph's (or vector of pointers to them?)? except also needs some special case term that refers to the turns (can I think make a struct inherit so use polymorphism with a separate class for turns?)
    //  then constructing a primitive's boundary string: just iterate through half-edges and sort by direction
    //      more complex graph's strings given by glue operations of the primitive's strings
    std::vector<std::variant<BoundaryHE, Turn>> boundaryString;

    std::vector<glm::vec3> face_colors;

    bool isSingleEdge = false;
    // TODO we might want to figure out if any of these should be private? but as it is it's more convenient this way since we do access them from other places some. boundary string could probably be private

    // TODO see note below should probably return separate graph
    bool applyReplacementRule(const Graph& from, const Graph& to);

    // TODO might want to change how this is passed in later
    // TODO note at the moment this is directed first->second but should actually do both, not sure if want to do within this function (just the same thing but try both first,second and second,first in a random order) or construct that elsewhere then pass in
    // TODO should probably change to return a separate graph rather than edit the calling one now that the copier works right
    bool applyRandomReplacementRule(const std::vector<std::pair<Graph, Graph> > &grammar, bool bidirectional=true, bool skipStarters=false);

    // TODO probably need to add special case graph: 0 primitives, just an edge (or rather two half edges in opposite direction)
    //  (generation 0 of hierarchy)

    std::vector<Graph> branchGlue(const Graph& other) const;
    std::vector<Graph> loopGlue() const;

    void cancelTurns();

    // TODO not sure really if this makes more sense as an overriding operator or a separate sort of function. checking for graph isomorphism
    //  if this is == should probably also check face_colors is the same?
    // bool operator==(const Graph& other);
    bool hasSubgraph(const Graph& other) const;
    bool isIsomorphicTo(const Graph& other) const;
    bool sameBoundaryString(const Graph& other) const;

    static std::vector<std::pair<Graph, Graph> > generateRules(const std::vector<Primitive *> &primitives, std::vector<glm::vec3> faceColors);
    std::vector<std::pair<Graph, Graph>> generateRules() const;
    // std::vector<Graph> splice();


    // just for testing, since in main algorithm we want to interject in the hierarchy construction rather than fully construct
    static std::vector<Graph> generateHierarchy(const std::vector<Primitive*> &primitives, int steps, std::vector<glm::vec3> faceColors);
    std::vector<Graph> generateHierarchy(int steps) const;


    std::vector<glm::vec3> samplePositions(float minEdgeLength = 1, float maxEdgeLength = 20, float minPosition = -10, float maxPosition = 10) const;
    // TODO probably store result between iterations and update

};

#endif // GRAPH_H
