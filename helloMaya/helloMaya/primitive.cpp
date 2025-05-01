#include "primitive.h"

unsigned int Primitive::lastID{0};

Primitive::Primitive()
    :ID(lastID++)
{
    setText(QString::number(ID));
}

void Primitive::resetIDs(){
    Primitive::lastID = 0;
}

bool Primitive::sameHEs(const Primitive &p) const
{
    // return std::equal(this->halfEdgeData.cbegin(),
    //                   this->halfEdgeData.cend(),
    //                   p.halfEdgeData.cbegin(),
    //                   [](const auto& a, const auto& b) {
    //                       return a.first == b.first;
    //                   });

    if (this->halfEdges.size() != p.halfEdges.size()) {
        return false;
    }
    // note doesn't check that there's the same number of each HE but that shouldn't come up in any properly constructed primitive anyway so not really a concern
    for (const HalfEdgeGraph &he : this->halfEdges) {
        bool noMatch = true;
        for (const HalfEdgeGraph &he2 : p.halfEdges) {
            if (he == he2) {
                noMatch = false;
                break;
            }
        }
        if (noMatch) {
            return false;
        }
    }
    return true;
}
