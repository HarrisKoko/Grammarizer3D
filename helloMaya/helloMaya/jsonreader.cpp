#include "jsonreader.h"

JSONReader::JSONReader() {}

bool JSONReader::LoadGeometryFromFile(const QString &local_path,
                              std::vector<glm::vec3> *oVertPos,
                              std::vector<std::pair<std::vector<int>, int> > *oFaces,
                              std::vector<glm::vec3> *oFaceColors)
{
    QFile file(local_path);
    if (file.open(QIODevice::ReadOnly)) {
        QByteArray rawData = file.readAll();
        QJsonDocument doc(QJsonDocument::fromJson(rawData));
        QJsonObject json = doc.object();
        QJsonArray vertices = json["vertices"].toArray();
        QJsonArray faces = json["faces"].toArray();
        QJsonArray faceColors = json["faceColors"].toArray();
        // TODO maybe add checks for if not present but I think actually kinda fine as is?
        // since e.g. toArray() gives an empty array when the member isn't present that should just behave fine


        QJsonArray vArr;
        for (const QJsonValue &v : vertices) {
            vArr = v.toArray();
            oVertPos->push_back(glm::vec3(vArr.at(0).toDouble(),
                                          vArr.at(1).toDouble(),
                                          vArr.at(2).toDouble()));
        }

        QJsonObject fData;
        QJsonArray iArr;
        for (const QJsonValue &f : faces) {
            fData = f.toObject();
            std::pair<std::vector<int>, int> newFace;
            newFace.second = fData["label"].toInt(-1);
            iArr = fData["indices"].toArray();
            for (const QJsonValue &i : std::as_const(iArr)) {
                newFace.first.push_back(i.toInt());
            }
            oFaces->push_back(newFace);
        }

        QJsonArray cArr;
        for (const QJsonValue &fc : faceColors) {
            cArr = fc.toArray();
            oFaceColors->push_back(glm::vec3(cArr.at(0).toDouble(),
                                             cArr.at(1).toDouble(),
                                             cArr.at(2).toDouble()));
        }

        return true;
    }
    return false;
}

void JSONReader::jsonToGraph(Graph* result, QJsonObject json)
{
    QJsonArray prims = json["primitives"].toArray();

    for (const QJsonValue& p : prims) {
        uPtr<Primitive> pPtr = mkU<Primitive>();
        QJsonObject pObj = p.toObject();
        QJsonArray hes = pObj["halfEdges"].toArray();
        for (const QJsonValue& he : hes) {
            QJsonObject heObj = he.toObject();
            HalfEdgeGraph heResult;
            heResult.angle = heObj["angle"].toDouble();
            QJsonArray heFaces = heObj["faces"].toArray();
            heResult.faces[0] = heFaces.at(0).toInt();
            heResult.faces[1] = heFaces.at(1).toInt();
            pPtr->halfEdges.push_back(heResult);
        }
        result->primitives.push_back(std::move(pPtr));
    }

    for (unsigned int i = 0; i < prims.size(); ++i) {
        QJsonObject pObj = prims.at(i).toObject();
        Primitive* pPtr = result->primitives.at(i).get();

        QJsonArray connectionIndices = pObj["connections"].toArray();
        for (const QJsonValue& ci : connectionIndices) {
            int val = ci.toInt();
            if (val < 0) {
                pPtr->connections.push_back(nullptr);
            } else {
                pPtr->connections.push_back(result->primitives.at(val).get());
            }
        }
    }

    QJsonArray boundary = json["boundaryString"].toArray();
    for (const QJsonValue& elem : boundary) {
        QJsonObject elemObj = elem.toObject();
        QString elemType = elemObj["type"].toString();
        if (elemType == "turn") {
            Turn val = Turn(elemObj["value"].toInt());
            result->boundaryString.push_back(val);
            // TODO check that works
        } else if (elemType == "halfEdge") {
            QJsonObject heObj = elemObj["h"].toObject();
            BoundaryHE b;
            b.primIndex = elemObj["primIndex"].toInt();
            QJsonArray heFaces = heObj["faces"].toArray();
            b.h.angle = heObj["angle"].toDouble();
            b.h.faces[0] = heFaces[0].toInt();
            b.h.faces[1] = heFaces[1].toInt();
            result->boundaryString.push_back(b);
        } else {
            // invalid TODO dunno how/if wanna handle
        }
    }

    result->isSingleEdge = json["isSingleEdge"].toBool(false);
}

Graph JSONReader::LoadGraphFromFile(const QString &local_path)
{
    // TODO not sure if makes more sense to pass out by pointer argument rather than as return type?
    QFile file(local_path);
    Graph result;
    if (file.open(QIODevice::ReadOnly)) {
        QByteArray rawData = file.readAll();
        QJsonDocument doc(QJsonDocument::fromJson(rawData));
        QJsonObject json = doc.object();
        jsonToGraph(&result, json);
    }
    return result;
}

// TODO could save positional data too? for saving during grammar applying steps
// TODO rn:
//  add GUI button; try to add positional data to graph class but maybe not to this part yet; try to add intersection constraint
//      done;   done (but not this part yet);   done for random case but not done for eigen version
bool JSONReader::graphToJSON(QJsonObject* json, const Graph* graph) {
    QJsonArray prims;
    for (unsigned int i = 0; i < graph->primitives.size(); ++i) {

        QJsonArray halfEdgesArr;
        QJsonArray connectionsArr;
        // TODO need to store connections as integer indices I think cuz can't do pointers per se in JSON. not sure if there's a more efficient way to do but I think just iterate over all the prims every time and check which it equals?
        //  not sure if should first make a separate vector of raw pointers or just use .get every time? I think in C++23 on latter would be better (since constexpr anyway)? efficiency doesn't really matter I guess though
        Primitive* curPrim = graph->primitives.at(i).get();
        for (const HalfEdgeGraph &curHalfEdge : curPrim->halfEdges) {
            QJsonObject heObj({{"angle", curHalfEdge.angle}, {"faces", QJsonArray{curHalfEdge.faces[0], curHalfEdge.faces[1]}}});
            halfEdgesArr.push_back(heObj);
        }
        for (const Primitive* connect : curPrim->connections) {
            if (connect == nullptr) {
                connectionsArr.push_back(-1);
            } else {
                bool added = false;
                for (unsigned int j = 0; j < graph->primitives.size(); ++j) {
                    if (connect == graph->primitives.at(j).get()) {
                        connectionsArr.push_back(QJsonValue(int(j)));
                        added = true;
                        break;
                    }
                }
                if (!added) {
                    return false; // failure to construct JSON so returning false
                }
            }
        }
        QJsonObject primObj({{"halfEdges", halfEdgesArr}, {"connections", connectionsArr}});
        prims.push_back(primObj);
    }
    json->insert("primitives", prims);


    QJsonArray boundary;
    for (unsigned int i = 0; i < graph->boundaryString.size(); ++i) {
        std::variant<BoundaryHE, Turn> elem = graph->boundaryString.at(i);
        Turn* t = std::get_if<Turn>(&elem);
        if (t != nullptr) {
            // Turn
            QJsonObject turnObj({{"type", "turn"}, {"value", *t}});
            boundary.push_back(turnObj);
        } else {
            // Half edge
            BoundaryHE* he = std::get_if<BoundaryHE>(&elem);
            if (he == nullptr) {
                return false; // invalid boundary string
            }
            QJsonObject heObj({{"angle", he->h.angle}, {"faces", QJsonArray{he->h.faces[0], he->h.faces[1]}}});
            QJsonObject elemObj({{"type", "halfEdge"}, {"h", heObj}, {"primIndex", int(he->primIndex)}});
            boundary.push_back(elemObj);
        }
    }
    json->insert("boundaryString", boundary);

    json->insert("isSingleEdge", graph->isSingleEdge);
    return true;
}

bool JSONReader::WriteGraphToFile(const QString &local_path, const Graph *graph)
{
    // TODO will need to make variant that takes several graphs later for writing whole grammars (array of arrays of 2 graphs); basically just this same code but an extra few layers after opening file
    QFile file(local_path);
    if (file.open(QIODevice::WriteOnly)) {
        // QJsonDocument doc{};
        // TODO
        // QJsonObject json = doc.object();
        QJsonObject json;
        if (!graphToJSON(&json, graph)) {
            return false;
        }

        QJsonDocument doc(json);
        file.write(doc.toJson());
        return true;
    }
    return false;
}

bool JSONReader::WriteGrammarToFile(const QString &local_path, const std::vector<std::array<Graph *, 2> > &graphGrammar)
{
    //TODO UPDATE
    // TODO probably need graphs in grammar to contain boundary strings? or at least some other sort of correspondence to let know which boundary edges correspond to which others in pair graph
    // TODO Need to test this and write reader but might figure out more first how grammar will be used in rest of code
    QFile file(local_path);
    if (file.open(QIODevice::WriteOnly)) {
        QJsonObject json;
        QJsonArray rules;
        for (const std::array<Graph*, 2> &graphPair : graphGrammar) {
            QJsonArray pair;
            for (const Graph* graph : graphPair) {
                QJsonObject graphObj;
                if (!graphToJSON(&graphObj, graph)) {
                    return false;
                }
                pair.push_back(graphObj);
            }
            rules.push_back(pair);
        }
        json.insert("rules", rules);
        QJsonDocument doc(json);
        file.write(doc.toJson());
        return true;
    }
    return false;
}


std::vector<std::pair<Graph,Graph>> JSONReader::LoadGrammarFromFile(const QString &local_path)
{
    // TODO not sure exact best return type? probably should just make a grammar class that stores this
    // NOTE currently this uses arrays of 2 graphs while other uses arrays of 2 graph pointers so inconsistent
    //      now uses a pair of graphs which still isn't same
    // TODO now uses same type as the random rule applying function but I HAVE NO IDEA if that's a reasonable way to pass it out (seems like might have some issues with clearing memory just intuitively speaking? but not sure)
    QFile file(local_path);
    std::vector<std::pair<Graph,Graph>> result;
    // std::vector<std::array<Graph,2>> result;
    // Graph result;
    if (file.open(QIODevice::ReadOnly)) {
        QByteArray rawData = file.readAll();
        QJsonDocument doc(QJsonDocument::fromJson(rawData));
        // QJsonObject json = doc.object();
        // jsonToGraph(&result, json);
        QJsonObject json = doc.object();
        QJsonArray rules = json["rules"].toArray();
        for (const QJsonValue &pair : rules) {
            QJsonArray pairArr = pair.toArray();
            // std::array<Graph,2> graphPair;
            // for (unsigned int i = 0; i < 2; ++i) {
                // for (const QJsonValue &graphVal : pairArr) {

            // Graph g0, g1;
            QJsonObject graphObj0 = pairArr.at(0).toObject();
            QJsonObject graphObj1 = pairArr.at(1).toObject();
            std::pair<Graph,Graph> gPair;
            jsonToGraph(&gPair.first, graphObj0);
            jsonToGraph(&gPair.second, graphObj1);

            // }
            result.push_back(gPair);
        }
    }
    return result;
}


