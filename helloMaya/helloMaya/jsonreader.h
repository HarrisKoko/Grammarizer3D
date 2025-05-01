#ifndef JSONREADER_H
#define JSONREADER_H

#include "graph.h"
#include "utils.h"
#include <QJsonObject>
#include <QJsonArray>
#include <QJsonDocument>
#include <QFile>
#include <QList>

class JSONReader
{
public:
    JSONReader();

    static bool LoadGeometryFromFile(const QString &local_path, std::vector<glm::vec3>* oVertPos, std::vector<std::pair<std::vector<int>,int>>* oFaces, std::vector<glm::vec3>* oFaceColors);
    static void extracted(Graph &result, QJsonArray &boundary);
    static Graph LoadGraphFromFile(const QString &local_path);

    static bool WriteGraphToFile(const QString &local_path, const Graph *graph);
    static bool WriteGrammarToFile(const QString &local_path, const std::vector<std::array<Graph*,2>> &graphGrammar);
    static std::vector<std::pair<Graph, Graph>> LoadGrammarFromFile(const QString &local_path);
private:
    static void jsonToGraph(Graph *result, QJsonObject json);
    static bool graphToJSON(QJsonObject* json, const Graph* graph);
};

#endif // JSONREADER_H
