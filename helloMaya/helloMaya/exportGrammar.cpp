#include "exportGrammar.h"
#include "LoadMesh.h"
#include "generateGrammar.h"
#include <QDir>
#include <QFileInfo>



void* ExportGrammarCmd::creator() {

	return new ExportGrammarCmd();
}


MStatus ExportGrammarCmd::doIt(const MArgList& args) {
    // Expect full file path (including name) from MEL
    if (args.length() < 1) {
        MGlobal::displayError("ExportGrammarCmd requires a file path argument.");
        return MS::kFailure;
    }

    MString filePathArg = args.asString(0);
    QString fullPath = QString::fromUtf8(filePathArg.asChar());

    QFileInfo fileInfo(fullPath);
    QDir dir = fileInfo.dir();
    if (!dir.exists()) {
        dir.mkpath(".");
    }

    // Prepare grammar data
    std::vector<std::array<Graph*, 2>> gram;
    for (auto& [first, second] : GenerateGrammarCmd::rules) {
        gram.push_back({ &first, &second });
    }

    // Write to file
    JSONReader::WriteGrammarToFile(fullPath, gram);

    MGlobal::displayInfo("Grammar exported to: " + MString(fullPath.toUtf8().constData()));
    return MS::kSuccess;
}
