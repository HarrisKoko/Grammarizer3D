#include "hello_maya.h"
#include <maya/MFnPlugin.h>
// define EXPORT for exporting dll functions
#define EXPORT _declspec(dllexport)
// Maya Plugin creator function
void* helloMaya::creator()
{
	return new helloMaya;
}
// Plugin doIt function
MStatus helloMaya::doIt(const MArgList& argList)
{
	MStatus status;
	MGlobal::displayInfo("Hello World!");
	// <<<your code goes here>>>
	MString name = "";
	MString id = "";
	for (int i = 0;i < argList.length();i++) {
		if (argList.asString(i) == "-name" || argList.asString(i) == "-n") {
			name = argList.asString(i + 1);
		}
		if (argList.asString(i) == "-id" || argList.asString(i) == "-i") {
			id = argList.asString(i + 1);
		}
	}
	MString str = "confirmDialog -title \"Hello Maya\" -message \"Name: " +
		name +
		"                     \\nID: " +
		id +
		"\" -button \"OK\";";
	MGlobal::executeCommand(str);

	return status;
}

EXPORT MStatus initializePlugin(MObject obj)
{
	MStatus status;
	MFnPlugin plugin(obj, "CIS660", "1.0", "Any");
	status = plugin.registerCommand("helloMaya", helloMaya::creator);
	if (!status)
		status.perror("registerCommand failed");
	return status;
}
// Cleanup Plugin upon unloading
EXPORT MStatus uninitializePlugin(MObject obj)
{
	MStatus status;
	MFnPlugin plugin(obj);
	status = plugin.deregisterCommand("helloMaya");
	if (!status)
		status.perror("deregisterCommand failed");
	return status;
}

