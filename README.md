# Grammarizer3D
This repository is for the 3D version of the grammarizer tool for Maya. This tool allows artists to use their models in Maya as a basis for a procedural grammar. They can then use the plugin to generate proceudral variations of their model. 

# Running this code
For this code to run properly, you will need to install Qt here https://doc.qt.io/qt-6/get-and-install-qt.html. Upon running the installer, select Qt 6.9 msvc 2022 64 bit and qQ 6.9 additional libraries -> multimedia. Once this has installed, you must add a system environment variable to your path for the installed msvc bin folder (Ex. C:\Qt\6.9.0\msvc2022_64\bin). This allows the Maya plugin (.mll file) to access Qt at runtime. It has been tested primarily with Maya 2022 which can be found here https://www.autodesk.com/products/maya/overview.

Once the required dependencies have been downloaded, please ensure your visual studio solution file is updated with the correct paths. The first section that needs updating is the C++ Additional Includes which needs:
1. Maya 2022 include (Ex. C:\Program Files\Autodesk\Maya2022\include)
2. Grammarizer include folder (Ex. C:\Users\Harris\Downloads\Grammarizer3D\helloMaya\include)
3. QtGui, QtWidgets, Qtcore, msvc include (Ex. C:\Qt\6.9.0\msvc2022_64\include\QtCore)
4. MSVC include (Ex. C:\Qt\6.9.0\msvc2022_64\include)

Then, the Linker's Additional Library Directories needs:
1. MSVC Library (Ex. C:\Qt\6.9.0\msvc2022_64\lib)
2. Maya 2022 Library (Ex. C:\Program Files\Autodesk\Maya2022\lib)

Once all of this has been completed, the code can be built. Upon a succesful build, the .mll plugin file is sent to helloMaya/x64/debug/helloMaya.mll or helloMaya/x64/release/helloMaya.mll. This can now be imported into Maya using the plugin window. 
