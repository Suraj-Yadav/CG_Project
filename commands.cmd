cd /D C:\ExtLibs\CGAL-4.9

cd /D G:\work\CG_Project\build
"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" amd64
MSBuild /p:Configuration=Release test3D.vcxproj

cmake -G "Visual Studio 14 2015 Win64"

forfiles /P ..\input_3D /M *.xyz /c "cmd /C ..\build\Release\surface3D.exe @file ..\build\output.txt" > log.txt 2>&1

RMDIR /S /Q .