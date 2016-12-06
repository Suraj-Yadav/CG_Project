cd /D C:\ExtLibs\CGAL-4.9
cd /D G:\work\CG_Project\build
cmake -G "Visual Studio 14 2015 Win64"
MSBuild /p:Configuration=Release 
RMDIR /S /Q .