#include<bits/stdc++.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Union_find.h>

using namespace std;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point2D;
typedef Kernel::Vector_2 Vector2D;
typedef Kernel::Segment_2 Segment2D;
typedef Kernel::Triangle_2 Triangle2D;
typedef CGAL::Delaunay_triangulation_2<Kernel> DT2;

int main(int argc, char *argv[]){

    if (argc != 3) {
		cerr << "Invalid Argument\n";
		return 1;
	}

	std::ifstream inputFile(argv[1]);
	std::ofstream outputFile(argv[2]);

    vector<Point2D> a;
    map<Point2D, int> index;
    Point2D p;
    int n = 0;

    while(inputFile >> p){
        // ignore whatever comes after x and y
		inputFile.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
		a.push_back(p);
        index[p] = n;
        n++;
    }

    DT2 delaunay;

    delaunay.insert(a.begin(),a.end());

    vector<pair<int,int>> allEdges;

    for(auto edgeItr = delaunay.finite_edges_begin(); edgeItr != delaunay.finite_edges_end(); edgeItr++){
        Segment2D seg = delaunay.segment(*edgeItr);
        allEdges.push_back(make_pair(index[seg[0]],index[seg[1]]));
    }

    outputFile << a.size() <<"\n";
    for(auto x: a){
        outputFile << x << " 0.0\n";
    }

    outputFile << allEdges.size() <<"\n";
    for(auto x: allEdges){
        outputFile << x.first << " " << x.second << "\n";
    }

    outputFile << 0 <<"\n";

    return 0;
}

