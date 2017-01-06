#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "util.h"

using namespace std;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Triangle_3 Triangle3D;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;


int main(int argc, char *argv[]) {
	if (argc != 3) {
		cerr << "Invalid Argument\n";
		return 1;
	}

	ifstream inputFile(argv[1]);
	ofstream outputFile(argv[2]);
	DT3 dt;

	vector<Point3D> points;

	if (!CGAL::read_xyz_points(inputFile, back_inserter(points))) { // output iterator over points
		cerr << "Error: cannot read file.";
		return 1;
	}

	sortAndRemoveDuplicate(points);

	dt.insert(points.begin(), points.end());

	if (!dt.is_valid(true)) {
		cerr << "Error: fail to build a Delaunay triangulation.\n";
		return 1;
	}
	if (dt.dimension() != 3) {
		cerr << "Error: cannot built a 3D triangulation.\n Current dimension = " << dt.dimension() << "\n";
		return 1;
	}

	vector<Segment3D> edges;
	vector<Triangle3D> faces;

	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++)
		edges.push_back(dt.segment(*edgeItr));

	for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++)
		faces.push_back(dt.triangle(*faceItr));

	sort(edges.begin(), edges.end(), [](const Segment3D &a, const Segment3D &b) -> bool { return a.squared_length() < b.squared_length(); });
	sort(faces.begin(), faces.end(), [](const Triangle3D &a, const Triangle3D &b) -> bool { return a.squared_area() < b.squared_area(); });

	outputFile << points.size() << "\n";
	for (Point3D point : points) {
		outputFile << point << "\n";
	}

	outputFile << edges.size() << "\n";
	for (Segment3D edge : edges) {
		outputFile << getIndex(points, edge[0]) << " " << getIndex(points, edge[1]) << "\n";
	}

	outputFile << faces.size() << "\n";
	for (Triangle3D triangle : faces) {
		outputFile << getIndex(points, triangle[0]) << " "
				   << getIndex(points, triangle[1]) << " "
				   << getIndex(points, triangle[2]) << "\n";
	}

	return 0;
}