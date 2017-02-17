#define _USE_MATH_DEFINES
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Union_find.h>

#include "util.h"

// Alias for CGAL datatypes
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Vector_3 Vector3D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Triangle_3 Triangle3D;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;
const double inf = std::numeric_limits<double>::infinity();

int main(int argc, char *argv[]) {
	if (argc != 2) {
		cerr << "Invalid Argument\n";
		return 1;
	}

	println("Input File =", argv[1]);

	std::ifstream inputFile(argv[1]);

	vector<Point3D> pts;
	DT3 dt;

	auto start = std::chrono::high_resolution_clock::now();
	if (!CGAL::read_xyz_points(inputFile, back_inserter(pts))) { // output iterator over points
		cerr << "Error: cannot read file.";
		return 1;
	}
	sortAndRemoveDuplicate(pts);
	pts.shrink_to_fit();

	auto finish = std::chrono::high_resolution_clock::now();
	cout << pts.size() << " pts read in " << std::chrono::duration<double>(finish - start).count() << " secs\n";
	start = std::chrono::high_resolution_clock::now();
	dt.insert(pts.begin(), pts.end());
	if (!dt.is_valid(true)) {
		cerr << "Error: fail to build a Delaunay triangulation.\n";
		return 1;
	}
	if (dt.dimension() != 3) {
		cerr << "Error: cannot built a 3D triangulation.\n Current dimension = " << dt.dimension() << "\n";
		return 1;
	}
	finish = std::chrono::high_resolution_clock::now();
	cout << "Delaunay Triangulation created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	/*vector<ourFace> allFaces;
	for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
		Triangle3D tri = dt.triangle(*faceItr);
		allFaces.push_back(ourFace(getIndex(pts, tri[0]),
								   getIndex(pts, tri[1]),
								   getIndex(pts, tri[2])));
	}
	allFaces.shrink_to_fit();
	std::sort(allFaces.begin(), allFaces.end());*/

	vector<set<int>> adjList(pts.size());

	//~ for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
		//~ Segment3D seg = dt.segment(*edgeItr);
		//~ int u = getIndex(pts, seg[0]), v = getIndex(pts, seg[1]);
		//~ adjList[u].insert(v);
		//~ adjList[v].insert(u);
	//~ }
	
	for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
		Triangle3D tri = dt.triangle(*faceItr);
		
		int u = getIndex(pts, tri[0]);
		int v = getIndex(pts, tri[1]);
		int w = getIndex(pts, tri[2]);
		
		adjList[u].insert(v);
		adjList[v].insert(u);
		
		adjList[u].insert(w);
		adjList[w].insert(u);
		
		adjList[w].insert(v);
		adjList[v].insert(w);
								 
	}

	/*vector<DT3::Vertex_handle> ptsHandle;
	ptsHandle.resize(pts.size());
	for (auto itr = dt.finite_vertices_begin(); itr != dt.finite_vertices_end(); itr++) {
		ptsHandle[getIndex(pts, itr->point())] = itr;
	}

	DT3::Cell_handle cell;
	int u, v, w;*/

	cout << std::boolalpha;
	for (int i = 0; i < pts.size(); ++i) {
		for (int j = i + 1; j < pts.size(); ++j) {
			for (int k = j + 1; k < pts.size(); ++k) {
				cout << i << " " << j << " " << k<<" ";
				cout << (adjList[i].find(j) != adjList[i].end() &&
						 adjList[i].find(k) != adjList[i].end() &&
						 adjList[j].find(k) != adjList[j].end())
					 << "\n";
				// cout << std::binary_search(allFaces.begin(), allFaces.end(), ourFace(i, j, k)) << "\n";
				// cout << dt.is_facet(ptsHandle[i], ptsHandle[j], ptsHandle[k], cell, u, v, w) << "\n";
			}
		}
	}

	// vector<ourFace> allFaces;
	// for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
	// 	Triangle3D tri = dt.triangle(*faceItr);
	// 	allFaces.push_back(ourFace(getIndex(pts, tri[0]),
	// 							   getIndex(pts, tri[1]),
	// 							   getIndex(pts, tri[2])));
	// }
	// allFaces.shrink_to_fit();
	// std::sort( allFaces.begin(), allFaces.end());

	return 0;
}
