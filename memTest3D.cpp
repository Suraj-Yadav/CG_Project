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

class Trie {
  private:
	struct node {
		std::unordered_map<int, node*> mapping;
		std::unordered_map<int, bool> check;
	};
	node root;
	std::vector<int> count;
	void insert(node *x, int u, int v, int w) {
		if (x->check[u] == false){
			x->mapping[u] = new node;
			x->check[u]=true;
		}
		x = x->mapping[u];
		if (x->check[v] == false){
			x->check[v]=true;
			x->mapping[v] = new node;
		}
		x = x->mapping[v];
		x->check[w] = true;
		/*if (x->mapping[w] == nullptr)
			x->mapping[w] = new node;*/
	}
	bool find(node *x, int u, int v, int w) {
		if (x->check[u] == false)
			return false;
		x = x->mapping[u];
		if (x->check[v] == false)
			return false;
		x = x->mapping[v];
		if (x->check[w] == false)
			return false;
		return true;
	}
  public:
	Trie(const DT3 &dt, const vector<Point3D> &pts) {
		count.resize(dt.number_of_vertices(), 0);
		for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
			Triangle3D tri = dt.triangle(*faceItr);
			int u = getIndex(pts, tri[0]);
			int v = getIndex(pts, tri[1]);
			int w = getIndex(pts, tri[2]);
			count[u]++;
			count[v]++;
			count[w]++;
		}
		for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
			Triangle3D tri = dt.triangle(*faceItr);
			int u = getIndex(pts, tri[0]);
			int v = getIndex(pts, tri[1]);
			int w = getIndex(pts, tri[2]);
			if (u > v) std::swap(u, v);
			if (v > w) std::swap(v, w);
			if (u > v) std::swap(u, v);

			if (count[u] < count[v]) std::swap(u, v);
			if (count[v] < count[w]) std::swap(v, w);
			if (count[u] < count[v]) std::swap(u, v);
			insert(&root, u, v, w);
		}
	}
	bool find(int u, int v, int w) {
		if (u > v) std::swap(u, v);
			if (v > w) std::swap(v, w);
			if (u > v) std::swap(u, v);
		if (count[u] < count[v]) std::swap(u, v);
		if (count[v] < count[w]) std::swap(v, w);
		if (count[u] < count[v]) std::swap(u, v);
		return find(&root, u, v, w);
		// return find(&root, u, v, w)||find(&root, u, w, v)||find(&root, w, v, u)||find(&root, w, u, v)||find(&root, v, u, w)||find(&root, v, w, u);
	}
};

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

	// vector<ourFace> allFaces;
	// for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
	// 	Triangle3D tri = dt.triangle(*faceItr);
	// 	allFaces.push_back(ourFace(getIndex(pts, tri[0]),
	// 							   getIndex(pts, tri[1]),
	// 							   getIndex(pts, tri[2])));
	// }
	// allFaces.shrink_to_fit();
	// std::sort(allFaces.begin(), allFaces.end());

	//vector<set<int>> adjList(pts.size());

	//~ for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
	//~ Segment3D seg = dt.segment(*edgeItr);
	//~ int u = getIndex(pts, seg[0]), v = getIndex(pts, seg[1]);
	//~ adjList[u].insert(v);
	//~ adjList[v].insert(u);
	//~ }

	/*for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
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
	}*/

	Trie t(dt,pts);

	/*vector<DT3::Vertex_handle> ptsHandle;
	ptsHandle.resize(pts.size());
	for (auto itr = dt.finite_vertices_begin(); itr != dt.finite_vertices_end(); itr++) {
		ptsHandle[getIndex(pts, itr->point())] = itr;
	}

	DT3::Cell_handle cell;
	int u, v, w;*/

	cout << std::boolalpha;

	start = std::chrono::high_resolution_clock::now();
	for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
		Triangle3D tri = dt.triangle(*faceItr);
		/*allFaces.push_back(ourFace(getIndex(pts, tri[0]),
								   getIndex(pts, tri[1]),
								   getIndex(pts, tri[2])));*/
	}
	finish = std::chrono::high_resolution_clock::now();
	cout << "Search Complete in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

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
