#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Union_find.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <chrono>

using namespace std;

#define val(x) cout << #x "=" << x << "\n";

const double inf = std::numeric_limits<double>::infinity();

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Triangle_3 Triangle3D;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;

template <typename T>
size_t getIndex(const vector<T> &v, const T &elem) {
	return lower_bound(v.begin(), v.end(), elem) - v.begin();
}

vector<Triangle3D> get_All_Facets(const DT3 &dt) {
	vector<Triangle3D> allFacets;
	for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++)
		allFacets.push_back(dt.triangle(*faceItr));

	return allFacets;
}

vector<Segment3D> get_All_Edges(const DT3 &dt) {
	vector<Segment3D> allEdges;
	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++)
		allEdges.push_back(dt.segment(*edgeItr));

	return allEdges;
}

vector<Segment3D> get_Mst_Edges_Kruskal(const vector<Point3D> &points, const DT3 &dt) {
	vector<CGAL::Union_find<size_t>::handle> handle;
	CGAL::Union_find<size_t> uf;
	handle.reserve(points.size());

	for (size_t i = 0; i < points.size(); i++)
		handle.push_back(uf.make_set(i));

	vector<tuple<double, size_t, size_t>> allEdges;
	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
		auto edge = dt.segment(*edgeItr);
		auto u = getIndex(points, edge.start()),
			 v = getIndex(points, edge.end());
		allEdges.push_back(make_tuple(edge.squared_length(), u, v));
	}

	sort(allEdges.begin(),
		 allEdges.end(),
		 [](tuple<double, size_t, size_t> &a, tuple<double, size_t, size_t> &b) -> bool { return get<0>(a) < get<0>(b); });

	vector<Segment3D> mst;

	for (auto edge : allEdges) {
		auto u = get<1>(edge), v = get<2>(edge);
		if (!uf.same_set(handle[u], handle[v])) {
			mst.push_back(Segment3D(points[u], points[v]));
			uf.unify_sets(handle[u], handle[v]);
		}
	}

	return mst;
}

vector<Segment3D> get_Mst_Edges_Prim(const vector<Point3D> &points, const DT3 &dt) {
	vector<vector<pair<size_t, double>>> adjList;

	vector<Segment3D> mst;

	adjList.resize(points.size());
	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
		auto edge = dt.segment(*edgeItr);
		auto u = getIndex(points, edge.start()),
			 v = getIndex(points, edge.end());
		adjList[u].push_back({v, edge.squared_length()});
		adjList[v].push_back({u, edge.squared_length()});
	}

	vector<double> weight(points.size(), inf);
	vector<bool> inPQ(points.size(), true);
	vector<size_t> parent(points.size(), SIZE_MAX);
	set<pair<double, size_t>> pq;

	weight[0] = 0;
	for (size_t i = 0; i < points.size(); i++)
		pq.insert({weight[i], i});

	while (!pq.empty()) {
		size_t u = pq.begin()->second;
		pq.erase(pq.begin());
		inPQ[u] = false;
		if (parent[u] != SIZE_MAX)
			mst.push_back(Segment3D(points[u], points[parent[u]]));
		for (auto elem : adjList[u]) {
			size_t v = elem.first;
			double w = elem.second;
			if (inPQ[v] && weight[v] > w) {
				pq.erase(pq.find({weight[v], v}));
				weight[v] = w;
				pq.insert({w, v});
				parent[v] = u;
			}
		}
	}

	return mst;
}

int main(int argc, char *argv[]) {

	if (argc != 3) {
		cerr << "Invalid Argument\n";
		return 1;
	}

	ifstream inputFile(argv[1]);
	ofstream outputFile(argv[2]);
	DT3 dt;

	vector<Point3D> points;

	auto start = chrono::high_resolution_clock::now();
	if (!CGAL::read_xyz_points(inputFile, back_inserter(points))) { // output iterator over points
		cerr << "Error: cannot read file.";
		return 1;
	}

	sort(points.begin(), points.end());				  // vector may have repeated elements like 1 1 2 2 3 3 3 4 4 5 5 6 7
	auto last = unique(points.begin(), points.end()); // vector now holds {1 2 3 4 5 6 7 x x x x x x}, where 'x' is indeterminate
	points.erase(last, points.end());

	auto finish = chrono::high_resolution_clock::now();
	cout << points.size() << " points read in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	start = chrono::high_resolution_clock::now();

	dt.insert(points.begin(), points.end());

	if (!dt.is_valid(true)) {
		cerr << "Error: fail to build a Delaunay triangulation.\n";
		return 1;
	}
	if (dt.dimension() != 3) {
		cerr << "Error: cannot built a 3D triangulation.\n Current dimension = " << dt.dimension() << "\n";
		return 1;
	}
	finish = chrono::high_resolution_clock::now();
	cout << "Delaunay Triangulation created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	vector<Segment3D> edges;
	vector<Triangle3D> faces;

	start = chrono::high_resolution_clock::now();
	edges = get_Mst_Edges_Kruskal(points, dt);
	finish = chrono::high_resolution_clock::now();
	cout << "Kruskal MST created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	/*start = chrono::high_resolution_clock::now();
	faces = get_All_Facets(dt);
	finish = chrono::high_resolution_clock::now();
	cout << "Faces created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";*/

	//start = chrono::high_resolution_clock::now();
	//auto mst2 = get_Mst_Edges_Prim(points, dt);
	//finish = chrono::high_resolution_clock::now();
	//cout << "Prim MST created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	//if (mst1.size() != mst2.size()) {
	//	cerr << "Size Different : Kruskal = " << mst1.size() << ", Prim = " << mst2.size() << "\n";
	//	return 1;
	//}

	start = chrono::high_resolution_clock::now();
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

	finish = chrono::high_resolution_clock::now();
	cout << "Output created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	return 0;
}