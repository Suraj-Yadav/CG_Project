#define _USE_MATH_DEFINES
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Union_find.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <chrono>
#include <queue>
#include <cmath>
#include <map>

#include "util.h"

using namespace std;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Vector_3 Vector3D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Triangle_3 Triangle3D;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;

const double inf = std::numeric_limits<double>::infinity();

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

double getScore(const Point3D &a, const Point3D &b, const Point3D &c) {
	Vector3D AB(a, b), AC(a, c), BC(b, c);
	AB = AB / std::sqrt(AB.squared_length());
	AC = AC / std::sqrt(AC.squared_length());
	BC = BC / std::sqrt(BC.squared_length());
	double score = (AB * AC);
	score = std::min(AC * BC, score);
	score = std::min(-AB * BC, score);
	return score;
}

bool append(set<array<size_t, 3>> &trianglesCovered, array<size_t, 3> A) {
	array<size_t, 3> arr;
	arr[0] = std::min(std::min(A[0], A[1]), A[2]);
	arr[2] = std::max(std::max(A[0], A[1]), A[2]);
	arr[1] = A[0] + A[1] + A[2] - arr[0] - arr[2];
	auto ans = trianglesCovered.insert(arr);
	return ans.second;
}

bool validToAdd(const vector<Point3D> &points, const vector<size_t> &edgeDegree, const myPair<size_t> &edge, size_t newPoint) {
	bool ans;
	if (edgeDegree.size() == 2)
		ans = false;
	else if (edgeDegree.size() == 0)
		ans = true;
	else {
		Vector3D norm1 = CGAL::cross_product(points[edge[1]] - points[edge[0]], points[newPoint] - points[edge[1]]),
				 norm2 = CGAL::cross_product(points[edge[0]] - points[edge[1]], points[edgeDegree.back()] - points[edge[0]]);
		if (CGAL::angle(norm1, norm2) >= CGAL::RIGHT)
			ans = true;
		else
			ans = false;
	}
	if (!ans)
		print("Problem Adding ", newPoint, "to", edge[0], edge[1], ". Sizeof edgedegree=", edgeDegree.size());
	return ans;
}

template <typename T, typename T2>
size_t addToMap(map<T, size_t> &mapping, vector<T2> &vec, const T &object) {
	if (mapping.find(object) == mapping.end()) {
		mapping.insert({object, vec.size()});
		vec.emplace_back();
	}
	return mapping[object];
}

bool process2(const vector<Point3D> &points,
			  vector<set<size_t>> &adjList,
			  set<myTriple<size_t>> &trianglesCovered,
			  map<myPair<size_t>, size_t> &edgeDegreeMap,
			  vector<vector<size_t>> &edgeDegree,
			  myTriple<size_t> UVW) {
	print("Checking ", UVW[0], UVW[1], UVW[2]);
	
	auto u = UVW[0], v = UVW[1], w = UVW[2];

	if (trianglesCovered.find(UVW) != trianglesCovered.end()) {
		print(u, v, w, "is already present");
		return false;
	}

	double score = getScore(points[UVW[0]], points[UVW[1]], points[UVW[2]]);
	if (score < 0) {
		print("Score = ", score);
		return false;
	}
	myPair<size_t> UV = {UVW[0], UVW[1]}, UW = {UVW[0], UVW[2]}, WV = {UVW[1], UVW[2]};
	sortThis(UV), sortThis(UW), sortThis(WV);

	auto UVIndex = addToMap(edgeDegreeMap, edgeDegree, UV),
		 UWIndex = addToMap(edgeDegreeMap, edgeDegree, UW),
		 WVIndex = addToMap(edgeDegreeMap, edgeDegree, WV);


	if (!validToAdd(points, edgeDegree[UVIndex], UV, w)) return false;
	if (!validToAdd(points, edgeDegree[UWIndex], UW, v)) return false;
	if (!validToAdd(points, edgeDegree[WVIndex], WV, u)) return false;

	trianglesCovered.insert(UVW);

	edgeDegree[UVIndex].push_back(w);
	edgeDegree[UWIndex].push_back(v);
	edgeDegree[WVIndex].push_back(u);
	adjList[u].insert(v), adjList[u].insert(w);
	adjList[v].insert(w), adjList[v].insert(u);
	adjList[w].insert(u), adjList[w].insert(v);
	return true;
}

void process(const vector<Point3D> &points, vector<Triangle3D> &faces, vector<Segment3D> &edges) {
	faces.clear();
	vector<set<size_t>> adjList(points.size());
	vector<bool> visited(points.size(), false);

	map<myPair<size_t>, size_t> edgeDegreeMap;
	vector<vector<size_t>> edgeDegree;
	//
	// vector<array<size_t, 3>> triangles;
	set<myTriple<size_t>> trianglesCovered;

	vector<size_t> parent(points.size(), SIZE_MAX);

	for (auto edge : edges) {
		auto u = getIndex(points, edge.start()),
			v = getIndex(points, edge.end());
		adjList[u].insert(v);
		adjList[v].insert(u);
	}
	edges.clear();
	queue<size_t> S;
	S.push(0);
	while (S.size() > 0) {
		size_t u = S.front(), w = parent[u];
		S.pop();
		if (visited[u])
			continue;
		visited[u] = true;
		for (size_t v : adjList[u]) {
			if (visited[v])
				continue;
			S.push(v);
			parent[v] = u;
			if (u == 0)
				continue;
			myTriple<size_t> UVW = {u, v, w};
			sortThis(UVW);
			process2(points, adjList, trianglesCovered, edgeDegreeMap, edgeDegree, UVW);
			// if (score > 0)
			// faces.emplace_back(points[u], points[v], points[parent[u]]);

			edges.emplace_back(points[u], points[v]);
		}
	}

	bool flag = true;
	while (flag && edgeDegreeMap.size() < 6000) {
		flag = false;
		print("edgeDegreeMap.size=", edgeDegreeMap.size());
		for (auto elem : edgeDegreeMap) {
			if (edgeDegree[elem.second].size() >= 2)
				continue;
			size_t u = elem.first[0], v = elem.first[1];
			set<size_t> Union;
			Union.insert(adjList[u].begin(), adjList[u].end());
			Union.insert(adjList[v].begin(), adjList[v].end());

			for (size_t w : Union) {
				if (w == u || w == v)
					continue;
				myTriple<size_t> UVW = {u, v, w};
				sortThis(UVW);
				if (process2(points, adjList, trianglesCovered, edgeDegreeMap, edgeDegree, UVW))
					flag = true;

			}
		}
	}

	//set<pair<double, size_t>> pq;
	//for (size_t i = 0; i < points.size(); ++i) {
	//	for (size_t j = 0; j < adjList[i].size(); ++j) {
	//		for (size_t k = j + 1; k < adjList[i].size(); ++k) {
	//			size_t a = i, b = adjList[i][j], c = adjList[i][k];
	//			if (append(trianglesCovered, {a, b, c})) {
	//				pq.insert({getScore(points[a], points[b], points[c]), triangles.size()});
	//				triangles.push_back({a, b, c});
	//			}
	//		}
	//	}
	//}
	//while (pq.size() && faces.size() <= 3 * points.size()) {
	//	auto ptr = pq.end();
	//	ptr--;
	//	double score = ptr->first;
	//	size_t a = triangles[ptr->second][0],
	//		   b = triangles[ptr->second][1],
	//		   c = triangles[ptr->second][2];
	//	pq.erase(ptr);
	//	if (score > M_PI_2) {
	//		std::cout << "Skipped\n";
	//		continue;
	//	}
	//	faces.emplace_back(points[a], points[b], points[c]);

	//	for (size_t w : adjList[b])
	//		if (w != a && w != c && append(trianglesCovered, {b, c, w})) {
	//			pq.insert({getScore(points[b], points[c], points[w]), triangles.size()});
	//			triangles.push_back({b, c, w});
	//		}
	//	for (size_t w : adjList[c])
	//		if (w != a && w != b && append(trianglesCovered, {c, b, w})) {
	//			pq.insert({getScore(points[c], points[b], points[w]), triangles.size()});
	//			triangles.push_back({c, b, w});
	//		}
	//	if (pq.size() % 1000 == 0)
	//		val(pq.size());
	//}
	for (auto triple : trianglesCovered) {
		//std::cout << triple[0] << " " << triple[1] << " " << triple[2] << "\n";
		faces.emplace_back(points[triple[0]], points[triple[1]], points[triple[2]]);
	}
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

	/*double len = 0;
	for (auto edge : edges)
		len = std::max(len, edge.squared_length());

	edges = get_All_Edges(dt);*/

	start = chrono::high_resolution_clock::now();
	//faces = get_All_Facets(dt);
	process(points, faces, edges);
	finish = chrono::high_resolution_clock::now();
	cout << "Faces created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

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