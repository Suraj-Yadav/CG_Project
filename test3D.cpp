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

//using namespace std;
using std::vector;
using std::set;
using std::pair;
using std::array;
using std::map;
using std::cerr;
using std::cout;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Vector_3 Vector3D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Triangle_3 Triangle3D;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;

const double inf = std::numeric_limits<double>::infinity();

vector<myTriple<size_t>> get_All_Facets(const DT3 &dt, const vector<Point3D> &points) {
	vector<myTriple<size_t>> allFacets;
	for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
		Triangle3D tri = dt.triangle(*faceItr);
		allFacets.push_back({getIndex(points, tri[0]),
							 getIndex(points, tri[1]),
							 getIndex(points, tri[2])});
	}
	return allFacets;
}

vector<myPair<size_t>> get_All_Edges(const DT3 &dt, const vector<Point3D> &points) {
	vector<myPair<size_t>> allEdges;
	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
		Segment3D seg = dt.segment(*edgeItr);
		allEdges.push_back({getIndex(points, seg[0]),
							getIndex(points, seg[1])});
	}
	return allEdges;
}

vector<myPair<size_t>> get_Mst_Edges_Kruskal(const vector<Point3D> &points, const DT3 &dt) {
	vector<CGAL::Union_find<size_t>::handle> handle;
	CGAL::Union_find<size_t> uf;
	handle.reserve(points.size());

	for (size_t i = 0; i < points.size(); i++)
		handle.push_back(uf.make_set(i));

	vector<std::tuple<double, size_t, size_t>> allEdges;
	for (auto edge : get_All_Edges(dt, points)) {
		allEdges.push_back(std::make_tuple(CGAL::squared_distance(points[edge[0]], points[edge[1]]), edge[0], edge[1]));
	}

	sort(allEdges.begin(),
		 allEdges.end(),
		 [](std::tuple<double, size_t, size_t> &a, std::tuple<double, size_t, size_t> &b) -> bool { return get<0>(a) < get<0>(b); });

	vector<myPair<size_t>> mst;

	for (auto edge : allEdges) {
		auto u = get<1>(edge), v = get<2>(edge);
		if (!uf.same_set(handle[u], handle[v])) {
			mst.push_back({u, v});
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

double getTriangleScore(const Point3D &a, const Point3D &b, const Point3D &c) {
	Vector3D AB(a, b), AC(a, c), BC(b, c);
	return std::sqrt(AB.squared_length()) + std::sqrt(AC.squared_length()) + std::sqrt(BC.squared_length());
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
	//if (!ans)
		//print("Problem Adding ", newPoint, "to", edge[0], edge[1], ". Sizeof edgedegree=", edgeDegree.size());
	return ans;
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

	double score = getTriangleScore(points[UVW[0]], points[UVW[1]], points[UVW[2]]);
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

vector<set<size_t>> getAdjList(size_t pointsCount, const vector<myPair<size_t>> &edges) {
	vector<set<size_t>> adjList(pointsCount);
	for (auto edge : edges) {
		adjList[edge[0]].insert(edge[1]);
		adjList[edge[1]].insert(edge[0]);
	}
	return adjList;
}

void process(const vector<Point3D> &points, const DT3 &dt, vector<myTriple<size_t>> &faces, vector<myPair<size_t>> &edges) {
	faces.clear();
	vector<myPair<size_t>> allEdges = get_All_Edges(dt, points);
	vector<set<size_t>> adjList = getAdjList(points.size(), allEdges);
	
	//edges = get_Mst_Edges_Kruskal(points, dt);	//Edges are already of MST

	std::queue<size_t> Q;
	vector<vector<size_t>> edgeDegree(allEdges.size());
	
	set<myTriple<size_t>> trianglesCovered;
	set<myPair<size_t>> edgesCovered;

	for (size_t i = 0; i < edges.size(); i++) sortThis(edges[i]);
	for (size_t i = 0; i < allEdges.size(); i++) sortThis(allEdges[i]);
	
	std::sort(allEdges.begin(), allEdges.end());
	std::sort(edges.begin(), edges.end());

	for (size_t i = 0; i < edges.size(); i++) edgesCovered.insert(edges[i]);
	
	for (size_t i = 0; i < edges.size(); i++) Q.push(getIndex(allEdges,edges[i]));
	
	while (!Q.empty()) {
		auto edgeIndex = Q.front();
		Q.pop();
		auto edge = allEdges[edgeIndex];
		switch (edgeDegree[edgeIndex].size()) {
			case 0: {
				set<size_t> commonPoints;
				std::set_intersection(
					adjList[edge[0]].begin(), adjList[edge[0]].end(), 
					adjList[edge[1]].begin(), adjList[edge[1]].end(), 
					std::inserter(commonPoints, commonPoints.end()));
				size_t nextPoint = SIZE_MAX;
				double score = inf;
				myPair<size_t> newEdge1, newEdge2;
				size_t newEdgeIndex1, newEdgeIndex2;
				myTriple<size_t> tri;
				for (auto point : commonPoints) {
					tri = {edge[0], edge[1], point};
					newEdge1 = {edge[0], point}, newEdge2 = {edge[1], point};
					sortThis(tri), sortThis(newEdge1), sortThis(newEdge2);

					newEdgeIndex1 = getIndex(allEdges, newEdge1);
					newEdgeIndex2 = getIndex(allEdges, newEdge2);

					double tempScore = getTriangleScore(points[tri[0]], points[tri[1]], points[tri[2]]);
					//print(tempScore, " ");
					if (trianglesCovered.find(tri) == trianglesCovered.end() &&
						validToAdd(points, edgeDegree[newEdgeIndex1], newEdge1, edge[1]) &&
						validToAdd(points, edgeDegree[newEdgeIndex2], newEdge2, edge[0]) &&
						tempScore > 0 &&
						tempScore < score) {
						nextPoint = point;
						score = tempScore;
					}
				}
				//println(":", score);
				if (nextPoint == SIZE_MAX)
					break;
				tri = {edge[0], edge[1], nextPoint};
				newEdge1 = {edge[0], nextPoint}, newEdge2 = {edge[1], nextPoint};
				sortThis(tri), sortThis(newEdge1), sortThis(newEdge2);
				newEdgeIndex1 = getIndex(allEdges, newEdge1);
				newEdgeIndex2 = getIndex(allEdges, newEdge2);
				
				edgesCovered.insert(newEdge1);
				edgesCovered.insert(newEdge2);
				
				edgeDegree[newEdgeIndex1].push_back(edge[1]);
				edgeDegree[newEdgeIndex2].push_back(edge[0]);
				edgeDegree[edgeIndex].push_back(nextPoint);
				
				trianglesCovered.insert(tri);
				faces.push_back(tri);
			}
			case 1: {
			}
			case 2:
			default:
				continue;
		}
	}
	//edges.clear();
	//for (auto edge : edgesCovered)
		//edges.push_back(edge);
}

int main(int argc, char *argv[]) {

	if (argc != 3) {
		cerr << "Invalid Argument\n";
		return 1;
	}

	std::ifstream inputFile(argv[1]);
	std::ofstream outputFile(argv[2]);
	DT3 dt;

	vector<Point3D> points;

	auto start = std::chrono::high_resolution_clock::now();
	if (!CGAL::read_xyz_points(inputFile, back_inserter(points))) { // output iterator over points
		cerr << "Error: cannot read file.";
		return 1;
	}

	sort(points.begin(), points.end());				  // vector may have repeated elements like 1 1 2 2 3 3 3 4 4 5 5 6 7
	auto last = unique(points.begin(), points.end()); // vector now holds {1 2 3 4 5 6 7 x x x x x x}, where 'x' is indeterminate
	points.erase(last, points.end());

	auto finish = std::chrono::high_resolution_clock::now();
	cout << points.size() << " points read in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	start = std::chrono::high_resolution_clock::now();

	dt.insert(points.begin(), points.end());

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

	vector<myPair<size_t>> edges;
	vector<myTriple<size_t>> faces;

	start = std::chrono::high_resolution_clock::now();
	edges = get_Mst_Edges_Kruskal(points, dt);
	finish = std::chrono::high_resolution_clock::now();
	cout << "Kruskal MST created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	/*double len = 0;
	for (auto edge : edges)
		len = std::max(len, edge.squared_length());

	edges = get_All_Edges(dt);*/

	start = std::chrono::high_resolution_clock::now();
	//faces = get_All_Facets(dt);
	process(points, dt, faces, edges);
	finish = std::chrono::high_resolution_clock::now();
	cout << "Faces created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	//start = chrono::high_resolution_clock::now();
	//auto mst2 = get_Mst_Edges_Prim(points, dt);
	//finish = chrono::high_resolution_clock::now();
	//cout << "Prim MST created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	//if (mst1.size() != mst2.size()) {
	//	cerr << "Size Different : Kruskal = " << mst1.size() << ", Prim = " << mst2.size() << "\n";
	//	return 1;
	//}

	start = std::chrono::high_resolution_clock::now();
	outputFile << points.size() << "\n";
	for (Point3D point : points) {
		outputFile << point << "\n";
	}

	outputFile << edges.size() << "\n";
	for (auto edge : edges) {
		outputFile << edge[0] << " " << edge[1] << "\n";
	}

	outputFile << faces.size() << "\n";
	for (auto triangle : faces) {
		outputFile << triangle[0] << " " << triangle[1] << " " << triangle[2] << "\n";
	}

	finish = std::chrono::high_resolution_clock::now();
	cout << "Output created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	return 0;
}