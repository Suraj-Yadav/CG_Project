#define _USE_MATH_DEFINES
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
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
typedef Kernel::Point_2 Point2D;
typedef Kernel::Vector_2 Vector2D;
typedef Kernel::Segment_2 Segment2D;
typedef Kernel::Triangle_2 Triangle2D;
typedef CGAL::Delaunay_triangulation_2<Kernel> DT2;
const double inf = std::numeric_limits<double>::infinity();

vector<myPair<int>> get_All_Edges(const DT2 &dt, const vector<Point2D> &points) {
	vector<myPair<int>> allEdges;
	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
		Segment2D seg = dt.segment(*edgeItr);
		allEdges.push_back({getIndex(points, seg[0]),
							getIndex(points, seg[1])});
	}
	return allEdges;
}

vector<myPair<int>> get_Mst_Edges_Kruskal(const vector<Point2D> &points, const DT2 &dt) {
	vector<CGAL::Union_find<int>::handle> handle;
	CGAL::Union_find<int> uf;
	handle.reserve(points.size());

	for (int i = 0; i < points.size(); i++)
		handle.push_back(uf.make_set(i));

	vector<std::tuple<double, int, int>> allEdges;
	for (auto edge : get_All_Edges(dt, points)) {
		allEdges.push_back(std::make_tuple(CGAL::squared_distance(points[edge[0]], points[edge[1]]), edge[0], edge[1]));
	}

	sort(allEdges.begin(),
		 allEdges.end(),
		 [](std::tuple<double, int, int> &a, std::tuple<double, int, int> &b) -> bool { return get<0>(a) < get<0>(b); });

	vector<myPair<int>> mst;

	for (auto edge : allEdges) {
		auto u = get<1>(edge), v = get<2>(edge);
		if (!uf.same_set(handle[u], handle[v])) {
			mst.push_back({u, v});
			uf.unify_sets(handle[u], handle[v]);
		}
	}
	return mst;
}

vector<Segment2D> get_Mst_Edges_Prim(const vector<Point2D> &points, const DT2 &dt) {
	vector<vector<pair<int, double>>> adjList;

	vector<Segment2D> mst;

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
	vector<int> parent(points.size(), SIZE_MAX);
	set<pair<double, int>> pq;

	weight[0] = 0;
	for (int i = 0; i < points.size(); i++)
		pq.insert({weight[i], i});

	while (!pq.empty()) {
		int u = pq.begin()->second;
		pq.erase(pq.begin());
		inPQ[u] = false;
		if (parent[u] != SIZE_MAX)
			mst.push_back(Segment2D(points[u], points[parent[u]]));
		for (auto elem : adjList[u]) {
			int v = elem.first;
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

double getTriangleScore(const Point2D &a, const Point2D &b, const Point2D &c) {
	Vector2D AB(a, b), AC(a, c), BC(b, c);
	AB = AB / std::sqrt(AB.squared_length());
	AC = AC / std::sqrt(AC.squared_length());
	BC = BC / std::sqrt(BC.squared_length());
	double score = AB * AC;
	score = std::min(AC * BC, score);
	score = std::min(-AB * BC, score);
	return score;
}

bool append(set<array<int, 3>> &trianglesCovered, array<int, 3> A) {
	array<int, 3> arr;
	arr[0] = std::min(std::min(A[0], A[1]), A[2]);
	arr[2] = std::max(std::max(A[0], A[1]), A[2]);
	arr[1] = A[0] + A[1] + A[2] - arr[0] - arr[2];
	auto ans = trianglesCovered.insert(arr);
	return ans.second;
}

bool validToAdd(const vector<Point2D> &points, const vector<int> &edgeDegree, const myPair<int> &edge, int newPoint) {
	if (edgeDegree.size() == 2)
		return false;
	if (edgeDegree.size() == 0)
		return true;
	auto t1 = CGAL::orientation(points[edge[0]], points[edge[1]], points[edgeDegree.back()]),
		 t2 = CGAL::orientation(points[edge[0]], points[edge[1]], points[newPoint]);
	return t1 != t2;
}

/*void process(const vector<Point2D> &points, vector<Triangle2D> &faces, vector<Segment2D> &edges) {
	faces.clear();
	vector<vector<int>> adjList(points.size());
	vector<bool> visited(points.size(), false);

	map<pair<int, int>, vector<int>> edgeDegree;
	//
	//vector<array<int, 3>> triangles;
	//set<array<int, 3>> trianglesCovered;

	vector<int> parent(points.size(), -1);

	for (auto edge : edges) {
		auto u = getIndex(points, edge.start()),
			 v = getIndex(points, edge.end());
		adjList[u].push_back(v);
		adjList[v].push_back(u);
	}

	edges.clear();
	queue<int> S;
	S.push(0);
	while (S.size() > 0) {
		int u = S.front();
		S.pop();
		if (visited[u])
			continue;
		visited[u] = true;
		for (int v : adjList[u])
			if (!visited[v]) {
				S.push(v);
				parent[v] = u;
				if (parent[u] != -1) {
					double score = getScore(points[u], points[v], points[parent[u]]);
					if (validToAdd(points, edgeDegree[normalize_pair(make_pair(u, v))], normalize_pair(make_pair(u, v)), parent[u]) &&
						validToAdd(points, edgeDegree[normalize_pair(make_pair(u, parent[u]))], normalize_pair(make_pair(u, parent[u])), v) &&
						validToAdd(points, edgeDegree[normalize_pair(make_pair(parent[u], v))], normalize_pair(make_pair(parent[u], v)), u)) {
						//if (score > 0)
						faces.emplace_back(points[u], points[v], points[parent[u]]);
						edgeDegree[normalize_pair(make_pair(u, v))].push_back(parent[u]);
						edgeDegree[normalize_pair(make_pair(u, parent[u]))].push_back(v);
						edgeDegree[normalize_pair(make_pair(v, parent[u]))].push_back(u);
					}
				}
				edges.emplace_back(points[u], points[v]);
			}
	}

	/*set<pair<double, int>> pq;
	for (int i = 0; i < points.size(); ++i) {
		std::sort(adjList[i].begin(),
				  adjList[i].end(),
				  [i, points](int a, int &b) -> bool {
					  return turn(points[a], points[i], points[b]) <= 0;
				  });
		for (int j = 1; j < adjList[i].size(); ++j) {
			int a = i, b = adjList[i][j], c = adjList[i][j - 1];
			if (append(trianglesCovered, {a, b, c})) {
				auto score = getScore(points[a], points[b], points[c]);
				if (score > 0) {
					pq.insert({score, triangles.size()});
					triangles.push_back({a, b, c});
				}
			}
		}
		//	for (int k = j + 1; k < adjList[i].size(); ++k) {
		//		int a = i, b = adjList[i][j], c = adjList[i][k];
		//		if (append(trianglesCovered, {a, b, c})) {
		//			pq.insert({getScore(points[a], points[b], points[c]), triangles.size()});
		//			triangles.push_back({a, b, c});
		//		}
		//	}
		//}
	}
	//while (pq.size() && faces.size() <= 3 * points.size()) {
	//	auto ptr = pq.end();
	//	ptr--;
	//	double score = ptr->first;
	//	int a = triangles[ptr->second][0],
	//		   b = triangles[ptr->second][1],
	//		   c = triangles[ptr->second][2];
	//	pq.erase(ptr);
	//	//if (score < 0.0) {
	//		//std::cout << "Skipped\n";
	//		//continue;
	//	//}
	//	faces.emplace_back(points[a], points[b], points[c]);

	//	for (int w : adjList[b])
	//		if (w != a && w != c && append(trianglesCovered, {b, c, w})) {
	//			auto score = getScore(points[b], points[c], points[w]);
	//			if (score > 0) {
	//				pq.insert({score, triangles.size()});
	//				triangles.push_back({b, c, w});
	//			}
	//		}
	//	for (int w : adjList[c])
	//		if (w != a && w != b && append(trianglesCovered, {c, b, w})) {
	//			auto score = getScore(points[c], points[b], points[w]);
	//			if (score > 0) {
	//				pq.insert({score, triangles.size()});
	//				triangles.push_back({c, b, w});
	//			}
	//		}
	//	if (pq.size() % 1000 == 0)
	//		val(pq.size());
	//}
}*/
vector<set<int>> getAdjList(int pointsCount, const vector<myPair<int>> &edges) {
	vector<set<int>> adjList(pointsCount);
	for (auto edge : edges) {
		adjList[edge[0]].insert(edge[1]);
		adjList[edge[1]].insert(edge[0]);
	}
	return adjList;
}

void process(const vector<Point2D> &points, const DT2 &dt, vector<myTriple<int>> &faces, vector<myPair<int>> &edges) {
	faces.clear();
	vector<myPair<int>> allEdges = get_All_Edges(dt, points);
	vector<set<int>> adjList = getAdjList(points.size(), edges);

	//edges = get_Mst_Edges_Kruskal(points, dt);	//Edges are already of MST

	vector<vector<int>> edgeDegree(allEdges.size());

	set<myTriple<int>> trianglesCovered;
	set<myPair<int>> edgesCovered;

	set<std::tuple<double, myPair<int>, int>, std::greater<std::tuple<double, myPair<int>, int>>> pq;

	for (int i = 0; i < edges.size(); i++) sortThis(edges[i]);
	for (int i = 0; i < allEdges.size(); i++) sortThis(allEdges[i]);

	std::sort(allEdges.begin(), allEdges.end());
	std::sort(edges.begin(), edges.end());

	for (int i = 0; i < edges.size(); i++) edgesCovered.insert(edges[i]);

	for (auto edge : edges) {
		set<int> commonPoints;
		std::set_union(
			adjList[edge[0]].begin(), adjList[edge[0]].end(),
			adjList[edge[1]].begin(), adjList[edge[1]].end(),
			std::inserter(commonPoints, commonPoints.end()));
		commonPoints.erase(edge[0]);
		commonPoints.erase(edge[1]);
		for (auto point : commonPoints) {
			double tempScore = getTriangleScore(points[edge[0]], points[edge[1]], points[point]);
			pq.insert({tempScore, edge, point});
		}
	}

	while (!pq.empty()) {
		auto element = *pq.begin();
		pq.erase(pq.begin());
		myPair<int> edge = get<1>(element);
		int point = get<2>(element);
		int edgeIndex = getIndex(allEdges, edge);
		if (edgeDegree[edgeIndex].size() == 2)
			continue;

		myTriple<int> tri = {edge[0], edge[1], point};
		myPair<int> newEdge1 = {edge[0], point}, newEdge2 = {edge[1], point};
		sortThis(tri), sortThis(newEdge1), sortThis(newEdge2);
		int newEdgeIndex1 = getIndex(allEdges, newEdge1),
			newEdgeIndex2 = getIndex(allEdges, newEdge2);

		if (trianglesCovered.find(tri) == trianglesCovered.end() &&
			validToAdd(points, edgeDegree[newEdgeIndex1], newEdge1, edge[1]) &&
			validToAdd(points, edgeDegree[newEdgeIndex2], newEdge2, edge[0]) &&
			(edgeDegree[edgeIndex].size() == 0 || (edgeDegree[edgeIndex].size() == 1 && validToAdd(points, edgeDegree[edgeIndex], edge, point)))) {

			edgesCovered.insert(newEdge1);
			edgesCovered.insert(newEdge2);

			edgeDegree[newEdgeIndex1].push_back(edge[1]);
			edgeDegree[newEdgeIndex2].push_back(edge[0]);
			edgeDegree[edgeIndex].push_back(point);

			trianglesCovered.insert(tri);
			faces.push_back(tri);

			adjList[edge[0]].insert(edge[1]);
			adjList[edge[0]].insert(point);
			adjList[edge[1]].insert(edge[0]);
			adjList[edge[1]].insert(point);
			adjList[point].insert(edge[0]);
			adjList[point].insert(edge[1]);

			for (int i = 0; i < 3; i++) {
				int u = tri[i], v = tri[(i + 1) % 3];

				if (u > v) std::swap(u, v);

				set<int> commonPoints;
				std::set_union(
					adjList[u].begin(), adjList[u].end(),
					adjList[v].begin(), adjList[v].end(),
					std::inserter(commonPoints, commonPoints.end()));
				commonPoints.erase(u);
				commonPoints.erase(v);

				for (auto w : commonPoints) {
					myTriple<int> tempTri = {u, v, w};
					double tempScore = getTriangleScore(points[u], points[v], points[w]);
					if (trianglesCovered.find(tempTri) == trianglesCovered.end())
						pq.insert({tempScore,{u, v}, w});
				}
			}
		}
	}
	/*switch (edgeDegree[edgeIndex].size()) {
	case 0: {
	}
	case 1: {
	set<int> commonPoints;
	std::set_union(
	adjList[edge[0]].begin(), adjList[edge[0]].end(),
	adjList[edge[1]].begin(), adjList[edge[1]].end(),
	std::inserter(commonPoints, commonPoints.end()));
	commonPoints.erase(edge[0]);
	commonPoints.erase(edge[1]);
	int nextPoint = SIZE_MAX;
	double score = -2;
	myPair<int> newEdge1, newEdge2;
	int newEdgeIndex1, newEdgeIndex2;
	myTriple<int> tri;
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
	(edgeDegree[edgeIndex].size() == 0 || (edgeDegree[edgeIndex].size() == 1 && validToAdd(points, edgeDegree[edgeIndex], edge, point))) &&
	tempScore > -0.1 &&
	tempScore > score) {
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

	adjList[edge[0]].insert(edge[1]);
	adjList[edge[0]].insert(nextPoint);
	adjList[edge[1]].insert(edge[0]);
	adjList[edge[1]].insert(nextPoint);
	adjList[nextPoint].insert(edge[0]);
	adjList[nextPoint].insert(edge[1]);

	if (edgeDegree[edgeIndex].size() < 2)
	Q.push(edgeIndex);
	if (edgeDegree[newEdgeIndex1].size() < 2)
	Q.push(newEdgeIndex1);
	if (edgeDegree[newEdgeIndex2].size() < 2)
	Q.push(newEdgeIndex2);

	trianglesCovered.insert(tri);
	faces.push_back(tri);
	}
	case 2:
	default:
	continue;
	}
	}
	edges.clear();
	for (auto edge : edgesCovered)
	edges.push_back(edge);
	*/
}

/*void process(const vector<Point2D> &points, const DT2 &dt, vector<myTriple<int>> &faces, vector<myPair<int>> &edges) {
	faces.clear();
	vector<myPair<int>> allEdges = get_All_Edges(dt, points);
	vector<set<int>> adjList = getAdjList(points.size(), edges);

	//edges = get_Mst_Edges_Kruskal(points, dt);	//Edges are already of MST

	std::queue<int> Q;
	vector<vector<int>> edgeDegree(allEdges.size());

	set<myTriple<int>> trianglesCovered;
	set<myPair<int>> edgesCovered;

	for (int i = 0; i < edges.size(); i++) sortThis(edges[i]);
	for (int i = 0; i < allEdges.size(); i++) sortThis(allEdges[i]);

	std::sort(allEdges.begin(), allEdges.end());
	std::sort(edges.begin(), edges.end());

	for (int i = 0; i < edges.size(); i++) edgesCovered.insert(edges[i]);

	for (int i = 0; i < edges.size(); i++) Q.push(getIndex(allEdges, edges[i]));

	while (!Q.empty()) {
		auto edgeIndex = Q.front();
		Q.pop();
		auto edge = allEdges[edgeIndex];
		switch (edgeDegree[edgeIndex].size()) {
			case 0: {
			}
			case 1: {
				set<int> commonPoints;
				std::set_union(
					adjList[edge[0]].begin(), adjList[edge[0]].end(),
					adjList[edge[1]].begin(), adjList[edge[1]].end(),
					std::inserter(commonPoints, commonPoints.end()));
				commonPoints.erase(edge[0]);
				commonPoints.erase(edge[1]);
				int nextPoint = SIZE_MAX;
				double score = -2;
				myPair<int> newEdge1, newEdge2;
				int newEdgeIndex1, newEdgeIndex2;
				myTriple<int> tri;
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
						(edgeDegree[edgeIndex].size() == 0 || (edgeDegree[edgeIndex].size() == 1 && validToAdd(points, edgeDegree[edgeIndex], edge, point))) &&
						tempScore > -0.1 &&
						tempScore > score) {
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

				adjList[edge[0]].insert(edge[1]);
				adjList[edge[0]].insert(nextPoint);
				adjList[edge[1]].insert(edge[0]);
				adjList[edge[1]].insert(nextPoint);
				adjList[nextPoint].insert(edge[0]);
				adjList[nextPoint].insert(edge[1]);

				if (edgeDegree[edgeIndex].size() < 2)
					Q.push(edgeIndex);
				if (edgeDegree[newEdgeIndex1].size() < 2)
					Q.push(newEdgeIndex1);
				if (edgeDegree[newEdgeIndex2].size() < 2)
					Q.push(newEdgeIndex2);

				trianglesCovered.insert(tri);
				faces.push_back(tri);
			}
			case 2:
			default:
				continue;
		}
	}
	//edges.clear();
	//for (auto edge : edgesCovered)
	//edges.push_back(edge);
}*/

int main(int argc, char *argv[]) {

	if (argc != 3) {
		cerr << "Invalid Argument\n";
		return 1;
	}

	std::ifstream inputFile(argv[1]);
	std::ofstream outputFile(argv[2]);
	DT2 dt;

	vector<Point2D> points;
	Point2D p;

	auto start = std::chrono::high_resolution_clock::now();
	while (inputFile >> p) {
		// ignore whatever comes after x and y
		inputFile.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
		points.push_back(p);
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
	if (dt.dimension() != 2) {
		cerr << "Error: cannot built a 2D triangulation.\n Current dimension = " << dt.dimension() << "\n";
		return 1;
	}
	finish = std::chrono::high_resolution_clock::now();
	cout << "Delaunay Triangulation created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	vector<myPair<int>> edges;
	vector<myTriple<int>> faces;

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
	for (Point2D point : points) {
		outputFile << point << " 0.0\n";
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
