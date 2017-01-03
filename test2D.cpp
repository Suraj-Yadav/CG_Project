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

using namespace std;

#define val(x) cout << #x "=" << x << "\n";

const double inf = std::numeric_limits<double>::infinity();

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point2D;
typedef Kernel::Vector_2 Vector2D;
typedef Kernel::Segment_2 Segment2D;
typedef Kernel::Triangle_2 Triangle2D;
typedef CGAL::Delaunay_triangulation_2<Kernel> DT2;

template <typename T>
size_t getIndex(const vector<T> &v, const T &elem) {
	return lower_bound(v.begin(), v.end(), elem) - v.begin();
}

vector<Segment2D> get_All_Edges(const DT2 &dt) {
	vector<Segment2D> allEdges;
	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++)
		allEdges.push_back(dt.segment(*edgeItr));

	return allEdges;
}

vector<Segment2D> get_Mst_Edges_Kruskal(const vector<Point2D> &points, const DT2 &dt) {
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

	vector<Segment2D> mst;

	for (auto edge : allEdges) {
		auto u = get<1>(edge), v = get<2>(edge);
		if (!uf.same_set(handle[u], handle[v])) {
			mst.push_back(Segment2D(points[u], points[v]));
			uf.unify_sets(handle[u], handle[v]);
		}
	}

	return mst;
}

vector<Segment2D> get_Mst_Edges_Prim(const vector<Point2D> &points, const DT2 &dt) {
	vector<vector<pair<size_t, double>>> adjList;

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
			mst.push_back(Segment2D(points[u], points[parent[u]]));
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

double getScore(const Point2D &a, const Point2D &b, const Point2D &c) {
	Vector2D AB(a, b), AC(a, c), BC(b, c);
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

Kernel::FT cross(const Vector2D &a, const Vector2D &b) {
	return a.x() * b.y() - a.y() * b.x();
}

Kernel::FT turn(const Point2D &a, const Point2D &b, const Point2D &c) { //+ left , 0 collinear, - right
	return cross(b - a, c - b);
}

template <typename T>
pair<T, T> normalize_pair(const pair<T, T> &p) {
	if (p.first > p.second)
		return make_pair(p.second, p.first);
	else
		return make_pair(p.first, p.second);
}

bool validToAdd(const vector<Point2D> &points, const vector<size_t> &edgeDegree, const pair<size_t, size_t> edge, size_t newPoint) {
	if (edgeDegree.size() == 2)
		return false;
	if (edgeDegree.size() == 0)
		return true;
	auto t1 = turn(points[edge.first], points[edge.second], points[edgeDegree.back()]),
		 t2 = turn(points[edge.first], points[edge.second], points[newPoint]);
	return t1 * t2 < 0.0;
}

void process(const vector<Point2D> &points, vector<Triangle2D> &faces, vector<Segment2D> &edges) {
	faces.clear();
	vector<vector<size_t>> adjList(points.size());
	vector<bool> visited(points.size(), false);

	map<pair<size_t, size_t>, vector<size_t>> edgeDegree;
	//
	//vector<array<size_t, 3>> triangles;
	//set<array<size_t, 3>> trianglesCovered;

	vector<size_t> parent(points.size(), -1);

	for (auto edge : edges) {
		auto u = getIndex(points, edge.start()),
			 v = getIndex(points, edge.end());
		adjList[u].push_back(v);
		adjList[v].push_back(u);
	}

	edges.clear();
	queue<size_t> S;
	S.push(0);
	while (S.size() > 0) {
		size_t u = S.front();
		S.pop();
		if (visited[u])
			continue;
		visited[u] = true;
		for (size_t v : adjList[u])
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

	/*set<pair<double, size_t>> pq;
	for (size_t i = 0; i < points.size(); ++i) {
		std::sort(adjList[i].begin(),
				  adjList[i].end(),
				  [i, points](size_t a, size_t &b) -> bool {
					  return turn(points[a], points[i], points[b]) <= 0;
				  });
		for (size_t j = 1; j < adjList[i].size(); ++j) {
			size_t a = i, b = adjList[i][j], c = adjList[i][j - 1];
			if (append(trianglesCovered, {a, b, c})) {
				auto score = getScore(points[a], points[b], points[c]);
				if (score > 0) {
					pq.insert({score, triangles.size()});
					triangles.push_back({a, b, c});
				}
			}
		}
		//	for (size_t k = j + 1; k < adjList[i].size(); ++k) {
		//		size_t a = i, b = adjList[i][j], c = adjList[i][k];
		//		if (append(trianglesCovered, {a, b, c})) {
		//			pq.insert({getScore(points[a], points[b], points[c]), triangles.size()});
		//			triangles.push_back({a, b, c});
		//		}
		//	}
		//}
	}*/
	//while (pq.size() && faces.size() <= 3 * points.size()) {
	//	auto ptr = pq.end();
	//	ptr--;
	//	double score = ptr->first;
	//	size_t a = triangles[ptr->second][0],
	//		   b = triangles[ptr->second][1],
	//		   c = triangles[ptr->second][2];
	//	pq.erase(ptr);
	//	//if (score < 0.0) {
	//		//std::cout << "Skipped\n";
	//		//continue;
	//	//}
	//	faces.emplace_back(points[a], points[b], points[c]);

	//	for (size_t w : adjList[b])
	//		if (w != a && w != c && append(trianglesCovered, {b, c, w})) {
	//			auto score = getScore(points[b], points[c], points[w]);
	//			if (score > 0) {
	//				pq.insert({score, triangles.size()});
	//				triangles.push_back({b, c, w});
	//			}
	//		}
	//	for (size_t w : adjList[c])
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
}

int main(int argc, char *argv[]) {

	if (argc != 3) {
		cerr << "Invalid Argument\n";
		return 1;
	}

	ifstream inputFile(argv[1]);
	ofstream outputFile(argv[2]);
	DT2 dt;

	vector<Point2D> points;
	Point2D p;

	auto start = chrono::high_resolution_clock::now();
	while (inputFile >> p) {
		// ignore whatever comes after x and y
		inputFile.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
		points.push_back(p);
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
	if (dt.dimension() != 2) {
		cerr << "Error: cannot built a 2D triangulation.\n Current dimension = " << dt.dimension() << "\n";
		return 1;
	}
	finish = chrono::high_resolution_clock::now();
	cout << "Delaunay Triangulation created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	vector<Segment2D> edges;
	vector<Triangle2D> faces;

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
	for (Point2D point : points) {
		outputFile << point << " 0.0\n";
	}

	outputFile << edges.size() << "\n";
	for (Segment2D edge : edges) {
		outputFile << getIndex(points, edge[0]) << " " << getIndex(points, edge[1]) << "\n";
	}

	outputFile << faces.size() << "\n";
	for (Triangle2D triangle : faces) {
		outputFile << getIndex(points, triangle[0]) << " "
				   << getIndex(points, triangle[1]) << " "
				   << getIndex(points, triangle[2]) << "\n";
	}

	finish = chrono::high_resolution_clock::now();
	cout << "Output created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

	return 0;
}