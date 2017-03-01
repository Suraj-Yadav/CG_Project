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

// logFile for debugging
std::ofstream logFile;

vector<ourFace> get_All_Facets(const DT3 &dt, const vector<Point3D> &points) {
	vector<ourFace> allFacets;
	for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
		Triangle3D tri = dt.triangle(*faceItr);
		allFacets.push_back(ourFace(getIndex(points, tri[0]),
									getIndex(points, tri[1]),
									getIndex(points, tri[2])));
	}

	return allFacets;
}

vector<ourEdge> get_All_Edges(const DT3 &dt, const vector<Point3D> &points) {
	vector<ourEdge> allEdges;
	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
		Segment3D seg = dt.segment(*edgeItr);
		allEdges.push_back(ourEdge(getIndex(points, seg[0]),
								   getIndex(points, seg[1])));
	}

	return allEdges;
}

vector<ourEdge> get_Mst_Edges_Kruskal(const vector<Point3D> &points, const DT3 &dt) {
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
	vector<ourEdge> mst;
	for (auto edge : allEdges) {
		auto u = get<1>(edge), v = get<2>(edge);
		if (!uf.same_set(handle[u], handle[v])) {
			mst.push_back({u, v});
			uf.unify_sets(handle[u], handle[v]);
		}
	}

	return mst;
}

vector<ourEdge> get_Mst_Edges_Prim(const vector<Point3D> &points, const DT3 &dt) {
	vector<vector<pair<int, double>>> adjList;
	vector<ourEdge> mst;
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
	vector<int> parent(points.size(), -1);
	set<pair<double, int>> pq;
	weight[0] = 0;
	for (int i = 0; i < points.size(); i++)
		pq.insert({weight[i], i});
	//while (!pq.empty()) {
	for (int i = 0; i < points.size(); i++) {
		//int u = pq.begin()->second;
		int u = i;
		pq.erase(pq.begin());
		inPQ[u] = false;
		if (parent[u] != SIZE_MAX)
			mst.push_back({u, parent[u]});
		for (auto elem : adjList[u]) {
			int v = elem.first;
			double w = elem.second;
			if (inPQ[v] && weight[v] > w) {
				//pq.erase(pq.find({weight[v], v}));
				weight[v] = w;
				//pq.insert({w, v});
				parent[v] = u;
			}
		}
	}

	return mst;
}

double getTriangleScore(const Point3D &a, const Point3D &b, const Point3D &c) {
	Vector3D AB(a, b), AC(a, c), BC(b, c);
	//return std::sqrt(AB.squared_length()) + std::sqrt(AC.squared_length()) + std::sqrt(BC.squared_length());
	AB = AB / std::sqrt(AB.squared_length());
	AC = AC / std::sqrt(AC.squared_length());
	BC = BC / std::sqrt(BC.squared_length());
	double score = AB * AC;
	// score += std::acos(AC * BC);
	// score += std::acos(-AB * BC);
	score = std::min(AC * BC, score);
	score = std::min(-AB * BC, score);
	return score;
}

bool validToAdd(const vector<Point3D> &points, const vector<int> &edgeDegree, const ourEdge &edge, int newPoint, bool debug = false) {
	bool ans;
	double angle;
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
		Point3D A = points[edge[0]],
				B = points[edge[1]],
				C = points[edgeDegree.back()],
				D = points[newPoint];
		Vector3D AB(A, B), BD(B, D), BC(B, C);
		BD = BD - ((BD * AB) / AB.squared_length()) * AB;
		BC = BC - ((BC * AB) / AB.squared_length()) * AB;
		BD = BD / CGAL::sqrt(BD.squared_length());
		BC = BC / CGAL::sqrt(BC.squared_length());
		angle = BD * BC;
		angle = 180 * std::acos(angle) / M_PI;
		println( "Angle between the face (", edge[0], edge[1], edgeDegree.back(), ") and (", edge[0], edge[1], newPoint, ") is", angle);
	}

	if (!ans && debug) {
		println( "Problem Adding ", newPoint, "to", edge, ". edgedegree={");
		for (int p : edgeDegree)
			println( p);
		println( "}");
	}

	return ans;
}

vector<set<int>> getAdjList(int pointsCount, const vector<ourEdge> &edges) {
	vector<set<int>> adjList(pointsCount);
	for (auto edge : edges) {
		adjList[edge[0]].insert(edge[1]);
		adjList[edge[1]].insert(edge[0]);
	}

	return adjList;
}

int findCenterNode(vector<set<int>> adjList) {
	std::queue<int> Q;
	for (int i = 0; i < adjList.size(); i++)
		if (adjList[i].size() == 1)
			Q.push(i);
	int lastPushed = 0;
	while (!Q.empty()) {
		int u = Q.front();
		Q.pop();
		for (int v : adjList[u]) {
			adjList[v].erase(u);
			if (adjList[v].size() == 1) {
				Q.push(v);
				lastPushed = v;
			}
		}
	}

	return lastPushed;
}

int getShortestDistance(int a, int b, const vector<int> &distance, const vector<int> &parent) {
	int ans = 0;
	while (distance[a] > distance[b])
		ans++, a = parent[a];
	while (distance[a] < distance[b])
		ans++, b = parent[b];
	while (parent[a] != parent[b])
		ans += 2, a = parent[a], b = parent[b];
	return ans;
}

void processEdges(const vector<Point3D> &points, const DT3 &dt, vector<ourEdge> &edges, set<int> leftVerts, const set<ourEdge> &coveredEdges) {
	set<ourEdge> usedEdges;
	vector<ourEdge> allEdges;
	for (auto e : get_All_Edges(dt, points)) {
		if (leftVerts.find(e[0]) != leftVerts.end() &&
			leftVerts.find(e[1]) != leftVerts.end()) // &&
			//coveredEdges.find(e) == coveredEdges.end())
			allEdges.push_back(e);
	}

	for (ourEdge &e : allEdges)
		println( e);
	println( "###########################################################");
	vector<set<int>> adjList = getAdjList(points.size(), allEdges);
	for (int u : leftVerts) {
		int v = -1;
		double dist = inf;
		cout << "For " << u << " considering:";
		for (int tempV : adjList[u]) {
			if (dist > CGAL::squared_distance(points[u], points[tempV])) {
				v = tempV;
				dist = CGAL::squared_distance(points[u], points[tempV]);
				cout << v << " ";
			}
		}
		cout << "\n";

		if (v != -1 && usedEdges.find(ourEdge(u, v)) == usedEdges.end()) {
			usedEdges.insert(ourEdge(u, v));
			println( ourEdge(u, v));
		}
	}

	edges.clear();
	for (ourEdge e : usedEdges)
		edges.push_back(e);
}

int testEdges(Point3D a, Point3D b, Point3D c, double length) {
	int ans = 0;
	if (CGAL::squared_distance(a, b) <= length) ans++;
	if (CGAL::squared_distance(a, c) <= length) ans++;
	if (CGAL::squared_distance(c, b) <= length) ans++;
	return ans;
}

set<int> process(const vector<Point3D> &points, const DT3 &dt, vector<ourFace> &faces, vector<ourEdge> &edges) {
	faces.clear();
	vector<ourEdge> allEdges = get_All_Edges(dt, points);
	vector<set<int>> adjList = getAdjList(points.size(), edges); //make adj list from mst.
	vector<ourFace> allFaces = get_All_Facets(dt, points);
	// edges = get_Mst_Edges_Kruskal(points, dt);	//Edges are already of MST
	vector<vector<int>> edgeDegree(allEdges.size());
	set<ourFace> trianglesCovered;
	set<ourEdge> edgesCovered;
	double maxEdge = 0;
	for (ourEdge &e : edges)
		maxEdge = std::max(maxEdge, CGAL::squared_distance(points[e[0]], points[e[1]]));

	set<std::tuple<double, ourEdge, int>, std::greater<std::tuple<double, ourEdge, int>>> pq;

	std::sort(allEdges.begin(), allEdges.end());
	std::sort(allFaces.begin(), allFaces.end());
	std::sort(edges.begin(), edges.end());
	//for (ourEdge edge : allEdges)
	//println( edge);
	//for (int i = 0; i < allEdges.size(); i++)
	//println( i, ":", allEdges[i]);
	//println();
	//for (int i = 0; i < allFaces.size(); i++)
	//println( i, ":", allFaces[i]);
	//return;

	// cout<<"MST Edges\n";
	// 	for (auto &elem : edges) {
	// 		cout<<elem<<"\n";
	// 	}

	for (int i = 0; i < edges.size(); i++) edgesCovered.insert(edges[i]);
	for (auto &elem : edgeDegree) elem.clear();
	int edgeCondition = 2;
	for (auto edge : edges) {
		set<int> commonPoints;
		commonPoints.insert(adjList[edge[0]].begin(), adjList[edge[0]].end());
		commonPoints.insert(adjList[edge[1]].begin(), adjList[edge[1]].end());
		commonPoints.erase(edge[0]);
		commonPoints.erase(edge[1]);
		set<int> t, f;
		cout << edge << "\n";
		for (auto point : commonPoints) {
			// cout << point << "\n";
			double tempScore = getTriangleScore(points[edge[0]], points[edge[1]], points[point]);
			if (std::binary_search(allFaces.begin(), allFaces.end(), ourFace(edge[0], edge[1], point)) &&
				testEdges(points[edge[0]], points[edge[1]], points[point], maxEdge) >= edgeCondition) {
				pq.insert(std::make_tuple(tempScore, edge, point));
				// cout << "Initially added " << edge << " and " << point << "\n";
				t.insert(point);
			}
			else {
				f.insert(point);
			}
		}
		cout << "T:";
		for (auto &elem : t) {
			cout << elem << " ";
		}
		cout << "\nF:";
		for (auto &elem : f) {
			cout << elem << " ";
		}
		cout << "\n";
	}

	// for (auto &elem : pq) {
	// 	cout << "[" << get<0>(elem) << " " << get<1>(elem) << " " << get<2>(elem) << "],\n";
	// }
	// cout<<"\n";
	set<int> leftVerts;
	while (edgeCondition >= 2) {
		while (!pq.empty()) {
			/*for (auto &elem : pq) {
				cout << "[" << get<0>(elem) << " " << get<1>(elem) << " " << get<2>(elem) << "],\n";
			}
			cout << "\n";
			for (int i = 0; i < points.size(); ++i) {
				cout << i << ": [";
				for (auto elem : adjList[i]) {
					cout << elem << " ";
				}
				cout << "]\n";
			}
			for (int i = 0; i < edgeDegree.size(); ++i) {
				if (edgeDegree[i].size() >= 1) {
					cout << allEdges[i] << ": [";
					for (auto &p : edgeDegree[i])
						cout << p << " ";
					cout << "]\n";
				}
			}*/
			auto element = *pq.begin();
			pq.erase(pq.begin());
			ourEdge edge = get<1>(element);
			int point = get<2>(element);
			cout << get<0>(element) << " ";
			cout << edge << " " << point << "\n";
			int edgeIndex = getIndex(allEdges, edge);
			ourFace tri(edge[0], edge[1], point);
			println( "Considering", tri);
			if (edgeDegree[edgeIndex].size() == 2) {
				println( ": Rejected");
				println( val(edgeDegree[edgeIndex].size() == 2));
				continue;
			}

			ourEdge newEdge1(edge[0], point), newEdge2(edge[1], point);
			int newEdgeIndex1 = getIndex(allEdges, newEdge1),
				newEdgeIndex2 = getIndex(allEdges, newEdge2);
			if (trianglesCovered.find(tri) == trianglesCovered.end() &&
				validToAdd(points, edgeDegree[newEdgeIndex1], newEdge1, edge[1]) &&
				validToAdd(points, edgeDegree[newEdgeIndex2], newEdge2, edge[0]) &&
				testEdges(points[edge[0]], points[edge[1]], points[point], maxEdge) >= edgeCondition &&
				(edgeDegree[edgeIndex].size() == 0 || (edgeDegree[edgeIndex].size() == 1 && validToAdd(points, edgeDegree[edgeIndex], edge, point)))) {
				println( ": Accepted");
				cout << "Added\n";
				edgesCovered.insert(newEdge1);
				edgesCovered.insert(newEdge2);
				edgeDegree[newEdgeIndex1].push_back(edge[1]);
				edgeDegree[newEdgeIndex2].push_back(edge[0]);
				edgeDegree[edgeIndex].push_back(point);
				println( "Added", edge[1], "to", newEdge1, "which has index", newEdgeIndex1);
				println( "Added", edge[0], "to", newEdge2, "which has index", newEdgeIndex2);
				println( "Added", point, "to", edge, "which has index", edgeIndex);
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
					if (u == edge[0] && v == edge[1])
						continue;
					set<int> commonPoints;
					commonPoints.insert(adjList[u].begin(), adjList[u].end());
					commonPoints.insert(adjList[v].begin(), adjList[v].end());
					commonPoints.erase(u);
					commonPoints.erase(v);
					set<int> t, f;
					std::vector<std::pair<ourEdge, int>> slkjf;
					cout << ourEdge(u, v) << "\n";
					for (auto w : commonPoints) {
						ourFace tempTri(u, v, w);
						ourEdge tempEdge(u, v);
						// cout << w << "\n";
						double tempScore = getTriangleScore(points[u], points[v], points[w]);
						if (trianglesCovered.find(tempTri) == trianglesCovered.end() &&
							std::binary_search(allFaces.begin(), allFaces.end(), tempTri)) {
							pq.insert(std::make_tuple(tempScore, tempEdge, w));
							// cout << edge << " added " << tempEdge << " and " << w << "\n";
							slkjf.push_back(std::make_pair(tempEdge, w));
							t.insert(w);
						}
						else {
							f.insert(w);
						}
					}
					cout << "T:";
					for (auto &elem : t) {
						cout << elem << " ";
					}
					cout << "\nF:";
					for (auto &elem : f) {
						cout << elem << " ";
					}
					cout << "\n";
					for (auto elem : slkjf) {
						cout << edge << " added " << elem.first << " and " << elem.second << "\n";
					}
				}
			}

			else {
				println( ": Rejected");
				println( val(edge));
				println( val(newEdge1));
				println( val(newEdge2));
				println( val(trianglesCovered.find(tri) == trianglesCovered.end()));
				println( val(validToAdd(points, edgeDegree[newEdgeIndex1], newEdge1, edge[1], true)));
				println( val(validToAdd(points, edgeDegree[newEdgeIndex2], newEdge2, edge[0], true)));
				println( val(edgeDegree[edgeIndex].size() == 0));
				println( val(edgeDegree[edgeIndex].size() == 1 && validToAdd(points, edgeDegree[edgeIndex], edge, point, true)));
			}
		}

		// for (int i = 0; i < allEdges.size(); i++) {
		// 	if (edgeDegree[i].size() == 1) {
		// 		ourEdge edge = allEdges[i];
		// 		set<int> commonPoints;
		// 		commonPoints.insert(adjList[edge[0]].begin(), adjList[edge[0]].end());
		// 		commonPoints.insert(adjList[edge[1]].begin(), adjList[edge[1]].end());
		// 		commonPoints.erase(edge[0]);
		// 		commonPoints.erase(edge[1]);
		// 		for (auto point : commonPoints) {
		// 			double tempScore = getTriangleScore(points[edge[0]], points[edge[1]], points[point]);
		// 			if (std::binary_search(allFaces.begin(), allFaces.end(), ourFace(edge[0], edge[1], point)) &&
		// 				testEdges(points[edge[0]], points[edge[1]], points[point], maxEdge) >= edgeCondition)
		// 				pq.insert(std::make_tuple(tempScore, edge, point));
		// 		}
		// 	}
		// }

		edgeCondition--;

		leftVerts.clear();
		for (int i = 0; i < allEdges.size(); i++) {
			if (edgeDegree[i].size() == 1) {
				leftVerts.insert(allEdges[i][0]), leftVerts.insert(allEdges[i][1]);
			}
		}

		println( "Vertices with ", edgeCondition + 1);
		for (auto p : leftVerts)
			println( p);
		// println();
		vector<ourEdge> tempEdges;
		processEdges(points, dt, tempEdges, leftVerts, edgesCovered);
		for (auto edge : tempEdges) {
			set<int> commonPoints;
			commonPoints.insert(adjList[edge[0]].begin(), adjList[edge[0]].end());
			commonPoints.insert(adjList[edge[1]].begin(), adjList[edge[1]].end());
			commonPoints.erase(edge[0]);
			commonPoints.erase(edge[1]);
			for (auto point : commonPoints) {
				double tempScore = getTriangleScore(points[edge[0]], points[edge[1]], points[point]);
				if (std::binary_search(allFaces.begin(), allFaces.end(), ourFace(edge[0], edge[1], point)) &&
					testEdges(points[edge[0]], points[edge[1]], points[point], maxEdge) >= edgeCondition) {
					pq.insert(std::make_tuple(tempScore, edge, point));
					//cout << "Finally added " << edge << " and " << point << "\n";
				}
			}
		}
	}

	return leftVerts;
}

void writeOffFile(std::ofstream &outputFile, const vector<Point3D> &points, const vector<ourEdge> &edges, const vector<ourFace> &faces) {
	outputFile << "OFF\n";
	outputFile << points.size() << " " << faces.size() << "\n";
	for (const Point3D &p : points)
		outputFile << p << "\n";
	for (const ourFace &f : faces)
		outputFile << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
}

void writeForBlend(std::ofstream &outputFile, const vector<Point3D> &points, const vector<ourEdge> &edges, const vector<ourFace> &faces) {
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
}

int main(int argc, char *argv[]) {
	string logFilePath(argv[2]);
	// println(logFilePath.substr(0, logFilePath.find_last_of('\\')));

	if (argc != 3) {
		cerr << "Invalid Argument\n";
		return 1;
	}

	println("Input File =", argv[1]);
	println("Output File =", argv[2]);

	std::ifstream inputFile(argv[1]);
	std::ofstream outputFile(argv[2]);
	logFile.open("log.txt");
	logFile << std::boolalpha;

	DT3 dt;

	vector<Point3D> points;
	auto start = std::chrono::high_resolution_clock::now();
	if (!CGAL::read_xyz_points(inputFile, back_inserter(points))) { // output iterator over points
		cerr << "Error: cannot read file.";
		return 1;
	}
	sortAndRemoveDuplicate(points);

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
	// cout << points.size() << " " << std::chrono::duration<double>(finish - start).count() << " secs\n";
	vector<ourEdge> edges;
	vector<ourFace> faces;
	start = std::chrono::high_resolution_clock::now();
	//edges = get_Mst_Edges_Prim(points, dt);
	edges = get_Mst_Edges_Kruskal(points, dt);
	finish = std::chrono::high_resolution_clock::now();
	cout << "Kruskal MST created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";
	/*double len = 0;
	for (auto edge : edges)
		len = std::max(len, edge.squared_length());
	edges = get_All_Edges(dt);*/
	start = std::chrono::high_resolution_clock::now();
	set<int> p;
	//int lastSize = 0;
	// while (1) {
	p = process(points, dt, faces, edges);
	//println(p.size());
	// 	if (p.size() == lastSize || p.size() == 0)
	// 		break;
	// 	lastSize = p.size();
	// processEdges(points, dt, edges, p);
	// }

	finish = std::chrono::high_resolution_clock::now();
	cout << points.size() << " " << std::chrono::duration<double>(finish - start).count() << "\n";
	//start = chrono::high_resolution_clock::now();
	//auto mst2 = get_Mst_Edges_Prim(points, dt);
	//finish = chrono::high_resolution_clock::now();
	//cout << "Prim MST created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";
	//if (mst1.size() != mst2.size()) {
	//	cerr << "Size Different : Kruskal = " << mst1.size() << ", Prim = " << mst2.size() << "\n";
	//	return 1;
	//}

	start = std::chrono::high_resolution_clock::now();
	edges.clear();
	writeForBlend(outputFile, points, edges, faces);
	finish = std::chrono::high_resolution_clock::now();
	// cout << "Output created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";
	/*println( "___________________________MST edge lengths____________________");
	for (auto e : edges)
		println( CGAL::sqrt(CGAL::squared_distance(points[e[0]], points[e[1]])));
	println( "_______________________________________________________________");
	println( "_______________________Triangle Edge lengths___________________");
	for (auto f : faces) {
		println( CGAL::sqrt(CGAL::squared_distance(points[f[0]], points[f[1]])));
		println( CGAL::sqrt(CGAL::squared_distance(points[f[1]], points[f[2]])));
		println( CGAL::sqrt(CGAL::squared_distance(points[f[0]], points[f[2]])));
	}

	println( "_______________________________________________________________");
	println( "________________________Triangle Perimeter_____________________");
	for (auto f : faces) {
		println( CGAL::sqrt(CGAL::squared_distance(points[f[0]], points[f[1]])) + CGAL::sqrt(CGAL::squared_distance(points[f[1]], points[f[2]])) + CGAL::sqrt(CGAL::squared_distance(points[f[0]], points[f[2]])));
	}

	println( "_______________________________________________________________");*/
	return 0;
}
