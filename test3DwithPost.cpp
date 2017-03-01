#define _USE_MATH_DEFINES
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Union_find.h>

#include "util.h"

// Alias for CGAL datatypes
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Plane_3 Plane3D;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Vector_3 Vector3D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Triangle_3 Triangle3D;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;
const double inf = std::numeric_limits<double>::infinity();

vector<ourEdge> getMstEdges(const vector<Point3D> &points, const DT3 &dt) {
	vector<CGAL::Union_find<int>::handle> handle;
	CGAL::Union_find<int> uf;
	handle.reserve(points.size());
	for (int i = 0; i < points.size(); i++)
		handle.push_back(uf.make_set(i));
	vector<std::tuple<double, int, int>> allEdges;
	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
		Segment3D seg = dt.segment(*edgeItr);
		allEdges.push_back(std::make_tuple(CGAL::squared_distance(seg[0], seg[1]),
										   getIndex(points, seg[0]),
										   getIndex(points, seg[1])));
	}

	allEdges.shrink_to_fit();
	sort(allEdges.begin(), allEdges.end());

	vector<ourEdge> mst;

	for (auto edge : allEdges) {
		auto u = get<1>(edge), v = get<2>(edge);
		if (!uf.same_set(handle[u], handle[v])) {
			mst.push_back({u, v});
			uf.unify_sets(handle[u], handle[v]);
		}
	}
	std::sort(mst.begin(), mst.end());
	return mst;
}

vector<set<int>> getAdjList(int pointsCount, const vector<ourEdge> &edges) {
	vector<set<int>> adjList(pointsCount);
	for (auto edge : edges) {
		adjList[edge[0]].insert(edge[1]);
		adjList[edge[1]].insert(edge[0]);
	}

	return adjList;
}

class SurfaceReconstruct {
	vector<Point3D> pts;
	DT3 dt;
	bool valid;
	vector<ourFace> dtFaces;
	set<ourEdge> modelEdges;
	set<ourFace> modelFaces;
	map<ourEdge, vector<int>> edgeDegree;
	double maxEdge;
	double maxMinRatio;
	vector<vector<set<ourEdge>::iterator>> vertAdjEdges;

	bool isFaceDelaunay(const ourFace &face) {
		return std::binary_search(dtFaces.cbegin(), dtFaces.cend(), face);
	}

	double getTriScore(int u, int v, int w) {
		Point3D a = pts[u], b = pts[v], c = pts[w];
		Vector3D AB(a, b), AC(a, c), BC(b, c);
		AB = AB / std::sqrt(AB.squared_length());
		AC = AC / std::sqrt(AC.squared_length());
		BC = BC / std::sqrt(BC.squared_length());
		double score = AB * AC;
		score = std::min(AC * BC, score);
		score = std::min(-AB * BC, score);
		return score;
	}

	void addFaceToModel(ourFace face) {
		ourEdge e1(face[0], face[1]),
			e2(face[0], face[2]),
			e3(face[1], face[2]);

		int p1 = face[2], p2 = face[1], p3 = face[0];

		auto res = modelFaces.insert(face);

		if (!res.second)
			return;

		modelEdges.insert(e1);
		modelEdges.insert(e2);
		modelEdges.insert(e3);

		edgeDegree[e1].push_back(p1);
		edgeDegree[e2].push_back(p2);
		edgeDegree[e3].push_back(p3);

		vertAdjEdges[p1].push_back(modelEdges.find(e1));
		vertAdjEdges[p2].push_back(modelEdges.find(e2));
		vertAdjEdges[p3].push_back(modelEdges.find(e3));
	}

	int testEdges(Point3D a, Point3D b, Point3D c, double length) {
		int ans = 0;
		double maxLength = std::max(std::max(CGAL::squared_distance(a, b), CGAL::squared_distance(b, c)), CGAL::squared_distance(a, c)),
			   minLength = std::min(std::min(CGAL::squared_distance(a, b), CGAL::squared_distance(b, c)), CGAL::squared_distance(a, c));
		if (CGAL::squared_distance(a, b) <= length)
			ans++;
		if (CGAL::squared_distance(a, c) <= length)
			ans++;
		if (CGAL::squared_distance(c, b) <= length)
			ans++;
		if (maxLength < (4 * minLength * maxMinRatio)) {
			maxMinRatio = std::max(maxMinRatio, maxLength / minLength);
		}
		else {
			ans = -10;
		}
		return ans;
	}

	vector<pair<int, double>> getFacesFromEdge(const ourEdge &edge, const vector<set<int>> &currAdjList, int edgeCondition) {
		vector<pair<int, double>> elems;
		set<int> t, f;
		// cout << edge << "\n";
		for (auto &point : currAdjList[edge[0]]) {
			if (point != edge[1] &&
				isFaceDelaunay(ourFace(point, edge[0], edge[1])) &&
				modelFaces.find(ourFace(point, edge[0], edge[1])) == modelFaces.end()) { // &&
				// testEdges(pts[edge[0]], pts[edge[1]], pts[point], maxEdge) >= edgeCondition) {
				elems.push_back(std::make_pair(point, getTriScore(point, edge[0], edge[1])));
				t.insert(point);
			}
			else if (point != edge[1]) {
				f.insert(point);
			}
		}
		for (auto &point : currAdjList[edge[1]]) {
			if (point != edge[0] &&
				isFaceDelaunay(ourFace(point, edge[0], edge[1])) &&
				modelFaces.find(ourFace(point, edge[0], edge[1])) == modelFaces.end()) { // &&
				// testEdges(pts[edge[0]], pts[edge[1]], pts[point], maxEdge) >= edgeCondition) {
				elems.push_back(std::make_pair(point, getTriScore(point, edge[0], edge[1])));
				t.insert(point);
			}
			else if (point != edge[0]) {
				f.insert(point);
			}
		}
		// cout << "T:";
		// for (auto &elem : t) {
		// 	cout << elem << " ";
		// }
		// cout << "\nF:";
		// for (auto &elem : f) {
		// 	cout << elem << " ";
		// }
		// cout << "\n";
		sortAndRemoveDuplicate(elems);
		return elems;
	}

	bool isEdgeOverlap(const ourEdge &e) {
		//cout << "testing " << e << std::endl;
		ourEdge edge(0, 1);
		for (int i = 0; i < 2; ++i) {
			for (auto &elem : vertAdjEdges[e[i]]) {
				edge = *elem;
				//cout << "with " << edge << " " << std::endl;
				if (edge[0] == e[!i])
					continue;
				if (edge[1] == e[!i])
					continue;
				Point3D a = pts[edge[0]], b = pts[edge[1]], c, o = pts[e[i]];
				Plane3D plane(a, b, o);
				c = plane.projection(pts[e[!i]]);
				auto d1 = CGAL::cross_product(a - o, c - o) * CGAL::cross_product(a - o, b - o),
					 d2 = CGAL::cross_product(b - o, c - o) * CGAL::cross_product(b - o, a - o);
				if (d1 > 0 && d2 > 0) {
					// cout << e << " rejected by " << ourFace(edge[0], edge[1], e[i]) << "\n";
					return true;
				}
			}
		}
		//cout << "pass" << std::endl;
		return false;
	}

	bool validToAdd(ourEdge edge, int newPoint) {
		auto &edgeDeg = edgeDegree[edge];
		if (edgeDeg.size() == 2)
			return false;
		if (edgeDeg.size() == 0)
			return true;

		Vector3D norm1 = CGAL::normal(pts[edge[0]], pts[edge[1]], pts[newPoint]),
				 norm2 = CGAL::normal(pts[edge[1]], pts[edge[0]], pts[edgeDeg[0]]);
		norm1 = norm1 / std::sqrt(norm1.squared_length());
		norm2 = norm2 / std::sqrt(norm2.squared_length());
		return norm1 * norm2 >= 0.5;
	}

	set<ourEdge> fillHoles(set<int> &leftVerts, set<ourEdge> &allEdges) {
		set<ourEdge> newEdges;
		map<int, set<pair<double, int>>> tempAdjList;

		if (leftVerts.size() == 3) {
			int p[3], i = 0;
			for (auto &elem : leftVerts) {
				p[i++] = elem;
			}
			addFaceToModel({p[0], p[1], p[2]});
			return newEdges;
		}

		if (leftVerts.size() == 4) {
			int p[4], i = 0;
			for (auto &elem : leftVerts) {
				p[i++] = elem;
			}

			if (modelEdges.find({p[0], p[1]}) == modelEdges.end()) {
				addFaceToModel({p[0], p[1], p[2]});
				addFaceToModel({p[0], p[1], p[3]});
				return newEdges;
			}
			else if (modelEdges.find({p[0], p[2]}) == modelEdges.end()) {
				addFaceToModel({p[0], p[1], p[2]});
				addFaceToModel({p[0], p[2], p[3]});
				return newEdges;
			}
			else {
				addFaceToModel({p[0], p[1], p[3]});
				addFaceToModel({p[0], p[2], p[3]});
				return newEdges;
			}
		}

		set<std::pair<double, ourEdge>> edges;
		map<int, bool> covered;
		for (auto &elem : leftVerts) {
			covered[elem] = false;
		}

		for (auto itr = allEdges.begin(); itr != allEdges.end();) {
			int u = (*itr)[0], v = (*itr)[1];
			if (leftVerts.find(u) != leftVerts.end() &&
				leftVerts.find(v) != leftVerts.end() &&
				modelEdges.find(ourEdge(u, v)) == modelEdges.end()) {
				double sqDist = CGAL::squared_distance(this->pts[u], this->pts[v]);
				edges.insert({sqDist, *itr});
				tempAdjList[u].insert({sqDist, v});
				tempAdjList[v].insert({sqDist, u});
				++itr;
			}
			else
				++itr;
		}
		while (!edges.empty()) {
			ourEdge edge = edges.begin()->second;
			edges.erase(edges.begin());
			if (!covered[edge[0]] || !covered[edge[1]]) {
				newEdges.insert(edge);
				covered[edge[0]] = true;
				covered[edge[1]] = true;
			}
		}
		/*
		while (!leftVerts.empty()) {
			int u = *leftVerts.begin();
			int v = -1;
			double dist = inf;
			leftVerts.erase(leftVerts.begin());
			// cout << "U  : " << u << " : ";
			int i = tempAdjList[u].size() >> 1;
			for (auto tempP : tempAdjList[u]) {
				int tempV = tempP.second;
				if (dist > tempP.first && !isEdgeOverlap(ourEdge(u, tempV))) {
					v = tempV;
					// cout << "changing v to  " << tempV << " ";
					dist = tempP.first;
					break;
				}
				if (not--i)
					break;
			}

			// for (int tempV : tempAdjList[u]) {
			// 	if (dist > CGAL::squared_distance(pts[u], pts[tempV]) &&
			// 		!isEdgeOverlap(ourEdge(u, tempV))) {
			// 		v = tempV;
			// 		// cout << "changing v to  " << tempV << " ";
			// 		dist = CGAL::squared_distance(pts[u], pts[tempV]);
			// 	}
			// }
			// cout << "\n";
			if (v != -1 && tempAdjList[u].size() > 1) {
				// cout << "V  : " << v << '\n';
				newEdges.insert(ourEdge(u, v));
				allEdges.erase(ourEdge(u, v));
			}
		}*/
		return newEdges;
	}

	void reconstruct(vector<ourEdge> &initialEdges) {
		vector<set<int>> currAdjList = getAdjList(pts.size(), initialEdges);
		set<ourEdge> origEdges(initialEdges.begin(), initialEdges.end());
		maxEdge = 0;
		maxMinRatio = 1;
		for (auto &e : initialEdges)
			maxEdge = std::max(maxEdge, CGAL::squared_distance(pts[e[0]], pts[e[1]]));
		set<std::tuple<double, ourEdge, int>, std::greater<std::tuple<double, ourEdge, int>>> pq;
		modelEdges.insert(initialEdges.begin(), initialEdges.end());
		int edgeCondition = 2;
		// for (auto &elem : pq) {
		// 	cout << "[" << get<0>(elem) << " " << get<1>(elem) << " " << get<2>(elem) << "],\n";
		// }
		// cout << "\n";
		set<ourEdge> nextEdges;
		while (edgeCondition >= -1) {
			println("Iteration #=", 2 - edgeCondition + 1);
			nextEdges.insert(origEdges.begin(), origEdges.end());
			for (auto &edge : nextEdges) {
				bool inserted = false;
				for (auto &elem : getFacesFromEdge(edge, currAdjList, edgeCondition)) {
					pq.insert(std::make_tuple(elem.second, edge, elem.first));
					inserted = true;
					// cout << "Initially added " << edge << " and " << elem.first << std::endl;
				}
				if (inserted)
					origEdges.erase(edge);
			}

			while (!pq.empty()) {
				/*for (auto &elem : pq) {
					cout << "[" << get<0>(elem) << " " << get<1>(elem) << " " << get<2>(elem) << "]," << std::endl;
				}
				cout << std::endl;
				for (int i = 0; i < pts.size(); ++i) {
					cout << i << ": [";
					for (auto elem : currAdjList[i]) {
						cout << elem << " ";
					}
					cout << "]" << std::endl;
				}
				for (auto &elem : edgeDegree) {
					if (elem.second.size() >= 1) {
						cout << elem.first << ": [";
						for (auto &p : elem.second)
							cout << p << " ";
						cout << "]" << std::endl;
					}
				}*/
				ourEdge edge = get<1>(*pq.begin());
				int point = get<2>(*pq.begin());
				// cout << get<0>(*pq.begin()) << " " << edge << " " << point << std::endl;
				pq.erase(pq.begin());
				ourFace face(edge[0], edge[1], point);
				if (modelFaces.find(face) != modelFaces.end())
					continue;
				if (edgeDegree[edge].size() == 2)
					continue;
				if (isEdgeOverlap(edge)) {
					// cout << "FACE REJECTED\n";
					continue;
				}
				ourEdge newEdge1(edge[0], point), newEdge2(edge[1], point);
				// cerr << get<0>(*pq.begin()) << " " << edge << " " << point << "::";
				// cerr << (modelFaces.find(face) == modelFaces.end()) << " "
				// 	 << isFaceDelaunay(face) << " "
				// 	 << validToAdd(edge, point) << " "
				// 	 << validToAdd(newEdge1, edge[1]) << " "
				// 	 << validToAdd(newEdge2, edge[0]) << " "
				// 	 << (testEdges(pts[face[0]], pts[face[1]], pts[face[2]], maxEdge) >= edgeCondition) << "\n";
				if (modelFaces.find(face) == modelFaces.end() &&
					isFaceDelaunay(face) &&
					validToAdd(edge, point) &&
					validToAdd(newEdge1, edge[1]) &&
					validToAdd(newEdge2, edge[0]) &&
					testEdges(pts[face[0]], pts[face[1]], pts[face[2]], maxEdge) >= edgeCondition) {
					// cout << "Added" << std::endl;

					addFaceToModel(face);

					for (int i = 0; i < 3; i++) {
						currAdjList[face[i]].insert(face[(i + 1) % 3]);
						currAdjList[face[i]].insert(face[(i + 2) % 3]);
					}

					for (int i = 0; i < 3; i++) {
						int u = face[i], v = face[(i + 1) % 3];
						if (u > v)
							std::swap(u, v);
						if (u == edge[0] && v == edge[1])
							continue;
						for (auto &elem : getFacesFromEdge(ourEdge(u, v), currAdjList, 0)) {
							pq.insert(std::make_tuple(elem.second, ourEdge(u, v), elem.first));
							// cout << edge << " added " << ourEdge(u, v) << " and " << elem.first << "\n";
						}
					}
				}
			}

			edgeCondition--;

			map<int, bool> leftVerts;
			map<int, set<int>> tempAdjList;

			for (auto &elem : modelEdges) {
				// print(elem, ":");
				// for (auto &point : edgeDegree[elem]) {
				// 	print(point);
				// }
				// println();
				if (edgeDegree[elem].size() == 1) {
					leftVerts.insert({elem[0], false});
					leftVerts.insert({elem[1], false});

					// leftVerts.insert(elem[0]);
					// leftVerts.insert(elem[1]);

					tempAdjList[elem[0]].insert(elem[1]);
					tempAdjList[elem[1]].insert(elem[0]);
				}
			}
			// println();
			cout << "leftVerts=" << leftVerts.size() << "\n";
			// for (auto &elem : leftVerts) {
			// 	print("(", elem.first, elem.second, "),");
			// }
			// println();
			// for (auto &elem : tempAdjList) {
			// 	print(elem.first, ":");
			// 	for (auto &v : elem.second) {
			// 		print(v);
			// 	}
			// 	println();
			// }

			if (leftVerts.size() == 0)
				break;

			set<ourEdge> allEdges;

			//allEdges

			for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
				Segment3D seg = dt.segment(*edgeItr);
				int u = getIndex(pts, seg[0]), v = getIndex(pts, seg[1]);
				if (leftVerts.find(u) != leftVerts.end() &&
					leftVerts.find(v) != leftVerts.end() &&
					modelEdges.find(ourEdge(u, v)) == modelEdges.end()) {
					allEdges.insert(ourEdge(u, v));
				}
			}

			nextEdges.clear();

			for (auto &elem : leftVerts) {
				if (elem.second)
					continue;
				int u = elem.first;
				set<int> hole;
				std::queue<int> Q;
				Q.push(u);
				while (!Q.empty()) {
					int u = Q.front();
					Q.pop();
					leftVerts[u] = true;
					hole.insert(u);
					// cout << u << " ";
					for (auto &v : tempAdjList[u]) {
						if (!leftVerts[v]) {
							Q.push(v);
						}
					}
				}
				cout << "Considering Hole of size: " << hole.size() << " ";
				for (auto &elem : hole) {
					cout << elem << " ";
				}
				cout << "\n";
				auto E = fillHoles(hole, allEdges);
				cout << "Hole Edges=" << E.size() << "\n";
				nextEdges.insert(E.begin(), E.end());
			}

			cout << "nextEdges=" << nextEdges.size() << "\n";

			// modelEdges.clear();
			// break;
		}
		modelEdges.swap(nextEdges);
		cout << "MaxMinRatio=" << maxMinRatio << "\n";
	}

  public:
	SurfaceReconstruct(const char *inputFilePath) {
		std::ifstream inputFile(inputFilePath);

		auto start = std::chrono::high_resolution_clock::now();
		if (!CGAL::read_xyz_points(inputFile, back_inserter(pts))) { // output iterator over points
			cerr << "Error: cannot read file.";
			return;
		}
		sortAndRemoveDuplicate(pts);
		pts.shrink_to_fit();

		auto finish = std::chrono::high_resolution_clock::now();
		cout << pts.size() << " points read in " << std::chrono::duration<double>(finish - start).count() << " secs" << std::endl;

		start = std::chrono::high_resolution_clock::now();
		dt.insert(pts.begin(), pts.end());
		if (!dt.is_valid(true)) {
			cerr << "Error: fail to build a Delaunay triangulation." << std::endl;
			valid = false;
			return;
		}
		if (dt.dimension() != 3) {
			cerr << "Error: cannot built a 3D triangulation.\n Current dimension = " << dt.dimension() << std::endl;
			valid = false;
			return;
		}
		finish = std::chrono::high_resolution_clock::now();
		cout << "Delaunay Triangulation created in " << std::chrono::duration<double>(finish - start).count() << " secs" << std::endl;

		valid = true;

		vertAdjEdges.resize(pts.size());

		for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
			Triangle3D tri = dt.triangle(*faceItr);
			dtFaces.push_back(ourFace(getIndex(pts, tri[0]),
									  getIndex(pts, tri[1]),
									  getIndex(pts, tri[2])));
		}

		std::sort(dtFaces.begin(), dtFaces.end());

		start = std::chrono::high_resolution_clock::now();
		vector<ourEdge> initialEdges = getMstEdges(pts, dt);
		finish = std::chrono::high_resolution_clock::now();
		cout << "MST created in " << std::chrono::duration<double>(finish - start).count() << " secs" << std::endl;

		start = std::chrono::high_resolution_clock::now();
		// modelEdges.insert(initialEdges.begin(), initialEdges.end());
		reconstruct(initialEdges);
		finish = std::chrono::high_resolution_clock::now();
		cout << pts.size() << " Reconstructed in " << std::chrono::duration<double>(finish - start).count() << std::endl;
	}

	bool isValid() { return valid; }

	void writeModel(const char *outputFilePath) {
		std::ofstream outputFile(outputFilePath);
		outputFile << pts.size() << std::endl;
		for (Point3D point : pts) {
			outputFile << point << std::endl;
		}

		outputFile << modelEdges.size() << std::endl;
		for (auto edge : modelEdges) {
			outputFile << edge[0] << " " << edge[1] << std::endl;
		}

		outputFile << modelFaces.size() << std::endl;
		for (auto triangle : modelFaces) {
			outputFile << triangle[0] << " " << triangle[1] << " " << triangle[2] << std::endl;
		}
	}
};

int main(int argc, char *argv[]) {

	if (argc != 3) {
		cerr << "Invalid Argument" << std::endl;
		return 1;
	}

	println("Input File =", argv[1]);
	println("Output File =", argv[2]);

	cout << std::boolalpha;

	SurfaceReconstruct surface(argv[1]);
	if (surface.isValid()) {
		surface.writeModel(argv[2]);
	}
	else
		return 1;

	return 0;
}
