#define _USE_MATH_DEFINES
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/number_utils.h>
#include <CGAL/Union_find.h>

#include "util.h"

#define sqDist(x, y) CGAL::squared_distance(x, y)

// Alias for CGAL datatypes
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Plane_3 Plane3D;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Vector_3 Vector3D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Triangle_3 Triangle3D;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;
typedef CGAL::Surface_mesh<Point3D> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

const Kernel::FT inf = std::numeric_limits<Kernel::FT>::infinity();

vector<set<int>> getAdjList(int pointsCount, const vector<ourEdge> &edges) {
	vector<set<int>> adjList(pointsCount);
	for (auto &edge : edges) {
		adjList[edge[0]].insert(edge[1]);
		adjList[edge[1]].insert(edge[0]);
	}
	return adjList;
}

float progress = 0;

class progress_bar {
	std::atomic<bool> finished;
	std::thread t;

  public:
	progress_bar()
		: finished(false), t(std::ref(*this)) {
	}					// initiate the bar by starting a new thread
	void operator()() { // function executed by the thread
		int barWidth = 70;
		while (!finished) {
			std::this_thread::sleep_for(std::chrono::milliseconds(500));
			cerr << "[";
			int pos = barWidth * progress;
			for (int i = 0; i < barWidth; ++i) {
				if (i < pos)
					cerr << "=";
				else if (i == pos)
					cerr << ">";
				else
					cerr << " ";
			}
			cerr << "] " << int(progress * 100.0) << " %\r";
			cerr.flush();
		}
	}
	void terminate() { // tell the thread/bar to stop
		finished = true;
		if (t.joinable())
			t.join();
		int barWidth = 70;
		cerr << "[";
		int pos = barWidth;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos)
				cerr << "=";
			else if (i == pos)
				cerr << ">";
			else
				cerr << " ";
		}
		cerr << "] " << int(100.0) << " %\r";
		cerr.flush();
	}
};

class SurfaceReconstruct {
	vector<Point3D> pts;
	vector<Kernel::FT> maxLength;
	DT3 dt;
	bool valid;
	set<ourFace> possibleFaces;
	set<ourEdge> possibleEdges;
	set<ourEdge> modelEdges;
	set<ourFace> modelFaces;
	map<ourEdge, vector<int>> edgeDegree;
	Kernel::FT maxEdge;
	// Kernel::FT maxMinRatio;
	vector<vector<set<ourEdge>::iterator>> vertAdjEdges;

	vector<ourEdge> getMstEdges() {
		vector<CGAL::Union_find<int>::handle> handle;
		CGAL::Union_find<int> uf;
		handle.reserve(pts.size());
		for (int i = 0; i < pts.size(); i++)
			handle.push_back(uf.make_set(i));

		vector<std::pair<Kernel::FT, ourEdge>> allEdges;


		for (auto &edge : possibleEdges) {
			allEdges.push_back({sqDist(pts[edge[0]], pts[edge[1]]), edge});
		}

		allEdges.shrink_to_fit();
		sort(allEdges.begin(), allEdges.end());

		println(__LINE__, val(allEdges.size()));

		vector<ourEdge> mst;

		for (auto &elem : allEdges) {
			auto u = elem.second[0], v = elem.second[1];
			if (!uf.same_set(handle[u], handle[v])) {
				mst.push_back(elem.second);
				uf.unify_sets(handle[u], handle[v]);
				possibleEdges.erase(elem.second);
			}
		}
		std::sort(mst.begin(), mst.end());
		println(__LINE__, val(mst.size()));
		return mst;
	}

	bool isFaceDelaunay(const ourFace &face) {
		return possibleFaces.find(face) != possibleFaces.end();
	}

	Kernel::FT getTriScore(int u, int v, int w) {
		Point3D a = pts[u], b = pts[v], c = pts[w];
		Vector3D AB(a, b), AC(a, c), BC(b, c);
		AB = AB / CGAL::sqrt(AB.squared_length());
		AC = AC / CGAL::sqrt(AC.squared_length());
		BC = BC / CGAL::sqrt(BC.squared_length());
		Kernel::FT score = AB * AC;
		score = std::min(AC * BC, score);
		score = std::min(-AB * BC, score);
		return score;
	}

	void addFaceToModel(ourFace face) {
		ourEdge e1(face[0], face[1]), e2(face[0], face[2]), e3(face[1], face[2]);

		int p1 = face[2], p2 = face[1], p3 = face[0];

		auto res = modelFaces.insert(face);

		println("Added", face);

		if (!res.second)
			return;

		possibleFaces.erase(face);

		modelEdges.insert(e1);
		modelEdges.insert(e2);
		modelEdges.insert(e3);

		maxLength[p1] = std::max(std::max(sqDist(pts[e2[0]], pts[e2[1]]),
										  sqDist(pts[e3[0]], pts[e3[0]])),
								 maxLength[p1]);
		maxLength[p2] = std::max(std::max(sqDist(pts[e1[0]], pts[e1[1]]),
										  sqDist(pts[e3[0]], pts[e3[0]])),
								 maxLength[p2]);
		maxLength[p3] = std::max(std::max(sqDist(pts[e2[0]], pts[e2[1]]),
										  sqDist(pts[e1[0]], pts[e1[0]])),
								 maxLength[p3]);

		possibleEdges.erase(e1);
		possibleEdges.erase(e2);
		possibleEdges.erase(e3);

		edgeDegree[e1].push_back(p1);
		edgeDegree[e2].push_back(p2);
		edgeDegree[e3].push_back(p3);

		vertAdjEdges[p1].push_back(modelEdges.find(e1));
		vertAdjEdges[p2].push_back(modelEdges.find(e2));
		vertAdjEdges[p3].push_back(modelEdges.find(e3));
	}

	void removeFace(ourEdge edge, int remove) {
		ourEdge edge1(ourEdge(edge[0], remove)), edge0(ourEdge(edge[1], remove));

		modelFaces.erase(ourFace(edge[0], edge[1], remove));

		deleteFromVector(vertAdjEdges[remove], modelEdges.find(edge));
		deleteFromVector(vertAdjEdges[edge[0]], modelEdges.find(edge0));
		deleteFromVector(vertAdjEdges[edge[1]], modelEdges.find(edge1));

		deleteFromVector(edgeDegree[edge], remove);
		deleteFromVector(edgeDegree[edge1], edge[1]);
		deleteFromVector(edgeDegree[edge0], edge[0]);

		if (edgeDegree[edge1].size() == 0) {
			modelEdges.erase(edge1);
		}
		if (edgeDegree[edge0].size() == 0) {
			modelEdges.erase(edge0);
		}
	}

	int testEdges(Point3D a, Point3D b, Point3D c, Kernel::FT length) {
		int ans = 0;
		Kernel::FT maxiLength =
					   std::max(std::max(sqDist(a, b), sqDist(b, c)), sqDist(a, c)),
				   miniLength =
					   std::min(std::min(sqDist(a, b), sqDist(b, c)), sqDist(a, c));
		if (sqDist(a, b) <= length)
			ans++;
		if (sqDist(a, c) <= length)
			ans++;
		if (sqDist(c, b) <= length)
			ans++;
		// if (maxiLength < (4 * miniLength * maxMinRatio)) {
		// 	maxMinRatio = std::max(maxMinRatio, maxiLength / miniLength);
		// }
		// else {
		// 	ans = -10;
		// }
		return ans;
	}

	vector<pair<int, Kernel::FT>> getFacesFromEdge(const ourEdge &edge, const vector<set<int>> &currAdjList, int edgeCondition) {
		vector<pair<int, Kernel::FT>> elems;
		// set<int> t, f;
		// cout << edge << "\n";
		for (int i = 0; i < 2; ++i) {
			for (auto &point : currAdjList[edge[i]]) {
				if (point != edge[!i] &&
					isFaceDelaunay(ourFace(point, edge[i], edge[!i])) &&
					modelFaces.find(ourFace(point, edge[i], edge[!i])) == modelFaces.end()) { // &&
					// testEdges(pts[edge[0]], pts[edge[1]], pts[point], maxEdge) >=
					// edgeCondition) {
					elems.push_back(
						std::make_pair(point, getTriScore(point, edge[i], edge[!i])));
					// t.insert(point);
				}
				// else if (point != edge[1]) {
				// 	f.insert(point);
				// }
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

	bool faceEdgeOverlap(int A, int B, int O, int D) {
		if (A == D || B == D)
			return false;
		Point3D a = pts[A], b = pts[B], c, o = pts[O];
		Plane3D plane(a, b, o);
		c = plane.projection(pts[D]);

		Vector3D OA = a - o, OC = c - o, OB = b - o;
		OA = OA / CGAL::sqrt(OA.squared_length());
		OB = OB / CGAL::sqrt(OA.squared_length());
		OC = OC / CGAL::sqrt(OA.squared_length());

		auto d1 = CGAL::cross_product(OA, OC) * CGAL::cross_product(OA, OB),
			 d2 = CGAL::cross_product(OB, OC) * CGAL::cross_product(OB, OA);
		if (d1 > 1e-4 && d2 > 1e-4) {
			// cout << e << " rejected by " << ourFace(edge[0], edge[1], e[i]) << "\n";
			println(val(d1), val(d2));
			return true;
		}
		return false;
	}

	bool isLocked(ourFace face) {
		for (int i = 0; i < 3; ++i) {
			if (edgeDegree[ourEdge(face[i], face[(i + 1) % 3])].size() < 2)
				return false;
		}
		return true;
	}

	bool isEdgeOverlap(const ourEdge &e) {
		ourEdge edge(0, 1);
		for (int i = 0; i < 2; ++i) {
			for (auto &elem : vertAdjEdges[e[i]]) {
				edge = *elem;
				if (edge[0] == e[!i])
					continue;
				if (edge[1] == e[!i])
					continue;
				Point3D a = pts[edge[0]], b = pts[edge[1]], c, o = pts[e[i]];
				Plane3D plane(a, b, o);
				c = plane.projection(pts[e[!i]]);

				Vector3D OA = a - o, OC = c - o, OB = b - o;
				// OA = OA / CGAL::sqrt(OA.squared_length());
				// OB = OB / CGAL::sqrt(OA.squared_length());
				// OC = OC / CGAL::sqrt(OA.squared_length());

				auto d1 = CGAL::cross_product(OA, OC) * CGAL::cross_product(OA, OB),
					 d2 = CGAL::cross_product(OB, OC) * CGAL::cross_product(OB, OA);

				if (d1 > 1e-4 && d2 > 1e-4) {
					// cout << e << " rejected by " << ourFace(edge[0], edge[1], e[i]) <<
					// "\n";
					return true;
				}
			}
		}
		// cout << "pass" << std::endl;
		return false;
	}

	bool validToAdd(ourEdge edge, int newPoint, int &toBeRemoved) {
		auto &edgeDeg = edgeDegree[edge];
		if (edgeDeg.size() == 2) {
			println(edge[0], edge[1], edgeDeg[0], edgeDeg[1], newPoint);

			Kernel::FT a1, a2, a3;

			Vector3D norm1 = CGAL::normal(pts[edge[0]], pts[edge[1]], pts[edgeDeg[0]]),
					 norm2 = CGAL::normal(pts[edge[1]], pts[edge[0]], pts[edgeDeg[1]]);
			norm1 = norm1 / CGAL::sqrt(norm1.squared_length());
			norm2 = norm2 / CGAL::sqrt(norm2.squared_length());
			a1 = std::acos(norm1 * norm2);

			norm1 = CGAL::normal(pts[edge[0]], pts[edge[1]], pts[newPoint]);
			norm2 = CGAL::normal(pts[edge[1]], pts[edge[0]], pts[edgeDeg[0]]);
			norm1 = norm1 / CGAL::sqrt(norm1.squared_length());
			norm2 = norm2 / CGAL::sqrt(norm2.squared_length());
			a2 = std::acos(norm1 * norm2);

			norm1 = CGAL::normal(pts[edge[0]], pts[edge[1]], pts[newPoint]);
			norm2 = CGAL::normal(pts[edge[1]], pts[edge[0]], pts[edgeDeg[1]]);
			norm1 = norm1 / CGAL::sqrt(norm1.squared_length());
			norm2 = norm2 / CGAL::sqrt(norm2.squared_length());
			a3 = std::acos(norm1 * norm2);

			if (a2 < a3) {
				if (a1 > a2 && !isLocked(ourFace(edge[0], edge[1], edgeDeg[1]))) {
					toBeRemoved = edgeDeg[1];
					println("Angle", a1, "->", a2);
					println(edge, newPoint, toBeRemoved);
					println("Edge Degree=2, returning true");
					return true;
				}
			}
			else {
				if (a1 > a3 && !isLocked(ourFace(edge[0], edge[1], edgeDeg[0]))) {
					toBeRemoved = edgeDeg[0];
					println("Angle", a1, "->", a3);
					println(edge, newPoint, toBeRemoved);
					println("Edge Degree=2, returning true");
					return true;
				}
			}
			println("Angle", 180 * a1 / M_PI, "->", 180 * a2 / M_PI, ",", 180 * a3 / M_PI);
			println("Edge Degree=2, returning false");
			return false;
		}
		if (edgeDeg.size() == 0) {
			println("Edge Degree=0, returning true");
			return true;
		}

		Vector3D norm1 = CGAL::normal(pts[edge[0]], pts[edge[1]], pts[newPoint]),
				 norm2 = CGAL::normal(pts[edge[1]], pts[edge[0]], pts[edgeDeg[0]]);
		norm1 = norm1 / CGAL::sqrt(norm1.squared_length());
		norm2 = norm2 / CGAL::sqrt(norm2.squared_length());
		println("Edge Degree=1, returning ", (norm1 * norm2));
		return norm1 * norm2 >= 0.0;
		// return norm1 * norm2 >= 0.707106781187;
	}

	set<ourEdge> fillHoles(set<int> &leftVerts, map<int, set<int>> &tempAdjList) {
		set<ourEdge> newEdges;

		if (leftVerts.size() == 3) {
			int p[3], i = 0;
			for (auto &elem : leftVerts) {
				p[i++] = elem;
			}
			addFaceToModel({p[0], p[1], p[2]});
			return newEdges;
		}

		// if (leftVerts.size() == 4) {
		// 	int p[4], i = 0;
		// 	for (auto &elem : leftVerts) {
		// 		p[i++] = elem;
		// 	}

		// 	if (modelEdges.find({p[0], p[1]}) == modelEdges.end()) {
		// 		addFaceToModel({p[0], p[1], p[2]});
		// 		addFaceToModel({p[0], p[1], p[3]});
		// 	}
		// 	else if (modelEdges.find({p[0], p[2]}) == modelEdges.end()) {
		// 		addFaceToModel({p[0], p[1], p[2]});
		// 		addFaceToModel({p[0], p[2], p[3]});
		// 	}
		// 	else {
		// 		addFaceToModel({p[0], p[1], p[3]});
		// 		addFaceToModel({p[0], p[2], p[3]});
		// 	}
		// 	return newEdges;
		// }

		vector<std::pair<Kernel::FT, ourEdge>> allEdges;
		vector<bool> covered(pts.size(), false);
		//         Kernel::FT avgCount = 0, avgSum = 0, lastAvg = 0, lastEdge = 0;

		println("Lengths");
		for (auto &edge : possibleEdges) {
			int u = edge[0], v = edge[1];
			if (leftVerts.find(u) != leftVerts.end() &&
				leftVerts.find(v) != leftVerts.end() && !isEdgeOverlap(edge)) {
				allEdges.push_back({sqDist(pts[u], pts[v]), edge});
				println(edge);
			}
			// print(CGAL::sqrt(sqDist(pts[u], pts[v])));
		}

		sortAndRemoveDuplicate(allEdges);
		allEdges.shrink_to_fit();

		// for (auto &elem : allEdges) {
		// 	println(elem.second);
		// }
		// println();

		/*		if (allEdges.size() >= 2) {
newEdges.insert(allEdges[0].second);
covered[allEdges[0].second[0]] = true;
covered[allEdges[0].second[1]] = true;

newEdges.insert(allEdges[1].second);
covered[allEdges[1].second[0]] = true;
covered[allEdges[1].second[1]] = true;

lastAvg = avgSum = allEdges[1].first - allEdges[0].first;
lastEdge = allEdges[1].first;
avgCount++;
}*/
		// println(avgSum);
		for (int i = 0; i < allEdges.size(); ++i) {
			auto &elem = allEdges[i];
			auto u = elem.second[0], v = elem.second[1];
			print(u, v);
			print(!covered[u], !covered[v], elem.first <= 4 * maxLength[u],
				  elem.first <= 4 * maxLength[v]);
			if ((!covered[u] || !covered[v]) && elem.first <= 4 * maxLength[u] &&
				elem.first <= 4 * maxLength[v]) {
				// print(elem.first, elem.second, covered[u], covered[v]);
				newEdges.insert(elem.second);
				print("added");
				covered[u] = true;
				covered[v] = true;
				// for (auto &elem : tempAdjList[u]) {
				// 	covered[elem] = true;
				// }
				// for (auto &elem : tempAdjList[v]) {
				// 	covered[elem] = true;
				// }
				// print("||", covered[u], covered[v]);
				// avgSum += allEdges[i].first - lastEdge;
				// avgCount++;
				// print(allEdges[i].first - lastEdge, avgSum / avgCount, ((avgCount -
				// 1) * avgSum / avgCount - lastAvg) / lastAvg);
				// if ((((avgCount - 1) * avgSum / avgCount - lastAvg) / lastAvg) > 0.1)
				// 	break;
				// lastAvg = avgSum;
				// lastEdge = allEdges[i].first;
				// println();
			}
			println();
		}

		return newEdges;
	}

	void reconstruct(const vector<ourEdge> &initialEdges) {
		vector<set<int>> currAdjList = getAdjList(pts.size(), initialEdges);
		set<ourEdge> origEdges(initialEdges.begin(), initialEdges.end());

		maxEdge = 0;
		// maxMinRatio = 1;

		println(val(initialEdges.size()));

		for (auto &elem : initialEdges) {
			maxEdge = std::max(maxEdge, sqDist(pts[elem[0]], pts[elem[1]]));
		}

		maxEdge *= Kernel::FT(1.0001);

		cerr << CGAL::sqrt(maxEdge) << "\n";

		set<std::tuple<Kernel::FT, ourEdge, int>,
			std::greater<std::tuple<Kernel::FT, ourEdge, int>>>
			pq, nextPQ;
		// modelEdges.insert(initialEdges.begin(), initialEdges.end());
		int edgeCondition = 2;
		// for (auto &elem : pq) {
		// 	cout << "[" << get<0>(elem) << " " << get<1>(elem) << " " <<
		// get<2>(elem) << "],\n";
		// }
		// cout << "\n";
		int lastFaces = -1;
		set<ourEdge> nextEdges;
		int tempModelCount = 1;

		while (true) {
			// Validity Check
			for (auto &face : modelFaces) {
				if (modelEdges.find(ourEdge(face[0], face[1])) == modelEdges.end()) {
					println("Validity Check", face, 0, 1);
					exit(0);
				}
				if (modelEdges.find(ourEdge(face[1], face[2])) == modelEdges.end()) {
					println("Validity Check", face, 1, 2);
					exit(0);
				}
				if (modelEdges.find(ourEdge(face[0], face[2])) == modelEdges.end()) {
					println("Validity Check", face, 0, 2);
					exit(0);
				}
			}

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
			for (auto &elem : origEdges) {
				// print(elem, ":");
				// for (auto &point : edgeDegree[elem]) {
				// 	print(point);
				// }
				// println();
				if (edgeDegree[elem].size() <= 1) {
					leftVerts.insert({elem[0], false});
					leftVerts.insert({elem[1], false});
					tempAdjList[elem[0]].insert(elem[1]);
					tempAdjList[elem[1]].insert(elem[0]);
				}

				// leftVerts.insert(elem[0]);
				// leftVerts.insert(elem[1]);
			}
			// println();
			println(val(leftVerts.size()));
			// cout << "leftVerts=" << leftVerts.size() << "\n";
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

			print("Size Change", possibleEdges.size());
			if (leftVerts.size()) {
				for (auto itr = possibleEdges.begin(); itr != possibleEdges.end();) {
					int u = (*itr)[0], v = (*itr)[1];
					if (leftVerts.find(u) == leftVerts.end() ||
						leftVerts.find(v) == leftVerts.end())
						itr = possibleEdges.erase(itr);
					else
						++itr;
				}
			}
			println(possibleEdges.size());

			// allEdges

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
					if (leftVerts[u] == true)
						continue;
					leftVerts[u] = true;
					hole.insert(u);
					// print(u, ":");
					for (auto &v : tempAdjList[u]) {
						if (!leftVerts[v]) {
							// print(v);
							Q.push(v);
						}
					}
					// println();
				}
				print("Considering Hole of size:", hole.size(), ":");
				// for (auto &elem : hole) {
				// 	print(elem);
				// }
				println();
				// if (hole.size() == 2) {
				// 	int u = *hole.begin();
				// 	hole.erase(hole.begin());
				// 	int v = *hole.begin();
				// 	for (auto &elem : edgeDegree[ourEdge(u, v)]) {
				// 		cout << elem << " ";
				// 	}
				// 	cout << std::endl;
				// 	return;
				// }
				auto E = fillHoles(hole, tempAdjList);
				println("Hole Edges =", E.size());
				nextEdges.insert(E.begin(), E.end());
			}

			// println("nextEdges Start");
			// for (auto &elem : nextEdges) {
			// 	println(elem);
			// }
			// println("nextEdges End");

			println("Generated tempMod" + std::to_string(tempModelCount) + ".off");
			writeModel("tempMod" + std::to_string(tempModelCount++) + ".off");

			println("Iteration #=", 2 - edgeCondition + 1);
			for (auto &edge : nextEdges) {
				// bool inserted = false;
				for (auto &elem : getFacesFromEdge(edge, currAdjList, edgeCondition)) {
					pq.insert(std::make_tuple(elem.second, edge, elem.first));
					// inserted = true;
					// cout << "Initially added " << edge << " and " << elem.first << std::endl;
				}
			}
			for (auto &edge : origEdges) {
				// bool inserted = false;
				for (auto &elem : getFacesFromEdge(edge, currAdjList, edgeCondition)) {
					pq.insert(std::make_tuple(elem.second, edge, elem.first));
					// inserted = true;
					// cout << "Initially added " << edge << " and " << elem.first << std::endl;
				}
			}

			cerr << std::endl
				 << "Iteration #" << 2 - edgeCondition + 1 << std::endl;
			float total = pq.size();
			progress_bar bar;
			while (!pq.empty()) {
				progress = 1 - pq.size() / total;
				// if (pq.size() < total)
				// showProgress(1 - pq.size() / total);
				// for (auto &elem : pq) {
				// cout << "[" << get<0>(elem) << " " << get<1>(elem) << " " <<
				// get<2>(elem) << "]," << std::endl;
				// }
				// cout << std::endl;
				// for (int i = 0; i < pts.size(); ++i) {
				// 	cout << i << ": [";
				// 	for (auto elem : currAdjList[i]) {
				// 		cout << elem << " ";
				// 	}
				// 	cout << "]" << std::endl;
				// }
				// for (auto &elem : edgeDegree) {
				// 	if (elem.second.size() >= 1) {
				// 		cout << elem.first << ": [";
				// 		for (auto &p : elem.second)
				// 			cout << p << " ";
				// 		cout << "]" << std::endl;
				// 	}
				// }*/
				Kernel::FT score = get<0>(*pq.begin());
				ourEdge edge = get<1>(*pq.begin());
				int point = get<2>(*pq.begin());
				pq.erase(pq.begin());

				println(score, edge, point);

				ourFace face(edge[0], edge[1], point);

				if (score <= -0.9999)
					continue;
				if (modelFaces.find(face) != modelFaces.end())
					continue;
				// if (edgeDegree[edge].size() == 2) {
				// 	println("edgeDegree[edge].size() == 2");
				// 	continue;
				// }
				if (isEdgeOverlap(edge)) {
					println("isEdgeOverlap(edge)");
					continue;
				}
				bool passed = true;
				for (int i = 0; passed && i < 3; ++i) {
					for (auto &point : currAdjList[face[i]]) {
						if (faceEdgeOverlap(face[(i + 1) % 3], face[(i + 2) % 3], face[i],
											point)) {
							println(face, "rejected by", ourEdge(face[i], point));
							passed = false;
							break;
						}
					}
				}
				if (!passed)
					continue;
				ourEdge newEdge1(edge[0], point), newEdge2(edge[1], point);
				int remove = -1, remove1 = -1, remove2 = -1;
				//cout << get<0>(*pq.begin()) << " " << edge << " " << point << "::";
				println(
					(modelFaces.find(face) == modelFaces.end()),
					isFaceDelaunay(face),
					validToAdd(edge, point, remove),
					validToAdd(newEdge1, edge[1], remove1),
					validToAdd(newEdge2, edge[0], remove2),
					(remove1 != -1 || !isEdgeOverlap(newEdge1)),
					(remove2 != -1 || !isEdgeOverlap(newEdge2)),
					(testEdges(pts[face[0]], pts[face[1]], pts[face[2]], maxEdge) >= edgeCondition));
				if (modelFaces.find(face) == modelFaces.end() &&
					isFaceDelaunay(face) &&
					validToAdd(edge, point, remove) &&
					validToAdd(newEdge1, edge[1], remove1) &&
					validToAdd(newEdge2, edge[0], remove2) &&
					(remove1 != -1 || !isEdgeOverlap(newEdge1)) &&
					(remove2 != -1 || !isEdgeOverlap(newEdge2))) {

					if (testEdges(pts[face[0]], pts[face[1]], pts[face[2]], maxEdge) < edgeCondition) {
						nextPQ.insert(std::make_tuple(score, edge, point));
						continue;
					}

					if (remove != -1) {
						println("Remove", ourFace(edge[0], edge[1], remove), "due to", ourFace(edge[0], edge[1], point));
						removeFace(edge, remove);
					}
					if (remove1 != -1) {
						println("Remove", ourFace(newEdge1[0], newEdge1[1], remove1), "due to", ourFace(newEdge1[0], newEdge1[1], edge[1]));
						removeFace(newEdge1, remove1);
					}
					if (remove2 != -1) {
						println("Remove", ourFace(newEdge2[0], newEdge2[1], remove2), "due to", ourFace(newEdge2[0], newEdge2[1], edge[0]));
						removeFace(newEdge2, remove2);
					}

					addFaceToModel(face);
					origEdges.erase(edge);
					origEdges.erase(newEdge1);
					origEdges.erase(newEdge2);

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
							// cout << edge << " added " << ourEdge(u, v) << " and " <<
							// elem.first << "\n";
						}
					}
				}
			}

			for (auto &elem : nextPQ) {
				if (edgeDegree[get<1>(elem)].size() == 1)
					pq.insert(elem);
			}

			nextPQ.clear();

			bar.terminate();

			if (lastFaces == modelFaces.size())
				break;
			lastFaces = modelFaces.size();

			edgeCondition--;
			println("Generated tempMod" + std::to_string(tempModelCount) + ".off");
			writeModel("tempMod" + std::to_string(tempModelCount++) + ".off");

			// modelEdges.clear();
			// break;
		}
		// modelEdges.swap(nextEdges);
		// println(val(maxMinRatio));
	}

	vector<vector<int>> getAllCycles(map<int, set<int>> &tempAdjList,
									 set<int> &leftPts) {
		vector<vector<int>> allCycles;
		map<int, bool> inCycle;
		for (int point : leftPts) {
			if (!inCycle[point]) {
				pushCycles(point, tempAdjList, allCycles, inCycle);
			}
		}
		return allCycles;
	}

	void pushCycles(int point, map<int, set<int>> &tempAdjList,
					vector<vector<int>> &allCycles, map<int, bool> &inCycle) {

		stack<int> s;
		vector<int> eu; // for the euler's circuit
		int a = point;

		while (!s.empty() || tempAdjList[a].size() != 0) {
			if (tempAdjList[a].size()) {
				s.push(a);
				int t = *tempAdjList[a].begin();
				tempAdjList[a].erase(t);
				tempAdjList[t].erase(a);
				a = t;
			}
			else {
				eu.push_back(a);
				a = s.top();
				s.pop();
			}
		}

		eu.push_back(point);

		map<int, bool> visited;
		// print("Euler Cycle: ");
		// for (auto &a : eu) {
		// print(a);
		// }
		// cout << std::endl;

		for (int &a : eu) {
			if (!visited[a]) {
				s.push(a);
				visited[a] = true;
			}
			else {
				vector<int> cycle;
				cycle.push_back(a);
				inCycle[a] = true;
				while (s.top() != a) {
					cycle.push_back(s.top());
					inCycle[s.top()] = true;
					visited[s.top()] = false;
					// cout << "Popping " << s.top() << std::endl;
					s.pop();
				}
				allCycles.push_back(cycle);
			}
		}
	}

	void postProcess() {

		println("Doing Post Processing");

		set<int> leftVerts;
		map<int, set<int>> tempAdjList;

		for (auto &elem : modelEdges) {
			if (edgeDegree[elem].size() == 1) {
				leftVerts.insert(elem[0]);
				leftVerts.insert(elem[1]);
				tempAdjList[elem[0]].insert(elem[1]);
				tempAdjList[elem[1]].insert(elem[0]);
				println(elem);
			}
		}

		// for (auto elem : leftVerts) {
		// 	print(elem);
		// }
		println();

		vector<vector<int>> cycles = getAllCycles(tempAdjList, leftVerts);

		for (auto &cycle : cycles) {
			print("Cycle of Size:", cycle.size(), ":");
			// for (auto &p : cycle) {
			// 	print(p);
			// }
			println();

			if (cycle.size() == 3) {
				addFaceToModel({cycle[0], cycle[1], cycle[2]});
			}

			// if (cycle.size() == 4) {
			// 	if (modelEdges.find({cycle[0], cycle[1]}) == modelEdges.end()) {
			// 		addFaceToModel({cycle[0], cycle[1], cycle[2]});
			// 		addFaceToModel({cycle[0], cycle[1], cycle[3]});
			// 	}
			// 	else if (modelEdges.find({cycle[0], cycle[2]}) == modelEdges.end()) {
			// 		addFaceToModel({cycle[0], cycle[1], cycle[2]});
			// 		addFaceToModel({cycle[0], cycle[2], cycle[3]});
			// 	}
			// 	else {
			// 		addFaceToModel({cycle[0], cycle[1], cycle[3]});
			// 		addFaceToModel({cycle[0], cycle[2], cycle[3]});
			// 	}
			// }
		}
	}

  public:
	SurfaceReconstruct(const char *inputFilePath) {
		std::ifstream inputFile(inputFilePath);

		auto start = std::chrono::high_resolution_clock::now();
		if (!CGAL::read_xyz_points(
				inputFile, back_inserter(pts))) { // output iterator over points
			cerr << "Error: cannot read file.";
			return;
		}

		// for (auto &elem : pts) {
		// 	elem = Point3D(100.0 * elem.x(), 100.0 * elem.y(), 100.0 * elem.z());
		// }

		sortAndRemoveDuplicate(pts);
		pts.shrink_to_fit();

		maxLength.resize(pts.size(), 0);

		auto finish = std::chrono::high_resolution_clock::now();
		cout << pts.size() << " points read in "
			 << std::chrono::duration<double>(finish - start).count() << " secs"
			 << std::endl;

		start = std::chrono::high_resolution_clock::now();
		dt.insert(pts.begin(), pts.end());
		if (!dt.is_valid(true)) {
			cerr << "Error: fail to build a Delaunay triangulation." << std::endl;
			valid = false;
			return;
		}
		if (dt.dimension() != 3) {
			cerr << "Error: cannot built a 3D triangulation.\n Current dimension = "
				 << dt.dimension() << std::endl;
			valid = false;
			return;
		}
		finish = std::chrono::high_resolution_clock::now();
		cout << "Delaunay Triangulation created in "
			 << std::chrono::duration<double>(finish - start).count() << " secs"
			 << std::endl;

		valid = true;

		vertAdjEdges.resize(pts.size());

		for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
			Triangle3D tri = dt.triangle(*faceItr);
			possibleFaces.insert(ourFace(getIndex(pts, tri[0]),
										 getIndex(pts, tri[1]),
										 getIndex(pts, tri[2])));
		}

		for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
			Segment3D segment = dt.segment(*edgeItr);
			possibleEdges.insert(ourEdge(getIndex(pts, segment[0]),
										 getIndex(pts, segment[1])));
		}

		// for (auto &edge : possibleEdges) {
		// 	println(edge);
		// }

		start = std::chrono::high_resolution_clock::now();
		vector<ourEdge> initialEdges = getMstEdges();
		finish = std::chrono::high_resolution_clock::now();
		cout << "MST created in "
			 << std::chrono::duration<double>(finish - start).count() << " secs"
			 << std::endl;

		for (auto &edge : initialEdges) {
			maxLength[edge[0]] =
				std::max(maxLength[edge[0]], sqDist(pts[edge[0]], pts[edge[1]]));
			maxLength[edge[1]] =
				std::max(maxLength[edge[1]], sqDist(pts[edge[0]], pts[edge[1]]));
		}

		start = std::chrono::high_resolution_clock::now();
		// modelEdges.insert(initialEdges.begin(), initialEdges.end());
		reconstruct(initialEdges);

		postProcess();

		modelEdges.clear();
		// modelEdges.insert(initialEdges.begin(), initialEdges.end());

		finish = std::chrono::high_resolution_clock::now();
		cout << std::endl
			 << pts.size() << " Reconstructed in "
			 << std::chrono::duration<double>(finish - start).count() << std::endl;

		std::ofstream outputFile("MST.txt");
		outputFile << pts.size() << "\n";
		for (Point3D &point : pts) {
			outputFile << point << "\n";
		}
		outputFile << initialEdges.size() << "\n";
		for (auto &edge : initialEdges) {
			outputFile << edge[0] << " " << edge[1] << "\n";
		}
	}

	bool isValid() { return valid; }

	void writeModel(std::string outputFilePath) {
		std::ofstream outputFile(outputFilePath);
		outputFile << "OFF\n";
		outputFile << pts.size() << " " << modelFaces.size() << " 0\n";
		for (Point3D &point : pts) {
			outputFile << point << "\n";
		}
		for (auto &tri : modelFaces) {
			outputFile << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
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
