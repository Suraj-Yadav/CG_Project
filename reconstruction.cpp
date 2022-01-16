#define _USE_MATH_DEFINES
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Union_find.h>

#include "util.h"

#define sqDist(x, y) CGAL::squared_distance(x, y)

// Alias for CGAL dataTypes
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Plane_3 Plane3D;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Vector_3 Vector3D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Triangle_3 Triangle3D;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;

const Kernel::FT inf = std::numeric_limits<Kernel::FT>::infinity();

/// @file
/// @brief This file is marvelous.

/////////////////////////////////////////////////
/// \brief
/// Generates adjacency list from given Edges.
/// \param[in] pointsCount int  Number of points.
/// \param[in] edges const vector<ourEdge>& The edges.
/// \return vector<set<int>> Adjacency List
///
/////////////////////////////////////////////////
vector<set<int>> getAdjList(int pointsCount, const vector<ourEdge> &edges) {
	vector<set<int>> adjList(pointsCount);
	for (auto &edge : edges) {
		adjList[edge[0]].insert(edge[1]);
		adjList[edge[1]].insert(edge[0]);
	}
	return adjList;
}
float progress = 0;

/////////////////////////////////////////////////
/// \brief Implements a Progress bar which runs in a separate thread.
/////////////////////////////////////////////////
class progress_bar {
	std::atomic<bool> finished;
	std::thread t;

  public:
	progress_bar()
		: finished(false), t(std::ref(*this)) {
	}
	/////////////////////////////////////////////////
	/// \brief
	/// initiate the bar by starting a new thread
	/// \return void
	///
	/////////////////////////////////////////////////
	void operator()() {
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
	/////////////////////////////////////////////////
	/// \brief
	/// Stop the thread
	/// \return void
	///
	/////////////////////////////////////////////////
	void terminate() {
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
		cerr << "] " << int(100.0) << " %\n";
	}
};

/////////////////////////////////////////////////
/// \brief Contains the Reconstructed Surface
/////////////////////////////////////////////////
class SurfaceReconstruct {
	vector<Point3D> pts;								 ///< set of distinct points
	vector<Kernel::FT> maxLength;						 ///< stores length max incident edge on the vertex i
	DT3 dt;												 ///< Delaunay Triangulation
	bool valid;											 ///< Validity of Delaunay Triangulation/surface
	vector<ourEdge> initialEdges;						 ///< MST edges
	set<ourFace> possibleFaces;							 ///< All Delaunay Faces
	set<pair<Kernel::FT, ourEdge>> possibleEdges;		 ///< All Delaunay edges - (modelEdges + initialEdges) sorted by length
	set<ourEdge> modelEdges;							 ///< Current edges on the surface
	set<ourFace> modelFaces;							 ///< Current faces on the surface
	map<ourEdge, vector<int>> edgeDegree;				 ///< set of vertices that form face with edge e
	vector<set<int>> currAdjList;						 ///< Adjacency list of the current surface
	Kernel::FT maxEdge;									 ///< Length of maximum edge in MST
	vector<vector<set<ourEdge>::iterator>> vertAdjEdges; ///< set of edges that form face with vertex i

	/////////////////////////////////////////////////
	/// \brief
	/// Computes the MST edges from the possibleEdges
	/// \return vector<ourEdge> MST edges
	///
	/////////////////////////////////////////////////
	vector<ourEdge> getMstEdges() {
		vector<CGAL::Union_find<int>::handle> handle;
		CGAL::Union_find<int> uf;
		handle.reserve(pts.size());
		for (int i = 0; i < pts.size(); i++)
			handle.push_back(uf.make_set(i));

		vector<ourEdge> mst;

		for (auto &elem : possibleEdges) {
			int u = elem.second[0], v = elem.second[1];
			if (!uf.same_set(handle[u], handle[v])) {
				mst.push_back(elem.second);
				uf.unify_sets(handle[u], handle[v]);
			}
		}

		for (auto &elem : mst) {
			eraseFromPossibleEdges(elem);
		}

		std::sort(mst.begin(), mst.end());
		return mst;
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Determines whether  \p face is a Delaunay Face
	/// \param face const ourFace& Input Face
	/// \return bool
	///
	/////////////////////////////////////////////////
	inline bool isFaceDelaunay(const ourFace &face) {
		return possibleFaces.find(face) != possibleFaces.end();
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Returns the score (cosine of the maximum angle) of Triangle.
	/// \param u, v, w int indexes of vertex of Triangle
	/// \return Kernel::FT  score
	///
	/////////////////////////////////////////////////
	Kernel::FT getTriScore(int u, int v, int w) {
		Point3D a = pts[u], b = pts[v], c = pts[w];
		Vector3D AB(a, b), AC(a, c), BC(b, c);
		AB = AB / CGAL::sqrt(AB.squared_length());
		AC = AC / CGAL::sqrt(AC.squared_length());
		BC = BC / CGAL::sqrt(BC.squared_length());
		Kernel::FT score[3] = {AB * AC, AC * BC, -AB * BC};
		if (score[0] > score[1]) std::swap(score[0], score[1]);
		if (score[1] > score[2]) std::swap(score[1], score[2]);
		if (score[0] > score[1]) std::swap(score[0], score[1]);
		// if (score[0] > 0.0)
		// return (score[0] + score[1]) / 2.0;
		return score[0];
	}

	/////////////////////////////////////////////////
	/// \brief
	///
	/// \param u int
	/// \param v int
	/// \param w int
	/// \return Kernel::FT
	///
	/////////////////////////////////////////////////
	Kernel::FT getMaxAngle(int u, int v, int w) {
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

	/////////////////////////////////////////////////
	/// \brief
	/// Returns cosine of the angle \p uvw
	/// \param u, v, w int indexes of vertex of Triangle
	/// \return Kernel::FT
	///
	/////////////////////////////////////////////////
	Kernel::FT getAngle(int u, int v, int w) {
		Point3D a = pts[u], b = pts[v], c = pts[w];
		Vector3D AB(a, b), BC(b, c);
		AB = AB / CGAL::sqrt(AB.squared_length());
		BC = BC / CGAL::sqrt(BC.squared_length());
		return -AB * BC;
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Add \p face to the mesh and updates the auxiliary data structures of object.
	/// \return void
	///
	/////////////////////////////////////////////////
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

		double len[3] = {sqDist(pts[e1[0]], pts[e1[1]]), sqDist(pts[e2[0]], pts[e2[1]]), sqDist(pts[e3[0]], pts[e3[1]])};

		maxLength[p1] = std::max(std::max(len[1], len[2]), maxLength[p1]);
		maxLength[p2] = std::max(std::max(len[0], len[2]), maxLength[p2]);
		maxLength[p3] = std::max(std::max(len[1], len[0]), maxLength[p2]);

		if (!contains(initialEdges, e1))
			eraseFromPossibleEdges(e1);
		if (!contains(initialEdges, e2))
			eraseFromPossibleEdges(e2);
		if (!contains(initialEdges, e3))
			eraseFromPossibleEdges(e3);

		edgeDegree[e1].push_back(p1);
		edgeDegree[e2].push_back(p2);
		edgeDegree[e3].push_back(p3);

		vertAdjEdges[p1].push_back(modelEdges.find(e1));
		vertAdjEdges[p2].push_back(modelEdges.find(e2));
		vertAdjEdges[p3].push_back(modelEdges.find(e3));
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Remove \p face to the mesh and updates the auxiliary data structures of object.
	/// \return void
	///
	/////////////////////////////////////////////////
	void removeFace(ourEdge edge, int remove) {
		ourEdge edge1(ourEdge(edge[0], remove)), edge0(ourEdge(edge[1], remove));

		modelFaces.erase(ourFace(edge, remove));
		possibleFaces.insert(ourFace(edge, remove));

		deleteFromVector(vertAdjEdges[remove], modelEdges.find(edge));
		deleteFromVector(vertAdjEdges[edge[0]], modelEdges.find(edge0));
		deleteFromVector(vertAdjEdges[edge[1]], modelEdges.find(edge1));

		deleteFromVector(edgeDegree[edge], remove);
		deleteFromVector(edgeDegree[edge1], edge[1]);
		deleteFromVector(edgeDegree[edge0], edge[0]);

		if (edgeDegree[edge1].size() == 0) {
			modelEdges.erase(edge1);
			addToPossibleEdges(edge1);
		}
		if (edgeDegree[edge0].size() == 0) {
			modelEdges.erase(edge0);
			addToPossibleEdges(edge0);
		}
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Count the number of edges in Triangle \p abc with length <= \p length
	/// \param a,b,c Point3D Vertices of Triangle
	/// \param length Kernel::FT
	/// \return int
	///
	/////////////////////////////////////////////////
	static int testEdges(Point3D a, Point3D b, Point3D c, Kernel::FT length) {
		int ans = 0;
		if (sqDist(a, b) <= length)
			ans++;
		if (sqDist(a, c) <= length)
			ans++;
		if (sqDist(c, b) <= length)
			ans++;
		return ans;
	}

	void addToPossibleEdges(ourEdge &e) {
		possibleEdges.insert(make_pair(sqDist(pts[e[0]], pts[e[1]]),
									   e));
	}

	void eraseFromPossibleEdges(ourEdge &e) {
		possibleEdges.erase(make_pair(sqDist(pts[e[0]], pts[e[1]]),
									  e));
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Generate all possible Delaunay Triangles with their score adjacent to \p edge
	/// \param edge const ourEdge&
	/// \return vector<pair<int, Kernel::FT>>
	///
	/////////////////////////////////////////////////
	vector<pair<int, Kernel::FT>> getFacesFromEdge(const ourEdge &edge) {
		vector<pair<int, Kernel::FT>> elems;
		// set<int> t, f;
		// cout << edge << "\n";
		// println("Getting Adj Faces");
		for (int i = 0; i < 2; ++i) {
			for (auto &point : currAdjList[edge[i]]) {
				// println(point, ":", (point != edge[!i]),
				// 		isFaceDelaunay(ourFace(point, edge[i], edge[!i])),
				// 		(modelFaces.find(ourFace(point, edge[i], edge[!i])) == modelFaces.end()));
				if (point != edge[!i] &&
					isFaceDelaunay(ourFace(point, edge[i], edge[!i])) &&
					!contains(modelFaces, ourFace(point, edge[i], edge[!i]))) { // &&
																				// testEdges(pts[edge[0]], pts[edge[1]], pts[point], maxEdge) >=							  // edgeCondition) {
					auto score = getTriScore(point, edge[i], edge[!i]);
					if (score <= -0.8660254037844387)
						continue;
					// if (getMaxAngle(point, edge[i], edge[!i]) <= -0.999)
					// 	continue;
					// print("added");
					// println(edge);
					elems.push_back(std::make_pair(point, score));
					// t.insert(point);
				}
				// else if (point != edge[1]) {
				// 	f.insert(point);
				// }
				// println();
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
		// sortAndRemoveDuplicate(elems);
		return elems;
	}

	// bool faceEdgeOverlap(int A, int B, int O, int D) {
	// if (A == D || B == D)
	// return false;
	// Point3D a = pts[A], b = pts[B], c, o = pts[O];
	// Plane3D plane(a, b, o);
	// c = plane.projection(pts[D]);

	// Vector3D OA = a - o, OC = c - o, OB = b - o;
	// OA = OA / CGAL::sqrt(OA.squared_length());
	// OB = OB / CGAL::sqrt(OA.squared_length());
	// OC = OC / CGAL::sqrt(OA.squared_length());

	// auto d1 = CGAL::cross_product(OA, OC) * CGAL::cross_product(OA, OB),
	// d2 = CGAL::cross_product(OB, OC) * CGAL::cross_product(OB, OA);
	// if (d1 > 0 && d2 > 0) {
	// cout << e << " rejected by " << ourFace(edge[0], edge[1], e[i]) << "\n";
	// println(val(d1), val(d2));
	// return true;
	// }
	// return false;
	// }

	/////////////////////////////////////////////////
	/// \brief
	/// Check if \p face is surrounded by Triangles on all side
	/// \param face ourFace
	/// \return int
	///
	/////////////////////////////////////////////////
	int lockedLevel(ourFace face) {
		int count = 0;
		for (int i = 0; i < 3; ++i) {
			if (edgeDegree[ourEdge(face[i], face[(i + 1) % 3])].size() >= 2)
				count++;
		}
		return count;
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Check if Segment \p OC projects on Triangle \p ABO
	/// \param A, B, O, C int indexes
	/// \return bool
	///
	/////////////////////////////////////////////////
	bool testProjection(int A, int B, int O, int C) {
		Point3D a = pts[A], b = pts[B], o = pts[O];
		Plane3D plane(a, b, o);
		Point3D c = plane.projection(pts[C]);
		Vector3D OA = a - o, OC = c - o, OB = b - o;
		auto d1 = CGAL::cross_product(OA, OC) * CGAL::cross_product(OA, OB),
			 d2 = CGAL::cross_product(OB, OC) * CGAL::cross_product(OB, OA);
		return d1 > 0 && d2 > 0;
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Check if face formed by \p e and \p point overlaps with any face or edge at a distance <= 2.
	/// \return bool
	///
	/////////////////////////////////////////////////
	bool isFaceOverlap(const ourEdge &e, int point, set<ourFace> excludeSet = set<ourFace>(),
					   set<ourEdge> excludeEdge = set<ourEdge>()) {
		bool ret = isEdgeOverlap(e, excludeSet);
		if (ret)
			return true;
		set<int> testingEdges;
		for (auto &edgeItr : vertAdjEdges[point]) {
			ourEdge edge = *edgeItr;
			if (contains(excludeSet, ourFace(edge, point)))
				continue;
			if (edge[0] != e[0] && edge[0] != e[1])
				testingEdges.insert(edge[0]);
			if (edge[1] != e[0] && edge[1] != e[1])
				testingEdges.insert(edge[1]);
			// for (int i = 0; i < 2; ++i) {
			// 	if (testProjection(edge[0], edge[1], point, e[i])) {
			// 		println(e, "rejected by", edge);
			// 		return true;
			// 	}
			// }
		}
		for (auto &elem : currAdjList[point]) {
			if (elem == e[0] || elem == e[1])
				continue;
			if (contains(excludeEdge, ourEdge(elem, point)))
				continue;
			testingEdges.insert(elem);
		}
		for (auto &elem : testingEdges) {
			if (testProjection(e[0], e[1], point, elem)) {
				println(e, "rejected by", ourEdge(elem, point));
				return true;
			}
		}

		for (auto &elem : currAdjList[point]) {
			if (elem == e[0] || elem == e[1])
				continue;

			for (auto &edgeItr : vertAdjEdges[elem]) {

				ourEdge adjEdge = *edgeItr;

				if (adjEdge == ourEdge(point, e[0]) || adjEdge == ourEdge(point, e[1]))
					continue;

				if (point == adjEdge[0] || point == adjEdge[1])
					continue;
				if (e[0] == adjEdge[0] || e[0] == adjEdge[1] || e[1] == adjEdge[0] || e[1] == adjEdge[1])
					continue;
				for (int i = 0; i <= 1; i++) {
					if (testProjection(point, e[i], elem, adjEdge[0]) && testProjection(point, e[i], elem, adjEdge[1])) {
						println(e, point, " rejected by ", elem, adjEdge, " *** ");
						return true;
					}
				}
			}
		}

		return false;
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Check if \p e any adjacent face.
	/// \return bool
	///
	/////////////////////////////////////////////////
	bool isEdgeOverlap(const ourEdge &e, set<ourFace> excludeSet = set<ourFace>()) {
		ourEdge edge(0, 1);
		for (int i = 0; i < 2; ++i) {
			for (auto &elem : vertAdjEdges[e[i]]) {
				edge = *elem;
				if (edge[0] == e[!i])
					continue;
				if (edge[1] == e[!i])
					continue;
				if (excludeSet.find(ourFace(edge, e[i])) != excludeSet.end())
					continue;
				if (testProjection(edge[0], edge[1], e[i], e[!i])) {
					println(e, "rejected by", ourFace(edge[0], edge[1], e[i]));
					return true;
				}
			}
		}
		return false;
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Check if the face formed by \p edge and \p newPoint can be added to Surface, with or without removing an existing face.
	/// \param[in] edge ourEdge
	/// \param[in] newPoint int
	/// \param[in,out] toBeRemoved int&
	/// \return bool
	///
	/////////////////////////////////////////////////
	bool validToAdd(ourEdge edge, int newPoint, int &toBeRemoved) {
		auto &edgeDeg = edgeDegree[edge];
		println(edge, newPoint, edgeDeg.size());
		if (edgeDeg.size() == 2) {
			// println(edge, "[", edgeDeg[0], edgeDeg[1], "]", newPoint);

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
				if (
					// (
					// (std::abs(a1 - a2) < 0.017453292519943295) && //0.005 * (a1 + a2)) && //Error is twice of the value used
					lockedLevel(ourFace(edge[0], edge[1], edgeDeg[1])) < 3 &&
					getTriScore(edge[0], edge[1], newPoint) > getTriScore(edge[0], edge[1], edgeDeg[1]) &&
					((std::abs(a1 - a2) >= 0.017453292519943295) && (a1 > a2))) {
					toBeRemoved = edgeDeg[1];
					// println("Angle", a1, "->", a2);
					// println(edge, newPoint, val(toBeRemoved), getTriScore(edge[0], edge[1], newPoint), getTriScore(edge[0], edge[1], edgeDeg[1]));
					// println("Edge Degree=2, returning true");
					return true;
				}
			}
			else {
				if (
					// (
					// (std::abs(a1 - a3) < 0.017453292519943295) && // 0.005 * (a1 + a3)) && //Error is twice of the value used
					lockedLevel(ourFace(edge[0], edge[1], edgeDeg[0])) < 3 &&
					getTriScore(edge[0], edge[1], newPoint) > getTriScore(edge[0], edge[1], edgeDeg[0]) &&
					((std::abs(a1 - a3) >= 0.017453292519943295) && (a1 > a3))) {
					toBeRemoved = edgeDeg[0];
					// println("Angle", a1, "->", a3);
					// println(edge, newPoint, val(toBeRemoved), getTriScore(edge[0], edge[1], newPoint), getTriScore(edge[0], edge[1], edgeDeg[0]));
					// println("Edge Degree=2, returning true");
					return true;
				}
			}
			// println("Angle", 180 * a1 / M_PI, "->", 180 * a2 / M_PI, ",", 180 * a3 / M_PI);
			// println("Edge Degree=2, returning false");
			return false;
		}
		if (edgeDeg.size() == 0) {
			// println("Edge Degree=0, returning true");
			return true;
		}

		Vector3D norm1 = CGAL::normal(pts[edge[0]], pts[edge[1]], pts[newPoint]),
				 norm2 = CGAL::normal(pts[edge[1]], pts[edge[0]], pts[edgeDeg[0]]);
		norm1 = norm1 / CGAL::sqrt(norm1.squared_length());
		norm2 = norm2 / CGAL::sqrt(norm2.squared_length());
		// println("Edge Degree=1, returning ", (norm1 * norm2));
		return true;
		// return norm1 * norm2 >= 0.707106781187;
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Generate Stitch Edges
	/// \param leftVerts set<int>&
	/// \param tempAdjList map<int, set<int>>&
	/// \return set<ourEdge>
	///
	/////////////////////////////////////////////////
	set<ourEdge> fillHoles(set<int> &leftVerts, map<int, set<int>> &tempAdjList) {
		set<ourEdge> newEdges;

		if (leftVerts.size() == 3) {
			int p[3], i = 0;
			for (auto &elem : leftVerts) {
				p[i++] = elem;
			}
			int remove;
			if (validToAdd(ourEdge(p[1], p[2]), p[0], remove) &&
				validToAdd(ourEdge(p[0], p[2]), p[1], remove) &&
				validToAdd(ourEdge(p[0], p[1]), p[2], remove))
				addFaceToModel({p[0], p[1], p[2]});
			return newEdges;
		}

		vector<bool> covered(pts.size(), false);

		if (possibleEdges.size() < std::log2(possibleEdges.size()) * leftVerts.size() * leftVerts.size()) {
			for (auto &elem : possibleEdges) {
				ourEdge edge = elem.second;
				int u = edge[0], v = edge[1];
				if (contains(leftVerts, u) &&
					contains(leftVerts, v) &&
					!isEdgeOverlap(edge) &&
					(!covered[u] || !covered[v]) &&
					elem.first <= 4 * maxLength[u] &&
					elem.first <= 4 * maxLength[v]) {
					newEdges.insert(elem.second);
					// print("added");
					covered[u] = true;
					covered[v] = true;
					// println(edge);
				}
			}
		}
		else {
			vector<pair<Kernel::FT, ourEdge>> edges;
			for (auto &u : leftVerts) {
				for (auto &v : leftVerts) {
					if (u >= v)
						break;
					auto itr = possibleEdges.find(make_pair(sqDist(pts[u], pts[v]), ourEdge(u, v)));
					if (itr == possibleEdges.end())
						continue;
					edges.push_back(*itr);
				}
			}
			std::sort(edges.begin(), edges.end());
			for (auto &elem : edges) {
				ourEdge edge = elem.second;
				int u = edge[0], v = edge[1];
				if (!isEdgeOverlap(edge) &&
					(!covered[u] || !covered[v]) &&
					elem.first <= 4 * maxLength[u] &&
					elem.first <= 4 * maxLength[v]) {
					newEdges.insert(edge);
					// print("added");
					covered[u] = true;
					covered[v] = true;
					// println(edge);
				}
			}
		}
		return newEdges;
	}

	void addAdjacency(int u, int v) {
		currAdjList[u].insert(v);
		currAdjList[v].insert(u);
	}

	void removeAdjacency(int u, int v) {
		if (std::binary_search(initialEdges.begin(), initialEdges.end(), ourEdge(u, v)))
			return;
		currAdjList[u].erase(v);
		currAdjList[v].erase(u);
	}

	/////////////////////////////////////////////////
	/// \brief
	/// Main Reconstruction Function
	/// \return void
	///
	/////////////////////////////////////////////////
	void reconstruct() {

		currAdjList = getAdjList(pts.size(), initialEdges);
		set<ourEdge> origEdges(initialEdges.begin(), initialEdges.end());

		maxEdge = 0;

		println(val(initialEdges.size()));

		for (auto &elem : initialEdges) {
			maxEdge = std::max(maxEdge, sqDist(pts[elem[0]], pts[elem[1]]));
		}
		cerr << maxEdge << "\t" << CGAL::sqrt(maxEdge) << "\n";

		maxEdge *= Kernel::FT(1.0004);

		set<std::tuple<Kernel::FT, ourEdge, int>,
			std::greater<std::tuple<Kernel::FT, ourEdge, int>>>
			pq, nextPQ;

		int edgeCondition = 2;
		set<ourEdge> nextEdges;
		bool changed = true;
		int tempModelCount = 1;

		while (changed) {
			changed = false;

			map<int, bool> leftVerts;
			map<int, set<int>> tempAdjList;

			for (auto &elem : modelEdges) {

				if (edgeDegree[elem].size() == 1) {
					leftVerts.insert({elem[0], false});
					leftVerts.insert({elem[1], false});

					tempAdjList[elem[0]].insert(elem[1]);
					tempAdjList[elem[1]].insert(elem[0]);
				}
			}
			for (auto &elem : origEdges) {

				if (edgeDegree[elem].size() <= 1) {
					leftVerts.insert({elem[0], false});
					leftVerts.insert({elem[1], false});
					tempAdjList[elem[0]].insert(elem[1]);
					tempAdjList[elem[1]].insert(elem[0]);
				}
			}

			println(val(leftVerts.size()));

			print("Size Change", possibleEdges.size());
			if (leftVerts.size()) {
				for (auto itr = possibleEdges.begin(); itr != possibleEdges.end();) {
					int u = (*itr).second[0], v = (*itr).second[1];
					if (leftVerts.find(u) == leftVerts.end() ||
						leftVerts.find(v) == leftVerts.end())
						itr = possibleEdges.erase(itr);
					else
						++itr;
				}
			}
			println(possibleEdges.size());

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

#ifdef ENABLE_STEP_MODEL
			println("Generated tempMod" + std::to_string(tempModelCount) + ".off");
			writeModel("tempMod" + std::to_string(tempModelCount++) + ".off");
#endif

			println("Iteration #=", 2 - edgeCondition + 1);
			for (auto &edge : nextEdges) {
				// bool inserted = false;
				for (auto &elem : getFacesFromEdge(edge)) {
					pq.insert(std::make_tuple(elem.second, edge, elem.first));
					// inserted = true;
					println("Initially added by nextEdges", edge, elem.first);
				}
			}
			for (auto &edge : origEdges) {
				// bool inserted = false;
				for (auto &elem : getFacesFromEdge(edge)) {
					pq.insert(std::make_tuple(elem.second, edge, elem.first));
					// inserted = true;
					println("Initially added by origEdges", edge, elem.first);
				}
			}

			cerr << "Iteration #" << 2 - edgeCondition + 1 << std::endl;
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

				println("###################################################");
				println(score, edge, point);

				ourFace face(edge[0], edge[1], point);

				if (score <= -0.8660254037844387)
					continue;
				if (modelFaces.find(face) != modelFaces.end())
					continue;
				// if (currAdjList[edge[0]].find(edge[1]) == currAdjList[edge[0]].end() &&
				// 	currAdjList[edge[1]].find(edge[0]) == currAdjList[edge[1]].end())
				// continue;
				// if (edgeDegree[edge].size() == 2) {
				// 	println("edgeDegree[edge].size() == 2");
				// 	continue;
				// }
				// bool passed = true;
				// for (int i = 0; passed && i < 3; ++i) {
				// 	for (auto &point : currAdjList[face[i]]) {
				// 		if (faceEdgeOverlap(face[(i + 1) % 3], face[(i + 2) % 3], face[i],
				// 							point)) {
				// 			println(face, "rejected by", ourEdge(face[i], point));
				// 			passed = false;
				// 			break;
				// 		}
				// 	}
				// }
				// if (!passed)
				// 	continue;
				ourEdge newEdge1(edge[0], point), newEdge2(edge[1], point);
				int remove1 = -1, remove2 = -1, remove3 = -1;
				//cout << get<0>(*pq.begin()) << " " << edge << " " << point << "::";

				println(val(isFaceDelaunay(face)));
				println(val(validToAdd(edge, point, remove1)));
				println(val(validToAdd(newEdge1, edge[1], remove2)));
				println(val(validToAdd(newEdge2, edge[0], remove3)));
				if (isFaceDelaunay(face) &&
					validToAdd(edge, point, remove1) &&
					validToAdd(newEdge1, edge[1], remove2) &&
					validToAdd(newEdge2, edge[0], remove3)) {

					set<ourFace> excludeSet;
					set<ourEdge> excludeEdge;
					if (remove1 != -1) excludeSet.insert(ourFace(edge, remove1));
					if (remove2 != -1) excludeSet.insert(ourFace(newEdge1, remove2));
					if (remove3 != -1) excludeSet.insert(ourFace(newEdge2, remove3));

					for (auto &face : excludeSet) {
						if (!contains(initialEdges, ourEdge(face[0], face[1])))
							excludeEdge.insert(ourEdge(face[0], face[1]));
						if (!contains(initialEdges, ourEdge(face[0], face[2])))
							excludeEdge.insert(ourEdge(face[0], face[2]));
						if (!contains(initialEdges, ourEdge(face[1], face[2])))
							excludeEdge.insert(ourEdge(face[1], face[2]));
					}

					print("excludeSet:");
					for (auto &elem : excludeSet) {
						print(elem);
					}
					println();

					if (isFaceOverlap(edge, point, excludeSet, excludeEdge)) {
						println("isFaceOverlap1", edge, point);
						continue;
					}
					if (isFaceOverlap(newEdge1, edge[1], excludeSet, excludeEdge)) {
						println("isFaceOverlap2", newEdge1, edge[1]);
						continue;
					}
					if (isFaceOverlap(newEdge2, edge[0], excludeSet, excludeEdge)) {
						println("isFaceOverlap3", newEdge2, edge[0]);
						continue;
					}

					// (remove != -1 || !isFaceOverlap(edge, point, currAdjList)) &&
					// (remove1 != -1 || !isFaceOverlap(newEdge1, edge[1], currAdjList)) &&
					// (remove2 != -1 || !isFaceOverlap(newEdge2, edge[0], currAdjList)))

					if (testEdges(pts[face[0]], pts[face[1]], pts[face[2]], maxEdge) < edgeCondition) {
						nextPQ.insert(std::make_tuple(score, edge, point));
						println("testEdges Failed");
						println(sqDist(pts[face[0]], pts[face[1]]), sqDist(pts[face[0]], pts[face[2]]), sqDist(pts[face[1]], pts[face[2]]), maxEdge);

						continue;
					}

					if (remove1 != -1) {
						println("Remove", ourFace(edge, remove1), "due to", face);
						removeFace(edge, remove1);
						removeAdjacency(remove1, edge[0]);
						removeAdjacency(remove1, edge[1]);
					}
					if (remove2 != -1) {
						println("Remove", ourFace(newEdge1, remove2), "due to", face);
						removeFace(newEdge1, remove2);
						removeAdjacency(remove2, newEdge1[0]);
						removeAdjacency(remove2, newEdge1[1]);
					}
					if (remove3 != -1) {
						println("Remove", ourFace(newEdge2, remove3), "due to", face);
						removeFace(newEdge2, remove3);
						removeAdjacency(remove3, newEdge2[0]);
						removeAdjacency(remove3, newEdge2[1]);
					}

					changed = true;
					addFaceToModel(face);
					origEdges.erase(edge);
					origEdges.erase(newEdge1);
					origEdges.erase(newEdge2);

					for (int i = 0; i < 3; i++) {
						currAdjList[face[i]].insert(face[(i + 1) % 3]);
						currAdjList[face[i]].insert(face[(i + 2) % 3]);
					}

					for (int i = 0; i < 3; i++) {
						ourEdge e(face[i], face[(i + 1) % 3]);
						for (auto &elem : getFacesFromEdge(e)) {
							pq.insert(std::make_tuple(elem.second, e, elem.first));
							println(edge, point, "added", e, elem.first);
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
			edgeCondition--;

#ifdef ENABLE_STEP_MODEL
			println("Generated tempMod" + std::to_string(tempModelCount) + ".off");
			writeModel("tempMod" + std::to_string(tempModelCount++) + ".off");
#endif

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

	/////////////////////////////////////////////////
	/// \brief
	/// Perform Post Processing on th Surface
	/// \return void
	///
	/////////////////////////////////////////////////
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

		for (auto elem : leftVerts) {
			print(elem);
		}
		println();

		vector<vector<int>> cycles = getAllCycles(tempAdjList, leftVerts);

		for (auto &cycle : cycles) {
			print("Cycle of Size:", cycle.size(), ":");
#ifdef ENABLE_LOG
			for (auto &p : cycle) {
				print(p);
			}
#endif
			println();

			if (cycle.size() == 3) {
				addFaceToModel({cycle[0], cycle[1], cycle[2]});
				continue;
			}
			else {
				std::map<int, Kernel::FT> vals;
				std::map<int, std::list<int>::iterator> iters;
				std::list<int> l;
				set<ourFace> excludeSet;
				set<pair<Kernel::FT, int>, std::greater<pair<Kernel::FT, int>>> PQ;
				for (auto &p : cycle) {
					l.push_front(p);
					iters[p] = l.begin();
				}
				for (auto curr = l.begin(); curr != l.end(); ++curr) {
					auto prev = prevIterator(l, curr);
					auto next = nextIterator(l, curr);
					vals[*curr] = getAngle(*prev, *curr, *next);
					if (isEdgeOverlap(ourEdge(*prev, *next), excludeSet)) {
						vals[*curr] += -2;
					}
					PQ.insert(make_pair(vals[*curr], *curr));
				}
				while (!PQ.empty()) {
					auto score = PQ.begin()->first;
					auto curr = iters[PQ.begin()->second];
					PQ.erase(PQ.begin());
					println(__LINE__, score, *curr);
					auto prev = prevIterator(l, curr);
					auto next = nextIterator(l, curr);
					addFaceToModel(ourFace(*prev, *curr, *next));
					excludeSet.insert(ourFace(*prev, *curr, *next));
					if (l.size() == 3)
						break;
					auto prevprev = prevIterator(l, prev);
					auto nextnext = nextIterator(l, next);
					PQ.erase(make_pair(vals[*prev], *prev));
					PQ.erase(make_pair(vals[*next], *next));
					l.erase(curr);
					vals[*prev] = getAngle(*prevprev, *prev, *next);
					if (isEdgeOverlap(ourEdge(*prevprev, *next), excludeSet)) {
						vals[*prev] += -2;
					}
					// if (contains(modelFaces, ourFace(*prevprev, *prev, *next)) /*||
					// 	isFaceOverlap(ourEdge(*prevprev, *prev), *next) ||
					// 	isFaceOverlap(ourEdge(*prevprev, *next), *prev) ||
					// 	isFaceOverlap(ourEdge(*next, *curr), *prevprev)*/) {
					// 	vals[*prev] = -2;
					// }
					vals[*next] = getAngle(*prev, *next, *nextnext);
					if (isEdgeOverlap(ourEdge(*prev, *nextnext), excludeSet)) {
						vals[*next] += -2;
					}
					PQ.insert(make_pair(vals[*prev], *prev));
					PQ.insert(make_pair(vals[*next], *next));
				}
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
		println("Overhanging Faces");
		for (auto &face : modelFaces) {
			if (lockedLevel(face) < 3) {
				// println(face, lockedLevel(face));
				modelFaces.erase(face);
			}
		}
	}

  public:
	/////////////////////////////////////////////////
	/// \brief
	/// Opens the file \p inputFilePath and generate the Surface Reconstruction
	/// \param inputFilePath const char *
	/// \return
	///
	/////////////////////////////////////////////////
	SurfaceReconstruct(const char *inputFilePath) {
		std::ifstream inputFile(inputFilePath);

		auto start = std::chrono::high_resolution_clock::now();

		if (!CGAL::read_xyz_points(
				inputFile, back_inserter(pts))) { // output iterator over points
			cerr << "Error: cannot read file.";
			return;
		}

		sortAndRemoveDuplicate(pts);
		pts.shrink_to_fit();

		auto finish = std::chrono::high_resolution_clock::now();
		cout << pts.size() << " points read in "
			 << std::chrono::duration<double>(finish - start).count() << " secs"
			 << std::endl;

		// ###--> Create Delaunay triangulation and check Validity

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

		// ###--> DONE

		maxLength.resize(pts.size(), 0);
		vertAdjEdges.resize(pts.size());

		// ###--> Populate possibleFaces
		for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
			Triangle3D tri = dt.triangle(*faceItr);
			possibleFaces.insert(ourFace(getIndex(pts, tri[0]),
										 getIndex(pts, tri[1]),
										 getIndex(pts, tri[2])));
		}

		// ###--> Populate possibleEdges
		for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
			Segment3D segment = dt.segment(*edgeItr);
			ourEdge edge(getIndex(pts, segment[0]), getIndex(pts, segment[1]));
			possibleEdges.insert(make_pair(segment.squared_length(), edge));
		}

		// ###--> Create MST
		start = std::chrono::high_resolution_clock::now();
		initialEdges = getMstEdges();
		finish = std::chrono::high_resolution_clock::now();
		cout << "MST created in "
			 << std::chrono::duration<double>(finish - start).count() << " secs"
			 << std::endl;
		// ###--> DONE

		// ###--> Populate maxLength from initialEdges
		for (auto &edge : initialEdges) {
			maxLength[edge[0]] =
				std::max(maxLength[edge[0]], sqDist(pts[edge[0]], pts[edge[1]]));
			maxLength[edge[1]] =
				std::max(maxLength[edge[1]], sqDist(pts[edge[0]], pts[edge[1]]));
		}

		start = std::chrono::high_resolution_clock::now();

		modelFaces.clear();
		modelEdges.clear();

		sortAndRemoveDuplicate(initialEdges);

		// ###--> reconstruct surface
		reconstruct();

		// ###--> fill Holes if present
		postProcess();

		modelEdges.clear();

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

	/////////////////////////////////////////////////
	/// \brief
	/// Return if the Surface Reconstruction is valid or not. Reconstruction is Invalid only if a Delaunay can;t be created from the input Point Set.
	/// \return bool
	///
	/////////////////////////////////////////////////
	bool isValid() { return valid; }

	/////////////////////////////////////////////////
	/// \brief
	/// Write the Surface Constructed in \p outputFilePath in OFF File Format
	/// \param outputFilePath std::string
	/// \return void
	///
	/////////////////////////////////////////////////
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
