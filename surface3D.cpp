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

vector<ourEdge> get_All_Edges(const DT3 &dt, const vector<Point3D> &points) {
	vector<ourEdge> allEdges;
	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
		Segment3D seg = dt.segment(*edgeItr);
		allEdges.push_back(ourEdge(getIndex(points, seg[0]),
								   getIndex(points, seg[1])));
	}

	return allEdges;
}

/**********************************************************************************************/ /**
 * @fn	vector<ourEdge> getMstEdges(const vector<Point3D> &points, const DT3 &dt)
 *
 * @brief	Gets MST edges using Kruskal's Algorithm.
 *
 * @param	points	The input points.
 * @param	dt	  	The Delaunay triangulation of the points.
 *
 * @return	The MST edges.
 **************************************************************************************************/

vector<ourEdge> getMstEdges(const vector<Point3D> &pts, const vector<set<int>> &dtAdjList) {
	vector<CGAL::Union_find<int>::handle> handle;
	handle.reserve(pts.size());

	CGAL::Union_find<int> uf;

	for (int i = 0; i < pts.size(); i++)
		handle.push_back(uf.make_set(i));

	vector<std::tuple<double, int, int>> allEdges;

	for (int u = 0; u < pts.size(); ++u) {
		for (int v : dtAdjList[u])
			if (u < v) {
				allEdges.push_back(std::make_tuple(CGAL::squared_distance(pts[u], pts[v]), u, v));
			}
	}

	sort(allEdges.begin(), allEdges.end());

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

vector<set<int>> getAdjList(int pointsCount, const vector<ourEdge> &edges) {
	vector<set<int>> adjList(pointsCount);
	for (auto edge : edges) {
		adjList[edge[0]].insert(edge[1]);
		adjList[edge[1]].insert(edge[0]);
	}

	return adjList;
}

int testEdges(Point3D a, Point3D b, Point3D c, double length) {
	int ans = 0;
	if (CGAL::squared_distance(a, b) <= length)
		ans++;
	if (CGAL::squared_distance(a, c) <= length)
		ans++;
	if (CGAL::squared_distance(c, b) <= length)
		ans++;
	return ans;
}

class SurfaceReconstruct {
	/** @brief	The input Points. */
	vector<Point3D> pts;

	/** @brief	Delaunay triangulation. */
	DT3 dt;

	/** @brief	Handles for referencing points. */
	vector<DT3::Vertex_handle> ptsHandle;

	/** @brief	File for logs. */
	std::ofstream logFile;

	/** @brief	Validity of Model. */
	bool valid;

	/** @brief	Faces in Delaunay. */
	vector<ourFace> dtFaces;

	/** @brief	Edges of the Final Model. */
	set<ourEdge> modelEdges;

	/** @brief	Faces of the Final Model.  */
	set<ourFace> modelFaces;

	/** @brief	Adjacency List for Delaunay. */
	vector<set<int>> dtAdjList;

	/** @brief	Face Degree of Edges. */
	map<ourEdge, vector<int>> edgeDegree;

	/** @brief	The length of longest edge in MST. */
	double maxEdge;

	/**********************************************************************************************/ /**
	 * @fn	bool isFaceDelaunay(ourFace face)
	 *
	 * @brief	Query if 'face' is part of delaunay.
	 *
	 * @param	face	The face.
	 *
	 * @return	True if face is delaunay, false if not.
	 **************************************************************************************************/

	bool isFaceDelaunay(const ourFace &face) {

		// return dtAdjList[face[0]].find(face[1]) != dtAdjList[face[0]].end() &&
		// 	   dtAdjList[face[0]].find(face[2]) != dtAdjList[face[0]].end() &&
		// 	   dtAdjList[face[1]].find(face[2]) != dtAdjList[face[1]].end();
		//~ DT3::Cell_handle cell;
		//~ int u, v, w;
		//~ return dt.is_facet(ptsHandle[face[0]], ptsHandle[face[1]], ptsHandle[face[2]], cell, u, v, w);
		return std::binary_search(dtFaces.cbegin(), dtFaces.cend(), face);
		//return
	}

	/**********************************************************************************************/ /**
	 * @fn	bool validToAdd(const ourEdge &edge, int newPoint)
	 *
	 * @brief	Valid to add.
	 *
	 * @param	edge		The edge.
	 * @param	newPoint	The new point.
	 *
	 * @return	True if it succeeds, false if it fails.
	 **************************************************************************************************/

	bool validToAdd(const ourEdge &edge, int newPoint) {
		auto &edgeDeg = edgeDegree[edge];
		if (edgeDeg.size() == 2)
			return false;
		if (edgeDeg.size() == 0)
			return true;

		Vector3D norm1 = CGAL::cross_product(pts[edge[1]] - pts[edge[0]], pts[newPoint] - pts[edge[1]]),
				 norm2 = CGAL::cross_product(pts[edge[0]] - pts[edge[1]], pts[edgeDeg.back()] - pts[edge[0]]);
		return CGAL::angle(norm1, norm2) >= CGAL::RIGHT;

		/*Point3D A = points[edge[0]],
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
			fprintln(logFile, "Angle between the face (", edge[0], edge[1], edgeDegree.back(), ") and (", edge[0], edge[1], newPoint, ") is", angle);*/
		//}

		//if (!ans && debug) {
		//	fprint(logFile, "Problem Adding ", newPoint, "to", edge, ". edgedegree={");
		//	for (int p : edgeDegree)
		//		fprint(logFile, p);
		//	fprintln(logFile, "}");
		//}
	}

	/**********************************************************************************************/ /**
	 * @fn	vector<int> bfs(int start, map<int, bool> &verts, const vector<set<int>> &currAdjList)
	 *
	 * @brief	Find Connected Vertices present in verts using BFS.
	 *
	 * @param 		  	start	   	start.
	 * @param [in,out]	verts	   	vertices.
	 * @param 		  	currAdjList	Adjacency List.
	 *
	 * @return	A vector&lt;int&gt;
	 **************************************************************************************************/

	set<int> bfs(int start, map<int, bool> &verts, const vector<set<int>> &currAdjList) {
		set<int> connected;
		verts[start] = true;
		std::queue<int> Q;
		Q.push(start);
		while (!Q.empty()) {
			int u = Q.front();
			Q.pop();
			connected.insert(u);
			for (auto &v : currAdjList[u]) {
				cout<<u<<"->"<<v<<"\n";
				if (verts.find(v) != verts.end() && !verts[v]) {
					verts[v] = true;
					Q.push(v);
				}
			}
		}
		return connected;
	}

	/**********************************************************************************************/ /**
	 * @fn	double getScore(int u, int v, int w)
	 *
	 * @brief	Gets the score of the triangle formed by points u,v,w.
	 *
	 * @param	u	Point 1.
	 * @param	v	Point 2.
	 * @param	w	Point 3.
	 *
	 * @return	The score.
	 **************************************************************************************************/

	double getScore(int u, int v, int w) {
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

	/**********************************************************************************************/ /**
	 * @fn	vector<pair<double, int>> getFacesFromEdge(const ourEdge &edge, const vector<set<int>> &currAdjList)
	 *
	 * @brief	Gets Possible faces from edges.
	 *
	 * @param	edge	   	The edge.
	 * @param	currAdjList	current Adjacency List.
	 *
	 * @return	The faces from edges.
	 **************************************************************************************************/

	vector<pair<int, double>> getFacesFromEdge(const ourEdge &edge, const vector<set<int>> &currAdjList, int edgeCondition) {
		vector<pair<int, double>> elems;
		for (auto &point : currAdjList[edge[0]]) {
			if (edge[0] == 22 && edge[1] == 37 && point == 24)
				cout << (point != edge[1]) << " " << isFaceDelaunay(ourFace(point, edge[0], edge[1])) << " " << (modelFaces.find(ourFace(point, edge[0], edge[1])) == modelFaces.end()) << " " << (testEdges(pts[edge[0]], pts[edge[1]], pts[point], maxEdge) >= edgeCondition) << "\n";
			if (point != edge[1] &&
				isFaceDelaunay(ourFace(point, edge[0], edge[1])) &&
				modelFaces.find(ourFace(point, edge[0], edge[1])) == modelFaces.end() &&
				testEdges(pts[edge[0]], pts[edge[1]], pts[point], maxEdge) >= edgeCondition) {
				elems.push_back(std::make_pair(point, getScore(point, edge[0], edge[1])));
			}
		}
		for (auto &point : currAdjList[edge[1]]) {
			if (point != edge[0] &&
				isFaceDelaunay(ourFace(point, edge[0], edge[1])) &&
				modelFaces.find(ourFace(point, edge[0], edge[1])) == modelFaces.end() &&
				testEdges(pts[edge[0]], pts[edge[1]], pts[point], maxEdge) >= edgeCondition) {
				elems.push_back(std::make_pair(point, getScore(point, edge[0], edge[1])));
			}
		}
		sortAndRemoveDuplicate(elems);
		return elems;
	}

	set<ourEdge> processHole(const set<int> &hole, const vector<set<int>> &currAdjList) {
		set<ourEdge> newEdges;

		if (hole.size() <= 3) {
			ourEdge e(0, 1);
			auto itr = hole.begin();
			e[0] = *itr++;
			e[1] = *itr;
			newEdges.insert(e);
			return newEdges;
		}

		vector<ourEdge> allEdges;
		for (auto e : get_All_Edges(dt, pts)) {
			if (hole.find(e[0]) != hole.end() &&
				hole.find(e[1]) != hole.end() &&
				modelEdges.find(e) == modelEdges.end())
				allEdges.push_back(e);
		}

		vector<set<int>> adjList = getAdjList(pts.size(), allEdges);
		for (int u : hole) {
			int v = SIZE_MAX;
			double dist = inf;
			for (int tempV : adjList[u]) {
				if (dist > CGAL::squared_distance(pts[u], pts[tempV])) {
					v = tempV;
					dist = CGAL::squared_distance(pts[u], pts[tempV]);
				}
			}

			if (v != SIZE_MAX && newEdges.find(ourEdge(u, v)) == newEdges.end()) {
				newEdges.insert(ourEdge(u, v));
				// fprintln(logFile, ourEdge(u, v));
			}
		}

		/*for (auto &u : hole) {
			int minV = SIZE_MAX;
			double dist = inf;
			for (auto &v : hole) {
				if (u >= v)
					continue;
				if (dist > CGAL::squared_distance(pts[u], pts[v]) &&
					newEdges.find(ourEdge(u, v)) == newEdges.end() &&
					(hole.size() <= 3 || currAdjList[u].find(v) == currAdjList[u].end())) {
					minV = v;
					dist = CGAL::squared_distance(pts[u], pts[v]);
				}
			}
			if (minV != SIZE_MAX) {
				modelEdges.insert(ourEdge(u, minV));
				newEdges.insert(ourEdge(u, minV));
			}
		}*/
		return newEdges;
	}

	/**********************************************************************************************/ /**
	 * @fn	void reconstruct(vector<ourEdge> initialEdges)
	 *
	 * @brief	Construct a Surface Reconstruction using the input Edges as a base.
	 *
	 * @param	initialEdges	The initial edges.
	 **************************************************************************************************/

	void reconstruct(vector<ourEdge> &initialEdges) {
		vector<set<int>> currAdjList = getAdjList(pts.size(), initialEdges);
		maxEdge = 0;
		for (auto &e : initialEdges)
			maxEdge = std::max(maxEdge, CGAL::squared_distance(pts[e[0]], pts[e[1]]));

		set<std::tuple<double, ourEdge, int>, std::greater<std::tuple<double, ourEdge, int>>> pq;

		//std::priority_queue<std::tuple<double, ourEdge, int>> pq;

		std::sort(initialEdges.begin(), initialEdges.end());

		// cout<<"MST Edges\n";
		// for (auto &elem : initialEdges) {
		// 	cout<<elem<<"\n";
		// }

		modelEdges.insert(initialEdges.begin(), initialEdges.end());

		int edgeCondition = 2;

		for (auto &edge : initialEdges) {
			for (auto &elem : getFacesFromEdge(edge, currAdjList, edgeCondition)) {
				pq.insert(std::make_tuple(elem.second, edge, elem.first));
				//~ cout<<"Initially added "<<edge<<" and "<<elem.first<<"\n";
			}
		}

		while (edgeCondition >= 2) {
			while (!pq.empty()) {
				//~ for (auto &elem : pq) {
				//~ cout<<"["<<get<0>(elem)<<" "<<get<1>(elem)<<" "<<get<2>(elem)<<"],\n";
				//~ }
				//~ cout<<"\n";
				//~ for (int i = 0; i<pts.size(); ++i) {
				//~ cout<<i<<": [";
				//~ for (auto elem : currAdjList[i]) {
				//~ cout<<elem<<" ";
				//~ }
				//~ cout<<"]\n";
				//~ }
				//~ for (auto &elem : edgeDegree) {
				//~ if (elem.second.size() >= 1) {
				//~ cout<<elem.first<<": [";
				//~ for (auto &p : elem.second)
				//~ cout<<p<<" ";
				//~ cout<<"]\n";
				//~ }
				//~ }
				ourEdge edge = get<1>(*pq.begin());
				int point = get<2>(*pq.begin());
				cout << get<0>(*pq.begin()) << " " << edge << " " << point << "\n";
				pq.erase(pq.begin());
				ourFace face(edge[0], edge[1], point);
				if (edgeDegree[edge].size() == 2)
					continue;
				ourEdge newEdge1(edge[0], point), newEdge2(edge[1], point);
				// cout<<(modelFaces.find(face) == modelFaces.end())<<" "<<(isFaceDelaunay(face))<<" "<<(validToAdd(edge, point))<<" "<<(validToAdd(newEdge1, edge[1]))<<" "<<(validToAdd(newEdge2, edge[0]))<<" "<<(testEdges(pts[face[0]], pts[face[1]], pts[face[2]], maxEdge) >= edgeCondition)<<"\n";
				if (modelFaces.find(face) == modelFaces.end() &&
					isFaceDelaunay(face) &&
					validToAdd(edge, point) &&
					validToAdd(newEdge1, edge[1]) &&
					validToAdd(newEdge2, edge[0]) &&
					testEdges(pts[face[0]], pts[face[1]], pts[face[2]], maxEdge) >= edgeCondition) {
					//~ cout<<"Added\n";
					modelEdges.insert(newEdge1);
					modelEdges.insert(newEdge2);

					edgeDegree[edge].push_back(point);
					edgeDegree[newEdge1].push_back(edge[1]);
					edgeDegree[newEdge2].push_back(edge[0]);

					modelFaces.insert(face);

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
							//~ cout<<edge<<" added "<<ourEdge(u, v)<<" and "<<elem.first<<"\n";
						}
					}
				}
			}

			edgeCondition--;

			map<int, bool> leftVerts;
			set<ourEdge> nextEdges;
			vector<set<int>> tempAdjList(pts.size());

			for (auto &elem : modelEdges) {
				if (edgeDegree[elem].size() == 1) {
					leftVerts.insert({elem[0], false});
					leftVerts.insert({elem[1], false});
					tempAdjList[elem[0]].insert(elem[1]);
					tempAdjList[elem[1]].insert(elem[0]);
					// nextEdges.insert(elem);
				}
			}

			for (auto &elem : leftVerts) {
				if (!elem.second) {
					auto hole = bfs(elem.first, leftVerts, tempAdjList);
					auto Edges = processHole(hole, tempAdjList);
					nextEdges.insert(Edges.begin(), Edges.end());
				}
			}
			// modelEdges.clear();
			for (auto &edge : nextEdges) {
				for (auto &elem : getFacesFromEdge(edge, currAdjList, edgeCondition)) {
					pq.insert(std::make_tuple(elem.second, edge, elem.first));
					// cout<<"Finally added "<<edge<<" and "<<elem.first<<"\n";
				}
				// modelEdges.insert(edge);
			}
		}
		modelEdges.clear();
		while (!pq.empty()) {
			ourEdge edge = get<1>(*pq.begin());
			pq.erase(pq.begin());
			modelEdges.insert(edge);
		}
	}

  public:
	/**********************************************************************************************/ /**
	 * @fn	SurfaceReconstruct(const char *inputFilePath)
	 *
	 * @brief	Constructor.
	 *
	 * @param	inputFilePath	Full pathname of the input file.
	 **************************************************************************************************/

	SurfaceReconstruct(const char *inputFilePath) {

		std::ifstream inputFile(inputFilePath);
		logFile.open("log.txt");
		logFile << std::boolalpha;

		auto start = std::chrono::high_resolution_clock::now();
		if (!CGAL::read_xyz_points(inputFile, back_inserter(pts))) { // output iterator over points
			cerr << "Error: cannot read file.";
			return;
		}
		sortAndRemoveDuplicate(pts);
		pts.shrink_to_fit();

		auto finish = std::chrono::high_resolution_clock::now();
		cout << pts.size() << " points read in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

		start = std::chrono::high_resolution_clock::now();
		dt.insert(pts.begin(), pts.end());
		if (!dt.is_valid(true)) {
			cerr << "Error: fail to build a Delaunay triangulation.\n";
			valid = false;
			return;
		}
		if (dt.dimension() != 3) {
			cerr << "Error: cannot built a 3D triangulation.\n Current dimension = " << dt.dimension() << "\n";
			valid = false;
			return;
		}
		finish = std::chrono::high_resolution_clock::now();
		cout << "Delaunay Triangulation created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

		valid = true;

		ptsHandle.resize(pts.size());
		for (auto itr = dt.finite_vertices_begin(); itr != dt.finite_vertices_end(); itr++) {
			ptsHandle[getIndex(pts, itr->point())] = itr;
		}

		for (auto faceItr = dt.finite_facets_begin(); faceItr != dt.finite_facets_end(); faceItr++) {
			Triangle3D tri = dt.triangle(*faceItr);
			dtFaces.push_back(ourFace(getIndex(pts, tri[0]),
									  getIndex(pts, tri[1]),
									  getIndex(pts, tri[2])));
		}

		std::sort(dtFaces.begin(), dtFaces.end());

		dtAdjList = getAdjList(pts.size(), get_All_Edges(dt, pts));

		start = std::chrono::high_resolution_clock::now();
		vector<ourEdge> initialEdges = getMstEdges(pts, dtAdjList);
		finish = std::chrono::high_resolution_clock::now();
		cout << "MST created in " << std::chrono::duration<double>(finish - start).count() << " secs\n";

		start = std::chrono::high_resolution_clock::now();
		reconstruct(initialEdges);
		finish = std::chrono::high_resolution_clock::now();
		cout << pts.size() << " Reconstructed in " << std::chrono::duration<double>(finish - start).count() << "\n";
	}

	/**********************************************************************************************/ /**
	 * @fn	bool isValid()
	 *
	 * @brief	Query if this object is valid.
	 *
	 * @return	True if valid, false if not.
	 **************************************************************************************************/

	bool isValid() { return valid; }

	void writeModel(const char *outputFilePath) {
		std::ofstream outputFile(outputFilePath);
		outputFile << pts.size() << "\n";
		for (Point3D point : pts) {
			outputFile << point << "\n";
		}

		outputFile << modelEdges.size() << "\n";
		for (auto edge : modelEdges) {
			outputFile << edge[0] << " " << edge[1] << "\n";
		}

		outputFile << modelFaces.size() << "\n";
		for (auto triangle : modelFaces) {
			outputFile << triangle[0] << " " << triangle[1] << " " << triangle[2] << "\n";
		}
	}
};

int main(int argc, char *argv[]) {

	if (argc != 3) {
		cerr << "Invalid Argument\n";
		return 1;
	}

	println("Input File =", argv[1]);
	println("Output File =", argv[2]);

	cout << std::boolalpha;

	SurfaceReconstruct surface(argv[1]);
	if (surface.isValid()) {
		surface.writeModel(argv[2]);
	}

	return 0;
}
