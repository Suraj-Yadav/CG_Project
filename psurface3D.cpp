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
typedef Kernel::Plane_3 Plane3D;
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
	
	vector<vector<set<ourEdge>::iterator>> vertAdjEdges;

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
			fprintln(logFile, "Angle between the face (", edge[0], edge[1], edgeDegree.back(), ") and (", edge[0], edge[1], newPoint, ") is", angle);
		}

		//if (!ans && debug) {
		//	fprint(logFile, "Problem Adding ", newPoint, "to", edge, ". edgedegree={");
		//	for (int p : edgeDegree)
		//		fprint(logFile, p);
		//	fprintln(logFile, "}");
		//}*/
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

	bool isEdgeOverlap(const ourEdge &e) {
		//cout << "testing " << e << std::endl;
		ourEdge edge(0, 1);
		for (int i = 0; i < 2; ++i) {
			for (auto &elem : vertAdjEdges[e[i]]) {
				edge = *elem;
				//cout << "with " << edge << " " << std::endl;
				Point3D a = pts[edge[0]], b = pts[edge[1]], c, o = pts[e[i]];
				Plane3D plane(a, b, o);
				c = plane.projection(pts[e[!i]]);
				auto d1 = CGAL::cross_product(a - o, c - o) * CGAL::cross_product(a - o, b - o),
					 d2 = CGAL::cross_product(b - o, c - o) * CGAL::cross_product(b - o, a - o);
				if (d1 > 0 && d2 > 0)
					return true;
			}
		}
		//cout << "pass" << std::endl;
		return false;
	}

	bool isFaceOverlap(const ourFace &f){
		ourEdge e1(f[0],f[1]);
		ourEdge e2(f[2],f[1]);
		ourEdge e3(f[0],f[2]);
		return (isEdgeOverlap(e1) || isEdgeOverlap(e2) || isEdgeOverlap(e3));
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
			int v = INT_MAX;
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
	
	int processCycle(vector<int> &cycle){
		static int cycleNo = 1;
		int facesAdded = 0;
		// 1. get all triangles corresponding to the adjacent vertices in the cycle
		// 2. store them in a priority queue with the score
		// 3. get the best one, and insert the new triangles
		cout << "Cycle Number "<<cycleNo++<<std::endl;
		cout<<"Cycle Size :"<<cycle.size()<<std::endl;
		for(int c:cycle){
			cout<<c<<'\n';
		}
		cout<<"End of Cycle\n";

		set<std::tuple<double, pair<int,int>, int>, std::greater<std::tuple<double,pair<int,int>,int>>> pq;
		set<std::tuple<double, pair<int,int>, int>, std::greater<std::tuple<double,pair<int,int>,int>>> checkedFaces;
		cout<<"Reached "<< 0 <<std::endl;

		for(int i=0;i<cycle.size();i++){
			int j=(i+1) % cycle.size();
			int k=(i+2) % cycle.size();
			bool overlap = isEdgeOverlap(ourEdge(cycle[i],cycle[k]));
			if(!overlap)
				pq.insert(make_tuple(getScore(cycle[i],cycle[j],cycle[k]),make_pair(i,k),j));
			checkedFaces.insert(make_tuple(getScore(cycle[i],cycle[j],cycle[k]),make_pair(i,k),j));
		}

		cout<<"Reached "<< 1 << std::endl;

		while(!pq.empty()){
			facesAdded++;
			auto face = *pq.begin();
			int v1 = get<1>(face).first;
			int v3 = get<1>(face).second;
			int v2 = get<2>(face);
			pq.erase(face);
			ourEdge e1(cycle[v1],cycle[v2]);
			ourEdge e2(cycle[v1],cycle[v3]);
			ourEdge e3(cycle[v3],cycle[v2]);
			cout<<"Reached "<< 2 << std::endl;
			cout<<"Face:\n";
			cout<< ourFace(cycle[v1],cycle[v2],cycle[v3])<<'\n';
			cout<<"End of Face\n";
			if(cycle[v1]==cycle[v2] || cycle[v2]==cycle[v3] || cycle[v1]==cycle[v3])
				continue;
			//if(!isEdgeOverlap(ourEdge(cycle[v1],cycle[v3]))) {
			if(1){
				cout<<"Reached "<< 3 << std::endl;
				modelFaces.insert(ourFace(cycle[v1],cycle[v2],cycle[v3]));
				modelEdges.insert(e1);
				modelEdges.insert(e2);
				modelEdges.insert(e3);
				
				int v0 = (v1-1)%cycle.size();
				int v4 = (v3+1)%cycle.size();
				
				if(v0<0) 
					v0 = v0 + cycle.size();

				if(v4<0) 
					v4 = v4 + cycle.size();

				cout<<"Reached "<< 3.2 << std::endl;
				
				pq.erase(make_tuple(getScore(cycle[v0],cycle[v1],cycle[v2]),make_pair(v0,v2),v1));
				cout<<"Reached "<< 3.3 << std::endl;
				cout<< "Point Index "<<v2<<' '<<v3<<' '<<v4<<".\nPoints Size = "<<cycle.size()<<std::endl;

				cout<<"points : "<<cycle[v2]<<' '<<cycle[v3]<<' '<<cycle[v4]<<std::endl;

				pq.erase(make_tuple(getScore(cycle[v2],cycle[v3],cycle[v4]),make_pair(v2,v4),v3));

				cout<<"Reached "<< 3.4 << std::endl;
				auto newItem1 = make_tuple(getScore(cycle[v0],cycle[v1],cycle[v3]),make_pair(v0,v3),v1);
				auto newItem2 = make_tuple(getScore(cycle[v1],cycle[v3],cycle[v4]),make_pair(v1,v4),v3);
				cout<<"Reached "<< 4 << std::endl;

				// if newItem1 is not in checkedFaces and newItem1 is not isEdgeOverlap
					// insert in pq
					// do same for newItem2

				if(!isEdgeOverlap(ourEdge(cycle[v0],cycle[v3])) && 
					checkedFaces.find(newItem1)==checkedFaces.end()){
						pq.insert(newItem1);	
				}

				if(!isEdgeOverlap(ourEdge(cycle[v1],cycle[v4])) && 
					checkedFaces.find(newItem2)==checkedFaces.end()){
						pq.insert(newItem2);	
				}
				cout<<"Reached "<< 5 << std::endl;
				checkedFaces.insert(newItem1);
				checkedFaces.insert(newItem2);
				
			}
			cout<<"Reached "<< 6 << std::endl;
			pq.erase(face);
			cout<<"Reached "<< 7 << std::endl;
		}
		cout<<"Done"<<std::endl;
		cout<<"CheckedFaces\n";
		for(auto elem: checkedFaces){
			int v1 = get<1>(elem).first;
			int v3 = get<1>(elem).second;
			int v2 = get<2>(elem);
			cout<<cycle[v1]<<' '<<cycle[v2]<<' '<<cycle[v3]<<'\n';
		}
		cout<<"End of Checked Faces\n";
		return facesAdded;
	}
	
	void pushCycles(int point, map<int,set<int>> &tempAdjList, vector<vector<int>> &allCycles,map<int,bool> &inCycle){

		stack<int> s;
		vector<int> eu; // for the euler's circuit
		int a=point;
		
		while(!s.empty() || tempAdjList[a].size()!=0){
			if(tempAdjList[a].size()){
				s.push(a);
				int t=*tempAdjList[a].begin();
				tempAdjList[a].erase(t);
				tempAdjList[t].erase(a);
				a=t;
			}
			else{
				eu.push_back(a);
				a=s.top();
				s.pop();
			}
		}
		
		map<int,bool> visited;	
		
		for(int a: eu){
			if(!visited[a]){
				s.push(a);
				visited[a]=true;
			}
			else{
				vector<int> cycle;
				cycle.push_back(a);
				inCycle[a]=true;
				while(s.top() != a){
					cycle.push_back(s.top());
					inCycle[s.top()]=true;
					s.pop();
				}
				allCycles.push_back(cycle);
				
			}
		}

	}

	vector< vector<int> > getAllCycles(map<int,set<int>> &tempAdjList, set<int> &leftPts){
		vector< vector<int> > allCycles;
		map<int,bool> inCycle;
		for(int point: leftPts){
			if(!inCycle[point]){
				pushCycles(point,tempAdjList,allCycles,inCycle);
			}
		}
		return allCycles;
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
					modelEdges.insert(edge);

					edgeDegree[edge].push_back(point);
					edgeDegree[newEdge1].push_back(edge[1]);
					edgeDegree[newEdge2].push_back(edge[0]);

					modelFaces.insert(face);
					vertAdjEdges[edge[0]].push_back(modelEdges.find(newEdge2));
					vertAdjEdges[edge[1]].push_back(modelEdges.find(newEdge1));
					vertAdjEdges[point].push_back(modelEdges.find(edge));

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
			set<int> leftPts;
			
			set<ourEdge> nextEdges;
			map<int,set<int>> tempAdjList;

			for (auto &elem : modelEdges) {
				if (edgeDegree[elem].size() == 1) {
					leftVerts.insert({elem[0], false});
					leftVerts.insert({elem[1], false});
					
					leftPts.insert(elem[0]);
					leftPts.insert(elem[1]);
					
					tempAdjList[elem[0]].insert(elem[1]);
					tempAdjList[elem[1]].insert(elem[0]);
					
				}
			}
			cout<<"Left Vertices\n";
			for(int v:leftPts){
				cout<<v<<'\n';
			}
			cout<<"ENd of leftPts\n";
			vector<vector<int>> cycles = getAllCycles(tempAdjList,leftPts);
			cout << "All Cycles obtained\n";
			cout<< "Total Cycles: " << cycles.size() << std::endl;
			int facesAdded = 0;
			for(vector<int> cycle: cycles){
				cout<< cycle.size()<<std::endl;
				facesAdded +=processCycle(cycle);
			}
			cout<<"Total Faces Added "<<facesAdded<<'\n';
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
		vertAdjEdges.resize(pts.size());

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
