#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/IO/read_xyz_points.h>
#include <iostream>
#include <fstream>
#include <map>

#define val(x) std::cout << #x "=" << x << "\n";

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Segment_3 Segment3D;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;

int main(int argc, char *argv[]) {

	if (argc != 3) {
		std::cerr << "Invalid Argument\n";
		return 1;
	}

	std::ifstream inputFile(argv[1]);
	std::ofstream outputFile(argv[2]);

	std::vector<Point3D> points;
	DT3 dt;

	val(points.size());
	if (!CGAL::read_xyz_points(inputFile, std::back_inserter(points))) { // output iterator over points
		std::cerr << "Error: cannot read file.";
	}
	val(points.size());

	dt.insert(points.begin(), points.end());

	if (!dt.is_valid(true)) {
		std::cerr << "Error: fail to build a Delaunay triangulation.";
		return 0;
	}
	if (dt.dimension() != 3) {
		std::cerr << "Error: cannot built a 3D triangulation.";
		return 0;
	}

	std::vector<Segment3D> allEdges, mst;

	mst.clear();
	CGAL::Union_find<Point3D> uf;

	std::map<Point3D, CGAL::Union_find<Point3D>::handle> mapping;
	std::map<Point3D, unsigned> index;

	unsigned i = 0;
	for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
		mapping.insert({vit->point(), uf.make_set(vit->point())});
		if (index.find(vit->point()) == index.end())
			index.insert({vit->point(), i++});
	}

	// val(i);

	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++) {
		allEdges.push_back(dt.segment(*edgeItr));
	}

	// std::cout<<allEdges.size()<<"\n";
	// for(auto edge:allEdges)
		// std::cout << index[edge.start()] << " " << index[edge.end()] << "\n";

	std::sort(allEdges.begin(),
		allEdges.end(),
		[](Segment3D& a, Segment3D& b) -> bool { return a.squared_length() < b.squared_length(); });

	for (auto edge : allEdges) {
		if (mapping.find(edge.start()) == mapping.end()) {
			val(edge.start());
			return 0;
		}
		if (mapping.find(edge.end()) == mapping.end()) {
			val(edge.end());
			return 0;
		}
		if (!uf.same_set(mapping[edge.start()], mapping[edge.end()])) {
			mst.push_back(edge);
			uf.unify_sets(mapping[edge.start()], mapping[edge.end()]);
		}
	}

	// val(mapping.size());
	// val(allEdges.size());
	// val(mst.size());
	// val(points.size());

	/*if(mst.size()+1!=points.size())
		std::cerr<<"Error\n";
*/
	points.clear();
	points.resize(mapping.size());
	for (auto point : index) {
		points[point.second] = point.first;
	}
	outputFile << points.size() << "\n";
	for (auto point : points) {
		outputFile << point << "\n";
	}

	outputFile << mst.size() << "\n";

	for (auto edge : mst) {
		outputFile << index[edge.start()] << " " << index[edge.end()] << "\n";
	}

	return 0;
}