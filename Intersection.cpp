#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/intersections.h>
#include <CGAL/IO/read_xyz_points.h>

#include <bits/stdc++.h>

//CGAL::Exact_predicates_exact_constructions_kernel

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Line_3 Line3D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Intersect_3 Intersect3D;
typedef Kernel::Triangle_3 Triangle3D;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;

using std::cerr;
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
	// Point3D p[2];
	// Triangle3D tri1({-0.28358, -10.0897, -0.86783}, {-0.42034, -10.0688, -0.43415}, {0.34223, -9.9244, -0.3875});
	// Segment3D line({-0.28358, -10.0897, -0.86783}, {0.89782, -10.0743, 0.26342});
	// std::cin >> p[0] >> p[1];
	//Triangle3D tri({-1.0 , -1.0 , 0.0},{1.0 , -1.0 , 0.0},{0 , 0 , 0.0});
	//Segment3D line({0.0 , 0.0 , 0.0},{0.0 , -1.5 , 1.0});

	// auto d = tri1.supporting_plane().projection(line[1]);
	// Triangle3D tri2(line[0], line[1], d + (d - line[1]));

	// std::cout << p[0] << "\n";
	// std::cout << p[1] << "\n";

	// auto result = intersection(tri1, tri2);
	// if (result) {
	// if (const Segment3D *s = boost::get<Segment3D>(&*result)) {
	// std::cout << *s << std::endl;
	// }
	// else {
	// const Point3D *p = boost::get<Point3D>(&*result);
	// std::cout << *p << std::endl;
	// }
	// }
	// else {
	// std::cout << "Nope\n";
	// }

	std::ifstream inputFile(argv[1]);

	DT3 m_dt;

	std::list<Point3D> pts;
	if (!CGAL::read_xyz_points(inputFile,			  // input ifstream
							   back_inserter(pts))) { // output iterator over points
		// showError(QObject::tr("Error: cannot read file %1.").arg(filename));
		cerr << "Error: cannot read file.";
	}
	cout << argv[1] << "\n";
	cout << "Read " << pts.size() << " points ";
	/* Insert the points to build a Delaunay triangulation */
	/* Note: this function returns the number of inserted points;
      it is not guaranteed to insert the points following the order of iteraror. */
	m_dt.insert(pts.begin(), pts.end());
	cout << "Inserted ";
	/* Check the combinatorial validity of the triangulation */
	/* Note: when it is set to be true,
      messages describing the first invalidity encountered are printed. */
	if (!m_dt.is_valid()) { // default: false - verbosity off
		cout << "Error: fail to build a Delaunay triangulation." << std::endl;
		return -1;
	}
	/* Check the dimension */
	if (m_dt.dimension() != 3) {
		cout << "Error: cannot built a 3D triangulation.\n Current dimension = " << m_dt.dimension() << std::endl;
		return -1;
	}
	cout << "Good\n";
	return 0;
}