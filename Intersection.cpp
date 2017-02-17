#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Line_3 Line3D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Intersect_3 Intersect3D;
typedef Kernel::Triangle_3 Triangle3D;

int main() {
	Point3D p[2];
	Triangle3D tri1({-0.28358, -10.0897, -0.86783}, {-0.42034, -10.0688, -0.43415}, {0.34223, -9.9244, -0.3875});
	Segment3D line({-0.28358, -10.0897, -0.86783}, {0.89782, -10.0743, 0.26342});
	// std::cin >> p[0] >> p[1];
	//Triangle3D tri({-1.0 , -1.0 , 0.0},{1.0 , -1.0 , 0.0},{0 , 0 , 0.0});
	//Segment3D line({0.0 , 0.0 , 0.0},{0.0 , -1.5 , 1.0});

	auto d = tri1.supporting_plane().projection(line[1]);
	Triangle3D tri2(line[0], line[1], d+(d-line[1]));

	std::cout << p[0] << "\n";
	std::cout << p[1] << "\n";

	auto result = intersection(tri1, tri2);
	if (result) {
		if (const Segment3D *s = boost::get<Segment3D>(&*result)) {
			std::cout << *s << std::endl;
		}
		else {
			const Point3D *p = boost::get<Point3D>(&*result);
			std::cout << *p << std::endl;
		}
	}
	else {
		std::cout << "Nope\n";
	}
	return 0;
}