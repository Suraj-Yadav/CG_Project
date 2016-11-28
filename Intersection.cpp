#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
#include <iostream>

using Kernel = CGAL::Simple_cartesian<double>;
using Point2D = Kernel::Point_2;
using Segment2D = Kernel::Segment_2;
using Intersect2D = Kernel::Intersect_2;

int main() {
	Point2D p[4];
	std::cin >> p[0] >> p[1] >> p[2] >> p[3];
	Segment2D line1(p[0], p[1]), line2(p[2], p[3]);

	auto result = intersection(line1, line2);
	if (result) {
		if (const Segment2D* s = boost::get<Segment2D>(&*result)) {
			std::cout << *s << std::endl;
		}
		else {
			const Point2D* p = boost::get<Point2D>(&*result);
			std::cout << *p << std::endl;
		}
	}
	return 0;
}