#include "Scene.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Union_find.h>

using namespace std;

void Scene::readOFFPointsandFacets(const char* filename,
								   std::list<Point_3>& points) {
	ifstream fin;
	fin.open(filename);
	// Check whether the file is opened properly
	if (!fin) {
		showError(QObject::tr("Error: cannot open file %1 for reading.").arg(filename));
		return;
	}

	istream* pIn = &fin;

	// Use CGAL::File_scanner_OFF to read in data
	CGAL::File_scanner_OFF scanner(*pIn);
	if (!(*pIn)) {
		showError(QObject::tr("Input error: file %1 is not in OFF format.").arg(filename));
		return;
	}
	if (scanner.size_of_vertices() <= 0) {
		showError(QObject::tr("Input error: file %1 has no vertices.").arg(filename));
		return;
	}
	// Get points data from scanner
	double x, y, z;
	for (std::size_t i = 0; i < scanner.size_of_vertices(); ++i) {
		scanner.scan_vertex(x, y, z);
		Point_3 pt(x, y, z);
		points.push_back(pt);
		scanner.skip_to_next_vertex(i);
	} //end-for-vertex
}

void Scene::loadPointsOFF(const char* filename) {
	list<Point_3> pts;

	/* Read point data from file */
	/* 1. use CGAL::File_scanner_OFF to read in data --tested */
	readOFFPointsandFacets(filename, pts);

	/* 2. use CGAL::read_off_points to read in data -- tested */
	/* Note: read in points only, i.e. normals and faces are ignored */
	/* Note: this function can NOT omit comments (starting with '#') */
	//  ifstream fin;
	//  fin.open( filename );
	// check whether the file is opened properly
	//  if( !fin ) {
	//    showError( QObject::tr("Error: cannot open file %1 for reading.").arg(filename) );
	//    return;
	//  }
	//  if ( !CGAL::read_off_points( fin,  // inout ifstream
	//                               back_inserter(pts) ) ) {  // output iterator over points
	//    showError( QObject::tr("Error: cannot read file %1.").arg(filename) );
	//  }

	/* Insert the points to build a Delaunay triangulation */
	/* Note: this function returns the number of inserted points;
      it is not guaranteed to insert the points following the order of iteraror. */
	m_dt.insert(pts.begin(), pts.end());
	/* Check the combinatorial validity of the triangulation */
	/* Note: when it is set to be true,
      messages describing the first invalidity encountered are printed. */
		
		showError(QObject::tr("Error: fail to build a Delauna4y triangulation."));
	/* Check the dimension */
	if (m_dt.dimension() != 3)
		showError(QObject::tr("Error: cannot built a 3D triangulation."));
	/* Store the vertex handles into an array for future usage (move, delete, etc) */
	for (vertices_iterator vit = m_dt.finite_vertices_begin();
		 vit != m_dt.finite_vertices_end(); ++vit) {
		m_vhArray.push_back(vit);
	}

	updateViewer();
	updateMST();

	assert(m_dt.number_of_vertices() == (std::size_t)m_vhArray.size());
}

void Scene::updateViewer() {
	FT minX = std::numeric_limits<FT>::max(),
	   maxX = std::numeric_limits<FT>::lowest(),
	   minY = std::numeric_limits<FT>::max(),
	   maxY = std::numeric_limits<FT>::lowest(),
	   minZ = std::numeric_limits<FT>::max(),
	   maxZ = std::numeric_limits<FT>::lowest();

	for (auto vertHandle : m_vhArray) {
		auto p = vertHandle->point();
		minX = std::min(minX, p.x());
		maxX = std::max(maxX, p.x());
		minY = std::min(minY, p.y());
		maxY = std::max(maxY, p.y());
		minZ = std::min(minZ, p.z());
		maxZ = std::max(maxZ, p.z());
	}

	FT rangeX = maxX - minX,
	   rangeY = maxY - minY,
	   rangeZ = maxZ - minZ,
	   midX = (maxX + minX) / 2,
	   midY = (maxY + minY) / 2,
	   midZ = (maxZ + minZ) / 2;

	m_viewer->setSceneCenter(qglviewer::Vec(midX, midY, midZ));
	m_viewer->setSceneRadius(std::sqrt(rangeX * rangeX + rangeY * rangeY + rangeZ * rangeZ));
}

void Scene::updateMST() {
	std::vector<Segment_3> allEdges;

	m_mst.clear();
	CGAL::Union_find<Point_3> uf;

	std::map<Point_3, CGAL::Union_find<Point_3>::handle> mapping;

	for (auto vertHandle : m_vhArray)
		mapping.insert({vertHandle->point(), uf.make_set(vertHandle->point())});

	for (auto edgeItr = m_dt.finite_edges_begin(); edgeItr != m_dt.finite_edges_end(); edgeItr++)
		allEdges.push_back(m_dt.segment(*edgeItr));

	assert(uf.number_of_sets() == (std::size_t)m_vhArray.size());

	std::sort(allEdges.begin(),
			  allEdges.end(),
			  [](Segment_3& a, Segment_3& b) -> bool { return a.squared_length() < b.squared_length(); });

	for (auto edge : allEdges) {
		if (!uf.same_set(mapping[edge.start()], mapping[edge.end()])) {
			m_mst.push_back(edge);
			uf.unify_sets(mapping[edge.start()], mapping[edge.end()]);
		}
	}
}

void Scene::loadPointsXYZ(const char* filename) {
	ifstream fin;
	fin.open(filename);
	// Check whether the file is opened properly
	if (!fin) {
		showError(QObject::tr("Error: cannot open file %1 for reading.").arg(filename));
		return;
	}

	/* Use CGAL::read_xyz_points to read in data -- tested */
	/* Note: this function reads in points only (normals are ignored) */
	/* Note: this function can NOT omit comments (starting with '#') */
	list<Point_3> pts;
	if (!CGAL::read_xyz_points(fin,					  // input ifstream
							   back_inserter(pts))) { // output iterator over points
		showError(QObject::tr("Error: cannot read file %1.").arg(filename));
	}

	/* Insert the points to build a Delaunay triangulation */
	/* Note: this function returns the number of inserted points;
      it is not guaranteed to insert the points following the order of iteraror. */
	m_dt.insert(pts.begin(), pts.end());
	/* Check the combinatorial validity of the triangulation */
	/* Note: when it is set to be true,
      messages describing the first invalidity encountered are printed. */
	if (!m_dt.is_valid()) // default: false - verbosity off
		showError(QObject::tr("Error: fail to build a Delaunay triangulation."));
	/* Check the dimension */
	if (m_dt.dimension() != 3)
		showError(QObject::tr("Error: cannot build a 3D triangulation."));
	/* Store the vertex handles into an array for future usage (move, delete, etc) */
	for (vertices_iterator vit = m_dt.finite_vertices_begin();
		 vit != m_dt.finite_vertices_end(); ++vit) {
		m_vhArray.push_back(vit);
	}

	updateViewer();
	updateMST();

	assert(m_dt.number_of_vertices() == (std::size_t)m_vhArray.size());
}

void Scene::savePointsOFF(const char* filename) {
	ofstream fout;
	fout.open(filename);
	if (!fout) {
		showError(QObject::tr("Error: cannot open file %1 for writting.").arg(filename));
		return;
	}

	ostream* pOut = &fout;

	/* Use CGAL::File_writer_OFF to write points */
	// initialize header_OFF
	CGAL::File_header_OFF header(false, // true: binary output; false: ASCII
								 false, // true: no comments in file
								 false, // true: Geomview SKEL format
								 true); // true: verbosity on; false: verbosity off
	// a simpler way to initialize header_OFF
	//  CGAL::File_header_OFF header(true);  // true: verbosity on
	//                                       // (ASCII output, comments, no SKEL)
	CGAL::File_writer_OFF writer(header);
	// write header
	writer.write_header(*pOut,					   // output ostream
						m_dt.number_of_vertices(), // number of points/vertices
						0,						   // number of halfedges
						0,						   // number of facets
						false);					   // true: has normals
	// write points (get from point array)
	for (vertices_iterator vit = m_dt.finite_vertices_begin();
		 vit != m_dt.finite_vertices_end(); ++vit) {
		Point_3& p = vit->point();
		writer.write_vertex(p.x(), p.y(), p.z());
	}
	// write footer
	writer.write_footer();
}

void Scene::savePointsXYZ(const char* filename) {
	ofstream fout;
	fout.open(filename);
	// Check whether the file is opened properly
	if (!fout) {
		showError(QObject::tr("Error: cannot open file %1 for writting.").arg(filename));
		return;
	}

	/* Use CGAL::write_xyz_points to write out data */
	/* Note: this function writes out points only (normals are ignored) */
	if (!CGAL::write_xyz_points(fout,				  // output ofstream
								m_dt.points_begin(),  // first output point
								m_dt.points_end())) { // past-the-end output point
		showError(QObject::tr("Error: cannot read file %1.").arg(filename));
	}
}
