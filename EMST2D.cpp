#include <fstream>
#include <map>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Union_find.h>

#include <CGAL/point_generators_2.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/Qt/VoronoiGraphicsItem.h>
#include <CGAL/Qt/SegmentsGraphicsItem.h>

// for viewportsBbox
#include <CGAL/Qt/utility.h>

// the two base classes
#include "ui_EMST2D.h"
#include <CGAL/Qt/DemosMainWindow.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point2D = Kernel::Point_2;
using Segment2D = Kernel::Segment_2;
using Iso_rectangle2D = Kernel::Iso_rectangle_2;

using Delaunay = CGAL::Delaunay_triangulation_2<Kernel>;

class UnionFind {
	std::vector<unsigned> id;
	std::vector<unsigned> rank;

  public:
	UnionFind(unsigned N) {
		id.resize(N);
		rank.resize(N);
		for (unsigned i = 0; i < N; i++) {
			id[i] = i;
			rank[i] = 0;
		}
	}
	unsigned find(unsigned x) {
		if (id[x] != x)
			id[x] = find(id[x]);
		return id[x];
	}
	void Union(unsigned x, unsigned y) {
		unsigned first = find(x);
		unsigned second = find(y);
		if (first == second)
			return;
		if (rank[first] > rank[second]) {
			id[second] = first;
		}
		else if (rank[first] < rank[second]) {
			id[first] = second;
		}
		else if (second != first) {
			id[second] = first;
			rank[first]++;
		}
	}
};

class MainWindow : public CGAL::Qt::DemosMainWindow,
				   public Ui::EMST2DWindow {
	Q_OBJECT

  private:
	Delaunay dt;
	std::vector<Segment2D> mstSegments;

	QGraphicsScene scene;
	CGAL::Qt::TriangulationGraphicsItem<Delaunay> *dgi;
	CGAL::Qt::VoronoiGraphicsItem<Delaunay> *vgi;
	CGAL::Qt::SegmentsGraphicsItem<std::vector<Segment2D>> *sgi;

  public:
	MainWindow();

  public Q_SLOTS:
	void on_actionShowPoints_toggled(bool checked);
	void on_actionShowMST_toggled(bool checked);
	void on_actionShowDelaunay_toggled(bool checked);
	void on_actionShowVoronoi_toggled(bool checked);
	void on_actionLoadPoints_triggered();
	void on_actionSavePoints_triggered();
	void on_actionClear_triggered();
	void on_actionRecenter_triggered();
	virtual void open(QString fileName);

  Q_SIGNALS:
	void changed();
};

MainWindow::MainWindow()
	: DemosMainWindow() {
	setupUi(this);

	this->graphicsView->setAcceptDrops(false);

	// Add a GraphicItem for the Delaunay triangulation
	dgi = new CGAL::Qt::TriangulationGraphicsItem<Delaunay>(&dt);

	QObject::connect(this, SIGNAL(changed()),
					 dgi, SLOT(modelChanged()));

	dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
	scene.addItem(dgi);

	// Add a GraphicItem for the Voronoi diagram
	vgi = new CGAL::Qt::VoronoiGraphicsItem<Delaunay>(&dt);

	QObject::connect(this, SIGNAL(changed()),
					 vgi, SLOT(modelChanged()));

	vgi->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
	scene.addItem(vgi);
	vgi->hide();

	// Add a GraphicItem for the MST
	sgi = new CGAL::Qt::SegmentsGraphicsItem<std::vector<Segment2D>>(&mstSegments);

	QObject::connect(this, SIGNAL(changed()),
					 sgi, SLOT(modelChanged()));
	scene.addItem(sgi);
	sgi->hide();

	//
	// Manual handling of actions
	//

	QObject::connect(this->actionQuit, SIGNAL(triggered()),
					 this, SLOT(close()));

	// Check two actions
	// this->actionInsertPoint->setChecked(true);
	this->actionShowDelaunay->setChecked(true);
	this->actionShowPoints->setChecked(true);
	this->actionShowMST->setChecked(true);

	//
	// Setup the scene and the view
	//
	scene.setItemIndexMethod(QGraphicsScene::NoIndex);
	scene.setSceneRect(-100, -100, 100, 100);
	this->graphicsView->setScene(&scene);
	this->graphicsView->setMouseTracking(true);

	// Turn the vertical axis upside down
	this->graphicsView->scale(1, -1);

	// The navigation adds zooming and translation functionality to the
	// QGraphicsView
	this->addNavigation(this->graphicsView);

	this->setupStatusBar();
	this->setupOptionsMenu();
	this->addAboutDemo(":/cgal/help/about_Delaunay_triangulation_2.html");
	this->addAboutCGAL();
	this->setupExportSVG(actionExport_SVG, graphicsView);

	this->addRecentFiles(this->menuFile, this->actionQuit);
	connect(this, SIGNAL(openRecentFile(QString)),
			this, SLOT(open(QString)));
}

/*
 *  Qt Automatic Connections
 *  http://doc.qt.io/qt-5/designer-using-a-ui-file.html#automatic-connections
 *
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */
void MainWindow::on_actionShowPoints_toggled(bool checked) {
	dgi->setVisibleVertices(checked);
}

void MainWindow::on_actionShowDelaunay_toggled(bool checked) {
	dgi->setVisibleEdges(checked);
}

void MainWindow::on_actionShowVoronoi_toggled(bool checked) {
	vgi->setVisible(checked);
}

void MainWindow::on_actionShowMST_toggled(bool checked){
	sgi->setVisible(checked);
}

void MainWindow::on_actionClear_triggered() {
	dt.clear();
	Q_EMIT(changed());
}

void MainWindow::on_actionLoadPoints_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this,
													tr("Open Points file"),
													".");
	if (!fileName.isEmpty()) {
		open(fileName);
	}
}

void MainWindow::open(QString fileName) {
	// wait cursor
	QApplication::setOverrideCursor(Qt::WaitCursor);
	std::ifstream ifs(qPrintable(fileName));

	dt.clear();
	mstSegments.clear();

	Kernel::Point_2 p;
	std::vector<Kernel::Point_2> points;
	while (ifs >> p) {
		// ignore whatever comes after x and y
		ifs.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
		points.push_back(p);
	}
	dt.insert(points.begin(), points.end());
	
	std::vector<Segment2D> allEdges;

	std::map<Point2D, CGAL::Union_find<Point2D>::handle> mapping;
	CGAL::Union_find<Point2D> uf;

	for (auto vertItr = dt.finite_vertices_begin(); vertItr != dt.finite_vertices_end(); ++vertItr) {
		auto p = vertItr->handle()->point();
		mapping.insert({p, uf.make_set(p)});
	}

	for (auto edgeItr = dt.finite_edges_begin(); edgeItr != dt.finite_edges_end(); edgeItr++)
		allEdges.push_back(dt.segment(*edgeItr));
	
	std::sort(allEdges.begin(),
			  allEdges.end(),
			  [](Segment2D & a, Segment2D & b) -> bool { return a.squared_length() < b.squared_length(); });

	for (auto edge : allEdges) {
		if (!uf.same_set(mapping[edge.start()], mapping[edge.end()])) {
			mstSegments.push_back(edge);
			uf.unify_sets(mapping[edge.start()], mapping[edge.end()]);
		}
	}

	// default cursor
	QApplication::restoreOverrideCursor();
	this->addToRecentFiles(fileName);
	actionRecenter->trigger();
	Q_EMIT(changed());
}

void MainWindow::on_actionSavePoints_triggered() {
	QString fileName = QFileDialog::getSaveFileName(this,
													tr("Save points"),
													".");
	if (!fileName.isEmpty()) {
		std::ofstream ofs(qPrintable(fileName));
		for (Delaunay::Finite_vertices_iterator
				 vit = dt.finite_vertices_begin(),
				 end = dt.finite_vertices_end();
			 vit != end; ++vit) {
			ofs << vit->point() << std::endl;
		}
	}
}

void MainWindow::on_actionRecenter_triggered() {
	this->graphicsView->setSceneRect(dgi->boundingRect());
	this->graphicsView->fitInView(dgi->boundingRect(), Qt::KeepAspectRatio);
}

#include "EMST2D.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv) {
	QApplication app(argc, argv);

	app.setOrganizationDomain("geometryfactory.com");
	app.setOrganizationName("GeometryFactory");
	app.setApplicationName("Delaunay_triangulation_2 demo");

	// Import resources from libCGAL (QT5).
	CGAL_QT_INIT_RESOURCES;

	MainWindow mainWindow;
	mainWindow.show();

	QStringList args = app.arguments();
	args.removeAt(0);
	Q_FOREACH (QString filename, args) {
		mainWindow.open(filename);
	}

	return app.exec();
}
