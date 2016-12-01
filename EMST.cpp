#include <fstream>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

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

// for viewportsBbox
#include <CGAL/Qt/utility.h>

// the two base classes
#include "ui_EMST.h"
#include <CGAL/Qt/DemosMainWindow.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point2D = Kernel::Point_2;
using Iso_rectangle2D = Kernel::Iso_rectangle_2;

using Delaunay = CGAL::Delaunay_triangulation_2<Kernel>;

class MainWindow : public CGAL::Qt::DemosMainWindow,
				   public Ui::EMSTWindow {
	Q_OBJECT

  private:
	Delaunay dt;
	QGraphicsScene scene;
	CGAL::Qt::TriangulationGraphicsItem<Delaunay> *dgi;
	CGAL::Qt::VoronoiGraphicsItem<Delaunay> *vgi;

  public:
	MainWindow();

  public Q_SLOTS:

	void processInput(CGAL::Object o);
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

	//
	// Manual handling of actions
	//

	QObject::connect(this->actionQuit, SIGNAL(triggered()),
					 this, SLOT(close()));

	// Check two actions
	// this->actionInsertPoint->setChecked(true);
	this->actionShowDelaunay->setChecked(true);

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

void MainWindow::processInput(CGAL::Object o) {
	Point2D p;
	if (CGAL::assign(p, o)) {
		dt.insert(p);
	}
	Q_EMIT(changed());
}

/* 
 *  Qt Automatic Connections
 *  http://doc.qt.io/qt-5/designer-using-a-ui-file.html#automatic-connections
 * 
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */
void MainWindow::on_actionShowDelaunay_toggled(bool checked) {
	dgi->setVisibleEdges(checked);
}

void MainWindow::on_actionShowVoronoi_toggled(bool checked) {
	vgi->setVisible(checked);
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

	Kernel::Point_2 p;
	std::vector<Kernel::Point_2> points;
	while (ifs >> p) {
		// ignore whatever comes after x and y
		ifs.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
		points.push_back(p);
	}
	dt.insert(points.begin(), points.end());

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

#include "EMST.moc"
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
