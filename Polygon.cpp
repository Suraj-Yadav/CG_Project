#include <fstream>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

// Qt headers
#include <QFileDialog>
#include <QGraphicsLineItem>
#include <QInputDialog>
#include <QString>
#include <QtGui>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/PolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/PolygonGraphicsItem.h>
#include <CGAL/Qt/LineGraphicsItem.h>

// the two base classes
#include "ui_Polygon.h"
#include <CGAL/Qt/DemosMainWindow.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point2D = Kernel::Point_2;
using Segment2D = Kernel::Segment_2;
using Line2D = Kernel::Line_2;

// it must be a list for the partition

using Polygon2D = CGAL::Polygon_2<Kernel, std::list<Point2D>>;

class MainWindow : public CGAL::Qt::DemosMainWindow, public Ui::PolygonWindow {
	Q_OBJECT

  private:
	CGAL::Qt::Converter<Kernel> convert;
	Polygon2D poly, kgon;
	QGraphicsScene scene;

	CGAL::Qt::PolygonGraphicsItem<Polygon2D> *pgi;

	CGAL::Qt::GraphicsViewPolylineInput<Kernel> *pi;

	CGAL::Qt::LineGraphicsItem<Kernel> *lgi;

	CGAL::Qt::PolygonGraphicsItem<Polygon2D> *kgongi;

  public:
	MainWindow();

  public Q_SLOTS:

	void processInput(CGAL::Object o);

	void on_actionClear_triggered();

	void on_actionLoadPolygon_triggered();
	void on_actionSavePolygon_triggered();

	void on_actionRecenter_triggered();
	void on_actionCreateInputPolygon_toggled(bool);

	void clear();

	virtual void open(QString);
  Q_SIGNALS:
	void changed();
};

MainWindow::MainWindow() : DemosMainWindow() {
	setupUi(this);
	
	this->graphicsView->setAcceptDrops(false);

	// Add a GraphicItem for the Polygon2
	pgi = new CGAL::Qt::PolygonGraphicsItem<Polygon2D>(&poly);

	QObject::connect(this, SIGNAL(changed()), pgi, SLOT(modelChanged()));

	scene.addItem(pgi);

	kgongi = new CGAL::Qt::PolygonGraphicsItem<Polygon2D>(&kgon);
	kgongi->setEdgesPen(
		QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
	kgongi->hide();
	scene.addItem(kgongi);

	lgi = new CGAL::Qt::LineGraphicsItem<Kernel>();
	lgi->setPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
	lgi->hide();
	scene.addItem(lgi);
	assert(lgi->scene() == &scene);
	// Setup input handlers. They get events before the scene gets them
	pi = new CGAL::Qt::GraphicsViewPolylineInput<Kernel>(this, &scene, 0, true);

	this->actionCreateInputPolygon->setChecked(true);
	QObject::connect(pi, SIGNAL(generate(CGAL::Object)), this,
					 SLOT(processInput(CGAL::Object)));

	//
	// Manual handling of actions
	//
	QObject::connect(this->actionQuit, SIGNAL(triggered()), this, SLOT(close()));

	//
	// Setup the scene and the view
	//
	scene.setItemIndexMethod(QGraphicsScene::NoIndex);
	scene.setSceneRect(-100, -100, 100, 100);
	this->graphicsView->setScene(&scene);
	this->graphicsView->setMouseTracking(true);

	// Uncomment the following line to get antialiasing by default.
	actionUse_Antialiasing->setChecked(true);

	// Turn the vertical axis upside down
	this->graphicsView->scale(1, -1);

	// The navigation adds zooming and translation functionality to the
	// QGraphicsView
	this->addNavigation(this->graphicsView);

	this->setupStatusBar();
	this->setupOptionsMenu();
	this->addAboutDemo(":/cgal/help/about_Polygon.html");
	this->addAboutCGAL();
	this->setupExportSVG(action_Export_SVG, graphicsView);

	this->addRecentFiles(this->menuFile, this->actionQuit);
	connect(this, SIGNAL(openRecentFile(QString)), this, SLOT(open(QString)));
}

void MainWindow::processInput(CGAL::Object o) {
	this->actionCreateInputPolygon->setChecked(false);
	std::list<Point2D> points;
	if (CGAL::assign(points, o)) {
		if ((points.size() == 1) && poly.size() > 0) {
		}
		else {
			poly.clear();
			if (points.front() == points.back()) {
				points.pop_back();
			}
			poly.insert(poly.vertices_begin(), points.begin(), points.end());
		}
		Q_EMIT(changed());
	}
}

/*
 *	Qt Automatic Connections
 *	http://doc.qt.io/qt-5/designer-using-a-ui-file.html#automatic-connections
 *
 *	setupUi(this) generates connections to the slots named
 *	"on_<action_name>_<signal_name>"
 */

void MainWindow::on_actionClear_triggered() {
	poly.clear();
	clear();
	this->actionCreateInputPolygon->setChecked(true);
	Q_EMIT(changed());
}

void MainWindow::on_actionLoadPolygon_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open Polygon File"),
													".", tr("Any file (*.*)"));
	if (!fileName.isEmpty()) {
		open(fileName);
	}
}

void MainWindow::open(QString fileName) {
	this->actionCreateInputPolygon->setChecked(false);
	std::ifstream ifs(qPrintable(fileName));
	poly.clear();
	ifs >> poly;
	clear();

	this->addToRecentFiles(fileName);
	Q_EMIT(changed());
}

void MainWindow::on_actionSavePolygon_triggered() {
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Polygon"), ".",
													tr("Any files (*.*)"));
	if (!fileName.isEmpty()) {
		std::ofstream ofs(qPrintable(fileName));
		ofs << poly;
	}
}

void MainWindow::on_actionCreateInputPolygon_toggled(bool checked) {
	poly.clear();
	clear();
	if (checked) {
		scene.installEventFilter(pi);
	}
	else {
		scene.removeEventFilter(pi);
	}
	Q_EMIT(changed());
}

void MainWindow::on_actionRecenter_triggered() {
	this->graphicsView->setSceneRect(pgi->boundingRect());
	this->graphicsView->fitInView(pgi->boundingRect(), Qt::KeepAspectRatio);
}

void MainWindow::clear() {
	lgi->hide();
	kgongi->hide();
}

#include "Polygon.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv) {
	QApplication app(argc, argv);

	app.setOrganizationDomain("geometryfactory.com");
	app.setOrganizationName("GeometryFactory");
	app.setApplicationName("Polygon_2 demo");

	// Import resources from libCGAL (Qt5).
	// See http://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
	CGAL_QT_INIT_RESOURCES;

	Q_INIT_RESOURCE(Polygon); //Use the name of the file Polygon.qrc

	MainWindow mainWindow;
	mainWindow.show();
	return app.exec();
}
