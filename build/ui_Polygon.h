/********************************************************************************
** Form generated from reading UI file 'Polygon.ui'
**
** Created by: Qt User Interface Compiler version 5.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_POLYGON_H
#define UI_POLYGON_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_PolygonWindow
{
public:
    QAction *actionAbout;
    QAction *actionAboutCGAL;
    QAction *actionQuit;
    QAction *actionClear;
    QAction *actionLoadPolygon;
    QAction *actionSavePolygon;
    QAction *actionRecenter;
    QAction *actionCreateInputPolygon;
    QAction *action_Export_SVG;
    QWidget *centralwidget;
    QHBoxLayout *hboxLayout;
    QSplitter *splitter;
    QGraphicsView *graphicsView;
    QStatusBar *statusbar;
    QToolBar *fileToolBar;
    QToolBar *toolBar;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuTools;

    void setupUi(QMainWindow *PolygonWindow)
    {
        if (PolygonWindow->objectName().isEmpty())
            PolygonWindow->setObjectName(QStringLiteral("PolygonWindow"));
        PolygonWindow->resize(568, 325);
        QIcon icon;
        icon.addFile(QStringLiteral(":/cgal/logos/cgal_icon"), QSize(), QIcon::Normal, QIcon::Off);
        PolygonWindow->setWindowIcon(icon);
        actionAbout = new QAction(PolygonWindow);
        actionAbout->setObjectName(QStringLiteral("actionAbout"));
        actionAboutCGAL = new QAction(PolygonWindow);
        actionAboutCGAL->setObjectName(QStringLiteral("actionAboutCGAL"));
        actionQuit = new QAction(PolygonWindow);
        actionQuit->setObjectName(QStringLiteral("actionQuit"));
        actionClear = new QAction(PolygonWindow);
        actionClear->setObjectName(QStringLiteral("actionClear"));
        QIcon icon1;
        icon1.addFile(QStringLiteral(":/cgal/fileToolbar/fileNew.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionClear->setIcon(icon1);
        actionLoadPolygon = new QAction(PolygonWindow);
        actionLoadPolygon->setObjectName(QStringLiteral("actionLoadPolygon"));
        QIcon icon2;
        icon2.addFile(QStringLiteral(":/cgal/fileToolbar/fileOpen.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLoadPolygon->setIcon(icon2);
        actionSavePolygon = new QAction(PolygonWindow);
        actionSavePolygon->setObjectName(QStringLiteral("actionSavePolygon"));
        QIcon icon3;
        icon3.addFile(QStringLiteral(":/cgal/fileToolbar/fileSave.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSavePolygon->setIcon(icon3);
        actionRecenter = new QAction(PolygonWindow);
        actionRecenter->setObjectName(QStringLiteral("actionRecenter"));
        QIcon icon4;
        icon4.addFile(QStringLiteral(":/cgal/Input/zoom-best-fit"), QSize(), QIcon::Normal, QIcon::Off);
        actionRecenter->setIcon(icon4);
        actionCreateInputPolygon = new QAction(PolygonWindow);
        actionCreateInputPolygon->setObjectName(QStringLiteral("actionCreateInputPolygon"));
        actionCreateInputPolygon->setCheckable(true);
        action_Export_SVG = new QAction(PolygonWindow);
        action_Export_SVG->setObjectName(QStringLiteral("action_Export_SVG"));
        centralwidget = new QWidget(PolygonWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        hboxLayout = new QHBoxLayout(centralwidget);
        hboxLayout->setObjectName(QStringLiteral("hboxLayout"));
        splitter = new QSplitter(centralwidget);
        splitter->setObjectName(QStringLiteral("splitter"));
        splitter->setOrientation(Qt::Horizontal);
        graphicsView = new QGraphicsView(splitter);
        graphicsView->setObjectName(QStringLiteral("graphicsView"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(2);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(graphicsView->sizePolicy().hasHeightForWidth());
        graphicsView->setSizePolicy(sizePolicy);
        graphicsView->setFocusPolicy(Qt::StrongFocus);
        graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
        graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
        graphicsView->setTransformationAnchor(QGraphicsView::NoAnchor);
        splitter->addWidget(graphicsView);

        hboxLayout->addWidget(splitter);

        PolygonWindow->setCentralWidget(centralwidget);
        statusbar = new QStatusBar(PolygonWindow);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        PolygonWindow->setStatusBar(statusbar);
        fileToolBar = new QToolBar(PolygonWindow);
        fileToolBar->setObjectName(QStringLiteral("fileToolBar"));
        PolygonWindow->addToolBar(Qt::TopToolBarArea, fileToolBar);
        toolBar = new QToolBar(PolygonWindow);
        toolBar->setObjectName(QStringLiteral("toolBar"));
        PolygonWindow->addToolBar(Qt::TopToolBarArea, toolBar);
        menubar = new QMenuBar(PolygonWindow);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 568, 21));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        menuTools = new QMenu(menubar);
        menuTools->setObjectName(QStringLiteral("menuTools"));
        PolygonWindow->setMenuBar(menubar);

        fileToolBar->addAction(actionClear);
        fileToolBar->addAction(actionLoadPolygon);
        fileToolBar->addAction(actionSavePolygon);
        toolBar->addAction(actionRecenter);
        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuTools->menuAction());
        menuFile->addSeparator();
        menuFile->addAction(actionClear);
        menuFile->addAction(actionLoadPolygon);
        menuFile->addAction(actionSavePolygon);
        menuFile->addAction(action_Export_SVG);
        menuFile->addSeparator();
        menuFile->addAction(actionQuit);
        menuTools->addSeparator();
        menuTools->addAction(actionRecenter);

        retranslateUi(PolygonWindow);

        QMetaObject::connectSlotsByName(PolygonWindow);
    } // setupUi

    void retranslateUi(QMainWindow *PolygonWindow)
    {
        PolygonWindow->setWindowTitle(QApplication::translate("PolygonWindow", "CGAL 2D Polygon", 0));
        actionAbout->setText(QApplication::translate("PolygonWindow", "&About", 0));
        actionAboutCGAL->setText(QApplication::translate("PolygonWindow", "About &CGAL", 0));
        actionQuit->setText(QApplication::translate("PolygonWindow", "&Quit", 0));
        actionQuit->setShortcut(QApplication::translate("PolygonWindow", "Ctrl+Q", 0));
        actionClear->setText(QApplication::translate("PolygonWindow", "&Clear", 0));
        actionClear->setShortcut(QApplication::translate("PolygonWindow", "Ctrl+C", 0));
        actionLoadPolygon->setText(QApplication::translate("PolygonWindow", "&Load Polygon", 0));
        actionLoadPolygon->setShortcut(QApplication::translate("PolygonWindow", "Ctrl+L", 0));
        actionSavePolygon->setText(QApplication::translate("PolygonWindow", "&Save Polygon", 0));
        actionSavePolygon->setShortcut(QApplication::translate("PolygonWindow", "Ctrl+S", 0));
        actionRecenter->setText(QApplication::translate("PolygonWindow", "Re&center the viewport", 0));
        actionRecenter->setShortcut(QApplication::translate("PolygonWindow", "Ctrl+R", 0));
        actionCreateInputPolygon->setText(QApplication::translate("PolygonWindow", "Create Input Polygon", 0));
        action_Export_SVG->setText(QApplication::translate("PolygonWindow", "&Export SVG...", 0));
        fileToolBar->setWindowTitle(QApplication::translate("PolygonWindow", "File Tools", 0));
        toolBar->setWindowTitle(QApplication::translate("PolygonWindow", "Visualization Tools", 0));
        menuFile->setTitle(QApplication::translate("PolygonWindow", "&File", 0));
        menuTools->setTitle(QApplication::translate("PolygonWindow", "&Algorithms", 0));
    } // retranslateUi

};

namespace Ui {
    class PolygonWindow: public Ui_PolygonWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_POLYGON_H
