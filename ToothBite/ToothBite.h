#pragma once

#include <QInputDialog>
#include <QtWidgets/QMainWindow>
#include "MyInteractorStyleTrackballActor.h"
#include "include.h"
#include "ui_ToothBite.h"

#define MAX_MODEL_NUM 10
class ToothBite : public QMainWindow {
  Q_OBJECT

 public:
  ToothBite(QWidget *parent = nullptr);
  ~ToothBite();

 private:
  Ui::ToothBiteClass ui;
  vtkSmartPointer<vtkImageActor> imgActor;
  vtkSmartPointer<vtkScalarBarActor> scalarBar;
  vtkSmartPointer<vtkAxesActor> axesActor;
  vtkSmartPointer<vtkAppendPolyData> m_UpPolyData;
  vtkSmartPointer<vtkAppendPolyData> m_LowPolyData;
  vtkSmartPointer<InteractorStyleTrackballActor> M_InteractorStyle;
  vtkSmartPointer<vtkEventQtSlotConnect> vtkQTconnect;
  vtkSmartPointer<vtkRenderer> renderer;
  vtkSmartPointer<vtkActor> m_actors[MAX_MODEL_NUM];
  vtkSmartPointer<vtkActor> intersectionActor;
  vtkSmartPointer<vtkLookupTable> seriesLut;
  vtkSmartPointer<vtkLookupTable> seriesLut2;
  vtkSmartPointer<vtkPolyDataMapper> intersectionMapper;
  vtkSmartPointer<vtkActor> intersectionActor2;
  vtkSmartPointer<vtkPolyDataMapper> intersectionMapper2;
  vtkSmartPointer<vtkAssembly> assembly1;
  vtkSmartPointer<vtkAssembly> assembly2;
  vtkSmartPointer<vtkPolyData> intersectionPolyData;
  vtkSmartPointer<vtkPolyDataMapper> m_mapper[10];
  vtkSmartPointer<vtkSTLReader> m_STLreader[10];
  int ActorNoteArray[MAX_MODEL_NUM];
  int actor_index = 0;

  bool isOkToAna = false;
  bool isOkToAna2 = false;
  bool isSimpler = false;
  int ar[11]{0};
  vtkSmartPointer<vtkActor> m_upModel, m_lowModel;
  vtkSmartPointer<vtkPolyDataMapper> m_upModelMapper, m_lowModelMapper;
  vtkSmartPointer<vtkSTLReader> m_upModelReader, m_lowModelReader;
  vtkSmartPointer<vtkPolyData> m_upModelPoly{nullptr}, m_lowModelPoly{nullptr};

  // functions
  void fun_ReadSTLFile(void);      // read file
  void fun_ClearSTLFile(void);     // clear files
  void fun_ReadAnaFile(void);      // read files for ana
  void fun_ImpactAnalysis(void);   // analysis impact situation
  void fun_ImpactAnalysis2(void);  // analysis impact situation
  void fun_ChangeBackgroudColor(
      QColor color);                // change show window background color
  void fun_ChangeLight(void);       // change light
  void fun_ChangeModelColor(void);  // change color of models
  void fun_ResetWidget(void);       // reset widget
  void fun_bind(void);
  void fun_Fillholes(vtkSmartPointer<vtkPolyData> pPolydata, bool pIsUp);
  void fun_GetUpPolyData();
  void fun_GetLowPolyData();

 public slots:
  void MykeyPressEvent();
 private slots:
  void on_pushButton_clicked();
};
