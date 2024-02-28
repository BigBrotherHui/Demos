#include "PluginExample.h"
#include "widget.h"
#include <iostream>
#include "mesh_processing.h"
#include "ProjectionWidget.h"

#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkAppendFilter.h>
#include <vtkProperty.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellDataToPointData.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkContourFilter.h>
#include <vtkAutoInit.h>
#include <vtkStripper.h>
#include <vtkCellArrayIterator.h>
#include <cstdlib>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkLabeledDataMapper.h>
#include <vtkPolyData.h>
#include <vtkTextProperty.h>
#include <vtkActor2D.h>
#include <vtkPolyDataMapper.h>
#include <vtkDelaunay3D.h>
#include <vtkCubeSource.h>
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);
int PluginExample::typeId = qRegisterMetaType<PluginExample*>("PluginExample");
PluginExample::PluginExample()
{

}

PluginExample::~PluginExample()
{
}

QWidget* PluginExample::createWidget()
{


        //vtkNew<vtkUnstructuredGrid> grid;   /* 创建一个非结构化网格对象 */
        ///*=============================网格节点=============================*/
        //vtkNew<vtkPoints> meshPoints; /* 用于存储网格节点 */
        //meshPoints->InsertNextPoint(20, 0, 0);
        //meshPoints->InsertNextPoint(17.320508109792, 10.0000000321997, 0);
        //meshPoints->InsertNextPoint(10.0000000084796, 17.3205081212641, 0);
        //meshPoints->InsertNextPoint(0, 20, 0);
        //meshPoints->InsertNextPoint(-10.0000000245946, 17.3205081142845, 0);
        //meshPoints->InsertNextPoint(-17.320508106976, 10.0000000330863, 0);
        //meshPoints->InsertNextPoint(-20, 0, 0);
        //meshPoints->InsertNextPoint(-17.3205081361632, -9.99999998652343, 0);
        //meshPoints->InsertNextPoint(-10.0000000640139, -17.3205080893384, 0);
        //meshPoints->InsertNextPoint(0, -20, 0);
        //meshPoints->InsertNextPoint(9.99999997546481, -17.3205081119126, 0);
        //meshPoints->InsertNextPoint(17.3205080670077, -10.000000025522, 0);
        //meshPoints->InsertNextPoint(-7.15966651257561, -7.12683602128246, 0);
        //meshPoints->InsertNextPoint(-2.63739358105485, -12.9978938464381, 0);
        //meshPoints->InsertNextPoint(-9.10592684278358, -0.737849316405346, 0);
        //meshPoints->InsertNextPoint(6.59947131139743, -9.08284206416132, 0);
        //meshPoints->InsertNextPoint(-0.137734590259206, -2.71811814272003, 0);
        //meshPoints->InsertNextPoint(-7.75933677195531, 5.80947235781205, 0);
        //meshPoints->InsertNextPoint(9.11726442880688, -1.21602962482667, 0);
        //meshPoints->InsertNextPoint(-0.484322809037511, 4.65580556604978, 0);
        //meshPoints->InsertNextPoint(-3.48354780828141, 10.4142340529639, 0);
        //meshPoints->InsertNextPoint(7.28614335171517, 6.31993800133739, 0);
        //meshPoints->InsertNextPoint(3.36376652373578, 13.053010718073, 0);
        //grid->SetPoints(meshPoints);
        ///*=================================================================*/

        ///*=============================设置单元=============================*/
        //vtkNew<vtkIdList> pointsId;
        //pointsId->SetNumberOfIds(3);
        //pointsId->SetId(0, 9);
        //pointsId->SetId(1, 13);
        //pointsId->SetId(2, 8);
        //grid->InsertNextCell(VTK_TRIANGLE, 3, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(3);
        //pointsId->SetId(0, 14);
        //pointsId->SetId(1, 12);
        //pointsId->SetId(2, 16);
        //grid->InsertNextCell(VTK_TRIANGLE, 3, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(3);
        //pointsId->SetId(0, 11);
        //pointsId->SetId(1, 15);
        //pointsId->SetId(2, 10);
        //grid->InsertNextCell(VTK_TRIANGLE, 3, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(3);
        //pointsId->SetId(0, 15);
        //pointsId->SetId(1, 18);
        //pointsId->SetId(2, 16);
        //grid->InsertNextCell(VTK_TRIANGLE, 3, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(3);
        //pointsId->SetId(0, 17);
        //pointsId->SetId(1, 19);
        //pointsId->SetId(2, 20);
        //grid->InsertNextCell(VTK_TRIANGLE, 3, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(3);
        //pointsId->SetId(0, 3);
        //pointsId->SetId(1, 22);
        //pointsId->SetId(2, 2);
        //grid->InsertNextCell(VTK_TRIANGLE, 3, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 8);
        //pointsId->SetId(1, 13);
        //pointsId->SetId(2, 12);
        //pointsId->SetId(3, 7);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 12);
        //pointsId->SetId(1, 14);
        //pointsId->SetId(2, 6);
        //pointsId->SetId(3, 7);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 13);
        //pointsId->SetId(1, 9);
        //pointsId->SetId(2, 10);
        //pointsId->SetId(3, 15);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 6);
        //pointsId->SetId(1, 14);
        //pointsId->SetId(2, 17);
        //pointsId->SetId(3, 5);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 16);
        //pointsId->SetId(1, 12);
        //pointsId->SetId(2, 13);
        //pointsId->SetId(3, 15);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 14);
        //pointsId->SetId(1, 16);
        //pointsId->SetId(2, 19);
        //pointsId->SetId(3, 17);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 17);
        //pointsId->SetId(1, 20);
        //pointsId->SetId(2, 4);
        //pointsId->SetId(3, 5);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 16);
        //pointsId->SetId(1, 18);
        //pointsId->SetId(2, 21);
        //pointsId->SetId(3, 19);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 15);
        //pointsId->SetId(1, 11);
        //pointsId->SetId(2, 0);
        //pointsId->SetId(3, 18);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 20);
        //pointsId->SetId(1, 19);
        //pointsId->SetId(2, 21);
        //pointsId->SetId(3, 22);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 22);
        //pointsId->SetId(1, 3);
        //pointsId->SetId(2, 4);
        //pointsId->SetId(3, 20);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 21);
        //pointsId->SetId(1, 18);
        //pointsId->SetId(2, 0);
        //pointsId->SetId(3, 1);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));

        //pointsId->Initialize();
        //pointsId->SetNumberOfIds(4);
        //pointsId->SetId(0, 2);
        //pointsId->SetId(1, 22);
        //pointsId->SetId(2, 21);
        //pointsId->SetId(3, 1);
        //grid->InsertNextCell(VTK_QUAD, 4, pointsId->GetPointer(0));
        ///*=================================================================*/

        ///*==============================设置单元颜色========================*/
        //vtkNew<vtkFloatArray> color;                        /* 用于设置每个单元的随机颜色 */
        //for (int i = 0; i < grid->GetNumberOfCells(); ++i)
        //{
        //    color->InsertNextValue(/*std::rand()*/i);
        //}
        //grid->GetCellData()->SetScalars(color);             /* 设置单元格的颜色 */
        ///*=================================================================*/

        ///*==============================单元——>节点========================*/
        //double rangeMin = grid->GetScalarRange()[0];
        //double rangeMax = grid->GetScalarRange()[1];

        //vtkNew<vtkCellDataToPointData> theCellDataToPointData;
        //theCellDataToPointData->SetInputData(grid);
        //theCellDataToPointData->PassCellDataOn();
        //theCellDataToPointData->Update();
        ///*=================================================================*/

        ///*============================等值线================================*/
        //vtkNew<vtkContourFilter> theContourFilter;
        //theContourFilter->SetInputData(theCellDataToPointData->GetOutput());
        //theContourFilter->GenerateValues(5, rangeMin, rangeMax);
        //theContourFilter->Update();
        ///*=================================================================*/

        ///*============================等值线值标记===========================*/
        //// 用于将离散的三角面片拼接为连续的等值面
        //vtkNew<vtkStripper> theStripper;
        //theStripper->SetInputData(theContourFilter->GetOutput());
        //theStripper->Update();
        //// 等值线数目
        //vtkIdType lines = theStripper->GetOutput()->GetNumberOfLines();
        //vtkPoints* points = theStripper->GetOutput()->GetPoints();
        //vtkCellArray* cells = theStripper->GetOutput()->GetLines();
        //vtkDataArray* scalars = theStripper->GetOutput()->GetPointData()->GetScalars();
        //vtkNew<vtkPolyData> labelPolyData;
        //vtkNew<vtkPoints> labelPoints;
        //vtkNew<vtkDoubleArray> labelScalars;
        //labelScalars->SetNumberOfComponents(1);
        //labelScalars->SetName("IsoValues");

        //// 创建单元数组迭代器
        //auto cellIter = vtk::TakeSmartPointer(cells->NewIterator());
        //for (cellIter->GoToFirstCell(); !cellIter->IsDoneWithTraversal(); cellIter->GoToNextCell())
        //{
        //    vtkIdList* cell = cellIter->GetCurrentCell();
        //    const vtkIdType samplePtIdx = static_cast<vtkIdType>(vtkMath::Random(0, cell->GetNumberOfIds()));
        //    vtkIdType midPointId = cell->GetId(samplePtIdx);

        //    double midPoint[3];
        //    points->GetPoint(midPointId, midPoint);
        //    labelPoints->InsertNextPoint(midPoint);
        //    labelScalars->InsertNextTuple1(scalars->GetTuple1(midPointId));
        //}
        //labelPolyData->SetPoints(labelPoints);

        //labelPolyData->GetPointData()->SetScalars(labelScalars);
        ///*=================================================================*/

        //vtkNew<vtkDataSetMapper> contourMapper;
        //contourMapper->SetInputData(theContourFilter->GetOutput());
        //vtkNew<vtkDataSetMapper> mapper;
        //mapper->SetInputData(theCellDataToPointData->GetOutput());
        //mapper->SetScalarModeToUsePointData();
        //mapper->SetScalarRange(rangeMin, rangeMax);

        //vtkNew<vtkActor> contourActor;
        //contourActor->SetMapper(contourMapper);

		vtkNew<vtkPolyData> pp;
		vtkNew<vtkCubeSource> cube;
        cube->Update();
        vtkNew<vtkPoints> pts;
        pts->DeepCopy(cube->GetOutput()->GetPoints());
        pp->SetPoints(pts);
        vtkNew<vtkDelaunay3D> dela;
        dela->SetInputData(pp);
        dela->Update();
        vtkNew<vtkActor> actor;
        //vtkNew<vtkPolyDataMapper> mapper;
        vtkNew<vtkDataSetMapper> mapper;
        mapper->SetInputData(dela->GetOutput());
        actor->SetMapper(mapper);

        //vtkNew<vtkLabeledDataMapper> labelMapper;
        //labelMapper->SetFieldDataName("IsoValues");
        //labelMapper->SetInputData(labelPolyData);
        //labelMapper->SetLabelModeToLabelScalars();
        //labelMapper->SetLabelFormat("%6.2f");
        //labelMapper->GetLabelTextProperty()->SetColor(0, 0, 1);

        //vtkNew<vtkActor2D> labelActor;
        //labelActor->SetMapper(labelMapper);

        vtkNew<vtkRenderer> renderer;
        renderer->AddActor(actor);
        ////renderer->AddActor(labelActor);
        //renderer->AddActor(contourActor);

        vtkNew<vtkRenderWindowInteractor> windowInteractor;
        vtkNew<vtkRenderWindow> renderWindow;

        windowInteractor->SetRenderWindow(renderWindow);
        renderWindow->AddRenderer(renderer);


        renderWindow->Render();

        windowInteractor->Start();

	if(!m_widget)
		m_widget = new Widget;
	/*Eigen::MatrixX3d points;
	points.resize(8, 3);
	points.row(0) << 1, 0, 0;
	points.row(1) << 11, 0, 0;
	points.row(2) << 1, 10, 0;
	points.row(3) << 11, 10, 0;
	points.row(4) << 1, 0, 10;
	points.row(5) << 11, 0, 10;
	points.row(6) << 1, 10, 10;
	points.row(7) << 11, 10, 10;
	static_cast<ProjectionWidget*>(m_widget)->drawRect(points);*/
	return m_widget;
}
