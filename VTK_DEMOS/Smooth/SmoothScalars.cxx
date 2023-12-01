#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeProperty.h>
//#include <vtkVolumeRayCastCompositeFunction.h>
//#include <vtkVolumeRayCastMIPFunction.h>
//#include <vtkVolumeRayCastMapper.h>
#include <vtkVolume.h>
//#include <vtkVolumeTextureMapper2D.h>
#include <vtkImageReader.h>
#include <vtkImageImport.h>
#include <vtkImageCast.h>
//#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include "vtkPointData.h"
#include <vtkPointPicker.h>
#include <vtkCommand.h>
#include <vtkMarchingCubes.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataMapper.h>
#include <vtkCamera.h>
#include <vtkMarchingCubes.h>
#include <vtkVectorNorm.h>
#include <vtkDataSetMapper.h>
#include <vtkImageToPolyDataFilter.h>
#include <vtkPolyDataReader.h>
#include <vtkLookupTable.h>
#include <vtkSphereSource.h>
#include <vtkCallbackCommand.h>
#include <vtkProperty.h>
#include <vtkImagePlaneWidget.h>
#include <vtkImageActor.h>

#include <vtkImageMapToColors.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleFlight.h>
#include <vtkPropPicker.h>
#include <vtkBoxWidget.h>
#include <vtkPlanes.h>
#include <vtkRendererCollection.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkImageThreshold.h>
#include <vtkWindowToImageFilter.h>
//#include <vtkAVIWriter.h>
#include <vtkSmartPointer.h>
#include <vtkOBJExporter.h>
#include <vtkMaskPoints.h>
#include <vtkConeSource.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>
#include <vtkThreshold.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkCellLocator.h>
#include <vtkCellArray.h>
#include "vtkGenericCell.h"
#include "vtkOrientationMarkerWidget.h"
#include "vtkAxesActor.h"
#include "vtkAnnotatedCubeActor.h"
#include "vtkPropAssembly.h"
#include "vtkPlaneWidget.h"
#include "vtkCellPicker.h"
#include "vtkLandmarkTransform.h"
#include <vtkMatrix4x4.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkTriangleFilter.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>


/*
* Please change this distance to suit your application 
* You can obtain this from the bounding box of your mesh when it is loaded on Paraview
* Calculate the diagonal of the bounding box using the Euclidean distance between 3D points 
*/
#define MAXDIST 100.0

int main( int argc, char * argv[] )
{
  // Input file
  char * input_file = "D:/Demos/stl/out.stl";
  // Output file
  char * output_file = "D:/Demos/stl/out1.stl";
  // Iterations
  vtkIdType itr = 20;

 // int probe = atoi ( argv[4] );
  // Reader
  vtkSTLReader* reader = vtkSTLReader::New();
  reader->SetFileName( input_file );
  reader->Update();
  // Triangle Filter
  vtkTriangleFilter * triangle = vtkTriangleFilter::New();
  triangle->PassVertsOn();
  triangle->PassLinesOn();
  triangle->SetInputConnection( reader->GetOutputPort() );
  triangle->Update();

  vtkPolyData * tree        = triangle->GetOutput();
  vtkPoints   * tree_points = tree->GetPoints();
  vtkIdType     num_points  = tree_points->GetNumberOfPoints();
  vtkIdList   * point_cells = vtkIdList::New();
  vtkIdList   * cell_points = vtkIdList::New();
  vtkIdType point_id;

    vtkPolyData * tree_out = vtkPolyData::New();
    tree_out->DeepCopy( tree );

  // This will hold a temporary copy of the points.
   vtkDoubleArray *tmp_DA = vtkDoubleArray::New();
   tmp_DA->SetNumberOfTuples(num_points);
   tmp_DA->SetNumberOfComponents(1);

  // Iterations
  tree->BuildLinks();
  vtkIdType i;
  for ( i = 0; i < itr; ++i ) {
    // Try to smooth the point scalars.
    for ( point_id = 0; point_id < num_points; ++point_id ) {
        double point_scalar = tree_out->GetPointData()->GetScalars()->GetTuple1(point_id);
        double *point_position = new double[3];
        point_position = tree->GetPoint(point_id);

        double average_scalar = point_scalar;
        double den_fact = 1;

        tree->GetPointCells( point_id, point_cells );
        vtkIdType num_cells = point_cells->GetNumberOfIds();

        vtkIdType cell;
        for (cell = 0; (cell < num_cells); ++cell ) {
            vtkIdType neighbor_cell_id;
            neighbor_cell_id = point_cells->GetId( cell );
            tree->GetCellPoints( neighbor_cell_id, cell_points );
            vtkIdType num_cell_points = cell_points->GetNumberOfIds();
            // Loop over the cell points.
            vtkIdType neighbor_point;
            for (neighbor_point = 0;neighbor_point < num_cell_points;++neighbor_point ) {
                // Get the neighbor point id.
                vtkIdType neighbor_point_id = cell_points->GetId(neighbor_point);
                double *neighbor_position = new double[3];
                neighbor_position = tree->GetPoint(neighbor_point_id);

                double a = neighbor_position[0] - point_position[0];
                double b = neighbor_position[1] - point_position[1];
                double c = neighbor_position[2] - point_position[2];
                double distance = sqrt(a*a + b*b + c*c);
                double weight = (1-distance/MAXDIST) * int(distance<MAXDIST);

                // Get the neighbor point position.
                double neighbor_scalar;
                neighbor_scalar = tree_out->GetPointData()->GetScalars()->GetTuple1(neighbor_point_id);
                // Add it to the average.
                average_scalar += neighbor_scalar * weight;
                den_fact += weight;
                } // End neighbor points loop.
            } // End of cells loop
        average_scalar /= den_fact;
        tmp_DA->SetTuple1(point_id, average_scalar );
        } // End of smooth loop.
//    tree_points->DeepCopy( tmp_points );
    tree_out->GetPointData()->SetScalars(tmp_DA);
    //tree_out->Update();
    } // End of itr loop

  // Eventually Write down the new polydata
  vtkSTLWriter* writer = vtkSTLWriter::New();
  writer->SetFileName(output_file);
  writer->SetInputData(reader->GetOutput());
  writer->Write();
  std::cout << "*************" << std::endl;
  return 0;
}
