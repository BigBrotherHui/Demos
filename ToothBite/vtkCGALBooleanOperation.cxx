#include "vtkCGALBooleanOperation.h"

// VTK related includes
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

// CGAL related includes
#include <CGAL/Polygon_mesh_processing/corefinement.h>

vtkStandardNewMacro(vtkCGALBooleanOperation);

namespace pmp = CGAL::Polygon_mesh_processing;

//------------------------------------------------------------------------------
vtkCGALBooleanOperation::vtkCGALBooleanOperation()
{
  this->SetNumberOfInputPorts(2);
}

//------------------------------------------------------------------------------
void vtkCGALBooleanOperation::PrintSelf(ostream& os, vtkIndent indent)
{
  switch (this->OperationType)
  {
    case vtkCGALBooleanOperation::Operation___DIFFERENCE:
      os << indent << "OperationType: "
         << "Difference" << std::endl;
      break;
    case vtkCGALBooleanOperation::Operation___INTERSECTION:
      os << indent << "OperationType: "
         << "Intersection" << std::endl;
      break;
    case vtkCGALBooleanOperation::Operation___UNION:
      os << indent << "OperationType: "
         << "Union" << std::endl;
      break;
    default:
      os << indent << "OperationType: "
         << "Unknown" << std::endl;
      break;
  }

  this->Superclass::PrintSelf(os, indent);
}

//------------------------------------------------------------------------------
void vtkCGALBooleanOperation::SetSourceConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}

//------------------------------------------------------------------------------
int vtkCGALBooleanOperation::RequestData(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // Get the input and source data object
  vtkPolyData* inputData  = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData* sourceData = vtkPolyData::GetData(inputVector[1]);
  vtkPolyData* output     = vtkPolyData::GetData(outputVector);

  if (!inputData || !sourceData)
  {
    vtkErrorMacro("Missing input or source.");
  }

  // Create the surface meshes for CGAL
  // ----------------------------------

  std::unique_ptr<CGAL_Mesh> cgalInputMesh  = this->toCGAL(inputData);
  std::unique_ptr<CGAL_Mesh> cgalSourceMesh = this->toCGAL(sourceData);
  std::unique_ptr<CGAL_Mesh> cgalOutMesh    = this->toCGAL(output);

  // CGAL Processing
  // ---------------

  bool res = true;
  try
  {
    // Preprocess
    if (!CGAL::Polygon_mesh_processing::does_bound_a_volume(cgalInputMesh->surface))
    {
      pmp::orient_to_bound_a_volume(cgalInputMesh->surface);
    }
    if (!CGAL::Polygon_mesh_processing::does_bound_a_volume(cgalSourceMesh->surface))
    {
      pmp::orient_to_bound_a_volume(cgalSourceMesh->surface);
    }

    // Main process
    switch (this->OperationType)
    {
      case vtkCGALBooleanOperation::Operation___DIFFERENCE:
        res = pmp::corefine_and_compute_difference(cgalInputMesh->surface, cgalSourceMesh->surface,
          cgalOutMesh->surface, pmp::parameters::all_default(), pmp::parameters::all_default());
        break;
      case vtkCGALBooleanOperation::Operation___INTERSECTION:
        res =
          pmp::corefine_and_compute_intersection(cgalInputMesh->surface, cgalSourceMesh->surface,
            cgalOutMesh->surface, pmp::parameters::all_default(), pmp::parameters::all_default());
        break;
      case vtkCGALBooleanOperation::Operation___UNION:
        res = pmp::corefine_and_compute_union(cgalInputMesh->surface, cgalSourceMesh->surface,
          cgalOutMesh->surface, pmp::parameters::all_default(), pmp::parameters::all_default());
        break;
      default:
        vtkErrorMacro("Unknown boolean operation!");
        res = false;
        break;
    }
  }
  catch (std::exception& e)
  {
    vtkErrorMacro("CGAL Exception: " << e.what());
    return 0;
  }

  if (!res)
  {
    // help user know the issue with their data.
    // Most of these checks are done after the processing for performance reasons
    std::cerr << "Boolean operation failed. Checking precondition:" << std::endl;
    bool se1 = CGAL::Polygon_mesh_processing::does_self_intersect(cgalInputMesh->surface);
    bool se2 = CGAL::Polygon_mesh_processing::does_self_intersect(cgalSourceMesh->surface);
    std::cerr << "Input self intersect: " << se1 << std::endl;
    std::cerr << "Source self intersect: " << se2 << std::endl;
    bool bv1 = CGAL::Polygon_mesh_processing::does_bound_a_volume(cgalInputMesh->surface);
    bool bv2 = CGAL::Polygon_mesh_processing::does_bound_a_volume(cgalSourceMesh->surface);
    std::cerr << "Input bounds a volume: " << bv1 << std::endl;
    std::cerr << "Source bounds a volume: " << bv2 << std::endl;

    return 0;
  }

  // VTK Output
  // ----------

  output->ShallowCopy(this->toVTK(cgalOutMesh.get()));

  this->interpolateAttributes(inputData, output);

  return 1;
}
