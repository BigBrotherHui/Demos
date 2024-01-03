/**
 * @class   vtkCGALBooleanOperation
 * @brief   Performs a boolean operation between two vtkPolyData objects
 *
 * vtkCGALBooleanOperation is a filter allowing to perform a boolean operation
 * between two closed, triangulated polygonal meshes.
 * These operations include union, intersection, and difference.
 * The resulting mesh is closed.
 */

#ifndef vtkCGALBooleanOperation_h
#define vtkCGALBooleanOperation_h
#include "vtkCGALPolyDataAlgorithm.h"

//#include "vtkCGALPMPModule.h" // For export macro

class vtkCGALBooleanOperation : public vtkCGALPolyDataAlgorithm
{
public:
  static vtkCGALBooleanOperation* New();
  vtkTypeMacro(vtkCGALBooleanOperation, vtkCGALPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * List of available boolean operations between the volume bounded by the input
   * and the volume bounded by the source. The result is a closed mesh.
   **/
  enum Operation___
  {
    Operation___DIFFERENCE = 0,
    Operation___INTERSECTION,
    Operation___UNION
  };

  ///@{
  /**
   * Get/set the type of boolean operation.
   * Default is DIFFERENCE.
   **/
  vtkGetMacro(OperationType, int);
  vtkSetClampMacro(
    OperationType, int, vtkCGALBooleanOperation::Operation___DIFFERENCE, vtkCGALBooleanOperation::Operation___UNION);
  ///@}

  /**
   * Set input connection for the second vtkPolyData.
   **/
  void SetSourceConnection(vtkAlgorithmOutput* algOutput);

protected:
  vtkCGALBooleanOperation();
  ~vtkCGALBooleanOperation() override = default;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  int OperationType = vtkCGALBooleanOperation::Operation___DIFFERENCE;

private:
  vtkCGALBooleanOperation(const vtkCGALBooleanOperation&) = delete;
  void operator=(const vtkCGALBooleanOperation&)          = delete;
};

#endif
