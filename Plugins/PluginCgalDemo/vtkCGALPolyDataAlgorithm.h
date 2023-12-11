/**
 * @class   vtkCGALPolyDataAlgorithm
 * @brief   virtual superclass of all polydata algorithms using CGAL.
 *
 * vtkCGALPolyDataAlgorithm is the superclass of all polydata algorithms
 * using CGAL. It defines useful methods to translate from VTK datasets
 * to CGAL surfaces back and forth.
 * This class also ensure face orientation consistency during the VTK to CGAL
 * export.
 * All filters inheriting from this class expect triangulated polydata as input,
 * expect when stated otherwise.
 */

#ifndef vtkCGALPolyDataAlgorithm_h
#define vtkCGALPolyDataAlgorithm_h

// VTK includes
#include <vtkPolyDataAlgorithm.h>
// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

using CGAL_Kernel  = CGAL::Simple_cartesian<double>;
using CGAL_Surface = CGAL::Surface_mesh<CGAL_Kernel::Point_3>;
using Graph_Verts  = boost::graph_traits<CGAL_Surface>::vertex_descriptor;
using Graph_Faces  = boost::graph_traits<CGAL_Surface>::face_descriptor;
using Graph_Coord  = boost::property_map<CGAL_Surface, CGAL::vertex_point_t>::type;
//#include "vtkCGALAlgorithmModule.h" // For export macro

// Container for CGAL related info
struct CGAL_Mesh
{
  CGAL_Surface surface;
  Graph_Coord  coords;

  CGAL_Mesh() { coords = get(CGAL::vertex_point, surface); }
};

// Filter
class vtkCGALPolyDataAlgorithm : public vtkPolyDataAlgorithm
{

public:
  static vtkCGALPolyDataAlgorithm* New();
  vtkTypeMacro(vtkCGALPolyDataAlgorithm, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //@{
  /**
   * Choose if the result mesh should have the
   * point / cell data attributes of the input.
   * If so, a vtkProbeFilter is called in order
   * to interpolate values to the new mesh.
   * Default is true
   **/
  vtkGetMacro(UpdateAttributes, bool);
  vtkSetMacro(UpdateAttributes, bool);
  vtkBooleanMacro(UpdateAttributes, bool);
  /**
 * Convert a vtkPolyData to a CGAL mesh.
 * This method fills the surface and coords
 * in the CGAL_Mesh data object
 */
  static std::unique_ptr<CGAL_Mesh> toCGAL(vtkPolyData* vtkMesh);

  /**
   * Convert a CGAL mesh to a vtk polydata.
   */
  static vtkSmartPointer<vtkPolyData> toVTK(CGAL_Mesh* cgalMesh);
  //@}

protected:
  vtkCGALPolyDataAlgorithm()           = default;
  ~vtkCGALPolyDataAlgorithm() override = default;
  /**
   * interpolate attributes of input onto the new VTK mesh
   * if UpdateAttributes is true
   */
  bool interpolateAttributes(vtkPolyData* input, vtkPolyData* vtkMesh);

  /**
   * Copy the attributes of input onto vtkMesh
   * if UpdateAttributes is true.
   */
  bool copyAttributes(vtkPolyData* input, vtkPolyData* vtkMesh);

  // Fields

  bool UpdateAttributes = true;

private:
  vtkCGALPolyDataAlgorithm(const vtkCGALPolyDataAlgorithm&) = delete;
  void operator=(const vtkCGALPolyDataAlgorithm&)           = delete;
};

#endif
