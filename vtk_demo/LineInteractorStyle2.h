#pragma once
#include <vtkInteractorStyleImage.h>
#include <itkImage.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkImageActor.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkImageActor.h>
class LineInteractorStyle2 : public vtkInteractorStyleImage
{
public:
	static LineInteractorStyle2* New();
	vtkTypeMacro(LineInteractorStyle2, vtkInteractorStyleImage);

	LineInteractorStyle2();
	~LineInteractorStyle2();
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkImageActor> actorImage;
protected:
	void OnMouseMove() override;
private:

};
