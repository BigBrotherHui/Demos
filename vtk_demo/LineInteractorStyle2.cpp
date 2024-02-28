
#include "LineInteractorStyle2.h"
#include <vtkImageActor.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkImageCast.h>
#include <vtkImageShiftScale.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <itkImageToVTKImageFilter.h>
#include <itkVTKImageToImageFilter.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkImageData.h>
#include <itkConnectedThresholdImageFilter.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkImageMapper3D.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformation.h>
#include <vtkTransform.h>
#include <vtkImageFlip.h>
#include <QDebug>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkPointPicker.h>
vtkStandardNewMacro(LineInteractorStyle2);

LineInteractorStyle2::LineInteractorStyle2()
{

}

LineInteractorStyle2::~LineInteractorStyle2()
{
}

void LineInteractorStyle2::OnMouseMove()
{
	int *screenPos=this->Interactor->GetEventPosition();
	vtkRenderer* renderer = this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
	vtkPointPicker *picker=vtkPointPicker::SafeDownCast(this->Interactor->GetPicker());
	if(picker)
	{
		picker->Pick(screenPos[0], screenPos[1],0,renderer);
		double *worldPos=picker->GetPickPosition();
		this->actor->SetPosition(worldPos[0], worldPos[1], 0);
		this->Interactor->Render();
	}
	return Superclass::OnMouseMove();
}
