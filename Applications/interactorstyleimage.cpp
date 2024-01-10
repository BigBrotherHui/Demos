#include "interactorstyleimage.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPointPicker.h"
#include "vtkRenderWindow.h"
#include "vtkRendererCollection.h"
#include "vtkAssemblyPath.h"
#include "vtkImageActor.h"
#include <QApplication>
#include <array>
#include <QDebug>
vtkStandardNewMacro(InteractorStyleImage);
void InteractorStyleImage::OnLeftButtonDown()
{
    if(!mIsEnable)
        return vtkInteractorStyleImage::OnLeftButtonDown();
    int *pixelPos=this->Interactor->GetEventPosition();
    this->Interactor->GetPicker()->Pick(
        pixelPos[0],pixelPos[1], 0,
        this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()
        );
    auto picker=this->Interactor->GetPicker();
    double picked[3];
    picker->GetPickPosition(picked);
    if(picker->IsA("vtkPointPicker"))
    {
        auto pointPicker=vtkPointPicker::SafeDownCast(picker);
        if(pointPicker->GetPointId()!=-1 && nullptr!=mImageData)
        {
            double *point =mImageData->GetPoint(pointPicker->GetPointId());
            double ijk[3];
            mImageData->TransformPhysicalPointToContinuousIndex(point,ijk);
            int *dim=mImageData->GetDimensions();
            bool isWidthInBound=ijk[0]>4 && ijk[0]<dim[0]-4;
            bool isHeightInBound=ijk[1]>4 && ijk[1]<dim[1]-4;
            if(isWidthInBound && isHeightInBound)
            {
                selectedPoints.push_back({ijk[0]+1,ijk[1]+1,ijk[2]});
                selectedPoints.push_back({ijk[0]+2,ijk[1]+2,ijk[2]});
                selectedPoints.push_back({ijk[0]+3,ijk[1]+3,ijk[2]});
                selectedPoints.push_back({ijk[0]+4,ijk[1]+4,ijk[2]});
                selectedPoints.push_back({ijk[0]-1,ijk[1]-1,ijk[2]});
                selectedPoints.push_back({ijk[0]-2,ijk[1]-2,ijk[2]});
                selectedPoints.push_back({ijk[0]-3,ijk[1]-3,ijk[2]});
                selectedPoints.push_back({ijk[0]-4,ijk[1]-4,ijk[2]});
                selectedPoints.push_back({ijk[0]+1,ijk[1]-1,ijk[2]});
                selectedPoints.push_back({ijk[0]+2,ijk[1]-2,ijk[2]});
                selectedPoints.push_back({ijk[0]+3,ijk[1]-3,ijk[2]});
                selectedPoints.push_back({ijk[0]+4,ijk[1]-4,ijk[2]});
                selectedPoints.push_back({ijk[0]-1,ijk[1]+1,ijk[2]});
                selectedPoints.push_back({ijk[0]-2,ijk[1]+2,ijk[2]});
                selectedPoints.push_back({ijk[0]-3,ijk[1]+3,ijk[2]});
                selectedPoints.push_back({ijk[0]-4,ijk[1]+4,ijk[2]});
            }
            selectedPoints.push_back({ijk[0],ijk[1],ijk[2]});
            mPickedPoints.push_back({ijk[0],ijk[1],ijk[2]});
            markPixel();
        }
    }
    return vtkInteractorStyleImage::OnLeftButtonDown();
}

void InteractorStyleImage::setImageData(vtkImageData *pImageData)
{
    this->mImageData=pImageData;
    if(!mUnMarkedImageData)
    {
        mUnMarkedImageData=vtkSmartPointer<vtkImageData>::New();
    }
    mUnMarkedImageData->DeepCopy(pImageData);
}

void InteractorStyleImage::setPickEnable(bool pIsEnable)
{
    mIsEnable=pIsEnable;
}

void InteractorStyleImage::markPixel()
{
    int *dim=mImageData->GetDimensions();
    for(size_t i=0;i<selectedPoints.size();i++)
    {
        std::array<double,3> tmp=selectedPoints.at(i);
        double pixelIndex[3]{tmp[0],tmp[1],tmp[2]};
        bool isWidthInBound=pixelIndex[0]>=0 && pixelIndex[0]<dim[0];
        bool isHeightInBound=pixelIndex[1]>=0 && pixelIndex[1]<dim[1];
        if(isWidthInBound && isHeightInBound)
            markPixel(pixelIndex);
    }
}

void InteractorStyleImage::clearPoints()
{
    if(!mImageData || !mUnMarkedImageData)
        return;
    selectedPoints.clear();
    mPickedPoints.clear();
    this->mImageData=mUnMarkedImageData;
    auto pointPicker=vtkPointPicker::SafeDownCast(this->Interactor->GetPicker());
    vtkAssemblyPath *path=pointPicker->GetPath();
    if(path!=nullptr){
        vtkImageActor *pickedActor =  static_cast<vtkImageActor*>(path->GetLastNode()->GetViewProp());
        pickedActor->SetInputData(mImageData);
        this->Interactor->Render();
    }
}

void InteractorStyleImage::lastPickedPoint(double *pt)
{
    int size=mPickedPoints.size();
    if(size==0 || !pt)
        return;
    std::array<double,3> ary=mPickedPoints.at(size-1);
    for(int i=0;i<3;i++)
        pt[i]=ary[i];
}

vtkImageData *InteractorStyleImage::imageData()
{
    return mImageData;
}

void InteractorStyleImage::markPixel(double *pPixelIndex)
{
    if(nullptr==mImageData)
        return;
    float *pixel=static_cast<float *>(mImageData->GetScalarPointer(pPixelIndex[0],pPixelIndex[1],pPixelIndex[2]));
    pixel[0]=4095;
    pixel[1]=0;
    pixel[2]=0;
    mImageData->Modified();
}
