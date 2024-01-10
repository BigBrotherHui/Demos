#ifndef INTERACTORSTYLEIMAGE_H
#define INTERACTORSTYLEIMAGE_H

#include "vtkInteractorStyleImage.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include <vector>
#include <QVTKOpenGLNativeWidget.h>
class InteractorStyleImage : public vtkInteractorStyleImage
{
public:
    static InteractorStyleImage* New();
    vtkTypeMacro(InteractorStyleImage, vtkInteractorStyleImage);
    virtual void OnLeftButtonDown() override;
    void setImageData(vtkImageData *pImageData);
    void setPickEnable(bool pIsEnable=true);
    void markPixel();
    void clearPoints();
    void lastPickedPoint(double *pt);
    vtkImageData *imageData();
protected:
    void markPixel(double *pPixelIndex);
private:
    std::vector<std::array<double,3>> selectedPoints;
    std::vector<std::array<double,3>> mPickedPoints;
    vtkSmartPointer<vtkImageData> mImageData{nullptr};
    vtkSmartPointer<vtkImageData> mUnMarkedImageData{nullptr};
    bool mIsEnable{false};
};

#endif // INTERACTORSTYLEIMAGE_H
