#ifndef THAEmptyImageData_H
#define THAEmptyImageData_H
#pragma execution_character_set("utf-8")
#include "vtkObject.h"
#include "vtkNew.h"
#include "vtkImageData.h"
#include "vtkInformation.h"


class EmptyImageData  : public vtkObject
{
public:
    static EmptyImageData* New();
    ~EmptyImageData();
    vtkNew<vtkImageData> mEmpty;
    vtkNew<vtkInformation> mInfomation;
    vtkImageData* f_GetEmptyImageData();
protected:
    EmptyImageData();
};

#endif // THAEmptyImageData_H
