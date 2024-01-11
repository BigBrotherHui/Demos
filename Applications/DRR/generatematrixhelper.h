#ifndef GENERATEMATRIXHELPER_H
#define GENERATEMATRIXHELPER_H
#include "itkTransform.h"
#include "itkVector.h"
#include "itkEuler3DTransform.h"
class vtkMatrix4x4;
class GenerateMatrixHelper
{
public:
    using TransformType = itk::Euler3DTransform<float>;
    using TransformPointer = typename TransformType::Pointer;
    GenerateMatrixHelper();
    ~GenerateMatrixHelper();
    void getTransformMatrix(bool isFront,float *);
    void setCenterOffset(float *centeroffset);
    void setFocalDistance(float focaldistance);
    void setImageSize(int *imagesize);
    void setPixelSpacing(float *pixelspacing);
    void setMatrixS1toS2(float *matrix);
    void generateEulerTransform(double rotateX,double rotateY,double rotateZ,double transX,double transY,double transZ);
protected:
    void generateTransformMatrix(float *);
private:
    int imageSize[2]{512,512};
    double pixelSpacing[2]{0.5,0.5};
    bool isFront{true};
    float centerOffset[3]{0,0,0};
    float mtS1ToS2[16];
    float focalDistance{1000};
    TransformPointer m_Transform{nullptr};
};
#endif // GENERATEMATRIXHELPER_H
