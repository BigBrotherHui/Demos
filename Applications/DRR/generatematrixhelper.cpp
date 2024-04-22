#pragma execution_character_set("utf-8")
#include "generatematrixhelper.h"
#include <vtkMatrix4x4.h>
#include <QDebug>
void copyArray(float* src, float* dst, size_t size)
{
    for (size_t i = 0; i < size; i++)
    {
        dst[i] = src[i];
    }
}

void  GenerateMatrixHelper::generateTransformMatrix(float* M)
{
    float cs = imageSize[0] / 2 + 0.001; //竖直改变，取projM 中间值
    float ls = imageSize[1] / 2 + 0.001; //水平改变，取projN 中间值
    float offset = 0.001;
    cs = cs + offset;
    ls = ls + offset;

    float sp = pixelSpacing[0]; //平板像素分辨率，默认体数据分辨率为1
    float d = focalDistance;

    float P[3][3];//存储临时的三维边界点
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            P[i][k] = 0.0;
        }
    }
    P[0][0] = cs / d;
    P[0][1] = 1 / sp;
    P[0][2] = 0;
    P[1][0] = ls / d;
    P[1][1] = 0;
    P[1][2] = 1 / sp;
    P[2][0] = 1 / d;
    P[2][1] = 0.0;
    P[2][2] = 0.0;

    float Rs1tos2[4][4]{};
    memset(Rs1tos2, 0, sizeof(float) * 16);
    if (!isFront) //判定正侧位
    {
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                int ind = i * 4 + j;
                Rs1tos2[i][j] = mtS1ToS2[ind];
            }

        }
    }
    else
    {
        Rs1tos2[0][0] = 1;
        Rs1tos2[1][1] = 1;
        Rs1tos2[2][2] = 1;
        Rs1tos2[3][3] = 1;
    }

    TransformType::MatrixType R = m_Transform->GetMatrix();  //旋转变换矩阵
    float R2[4][4]{};
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
            R2[i][j] = R[i][j];
    }
    R2[3][3] = 1;

    float T[4][4]{};//分割结果相对CT中心偏移位置
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 4; k++) {
            T[i][k] = 0.0;
        }
    }
    T[0][0] = 1;
    T[1][1] = 1;
    T[2][2] = 1;
    T[3][3] = 1;
    T[0][3] = centerOffset[0];
    T[1][3] = centerOffset[1];
    T[2][3] = centerOffset[2];


    // M = P * T*R;
    //R/R2是旋转矩阵
    //T是中心偏移矩阵
    //P是
    float R3[4][4];
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            R3[i][j] = 0;
            for (int k = 0; k < 4; k++)
                R3[i][j] += R2[i][k] * T[k][j];
        }
    }

    auto par = m_Transform->GetParameters();  //位移变换矩阵
    float T2[4][4]{};
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 4; k++) {
            T2[i][k] = 0.0;
        }
    }
    T2[0][0] = 1;
    T2[1][1] = 1;
    T2[2][2] = 1;
    T2[3][3] = 1;
    T2[0][3] = par[3];
    T2[1][3] = par[4];
    T2[2][3] = par[5];

    // M = P * T*R;
    float M1[4][4];
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            M1[i][j] = 0;
            for (int k = 0; k < 4; k++)
                M1[i][j] += T2[i][k] * R3[k][j];
        }
    }

    float M2[4][4];
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            M2[i][j] = 0;
            for (int k = 0; k < 4; k++)
                M2[i][j] += Rs1tos2[i][k] * M1[k][j];
        }
    }

    int ind = 0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            ind = i * 4 + j;
            M[ind] = 0;
            for (int k = 0; k < 3; k++)
                M[ind] += P[i][k] * M2[k][j];
        }
    }
}

GenerateMatrixHelper::~GenerateMatrixHelper()
{

}

void GenerateMatrixHelper::getTransformMatrix(bool isfront, float* M)
{
    this->isFront = isfront;
    generateTransformMatrix(M);
}

void GenerateMatrixHelper::setCenterOffset(float* centeroffset)
{
    for (size_t i = 0; i < 3; i++)
    {
        centerOffset[i] = centeroffset[i];
    }
}

void GenerateMatrixHelper::setFocalDistance(float focaldistance)
{
    focalDistance = focaldistance;
}

void GenerateMatrixHelper::setImageSize(int* imagesize)
{
    imageSize[0] = imagesize[0];
    imageSize[1] = imagesize[1];
}

void GenerateMatrixHelper::setPixelSpacing(float* pixelspacing)
{
    pixelSpacing[0] = pixelspacing[0];
    pixelSpacing[1] = pixelspacing[1];
}

void GenerateMatrixHelper::setMatrixS1toS2(float* matrix)
{
    memcpy(mtS1ToS2, matrix, sizeof(float) * 16);
}

void GenerateMatrixHelper::generateEulerTransform(double rotateX, double rotateY, double rotateZ, double transX, double transY, double transZ)
{
    itk::Vector<double, 3> vecTranslation;
    vecTranslation.SetElement(0, transX);
    vecTranslation.SetElement(1, transY);
    vecTranslation.SetElement(2, transZ);
    m_Transform->SetTranslation(vecTranslation);
    const double dtr = (atan(1.0) * 4.0) / 180.0;
    m_Transform->SetRotation(dtr * rotateX, dtr * rotateY, dtr * rotateZ);
}

GenerateMatrixHelper::GenerateMatrixHelper()
{
    m_Transform = itk::Euler3DTransform<float>::New();
    m_Transform->SetComputeZYX(1);
}
