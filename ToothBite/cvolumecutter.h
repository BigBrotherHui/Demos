///  @file CVolumeCutter.h
///  @brief 对vtkImageData进行截取
///  @author project4gogo@163.com
///  @date 2018
///
///

#ifndef CVOLUMECUTTER_H
#define CVOLUMECUTTER_H

#include <QtGlobal>

class vtkPolyData;
class vtkLinearExtrusionFilter;
class vtkTransformPolyDataFilter;
class vtkPolyDataToImageStencil;
class vtkImageStencil;
class vtkExtractVOI;
class vtkImageData;
class vtkImageAlgorithm;


class CVolumeCutter
{
public:
    enum enAxis{AXIS_X, AXIS_Y, AXIS_Z, AXIS_XYZ};
public:
    //static CVolumeCutter* GetInstace();
    CVolumeCutter();
    ~CVolumeCutter();

    void SetThickness(double thickness = 5);
    void SetBackgroundValue(qint32 bgValue = 0);
    void SetCutAxix(enAxis axis = AXIS_X);
    void SetCutAxix(double direction[3]);
    void SetCutAxix(double dx, double dy, double dz);

    //平移
    void SetTranslate(double dx, double dy, double dz);

    void GetCutAreaBounds(double bounds[6]);

    //保留vtkPolyData包围的所有点
    vtkImageData *GetCutResult(vtkImageData *pSrc, vtkPolyData *pCutPolyData);


private:
    void updateSrcSource(vtkImageData *pImageData);

private:
    vtkPolyData*                 m_pMaskPolyData = nullptr;
    vtkLinearExtrusionFilter*    m_pMaskExtrusionFilter = nullptr;
    vtkTransformPolyDataFilter*  m_pMaskTransformPolyDatafilter = nullptr;
    vtkPolyDataToImageStencil*   m_pMaskPolydataToStencil = nullptr;
    vtkImageStencil*             m_pMaskImageStencil = nullptr;
    vtkExtractVOI*               m_pMaskExtractVoi = nullptr;

    vtkImageData*                m_pMaskOutput = nullptr;


    double m_iThickness = 5;
    qint32 m_iBackgroundValue = 0;
    enAxis m_enAxis = AXIS_X;
};

#endif // CVOLUMECUTTER_H
