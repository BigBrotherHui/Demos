#include "cvolumecutter.h"

//vtk .h file
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkTransform.h>
#include <vtkExtractVOI.h>
#include <vtkImageData.h>
#include <vtkImageAlgorithm.h>

#define RELEASE_VTK(p) if(nullptr != p) { p->Delete(); p = nullptr;}


CVolumeCutter::CVolumeCutter()
{
    if(nullptr == m_pMaskPolyData)
    {
        m_pMaskPolyData = vtkPolyData::New();
    }


    if(nullptr == m_pMaskTransformPolyDatafilter)
    {
        m_pMaskTransformPolyDatafilter = vtkTransformPolyDataFilter::New();
        m_pMaskTransformPolyDatafilter->SetInputData(m_pMaskPolyData);
//        vtkSmartPointer<vtkTransform> transform = vtkTransform::New();
//        transform->Translate(0, 0, 0);
//        m_pMaskTransformPolyDatafilter->SetTransform(transform);
    }


    if(nullptr == m_pMaskExtrusionFilter)
    {
        m_pMaskExtrusionFilter = vtkLinearExtrusionFilter::New();
        m_pMaskExtrusionFilter->SetInputConnection(m_pMaskTransformPolyDatafilter->GetOutputPort());
        m_pMaskExtrusionFilter->SetScaleFactor(m_iThickness);
        m_pMaskExtrusionFilter->SetExtrusionTypeToNormalExtrusion();
        m_pMaskExtrusionFilter->SetVector(1, 0, 0);
    }

    if(nullptr == m_pMaskExtractVoi)
    {
        m_pMaskExtractVoi = vtkExtractVOI::New();
        m_pMaskExtractVoi->SetSampleRate(1, 1, 1);
        m_pMaskExtractVoi->ReleaseDataFlagOff();
    }


    if(nullptr == m_pMaskPolydataToStencil)
    {
        m_pMaskPolydataToStencil = vtkPolyDataToImageStencil::New();
        m_pMaskPolydataToStencil->SetInputConnection( m_pMaskExtrusionFilter->GetOutputPort() );
        m_pMaskPolydataToStencil->SetInformationInput( m_pMaskExtractVoi->GetOutput() );
    }

    if(nullptr == m_pMaskImageStencil)
    {
        m_pMaskImageStencil = vtkImageStencil::New();
        m_pMaskImageStencil->SetInputConnection(m_pMaskExtractVoi->GetOutputPort());
        m_pMaskImageStencil->SetStencilConnection(m_pMaskPolydataToStencil->GetOutputPort());
        m_pMaskImageStencil->ReverseStencilOff();
        m_pMaskImageStencil->SetBackgroundValue(m_iBackgroundValue);
    }

}

CVolumeCutter::~CVolumeCutter()
{
    RELEASE_VTK(m_pMaskPolyData);
    RELEASE_VTK(m_pMaskExtrusionFilter);
    RELEASE_VTK(m_pMaskTransformPolyDatafilter);
    RELEASE_VTK(m_pMaskPolydataToStencil);
    RELEASE_VTK(m_pMaskImageStencil);
    RELEASE_VTK(m_pMaskExtractVoi);
    RELEASE_VTK(m_pMaskOutput);
}

void CVolumeCutter::SetThickness(double thickness)
{
    m_iThickness = thickness;
    m_pMaskExtrusionFilter->SetScaleFactor(m_iThickness);
}

void CVolumeCutter::SetBackgroundValue(qint32 bgValue)
{
    m_iBackgroundValue = bgValue;
    m_pMaskImageStencil->SetBackgroundValue(m_iBackgroundValue);
}

void CVolumeCutter::SetCutAxix(CVolumeCutter::enAxis axis)
{
    m_enAxis = axis;
    m_pMaskExtrusionFilter->SetExtrusionTypeToVectorExtrusion();
    switch(m_enAxis)
    {
    case AXIS_X:
        m_pMaskExtrusionFilter->SetVector(1, 0, 0);
        break;
    case AXIS_Y:
        m_pMaskExtrusionFilter->SetVector(0, 1, 0);
        break;
    case AXIS_Z:
        m_pMaskExtrusionFilter->SetVector(0, 0, 1);
        break;
    default:
        //m_pMaskExtrusionFilter->SetExtrusionTypeToNormalExtrusion();
        m_pMaskExtrusionFilter->SetExtrusionTypeToPointExtrusion();
        break;
    }
}

void CVolumeCutter::SetCutAxix(double direction[])
{
    SetCutAxix(direction[0], direction[1], direction[2]);
}

void CVolumeCutter::SetCutAxix(double dx, double dy, double dz)
{
    m_pMaskExtrusionFilter->SetExtrusionTypeToVectorExtrusion();
    m_pMaskExtrusionFilter->SetVector(dx, dy, dz);
}

void CVolumeCutter::SetTranslate(double dx, double dy, double dz)
{
    vtkSmartPointer<vtkTransform> transform = vtkTransform::New();
    transform->Translate(dx, dy, dz);
    m_pMaskTransformPolyDatafilter->SetTransform(transform);

}


vtkImageData *CVolumeCutter::GetCutResult(vtkImageData *pSrc, vtkPolyData *pCutPolyData)
{
    vtkImageData *pResult = nullptr;
    if(nullptr != pSrc && nullptr != pCutPolyData)
    {
        updateSrcSource(pSrc);
        m_pMaskPolyData->DeepCopy(pCutPolyData);
        m_pMaskImageStencil->Update();
        pResult = m_pMaskImageStencil->GetOutput();
    }
    return pResult;
}



void CVolumeCutter::updateSrcSource(vtkImageData *pImageData)
{
    if(NULL != pImageData)
    {
        m_pMaskExtractVoi->SetVOI(pImageData->GetExtent());
        m_pMaskExtractVoi->SetSampleRate(1, 1, 1);
        m_pMaskExtractVoi->SetInputData(pImageData);
        m_pMaskExtractVoi->Update();
    }
}

void CVolumeCutter::GetCutAreaBounds(double bounds[6])
{
    m_pMaskExtrusionFilter->GetOutput()->GetBounds(bounds);
}
