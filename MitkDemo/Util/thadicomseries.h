#pragma execution_character_set("utf-8")
#ifndef THADicomSeries_H
#define THADicomSeries_H

#include <vtkStringArray.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>
#include <vtkNew.h>
#include <vtkImageChangeInformation.h>
#include <vtkDICOMMetaData.h>
#include "ustatus.h"
#include  "udicomtag.h"
class THADicomSeries
{
public:
    THADicomSeries();
    ~THADicomSeries();

	//vtkSmartPointer<vtkDICOMReader> mDICOMReader ;
	vtkNew<vtkDICOMReader> mDICOMReader;
    vtkImageData* mImageDate;
    vtkDICOMReader* f_GetDICOMReader();
    void f_GetDicomSeries(vtkStringArray* pDicomSeriesFiles);
    vtkImageData* f_GetImageData();
    vtkSmartPointer<vtkDICOMMetaData> mDicomMetaData;
    vtkSmartPointer<vtkIntArray> fileIndexArray;
    vtkSmartPointer<vtkIntArray> fameIndexArray;

};

#endif // THADicomSeries_H
