#include "THADicomSeries.h"

#include <QDebug>

THADicomSeries::THADicomSeries()
{

}
THADicomSeries::~THADicomSeries()
{

    mDICOMReader->Register(0);
    //qDebug() << mDICOMReader;
    mDICOMReader->Delete();
    //qDebug() << "THADicomSeries::~THADicomSeries";
}
vtkImageData* THADicomSeries::f_GetImageData()
{

    return mDICOMReader->GetOutput();
}

vtkDICOMReader* THADicomSeries::f_GetDICOMReader()
{
    return mDICOMReader;
}

void THADicomSeries::f_GetDicomSeries(vtkStringArray* pDicomSeriesFiles)
{
    // 处理
    int imageDims[3];
	uStatus::mDICOMReader = nullptr;

    mDICOMReader->SetFileNames(pDicomSeriesFiles);
    mDICOMReader->SetMemoryRowOrderToFileNative();
    mDICOMReader->SortingOn();
    mDICOMReader->Update();

    uStatus::mDICOMReader = mDICOMReader;
    mDICOMReader->GetOutput()->GetDimensions(uStatus::Dims);
    qDebug() << "dicom dimensions:" << uStatus::Dims[0] << uStatus::Dims[1]  <<uStatus::Dims[2]  ;


    mDICOMReader->GetOutput()->GetExtent(uStatus::mExtent);

    mDICOMReader->GetOutput()->GetOrigin(uStatus::mOrigin);


    mDICOMReader->GetOutput()->GetSpacing(uStatus::mSpacing);
    //p_uStatus->f_SetSpacing(mSpacing);
    qDebug() << "dicom thickness: "  << QString::number(uStatus::mSpacing[0] , 'f', 6)<< "," << QString::number(uStatus::mSpacing[1], 'f', 6)<< "," << QString::number(uStatus::mSpacing[2], 'f', 6);

    uStatus::mCenter[0] = uStatus::mOrigin[0] + uStatus::mSpacing[0] * 0.5 * (uStatus::mExtent[0] + uStatus::mExtent[1]);
    uStatus::mCenter[1] = uStatus::mOrigin[1] + uStatus::mSpacing[1] * 0.5 * (uStatus::mExtent[2] + uStatus::mExtent[3]);
    uStatus::mCenter[2] = uStatus::mOrigin[2] + uStatus::mSpacing[2] * 0.5 * (uStatus::mExtent[4] + uStatus::mExtent[5]);

    mDicomMetaData = mDICOMReader->GetMetaData();
    fileIndexArray = mDICOMReader->GetFileIndexArray() ;
    fameIndexArray = mDICOMReader->GetFrameIndexArray() ;

    vtkDICOMValue mWindowWidth = mDicomMetaData->GetAttributeValue(0, DC::WindowWidth);
    vtkDICOMValue mWindowCenter = mDicomMetaData->GetAttributeValue(0, DC::WindowCenter);

    if(mWindowWidth.IsValid() && mWindowCenter.IsValid()){
        uStatus::mWW = mWindowWidth.GetInt(0);
        uStatus::mWL = mWindowCenter.GetInt(0);
    }
    else
    {
        short minPixel;
        short maxPixel;
        bool start = true;
        for(int x=uStatus::mExtent[0];x<=uStatus::mExtent[1];x++){
            for(int y=uStatus::mExtent[2];y<=uStatus::mExtent[3];y++){
                int z=uStatus::mExtent[4];
                short* scalar = (short *)mDICOMReader->GetOutput()->GetScalarPointer(x, y, z) ;
                if(start || *scalar>maxPixel){
                    maxPixel = *scalar;
                }
                if(start || *scalar<minPixel){
                    minPixel = *scalar;
                }
                start=false;
            }
        }
        //qDebug() << "f_Load_Files : minPixel: " <<  minPixel << " maxPixel: " << maxPixel;
        uStatus::mWW = maxPixel - minPixel;
        uStatus::mWL = minPixel + uStatus::mWW / 2;
        uStatus::mDefaultWW = maxPixel - minPixel;
        uStatus::mDefaultWL = minPixel + uStatus::mWW / 2;
    }

}
