#pragma execution_character_set("utf-8")
#include "EmptyImageData.h"
#include "vtkObjectFactory.h"


vtkStandardNewMacro(EmptyImageData);

EmptyImageData::EmptyImageData()
{

}

EmptyImageData::~EmptyImageData()
{
//    mInfomation->Delete();
//    mEmpty->Delete();
}


vtkImageData* EmptyImageData::f_GetEmptyImageData()
{

    mEmpty->SetSpacing(0.5,0.5,1.0);
    mEmpty->SetExtent(0 , 511  , 0 , 511, 0 , 0);
    mEmpty->SetOrigin(0 , 0 , 0);
    mEmpty->SetScalarType(VTK_SHORT,mInfomation);
    mEmpty->AllocateScalars(VTK_SHORT,2);

    for(int i=0; i<512; i++)
    {
        for(int j=0; j<512; j++)
        {
                int *pixel_mask = (int*) mEmpty->GetScalarPointer (j,i,0);
                *pixel_mask = -1024;
        }
    }

    return mEmpty;
}
