#ifndef DicomTag_H
#define DicomTag_H

#include <vtkStringArray.h>
#include <vtkDICOMReader.h>
#include <vtkDICOMValue.h>
#include <vtkDICOMMetaData.h>
#include <vtkIntArray.h>
#include <QString>
#include <QDebug>
#include <vtkNew.h>
#include "ustatus.h"

class DicomTag
{
public:
    DicomTag();
    ~DicomTag();
    vtkSmartPointer<vtkDICOMMetaData> mDicomMetaData;
    vtkSmartPointer<vtkIntArray> fileIndexArray;
    vtkSmartPointer<vtkIntArray> fameIndexArray;
    int sliceIndex = 0 ;
    int fileIndex  = 0 ;
    int frameIndex = 0 ;

    vtkNew<vtkDICOMReader> mDICOMReader;
    void f_OpenFileArray(vtkStringArray* pDicomSeriesFiles);
    void f_SetSystemData();
    //////////
    std::string f_GetPatientName(int pCompIdx=0);
    std::string f_GetSOPUID(int pCompIdx=0);
    std::string f_GetPatientID(int pCompIdx=0);
    std::string f_GetPatientBirthDay(int pCompIdx=0);
    std::string f_GetPatientSex(int pCompIdx=0);
    std::string f_GetContentDate(int pCompIdx=0);
    std::string f_GetContentTime(int pCompIdx=0);
    std::string f_GetPatientBirthTime(int pCompIdx=0);
    std::string f_GetPatientWeight(int pCompIdx=0);
    std::string f_GetPregnancyStatus(int pCompIdx=0);
    std::string f_GetSliceThickness(int pCompIdx=0);
    std::string f_GetFileMetaInformationGroupLength(int pCompIdx=0);
    std::string f_GetFileMetaInformationVersion(int pCompIdx=0);
    std::string f_GetMediaStorageSOPClassUID(int pCompIdx=0);
    std::string f_GetMediaStorageSOPInstanceUID(int pCompIdx=0);
    std::string f_GetTransferSyntaxUID(int pCompIdx=0);
    std::string f_GetImplementationClassUID(int pCompIdx=0);
//    std::string f_GetRTVCommunicationSOPClassUID(int pCompIdx=0);
//    std::string f_GetRTVCommunicationSOPInstanceUID(int pCompIdx=0);
    std::string f_GetPrivateInformationCreatorUID(int pCompIdx=0);
    std::string f_GetPrivateRecordUID(int pCompIdx=0);
    std::string f_GetInstanceCreatorUID(int pCompIdx=0);
    std::string f_GetSOPClassUID(int pCompIdx=0);
    std::string f_GetOriginalSpecializedSOPClassUID(int pCompIdx=0);
    std::string f_GetCodingSchemeUID(int pCompIdx=0);
    std::string f_GetContextGroupExtensionCreatorUID(int pCompIdx=0);
    std::string f_GetImageType(int pCompIdx=0);
    std::string f_GetInstanceCreationDate(int pCompIdx=0);
    std::string f_GetInstanceCreationTime(int pCompIdx=0);
    std::string f_GetStudyDate(int pCompIdx=0);
    std::string f_GetAcquisitionDate(int pCompIdx=0);
    std::string f_GetAcquisitionDateTime(int pCompIdx=0);
    std::string f_GetStudyTime(int pCompIdx=0);
    std::string f_GetAcquisitionTime(int pCompIdx=0);
    std::string f_GetModality(int pCompIdx=0);
    std::string f_GetReferencedSOPClassUID(int pCompIdx=0);
    std::string f_GetReferencedSOPInstanceUID(int pCompIdx=0);
    std::string f_GetRotationDirection(int pCompIdx=0);
    std::string f_GetStudyInstanceUID(int pCompIdx=0);
    std::string f_GetSeriesInstanceUID(int pCompIdx=0);
    std::string f_GetFrameOfReferenceUID(int pCompIdx=0);
    std::string f_GetRows(int pCompIdx=0);
    std::string f_GetColumns(int pCompIdx=0);
    std::string f_GetBitsAllocated(int pCompIdx=0);
    std::string f_GetBitsStored(int pCompIdx=0);
    std::string f_GetHighBit(int pCompIdx=0);
    std::string f_GetAccessionNumber(int pCompIdx=0);
    std::string f_GetReferencedImageSequence(int pCompIdx=0);
    std::string f_GetItem(int pCompIdx=0);
    std::string f_GetXRayTubeCurrent(int pCompIdx=0);
    std::string f_GetAcquisitionNumber(int pCompIdx=0);
    std::string f_GetPositionReferenceIndicator(int pCompIdx=0);
    std::string f_GetPreMedication(int pCompIdx=0);
    std::string f_GetReferringPhysicianName(int pCompIdx=0);
    std::string f_GetRequestedContrastAgent(int pCompIdx=0);
    std::string f_GetSeriesNumber(int pCompIdx=0);
    std::string f_GetInstanceNumber(int pCompIdx=0);
    std::string f_GetInstitutionName(int pCompIdx=0);
    std::string f_GetStudyDescription(int pCompIdx=0);
    std::string f_GetSeriesDescription(int pCompIdx=0);
    std::string f_GetImageComments(int pCompIdx=0);
    std::string f_GetWindowCenter(int pCompIdx=0);
    std::string f_GetWindowWidth(int pCompIdx=0);
    std::string f_GetSliceLocation(int pCompIdx=0);
    std::string f_GetPatientPosition(int pCompIdx=0);
    std::string f_GetRescaleIntercept(int pCompIdx=0);
    std::string f_GetRescaleSlope(int pCompIdx=0);
    int f_GetCount();
    /////////
    /// \brief f_GetPatientSex
    /// \param pCompIdx
    /// \return
    ///
    std::string *mPatientName;
//    std::string *mPatientName;
//    std::string *mPatientName;
//    std::string *mPatientName;
//    std::string *mPatientName;
//    std::string *mPatientName;

    vtkDICOMValue f_GetdicomValue(int pCompIdx,vtkDICOMTag p);

};

#endif // DicomTag_H
