#pragma execution_character_set("utf-8")
#include "DicomTag.h"
#include <QJsonArray>
#include <QJsonObject>
#include <QTabWidget>

DicomTag::DicomTag()
{
    //qDebug() << "DicomTag : DicomTag : ";
}
DicomTag::~DicomTag()
{
//    mDICOMReader->Delete();
}
void DicomTag::f_OpenFileArray(vtkStringArray* pDicomSeriesFiles)
{
    mDICOMReader->SetFileNames(pDicomSeriesFiles);
    mDICOMReader->SetMemoryRowOrderToFileNative();
    mDICOMReader->SortingOn();
    mDICOMReader->UpdateInformation() ;
    mDICOMReader->Update();

    DicomTag::mDicomMetaData = DicomTag::mDICOMReader->GetMetaData();
    DicomTag::fileIndexArray = DicomTag::mDicomMetaData->GetFileIndexArray() ;
    DicomTag::fameIndexArray = DicomTag::mDicomMetaData->GetFrameIndexArray() ;
}
void DicomTag::f_SetSystemData()
{
    //qDebug() << 10;
    DicomTag::mDicomMetaData = uStatus::mDICOMReader->GetMetaData();
    //qDebug() << 20;
    DicomTag::fileIndexArray = DicomTag::mDicomMetaData->GetFileIndexArray() ;
    //qDebug() << 30;
    DicomTag::fameIndexArray = DicomTag::mDicomMetaData->GetFrameIndexArray() ;
}
int DicomTag::f_GetCount()
{
    return -1;
}

std::string DicomTag::f_GetPatientName(int ppCompIdx)
{

    std::string mPatientName;
    std::string StrPatientName;
    vtkDICOMValue pv =f_GetdicomValue(ppCompIdx,DC::PatientName);
    //pv.GetValues(mPatientName , 0) ;
    //    int cnt = pv.GetNumberOfValues ();
    StrPatientName=pv.GetString(0);
    QString qStr=QString::fromStdString(StrPatientName);
    qDebug() << "dicom:patient name" << qStr;
    return StrPatientName;
}

std::string DicomTag::f_GetSOPUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SOPInstanceUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: SOPInstanceUID: " << qStr;
    return mStr;
}

std::string DicomTag::f_GetPatientID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PatientID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: PatientID: " << qStr;
    return mStr;
}

std::string DicomTag::f_GetPatientBirthDay(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PatientBirthDate);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: patient birth date: " << qStr;
    return mStr;
}
std::string DicomTag::f_GetPatientSex(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PatientSex);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: patient gender: " << qStr;
    return mStr;
}

//影像拍摄的日期.
std::string DicomTag::f_GetContentDate(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ContentDate);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: images create date: " << qStr;
    return mStr;
}

//影像拍摄的时间
std::string DicomTag::f_GetContentTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ContentTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: images create date: " << qStr;
    return mStr;
}

//患者出生时间
std::string DicomTag::f_GetPatientBirthTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PatientBirthTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: patient birth time: " << qStr;
    return mStr;
}

//患者体重
std::string DicomTag::f_GetPatientWeight(int pCompIdx)
{
    vtkDICOMValue pv = f_GetdicomValue(pCompIdx,DC ::PatientWeight);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: patient weight: " << qStr;
    return mStr;
}

// 怀孕状态
std::string DicomTag::f_GetPregnancyStatus(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PregnancyStatus);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: patient pregnancy status: " << qStr;
    return mStr;
}

//层厚
std::string DicomTag::f_GetSliceThickness(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SliceThickness);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: thickness: " << qStr;
    return mStr;
}


//FileMetaInformationGroupLength
std::string DicomTag::f_GetFileMetaInformationGroupLength(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::FileMetaInformationGroupLength);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: FileMetaInformationGroupLength: " << qStr;
    return mStr;
}

//FileMetaInformationVersion
std::string DicomTag::f_GetFileMetaInformationVersion(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::FileMetaInformationVersion);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: FileMetaInformationVersion: " << qStr;
    return mStr;
}

//MediaStorageSOPClassUID
std::string DicomTag::f_GetMediaStorageSOPClassUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::MediaStorageSOPClassUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: MediaStorageSOPClassUID: " << qStr;
    return mStr;
}
//MediaStorageSOPInstanceUID
std::string DicomTag::f_GetMediaStorageSOPInstanceUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::MediaStorageSOPInstanceUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: MediaStorageSOPInstanceUID: " << qStr;
    return mStr;
}

//TransferSyntaxUID
std::string DicomTag::f_GetTransferSyntaxUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::TransferSyntaxUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: TransferSyntaxUID: " << qStr;
    return mStr;
}

//ImplementationClassUID
std::string DicomTag::f_GetImplementationClassUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ImplementationClassUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: ImplementationClassUID: " << qStr;
    return mStr;
}

////RTVCommunicationSOPClassUID
//std::string DicomTag::f_GetRTVCommunicationSOPClassUID(int pCompIdx)
//{
//    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RTVCommunicationSOPClassUID);
//    std::string mStr;
//    mStr=pv.GetString(0);
//    QString qStr=QString::fromStdString(mStr);
//    qDebug() << "DicomTag: RTVCommunicationSOPClassUID: " << qStr;
//    return mStr;
//}

////RTVCommunicationSOPInstanceUID
//std::string DicomTag::f_GetRTVCommunicationSOPInstanceUID(int pCompIdx)
//{
//    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RTVCommunicationSOPInstanceUID);
//    std::string mStr;
//    mStr=pv.GetString(0);
//    QString qStr=QString::fromStdString(mStr);
//    qDebug() << "DicomTag: RTVCommunicationSOPInstanceUID: " << qStr;
//    return mStr;
//}

//PrivateInformationCreatorUID
std::string DicomTag::f_GetPrivateInformationCreatorUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PrivateInformationCreatorUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: PrivateInformationCreatorUID: " << qStr;
    return mStr;
}

//PrivateRecordUID
std::string DicomTag::f_GetPrivateRecordUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PrivateRecordUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: PrivateRecordUID: " << qStr;
    return mStr;
}

//InstanceCreatorUID
std::string DicomTag::f_GetInstanceCreatorUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::InstanceCreatorUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: InstanceCreatorUID: " << qStr;
    return mStr;
}

//SOPClassUID
std::string DicomTag::f_GetSOPClassUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SOPClassUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: SOPClassUID: " << qStr;
    return mStr;
}

//OriginalSpecializedSOPClassUID
std::string DicomTag::f_GetOriginalSpecializedSOPClassUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::OriginalSpecializedSOPClassUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: OriginalSpecializedSOPClassUID: " << qStr;
    return mStr;
}

//CodingSchemeUID
std::string DicomTag::f_GetCodingSchemeUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::CodingSchemeUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: CodingSchemeUID: " << qStr;
    return mStr;
}
//ContextGroupExtensionCreatorUID
std::string DicomTag::f_GetContextGroupExtensionCreatorUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ContextGroupExtensionCreatorUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: ContextGroupExtensionCreatorUID: " << qStr;
    return mStr;
}

//ImageType
std::string DicomTag::f_GetImageType(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ImageType);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: ImageType: " << qStr;
    return mStr;
}

//InstanceCreationDate
std::string DicomTag::f_GetInstanceCreationDate(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::InstanceCreationDate);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: InstanceCreationDate: " << qStr;
    return mStr;
}

//InstanceCreationTime
std::string DicomTag::f_GetInstanceCreationTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::InstanceCreationTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: InstanceCreationTime: " << qStr;
    return mStr;
}

//StudyDate
std::string DicomTag::f_GetStudyDate(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::StudyDate);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: StudyDate: " << qStr;
    return mStr;
}

//AcquisitionDate
std::string DicomTag::f_GetAcquisitionDate(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::AcquisitionDate);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: AcquisitionDate: " << qStr;
    return mStr;
}

//AcquisitionDateTime
std::string DicomTag::f_GetAcquisitionDateTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::AcquisitionDateTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: AcquisitionDateTime: " << qStr;
    return mStr;
}

//StudyTime
std::string DicomTag::f_GetStudyTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::StudyTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: StudyTime: " << qStr;
    return mStr;
}

//AcquisitionTime
std::string DicomTag::f_GetAcquisitionTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::AcquisitionTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: AcquisitionTime: " << qStr;
    return mStr;
}

//Modality
std::string DicomTag::f_GetModality(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Modality);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: Modality: " << qStr;
    return mStr;
}

//ReferencedSOPClassUID
std::string DicomTag::f_GetReferencedSOPClassUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Modality);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: ReferencedSOPClassUID: " << qStr;
    return mStr;
}

//ReferencedSOPInstanceUID
std::string DicomTag::f_GetReferencedSOPInstanceUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Modality);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: ReferencedSOPInstanceUID: " << qStr;
    return mStr;
}

//RotationDirection
std::string DicomTag::f_GetRotationDirection(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RotationDirection);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: RotationDirection: " << qStr;
    return mStr;
}

//StudyInstanceUID
std::string DicomTag::f_GetStudyInstanceUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::StudyInstanceUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: StudyInstanceUID: " << qStr;
    return mStr;
}

//SeriesInstanceUID
std::string DicomTag::f_GetSeriesInstanceUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SeriesInstanceUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: SeriesInstanceUID: " << qStr;
    return mStr;
}

//FrameOfReferenceUID
std::string DicomTag::f_GetFrameOfReferenceUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::FrameOfReferenceUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: FrameOfReferenceUID: " << qStr;
    return mStr;
}

//Rows
std::string DicomTag::f_GetRows(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Rows);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: rows: " << qStr;
    return mStr;
}

//Columns
std::string DicomTag::f_GetColumns(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Columns);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: columns: " << qStr;
    return mStr;
}

//BitsAllocated
std::string DicomTag::f_GetBitsAllocated(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::BitsAllocated);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: BitsAllocated: " << qStr;
    return mStr;
}

//BitsStored
std::string DicomTag::f_GetBitsStored(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::BitsStored);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: BitsStored: " << qStr;
    return mStr;
}

//HighBit
std::string DicomTag::f_GetHighBit(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::HighBit);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: HighBit: " << qStr;
    return mStr;
}

//AccessionNumber
std::string DicomTag::f_GetAccessionNumber(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::AccessionNumber);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: AccessionNumber: " << qStr;
    return mStr;
}
//ReferencedImageSequence
std::string DicomTag::f_GetReferencedImageSequence(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ReferencedImageSequence);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: ReferencedImageSequence: " << qStr;
    return mStr;
}
//Item
std::string DicomTag::f_GetItem(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Item);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: Item: " << qStr;
    return mStr;
}
//XRayTubeCurrent
std::string DicomTag::f_GetXRayTubeCurrent(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::XRayTubeCurrent);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: XRayTubeCurrent: " << qStr;
    return mStr;
}
//AcquisitionNumber
std::string DicomTag::f_GetAcquisitionNumber(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::AcquisitionNumber);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: AcquisitionNumber: " << qStr;
    return mStr;
}

//PositionReferenceIndicator
std::string DicomTag::f_GetPositionReferenceIndicator(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PositionReferenceIndicator);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: PositionReferenceIndicator: " << qStr;
    return mStr;
}
//PreMedication
std::string DicomTag::f_GetPreMedication(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PreMedication);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: PreMedication: " << qStr;
    return mStr;
}
//RequestedContrastAgent
std::string DicomTag::f_GetRequestedContrastAgent(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RequestedContrastAgent);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: RequestedContrastAgent: " << qStr;
    return mStr;
}
//ReferringPhysicianName
std::string DicomTag::f_GetReferringPhysicianName(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ReferringPhysicianName);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: ReferringPhysicianName: " << qStr;
    return mStr;
}
//SeriesNumber
std::string DicomTag::f_GetSeriesNumber(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SeriesNumber);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: SeriesNumber: " << qStr;

    return mStr;
}

//InstanceNumber
std::string DicomTag::f_GetInstanceNumber(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::InstanceNumber);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: InstanceNumber: " << qStr;

    return mStr;
}


//InstitutionName
std::string DicomTag::f_GetInstitutionName(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::InstitutionName);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: InstitutionName: " << qStr;

    return mStr;
}

//StudyDescription
std::string DicomTag::f_GetStudyDescription(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::StudyDescription);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: StudyDescription: " << qStr;

    return mStr;
}

//SeriesDescription
std::string DicomTag::f_GetSeriesDescription(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SeriesDescription);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: SeriesDescription: " << qStr;

    return mStr;
}


//ImageComments
std::string DicomTag::f_GetImageComments(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ImageComments);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: ImageComments: " << qStr;

    return mStr;
}

//WindowCenter
std::string DicomTag::f_GetWindowCenter(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::WindowCenter);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: WindowCenter: " << qStr;

    return mStr;
}

//WindowWidth
std::string DicomTag::f_GetWindowWidth(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::WindowWidth);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: WindowWidth: " << qStr;

    return mStr;
}

//SliceLocation
std::string DicomTag::f_GetSliceLocation(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SliceLocation);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: SliceLocation: " << qStr;

    return mStr;
}

//SliceLocation
std::string DicomTag::f_GetPatientPosition(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PatientPosition);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: PatientPosition: " << qStr;

    return mStr;
}

//RescaleIntercept
std::string DicomTag::f_GetRescaleIntercept(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RescaleIntercept);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: RescaleIntercept: " << qStr;

    return mStr;
}

//RescaleSlope
std::string DicomTag::f_GetRescaleSlope(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RescaleSlope);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "DicomTag: RescaleSlope: " << qStr;

    return mStr;
}

vtkDICOMValue DicomTag::f_GetdicomValue(int pCompIdx,vtkDICOMTag p)
{
    DicomTag::fileIndex  = DicomTag::fileIndexArray->GetComponent (DicomTag::sliceIndex,pCompIdx) ;
    DicomTag::frameIndex = DicomTag::fameIndexArray->GetComponent (DicomTag::frameIndex,pCompIdx) ;
    vtkDICOMValue pv = DicomTag::mDicomMetaData->GetAttributeValue (DicomTag::fileIndex,
                                                                     DicomTag::frameIndex,p);
    return  pv;
}

