#include "udicomtag.h"
#include <QJsonArray>
#include <QJsonObject>
#include <QTabWidget>

uDicomTag::uDicomTag()
{
    //qDebug() << "uDicomTag : uDicomTag : ";
}
uDicomTag::~uDicomTag()
{
//    mDICOMReader->Delete();
}
void uDicomTag::f_OpenFileArray(vtkStringArray* pDicomSeriesFiles)
{
    mDICOMReader->SetFileNames(pDicomSeriesFiles);
    mDICOMReader->SetMemoryRowOrderToFileNative();
    mDICOMReader->SortingOn();
    mDICOMReader->UpdateInformation() ;
    mDICOMReader->Update();

    uDicomTag::mDicomMetaData = uDicomTag::mDICOMReader->GetMetaData();
    uDicomTag::fileIndexArray = uDicomTag::mDICOMReader->GetFileIndexArray() ;
    uDicomTag::fameIndexArray = uDicomTag::mDICOMReader->GetFrameIndexArray() ;
}
void uDicomTag::f_SetSystemData()
{
    //qDebug() << 10;
    uDicomTag::mDicomMetaData = mDICOMReader->GetMetaData();
    //qDebug() << 20;
    uDicomTag::fileIndexArray = mDICOMReader->GetFileIndexArray() ;
    //qDebug() << 30;
    uDicomTag::fameIndexArray = mDICOMReader->GetFrameIndexArray() ;
}
int uDicomTag::f_GetCount()
{
    return -1;
}

std::string uDicomTag::f_GetPatientName(int ppCompIdx)
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

std::string uDicomTag::f_GetSOPUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SOPInstanceUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: SOPInstanceUID: " << qStr;
    return mStr;
}

std::string uDicomTag::f_GetPatientID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PatientID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: PatientID: " << qStr;
    return mStr;
}

std::string uDicomTag::f_GetPatientBirthDay(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PatientBirthDate);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: patient birth date: " << qStr;
    return mStr;
}
std::string uDicomTag::f_GetPatientSex(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PatientSex);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: patient gender: " << qStr;
    return mStr;
}

//影像拍摄的日期.
std::string uDicomTag::f_GetContentDate(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ContentDate);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: images create date: " << qStr;
    return mStr;
}

//影像拍摄的时间
std::string uDicomTag::f_GetContentTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ContentTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: images create date: " << qStr;
    return mStr;
}

//患者出生时间
std::string uDicomTag::f_GetPatientBirthTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PatientBirthTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: patient birth time: " << qStr;
    return mStr;
}

//患者体重
std::string uDicomTag::f_GetPatientWeight(int pCompIdx)
{
    vtkDICOMValue pv = f_GetdicomValue(pCompIdx,DC ::PatientWeight);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: patient weight: " << qStr;
    return mStr;
}

// 怀孕状态
std::string uDicomTag::f_GetPregnancyStatus(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PregnancyStatus);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: patient pregnancy status: " << qStr;
    return mStr;
}

//层厚
std::string uDicomTag::f_GetSliceThickness(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SliceThickness);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: thickness: " << qStr;
    return mStr;
}


//FileMetaInformationGroupLength
std::string uDicomTag::f_GetFileMetaInformationGroupLength(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::FileMetaInformationGroupLength);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: FileMetaInformationGroupLength: " << qStr;
    return mStr;
}

//FileMetaInformationVersion
std::string uDicomTag::f_GetFileMetaInformationVersion(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::FileMetaInformationVersion);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: FileMetaInformationVersion: " << qStr;
    return mStr;
}

//MediaStorageSOPClassUID
std::string uDicomTag::f_GetMediaStorageSOPClassUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::MediaStorageSOPClassUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: MediaStorageSOPClassUID: " << qStr;
    return mStr;
}
//MediaStorageSOPInstanceUID
std::string uDicomTag::f_GetMediaStorageSOPInstanceUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::MediaStorageSOPInstanceUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: MediaStorageSOPInstanceUID: " << qStr;
    return mStr;
}

//TransferSyntaxUID
std::string uDicomTag::f_GetTransferSyntaxUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::TransferSyntaxUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: TransferSyntaxUID: " << qStr;
    return mStr;
}

//ImplementationClassUID
std::string uDicomTag::f_GetImplementationClassUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ImplementationClassUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: ImplementationClassUID: " << qStr;
    return mStr;
}

////RTVCommunicationSOPClassUID
//std::string uDicomTag::f_GetRTVCommunicationSOPClassUID(int pCompIdx)
//{
//    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RTVCommunicationSOPClassUID);
//    std::string mStr;
//    mStr=pv.GetString(0);
//    QString qStr=QString::fromStdString(mStr);
//    qDebug() << "uDicomTag: RTVCommunicationSOPClassUID: " << qStr;
//    return mStr;
//}

////RTVCommunicationSOPInstanceUID
//std::string uDicomTag::f_GetRTVCommunicationSOPInstanceUID(int pCompIdx)
//{
//    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RTVCommunicationSOPInstanceUID);
//    std::string mStr;
//    mStr=pv.GetString(0);
//    QString qStr=QString::fromStdString(mStr);
//    qDebug() << "uDicomTag: RTVCommunicationSOPInstanceUID: " << qStr;
//    return mStr;
//}

//PrivateInformationCreatorUID
std::string uDicomTag::f_GetPrivateInformationCreatorUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PrivateInformationCreatorUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: PrivateInformationCreatorUID: " << qStr;
    return mStr;
}

//PrivateRecordUID
std::string uDicomTag::f_GetPrivateRecordUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PrivateRecordUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: PrivateRecordUID: " << qStr;
    return mStr;
}

//InstanceCreatorUID
std::string uDicomTag::f_GetInstanceCreatorUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::InstanceCreatorUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: InstanceCreatorUID: " << qStr;
    return mStr;
}

//SOPClassUID
std::string uDicomTag::f_GetSOPClassUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SOPClassUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: SOPClassUID: " << qStr;
    return mStr;
}

//OriginalSpecializedSOPClassUID
std::string uDicomTag::f_GetOriginalSpecializedSOPClassUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::OriginalSpecializedSOPClassUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: OriginalSpecializedSOPClassUID: " << qStr;
    return mStr;
}

//CodingSchemeUID
std::string uDicomTag::f_GetCodingSchemeUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::CodingSchemeUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: CodingSchemeUID: " << qStr;
    return mStr;
}
//ContextGroupExtensionCreatorUID
std::string uDicomTag::f_GetContextGroupExtensionCreatorUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ContextGroupExtensionCreatorUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: ContextGroupExtensionCreatorUID: " << qStr;
    return mStr;
}

//ImageType
std::string uDicomTag::f_GetImageType(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ImageType);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: ImageType: " << qStr;
    return mStr;
}

//InstanceCreationDate
std::string uDicomTag::f_GetInstanceCreationDate(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::InstanceCreationDate);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: InstanceCreationDate: " << qStr;
    return mStr;
}

//InstanceCreationTime
std::string uDicomTag::f_GetInstanceCreationTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::InstanceCreationTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: InstanceCreationTime: " << qStr;
    return mStr;
}

//StudyDate
std::string uDicomTag::f_GetStudyDate(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::StudyDate);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: StudyDate: " << qStr;
    return mStr;
}

//AcquisitionDate
std::string uDicomTag::f_GetAcquisitionDate(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::AcquisitionDate);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: AcquisitionDate: " << qStr;
    return mStr;
}

//AcquisitionDateTime
std::string uDicomTag::f_GetAcquisitionDateTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::AcquisitionDateTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: AcquisitionDateTime: " << qStr;
    return mStr;
}

//StudyTime
std::string uDicomTag::f_GetStudyTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::StudyTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: StudyTime: " << qStr;
    return mStr;
}

//AcquisitionTime
std::string uDicomTag::f_GetAcquisitionTime(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::AcquisitionTime);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: AcquisitionTime: " << qStr;
    return mStr;
}

//Modality
std::string uDicomTag::f_GetModality(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Modality);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: Modality: " << qStr;
    return mStr;
}

//ReferencedSOPClassUID
std::string uDicomTag::f_GetReferencedSOPClassUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Modality);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: ReferencedSOPClassUID: " << qStr;
    return mStr;
}

//ReferencedSOPInstanceUID
std::string uDicomTag::f_GetReferencedSOPInstanceUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Modality);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: ReferencedSOPInstanceUID: " << qStr;
    return mStr;
}

//RotationDirection
std::string uDicomTag::f_GetRotationDirection(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RotationDirection);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: RotationDirection: " << qStr;
    return mStr;
}

//StudyInstanceUID
std::string uDicomTag::f_GetStudyInstanceUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::StudyInstanceUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: StudyInstanceUID: " << qStr;
    return mStr;
}

//SeriesInstanceUID
std::string uDicomTag::f_GetSeriesInstanceUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SeriesInstanceUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: SeriesInstanceUID: " << qStr;
    return mStr;
}

//FrameOfReferenceUID
std::string uDicomTag::f_GetFrameOfReferenceUID(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::FrameOfReferenceUID);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: FrameOfReferenceUID: " << qStr;
    return mStr;
}

//Rows
std::string uDicomTag::f_GetRows(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Rows);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: rows: " << qStr;
    return mStr;
}

//Columns
std::string uDicomTag::f_GetColumns(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Columns);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    qDebug() << "dicom: columns: " << qStr;
    return mStr;
}

//BitsAllocated
std::string uDicomTag::f_GetBitsAllocated(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::BitsAllocated);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: BitsAllocated: " << qStr;
    return mStr;
}

//BitsStored
std::string uDicomTag::f_GetBitsStored(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::BitsStored);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: BitsStored: " << qStr;
    return mStr;
}

//HighBit
std::string uDicomTag::f_GetHighBit(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::HighBit);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: HighBit: " << qStr;
    return mStr;
}

//AccessionNumber
std::string uDicomTag::f_GetAccessionNumber(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::AccessionNumber);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: AccessionNumber: " << qStr;
    return mStr;
}
//ReferencedImageSequence
std::string uDicomTag::f_GetReferencedImageSequence(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ReferencedImageSequence);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: ReferencedImageSequence: " << qStr;
    return mStr;
}
//Item
std::string uDicomTag::f_GetItem(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::Item);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: Item: " << qStr;
    return mStr;
}
//XRayTubeCurrent
std::string uDicomTag::f_GetXRayTubeCurrent(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::XRayTubeCurrent);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: XRayTubeCurrent: " << qStr;
    return mStr;
}
//AcquisitionNumber
std::string uDicomTag::f_GetAcquisitionNumber(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::AcquisitionNumber);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: AcquisitionNumber: " << qStr;
    return mStr;
}

//PositionReferenceIndicator
std::string uDicomTag::f_GetPositionReferenceIndicator(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PositionReferenceIndicator);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: PositionReferenceIndicator: " << qStr;
    return mStr;
}
//PreMedication
std::string uDicomTag::f_GetPreMedication(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PreMedication);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: PreMedication: " << qStr;
    return mStr;
}
//RequestedContrastAgent
std::string uDicomTag::f_GetRequestedContrastAgent(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RequestedContrastAgent);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: RequestedContrastAgent: " << qStr;
    return mStr;
}
//ReferringPhysicianName
std::string uDicomTag::f_GetReferringPhysicianName(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ReferringPhysicianName);
    std::string mStr;
    mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: ReferringPhysicianName: " << qStr;
    return mStr;
}
//SeriesNumber
std::string uDicomTag::f_GetSeriesNumber(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SeriesNumber);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: SeriesNumber: " << qStr;

    return mStr;
}

//InstanceNumber
std::string uDicomTag::f_GetInstanceNumber(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::InstanceNumber);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: InstanceNumber: " << qStr;

    return mStr;
}


//InstitutionName
std::string uDicomTag::f_GetInstitutionName(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::InstitutionName);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: InstitutionName: " << qStr;

    return mStr;
}

//StudyDescription
std::string uDicomTag::f_GetStudyDescription(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::StudyDescription);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: StudyDescription: " << qStr;

    return mStr;
}

//SeriesDescription
std::string uDicomTag::f_GetSeriesDescription(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SeriesDescription);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: SeriesDescription: " << qStr;

    return mStr;
}


//ImageComments
std::string uDicomTag::f_GetImageComments(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::ImageComments);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: ImageComments: " << qStr;

    return mStr;
}

//WindowCenter
std::string uDicomTag::f_GetWindowCenter(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::WindowCenter);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: WindowCenter: " << qStr;

    return mStr;
}

//WindowWidth
std::string uDicomTag::f_GetWindowWidth(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::WindowWidth);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: WindowWidth: " << qStr;

    return mStr;
}

//SliceLocation
std::string uDicomTag::f_GetSliceLocation(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::SliceLocation);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: SliceLocation: " << qStr;

    return mStr;
}

//SliceLocation
std::string uDicomTag::f_GetPatientPosition(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::PatientPosition);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: PatientPosition: " << qStr;

    return mStr;
}

//RescaleIntercept
std::string uDicomTag::f_GetRescaleIntercept(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RescaleIntercept);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: RescaleIntercept: " << qStr;

    return mStr;
}

//RescaleSlope
std::string uDicomTag::f_GetRescaleSlope(int pCompIdx)
{
    vtkDICOMValue pv =f_GetdicomValue(pCompIdx,DC::RescaleSlope);
    std::string mStr=pv.GetString(0);
    QString qStr=QString::fromStdString(mStr);
    //qDebug() << "uDicomTag: RescaleSlope: " << qStr;

    return mStr;
}

vtkDICOMValue uDicomTag::f_GetdicomValue(int pCompIdx,vtkDICOMTag p)
{
    uDicomTag::fileIndex  = uDicomTag::fileIndexArray->GetComponent (uDicomTag::sliceIndex,pCompIdx) ;
    uDicomTag::frameIndex = uDicomTag::fameIndexArray->GetComponent (uDicomTag::frameIndex,pCompIdx) ;
    vtkDICOMValue pv = uDicomTag::mDicomMetaData->GetAttributeValue (uDicomTag::fileIndex,
                                                                     uDicomTag::frameIndex,p);
    return  pv;
}

vtkImageData* uDicomTag::f_GetImageData(){
    return mDICOMReader->GetOutput();
}
