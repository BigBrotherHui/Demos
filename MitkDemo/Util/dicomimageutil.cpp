#pragma execution_character_set("utf-8")
#include "dicomimageutil.h"
#include <QFuture>
#include <QFutureWatcher>
#include "QMessageBox.h"
#include <QtConcurrent>
#include <fstream>
#include <iostream>

#include "itkBinaryDilateImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCoordinate.h"
#include "vtkDICOMDirectory.h"
#include "vtkDICOMReader.h"
#include "vtkImageActor.h"
#include "vtkImageDilateErode3D.h"
#include "vtkImageMapper3D.h"
#include "vtkImageProperty.h"
#include "vtkImageReader.h"
#include "vtkImageStencil.h"
#include "vtkImageStencilData.h"
#include "vtkMetaImageWriter.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkProperty.h"
#include "vtkRendererCollection.h"
#include "vtkSTLReader.h"
#include "vtkInteractorStyleImage.h"

#include "progressutil.h"
#include "systemutil.h"
#include "EmptyImageData.h"
//#include "exception.h"

#include <QApplication>
#include <QDebug>
#include <QThread>
#include <QTime>

#ifndef tr
#define tr QObject::tr
#endif

DicomImageUtil::DicomImageUtil(QObject* parent) : QObject(parent) {}

DicomImageUtil::~DicomImageUtil()
{
    qDebug() << "DicomImageUtil::~DicomImageUtil";
    if(mMyVTK)
    {
        mMyVTK->deleteLater();
    }
    qDebug() << "DicomImageUtil::~DicomImageUtil finish";
}

void DicomImageUtil::Init() {
  qDebug() << "DicomImageUtil::Init";
  mRenderWindow = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
  mMyVTK = new QVTKOpenGLNativeWidget();
  mMyVTK->setRenderWindow(mRenderWindow);
  qDebug() << "DicomImageUtil::Init finished";
}

void DicomImageUtil::SetContainer(QBoxLayout* pContainer, int pIndex) {
  qDebug() << "DicomImageUtil::SetContainer";
  if (mMyVTK == nullptr) {
    Init();
  }
  pContainer->removeItem(pContainer->takeAt(pIndex));
  pContainer->insertWidget(pIndex, mMyVTK);
  qDebug() << "DicomImageUtil::SetContainer finished";
}

void DicomImageUtil::initImageView() {
  if (mRenderWindow == nullptr) {
    Init();
  }
  mImageView = vtkSmartPointer<vtkImageViewer2>::New();
  mImageView->SetupInteractor(mRenderWindow->GetInteractor());
  mImageView->GetRenderer()->GetActiveCamera()->SetParallelProjection(1);
  mImageView->SetRenderWindow(mRenderWindow);
  mImageView->GetRenderer()->GetActiveCamera()->SetViewUp(0, -1, 0);
  mImageView->GetRenderer()->GetActiveCamera()->SetPosition(0, 0, -1);
}
void DicomImageUtil::renderImage(vtkImageData* pImage) {
  qDebug() << "DicomImageUtil::renderImage";
  initImageView();
  mImageView->SetInputData(pImage);
  if (!mEmpty) {
    getDefaultWWWL();
    resetWWWL(false);
    setSlice(0, false);
  }
  mDefaultZoom = 0;
  render();
  resetZoom(true);
  qDebug() << "DicomImageUtil::renderImage finished";
}

void DicomImageUtil::SetInteractorStyle(
    vtkInteractorObserver* pInteractorStyle) {
  mRenderWindow->GetInteractor()->SetInteractorStyle(pInteractorStyle);
}

vtkRenderWindow* DicomImageUtil::GetRenderWindow() const {
  return mRenderWindow;
}

vtkRenderer* DicomImageUtil::GetRenderer() {
  if (mImageView == nullptr) {
    initImageView();
  }
  return mImageView->GetRenderer();
}

vtkRenderWindowInteractor* DicomImageUtil::GetInteractor() const {
  return mRenderWindow->GetInteractor();
}

void DicomImageUtil::render() {
  qDebug() << "DicomImageUtil::render";
  mImageView->GetRenderer()->ResetCamera();
  mImageView->Render();
  mRenderWindow->Render();
  mImageView->GetRenderer()->GetActiveCamera()->GetFocalPoint(
      mCameraFocalPoint);
  mImageView->GetRenderer()->GetActiveCamera()->GetPosition(mCameraPosition);
  qDebug() << "DicomImageUtil::render finished";
}

int* DicomImageUtil::getDisplayPixelSize() {
    qDebug() << "DicomImageUtil::getDisplayPixelSize";
    //    if(mEmpty){
    //        qDebug() << "DicomImageUtil.getDisplayPixelSize : empty image";
    //        return nullptr;
    //    }
    double bounds[6];
    auto imageActor = mActor;
    if(imageActor==nullptr && mImageView!=nullptr)
    {
        imageActor = mImageView->GetImageActor();
    }
    if(imageActor==nullptr)
    {
        qDebug() << "DicomImageUtil.getDisplayPixelSize : no actor";
        return nullptr;
    }
    auto renderer = mRenderer;
    if(renderer==nullptr && mImageView!=nullptr)
    {
        renderer = mImageView->GetRenderer();
    }
    if(renderer==nullptr)
    {
        qDebug() << "DicomImageUtil.getDisplayPixelSize : no renderer";
        return nullptr;
    }
    imageActor->GetBounds(bounds);
    qDebug() << "DicomImageUtil.getDisplayPixelSize : Bounds:" << bounds[0] << " "
             << bounds[1] << " " << bounds[2] << " " << bounds[3] << " "
             << bounds[4] << " " << bounds[5] << " ";
    vtkNew<vtkCoordinate> coordinate;
    coordinate->SetCoordinateSystemToWorld();
    coordinate->SetValue(bounds[1], bounds[3], bounds[5]);
    int* point = coordinate->GetComputedDisplayValue(renderer);
    qDebug() << "DicomImageUtil.getDisplayPixelSize : Compute display value 135"
             << point[0] << " " << point[1] << " " << point[2] << " ";
    vtkNew<vtkCoordinate> coordinate1;
    coordinate1->SetCoordinateSystemToWorld();
    coordinate1->SetValue(bounds[0], bounds[2], bounds[4]);
    int* point1 = coordinate1->GetComputedDisplayValue(renderer);
    qDebug() << "DicomImageUtil.getDisplayPixelSize : Compute display value 024"
             << point1[0] << " " << point1[1] << " " << point1[2] << " ";
    int* displaySize =
        new int[2]{labs(point[0] - point1[0]), labs(point[1] - point1[1])};
    qDebug() << "DicomImageUtil::getDisplayPixelSize finished : Display size:"
             << displaySize[0] << " " << displaySize[1];
    return displaySize;
}

double DicomImageUtil::getZoom() {
  qDebug() << "DicomImageUtil::getZoom";
  if (mEmpty) {
    qDebug() << "DicomImageUtil.getZoom : empty image";
    return 0;
  }
  int* displaySize = getDisplayPixelSize();
  int realSize[2];
  int orientation = mImageView->GetSliceOrientation();
  if (orientation == vtkImageViewer2::SLICE_ORIENTATION_YZ) {
    realSize[0] = (mExtent[3] - mExtent[2]);
    realSize[1] = (mExtent[5] - mExtent[4]);
  } else {
    realSize[0] = (mExtent[1] - mExtent[0]);
    if (orientation == vtkImageViewer2::SLICE_ORIENTATION_XY) {
      realSize[1] = (mExtent[3] - mExtent[2]);
    } else {
      realSize[1] = (mExtent[5] - mExtent[4]);
    }
  }
  qDebug() << "DicomImageUtil::getZoom finished"
           << "Extents:" << mExtent[0] << " " << mExtent[1] << " " << mExtent[2]
           << " " << mExtent[3] << " " << mExtent[4] << " " << mExtent[5] << " "
           << "Real size:" << realSize[0] << " " << realSize[1];
  return (displaySize[0] + 0.0) / (realSize[0] + 0.0);
}

void DicomImageUtil::getDefaultWWWL(){
    qDebug() << "DicomImageUtil::getDefaultWWWL";
    vtkDICOMValue windowWidth = mDicomMetaData->GetAttributeValue(0, DC::WindowWidth);
    vtkDICOMValue windowCenter = mDicomMetaData->GetAttributeValue(0, DC::WindowCenter);

    if(windowWidth.IsValid() && windowCenter.IsValid()){
        mDefaultWW = windowWidth.GetInt(0);
        mDefaultWL = windowCenter.GetInt(0);
    }else{
        mDefaultWW = 1500;
        mDefaultWL = 300;
        /*
        short minPixel;
        short maxPixel;
        bool start = true;
        for(int x=mExtent[0];x<=mExtent[1];x++){
            for(int y=mExtent[2];y<=mExtent[3];y++){
                int z=mExtent[4];
                short* scalar = (short *)mImageData->GetScalarPointer(x, y, z) ;
                if(start || *scalar>maxPixel){
                    maxPixel = *scalar;
                }
                if(start || *scalar<minPixel){
                    minPixel = *scalar;
                }
                start=false;
            }
        }
        qDebug() << "DicomImageUtil.getDefaultWWWL : minPixel: " <<  minPixel << " maxPixel: " << maxPixel;
        mDefaultWW = maxPixel - minPixel;
        mDefaultWL = minPixel + mDefaultWW / 2;
        */
    }
    qDebug() << "DicomImageUtil::getDefaultWWWL finished: default WW: " <<  mDefaultWW << " default WL: " << mDefaultWL;
}

void DicomImageUtil::setWWWL(int pWW, int pWL, bool pRender){
    qDebug() << "DicomImageUtil::setWWWL: ww="<<pWW<<" wl="<<pWL;
    mWW = pWW;
    mWL = pWL;
    mImageView->SetColorLevel(mWL);
    mImageView->SetColorWindow(mWW);
    if(pRender){
        mImageView->Render();
    }
    qDebug() << "DicomImageUtil::setWWWL finished";
}

void DicomImageUtil::SetWWWL(int pWW, int pWL){
    setWWWL(pWW, pWL, true);
}

void DicomImageUtil::resetWWWL(bool pRender){
    qDebug() << "DicomImageUtil : rset ww&wl";
    setWWWL(mDefaultWW, mDefaultWL, pRender);
}

void DicomImageUtil::ResetWWWL(){
    resetWWWL(true);
}

void DicomImageUtil::zoom(double pRate, bool pRender){
    qDebug() << "DicomImageUtil::zoom: rate="<<pRate;
    vtkRenderer *renderer = mRenderer;
    if(renderer==nullptr && mImageView!=nullptr)
    {
        renderer = mImageView->GetRenderer();
    }
    if(renderer==nullptr)
    {
        qDebug() << "DicomImageUtil::zoom renderer is null";
        return;
    }
    qDebug() << "DicomImageUtil::zoom 1";
    double zoomRate;
    if(pRate<=0)
    {
        qDebug() << "DicomImageUtil::zoom: zoom to full:"<<mMyVTK->width()<<mMyVTK->height();
        int* displaySize = getDisplayPixelSize();
        double rate1 = mMyVTK->width()/(displaySize[0]+0.0);
        double rate2 = mMyVTK->height()/(displaySize[1]+0.0);
        if(rate1>rate2)
        {
            zoomRate = rate2;
        }
        else
        {
            zoomRate = rate1;
        }
    }
    else
    {
        double currentRate = getZoom();
        qDebug() << "DicomImageUtil::zoom: Current zoom rate=" << currentRate;
        zoomRate = pRate/currentRate;
    }
    qDebug() << "DicomImageUtil::zoom: zoomRate="<<zoomRate;
    vtkCamera *camera =renderer->GetActiveCamera();
    if(camera==nullptr)
    {
        qDebug() << "DicomImageUtil::zoom camera is null";
        return;
    }
    camera->Zoom(zoomRate);
    if(pRender){
//        mImageView->Render();
//        renderer->Render();
        mRenderWindow->Render();
    }
    qDebug() << "DicomImageUtil::zoom finished";
}

void DicomImageUtil::Zoom(double pRate){
    zoom(pRate, true);
}

void DicomImageUtil::resetZoom(bool pRender){
    qDebug() << "DicomImageUtil : reset zoom";
    zoom(mDefaultZoom, pRender);
}

void DicomImageUtil::ResetZoom(){
    resetZoom(true);
}

void DicomImageUtil::Reset(){
    qDebug() << "DicomImageUtil::Reset";
    mImageView->GetRenderer()->GetActiveCamera()->SetFocalPoint(mCameraFocalPoint);
    mImageView->GetRenderer()->GetActiveCamera()->SetPosition(mCameraPosition);
    resetWWWL(false);
    resetZoom(false);
    mRenderWindow->Render();
    qDebug() << "DicomImageUtil::Reset finished";
}

void DicomImageUtil::RequestUpdate()
{
    render();
}

double *DicomImageUtil::GetOrigin()
{
    return mOrigin;
}

double *DicomImageUtil::GetCenter()
{
    return mCenter;
}

double *DicomImageUtil::GetSpacing()
{
    return mSpacing;
}

int *DicomImageUtil::GetExtent()
{
    return mExtent;
}

int *DicomImageUtil::GetDims()
{
    return mDims;
}

void DicomImageUtil::LoadEmptyDicomImage(){
    qDebug() << "DicomImageUtil::LoadEmptyDicomImage";
    mEmpty=true;
    vtkSmartPointer< EmptyImageData > emptyImageData = vtkSmartPointer< EmptyImageData >(EmptyImageData::New());
    renderImage(emptyImageData->f_GetEmptyImageData());
    qDebug() << "DicomImageUtil::LoadEmptyDicomImage finished";
}

vtkStringArray* DicomImageUtil::GetFiles(QString pFolder, QString &pErrorMessage, int pSeriesIndex)
{
    pFolder = SystemUtil::GetFullDir(pFolder);
    qDebug() << "DicomImageUtil::GetFiles : get files from " << pFolder;
    if(!SystemUtil::CheckDir(pFolder, tr("图像文件夹"), false)){
        pErrorMessage = tr("加载图像失败，图像文件夹不存在！");
        return nullptr;
    }
    vtkNew<vtkDICOMDirectory> reader;
    reader->SetDirectoryName(pFolder.toStdString().c_str());
    reader->SetScanDepth(3);
    qDebug() << "DicomImageUtil::GetFiles: update reader";
    reader->Update();
    if(reader->GetNumberOfSeries()<=0){
        pErrorMessage = tr("加载图像失败，文件夹中没有找到图像！");
        return nullptr;
    }
    vtkStringArray* files = vtkStringArray::New();
    qDebug() << "DicomImageUtil::GetFiles: deep copy files";
    files->DeepCopy(reader->GetFileNamesForSeries(pSeriesIndex));
    qDebug() << "DicomImageUtil::GetFiles: return files";
    return files;
}

bool DicomImageUtil::IsInvert(vtkStringArray *pDicomSeriesFiles)
{
    qDebug()<<"IsInvert start";
    if(pDicomSeriesFiles==nullptr || pDicomSeriesFiles->GetSize()<=1)
    {
        return false;
    }
    readMetaData(pDicomSeriesFiles);
    vtkDICOMValue imageOrientation = mDicomMetaData->GetAttributeValue (0,DC::ImageOrientationPatient);
    vtkDICOMValue imagePostion = mDicomMetaData->GetAttributeValue (0,DC::ImagePositionPatient);
    vtkDICOMValue imagePostion1 = mDicomMetaData->GetAttributeValue (1,DC::ImagePositionPatient);
    vtkDICOMValue instanceNumber = mDicomMetaData->GetAttributeValue(0, DC::InstanceNumber);
    vtkDICOMValue instanceNumber1 = mDicomMetaData->GetAttributeValue(1, DC::InstanceNumber);
    double p1,p2;
    if(imageOrientation.GetDouble(1)!=0&&imageOrientation.GetDouble(5)!=0){
        p1=imagePostion.GetDouble(0);
        p2=imagePostion1.GetDouble(0);
    }else if(imageOrientation.GetDouble(0)!=0&&imageOrientation.GetDouble(5)!=0){
        p1=imagePostion.GetDouble(1);
        p2=imagePostion1.GetDouble(1);
    }else{
        p1=imagePostion.GetDouble(2);
        p2=imagePostion1.GetDouble(2);
    }
    int in1 = instanceNumber.GetInt(0), in2 = instanceNumber1.GetInt(0);
    qDebug()<<"IsInvert start";
    return (in1>in2&&p1>p2) || (in1<in2&&p1<p2);
}

bool DicomImageUtil::CheckSeries(vtkStringArray *pDicomSeriesFiles)
{
    qDebug()<<"DicomImageUtil::CheckSeries";
    //读取dicom元信息
    readMetaData(pDicomSeriesFiles);
    auto instanceCount = mDicomMetaData->GetNumberOfInstances();
    if(instanceCount<=0)
    {
        return false;
    }
    //获取层厚
    auto sliceThickness = mDicomMetaData->GetAttributeValue (0,DC::SliceThickness).GetDouble(0);
    //遍历判断是否缺失数据
    auto sliceLocation = mDicomMetaData->GetAttributeValue (0,DC::SliceLocation).GetDouble(0);
    for(int i=1;i<instanceCount;i++)
    {
        auto sliceLocation1 = mDicomMetaData->GetAttributeValue(i,DC::SliceLocation).GetDouble(0);
        //判断是否缺帧
        if(abs(sliceLocation1 - sliceLocation)>sliceThickness*1.5)
        {
            return false;
        }
        sliceLocation = sliceLocation1;
    }
    qDebug()<<"DicomImageUtil::CheckSeries finish";
    return true;
}

void DicomImageUtil::readDicom(vtkStringArray *pDicomSeriesFiles)
{
    qDebug() << "DicomImageUtil::loadDicom";
    // 处理
    mDICOMReader->SetFileNames(pDicomSeriesFiles);
    mDICOMReader->SetMemoryRowOrderToFileNative();
    mDICOMReader->SortingOn();
    mDICOMReader->Update();
    qDebug() << "DicomImageUtil::loadDicom finish";
}

void DicomImageUtil::readMetaData(vtkStringArray *pDicomSeriesFiles)
{
    qDebug() << "DicomImageUtil::readMetaData";
    readDicom(pDicomSeriesFiles);
    mDicomMetaData = mDICOMReader->GetMetaData();
    qDebug() << "DicomImageUtil::readMetaData finish";
}

void DicomImageUtil::GetDicomSeries(vtkStringArray *pDicomSeriesFiles)
{
    qDebug() << "DicomImageUtil::GetDicomSeries";
    readMetaData(pDicomSeriesFiles);
    mDICOMReader->GetOutput()->GetDimensions(mDims);
    mDICOMReader->GetOutput()->GetExtent(mExtent);
    mDICOMReader->GetOutput()->GetOrigin(mOrigin);
    mDICOMReader->GetOutput()->GetSpacing(mSpacing);
    mCenter[0] = mOrigin[0] + mSpacing[0] * 0.5 * (mExtent[0] + mExtent[1]);
    mCenter[1] = mOrigin[1] + mSpacing[1] * 0.5 * (mExtent[2] + mExtent[3]);
    mCenter[2] = mOrigin[2] + mSpacing[2] * 0.5 * (mExtent[4] + mExtent[5]);
    //read ww&wc
    vtkDICOMValue windowWidth = mDicomMetaData->GetAttributeValue(0, DC::WindowWidth);
    vtkDICOMValue windowCenter = mDicomMetaData->GetAttributeValue(0, DC::WindowCenter);
    if(windowWidth.IsValid() && windowCenter.IsValid()){
        mWW = windowWidth.GetInt(0);
        mWL = windowCenter.GetInt(0);
    }
    else
    {
        short minPixel;
        short maxPixel;
        bool start = true;
        for(int x=mExtent[0];x<=mExtent[1];x++){
            for(int y=mExtent[2];y<=mExtent[3];y++){
                int z=mExtent[4];
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
        mWW = maxPixel - minPixel;
        mWL = minPixel + mWW / 2;
        mDefaultWW = maxPixel - minPixel;
        mDefaultWL = minPixel + mWW / 2;
    }
    qDebug() << "DicomImageUtil::GetDicomSeries finish";
}

int DicomImageUtil::LoadDicomImage(QString pFolder){
    return doLoadDicomImage([=]{
        QString errorMessage;
        return GetFiles(pFolder, errorMessage);
    });
}

int DicomImageUtil::LoadSingleDicomImage(QString pFilePath){
    return doLoadDicomImage([=]{
        vtkStringArray* files = vtkStringArray::New();
        files->InsertValue(0, pFilePath.toLocal8Bit());
        return files;
    });
}

int DicomImageUtil::LoadDicomImage(vtkStringArray* pFiles){
    return doLoadDicomImage([=]{
        return pFiles;
    });
}
int DicomImageUtil::doLoadDicomImage(std::function<vtkStringArray* ()> pCallback){
//    MY_TRY_CATCH("SurgeryNavigation::onToolButtonClick",
    qDebug() << "DicomImageUtil::LoadDicomImage";
    ProgressUtil progressUtil;
    progressUtil.Start(tr("预览图像"));
    progressUtil.SetTips(tr("正在读取图像"));
    vtkStringArray* files=nullptr;
    progressUtil.AutoStep([&]()mutable{
        files = pCallback();
    }, 30, tr("正在加载图像"));
    if(files==nullptr){
        progressUtil.Finish();
        return -1;
    }
    qDebug() << "DicomImageUtil::LoadDicomImage : load image from files, image count="<<files->GetNumberOfValues();
    QApplication::processEvents();
    mEmpty=false;
    progressUtil.AutoStep([=]()mutable{
        GetDicomSeries(files);
    }, 50);
    progressUtil.AutoStep([=]()mutable{
        mImageData=mDICOMReader->GetOutput();
    }, 70);
    progressUtil.SetTips(tr("正在显示图像"));
    renderImage(mImageData);
    progressUtil.Finish();
    qDebug() << "DicomImageUtil::LoadDicomImage finished";
    return mDims[2];
//    )
    return -1;
}

int DicomImageUtil::LoadRawImage(QString pRawDataFile, int *pDims, double pOpacity, int pDataScalarType)
{
    qDebug()<<"DicomImageUtil::LoadRawImage"<<pRawDataFile;
    if(pRawDataFile.isEmpty())
    {
        return -1;
    }
    //raw data reader
//    if(mImageView==nullptr)
//    {
//        initImageView();
//    }
    vtkSmartPointer<vtkImageReader> reader = vtkSmartPointer<vtkImageReader>::New();
    reader->SetFileName(pRawDataFile.toStdString().data());
    reader->SetFileDimensionality(3);
    reader->SetDataScalarType(pDataScalarType);
    reader->SetDataExtent(0, pDims[0]-1, 0, pDims[1]-1, 0, 0);
    reader->SetNumberOfScalarComponents(1);
    reader->Update();
    if(mImageBlend==nullptr)
    {
        mImageBlend = vtkSmartPointer< vtkImageBlend>::New();
    }
//    mImageBlend->AddInputConnection(reader->GetOutputPort());
    mImageBlend->AddInputData(reader->GetOutput());
    if(mRenderer==nullptr)
    {
        mRenderer = vtkSmartPointer<vtkRenderer>::New();
        mRenderWindow->AddRenderer(mRenderer);
    }
    if(mActor==nullptr)
    {
        mActor = vtkSmartPointer<vtkImageActor>::New();
        mActor->GetProperty()->SetColorLevel(2047);
        mActor->GetProperty()->SetColorWindow(4096);
        mRenderer->AddActor(mActor);
    }
    mActor->GetMapper()->SetInputConnection(mImageBlend->GetOutputPort());
//    mRenderer->Render();

//    mImageView->SetColorLevel(2047);
//    mImageView->SetColorWindow(4096);
//    mImageView->SetInputData(mImageBlend->GetOutput());
    auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    auto style = vtkSmartPointer<vtkInteractorStyleImage>::New();
    interactor->SetInteractorStyle(style);
    interactor->SetRenderWindow(mRenderWindow);
    interactor->Initialize();
    int num = mImageBlend->GetNumberOfInputs();
    mImageBlend->SetOpacity(num, pOpacity);
    mImageBlend->Update();
//    mRenderWindow->Render();
    qDebug()<<"DicomImageUtil::LoadRawImage finish"<<num<<pOpacity;
    return num;
}

void DicomImageUtil::SetOpacity(int pNum, double pOpacity){
    qDebug() << "DicomImageUtil : SetOpacity"<<pNum<<pOpacity;
    if(mImageBlend!=nullptr)
    {
        mImageBlend->SetOpacity(pNum, pOpacity);
        mImageBlend->Update();
//        mRenderer->Render();
//        mRenderWindow->Render();
    }
}

void DicomImageUtil::ClearRawImage(){
    qDebug() << "DicomImageUtil : ClearRawImage";
    if(mImageBlend!=nullptr)
    {
        mImageBlend->RemoveAllInputs();
        mImageBlend->Update();
        mImageBlend = nullptr;
//        mRenderer->Render();
//        mRenderWindow->Render();
//        qDebug()<<"num = "<< mImageBlend->GetNumberOfInputs();
    }
    qDebug() << "DicomImageUtil : ClearRawImage finish";
}

void DicomImageUtil::setSlice(int pSlice, bool pRender){
    qDebug() << "DicomImageUtil : set slice to "<<pSlice;
    mImageView->SetSlice(pSlice);
    if(pRender){
        mImageView->Render();
        mRenderWindow->Render();
    }
}

void DicomImageUtil::SetSlice(int pSlice){
    setSlice(pSlice, true);
}

mitk::Vector3D DicomImageUtil::GetSpacing(vtkStringArray *pFiles)
{
    qDebug()<<"DicomImageUtil::GetSpacing";
    vtkNew<vtkDICOMReader> reader;
    reader->SetFileNames(pFiles);
    reader->Update();
    double spacing[3];
    reader->GetOutput()->GetSpacing(spacing);
    mitk::Vector3D spacing3D;
    spacing3D.FillVector(spacing);
    qDebug()<<"DicomImageUtil::GetSpacing finish";
    return spacing3D;
}

bool DicomImageUtil::StlToMask(QString pStlFile, QString pInfoFile, QString pRawFile, double *pImageCenter, double pSpacing, int pInValue, int pOutValue, int pDilateValue)
{
    return StlToMask(pStlFile, pInfoFile, pRawFile, pImageCenter, new double[3]{pSpacing, pSpacing, pSpacing}, pInValue, pOutValue, pDilateValue);
}

bool DicomImageUtil::StlToMask(QString pStlFile, QString pInfoFile, QString pRawFile, double* pImageCenter,
                               double* pSpacing, int pForegroundValue, int pBackgroundValue, int pDilateValue)
{
    qDebug()<<QTime::currentTime()<<"DicomImageUtil::StlToMask"<<pStlFile<<pInfoFile<<pRawFile;
    //read poly data from stl
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(pStlFile.toStdString().data());
    reader->Update();
    auto polyData = reader->GetOutput();
    //convert poly data to image data
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    double* bounds = polyData->GetBounds();
    imageData->SetSpacing(pSpacing);
    // compute dimensions,extents,origin
    int dim[3];
    int extent[6];
    double origin[3];
    for (int i = 0; i < 3; i++)
    {
        dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / pSpacing[i]));
        extent[i*2] = 0;
        extent[i*2 + 1] = dim[i] - 1;
        origin[i] = bounds[i*2] + pSpacing[i] / 2;
    }
    imageData->SetDimensions(dim);
    imageData->SetExtent(extent);

    imageData->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
    whiteImage->SetScalarTypeToShort();
    whiteImage->AllocateScalars();
#else
    imageData->AllocateScalars(VTK_SHORT, 1);
#endif
    // fill the image with foreground voxels:
    vtkIdType count = imageData->GetNumberOfPoints();
    for (vtkIdType i = 0; i < count; ++i) {
      imageData->GetPointData()->GetScalars()->SetTuple1(i, pForegroundValue);
    }

    // polygonal data --> image stencil:
    vtkSmartPointer<vtkPolyDataToImageStencil> pdtoImageStencil =
        vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    pdtoImageStencil->SetInputData(polyData);
    pdtoImageStencil->SetOutputOrigin(origin);
    pdtoImageStencil->SetOutputSpacing(pSpacing);
    pdtoImageStencil->SetOutputWholeExtent(extent);
    pdtoImageStencil->Update();

    //cut the corresponding white image and set the background:
    vtkSmartPointer<vtkImageStencil> imageStencil = vtkSmartPointer<vtkImageStencil>::New();
    imageStencil->SetInputData(imageData);
    imageStencil->SetStencilConnection(pdtoImageStencil->GetOutputPort());
    imageStencil->ReverseStencilOff();
    imageStencil->SetBackgroundValue(pBackgroundValue);
    imageStencil->Update();
    imageData->DeepCopy(imageStencil->GetOutput());
    if(pDilateValue>0)
    {
        imageData=DilateImage(imageData, pDilateValue, pForegroundValue, pBackgroundValue);
    }
    //save image data to file
    vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
    writer->SetFileName(pInfoFile.toStdString().data());
    writer->SetRAWFileName(pRawFile.toStdString().data());
    writer->SetInputData(imageData);
    writer->SetCompression(false);
    writer->Write();
    //add raw center position via image center
    double* center = polyData->GetCenter();
    double* position = new double[3];
    for(int i=0;i<3;i++)
    {
        position[i] = center[i] - pImageCenter[i];
    }
    std::ofstream file;
    file.open(pInfoFile.toStdString().data(), std::ios::out|std::ios::app);
    file<<"CenterToImage =";
    for(int i=0;i<3;i++)
    {
        file<<" "<<position[i];
    }
    file<<std::endl;
    file.close();
    qDebug()<<QTime::currentTime()<<"DicomImageUtil::StlToMask finish";
    return true;
}

bool DicomImageUtil::isInside(vtkImageStencilData* pStencilImage, int x, int y, int z, int pDilateValue)
{
    if(pStencilImage->IsInside(x, y, z))
    {
        return true;
    }
    for(int x1 = x - pDilateValue;x1<x+pDilateValue;x1++)
    {
        for(int y1 = y - pDilateValue;y1<y+pDilateValue;y1++)
        {
            for(int z1 = z - pDilateValue;z1<z+pDilateValue;z1++)
            {
                if(sqrt(pow(x1-x,2)+pow(y1-y,2)+pow(z1-z,2))>pDilateValue)
                {
                    continue;
                }
                if(pStencilImage->IsInside(x1, y1, z1))
                {
                    return true;
                }
            }
        }
    }
    return false;
}

vtkSmartPointer<vtkImageData> DicomImageUtil::maskImage(vtkSmartPointer<vtkImageData> pImageData, vtkSmartPointer<vtkImageData> pNewImageData, vtkSmartPointer<vtkPolyData> pPolyData, double* pImageOrigin, double *pImageCenter, int pBackgroundValue, int pDilateValue)
{
    qDebug()<<"DicomImageUtil::maskImage"<<pBackgroundValue<<pDilateValue;
    //获取图像信息
    auto spacing = pImageData->GetSpacing();
    auto imageExtent = pImageData->GetExtent();
    auto newImageOrigin = pNewImageData->GetOrigin();
    auto imageOrigin = pImageOrigin;
    //polygonal data --> image stencil:
    //获取stl信息
    double* stlBounds = pPolyData->GetBounds();
    int stlDim[3],stlExtent[6],offset[3],offsetNew[3];
    int flag[3];
    double stlOrigin[3];
    for(int i=0;i<3;i++)
    {
        stlDim[i] = ceil((stlBounds[i*2+1]-stlBounds[i*2]) / spacing[i]);
        if(pImageOrigin[i]>pImageCenter[i])
        {
            flag[i] = -1;
        }
        else
        {
            flag[i] = 1;
        }
    }
    if(stlDim[1]>stlDim[0])
    {
        stlDim[0] = stlDim[1];
    }
    else
    {
        stlDim[1] = stlDim[0];
    }
    for(int i=0;i<3;i++)
    {
        stlOrigin[i] = stlBounds[i*2] + spacing[i] / 2;
        stlExtent[i*2] = 0;
        stlExtent[i*2+1] = stlDim[i] - 1;
        offset[i] = ceil((stlOrigin[i] - imageOrigin[i]) / spacing[i]);
        offsetNew[i] = ceil((stlOrigin[i] - newImageOrigin[i]) / spacing[i]);
    }
    //convet poly data to image stencil
    vtkSmartPointer<vtkPolyDataToImageStencil> polyDatatoImageStencil = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    polyDatatoImageStencil->SetInputData(pPolyData);
    polyDatatoImageStencil->SetOutputOrigin(stlOrigin);
    polyDatatoImageStencil->SetOutputSpacing(spacing);
    polyDatatoImageStencil->SetOutputWholeExtent(stlExtent);
    polyDatatoImageStencil->Update();
    auto stencilImage = polyDatatoImageStencil->GetOutput();
    vtkSmartPointer<vtkImageData> newImageData;
    if(pDilateValue<0)
    {
        //cut image by extent
        int cutExtent[6];
        for(int i=0;i<3;i++)
        {
            if(flag[i]==1)
            {
                cutExtent[i*2] = offset[i];
                cutExtent[i*2+1] = offset[i] + stlDim[i] - 1;
            }
            else
            {
                cutExtent[i*2] = -offset[i] - stlDim[i] + 1;
                cutExtent[i*2+1] = -offset[i];
            }
        }
        newImageData = CutImage(pImageData, cutExtent);
        newImageData->SetOrigin(stlOrigin);
//        //cut the corresponding white image and set the background:
//        vtkSmartPointer<vtkImageStencil> imageStencil = vtkSmartPointer<vtkImageStencil>::New();
//        imageStencil->SetInputData(newImageData);
//        imageStencil->SetStencilData(stencilImage);
//        imageStencil->ReverseStencilOff();
//        imageStencil->SetBackgroundValue(pBackgroundValue);
//        imageStencil->Update();
//        newImageData->DeepCopy(imageStencil->GetOutput());
//        newImageData->SetOrigin(stlOrigin);
    }
    else
    {
        //cut image by stencil with dilate
        auto newExtent = new int[6];
        auto newOrigin = new double[6];
        for(int i=0;i<3;i++)
        {
            newExtent[2*i] = stlExtent[2*i];
            newExtent[2*i+1] = stlExtent[2*i+1]+2*pDilateValue;
            newOrigin[i] = stlOrigin[i] - pDilateValue*spacing[i];
        }
        newImageData = vtkSmartPointer<vtkImageData>::New();
        newImageData->SetSpacing(spacing);
        newImageData->SetExtent(newExtent);
        newImageData->SetOrigin(newOrigin);
        auto scalarType = pImageData->GetScalarType();
        auto numComponents = pImageData->GetNumberOfScalarComponents();
        newImageData->AllocateScalars(scalarType, numComponents);
        //使用多线程并发进行处理，提高处理速度
        qDebug()<<"DicomImageUtil::maskImage create process";
        mProcessCount = 0;
        QFutureWatcher<void> *mWatcher;
        mWatcher = new QFutureWatcher<void>();
        for(int x=stlExtent[0] - pDilateValue;x<=stlExtent[1] + pDilateValue;x++)
        {
            mMutex.lock();
            mProcessCount++;
            mMutex.unlock();
            QFuture<void> mFuture = QtConcurrent::run([=](int x) mutable{
                try{
                    for(int y=stlExtent[2] - pDilateValue;y<=stlExtent[3] + pDilateValue;y++)
                    {
                        for(int z=stlExtent[4] - pDilateValue;z<=stlExtent[5] + pDilateValue;z++)
                        {
                            int cood[3] = {x, y ,z};
                            int newCood[3], imageCood[3],newImageCood[3];
                            for(int i=0;i<3;i++)
                            {
                                newCood[i] = cood[i] + pDilateValue;
                                newImageCood[i] = cood[i] + offsetNew[i];
                                imageCood[i] = cood[i] + offset[i];
                                if(flag[i]==-1)
                                {
                                    imageCood[i] = - imageCood[i];
                                }
                            }
                            short *newPixelValue = (short *)newImageData->GetScalarPointer(newCood);
                            if(newPixelValue==nullptr)
                            {
                                qDebug()<<"DicomImageUtil::maskImage new pixel is null"<<newCood[0]<<newCood[1]<<newCood[2];
                                continue;;
                            }
                            if(imageCood[0]<imageExtent[0] || imageCood[0]>imageExtent[1] ||
                                    imageCood[1]<imageExtent[2] || imageCood[1]>imageExtent[3] ||
                                    imageCood[2]<imageExtent[4] || imageCood[2]>imageExtent[5])
                            {
                                *newPixelValue = pBackgroundValue;
                            }
                            else if(isInside(stencilImage, x, y, z, pDilateValue))
                            {
                                short *pixelValue = (short *)pImageData->GetScalarPointer(imageCood);
                                if(pixelValue==nullptr)
                                {
                                    qDebug()<<"DicomImageUtil::maskImage pixel is null"<<imageCood[0]<<imageCood[1]<<imageCood[2];
                                    continue;;
                                }
                                *newPixelValue = *pixelValue;
                                short *pixelValueNew = (short *)pNewImageData->GetScalarPointer(newImageCood);
                                if(pixelValueNew==nullptr)
                                {
                                    qDebug()<<"DicomImageUtil::maskImage new image pixel is null"<<newImageCood[0]<<newImageCood[1]<<newImageCood[2];
                                    continue;;
                                }
                                *pixelValueNew = *pixelValue;
                            }
                            else
                            {
                                *newPixelValue = pBackgroundValue;
                            }
                        }
                    }
                }  catch (std::exception e) {
                    qDebug()<<"exception:"<<e.what();
                }  catch (...) {
                    qDebug()<<"exception...";
                }
                mMutex.lock();
                mProcessCount--;
                mMutex.unlock();
            }, x);
            mWatcher->setFuture(mFuture);
        }
        qDebug()<<"DicomImageUtil::maskImage check process";
        while(mProcessCount>0){
            int waitMs = 1000;
            QTime dieTime = QTime::currentTime().addMSecs(waitMs);
            while( QTime::currentTime() < dieTime )
            {
                QApplication::processEvents(QEventLoop::AllEvents, waitMs);
            }
        }
    }
    qDebug()<<"DicomImageUtil::maskImage finish";
    return newImageData;
}

bool DicomImageUtil::ImageMask(vtkSmartPointer<vtkImageData> pImageData, QStringList pStlFiles, double* pImageOrigin,
                               double* pImageCenter, QStringList pInfoFiles, QStringList pRawFiles,
                               QString pInfoFile, QString pRawFile, int pBackgroundValue, int pDilateValue,
                               std::function<void (int pStep)> pCallback)
{
    qDebug()<<QTime::currentTime()<<"DicomImageUtil::ImageMask";
    QElapsedTimer elapsedTimer;
    elapsedTimer.start();
    //TODO 有个奇怪的现象，在这里用pImageData的GetOrigin和GetCenter返回的结果都是不对，因此需要作为函数参数传递进来
    auto spacing = pImageData->GetSpacing();
    double newBounds[6];
    int stepReadStl = 20, stepInitImage = 10, stepMask = 60;
    int step = 0;
    QList<vtkSmartPointer<vtkPolyData>> polyDatas;
    for(int i=0;i<pStlFiles.size();i++)
    {
        //read poly data from stl
        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(pStlFiles.at(i).toStdString().data());
        reader->Update();
        auto polyData = reader->GetOutput();
        polyDatas.append(polyData);
        double* stlBounds = polyData->GetBounds();
        if(i==0)
        {
            for(int j=0;j<6;j++)
            {
                newBounds[j] = stlBounds[j];
            }
        }
        else
        {
            for(int j=0;j<3;j++)
            {
                if(stlBounds[j*2]<newBounds[j*2])
                {
                    newBounds[j*2] = stlBounds[j*2];
                }
                if(stlBounds[j*2+1]>newBounds[j*2+1])
                {
                    newBounds[j*2+1] = stlBounds[j*2+1];
                }
            }
        }
    }
    step += stepReadStl;
    if(pCallback)
    {
        pCallback(step);
    }
    double newOrigin[3];
    int newDim[3],newExtent[6];
    for(int i=0;i<3;i++)
    {
        newDim[i] = ceil((newBounds[i*2+1]-newBounds[i*2]) / spacing[i]);
        if(pDilateValue>=0)
        {
            newDim[i] += pDilateValue*2;
        }
    }
    if(newDim[1]>newDim[0])
    {
        newDim[0] = newDim[1];
    }
    else
    {
        newDim[1] = newDim[0];
    }
    for(int i=0;i<3;i++)
    {
        newOrigin[i] = newBounds[i*2] + spacing[i]/2;
        if(pDilateValue>=0)
        {
            newOrigin[i] += - pDilateValue*spacing[i];
            newExtent[i*2] = 0;
            newExtent[i*2 + 1] = newDim[i] - 1;
        }
        else
        {
            newExtent[i*2] = (newOrigin[i] - pImageOrigin[i])/spacing[i];
            newExtent[i*2 + 1] = newExtent[i*2] + newDim[i] - 1;
        }
    }
    vtkSmartPointer<vtkImageData> newImageData;
    if(pDilateValue<0)
    {
        newImageData = CutImage(pImageData, newExtent);
        newImageData->SetOrigin(newOrigin);
    }
    else
    {
        newImageData = vtkSmartPointer<vtkImageData>::New();
        newImageData->SetSpacing(spacing);
        newImageData->SetExtent(newExtent);
        newImageData->SetOrigin(newOrigin);
        auto scalarType = pImageData->GetScalarType();
        auto numComponents = pImageData->GetNumberOfScalarComponents();
        newImageData->AllocateScalars(scalarType, numComponents);
        for(int i=newExtent[0];i<=newExtent[1];i++)
        {
            for(int j=newExtent[2];j<=newExtent[3];j++)
            {
                for(int k=newExtent[4];k<=newExtent[5];k++)
                {
                    short *pixel = (short*) newImageData->GetScalarPointer (i,j,k);
                    *pixel = pBackgroundValue;
                }
            }
        }
    }
    step += stepInitImage;
    if(pCallback)
    {
        pCallback(step);
    }
    auto size = polyDatas.size();
    auto stepPerImage = stepMask / size;
    auto stepPerMaskImage = stepPerImage * 4 / 5;
    auto stepLast = stepPerImage - stepPerMaskImage;
    for(int i=0;i<size;i++)
    {
//        //read poly data from stl
//        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
//        reader->SetFileName(pStlFiles.at(i).toStdString().data());
//        reader->Update();
//        auto polyData = reader->GetOutput();
        auto polyData = polyDatas.at(i);
        auto polyCenter = polyData->GetCenter();
        //mask image with poly data
        auto imageData = maskImage(pImageData, newImageData, polyData, pImageOrigin, pImageCenter, pBackgroundValue, pDilateValue);
        if(imageData == nullptr)
        {
            return false;
        }
        step += stepPerMaskImage;
        if(pCallback)
        {
            pCallback(step);
        }
        //save image data to file
        vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
        writer->SetFileName(pInfoFiles.at(i).toStdString().data());
        writer->SetRAWFileName(pRawFiles.at(i).toStdString().data());
        writer->SetCompression(false);
        writer->SetInputData(imageData);
        writer->Write();
        //add raw center position via image center
        double* center = imageData->GetCenter();
        double* imageCenter = pImageCenter;
        double* position = new double[3];
        double* positionPoly = new double[3];
        for(int i=0;i<3;i++)
        {
            position[i] = center[i] - imageCenter[i];
            positionPoly[i] = polyCenter[i] - imageCenter[i];
        }
        std::ofstream file;
        file.open(pInfoFiles.at(i).toStdString().data(), std::ios::out|std::ios::app);
        file<<"CenterToImage =";
        for(int i=0;i<3;i++)
        {
            file<<" "<<position[i];
        }
        file<<std::endl;
        file<<"SegmentCenterToImage =";
        for(int i=0;i<3;i++)
        {
            file<<" "<<positionPoly[i];
        }
        file<<std::endl;
        file.close();
        step += stepLast;
        if(pCallback)
        {
            pCallback(step);
        }
    }
    //save image data to file
    vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
    writer->SetFileName(pInfoFile.toStdString().data());
    writer->SetRAWFileName(pRawFile.toStdString().data());
    writer->SetCompression(false);
    writer->SetInputData(newImageData);
    writer->Write();
    double* center = newImageData->GetCenter();
    double* imageCenter = pImageCenter;
    double* position = new double[3];
    for(int i=0;i<3;i++)
    {
        position[i] = center[i] - imageCenter[i];
    }
    std::ofstream file;
    file.open(pInfoFile.toStdString().data(), std::ios::out|std::ios::app);
    file<<"CenterToImage =";
    for(int i=0;i<3;i++)
    {
        file<<" "<<position[i];
    }
    file<<std::endl;
    file.close();
    step = 100;
    if(pCallback)
    {
        pCallback(step);
    }

    qDebug()<<"DicomImageUtil::ImageMask finish"<<elapsedTimer.elapsed()<<"ms";
    return true;
}

vtkSmartPointer<vtkImageData> DicomImageUtil::CutImage(vtkSmartPointer<vtkImageData> pImageData, int pExtent[6])
{
    qDebug()<<QTime::currentTime()<<"DicomImageUtil::CutImage"<<pExtent[0]<<pExtent[1]<<pExtent[2]<<pExtent[3]<<pExtent[4]<<pExtent[5];
    vtkSmartPointer<vtkImageData> newImageData = vtkSmartPointer<vtkImageData>::New();
    newImageData->SetSpacing(pImageData->GetSpacing());
    newImageData->SetExtent(0, pExtent[1] - pExtent[0],0, pExtent[3] - pExtent[2],0, pExtent[5] - pExtent[4]);
    newImageData->SetOrigin(0,0,0);
    auto scalarType = pImageData->GetScalarType();
    auto numComponents = pImageData->GetNumberOfScalarComponents();
    newImageData->AllocateScalars(scalarType, numComponents);
    int newX=0;
    for(int x=pExtent[0];x<=pExtent[1];x++,newX++)
    {
        int newY=0;
        for(int y=pExtent[2];y<=pExtent[3];y++,newY++)
        {
            int newZ=0;
            for(int z=pExtent[4];z<=pExtent[5];z++,newZ++)
            {
                short *pixelValue = (short *)pImageData->GetScalarPointer(x, y, z);
                short *newPixelValue = (short *)newImageData->GetScalarPointer(newX, newY, newZ);
                *newPixelValue = *pixelValue;
            }
        }
    }
    qDebug()<<QTime::currentTime()<<"DicomImageUtil::CutImage finish";
    return newImageData;
}

vtkSmartPointer<vtkImageData> DicomImageUtil::DilateImage(vtkSmartPointer<vtkImageData> pImageData, int pDilateValue, int pForegroundValue, int pBackgroundValue)
{
    qDebug()<<QTime::currentTime()<<"DicomImageUti::DilateImage"<<pDilateValue;
    //由于膨胀时图像边界不会扩张，因此，在调用膨胀处理时，先将图像边界进行扩张
    auto extent = pImageData->GetExtent();
    auto origin = pImageData->GetOrigin();
    auto spacing = pImageData->GetSpacing();
    auto newExtent = new int[6];
    auto newOrigin = new double[3];
    for(int i=0;i<3;i++)
    {
        newExtent[2*i] = extent[2*i];
        newExtent[2*i+1] = extent[2*i+1]+2*pDilateValue;
        newOrigin[i] = origin[i] - pDilateValue*spacing[i];
    }
    vtkSmartPointer<vtkImageData> newImageData = vtkSmartPointer<vtkImageData>::New();
    newImageData->SetSpacing(spacing);
    newImageData->SetExtent(newExtent);
    newImageData->SetOrigin(newOrigin);
    auto scalarType = pImageData->GetScalarType();
    auto numComponents = pImageData->GetNumberOfScalarComponents();
    newImageData->AllocateScalars(scalarType, numComponents);
    for(int x=extent[0] - pDilateValue;x<=extent[1] + pDilateValue;x++)
    {
        for(int y=extent[2] - pDilateValue;y<=extent[3] + pDilateValue;y++)
        {
            for(int z=extent[4] - pDilateValue;z<=extent[5] + pDilateValue;z++)
            {
                auto newX = x + pDilateValue, newY = y + pDilateValue, newZ = z + pDilateValue;
                short *newPixelValue = (short *)newImageData->GetScalarPointer(newX, newY, newZ);
                if(x<extent[0] || y<extent[2] || z<extent[4] || x>=extent[1] || y>=extent[3] || z>=extent[5])
                {
                    *newPixelValue = pBackgroundValue;
                }
                else
                {
                    short *pixelValue = (short *)pImageData->GetScalarPointer(x, y, z);
                    *newPixelValue = *pixelValue;
                }
            }
        }
    }
    //TODO 由于测试了vtk的膨胀类没有效果，暂时使用itk的膨胀类
    //转换为itkImage
    const int dimension = 3;
    using ImageType = itk::Image<signed short, dimension>;
    using VTKToITK = itk::VTKImageToImageFilter<ImageType>;
    VTKToITK::Pointer vtkToItk = VTKToITK::New();
    vtkToItk->SetInput(newImageData);
    vtkToItk->Update();
    auto itkImage = vtkToItk->GetOutput();
    //使用itk方法进行膨胀
    using StructuringElementType = itk::FlatStructuringElement<dimension>;
    StructuringElementType::RadiusType radius;
    radius.Fill(pDilateValue);
    StructuringElementType structuringElement = StructuringElementType::Ball(radius);
    using BinaryDilateImageFilterType = itk::BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType>;
    BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(itkImage);
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->SetForegroundValue(pForegroundValue);
    //转换为vtkImage
    using ITKToVTK = itk::ImageToVTKImageFilter<ImageType>;
    ITKToVTK::Pointer itkToVtk = ITKToVTK::New();
    itkToVtk->SetInput(dilateFilter->GetOutput());
    itkToVtk->Update();
    newImageData = vtkSmartPointer<vtkImageData>::New();
    newImageData->DeepCopy(itkToVtk->GetOutput());
    qDebug()<<QTime::currentTime()<<"DicomImageUtil::DilateImage finish";
    return newImageData;
}
