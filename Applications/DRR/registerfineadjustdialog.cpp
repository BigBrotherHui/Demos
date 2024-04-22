#pragma execution_character_set("utf-8")
#include "registerfineadjustdialog.h"
#include "ui_registerfineadjustdialog.h"
#include "generatematrixhelper.h"
#include "qmessagebox.h"
//itk
#include "itkImageToVTKImageFilter.h"
#include "itkOpenCVImageBridge.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "itkMatrix.h"
#include "itkRGBPixel.h"
#include "itkImportImageFilter.h"
#include "itkFlipImageFilter.h"
//vtk
#include "vtkImageBlend.h"
#include "QVTKOpenGLNativeWidget.h"
#include "vtkImageMapper3D.h"
#include "vtkImageProperty.h"
#include "vtkGenericRenderWindowInteractor.h"
#include "vtkPropPicker.h"

#include "vtkPointPicker.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkImageCast.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkInteractorStyleImage.h"
#include <vtkImageFlip.h>
#include <vtkCamera.h>
//Qt
#include <QDebug>
#include <QFileDialog>
#include <QJsonDocument>
#include <QJsonArray>
#include <QAction>
#include <QMenu>
#include <QPainter>
#include <QPolygon>
#include <QtConcurrent>
#include <QKeyEvent>
#include "opencv2/opencv.hpp"
#include <vtkImageWriter.h>
#include <vtkInformation.h>
#include "vtkMetaImageWriter.h"
#include <itkGDCMSeriesFileNames.h>
#include <vtkImageShiftScale.h>


RegisterFineAdjustDialog::RegisterFineAdjustDialog(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::RegisterFineAdjustDialog)
{
    ui->setupUi(this);
    setWindowFlag(Qt::Window);
    MatrixHelper=new GenerateMatrixHelper;
//    ui->comboBox_vetebra->setCheckExclusive(true);
    take_over();

    mActorMergedFront=vtkSmartPointer<vtkImageActor>::New();
    mActorMergedSide=vtkSmartPointer<vtkImageActor>::New();
    mRenderwindowMergedFront=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    mRenderwindowMergedSide=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    mRenderwindowProjectionFront=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    mRenderwindowProjectionSide=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();

    wgMergedFront=new QVTKOpenGLNativeWidget(this);
    wgMergedSide=new QVTKOpenGLNativeWidget(this);

    QHBoxLayout *left=new QHBoxLayout;
    left->addWidget(wgMergedFront);
    QHBoxLayout *right=new QHBoxLayout;
    right->addWidget(wgMergedSide);
    ui->widget_front->setLayout(left);
    ui->widget_side->setLayout(right);

    vtkNew<vtkPointPicker> pickerFront,pickerSide;
    rendererFront=vtkSmartPointer<vtkRenderer>::New();
    rendererSide=vtkSmartPointer<vtkRenderer>::New();
    rendererFront->AddActor(mActorMergedFront);
    rendererSide->AddActor(mActorMergedSide);

    interatorFront=vtkSmartPointer<vtkGenericRenderWindowInteractor>::New();
    interatorFront->SetPicker(vtkPointPicker::New());
    mRenderwindowMergedFront->SetInteractor(interatorFront);
    wgMergedFront->setRenderWindow(mRenderwindowMergedFront);

    interatorSide=vtkSmartPointer<vtkGenericRenderWindowInteractor>::New();
    interatorSide->SetPicker(vtkPointPicker::New());
    mRenderwindowMergedSide->SetInteractor(interatorSide);
    wgMergedSide->setRenderWindow(mRenderwindowMergedSide);

    //duplicatorSide = DuplicatorType::New();
    //duplicatorFront = DuplicatorType::New();

    //ui->widget_titleBar->setTitle("手动配准");
}

RegisterFineAdjustDialog::~RegisterFineAdjustDialog()
{
    delete ui;
    if(siddonGPUFront)
    {
        delete siddonGPUFront;
        siddonGPUFront=nullptr;
    }
    //if(siddonGPUSide)
    //{
    //    delete siddonGPUSide;
    //    siddonGPUSide=nullptr;
    //}
    if(MatrixHelper)
    {
        delete MatrixHelper;
        MatrixHelper=nullptr;
    }
    if(resultFront)
    {
        delete[] resultFront;
        resultFront=nullptr;
    }
    if(resultSide)
    {
        delete[] resultSide;
        resultSide=nullptr;
    }

}

void RegisterFineAdjustDialog::itkReadMhdFile(QString path)
{
    itk::ObjectFactoryBase::RegisterFactory(itk::MetaImageIOFactory::New());
    RawReaderType::Pointer reader = RawReaderType::New();
    reader->SetFileName(path.toUtf8().data());
    reader->Update();
    m_rawImage=reader->GetOutput();
    m_rawImage=imageinterpolation(m_rawImage,m_rawImage->GetSpacing()[0]);
    //若病案类型是模型，
    //int caseType=gCurrentCase.CaseType;
    CT3DNormalization(m_rawImage,100);
}

itk::Image<float,3>::Pointer RegisterFineAdjustDialog::imageinterpolation(RawImageType::Pointer imginput, double isoSpacing)
{
    using ImageType=RawImageType;
    const int Dimension=3;
    using ResampleFilterType4Interp =
        itk::ResampleImageFilter<ImageType, ImageType>;
    ResampleFilterType4Interp::Pointer resamplerInterp = ResampleFilterType4Interp::New();
    using TransformType4Interp = itk::IdentityTransform<double, Dimension>;
    TransformType4Interp::Pointer transform = TransformType4Interp::New();
    transform->SetIdentity();
    resamplerInterp->SetTransform(transform);
    using InterpolatorType4Interp =
        itk::LinearInterpolateImageFunction<ImageType, double>;

    InterpolatorType4Interp::Pointer interpolator = InterpolatorType4Interp::New();
    resamplerInterp->SetInterpolator(interpolator);
    resamplerInterp->SetDefaultPixelValue(0);

    const ImageType::SpacingType & inputSpacing = imginput->GetSpacing();
    ImageType::SpacingType spacing;
    spacing[0] = isoSpacing;
    spacing[1] = isoSpacing;
    spacing[2] = isoSpacing;
    resamplerInterp->SetOutputSpacing(spacing);
    resamplerInterp->SetOutputOrigin(imginput->GetOrigin());
    resamplerInterp->SetOutputDirection(imginput->GetDirection());

    ImageType::SizeType inputSize =
        imginput->GetLargestPossibleRegion().GetSize();

    using SizeValueType = ImageType::SizeType::SizeValueType;
    const double dx = inputSize[0] * inputSpacing[0] / isoSpacing;
    const double dy = inputSize[1] * inputSpacing[1] / isoSpacing;
    const double dz = std::floor((inputSize[2]) * inputSpacing[2] / isoSpacing);
    ImageType::SizeType size;

    size[0] = static_cast<SizeValueType>(dx);
    size[1] = static_cast<SizeValueType>(dy);
    size[2] = static_cast<SizeValueType>(dz);

    resamplerInterp->SetSize(size);
    resamplerInterp->SetInput(imginput);

    try
    {
        resamplerInterp->Update();
    }
    catch (itk::ExceptionObject & ex)
    {
        std::cout << ex << std::endl;
    }
    return resamplerInterp->GetOutput();;
}

void RegisterFineAdjustDialog::CT3DNormalization(RawImageType::Pointer img, double dThreshold)
{
    using ImageType=RawImageType;
    ImageType::SizeType inputSizess =
        img->GetLargestPossibleRegion().GetSize();
    float * buffer = img->GetBufferPointer();
    int len = inputSizess[0] * inputSizess[1] * inputSizess[2];
    for (int i = 0; i < len; ++i)
    {
        if (buffer[i] > dThreshold)
        {
            buffer[i] = buffer[i] / 1000 - dThreshold / 1000;
        }
        else
        {
            buffer[i] = 0;
        }

    }
}

void RegisterFineAdjustDialog::itkReadXRayImageFront(QString path)
{
    XRayReaderType::Pointer reader = XRayReaderType::New();
    ImageIOType::Pointer gdcmIO = ImageIOType::New();
    reader->SetImageIO(gdcmIO);
    reader->SetFileName(path.toStdString());
    reader->Update();
    typedef itk::FlipImageFilter<XRayImageType> FlipfilterType;
    typedef FlipfilterType::FlipAxesArrayType FlipAxesArrayType;
    FlipfilterType::Pointer filter=FlipfilterType::New();
    FlipAxesArrayType flipArray;
    flipArray[0]=0;
    flipArray[1]=1;
    filter->SetFlipAxes(flipArray);
    filter->SetInput(reader->GetOutput());
    filter->Update();
    m_frontImage=filter->GetOutput();
}

void RegisterFineAdjustDialog::itkReadXRayImageSide(QString path)
{
    XRayReaderType::Pointer reader = XRayReaderType::New();
    ImageIOType::Pointer gdcmIO = ImageIOType::New();
    reader->SetImageIO(gdcmIO);
    reader->SetFileName(path.toStdString());
    reader->Update();
    typedef itk::FlipImageFilter<XRayImageType> FlipfilterType;
    typedef FlipfilterType::FlipAxesArrayType FlipAxesArrayType;
    FlipfilterType::Pointer filter=FlipfilterType::New();
    FlipAxesArrayType flipArray;
    flipArray[0]=0;
    flipArray[1]=1;
    filter->SetFlipAxes(flipArray);
    filter->SetInput(reader->GetOutput());
    filter->Update();
    m_sideImage=filter->GetOutput();
}

void RegisterFineAdjustDialog::ConvertedInverse(itk::Image<float,2>::Pointer img)
{
    using ImageType=itk::Image<float,2>;
    ImageType::SizeType inputSizess2 =
        img->GetLargestPossibleRegion().GetSize();
    float * buffer3 = img->GetBufferPointer();
    int k, f;
    float temp;
    int len2 = inputSizess2[0] * inputSizess2[1] * inputSizess2[2];

    for (int i = 0; i < inputSizess2[0]; i++) {
        for (int j = 0; j < inputSizess2[0] / 2; j++) {//注意此处循环限制语句j<i
            k = j * inputSizess2[0] + i;
            f = (inputSizess2[0] - j - 1) * inputSizess2[0] + i;
            temp = buffer3[k];
            buffer3[k] = buffer3[f];
            buffer3[f] = temp;
        }
    }
    for (int i = 0; i < inputSizess2[0]; i++) {
        for (int j = 0; j < i; j++) {//注意此处循环限制语句j<i
            k = i * inputSizess2[0] + j;
            f = j * inputSizess2[0] + i;
            temp = buffer3[k];
            buffer3[k] = buffer3[f];
            buffer3[f] = temp;
        }
    }
}

vtkSmartPointer<vtkImageData> RegisterFineAdjustDialog::mergeTwoImages(vtkImageData *origin, vtkImageData *processed,double v)
{
    vtkNew<vtkImageData> imgData;
    vtkSmartPointer<vtkImageBlend> blend = vtkSmartPointer<vtkImageBlend>::New();
    blend->AddInputData(origin);
    blend->AddInputData(processed);
    blend->SetOpacity(0, v);
    blend->SetOpacity(1, 1-v);
    blend->Update();
    imgData->DeepCopy(blend->GetOutput());
    return imgData;
}

itk::Image<float,2>::Pointer RegisterFineAdjustDialog::mergeTwoImages(XRayImageType::Pointer pimg1, XRayImageType::Pointer pimg2,double v)
{
    vtkNew<vtkImageData> vimg1;
    vtkNew<vtkImageData> vimg2;
    vtkNew<vtkImageData> vimgres;

    itk::ImageToVTKImageFilter<itk::Image<float,2>>::Pointer filter=itk::ImageToVTKImageFilter<itk::Image<float,2>>::New();
    filter->SetInput(pimg1);
    filter->Update();
    vimg1->DeepCopy(filter->GetOutput());
    filter->SetInput(pimg2);
    filter->Update();
    vimg2->DeepCopy(filter->GetOutput());
    vimgres->DeepCopy(mergeTwoImages(vimg1,vimg2,v));
    itk::VTKImageToImageFilter<itk::Image<float,2>>::Pointer filter1=itk::VTKImageToImageFilter<itk::Image<float,2>>::New();
    filter1->SetInput(vimgres);
    filter1->Update();
    return filter1->GetOutput();
}

void RegisterFineAdjustDialog::readMetaHeaderFile(QString path)
{
    if(path.isEmpty())
        return;
    //读取CenterToImage字段的值
    QFile file(path);
    if(file.open(QFile::ReadOnly))
    {
        while(!file.atEnd())
        {
            QString line=file.readLine();
            QStringList foobar=line.simplified().split(" ");
            if(foobar.indexOf("CenterToImage")!=-1 && foobar.size()>=5)
            {
                m_CenterToImage[0]=foobar[2].toFloat();
                m_CenterToImage[1]=foobar[3].toFloat();
                m_CenterToImage[2]=foobar[4].toFloat();
            }
        }
    }
}

void RegisterFineAdjustDialog::readMatrixInternal(QString path)
{
    QDir dir(path);
    QStringList filepaths=dir.entryList();
    if(filepaths.indexOf("s1tos2.txt")==-1 || filepaths.indexOf("s1torobot.txt")==-1)
        return;
    QString s1tos2FilePath=path+"/s1tos2.txt",s1torobotFilePath=path+"/s1torobot.txt";
    auto getMatrix=[&](QString path,float *out)
    {
        QFile file(path);
        bool isOpen=file.open(QIODevice::ReadOnly);
        if(!isOpen)
            return;
        QString str_matrix=file.readAll();
        str_matrix=str_matrix.simplified();
        QStringList list_element=str_matrix.split(" ");
        int elementSize=list_element.size();
        if(elementSize<16)
            return;
        for(int i=0;i<16;i++)
        {
            out[i]=list_element[i].toFloat();
        }
    };
    getMatrix(s1tos2FilePath,s1tos2);
    getMatrix(s1torobotFilePath,s1torobot);
}

void RegisterFineAdjustDialog::inversePixel(itk::Image<float,2>::Pointer pImage)
{
    auto imageSize=pImage->GetLargestPossibleRegion().GetSize();
    int size=imageSize[0]*imageSize[1];
    float *buffer=pImage->GetBufferPointer();
    for(int i=0;i<size;i++)
    {
        buffer[i]=4095-buffer[i];
    }
}

void RegisterFineAdjustDialog::showImagesVtk(itk::Image<float,2>::Pointer imageFrontProjection, itk::Image<float,2>::Pointer imageSideProjection,double pXRayOpacity,bool pIsRgb)
{
    auto castFunc1=[&](RGBImageType::Pointer img)
    {
        vtkNew<vtkImageData> imgData;
        itk::ImageToVTKImageFilter<RGBImageType>::Pointer filter=itk::ImageToVTKImageFilter<RGBImageType>::New();
        filter->SetInput(img);
        filter->Update();
        imgData->DeepCopy(filter->GetOutput());
        return imgData;
    };

    auto castFunc2=[&](GrayImageType::Pointer img)
    {
        vtkNew<vtkImageData> imgData;
        itk::ImageToVTKImageFilter<GrayImageType>::Pointer filter=itk::ImageToVTKImageFilter<GrayImageType>::New();
        filter->SetInput(img);
        filter->Update();
        imgData->DeepCopy(filter->GetOutput());
        return imgData;
    };


    if(pIsRgb)
    {
        RGBImageType::Pointer rgbImageFront=convertToRgb(imageFrontProjection,m_frontImage,pXRayOpacity);
        vtkSmartPointer<vtkImageData> vRgbImageFront=castFunc1(rgbImageFront);
        //RGBImageType::Pointer rgbImageSide=convertToRgb(imageSideProjection,m_sideImage,pXRayOpacity);
        //vtkSmartPointer<vtkImageData> vRgbImageSide=castFunc1(rgbImageSide);
        showImages(vRgbImageFront,nullptr);
    }
    else
    {
        vtkSmartPointer<vtkImageData> originFront=castFunc2(m_frontImage);
        vtkSmartPointer<vtkImageData> projectionFront=castFunc2(imageFrontProjection);
        vtkSmartPointer<vtkImageData> mergedFront=mergeTwoImages(originFront,projectionFront,pXRayOpacity);
        /*vtkSmartPointer<vtkImageData> originSide=castFunc2(m_sideImage);
        vtkSmartPointer<vtkImageData> projectionSide=castFunc2(imageSideProjection);
        vtkSmartPointer<vtkImageData> mergedSide=mergeTwoImages(originSide,projectionSide,pXRayOpacity);*/
        showImages(mergedFront,nullptr);
    }
}

void RegisterFineAdjustDialog::showImages(vtkImageData *pMergedFront,vtkImageData *pMergedSide)
{
    if(!mIsWindowInitialized)
    {
        mIsWindowInitialized=true;
        interatorFront->Initialize();
        mRenderwindowMergedFront->AddRenderer(rendererFront);

        //interatorSide->Initialize();
        //mRenderwindowMergedSide->AddRenderer(rendererSide);
    }

    mActorMergedFront->GetMapper()->SetInputData(pMergedFront);
    //mActorMergedSide->GetMapper()->SetInputData(pMergedSide);


	wgMergedFront->renderWindow()->Render();
    //wgMergedSide->renderWindow()->Render();
}

void RegisterFineAdjustDialog::writerToFile(float *source, QString path,int size)
{
    size=imageSizeFront[0]*imageSizeFront[1];
    QFile file(path);
    if(file.open(QIODevice::WriteOnly))
    {
        for(int i=0;i<size;i++)
        {
            file.write(QString::number(source[i],'f',6).toUtf8()+" ");
        }
    }
    file.close();
}

void RegisterFineAdjustDialog::WriterImageHomogenization(itk::Image<float,2>::Pointer &image)
{
    using ImageType=itk::Image<float,2>;
    using OutputImageType = itk::Image<float, 2>;
    using RescaleFilterType = itk::RescaleIntensityImageFilter<ImageType, OutputImageType>;
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(image);
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(4095);
    rescaleFilter->Update();
    image=rescaleFilter->GetOutput();
}

void RegisterFineAdjustDialog::PatientourceToRobotMatrix(double *par, float *ds1torobot,  float *result)
{
    using TransformType = itk::Euler3DTransform<double>;
    TransformType::Pointer transform = TransformType::New();
    transform->SetComputeZYX(true);
    double dtr = (atan(1.0) * 4.0) / 180.0;
    transform->SetRotation(dtr * par[0], dtr * par[1], dtr * par[2]);
    TransformType::MatrixType R = transform->GetMatrix();
    vtkNew<vtkMatrix4x4> MatrixEuler;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
            MatrixEuler->SetElement(i,j,R[i][j]);
    }
    MatrixEuler->SetElement(0,3,par[3]);
    MatrixEuler->SetElement(1,3,par[4]);
    MatrixEuler->SetElement(2,3,par[5]);

    vtkNew<vtkMatrix4x4> matrixRs1ToRobot;
    for(int i=0;i<16;i++)
        matrixRs1ToRobot->SetElement(i/4,i%4,ds1torobot[i]);
    vtkNew<vtkMatrix4x4> MatrixDst;
    vtkMatrix4x4::Multiply4x4(matrixRs1ToRobot,MatrixEuler,MatrixDst);
    for(int i=0;i<16;i++)
        result[i]=MatrixDst->GetElement(i/4,i%4);
}

itk::Image<itk::RGBPixel<float>>::Pointer RegisterFineAdjustDialog::convertToRgb(itk::Image<float,2>::Pointer pImage1, itk::Image<float,2>::Pointer pImage2,double pXRayOpacity)
{
    using PixelType=itk::RGBPixel<float>;
    using RGBImageType=itk::Image<PixelType,2>;
    typedef itk::ImportImageFilter<PixelType, 2 >   ImportFilterType;
    RGBImageType::Pointer rgbImage=RGBImageType::New();
    ImportFilterType::Pointer importFilter = ImportFilterType::New();
    ImportFilterType::SizeType  size;
    auto sourceImageSize=pImage1->GetLargestPossibleRegion().GetSize();
    size[0]  = sourceImageSize[0];
    size[1]  = sourceImageSize[1];
    ImportFilterType::IndexType start;
    start.Fill(0);
    ImportFilterType::RegionType region;
    region.SetIndex(start);
    region.SetSize( size);
    importFilter->SetRegion(region);
    importFilter->SetOrigin(pImage1->GetOrigin());
    importFilter->SetSpacing(pImage1->GetSpacing());
    int numberOfPixels=size[0]*size[1];
    PixelType * localBuffer = new PixelType[numberOfPixels];
    bool importImageFilterWillOwnTheBuffer = true;
    importFilter->SetImportPointer( localBuffer, numberOfPixels,
                                   importImageFilterWillOwnTheBuffer );
    importFilter->Update();
    rgbImage= importFilter->GetOutput();
    float *buf1=pImage1->GetBufferPointer();
    float *buf2=pImage2->GetBufferPointer();
    for(int i=0;i<size[0];i++)
    {
        for(int j=0;j<size[1];j++)
        {
            RGBImageType::IndexType index{j,i};
            PixelType &pixel=rgbImage->GetPixel(index);
            pixel.SetRed((1-pXRayOpacity)*buf1[i*size[0]+j]);
            pixel.SetGreen(pXRayOpacity*buf2[i*size[0]+j]);
        }
    }
    return rgbImage;
}

void RegisterFineAdjustDialog::readMatrix(QString path, float *matrix)
{
    QFile file(path);
    if(!file.open(QIODevice::ReadOnly))
        return;
    QList<QByteArray> elements=file.readAll().simplified().split(' ');
    for(int i=0;i<elements.size();i++)
    {
        matrix[i]=elements.at(i).toFloat();
    }
}

void RegisterFineAdjustDialog::initAdjustInfo()
{
    qDebug()<<"RegisterFineAdjustDialog::initAdjustInfo start";
    auto findObjectAdjustInfoFunc=[&](QString awlName) -> AdjustInfo &
    {
        for(auto &pair : mAdjustInfo)
        {
            auto &info=pair.second;
            if(info.awlName==awlName)
                return info;
        }
    };
    if(mAwlList.indexOf("All")==-1)
        mAwlList.push_back("All");
    for(auto awlName : mAwlList)
    {
        QString name=awlName+"Result";
        QString prefix=mJsonObject.value("path").toString();
        QString resultPath;
        if(awlName=="All")
            resultPath="Result.txt";
        else
           resultPath=mJsonObject.value(name).toString();
        resultPath.prepend(prefix+"/");

        float result[16];
        readMatrix(resultPath,result);

        vtkNew<vtkMatrix4x4> matrixS1toRobot;
        for(int i=0;i<16;i++)
            matrixS1toRobot->SetElement(i/4,i%4,s1torobot[i]);
        vtkSmartPointer<vtkMatrix4x4> matrixRobotToS1=matrixInverse(matrixS1toRobot);
        vtkNew<vtkMatrix4x4> MatrixInitial,MatrixResult;
        for(int i=0;i<16;i++)
            MatrixResult->SetElement(i/4,i%4,result[i]);
        vtkMatrix4x4::Multiply4x4(matrixRobotToS1,MatrixResult,MatrixInitial);
        double sy = std::sqrt(MatrixInitial->GetElement(0, 0) * MatrixInitial->GetElement(0, 0) +MatrixInitial->GetElement(1, 0) * MatrixInitial->GetElement(1, 0));
        const double dtr = (atan(1.0) * 4.0) / 180.0;
        double rotateX, rotateY, rotateZ,translateX,translateY,translateZ;
        if (sy > 1e-6)
        {
            rotateX = atan2(MatrixInitial->GetElement(2, 1), MatrixInitial->GetElement(2, 2))/ dtr;
            rotateY = atan2(-MatrixInitial->GetElement(2, 0), sy)/ dtr;
            rotateZ = atan2(MatrixInitial->GetElement(1, 0),MatrixInitial->GetElement(0, 0))/ dtr;
        }
        else
        {
            rotateX = atan2(-MatrixInitial->GetElement(1, 2),MatrixInitial->GetElement(1, 1))/ dtr;
            rotateY = atan2(-MatrixInitial->GetElement(2, 0), sy)/ dtr;
            rotateZ = 0/ dtr;
        }
        translateX=MatrixInitial->GetElement(0,3);
        translateY=MatrixInitial->GetElement(1,3);
        translateZ=MatrixInitial->GetElement(2,3);

        AdjustInfo &info=findObjectAdjustInfoFunc(awlName);
        info.rotateX=rotateX;
        info.rotateY=rotateY;
        info.rotateZ=rotateZ;
        info.translateX=translateX;
        info.translateY=translateY;
        info.translateZ=translateZ;
    }
    qDebug()<<"RegisterFineAdjustDialog::initAdjustInfo end";
}

vtkSmartPointer<vtkMatrix4x4> RegisterFineAdjustDialog::matrixInverse(vtkMatrix4x4 *r)
{
    vtkNew<vtkMatrix4x4> RM;
    vtkMatrix4x4::Transpose(r,RM);
    for(int i=0;i<3;i++)
    {
        RM->SetElement(3,i,0);
    }
    vtkNew<vtkMatrix4x4> TM;
    TM->SetElement(0,3,-r->GetElement(0,3));
    TM->SetElement(1,3,-r->GetElement(1,3));
    TM->SetElement(2,3,-r->GetElement(2,3));
    vtkNew<vtkMatrix4x4> result;
    vtkMatrix4x4::Multiply4x4(RM,TM,result);
    return result;
}

void RegisterFineAdjustDialog::showEvent(QShowEvent *event)
{
    mLastPath="";
    mIsLastValueChanged=false;
    if(mItems.size()>0)
    {
//        ui->comboBox_vetebra->setCurrentIndex(0);
//        ui->comboBox_vetebra->setItemCheckState(0,Qt::Checked);
        mItems.first()->click();
//        checkChanged();
    }
    return QWidget::showEvent(event);
}

//void RegisterFineAdjustDialog::paintEvent(QPaintEvent *e)
//{
//    QPainter painter(this);
//    QPen pen(QColor("#2F3443"),2);
//    painter.setPen(pen);
//    painter.drawRoundRect(rect(),2,2);
//    return QWidget::paintEvent(e);
//}

void RegisterFineAdjustDialog::show()
{
    QEventLoop loop;
    QWidget::show();
    loop.exec();
}

void RegisterFineAdjustDialog::addMenu(QVTKOpenGLNativeWidget *pTarget)
{
    //新建了5个椎节、2个方位
    QAction* L1 = new QAction(tr("L1"), pTarget);
    QAction* L2 = new QAction(tr("L2"), pTarget);
    QAction* L3 = new QAction(tr("L3"), pTarget);
    QAction* L4 = new QAction(tr("L4"), pTarget);
    QAction* L5 = new QAction(tr("L5"),pTarget);
    QAction* Left1 = new QAction(tr("Left"), L1);
    QAction* Right1 = new QAction(tr("Right"),L1);
    QAction* Left2 = new QAction(tr("Left"), L2);
    QAction* Right2 = new QAction(tr("Right"),L2);
    QAction* Left3 = new QAction(tr("Left"), L3);
    QAction* Right3 = new QAction(tr("Right"),L3);
    QAction* Left4 = new QAction(tr("Left"), L4);
    QAction* Right4 = new QAction(tr("Right"),L4);
    QAction* Left5 = new QAction(tr("Left"), L5);
    QAction* Right5 = new QAction(tr("Right"),L5);

    QMenu *menuL1=new QMenu(pTarget);
    menuL1->addAction(Left1);
    menuL1->addAction(Right1);

    QMenu *menuL2=new QMenu(pTarget);
    menuL2->addAction(Left2);
    menuL2->addAction(Right2);

    QMenu *menuL3=new QMenu(pTarget);
    menuL3->addAction(Left3);
    menuL3->addAction(Right3);

    QMenu *menuL4=new QMenu(pTarget);
    menuL4->addAction(Left4);
    menuL4->addAction(Right4);

    QMenu *menuL5=new QMenu(pTarget);
    menuL5->addAction(Left5);
    menuL5->addAction(Right5);

    L1->setMenu(menuL1);
    L2->setMenu(menuL2);
    L3->setMenu(menuL3);
    L4->setMenu(menuL4);
    L5->setMenu(menuL5);

    pTarget->addAction(L1);
    pTarget->addAction(L2);
    pTarget->addAction(L3);
    pTarget->addAction(L4);
    pTarget->addAction(L5);
}

void RegisterFineAdjustDialog::take_over()
{
    QList<QDoubleSpinBox *> spinboxs;
    spinboxs<<ui->spinbox_rotateX<<ui->spinbox_rotateY<<ui->spinbox_rotateZ
             <<ui->spinbox_translateX<<ui->spinbox_translateY<<ui->spinbox_translateZ;
    for(auto &spinbox : spinboxs)
    {
        connect(spinbox,static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),[&]
                {
                    mIsLastValueChanged=true;
                    adjust();
                });
    }
//    connect(ui->comboBox_vetebra,&CheckComboBox::signal_checkStateChanged,this,&RegisterFineAdjustDialog::checkChanged);
    connect(ui->horizontalSlider_Opacity,&QSlider::valueChanged,this,&RegisterFineAdjustDialog::adjustByOpacity);
    
    connect(ui->horizontalSlider_transX,&QSlider::valueChanged,ui->spinbox_translateX,&QDoubleSpinBox::setValue);
    connect(ui->horizontalSlider_transY,&QSlider::valueChanged,ui->spinbox_translateY,&QDoubleSpinBox::setValue);
    connect(ui->horizontalSlider_transZ,&QSlider::valueChanged,ui->spinbox_translateZ,&QDoubleSpinBox::setValue);
    connect(ui->horizontalSlider_rotateX,&QSlider::valueChanged,ui->spinbox_rotateX,&QDoubleSpinBox::setValue);
    connect(ui->horizontalSlider_rotateY,&QSlider::valueChanged,ui->spinbox_rotateY,&QDoubleSpinBox::setValue);
    connect(ui->horizontalSlider_rotateZ,&QSlider::valueChanged,ui->spinbox_rotateZ,&QDoubleSpinBox::setValue);

    connect(ui->spinbox_rotateX,static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),[&](double pValue)
            {
                ui->horizontalSlider_rotateX->blockSignals(true);
                ui->horizontalSlider_rotateX->setValue(pValue);
                ui->horizontalSlider_rotateX->blockSignals(false);
            });
    connect(ui->spinbox_rotateY,static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),[&](double pValue)
            {
                ui->horizontalSlider_rotateY->blockSignals(true);
                ui->horizontalSlider_rotateY->setValue(pValue);
                ui->horizontalSlider_rotateY->blockSignals(false);
            });

    connect(ui->spinbox_rotateZ,static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),[&](double pValue)
            {
                ui->horizontalSlider_rotateZ->blockSignals(true);
                ui->horizontalSlider_rotateZ->setValue(pValue);
                ui->horizontalSlider_rotateZ->blockSignals(false);
            });

    connect(ui->spinbox_translateX,static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),[&](double pValue)
            {
                ui->horizontalSlider_transX->blockSignals(true);
                ui->horizontalSlider_transX->setValue(pValue);
                ui->horizontalSlider_transX->blockSignals(false);
            });

    connect(ui->spinbox_translateY,static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),[&](double pValue)
            {
                ui->horizontalSlider_transY->blockSignals(true);
                ui->horizontalSlider_transY->setValue(pValue);
                ui->horizontalSlider_transY->blockSignals(false);
            });

    connect(ui->spinbox_translateZ,static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),[&](double pValue)
            {
                ui->horizontalSlider_transZ->blockSignals(true);
                ui->horizontalSlider_transZ->setValue(pValue);
                ui->horizontalSlider_transZ->blockSignals(false);
            });

    
    connect(ui->pushButton_importImage,&QPushButton::clicked,this,&RegisterFineAdjustDialog::slot_importImage_clicked);
}

void RegisterFineAdjustDialog::adjust()
{
    //获取当前6个自由度的值
    double transX,transY,transZ,rotateX,rotateY,rotateZ;
    transX=ui->spinbox_translateX->value();
    transY=ui->spinbox_translateY->value();
    transZ=ui->spinbox_translateZ->value();
    rotateX=ui->spinbox_rotateX->value();
    rotateY=ui->spinbox_rotateY->value();
    rotateZ=ui->spinbox_rotateZ->value();
    MatrixHelper->generateEulerTransform(rotateX,rotateY,rotateZ,transX,transY,transZ);
    float transformMatrixFront[12];
    MatrixHelper->getTransformMatrix(true,transformMatrixFront);
    int ImageSize=imageSizeFront[0]*imageSizeFront[1];


    if(!resultFront)
        resultFront=new float[ImageSize]{0.0};

    float tmp[12]{ 0 };
    siddonGPUFront->SetTransformMatrix(tmp);


    
    siddonGPUFront->Run(transformMatrixFront,resultFront);
    //for(int i=0;i< ImageSize;++i)
    //{
    //    std::cout << i << " " << resultFront[i] << std::endl;
    //}
    vtkSmartPointer<vtkImageData> ret = vtkSmartPointer<vtkImageData>::New();
    vtkSmartPointer<vtkInformation> in = vtkSmartPointer<vtkInformation>::New();
    int mDimensions[3]{ imageSizeFront[0],imageSizeFront[1],1};
    int mExtent[6]{ 0,0,0,0,0,0 };
    mExtent[1] = mDimensions[0];
    mExtent[3] = mDimensions[1];
    mExtent[5] = mDimensions[2];
    ret->SetScalarType(VTK_FLOAT, in);
    ret->SetSpacing(.5, .5, 1);
    ret->SetOrigin(0,0,0);
    ret->SetExtent(mExtent);
    ret->SetDimensions(mDimensions);
    ret->SetNumberOfScalarComponents(1, in);
    ret->AllocateScalars(VTK_FLOAT, 1);
    memcpy(ret->GetScalarPointer(),resultFront,ImageSize*sizeof (float));
    ret->Modified();
    float* buffer = (float *)ret->GetScalarPointer();
    float min = ret->GetScalarRange()[0];
    float max = ret->GetScalarRange()[1];
    for (int i = 0; i < ImageSize; i++)
    {
        buffer[i] = 4096*(buffer[i]-min)/max-min;
    }
    ret->Modified();
    double originOpacity=ui->horizontalSlider_Opacity->value()/100.;
    mActorMergedFront->GetProperty()->SetColorLevel(2047);
    mActorMergedFront->GetProperty()->SetColorWindow(4095);
    showImages(ret,nullptr);
}

void RegisterFineAdjustDialog::adjustByOpacity()
{
    itk::Image<PixelType, 2>::Pointer processedImageFront{nullptr},processedImageSide{nullptr};
    for(int i=0;i<mItems.size();i++)
    {
        if(mItems[i]->isChecked())
        {
            std::array<itk::Image<PixelType, 2>::Pointer,2> imgs=mImageMap[mItems[i]];
            if(imgs.size()==0)
                return;
            processedImageFront=imgs[0];
            processedImageSide=imgs[1];
            break;
        }
    }
    double originOpacity=ui->horizontalSlider_Opacity->value()/100.;
    
}

void RegisterFineAdjustDialog::clearPoints()
{
    int levelFront=mActorMergedFront->GetProperty()->GetColorLevel();
    int windowFront=mActorMergedFront->GetProperty()->GetColorWindow();
    int levelSide=mActorMergedSide->GetProperty()->GetColorLevel();
    int windowSide=mActorMergedSide->GetProperty()->GetColorWindow();

    adjust();
    mActorMergedFront->GetProperty()->SetColorLevel(levelFront);
    mActorMergedFront->GetProperty()->SetColorWindow(windowFront);
    mActorMergedSide->GetProperty()->SetColorLevel(levelSide);
    mActorMergedSide->GetProperty()->SetColorWindow(windowSide);
    mRenderwindowMergedFront->Render();
    mRenderwindowMergedSide->Render();
}

void RegisterFineAdjustDialog::checkChanged()
{
    QStringList items=/*ui->comboBox_vetebra->*/checkedItems();
    if(items.size()==0)
        return;
    QString checkedItem=items[0];
    QString path=pair.at(checkedItem);
    if(!mLastPath.isEmpty() && mIsLastValueChanged)
    {
        AdjustInfo &info=mAdjustInfo[mLastPath];
        double transX,transY,transZ,rotateX,rotateY,rotateZ;
        transX=ui->spinbox_translateX->value();
        transY=ui->spinbox_translateY->value();
        transZ=ui->spinbox_translateZ->value();
        rotateX=ui->spinbox_rotateX->value();
        rotateY=ui->spinbox_rotateY->value();
        rotateZ=ui->spinbox_rotateZ->value();
        info.rotateX=rotateX;
        info.rotateY=rotateY;
        info.rotateZ=rotateZ;
        info.translateX=transX;
        info.translateY=transY;
        info.translateZ=transZ;
    }

    mLastPath=path;

    mIsLastValueChanged=false;
    //读取CenterToImage字段的值
    readMetaHeaderFile(path);
    MatrixHelper->setCenterOffset(m_CenterToImage);
    //读取椎节
    itkReadMhdFile(path);

    itk::Image<float,3>::SizeType size=m_rawImage->GetLargestPossibleRegion().GetSize();
    itk::Image<float,3>::SpacingType spacing=m_rawImage->GetSpacing();
    for(int i=0;i<3;i++)
    {
        _3dpixelNum[i]=size[i];
        _3dpixelSpacing[i]=spacing[0];
    }
    MatrixHelper->setPixelSpacing(spacingFront);

    ui->spinbox_translateX->blockSignals(true);
    ui->spinbox_translateY->blockSignals(true);
    ui->spinbox_translateZ->blockSignals(true);
    ui->spinbox_rotateX->blockSignals(true);
    ui->spinbox_rotateY->blockSignals(true);
    ui->spinbox_rotateZ->blockSignals(true);
    AdjustInfo info=mAdjustInfo[mLastPath];
    ui->horizontalSlider_transX->setValue(info.translateX);
    ui->horizontalSlider_transY->setValue(info.translateY);
    ui->horizontalSlider_transZ->setValue(info.translateZ);
    ui->horizontalSlider_rotateX->setValue(info.rotateX);
    ui->horizontalSlider_rotateY->setValue(info.rotateY);
    ui->horizontalSlider_rotateZ->setValue(info.rotateZ);

    ui->spinbox_translateX->setValue(info.translateX);
    ui->spinbox_translateY->setValue(info.translateY);
    ui->spinbox_translateZ->setValue(info.translateZ);
    ui->spinbox_rotateX->setValue(info.rotateX);
    ui->spinbox_rotateY->setValue(info.rotateY);
    ui->spinbox_rotateZ->setValue(info.rotateZ);

    ui->spinbox_translateX->blockSignals(false);
    ui->spinbox_translateY->blockSignals(false);
    ui->spinbox_translateZ->blockSignals(false);
    ui->spinbox_rotateX->blockSignals(false);
    ui->spinbox_rotateY->blockSignals(false);
    ui->spinbox_rotateZ->blockSignals(false);

    if(siddonGPUFront)
    {
        delete siddonGPUFront;
        siddonGPUFront=new SiddonGPU;
    }
    /*if(siddonGPUSide)
    {
        delete siddonGPUSide;
        siddonGPUSide=new SiddonGPU;
    }*/
    if(!siddonGPUFront)
        siddonGPUFront=new SiddonGPU;
    //if(!siddonGPUSide)
    //    siddonGPUSide=new SiddonGPU;
    float tmp[12]{0};
    siddonGPUFront->SetTransformMatrix(tmp);
    //siddonGPUSide->SetTransformMatrix(tmp);
    float *img=m_rawImage->GetBufferPointer();
    siddonGPUFront->SetImg3d(img,_3dpixelSpacing,_3dpixelNum);
    //siddonGPUSide->SetImg3d(img,_3dpixelSpacing,_3dpixelNum);
    siddonGPUFront->SetImg2dParameter(spacingFront,imageSizeFront);
    //siddonGPUSide->SetImg2dParameter(spacingSide,imageSizeSide,mask);
    adjust();
    if(!mIsBoneAwlSelected)
    {
        mIsBoneAwlSelected=true;
//        addMenu(wgMergedFront);
//        addMenu(wgMergedSide);
    }
}

void RegisterFineAdjustDialog::showRgb()
{
    adjust();
    rendererFront->ResetCamera();
    rendererSide->ResetCamera();
    rendererFront->Render();
    rendererSide->Render();
}

void RegisterFineAdjustDialog::SaveDrrToImg(const itk::Image<float, 2>::Pointer pImg, QString path)
{
    using OutputImageType = itk::Image<unsigned short, 2>;
    using RescaleFilterType = itk::RescaleIntensityImageFilter<itk::Image<float, 2>, OutputImageType>;
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(pImg);
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(4095);
    rescaleFilter->Update();
    OutputImageType::Pointer img = rescaleFilter->GetOutput();
    itk::Image<float, 2>::SizeType inputSizess2 =
        img->GetLargestPossibleRegion().GetSize();
    unsigned short * buffer = img->GetBufferPointer();
    int len = inputSizess2[0] * inputSizess2[1];
    FILE *fpwrt = NULL;
    const char* file_c = path.toStdString().c_str();
    fopen_s(&fpwrt, file_c, "wb+");
    if (fpwrt == NULL)
    {
        qDebug() << "写入文件创建错误";
    }
    fwrite(buffer, sizeof(unsigned short), len, fpwrt);
    fclose(fpwrt);
}

void RegisterFineAdjustDialog::slot_importImage_clicked()
{
    QString path = QFileDialog::getOpenFileName(this, "选择图像序列", "D:/Images/");
    if (path.isEmpty())
        return;
    /*using SeriesFileNamesType = itk::GDCMSeriesFileNames;
    using ReaderType = itk::ImageSeriesReader<RawImageType>;
    using ImageIOType = itk::GDCMImageIO;
    ImageIOType::Pointer m_gdcmImageIO = ImageIOType::New();
    SeriesFileNamesType::Pointer m_gdcmSeriesFileNames = SeriesFileNamesType::New();
    m_gdcmSeriesFileNames->SetUseSeriesDetails(true);
    m_gdcmSeriesFileNames->SetDirectory(dir.toStdString());
    using SeriesIdContainer = std::vector<std::string>;
    const SeriesIdContainer& seriesUID = m_gdcmSeriesFileNames->GetSeriesUIDs();
    const std::string& seriesUIDStr = seriesUID.at(0);
    using FileNamesContainer = std::vector<std::string>;
    FileNamesContainer fileNames = m_gdcmSeriesFileNames->GetFileNames(seriesUIDStr);
    ReaderType::Pointer m_reader= ReaderType::New();
    m_reader->SetImageIO(m_gdcmImageIO);
    m_reader->SetFileNames(fileNames);
    m_reader->ForceOrthogonalDirectionOff();*/
    itk::ObjectFactoryBase::RegisterFactory(itk::MetaImageIOFactory::New());
    typedef float PixelType;
    typedef itk::Image<PixelType, 3> RawImageType;
    typedef itk::ImageFileReader<RawImageType> RawReaderType;
    RawReaderType::Pointer reader = RawReaderType::New();
    reader->SetFileName(path.toUtf8().data());
    try {
        reader->Update();
        m_rawImage = reader->GetOutput();
        readMetaHeaderFile(path);
        MatrixHelper->setCenterOffset(m_CenterToImage);
        m_rawImage = imageinterpolation(m_rawImage, m_rawImage->GetSpacing()[0]);
        CT3DNormalization(m_rawImage, 100);
        itk::Image<float, 3>::SizeType size = m_rawImage->GetLargestPossibleRegion().GetSize();
        itk::Image<float, 3>::SpacingType spacing = m_rawImage->GetSpacing();
        for (int i = 0; i < 3; i++)
        {
            _3dpixelNum[i] = size[i];
            _3dpixelSpacing[i] = spacing[0];
        }
        MatrixHelper->setPixelSpacing(spacingFront);
        if (!siddonGPUFront)
            siddonGPUFront = new SiddonGPU;
        float tmp[12]{ 0 };
        siddonGPUFront->SetTransformMatrix(tmp);
        float* img = m_rawImage->GetBufferPointer();
        siddonGPUFront->SetImg3d(img, _3dpixelSpacing, _3dpixelNum);
        siddonGPUFront->SetImg2dParameter(spacingFront, imageSizeFront);
        QMessageBox::warning(this, "提示", "导入影像成功！");
    }
    catch (itk::ExceptionObject& err)
    {
        std::cerr << "ERROR: ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return;
    }    
}

void RegisterFineAdjustDialog::insertItems(QStringList items)
{
    if(items.size()==0)
        return;
    mItems.clear();
    mImageMap.clear();
    QWidget *container=new QWidget(this);
    container->setStyleSheet("background:#1D212C;border:none;");
    QVBoxLayout *l=new QVBoxLayout(container);
    l->setAlignment(Qt::AlignCenter);
    l->setSpacing(10);
    l->setContentsMargins(0,0,60,0);
    for(int i=0;i<items.size();i++)
    {
        QRadioButton *btn=new QRadioButton(container);
        mImageMap.emplace(btn,std::array<itk::Image<PixelType, 2>::Pointer,2>());
        connect(btn,&QRadioButton::clicked,this,[&]
        {
            QRadioButton *btn=static_cast<QRadioButton *>(sender());
            if(btn->isChecked())
            {
                checkChanged();
            }
        });
        mItems.push_back(btn);
        btn->setText(items[i]);
        l->addWidget(btn);
    }
    l->addStretch();
    ui->scrollArea->setWidget(container);
}

QStringList RegisterFineAdjustDialog::checkedItems()
{
    QStringList rst;
    for(auto ckbox:mItems)
    {
        if(ckbox->isChecked())
        {
            rst.push_back(ckbox->text());
            break;
        }
    }
    return rst;
}

void RegisterFineAdjustDialog::keyPressEvent(QKeyEvent *event)
{
    if(event->key()==Qt::Key_C && event->modifiers()==Qt::ControlModifier)
    {
        //复制当前参数
        mInfo[0]=ui->spinbox_translateX->value();
        mInfo[1]=ui->spinbox_translateY->value();
        mInfo[2]=ui->spinbox_translateZ->value();
        mInfo[3]=ui->spinbox_rotateX->value();
        mInfo[4]=ui->spinbox_rotateY->value();
        mInfo[5]=ui->spinbox_rotateZ->value();
    }
    else if(event->key()==Qt::Key_V)
    {
        //粘贴已复制的参数
/*        ui->horizontalSlider_transX->setValue(mInfo[0]);
        ui->horizontalSlider_transY->setValue(mInfo[1]);
        ui->horizontalSlider_transZ->setValue(mInfo[2]);
        ui->horizontalSlider_rotateX->setValue(mInfo[3]);
        ui->horizontalSlider_rotateY->setValue(mInfo[4]);
        ui->horizontalSlider_rotateZ->setValue(mInfo[5]);
*/
        ui->spinbox_translateX->setValue(mInfo[0]);
        ui->spinbox_translateY->setValue(mInfo[1]);
        ui->spinbox_translateZ->setValue(mInfo[2]);
        ui->spinbox_rotateX->setValue(mInfo[3]);
        ui->spinbox_rotateY->setValue(mInfo[4]);
        ui->spinbox_rotateZ->setValue(mInfo[5]);
    }
}
