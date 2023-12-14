#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QHBoxLayout>
#include <vtkOutputWindow.h>
#include <QmitkStdMultiWidget.h>
#include <mitkIOUtil.h>
#include <QmitkRenderWindow.h>
#include <QmitkRenderWindowWidget.h>
#include <mitkDataNode.h>
#include <mitkPointSet.h>
#include <mitkPointSetDataInteractor.h>
#include <mitkNodePredicateDataType.h>
#include <mitkImage.h>
#include <QmitkRegisterClasses.h>
#include <QmitkLevelWindowWidget.h>
//#include <mitkDICOMEnums.h>
//#include <mitkDICOMFilesHelper.h>
//#include <mitkDICOMITKSeriesGDCMReader.h>
//#include <mitkDICOMTagsOfInterestHelper.h>
//#include <mitkDICOMProperty.h>
//#include <mitkDICOMDCMTKTagScanner.h>
//#include <mitkPropertyNameHelper.h>
//#include <mitkPropertyKeyPath.h>
//#include <mitkDICOMIOMetaInformationPropertyConstants.h>
//#include <mitkSplineVtkMapper3D.h>

#include <itkImageSeriesReader.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <mitkITKImageImport.h>

MainWindow::MainWindow(QWidget *parent)
    : WidgetBase(parent)
    , ui(new Ui::MainWindow)
{
    QmitkRegisterClasses();
    vtkOutputWindow::GlobalWarningDisplayOff();
    ui->setupUi(this);
    ds1=mitk::StandaloneDataStorage::New();
    /*ds2=mitk::StandaloneDataStorage::New();
    ds3=mitk::StandaloneDataStorage::New();*/
    w1=new QmitkStdMultiWidget(this);
    /*w2=new QmitkStdMultiWidget(this);
    rw=new QmitkRenderWindow(this);*/
    QHBoxLayout *l=new QHBoxLayout(ui->widget);
    l->addWidget(w1);
    QmitkLevelWindowWidget* level_window_widget = new QmitkLevelWindowWidget(this); 
    l->addWidget(level_window_widget);
    level_window_widget->SetDataStorage(ds1);
    /*l->addWidget(w2);
    l->addWidget(rw);*/
    l->setStretch(0,100);
    l->setStretch(1, 1);
    /*l->setStretch(2, 100);
    l->setStretch(3, 100);*/
    w1->SetDataStorage(ds1);
    //w2->SetDataStorage(ds2);
    //rw->GetRenderer()->SetDataStorage(ds3);
    //rw->GetSliceNavigationController()->SetDefaultViewDirection(mitk::AnatomicalPlane::Coronal);
    //rw->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard3D);
    w1->InitializeMultiWidget();
    //w2->InitializeMultiWidget();
    w1->AddPlanesToDataStorage();
    //w2->AddPlanesToDataStorage();
}

MainWindow::~MainWindow()
{
    delete ui;
}

//std::string GenerateNameFromDICOMProperties(const mitk::IPropertyProvider* provider)
//{
//    std::string nodeName = mitk::DataNode::NO_NAME_VALUE();
//
//    auto studyProp = provider->GetConstProperty(mitk::GeneratePropertyNameForDICOMTag(0x0008, 0x1030).c_str());
//    if (studyProp.IsNotNull())
//    {
//        nodeName = studyProp->GetValueAsString();
//    }
//
//    auto seriesProp = provider->GetConstProperty(mitk::GeneratePropertyNameForDICOMTag(0x0008, 0x103E).c_str());
//
//    if (seriesProp.IsNotNull())
//    {
//        if (studyProp.IsNotNull())
//        {
//            nodeName += " / ";
//        }
//        else
//        {
//            nodeName = "";
//
//        }
//        nodeName += seriesProp->GetValueAsString();
//    }
//
//    return nodeName;
//};

void MainWindow::on_pushButton_loaddata_clicked()
{
    std::string dicomSeriesPath="D:/Images/THA Datas/HIP_MODEL/ct/raw data";
    using TPixel = signed short;
    const unsigned int DIM3 = 3;
    using TImage = itk::Image<TPixel, DIM3>;
    using TImagePtr = TImage::Pointer;
    using TPoint = TImage::IndexType;
    using ImageType = TImage;
    using ReaderType = itk::ImageSeriesReader<ImageType>;
    using ImageIOType = itk::GDCMImageIO;
    using SeriesFileNamesType = itk::GDCMSeriesFileNames;
    using FileNamesContainer = std::vector<std::string>;
    using SeriesIdContainer = std::vector<std::string>;
    using SeriesIdContainer = std::vector<std::string>;

    ImageIOType::Pointer m_gdcmImageIO;
    SeriesFileNamesType::Pointer m_gdcmSeriesFileNames = SeriesFileNamesType::New();
    ReaderType::Pointer m_gdcmReader;
    m_gdcmReader = ReaderType::New();
    m_gdcmImageIO = ImageIOType::New();

    m_gdcmSeriesFileNames->SetUseSeriesDetails(true);
    m_gdcmSeriesFileNames->SetDirectory(dicomSeriesPath.c_str());
    const SeriesIdContainer& seriesUID = m_gdcmSeriesFileNames->GetSeriesUIDs();
    FileNamesContainer fileNames = m_gdcmSeriesFileNames->GetFileNames(seriesUID[0]);
    m_gdcmReader->SetImageIO(m_gdcmImageIO);
    m_gdcmReader->SetFileNames(fileNames);
    m_gdcmReader->ForceOrthogonalDirectionOff();
    try {
        // Read the files
        m_gdcmReader->Update();
        // Store the itk::Image pointer
        mitk::Image::Pointer mitkImage = mitk::Image::New();
        mitk::GrabItkImageMemory(m_gdcmReader->GetOutput(), mitkImage);
        mitk::DataNode::Pointer imgNode = mitk::DataNode::New();
        imgNode->SetData(mitkImage);
        ds1->Add(imgNode);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return;
    }
    //mitk::IOUtil::Load(dicomSeriesPath,*ds1);
    //mitk::IOUtil::Load("D:/Images/Pelvis_Right.stl", *ds2);
    /*auto subset=ds1->GetSubset(mitk::NodePredicateDataType::New("Image"));
    for(auto iter : *subset){
        iter->SetProperty("volumerendering",mitk::BoolProperty::New(1));
    }
    w1->ResetCrosshair();
    ds1->GetNamedNode("stdmulti.widget0.plane")->GetData()->Print(std::cout);
    w1->RequestUpdateAll();*/

    /*mitk::StringList relevantFiles = mitk::GetDICOMFilesInSameDirectory(dicomSeriesPath);
    mitk::DICOMITKSeriesGDCMReader::Pointer reader = mitk::DICOMITKSeriesGDCMReader::New();
    std::vector< std::string > m_ReadFiles;
    const unsigned int ntotalfiles = relevantFiles.size();
    for (unsigned int i = 0; i < ntotalfiles; i++)
    {
        m_ReadFiles.push_back(relevantFiles.at(i));
    }
    reader->SetAdditionalTagsOfInterest(mitk::GetCurrentDICOMTagsOfInterest());
    reader->SetTagLookupTableToPropertyFunctor(mitk::GetDICOMPropertyForDICOMValuesFunctor);
    reader->SetInputFiles(relevantFiles);

    mitk::DICOMDCMTKTagScanner::Pointer scanner = mitk::DICOMDCMTKTagScanner::New();
    scanner->AddTagPaths(reader->GetTagsOfInterest());
    scanner->SetInputFiles(relevantFiles);
    scanner->Scan();

    reader->SetTagCache(scanner->GetScanCache());
    reader->AnalyzeInputFiles();
    reader->LoadImages();
    std::vector<mitk::BaseData::Pointer> result;
    for (unsigned int i = 0; i < reader->GetNumberOfOutputs(); ++i)
    {
        const mitk::DICOMImageBlockDescriptor& desc = reader->GetOutput(i);
        mitk::BaseData::Pointer data = desc.GetMitkImage().GetPointer();

        std::string nodeName = GenerateNameFromDICOMProperties(&desc);

        mitk::StringProperty::Pointer nameProp = mitk::StringProperty::New(nodeName);
        data->SetProperty("name", nameProp);
        data->GetGeometry()->Print(std::cout);
        data->SetProperty(mitk::PropertyKeyPathToPropertyName(
            mitk::DICOMIOMetaInformationPropertyConstants::READER_CONFIGURATION()), mitk::StringProperty::New(reader->GetConfigurationLabel()));
        result.push_back(data);
    }
    for(int i=0;i<result.size();++i)
    {
        mitk::DataNode::Pointer dt = mitk::DataNode::New();
        dt->SetData(result[i]);
        ds1->Add(dt);
    }*/
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(ds1);
    w1->ResetCrosshair();
}


void MainWindow::on_pushButton_addpoint_clicked()
{
    //mitk::DataNode::Pointer dt=mitk::DataNode::New();
    //mitk::PointSet::Pointer ps=mitk::PointSet::New();
    //dt->SetData(ps);
    ////mitk::SplineVtkMapper3D::Pointer mapper = mitk::SplineVtkMapper3D::New();
    ////dt->SetMapper(mitk::BaseRenderer::Standard3D, mapper);
    //dt->SetProperty("pointsize",mitk::FloatProperty::New(10.));
    //mitk::PointSetDataInteractor::Pointer inter=mitk::PointSetDataInteractor::New();
    //inter->SetMaxPoints(10);
    //inter->LoadStateMachine("PointSet.xml");
    //inter->SetEventConfig("PointSetConfig.xml");
    //inter->SetDataNode(dt);
    //ds1->Add(dt);
    //ds2->Add(dt);
}

void MainWindow::on_pushButton_sharedTo2_clicked()
{
    //w2->SetDataStorage(ds2);
    //w2->ResetCrosshair();
}

void MainWindow::on_pushButton_sharedTo3_clicked()
{
    //w2->SetDataStorage(ds1);
    //w2->ResetCrosshair();
    /*rw->GetRenderer()->SetDataStorage(ds1);
    mitk::RenderingManager::GetInstance()->InitializeViewByBoundingObjects(rw->GetVtkRenderWindow(), ds1);*/
}
