#include "widget.h"
#include "ui_widget.h"
#include <QDir>
#include <QDebug>
#include <QFileDialog>
#include <itkGDCMSeriesFileNames.h>
#include <itkImageSeriesReader.h>
#include <itkGDCMImageIO.h>
#include "drr.hpp"

#include <itkImageToVTKImageFilter.h>
#include <vtkImageActor.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkGDCMImageIO.h"

#include <vtkSmartPointer.h>
#include <vtkImageActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkOutputWindow.h>
#include <vtkConeSource.h>
#include <vtkCubeSource.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <mitkITKImageImport.h>
#include <mitkSurface.h>
#include <mitkNodePredicateDataType.h>
#include <mitkTransferFunctionProperty.h>
#include <QmitkStdMultiWidget.h>
#include <mitkIOUtil.h>
#include <mitkDisplayActionEventHandlerStd.h>
#include <itkImageDuplicator.h>
#include <itkPNGImageIO.h>

Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)
{
	
    ui->setupUi(this);
	m_DisplayActionEventBroadcast = mitk::DisplayActionEventBroadcast::New();
	m_DisplayActionEventBroadcast->LoadStateMachine("DisplayInteraction.xml");
	m_DisplayActionEventBroadcast->SetEventConfig("DisplayConfigMITKBase.xml");
	//m_DisplayActionEventBroadcast->AddEventConfig("DisplayConfigCrosshair.xml");
	m_DisplayActionEventHandler = std::make_unique<mitk::DisplayActionEventHandlerStd>();
	m_DisplayActionEventHandler->SetObservableBroadcast(m_DisplayActionEventBroadcast);
	m_DisplayActionEventHandler->InitActions();
	vtkOutputWindow::GlobalWarningDisplayOff();
	m_renderwindow = new QmitkRenderWindow(this);
	m_renderwindow2d = new QmitkRenderWindow(this);
	m_data_storage_ = mitk::StandaloneDataStorage::New();
	m_data_storage_2d= mitk::StandaloneDataStorage::New();
	m_renderwindow->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard3D);
	m_renderwindow->GetRenderer()->SetDataStorage(m_data_storage_);
	m_renderwindow2d->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard2D);
	m_renderwindow2d->GetSliceNavigationController()->SetDefaultViewDirection(mitk::AnatomicalPlane::Axial);
	m_renderwindow2d->GetRenderer()->SetDataStorage(m_data_storage_2d);
    QList<int> sizes;
	sizes << 900 << 100;
    ui->splitter->setSizes(sizes);
    QHBoxLayout* l = new QHBoxLayout(ui->widget_window);
    l->addWidget(m_renderwindow);
	m_lw= new QmitkLevelWindowWidget(this);
	m_lw->SetDataStorage(m_data_storage_);
	m_lw->GetManager()->LevelWindowChanged.AddListener(mitk::MessageDelegate1<Widget, const mitk::LevelWindow&>(this, &Widget::levelWindowChanged));
	l->addWidget(m_lw);
	l->addWidget(m_renderwindow2d);
	l->setStretchFactor(m_renderwindow, 5);
	l->setStretchFactor(m_lw, 1);
	l->setStretchFactor(m_renderwindow2d, 5);
}

Widget::~Widget()
{
    delete ui;
}

void Widget::levelWindowChanged(const mitk::LevelWindow& levelWindow)
{
	auto subset = m_data_storage_->GetSubset(mitk::NodePredicateDataType::New("Image"));
	for(auto pNode : *subset)
	{
		bool v;
		pNode->GetBoolProperty("volumerendering", v);
		if(!v)
		{
			continue;
		}
		double opacity = 1;
		auto levelWindow = m_lw->GetManager()->GetLevelWindow();
		auto level = levelWindow.GetLevel();
		auto window = levelWindow.GetWindow();
		auto min = level - window / 2;
		auto perWindow = window / 8;
		pNode->SetProperty("volumerendering.usegpu", mitk::BoolProperty::New(false));
		pNode->SetProperty("volumerendering.usemip", mitk::BoolProperty::New(false));
		pNode->SetProperty("volumerendering.blendmode", mitk::IntProperty::New(0));  //使用混合形式的投影
		auto transferFunction = mitk::TransferFunction::New();
		double scalarOpacities[8] = { 0, 0.001846, 0.024414, 0.100113, 0.467300, 0.711914, 0.915909, 1 };
		transferFunction->GetScalarOpacityFunction()->RemoveAllPoints();
		for (int i = 0; i < 8; i++)
		{
			auto x = min + i * perWindow;
			if (i >= 4)
			{
				x += perWindow;
			}
			auto y = scalarOpacities[i] * opacity;
			transferFunction->GetScalarOpacityFunction()->AddPoint(x, y);
		}
		transferFunction->GetScalarOpacityFunction()->Modified();

		transferFunction->GetGradientOpacityFunction()->RemoveAllPoints();
		transferFunction->GetGradientOpacityFunction()->AddPoint(0, 1.000000 * opacity);
		transferFunction->GetGradientOpacityFunction()->Modified();

		transferFunction->GetColorTransferFunction()->RemoveAllPoints();
		transferFunction->GetColorTransferFunction()->AddRGBPoint(312.382940, 1.000000, 0.564706, 0.274510);
		transferFunction->GetColorTransferFunction()->AddRGBPoint(455.103448, 1.000000, 0.945098, 0.768627);
		transferFunction->GetColorTransferFunction()->AddRGBPoint(623.773140, 1.000000, 0.800000, 0.333333);
		transferFunction->GetColorTransferFunction()->AddRGBPoint(796.767695, 1.000000, 0.901961, 0.815686);
		transferFunction->GetColorTransferFunction()->AddRGBPoint(930.838475, 1.000000, 1.000000, 1.000000);
		transferFunction->GetColorTransferFunction()->AddRGBPoint(1073.558984, 1.000000, 0.839216, 0.423529);
		transferFunction->GetColorTransferFunction()->AddRGBPoint(1220.604356, 1.000000, 0.772549, 0.490196);
		transferFunction->GetColorTransferFunction()->Modified();
		pNode->SetProperty("TransferFunction", mitk::TransferFunctionProperty::New(transferFunction.GetPointer()));
	}
	mitk::RenderingManager::GetInstance()->RequestUpdate(m_renderwindow->GetVtkRenderWindow());
}

void Widget::on_pushButton_importImage_clicked()
{
    QString dir=QFileDialog::getExistingDirectory(this, "选择图像序列", "D:/Images/");
    if (dir.isEmpty())
        return;
	using SeriesFileNamesType = itk::GDCMSeriesFileNames;
	using ImageType = itk::Image<short, 3>;
	using ReaderType = itk::ImageSeriesReader<ImageType>;
	using ImageIOType = itk::GDCMImageIO;
	ImageIOType::Pointer m_gdcmImageIO = ImageIOType::New();
	SeriesFileNamesType::Pointer m_gdcmSeriesFileNames = SeriesFileNamesType::New();
	m_gdcmSeriesFileNames->SetUseSeriesDetails(true);
	// m_gdcmSeriesFileNames->AddSeriesRestriction("0008|0021");
	m_gdcmSeriesFileNames->SetDirectory(dir.toStdString());
	using SeriesIdContainer = std::vector<std::string>;
	const SeriesIdContainer& seriesUID = m_gdcmSeriesFileNames->GetSeriesUIDs();
	const std::string& seriesUIDStr = seriesUID.at(0);
	using FileNamesContainer = std::vector<std::string>;
	FileNamesContainer fileNames = m_gdcmSeriesFileNames->GetFileNames(seriesUIDStr);
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetImageIO(m_gdcmImageIO);
	reader->SetFileNames(fileNames);
	reader->ForceOrthogonalDirectionOff();
	try {
		reader->Update();
		using DuplicatorType = itk::ImageDuplicator<ImageType>;
		auto duplicator = DuplicatorType::New();
		duplicator->SetInputImage(reader->GetOutput());
		duplicator->Update();
		m_image = duplicator->GetOutput();

		ImageType::Pointer clonedImage = duplicator->GetOutput();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cerr << "ERROR: ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
	// image construct of mitk
	mitk::Image::Pointer mitkImage = mitk::Image::New();
	mitk::GrabItkImageMemory(reader->GetOutput(), mitkImage);
	mitk::DataNode::Pointer dt = mitk::DataNode::New();
	dt->SetName("image");
	dt->SetBoolProperty("volumerendering", 1);
	dt->SetData(mitkImage);
	m_data_storage_->Add(dt);

	float scd = 200;
	drr(reader->GetOutput(), "drr.png", 200);

	mitk::IOUtil::Load("drr.png", *m_data_storage_2d);
	mitk::RenderingManager::GetInstance()->InitializeViewByBoundingObjects(m_renderwindow2d->GetVtkRenderWindow(), m_data_storage_2d);
	/*using PixelType = unsigned char;
	const unsigned int  Dimension = 2;
	typedef itk::Image<PixelType, Dimension> DRRImageType;
	typedef itk::ImageFileReader<DRRImageType> DRRReaderType;*/
	//DRRReaderType::Pointer drrreader = DRRReaderType::New();
	//drrreader->SetFileName("drr.dcm");
	//using ImageIOType = itk::GDCMImageIO;
	//ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
	//drrreader->SetImageIO(gdcmImageIO);
	//try
	//{
	//	drrreader->Update();
	//}
	//catch (itk::ExceptionObject& e)
	//{
	//	std::cerr << "exception in file reader" << std::endl;
	//	std::cerr << e << std::endl;
	//	return;
	//}
	//mitk::Image::Pointer mitkImagedrr = mitk::Image::New();
	//mitk::GrabItkImageMemory(drrreader->GetOutput(), mitkImagedrr);
	//mitk::DataNode::Pointer dtdrr = mitk::DataNode::New();
	//dtdrr->SetName("drr");
	//dtdrr->SetBoolProperty("volumerendering", 1);
	//dtdrr->SetData(mitkImagedrr);
	//m_data_storage_->Add(dtdrr);

	//typedef itk::ImageToVTKImageFilter<DRRImageType> ConnectorType;
	//ConnectorType::Pointer connector = ConnectorType::New();
	//connector->SetInput(drrreader->GetOutput());
	//try
	//{
	//	connector->Update();
	//}
	//catch (itk::ExceptionObject& e)
	//{
	//	std::cerr << "exception in file reader" << std::endl;
	//	std::cerr << e << std::endl;
	//	return;
	//}


	vtkNew<vtkConeSource> camCS;
	camCS->SetHeight(1.5);
	camCS->SetResolution(12);
	camCS->SetRadius(0.4);

	vtkNew<vtkCubeSource> camCBS;
	camCBS->SetXLength(1.5);
	camCBS->SetZLength(0.8);
	camCBS->SetCenter(0.4, 0, 0);

	vtkNew<vtkAppendPolyData> camAPD;
	camAPD->AddInputConnection(camCS->GetOutputPort());
	camAPD->AddInputConnection(camCBS->GetOutputPort());
	camAPD->Update();

	mitk::Surface::Pointer sur = mitk::Surface::New();
	sur->SetVtkPolyData(camAPD->GetOutput());
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->Translate(0, 0, scd);
	transform->RotateY(-90);
	transform->Scale(50, 50, 50);
	sur->GetGeometry()->SetIndexToWorldTransformByVtkMatrix(transform->GetMatrix());
	mitk::DataNode::Pointer surNode = mitk::DataNode::New();
	surNode->SetData(sur);
	surNode->SetName("sur");
	m_data_storage_->Add(surNode);
	
	mitk::RenderingManager::GetInstance()->InitializeViewByBoundingObjects(m_renderwindow->GetVtkRenderWindow(), m_data_storage_);
}

void Widget::on_horizontalSlider_translate_z_valueChanged(int v)
{
	drr(m_image, "drr.png", v,0,0,-90);
	//m_data_storage_->Remove(m_data_storage_->GetNamedNode("drr"));
	//mitk::IOUtil::Load("drr.png", *m_data_storage_2d);
	using PixelType = unsigned char;
	const unsigned int  Dimension = 2;
	typedef itk::Image<PixelType, Dimension> DRRImageType;
	typedef itk::ImageFileReader<DRRImageType> DRRReaderType;
	DRRReaderType::Pointer drrreader = DRRReaderType::New();
	drrreader->SetImageIO(itk::PNGImageIO::New());
	drrreader->SetFileName("drr.png");
	try
	{
		drrreader->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << "exception in file reader" << std::endl;
		std::cerr << e << std::endl;
		return;
	}
	m_data_storage_2d->GetNamedObject<mitk::Image>("drr")->SetImportChannel(drrreader->GetOutput()->GetBufferPointer());
	
	mitk::Surface::Pointer sur =m_data_storage_->GetNamedObject<mitk::Surface>("sur");
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->Translate(0, 0, v);
	transform->RotateZ(-90);
	transform->Scale(50, 50, 50);
	sur->GetGeometry()->SetIndexToWorldTransformByVtkMatrix(transform->GetMatrix());
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}
