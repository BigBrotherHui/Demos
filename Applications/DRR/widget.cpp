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
#include <itkImageDuplicator.h>
#include <vtkColorTransferFunction.h>
#include <itkPNGImageIO.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkCamera.h>
#include <vtkFrustumSource.h>
#include <vtkPlaneSource.h>
#include <vtkTexture.h>
#include <vtkPiecewiseFunction.h>
#include <vtkNamedColors.h>
#include <vtkPlanes.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <itkMetaImageIOFactory.h>
#include <vtkInteractorStyleImage.h>
#include <vtkVolumeProperty.h>
#include <vtkSmartVolumeMapper.h>
Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)
{
	
    ui->setupUi(this);
	QHBoxLayout* layout = new QHBoxLayout(ui->widget_window);
	m_renderwindow2dside = new QVTKOpenGLNativeWidget(this);
	m_renderwindow2dside->setRenderWindow(rwside);
	rwside->AddRenderer(rside);
	rside->AddActor(actorside);
	layout->addWidget(m_renderwindow2dside);
	m_renderwindow2dside->interactor()->SetInteractorStyle(vtkInteractorStyleImage::New());
	m_renderwindow2dfront = new QVTKOpenGLNativeWidget(this);
	m_renderwindow2dfront->setRenderWindow(rwfront);
	rwfront->AddRenderer(rfront);
	rfront->AddActor(actorfront);	
	layout->addWidget(m_renderwindow2dfront);
	m_renderwindow2dfront->interactor()->SetInteractorStyle(vtkInteractorStyleImage::New());
}

Widget::~Widget()
{
    delete ui;
}

void Widget::requestUpdateAll()
{
	rwfront->Render();
	rwside->Render();
}

void Widget::setImage(bool front, itk::Image<unsigned char, 3>::Pointer img)
{
	itk::ImageToVTKImageFilter<itk::Image<unsigned char, 3>>::Pointer cast = itk::ImageToVTKImageFilter<itk::Image<unsigned char, 3>>::New();
	cast->SetInput(img);
	cast->Update();
	vtkImageData *vimg=cast->GetOutput();
	if (front)
		actorfront->SetInputData(vimg);
	else
		actorside->SetInputData(vimg);
	requestUpdateAll();
}

vtkSmartPointer<vtkPolyData> Widget::transformPolyData(vtkSmartPointer<vtkMatrix4x4> mt, vtkSmartPointer<vtkPolyData> p)
{
	if (!p)
		return nullptr;
	if (!mt)
		mt = vtkSmartPointer<vtkMatrix4x4>::New();
	vtkSmartPointer<vtkPolyData> ret = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->SetMatrix(mt);
	vtkSmartPointer<vtkTransformPolyDataFilter> fi = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	fi->SetTransform(transform);
	fi->SetInputData(p);
	fi->Update();
	ret->DeepCopy(fi->GetOutput());
	return ret;
}

void Widget::on_pushButton_importImage_clicked()
{
    QString dir=QFileDialog::getExistingDirectory(this, "选择图像", "D:/AA/");
    if (dir.isEmpty())
        return;
	using SeriesFileNamesType = itk::GDCMSeriesFileNames;
	using ImageType = itk::Image<float, 3>;
	using ReaderType = itk::ImageSeriesReader<ImageType>;
	using ImageIOType = itk::GDCMImageIO;
	ImageIOType::Pointer m_gdcmImageIO = ImageIOType::New();
	SeriesFileNamesType::Pointer m_gdcmSeriesFileNames = SeriesFileNamesType::New();
	m_gdcmSeriesFileNames->SetDirectory(dir.toStdString());
	using SeriesIdContainer = std::vector<std::string>;
	const SeriesIdContainer& seriesUID = m_gdcmSeriesFileNames->GetSeriesUIDs();
	const std::string& seriesUIDStr = seriesUID.at(0);
	using FileNamesContainer = std::vector<std::string>;
	FileNamesContainer fileNames = m_gdcmSeriesFileNames->GetFileNames(seriesUIDStr);
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetImageIO(m_gdcmImageIO);
	reader->SetFileNames(fileNames);
	try {
		reader->Update();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cerr << "ERROR: ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
	m_image = reader->GetOutput();
	
	auto frontImage = drr(m_image, "",ui->horizontalSlider_scd->value(), ui->horizontalSlider_rotate_x->value(), ui->horizontalSlider_rotate_y->value(), ui->horizontalSlider_rotate_z->value(),
		ui->horizontalSlider_translate_x->value(), ui->horizontalSlider_translate_y->value(), ui->horizontalSlider_translate_z->value());
	setImage(1, frontImage);
	requestUpdateAll();
}

void Widget::on_horizontalSlider_scd_valueChanged(int v)
{
	auto frontImage = drr(m_image, "", v, ui->horizontalSlider_rotate_x->value(), ui->horizontalSlider_rotate_y->value(), ui->horizontalSlider_rotate_z->value(),
		ui->horizontalSlider_translate_x->value(), ui->horizontalSlider_translate_y->value(), ui->horizontalSlider_translate_z->value());
	setImage(1, frontImage);
	requestUpdateAll();
}

void Widget::on_horizontalSlider_translate_x_valueChanged(int v)
{
	auto frontImage = drr(m_image, "",  ui->horizontalSlider_scd->value(), ui->horizontalSlider_rotate_x->value(), ui->horizontalSlider_rotate_y->value(), ui->horizontalSlider_rotate_z->value(),
		v,ui->horizontalSlider_translate_y->value(),ui->horizontalSlider_translate_z->value());
	setImage(1, frontImage);
	requestUpdateAll();
}

void Widget::on_horizontalSlider_translate_y_valueChanged(int v)
{
	auto frontImage = drr(m_image, "", ui->horizontalSlider_scd->value(), ui->horizontalSlider_rotate_x->value(), ui->horizontalSlider_rotate_y->value(), ui->horizontalSlider_rotate_z->value(),
		ui->horizontalSlider_translate_x->value(), v,  ui->horizontalSlider_translate_z->value());
	setImage(1, frontImage);
	requestUpdateAll();
}

void Widget::on_horizontalSlider_translate_z_valueChanged(int v)
{
	auto frontimage = drr(m_image, "", ui->horizontalSlider_scd->value(), ui->horizontalSlider_rotate_x->value(), ui->horizontalSlider_rotate_y->value(), ui->horizontalSlider_rotate_z->value(),
		ui->horizontalSlider_translate_x->value(), ui->horizontalSlider_translate_y->value(),v );
	auto sideimage = drr(m_image, "drr.png", v, 0, 0, -90, 0);
	
	setImage(1, frontimage);
	setImage(0, sideimage);
	
	requestUpdateAll();
}

void Widget::on_horizontalSlider_rotate_x_valueChanged(int v)
{
	auto frontImage = drr(m_image, "", ui->horizontalSlider_scd->value(), v, ui->horizontalSlider_rotate_y->value(), ui->horizontalSlider_rotate_z->value(),ui->horizontalSlider_translate_x->value()
		, ui->horizontalSlider_translate_y->value(), ui->horizontalSlider_translate_z->value());
	setImage(1, frontImage);
	requestUpdateAll();
}

void Widget::on_horizontalSlider_rotate_y_valueChanged(int v)
{
	auto frontImage = drr(m_image, "", ui->horizontalSlider_scd->value(), ui->horizontalSlider_rotate_x->value(), v,  ui->horizontalSlider_rotate_z->value(), ui->horizontalSlider_translate_x->value()
		, ui->horizontalSlider_translate_y->value(), ui->horizontalSlider_translate_z->value());
	setImage(1, frontImage);
	requestUpdateAll();
}

void Widget::on_horizontalSlider_rotate_z_valueChanged(int v)
{
	auto frontImage = drr(m_image, "", ui->horizontalSlider_scd->value(), ui->horizontalSlider_rotate_x->value(), ui->horizontalSlider_rotate_y->value(),v , ui->horizontalSlider_translate_x->value()
		, ui->horizontalSlider_translate_y->value(), ui->horizontalSlider_translate_z->value());
	setImage(1, frontImage);
	requestUpdateAll();
}
