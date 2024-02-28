#include "mainwindow.h"
#include "ui_mainwindow.h"
//Qt
#include <QProgressDialog>
#include <QFileDialog>
#include <QMessageBox>
#include <QLabel>
#include <QListWidgetItem>
#include <QDebug>
#include <QTextCodec>

#include <QMouseEvent>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkPointPicker.h>
#include <itkNiftiImageIO.h>
#include <itkImageFileReader.h>
#include <itkImageToVTKImageFilter.h>
#include <vtkImageFlip.h>
#include <vtkImageCast.h>
#include <vtkImageProperty.h>
#include <QTimer>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    mrenderwindow = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    mrenderer = vtkSmartPointer<vtkRenderer>::New();
    mrenderwindow->AddRenderer(mrenderer);
    ui->widget->setRenderWindow(mrenderwindow);
    mstyle = vtkSmartPointer<LineInteractorStyle2>::New();
    ui->widget->GetInteractor()->SetInteractorStyle(mstyle);
	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
	plane->SetCenter(0, 0, 0);
	plane->SetNormal(0.0, 0, 1);
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(plane->GetOutputPort());

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(1.0, 0, 0);
	//actor->SetScale(5 * space[0], 5 * space[1], 1);//space是像素间隔
	actor->SetPickable(false);
	actor->GetProperty()->SetRepresentationToWireframe();
	actor->GetProperty()->SetEdgeColor(1.0, 1.0, 0.0);
	actor->GetProperty()->SetEdgeVisibility(true);
	actor->GetProperty()->SetLineWidth(2.0);
	actor->GetProperty()->SetRenderLinesAsTubes(true);

	mrenderer->AddActor(actor);
    mstyle->actor = actor;
	vtkSmartPointer<vtkPointPicker> picker = vtkSmartPointer<vtkPointPicker>::New();
	picker->SetTolerance(1e-10);
	ui->widget->GetInteractor()->SetPicker(picker);
    mimageactor = vtkSmartPointer<vtkImageActor>::New();
    mstyle->actorImage = mimageactor;
    mimageactor->SetPickable(1);
    mrenderer->AddActor(mimageactor);
    mimageactor2 = vtkSmartPointer<vtkImageActor>::New();
    mimageactor2->SetPickable(0);
    mimageactor2->SetOpacity(.2);
    mrenderer->AddActor(mimageactor2);
    actor->SetVisibility(0);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_openFile_clicked()
{
    QString path=QFileDialog::getOpenFileName(this, "open file", "D:/Images","*.nii");
    if (path.isEmpty())
        return;
    typedef itk::NiftiImageIO   ImageIOType;
    ImageIOType::Pointer NiftiImageIO = ImageIOType::New();
    const int Dimension = 3;
    typedef float PixelType;
    typedef itk::Image<PixelType, Dimension> ImageType;
    itk::ImageFileReader<ImageType>::Pointer reader = itk::ImageFileReader<ImageType>::New();
    reader->SetFileName(path.toStdString());
    reader->SetImageIO(ImageIOType::New());
    try
    {
        reader->Update();
    }
    catch (const itk::ExceptionObject& ex) {
        std::cout << ex << std::endl;
        return;
    }
    itk::ImageToVTKImageFilter<ImageType>::Pointer importer = itk::ImageToVTKImageFilter<ImageType>::New();
    importer->SetInput(reader->GetOutput());
    try
    {
        importer->Update();
    }
    catch (const itk::ExceptionObject& ex)
    {
        std::cout << ex << std::endl;
        return;
    }
    vtkImageData* vtkimage = importer->GetOutput();

    //由于nii坐标系问题，会出现图像向上去，需要反转y轴
    //vtkNew<vtkImageFlip> imageFlip;
    //imageFlip->SetInputData(vtkimage);
    //imageFlip->SetFilteredAxes(1);
    //imageFlip->Update();

    //vtkNew<vtkImageFlip> imageFlip1;
    //imageFlip1->SetInputConnection(imageFlip->GetOutputPort());
    //imageFlip1->SetFilteredAxes(0);
    //imageFlip1->Update();

    //vtkNew<vtkImageCast> imageCast;
    //imageCast->SetInputConnection(imageFlip1->GetOutputPort());
    //imageCast->SetOutputScalarTypeToShort();
    //imageCast->Update();

    mimagedata = vtkSmartPointer<vtkImageData>::New();
    mimagedata->DeepCopy(vtkimage);
    mimageactor->SetInputData(mimagedata);

    mrenderwindow->Render();
    mrenderer->ResetCameraClippingRange();

}

void MainWindow::on_pushButton_openFile_2_clicked()
{
    QString path = QFileDialog::getOpenFileName(this, "open file", "D:/Images", "*.nii");
    if (path.isEmpty())
        return;
    typedef itk::NiftiImageIO   ImageIOType;
    ImageIOType::Pointer NiftiImageIO = ImageIOType::New();
    const int Dimension = 3;
    typedef float PixelType;
    typedef itk::Image<PixelType, Dimension> ImageType;
    itk::ImageFileReader<ImageType>::Pointer reader = itk::ImageFileReader<ImageType>::New();
    reader->SetFileName(path.toStdString());
    reader->SetImageIO(ImageIOType::New());
    try
    {
        reader->Update();
    }
    catch (const itk::ExceptionObject& ex) {
        std::cout << ex << std::endl;
        return;
    }
    itk::ImageToVTKImageFilter<ImageType>::Pointer importer = itk::ImageToVTKImageFilter<ImageType>::New();
    importer->SetInput(reader->GetOutput());
    try
    {
        importer->Update();
    }
    catch (const itk::ExceptionObject& ex)
    {
        std::cout << ex << std::endl;
        return;
    }
    vtkImageData* vtkimage = importer->GetOutput();

    mimagedata2 = vtkSmartPointer<vtkImageData>::New();
    mimagedata2->DeepCopy(vtkimage);
    mimageactor2->SetInputData(mimagedata2);
    mrenderwindow->Render();    
}

void MainWindow::on_pushButton_show_clicked()
{
    mstyle->actor->SetVisibility(1);
    mrenderwindow->Render();
}
