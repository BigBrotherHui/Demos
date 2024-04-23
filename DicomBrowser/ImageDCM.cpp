#include "ImageDCM.h"

ImageDCM::ImageDCM()
{
	image = IntImageType::New();
	gdcmImageIO = ImageIOType::New();
	reader = ReaderType::New();
	reader->SetImageIO(gdcmImageIO);
	renderer = vtkSmartPointer<vtkRenderer>::New();
	viewer = vtkSmartPointer<vtkImageViewer2>::New();
}

ImageDCM::ImageDCM(QString fileName, QVTKOpenGLNativeWidget*vtkWidget)
{
	ImageDCM();
	this->fileName = fileName;
	this->qvtkWidget = vtkWidget;
	this->qvtkWidget->renderWindow()->AddRenderer(renderer);
}

ImageDCM::~ImageDCM()
{
}

void ImageDCM::openDCM(QString filePath, QVTKOpenGLNativeWidget*vtkWidget)
{
	bool flag=setFileName(filePath);
	setVTKWidget(vtkWidget);
	this->qvtkWidget->renderWindow()->AddRenderer(renderer);//设置qvtkWidget的渲染器

	if (flag)
	{
		//itk获取图像,vtk显示在qvtkWidget里
		ITKFilter::Painting(image, qvtkWidget);
	}
}

bool ImageDCM::setFileName(QString fileName)
{
	this->fileName = fileName;
	if (!fileName.isEmpty())
	{
		//QString转std::string
		QByteArray b;
		b.append(fileName);
		char *str = b.data();

		reader->SetFileName(str);
		reader->Update();
		image = reader->GetOutput();
		return true;
	}
	else
	{
		QMessageBox::warning(NULL, "Error", "图片地址无效", QMessageBox::Ok);
		return false;
	}

}

void ImageDCM::setVTKWidget(QVTKOpenGLNativeWidget* vtkWidget)
{
	this->qvtkWidget = vtkWidget;
}
