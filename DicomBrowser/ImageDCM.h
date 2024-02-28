#ifndef IMAGEDCM_H
#define IMAGEDCM_H

#include "ITKFilter.h"

class ImageDCM
{
public:
	IntImageType::Pointer image;
	ImageDCM();
	ImageDCM(QString fileName, QVTKOpenGLNativeWidget*vtkWidget);
	virtual ~ImageDCM();
	bool setFileName(QString fileName);
	void setVTKWidget(QVTKOpenGLNativeWidget*vtkWidget);
	void openDCM(QString fileName, QVTKOpenGLNativeWidget*vtkWidget);

private:
	QString fileName;
	QVTKOpenGLNativeWidget*qvtkWidget;//��ʾDCMͼ���QVTKWidget�ؼ�

	ImageIOType::Pointer gdcmImageIO;
	ReaderType::Pointer reader;

	/** The VTK image to display in this window */
	vtkSmartPointer<vtkImageViewer2> viewer;
	vtkSmartPointer<vtkRenderer> renderer;
};
#endif
