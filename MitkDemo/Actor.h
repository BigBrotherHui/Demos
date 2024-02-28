#pragma once

#pragma execution_character_set("utf-8")
#include <QString>
#include <QSettings>
#include <QFile>
#include "QCryptographicHash"
#include "vtkObject.h"
#include "mitkDataStorage.h"
#include "mitkDataNode.h"
#include "vtkSmartPointer.h"
#include <vtkMatrix4x4.h>
#include <QtCore/QStringList>
#include <itkImage.h>
#include <itkImageSeriesReader.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <mitkITKImageImport.h>

class Actor : public vtkObject
{
public:
	static Actor* New();
	virtual ~Actor();
	vtkTypeMacro(Actor, vtkObject);

	QString mObjectName = "";
	float mOpacity = 1.0;
	int mLayer = 0;
		virtual int f_LoadActor();

	mitk::DataNode::Pointer outputNode{nullptr};
	mitk::DataStorage *f_GetDataStorage();
	void f_SetDataStorage(mitk::DataStorage *pDataStorage);
	virtual mitk::DataNode* f_GetDataNode();

	virtual vtkMatrix4x4 *f_GetOrgMatrix();

private:
	Actor(const Actor&) = delete;
	void operator=(const Actor&) = delete;
	mitk::DataStorage *m_GetDataStorage = nullptr;
	vtkMatrix4x4 *mGetOrgMatrix{nullptr};

protected:
	Actor();
private:
	QString mFileName;
	using PixelType = signed short;
	static constexpr unsigned int Dimension = 3;
	using ImageType = itk::Image<PixelType, Dimension>;
	using ReaderType = itk::ImageSeriesReader<ImageType>;
	using ImageIOType = itk::GDCMImageIO;
	using NamesGeneratorType = itk::GDCMSeriesFileNames;
	itk::ImageSeriesReader<ImageType>::Pointer mITKreader{ nullptr };
	ImageIOType::Pointer  gdcmIO{ nullptr };
};
