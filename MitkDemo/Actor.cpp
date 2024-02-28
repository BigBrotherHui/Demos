#include "vtkObjectFactory.h"
#include <QDebug>
#include <QFileInfo>
#include <QDir>
#include "Actor.h"
#include <QMessageBox>

#include "ustatus.h"
#include <mitkProperties.h>

vtkStandardNewMacro(Actor);


Actor::Actor()
{
	mGetOrgMatrix = vtkMatrix4x4::New();
}

Actor::~Actor()
{
	if(mGetOrgMatrix)
	{
		mGetOrgMatrix->Delete();
		mGetOrgMatrix=nullptr;
	}
}

mitk::DataNode* Actor::f_GetDataNode()
{
	return nullptr;
}

vtkMatrix4x4* Actor::f_GetOrgMatrix()
{
	return mGetOrgMatrix;
}

mitk::DataStorage * Actor::f_GetDataStorage()
{
	return m_GetDataStorage;
}

void Actor::f_SetDataStorage(mitk::DataStorage *pDataStorage)
{
	m_GetDataStorage = pDataStorage;
}


int Actor::f_LoadActor()
{
	if (!outputNode)
		outputNode = mitk::DataNode::New();

	mFileName = uStatus::mMainDataDir + "/" + uStatus::mCurrentPatientID + uStatus::mDicomOrgDir;

	//判断是否存在文件
	std::vector< std::string >  mPaths;
	QDir dir(mFileName);
	dir.setFilter(QDir::Files | QDir::NoDotAndDotDot);
	QFileInfoList list = dir.entryInfoList();

	mObjectName = "DicomImage";
	
	mitk::DataNode* mTmpNode = f_GetDataNode();
	if (mTmpNode != nullptr)
	{
		f_GetDataStorage()->Remove(mTmpNode);
	}
	for (int i = 0; i < list.size(); i++)
	{
		QFileInfo mFileInfo = list.at(i);
		mPaths.push_back(mFileInfo.filePath().toStdString());
	}


	// 无法加载数据
	if (mPaths.size() == 0)
	{
		QString msg = "无法加载指定数据在 \r\n" + mFileName;
		QMessageBox::critical(nullptr, "错误", msg);
		return -1;
	}

	if (!mITKreader)
		mITKreader = ReaderType::New();
	if (!gdcmIO)
		gdcmIO = ImageIOType::New();
	mITKreader->SetImageIO(gdcmIO);
	mITKreader->SetFileNames(mPaths);
	mITKreader->Update();
	mLayer = 0;
	outputNode->SetProperty("name", mitk::StringProperty::New(mObjectName.toStdString().c_str()));
	outputNode->SetProperty("dicom", mitk::BoolProperty::New(true));
	outputNode->SetProperty("layer", mitk::IntProperty::New(mLayer));
	outputNode->SetProperty("volumerendering", mitk::BoolProperty::New(true));
	mitk::Image::Pointer image = mitk::ImportItkImage(mITKreader->GetOutput());
	outputNode->SetData(image);

	outputNode->Update();
	f_GetDataStorage()->Add(outputNode);

	if (outputNode && mGetOrgMatrix)
	{
		mGetOrgMatrix->DeepCopy(outputNode->GetData()->GetGeometry()->GetVtkMatrix());
		
	}
	return 0;
}