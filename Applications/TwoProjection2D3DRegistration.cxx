#ifndef GLOG_NO_ABBREVIATED_SEVERITIES
#define GLOG_NO_ABBREVIATED_SEVERITIES // 如果不加这个宏定义代码就会报错
#endif
#pragma comment(lib, "glog.lib")
#include <glog\logging.h>
#include "itkTwoProjectionImageRegistrationMethod.h"
#include "itkEuler3DTransform.h"
#include "itkNormalizedCorrelationTwoImageToOneImageMetric.h"
#include "itkSiddonJacobsRayCastInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkCommand.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkGiplImageIOFactory.h"
#include "itkSize.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkImage.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include<cmath>
#include "itkJPEGImageIOFactory.h"
#include "itkNIFTIImageIO.h"
#include "itkNIFTIImageIOFactory.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryBallStructuringElement.h"//基本球形
#include "itkBinaryDilateImageFilter.h"//二值膨胀
#include <json/json.h>
//#include "itkPowellOptimizertest.h"
#include "itkPowellOptimizer.h"
// mhd格式图像
#include <itkMetaImageIOFactory.h>
#include <time.h>
// First we define the command class to allow us to monitor the registration.
using InternalPixelType = float;
constexpr unsigned int Dimension = 3;
using ImageType = itk::Image<InternalPixelType, Dimension>;

class CommandIterationUpdate : public itk::Command
{
public:
	using Self = CommandIterationUpdate;
	using Superclass = itk::Command;
	using Pointer = itk::SmartPointer<Self>;
	itkNewMacro(Self);

protected:
	CommandIterationUpdate() = default;

public:
	using OptimizerType = itk::PowellOptimizer;
	using OptimizerPointer = const OptimizerType *;

	void
		Execute(itk::Object * caller, const itk::EventObject & event) override
	{
		Execute((const itk::Object *)caller, event);
	}

	void
		Execute(const itk::Object * object, const itk::EventObject & event) override
	{
		auto optimizer = dynamic_cast<OptimizerPointer>(object);
		if (typeid(event) != typeid(itk::IterationEvent))
		{
			return;
		}
	}
};


bool isFileExists_stat(std::string& name) {
	struct stat buffer;
	bool bExit = (stat(name.c_str(), &buffer) == 0);
	if (!bExit)
	{
		LOG(ERROR) << name << " doesn't exist";
		std::cout << name << " doesn't exist" << std::endl;
	}
	return bExit;
}

bool downloadTXT(std::string file, double* p, int num)
{

	if (!isFileExists_stat(file))
		return false;
	FILE* fp;
	int i;
	const char* file_c = file.c_str();
	fp = fopen(file_c, "rt");
	if (fp == NULL)
	{
		printf("%s", "File does not exist!");
		exit(-1);
	}
	for (i = 0; i < num; i++)
	{
		fscanf(fp, "%lf", &p[i]);
	}
	fclose(fp);
	return true;
}

void writerTXT(std::string file, double data[][4])
{
	const char* file_c = file.c_str();
	int cnt = 0;
	FILE *fp = fopen(file_c, "wb");

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; ++j)
		{
			fprintf(fp, "%lf ", data[i][j]);
		}
		fprintf(fp, "%\n ");
	}
	fclose(fp);
}


void writerimage(ImageType::Pointer image, std::string str)
{
	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	itk::JPEGImageIOFactory::RegisterOneFactory();
	writer->SetFileName(str);
	writer->SetInput(image);
	try {
		writer->Update();
	}
	catch (itk::ExceptionObject &e) {
		std::cerr << e << std::endl;
	}
}

void WriterImageHomogenization(ImageType::Pointer image, std::string st)
{
	using OutputImageType = itk::Image<unsigned short, 3>;
	using RescaleFilterType = itk::RescaleIntensityImageFilter<ImageType, OutputImageType>;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(image);
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(4095);
	rescaleFilter->Update();
	OutputImageType::Pointer img = rescaleFilter->GetOutput();
	ImageType::SizeType inputSizess2 =
		img->GetLargestPossibleRegion().GetSize();
	unsigned short * buffer = img->GetBufferPointer();
	int len = inputSizess2[0] * inputSizess2[1] * inputSizess2[2];

	FILE *fpwrt = NULL;
	const char* file_c = st.c_str();
	fopen_s(&fpwrt, file_c, "wb+");
	if (fpwrt == NULL)
	{
		std::cout << "写入文件创建错误" << std::endl;
	}
	fwrite(buffer, sizeof(unsigned short), len, fpwrt);
	fclose(fpwrt);
}

void PatientourceToRobotMatrix(double * par, double* ds1torobot, double M1[][4])
{
	using TransformType = itk::Euler3DTransform<double>;
	TransformType::Pointer    transform = TransformType::New();
	transform->SetComputeZYX(true);
	double dtr = (atan(1.0) * 4.0) / 180.0;
	transform->SetRotation(dtr * par[0], dtr * par[1], dtr * par[2]);
	TransformType::MatrixType R = transform->GetMatrix();
	double m[4][4];

	for (int i = 0; i < 3; i++) {
		m[3][i] = 0.0;

	}

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
			m[i][j] = R[i][j];
	}
	m[3][3] = 1;
	m[0][3] = par[3];
	m[1][3] = par[4];
	m[2][3] = par[5];

	float Rs1toRB[4][4]{};
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			int ind = i * 4 + j;
			Rs1toRB[i][j] = ds1torobot[ind];
		}

	}


	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			M1[i][j] = 0;
			for (int k = 0; k < 4; k++)
				M1[i][j] += Rs1toRB[i][k] * m[k][j];
		}
	}
}


void saveimg(ImageType::Pointer img, const char* file)
{
	ImageType::SizeType inputSizes =
		img->GetLargestPossibleRegion().GetSize();
	float * buffer = img->GetBufferPointer();
	int len = inputSizes[0] * inputSizes[1] * inputSizes[2];

	savedata(buffer, len, file);
}


void Converted(ImageType::Pointer img)
{
	ImageType::SizeType inputSizess2 =
		img->GetLargestPossibleRegion().GetSize();
	float* buffer3 = img->GetBufferPointer();
	int k, f;
	float temp;
	int len2 = inputSizess2[0] * inputSizess2[1] * inputSizess2[2];
	for (int i = 0; i < inputSizess2[0]; i++) {
		for (int j = 0; j < i; j++) {//注意此处循环限制语句j<i
			k = i * inputSizess2[0] + j;
			f = j * inputSizess2[0] + i;
			temp = buffer3[k];
			buffer3[k] = buffer3[f];
			buffer3[f] = temp;
		}
	}

	for (int i = 0; i < inputSizess2[0]; i++) {
		for (int j = 0; j < inputSizess2[0] / 2; j++) {//注意此处循环限制语句j<i
			k = j * inputSizess2[0] + i;
			f = (inputSizess2[0] - j - 1) * inputSizess2[0] + i;
			temp = buffer3[k];
			buffer3[k] = buffer3[f];
			buffer3[f] = temp;
		}
	}
}

void DRRMask(float* rst1, int len)
{
	for (int i = 0; i < len; ++i)
	{
		if (std::abs(rst1[i]) < 0.0001)
		{
			rst1[i] = 0;
		}
		else
		{
			rst1[i] = 1;
		}
	}
}

ImageType::Pointer readMask(std::string filename, ImageType::SizeType size)
{
	ImageType::Pointer image = ImageType::New();
	ImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;  // first index on Y 图像Y维最初的像素值
	start[2] = 0;  // first index on Z 图像Z维最初的像素值
	ImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);
	image->SetRegions(region);
	image->Allocate();

	float * buffer2 = image->GetBufferPointer();
	int len = size[0] * size[1] * size[2];

	float * rst1 = new float[len];
	readdata(filename, len, rst1);
	memcpy(buffer2, rst1, len * sizeof(float));
	return image;
}

ImageType::Pointer imageinterpolation(ImageType::Pointer imginput,double isoSpacing)
{
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

void ConvertedInverse(ImageType::Pointer img)
{
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



void SubtractImage(ImageType::Pointer img)
{
	ImageType::SizeType inputSizess2 =
		img->GetLargestPossibleRegion().GetSize();
	float * buffer = img->GetBufferPointer();
	int len = inputSizess2[0] * inputSizess2[1] * inputSizess2[2];
	for (int i = 0; i < len; ++i)
	{
		buffer[i] = 255 - buffer[i];
	}
}


void CT3DNormalization(ImageType::Pointer img,double dThreshold)
{
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


std::vector<std::string> read_txt(const std::string &path)
{
	std::ifstream ifile(path);//读取文件
	std::streampos len = ifile.tellg();//获取文件长度
	std::vector<std::string> data(len);//保存文件
	std::string str1; //temp
	if (ifile) //文件打开是否成功
	{
		while (std::getline(ifile, str1)) {
			//std::cout << str1 << std::endl;
			data.push_back(str1);
		}
		ifile.close();
		return data;
	}
	else
	{
		std::cout << path +" doesn't exist" << std::endl;
		return {};
	}
}

std::vector<float> ReadCenterOff(const std::string &path)
{
	std::vector<std::string> data = read_txt(path);
	std::string sCenterOff = "CenterToImage = ";
	std::vector<float> vN;
	vN.clear();
	for (std::string st : data)
	{
		std::string::size_type position;
		position = st.find(sCenterOff);
		
		if (st.npos != position&& 0== position)
		{
			st.erase(st.begin(), st.begin() + sCenterOff.size()); 
			std::size_t size_t(1);
			float num_float;
			while (!st.empty())
			{
				num_float = std::stof(st, &size_t);
				vN.push_back(num_float);
				st.erase(st.begin(), st.begin() + size_t);
			}
		}
	}

	return vN;
}

ImageType::Pointer  ReadImageMask(std::string file)
{
	typedef itk::ImageFileReader<ImageType> ReaderType;
	itk::ObjectFactoryBase::RegisterFactory(itk::MetaImageIOFactory::New());
	ReaderType::Pointer maskreader = ReaderType::New();
	maskreader->SetFileName(file);  //当前文件夹下需要有3.mhd  3.raw两个文件
	maskreader->Update();
	return maskreader->GetOutput();
}

ImageType::Pointer  ReadDicom(std::string file)
{
	using ImageIOType = itk::GDCMImageIO;
	ImageIOType::Pointer dicomIO = ImageIOType::New();
	using ImageReaderType2D = itk::ImageFileReader<ImageType>;
	ImageReaderType2D::Pointer imageReader2D = ImageReaderType2D::New();
	imageReader2D->SetFileName(file);
	imageReader2D->SetImageIO(dicomIO);
	imageReader2D->Update();
	return  imageReader2D->GetOutput();
}


ImageType::Pointer creatBlankImage(ImageType::SizeType  tempsize)
{
	ImageType::Pointer imagetemp = ImageType::New();
	ImageType::IndexType start;
	start[0] = 0;  // first index on X 图像X维最初的像素值
	start[1] = 0;  // first index on Y 图像Y维最初的像素值
	start[2] = 0;  // first index on Z 图像Z维最初的像素值
	ImageType::RegionType region;
	region.SetSize(tempsize);
	region.SetIndex(start);
	imagetemp->SetRegions(region);
	imagetemp->Allocate();
	return imagetemp;
}

bool ReadJson(char *file, Json::Value &root)
{
	Json::Reader reader;
	std::ifstream srcFile(file, std::ios::binary);
	if (!srcFile.is_open())
	{
		//LOG(ERROR) << "未找到配置文件: " << file;
		return false;
	}
	reader.parse(srcFile, root);
	return true;
}


int registration(Json::Value& root)
{
	bool bModel(true);
	if (root["isModel"].isBool())
	{
		bModel = root["isModel"].asBool();
	}
	
	std::string apFile = root["apFile"].asString();
	std::string file3D = root["imagePath"].asString();
	std::string oblFile = root["oblFile"].asString();
	std::string temppath = root["path"].asString();
	std::string allmask = root["MaskInfo"].asString();
	int lumbarlen = root["list"].size();
	std::string* lumbar = new std::string[lumbarlen];
	for (int i = 0; i < lumbarlen; ++i)
	{
		lumbar[i] = root["list"][i].asString();
		if( !isFileExists_stat(root[lumbar[i] + "MaskInfo"].asString()))
		{
			return EXIT_FAILURE;
		}
	}
	std::string focal = temppath + "/s1tos2.txt";
	std::string fapmask = temppath + "/topMask.jpg";
	std::string foblmask = temppath + "/sideMask.jpg";
	std::string fInitialValue = temppath + "/InitialValue.txt";
	std::string fs1torobot = temppath + "/s1torobot.txt";
	double mtmp[16];
	if (!downloadTXT(focal, mtmp, 16)) //读取变换矩阵
		return EXIT_FAILURE;

	double ds1torobot[16];
	if (!downloadTXT(fs1torobot, ds1torobot, 16))//读取变换矩阵
		return EXIT_FAILURE;

	double InitialValue[6];
	if (!downloadTXT(fInitialValue, InitialValue, 6))//读取变换矩阵
		return EXIT_FAILURE;

	LOG(INFO) << "完成数据准备";
	
	// The following lines define each of the components used in the
	// registration: The transform, optimizer, metric, interpolator and
	// the registration method itself.
	//constexpr unsigned int Dimension = 3;
	using InternalImageType = itk::Image<InternalPixelType, Dimension>;

	using TransformType = itk::Euler3DTransform<double>;

	using OptimizerType = itk::PowellOptimizer;

	// using MetricType = itk::GradientDifferenceTwoImageToOneImageMetric<
	using MetricType = itk::NormalizedCorrelationTwoImageToOneImageMetric<InternalImageType, InternalImageType>;

	using InterpolatorType = itk::SiddonJacobsRayCastInterpolateImageFunction<InternalImageType, double>;

	using RegistrationType = itk::TwoProjectionImageRegistrationMethod<InternalImageType, InternalImageType>;

	//  The 2- and 3-D images are read from files,
	using ImageIOType = itk::GDCMImageIO;
	ImageIOType::Pointer dicomIO = ImageIOType::New();
	using ImageReaderType2D = itk::ImageFileReader<InternalImageType>;
	ImageReaderType2D::Pointer imageReader2D1 = ImageReaderType2D::New();
	ImageReaderType2D::Pointer imageReader2D2 = ImageReaderType2D::New();

	ImageReaderType2D::Pointer maskReader2D1 = ImageReaderType2D::New();
	ImageReaderType2D::Pointer maskReader2D2 = ImageReaderType2D::New();

	imageReader2D1->SetFileName(apFile);
	imageReader2D2->SetFileName(oblFile);

	imageReader2D1->SetImageIO(dicomIO);
	imageReader2D2->SetImageIO(dicomIO);
	imageReader2D1->Update();
	imageReader2D2->Update();
	InternalImageType::Pointer imageCArmAP = imageReader2D1->GetOutput();
	InternalImageType::Pointer imageCArmOBL = imageReader2D2->GetOutput();
	//saveimg(imageCArmOBL, "img.raw");
	Converted(imageCArmAP);
	Converted(imageCArmOBL);
	
	maskReader2D1->SetFileName(fapmask);   //读取 c臂 mask
	maskReader2D2->SetFileName(foblmask);
	itk::JPEGImageIOFactory::RegisterOneFactory();
	maskReader2D1->Update();
	maskReader2D2->Update();
	InternalImageType::Pointer mask1 = maskReader2D1->GetOutput();
	InternalImageType::Pointer mask2 = maskReader2D2->GetOutput();
	Converted(mask1);
	Converted(mask2);
	const InternalImageType::SpacingType tCArmSpacingType = imageCArmAP->GetSpacing();
	float fCArmSpacing = tCArmSpacingType[0];
	ImageType::Pointer image3D = ReadImageMask(allmask);
	std::vector<float> vN = ReadCenterOff(allmask);
	const ImageType::SpacingType  inputSpacing = image3D->GetSpacing();
	ImageType::Pointer image3DIn = imageinterpolation(image3D, inputSpacing[0]*2);
	double lthreshold(-1000);
	if (!bModel)
		lthreshold = 100;
	CT3DNormalization(image3DIn, lthreshold);
	// To simply Siddon-Jacob's fast ray-tracing algorithm, we force the origin of the CT image
	// to be (0,0,0). Because we align the CT isocenter with the central axis, the projection
	// geometry is fully defined. The origin of the CT image becomes irrelavent.
	ImageType::PointType image3DOrigin;
	image3DOrigin[0] = 0.0;
	image3DOrigin[1] = 0.0;
	image3DOrigin[2] = 0.0;
	image3DIn->SetOrigin(image3DOrigin);

	using ParametersType = RegistrationType::ParametersType;
	ParametersType finalParameters;

	int numberOfLevels(3);
	using RecursiveMultiResolutionPyramidImageFilterType =
		itk::RecursiveMultiResolutionPyramidImageFilter<InternalImageType, InternalImageType>;
	RecursiveMultiResolutionPyramidImageFilterType::Pointer recursiveMultiResolutionPyramidImageFilter1 =
		RecursiveMultiResolutionPyramidImageFilterType::New();
	RecursiveMultiResolutionPyramidImageFilterType::Pointer recursiveMultiResolutionPyramidImageFilter2 =
		RecursiveMultiResolutionPyramidImageFilterType::New();
	recursiveMultiResolutionPyramidImageFilter1->SetInput(imageCArmAP);
	recursiveMultiResolutionPyramidImageFilter1->SetNumberOfLevels(numberOfLevels);
	recursiveMultiResolutionPyramidImageFilter1->Update();
	recursiveMultiResolutionPyramidImageFilter2->SetInput(imageCArmOBL);
	recursiveMultiResolutionPyramidImageFilter2->SetNumberOfLevels(numberOfLevels);
	recursiveMultiResolutionPyramidImageFilter2->Update();
	RecursiveMultiResolutionPyramidImageFilterType::Pointer recursiveMultiResolutionPyramidMaskFilter1 =
		RecursiveMultiResolutionPyramidImageFilterType::New();
	RecursiveMultiResolutionPyramidImageFilterType::Pointer recursiveMultiResolutionPyramidMaskFilter2 =
		RecursiveMultiResolutionPyramidImageFilterType::New();
	recursiveMultiResolutionPyramidMaskFilter1->SetInput(mask1);
	recursiveMultiResolutionPyramidMaskFilter1->SetNumberOfLevels(numberOfLevels);
	recursiveMultiResolutionPyramidMaskFilter1->Update();
	recursiveMultiResolutionPyramidMaskFilter2->SetInput(mask2);
	recursiveMultiResolutionPyramidMaskFilter2->SetNumberOfLevels(numberOfLevels);
	recursiveMultiResolutionPyramidMaskFilter2->Update();



	//InitialValue[0] = -9.01138; // Convert radian to degree
	//InitialValue[1] = 1.35299;
	//InitialValue[2] = 91.7895;
	//InitialValue[3] = 812.054;
	//InitialValue[4] = 1.31795;
	//InitialValue[5] = -46.1212;

	int CoarseRegistration(3);
	for (int i = 0; i < CoarseRegistration + lumbarlen; ++i)
	{
		LOG(INFO) << "开始第"<<i+1<<"次配准";
		std::clock_t  tStart, tEnd;
		tStart = clock();
		int CurrentLumbar(0); //当前椎节
		int iLevels(0);
		if (0 == i)
		{
			iLevels = 0;
		}
		else if (CoarseRegistration > i)
		{
			iLevels = i - 1;
			//image3DIn = imageinterpolation(image3D, inputSpacing[0] * std::pow(2.0,2-iLevels));
			image3DIn = imageinterpolation(image3D, inputSpacing[0] * 2);
			CT3DNormalization(image3DIn, lthreshold);
		}
		else
		{
			iLevels = 2;
			CurrentLumbar = i - CoarseRegistration;

			std::string datafile = root[lumbar[CurrentLumbar] + "MaskInfo"].asString();
			image3DIn = ReadImageMask(datafile);
			const ImageType::SpacingType Spacing = image3DIn->GetSpacing();
			vN = ReadCenterOff(datafile);
			image3DIn = imageinterpolation(image3DIn, Spacing[0]);
			CT3DNormalization(image3DIn, lthreshold);
		}
		using RescaleFilterType = itk::RescaleIntensityImageFilter<InternalImageType, InternalImageType>;
		RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		rescaleFilter->SetInput(recursiveMultiResolutionPyramidImageFilter1->GetOutput(iLevels));
		rescaleFilter->SetOutputMinimum(0);
		rescaleFilter->SetOutputMaximum(255);
		rescaleFilter->Update();
		InternalImageType::Pointer img1 = rescaleFilter->GetOutput();
		SubtractImage(img1);

		RescaleFilterType::Pointer rescaleFilter2 = RescaleFilterType::New();
		rescaleFilter2->SetInput(recursiveMultiResolutionPyramidImageFilter2->GetOutput(iLevels));
		rescaleFilter2->SetOutputMinimum(0);
		rescaleFilter2->SetOutputMaximum(255);
		rescaleFilter2->Update();
		InternalImageType::Pointer img2 = rescaleFilter2->GetOutput();
		SubtractImage(img2);


		InternalImageType::Pointer masktemp1 = recursiveMultiResolutionPyramidMaskFilter1->GetOutput(iLevels);
		InternalImageType::Pointer masktemp2 = recursiveMultiResolutionPyramidMaskFilter2->GetOutput(iLevels);

		MetricType::Pointer       metric = MetricType::New();
		TransformType::Pointer    transform = TransformType::New();
		OptimizerType::Pointer    optimizer = OptimizerType::New();
		InterpolatorType::Pointer interpolator1 = InterpolatorType::New();
		InterpolatorType::Pointer interpolator2 = InterpolatorType::New();
		RegistrationType::Pointer registration = RegistrationType::New();
		metric->ComputeGradientOff();
		metric->SetSubtractMean(true);
		if (0 == i)
		{
			metric->SetActive(true, true);
		}
		else
		{
			metric->SetActive(true, true);
		}
		// and passed to the registration method:
		registration->SetMetric(metric);
		registration->SetOptimizer(optimizer);
		registration->SetTransform(transform);
		registration->SetInterpolator1(interpolator1);
		registration->SetInterpolator2(interpolator2);

		using FlipFilterType = itk::FlipImageFilter<InternalImageType>;
		FlipFilterType::Pointer flipFilter1 = FlipFilterType::New();
		FlipFilterType::Pointer flipFilter2 = FlipFilterType::New();

		using FlipAxesArrayType = FlipFilterType::FlipAxesArrayType;
		FlipAxesArrayType flipArray;
		flipArray[0] = false;
		flipArray[1] = false;
		flipArray[2] = false;

		flipFilter1->SetFlipAxes(flipArray);
		flipFilter2->SetFlipAxes(flipArray);

		flipFilter1->SetInput(img1);
		flipFilter2->SetInput(img2);

		// The input 2D images may have 16 bits. We rescale the pixel value to between 0-255.
		using Input2DRescaleFilterType = itk::RescaleIntensityImageFilter<InternalImageType, InternalImageType>;

		Input2DRescaleFilterType::Pointer rescaler2D1 = Input2DRescaleFilterType::New();
		rescaler2D1->SetOutputMinimum(0);
		rescaler2D1->SetOutputMaximum(255);
		rescaler2D1->SetInput(flipFilter1->GetOutput());
		Input2DRescaleFilterType::Pointer rescaler2D2 = Input2DRescaleFilterType::New();
		rescaler2D2->SetOutputMinimum(0);
		rescaler2D2->SetOutputMaximum(255);
		rescaler2D2->SetInput(flipFilter2->GetOutput());

		//  The 3D CT dataset is casted to the internal image type using
		//  {CastImageFilters}.

		using CastFilterType3D = itk::CastImageFilter<ImageType, InternalImageType>;

		CastFilterType3D::Pointer caster3D = CastFilterType3D::New();
		caster3D->SetInput(image3DIn);
		rescaler2D1->Update();
		rescaler2D2->Update();
		caster3D->Update();

		// registration
		ImageType::SizeType  tempsize = img1->GetLargestPossibleRegion().GetSize();
		int lenfix = tempsize[0] * tempsize[1];
		float* drrmask1 = new float[lenfix];
		float* drrmask2 = new float[lenfix];
		for (int i = 0; i < lenfix; ++i)
		{
			drrmask1[i] = 1;
			drrmask2[i] = 1;
		}
		registration->SetFixedImage1(rescaler2D1->GetOutput(), masktemp1);
		registration->SetFixedImage2(rescaler2D2->GetOutput(), masktemp2);
		registration->SetMovingImage(caster3D->GetOutput());

		interpolator1->SetDrrMask(drrmask1);
		interpolator2->SetDrrMask(drrmask2);
		//savedata(drrmask1, lenfix, "mask1.raw");
		// Initialise the transform
		// ~~~~~~~~~~~~~~~~~~~~~~~~

		// Set the order of the computation. Default ZXY
		transform->SetComputeZYX(true);

		// The transform is initialised with the translation [tx,ty,tz] and
		// rotation [rx,ry,rz] specified on the command line

		TransformType::OutputVectorType translation;
		translation[0] = InitialValue[3];
		translation[1] = InitialValue[4];
		translation[2] = InitialValue[5];
		transform->SetTranslation(translation);

		// constant for converting degrees to radians
		const double dtr = (atan(1.0) * 4.0) / 180.0;
		transform->SetRotation(dtr * InitialValue[0], dtr * InitialValue[1], dtr * InitialValue[2]);

		// The centre of rotation is set by default to the centre of the 3D
		// volume but can be offset from this position using a command
		// line specified translation [cx,cy,cz]

		ImageType::PointType       origin3D = image3DIn->GetOrigin();
		const itk::Vector<double, 3> resolution3D = image3DIn->GetSpacing();

		using ImageRegionType3D = ImageType::RegionType;
		using SizeType3D = ImageRegionType3D::SizeType;

		ImageRegionType3D region3D = caster3D->GetOutput()->GetBufferedRegion();
		SizeType3D        size3D = region3D.GetSize();

		TransformType::InputPointType isocenter;

		// Set the center of the image as the isocenter.
		isocenter[0] = origin3D[0] + resolution3D[0] * static_cast<double>(size3D[0]) / 2.0;
		isocenter[1] = origin3D[1] + resolution3D[1] * static_cast<double>(size3D[1]) / 2.0;
		isocenter[2] = origin3D[2] + resolution3D[2] * static_cast<double>(size3D[2]) / 2.0;
		
		transform->SetCenter(isocenter);

		// Note: Two 2D images may have different image sizes and pixel dimensions, although
		// scd are the same.

		const itk::Vector<double, 3> resolution2D1 = rescaler2D1->GetOutput()->GetSpacing();
		const itk::Vector<double, 3> resolution2D2 = rescaler2D1->GetOutput()->GetSpacing();

		using ImageRegionType2D = InternalImageType::RegionType;
		using SizeType2D = ImageRegionType2D::SizeType;

		ImageRegionType2D region2D1 = rescaler2D1->GetOutput()->GetBufferedRegion();
		ImageRegionType2D region2D2 = rescaler2D2->GetOutput()->GetBufferedRegion();
		SizeType2D        size2D1 = region2D1.GetSize();
		SizeType2D        size2D2 = region2D2.GetSize();

		registration->SetFixedImageRegion1(rescaler2D1->GetOutput()->GetBufferedRegion());
		registration->SetFixedImageRegion2(rescaler2D2->GetOutput()->GetBufferedRegion());

		// Initialize the ray cast interpolator

		// 2D Image 1
		interpolator1->SetTransform(transform);
		interpolator1->Initialize();

		// 2D Image 2
		interpolator2->SetTransform(transform);
		interpolator2->Initialize();

		interpolator1->SetPixelspacing(resolution2D1);
		interpolator1->SetImgSize(size2D1);
		interpolator2->SetPixelspacing(resolution2D2);
		interpolator2->SetImgSize(size2D1);

		//float FocalDistance = 5069.9*fCArmSpacing;
		float FocalDistance = 5042.04*fCArmSpacing;
		interpolator1->SetFocalDistance(FocalDistance);
		interpolator2->SetFocalDistance(FocalDistance);

		interpolator1->SetbObl(false);
		interpolator2->SetbObl(true);

		interpolator2->SetToS1M(mtmp);

		interpolator1->SetCenterOff(vN);
		interpolator2->SetCenterOff(vN);
		// Set up the transform and start position
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// The registration start position is intialised using the
		// transformation parameters.
		registration->SetInitialTransformParameters(transform->GetParameters());
		metric->SetInitialTransform(transform->GetParameters());
		// The optimizer weightings are set such that one degree equates to
		// one millimeter.
		itk::Optimizer::ScalesType weightings(transform->GetNumberOfParameters());
		if (CoarseRegistration  <= i)
		{
			optimizer->SetMaximize(true);  // for GradientDifferenceTwoImageToOneImageMetric
			optimizer->SetMaximumIteration(11);
			optimizer->SetMaximumLineIteration(9); // for Powell's method
			optimizer->SetStepLength(1);
			optimizer->SetStepTolerance(0.0002);
			optimizer->SetValueTolerance(0.000001);
			weightings[0] = 1. / dtr*1 ;
			weightings[1] = 1. / dtr*1;
			weightings[2] = 1. / dtr*1;
			weightings[3] = 1.;
			weightings[4] = 1.;
			weightings[5] = 1.;
		}
		else
		{
			optimizer->SetMaximize(true);  // for GradientDifferenceTwoImageToOneImageMetric
			optimizer->SetMaximumIteration(11);
			optimizer->SetMaximumLineIteration(9); // for Powell's method
			optimizer->SetStepLength(15 / (iLevels*0.5 + 1));
			optimizer->SetStepTolerance(0.002);
			optimizer->SetValueTolerance(0.0001 / (iLevels * 10 + 1));
			weightings[0] = 1. / dtr * 15;
			weightings[1] = 1. / dtr * 15;
			weightings[2] = 1. / dtr * 15;
			weightings[3] = 1.;
			weightings[4] = 1.;
			weightings[5] = 1.;
		}

		optimizer->SetScales(weightings);
		// Create the observers
		// ~~~~~~~~~~~~~~~~~~~~
		CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
		optimizer->AddObserver(itk::IterationEvent(), observer);
		try
		{
			registration->StartRegistration();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cout << "ExceptionObject caught !" << std::endl;
			LOG(ERROR) << "第" << i + 1 << "次配准错误"<< err;
			std::cout << err << std::endl;
			return -1;
		}
		
		//单椎节时使用整体为初始值，不更新
		finalParameters = registration->GetLastTransformParameters();
		//if (CoarseRegistration > i)
		//if (0 == i)
		//{
		//	InitialValue[3] = finalParameters[3];
		//}
		//else
		{
			InitialValue[0] = finalParameters[0] / dtr; // Convert radian to degree
			InitialValue[1] = finalParameters[1] / dtr;
			InitialValue[2] = finalParameters[2] / dtr;
			InitialValue[3] = finalParameters[3];
			InitialValue[4] = finalParameters[4];
			InitialValue[5] = finalParameters[5];
		}
		double FinalValue[6];
		FinalValue[0] = finalParameters[0] / dtr; // Convert radian to degree
		FinalValue[1] = finalParameters[1] / dtr;
		FinalValue[2] = finalParameters[2] / dtr;
		FinalValue[3] = finalParameters[3];
		FinalValue[4] = finalParameters[4];
		FinalValue[5] = finalParameters[5];

		const int numberOfIterations = optimizer->GetCurrentIteration();
		const double bestValue = optimizer->GetValue();
		tEnd = clock();
		LOG(INFO) << " Result = " << finalParameters[0] / dtr<<"  "<<finalParameters[1] / dtr<< "  " << finalParameters[2] / dtr << "  " << finalParameters[3] << "  " << finalParameters[4] << "  " << finalParameters[5] << std::endl;
		LOG(INFO) << " Number Of Iterations = " << numberOfIterations << " Metric value  = " << bestValue << " Time  = " << (tEnd - tStart)/1000.0 << std::endl;
		LOG(INFO) << "完成第" << i + 1 << "次配准";
		if (CoarseRegistration-1 <= i)
		{
			if (CoarseRegistration - 1 == i)
			{
				std::string tempLumbar = lumbar[CurrentLumbar];
				ConvertedInverse(img1);
				ConvertedInverse(img2);
				WriterImageHomogenization(img1, temppath + "/" + "ImageAP.raw");
				WriterImageHomogenization(img2, temppath + "/" + "ImageOBL.raw");
				Converted(img1);
				Converted(img2);
				
				ImageType::Pointer imagetemp = creatBlankImage(tempsize);
				float* rst1 = new float[lenfix];
				float * buffertemp = imagetemp->GetBufferPointer();
				interpolator1->DRRCudaRun(rst1);
				memcpy(buffertemp, rst1, lenfix * sizeof(float));
				ConvertedInverse(imagetemp);
				WriterImageHomogenization(imagetemp, temppath + "/" + "XrayAP.raw");

				interpolator2->DRRCudaRun(rst1);
				memcpy(buffertemp, rst1, lenfix * sizeof(float));
				ConvertedInverse(imagetemp);
				WriterImageHomogenization(imagetemp, temppath + "/" +"XrayOBL.raw");

				double ptr[4][4];
				PatientourceToRobotMatrix(FinalValue, ds1torobot, ptr);
				std::string fresoult = temppath + "/" + "Result.txt";
				writerTXT(fresoult, ptr);
			}
			else
			{
				std::string tempLumbar = lumbar[CurrentLumbar];
				ConvertedInverse(img1);
				ConvertedInverse(img2);
				WriterImageHomogenization(img1, temppath + "/" + root[tempLumbar + "ImageAP"].asString());
				WriterImageHomogenization(img2, temppath + "/" + root[tempLumbar + "ImageOBL"].asString());
				Converted(img1);
				Converted(img2);

				ImageType::Pointer imagetemp = creatBlankImage(tempsize);
				float* rst1 = new float[lenfix];
				float * buffertemp = imagetemp->GetBufferPointer();
				interpolator1->DRRCudaRun(rst1);
				memcpy(buffertemp, rst1, lenfix * sizeof(float));
				ConvertedInverse(imagetemp);
				WriterImageHomogenization(imagetemp, temppath + "/" + root[tempLumbar + "XrayAP"].asString());

				interpolator2->DRRCudaRun(rst1);
				memcpy(buffertemp, rst1, lenfix * sizeof(float));
				ConvertedInverse(imagetemp);
				WriterImageHomogenization(imagetemp, temppath + "/" + root[tempLumbar + "XrayOBL"].asString());

				double ptr[4][4];
				PatientourceToRobotMatrix(FinalValue, ds1torobot, ptr);
				std::string fresoult = temppath + "/" + root[tempLumbar + "Result"].asString();
				writerTXT(fresoult, ptr);
			}
			
		}

	}
	LOG(INFO) << "完成全部配准";
	return EXIT_SUCCESS;
}

int _main(int argc, char * argv[])
{
	char * file = argv[1];

	Json::Value root;
	if (!ReadJson(file, root))
	{
		std::cout << "未找到配置文件: " << file;
		return false;
	}
	google::InitGoogleLogging("test1");//使用glog之前必须先初始化库，仅需执行一次，括号内为程序名
	FLAGS_alsologtostderr = true;//是否将日志输出到文件和stderr
	FLAGS_colorlogtostderr = true;//是否启用不同颜色显示
	std::string logfile = root["logPath"].asString();
	google::SetLogDestination(google::GLOG_INFO, (logfile + "/INFO_").c_str());//INFO级别的日志都存放到logs目录下且前缀为INFO_
	try
	{
		if (registration(root))
		{
			LOG(ERROR) << "配准失败";
			return false;
		}
	}
	catch (const itk::ExceptionObject & excp)
	{
		LOG(ERROR) << "配准失败";
		return false;
	}

	google::ShutdownGoogleLogging();
	return true;
}
