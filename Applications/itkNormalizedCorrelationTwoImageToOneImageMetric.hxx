/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkNormalizedCorrelationTwoImageToOneImageMetric_hxx
#define itkNormalizedCorrelationTwoImageToOneImageMetric_hxx

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageDuplicator.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkJoinImageFilter.h"
#include "itkImageToHistogramFilter.h"
#include <time.h>
#include <string>

void savedata(const float* rst1, int len, std::string  st)
{
	FILE *fpwrt = NULL;
	const char* file_c = st.c_str();
	fopen_s(&fpwrt, file_c, "wb+");
	if (fpwrt == NULL)
	{
		std::cout << "写入文件创建错误" << std::endl;
	}
	fwrite(rst1, sizeof(float), len, fpwrt);
	fclose(fpwrt);
}



void readdata(std::string st, int len, float* rst1)
{
	FILE *fp = NULL;
	const char* file_c = st.c_str();
	if ((fp = fopen(file_c, "rb")) == NULL)
	{
		printf("can not open the raw image");
	}
	else
	{
		printf("read OK\n");
	}
	fread(rst1, sizeof(float), len, fp);
}


double CalcCnn(const float* ptr, float* rst1, const float* FixMask ,int len, bool m_SubtractMean)
{
	double fixedValue;
	double movingValue;
	int NumberOfPixelsCounted = 0;
	double sff(0.0), smm(0.0), sfm(0.0), sf(0.0), sm(0.0);
	for (int i = 0; i < len; ++i)
	{

		fixedValue = ptr[i];
		movingValue = rst1[i];

		if (FixMask[i] < 1&& movingValue>0.0001)
		{
			sff += fixedValue * fixedValue;
			smm += movingValue * movingValue;
			sfm += fixedValue * movingValue;
			if (m_SubtractMean)
			{
				sf += fixedValue;
				sm += movingValue;
			}
			NumberOfPixelsCounted++;
		}
		//NumberOfPixelsCounted++;

	}

	if (m_SubtractMean && NumberOfPixelsCounted > 0)
	{
		sff -= (sf * sf / NumberOfPixelsCounted);
		smm -= (sm * sm / NumberOfPixelsCounted);
		sfm -= (sf * sm / NumberOfPixelsCounted);
	}

	double denom1 = -1.0 * sqrt(sff * smm);
	double measure1(0.0);
	if (NumberOfPixelsCounted > 0 && denom1 != 0.0)
	{
		measure1 = sfm / denom1;
	}
	else
	{
		measure1 = 0.0;
	}
	return measure1;
}


namespace itk
{

	template <typename TFixedImage, typename TMovingImage>
	NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage,
		TMovingImage>::NormalizedCorrelationTwoImageToOneImageMetric()
	{
		m_SubtractMean = false;
		m_bInterpolator1 = true;
		m_bInterpolator2 = true;
		ii = 0;
	}

	template <typename TFixedImage, typename TMovingImage>
	NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage,
		TMovingImage>::~NormalizedCorrelationTwoImageToOneImageMetric()
	{
		std::cout << "执行次数 " << ii << std::endl;
	}


	template <typename TFixedImage, typename TMovingImage>
	bool
		NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage, TMovingImage>::CalculateMI(FixedImageConstPointer image1,
			FixedImageConstPointer image2,
			int      HistogramBins,
			double   MarginalScale,
			double & Entropy) const
	{
		using JoinFilterType = itk::JoinImageFilter<FixedImageType, FixedImageType>;

		auto joinFilter = JoinFilterType::New();

		joinFilter->SetInput1(image1);
		joinFilter->SetInput2(image2);

		try
		{
			joinFilter->Update();
		}
		catch (const itk::ExceptionObject & excp)
		{
			std::cerr << excp << std::endl;
			return EXIT_FAILURE;
		}
		using VectorImageType = JoinFilterType::OutputImageType;

		using HistogramFilterType = itk::Statistics::ImageToHistogramFilter<VectorImageType>;

		auto histogramFilter = HistogramFilterType::New();

		histogramFilter->SetInput(joinFilter->GetOutput());

		histogramFilter->SetMarginalScale(MarginalScale);

		using HistogramSizeType = HistogramFilterType::HistogramSizeType;

		HistogramSizeType size(2);

		size[0] = 255; // number of bins for the first  channel
		size[1] = 255; // number of bins for the second channel

		histogramFilter->SetHistogramSize(size);
		// Software Guide : EndCodeSnippet
		using HistogramMeasurementVectorType = HistogramFilterType::HistogramMeasurementVectorType;

		HistogramMeasurementVectorType binMinimum(3);
		HistogramMeasurementVectorType binMaximum(3);

		binMinimum[0] = -0.5;
		binMinimum[1] = -0.5;
		binMinimum[2] = -0.5;

		binMaximum[0] = 255.5;
		binMaximum[1] = 255.5;
		binMaximum[2] = 255.5;

		histogramFilter->SetHistogramBinMinimum(binMinimum);
		histogramFilter->SetHistogramBinMaximum(binMaximum);

		histogramFilter->Update();
		using HistogramType = HistogramFilterType::HistogramType;
		const HistogramType * histogram = histogramFilter->GetOutput();

		HistogramType::ConstIterator itr = histogram->Begin();
		HistogramType::ConstIterator end = histogram->End();

		const double Sum = histogram->GetTotalFrequency();
		double JointEntropy = 0.0;

		while (itr != end)
		{
			const double count = itr.GetFrequency();
			if (count > 0.0)
			{
				const double probability = count / Sum;
				JointEntropy += -probability * std::log(probability) / std::log(2.0);
			}
			++itr;
		}

		size[0] = 255; // number of bins for the first  channel
		size[1] = 0;   // number of bins for the second channel

		histogramFilter->SetHistogramSize(size);
		histogramFilter->Update();

		itr = histogram->Begin();
		end = histogram->End();

		double Entropy1 = 0.0;

		while (itr != end)
		{
			const double count = itr.GetFrequency();
			if (count > 0.0)
			{
				const double probability = count / Sum;
				Entropy1 += -probability * std::log(probability) / std::log(2.0);
			}
			++itr;
		}

		size[0] = 0;   // number of bins for the first channel
		size[1] = 255; // number of bins for the second channel

		histogramFilter->SetHistogramSize(size);
		histogramFilter->Update();
		itr = histogram->Begin();
		end = histogram->End();

		double Entropy2 = 0.0;

		while (itr != end)
		{
			const double count = itr.GetFrequency();
			if (count > 0.0)
			{
				const double probability = count / Sum;
				Entropy2 += -probability * std::log(probability) / std::log(2.0);
			}
			++itr;
		}

		Entropy = Entropy1 + Entropy2 - JointEntropy;
		return true;
	}



	template <typename TFixedImage, typename TMovingImage>
	typename NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage, TMovingImage>::MeasureType
		NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage, TMovingImage>::GetValue(
			const TransformParametersType & parameters) const
	{
		FixedImageConstPointer fixedImage1 = this->m_FixedImage1;

		if (!fixedImage1)
		{
			itkExceptionMacro(<< "Fixed image1 has not been assigned");
		}

		FixedImageConstPointer fixedImage2 = this->m_FixedImage2;

		if (!fixedImage2)
		{
			itkExceptionMacro(<< "Fixed image2 has not been assigned");
		}

		using FixedIteratorType = itk::ImageRegionConstIteratorWithIndex<FixedImageType>;

		// Calculate the measure value between fixed image 1 and the moving image

		MeasureType measure1(0.0);
		MeasureType measure2(0.0);
		MeasureType measuresum(0.0);

		this->m_NumberOfPixelsCounted = 0;
		this->SetTransformParameters(parameters);

		const double dtr = (atan(1.0) * 4.0) / 180.0;
		//std::cout << parameters[0] / dtr << " " << parameters[1] / dtr << " " << parameters[2] / dtr << " "
		//    << parameters[3] << " " << parameters[4] << " " << parameters[5] << " " << std::endl;

		//超出范围为0
		double AngleRange = 30 * dtr;
		double shiftRange = 200;
		for (int i = 0; i < 6; ++i)
		{
			double  a = parameters[i];
			if (i < 3)
			{

				if ((parameters[i] > (m_InitialTransform[i] + AngleRange)) || (parameters[i] < (m_InitialTransform[i] - AngleRange)))
				{
					return 0;
				}
			}
			else
			{
				if ((parameters[i] > (m_InitialTransform[i] + shiftRange)) || (parameters[i] < (m_InitialTransform[i] - shiftRange)))
				{
					return 0;
				}
			}
		}

		ImageType::SizeType imgSize = m_FixedImage1->GetLargestPossibleRegion().GetSize();
		int len = imgSize[0] * imgSize[1];
		float* rst1 = new float[len];
		float* rst2 = new float[len];
		const float * buffer1 = m_FixedImage1->GetBufferPointer();
		const float * buffer2 = m_FixedImage2->GetBufferPointer();

		const float * fixmask1 = m_FixedMask1->GetBufferPointer();
		const float * fixmask2 = m_FixedMask2->GetBufferPointer();

		std::clock_t  t1, t2;
		t1 = clock();
		if (true == m_bInterpolator1)
		{
			this->m_Interpolator1->DRRCudaRun(rst1);
			measure1 = CalcCnn(buffer1, rst1, fixmask1,len, m_SubtractMean);
		}

		if (true == m_bInterpolator2)
		{
			this->m_Interpolator2->DRRCudaRun(rst2);
			measure2 = CalcCnn(buffer2, rst2, fixmask2,len, m_SubtractMean);
		}
		//t2 = clock();
		//std::cout << t2 - t1 << std::endl;
		//std::cout << measure1 <<"   "<< measure2 << std::endl;

		//savedata(buffer1, len, "fix1.raw");
		//savedata(buffer2, len, "fix2.raw");
		//savedata(rst1, len, "drr1.raw");
		//savedata(rst2, len, "drr2.raw");

		//savedata(fixmask1, len, "fixmask1.raw");
		//savedata(fixmask2, len, "fixmask2.raw");

		delete rst1;
		delete rst2;
		measuresum = -(measure1 + measure2) / 2.0;
		//std::cout << measuresum << std::endl;
		//ii =ii+1;
		++ii;
		return measuresum;
	}


	template <typename TFixedImage, typename TMovingImage>
	void
		NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage, TMovingImage>::GetDerivative(
			const TransformParametersType & itkNotUsed(parameters),
			DerivativeType &                itkNotUsed(derivative)) const
	{
		// under construction
	}


	template <typename TFixedImage, typename TMovingImage>
	void
		NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage, TMovingImage>::GetValueAndDerivative(
			const TransformParametersType & itkNotUsed(parameters),
			MeasureType &                   itkNotUsed(value),
			DerivativeType &                itkNotUsed(derivative)) const
	{
		// under construction
	}


	template <typename TFixedImage, typename TMovingImage>
	void
		NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage, TMovingImage>::PrintSelf(std::ostream & os,
			Indent         indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "SubtractMean: " << m_SubtractMean << std::endl;
	}

} // end namespace itk


#endif
