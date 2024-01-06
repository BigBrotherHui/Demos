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
#include "itkPointSet.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "itkCenteredTransformInitializer.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkEuler3DTransform.h"
#ifndef itkTwoImageToOneImageMetric_hxx
#define itkTwoImageToOneImageMetric_hxx


#include <time.h>

void matrix_multiply(float* A, float* B)
{
	//int m1 = sizeof(A) / sizeof(A[0]);
	//int n1 = sizeof(A[0]);

	//int m2 = sizeof(B) / sizeof(B[0]);
	//int n2 = sizeof(B[0]);

	//for (int i = 0; i < M1; i++)
	//{
	//	for (int j = 0; j < 3; j++)
	//	{
	//		c[i][j] = 0;
	//		for (int k = 0; k < 2; k++)
	//			c[i][j] += a[i][k] * b[k][j];
	//	}
	//}
}

//void getmatrix(float cs, float ls, float sp, float d, float xs, float tx, float ty, float tz, float Rx, float Ry, float Rz, float* M)
//{
//	float offset = 0.001;
//	cs = cs + offset;
//	ls = ls + offset;
//	//using MatrixType = itk::Matrix<float, 3, 3>;
//	//MatrixType P;
//	//P(0, 0) = cs/d;
//	//P(0, 1) = 1 / sp;
//	//P(0, 2) = 0;
//	//P(1, 0) = ls / d;
//	//P(1, 1) = 0;
//	//P(1, 2) = 1 / sp;
//	//P(2, 0) = 1 / d;
//	//P(2, 1) = 0.0;
//	//P(2, 2) = 0.0;
//
//	float P[3][3];//存储临时的三维边界点
//	for (int i = 0; i < 3; i++) {
//		for (int k = 0; k < 3; k++) {
//			P[i][k] = 0.0;
//		}
//	}
//	P[0][0] = cs / d;
//	P[0][1] = 1 / sp;
//	P[0][2] = 0;
//	P[1][0] = ls / d;
//	P[1][1] = 0;
//	P[1][2] = 1 / sp;
//	P[2][0] = 1 / d;
//	P[2][1] = 0.0;
//	P[2][2] = 0.0;
//
//	using TransformType = itk::Euler3DTransform<double>;
//	TransformType::Pointer    transform = TransformType::New();
//	transform->SetComputeZYX(true);
//	TransformType::OutputVectorType translation;
//
//	translation[0] = tx;
//	translation[1] = ty;
//	translation[2] = tz;
//
//	transform->SetTranslation(translation);
//	const double dtr = (atan(1.0) * 4.0) / 180.0;
//	transform->SetRotation(dtr * Rx, dtr * Ry, dtr * Rz);
//
//	TransformType::MatrixType R = transform->GetMatrix();
//
//
//
//	//using MatrixRotationType = itk::Matrix<float, 3, 4>;
//	//MatrixRotationType T ;
//	//T(0, 3) = -tx;
//	//T(1, 3) = -ty;
//	//T(2, 3) = -tz;
//	//using MatrixMType = itk::Matrix<float, 3, 4>;
//
//	float T[3][4];//存储临时的三维边界点
//
//	for (int i = 0; i < 3; i++) {
//		for (int k = 0; k < 4; k++) {
//			T[i][k] = 0.0;
//		}
//	}
//	T[0][0] = 1;
//	T[1][1] = 1;
//	T[2][2] = 1;
//	T[0][3] = -tx;
//	T[1][3] = -ty;
//	T[2][3] = -tz;
//
//
//	// M = P * R*T;
//	float M1[3][4];
//	for (int i = 0; i < 3; i++)
//	{
//		for (int j = 0; j < 4; j++)
//		{
//			M1[i][j] = 0;
//			for (int k = 0; k < 3; k++)
//				M1[i][j] += R[i][k] * T[k][j];
//		}
//	}
//
//	int ind = 0;
//	for (int i = 0; i < 3; i++)
//	{
//		for (int j = 0; j < 4; j++)
//		{
//			ind = i * 4 + j;
//			M[ind] = 0;
//			for (int k = 0; k < 3; k++)
//				M[ind] += P[i][k] * M1[k][j];
//		}
//	}
//
//
//}

namespace itk
{

template <typename TFixedImage, typename TMovingImage>
TwoImageToOneImageMetric<TFixedImage, TMovingImage>::TwoImageToOneImageMetric()
{
  m_FixedImage1 = nullptr;     // has to be provided by the user.
  m_FixedImage2 = nullptr;     // has to be provided by the user.
  m_MovingImage = nullptr;     // has to be provided by the user.
  m_Transform = nullptr;       // has to be provided by the user.
  m_Interpolator1 = nullptr;   // has to be provided by the user.
  m_Interpolator2 = nullptr;   // has to be provided by the user.
  m_GradientImage = nullptr;   // will receive the output of the filter;
  m_ComputeGradient = true;    // metric computes gradient by default
  m_NumberOfPixelsCounted = 0; // initialize to zero
  m_GradientImage = nullptr;   // computed at initialization
}


/*
 * Set the parameters that define a unique transform
 */
template <typename TFixedImage, typename TMovingImage>
void
TwoImageToOneImageMetric<TFixedImage, TMovingImage>::SetTransformParameters(const ParametersType & parameters) const
{
  if (!m_Transform)
  {
    itkExceptionMacro(<< "Transform has not been assigned");
  }
  m_Transform->SetParameters(parameters);
}


template <typename TFixedImage, typename TMovingImage>
void
TwoImageToOneImageMetric<TFixedImage, TMovingImage>::Initialize()
{

  if (!m_Transform)
  {
    itkExceptionMacro(<< "Transform is not present");
  }

  if (!m_Interpolator1)
  {
    itkExceptionMacro(<< "Interpolator1 is not present");
  }
  if (!m_Interpolator2)
  {
    itkExceptionMacro(<< "Interpolator2 is not present");
  }

  if (!m_MovingImage)
  {
    itkExceptionMacro(<< "MovingImage is not present");
  }

  if (!m_FixedImage1)
  {
    itkExceptionMacro(<< "FixedImage1 is not present");
  }

  if (!m_FixedImage2)
  {
    itkExceptionMacro(<< "FixedImage2 is not present");
  }

  if (m_FixedImageRegion1.GetNumberOfPixels() == 0)
  {
    itkExceptionMacro(<< "FixedImageRegion1 is empty");
  }

  if (m_FixedImageRegion2.GetNumberOfPixels() == 0)
  {
    itkExceptionMacro(<< "FixedImageRegion2 is empty");
  }

  // If the image is provided by a source, update the source.
  if (m_MovingImage->GetSource())
  {
    m_MovingImage->GetSource()->Update();
  }

  //If the image is provided by a source, update the source.
  if (m_FixedImage1->GetSource())
  {
    m_FixedImage1->GetSource()->Update();
  }

  if (m_FixedImage2->GetSource())
  {
    m_FixedImage2->GetSource()->Update();
  }

  // Make sure the FixedImageRegion is within the FixedImage buffered region
  if (!m_FixedImageRegion1.Crop(m_FixedImage1->GetBufferedRegion()))
  {
    itkExceptionMacro(<< "FixedImageRegion1 does not overlap the fixed image buffered region");
  }

  if (!m_FixedImageRegion2.Crop(m_FixedImage2->GetBufferedRegion()))
  {
    itkExceptionMacro(<< "FixedImageRegion2 does not overlap the fixed image buffered region");
  }

  m_Interpolator1->SetInputImage(m_MovingImage);
  m_Interpolator2->SetInputImage(m_MovingImage);

  ImageType::SizeType imgSize1 = m_FixedImage1->GetLargestPossibleRegion().GetSize();
  ImageType::SpacingType spacing1 = m_FixedImage1->GetSpacing();
  int img2dPn[2];
  img2dPn[0] = imgSize1[0];
  img2dPn[1] = imgSize1[1];

  float img2dPs[2];
  img2dPs[0] = spacing1[0];
  img2dPs[1] = spacing1[1];

  m_Interpolator1->PrepareMemory4Cuda();
  m_Interpolator2->PrepareMemory4Cuda();

  if (m_ComputeGradient)
  {

    GradientImageFilterPointer gradientFilter = GradientImageFilterType::New();

    gradientFilter->SetInput(m_MovingImage);

    const typename MovingImageType::SpacingType & spacing = m_MovingImage->GetSpacing();
    double                                        maximumSpacing = 0.0;
    for (unsigned int i = 0; i < MovingImageDimension; i++)
    {
      if (spacing[i] > maximumSpacing)
      {
        maximumSpacing = spacing[i];
      }
    }
    gradientFilter->SetSigma(maximumSpacing);
    gradientFilter->SetNormalizeAcrossScale(true);

    gradientFilter->Update();

    m_GradientImage = gradientFilter->GetOutput();
  }

  // If there are any observers on the metric, call them to give the
  // user code a chance to set parameters on the metric
  this->InvokeEvent(InitializeEvent());
}


template <typename TFixedImage, typename TMovingImage>
void
TwoImageToOneImageMetric<TFixedImage, TMovingImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "ComputeGradient: " << static_cast<typename NumericTraits<bool>::PrintType>(m_ComputeGradient)
     << std::endl;
  os << indent << "Moving Image: " << m_MovingImage.GetPointer() << std::endl;
  os << indent << "Fixed  Image 1: " << m_FixedImage1.GetPointer() << std::endl;
  os << indent << "Fixed  Image 2: " << m_FixedImage2.GetPointer() << std::endl;
  os << indent << "Gradient Image: " << m_GradientImage.GetPointer() << std::endl;
  os << indent << "Transform:    " << m_Transform.GetPointer() << std::endl;
  os << indent << "Interpolator 1: " << m_Interpolator1.GetPointer() << std::endl;
  os << indent << "Interpolator 2: " << m_Interpolator2.GetPointer() << std::endl;
  os << indent << "FixedImageRegion 1: " << m_FixedImageRegion1 << std::endl;
  os << indent << "FixedImageRegion 2: " << m_FixedImageRegion2 << std::endl;
  os << indent << "Moving Image Mask: " << m_MovingImageMask.GetPointer() << std::endl;
  os << indent << "Fixed Image Mask 1: " << m_FixedImageMask1.GetPointer() << std::endl;
  os << indent << "Fixed Image Mask 2: " << m_FixedImageMask2.GetPointer() << std::endl;
  os << indent << "Number of Pixels Counted: " << m_NumberOfPixelsCounted << std::endl;
}


} // end namespace itk

#endif
