/*=========================================================================
  Program:   vv                     http://www.creatis.insa-lyon.fr/rio/vv

  Authors belong to:
  - University of LYON              http://www.universite-lyon.fr/
  - Léon Bérard cancer center       http://www.centreleonberard.fr
  - CREATIS CNRS laboratory         http://www.creatis.insa-lyon.fr

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the copyright notices for more information.

  It is distributed under dual licence

  - BSD        See included LICENSE.txt file
  - CeCILL-B   http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
===========================================================================**/
#ifndef clitkSetBackgroundGenericFilter_txx
#define clitkSetBackgroundGenericFilter_txx

namespace clitk
{

//-------------------------------------------------------------------
// Update with the number of dimensions
//-------------------------------------------------------------------
template<unsigned int Dimension>
void
SetBackgroundGenericFilter::UpdateWithDim(std::string PixelType)
{
  if (m_Verbose) std::cout << "Image was detected to be "<<Dimension<<"D and "<< PixelType<<"..."<<std::endl;

  if(PixelType == "short") {
    if (m_Verbose) std::cout << "Launching filter in "<< Dimension <<"D and signed short..." << std::endl;
    UpdateWithDimAndPixelType<Dimension, signed short>();
  }
  //    else if(PixelType == "unsigned_short"){
  //       if (m_Verbose) std::cout  << "Launching filter in "<< Dimension <<"D and unsigned_short..." << std::endl;
  //       UpdateWithDimAndPixelType<Dimension, unsigned short>();
  //     }

  else if (PixelType == "unsigned_char") {
    if (m_Verbose) std::cout  << "Launching filter in "<< Dimension <<"D and unsigned_char..." << std::endl;
    UpdateWithDimAndPixelType<Dimension, unsigned char>();
  }

  //     else if (PixelType == "char"){
  //       if (m_Verbose) std::cout  << "Launching filter in "<< Dimension <<"D and signed_char..." << std::endl;
  //       UpdateWithDimAndPixelType<Dimension, signed char>();
  //     }
  else {
    if (m_Verbose) std::cout  << "Launching filter in "<< Dimension <<"D and float..." << std::endl;
    UpdateWithDimAndPixelType<Dimension, float>();
  }
}

//-------------------------------------------------------------------
// Update with the number of dimensions and the pixeltype
//-------------------------------------------------------------------
template <unsigned int Dimension, class  PixelType>
void
SetBackgroundGenericFilter::UpdateWithDimAndPixelType()
{

  // ImageTypes
  typedef itk::Image<PixelType, Dimension> InputImageType;
  typedef itk::Image<unsigned char, Dimension> MaskImageType;

  // Read the input
  typedef itk::ImageFileReader<InputImageType> InputReaderType;
  typename InputReaderType::Pointer reader = InputReaderType::New();
  reader->SetFileName( m_InputFileName);
  typename InputImageType::Pointer input= reader->GetOutput();

  // Read the mask
  typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
  typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
  maskReader->SetFileName( m_ArgsInfo.mask_arg);
  typename MaskImageType::Pointer mask= maskReader->GetOutput();

  // Filter setting background
  typedef clitk::SetBackgroundImageFilter<InputImageType,MaskImageType, InputImageType> SetBackgroundFilterType;
  typename SetBackgroundFilterType::Pointer setBackgroundFilter = SetBackgroundFilterType::New();
  setBackgroundFilter->SetInput(input);
  setBackgroundFilter->SetInput2(mask);
  if(m_ArgsInfo.fg_flag)  setBackgroundFilter->SetForeground(m_ArgsInfo.fg_flag);
  if(m_ArgsInfo.maskValue_given)  setBackgroundFilter->SetMaskValue(m_ArgsInfo.maskValue_arg);
  if(m_ArgsInfo.outsideValue_given)  setBackgroundFilter->SetOutsideValue(m_ArgsInfo.outsideValue_arg);
  setBackgroundFilter->Update();
  typename InputImageType::Pointer output =setBackgroundFilter->GetOutput();

  // Output
  typedef itk::ImageFileWriter<InputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(m_ArgsInfo.output_arg);
  writer->SetInput(output);
  writer->Update();

}


}//end clitk

#endif //#define clitkSetBackgroundGenericFilter_txx
