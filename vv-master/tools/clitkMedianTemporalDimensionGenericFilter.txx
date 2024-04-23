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
#ifndef clitkMedianTemporalDimensionGenericFilter_txx
#define clitkMedianTemporalDimensionGenericFilter_txx

/* =================================================
 * @file   clitkMedianTemporalDimensionGenericFilter.txx
 * @author
 * @date
 *
 * @brief
 *
 ===================================================*/

#include <algorithm>

namespace clitk
{

  //-----------------------------------------------------------
  // Constructor
  //-----------------------------------------------------------
  template<class args_info_type>
    MedianTemporalDimensionGenericFilter<args_info_type>::MedianTemporalDimensionGenericFilter()
    {
      m_Verbose=false;
      m_InputFileName="";
    }


  //-----------------------------------------------------------
  // Update
  //-----------------------------------------------------------
  template<class args_info_type>
    void MedianTemporalDimensionGenericFilter<args_info_type>::Update()
    {
      // Read the Dimension and PixelType
      int Dimension, Components;
      std::string PixelType;
      ReadImageDimensionAndPixelType(m_InputFileName, Dimension, PixelType, Components);


      // Call UpdateWithDim
      if (m_ArgsInfo.input_given>1) Dimension+=1;
      if(Dimension==2) UpdateWithDim<2>(PixelType, Components);
      else if(Dimension==3)UpdateWithDim<3>(PixelType, Components);
      else if (Dimension==4)UpdateWithDim<4>(PixelType, Components);
      else {
        std::cout<<"Error, Only for 2, 3 or 4D!!!"<<std::endl ;
        return;
      }
    }

  //-------------------------------------------------------------------
  // Update with the number of dimensions
  //-------------------------------------------------------------------
  template<class args_info_type>
    template<unsigned int Dimension>
    void
    MedianTemporalDimensionGenericFilter<args_info_type>::UpdateWithDim(const std::string PixelType, const int Components)
    {
      if (m_Verbose) std::cout << "Image was detected to be "<<Dimension<<"D with "<<Components<<" component(s) of "<<  PixelType<<"..."<<std::endl;

      if (Components==1) {
        if(PixelType == "short") {
          if (m_Verbose) std::cout << "Launching filter in "<< Dimension <<"D and signed short..." << std::endl;
          UpdateWithDimAndPixelType<Dimension, signed short>();
        }
        else if (PixelType == "unsigned_char") {
          if (m_Verbose) std::cout  << "Launching filter in "<< Dimension <<"D and unsigned_char..." << std::endl;
          UpdateWithDimAndPixelType<Dimension, unsigned char>();
        }
        else if (PixelType == "float") {
          if (m_Verbose) std::cout  << "Launching filter in "<< Dimension <<"D and float..." << std::endl;
          UpdateWithDimAndPixelType<Dimension, float>();
        }
        else if (PixelType == "double") {
          if (m_Verbose) std::cout  << "Launching filter in "<< Dimension <<"D and double..." << std::endl;
          UpdateWithDimAndPixelType<Dimension, double>();
        }
        else std::cerr<<"Number of components is "<<Components<<", not supported!"<<std::endl;
      }
    }


  //-------------------------------------------------------------------
  // Update with the number of dimensions and the pixeltype
  //-------------------------------------------------------------------
  template<class args_info_type>
    template <unsigned int Dimension, class  PixelType>
    void
    MedianTemporalDimensionGenericFilter<args_info_type>::UpdateWithDimAndPixelType()
    {

      // ImageTypes
      typedef itk::Image<PixelType, Dimension> InputImageType;
      typedef itk::Image<PixelType, Dimension-1> OutputImageType;
      typename InputImageType::Pointer input;

      if (m_ArgsInfo.input_given ==1 ) {
        // Read the input
        typedef itk::ImageFileReader<InputImageType> InputReaderType;
        typename InputReaderType::Pointer reader = InputReaderType::New();
        reader->SetFileName( m_InputFileName);
        reader->Update();
        input= reader->GetOutput();
      }

      else {
        // Read and join multiple inputs
        if (m_Verbose) std::cout<<m_ArgsInfo.input_given<<" inputs given..."<<std::endl;
        std::vector<std::string> filenames;
        for(unsigned int i=0; i<m_ArgsInfo.input_given; i++) {
          if (m_Verbose) std::cout<<m_ArgsInfo.input_arg[i]<<std::endl;
          filenames.push_back(m_ArgsInfo.input_arg[i]);
        }

        typedef itk::ImageSeriesReader<InputImageType> ImageReaderType;
        typename  ImageReaderType::Pointer reader= ImageReaderType::New();
        reader->SetFileNames(filenames);
        reader->Update();
        input =reader->GetOutput();
      }


      // Output properties
      typename OutputImageType::RegionType region;
      typename OutputImageType::RegionType::SizeType size;
      typename OutputImageType::IndexType index;
      typename OutputImageType::SpacingType spacing;
      typename OutputImageType::PointType origin;
      typename InputImageType::RegionType region4D=input->GetLargestPossibleRegion();
      typename InputImageType::RegionType::SizeType size4D=region4D.GetSize();
      typename InputImageType::IndexType index4D=region4D.GetIndex();
      typename InputImageType::SpacingType spacing4D=input->GetSpacing();
      typename InputImageType::PointType origin4D=input->GetOrigin();

      for (unsigned int i=0; i< Dimension-1; i++) {
        size[i]=size4D[i];
        index[i]=index4D[i];
        spacing[i]=spacing4D[i];
        origin[i]=origin4D[i];
      }
      region.SetSize(size);
      region.SetIndex(index);
      typename OutputImageType::Pointer output= OutputImageType::New();
      output->SetRegions(region);
      output->SetSpacing(spacing);
      output->SetOrigin(origin);
      output->Allocate();


      // Region iterators
      typedef itk::ImageRegionIterator<InputImageType> IteratorType;
      std::vector<IteratorType> iterators(size4D[Dimension-1]);
      for (unsigned int i=0; i< size4D[Dimension-1]; i++) {
        typename InputImageType::RegionType regionIt=region4D;
        typename InputImageType::RegionType::SizeType sizeIt=regionIt.GetSize();
        sizeIt[Dimension-1]=1;
        regionIt.SetSize(sizeIt);
        typename InputImageType::IndexType indexIt=regionIt.GetIndex();
        indexIt[Dimension-1]=i;
        regionIt.SetIndex(indexIt);
        iterators[i]=IteratorType(input, regionIt);
      }

      typedef itk::ImageRegionIterator<OutputImageType> OutputIteratorType;
      OutputIteratorType avIt(output, output->GetLargestPossibleRegion());

      // Take the median
      double value;
      typename InputImageType::IndexValueType temporal_dimension = size4D[Dimension-1];
      std::vector<PixelType> temp(temporal_dimension);
      while (!(iterators[0]).IsAtEnd()) {
        value=0.;
        for (unsigned int i=0; i<temporal_dimension; i++) {
          temp[i] = iterators[i].Get();
          ++(iterators[i]);
        }
        if (temporal_dimension & 1) {
          nth_element(temp.begin(), temp.begin() + temporal_dimension/2,temp.end());
          value = temp[temporal_dimension/2];
        } else {
          nth_element(temp.begin(), temp.begin() + temporal_dimension/2 - 1, temp.end());
          value = temp[temporal_dimension/2 - 1];
          nth_element(temp.begin(), temp.begin() + temporal_dimension/2, temp.end());
          value += temp[temporal_dimension/2];
          value /= 2;
        }
        avIt.Set(value);
        ++avIt;
      }

  // Output
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(m_ArgsInfo.output_arg);
  writer->SetInput(output);
  writer->Update();

}


}//end clitk

#endif //#define clitkMedianTemporalDimensionGenericFilter_txx
