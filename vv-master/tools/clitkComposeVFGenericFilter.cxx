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
#ifndef CLITKCOMPOSEVFGENERICFILTER_CXX
#define CLITKCOMPOSEVFGENERICFILTER_CXX
#include "clitkComposeVFGenericFilter.h"


namespace clitk {

  clitk::ComposeVFGenericFilter::ComposeVFGenericFilter()
  {
    m_Verbose=false;
    m_Type = 0;
  }


  void clitk::ComposeVFGenericFilter::Update()
  {
    //Get the image Dimension and PixelType
    int Dimension;
    std::string PixelType;

    clitk::ReadImageDimensionAndPixelType(m_InputName1, Dimension, PixelType);

    if(Dimension==2) UpdateWithDim<2>(PixelType);
    else if(Dimension==3) UpdateWithDim<3>(PixelType);
    else 
      {
	std::cout<<"Error, Only for 2 and 3 Dimensions!!!"<<std::endl ;
	return;
      }

  }
} //end namespace

#endif
