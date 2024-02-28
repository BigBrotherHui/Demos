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
===========================================================================*/
#ifndef clitkExtractPatientGenericFilter_cxx
#define clitkExtractPatientGenericFilter_cxx

/* =================================================
 * @file   clitkExtractPatientGenericFilter.cxx
 * @author 
 * @date   
 * 
 * @brief 
 * 
 ===================================================*/

#include "clitkExtractPatientGenericFilter.h"


namespace clitk
{


  //-----------------------------------------------------------
  // Constructor
  //-----------------------------------------------------------
  ExtractPatientGenericFilter:::ExtractPatientGenericFilter()
  {
    m_Verbose=false;
    m_InputFileName="";
  }


  //-----------------------------------------------------------
  // Update
  //-----------------------------------------------------------
  void ExtractPatientGenericFilter::Update()
  {
    // Read the Dimension and PixelType
    int Dimension;
    std::string PixelType;
    ReadImageDimensionAndPixelType(m_InputFileName, Dimension, PixelType);

    
    // Call UpdateWithDim
    if(Dimension==2) UpdateWithDim<2>(PixelType);
    else if(Dimension==3) UpdateWithDim<3>(PixelType);
    // else if (Dimension==4)UpdateWithDim<4>(PixelType); 
    else 
      {
	std::cout<<"Error, Only for 2 or 3  Dimensions!!!"<<std::endl ;
	return;
      }
  }


} //end clitk

#endif  //#define clitkExtractPatientGenericFilter_cxx
