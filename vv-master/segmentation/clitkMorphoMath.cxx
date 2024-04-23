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
#ifndef clitkMorphoMath_cxx
#define clitkMorphoMath_cxx

// clitk include
#include "clitkMorphoMath_ggo.h"
#include "clitkMorphoMathGenericFilter.h"

//--------------------------------------------------------------------
int main(int argc, char * argv[]) {

  // Init command line
  GGO(clitkMorphoMath,args_info);
  CLITK_INIT;

  // Creation of a generic filter
  typedef clitk::MorphoMathGenericFilter<args_info_clitkMorphoMath> FilterType;
  FilterType::Pointer filter = FilterType::New();

  // Passing the arguments to the generic filter
  filter->SetArgsInfo(args_info);
 
  CLITK_TRY_CATCH_EXIT(filter->Update());

  return EXIT_SUCCESS;
} // This is the end, my friend
//--------------------------------------------------------------------

#endif
