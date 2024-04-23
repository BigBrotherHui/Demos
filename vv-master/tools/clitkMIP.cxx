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

// clitk
#include "clitkMIP_ggo.h"
#include "clitkMIPGenericFilter.h"

//--------------------------------------------------------------------
int main(int argc, char * argv[])
{

  // Init command line
  GGO(clitkMIP, args_info);
  CLITK_INIT;

  // Filter
  typedef clitk::MIPGenericFilter FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter->SetArgsInfo(args_info);

  try {
    filter->Update();
  } catch(std::runtime_error e) {
    std::cout << e.what() << std::endl;
  }

  return EXIT_SUCCESS;
} // This is the end, my friend
//--------------------------------------------------------------------
