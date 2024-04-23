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

//------------------------------------------------------------------------------
template<class ToolType>
void vvToolCreator<ToolType>::InsertToolInMenu(vvMainWindowBase * m)
{
  mMainWindow = m;

  // Default Initialization
  mToolMenuName = mToolName;
  mToolIconFilename = "noicon";
  mToolTip = mToolName;

  // User Tool Initialization
  ToolType::Initialize();

  // Common Initialization (insertion into menu)
  vvToolCreatorBase::InsertToolInMenu(mMainWindow);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<class ToolType>
vvToolCreator<ToolType>* & vvToolCreator<ToolType>::GetInstance()
{
  if(!mSingleton)
    mSingleton = new vvToolCreator<ToolType>;
  return mSingleton;
}
//------------------------------------------------------------------------------
