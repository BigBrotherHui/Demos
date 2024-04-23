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

#include "vvToolManager.h"
#include "vvToolCreatorBase.h"
#include "vvMainWindowBase.h"
#include <QAction>
//------------------------------------------------------------------------------
/// Unique static instance
vvToolManager* vvToolManager::mSingleton=0;
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
vvToolManager * vvToolManager::GetInstance()
{
  if (mSingleton == 0) {
    mSingleton = new vvToolManager;
  }
  return mSingleton;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void vvToolManager::AddTool(vvToolCreatorBase * v)
{
  //std::cout << "Adding the tool <" << v->mToolName.toStdString() << ">." << std::endl;
  GetInstance()->mListOfTools.push_back(v);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void vvToolManager::InsertToolsInMenu(vvMainWindowBase * m)
{
  for(unsigned int i=0; i<GetInstance()->mListOfTools.size(); i++) {
    GetInstance()->mListOfTools[i]->InsertToolInMenu(m);
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void vvToolManager::EnableToolsInMenu(vvMainWindowBase * /* m */, bool enable){
  std::vector<vvToolCreatorBase *>::iterator it;
  for(it=GetInstance()->mListOfTools.begin(); it!=GetInstance()->mListOfTools.end(); ++it){
    if((*it)->mAction){
      (*it)->mAction->setEnabled(enable);
    }
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
vvToolCreatorBase * vvToolManager::GetToolCreatorFromName(QString toolTypeName)
{
  std::vector<vvToolCreatorBase *> & v = vvToolManager::GetInstance()->GetListOfTools();
  int index=-1;
  for(uint i=0; i<v.size(); i++) {
    // DD(v[i]->mToolName.toStdString());
    if (v[i]->mToolName == toolTypeName) {
      index = i;
    }
  }
  if (index == -1) {
    std::cerr << "Error, ToolCreator named '" << toolTypeName.toStdString() 
              << "' does not exist. Abort." << std::endl;
    return NULL;
  }
  return v[index];
}
//------------------------------------------------------------------------------


