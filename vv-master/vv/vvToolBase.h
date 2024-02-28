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

#ifndef VVTOOLBASE_H
#define VVTOOLBASE_H

#include "vvToolBaseBase.h"
#include "vvToolCreator.h"
#include <QIcon>

//------------------------------------------------------------------------------
template<class ToolType>
class vvToolBase : public vvToolBaseBase {
public:
  vvToolBase(vvMainWindowBase * m);
  virtual ~vvToolBase() {};
  
  static void Initialize();  // can't be virtual, must be overwritten

  static void SetToolName(QString n) { vvToolCreator<ToolType>::GetInstance()->mToolName = n; }
  static void SetToolMenuName(QString n) { vvToolCreator<ToolType>::GetInstance()->mToolMenuName = n; }
  static void SetToolIconFilename(QString n) { vvToolCreator<ToolType>::GetInstance()->mToolIconFilename = n; }
  static void SetToolTip(QString n) { vvToolCreator<ToolType>::GetInstance()->mToolTip = n; }
  static void SetToolExperimental(bool exp) { vvToolCreator<ToolType>::GetInstance()->mExperimental = exp; }
  static void SetToolInMenu(std::string menu) { vvToolCreator<ToolType>::GetInstance()->SetMenuName(menu); }

  static QIcon GetToolIcon() { return QIcon(vvToolCreator<ToolType>::GetInstance()->mToolIconFilename); }
  static QString & GetToolName() { return vvToolCreator<ToolType>::GetInstance()->mToolName; }

  vvSlicerManager * AddImage(vvImage::Pointer image,std::string filename) {
    return CREATOR(ToolType)->GetMainWindow()->AddImage(image,filename); 
  }

  void setSender(QObject *s) { mSender = s; }

 protected:
  QObject * mSender;
};
//------------------------------------------------------------------------------

#include "vvToolBase.txx"

#endif

