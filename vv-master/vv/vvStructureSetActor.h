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

#ifndef VVSTRUCTURESETACTOR_H
#define VVSTRUCTURESETACTOR_H

#include "clitkDicomRT_StructureSet.h"
#include <QObject>
#include <QSharedPointer>

class vvSlicerManager;
class vvROIActor;

//------------------------------------------------------------------------------
class vvStructureSetActor: public QObject {
  Q_OBJECT
    public:
  vvStructureSetActor();
  ~vvStructureSetActor();

  void SetStructureSet(clitk::DicomRT_StructureSet * s);
  void SetSlicerManager(vvSlicerManager * s);
  void CreateNewROIActor(int n, bool modeBG=true);
  vvROIActor* GetROIActor(int n);
  int GetNumberOfROIs();
  std::vector< QSharedPointer<vvROIActor> > & GetROIList();
  
 protected:
  clitk::DicomRT_StructureSet * mStructureSet;
  vvSlicerManager * mSlicerManager;
  std::vector< QSharedPointer<vvROIActor> > mROIActors;
  std::map<int, int> mMapROIIndex;

}; // end class vvStructureSetActor
//------------------------------------------------------------------------------

#endif

