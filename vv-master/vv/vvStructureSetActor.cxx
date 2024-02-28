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

#include "vvStructureSetActor.h"
#include "vvROIActor.h"

//------------------------------------------------------------------------------
vvStructureSetActor::vvStructureSetActor()
{

}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
vvStructureSetActor::~vvStructureSetActor()
{

}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
int vvStructureSetActor::GetNumberOfROIs() 
{ 
  return mROIActors.size(); 
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
std::vector< QSharedPointer<vvROIActor> > & vvStructureSetActor::GetROIList() 
{ 
  return mROIActors; 
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void vvStructureSetActor::SetStructureSet(clitk::DicomRT_StructureSet * s)
{
  mStructureSet = s;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void vvStructureSetActor::SetSlicerManager(vvSlicerManager * s)
{
  mSlicerManager = s;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
vvROIActor * vvStructureSetActor::GetROIActor(int n)
{
  if (mMapROIIndex.find(n) == mMapROIIndex.end()) {
    std::cerr << "No ROI number " << n << std::endl;
    return NULL;
  }
  return mROIActors[mMapROIIndex[n]].data();
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void vvStructureSetActor::CreateNewROIActor(int n, bool modeBG)
{
  // Check
  clitk::DicomRT_ROI * roi = mStructureSet->GetROIFromROINumber(n);
  if (roi == NULL) {
    std::cerr << "Error. No ROI number " << n << std::endl;
    exit(0);
  }

  // If already exist : delete it
  int old = -1;
  if (mMapROIIndex.find(n) != mMapROIIndex.end())
    old = mMapROIIndex[n];

  // Add ROI Actors
  QSharedPointer<vvROIActor> actor = QSharedPointer<vvROIActor>(new vvROIActor);
  if (old == -1)
    mROIActors.push_back(actor);
  else
    mROIActors[old] = actor;
  actor->SetBGMode(modeBG);
  actor->SetROI(roi);
  actor->SetSlicerManager(mSlicerManager);
  actor->Initialize(n+1); // depth is n+1 to start at 1
  mMapROIIndex[n] = mROIActors.size()-1;
}
//------------------------------------------------------------------------------


