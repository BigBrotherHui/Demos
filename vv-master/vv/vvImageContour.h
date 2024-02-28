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
#ifndef VVIMAGECONTOUR_H
#define VVIMAGECONTOUR_H

#include "clitkCommon.h"
#include "vvSlicer.h"

class vtkImageClip;
class vtkMarchingSquares;
class vtkActor;
class vvImage;

//------------------------------------------------------------------------------
class vvImageContour: public itk::LightObject
{
  //  Q_OBJECT
public:
  typedef vvImageContour Self;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro(Self);

  void SetSlicer(vvSlicer * slicer);
  void Update(double value);
  void HideActors();
  void ShowActors();
  void SetColor(double r, double g, double b);
  void SetLineWidth(double w);
  void SetImage(vvImage::Pointer image);
  void SetPreserveMemoryModeEnabled(bool b);
  void SetDepth(double d);
  void RemoveActors();

protected:
  vvSlicer * mSlicer;
  int mSlice;
  int mTSlice;
  double mValue;
  int mPreviousTSlice;
  double mPreviousValue;
  bool mHiddenImageIsUsed;
  vvImage::Pointer mHiddenImage;
  bool mDisplayModeIsPreserveMemory;
  double mDepth;

  // For preserveMemory mode
  std::vector<vtkSmartPointer<vtkActor> > mSquaresActorList;
  std::vector<vtkSmartPointer<vtkImageClip> > mClipperList;
  std::vector<vtkSmartPointer<vtkMarchingSquares> > mSquaresList;
  std::vector<vtkSmartPointer<vtkPolyDataMapper> > mSquaresMapperList;

  // For fast cache mode
  int mPreviousSlice;
  int mPreviousOrientation;
  std::vector<std::vector<vtkActor*> > mListOfCachedContourActors;

  // Functions
  void InitializeCacheMode();
  void UpdateWithPreserveMemoryMode();
  void UpdateWithFastCacheMode();
  void CreateNewActor(int numImage);
  void UpdateActor(vtkActor * actor, vtkPolyDataMapper * mapper, vtkMarchingSquares * squares, vtkImageClip * clipper,
                   double threshold, int orientation, int slice);
  void CreateActor(int orientation, int slice);
  int ComputeCurrentOrientation();

private:
  vvImageContour();
  ~vvImageContour();
    int mPreviousTslice;
}; // end class vvImageContour
//------------------------------------------------------------------------------

#endif

