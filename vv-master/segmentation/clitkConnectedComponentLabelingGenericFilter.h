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

#ifndef CLITKCONNECTEDCOMPONENTLABELINGSGENERICFILTER_H
#define CLITKCONNECTEDCOMPONENTLABELINGSGENERICFILTER_H

#include "clitkIO.h"
#include "clitkImageToImageGenericFilter.h"
#include "itkMacro.h"

//--------------------------------------------------------------------
namespace clitk 
{
  
  template<class ArgsInfoType>
  class ITK_EXPORT ConnectedComponentLabelingGenericFilter: 
    public ImageToImageGenericFilter<ConnectedComponentLabelingGenericFilter<ArgsInfoType> >
  {
    
  public:
    //--------------------------------------------------------------------
    ConnectedComponentLabelingGenericFilter();

    //--------------------------------------------------------------------
    typedef ConnectedComponentLabelingGenericFilter Self;
    typedef itk::SmartPointer<Self>                 Pointer;
    typedef itk::SmartPointer<const Self>           ConstPointer;

    //--------------------------------------------------------------------
    itkNewMacro(Self);  
    itkTypeMacro(ConnectedComponentLabelingGenericFilter, LightObject);

    //--------------------------------------------------------------------
    void SetArgsInfo(const ArgsInfoType & a);
    std::vector<float> GetSizeOfObjectsInPixels() const { return m_SizeOfObjectsInPixels; }
    std::vector<float> GetSizeOfObjectsInPhysicalUnits() const { return m_SizeOfObjectsInPhysicalUnits; }
    itkGetConstMacro(OriginalNumberOfObjects, unsigned long int);

    //--------------------------------------------------------------------
    // Main function called each time the filter is updated
    template<class ImageType>  
    void UpdateWithInputImageType();

  protected:
    virtual void Modified() const {} // Need for using itkMacros
    template<unsigned int Dim> void InitializeImageType();
    ArgsInfoType mArgsInfo;
    std::vector<float> m_SizeOfObjectsInPixels;
    std::vector<float> m_SizeOfObjectsInPhysicalUnits;
    unsigned long int m_OriginalNumberOfObjects;
    
  }; // end class
  //--------------------------------------------------------------------
    
} // end namespace clitk

#ifndef ITK_MANUAL_INSTANTIATION
#include "clitkConnectedComponentLabelingGenericFilter.txx"
#endif

#endif // #define CLITKCONNECTEDCOMPONENTLABELINGSGENERICFILTER_H
