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
#ifndef clitkDecomposeThroughErosionImageFilter_h
#define clitkDecomposeThroughErosionImageFilter_h

/* =================================================
 * @file   clitkDecomposeThroughErosionImageFilter.h
 * @author 
 * @date   
 * 
 * @brief 
 * 
 ===================================================*/


// clitk include
#include "clitkIO.h"
#include "clitkCommon.h"
#include "clitkSetBackgroundImageFilter.h"

//itk include
#include "itkImageToImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkRelabelComponentImageFilter.h"

namespace clitk 
{

  template <class InputImageType, class OutputImageType>
  class ITK_EXPORT DecomposeThroughErosionImageFilter :
    public itk::ImageToImageFilter<InputImageType, OutputImageType>
  {
  public:
    //----------------------------------------
    // ITK
    //----------------------------------------
    typedef DecomposeThroughErosionImageFilter                                                 Self;
    typedef itk::ImageToImageFilter<InputImageType, OutputImageType>  Superclass;
    typedef itk::SmartPointer<Self>                                   Pointer;
    typedef itk::SmartPointer<const Self>                             ConstPointer;
   
    // Method for creation through the object factory
    itkNewMacro(Self);  

    // Run-time type information (and related methods)
    itkTypeMacro( DecomposeThroughErosionImageFilter, ImageToImageFilter );

    /** Dimension of the domain space. */
    itkStaticConstMacro(InputImageDimension, unsigned int, Superclass::InputImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int, Superclass::OutputImageDimension);

    //----------------------------------------
    // Typedefs
    //----------------------------------------
    typedef typename OutputImageType::RegionType OutputImageRegionType;
    typedef int InternalPixelType;
    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename OutputImageType::PixelType OutputPixelType;
    typedef typename InputImageType::SizeType SizeType;

    //----------------------------------------
    // Set & Get
    //----------------------------------------    
    itkBooleanMacro(Verbose);
    itkSetMacro( Verbose, bool);
    itkGetConstReferenceMacro( Verbose, bool);
    itkBooleanMacro(FullyConnected);
    itkSetMacro( FullyConnected, bool);
    itkGetConstReferenceMacro( FullyConnected, bool);
    itkSetMacro( Lower, InputPixelType);
    itkGetConstMacro( Lower, InputPixelType);
    itkSetMacro( Upper, InputPixelType);
    itkGetConstMacro( Upper, InputPixelType);
    itkSetMacro( Inside, InternalPixelType);
    itkGetConstMacro( Inside, InternalPixelType);  
    itkSetMacro( Outside, InternalPixelType);
    itkGetConstMacro( Outside, InternalPixelType);  
    itkSetMacro( ErosionPaddingValue, OutputPixelType);
    itkGetConstMacro( ErosionPaddingValue, OutputPixelType);  
    void SetRadius ( const SizeType& s) { m_Radius=s; this->Modified();}
    SizeType GetRadius(void){return m_Radius;}
    void SetRadius(const int r) { for(uint i=0; i<InputImageDimension; i++) m_Radius[i] = r; SetRadius(m_Radius); }
    itkSetMacro( NumberOfNewLabels, unsigned int);
    itkGetConstMacro( NumberOfNewLabels, unsigned int);
    itkSetMacro( MinimumObjectSize, unsigned int);
    itkGetConstMacro( MinimumObjectSize, unsigned int);
    itkSetMacro( MinimumNumberOfIterations, unsigned int);
    itkGetConstMacro( MinimumNumberOfIterations, unsigned int);

  protected:

    //----------------------------------------  
    // Constructor & Destructor
    //----------------------------------------  
    DecomposeThroughErosionImageFilter();
    ~DecomposeThroughErosionImageFilter() {};

    //----------------------------------------  
    // Update
    //----------------------------------------  
    // Generate Data
    void GenerateData(void);
    void AllocateOutput(){;}

    //     // Threaded Generate Data
    //     void BeforeThreadedGenerateData(void );
    //     void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, int threadId );
    //     void AfterThreadedGenerateData(void );
    //     // Override defaults
    //     virtual void GenerateInputRequestedRegion();
    //     virtual void GenerateOutputInformation (void);
    //     virtual void EnlargeOutputRequestedRegion(DataObject *data);
    //     void AllocateOutputs();
    //----------------------------------------  
    // Data members
    //----------------------------------------
    bool m_Verbose;
    bool m_FullyConnected;
    InputPixelType m_Lower;
    InputPixelType m_Upper;
    OutputPixelType m_ErosionPaddingValue;
    InputPixelType m_Inside;
    InputPixelType m_Outside;
    SizeType m_Radius;
    unsigned int m_NumberOfNewLabels;
    unsigned int m_MinimumObjectSize;
    unsigned int m_MinimumNumberOfIterations;

  };


} // end namespace clitk

#ifndef ITK_MANUAL_INSTANTIATION
#include "clitkDecomposeThroughErosionImageFilter.txx"
#endif

#endif // #define clitkDecomposeThroughErosionImageFilter_h


