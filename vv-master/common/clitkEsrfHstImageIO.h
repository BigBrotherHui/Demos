/*=========================================================================
  Program:   vv                     http://www.creatis.insa-lyon.fr/rio/vv

  Authors belong to:
  - University of LYON              http://www.universite-lyon.fr/
  - L�on B�rard cancer center       http://www.centreleonberard.fr
  - CREATIS CNRS laboratory         http://www.creatis.insa-lyon.fr

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the copyright notices for more information.

  It is distributed under dual licence

  - BSD        See included LICENSE.txt file
  - CeCILL-B   http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
===========================================================================**/

#ifndef __clitkEsrfHstImageIO_h
#define __clitkEsrfHstImageIO_h

#include <itkImageIOBase.h>
#include <fstream>
#include <string.h>

namespace clitk
{

//====================================================================
// Class for reading Esrf Hst Image file format
class EsrfHstImageIO: public itk::ImageIOBase
{
public:
  /** Standard class typedefs. */
  typedef EsrfHstImageIO          Self;
  typedef itk::ImageIOBase        Superclass;
  typedef itk::SmartPointer<Self> Pointer;

  EsrfHstImageIO():Superclass() { }

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(EsrfHstImageIO, ImageIOBase);

  /*-------- This part of the interface deals with reading data. ------ */
  virtual void ReadImageInformation();
  virtual bool CanReadFile( const char* FileNameToRead );
  virtual void Read(void * buffer);

  /*-------- This part of the interfaces deals with writing data. ----- */
  virtual void WriteImageInformation(bool keepOfStream);
  virtual void WriteImageInformation() { WriteImageInformation(false); }
  virtual bool CanWriteFile(const char* filename);
  virtual void Write(const void* buffer);

protected:
  std::string m_XmlFileName;
  std::string m_RawFileName;
};

} // end namespace

#endif

