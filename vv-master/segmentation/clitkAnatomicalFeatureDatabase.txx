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


//--------------------------------------------------------------------
template<class ImageType>
typename ImageType::Pointer AnatomicalFeatureDatabase::
GetImage(std::string tag, bool reload)
{
  if (m_MapOfTag.find(tag) == m_MapOfTag.end()) {
    clitkExceptionMacro("Could not find the tag <" << tag << "> of type 'Image' in the DB ('"
                        << GetFilename() << "')");
  }
  else {
    typename ImageType::Pointer image;
    if ((!reload) && (m_MapOfImage[tag])) {
      image = static_cast<ImageType *>(m_MapOfImage[tag]);
    }
    else {
      std::string s = m_MapOfTag[tag];
      // Read the file
      if (s[0] != '/')
        image = readImage<ImageType>(GetPath() + "/" + s);
      else
        image = readImage<ImageType>(s);
      // I add a reference count because the cache is not a smartpointer
      image->SetReferenceCount(image->GetReferenceCount()+1);
      // Insert into the cache
      m_MapOfImage.erase(tag); 
      m_MapOfImage[tag] = &(*image); // pointer
    }
    return image;
  }
}
//--------------------------------------------------------------------


//--------------------------------------------------------------------
template<class ImageType>
bool AnatomicalFeatureDatabase::
CheckImage(std::string tag)
{
  if (m_MapOfTag.find(tag) == m_MapOfTag.end()) {
    return false;
  }
  else {
    typename ImageType::Pointer image;
    if (m_MapOfImage[tag]) {
      image = static_cast<ImageType *>(m_MapOfImage[tag]);
    }
    else {
      std::string s = m_MapOfTag[tag];
      // Read the file
      std::string n = GetPath()+"/"+s;
      //      image = readImage<ImageType>();
      typedef itk::ImageFileReader<ImageType> ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName(n.c_str());
      try {
        reader->Update();
      } catch(itk::ExceptionObject & err) {
        return false;
      }
    }
  }
  return true;
}
//--------------------------------------------------------------------


//--------------------------------------------------------------------
template<class ImageType>
void AnatomicalFeatureDatabase::
SetImage(TagType tag, std::string f, typename ImageType::Pointer image, bool write)
{
  SetImageFilename(tag, f);
  m_MapOfImage[tag] = &(*image);
  // I add a reference count because the cache is not a smartpointer
  image->SetReferenceCount(image->GetReferenceCount()+1);
  if (write) {
    writeImage<ImageType>(image, f);
  }
}
//--------------------------------------------------------------------


//--------------------------------------------------------------------
template<class ImageType>
void AnatomicalFeatureDatabase::
ReleaseImage(std::string tag)
{
  if (m_MapOfTag.find(tag) == m_MapOfTag.end()) {
    clitkExceptionMacro("Could not find the tag <" << tag << "> of type Image Filename in the DB");
  }
  else {
    typename ImageType::Pointer image = GetImage<ImageType>(tag);
    image->SetReferenceCount(image->GetReferenceCount()-1);
    m_MapOfImage.erase(tag);
  }
}
//--------------------------------------------------------------------
