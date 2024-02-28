#ifndef DICOMIMAGEUTIL_H
#define DICOMIMAGEUTIL_H

#include <QObject>
#include <QBoxLayout>
#include <QThreadPool>
#include <QMutex>
#include "vtkSmartPointer.h"
#include "vtkRenderer.h"
#include "vtkPolyData.h"
#include "QVTKOpenGLNativeWidget.h"
#include "vtkDICOMMetaData.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkImageData.h"
#include "vtkDICOMReader.h"
#include "vtkImageViewer2.h"
#include "vtkStringArray.h"
#include "vtkImageBlend.h"
#include "vtkImageActor.h"
#include "mitkVector.h"

class DicomImageUtil:public QObject
{
    Q_OBJECT
public:

    DicomImageUtil(QObject* parent=nullptr);
    ~DicomImageUtil();
    void Init();
    void SetContainer(QBoxLayout* pContainer, int pIndex=0);
    void LoadEmptyDicomImage();
    vtkStringArray* GetFiles(QString pFolder, QString &pErrMessage, int pSeriesIndex=0);
    void GetDicomSeries(vtkStringArray* pDicomSeriesFiles);
    int LoadDicomImage(QString pFolder);
    int LoadSingleDicomImage(QString pFilePath);
    int LoadDicomImage(vtkStringArray* pFiles);
    int LoadRawImage(QString pRawDataFile, int* pDims= new int[2]{512,512}, double pOpacity=1,
                     int pDataScalarType=VTK_UNSIGNED_SHORT);
    void ClearRawImage();
    bool StlToMask(QString pStlFile, QString pInfoFile, QString pRawFile, double* pImageCenter,
                          double pSpacing=1, int pForegroundValue=255, int pBackgroundValue=0, int pDilateValue=0);
    bool StlToMask(QString pStlFile, QString pInfoFile, QString pRawFile, double* pImageCenter,
                          double* pSpacing, int pForegroundValue=255, int pBackgroundValue=0, int pDilateValue=0);
    bool ImageMask(vtkSmartPointer<vtkImageData> imageData, QStringList pStlFiles, double* pImageOrigin,
                          double* pImageCenter, QStringList pInfoFiles, QStringList pRawFiles,
                          QString pInfoFile, QString pRawFile, int pBackgroundValue=0, int pDilateValue=0, std::function<void (int pStep)> pCallback=nullptr);
    vtkSmartPointer<vtkImageData> CutImage(vtkSmartPointer<vtkImageData> pImageData, int pExtent[6]);
    vtkSmartPointer<vtkImageData> DilateImage(vtkSmartPointer<vtkImageData> pImageData, int pDilateValue, int pForegroundValue=255, int pBackgroundValue=0);
    void SetSlice(int pSlice);
    void SetWWWL(int pWW, int pWL);
    void ResetWWWL();
    void Zoom(double pRate);
    void ResetZoom();
    void Reset();
    void RequestUpdate();
    double* GetOrigin();
    double* GetCenter();
    double* GetSpacing();
    int* GetExtent();
    int* GetDims();
    bool IsInvert(vtkStringArray *pDicomSeriesFiles);
    bool CheckSeries(vtkStringArray *pDicomSeriesFiles);
    mitk::Vector3D GetSpacing(vtkStringArray* pFiles);
    void SetOpacity(int pNum, double pOpacity);
    void SetInteractorStyle(vtkInteractorObserver *pInteractorStyle);
    vtkRenderWindow* GetRenderWindow() const;
    vtkRenderer* GetRenderer();
    vtkRenderWindowInteractor* GetInteractor() const;

private:
    bool mEmpty;
    vtkSmartPointer< vtkGenericOpenGLRenderWindow > mRenderWindow=nullptr;
    vtkSmartPointer< vtkImageViewer2> mImageView=nullptr ;
    vtkSmartPointer< vtkImageBlend> mImageBlend=nullptr;
    vtkSmartPointer<vtkRenderer> mRenderer=nullptr;
    vtkSmartPointer<vtkImageActor> mActor=nullptr;
    QVTKOpenGLNativeWidget* mMyVTK=nullptr;
    vtkSmartPointer<vtkDICOMMetaData> mDicomMetaData;
    vtkSmartPointer<vtkImageData> mImageData;
//    vtkSmartPointer<vtkIntArray> mFileIndexArray;
//    vtkSmartPointer<vtkIntArray> mFrameIndexArray;
    vtkNew<vtkDICOMReader> mDICOMReader ;
    int mExtent[6];
    int mDims[3];
    double mOrigin[3];
    double mSpacing[3];
    double mCenter[3];
    double mDefaultZoom = 0;
    double mCameraFocalPoint[3];
    double mCameraPosition[3];
    int mDefaultWW;
    int mDefaultWL;
    int mMaxWW;
    int mMaxWL;
    int mMinWW;
    int mMinWL;
    int mWW;
    int mWL;
    QMutex mMutex;
    int mProcessCount;

    void renderImage(vtkImageData* pImage);
    void render();
    void getDefaultWWWL();
    void setSlice(int pSlice, bool pRender);
    void setWWWL(int pWW, int pWL, bool pRender);
    void resetWWWL(bool pRender);
    void zoom(double pRate, bool pRender);
    void resetZoom(bool pRender);
    int* getDisplayPixelSize();
    double getZoom();
    void initImageView();
    vtkSmartPointer<vtkImageData> maskImage(vtkSmartPointer<vtkImageData> pImageData, vtkSmartPointer<vtkImageData> pNewImageData, vtkSmartPointer<vtkPolyData> pPolyData, double *pImageOrigin, double *pImageCenter, int pBackgroundValue, int pDilateValue);
    bool isInside(vtkImageStencilData *pStencilImage, int x, int y, int z, int pDilateValue);
    void readMetaData(vtkStringArray *pDicomSeriesFiles);
    void readDicom(vtkStringArray *pDicomSeriesFiles);
    int doLoadDicomImage(std::function<vtkStringArray *()> pCallback);
};

#endif // DICOMIMAGEUTIL_H
