#ifndef REGISTERFINEADJUSTDIALOG_H
#define REGISTERFINEADJUSTDIALOG_H


//#include "kernel.cuh"
//#include "register.h"
//itk
#include "itkImageFileReader.h"
#include "itkMetaImageIOFactory.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkImageDuplicator.h"

//vtk
#include "vtkImageActor.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkGenericRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
//Qt
#include <QWidget>
#include <QRadioButton>
#include <QJsonObject>
#include <array>
#include "kernel.cuh"

namespace Ui {
class RegisterFineAdjustDialog;
}
class QVTKOpenGLNativeWidget;
class GenerateMatrixHelper;
class InteractorStyleImage;
class QCheckBox;

class RegisterFineAdjustDialog : public QWidget
{
    Q_OBJECT

public:
    typedef struct AdjustInfo
    {
        double rotateX;
        double rotateY;
        double rotateZ;
        double translateX;
        double translateY;
        double translateZ;
        QString awlName;
        AdjustInfo(QString pAwlName) : awlName(pAwlName)
        {
            rotateX=0;
            rotateY=0;
            rotateZ=90;
            translateX=700;
            translateY=0;
            translateZ=0;
        }
        AdjustInfo()
        {
            rotateX=0;
            rotateY=0;
            rotateZ=90;
            translateX=700;
            translateY=0;
            translateZ=0;
        }
    }AdjustInfo;
    typedef float PixelType;
    typedef itk::Image<PixelType, 3> RawImageType;
    typedef itk::Image<PixelType, 2> XRayImageType;
    typedef itk::ImageFileReader<RawImageType> RawReaderType;
    typedef itk::ImageSeriesReader<XRayImageType> XRayReaderType ;
    typedef itk::Point<float, 3> PointType;
    typedef itk::GDCMImageIO ImageIOType;
    typedef itk::ImageDuplicator< XRayImageType > DuplicatorType;
    using RGBImageType= itk::Image<itk::RGBPixel<float>,2>;
    using GrayImageType=itk::Image<float,2>;
    explicit RegisterFineAdjustDialog(QWidget *parent = nullptr);
    ~RegisterFineAdjustDialog();
    //void setRegister(Register *pRegister);
    //读取.raw文件信息
    vtkSmartPointer<vtkImageData> readBoneawl(QString path);

    //读取.raw文件信息
    void itkReadMhdFile(QString path);
    //对raw文件插值
    RawImageType::Pointer imageinterpolation(RawImageType::Pointer imginput,double isoSpacing);
    //对raw归一化
    void CT3DNormalization(RawImageType::Pointer img,double dThreshold);
    //读取正/侧位图
    void itkReadXRayImageFront(QString path);
    void itkReadXRayImageSide(QString path);

    void ConvertedInverse(itk::Image<float,2>::Pointer img);

    //执行配准调整
    void adjust();
    void adjustByOpacity();
    //清除已选点
    void clearPoints();
    void show();
    //添加上下文菜单
    void addMenu(QVTKOpenGLNativeWidget *pTarget);
protected:
    vtkSmartPointer<vtkImageData> mergeTwoImages(vtkImageData *pimg1,vtkImageData *pimg2,double v);
    XRayImageType::Pointer mergeTwoImages(XRayImageType::Pointer pimg1,XRayImageType::Pointer pimg2,double v);

    //读取mhd文件信息，读取CenterToImage字段
    void readMetaHeaderFile(QString path);

    //读取s1torobot和s1tos2
    void readMatrixInternal(QString path);

    //反转图像像素
    void inversePixel(itk::Image<float,2>::Pointer pImage);

    void showImagesVtk(itk::Image<float,2>::Pointer imageFrontProjection,itk::Image<float,2>::Pointer imageSideProjection,double pOriginOpacity,bool pIsRgb=false);
    void showImages(vtkImageData *pMergedFront,vtkImageData *pMergedSide);
    void insertItems(QStringList items);
    QStringList checkedItems();
protected:
    virtual void keyPressEvent(QKeyEvent *event) override;

    void writerToFile(float *source,QString path,int size=0);
    void take_over();
    void WriterImageHomogenization(itk::Image<float,2>::Pointer &image);

   void PatientourceToRobotMatrix(double * par, float* ds1torobot, float *result);

   //分别将image1,image2添加到rgb图中的红绿通道中，并返回图像
   itk::Image<itk::RGBPixel<float>>::Pointer convertToRgb(itk::Image<float,2>::Pointer pImage1,itk::Image<float,2>::Pointer pImage2,double pXRayOpacity);
   void readMatrix(QString path,float *matrix);
   void initAdjustInfo();
   vtkSmartPointer<vtkMatrix4x4> matrixInverse(vtkMatrix4x4 *mt);
   //void paintEvent(QPaintEvent *e) override;
   void showEvent(QShowEvent *event) override;
   void checkChanged();
   void showRgb();
   void SaveDrrToImg(const itk::Image<PixelType, 2>::Pointer pImg,QString path);
private slots:
    void slot_importImage_clicked();
private:
    //Register *mRegister{nullptr};
    //椎节选择完毕
    bool mIsBoneAwlSelected{false};
    QString mRegistrationFile;
    //当前选中的点是否标记
    bool mIsCurrentPosMarked{false};
    //标记的椎节位置
    double mBoneAwlPos[2][2];
    //DuplicatorType::Pointer duplicatorFront,duplicatorSide;
    bool mIsImageInversed{false};
    bool mIsWindowInitialized{false};
    vtkSmartPointer<vtkRenderer> rendererFront,rendererSide;
    vtkSmartPointer<vtkGenericRenderWindowInteractor> interatorFront,interatorSide;
    QVTKOpenGLNativeWidget *wgMergedFront{nullptr},*wgMergedSide{nullptr};
    bool mIsRegistrationInfoImported{false};
    QJsonObject mJsonObject;
    std::map<QString,AdjustInfo> mAdjustInfo;
    bool mIsLastValueChanged{false};
    QString mLastPath;
    QStringList mAwlList;
    vtkSmartPointer<vtkImageActor> mActorMergedFront,mActorMergedSide;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> mRenderwindowMergedFront,mRenderwindowMergedSide;
    vtkSmartPointer<vtkImageActor> mActorProjectionFront,mActorProjectionSide;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> mRenderwindowProjectionFront,mRenderwindowProjectionSide;
    Ui::RegisterFineAdjustDialog *ui;
    RawImageType::Pointer m_rawImage;
    XRayImageType::Pointer m_frontImage,m_sideImage;
    SiddonGPU *siddonGPUFront{nullptr};//, *siddonGPUSide{ nullptr }
    //椎节中心与CT中心的距离
    float m_CenterToImage[3];
    //焦距
    float focalDistance;

    //三维节点的dimensions
    int _3dpixelNum[3];
    //三维节点的Spacing
    float _3dpixelSpacing[3];

    //侧位图的dimensions
    int imageSizeSide[2]{0};
    //侧位图的spacing
    float spacingSide[2]{0.0};

    //正位图的dimensions
    int imageSizeFront[2]{512,512};
    //正位图的spacing
    float spacingFront[2]{0.5,0.5};
    double originFront[3]{0.0};
    double originSide[3]{0.0};

    //椎节按钮与其.mhd,.raw文件的映射关系
    std::map<QString,QString> pair;
    //侧位片
    float s1tos2[16];
    //点光源到CT坐标系的转换矩阵
    float s1torobot[16];

    GenerateMatrixHelper *MatrixHelper{nullptr};

    float *resultFront{nullptr};
    float *resultSide{nullptr};


    //存储配准矩阵的结果的文件路径
    QString resultTransformFile;

    QStringList mVertebraIds;
    QList<QRadioButton *> mItems;
    std::map<QRadioButton*,std::array<itk::Image<PixelType, 2>::Pointer,2>> mImageMap;
    double mInfo[6];
    DuplicatorType::Pointer duplicatorFront;
};

#endif // REGISTERFINEADJUSTDIALOG_H
