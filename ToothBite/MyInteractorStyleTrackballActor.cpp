#include "MyInteractorStyleTrackballActor.h"

#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCellPicker.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkProp3D.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkTransform.h"

#include <vtkTransformPolyDataFilter.h>
#include <QDebug>
#include "vtkRenderWindow.h"
#include "vtkRendererCollection.h"

vtkStandardNewMacro(InteractorStyleTrackballActor);  // 不

InteractorStyleTrackballActor::InteractorStyleTrackballActor()
    : m_isDoTranslation(true) {  // 初始化函数
  rotateAngle = 1;               // 每次旋转的旋转角度
  translateDistance = 1;         // 每次移动移动的距离
  this->MotionFactor = 10.0;
  this->InteractionProp =
      nullptr;  // InteractionProp是一个交互对象类，用这个类来指向我们左键点击选中的模型
  this->InteractionPicker = vtkCellPicker::New();
  this->InteractionPicker->SetTolerance(
      0.001);  // 设置鼠标点击选取模型的误差，即鼠标点击位置和模型边界之间的距离在0.001以内都认为是点到了模型
}

InteractorStyleTrackballActor::~InteractorStyleTrackballActor()  // 析构函数
{
  this->InteractionPicker->Delete();
}

void InteractorStyleTrackballActor::Rotate(
    const double x, const double y,
    const double z)  // 旋转模型，x,y,z分别是模型要沿x,y,z三个轴旋转的角度
{
  if (CurrentRenderer == nullptr ||
      InteractionProp == nullptr)  // 如果没选到任何模型或者没有舞台就直接返回
  {
    return;
  }

  vtkRenderWindowInteractor* rwi = this->Interactor;  // 拿到当前用的交互器
  vtkCamera* cam = CurrentRenderer->GetActiveCamera();  // 拿到当前用的的摄像机

  if (InteractionProp->GetUserMatrix() !=
      nullptr)  // 如果当前对象有人为给他设置的变换矩阵，就进行以下变换，这里你不用管，我们的项目不会用到这一块
  {
    vtkTransform* t = vtkTransform::New();
    t->PostMultiply();
    t->SetMatrix(InteractionProp->GetUserMatrix());
    t->RotateWXYZ(rotateAngle, x, y, z);
    InteractionProp->GetUserMatrix()->DeepCopy(t->GetMatrix());
    t->Delete();
  } else  // 我们的旋转程序
  {
    VTK_CREATE(
        vtkActorCollection,
        m_Collection);  // vtkActorCollection类是一个收集类，作用是把几个actor合在一起
    VTK_CREATE(vtkActor, test);
    InteractionProp->GetActors(m_Collection);
    m_Collection->InitTraversal();
    test =
        m_Collection->GetNextActor();  // 创建一个临时actor，指向我们选择的actor
    VTK_CREATE(vtkTransform, pTransform);  // 创建一个变换矩阵
    // pTransform->RotateWXYZ(rotateAngle, x, y, z);   //
    // 绕x,y,z轴旋转rotateAngle角度，如果xyz是001，那就是绕着z轴旋转rotateAngle度
    // VTK_CREATE(vtkTransformPolyDataFilter, pfilter);// 旋转过滤器
    // pfilter->SetInputData(test->GetMapper()->GetInput());//
    // 把选中的模型给过滤器 pfilter->SetTransform(pTransform);              //
    // 设置之前设置的旋转模式 pfilter->Update();
    // test->GetMapper()->SetInputConnection(pfilter->GetOutputPort());

    double* center = test->GetCenter();  //{ 0 };
    double* bounds = test->GetBounds();
    center[0] = 0;
    center[1] = 0;
    center[2] = 0;
    double newCenter[3] = {0};
    newCenter[0] = center[0] - ((bounds[1] - bounds[0]) / 2);
    newCenter[1] = center[1] - ((bounds[3] - bounds[2]) / 2);
    newCenter[2] = center[2] - ((bounds[5] - bounds[4]) / 2);

    pTransform->Translate(newCenter[0], newCenter[1], newCenter[2]);
    if (x > 0) {
      pTransform->RotateX(rotateAngle);
    } else if (x < 0) {
      pTransform->RotateX(-rotateAngle);
    } else if (y > 0) {
      pTransform->RotateY(rotateAngle);
    } else if (y < 0) {
      pTransform->RotateY(-rotateAngle);
    } else if (z > 0) {
      pTransform->RotateZ(rotateAngle);
    } else if (z < 0) {
      pTransform->RotateZ(-rotateAngle);
    }
    pTransform->Translate(-(newCenter[0]), -(newCenter[1]), -(newCenter[2]));
    pTransform->Update();

    VTK_CREATE(vtkTransformPolyDataFilter, pfilter);  // 旋转过滤器
    pfilter->SetInputData(
        test->GetMapper()->GetInput());  // 把选中的模型给过滤器
    pfilter->SetTransform(pTransform);   // 设置之前设置的旋转模式
    pfilter->Update();
    test->GetMapper()->SetInputConnection(pfilter->GetOutputPort());
  }

  if (AutoAdjustCameraClippingRange)  // 无用
  {
    CurrentRenderer->ResetCameraClippingRange();
  }
  rwi->Render();
}

// 平滑移动，就是平移
void InteractorStyleTrackballActor::Pan(const double x, const double y,
                                        const double z) {
  if (CurrentRenderer == nullptr || InteractionProp == nullptr) {
    return;
  }

  vtkRenderWindowInteractor* rwi = Interactor;

  // Use initial center as the origin from which to pan

  double* obj_center =
      InteractionProp
          ->GetCenter();  // 先获得当前选中模型的中心，这里的中心坐标是世界坐标系下的坐标，就是xyz坐标系下面的

  double old_pick_point[4], motion_vector[3];

  double new_location
      [3];  // 新的位置，就是平移之后的位置，是原来的中心位置x,y,z三个方向的平移距离，比如说xyz是110，那就是往x正方向和y正方向移动1，z不变
  new_location[0] = obj_center[0] + x;
  new_location[1] = obj_center[1] + y;
  new_location[2] = obj_center[2] + z;

  motion_vector[0] =
      new_location[0] -
      obj_center
          [0];  // motion_vector是一个矢量，代表了原来的位置到新的位置移动的矢量，
  motion_vector[1] = new_location[1] - obj_center[1];
  motion_vector[2] = new_location[2] - obj_center[2];

  if (InteractionProp->GetUserMatrix() != nullptr)  // 和旋转一样，不用管
  {
    qDebug() << "Pan...";
    vtkTransform* t = vtkTransform::New();
    t->PostMultiply();
    t->SetMatrix(InteractionProp->GetUserMatrix());
    t->Translate(motion_vector[0], motion_vector[1], motion_vector[2]);
    vtkSmartPointer<vtkTransformPolyDataFilter> filter =
        vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    vtkSmartPointer<vtkTransform> transform =
        vtkSmartPointer<vtkTransform>::New();
    transform->SetMatrix(t->GetMatrix());
    filter->SetTransform(transform);
    filter->SetInputData(
        dynamic_cast<vtkActor*>(InteractionProp)->GetMapper()->GetInput());
    filter->Update();
    dynamic_cast<vtkActor*>(InteractionProp)
        ->GetMapper()
        ->GetInput()
        ->DeepCopy(filter->GetOutput());
    //		InteractionProp->GetUserMatrix()->DeepCopy(t->GetMatrix());
    t->Delete();
  } else {  // 让选中的模型加上要移动的距离，就是矢量的xyz
    qDebug() << "Add position";
    vtkSmartPointer<vtkTransformPolyDataFilter> filter =
        vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    vtkSmartPointer<vtkTransform> transform =
        vtkSmartPointer<vtkTransform>::New();
    transform->Translate(motion_vector[0], motion_vector[1], motion_vector[2]);
    filter->SetTransform(transform);
    filter->SetInputData(
        dynamic_cast<vtkActor*>(InteractionProp)->GetMapper()->GetInput());
    filter->Update();
    dynamic_cast<vtkActor*>(InteractionProp)
        ->GetMapper()
        ->GetInput()
        ->DeepCopy(filter->GetOutput());
    //    InteractionProp->AddPosition(motion_vector[0], motion_vector[1],
    //                                 motion_vector[2]);
  }

  if (AutoAdjustCameraClippingRange) {
    CurrentRenderer->ResetCameraClippingRange();
  }

  rwi->Render();
}

void InteractorStyleTrackballActor::OnChar()  // 按下不同按键该如何做函数
{
  if (m_isDoTranslation)  //平移
  {
    switch (this->Interactor->GetKeyCode()) {
      case 'w':  //上
        Pan(0, 0, translateDistance);
        break;
      case 'a':  //左
        Pan(0, -translateDistance, 0);
        break;
      case 's':  //下
        Pan(0, 0, -translateDistance);
        break;
      case 'd':  //右
        Pan(0, translateDistance, 0);
        break;
      case 'q':  //前
        Pan(-translateDistance, 0, 0);
        break;
      case 'e':  //后
        Pan(translateDistance, 0, 0);
        break;
      default:
        break;
    }
  }

  else  //旋转
  {
    double x, y, z;
    x = Interactor->GetEventPosition()[0];
    y = Interactor->GetEventPosition()[1];
    z = Interactor->GetEventPosition()[2];
    // cout << x << " ," << y << " ," << z << endl;
    switch (Interactor->GetKeyCode()) {
      case 'w':  //上
        Rotate(0, -y, 0);
        break;
      case 's':  //下
        Rotate(0, y, 0);
        break;
      case 'a':  //左
        Rotate(x, 0, 0);
        break;
      case 'd':  //右
        Rotate(-x, 0, 0);
        break;
      case 'q':  //前
        Rotate(0, 0, -z);
        break;
      case 'e':  //后
        Rotate(0, 0, z);
        break;
      default:
        break;
    }
    // switch (Interactor->GetKeyCode())
    //{
    // case 'w'://上
    //	Rotate(0, -0.1, 0);
    //	break;
    // case 'a'://左
    //	Rotate(0.1, 0, 0);
    //	break;
    // case 's'://下
    //	Rotate(0, 0.1, 0);
    //	break;
    // case 'd'://右
    //	Rotate(-0.1, 0, 0);
    //	break;
    // case 'q'://前
    //	Rotate(0, 0, -0.1);
    //	break;
    // case 'e'://后
    //	Rotate(0, 0, 0.1);
    //	break;
    // default:
    //	break;
    //}
  }
}

void InteractorStyleTrackballActor::OnKeyPress() {
  auto keyChar = this->Interactor->GetKeySym();
  const char* Control_L =
      "Control_L";  // 如果按下左下角的ctrl键，就在旋转和平移这两个模式中来回切换
  auto isEqual = strcmp(keyChar, Control_L);
  if (isEqual == 0) {
    m_isDoTranslation = !m_isDoTranslation;
    if (m_isDoTranslation)
      std::cout << "translation" << std::endl;
    else
      std::cout << "rotation" << std::endl;
  }
}

//----------------------------------------------------------------------------
void InteractorStyleTrackballActor::PrintSelf(ostream& os,
                                              vtkIndent indent)  //不用管
{
  Superclass::PrintSelf(os, indent);
}

void InteractorStyleTrackballActor::OnLeftButtonDown()  // 左键按下
{
  int x = Interactor->GetEventPosition()[0];  // 获得按下的按键位置的xy值
  int y = Interactor->GetEventPosition()[1];

  FindPokedRenderer(x, y);  // 先看看是按在哪一个舞台上，因为可能会有多个舞台
  FindPickedActor(x, y);  // 找到按下的是那一个模型

  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

//----------------------------------------------------------------------------
void InteractorStyleTrackballActor::FindPickedActor(
    int x, int y)  // 找到按下的是哪一个模型
{
  InteractionPicker->Pick(x, y, 0.0, CurrentRenderer);
  vtkProp* prop =
      InteractionPicker
          ->GetViewProp();  // 这两句就是说按照上面获得的x,y,来确定模型
  if (prop != nullptr)
    InteractionProp = vtkProp3D::SafeDownCast(prop);
  else
    InteractionProp = nullptr;
}

//----------------------------------------------------------------------------
void InteractorStyleTrackballActor::Prop3DTransform(
    vtkProp3D* prop3D,  // 不用管
    double* boxCenter, int numRotation, double** rotate, double* scale) {
  vtkMatrix4x4* oldMatrix = vtkMatrix4x4::New();
  prop3D->GetMatrix(oldMatrix);

  double orig[3];
  prop3D->GetOrigin(orig);

  vtkTransform* newTransform = vtkTransform::New();
  newTransform->PostMultiply();
  if (prop3D->GetUserMatrix() != nullptr)
    newTransform->SetMatrix(prop3D->GetUserMatrix());
  else
    newTransform->SetMatrix(oldMatrix);

  newTransform->Translate(-(boxCenter[0]), -(boxCenter[1]), -(boxCenter[2]));

  for (int i = 0; i < numRotation; i++)
    newTransform->RotateWXYZ(rotate[i][0], rotate[i][1], rotate[i][2],
                             rotate[i][3]);

  if ((scale[0] * scale[1] * scale[2]) != 0.0)
    newTransform->Scale(scale[0], scale[1], scale[2]);

  newTransform->Translate(boxCenter[0], boxCenter[1], boxCenter[2]);

  // now try to get the composite of translate, rotate, and scale
  newTransform->Translate(-(orig[0]), -(orig[1]), -(orig[2]));
  newTransform->PreMultiply();
  newTransform->Translate(orig[0], orig[1], orig[2]);

  if (prop3D->GetUserMatrix() != nullptr)
    newTransform->GetMatrix(prop3D->GetUserMatrix());
  else {
    prop3D->SetPosition(newTransform->GetPosition());
    prop3D->SetScale(newTransform->GetScale());
    prop3D->SetOrientation(newTransform->GetOrientation());
  }
  oldMatrix->Delete();
  newTransform->Delete();
}
