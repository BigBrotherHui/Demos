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

vtkStandardNewMacro(InteractorStyleTrackballActor);  // ��

InteractorStyleTrackballActor::InteractorStyleTrackballActor()
    : m_isDoTranslation(true) {  // ��ʼ������
  rotateAngle = 1;               // ÿ����ת����ת�Ƕ�
  translateDistance = 1;         // ÿ���ƶ��ƶ��ľ���
  this->MotionFactor = 10.0;
  this->InteractionProp =
      nullptr;  // InteractionProp��һ�����������࣬���������ָ������������ѡ�е�ģ��
  this->InteractionPicker = vtkCellPicker::New();
  this->InteractionPicker->SetTolerance(
      0.001);  // ���������ѡȡģ�͵����������λ�ú�ģ�ͱ߽�֮��ľ�����0.001���ڶ���Ϊ�ǵ㵽��ģ��
}

InteractorStyleTrackballActor::~InteractorStyleTrackballActor()  // ��������
{
  this->InteractionPicker->Delete();
}

void InteractorStyleTrackballActor::Rotate(
    const double x, const double y,
    const double z)  // ��תģ�ͣ�x,y,z�ֱ���ģ��Ҫ��x,y,z��������ת�ĽǶ�
{
  if (CurrentRenderer == nullptr ||
      InteractionProp == nullptr)  // ���ûѡ���κ�ģ�ͻ���û����̨��ֱ�ӷ���
  {
    return;
  }

  vtkRenderWindowInteractor* rwi = this->Interactor;  // �õ���ǰ�õĽ�����
  vtkCamera* cam = CurrentRenderer->GetActiveCamera();  // �õ���ǰ�õĵ������

  if (InteractionProp->GetUserMatrix() !=
      nullptr)  // �����ǰ��������Ϊ�������õı任���󣬾ͽ������±任�������㲻�ùܣ����ǵ���Ŀ�����õ���һ��
  {
    vtkTransform* t = vtkTransform::New();
    t->PostMultiply();
    t->SetMatrix(InteractionProp->GetUserMatrix());
    t->RotateWXYZ(rotateAngle, x, y, z);
    InteractionProp->GetUserMatrix()->DeepCopy(t->GetMatrix());
    t->Delete();
  } else  // ���ǵ���ת����
  {
    VTK_CREATE(
        vtkActorCollection,
        m_Collection);  // vtkActorCollection����һ���ռ��࣬�����ǰѼ���actor����һ��
    VTK_CREATE(vtkActor, test);
    InteractionProp->GetActors(m_Collection);
    m_Collection->InitTraversal();
    test =
        m_Collection->GetNextActor();  // ����һ����ʱactor��ָ������ѡ���actor
    VTK_CREATE(vtkTransform, pTransform);  // ����һ���任����
    // pTransform->RotateWXYZ(rotateAngle, x, y, z);   //
    // ��x,y,z����תrotateAngle�Ƕȣ����xyz��001���Ǿ�������z����תrotateAngle��
    // VTK_CREATE(vtkTransformPolyDataFilter, pfilter);// ��ת������
    // pfilter->SetInputData(test->GetMapper()->GetInput());//
    // ��ѡ�е�ģ�͸������� pfilter->SetTransform(pTransform);              //
    // ����֮ǰ���õ���תģʽ pfilter->Update();
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

    VTK_CREATE(vtkTransformPolyDataFilter, pfilter);  // ��ת������
    pfilter->SetInputData(
        test->GetMapper()->GetInput());  // ��ѡ�е�ģ�͸�������
    pfilter->SetTransform(pTransform);   // ����֮ǰ���õ���תģʽ
    pfilter->Update();
    test->GetMapper()->SetInputConnection(pfilter->GetOutputPort());
  }

  if (AutoAdjustCameraClippingRange)  // ����
  {
    CurrentRenderer->ResetCameraClippingRange();
  }
  rwi->Render();
}

// ƽ���ƶ�������ƽ��
void InteractorStyleTrackballActor::Pan(const double x, const double y,
                                        const double z) {
  if (CurrentRenderer == nullptr || InteractionProp == nullptr) {
    return;
  }

  vtkRenderWindowInteractor* rwi = Interactor;

  // Use initial center as the origin from which to pan

  double* obj_center =
      InteractionProp
          ->GetCenter();  // �Ȼ�õ�ǰѡ��ģ�͵����ģ������������������������ϵ�µ����꣬����xyz����ϵ�����

  double old_pick_point[4], motion_vector[3];

  double new_location
      [3];  // �µ�λ�ã�����ƽ��֮���λ�ã���ԭ��������λ��x,y,z���������ƽ�ƾ��룬����˵xyz��110���Ǿ�����x�������y�������ƶ�1��z����
  new_location[0] = obj_center[0] + x;
  new_location[1] = obj_center[1] + y;
  new_location[2] = obj_center[2] + z;

  motion_vector[0] =
      new_location[0] -
      obj_center
          [0];  // motion_vector��һ��ʸ����������ԭ����λ�õ��µ�λ���ƶ���ʸ����
  motion_vector[1] = new_location[1] - obj_center[1];
  motion_vector[2] = new_location[2] - obj_center[2];

  if (InteractionProp->GetUserMatrix() != nullptr)  // ����תһ�������ù�
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
  } else {  // ��ѡ�е�ģ�ͼ���Ҫ�ƶ��ľ��룬����ʸ����xyz
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

void InteractorStyleTrackballActor::OnChar()  // ���²�ͬ���������������
{
  if (m_isDoTranslation)  //ƽ��
  {
    switch (this->Interactor->GetKeyCode()) {
      case 'w':  //��
        Pan(0, 0, translateDistance);
        break;
      case 'a':  //��
        Pan(0, -translateDistance, 0);
        break;
      case 's':  //��
        Pan(0, 0, -translateDistance);
        break;
      case 'd':  //��
        Pan(0, translateDistance, 0);
        break;
      case 'q':  //ǰ
        Pan(-translateDistance, 0, 0);
        break;
      case 'e':  //��
        Pan(translateDistance, 0, 0);
        break;
      default:
        break;
    }
  }

  else  //��ת
  {
    double x, y, z;
    x = Interactor->GetEventPosition()[0];
    y = Interactor->GetEventPosition()[1];
    z = Interactor->GetEventPosition()[2];
    // cout << x << " ," << y << " ," << z << endl;
    switch (Interactor->GetKeyCode()) {
      case 'w':  //��
        Rotate(0, -y, 0);
        break;
      case 's':  //��
        Rotate(0, y, 0);
        break;
      case 'a':  //��
        Rotate(x, 0, 0);
        break;
      case 'd':  //��
        Rotate(-x, 0, 0);
        break;
      case 'q':  //ǰ
        Rotate(0, 0, -z);
        break;
      case 'e':  //��
        Rotate(0, 0, z);
        break;
      default:
        break;
    }
    // switch (Interactor->GetKeyCode())
    //{
    // case 'w'://��
    //	Rotate(0, -0.1, 0);
    //	break;
    // case 'a'://��
    //	Rotate(0.1, 0, 0);
    //	break;
    // case 's'://��
    //	Rotate(0, 0.1, 0);
    //	break;
    // case 'd'://��
    //	Rotate(-0.1, 0, 0);
    //	break;
    // case 'q'://ǰ
    //	Rotate(0, 0, -0.1);
    //	break;
    // case 'e'://��
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
      "Control_L";  // ����������½ǵ�ctrl����������ת��ƽ��������ģʽ�������л�
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
                                              vtkIndent indent)  //���ù�
{
  Superclass::PrintSelf(os, indent);
}

void InteractorStyleTrackballActor::OnLeftButtonDown()  // �������
{
  int x = Interactor->GetEventPosition()[0];  // ��ð��µİ���λ�õ�xyֵ
  int y = Interactor->GetEventPosition()[1];

  FindPokedRenderer(x, y);  // �ȿ����ǰ�����һ����̨�ϣ���Ϊ���ܻ��ж����̨
  FindPickedActor(x, y);  // �ҵ����µ�����һ��ģ��

  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

//----------------------------------------------------------------------------
void InteractorStyleTrackballActor::FindPickedActor(
    int x, int y)  // �ҵ����µ�����һ��ģ��
{
  InteractionPicker->Pick(x, y, 0.0, CurrentRenderer);
  vtkProp* prop =
      InteractionPicker
          ->GetViewProp();  // ���������˵���������õ�x,y,��ȷ��ģ��
  if (prop != nullptr)
    InteractionProp = vtkProp3D::SafeDownCast(prop);
  else
    InteractionProp = nullptr;
}

//----------------------------------------------------------------------------
void InteractorStyleTrackballActor::Prop3DTransform(
    vtkProp3D* prop3D,  // ���ù�
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
