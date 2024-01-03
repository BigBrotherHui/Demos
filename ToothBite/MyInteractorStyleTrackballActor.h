#ifndef MYINTERACTORSTYLETRACKBALLACTOR_H
#define MYINTERACTORSTYLETRACKBALLACTOR_H

#include <vtkInteractionStyleModule.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyle.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkSelectEnclosedPoints.h>

#include <vtkRenderWindowInteractor.h>
#include <vtkPropCollection.h>
#include <vtkProp3D.h>
#include <vtkActor.h>
#include <vtkMapper.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkAutoInit.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#define VTK_CREATE(type, name) \
vtkSmartPointer<type> name = vtkSmartPointer<type>::New();    // ����һ�� VTK ����ָ�룬ʹ�� New() ����������󲢸�ֵ������������ָ��

class vtkCellPicker;                           // ǰ�����������������ڴ˴������õ���δ���壬��������Ӧ���ں���Ĵ�����

class InteractorStyleTrackballActor : public vtkInteractorStyleTrackballCamera   // �Լ�����Ľ�������𣬼̳��ڹٷ��Ľ�������ʽ vtkInteractorStyleTrackballCamera
{
public:
	static InteractorStyleTrackballActor* New();  // ��̬�������� InteractorStyleTrackballActor ����ָ��
	vtkTypeMacro(InteractorStyleTrackballActor, vtkInteractorStyleTrackballCamera); // �ƺ��Ǻ궨�壬�����˸����һЩ����
	void PrintSelf(ostream& os, vtkIndent indent) override;   // ��ӡ������Ϣ���ڴ�ʡ��

	void OnLeftButtonDown() override; // ������µĴ�������override ��ʾ���ǻ���ͬ������

	void Rotate(const double x, const double y, const double z);  // ��ת����ķ���������ŷ����

	void Pan(const double x, const double y, const double z); // ƽ������ķ�����������������ϵ���ƶ�����Ԫ��

	void OnChar() override;            // �������ַ��ķ���
	void OnKeyPress() override;        // �������ķ���

protected:
	InteractorStyleTrackballActor();   // ���캯����������
	~InteractorStyleTrackballActor() override;  // ������������ͬ��������

	void FindPickedActor(int x, int y);    // ���������λ���ҵ���ѡ��� Actor

	void Prop3DTransform(vtkProp3D* prop3D,
		double* boxCenter,
		int NumRotation,
		double** rotate,
		double* scale);   // ����ִ����ά����ı任

	double MotionFactor;      // �����ƶ�����ת������

	vtkProp3D* InteractionProp;    // �洢�û������е�һ�� prop������ actor �����壩��
	vtkCellPicker* InteractionPicker;  // �洢���ڴӾ�������ʰȡ�����ɵ�ѡȡ·�����Ա��ڴ����д����û�������

private:
	float rotateAngle, translateDistance;  // ��ת�ǶȺ�ƽ�ƾ���
	InteractorStyleTrackballActor(const InteractorStyleTrackballActor&) = delete;  // ��ֹ���ƹ��캯�����ڲ�ʵ��
	void operator=(const InteractorStyleTrackballActor&) = delete;

	bool m_isDoTranslation;

	//signals:
	//	void modelTransformSignal();

};

#endif


