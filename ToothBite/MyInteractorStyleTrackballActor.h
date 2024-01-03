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
vtkSmartPointer<type> name = vtkSmartPointer<type>::New();    // 创建一个 VTK 智能指针，使用 New() 方法构造对象并赋值给创建的智能指针

class vtkCellPicker;                           // 前向声明，表明该类在此处被引用但尚未定义，真正定义应该在后面的代码中

class InteractorStyleTrackballActor : public vtkInteractorStyleTrackballCamera   // 自己定义的交互器类别，继承于官方的交互器形式 vtkInteractorStyleTrackballCamera
{
public:
	static InteractorStyleTrackballActor* New();  // 静态函数返回 InteractorStyleTrackballActor 对象指针
	vtkTypeMacro(InteractorStyleTrackballActor, vtkInteractorStyleTrackballCamera); // 似乎是宏定义，定义了该类的一些属性
	void PrintSelf(ostream& os, vtkIndent indent) override;   // 打印自身信息，在此省略

	void OnLeftButtonDown() override; // 左键按下的处理函数，override 表示覆盖基类同名函数

	void Rotate(const double x, const double y, const double z);  // 旋转物体的方法，给出欧拉角

	void Pan(const double x, const double y, const double z); // 平移物体的方法，给出世界坐标系下移动量三元组

	void OnChar() override;            // 处理按键字符的方法
	void OnKeyPress() override;        // 处理按键的方法

protected:
	InteractorStyleTrackballActor();   // 构造函数，被保护
	~InteractorStyleTrackballActor() override;  // 虚析构函数，同样被保护

	void FindPickedActor(int x, int y);    // 根据鼠标点击位置找到被选择的 Actor

	void Prop3DTransform(vtkProp3D* prop3D,
		double* boxCenter,
		int NumRotation,
		double** rotate,
		double* scale);   // 用于执行三维物体的变换

	double MotionFactor;      // 控制移动和旋转的因子

	vtkProp3D* InteractionProp;    // 存储用户交互中的一个 prop（比如 actor 或物体）。
	vtkCellPicker* InteractionPicker;  // 存储用于从绝对坐标拾取器生成的选取路径，以便在窗口中处理用户交互。

private:
	float rotateAngle, translateDistance;  // 旋转角度和平移距离
	InteractorStyleTrackballActor(const InteractorStyleTrackballActor&) = delete;  // 禁止复制构造函数的内部实现
	void operator=(const InteractorStyleTrackballActor&) = delete;

	bool m_isDoTranslation;

	//signals:
	//	void modelTransformSignal();

};

#endif


