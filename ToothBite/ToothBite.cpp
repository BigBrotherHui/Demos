#include "ToothBite.h"
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <vtkAutoInit.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkMassProperties.h>
#include "vtkOFFReader.h"
#include "vtkOFFWriter.h"
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <Eigen/Eigen>
#include <QDebug>
//#undef IGL_STATIC_LIBRARY
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/opengl/glfw/Viewer.h>
#pragma execution_character_set("utf-8")
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

/*   写在前面：VTK的一般流程
                可以把vtk的使用流程当成一个大剧院演戏
                如果剧院要演习，首先我们要有一个舞台（vtkRendererWindow）。这个舞台现在是空荡荡的
                所以我们要布景，给舞台上摆放上布景的道具，形成一个场景，也叫渲染器（renderer）
                我们还需要有演员，演员需要化妆对吧。
                所以我们先去找一个素人（vtkSTLReader，读取一个素人进来）
                接着给这个素人化妆
                （vtkPolydataMapper,把vtkSTLReader读取进来的素人（vtkSTLReader-getoutput()）,
                        传给化妆师mapper->setinput(（vtkSTLReader-getoutput())）
                此时我们就可以根据剧本创建一个角色了（vtkActor），我们让化妆师画好的素人来扮演角色：
                vtkActor->setmapper(mapper);给角色指定由哪个化妆师画出来的素人演
                最后让演员站到舞台上：
                renderer->addactor(actor);
                开始表演：
                renderer->render();
                此外，还可以设置交互器来让观众可以和舞台交流，不同的交互器就有不同的交互功能，比如说灯牌只能让演员看见，鲜花却可以直接送
                交互器是interactor
*/

//定义了toothbite应用程序的用户界面，并创建了一个包含三个子菜单项的菜单栏:“打开STL文件”、“清除文件”和“将咬合区域保存为STL”。
ToothBite::ToothBite(QWidget* parent) : QMainWindow(parent) {
  ui.setupUi(this);  //初始化窗口
  setWindowIcon(QIcon("../icons/tubiao.ico"));
  vtkOutputWindow::SetGlobalWarningDisplay(0);  // 不显示警告窗口

  /*  程序说明
  这个部分是在创建菜单栏，就是最上面的“文件”，“工具”，“帮助”三个选项，一个QMenu对应一个大类
  每一个大类下面还有细分的功能，比如说“文件”，下面就细分为”打开STL文件“，”清空文件“，”保存咬合区“等
  可以参考这个演示https://blog.csdn.net/ppss177/article/details/107627617
  */
  // “文件“部分添加子功能和设置
  QMenu* menu_File = ui.menuBar->addMenu(QString(
      "文件"));  // 创建一个”文件“菜单选项，ui.menuBar中ui是整个界面，menuBar是属于ui的菜单栏（就是菜单栏那一整个长条），addMenu就是在长条里面加入一个”文件“菜单选项
  //使用QAction类创建子功能按钮，并为它们设置图标
  QAction* act_ModelOpen = menu_File->addAction(QString(
      "打开STL文件"));  // 给”文件“菜单选项下面增加一个”动作“，也就是一个名为act_ModelOpen的按钮，
  act_ModelOpen->setIcon(QIcon("../icons/files.png"));  // 这是在设置按钮的图标
  menu_File->addSeparator();  // 这句话是加入一个分隔线，就是一条横线
  QAction* act_Clear = menu_File->addAction(QString("清空文件"));
  act_Clear->setIcon(QIcon("../icons/clear.png"));
  QAction* act_Save = menu_File->addAction(QString("保存咬合区为stl"));
  act_Save->setIcon(QIcon("../icons/save.png"));

  // connect很重要，他的作用就是用来绑定时间，qt里面的术语叫做信号和槽，他的作用就是把信号和如何处理这个信号的函数（槽函数）联系起来
  // 比如说 connect(act_ModelOpen, &QAction::triggered, this,
  // &ToothBite::fun_ReadSTLFile);这一句
  // 就是检测到你按下act_ModelOpen这个按钮所产生的信号后（&QAction::triggered），就立刻让程序执行fun_ReadSTLFile这个函数
  // 表现出来就是你按下这个按钮，就弹出读取文件的窗口
  connect(act_ModelOpen, &QAction::triggered, this,
          &ToothBite::fun_ReadSTLFile);
  connect(act_Clear, &QAction::triggered, this, &ToothBite::fun_ClearSTLFile);
  connect(act_Save, &QAction::triggered, this, [=]() {
    // 这是lamda表达式，就是识别到那个信号后就执行下面的函数
    if (intersectionPolyData->GetNumberOfCells() >
        0)  // 先判断intersectionPolyData（这就是咬合产生的咬合区域）里面的单元数量是不是0，如果是0，就说明没有碰撞，没有产生碰撞区域
    {
      // 先弹出一个窗口让用户自己选择要保存的地方
      QString strPath = QFileDialog::getSaveFileName(
          NULL, QString("open"), QString(" "), QString("STL(*.STL);"));
      // 创建一个vtkSTLWriter类，这个类是用来保存stl的
      vtkSmartPointer<vtkSTLWriter> writer =
          vtkSmartPointer<vtkSTLWriter>::New();
      writer->SetFileName(
          QString(strPath).toLocal8Bit());  // 设置路径为用户所选择的
      writer->SetInputData(intersectionPolyData);  // 把咬合区域的数据传入进来
      writer->Update();                            // 更新一下
      writer->Write();                             // 最后保存
    }
  });
  // “工具“部分添加子功能和设置，具体如同上面”文件“菜单选项
  QMenu* menu_Tools = ui.menuBar->addMenu(QString("工具"));
  QAction* act_anaModelImport = menu_Tools->addAction(QString("牙颌导入"));
  act_anaModelImport->setIcon(QIcon("../icons/ana_import.png"));
  QAction* act_Fillhode = menu_Tools->addAction("补洞（去除边界边）");
  connect(act_Fillhode, &QAction::triggered, this, [&] {
    if (m_upModelPoly) {
      fun_Fillholes(m_upModelPoly, true);
    }
    if (m_lowModelPoly) {
      fun_Fillholes(m_lowModelPoly, false);
    }
    ui.qvtkWidget->renderWindow()->Render();
  });

  QAction* act_BiteAna = menu_Tools->addAction(QString("咬合分析"));
  act_BiteAna->setIcon(QIcon("../icons/analysis.png"));
  connect(act_BiteAna, &QAction::triggered, this,
          &ToothBite::fun_ImpactAnalysis);

  QAction* act_BiteAna2 = menu_Tools->addAction(QString("咬合纸分析"));
  act_BiteAna2->setIcon(QIcon("../icons/analysis.png"));
  QAction* act_EndAna = menu_Tools->addAction(QString("结束分析"));
  act_EndAna->setIcon(QIcon("../icons/end_analysis.png"));
  menu_Tools->addSeparator();
  QAction* act_Reset = menu_Tools->addAction(QString("重置窗口"));
  act_Reset->setIcon(QIcon("../icons/reset.png"));
  menu_Tools->addSeparator();
  QAction* act_Flag = menu_Tools->addAction(QString("简化模型"));
  act_Flag->setIcon(QIcon("../icons/flag.png"));
  QAction* act_Bind = menu_Tools->addAction(QString("绑定模型"));
  act_Bind->setIcon(QIcon("../icons/bind.png"));
  connect(act_Flag, &QAction::triggered, this,
          [=]()  //简化模型
          {
            if (isSimpler == false) {
              isSimpler = true;
              act_Flag->setIcon(QIcon("../icons/flag_down.png"));
            } else {
              act_Flag->setIcon(QIcon("../icons/flag.png"));
              isSimpler = false;
            }
          });
  connect(act_anaModelImport, &QAction::triggered, this,
          &ToothBite::fun_ReadAnaFile);

  /*{
          std::thread t(&ToothBite::fun_ImpactAnalysis, this);
          t.detach();
  });*/
  connect(act_BiteAna2, &QAction::triggered, this,
          &ToothBite::fun_ImpactAnalysis2);
  connect(act_Bind, &QAction::triggered, this, &ToothBite::fun_bind);
  connect(act_Reset, &QAction::triggered, this,
          [=]() {  // 重置我们的窗口，就是说让摄像机重新对准舞台的某一个地方
            renderer->GetActiveCamera()->SetFocalPoint(-150.0, 50.0,
                                                       .0);  // 对准的位置坐标
            renderer->GetActiveCamera()->SetPosition(300., 100.,
                                                     150.);  // 摄像机摆放的位置
            renderer->GetActiveCamera()->SetViewUp(
                0, 0, 1);  // 摄像师站着拍，而不是躺着拍
            renderer->GetActiveCamera()->SetParallelProjection(1);  // 正交投影
            renderer->GetActiveCamera()->SetParallelScale(220);
            ui.qvtkWidget->renderWindow()->Render();  // 交互器开始工作
            renderer->Render();                       // 开始演出
          });
  connect(act_EndAna, &QAction::triggered, this, [=]() {  // 结束分析，清屏
    if (intersectionActor != NULL)  // 如果存在咬合产生的部分（演员）
    {
      renderer->RemoveActor(intersectionActor);  // 演员从舞台离开
      renderer->Render();
      renderer->RemoveActor2D(scalarBar);
      ui.qvtkWidget->renderWindow()->Render();
    }
  });
  // “菜单“部分添加子功能和设置
  QMenu* menu_Help =
      ui.menuBar->addMenu(QString("帮助"));  // 添加菜单栏的菜单——帮助
  QAction* act_Help_use = menu_Help->addAction(QString("如何使用..."));
  act_Help_use->setIcon(QIcon("../icons/how.png"));
  connect(act_Help_use, &QAction::triggered, [=]() {
    QMessageBox::information(this, QString("使用流程"),
                             QString("请阅读说明书"));
  });
  QAction* act_Help_button = menu_Help->addAction(QString("按键说明..."));
  act_Help_button->setIcon(QIcon("../icons/keypress.png"));
  connect(act_Help_button, &QAction::triggered, [=]() {
    QMessageBox::information(
        this, QString("按键功能"),
        QString("r/R:\t\t重置显示角度和画面\nConrol_left:\t进入旋转/"
                "移动模式\ne/E:\t\t让模型沿着x轴正方向移动\nq/"
                "Q:\t\t让模型沿着x轴负方向移动\nWASD:\t\t模型上下左右移动"));
  });
  QAction* act_Settings = menu_Help->addAction(QString("设置"));
  act_Settings->setIcon(QIcon("../icons/setting.png"));
  connect(act_Settings, &QAction::triggered, [=]() {  // 设置背景颜色
    QColor color = QColorDialog::
        getColor();  // 弹出一个窗口让用户选颜色，并且记录下用户选的什么（color）
    int r, g, b;
    color.getRgb(&r, &g, &b);  // 拿到用户选择的颜色的r,g,b分量
    if (r == 0 && g == 0 && b == 0)  // 如果都是0，那就设置为默认的颜色
      color.setRgb(85, 85, 127);
    fun_ChangeBackgroudColor(color);
  });

  // tool
  // bar，工具栏，这是把上面的那些功能都单独方程一排，方便用户快速使用，mainToolBar是那些图标所在的一整个长条，每次addAction，就会往里面增加一个图标
  ui.mainToolBar->setStyleSheet("spacing: 5px;");
  ui.mainToolBar->addAction(act_ModelOpen);
  ui.mainToolBar->addAction(act_Clear);
  ui.mainToolBar->addAction(act_anaModelImport);
  ui.mainToolBar->addAction(act_BiteAna);
  ui.mainToolBar->addAction(act_BiteAna2);
  ui.mainToolBar->addAction(act_EndAna);
  ui.mainToolBar->addAction(act_Settings);
  ui.mainToolBar->addAction(act_Reset);
  ui.mainToolBar->addAction(act_Flag);
  ui.mainToolBar->addAction(act_Save);
  ui.mainToolBar->addAction(act_Bind);

  // show
  imgActor = vtkSmartPointer<vtkImageActor>::New();
  vtkQTconnect = vtkSmartPointer<vtkEventQtSlotConnect>::
      New();  // 这个是让vtk和qt能够联系起来的工具，如果不加入这个工具，那么在vtk里面产生的信号qt就无法识别到，比如说在qvtkwidet里面按下了w按键，我们就没办法在qt程序里面识别到，从而无法操作让模型运动
  renderer = vtkSmartPointer<vtkRenderer>::New();  // 创建”布景“舞台
  vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderwindow =
      vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
  ui.qvtkWidget->setRenderWindow(renderwindow);
  ui.qvtkWidget->renderWindow()->AddRenderer(
      renderer);  // GetRenderWindow（）会返回一个vtkRenderWindow类，就是那个没有布景的大舞台，我们给他加上布景renderer，舞台就有了布景
  renderer->SetBackground(85.0 / 255, 85.0 / 255,
                          127.0 / 255);  // 设置布景的颜色
  renderer->GetActiveCamera()->SetFocalPoint(
      -0, 0,
      .0);  // GetActiveCamera()会返回当前使用的摄像机类vtkCamera，这一句是让摄像机对准某个点拍摄
  renderer->GetActiveCamera()->SetPosition(300., 100.,
                                           150.);  // 设置摄像机的位置
  renderer->GetActiveCamera()->SetViewUp(
      0, 0,
      1);  // 设置摄像师的脑袋方向，这里是(0,0,1)，就是说摄像师脑袋沿着z轴向上，相当于摄像师站着拍，如果把最后一个设置为-1，摄像师就会倒过来拍
  renderer->GetActiveCamera()->SetParallelProjection(1);  // 正交投影
  renderer->GetActiveCamera()->SetParallelScale(220);
  //  renderer->Render();

  vtkSmartPointer<vtkTransform> axesTransformer =
      vtkSmartPointer<vtkTransform>::New();  // 创建一个移动方法
  axesTransformer->Translate(
      -20, -20, -20);  // 这个移动方法是让某个东西移动到(-20，-20，-20)的位置
  axesActor = vtkSmartPointer<vtkAxesActor>::New();  // 创建一个坐标轴
  axesActor->SetUserTransform(
      axesTransformer);  // 让坐标轴按照刚才设置的移动方法移动
  axesActor->SetTotalLength(50, 50, 50);  // 让坐标轴轴的长度是50，50，50
  axesActor->SetPickable(0);              // 不允许用户移动坐标轴
  axesActor->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(
      1, 0, 0);  // x颜色
  axesActor->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(
      0, 1, 0);  // y颜色
  axesActor->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(
      0, 0, 1);  // z颜色
  axesActor->GetXAxisCaptionActor2D()
      ->GetTextActor()
      ->SetTextScaleModeToNone();  //设置x字体大小
  axesActor->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(
      10);
  axesActor->GetYAxisCaptionActor2D()
      ->GetTextActor()
      ->SetTextScaleModeToNone();  //设置y字体大小
  axesActor->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(
      10);
  axesActor->GetZAxisCaptionActor2D()
      ->GetTextActor()
      ->SetTextScaleModeToNone();  //设置z字体大小
  axesActor->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(
      10);
  renderer->AddActor(axesActor);  // 让布景舞台上出现坐标轴这个角色

  M_InteractorStyle = vtkSmartPointer<InteractorStyleTrackballActor>::
      New();  // 创建一个交互器类型，是灯牌还是鲜花，这里创建的类型是我们自己写的一个类型，叫做InteractorStyleTrackballActor，在MyInteractorStyleTrackballActor.h中进行了声明
  ui.qvtkWidget->interactor()->SetInteractorStyle(
      M_InteractorStyle);  // 确定之前创建的交互器就用我们自己创建的类型
  ui.qvtkWidget->GetRenderWindow()
      ->GetInteractor()
      ->Initialize();  // 交互器初始化
  vtkQTconnect->Connect(
      ui.qvtkWidget->GetRenderWindow()
          ->GetInteractor(),  // 这里就是之前提到的vtk和qt之间分享信号的方法，让在vtk界面里产生的vtkCommand::KeyPressEvent()（键盘按下信号）也能被qt识别到，并且用MykeyPressEvent这个函数去处理它
      vtkCommand::KeyPressEvent, this, SLOT(MykeyPressEvent()));
}
// 类的析构函数
ToothBite::~ToothBite() {}

// fun_ClearSTLFile函数的作用是将所有已经加载的模型从场景中移除，清空屏幕，以便打开新的模型文件。
/*void ToothBite::fun_ClearSTLFile(void)    // 清空屏幕
{
        for (int i = 0; i < MAX_MODEL_NUM; i++)
        {
                if (ActorNoteArray[i] == 1)    //
如果第i个模型的标志为1，说明这个模型是读取进来的，不为1说明没有读取第i个模型
                {
                        renderer->RemoveActor(m_actors[i]);  // 角色退场
                        ActorNoteArray[i] = 0;           // 标志位重新设为0
                }
        }
        if (m_upModel != 0)     // 如果上牙和下牙模型也读取了的，就也退场
        {
                renderer->RemoveActor(m_upModel);
                if (m_lowModel != 0)
                        renderer->RemoveActor(m_lowModel);
        }
        isOkToAna = false;     // 模型均已退场，无法进行分析
        renderer->Render();
        renderer->RemoveActor2D(scalarBar);
        ui.qvtkWidget->GetRenderWindow()->GetInteractor()->Render();
}*/

void ToothBite::fun_ClearSTLFile(void)  // 清空屏幕
{
  for (int i = 0; i < MAX_MODEL_NUM; i++) {
    if (ActorNoteArray[i] ==
        1)  // 如果第i个模型的标志为1，说明这个模型是读取进来的，不为1说明没有读取第i个模型
    {
      renderer->RemoveActor(m_actors[i]);  // 角色退场
      ActorNoteArray[i] = 0;               // 标志位重新设为0
    }
  }

  // 清除上牙和下牙模型
  if (m_upModel != 0) {
    renderer->RemoveActor(m_upModel);
    m_upModel = 0;
  }
  if (m_lowModel != 0) {
    renderer->RemoveActor(m_lowModel);
    m_lowModel = 0;
  }

  // 清除所有其他加载的模型
  renderer->RemoveAllViewProps();
  renderer->AddActor(axesActor);

  isOkToAna = false;  // 模型均已退场，无法进行分析
  renderer->Render();
  renderer->RemoveActor2D(scalarBar);
  ui.qvtkWidget->GetRenderWindow()->GetInteractor()->Render();
}

void ToothBite::fun_ReadSTLFile(
    void)  // 读取 STL 模型文件的函数。函数的主要作用是将用户选择的 STL
           // 模型文件读取进来，并将其添加到渲染器中进行显示。
{
  QStringList filePath = QFileDialog::getOpenFileNames(
      this, "open", "  ",
      "STL(*.STL);");  // 弹出一个窗口让用户选择文件，并且把这个路径保存在filePath变量中，可以一次读取多个模型
  if (filePath.isEmpty())  // 如果用户没有读取，就直接返回
    return;
  else {
    for (int i = 0; i < filePath.count(); i++)  // 挨个读取刚才选的文件
    {
      if (actor_index <
          MAX_MODEL_NUM)  // 小于设置的最大文件数就读入，否则就告诉用户超员了
      {
        m_actors[actor_index] = vtkSmartPointer<vtkActor>::New();  // 创建角色
        m_mapper[actor_index] =
            vtkSmartPointer<vtkPolyDataMapper>::New();  // 创建化妆师
        m_STLreader[actor_index] =
            vtkSmartPointer<vtkSTLReader>::New();  // 创建stl读取器

        m_STLreader[actor_index]->SetFileName(
            filePath[i]
                .toLocal8Bit()
                .data());  // 用stl读取器来读取stl模型，就是读取素人
        m_STLreader[actor_index]
            ->Update();  // 更新一下，很重要，不然不会真正读取进来

        m_mapper[actor_index]->SetInputConnection(
            m_STLreader[actor_index]->GetOutputPort());  // 给素人指定化妆师

        m_actors[actor_index]->SetMapper(
            m_mapper
                [actor_index]);  // 把化妆师分配给角色，这个化妆师画出来的人就只演这一个角色
        m_actors[actor_index]->GetProperty()->SetColor(0.6, 0.6,
                                                       0.6);  // 设置角色的颜色

        renderer->AddActor(m_actors[actor_index]);
        ActorNoteArray[actor_index] = 1;  // 第i个模型读取成功，标志位设置为1

        if (isSimpler)  // 如果要简化模型，isSimpler也是标志位，只有按下简化模型按钮，才会设置为1
        {
          vtkSmartPointer<vtkQuadricDecimation> quadric =
              vtkSmartPointer<vtkQuadricDecimation>::New();  // 创建简化器
          quadric->SetInputConnection(
              m_STLreader[actor_index]
                  ->GetOutputPort());  // 把读取到的素人给简化器
          quadric->SetTargetReduction(0.9);  // 设置简化比例为90%
          quadric->Update();                 // 更新一下

          QString strPath = QFileDialog::getSaveFileName(
              NULL, QString("open"), QString(" "),
              QString("STL(*.STL);"));  // qt的类弹窗口让用户选择要保存的路径
          vtkSmartPointer<vtkSTLWriter> writer =
              vtkSmartPointer<vtkSTLWriter>::New();  // 创建stl保存器
          writer->SetFileName(
              QString(strPath).toLocal8Bit());  // 设置要保存的路径和名字
          writer->SetInputData(quadric->GetOutput());  // 把简化后的素人给保存期
          writer->Update();
          writer->Write();
        }

        actor_index++;  // 已经完成一个模型的读取了，就索引小标加1
      } else            // 如果超过最大值，就弹出警告
        QMessageBox::about(NULL, QString("警告"),
                           QString("读取模型数量超出上限"));
    }
    renderer->Render();
    renderer->GetActiveCamera()->SetFocalPoint(-100., 10., -20.);
    renderer->GetActiveCamera()->SetPosition(200., 50., 150.);
    renderer->GetActiveCamera()->SetViewUp(0, 0, 1);
    renderer->GetActiveCamera()->SetParallelProjection(1);
    renderer->GetActiveCamera()->SetParallelScale(400);
    ui.qvtkWidget->GetRenderWindow()->GetInteractor()->Render();
  }
}

void ToothBite::fun_ReadAnaFile(void)  // 读取上牙和下牙的模型
{
  QString filePath = QFileDialog::getOpenFileName(this, "选择上颌模型", "  ",
                                                  "STL(*.STL);");  // 同之前的
  if (filePath.isEmpty()) {
    isOkToAna = false;
    isOkToAna2 = false;
    return;
  } else {  // 读一个上牙模型进来
    m_upModel = vtkSmartPointer<vtkActor>::New();
    m_upModelMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    m_upModelReader = vtkSmartPointer<vtkSTLReader>::New();
    m_upModelPoly = vtkSmartPointer<vtkPolyData>::
        New();  // 叫做vtkPolyData,他的作用就是获得一个模型的几何数据，后面要用的
    m_upModelReader->SetFileName(
        filePath.toLocal8Bit().data());  //设置STL文件的路径
    m_upModelReader->Update();
    m_upModelPoly =
        m_upModelReader
            ->GetOutput();  //获取m_upModelReader对象读取的STL文件的vtkPolyData对象，并将其存储到m_upModelPoly对象中
    // fun_Fillholes(m_upModelPoly,true);
    //使用GetOutputPort()方法将m_upModelReader对象的输出端口连接到m_upModelMapper对象的输入端口上
    m_upModelMapper->SetInputConnection(m_upModelReader->GetOutputPort());
    m_upModel->SetMapper(m_upModelMapper);

    // GetProperty()方法获取m_upModel的属性对象，并使用SetColor方法设置m_upModel的颜色，使用SetOpacity方法设置m_upModel的透明度。
    m_upModel->GetProperty()->SetColor(0.6, 0.6, 0.6);
    m_upModel->GetProperty()->SetOpacity(1.0);  // 这是在设置透明度

    renderer->AddActor(m_upModel);

    renderer->Render();
    renderer->GetActiveCamera()->SetFocalPoint(-100., 10., -20.);
    renderer->GetActiveCamera()->SetPosition(200., 50., 150.);
    renderer->GetActiveCamera()->SetViewUp(0, 0, 1);
    renderer->GetActiveCamera()->SetParallelProjection(1);
    renderer->GetActiveCamera()->SetParallelScale(400);
    ui.qvtkWidget->GetRenderWindow()->GetInteractor()->Render();
    isOkToAna2 = true;
  }
  filePath = QFileDialog::getOpenFileName(this, "选择下颌模型", "  ",
                                          "STL(*.STL);");  // 同之前的
  if (filePath.isEmpty())
    return;
  else {
    m_lowModel = vtkSmartPointer<vtkActor>::New();
    m_lowModelMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    m_lowModelReader = vtkSmartPointer<vtkSTLReader>::New();
    m_lowModelPoly = vtkSmartPointer<vtkPolyData>::New();
    m_lowModelReader->SetFileName(filePath.toLocal8Bit().data());
    m_lowModelReader->Update();
    m_lowModelPoly = m_lowModelReader->GetOutput();
    // fun_Fillholes(m_lowModelPoly,false);
    m_lowModelMapper->SetInputConnection(m_lowModelReader->GetOutputPort());
    m_lowModel->SetMapper(m_lowModelMapper);
    m_lowModel->AddPosition(0, 0, 0);
    m_lowModel->GetProperty()->SetOpacity(1.0);
    m_lowModel->GetProperty()->SetColor(0.6, 0.6, 0.6);
    renderer->AddActor(m_lowModel);

    renderer->Render();
    renderer->GetActiveCamera()->SetFocalPoint(-100., 10., -20.);
    renderer->GetActiveCamera()->SetPosition(200., 50., 150.);
    renderer->GetActiveCamera()->SetViewUp(0, 0, 1);
    renderer->GetActiveCamera()->SetParallelProjection(1);
    renderer->GetActiveCamera()->SetParallelScale(400);
    ui.qvtkWidget->GetRenderWindow()->GetInteractor()->Render();
    isOkToAna = true;
  }
}

void ToothBite::fun_ImpactAnalysis(void)  //咬合分析
{
  if (isOkToAna != false) {
    //调用 fun_GetUpPolyData 和 fun_GetLowPolyData
    //函数，分别获取上下牙的几何数据。
    //    fun_GetUpPolyData();   // 先获得上牙的几何数据
    //    fun_GetLowPolyData();  // 再获得下牙的几何数据
    //    vtkSmartPointer<vtkCleanPolyData> upperJawCleaner =
    //        vtkSmartPointer<vtkCleanPolyData>::
    //            New();  // 创建一个vtkCleanPolyData类，对几何数据进行清理
    //    //
    //    vtkCleanPolyData类可以清除几何数据中的重复点、共面面片和孔洞等问题，以确保数据的完整性和一致性。
    //    upperJawCleaner->SetInputData(
    //        m_UpPolyData->GetOutput());  // 把上牙的几何数据传给cleanPolyData
    //    upperJawCleaner->Update();
    //    vtkSmartPointer<vtkTriangleFilter> upfilter =
    //        vtkSmartPointer<vtkTriangleFilter>::
    //            New();  //
    //            创建一个vtkTriangleFilter类，这个是一个三角面片过滤器类
    //    //
    //    vtkTriangleFilter类可以将几何数据中的四边形面片转化为三角形面片，以便后续处理
    //    upfilter->SetInputConnection(upperJawCleaner->GetOutputPort());
    //    upfilter->Update();

    //    vtkSmartPointer<vtkCleanPolyData> lowerJawCleaner =
    //        vtkSmartPointer<vtkCleanPolyData>::New();
    //    lowerJawCleaner->SetInputData(m_LowPolyData->GetOutput());
    //    lowerJawCleaner->Update();
    //    vtkSmartPointer<vtkTriangleFilter> lowfilter =
    //        vtkSmartPointer<vtkTriangleFilter>::New();
    //    lowfilter->SetInputConnection(lowerJawCleaner->GetOutputPort());
    //    lowfilter->Update();

    //    //
    //    使用vtkBooleanOperationPolyDataFilter求交集，这个类的原理就是当你给定两个封闭的模型后，这个类会便利两个模型里的全部点，然后判断该点是否在另一个模型内部，如果在就是相交点，然后把所有的相交点提取出来，形成一个新的模型，就是我们要的
    //    qDebug() << "start processing";
    //    vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation =
    //        vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
    //    booleanOperation
    //        ->SetOperationToIntersection();  //
    //        设置boolean类的算法是求解两个集合的交集，也就是说两个模型的碰撞区域
    //    booleanOperation->SetInputData(
    //        0,
    //        upfilter->GetOutput());  // 第一个集合是上牙经过梳理过的数据
    //    booleanOperation->SetInputData(
    //        1, lowfilter->GetOutput());  // 第二个集合是下牙经过梳理过的数据
    //    booleanOperation->Update();
    //    qDebug() << "end processing";
    vtkNew<vtkOFFWriter> writer;
    writer->SetInputData(m_upModelPoly);
    writer->SetFileName("1.off");
    writer->Update();
    writer->Write();
    vtkNew<vtkOFFWriter> writer2;
    writer2->SetInputData(m_lowModelPoly);
    writer2->SetFileName("2.off");
    writer2->Update();
    writer2->Write();
    Eigen::MatrixXd VA, VB, VC;
    Eigen::VectorXi J, I;
    Eigen::MatrixXi FA, FB, FC;
    igl::MeshBooleanType boolean_type(igl::MESH_BOOLEAN_TYPE_INTERSECT);
    igl::readOFF("1.off", VA, FA);
    igl::readOFF("2.off", VB, FB);
    igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, boolean_type, VC, FC, J);
    igl::writeOFF("3.off", VC, FC);
    vtkNew<vtkOFFReader> offreader;
    offreader->SetFileName("3.off");
    offreader->Update();
    intersectionPolyData = vtkSmartPointer<
        vtkPolyData>::New();  // 给我们咬合产生的交集区域创建一个对象来保存
    intersectionPolyData->DeepCopy(
        offreader->GetOutput());  // 把boolean类的输出，复制到上一句创建的类中
    // 判断交集的vtkPolyData对象是否为空，如果不为空，则说明上下牙在该位置存在碰撞，需要进行处理
    if (intersectionPolyData->GetNumberOfCells() > 0) {
      int numPts = intersectionPolyData->GetPoints()
                       ->GetNumberOfPoints();  // 获得咬合区域的点的数量
      vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::
          New();  // vtkFloatArray类用于表示和处理各种浮点数数据，如标量、向量、张量等。
      scalars->SetNumberOfValues(
          numPts);  // 建立一个标量数组，数组长度就是咬合区域点的数量，这意味着为咬合区域每一个点分配一个额外的属性，这个属性就是他的颜色
      intersectionPolyData->GetPointData()->SetScalars(
          scalars);  // 给咬合区域分配刚才建立的标量数组
      for (int i = 0; i < numPts; ++i) {
        //对咬合区域的每个点进行颜色标量的赋值，颜色值代表该点的相对咬合程度。(将咬合区域中每个点的相对位置关系转换为颜色标量值，并将这些标量值存储在标量数组中)
        // for循环是对每一个点的颜色标量进行赋值，具体思想是获得当前点的z位置，并且获得咬合区域的最大Z和最小Z值
        //使用intersectionPolyData->GetBounds()方法获取咬合区域的边界框信息，具体是获取咬合区域在x、y、z三个方向上的最小值和最大值。
        double min = intersectionPolyData->GetBounds()[4];
        double max = intersectionPolyData->GetBounds()[5];
        double pos = intersectionPolyData->GetPoints()->GetPoint(i)[2];
        float val =
            (pos - min) /
            (max -
             min);  // 算出相对咬合值.得到一个0~1之间的浮点数val，表示该点在咬合区域中的相对位置关系，值越接近1表示该点越靠近咬合区域的最高点，值越接近0表示该点越靠近咬合区域的最低点。
        scalars->SetValue(
            i,
            val);  // 为标量数组每一个元素赋值，刚才只是建立，是空的。现在才是赋值
        ar[static_cast<int>(val * 10)]++;
      }
      VTK_CREATE(vtkMassProperties, mass);
      mass->SetInputData(intersectionPolyData);
      mass->Update();
      double area = mass->GetSurfaceArea();
      double vol = mass->GetVolume();
      QStringList ret;
      double height = intersectionPolyData->GetBounds()[5] -
                      intersectionPolyData->GetBounds()[4];
      ret.append(
          QString("咬合高度差：%1毫米").arg(QString::number(height, 'f', 2)));
      ret.append(QString("咬合区域面积：%1平方毫米")
                     .arg(QString::number(area, 'f', 2)));
      ret.append(QString("咬合区域体积：%1立方毫米")
                     .arg(QString::number(vol, 'f', 2)));
      for (int i = 0; i < 10; i++) {
        ret.append(
            QString("咬合区域在距咬合最低点高度为%1毫米处，占比%2%")
                .arg(QString::number((i + 1) * 1.0 / 10 * height, 'f', 2))
                .arg(ar[i] * 1.0 / numPts)
                .toStdString()
                .c_str());
      }
      ui.plainTextEdit->setPlainText(ret.join("\n"));
      // 创建颜色查找表，其意义在于颜色映射，详情可以直接百度LUT技术(Look Up
      // Table)
      seriesLut = vtkSmartPointer<vtkLookupTable>::New();
      // vtkLookupTable是VTK库中的一个颜色查找表类，可以用于将标量值映射为对应的颜色值，从而实现颜色编码和可视化等功能。
      seriesLut->SetNumberOfColors(numPts);  // 指定颜色查找表中有多少种颜色
      seriesLut->SetHueRange(
          0.67,
          0.0);  // 设定HSV颜色范围，色调H取值范围为0°～360°，从红色开始按逆时针方向计算，红色为0°/0.0，绿色为120°/0.34,蓝色为240°/0.67
      seriesLut->Build();

      seriesLut2 = vtkSmartPointer<vtkLookupTable>::New();
      seriesLut2->SetNumberOfColors(numPts);  // 指定颜色查找表中有多少种颜色
      seriesLut2->SetHueRange(
          0.0,
          0.67);  // 设定HSV颜色范围，色调H取值范围为0°～360°，从红色开始按逆时针方向计算，红色为0°/0.0，绿色为120°/0.34,蓝色为240°/0.67
      seriesLut2->Build();

      //将咬合区域的几何信息和颜色信息映射为vtkPolyDataMapper对象intersectionMapper
      intersectionMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      intersectionMapper->SetInputData(
          intersectionPolyData);  // 用于指定mapper的输入数据，即咬合区域的vtkPolyData对象intersectionPolyData。
      intersectionMapper
          ->ScalarVisibilityOn();  // 开启标量显示,即将咬合区域中每个点的标量值（对应着颜色查找表中的颜色）映射到对应的颜色上，实现颜色编码和可视化。
      intersectionMapper->SetColorModeToDefault();  // 颜色模式自动
      intersectionMapper->SetScalarRange(0, 1);     // 标量范围0-1
      intersectionMapper->SetLookupTable(
          seriesLut);  // 给mapper设置颜色查找表为刚才设置的seriesLut

      //将咬合区域的几何信息和颜色信息封装在vtkActor对象intersectionActor中，使得其能够在VTK渲染窗口中显示出来，实现咬合区域的可视化效果
      intersectionActor =
          vtkSmartPointer<vtkActor>::New();  //创建了一个 vtkActor 类型的变量
                                             // intersectionActor
      intersectionActor->SetMapper(
          intersectionMapper);  //将 intersectionMapper 设置为其 mapper
      intersectionActor->GetProperty()
          ->SetInterpolationToFlat();  // SetInterpolationToFlat()方法，将插值方式设置为平面插值，使得咬合区域的表面显示更加平滑和自然

      // added 2023.5
      intersectionMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
      intersectionMapper2->SetInputData(
          intersectionPolyData);                     // 咬合区域的mapper
      intersectionMapper2->ScalarVisibilityOn();     // 开启标量显示
      intersectionMapper2->SetColorModeToDefault();  // 颜色模式自动
      intersectionMapper2->SetScalarRange(0, 1);     // 标量范围0-1
      intersectionMapper2->SetLookupTable(
          seriesLut2);  // 给mapper设置颜色查找表为我们刚才设置的

      intersectionActor2 = vtkSmartPointer<vtkActor>::New();
      intersectionActor2->SetMapper(
          intersectionMapper2);  // 让咬合区域素人去对应的化妆师
      intersectionActor2->GetProperty()->SetInterpolationToFlat();

      //创建了一个 vtkScalarBarActor 类型的变量
      // scalarBar，用于创建图例。将该图例的颜色查找表设置为
      // intersectionMapper
      //的颜色查找表，标题设置为 "Bite Level"，并显示 4 个数字标签。
      scalarBar = vtkSmartPointer<
          vtkScalarBarActor>::New();  // 这是在创建图例，就是先是在右边的长条
      scalarBar->SetLookupTable(
          intersectionMapper->GetLookupTable());  // 设置LUT
      scalarBar->SetTitle("bite level");
      scalarBar->SetNumberOfLabels(4);  // 显示几个数字

      //设置 intersectionActor 的透明度为 0.8，并将其添加到渲染器 renderer
      //中。同时，将 scalarBar 添加到 2D
      //渲染器中，以便在渲染窗口的右侧显示图例。
      intersectionActor->GetProperty()->SetOpacity(0.8);
      renderer->AddActor(intersectionActor);
      renderer->AddActor2D(scalarBar);
    } else
      QMessageBox::about(NULL, QString("提示"), QString("当前模型无碰撞"));

    ui.qvtkWidget->renderWindow()->GetInteractor()->Render();
  } else
    QMessageBox::about(NULL, QString("提示"),
                       QString("当前未导入上下牙的模型，无法分析"));
}

void ToothBite::fun_bind(void)  //咬合区域绑定功能
{
  assembly1 = vtkSmartPointer<vtkAssembly>::New();
  assembly2 = vtkSmartPointer<vtkAssembly>::New();

  assembly1->AddPart(m_lowModel);
  assembly1->AddPart(intersectionActor);

  assembly2->AddPart(intersectionActor2);
  assembly2->AddPart(m_upModel);

  renderer->AddActor(assembly1);
  renderer->AddActor(assembly2);
  renderer->RemoveActor(m_upModel);
  renderer->RemoveActor(m_lowModel);
  renderer->RemoveActor(intersectionActor);
  renderer->RemoveActor(intersectionActor2);
}

vtkPolyData* GenerateHoleCover(vtkPolyData* edges, bool pIsUp) {
  vtkPolyData* cover = vtkPolyData::New();
  vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkFloatArray> scalars =
      vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCleanPolyData> sur_filt =
      vtkSmartPointer<vtkCleanPolyData>::New();
  sur_filt->SetInputData(edges);
  sur_filt->Update();
  points->DeepCopy(sur_filt->GetOutput()->GetPoints());
  double centr[3] = {0, 0, 0};
  double h;
  if (pIsUp)
    h = 0;
  else
    h = 1e10;
  for (int i = 0; i < points->GetNumberOfPoints(); i++) {
    for (int j = 0; j < 2; j++) {
      centr[j] += points->GetPoint(i)[j];
    }
    if (pIsUp) {
      if (points->GetPoint(i)[2] > h) h = points->GetPoint(i)[2];
    } else {
      if (points->GetPoint(i)[2] < h) h = points->GetPoint(i)[2];
    }
  }
  for (int j = 0; j < 2; j++) {
    centr[j] /= points->GetNumberOfPoints();
  }
  centr[2] = h;
  vtkIdType cnt_pt = points->InsertNextPoint(centr);
  for (int i = 0; i < sur_filt->GetOutput()->GetNumberOfCells(); i++) {
    vtkCell* cell = sur_filt->GetOutput()->GetCell(i);
    vtkIdType pts3[3];
    pts3[0] = cell->GetPointId(0);
    pts3[1] = cell->GetPointId(1);
    pts3[2] = cnt_pt;
    polys->InsertNextCell(3, pts3);
  }
  cover->SetPoints(points);
  cover->SetPolys(polys);
  return cover;
}

void ToothBite::fun_Fillholes(vtkSmartPointer<vtkPolyData> pPolydata,
                              bool pIsUp) {
  vtkSmartPointer<vtkFeatureEdges> fe = vtkSmartPointer<vtkFeatureEdges>::New();
  fe->SetInputData(pPolydata);
  fe->BoundaryEdgesOn();
  fe->NonManifoldEdgesOff();
  fe->FeatureEdgesOff();
  fe->ManifoldEdgesOff();
  fe->Update();
  vtkSmartPointer<vtkPolyDataConnectivityFilter> connect =
      vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  connect->SetInputData(fe->GetOutput());
  connect->Update();
  const int ncontours = connect->GetNumberOfExtractedRegions();
  vtkSmartPointer<vtkAppendPolyData> append =
      vtkSmartPointer<vtkAppendPolyData>::New();
  append->AddInputData(pPolydata);
  for (int i = 0; i < ncontours; i++) {
    connect->AddSpecifiedRegion(i);
    connect->SetExtractionModeToSpecifiedRegions();
    connect->Update();
    vtkPolyData* edges = connect->GetOutput();
    vtkPolyData* cover = GenerateHoleCover(edges, pIsUp);
    append->AddInputData(cover);
    cover->Delete();
    connect->DeleteSpecifiedRegion(i);
  }
  append->Update();
  vtkSmartPointer<vtkCleanPolyData> cleaner =
      vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputData(append->GetOutput());
  cleaner->Update();
  vtkNew<vtkSmoothPolyDataFilter> tri;
  tri->SetInputData(cleaner->GetOutput());
  tri->Update();
  pPolydata->DeepCopy(tri->GetOutput());
}

void ToothBite::fun_ImpactAnalysis2(void)  // 咬合纸
//函数的主要功能是根据用户输入的咬合深度，
//从已导入的上牙模型中切割出与下牙咬合的部分，
//并将切割出的部分转换为二值图像进行显示。
{
  if (isOkToAna2 !=
      false)  // 检查是否已经导入上牙的模型;`isOkToAna2`变量为false表示未导入
  {
    if (imgActor != NULL)  // 如果之前有显示咬合纸结果，则先移除
    {
      renderer->RemoveActor(imgActor);
    }

    // 获取用户输入的咬合深度
    //用户输入的咬合深度会保存在变量`distance1`中
    double distance1 = QInputDialog::getDouble(
        NULL, "提示", "请输入想要咬合的深度：", 0.00, 0.00, 100.0, 3);
    // 先获得上牙的几何数据
    fun_GetUpPolyData();

    // 对上牙的几何数据进行清理，为后续分析做准备
    // vtkCleanPolyData类：该类用于对几何数据进行清理操作，以便后续分析。
    vtkSmartPointer<vtkCleanPolyData> upperJawCleaner = vtkSmartPointer<
        vtkCleanPolyData>::
        New();  // 创建一个vtkCleanPolyData类，这个类用来对几何数据进行清理，为啥要清理能，因为某些模型内部的数据可能有缺陷，会影响后续的分析
    upperJawCleaner->SetInputData(
        m_UpPolyData->GetOutput());  // 把上牙的几何数据传给cleanPolyData
    upperJawCleaner->Update();

    // 创建一个三角面片过滤器类，对数据结构进行重新梳理，为后续分析做准备
    // vtkTriangleFilter类是一个三角面片过滤器，它可以将几何数据转换为三角面片结构，便于后续处理和分析。
    vtkSmartPointer<vtkTriangleFilter> upfilter = vtkSmartPointer<
        vtkTriangleFilter>::
        New();  // 创建一个vtkTriangleFilter类，这个是一个三角面片过滤器类，之所以需要这个是因为经过上面的清理后可能会有一些数据结构被打乱，所以需要把所有点重新梳理一下
    upfilter->SetInputConnection(upperJawCleaner->GetOutputPort());
    upfilter->Update();
    //经过前面的清洗操作后，可能会有一些数据结构被打乱，导致几何数据的拓扑结构发生变化。因此，需要使用三角面片过滤器将所有点重新梳理一下，以便后续处理和分析。通过使用三角面片过滤器，可以将几何数据转换为三角面片结构，保证几何数据的拓扑结构正确无误。

    // 计算咬合平面的位置:根据上牙模型的边界，计算出咬合平面的位置。
    double xx = (m_upModel->GetBounds()[0] + m_upModel->GetBounds()[1]) / 2;
    double yy = (m_upModel->GetBounds()[2] + m_upModel->GetBounds()[3]) / 2;
    // double xx = (m_upModel->GetBounds()[0] + m_upModel->GetBounds()[1]) /
    // 2;
    double zz = m_upModel->GetBounds()[4] + distance1;
    double slicePlaneLocation[3] = {
        xx, yy,
        zz};  //在这里，切片平面的位置由模型的中心点坐标决定，即将模型的中心点作为切片平面的位置。

    // vtkPlane类：该类用于创建一个平面切割器，并设置切割平面的位置和法向量。
    // vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    VTK_CREATE(vtkPlane, plane);
    plane->SetOrigin(slicePlaneLocation);  // z最高为129
    plane->SetNormal(0, 0, 1);
    // Create cutter

    // 创建一个切割器，并设置切割函数和输入数据
    VTK_CREATE(vtkCutter, cutter);

    // SetCutFunction方法用于设置切割函数，它的参数是一个可调用对象，表示切割函数的计算方式。
    // 在这里，使用之前创建的vtkPlane对象plane作为切割函数，表示按照切片平面进行切割。
    cutter->SetCutFunction(plane);

    // cutter->SetInputData(m_STLreader[0]->GetOutput());  // 不走内存的方式
    cutter->SetInputData(upfilter->GetOutput());  // 不走内存的方式
    cutter->SetSortByToSortByValue();
    // SetSortByToSortByValue方法用于设置排序方式，将排序方式设置为按照数值排序。
    //这种排序方式可以在一定程度上提高渲染效率，使得切割结果更加清晰和平滑。

    // 获取线条
    //创建了一个vtkStripper类对象stripper，并将切割器的输出数据传递给它进行线条提取。
    //线条提取后得到的数据就是咬合圆圈。
    VTK_CREATE(vtkStripper, stripper);
    stripper->SetInputConnection(cutter->GetOutputPort());  // valid circle
    stripper->JoinContiguousSegmentsOn();
    stripper->Update();

    // that's our circle
    auto circle = stripper->GetOutput();  // 获取到的线条就是咬合圆圈

    // 准备二值化图像的体素网格 prepare the binary image's voxel grid
    //二值化的目的是将三维模型的几何信息转换为二元信息，以方便进行图像处理和分析。
    VTK_CREATE(vtkImageData, whiteImage);
    double bounds[6], spacing[3],
        spacing_value = 0.01;  // desired volume spacing
    cutter->GetOutput()->GetBounds(
        bounds);  // 获取切割器cutter的输出的边界信息，并将其存储在bounds数组中
    spacing[0] = spacing_value;
    spacing[1] = spacing_value;
    spacing[2] = spacing_value;
    whiteImage->SetSpacing(spacing);
    // 计算尺寸(dimensions)
    int dim[3];
    for (int i = 0; i < 3; i++) {
      dim[i] = static_cast<int>(
                   ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i])) +
               1;
      if (dim[i] < 1) dim[i] = 1;
    }
    whiteImage->SetDimensions(dim);
    whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);
    double origin[3];
    origin[0] = bounds[0];  // + spacing[0] / 2;
    origin[1] = bounds[2];  // + spacing[1] / 2;
    origin[2] = bounds[4];  // + spacing[2] / 2;

    //将体素网格的原点设置为bounds数组中的最小坐标值，并使用SetOrigin方法将其存储
    //这里的原点表示体素网格的最小坐标值。
    whiteImage->SetOrigin(origin);

    //使用AllocateScalars方法为二值化图像的体素网格分配内存，并将数据类型设置为VTK_UNSIGNED_CHAR，
    //表示每个体素的值为无符号字符型，占用一个字节。
    whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

    // 对咬合圆圈进行扫掠操作，得到二维图像
    // 扫掠 sweep polygonal data (this is the important thing with contours!)
    //使用vtkLinearExtrusionFilter类将咬合圆圈沿着指定的方向进行拉伸，从而生成一个圆柱体
    VTK_CREATE(vtkLinearExtrusionFilter, extruder);
    extruder->SetInputData(circle);
    extruder->SetScaleFactor(1.);                   //表示拉伸的比例为1
    extruder->SetExtrusionTypeToVectorExtrusion();  //将拉伸类型设置为向量拉伸
    extruder->SetVector(
        0, 0,
        100);  // 拉伸方向和长度，原本是1,这里表示沿着z轴方向拉伸100个单位长度
    extruder->Update();

    // 将多边形数据转换为图像模具，以便进行二值化
    // 多边形数据->图像模具polygonal data --> image stencil:
    VTK_CREATE(
        vtkPolyDataToImageStencil,
        pol2stenc);  //使用vtkPolyDataToImageStencil类将圆柱体数据转换为二值化图像的掩模
    pol2stenc->SetTolerance(
        1);  // important if extruder->SetVector(0, 0, 1) !!!
    pol2stenc->SetInputConnection(
        extruder
            ->GetOutputPort());  //将圆柱体的输出数据连接到pol2stenc的输入端口上

    //使用SetOutputOrigin和SetOutputSpacing方法分别设置掩模的原点和体素间距。
    //这里的origin和spacing分别表示圆柱体的原点和体素间距
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->Update();

    // 对白色图像进行剪切，并设置背景
    // 剪切相应的白色图像并设置背景：cut the corresponding white image and set
    // the background:
    //使用vtkImageStencil类将二值化图像掩模应用到白色图像上，将圆柱体的内部区域标记为1，圆柱体的外部区域标记为0。
    //这个处理过程可以用于将圆柱体数据映射到体素网格上，从而生成咬合面的体素网格。
    VTK_CREATE(vtkImageStencil, imgstenc);
    imgstenc->SetInputData(whiteImage);  // 用whiteimage会是灰的
    imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
    imgstenc
        ->ReverseStencilOn();  //将掩模翻转，使圆柱体内部的区域被标记为1，圆柱体外部的区域被标记为0。
    imgstenc->SetBackgroundValue(
        255);  //将背景值设置为255，即将未标记的区域设置为白色
    imgstenc->Update();  //更新体素网格化数据。

    vtkSmartPointer<vtkImageData> img = imgstenc->GetOutput();
    int imgDimension[3];
    img->GetDimensions(imgDimension);

    // 将二值图像添加到渲染器中进行显示
    imgActor->SetInputData(img);
    renderer->AddActor(imgActor);

  } else
    // 如果没有导入上牙的模型，则弹出提示框
    QMessageBox::about(NULL, QString("提示"),
                       QString("当前未导入上牙的模型，无法进行咬合纸分析"));
}

//改变渲染窗口的背景颜色
void ToothBite::fun_ChangeBackgroudColor(QColor color)  // 改变颜色
{
  int r, g, b;
  color.getRgb(&r, &g, &b);
  std::cout << r << "  " << g << "  " << b << "  " << endl;
  renderer->SetBackground(r / 255.0, g / 255., b / 255.);
}

void ToothBite::fun_GetUpPolyData() {  // 获得上牙的完全几何数据，分析的时候要用
  if (m_upModel != NULL) {
    m_UpPolyData = vtkSmartPointer<
        vtkAppendPolyData>::New();  // 创建一个几何数据集合，这会还是空的
    // m_UpPolyData->RemoveAllInputs();
    double pPos[3];
    m_upModel->GetPosition(pPos);  // 先获得当前模型的中心点位置
    vtkSmartPointer<vtkPoints> new_points =
        vtkSmartPointer<vtkPoints>::New();  // 创建一些点，这会还是空的

    // 三角面片数据的数据，获得了几何数据
    vtkSmartPointer<vtkPolyData> polydata =
        ((vtkPolyDataMapper*)m_upModel->GetMapper())->GetInput();
    new_points =
        polydata
            ->GetPoints();  // GetPoints()方法可以获取vtkPolyData对象中所有的点的坐标信息，并将其存储到一个vtkPoints对象中
    vtkIdType iSumPts = new_points->GetNumberOfPoints();  //获得三角片点集
    // vtkIdType是一种用于表示vtk数据对象中元素数量的数据类型，例如vtkPoints对象中的点数量
    // GetNumberOfPoints()方法是vtkPoints对象的一个成员函数，用于获取该对象中包含的点的数量。
    for (int j = 0; j < iSumPts; j++) {
      double tmp[3];
      new_points->GetPoint(j,
                           tmp);  // 获取该点的三维坐标，将其存储到临时暑促tmp中
      for (int k = 0; k < 3; k++) {
        tmp[k] += pPos[k];  // 将该点坐标加上模型的中心点坐标
      }
      new_points->SetPoint(j, tmp);  // 设置新的点坐标
    }
    //获取顶点数量的作用是，在遍历每个顶点时，我们需要知道总共有多少个顶点，以便确定循环的范围。这样可以确保我们遍历到所有的顶点，并对它们进行适当的处理。
    m_UpPolyData->AddInputData(polydata);  // 将三角面片数据添加到几何数据集合中
    m_UpPolyData->Update();  // 更新几何数据集合
  }
}
void ToothBite::fun_GetLowPolyData()  // 获得下牙的完全几何数据
{
  if (m_lowModel != NULL) {
    m_LowPolyData =
        vtkSmartPointer<vtkAppendPolyData>::New();  // 创建 vtkAppendPolyData
                                                    // 变量 m_LowPolyData
    double pPos[3];  // 创建数组 pPos 用于存储模型位置信息
    m_lowModel->GetPosition(pPos);  // 获取模型位置信息，存储在 pPos 中
    vtkSmartPointer<vtkPoints> new_points =
        vtkSmartPointer<vtkPoints>::New();  // 创建 vtkPoints 变量
                                            // new_points，用于存储点集数据
    vtkSmartPointer<vtkPolyData> polydata =
        ((vtkPolyDataMapper*)m_lowModel->GetMapper())
            ->GetInput();                // 获取模型几何数据
    new_points = polydata->GetPoints();  // 获取几何数据中的点集
    vtkIdType iSumPts = new_points->GetNumberOfPoints();  // //获得三角片点集
    for (int j = 0; j < iSumPts; j++)  // 遍历所有点，更新点集坐标信息
    {
      double tmp[3];  // 创建数组 tmp 存储当前点的坐标信息
      new_points->GetPoint(j, tmp);  // 获取当前点的坐标信息
      for (int k = 0; k < 3; k++)    // 更新当前点的坐标信息
      {
        tmp[k] += pPos[k];
      }
      new_points->SetPoint(j, tmp);  // 更新点集中的点坐标信息
    }
    m_LowPolyData->AddInputData(
        polydata);  // 将低模型几何数据添加到 m_LowPolyData 中
    m_LowPolyData->Update();  // 更新数据
  }
}
void ToothBite::
    MykeyPressEvent()  // 键盘按下，我们要额外处理的部分，就是按下r键，让屏幕显示重新聚焦
{
  int cases;
  switch ((int)ui.qvtkWidget->interactor()
              ->GetKeyCode())  // 从交互器获取按下的按键值
  {
    case 114:  // r 的keycode，按下r键
      renderer->GetActiveCamera()->SetFocalPoint(-300.0, 10.0, -20);
      renderer->GetActiveCamera()->SetPosition(200., 50., 150.);
      renderer->GetActiveCamera()->SetViewUp(0, 0, 1);
      renderer->GetActiveCamera()->SetParallelProjection(1);
      renderer->GetActiveCamera()->SetParallelScale(400);
      renderer->GetActiveCamera()->ComputeViewPlaneNormal();
      renderer->ResetCamera();

      ui.qvtkWidget->interactor()->Render();
      break;
    default:
      break;
  }
}
double calculateStandardDeviation(double x, double y, double z) {
  double mean = (x + y + z) / 3.0;
  double variance = ((x - mean) * (x - mean) + (y - mean) * (y - mean) +
                     (z - mean) * (z - mean)) /
                    3.0;
  double standardDeviation = std::sqrt(variance);
  return standardDeviation;
}
void ToothBite::on_pushButton_clicked() {
  double v1 = ui.lineEdit->text().toDouble();
  double v2 = ui.lineEdit_2->text().toDouble();
  double v3 = ui.lineEdit_3->text().toDouble();
  ui.plainTextEdit->appendPlainText(
      QString("所选3个高度的标准差为：%1")
          .arg(
              QString::number(calculateStandardDeviation(v1, v2, v3), 'f', 2)));
}
