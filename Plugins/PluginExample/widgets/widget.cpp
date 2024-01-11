#include "widget.h"
#include "ui_widget.h"
#include <QDir>
#include <QDebug>
#include <iostream>
#include <QmitkRenderWindow.h>
#include <mitkPointSet.h>
#include <mitkIOUtil.h>
#include <mitkDisplayActionEventHandlerStd.h>
#include "cereal/archives/json.hpp"
#include "cereal/archives/xml.hpp"
#include <fstream>
#include <random>
#include <vtkSTLReader.h>
#include <vtkPlaneSource.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkCollisionDetectionFilter.h>
#include <vtkSphereSource.h>
#include <mitkVtkScalarModeProperty.h>
#include <vtkPolyLine.h>
#include <vtkStripper.h>
#include <vtkPolygon.h>
#include <vtkTriangleFilter.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkCleanPolyData.h>
using namespace std;
namespace
{
    vector<Eigen::Vector3d>* curveDataPoints{ nullptr };
    vector<vector<Eigen::Vector3d>>* surfaceDataPoints{ nullptr };
    double radius = 0.5;
    int pType = 0, kType = 0, bType = 0, uType = 0, vType = 0;
    Fitter fitter(pType, kType);
    BSpline* bc = nullptr;
    BSplineSurface* bs = nullptr;
    int dataNum = 20, bcP = 3, bsP = 3, bsQ = 3, bcH = 3, bsE = 10, bsF = 10;
    int gType = 1; // 0 == curve, 1 == surface
    int fType = 0; // 0 == interpolation, 1 == approximation
    double error = 0;
    double err = 1e-8;
    const int sampleNum = 100;
    const double PI = acos(-1);
    const double step = 0.01;

}
using namespace std;
class PT
{
public:
    PT(int tx, int ty) : x(tx), y(ty) {}
    PT(){}
    int x, y;
    friend std::ostream & operator<<(std::ostream & os, const PT& mr);
    friend QDebug operator<<(QDebug dbg, const PT& p);
private:
    friend class cereal::access;  // 声明cereal::access为友元类，以便访问私有成员（一定要有，否则会编译报错）
    // define a serialization function
    template <class Archive>
    void serialize(Archive& ar)
    {
        ar(CEREAL_NVP(x), CEREAL_NVP(y));
    }
};
std::ostream& operator<<(std::ostream& os, const PT& mr)
{
    os << "PT(" << mr.x << ", " << mr.y << ")\n";
    return os;
}
QDebug operator<<(QDebug dbg, const PT& p)
{
    dbg << "qDebug()<<" << p.x << p.y;
    return dbg;
}
void generateCircle(int n) {
    delete curveDataPoints;
    curveDataPoints = new vector<Eigen::Vector3d>;
    double degreeStep = 360. / n;
    for (int i = 0; i <= n; ++i) {
        double radians = 3.1415926/180*(degreeStep * i);
        curveDataPoints->emplace_back(radius * cos(radians), radius * sin(radians), 0);
    }
}

void generateCylinder(int n) {
    delete surfaceDataPoints;
    surfaceDataPoints = new vector<vector<Eigen::Vector3d>>;
    double degreeStep = 360. / n, heightStep = 1. / n;
    for (int i = 0; i <= n; ++i) {
        vector<Eigen::Vector3d> row;
        double radian = 3.1415926 / 180*(degreeStep * i);
        for (int j = 0; j <= n; ++j) {
            double height = -0.5 + heightStep * j;
            row.emplace_back(radius * sin(radian), radius * cos(radian), height);
        }
        surfaceDataPoints->emplace_back(row);
    }
}
void generateSphere(int n) {
    surfaceDataPoints = new vector<vector<Eigen::Vector3d>>;
    double thetaStep = 360 / (1. * n);
    double phiStep = 360 / (1. * n);
    for (int i = 0; i <= n; ++i) {
        vector<Eigen::Vector3d> row;
        double phi = 3.1415926 / 180 * (phiStep * i);
        for (int j = 0; j <= n; ++j) {
            double theta = 3.1415926 / 180 * (thetaStep * j);
            row.emplace_back(radius * sin(phi) * cos(theta), radius * sin(phi) * sin(theta), radius * cos(phi));
        }
        surfaceDataPoints->emplace_back(row);
    }
}
void bspCurve(int p, const std::vector<Eigen::Vector3d>& dataPoints, int bType) {
    std::vector<Eigen::Vector3d> controlPoints;
    if (fType == 0) bc = fitter.interpolateCurve(p, dataPoints, controlPoints, bType);
    else bc = fitter.approximateCurve(p, bcH, dataPoints, controlPoints);
    error = 0;
    std::default_random_engine random(static_cast<int>(time(nullptr)));
    uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
    for (int i = 0; i < sampleNum; ++i) {
        double u = uniformDistribution(random);
        auto samplePoint = Eigen::Vector3d(radius * cos(u * 2 * PI), radius * sin(u * 2 * PI), 0);
        error += (samplePoint - (*bc)(controlPoints, u)).norm();
    }
    error /= radius * sampleNum;

    int dataNum = static_cast<int>(dataPoints.size());
    std::vector<Eigen::Vector3d> data(dataNum);
    for (int i = 0; i < dataNum; i++) {
        data[i] = dataPoints[i];
    }
    std::vector<int> pointIndex(dataNum);
    for (int i = 0; i < dataNum; i++) {
        pointIndex[i] = i;
    }
    int h = static_cast<int>(controlPoints.size()) - 1;

    std::vector<int> controlPointIndex;
    for (int i = 0; i <= h; i++) {
        controlPointIndex.push_back(i);
    }
    std::vector<int> lineIndex;
    for (int i = 0; i <= h - 1; i++) {
        lineIndex.push_back(i);
        lineIndex.push_back(i + 1);
    }
    std::vector<Eigen::Vector3d> curvePoints;
    for (double u = 0; u < 1 + err; u += step) {
        curvePoints.emplace_back((*bc)(controlPoints, u));
    }
    std::vector<int> curveIndex;
    for (int i = 0; i < curvePoints.size() - 1; i++) {
        curveIndex.push_back(i);
        curveIndex.push_back(i + 1);
    }
}
void bspSurface(int p,
    int q,
    const vector<vector<Eigen::Vector3d>>& dataPoints,
    int uType,
    int vType) {
    std::vector<vector<Eigen::Vector3d>> controlPoints;
    if (fType == 0) {
        bs = fitter.interpolateSurface(p, q, dataPoints, controlPoints, uType, vType);
    }
    else {
        bs = fitter.approximateSurface(p, q, bsE, bsF, dataPoints, controlPoints);
    }

    error = 0;
    std::default_random_engine random(static_cast<int>(time(nullptr)));
    uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
    for (int i = 0; i < sampleNum; ++i) {
        for (int j = 0; j < sampleNum; ++j) {
            double u = uniformDistribution(random), v = uniformDistribution(random);
            double uradian = 2 * PI * u, vradian = 2 * PI * v;
            auto p1 = Eigen::Vector3d(radius * sin(uradian) * cos(vradian),
                radius * sin(uradian) * sin(vradian),
                radius * cos(uradian));
            auto p2 = (*bs)(controlPoints, u, v);
            error += static_cast<double>((p1 - p2).norm());
        }
    }
    error /= radius * sampleNum * sampleNum;

    int m = static_cast<int>(dataPoints.size()), n = static_cast<int>(dataPoints[0].size());
    std::vector<Eigen::Vector3d> data(m * n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            data[i * n + j] = dataPoints[i][j];
        }
    }
    std::vector<int> pointIndex(m * n);
    for (int i = 0; i < m * n; i++) {
        pointIndex[i] = i;
    }
    std::vector<int> controlPointIndex;
    for (int i = 0; i < m * n; i++) {
        controlPointIndex.push_back(i);
    }
    std::vector<int> lineIndex;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n - 1; ++j) {
            lineIndex.push_back(i * m + j);
            lineIndex.push_back(i * m + j + 1);
        }
    }
    for (int i = 0; i < m - 1; ++i) {
        for (int j = 0; j < n; ++j) {
            lineIndex.push_back(i * m + j);
            lineIndex.push_back((i + 1) * m + j);
        }
    }
    std::vector<Eigen::Vector3d> surfacePoints;
    int count = 0;
    for (double u = 0; u < 1 + err; u += step) {
        for (double v = 0; v < 1 + err; v += step) {
            surfacePoints.emplace_back((*bs)(controlPoints, u, v));
        }
        count++;
    }
    std::vector<int> surfaceIndex;
    for (int i = 0; i < count - 1; ++i) {
        for (int j = 0; j < count - 1; ++j) {
            int idx = i * count + j;
            surfaceIndex.push_back(idx);
            surfaceIndex.push_back(idx + count);
            surfaceIndex.push_back(idx + 1);
            surfaceIndex.push_back(idx + count);
            surfaceIndex.push_back(idx + count + 1);
            surfaceIndex.push_back(idx + 1);
        }
    }
}
Widget::Widget(QWidget *parent)
    : WidgetBase(parent)
{
    PT a;
    a.x = 100;
    a.y = 121;
    {
	    std::ofstream os("PT.json",std::ios::trunc);
	    cereal::JSONOutputArchive ar(os);
	    ar(cereal::make_nvp("PT",a));
        os.close();
        os.flush();
    }
    try
    {
        std::ifstream is("PT.json");
        cereal::JSONInputArchive iarchive(is);
        iarchive(cereal::make_nvp("PT",a));
        is.close();
        cout << a;
        qDebug() << a;
    }catch(exception &e)
    {
        qDebug() << "error:" << QString::fromStdString(e.what());
    }
    QVBoxLayout* l = new QVBoxLayout(this);
    m_rw = new QmitkRenderWindow(this);
    m_rw->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard3D);
    //m_rw->GetRenderer()->GetSliceNavigationController()->SetDefaultViewDirection(mitk::AnatomicalPlane::Coronal);
    m_ds = mitk::StandaloneDataStorage::New();
    m_rw->GetRenderer()->SetDataStorage(m_ds);
    m_DisplayActionEventBroadcast = mitk::DisplayActionEventBroadcast::New();
    m_DisplayActionEventBroadcast->LoadStateMachine("DisplayInteraction.xml");
    m_DisplayActionEventBroadcast->SetEventConfig("DisplayConfigMITKBase.xml");
    m_DisplayActionEventBroadcast->AddEventConfig("DisplayConfigCrosshair.xml");
    m_DisplayActionEventHandler = std::make_unique<mitk::DisplayActionEventHandlerStd>();
    m_DisplayActionEventHandler->SetObservableBroadcast(m_DisplayActionEventBroadcast);
    m_DisplayActionEventHandler->InitActions();
    l->addWidget(m_rw);
    l->setContentsMargins(0, 0, 0, 0);

    generateCircle(dataNum);
    generateSphere(dataNum);
    bspCurve(bcP, *curveDataPoints, bType);
    bspSurface(bsP, bsQ, *surfaceDataPoints, uType, vType);


    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("D:\\kasystem\\build\\bin\\x64\\Release\\prosthesisdata\\MeshSTL\\Femur_left.stl");
    reader->Update();

    mitk::Surface::Pointer sur = mitk::Surface::New();
    sur->SetVtkPolyData(reader->GetOutput());

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(sur);
    //node->SetBoolProperty("scalar visibility", 1);
    //mitk::VtkScalarModeProperty::Pointer scalarMode2 = mitk::VtkScalarModeProperty::New();
    //scalarMode2->SetScalarModeToCellData();
    //node->SetProperty("scalar mode", scalarMode2);
    node->SetName("femur");
    node->SetProperty("layer", mitk::IntProperty::New(9));
    //addNode(node);

    vtkNew<vtkPlaneSource> sp;
    Eigen::Vector3d tcutpoint{reader->GetOutput()->GetCenter()};
    sp->SetOrigin(Eigen::Vector3d(tcutpoint + 200 * Eigen::Vector3d::UnitX()).data());
    sp->SetPoint1(Eigen::Vector3d(tcutpoint - 200 * Eigen::Vector3d::UnitX()).data());
    sp->SetPoint2(Eigen::Vector3d(tcutpoint - 100 * Eigen::Vector3d::UnitY()).data());
    //sp->SetXResolution(100);
    //sp->SetYResolution(100);
    sp->Update();
   /* vtkSmartPointer<vtkLinearExtrusionFilter> ll = vtkSmartPointer<vtkLinearExtrusionFilter>::New();
    ll->SetInputData(sp->GetOutput());
    ll->SetExtrusionTypeToNormalExtrusion();
    ll->SetVector(0, 0, 1);
    ll->Update();*/
    vtkSmartPointer<vtkPlane> plane1 = vtkSmartPointer<vtkPlane>::New();
    plane1->SetOrigin(Eigen::Vector3d(tcutpoint + 200 * Eigen::Vector3d::UnitX()).data());
    plane1->SetNormal(sp->GetNormal());
    vtkSmartPointer<vtkCutter> c1 = vtkSmartPointer<vtkCutter>::New();
    c1->SetInputData(reader->GetOutput());
    c1->SetCutFunction(plane1);
    c1->Update();

    mitk::Surface::Pointer sur2 = mitk::Surface::New();
    sur2->SetVtkPolyData(c1->GetOutput());
    //sur2->SetVtkPolyData(ll->GetOutput());
    mitk::DataNode::Pointer node2 = mitk::DataNode::New();
    node2->SetData(sur2);
    node2->SetProperty("layer", mitk::IntProperty::New(10));
    node2->SetName("plane");
    node2->SetColor(0, 1, 0);
    node2->SetFloatProperty("material.wireframeLineWidth", 10);
    //node2->SetBoolProperty("scalar visibility", 1);
    //mitk::VtkScalarModeProperty::Pointer scalarMode = mitk::VtkScalarModeProperty::New();
    //scalarMode->SetScalarModeToCellData();
    //node2->SetProperty("scalar mode", scalarMode);
    addNode(node2);

 //   vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
 //   cutter->SetInputData(inputdata);
 //   cutter->SetCutFunction(plane);
 //   cutter->Update();
 //   vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
 //   stripper->SetInputData(cutter->GetOutput());
 //   stripper->Update();
 //   vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
 //   polydata->DeepCopy(stripper->GetOutput());
 //   polydata->SetLines(stripper->GetOutput()->GetLines());
 //   vtkSmartPointer<vtkTriangleFilter > triangle = vtkSmartPointer<vtkTriangleFilter>::New();
 //   triangle->SetInputData(polydata);
 //   triangle->Update();

 //   //如果为非封闭体，就需要判断stripper结果GetLines()中单元的属性，如果单元的首节点与尾节点的序号一致，
	////才认为是封闭的多边形，否则就是非封闭的多义线，不能直接加入到Polys中，仍然要加载在Lines单元中。
 //   vtkSmartPointer<vtkCellArray> polyAry = vtkSmartPointer<vtkCellArray>::New();
 //   vtkSmartPointer<vtkCellArray> lineAry = vtkSmartPointer<vtkCellArray>::New();
 //   for (int i = 0; i < stripper->GetOutput()->GetNumberOfCells(); i++)
 //   {
 //       vtkCell *cell = stripper->GetOutput()->GetCell(i);
 //       if (cell->GetNumberOfPoints() > 0)
 //       {
 //           int iId = cell->GetPointId(0);
 //           int iId1 = cell->GetPointId(cell->GetNumberOfPoints() - 1);
 //           if (iId != iId1)
 //           {
 //               lineAry->InsertNextCell(cell);
 //           }
 //           else
 //           {
 //               polyAry->InsertNextCell(cell);
 //           }
 //       }
 //   }
 //   polydata->SetLines(lineAry);
 //   polydata->SetPolys(polyAry);
}

Widget::~Widget()
{
}

void Widget::addNode(mitk::DataNode::Pointer dt)
{
    if(!m_ds->Exists(dt))
    {
        m_ds->Add(dt);
        mitk::RenderingManager::GetInstance()->InitializeViewByBoundingObjects(m_rw->GetVtkRenderWindow(), m_ds);
    }
}

void Widget::addNodeByPoint(double* pt,std::string name)
{
    mitk::DataNode::Pointer dt = mitk::DataNode::New();
    mitk::PointSet::Pointer ps = mitk::PointSet::New();
    dt->SetData(ps);
    ps->SetPoint(0, mitk::Point3D(pt));
    dt->SetName(name);
    addNode(dt);
}

mitk::DataNode::Pointer Widget::getNode(std::string name)
{
    return m_ds->GetNamedNode(name);
}

void Widget::getPointByName(double *pt,std::string name)
{
    auto ps=dynamic_cast<mitk::PointSet *>(getNode(name)->GetData());
    if (!ps || ps->GetSize()==0)
       return;
    for(int i=0;i<3;++i)
    {
        pt[i] = ps->GetPoint(0)[i];
    }
}
