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
    if (fType == 0) bs = fitter.interpolateSurface(p, q, dataPoints, controlPoints, uType, vType);
    else bs = fitter.approximateSurface(p, q, bsE, bsF, dataPoints, controlPoints);

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
    }
    try
    {
        std::ifstream is("PT.json");
        cereal::JSONInputArchive iarchive(is);
        iarchive(cereal::make_nvp("PT",a));
        cout << a;
        qDebug() << a;
    }catch(exception &e)
    {
        qDebug() << "error:" << QString::fromStdString(e.what());
    }
    

    QVBoxLayout* l = new QVBoxLayout(this);
    m_rw = new QmitkRenderWindow(this);
    m_rw->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard2D);
    m_rw->GetRenderer()->GetSliceNavigationController()->SetDefaultViewDirection(mitk::AnatomicalPlane::Coronal);
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

    //qDebug() << "111";
    //generateCircle(dataNum);
    //qDebug() << "112";
    //bspCurve(bcP, *curveDataPoints, bType);
    //qDebug() << "113";
    //bspSurface(bsP, bsQ, *surfaceDataPoints, uType, vType);
    //qDebug() << "114";
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
