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
#include <vtkPointData.h>
#include <vtkContourFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkSTLWriter.h>
#include <vtkHull.h>
#include <vtkPointsProjectedHull.h>
#include <vtkTubeFilter.h>
#include <vtkGlyph3D.h>
#include <vtkSphereSource.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkTriangle.h>
#include <vtkAppendPolyData.h>
#include <vtkDelaunay2D.h>
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

vtkPolyData* GenerateHoleCover(vtkPolyData* edges) {
    vtkPolyData* cover = vtkPolyData::New();
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCleanPolyData> sur_filt =
        vtkSmartPointer<vtkCleanPolyData>::New();
    sur_filt->SetInputData(edges);
    sur_filt->Update();
    points->DeepCopy(sur_filt->GetOutput()->GetPoints());
    double center[3] = { 0, 0, 0 };
    double h= 0;
    for (int i = 0; i < points->GetNumberOfPoints(); i++) {
        for (int j = 0; j < 2; j++) {
            center[j] += points->GetPoint(i)[j];
        }
    }
    for (int j = 0; j < 2; j++) {
        center[j] /= points->GetNumberOfPoints();
    }
    center[2] = h;
    vtkIdType cnt_pt = points->InsertNextPoint(center);
    int total = sur_filt->GetOutput()->GetNumberOfPoints();
    for (int i = 0; i < total-1; i++) {
        vtkIdType pts3[3];
        pts3[0] = i;
        pts3[1] = i+1;
        pts3[2] = cnt_pt;
        polys->InsertNextCell(3, pts3);
    }
    
    vtkIdType pts3[3];
    pts3[0] = total-1;
    pts3[1] = 0;
    pts3[2] = cnt_pt;
    polys->InsertNextCell(3, pts3);
    cover->SetPoints(points);
    cover->SetPolys(polys);
    return cover;
}

vtkSmartPointer<vtkPolyData> transformPolyData(vtkMatrix4x4* mt, vtkPolyData* p)
{
    vtkSmartPointer<vtkPolyData> ret = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->SetMatrix(mt);
    vtkSmartPointer<vtkTransformPolyDataFilter> fi = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    fi->SetTransform(transform);
    fi->SetInputData(p);
    fi->Update();
    ret->DeepCopy(fi->GetOutput());
    return ret;
}

vtkSmartPointer<vtkPolyData> generateWall(vtkPolyData *p,Eigen::Vector3d normal, Eigen::Vector3d origin,Eigen::Vector3d exdir)
{
    /*Eigen::Matrix3d rotate = Eigen::Quaterniond::FromTwoVectors(normal, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    Eigen::Matrix4d mt = Eigen::Matrix4d::Identity();
    mt.block<3, 3>(0, 0) = rotate;

    auto vmt = vtkSmartPointer<vtkMatrix4x4>::New();
    memcpy(vmt->GetData(), mt.data(), sizeof(double) * 16);
    vmt->Transpose();*/

    //vtkSmartPointer<vtkPolyData> outp = transformPolyData(vmt, p);
    auto newpts = vtkSmartPointer<vtkPoints>::New();
    auto plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(origin.data());
    plane->SetNormal(normal.data());
    for (int i = 0; i < p->GetNumberOfPoints(); ++i)
    {
        double pt[3];
        plane->ProjectPoint(p->GetPoint(i), pt);
        newpts->InsertNextPoint(pt[0],pt[1],pt[2]);
    }
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(newpts);
    vtkSmartPointer<vtkDelaunay2D> delaunay1 = vtkSmartPointer<vtkDelaunay2D>::New();
    delaunay1->SetInputData(polyData);
    delaunay1->Update();
    //提取边界边即可提取轮廓
    auto ex = vtkSmartPointer<vtkLinearExtrusionFilter>::New();
    ex->SetInputData(delaunay1->GetOutput());
    ex->SetExtrusionTypeToVectorExtrusion();
    ex->SetVector(exdir.data());
    ex->SetScaleFactor(30);
    ex->Update();
    vtkSmartPointer<vtkPolyData> op = vtkSmartPointer<vtkPolyData>::New();
    op->DeepCopy(ex->GetOutput());
    return op;
}

void Widget::addPolyData(vtkSmartPointer<vtkPolyData> p,double opa,double *color)
{
    mitk::Surface::Pointer sur = mitk::Surface::New();
    sur->SetVtkPolyData(p);
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(sur);
    node->SetOpacity(opa);
    node->SetColor(color[0],color[1],color[2]);
    addNode(node);
}

Widget::Widget(QWidget *parent)
    : WidgetBase(parent),ui(new Ui::Widget)
{
    ui->setupUi(this);
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
    //QVBoxLayout* l = new QVBoxLayout(ui->widget);
    m_rw = new QmitkRenderWindow();
    m_rw->show();
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
    /*l->addWidget(m_rw);
    l->setContentsMargins(0, 0, 0, 0);*/

    /*generateCircle(dataNum);
    generateSphere(dataNum);
    bspCurve(bcP, *curveDataPoints, bType);
    bspSurface(bsP, bsQ, *surfaceDataPoints, uType, vType);*/


    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:/XJT_KASystem/Cfg/ProsthesisData/科仪邦恩_NA_Femoral_UN_5_LLRM.STL");
    reader->Update();

    Eigen::Vector3d normal{0,-0.707,0.707};
    Eigen::Vector3d origin{ 0,0,0 };
    Eigen::Vector3d pt2{ 3.15 ,18.45 ,6.59 };
    double len = pt2.dot(normal);
    origin = Eigen::Vector3d(origin + normal * len);

    vtkSmartPointer<vtkPolyData> w1 = generateWall(reader->GetOutput(), normal, origin, -normal);
    double color[3]{ 0,1,0 };
    addPolyData(w1,.2,color);
    Eigen::Vector3d normal2{ 0,0,1 };
    vtkSmartPointer<vtkPolyData> w2 = generateWall(reader->GetOutput(), normal2, Eigen::Vector3d::Zero(),-normal2);
    double color2[3]{ 1,0,0 };
    addPolyData(w2, .2,color2);

    Eigen::Vector3d normal3{ 0,-0.97 ,-0.26 };
    Eigen::Vector3d pt3{ 3.15 ,25.03 ,13.17 };
    double len3 = pt3.dot(normal3);
    Eigen::Vector3d origin3 = Eigen::Vector3d(normal3 * len3);
    vtkSmartPointer<vtkPolyData> w3 = generateWall(reader->GetOutput(), normal3, origin3, -normal3);
    double color3[3]{ 0,0,1 };
    addPolyData(w3, .2, color3);


    Eigen::Matrix3d rotate = Eigen::Quaterniond::FromTwoVectors(normal, Eigen::Vector3d::UnitZ()).toRotationMatrix();

    //Eigen::Matrix4d mt=Eigen::Matrix4d::Identity();
    //mt.block<3, 3>(0, 0) = rotate;
    ////mt.block<3, 1>(0, 2) = Eigen::Vector3d::Zero();
    //auto vmt = vtkSmartPointer<vtkMatrix4x4>::New();
    //memcpy(vmt->GetData(), mt.data(), sizeof(double) * 16);
    //vmt->Transpose();


    /*vtkSmartPointer<vtkPolyData> outp=transformPolyData(vmt, reader->GetOutput());
    for(int i=0;i<outp->GetNumberOfPoints();++i)
    {
        outp->GetPoints()->SetPoint(i, outp->GetPoint(i)[0], outp->GetPoint(i)[1], len);
    }*/
    //vtkNew<vtkPointsProjectedHull> points;
    //points->DeepCopy(outp->GetPoints());
    //int xSize = points->GetSizeCCWHullZ();

    //std::unique_ptr<double[]> pts{ new double[xSize * 2] };

    //points->GetCCWHullZ(pts.get(), xSize);

    //mitk::PointSet::Pointer ps = mitk::PointSet::New();
    //vtkNew<vtkPoints> zHullPoints;
    //for (int i = 0; i < xSize; i++)
    //{
    //    double xval = pts[2 * i];
    //    double yval = pts[2 * i + 1];

    //    zHullPoints->InsertNextPoint(xval, yval,0);
    //    double p[3]{ xval, yval,0 };
    //    mitk::Point3D pt(p);
    //    ps->InsertPoint(i,pt);
    //}
    //// Insert the first point again to close the loop.
    //zHullPoints->InsertNextPoint(pts[0], pts[1],0);

    //// Display the x hull.
    //vtkNew<vtkPolyLine> zPolyLine;
    //zPolyLine->GetPointIds()->SetNumberOfIds(zHullPoints->GetNumberOfPoints());

    //for (vtkIdType i = 0; i < zHullPoints->GetNumberOfPoints(); i++)
    //{
    //    zPolyLine->GetPointIds()->SetId(i, i);
    //}

    //// Create a cell array to store the lines in and add the lines to it.
    //vtkNew<vtkCellArray> cells;
    //cells->InsertNextCell(zPolyLine);

    //// Create a polydata to store everything in.
    //vtkSmartPointer<vtkPolyData> polyData= vtkSmartPointer<vtkPolyData>::New();

    //// Add the points to the dataset.
    //polyData->SetPoints(zHullPoints);

    //// Add the lines to the dataset.
    //polyData->SetLines(cells);
    /*auto points = outp->GetPoints();
    auto newpts = vtkSmartPointer<vtkPoints>::New();
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    newpts->DeepCopy(points);
    polyData->SetPoints(newpts);
    vtkSmartPointer<vtkDelaunay2D> delaunay1 = vtkSmartPointer<vtkDelaunay2D>::New();
    delaunay1->SetInputData(polyData);
    delaunay1->Update();*/
    //auto ap = vtkSmartPointer<vtkAppendPolyData>::New();
    //ap->AddInputData(polyData);
    //ap->AddInputData(GenerateHoleCover(polyData));
    //ap->Update();


    //auto tri = vtkSmartPointer<vtkTriangleFilter>::New();
    //tri->SetInputData(polyData);
    //tri->Update();
    /*auto ex = vtkSmartPointer<vtkLinearExtrusionFilter>::New();
    ex->SetInputData(delaunay1->GetOutput());
    ex->SetExtrusionTypeToVectorExtrusion();
    ex->SetVector(0, 0, 1);
    ex->SetScaleFactor(30);
    ex->Update();*/
    //auto writer = vtkSmartPointer<vtkSTLWriter>::New();
    //writer->SetInputData(polyData);
    //writer->SetFileName("oooooooo.stl");
    //writer->Write();
    //qDebug() << polyData->GetNumberOfPoints() << "*********";
    
    /*mitk::Surface::Pointer sur = mitk::Surface::New();
    sur->SetVtkPolyData(w1);

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(sur);
    node->SetOpacity(.2);*/
    //node->SetBoolProperty("scalar visibility", 1);
    //mitk::VtkScalarModeProperty::Pointer scalarMode2 = mitk::VtkScalarModeProperty::New();
    //scalarMode2->SetScalarModeToCellData();
    //node->SetProperty("scalar mode", scalarMode2);
    //node->SetName("pro"); node->SetColor(0, 1, 0);
    //node->SetOpacity(.8);
    //node->SetProperty("layer", mitk::IntProperty::New(9));
    //addNode(node);

    mitk::Surface::Pointer sur2 = mitk::Surface::New();
    sur2->SetVtkPolyData(reader->GetOutput());
    mitk::DataNode::Pointer node2 = mitk::DataNode::New();
    node2->SetData(sur2);
    node2->SetOpacity(.5);
    addNode(node2);
    //Eigen::Matrix4d mt2 = Eigen::Matrix4d::Identity();
    //mt2.block<3, 3>(0, 0) = rotate;
    //auto vmt2 = vtkSmartPointer<vtkMatrix4x4>::New();
    //memcpy(vmt2->GetData(), mt2.data(), sizeof(double) * 16);
    //vmt2->Transpose();
    //sur2->GetGeometry()->SetIndexToWorldTransformByVtkMatrix(vmt2);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    /*vtkSmartPointer<vtkSTLReader> reader2 = vtkSmartPointer<vtkSTLReader>::New();
    reader2->SetFileName("D:\\kasystem\\build\\bin\\x64\\Release\\prosthesisdata\\companyA\\BrandTHA\\Cup_54.STL");
    reader2->Update();
    mitk::Surface::Pointer sur2 = mitk::Surface::New();
    sur2->SetVtkPolyData(reader2->GetOutput());
    mitk::DataNode::Pointer node2 = mitk::DataNode::New();
    node2->SetData(sur2);*/
    //node->SetBoolProperty("scalar visibility", 1);
    //mitk::VtkScalarModeProperty::Pointer scalarMode2 = mitk::VtkScalarModeProperty::New();
    //scalarMode2->SetScalarModeToCellData();
    //node->SetProperty("scalar mode", scalarMode2);
    //node2->SetName("stem2");
    //node->SetProperty("layer", mitk::IntProperty::New(9));
    //addNode(node2);

    //vtkNew<vtkPlaneSource> sp;
    //Eigen::Vector3d tcutpoint{reader->GetOutput()->GetCenter()};
    //sp->SetOrigin(Eigen::Vector3d(tcutpoint + 200 * Eigen::Vector3d::UnitX()).data());
    //sp->SetPoint1(Eigen::Vector3d(tcutpoint - 200 * Eigen::Vector3d::UnitX()).data());
    //sp->SetPoint2(Eigen::Vector3d(tcutpoint - 100 * Eigen::Vector3d::UnitY()).data());
    //sp->SetXResolution(100);
    //sp->SetYResolution(100);
    //sp->Update();
   /* vtkSmartPointer<vtkLinearExtrusionFilter> ll = vtkSmartPointer<vtkLinearExtrusionFilter>::New();
    ll->SetInputData(sp->GetOutput());
    ll->SetExtrusionTypeToNormalExtrusion();
    ll->SetVector(0, 0, 1);
    ll->Update();*/
    //vtkSmartPointer<vtkPlane> plane1 = vtkSmartPointer<vtkPlane>::New();
    //plane1->SetOrigin(Eigen::Vector3d(tcutpoint + 200 * Eigen::Vector3d::UnitX()).data());
    //plane1->SetNormal(sp->GetNormal());
    //vtkSmartPointer<vtkCutter> c1 = vtkSmartPointer<vtkCutter>::New();
    //c1->SetInputData(reader->GetOutput());
    //c1->SetCutFunction(plane1);
    //c1->Update();

    /*mitk::Surface::Pointer sur2 = mitk::Surface::New();
    sur2->SetVtkPolyData(c1->GetOutput());*/
    //sur2->SetVtkPolyData(ll->GetOutput());
    /*mitk::DataNode::Pointer node2 = mitk::DataNode::New();
    node2->SetData(sur2);
    node2->SetProperty("layer", mitk::IntProperty::New(10));
    node2->SetName("plane");
    node2->SetColor(0, 1, 0);
    node2->SetFloatProperty("material.wireframeLineWidth", 10);*/
    //node2->SetBoolProperty("scalar visibility", 1);
    //mitk::VtkScalarModeProperty::Pointer scalarMode = mitk::VtkScalarModeProperty::New();
    //scalarMode->SetScalarModeToCellData();
    //node2->SetProperty("scalar mode", scalarMode);
    //addNode(node2);

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
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto polydata = vtkSmartPointer<vtkPolyData>::New();
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    for (int i = 0; i < 10; ++i) {
        points->InsertNextPoint(i * 3, i * 3, i * 3);
        if (i % 2) {
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            line->GetPointIds()->SetId(0, i - 1);
            line->GetPointIds()->SetId(1, i);
            cells->InsertNextCell(line);
        }
    }
    polydata->SetPoints(points);
    polydata->SetLines(cells);
    mitk::DataNode::Pointer dd = mitk::DataNode::New();
    mitk::Surface::Pointer pp = mitk::Surface::New();
    pp->SetVtkPolyData(polydata);
    dd->SetData(pp);
    auto rep = mitk::VtkRepresentationProperty::New();
    rep->SetRepresentationToWireframe();
    dd->SetProperty("material.representation", rep);
    dd->SetProperty("material.wireframeLineWidth",mitk::FloatProperty::New(10));
    addNode(dd);
}

void Widget::on_pushButton_clicked()
{
    double x = ui->lineEditx->text().toDouble();
    double y = ui->lineEdity->text().toDouble();
    double z = ui->lineEditz->text().toDouble();
    double p1[3]{ 40.96,0,-40.95 };
    Eigen::Matrix3d rotate=Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(p1).normalized(),Eigen::Vector3d(x,y,z).normalized()).matrix();
    Eigen::Matrix4d mt;
    mt.block<3, 3>(0, 0) = rotate;
    vtkSmartPointer<vtkMatrix4x4> vmt = vtkSmartPointer<vtkMatrix4x4>::New();
    memcpy(vmt->GetData(), mt.data(), sizeof(double) * 16);
    vmt->Transpose();
    m_ds->GetNamedObject<mitk::Surface>("stem")->GetGeometry()->SetIndexToWorldTransformByVtkMatrix(vmt);
    mitk::RenderingManager::GetInstance()->RequestUpdate(m_rw->renderWindow());
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