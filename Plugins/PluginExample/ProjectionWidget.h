#ifndef PROJECTIONWIDGET_H
#define PROJECTIONWIDGET_H

#include <QWidget>
#include <Eigen/Eigen>

QT_BEGIN_NAMESPACE
namespace Ui { class ProjectionWidget; }
QT_END_NAMESPACE

class ProjectionWidget : public QWidget
{
    Q_OBJECT

public:
    ProjectionWidget(QWidget *parent = nullptr);
    ~ProjectionWidget();
    void drawRect(Eigen::MatrixX3d points);
    Eigen::Vector3d WorldToView(const Eigen::Vector3d &worldpt);
protected:
    void paintEvent(QPaintEvent* event) override;
    Eigen::Matrix4d getPerspectiveMatrix();
private:
    Ui::ProjectionWidget *ui;
    //假设视角轴为-Z轴
    //注意相机坐标系为左手系，与世界坐标系不同
    double fov{ 45 };
    double aspect{16./9};//裁剪面宽高比
    double h_b;//近裁剪面高 h_b=-h_t
    double w_r;//右侧视平面宽 w_r=h_t*aspect
    double w_l;//左侧视平面宽w_l=-w_
    double h_t;//近裁剪面高 h_t=near*tan(fov/2)
    double focal_near{0.1};//近裁剪面距离
    double focal_far{100};//远裁剪面距离
    Eigen::MatrixX3d m_points;
};
#endif // PROJECTIONWIDGET_H
