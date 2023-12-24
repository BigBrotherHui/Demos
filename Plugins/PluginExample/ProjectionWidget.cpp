#include "ProjectionWidget.h"
#include "ui_ProjectionWidget.h"
#include <QPainter>
#include <QDebug>
#include <iostream>
using std::cout;
ProjectionWidget::ProjectionWidget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::ProjectionWidget)
{
    ui->setupUi(this);
}

ProjectionWidget::~ProjectionWidget()
{
    delete ui;
}

void ProjectionWidget::drawRect(Eigen::MatrixX3d points)
{
    m_points.resize(8, 3);
    for (int i = 0; i < points.rows(); ++i)
    {
        m_points.row(i) = WorldToView(points.row(i));
    }
}

Eigen::Vector3d ProjectionWidget::WorldToView(const Eigen::Vector3d& worldpt)
{
    Eigen::Vector4d tmpW;
    tmpW << worldpt, 1;
    Eigen::Vector3d ret=(getPerspectiveMatrix() * tmpW).block<3, 1>(0, 3);
    ret /= tmpW[3];
    return ret;
}

void ProjectionWidget::paintEvent(QPaintEvent* event)
{
    QPainter painter(this);
    painter.fillRect(rect(), Qt::black);
    QPen pen(Qt::white, 5);
    painter.setPen(pen);
    //m_points为经过投影变换计算后得到的归一化(-1,1)的坐标
    for(int i=0;i< m_points.rows();++i)
    {
        std::cout << width() << " " << height() << std::endl;
        Eigen::Vector2d sc;
        sc[0] = (m_points.row(i)[0] + 1.0) * 0.5 * width();
        sc[1] = (m_points.row(i)[1] + 1.0) * 0.5 * height();
        painter.drawPoint(sc[0], sc[1]);
    }
	QWidget::paintEvent(event);
}

Eigen::Matrix4d ProjectionWidget::getPerspectiveMatrix()
{
    Eigen::Matrix4d projectionMatrix;
    h_t = focal_near * tan(fov * acos(-1) * 180. / 2);
    h_b = -h_t;
    w_r = aspect * h_t;
    w_l = -w_r;
    projectionMatrix << 2 * focal_near / (w_r - w_l), 0, 0, 0, 0, 2 * focal_near / (h_t - h_b), 0, 0
        , (w_r + w_l) / (w_r - w_l), (h_t + h_b) / (h_t - h_b), -(focal_far + focal_near) / (focal_far - focal_near), -1,
        0, 0, 2 * focal_far * focal_near / (focal_far + focal_near), 0;//貌似该矩阵有问题，无法计算出正确的结果
    projectionMatrix.transposeInPlace();
    return projectionMatrix;
}

