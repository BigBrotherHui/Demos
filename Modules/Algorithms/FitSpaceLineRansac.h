#pragma once

#include<iostream>
#include<vector>
#include<ctime>
using namespace std;

#include"Eigen/Core"
#include"Eigen/Geometry"

/*******************************************************************************************************
** 功能描述：鲁棒性空间3D点直线拟合
** 1、基于RANSAC算法的3D样本点集筛选
** 2、基于SVD分解的直线拟合
** 3、本模块依赖Eigen算法库，使用前需配置Eigen
********************************************************************************************************
** 模块名：基于RANSAC和SVD分解的3D点直线拟合
** 输  入：_pointsMatrix               3D点集
**         _threshold                  若点到直线距离为d；如果d < sqrt(threshold), 则该点为内点，反之为外点
**           _maxIterNum                 最大迭代次数
**           _delta                      随机选取两点距离必须大于 delta
** 输  出：_V                          直线方向向量
**        (_centerX _centerY _centerZ) 3D样本点集质心坐标
** 备 注: 当输入3D点数量小于4时，成员函数getModel() 输出为0向量
*******************************************************************************************************/
class FitSpaceLineRansac
{
public:
    FitSpaceLineRansac(const Eigen::MatrixXd& pointsMatrix, const double& threshold, const int& maxIterNum, const double&delta);
    ~FitSpaceLineRansac();
    // 执行拟合流程
    void compute();
    // 获取拟合结果
    Eigen::Matrix<double, 1, 6>  getModel() const;

private:
    // 获取随机数
    int randIndex(const int& min, const int& max);
    // 基于RANSAC算法的3D样本点集筛选
    int selectInliersRansac(const Eigen::MatrixXd& pointsMatrix, const double& threshold, const int& maxIterNum, const double&delta, vector<int>& inLierVec);
    // 基于SVD分解的直线拟合
    void spacialLineFitting(const Eigen::MatrixXd& pointsMatrix, Eigen::Matrix3d& V, double& centerX, double& centerY, double& centerZ);

private:
    // 输入
    Eigen::MatrixXd _pointsMatrix;       // 3D点样本集
    double          _threshold;          // 若点到直线距离为d； d < sqrt(_threshold),则该点为内点
    int             _maxIterNum;         // 最大迭代次数
    double          _delta;              // 随机选取两点距离必须大于 _delta

    vector<int>     _inLierVec;          // 内点集合inlierPointsMatrix在_pointsMatrix中的索引集
    Eigen::MatrixXd _inlierPointsMatrix; // 内点集合

    // 输出直线方程参数，直线方程:(x - x0)/a = (y - y0)/b = (z - z0)/c
    Eigen::Matrix3d _V;                  // V.col(0) = [a b c]
    double          _centerX;             // x0
    double          _centerY;            // y0
    double          _centerZ;            // z0
};
