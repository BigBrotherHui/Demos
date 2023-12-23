#include "FitSpaceLineRansac.h"

FitSpaceLineRansac::FitSpaceLineRansac(const Eigen::MatrixXd& pointsMatrix, const double& threshold, const int& maxIterNum, const double&delta) :
                                       _pointsMatrix(pointsMatrix),
                                       _threshold(threshold),
                                       _maxIterNum(maxIterNum),
                                       _delta(delta)
{}

FitSpaceLineRansac::~FitSpaceLineRansac()
{}

// ********************************************************************* //
// 功能：执行拟合流程
// 输入：空
// 输出：空
// 备注：无
// ********************************************************************* //
void FitSpaceLineRansac::compute()
{
    // <1> 选内点集
    if (this->selectInliersRansac(this->_pointsMatrix, this->_threshold, this->_maxIterNum, this->_delta, this->_inLierVec) != 0)
        return;

    // <2> 深拷贝数据
    this->_inlierPointsMatrix.resize(this->_inLierVec.size(), 3);
    for (int i = 0; i < this->_inLierVec.size(); i++)
    {
        this->_inlierPointsMatrix.row(i) = this->_pointsMatrix.row(this->_inLierVec[i]);
    }

    // <3> 拟合
    this->_centerX = 0.0;
    this->_centerY = 0.0;
    this->_centerZ = 0.0;
    this->spacialLineFitting(this->_inlierPointsMatrix, this->_V, this->_centerX, this->_centerY, this->_centerZ);
}

// ********************************************************************* //
// 功能：获取拟合结果
// 输入：空
// 输出：直线方程参数 [a b c x0 y0 z0]
// 备注：直线方程 (x - x0)/a = (y - y0)/b = (z - z0)/c
// ********************************************************************* //
Eigen::Matrix<double, 1, 6> FitSpaceLineRansac::getModel() const
{
    Eigen::Matrix<double, 1, 6> model;
    model <<(this->_V)(0, 0), // a
            (this->_V)(1, 0), // b
            (this->_V)(2, 0), // c
             this->_centerX,  // x0
             this->_centerY,  // y0
             this->_centerZ;  // z0
    return model;
}

// ********************************************************************* //
// 功能：获取随机值
// 输入：min 区间下限、max 区间上限
// 输出：空
// 返回：随机值，取值范围 [min, max)
// 备注：无
// ********************************************************************* //
int FitSpaceLineRansac::randIndex(const int & min, const int & max)
{
    return rand() % (max - min) + min;
}

// ********************************************************************* //
// 功能：基于RANSAC算法的3D样本点集筛选
// 输入：pointsMatrix 3D点集
//       threshold    若点到直线距离为d；如果d < sqrt(threshold), 则该点为内点，反之为外点
//         maxIterNum   最大迭代次数
//         delta        随机选取两点距离必须大于 delta
// 输出：inLierVec    内点集合在pointsMatrix中的索引集
// 返回：-1           异常
//        0           正常
// 备注：3D点集至少包含4个点
// ********************************************************************* //
int FitSpaceLineRansac::selectInliersRansac(const Eigen::MatrixXd& pointsMatrix, const double& threshold, const int& maxIterNum, const double&delta, vector<int>& inLierVec)
{
    if (pointsMatrix.data() == nullptr || pointsMatrix.rows() < 4 ||pointsMatrix.cols() != 3)
        return -1;

    srand((unsigned)time(nullptr)); // seed of time
    int preTotal = 0;               // 上一次迭代内点总数
    int total = 0;                  // 当前迭代内点总数

    for (int i = 0; i < maxIterNum; i++)
    {
        // <1> 随机取两点p0、p1；p0与p1连线记为L0
        Eigen::Matrix<double, 2, 3> sample;
        sample.row(0) = pointsMatrix.row(randIndex(0, pointsMatrix.rows())); // p0(x0 y0 z0)
        sample.row(1) = pointsMatrix.row(randIndex(0, pointsMatrix.rows())); // p1(x1 y1 z1)

        // <2> 确保p0 p1距离适当
        Eigen::Vector2d x, y, z;
        x = sample.col(0);
        y = sample.col(1);
        z = sample.col(2);
        Eigen::Vector3d v01;
        v01 << x(0) - x(1), y(0) - y(1), z(0) - z(1);// 向量v01 = p0 - p1
        if (v01.norm() < delta) // p0 p1太近？
            continue;
        v01.normalize(); // 归一化

        // <3> 筛选点
        vector<int> idxVec;
        for (int i = 0; i < pointsMatrix.rows(); i++)
        {
            // pi(xi yi zi)
            double xi = pointsMatrix.row(i).x();
            double yi = pointsMatrix.row(i).y();
            double zi = pointsMatrix.row(i).z();

            // vi0 = pi - p0
            Eigen::Vector3d vi0(xi - x(0), yi - y(0), zi - z(0));

            // vi0.cross(v01) : pi 到 直线L0 的距离
            // norm : 平方根
            if ((vi0.cross(v01)).norm() < sqrt(threshold))
            {
                idxVec.push_back(i);
            }
        }
        // <4> 更新
        total = (int)idxVec.size();
        if (total > preTotal)
        {
            preTotal = total;
            inLierVec = idxVec;
        }

    }

    return 0;
}

// ********************************************************************* //
// 功能：基于SVD分解的直线拟合
// 输入：pointsMatrix              3D点样本集
// 输出：V                         直线方向向量
//       (centerX centerY centerZ) 3D样本点集质心坐标
// 备注：
// ********************************************************************* //
void FitSpaceLineRansac::spacialLineFitting(const Eigen::MatrixXd& pointsMatrix, Eigen::Matrix3d& V, double& centerX, double& centerY, double& centerZ)
{
    int r = (int)pointsMatrix.rows();
    Eigen::MatrixXd sumMatrixXd(1, 3);
    sumMatrixXd = pointsMatrix.colwise().sum();//求x，y，z和
    // <1>、求质心
    centerX = sumMatrixXd(0, 0) / r;
    centerY = sumMatrixXd(0, 1) / r;
    centerZ = sumMatrixXd(0, 2) / r;
    Eigen::MatrixXd centerdMatrix(r, 3);
    // <2>、去质心
    for (int i = 0; i < r; i++)
    {
        centerdMatrix.row(i)[0] = pointsMatrix.row(i)[0] - centerX;
        centerdMatrix.row(i)[1] = pointsMatrix.row(i)[1] - centerY;
        centerdMatrix.row(i)[2] = pointsMatrix.row(i)[2] - centerZ;
    }
    // <3>、SVD分解
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(centerdMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    V = svd.matrixV();
}