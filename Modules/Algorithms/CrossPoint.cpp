#include"CrossPoint.h"

// **************************************************************************************************************

// 功能：空间直线拟合

// 输入：pointsMatrix：3D点样本集

// 输出：V：直线方向向量；(centerX centerY centerZ):3D样本点集质心

// return: void

// **************************************************************************************************************
void spacialLineFitting(Eigen::MatrixXd& pointsMatrix, Eigen::Matrix3d& V, double& centerX, double& centerY, double& centerZ)
{
    int r = pointsMatrix.rows();
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
//字符转double
template <class Type>
Type stringToNum(string& str)
{
    istringstream iss(str);
    Type num;
    iss >> num;
    return num;
}
//读取csv文件点
void readPoint(string strPath, vector<vector<double>>& pointsVec)
{
    std::ifstream _csvInput;
    _csvInput.open(strPath, std::ios::in);
    std::string _Oneline;
    std::vector<double> _lineOfstr;
    vector<vector<double>> verticalLinePoints;
    while (std::getline(_csvInput, _Oneline))  // 整行读取，换行符“\n”区分，遇到文件尾标志eof终止读取
    {
        std::istringstream _Readstr(_Oneline); // 将整行字符串_Oneline读入到字符串流istringstream中
        std::string _partOfstr;
        double pointCorr;
        int i = 0;
        while (std::getline(_Readstr, _partOfstr, ',')) // 将字符串流_Readstr中的字符读入到_partOfstr字符串中，以逗号为分隔符
        {
            if (i < 3)
            {
                pointCorr = stringToNum<double>(_partOfstr);
                _lineOfstr.push_back(pointCorr); // 将刚刚读取的字符串添加到向量_lineOfstr中,
            }
            else if (i == 3)
            {
                verticalLinePoints.push_back(_lineOfstr);
                _lineOfstr.clear();
                _Readstr.clear();
                break;
            }
            ++i;
        }
    }
    pointsVec = verticalLinePoints;
}


// **************************************************************************************************************

// 功能：计算公垂线L与两直线（L0、L1）交点

// 输入：dir_V:竖向直线（L0）方向向量； center_V：竖向3D点质心
//       dir_H:横向直线（L1）方向向量； center_H：横向3D点质心

// 输出：p0: L0与公垂线（L）的交点； p1: L1与公垂线（L）的交点；

// return: -1: L0、L1平行（重合）
//          0: L0、L1不平行（不重合）

// **************************************************************************************************************
int Intersection3DPoint(Eigen::Matrix3d dir_V, vector<double> center_V, Eigen::Matrix3d dir_H, vector<double> center_H, vector<double>& p0, vector<double>& p1)
{
    //直线1： (x - x0)/a0 = (y - y0)/b0 = (z - z0)/c0
    double a0 = dir_V(0, 0), b0 = dir_V(1, 0), c0 = dir_V(2, 0); // (a0 b0 c0) 是竖向直线L0的方向向量
    double x0 = center_V[0], y0 = center_V[1], z0 = center_V[2]; // (x0 y0 z0) 是竖向3D样本点的质心
    //直线2： (x - x1)/a1 = (y - y1)/b1 = (z - z1)/c1
    double a1 = dir_H(0, 0), b1 = dir_H(1, 0), c1 = dir_H(2, 0);
    double x1 = center_H[0], y1 = center_H[1], z1 = center_H[2];
    Eigen::Matrix<double, 1, 3> dir_V_, dir_H_;

    // 方向向量零分量补偿
    if (a0 == 0) // abs(a0 - 1e-5)
        a0 += 1e-5;
    if (b0 == 0)
        b0 += 1e-5;
    if (c0 == 0)
        c0 += 1e-5;

    if (a1 == 0)
        a1 += 1e-5;
    if (b1 == 0)
        b1 += 1e-5;
    if (c1 == 0)
        c1 += 1e-5;

    //cout << a0 / a1 << endl;
    //cout << b0 / b1 << endl;
    //cout << c0 / c1 << endl;
    // 判断平行(重合)
    if (a0 / a1 == b0 / b1 == c0 / c1)
        return -1;

    dir_V_(0, 0) = a0;
    dir_V_(0, 1) = b0;
    dir_V_(0, 2) = c0;
    dir_H_(0, 0) = a1;
    dir_H_(0, 1) = b1;
    dir_H_(0, 2) = c1;
    // <1>、求公垂线方向
    Eigen::MatrixXd crossM = dir_H_.cross(dir_V_); // (a0 b0 c0) 与 (a1 b1 c1) 作叉积 -> 公垂向量(a2 b2 c2)
    double a2 = crossM(0, 0);
    double b2 = crossM(0, 1);
    double c2 = crossM(0, 2);
    // <2>、求交点p0、p1(推导见 readme.txt)
    Eigen::Matrix<double, 2, 2> A;
    A(0, 0) = a1 * b2 - a2 * b1;
    A(0, 1) = a2 * b0 - a0 * b2;
    A(1, 0) = b1 * c2 - b2 * c1;
    A(1, 1) = b2 * c0 - b0 * c2;
    Eigen::Matrix<double, 2, 1> B;
    B(0, 0) = a2 * (y1 - y0) + b2 * (x0 - x1);
    B(1, 0) = b2 * (z1 - z0) + c2 * (y0 - y1);
    Eigen::Matrix<double, 2, 1> X;
    X = A.inverse() * B; // 解方程 A*X = B
    double alpha = X(1, 0);
    double beta = X(0, 0);
    /*p0.resize(3);
    p1.resize(3);*/
    p0[0] = alpha * a0 + x0;
    p0[1] = alpha * b0 + y0;
    p0[2] = alpha * c0 + z0;
    p1[0] = beta * a1 + x1;
    p1[1] = beta * b1 + y1;
    p1[2] = beta * c1 + z1;
    return 0;
}