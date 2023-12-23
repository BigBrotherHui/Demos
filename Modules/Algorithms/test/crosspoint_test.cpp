#include"../CrossPoint.h"
#include <iostream>
using namespace std;
/**************************************************************************************

//function: 2D平面上，求解两条直线交点

//input:    line1: 直线方程:a0 * x + b0 * y + c0 = 0
            line2: 直线方程:a1 * x + b1 * y + c1 = 0

 //output:   point: 交点

 //return:   -1: line1平行于line2
             0: 有交点值返回

**************************************************************************************/

int CalculateIntersectionPoint(const Eigen::Vector3d & line1, const Eigen::Vector3d& line2, Eigen::Vector2d &point)
{
	// ...
	double a0 = line1[0];
	double b0 = line1[1];
	double c0 = line1[2];

	double a1 = line2[0];
	double b1 = line2[1];
	double c1 = line2[2];

	double D = a0 * b1 - a1 * b0; // ...
	if (abs(D) < 1e-5)
		return -1;
	point[0] = (b0 * c1 - b1 * c0) / D;
	point[1] = (a1 * c0 - a0 * c1) / D;
	return 0;
}

//int main()
//{
//	cv::Vec3d line1(0, 1, -4);
//	cv::Vec3d line2(1, -2, 2);
//	cv::Point2d point;
//	if (CalculateIntersectionPoint(line1, line2, point) == 0)
//	    {
//	        cout << point.x << "," << point.y;
//	    }
//	return 1;
//}
int test2(Eigen::MatrixXd& pointsMatrix_H,
    Eigen::MatrixXd& pointsMatrix_V,
    vector<double>& p0,
    vector<double>& p1)
{
    // <1>
    Eigen::Matrix3d V_H;
    double centerX_H = 0.0, centerY_H = 0.0, centerZ_H = 0.0;
    spacialLineFitting(pointsMatrix_H, V_H, centerX_H, centerY_H, centerZ_H); // 拟合
    //cout << V_H(0, 0) << V_H(1, 0)  << V_H(2, 0) << endl; 方向
    vector<double> center_H; // 质心
    center_H.reserve(3);
    center_H.push_back(centerX_H);
    center_H.push_back(centerY_H);
    center_H.push_back(centerZ_H);
    // <2>
    Eigen::Matrix3d V_V;
    double centerX_V = 0.0, centerY_V = 0.0, centerZ_V = 0.0;
    spacialLineFitting(pointsMatrix_V, V_V, centerX_V, centerY_V, centerZ_V);
    //cout << V_V(0, 0) << V_V(1, 0) << V_V(2, 0) << endl;
    vector<double> center_V;
    center_V.reserve(3);
    center_V.push_back(centerX_V);
    center_V.push_back(centerY_V);
    center_V.push_back(centerZ_V);
    // <3>
    return Intersection3DPoint(V_V, center_V, V_H, center_H, p0, p1);
}

int t()
{
    cout << "p0、p1分别是公垂线L与直线L0、L1的交点" << endl;

    vector<double> p0, p1;//与两直线交点
    p0.resize(3);
    p1.resize(3);

    // 情况1:实测数据
    /*if (test1(p0, p1) == 0)
    {
        cout << "p0(" << p0[0] << "," << p0[1] << "," << p0[2] << ")" << endl;
        cout << "p1(" << p1[0] << "," << p1[1] << "," << p1[2] << ")" << endl;
        cout << "************************" << endl;
    }*/


    // 情况2:方向向量有零分量
    Eigen::MatrixXd pointsMatrix_H(3, 3); // 3个3D点，水平
    pointsMatrix_H << 0, 0, 1, 2, 0, 1, 3, 0, 1; // (0, 0, 1), (2, 0, 1), (3, 0, 1)
    Eigen::MatrixXd pointsMatrix_V(3, 3);
    pointsMatrix_V << 0, 0, 4, 0, 2, 4, 0, 9, 4; // (0, 0, 4), (0, 2, 4), (0, 9, 4)

    if (test2(pointsMatrix_H, pointsMatrix_V, p0, p1) == 0)
    {
        cout << "p0(" << p0[0] << "," << p0[1] << "," << p0[2] << ")" << endl; //  理论：(0 0 4)
        cout << "p1(" << p1[0] << "," << p1[1] << "," << p1[2] << ")" << endl; //  理论：(0 0 1)
        cout << "************************" << endl;
    }


    // 情况3:两直线有交点
    pointsMatrix_H << 0, 0, 0, 0, 0, 2, 0, 0, 7; // (0, 0, 0), (0, 0, 2), (0, 0, 7)
    pointsMatrix_V << 0, 0, 3, 3, 0, 3, 4, 0, 3; // (0, 0, 3), (3, 0, 3), (4, 0, 3)

    if (test2(pointsMatrix_H, pointsMatrix_V, p0, p1) == 0)
    {
        cout << "p0(" << p0[0] << "," << p0[1] << "," << p0[2] << ")" << endl; // 理论：(0 0 3)
        cout << "p1(" << p1[0] << "," << p1[1] << "," << p1[2] << ")" << endl; // 理论：(0 0 3)
        cout << "************************" << endl;
    }

    // 情况4:平行(重合)(既平行，又是零向量，所以要先判断是否平行)
    pointsMatrix_H << 1.33, 0, 3.21, 1.33, 0, 3.31, 1.33, 0, 4.81; // (1.33, 0, 3.21), (1.33, 0, 3.31), (1.33, 0, 4.81)
    pointsMatrix_V << 0, 2.3, 3.1, 0, 2.3, 3.2, 0, 2.3, 1.222; // (0, 2.3, 3.1), (0, 2.3, 3.2), (0, 2.3, 1.222)

    if (test2(pointsMatrix_H, pointsMatrix_V, p0, p1) == 0)
    {
        cout << "p0(" << p0[0] << "," << p0[1] << "," << p0[2] << ")" << endl; // 理论：由于平行，直接return
        cout << "p1(" << p1[0] << "," << p1[1] << "," << p1[2] << ")" << endl; //
        cout << "************************" << endl;
    }




    system("pause");
    return 0;
}