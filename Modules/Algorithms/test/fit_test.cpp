// FitSpaceLineRansacTest.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include"../FitSpaceLineRansac.h"

int main()
{
    Eigen::MatrixXd pointsMatrix(12, 3); // 3D点数量最少四个点
    pointsMatrix <<
        // 内点
        0, 0, 0,
        1, 1, 1,
        2, 2, 2,
        3, 3, 3,
        4, 4, 4,
        5, 5, 5,
        6, 6, 6,
        7, 7, 7,
        //外点
        3, 4, 5,
        4, 5, 6,
        5, 4, 3,
        4, 4, 7;

    FitSpaceLineRansac fitLine(pointsMatrix, 0.0000005, 50, 0.01);
    fitLine.compute();
    cout << fitLine.getModel() << endl;

    // TODO
    return -1;
}