// FitSpaceLineRansacTest.cpp : ���ļ����� "main" ����������ִ�н��ڴ˴���ʼ��������
//

#include"../FitSpaceLineRansac.h"

int main()
{
    Eigen::MatrixXd pointsMatrix(12, 3); // 3D�����������ĸ���
    pointsMatrix <<
        // �ڵ�
        0, 0, 0,
        1, 1, 1,
        2, 2, 2,
        3, 3, 3,
        4, 4, 4,
        5, 5, 5,
        6, 6, 6,
        7, 7, 7,
        //���
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