#include"CrossPoint.h"

// **************************************************************************************************************

// ���ܣ��ռ�ֱ�����

// ���룺pointsMatrix��3D��������

// �����V��ֱ�߷���������(centerX centerY centerZ):3D�����㼯����

// return: void

// **************************************************************************************************************
void spacialLineFitting(Eigen::MatrixXd& pointsMatrix, Eigen::Matrix3d& V, double& centerX, double& centerY, double& centerZ)
{
    int r = pointsMatrix.rows();
    Eigen::MatrixXd sumMatrixXd(1, 3);
    sumMatrixXd = pointsMatrix.colwise().sum();//��x��y��z��
    // <1>��������
    centerX = sumMatrixXd(0, 0) / r;
    centerY = sumMatrixXd(0, 1) / r;
    centerZ = sumMatrixXd(0, 2) / r;
    Eigen::MatrixXd centerdMatrix(r, 3);
    // <2>��ȥ����
    for (int i = 0; i < r; i++)
    {
        centerdMatrix.row(i)[0] = pointsMatrix.row(i)[0] - centerX;
        centerdMatrix.row(i)[1] = pointsMatrix.row(i)[1] - centerY;
        centerdMatrix.row(i)[2] = pointsMatrix.row(i)[2] - centerZ;
    }
    // <3>��SVD�ֽ�
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(centerdMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    V = svd.matrixV();
}
//�ַ�תdouble
template <class Type>
Type stringToNum(string& str)
{
    istringstream iss(str);
    Type num;
    iss >> num;
    return num;
}
//��ȡcsv�ļ���
void readPoint(string strPath, vector<vector<double>>& pointsVec)
{
    std::ifstream _csvInput;
    _csvInput.open(strPath, std::ios::in);
    std::string _Oneline;
    std::vector<double> _lineOfstr;
    vector<vector<double>> verticalLinePoints;
    while (std::getline(_csvInput, _Oneline))  // ���ж�ȡ�����з���\n�����֣������ļ�β��־eof��ֹ��ȡ
    {
        std::istringstream _Readstr(_Oneline); // �������ַ���_Oneline���뵽�ַ�����istringstream��
        std::string _partOfstr;
        double pointCorr;
        int i = 0;
        while (std::getline(_Readstr, _partOfstr, ',')) // ���ַ�����_Readstr�е��ַ����뵽_partOfstr�ַ����У��Զ���Ϊ�ָ���
        {
            if (i < 3)
            {
                pointCorr = stringToNum<double>(_partOfstr);
                _lineOfstr.push_back(pointCorr); // ���ոն�ȡ���ַ�����ӵ�����_lineOfstr��,
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

// ���ܣ����㹫����L����ֱ�ߣ�L0��L1������

// ���룺dir_V:����ֱ�ߣ�L0������������ center_V������3D������
//       dir_H:����ֱ�ߣ�L1������������ center_H������3D������

// �����p0: L0�빫���ߣ�L���Ľ��㣻 p1: L1�빫���ߣ�L���Ľ��㣻

// return: -1: L0��L1ƽ�У��غϣ�
//          0: L0��L1��ƽ�У����غϣ�

// **************************************************************************************************************
int Intersection3DPoint(Eigen::Matrix3d dir_V, vector<double> center_V, Eigen::Matrix3d dir_H, vector<double> center_H, vector<double>& p0, vector<double>& p1)
{
    //ֱ��1�� (x - x0)/a0 = (y - y0)/b0 = (z - z0)/c0
    double a0 = dir_V(0, 0), b0 = dir_V(1, 0), c0 = dir_V(2, 0); // (a0 b0 c0) ������ֱ��L0�ķ�������
    double x0 = center_V[0], y0 = center_V[1], z0 = center_V[2]; // (x0 y0 z0) ������3D�����������
    //ֱ��2�� (x - x1)/a1 = (y - y1)/b1 = (z - z1)/c1
    double a1 = dir_H(0, 0), b1 = dir_H(1, 0), c1 = dir_H(2, 0);
    double x1 = center_H[0], y1 = center_H[1], z1 = center_H[2];
    Eigen::Matrix<double, 1, 3> dir_V_, dir_H_;

    // �����������������
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
    // �ж�ƽ��(�غ�)
    if (a0 / a1 == b0 / b1 == c0 / c1)
        return -1;

    dir_V_(0, 0) = a0;
    dir_V_(0, 1) = b0;
    dir_V_(0, 2) = c0;
    dir_H_(0, 0) = a1;
    dir_H_(0, 1) = b1;
    dir_H_(0, 2) = c1;
    // <1>���󹫴��߷���
    Eigen::MatrixXd crossM = dir_H_.cross(dir_V_); // (a0 b0 c0) �� (a1 b1 c1) ����� -> ��������(a2 b2 c2)
    double a2 = crossM(0, 0);
    double b2 = crossM(0, 1);
    double c2 = crossM(0, 2);
    // <2>���󽻵�p0��p1(�Ƶ��� readme.txt)
    Eigen::Matrix<double, 2, 2> A;
    A(0, 0) = a1 * b2 - a2 * b1;
    A(0, 1) = a2 * b0 - a0 * b2;
    A(1, 0) = b1 * c2 - b2 * c1;
    A(1, 1) = b2 * c0 - b0 * c2;
    Eigen::Matrix<double, 2, 1> B;
    B(0, 0) = a2 * (y1 - y0) + b2 * (x0 - x1);
    B(1, 0) = b2 * (z1 - z0) + c2 * (y0 - y1);
    Eigen::Matrix<double, 2, 1> X;
    X = A.inverse() * B; // �ⷽ�� A*X = B
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