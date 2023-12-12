#include "BasicAlgorithm.h"

static const double M_PI=3.14159265358979323846;

Eigen::Quaterniond BasicAlgorithm::Quaternion_Slerp(Eigen::Quaterniond& start_q, Eigen::Quaterniond& end_q, double t)
{
    //Eigen::Quaterniond q = q0.slerp(t0,q1);可以使用Eigen默认的slerp方法
    Eigen::Quaterniond slerp_q;

    double cos_angle = start_q.x() * end_q.x()
        + start_q.y() * end_q.y()
        + start_q.z() * end_q.z()
        + start_q.w() * end_q.w();

    // 如果四元数点积的结果是负值（夹角大于90°），那么后面的插值就会在4D球面上绕远路。为了解决这个问题，
	//先测试点积的结果，当结果是负值时，将2个四元数的其中一个取反（并不会改变它代表的朝向）。
	//而经过这一步操作，可以保证这个旋转走的是最短路径。四元数的三个虚变量对应旋转轴，
	//因此两个相反向量对应的是同一旋转轴，因此两个相反的四元数表示的旋转是一样的。
    if (cos_angle < 0) {
        end_q.x() = -end_q.x();
        end_q.y() = -end_q.y();
        end_q.z() = -end_q.z();
        end_q.w() = -end_q.w();
        cos_angle = -cos_angle;
    }

    double ratio_A, ratio_B;
    //当p和q的夹角差非常小时会导致sinθ → 0，这时除法可能会出现问题。为了避免这样的问题，
	//当θ非常小时可以使用简单的线性插值代替（θ → 0 时，sinθ ≈ θ，
	//因此方程退化为线性方程：slerp(p,q,t)=(1−t)p+tq)slerp(p,q,t)=(1-t)p+tq)slerp(p,q,t)=(1−t)p+tq)
    if (cos_angle > 0.99995f) {
        ratio_A = 1.0f - t;
        ratio_B = t;
    }
    else {
        double sin_angle = sqrt(1.0f - cos_angle * cos_angle);
        double angle = atan2(sin_angle, cos_angle);
        ratio_A = sin((1.0f - t) * angle) / sin_angle;
        ratio_B = sin(t * angle) / sin_angle;
    }

    slerp_q.x() = ratio_A * start_q.x() + ratio_B * end_q.x();
    slerp_q.y() = ratio_A * start_q.y() + ratio_B * end_q.y();
    slerp_q.z() = ratio_A * start_q.z() + ratio_B * end_q.z();
    slerp_q.w() = ratio_A * start_q.w() + ratio_B * end_q.w();

    return slerp_q.normalized();
}

Eigen::Matrix4d BasicAlgorithm::RotateAroundAxisAtFixedPoint(const Eigen::Vector3d& point, const Eigen::Vector3d& dir, double angle)
{
    // 平移到指定点p，绕轴旋转，再平移回原点
    Eigen::Vector3d normdir = dir.normalized();

    // 角度转为弧度
    angle = angle * M_PI / 180.0;

    Eigen::Affine3d A;
    A = Eigen::Translation3d(point) * Eigen::AngleAxisd(angle, normdir) * Eigen::Translation3d(-point);

    Eigen::Matrix4d mat4 = A.matrix();

    return mat4;
}

Eigen::Matrix3d BasicAlgorithm::GenerateRotationMatrix(const Eigen::Vector3d& vectorBefore, const Eigen::Vector3d& vectorAfter)
{
    Eigen::Vector3d vector = vectorBefore.cross(vectorAfter);
    double angle = AngleBetween2Vector(vectorBefore, vectorAfter, true);
    Eigen::AngleAxisd rotationVetor(angle, vector.normalized());
    return rotationVetor.toRotationMatrix();
}

double BasicAlgorithm::GetPointToLineDistance(double* point, double* linePoint_0, double* linePoint_1)
{
	Eigen::Vector3d lineVector;
	lineVector[0] = linePoint_0[0] - linePoint_1[0];
	lineVector[1] = linePoint_0[1] - linePoint_1[1];
	lineVector[2] = linePoint_0[2] - linePoint_1[2];

	Eigen::Vector3d tmpVector;
	tmpVector[0] = linePoint_0[0] - point[0];
	tmpVector[1] = linePoint_0[1] - point[1];
	tmpVector[2] = linePoint_0[2] - point[2];

	double projection = lineVector.dot(tmpVector) / lineVector.norm();

	double distance = sqrt(pow(tmpVector.norm(), 2) - pow(projection, 2));

	return distance;
}

Eigen::Matrix4d BasicAlgorithm::ConvertCoordstoTransform(Eigen::Vector3d& o, Eigen::Vector3d& x, Eigen::Vector3d& y,
	Eigen::Vector3d& z)
{
	Eigen::Matrix4d transform;
	transform <<
		x(0), y(0), z(0), o(0),
		x(1), y(1), z(1), o(1),
		x(2), y(2), z(2), o(2),
		0, 0, 0, 1;
	return transform;
}

void BasicAlgorithm::ProjectToPlane(const double x[3], const double origin[3], const double normal[3], double xproj[3])
{
	double t, xo[3];

	xo[0] = x[0] - origin[0];
	xo[1] = x[1] - origin[1];
	xo[2] = x[2] - origin[2];

	Eigen::Vector3d nl(normal);
	nl.normalize();
	Eigen::Vector3d xv(xo);

	t = xv.dot(nl);

	xproj[0] = x[0] - t * nl[0];
	xproj[1] = x[1] - t * nl[1];
	xproj[2] = x[2] - t * nl[2];
}

double BasicAlgorithm::AngleBetween2Vector(const double vec1[3], const double vec2[3], bool radian)
{
	Eigen::Vector3d v1(vec1);
	Eigen::Vector3d v2(vec2);
	auto cross = v1.cross(v2);

	if (radian)
	{
		return atan2(cross.norm(), v1.dot(v2));
	}
	else
	{
		auto DegreesFromRadians = [](double x) -> double { return x * 57.2957795131; };
		double degreeAngle = DegreesFromRadians(atan2(cross.norm(), v1.dot(v2)));
		return degreeAngle;
	}
}

double BasicAlgorithm::DistanceBetweenTwoPoints(const double p1[3], const double p2[3])
{
	return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
}

double BasicAlgorithm::DistanceFromPointToLineWithDirectionVector(const double PointOutsideLine[3],
	const double PointOnLine[3], const double direction[3])
{
	Eigen::Vector3d xp1, p1p2;
	double proj;

	for (int i = 0; i < 3; i++)
	{
		xp1[i] = PointOutsideLine[i] - PointOnLine[i];
		p1p2[i] = direction[i];
	}

	if (p1p2.norm() != 0.0)
	{
		p1p2.normalize();
	}
	else
	{
		return xp1.norm();
	}

	//calculate the dot product
	proj = xp1.dot(p1p2);
	return sqrt(xp1.dot(xp1) - proj * proj);
}

double BasicAlgorithm::DistanceFromPointToPlane(const double x[3], const double p0[3], const double n[3])
{
	return abs(n[0] * (x[0] - p0[0]) + n[1] * (x[1] - p0[1]) + n[2] * (x[2] - p0[2])) / sqrt(pow(n[0], 2) + pow(n[1], 2) + pow(n[2], 2));
}

double BasicAlgorithm::AngleBetween2Vector(const Eigen::Vector3d vec1, const Eigen::Vector3d vec2, bool radian)
{
	Eigen::Vector3d v1(vec1);
	Eigen::Vector3d v2(vec2);
	auto cross = v1.cross(v2);

	if (radian)
	{
		return atan2(cross.norm(), v1.dot(v2));
	}
	else
	{
		auto DegreesFromRadians = [](double x) -> double { return x * 57.2957795131; };
		double degreeAngle = DegreesFromRadians(atan2(cross.norm(), v1.dot(v2)));
		return degreeAngle;
	}
}

double BasicAlgorithm::AngleBetween2Vector(const double vec1[3], const double vec2[3], const double normal[3],
	bool radian)
{
	Eigen::Vector3d v1(vec1);
	Eigen::Vector3d v2(vec2);
	Eigen::Vector3d norm(normal);
	auto cross = v1.cross(v2);

	if (radian)
	{
		if (cross.dot(norm) >= 0)
		{
			return atan2(cross.norm(), v1.dot(v2));
		}
		else
		{
			return -atan2(cross.norm(), v1.dot(v2));
		}
	}
	else
	{
		if (cross.dot(norm) >= 0)
		{
			auto DegreesFromRadians = [](double x) -> double { return x * 57.2957795131; };
			double degreeAngle = DegreesFromRadians(atan2(cross.norm(), v1.dot(v2)));
			return degreeAngle;
		}
		else
		{
			auto DegreesFromRadians = [](double x) -> double { return x * 57.2957795131; };
			double degreeAngle = DegreesFromRadians(-atan2(cross.norm(), v1.dot(v2)));
			return degreeAngle;
		}
	}
}

double BasicAlgorithm::AngleBetween2Vector(const Eigen::Vector3d vec1, const Eigen::Vector3d vec2,
	const Eigen::Vector3d normal, bool radian)
{
	auto cross = vec1.cross(vec2);

	if (radian)
	{
		if (cross.dot(normal) >= 0)
		{
			return atan2(cross.norm(), vec1.dot(vec2));
		}
		else
		{
			return -atan2(cross.norm(), vec1.dot(vec2));
		}
	}
	else
	{
		if (cross.dot(normal) >= 0)
		{
			auto DegreesFromRadians = [](double x) -> double { return x * 57.2957795131; };
			double degreeAngle = DegreesFromRadians(atan2(cross.norm(), vec1.dot(vec2)));
			return degreeAngle;
		}
		else
		{
			auto DegreesFromRadians = [](double x) -> double { return x * 57.2957795131; };
			double degreeAngle = DegreesFromRadians(-atan2(cross.norm(), vec1.dot(vec2)));
			return degreeAngle;
		}
	}
}

double BasicAlgorithm::AngleBetween2Line(const double* p11, const double* p12, const double* p21, const double* p22,
	bool radian)
{
	double vec1[3] = { p12[0] - p11[0], p12[1] - p11[1], p12[2] - p11[2] };
	double vec2[3] = { p22[0] - p21[0], p22[1] - p21[1], p22[2] - p21[2] };

	double angle = AngleBetween2Vector(vec1, vec2, radian);

	if (angle > 90)
	{
		return 180 - angle;
	}
	else
	{
		return angle;
	}
}

double BasicAlgorithm::AngleBetweenLineAndPlane(const double* p1, const double* p2, const double* normal, bool radian)
{
	auto RadiansFromDegrees = [](double x) -> double { return x * 0.017453292; };
	Eigen::Vector3d LineVector(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
	Eigen::Vector3d PlaneNormal(normal);

	//平行的时候，向量方向相同则为0， 相反则为180
	double degreeAngle = AngleBetween2Vector(LineVector, PlaneNormal, false);
	double result;
	double diff = fabs(degreeAngle - 90);

	if (diff < 1e-6)
	{
		//直线与平面平行
		return 0;
	}
	else
	{
		if (degreeAngle > 90)
		{
			result = degreeAngle - 90;
		}
		else
		{
			result = 90 - degreeAngle;
		}
	}

	if (radian)
	{
		return RadiansFromDegrees(result);
	}

	return result;
}

double BasicAlgorithm::AngleBetweenLineAndPlane(const Eigen::Vector3d vec, const Eigen::Vector3d normal, bool radian)
{
	auto RadiansFromDegrees = [](double x) -> double { return x * 0.017453292; };

	//平行的时候，向量方向相同则为0， 相反则为180
	double degreeAngle = AngleBetween2Vector(vec, normal, false);
	double result;
	double diff = fabs(degreeAngle - 90);

	if (diff < 1e-6)
	{
		//直线与平面平行
		return 0;
	}
	else
	{
		if (degreeAngle > 90)
		{
			result = degreeAngle - 90;
		}
		else
		{
			result = 90 - degreeAngle;
		}
	}

	if (radian)
	{
		return RadiansFromDegrees(result);
	}

	return result;
}
void BasicAlgorithm::fit_circle_2d(std::vector<double>& inp_x, std::vector<double>& inp_y, double& outp_Cx,
    double& outp_Cy, double& outp_R)
{
    if (inp_x.size() != inp_y.size())
    {
        std::cout << "fit_circle err: input x y must have same length" << std::endl;
        return;
    }

    //建立矩阵
    int length = inp_x.size();
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::Vector3d C;
    A.resize(length, 3);
    b.resize(length);
    C.resize(3);

    //input A b
    for (int i = 0; i < length; i++)
    {
        A(i, 0) = inp_x[i];
        A(i, 1) = inp_y[i];
        A(i, 2) = 1;

        b(i) = pow(inp_x[i], 2) + pow(inp_y[i], 2);
    }

    C = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    outp_Cx = C[0] / 2;
    outp_Cy = C[1] / 2;
    outp_R = sqrt(C[2] + pow(outp_Cx, 2) + pow(outp_Cy, 2));
}
bool BasicAlgorithm::fit_circle_3d(std::vector<double>& inp_pointset, std::array<double, 3>& outp_center, double& outp_radius,
    std::array<double, 3>& outp_normal)
{
    auto len = inp_pointset.size() / 3;

    if (len < 3)
    {
        return false;
    }

    //Store input point set into a 3xN eigen Matrix, N is the size of input point set
    Eigen::MatrixXd pointSet;
    pointSet.resize(3, len);
    pointSet.setZero();

    for (int i = 0; i < len; i++)
    {
        pointSet(0, i) = inp_pointset[3 * i + 0];
        pointSet(1, i) = inp_pointset[3 * i + 1];
        pointSet(2, i) = inp_pointset[3 * i + 2];
    }

    //Centralize data
    auto p_mean = pointSet.rowwise().mean().eval();
    pointSet = pointSet.colwise() - p_mean;


    //Using SVD to get local coordinate, Which can describe the attitude of the point set
    //Rg2l = U
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(pointSet, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::Matrix3d Rg2l = svd.matrixU();
    //auto V = svd.matrixV();
    //auto A = svd.singularValues();


    /*Eigen::Vector3d local_x = U.col(0);
    Eigen::Vector3d local_y = U.col(1);
    Eigen::Vector3d local_z = U.col(2);

    Eigen::Matrix3d Rg2l;
    Rg2l << local_x[0], local_y[0], local_z[0],
            local_x[1], local_y[1], local_z[1],
            local_x[2], local_y[2], local_z[2];*/


            //calculate point position in local coordinates; P_local = Rl2g * P_global = Rg2l.transpose() * P_global
    pointSet = Rg2l.transpose().eval() * pointSet;


    //Using the points in local coordinates xy plane to fit the circle, ignore the z-axis information
    std::vector<double> x;
    std::vector<double> y;

    for (int i = 0; i < pointSet.cols(); i++)
    {
        x.push_back(pointSet.col(i)[0]);
        y.push_back(pointSet.col(i)[1]);
    }

    double cx{ 0 }, cy{ 0 }, r{ 0 };
    fit_circle_2d(x, y, cx, cy, r);

    //Convert the circle center back to the global coordinates;
    Eigen::Vector3d center;
    center << cx, cy, 0;
    center = Rg2l * center;

    center = center + p_mean;
    //output
    outp_center[0] = center[0];
    outp_center[1] = center[1];
    outp_center[2] = center[2];

    outp_radius = r;

    outp_normal[0] = Rg2l(0, 2);
    outp_normal[1] = Rg2l(1, 2);
    outp_normal[2] = Rg2l(2, 2);

    return true;
}

bool BasicAlgorithm::fit_sphere(std::vector<double>& inp_x, std::vector<double>& inp_y, std::vector<double>& inp_z,
    double& outp_cx, double& outp_cy, double& outp_cz, double& outp_R)
{
    if (inp_x.size() != inp_y.size() || inp_x.size() != inp_z.size() || inp_y.size() != inp_z.size())
    {
        std::cout << "fit_circle err: input x y must have same length" << std::endl;
        return false;
    }

    //建立矩阵
    int length = inp_x.size();
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::Vector4d C;
    A.resize(length, 4);
    b.resize(length);
    C.resize(4);

    //input A b
    for (int i = 0; i < length; i++)
    {
        A(i, 0) = inp_x[i];
        A(i, 1) = inp_y[i];
        A(i, 2) = inp_z[i];
        A(i, 3) = 1;

        b(i) = pow(inp_x[i], 2) + pow(inp_y[i], 2) + pow(inp_z[i], 2);
    }

    C = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    outp_cx = C[0] / 2;
    outp_cy = C[1] / 2;
    outp_cz = C[2] / 2;
    outp_R = sqrt(C[3] + pow(outp_cx, 2) + pow(outp_cy, 2) + pow(outp_cz, 2));

    return true;
}
bool BasicAlgorithm::fit_sphere(std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, double& outp_r)
{
    auto length = inp_pSet.size() / 3;

    if (inp_pSet.size() % 3 != 0)
    {
        std::cout << "fit_sphere err: input point set must be Multiple of 3" << std::endl;
        return false;
    }

    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::Vector4d C;
    A.resize(length, 4);
    b.resize(length);
    C.resize(4);

    //input A b
    for (int i = 0; i < length; i++)
    {
        A(i, 0) = inp_pSet[0 + 3 * i];
        A(i, 1) = inp_pSet[1 + 3 * i];
        A(i, 2) = inp_pSet[2 + 3 * i];
        A(i, 3) = 1;

        b(i) = pow(inp_pSet[0 + 3 * i], 2) + pow(inp_pSet[1 + 3 * i], 2) + pow(inp_pSet[2 + 3 * i], 2);
    }

    C = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    outp_center[0] = C[0] / 2;
    outp_center[1] = C[1] / 2;
    outp_center[2] = C[2] / 2;

    outp_r = sqrt(C[3] + pow(outp_center[0], 2) + pow(outp_center[1], 2) + pow(outp_center[2], 2));

    return true;
}
bool BasicAlgorithm::fit_sphere(std::vector<std::array<double, 3>>& inp_pSet, std::array<double, 3>& outp_center, double& outp_r)
{
    auto length = inp_pSet.size();

    if (inp_pSet.size() < 4)
    {
        std::cout << "fit_sphere err: not enough points" << std::endl;
        return false;
    }

    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::Vector4d C;
    A.resize(length, 4);
    b.resize(length);
    C.resize(4);

    //input A b
    for (int i = 0; i < length; i++)
    {
        A(i, 0) = inp_pSet[i][0];
        A(i, 1) = inp_pSet[i][1];
        A(i, 2) = inp_pSet[i][2];
        A(i, 3) = 1;

        b(i) = pow(inp_pSet[i][0], 2) + pow(inp_pSet[i][1], 2) + pow(inp_pSet[i][2], 2);
    }

    C = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    outp_center[0] = C[0] / 2;
    outp_center[1] = C[1] / 2;
    outp_center[2] = C[2] / 2;

    outp_r = sqrt(C[3] + pow(outp_center[0], 2) + pow(outp_center[1], 2) + pow(outp_center[2], 2));

    return true;
}

bool BasicAlgorithm::fit_sphere_fixR(const std::vector<double>& inp_x, const std::vector<double>& inp_y, const std::vector<double>& inp_z, const double inp_r,
    double& outp_cx, double& outp_cy, double& outp_cz)
{
    if (inp_x.size() != inp_y.size() || inp_x.size() != inp_z.size() || inp_y.size() != inp_z.size())
    {
        std::cout << "fit_circle err: input x y must have same length" << std::endl;
        return false;
    }

    int length = inp_x.size();
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::Vector4d C;
    A.resize(length, 4);
    b.resize(length);
    C.resize(4);

    //input A b
    for (int i = 0; i < length; i++)
    {
        A(i, 0) = inp_x[i];
        A(i, 1) = inp_y[i];
        A(i, 2) = inp_z[i];
        A(i, 3) = 1;

        b(i) = pow(inp_x[i], 2) + pow(inp_y[i], 2) + pow(inp_z[i], 2) - pow(inp_r, 2);
    }

    C = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    outp_cx = C[0] / 2;
    outp_cy = C[1] / 2;
    outp_cz = C[2] / 2;

    return true;
}

bool BasicAlgorithm::fit_plane(const std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, std::array<double, 3>& outp_normal)
{
    auto len = inp_pSet.size() / 3;

    if (len < 3)
    {
        return false;
    }

    //Store input point set into a 3xN eigen Matrix, N is the size of input point set
    Eigen::MatrixXd pointSet;
    pointSet.resize(3, len);
    pointSet.setZero();

    for (int i = 0; i < len; i++)
    {
        pointSet(0, i) = inp_pSet[3 * i + 0];
        pointSet(1, i) = inp_pSet[3 * i + 1];
        pointSet(2, i) = inp_pSet[3 * i + 2];
    }

    //Centralize data
    auto p_mean = pointSet.rowwise().mean().eval();
    pointSet = pointSet.colwise() - p_mean;


    //Using SVD to get local coordinate, Which can describe the attitude of the point set
    //Rg2l = U
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(pointSet, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::Matrix3d Rg2l = svd.matrixU();

    outp_center[0] = p_mean[0];
    outp_center[1] = p_mean[1];
    outp_center[2] = p_mean[2];

    outp_normal[0] = Rg2l(0, 2);
    outp_normal[1] = Rg2l(1, 2);
    outp_normal[2] = Rg2l(2, 2);

    return true;
}
bool BasicAlgorithm::fit_rectangle(const std::vector<double>& inp_pSet,
    std::array<double, 3>& outp_center,
    std::array<double, 3>& outp_normal,
    std::array<double, 3>& outp_x,
    std::array<double, 3>& outp_y,
    double& length,
    double& width)
{
    auto len = inp_pSet.size() / 3;

    if (len < 3)
    {
        return false;
    }

    //Store input point set into a 3xN eigen Matrix, N is the size of input point set
    Eigen::MatrixXd pointSet;
    pointSet.resize(3, len);
    pointSet.setZero();

    for (int i = 0; i < len; i++)
    {
        pointSet(0, i) = inp_pSet[3 * i + 0];
        pointSet(1, i) = inp_pSet[3 * i + 1];
        pointSet(2, i) = inp_pSet[3 * i + 2];
    }

    //Centralize data
    auto p_mean = pointSet.rowwise().mean().eval();
    pointSet = pointSet.colwise() - p_mean;


    //Using SVD to get local coordinate, Which can describe the attitude of the point set
    //Rg2l = U
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(pointSet, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::Matrix3d Rg2l = svd.matrixU();

    outp_center[0] = p_mean[0];
    outp_center[1] = p_mean[1];
    outp_center[2] = p_mean[2];

    outp_normal[0] = Rg2l(0, 2);
    outp_normal[1] = Rg2l(1, 2);
    outp_normal[2] = Rg2l(2, 2);

    outp_x[0] = Rg2l(0, 0);
    outp_x[1] = Rg2l(1, 0);
    outp_x[2] = Rg2l(2, 0);

    outp_y[0] = Rg2l(0, 1);
    outp_y[1] = Rg2l(1, 1);
    outp_y[2] = Rg2l(2, 1);

    //calculate point position in local coordinates; P_local = Rl2g * P_global = Rg2l.transpose() * P_global
    pointSet = Rg2l.transpose().eval() * pointSet;


    //Using the points in local coordinates xy plane to cal width and length, ignore the z-axis information
    //std::vector<double> x;
    //std::vector<double> y;

    for (int i = 0; i < pointSet.cols(); i++)
    {
        //x.push_back(pointSet.col(i)[0]);
        //y.push_back(pointSet.col(i)[1]);
        length = pointSet.row(0).maxCoeff() - pointSet.row(0).minCoeff();
        width = pointSet.row(1).maxCoeff() - pointSet.row(1).minCoeff();
    }

    return true;
}
