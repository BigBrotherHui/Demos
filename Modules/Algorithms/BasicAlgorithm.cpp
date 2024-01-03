#include "BasicAlgorithm.h"

static const double M_PI=3.14159265358979323846;

bool BasicAlgorithm::fit_plane(const Eigen::MatrixX3d& plane_pts, Eigen::Vector4d &out)
{
    if(plane_pts.rows()<3)
    {
        std::cout << "at least input 3 points" << std::endl;
        return false;
    }
    Eigen::Vector3d center = plane_pts.colwise().mean().eval();
    Eigen::MatrixX3d tmp = plane_pts;
    tmp.rowwise() -= center.transpose();

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(tmp, Eigen::ComputeThinV);
    out.head(3) = svd.matrixV().block<3, 1>(0, 2);
    out[3] = -(out[0] * center[0] + out[1] * center[1] + out[2] * center[2]);
    return true;
}

void BasicAlgorithm::fit_sphere(const MatrixXd& A, Eigen::Vector3d& center, double& radius, double& rms)
{
    MatrixXd B = A;
    MatrixXd id(A.rows(), 1);
    for (int i = 0; i < A.rows(); i++) {
        id(i) = 1;
    }
    B.resize(A.rows(), 4);
    B.col(B.cols() - 1) = id;
    MatrixXd f;
    f = A.rowwise().squaredNorm();
    MatrixXd sol = (((B.transpose() * B).inverse()) * B.transpose()) * f;
    cout << "The least-squares solution using normal equations is:\n"
        << sol << endl;

    double x0 = sol(0) / 2;
    double y0 = sol(1) / 2;
    double z0 = sol(2) / 2;
    radius = sqrt(sol(3) + pow(x0, 2) + pow(y0, 2) + pow(z0, 2));

    RowVector3d rcenter;
    rcenter << x0, y0, z0;
    center << x0, y0, z0;

    cout << "Coordinates of the center is : (" << x0 << ", " << y0 << ", " << z0 << ")" << endl;
    cout << "Radius of the sphere is : " << radius << endl;

    //Cost of fitting
    MatrixXd d;
    MatrixXd del = A.rowwise() - rcenter;
    d = del.rowwise().norm();
    MatrixXd diff = d.array() - radius;
    MatrixXd result = diff.colwise().squaredNorm();
    double tmprms = result(0) / del.rows();
    rms = sqrt(tmprms);

    cout << "The cost of fitting is : " << rms << endl;
}

void BasicAlgorithm::polyfit(const std::vector<double>& v,const std::vector<double>& t, std::vector<double>& coeff, int order)
{
    // Create Matrix Placeholder of size n x k, n= number of datapoints, k = order of polynomial, for exame k = 3 for cubic polynomial
    Eigen::MatrixXd T(t.size(), order + 1);
    Eigen::VectorXd V = Eigen::VectorXd::Map(&v.front(), v.size());
    //std::cout<<"ceshi"<<std::endl;
    //std::cout<<V<<std::endl;
    Eigen::VectorXd result;

    // check to make sure inputs are correct
    assert(t.size() == v.size());
    assert(t.size() >= order + 1);
    // Populate the matrix
    for (size_t i = 0; i < t.size(); ++i)
    {
        for (size_t j = 0; j < order + 1; ++j)
        {
            T(i, j) = pow(t.at(i), j);
        }
    }
    std::cout << T << std::endl;

    // Solve for linear least square fit
    result = T.householderQr().solve(V);
    coeff.resize(order + 1);
    for (int k = 0; k < order + 1; k++)
    {
        coeff[k] = result[k];
    }
}

void BasicAlgorithm::polyfit_example()
{
    // time value
    std::vector<double> time = { 0, 0.0192341804504395, 0.0394501686096191,  0.059575080871582, 0.0790810585021973, 0.0792751312255859, 0.0987141132354736,  0.119336366653442,  0.138712167739868,  0.159000158309937,  0.178890228271484,   0.19960618019104,  0.219112157821655,   0.23919415473938,  0.259442090988159,  0.279186248779297,  0.299112319946289,  0.319219350814819,  0.339494228363037,  0.339675188064575,  0.359552145004272,   0.37941837310791,  0.399189233779907,  0.419828176498413,  0.439810276031494,  0.459331274032593,  0.479461193084717,  0.499663114547729,  0.519809246063232,  0.539092063903809,  0.559118270874023,  0.579315185546875,  0.598889112472534,  0.619685173034668,  0.638863086700439,  0.639052152633667,  0.658920288085938,  0.679149150848389,  0.699787139892578,   0.71905517578125,   0.73898720741272,  0.739143371582031,  0.758654117584229,  0.779210329055786,  0.799195289611816,  0.819046258926392,  0.839539289474487,   0.85923433303833,   0.87903618812561,  0.899263143539429,  0.919251203536987,  0.939138174057007,  0.959244251251221,  0.979074239730835,  0.998935222625732,   1.01904726028442,    1.0387852191925,   1.03895926475525,   1.05906510353088,   1.07873225212097,   1.09908628463745,   1.11907029151917,   1.13899827003479,   1.15879201889038 };
    // velocity value
    std::vector<double> velocity = { 1.8, 1.86, 2.03, 2.08, 2.14, 2.14, 2.25, 2.36, 2.42, 2.59,  2.7, 2.81, 2.87, 3.04, 3.15, 3.26, 3.32, 3.43, 3.54, 3.54,  3.6, 3.71, 3.83, 3.94, 4.11, 4.22, 4.33, 4.44, 4.56, 4.67, 4.78, 4.84, 4.84, 4.89, 4.89, 4.89, 4.95, 5.01, 5.06, 5.06, 5.06, 5.06, 5.01, 5.06, 5.12, 5.18, 5.18, 5.23, 5.23, 5.23, 5.29, 5.34, 5.29,  5.4,  5.4, 5.46, 5.51, 5.51, 5.51, 5.46,  5.4, 5.34, 5.34, 5.34 };

    // placeholder for storing polynomial coefficient
    std::vector<double> coeff;
    polyfit(time, velocity, coeff, 3);

    std::vector<double> fitted_velocity;
    std::cout << "Printing fitted values" << std::endl;
    for (int p = 0; p < time.size(); ++p)
    {
        double vfitted = coeff[0] + coeff[1] * time.at(p) + coeff[2] * (pow(time.at(p), 2)) + coeff[3] * (pow(time.at(p), 3));
        std::cout << vfitted << ", ";
        fitted_velocity.push_back(vfitted);
    }
    std::cout << std::endl;
}

Eigen::Matrix4d BasicAlgorithm::eulerAngle_translate_ToMatrix(double zAngle,double yAngle,double xAngle,double transZ,double transY,double transX)
{
    Eigen::Quaterniond quaternion3 = Eigen::AngleAxisd(zAngle, Eigen::Vector3d::UnitZ()) *
                  Eigen::AngleAxisd(yAngle, Eigen::Vector3d::UnitY()) *
                  Eigen::AngleAxisd(xAngle, Eigen::Vector3d::UnitX());
    Eigen::Vector3d translate(transX, transY, transZ);
    Eigen::Isometry3d T;
    T.setIdentity();
    T.rotate(quaternion3.toRotationMatrix());
    T.pretranslate(translate);
    return T.matrix();
}

Eigen::Quaterniond BasicAlgorithm::Quaternion_Slerp(Eigen::Quaterniond &start_q, Eigen::Quaterniond &end_q, double t)
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
    //第二种写法：旋转中心不变，旋转后用初始位置减去旋转后旋转中心位置
    /*Eigen::Vector3d x_world;
    x_world << 1, 0, 0;
    Eigen::Vector3d x_pelvis;
    x_pelvis = ASIS_L - ASIS_R;
    Eigen::Matrix3d rot = Eigen::Quaterniond().FromTwoVectors(x_pelvis, x_world).matrix();

    Eigen::Isometry3d T;
    T.setIdentity();
    T.rotate(rot);

    Eigen::Vector3d ASIS_R_roted = T * ASIS_R;
    Eigen::Vector3d t = ASIS_R - ASIS_R_roted;

    T.pretranslate(t);

    return T.matrix();*/
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
