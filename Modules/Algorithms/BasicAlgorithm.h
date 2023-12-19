#include <iostream>
#include <Eigen/Eigen>
using namespace Eigen;
using namespace std;
class BasicAlgorithm
{
public:
	//A为N行3列的点
	static void fit_sphere(const MatrixXd& A, Eigen::Vector3d& center, double& radius, double& rms);
	//n次多项式拟合，order为多项式次数，coeff为多项式系数
	static void polyfit(const std::vector<double>& v,const std::vector<double>& t, std::vector<double>& coeff, int order);
	void polyfit_example();
	//default by zyx
	Eigen::Matrix4d eulerAngle_translate_ToMatrix(double zAngle,double yAngle,double xAngle,double transZ,double transY,double transX);

	static Eigen::Quaterniond Quaternion_Slerp(Eigen::Quaterniond& start_q, Eigen::Quaterniond& end_q, double t);
	static Eigen::Matrix4d RotateAroundAxisAtFixedPoint(const Eigen::Vector3d& point, const Eigen::Vector3d& dir, double angle);
	static Eigen::Matrix3d GenerateRotationMatrix(const Eigen::Vector3d& vectorBefore, const Eigen::Vector3d& vectorAfter);
	static double GetPointToLineDistance(double* point, double* linePoint_0, double* linePoint_1);
	static Eigen::Matrix4d ConvertCoordstoTransform(Eigen::Vector3d& o, Eigen::Vector3d& x, Eigen::Vector3d& y, Eigen::Vector3d& z);
	static void ProjectToPlane(const double x[3], const double origin[3], const double normal[3], double xproj[3]);
	static double AngleBetween2Vector(const double vec1[3], const double vec2[3], bool radian);

	static double DistanceBetweenTwoPoints(const double p1[3], const double p2[3]);

	static double DistanceFromPointToLineWithDirectionVector(const double PointOutsideLine[3], const double PointOnLine[3], const double direction[3]);

	static double DistanceFromPointToPlane(const double x[3], const double p0[3], const double n[3]);

	static double AngleBetween2Vector(const Eigen::Vector3d vec1, const Eigen::Vector3d vec2, bool radian = false);

	static double AngleBetween2Vector(const double vec1[3], const double vec2[3], const double normal[3], bool radian = false);

	static double AngleBetween2Vector(const Eigen::Vector3d vec1, const Eigen::Vector3d vec2, const Eigen::Vector3d normal, bool radian = false);

	static double AngleBetween2Line(const double* p11, const double* p12, const double* p21, const double* p22, bool radian = false);

	static double AngleBetweenLineAndPlane(const double* p1, const double* p2, const double* normal, bool radian = false);

	static double AngleBetweenLineAndPlane(const Eigen::Vector3d vec, const Eigen::Vector3d normal, bool radian = false);

	static void fit_circle_2d(std::vector<double>& inp_x, std::vector<double>& inp_y, double& outp_Cx, double& outp_Cy, double& outp_R);

	static bool fit_circle_3d(std::vector<double>& inp_pointset, std::array<double, 3>& outp_center, double& outp_radius,
		std::array<double, 3>& outp_normal);

	static bool fit_sphere(std::vector<double>& inp_x, std::vector<double>& inp_y, std::vector<double>& inp_z, double& outp_cx, double& outp_cy, double& outp_cz, double& outp_R);

	static bool fit_sphere(std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, double& outp_r);

	static bool fit_sphere(std::vector<std::array<double, 3>>& inp_pSet, std::array<double, 3>& outp_center, double& outp_r);

	static bool fit_sphere_fixR(const std::vector<double>& inp_x, const std::vector<double>& inp_y, const std::vector<double>& inp_z, const double inp_r,
		double& outp_cx, double& outp_cy, double& outp_cz);

	static bool fit_plane(const std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, std::array<double, 3>& outp_normal);

	static bool fit_rectangle(const std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, std::array<double, 3>& outp_normal, std::array<double, 3>& outp_x,
		std::array<double, 3>& outp_y, double& length, double& width);
};