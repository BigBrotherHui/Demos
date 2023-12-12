#include <iostream>
#include <Eigen/Eigen>
class BasicAlgorithm
{
public:
	/**
	 * @brief 计算2个四元数的球面线性插值Spherical linear interpolation
	 * @param start_q起始位置
	 * @param end_q终止位置
	 * @param t某个时间的位置
	 * @return 插值矩阵
	 */
	static Eigen::Quaterniond Quaternion_Slerp(Eigen::Quaterniond& start_q, Eigen::Quaterniond& end_q, double t);
	static Eigen::Matrix4d RotateAroundAxisAtFixedPoint(const Eigen::Vector3d& point, const Eigen::Vector3d& dir, double angle);
	static Eigen::Matrix3d GenerateRotationMatrix(const Eigen::Vector3d& vectorBefore, const Eigen::Vector3d& vectorAfter);
	static double GetPointToLineDistance(double* point, double* linePoint_0, double* linePoint_1);
	static Eigen::Matrix4d ConvertCoordstoTransform(Eigen::Vector3d& o, Eigen::Vector3d& x, Eigen::Vector3d& y, Eigen::Vector3d& z);
	static void ProjectToPlane(const double x[3], const double origin[3], const double normal[3], double xproj[3]);
	static double AngleBetween2Vector(const double vec1[3], const double vec2[3], bool radian);
	/**
	 * @brief 计算两点距离
	 * @param double* 两点的点集
	 * @return 两点间距离
	 */
	static double DistanceBetweenTwoPoints(const double p1[3], const double p2[3]);
	/**
	 * @brief 计算点到直线的距离
	 * @param double* 点的坐标
	 * @param double* 直线上一个点的坐标
	 * @param double* 直线的方向向量(不一定为单位向量)
	 * @return 点到直线的距离
	 */
	static double DistanceFromPointToLineWithDirectionVector(const double PointOutsideLine[3], const double PointOnLine[3], const double direction[3]);
	/**
	 * @brief 计算点到平面的距离
	 * @param double* 点的坐标
	 * @param double* 平面上点的坐标
	 * @param double* 平面上的法线
	 * @return 点到直线的距离
	 */
	static double DistanceFromPointToPlane(const double x[3], const double p0[3], const double n[3]);
	/**
	 * @brief 计算两个矢量之间的角度
	 * @param Eigen::Vector3d 第一个向量
	 * @param Eigen::Vector3d 第二个向量
	 * @param bool 角度类型，默认值：false，返回度数；设置为true返回弧度类型
	 * @return double 点在平面上投影点的坐标
	 */
	static double AngleBetween2Vector(const Eigen::Vector3d vec1, const Eigen::Vector3d vec2, bool radian = false);
	/**
	 * @brief 计算两个矢量之间的角度
	 * @param double* 第一个向量
	 * @param double* 第二个向量
	 * @param double* 平面的正法线（平行于坐标轴）。
	 * @param bool 角度类型，默认值：false，返回度数；设置为true返回弧度类型
	 * @return double 返回的两个矢量之间的角度，由于有正法线，角度范围在-180度至180度之间，默认返回角度
	 */
	static double AngleBetween2Vector(const double vec1[3], const double vec2[3], const double normal[3], bool radian = false);
	/**
	 * @brief 计算两个矢量之间的角度
	 * @param Eigen::Vector3d 第一个向量
	 * @param Eigen::Vector3d 第二个向量
	 * @param double* 平面的正法线（平行于坐标轴）。
	 * @param bool 返回是否为弧度
	 * @return double 返回的两个矢量之间的角度，由于有正法线，默认返角度
	 */
	static double AngleBetween2Vector(const Eigen::Vector3d vec1, const Eigen::Vector3d vec2, const Eigen::Vector3d normal, bool radian = false);
	/**
	 * @brief 计算直线和直线之间的角度
	 * @param double* 直线一上的第一个点
	 * @param double* 直线一上的第二个点
	 * @param double* 直线二上的第一个点
	 * @param double* 直线二上的第二个点
	 * @param bool 返回是否为弧度
	 * @return double 返回的两个矢量之间的角度，由于有正法线，默认返角度
	 */
	static double AngleBetween2Line(const double* p11, const double* p12, const double* p21, const double* p22, bool radian = false);
	/**
	 * @brief 计算直线和平面之间的角度
	 * @param double* 直线一上的第一个点
	 * @param double* 直线一上的第二个点
	 * @param double* 平面的法线
	 * @param bool 返回是否为弧度
	 * @return double 返回直线和平面之间的角度
	 */
	static double AngleBetweenLineAndPlane(const double* p1, const double* p2, const double* normal, bool radian = false);
	/**
	 * @brief 计算直线和平面之间的角度
	 * @param Eigen::Vector3d 的方向向量
	 * @param Eigen::Vector3d 平面的法向量
	 * @param bool 返回是否为弧度
	 * @return double 返回直线和平面之间的角度
	 */
	static double AngleBetweenLineAndPlane(const Eigen::Vector3d vec, const Eigen::Vector3d normal, bool radian = false);
	/**
	 * @brief 在二维中用一组点拟合一个圆
	 * @param std::vector<double>& 点集x坐标集合
	 * @param std::vector<double>& 点集y坐标集合
	 * @param double& 圆心点x坐标
	 * @param double& 圆心点y坐标
	 */
	static void fit_circle_2d(std::vector<double>& inp_x, std::vector<double>& inp_y, double& outp_Cx, double& outp_Cy, double& outp_R);

	/**
	 * @brief 在三维中用一组点拟合一个圆
	 * @param std::vector<double>& 点集坐标
	 * @param std::vector<double>& 圆心
	 * @param double& 半径
	 * @param double& 法线
	 */
	static bool fit_circle_3d(std::vector<double>& inp_pointset, std::array<double, 3>& outp_center, double& outp_radius,
		std::array<double, 3>& outp_normal);
	/**
	 * @brief 在三维中用一组点拟合一个圆
	 * @param std::vector<double>& 点集坐标
	 * @param std::vector<double>& 点集圆心
	 * @param std::vector<double>& 点集半径
	 * @param double& 球心x坐标
	 * @param double& 球心y坐标
	 * @param double& 球心z坐标
	 * @param double& 球的半径
	 */
	static bool fit_sphere(std::vector<double>& inp_x, std::vector<double>& inp_y, std::vector<double>& inp_z, double& outp_cx, double& outp_cy, double& outp_cz, double& outp_R);
	/**
	 * @brief 在三维中用一组点拟合一个圆
	 * @param std::vector<double>&点集坐标，按x1,y1,z1,x2,y2,z2...顺序输入
	 * @param double& 球心坐标
	 * @param double& 球的半径
	 */
	static bool fit_sphere(std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, double& outp_r);
	/**
	 * @brief 在三维中用一组点拟合一个圆
	 * @param std::vector<std::array<double, 3>>&点坐标集合，
	 * @param double& 球心坐标
	 * @param double& 球的半径
	 */
	static bool fit_sphere(std::vector<std::array<double, 3>>& inp_pSet, std::array<double, 3>& outp_center, double& outp_r);

	/**
	 * @brief 将一组至少4个空间点拟合到一个固定半径的球体中。
	 * @param std::vector<double>& 点集x坐标
	 * @param std::vector<double>& 点集y坐标
	 * @param std::vector<double>& 点集z坐标
	 * @param std::vector<double>& 球体的半径
	 * @param double& 球心x坐标
	 * @param double& 球心y坐标
	 * @param double& 球心z坐标
	 */
	static bool fit_sphere_fixR(const std::vector<double>& inp_x, const std::vector<double>& inp_y, const std::vector<double>& inp_z, const double inp_r,
		double& outp_cx, double& outp_cy, double& outp_cz);

	/**
	 * @brief 点集拟合平面
	 * @param std::vector<double>& 点集坐标
	 * @param double& 平面中心
	 * @param double& 平面法线
	 */
	static bool fit_plane(const std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, std::array<double, 3>& outp_normal);

	/**
	 * @brief 点集拟合矩形
	 * @param std::vector<double>& 点集坐标
	 * @param double& 平面中心
	 * @param double& 平面法线
	 */
	static bool fit_rectangle(const std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, std::array<double, 3>& outp_normal, std::array<double, 3>& outp_x,
		std::array<double, 3>& outp_y, double& length, double& width);
};