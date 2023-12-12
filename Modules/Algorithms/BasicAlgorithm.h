#include <iostream>
#include <Eigen/Eigen>
class BasicAlgorithm
{
public:
	/**
	 * @brief ����2����Ԫ�����������Բ�ֵSpherical linear interpolation
	 * @param start_q��ʼλ��
	 * @param end_q��ֹλ��
	 * @param tĳ��ʱ���λ��
	 * @return ��ֵ����
	 */
	static Eigen::Quaterniond Quaternion_Slerp(Eigen::Quaterniond& start_q, Eigen::Quaterniond& end_q, double t);
	static Eigen::Matrix4d RotateAroundAxisAtFixedPoint(const Eigen::Vector3d& point, const Eigen::Vector3d& dir, double angle);
	static Eigen::Matrix3d GenerateRotationMatrix(const Eigen::Vector3d& vectorBefore, const Eigen::Vector3d& vectorAfter);
	static double GetPointToLineDistance(double* point, double* linePoint_0, double* linePoint_1);
	static Eigen::Matrix4d ConvertCoordstoTransform(Eigen::Vector3d& o, Eigen::Vector3d& x, Eigen::Vector3d& y, Eigen::Vector3d& z);
	static void ProjectToPlane(const double x[3], const double origin[3], const double normal[3], double xproj[3]);
	static double AngleBetween2Vector(const double vec1[3], const double vec2[3], bool radian);
	/**
	 * @brief �����������
	 * @param double* ����ĵ㼯
	 * @return ��������
	 */
	static double DistanceBetweenTwoPoints(const double p1[3], const double p2[3]);
	/**
	 * @brief ����㵽ֱ�ߵľ���
	 * @param double* �������
	 * @param double* ֱ����һ���������
	 * @param double* ֱ�ߵķ�������(��һ��Ϊ��λ����)
	 * @return �㵽ֱ�ߵľ���
	 */
	static double DistanceFromPointToLineWithDirectionVector(const double PointOutsideLine[3], const double PointOnLine[3], const double direction[3]);
	/**
	 * @brief ����㵽ƽ��ľ���
	 * @param double* �������
	 * @param double* ƽ���ϵ������
	 * @param double* ƽ���ϵķ���
	 * @return �㵽ֱ�ߵľ���
	 */
	static double DistanceFromPointToPlane(const double x[3], const double p0[3], const double n[3]);
	/**
	 * @brief ��������ʸ��֮��ĽǶ�
	 * @param Eigen::Vector3d ��һ������
	 * @param Eigen::Vector3d �ڶ�������
	 * @param bool �Ƕ����ͣ�Ĭ��ֵ��false�����ض���������Ϊtrue���ػ�������
	 * @return double ����ƽ����ͶӰ�������
	 */
	static double AngleBetween2Vector(const Eigen::Vector3d vec1, const Eigen::Vector3d vec2, bool radian = false);
	/**
	 * @brief ��������ʸ��֮��ĽǶ�
	 * @param double* ��һ������
	 * @param double* �ڶ�������
	 * @param double* ƽ��������ߣ�ƽ���������ᣩ��
	 * @param bool �Ƕ����ͣ�Ĭ��ֵ��false�����ض���������Ϊtrue���ػ�������
	 * @return double ���ص�����ʸ��֮��ĽǶȣ������������ߣ��Ƕȷ�Χ��-180����180��֮�䣬Ĭ�Ϸ��ؽǶ�
	 */
	static double AngleBetween2Vector(const double vec1[3], const double vec2[3], const double normal[3], bool radian = false);
	/**
	 * @brief ��������ʸ��֮��ĽǶ�
	 * @param Eigen::Vector3d ��һ������
	 * @param Eigen::Vector3d �ڶ�������
	 * @param double* ƽ��������ߣ�ƽ���������ᣩ��
	 * @param bool �����Ƿ�Ϊ����
	 * @return double ���ص�����ʸ��֮��ĽǶȣ������������ߣ�Ĭ�Ϸ��Ƕ�
	 */
	static double AngleBetween2Vector(const Eigen::Vector3d vec1, const Eigen::Vector3d vec2, const Eigen::Vector3d normal, bool radian = false);
	/**
	 * @brief ����ֱ�ߺ�ֱ��֮��ĽǶ�
	 * @param double* ֱ��һ�ϵĵ�һ����
	 * @param double* ֱ��һ�ϵĵڶ�����
	 * @param double* ֱ�߶��ϵĵ�һ����
	 * @param double* ֱ�߶��ϵĵڶ�����
	 * @param bool �����Ƿ�Ϊ����
	 * @return double ���ص�����ʸ��֮��ĽǶȣ������������ߣ�Ĭ�Ϸ��Ƕ�
	 */
	static double AngleBetween2Line(const double* p11, const double* p12, const double* p21, const double* p22, bool radian = false);
	/**
	 * @brief ����ֱ�ߺ�ƽ��֮��ĽǶ�
	 * @param double* ֱ��һ�ϵĵ�һ����
	 * @param double* ֱ��һ�ϵĵڶ�����
	 * @param double* ƽ��ķ���
	 * @param bool �����Ƿ�Ϊ����
	 * @return double ����ֱ�ߺ�ƽ��֮��ĽǶ�
	 */
	static double AngleBetweenLineAndPlane(const double* p1, const double* p2, const double* normal, bool radian = false);
	/**
	 * @brief ����ֱ�ߺ�ƽ��֮��ĽǶ�
	 * @param Eigen::Vector3d �ķ�������
	 * @param Eigen::Vector3d ƽ��ķ�����
	 * @param bool �����Ƿ�Ϊ����
	 * @return double ����ֱ�ߺ�ƽ��֮��ĽǶ�
	 */
	static double AngleBetweenLineAndPlane(const Eigen::Vector3d vec, const Eigen::Vector3d normal, bool radian = false);
	/**
	 * @brief �ڶ�ά����һ������һ��Բ
	 * @param std::vector<double>& �㼯x���꼯��
	 * @param std::vector<double>& �㼯y���꼯��
	 * @param double& Բ�ĵ�x����
	 * @param double& Բ�ĵ�y����
	 */
	static void fit_circle_2d(std::vector<double>& inp_x, std::vector<double>& inp_y, double& outp_Cx, double& outp_Cy, double& outp_R);

	/**
	 * @brief ����ά����һ������һ��Բ
	 * @param std::vector<double>& �㼯����
	 * @param std::vector<double>& Բ��
	 * @param double& �뾶
	 * @param double& ����
	 */
	static bool fit_circle_3d(std::vector<double>& inp_pointset, std::array<double, 3>& outp_center, double& outp_radius,
		std::array<double, 3>& outp_normal);
	/**
	 * @brief ����ά����һ������һ��Բ
	 * @param std::vector<double>& �㼯����
	 * @param std::vector<double>& �㼯Բ��
	 * @param std::vector<double>& �㼯�뾶
	 * @param double& ����x����
	 * @param double& ����y����
	 * @param double& ����z����
	 * @param double& ��İ뾶
	 */
	static bool fit_sphere(std::vector<double>& inp_x, std::vector<double>& inp_y, std::vector<double>& inp_z, double& outp_cx, double& outp_cy, double& outp_cz, double& outp_R);
	/**
	 * @brief ����ά����һ������һ��Բ
	 * @param std::vector<double>&�㼯���꣬��x1,y1,z1,x2,y2,z2...˳������
	 * @param double& ��������
	 * @param double& ��İ뾶
	 */
	static bool fit_sphere(std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, double& outp_r);
	/**
	 * @brief ����ά����һ������һ��Բ
	 * @param std::vector<std::array<double, 3>>&�����꼯�ϣ�
	 * @param double& ��������
	 * @param double& ��İ뾶
	 */
	static bool fit_sphere(std::vector<std::array<double, 3>>& inp_pSet, std::array<double, 3>& outp_center, double& outp_r);

	/**
	 * @brief ��һ������4���ռ����ϵ�һ���̶��뾶�������С�
	 * @param std::vector<double>& �㼯x����
	 * @param std::vector<double>& �㼯y����
	 * @param std::vector<double>& �㼯z����
	 * @param std::vector<double>& ����İ뾶
	 * @param double& ����x����
	 * @param double& ����y����
	 * @param double& ����z����
	 */
	static bool fit_sphere_fixR(const std::vector<double>& inp_x, const std::vector<double>& inp_y, const std::vector<double>& inp_z, const double inp_r,
		double& outp_cx, double& outp_cy, double& outp_cz);

	/**
	 * @brief �㼯���ƽ��
	 * @param std::vector<double>& �㼯����
	 * @param double& ƽ������
	 * @param double& ƽ�淨��
	 */
	static bool fit_plane(const std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, std::array<double, 3>& outp_normal);

	/**
	 * @brief �㼯��Ͼ���
	 * @param std::vector<double>& �㼯����
	 * @param double& ƽ������
	 * @param double& ƽ�淨��
	 */
	static bool fit_rectangle(const std::vector<double>& inp_pSet, std::array<double, 3>& outp_center, std::array<double, 3>& outp_normal, std::array<double, 3>& outp_x,
		std::array<double, 3>& outp_y, double& length, double& width);
};