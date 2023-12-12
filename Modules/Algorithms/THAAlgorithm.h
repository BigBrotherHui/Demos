#include <iostream>
#include <Eigen/Eigen>
#include "Definition.h"
#include "BasicAlgorithm.h"
class THAAlgorithm
{
	/**
	 * compute the Anteversion angle and the Inclination Angle
	 *
	 * calculation is based on a radio_graphic coordinate system: x left , y backward , z head.
	 *
	 * @param direction [Input]The direction of the cup normal or Acetabular grinding rod.
	 * @param ResultAnteversion [Output] The compute result of Anteversion.
	 * @param ResultInclination [Output] The compute result of Inclination.
	 * @param type [Input] different type of Anteversion/Inclination angle, default:Radiographic
	 */
	void CupAngle(double* direction, double& ResultAnteversion, double& ResultInclination,
		Defintion::ECupAngleType type = Defintion::ECupAngleType::RADIO_GRAPHIC);
	void cupAngleRadiographic(double* direction, double& ResultAnteversion, double& ResultInclination);//影像学
	void cupAngleOperative(double* direction, double& ResultAnteversion, double& ResultInclination);//手术学
	void cupAngleAnatomical(double* direction, double& ResultAnteversion, double& ResultInclination);//解剖学
	/**
	* \brief Compute the Femoral Version Angle.
	*
	* can use either local or global points, all points are on Femur.
	*	The native femoral version is the angle between the neck axis and epicondylar axis when
	* these 2 axes are projected on a plane perpendicular to the femur canal
	* a positive value (anteversion) is returned when the neck axis angled anteriorly relative to the
	* epicondylar axis.
	*
	* @param side [Input] define right or left femur.
	* @param FHC [Input] Femur Head Center.
	* @param FNC [Input] Femur Neck Center.
	* @param ME [Input] Medial Femoral Epicondyle.
	* @param LE [Input] Lateral Femoral Epicondyle.
	* @param DFCA [Input] Distal point of Femoral canal axis.
	* @param PFCA [Input] Proximal point of Femoral canal axis.
	* @return  Femoral Version Angle
	*/
	double FemoralVersionAngle(Defintion::ESide side, double* FHC, double* FNC, double* ME, double* LE, double* DFCA, double* PFCA);
	/**
	 * \brief Pelvic Tilt is the angle between APP and coronal plane. 
	 * when the line of the anterior superior iliac spine (X-axis) is virtually corrected
	 * to be parallel to the inner lateral (horizontal) direction.
	 * The Angle between the Y axis of the pelvis and the Y axis of the world coordinate system .
	 * “骨盆倾斜”。
	 * 在医学和运动学的背景下，指的是骨盆相对于垂直轴线的倾斜或旋转。
	 * 这个术语通常用于描述骨盆的前倾（anterior pelvic tilt）或后倾（posterior pelvic tilt），
	 * 这在解释姿势问题、运动损伤或治疗方案时非常有用
	 * \param pelvicYAxis Perpendicular to the plane defined by the ASIS and midLine point and pointing backward.
	 * \return Pelvic Tilt angle
	 */
	double PelvicTilt(double* pelvicYAxis);
	/**
	 * \brief hip length is the distance from the lesser trochanter to the ASIS line.
	 * when the femur is mechanically aligned and the maximum femural offset is parallel to the coronal plane.
	 *
	 * MUST use global points in aligned pose.
	 * \param LT Lesser Trochanter,Left or Right
	 * \param ASIS_L Anterior Superior lliac Spine��Left��
	 * \param ASIS_R Anterior Superior lliac Spine��Right��
	 * \return hip length(mm)
	 */
	double HipLength(double* LT, double* ASIS_L, double* ASIS_R);
};