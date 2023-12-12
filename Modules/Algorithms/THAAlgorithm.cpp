#include "THAAlgorithm.h"
using namespace Defintion;
void THAAlgorithm::CupAngle(double* direction, double& ResultAnteversion, double& ResultInclination,
	Defintion::ECupAngleType type)
{
	switch (type)
	{
	case ECupAngleType::RADIO_GRAPHIC:
		cupAngleRadiographic(direction, ResultAnteversion, ResultInclination);
		break;
	case ECupAngleType::OPERATIVE:
		cupAngleOperative(direction, ResultAnteversion, ResultInclination);
		break;
	case ECupAngleType::ANATOMICAL:
		cupAngleAnatomical(direction, ResultAnteversion, ResultInclination);
		break;
	default:
		cupAngleRadiographic(direction, ResultAnteversion, ResultInclination);
	}
}

void THAAlgorithm::cupAngleRadiographic(double* direction, double& ResultAnteversion, double& ResultInclination)
{
	double Coronalnormal[3] = { 0, -1, 0 };
	double D_value = 90.0;
	double CupLineTransform[3] = { direction[0], direction[1], direction[2] };
	double VersionRadians = BasicAlgorithm::AngleBetween2Vector(CupLineTransform, Coronalnormal, true);
	ResultAnteversion = (double)(VersionRadians * (180.0 / EIGEN_PI)) - D_value;
	if (ResultAnteversion < 0.0)
	{
		ResultAnteversion = -ResultAnteversion;
	}
	//std::cout << "ResultAnteversion:" << ResultAnteversion << std::endl;
	double CupLineProject[3];
	double origin[] = { 0, 0, 0 };
	//ͶӰ����״����
	BasicAlgorithm::ProjectToPlane(CupLineTransform, origin, Coronalnormal, CupLineProject);
	double z[3] = { 0, 0, 1 };
	double InclinationRadians = BasicAlgorithm::AngleBetween2Vector(CupLineProject, z, true);
	ResultInclination = (double)(InclinationRadians * (180.0 / EIGEN_PI));
	if (ResultInclination > 90.0)
	{
		ResultInclination = 180.0 - ResultInclination;
	}
	//std::cout << "ResultInclination:" << ResultInclination << std::endl;
}

void THAAlgorithm::cupAngleOperative(double* direction, double& ResultAnteversion, double& ResultInclination)
{
	double Sagittalnormal[3] = { -1, 0, 0 };
	double z[3] = { 0, 0, 1 };
	double D_value = 90.0;
	double CupLineTransform[3] = { direction[0], direction[1], direction[2] };
	double CupLineProject[3];
	double origin[] = { 0, 0, 0 };
	//ͶӰ��ʸ״����
	BasicAlgorithm::ProjectToPlane(CupLineTransform, origin, Sagittalnormal, CupLineProject);
	double VersionRadians = BasicAlgorithm::AngleBetween2Vector(z, CupLineProject, true);
	ResultAnteversion = (double)(VersionRadians * (180.0 / EIGEN_PI));
	if (ResultAnteversion > 90.0)
	{
		ResultAnteversion = 180.0 - ResultAnteversion;
	}
	//std::cout << "ResultAnteversion:" << ResultAnteversion << std::endl;

	double InclinationRadians = BasicAlgorithm::AngleBetween2Vector(CupLineTransform, Sagittalnormal, true);
	ResultInclination = (double)(InclinationRadians * (180.0 / EIGEN_PI)) - D_value;
	if (ResultInclination < 0.0)
	{
		ResultInclination = -ResultInclination;
	}
	//std::cout << "ResultInclination:" << ResultInclination << std::endl;
}

void THAAlgorithm::cupAngleAnatomical(double* direction, double& ResultAnteversion, double& ResultInclination)
{
	const double TransverseNormal[3] = { 0, 0, -1 };
	double z[3] = { 0, 0, 1 };
	double x[3] = { -1, 0, 0 };
	double CupLineTransform[3] = { direction[0], direction[1], direction[2] };
	double CupLineProject[3];
	double origin[] = { 0, 0, 0 };
	
	BasicAlgorithm::ProjectToPlane(CupLineTransform, origin, TransverseNormal, CupLineProject);
	double VersionRadians = BasicAlgorithm::AngleBetween2Vector(x, CupLineProject, true);
	ResultAnteversion = (double)(VersionRadians * (180.0 / EIGEN_PI));
	if (ResultAnteversion > 90.0)
	{
		ResultAnteversion = 180.0 - ResultAnteversion;
	}
	//std::cout << "ResultAnteversion:" << ResultAnteversion << std::endl;

	double InclinationRadians = BasicAlgorithm::AngleBetween2Vector(CupLineTransform, z, true);
	ResultInclination = (double)(InclinationRadians * (180.0 / EIGEN_PI));
	if (ResultInclination > 90.0)
	{
		ResultInclination = 180.0 - ResultInclination;
	}
	//std::cout << "ResultInclination:" << ResultInclination << std::endl;
}

double THAAlgorithm::FemoralVersionAngle(ESide side, double* FHC, double* FNC, double* ME, double* LE, double* DFCA, double* PFCA)
{
	// Compose neck axis: neck center + femurCOR
	Eigen::Vector3d neckAxis;
	neckAxis[0] = FHC[0] - FNC[0];
	neckAxis[1] = FHC[1] - FNC[1];
	neckAxis[2] = FHC[2] - FNC[2];
	// Compose epicondylar axis
	Eigen::Vector3d epicondylarAxis;
	epicondylarAxis[0] = ME[0] - LE[0];
	epicondylarAxis[1] = ME[1] - LE[1];
	epicondylarAxis[2] = ME[2] - LE[2];
	// Compose femur canal axis
	Eigen::Vector3d canalAxis;
	canalAxis[0] = PFCA[0] - DFCA[0];
	canalAxis[1] = PFCA[1] - DFCA[1];
	canalAxis[2] = PFCA[2] - DFCA[2];

	canalAxis.normalize();

	// Project the neckAxis onto the canal axis
	Eigen::Vector3d neckAxis_ontoCanalAxis = neckAxis.dot(canalAxis) * canalAxis;

	// neckAxis projection onto the perpendicular plane 
	Eigen::Vector3d neckAxis_ontoPlane = neckAxis - neckAxis_ontoCanalAxis;

	// Project the epicondylarAxis onto the canal axis
	Eigen::Vector3d epicondylarAxis_ontoCanalAxis = epicondylarAxis.dot(canalAxis) * canalAxis;

	// epicondylarAxis projection onto the perpendicular plane
	Eigen::Vector3d epicondylarAxis_ontoPlane = epicondylarAxis - epicondylarAxis_ontoCanalAxis;

	double femoralVersion = (180 / EIGEN_PI) * acos(epicondylarAxis_ontoPlane.dot(neckAxis_ontoPlane)
		/ (epicondylarAxis_ontoPlane.norm() * neckAxis_ontoPlane.norm()));

	// Determine anteversion or retroversion; if ante, assign femoralVersion as (+); if retro, assign as (-)
	double tmpValue = epicondylarAxis_ontoPlane.cross(neckAxis_ontoPlane).dot(canalAxis);

	if (side == ESide::left)
	{
		if (tmpValue < 0)
		{
			femoralVersion = -femoralVersion;
		}
	}
	else
	{
		if (tmpValue > 0)
		{
			femoralVersion = -femoralVersion;
		}
	}

	return femoralVersion;
}

double THAAlgorithm::PelvicTilt(double* pelvicYAxis)
{
	// if the patient direction meets the requirements, the supine pelvis tilt is angle between pelvicFrame's
	// y axis and worldFrame's y axis

	Eigen::Vector3d y_pelvicFrame;
	y_pelvicFrame << pelvicYAxis[0], pelvicYAxis[1], pelvicYAxis[2];

	Eigen::Vector3d y_worldFrame;
	y_worldFrame << 0, 1, 0;

	double tmpDotProduct = y_pelvicFrame.dot(y_worldFrame);
	Eigen::Vector3d tmpCrossProduct = y_pelvicFrame.cross(y_worldFrame);
	double angle = (180.0 / EIGEN_PI) * acos(tmpDotProduct);

	if (tmpCrossProduct[0] > 0)
	{
		angle = -angle;
	}
	return angle;
}

double THAAlgorithm::HipLength(double* LT, double* ASIS_L, double* ASIS_R)
{
	return ASIS_L[2] - LT[2];
}
