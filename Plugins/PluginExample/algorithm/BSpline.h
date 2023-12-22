#pragma once

#include <vector>
#include <Eigen/Eigen>
using namespace Eigen;

struct BSpline {
	BSpline(int _p, const std::vector<double>& _us, int _type = 1);
	std::vector<double> computeCoefficients(double u) const;
	Eigen::Vector3d operator()(const std::vector<Eigen::Vector3d>& controlPoints, double u) const;

	int n = 0, m = 0, p = 0;
	std::vector<double> us;
	int bType; // 0 = open, 1 = clamped, 2 = closed;
};

struct BSplineSurface {
	BSplineSurface(int _p,
	               int _q,
	               const std::vector<double>& us,
	               const std::vector<double>& vs,
	               int utype = 1,
	               int vtype = 1);
	Eigen::Vector3d operator()(const std::vector<std::vector<Eigen::Vector3d>>& controlPoints, double u, double v) const;

	BSpline ubsp, vbsp;
};
