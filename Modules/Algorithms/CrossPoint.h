#pragma once

#ifndef CROSSPOINT_H
#define CROSSPOINT_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include"Eigen/Dense"
using namespace std;

void spacialLineFitting(Eigen::MatrixXd& pointsMatrix, Eigen::Matrix3d& V, double& centerX, double& centerY, double& centerZ);

void readPoint(string strPath, vector<vector<double>>& pointsVec);

int Intersection3DPoint(Eigen::Matrix3d dir_V, vector<double> center_V, Eigen::Matrix3d dir_H, vector<double> center_H, vector<double>& p0, vector<double>& p1);

#endif // !CROSSPOINT_H