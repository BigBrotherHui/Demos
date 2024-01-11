#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <windows.h>
#include<algorithm>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<vector>

#include<time.h>
#include "kernel.cuh"
using namespace std;
//const int BlockSize = 512;
__device__ float dx = 0;
__device__ float dy = 0;
__device__ float dz = 0;//20


void savedata1(const short* rst1, int len, std::string  st)
{
	FILE *fpwrt = NULL;
	const char* file_c = st.c_str();
	fopen_s(&fpwrt, file_c, "wb+");
	if (fpwrt == NULL)
	{
		std::cout << "error write file" << std::endl;
	}
	fwrite(rst1, sizeof(short), len, fpwrt);
	fclose(fpwrt);
}

void savedata1(const float* rst1, int len, std::string  st)
{
	FILE* fpwrt = NULL;
	const char* file_c = st.c_str();
	fopen_s(&fpwrt, file_c, "wb+");
	if (fpwrt == NULL)
	{
		std::cout << "error write file" << std::endl;
	}
	fwrite(rst1, sizeof(float), len, fpwrt);
	fclose(fpwrt);
}

__global__ void CUDAprojection(float *image_3d, float *proj_m, float* _mask, float* GmtrcTrnsfrmtnMtrx, int imageM, int imageN, int imageH, float Pixeld, float Pixeldz, int projM, int projN, float Threshold) {
	//GmtrcTrnsfrmtnMtrx大小：float[12]							
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	int u = j / projN;
	int v = j % projN;
	proj_m[j] = 0;
	if (u > projM || v > projN)
		return ;


	//if (_mask[j] < 0.0001)
	//{
	//	proj_m[j] = 0;
	//	return;
	//}
	//int u = j % projN;
	//int v = j / projN;   //第一维度为u

	//if (v < 30 || u < 3) return;
	int flag = 0;
	//if(u<0||u>511||v<0||v>511)return;
	float point[3][4];//save tmp 3dim bound-point
	for (int i = 0; i < 3; i++) {
		for (int k = 0; k < 4; k++) {
			point[i][k] = 0.0;
		}
	}
	float p_xyz[6];//save final two point
	for (int i = 0; i < 6; i++) {
		p_xyz[i] = 0.0;
	}
	
	float py = 0;
	//image origin
	float xb = (0 - imageM / 2) * Pixeld + dx, yb = (0 - imageN / 2) * Pixeld + dy, zb = (0 - imageH / 2) * Pixeldz + dz - py;
	//image isocenter
	float xe = (imageM - 1 - imageM / 2) * Pixeld + dx, ye = (imageN - 1 - imageN / 2) * Pixeld + dy, ze = (imageH - 1 - imageH / 2) * Pixeldz + dz - py;

	//GmtrcTrnsfrmtnMtrx is transform matrix
	float t1 = GmtrcTrnsfrmtnMtrx[0] - u * GmtrcTrnsfrmtnMtrx[8];
	float t2 = GmtrcTrnsfrmtnMtrx[1] - u * GmtrcTrnsfrmtnMtrx[9];
	float t3 = GmtrcTrnsfrmtnMtrx[4] - v * GmtrcTrnsfrmtnMtrx[8];
	float t4 = GmtrcTrnsfrmtnMtrx[5] - v * GmtrcTrnsfrmtnMtrx[9];
	float a = (GmtrcTrnsfrmtnMtrx[10] * zb + GmtrcTrnsfrmtnMtrx[11]) * u - (GmtrcTrnsfrmtnMtrx[2] * zb + GmtrcTrnsfrmtnMtrx[3]);
	float b = (GmtrcTrnsfrmtnMtrx[10] * zb + GmtrcTrnsfrmtnMtrx[11]) * v - (GmtrcTrnsfrmtnMtrx[6] * zb + GmtrcTrnsfrmtnMtrx[7]);

	point[0][0] = (t4 * a - t2 * b) / (t1 * t4 - t2 * t3);
	point[1][0] = (t3 * a - t1 * b) / (t2 * t3 - t1 * t4);
	if (point[0][0] >= xb && point[0][0] <= xe && point[1][0] >= yb && point[1][0] <= ye) {
		flag++;
		//cout << "X:" << point[0][0] <<" "<< "Y:" << point[1][0] <<" "<< "Z:" << zb << endl;
		if (flag == 1) {
			p_xyz[0] = (point[0][0] - dx) / Pixeld + imageM / 2;
			p_xyz[1] = (point[1][0] - dy) / Pixeld + imageN / 2;
			p_xyz[2] = (zb - dz + py) / Pixeldz + imageH / 2;
			//proj_m[j] +=point[0][0];
		}
		if (flag == 2) {
			p_xyz[3] = (point[0][0] - dx) / Pixeld + imageM / 2;
			p_xyz[4] = (point[1][0] - dy) / Pixeld + imageN / 2;
			p_xyz[5] = (zb - dz + py) / Pixeldz + imageH / 2;
			//proj_m[j] +=point[0][0];
		}

	}


	t1 = GmtrcTrnsfrmtnMtrx[0] - u * GmtrcTrnsfrmtnMtrx[8];
	t2 = GmtrcTrnsfrmtnMtrx[1] - u * GmtrcTrnsfrmtnMtrx[9];
	t3 = GmtrcTrnsfrmtnMtrx[4] - v * GmtrcTrnsfrmtnMtrx[8];
	t4 = GmtrcTrnsfrmtnMtrx[5] - v * GmtrcTrnsfrmtnMtrx[9];
	a = (GmtrcTrnsfrmtnMtrx[10] * ze + GmtrcTrnsfrmtnMtrx[11]) * u - (GmtrcTrnsfrmtnMtrx[2] * ze + GmtrcTrnsfrmtnMtrx[3]);
	b = (GmtrcTrnsfrmtnMtrx[10] * ze + GmtrcTrnsfrmtnMtrx[11]) * v - (GmtrcTrnsfrmtnMtrx[6] * ze + GmtrcTrnsfrmtnMtrx[7]);
	point[0][1] = (t4 * a - t2 * b) / (t1 * t4 - t2 * t3);
	point[1][1] = (t3 * a - t1 * b) / (t2 * t3 - t1 * t4);
	if (point[0][1] >= xb && point[0][1] <= xe && point[1][1] >= yb && point[1][1] <= ye) {
		flag++;
		//cout << "X:" << point[0][1] <<" "<< "Y:" << point[1][1] <<" "<< "Z:" << ze << endl;
		if (flag == 1) {
			p_xyz[0] = (point[0][1] - dx) / Pixeld + imageM / 2;
			p_xyz[1] = (point[1][1] - dy) / Pixeld + imageN / 2;
			p_xyz[2] = (ze - dz + py) / Pixeldz + imageH / 2;
			//proj_m[j] +=point[0][1];
		}
		if (flag == 2) {
			p_xyz[3] = (point[0][1] - dx) / Pixeld + imageM / 2;
			p_xyz[4] = (point[1][1] - dy) / Pixeld + imageN / 2;
			p_xyz[5] = (ze - dz + py) / Pixeldz + imageH / 2;
			//proj_m[j] +=point[0][1];
		}
	}


	t1 = GmtrcTrnsfrmtnMtrx[0] - u * GmtrcTrnsfrmtnMtrx[8];
	t2 = GmtrcTrnsfrmtnMtrx[2] - u * GmtrcTrnsfrmtnMtrx[10];
	t3 = GmtrcTrnsfrmtnMtrx[4] - v * GmtrcTrnsfrmtnMtrx[8];
	t4 = GmtrcTrnsfrmtnMtrx[6] - v * GmtrcTrnsfrmtnMtrx[10];
	a = (GmtrcTrnsfrmtnMtrx[9] * yb + GmtrcTrnsfrmtnMtrx[11]) * u - (GmtrcTrnsfrmtnMtrx[1] * yb + GmtrcTrnsfrmtnMtrx[3]);
	b = (GmtrcTrnsfrmtnMtrx[9] * yb + GmtrcTrnsfrmtnMtrx[11]) * v - (GmtrcTrnsfrmtnMtrx[5] * yb + GmtrcTrnsfrmtnMtrx[7]);
	point[0][2] = (t4 * a - t2 * b) / (t1 * t4 - t2 * t3);
	point[2][0] = (t3 * a - t1 * b) / (t2 * t3 - t1 * t4);
	if (point[0][2] >= xb && point[0][2] <= xe && point[2][0] >= zb && point[2][0] <= ze) {
		flag++;
		//cout << "X:" << point[0][2] <<" "<< "Y:" << yb <<" "<< "Z:" << point[2][0] << endl;
		if (flag == 1) {
			p_xyz[0] = (point[0][2] - dx) / Pixeld + imageM / 2;
			p_xyz[1] = (yb - dy) / Pixeld + imageN / 2;
			p_xyz[2] = (point[2][0] - dz + py) / Pixeldz + imageH / 2;
			//proj_m[j] +=point[0][2];
		}
		if (flag == 2) {
			p_xyz[3] = (point[0][2] - dx) / Pixeld + imageM / 2;
			p_xyz[4] = (yb - dy) / Pixeld + imageN / 2;
			p_xyz[5] = (point[2][0] - dz + py) / Pixeldz + imageH / 2;
			//proj_m[j] +=point[0][2];
		}
	}


	t1 = GmtrcTrnsfrmtnMtrx[0] - u * GmtrcTrnsfrmtnMtrx[8];
	t2 = GmtrcTrnsfrmtnMtrx[2] - u * GmtrcTrnsfrmtnMtrx[10];
	t3 = GmtrcTrnsfrmtnMtrx[4] - v * GmtrcTrnsfrmtnMtrx[8];
	t4 = GmtrcTrnsfrmtnMtrx[6] - v * GmtrcTrnsfrmtnMtrx[10];
	a = (GmtrcTrnsfrmtnMtrx[9] * ye + GmtrcTrnsfrmtnMtrx[11]) * u - (GmtrcTrnsfrmtnMtrx[1] * ye + GmtrcTrnsfrmtnMtrx[3]);
	b = (GmtrcTrnsfrmtnMtrx[9] * ye + GmtrcTrnsfrmtnMtrx[11]) * v - (GmtrcTrnsfrmtnMtrx[5] * ye + GmtrcTrnsfrmtnMtrx[7]);
	point[0][3] = (t4 * a - t2 * b) / (t1 * t4 - t2 * t3);
	point[2][1] = (t3 * a - t1 * b) / (t2 * t3 - t1 * t4);
	if (point[0][3] >= xb && point[0][3] <= xe && point[2][1] >= zb && point[2][1] <= ze) {
		flag++;
		//cout << "X:" << point[0][3] <<" "<< "Y:" << ye <<" "<< "Z:" << point[2][1] << endl;
		if (flag == 1) {
			p_xyz[0] = (point[0][3] - dx) / Pixeld + imageM / 2;
			p_xyz[1] = (ye - dy) / Pixeld + imageN / 2;
			p_xyz[2] = (point[2][1] - dz + py) / Pixeldz + imageH / 2;
		}
		if (flag == 2) {
			p_xyz[3] = (point[0][3] - dx) / Pixeld + imageM / 2;
			p_xyz[4] = (ye - dy) / Pixeld + imageN / 2;
			p_xyz[5] = (point[2][1] - dz + py) / Pixeldz + imageH / 2;
		}
	}


	t1 = GmtrcTrnsfrmtnMtrx[1] - u * GmtrcTrnsfrmtnMtrx[9];
	t2 = GmtrcTrnsfrmtnMtrx[2] - u * GmtrcTrnsfrmtnMtrx[10];
	t3 = GmtrcTrnsfrmtnMtrx[5] - v * GmtrcTrnsfrmtnMtrx[9];
	t4 = GmtrcTrnsfrmtnMtrx[6] - v * GmtrcTrnsfrmtnMtrx[10];
	a = (GmtrcTrnsfrmtnMtrx[8] * xb + GmtrcTrnsfrmtnMtrx[11]) * u - (GmtrcTrnsfrmtnMtrx[0] * xb + GmtrcTrnsfrmtnMtrx[3]);
	b = (GmtrcTrnsfrmtnMtrx[8] * xb + GmtrcTrnsfrmtnMtrx[11]) * v - (GmtrcTrnsfrmtnMtrx[4] * xb + GmtrcTrnsfrmtnMtrx[7]);
	point[1][2] = (t4 * a - t2 * b) / (t1 * t4 - t2 * t3);
	point[2][2] = (t3 * a - t1 * b) / (t2 * t3 - t1 * t4);
	if (point[1][2] >= yb && point[1][2] <= ye && point[2][2] >= zb && point[2][2] <= ze) {
		flag++;
		//cout << "X:" << xb <<" "<< "Y:" << point[1][2] <<" "<< "Z:" << point[2][2] << endl;
		if (flag == 1) {
			p_xyz[0] = (xb - dx) / Pixeld + imageM / 2;
			p_xyz[1] = (point[1][2] - dy) / Pixeld + imageN / 2;
			p_xyz[2] = (point[2][2] - dz + py) / Pixeldz + imageH / 2;
		}
		if (flag == 2) {
			p_xyz[3] = (xb - dx) / Pixeld + imageM / 2;
			p_xyz[4] = (point[1][2] - dy) / Pixeld + imageN / 2;
			p_xyz[5] = (point[2][2] - dz + py) / Pixeldz + imageH / 2;
		}
	}


	t1 = GmtrcTrnsfrmtnMtrx[1] - u * GmtrcTrnsfrmtnMtrx[9];
	t2 = GmtrcTrnsfrmtnMtrx[2] - u * GmtrcTrnsfrmtnMtrx[10];
	t3 = GmtrcTrnsfrmtnMtrx[5] - v * GmtrcTrnsfrmtnMtrx[9];
	t4 = GmtrcTrnsfrmtnMtrx[6] - v * GmtrcTrnsfrmtnMtrx[10];
	a = (GmtrcTrnsfrmtnMtrx[8] * xe + GmtrcTrnsfrmtnMtrx[11]) * u - (GmtrcTrnsfrmtnMtrx[0] * xe + GmtrcTrnsfrmtnMtrx[3]);
	b = (GmtrcTrnsfrmtnMtrx[8] * xe + GmtrcTrnsfrmtnMtrx[11]) * v - (GmtrcTrnsfrmtnMtrx[4] * xe + GmtrcTrnsfrmtnMtrx[7]);
	point[1][3] = (t4 * a - t2 * b) / (t1 * t4 - t2 * t3);
	point[2][3] = (t3 * a - t1 * b) / (t2 * t3 - t1 * t4);
	if (point[1][3] >= yb && point[1][3] <= ye && point[2][3] >= zb && point[2][3] <= ze) {
		flag++;
		//cout << "X:" << xe <<" "<< "Y:" << point[1][3] <<" "<< "Z:" << point[2][3] << endl;
		if (flag == 1) {
			p_xyz[0] = (xe - dx) / Pixeld + imageM / 2;
			p_xyz[1] = (point[1][3] - dy) / Pixeld + imageN / 2;
			p_xyz[2] = (point[2][3] - dz + py) / Pixeldz + imageH / 2;
		}
		if (flag == 2) {
			p_xyz[3] = (xe - dx) / Pixeld + imageM / 2;
			p_xyz[4] = (point[1][3] - dy) / Pixeld + imageN / 2;
			p_xyz[5] = (point[2][3] - dz + py) / Pixeldz + imageH / 2;
		}
	}

	int i_min = floor(min(p_xyz[0], p_xyz[3]));
	int i_max = ceil(max(p_xyz[0], p_xyz[3]));
	int j_min = floor(min(p_xyz[1], p_xyz[4]));
	int j_max = ceil(max(p_xyz[1], p_xyz[4]));
	int k_min = floor(min(p_xyz[2], p_xyz[5]));
	int k_max = ceil(max(p_xyz[2], p_xyz[5]));
	int N = (i_max - i_min + 1) + (j_max - j_min + 1) + (k_max - k_min + 1);
	//matrix_jkl[0][max_N - 1] = N;
	int NX = i_max - i_min + 1;
	int NY = j_max - j_min + 1;
	int NZ = k_max - k_min + 1;
	const int dimx = 512;
	const int dimy = 512;
	const int dimz = 906;
	float alphax[dimx+1];
	float alphay[dimy+1];
	float alphaz[dimz+1];

	for (int i = i_min; i <= i_max; i++) {
		if (p_xyz[0] == p_xyz[3]) {
			alphax[i - i_min] = 1;
			break;
		}
		if (p_xyz[3] > p_xyz[0]) {
			alphax[i - i_min] = (i - p_xyz[0]) / (p_xyz[3] - p_xyz[0]);
		}
		else {
			alphax[i - i_min] = ((i_max - i + i_min) - p_xyz[0]) / (p_xyz[3] - p_xyz[0]);
		}

		if (alphax[i - i_min] < 0.0)
			alphax[i - i_min] = 0.0;
		else if (alphax[i - i_min] > 1)
			alphax[i - i_min] = 1.0;
	}

	for (int i = j_min; i <= j_max; i++) {
		if (p_xyz[4] == p_xyz[1]) {
			alphay[i - j_min] = 1;
			break;
		}
		if (p_xyz[4] > p_xyz[1]) {
			alphay[i - j_min] = (i - p_xyz[1]) / (p_xyz[4] - p_xyz[1]);
		}
		else {
			alphay[i - j_min] = ((j_max - i + j_min) - p_xyz[1]) / (p_xyz[4] - p_xyz[1]);
		}

		if (alphay[i - j_min] < 0.0)
			alphay[i - j_min] = 0.0;
		else if (alphay[i - j_min] > 1.0)
			alphay[i - j_min] = 1.0;

	}
	for (int i = k_min; i <= k_max; i++) {
		if (p_xyz[5] == p_xyz[2]) {
			alphaz[i - k_min] = 1;
			break;
		}
		if (p_xyz[5] > p_xyz[2]) {
			alphaz[i - k_min] = (i - p_xyz[2]) / (p_xyz[5] - p_xyz[2]);
		}
		else {
			alphaz[i - k_min] = ((k_max - i + k_min) - p_xyz[2]) / (p_xyz[5] - p_xyz[2]);
		}

		if (alphaz[i - k_min] < 0.0)
			alphaz[i - k_min] = 0.0;
		else if (alphaz[i - k_min] > 1.0)
			alphaz[i - k_min] = 1.0;
	}

	int NXY = NX + NY;
	float alphaxy[dimx + dimy+2];
	int ptrxy = 0, ptrx = 0, ptry = 0;
	while (ptrx < NX&&ptry < NY) {
		if (alphax[ptrx] <= alphay[ptry]) {
			alphaxy[ptrxy++] = alphax[ptrx++];
		}
		else {
			alphaxy[ptrxy++] = alphay[ptry++];
		}
	}
	while (ptrx < NX) alphaxy[ptrxy++] = alphax[ptrx++];
	while (ptry < NY) alphaxy[ptrxy++] = alphay[ptry++];

	//int NXYZ=NXY+NZ;
	float alpha[dimx + dimy + dimz+3];
	int ptr = 0, ptrz = 0;
	ptrxy = 0;
	while (ptrxy < NXY&&ptrz < NZ) {
		if (alphaxy[ptrxy] <= alphaz[ptrz]) {
			alpha[ptr++] = alphaxy[ptrxy++];
		}
		else {
			alpha[ptr++] = alphaz[ptrz++];
		}
	}
	while (ptrxy < NXY) alpha[ptr++] = alphaxy[ptrxy++];
	while (ptrz < NZ) alpha[ptr++] = alphaz[ptrz++];


	float L = sqrtf((p_xyz[3] - p_xyz[0]) * (p_xyz[3] - p_xyz[0]) + (p_xyz[4] - p_xyz[1]) * (p_xyz[4] - p_xyz[1])
		+ (p_xyz[5] - p_xyz[2]) * (p_xyz[5] - p_xyz[2]));
	int Maxlen = imageM * imageN * imageH-1;
	for (int i = 1; i < N; i++) {
		float a_mid = (alpha[i] + alpha[i - 1]) / 2;
		float w = (alpha[i] - alpha[i - 1]) * L;
		int x = floor(p_xyz[0] + a_mid * (p_xyz[3] - p_xyz[0]));
		int y = floor(p_xyz[1] + a_mid * (p_xyz[4] - p_xyz[1]));
		int z = floor(p_xyz[2] + a_mid * (p_xyz[5] - p_xyz[2]));

		int indx = z * imageM * imageN + y * imageN + x;
		if (indx < Maxlen)
		{
			proj_m[j] += image_3d[indx] * w;

		}
	}
	if (proj_m[j] < 0)
	{
		proj_m[j] = 0;
	}
}

bool Cprojection(float*image_3d, float *proj_m, float *_mask,float*GmtrcTrnsfrmtnMtrx, int imageM, int imageN, int imageH, float Pixeld, float Pixeldz, int projM, int projN, float Threshold)
{
	int BlockSize = 128;
	dim3 threads(BlockSize);
	int GridSize = projM * projN / BlockSize;
	dim3 blocks(GridSize);

	CUDAprojection <<< blocks, threads >>> (image_3d, proj_m, _mask, GmtrcTrnsfrmtnMtrx, imageM, imageN, imageH, Pixeld, Pixeldz,projM, projN, Threshold);
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) {
		printf("CUDA Error4: %s\n", cudaGetErrorString(err));
		// Possibly: exit(-1) if program cannot continue....
	}
	return true;
}

SiddonGPU::SiddonGPU() { 
	m_fThreshold = 0;
	m_fImg3d4Cuda = NULL;
	m_bPrepare3d = false;
}

SiddonGPU::~SiddonGPU()
{
	cudaFree(m_fImg3d4Cuda);
	cudaFree(m_fImg2d);
	cudaFree(m_fTransformMatrix);
	cudaFree(m_fImg2dMask4Cude);
	delete m_lImg3dPixelNumber;
	delete m_fImg3dPixelSpacing;
	if(_mask)
	{
		delete _mask;
	}
}



void SiddonGPU::SetImg3d(const float* _fimg3d, float* _PS, int* _PN)
{
	if (m_bPrepare3d)
		return;
	m_lImg3dPixelNumber = new int[3];
	memset(m_lImg3dPixelNumber, 0, 3 * sizeof(int));
	memcpy(m_lImg3dPixelNumber, _PN, 3 * sizeof(int));

	m_fImg3dPixelSpacing = new float[3];
	memset(m_fImg3dPixelSpacing, 0, 3 * sizeof(float));
	memcpy(m_fImg3dPixelSpacing, _PS, 3 * sizeof(float));

	int64_t len = m_lImg3dPixelNumber[0] * m_lImg3dPixelNumber[1] * m_lImg3dPixelNumber[2];
    m_fImg3d4Cuda = new float[len];
	cudaMalloc((void**)&m_fImg3d4Cuda, sizeof(short) * len);
	cudaMemcpy(m_fImg3d4Cuda, _fimg3d, sizeof(short) * len, cudaMemcpyHostToDevice);

	//savedata1(_fimg3d, len, "ct.raw");

	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) {
		printf("CUDA Error3: %s\n", cudaGetErrorString(error));
		// Possibly: exit(-1) if program cannot continue....
	}
	m_bPrepare3d = true;
}

void SiddonGPU::SetImg2dParameter(float* _PS, int* _PN)
{
	m_lImg2dPixelSpacing = new float[2];
	memset(m_lImg2dPixelSpacing, 0, 2 * sizeof(float));
	memcpy(m_lImg2dPixelSpacing, _PS, 2 * sizeof(float));

	m_lImg2dPixelNumber = new int[2];
	memset(m_lImg2dPixelNumber, 0, 2 * sizeof(int));
	memcpy(m_lImg2dPixelNumber, _PN, 2 * sizeof(int));

	int len = m_lImg2dPixelNumber[0] * m_lImg2dPixelNumber[1];
	m_fImg2d = new float[len];
	cudaMalloc((void**)&m_fImg2d, sizeof(float) * len);
	cudaMemset((void*)m_fImg2d, 0, sizeof(float) * len);

	m_fImg2dMask4Cude = new float[len];

	if (!_mask)
		_mask = new float[len];
	cudaMalloc((void**)&m_fImg2dMask4Cude, sizeof(float) * len);
	cudaMemcpy(m_fImg2dMask4Cude, _mask, sizeof(float) * len, cudaMemcpyHostToDevice);

	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) {
		printf("CUDA Error1: %s\n", cudaGetErrorString(error));
		// Possibly: exit(-1) if program cannot continue....
	}
}

void SiddonGPU::SetTransformMatrix(float* _fTransformMatrix)
{
	m_fTransformMatrix = new float[12]{0};
	cudaMalloc((void**)&m_fTransformMatrix, sizeof(float) * 12);
	cudaMemcpy(m_fTransformMatrix, _fTransformMatrix, sizeof(float) * 12, cudaMemcpyHostToDevice);
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) {
		printf("CUDA Error2: %s\n", cudaGetErrorString(error));
		// Possibly: exit(-1) if program cannot continue....
	}
}

bool SiddonGPU::Run(float* _fTransformMatrix,float* rst)
{
	//int BlockSize = 512;
	//dim3 threads(BlockSize);
	//int GridSize = m_lImg2dPixelNumber[0] * m_lImg2dPixelNumber[1] / BlockSize;
	//dim3 blocks(GridSize);
	cudaMemcpy(m_fTransformMatrix, _fTransformMatrix, sizeof(float) * 12, cudaMemcpyHostToDevice);
	//cudaError_t error = cudaGetLastError();
	//printf("CUDA error: %s\n", cudaGetErrorString(error));
	int len = m_lImg2dPixelNumber[0] * m_lImg2dPixelNumber[1];
	cudaMemset((void*)m_fImg2d, 0, sizeof(float) * len);
	Cprojection(m_fImg3d4Cuda, m_fImg2d, m_fImg2dMask4Cude, m_fTransformMatrix, m_lImg3dPixelNumber[1], m_lImg3dPixelNumber[0], m_lImg3dPixelNumber[2], m_fImg3dPixelSpacing[0], m_fImg3dPixelSpacing[2], m_lImg2dPixelNumber[0], m_lImg2dPixelNumber[1], m_fThreshold);
	cudaMemcpy(rst, m_fImg2d, sizeof(float) * len, cudaMemcpyDeviceToHost);
	return true;
}

