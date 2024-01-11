#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "DRR_GPU_Export.h"

class CUDA_EXPORT SiddonGPU
{
public:
	SiddonGPU();
	~SiddonGPU();

	void SetImg3d(const float* _fimg3d, float* _PS, int* _PN);
	void SetImg2dParameter(float* _PS, int* _PN);
	void SetTransformMatrix(float* _fTransformMatrix);
	bool Run(float* _fTransformMatrix, float* rst);
private:
	float* m_fImg3d;
	float* m_fImg2d;
	float* m_fTransformMatrix;

	float* m_fImg3dPixelSpacing;
	float* m_lImg2dPixelSpacing;
	int* m_lImg3dPixelNumber;
	int* m_lImg2dPixelNumber;

	float* m_fImg3d4Cuda;
	float* m_fImg2d4Cuda;
	float* m_fImg2dMask4Cude;
	float m_fThreshold;
	bool m_bPrepare3d;
	float* _mask{ nullptr };

};