#pragma once

#include "thrust/host_vector.h"
#include "thrust/device_vector.h"

class HeightField
{
	/*
	Heightfield approximation of a fluid surface
	See Sigraph 2007 course notes by Bridson and Muller-Fischer, Chapter 8:
	https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf
	*/
public:
	HeightField(int rows, int cols, float* heights);

	void step(const float& fElapsedTime);

	void setHeights(float* heights);
	void setHeight(const int& x, const int& y, const float& fHeight);
	void setDomain(bool* domain);
	void setDomainCell(const int& x, const int& y, const bool& b);
	void setFDamp(const float& damp) { fDamp = damp; }
	void zeroVelocities();
	void clearDomain();
	float getHeight(const int& x, const int& y);
	int getNCols() { return nCols; }
	int getNRows() { return nRows; }
	bool isInDomain(const int& x, const int& y);

private:
	int nRows;
	int nCols;
	float fDamp;

	thrust::host_vector<float> h_z;
	thrust::device_vector<float> d_z;
	thrust::host_vector<float> h_dz;
	thrust::device_vector<float> d_dz;
	thrust::device_vector<float> d_ddz;
	thrust::host_vector<bool> h_bDomain;
	thrust::device_vector<bool> d_bDomain;
};