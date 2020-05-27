#include "HeightField.h"

HeightField::HeightField(int rows, int cols)
{
	/*
		Creates a rectangular uniform heightfield with no internal restrictions on domain
	*/
	nRows = rows;
	nCols = cols;
	u = new float[nRows * nCols];
	v = new float[nRows * nCols];
	bDomain = new bool[nRows * nCols];
	for (int i = 0; i < (nRows * nCols); i++)
	{
		u[i] = 1.0f;
		v[i] = 0;
		bDomain[i] = true;
	}
}

HeightField::HeightField(int rows, int cols, float* heights)
{
	nRows = rows;
	nCols = cols;
	u = new float[nRows * nCols];
	v = new float[nRows * nCols];
	bDomain = new bool[nRows * nCols];
	for (int i = 0; i < (nRows * nCols); i++)
	{
		u[i] = heights[i];
		v[i] = 0;
		bDomain[i] = true;
	}
}

HeightField::HeightField(int rows, int cols, float* heights, bool* domain)
{
	nRows = rows;
	nCols = cols;
	u = new float[nRows * nCols];
	v = new float[nRows * nCols];
	bDomain = new bool[nRows * nCols];
	for (int i = 0; i < (nRows * nCols); i++)
	{
		u[i] = heights[i];
		v[i] = 0;
		bDomain[i] = domain[i];
	}
}

void HeightField::setHeights(float* heights)
{
	for (int i = 0; i < (nRows * nCols); i++) u[i] = heights[i];
}

void HeightField::step(const float& fElapsedTime, const float& fDamp)
{
	// Calculate new velocities
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			v[i * nCols + j] += getVelocityChange(j, i);
			v[i * nCols + j] *= fDamp; // dampen
		}
	}
	// Update heights based on new velocities
	for (int i = 0; i < (nRows * nCols); i++)
	{
		u[i] += v[i] * fElapsedTime;
	}
}

void HeightField::setHeight(const int& x, const int& y, const float& fHeight)
{
	u[y * nCols + x] = fHeight;
}

void HeightField::setDomain(bool* domain)
{
	for (int i = 0; i < nRows * nCols; i++) bDomain = domain;
}

void HeightField::setDomainCell(const int& x, const int& y, const bool& b)
{
	if (x >= 0 && x < nCols && y >= 0 && y < nRows) bDomain[y * nCols + x] = b;
}

void HeightField::zeroVelocities()
{
	for (int i = 0; i < (nRows * nCols); i++) v[i] = 0;
}

void HeightField::clearDomain()
{
	for (int i = 0; i < (nRows * nCols); i++) bDomain[i] = true;
}

float HeightField::getHeight(const int& x, const int& y)
{
	return u[y * nCols + x];
}

bool HeightField::isInDomain(const int& x, const int& y)
{
	if (x < 0 || x >= nCols || y < 0 || y >= nRows) return false;
	return bDomain[y * nCols + x];
}

float HeightField::getVelocityChange(const int& x, const int& y)
{
	float eastHeight, westHeight, northHeight, southHeight; // heights of neighbors

	// Mirrored boundary conditions
	if (y == 0 || !bDomain[(y - 1) * nCols + x]) northHeight = u[y * nCols + x];
	else northHeight = u[(y - 1) * nCols + x];
	if (y == nRows - 1 || !bDomain[(y + 1) * nCols + x]) southHeight = u[y * nCols + x];
	else southHeight = u[(y + 1) * nCols + x];
	if (x == 0 || !bDomain[y * nCols + (x - 1)]) westHeight = u[y * nCols + x];
	else westHeight = u[y * nCols + (x - 1)];
	if (x == nCols - 1 || !bDomain[y * nCols + (x + 1)]) eastHeight = u[y * nCols + x];
	else eastHeight = u[y * nCols + (x + 1)];

	float result = (northHeight + southHeight + westHeight + eastHeight) / 4.0f - u[y * nCols + x];

	return result;
}