#include "HeightField.h"

HeightField::HeightField(int rows, int cols)
{
	/*
		Creates a rectangular uniform heightfield with no internal restrictions on domain
	*/
	nRows = rows;
	nCols = cols;
	z = new float[nRows * nCols];
	dz = new float[nRows * nCols];
	bDomain = new bool[nRows * nCols];
	for (int i = 0; i < (nRows * nCols); i++)
	{
		z[i] = 1.0f;
		dz[i] = 0;
		bDomain[i] = true;
	}
}

HeightField::HeightField(int rows, int cols, float* heights)
{
	nRows = rows;
	nCols = cols;
	z = new float[nRows * nCols];
	dz = new float[nRows * nCols];
	bDomain = new bool[nRows * nCols];
	for (int i = 0; i < (nRows * nCols); i++)
	{
		z[i] = heights[i];
		dz[i] = 0;
		bDomain[i] = true;
	}
}

HeightField::HeightField(int rows, int cols, float* heights, bool* domain)
{
	nRows = rows;
	nCols = cols;
	z = new float[nRows * nCols];
	dz = new float[nRows * nCols];
	bDomain = new bool[nRows * nCols];
	for (int i = 0; i < (nRows * nCols); i++)
	{
		z[i] = heights[i];
		dz[i] = 0;
		bDomain[i] = domain[i];
	}
}

void HeightField::setHeights(float* heights)
{
	for (int i = 0; i < (nRows * nCols); i++) z[i] = heights[i];
}

void HeightField::step(const float& fElapsedTime, const float& fDamp)
{
	// Calculate new velocities
	for (int y = 0; y < nRows; y++)
	{
		for (int x = 0; x < nCols; x++)
		{
			dz[y * nCols + x] += getVelocityChange(x, y);
			dz[y * nCols + x] *= fDamp; // dampen
		}
	}
	// Update heights based on new velocities
	for (int i = 0; i < (nRows * nCols); i++)
	{
		z[i] += dz[i] * fElapsedTime;
	}
}

void HeightField::setHeight(const int& x, const int& y, const float& fHeight)
{
	z[y * nCols + x] = fHeight;
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
	for (int i = 0; i < (nRows * nCols); i++) dz[i] = 0;
}

void HeightField::clearDomain()
{
	for (int i = 0; i < (nRows * nCols); i++) bDomain[i] = true;
}

float HeightField::getHeight(const int& x, const int& y)
{
	return z[y * nCols + x];
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
	if (y == 0 || !bDomain[(y - 1) * nCols + x]) northHeight = z[y * nCols + x];
	else northHeight = z[(y - 1) * nCols + x];
	if (y == nRows - 1 || !bDomain[(y + 1) * nCols + x]) southHeight = z[y * nCols + x];
	else southHeight = z[(y + 1) * nCols + x];
	if (x == 0 || !bDomain[y * nCols + (x - 1)]) westHeight = z[y * nCols + x];
	else westHeight = z[y * nCols + (x - 1)];
	if (x == nCols - 1 || !bDomain[y * nCols + (x + 1)]) eastHeight = z[y * nCols + x];
	else eastHeight = z[y * nCols + (x + 1)];

	float result = (northHeight + southHeight + westHeight + eastHeight) / 4.0f - z[y * nCols + x];

	return result;
}