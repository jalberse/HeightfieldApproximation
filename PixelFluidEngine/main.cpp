#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

// TODO: 
// get it working just reading numbers
// Be able to render a 3D plane with pixelGameEngine
// Render heightmap with pixelGameEngine (3D)
// Be able to flag arbitrary cells as boundary cells and simulate arbitrary domains

class HeightField
{
public:
	float* u = nullptr; // height matrix, flattened

	HeightField(int rows, int cols, float* heights)
	{
		nRows = rows;
		nCols = cols;
		u = new float[nRows*nCols];
		v = new float[nRows*nCols];
		for (int i = 0; i < (nRows * nCols); i++)
		{
			u[i] = heights[i];
			v[i] = 0;
		}
	}

	//void step(const float& fElapsedTime)
	void step()
	{
		// Calculate new velocities based on surrounding, own heights
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				v[i * nCols + j] += getVelocityChange(i, j);
				v[i * nCols + j] *= 0.99f; // dampen
			}
		}
		// Update heights based on new velocities
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				u[i * nCols + j] += v[i * nCols + j];
				// u[i * nCols + j] += v[i * nCols + j] * fElapsedTime; 
			}
		}
	}

	float getHeight(const int& i, const int& j)
	{
		return u[i * nCols + j];
	}

	int getNCols()
	{
		return nCols;
	}

	int getNRows()
	{
		return nRows;
	}

private:
	int nRows;
	int nCols;
	float* v = nullptr; // velocity matrix, flattened

	float getVelocityChange(const int& i, const int& j)
	{
		
		float eastHeight, westHeight, northHeight, southHeight; // heights of neighbors

		// Mirrored boundary conditions
		if (i == 0) northHeight = u[j];
		else northHeight = u[(i - 1) * nCols + j];
		if (i == nRows - 1) southHeight = u[i * nCols + j];
		else southHeight = u[(i + 1) * nCols + j];
		if (j == 0) westHeight = u[i * nCols + j];
		else westHeight = u[i * nCols + (j - 1)];
		if (j == nCols - 1) eastHeight = u[i * nCols + j];
		else eastHeight = u[i * nCols + (j + 1)];

		return (northHeight + southHeight + westHeight + eastHeight) / 4.0f - u[i * nCols + j];
	}
};



// Override base class with your custom functionality
class PixelFluidEngine : public olc::PixelGameEngine
{
public:
	PixelFluidEngine()
	{
		// Name you application
		sAppName = "Pixel Fluid Engine";
	}

public:
	HeightField* hField = nullptr;
	float* initialHeights = nullptr;

	bool OnUserCreate() override
	{
		int nCols = 10;
		int nRows = 10;
		initialHeights = new float[nRows * nCols];
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				initialHeights[i * nCols + j] = (float)(i * (100 / nRows));
			}
		}
		hField = new HeightField(nRows, nCols, initialHeights);
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
				
		return true;
	}
};

int main()
{
	/*PixelFluidEngine demo;
	if (demo.Construct(256, 240, 4, 4))
		demo.Start();*/

	int nRows = 10;
	int nCols = 10;

	float* initialHeights = new float[nRows * nCols];
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			initialHeights[i * nCols + j] = (float)((i+1) * (100 / nRows));
		}
	}
	HeightField* hField = new HeightField(10, 10, initialHeights);

	for (int k = 0; k < 15; k++)
	{
		hField->step();
	}

	return 0;
}