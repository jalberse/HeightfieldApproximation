#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include <math.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

// TODO: 
// MAKE IT FASTERRRR IT SUCKS SO BAD RIGHT NOW!!!
//		I don't know how to make it better besides better bounds checking though. 
//		Maybe just won't do realtime since I won't be able to push it to GPU anyways
//		Other algorithms I definitely can't do realtime anyways, so... might as well just go compute->save->load->render route
//
//		BUT I'm also confused because it seems to be running very fast when I'm just like, stepping in main() without the olcEngine class?
//		Maybe we can write frames to a vector of frames? And we can scrub through?
// Toggle periodic boundary conditions vs mirror
// Render 2D by just coloring 0-255
// Be able to render a 3D plane with pixelGameEngine
// Render heightmap with pixelGameEngine (3D)
// Be able to flag arbitrary cells as boundary cells and simulate arbitrary domains

class HeightField
	/*
	Heightfield approximation of a fluid surface
	See Sigraph 2007 course notes by Bridson and Muller-Fischer, Chapter 8:
	https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf
	*/
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

	void step(const float& fElapsedTime)
	{
		// BEGIN MIRROR BOUDNARY CONDITION
		// Calculate new velocities - boundary cells
		// Top left corner [0,0]
		v[0] += (u[0] + u[(0 + 1) * nCols] + u[0] + u[1]) / 4.0f - u[0];
		v[0] *= 0.99f;
		// Top right corner [0,nCols - 1]
		v[nCols - 1] += (u[nCols - 1] + u[1 * nCols + (nCols - 1)] + u[nCols - 2] + u[nCols - 1]) / 4.0f - u[nCols - 1];
		v[nCols - 1] *= 0.99f;
		// Bottom left corner [nRows - 1, 0]
		v[(nRows - 1) * nCols] += (u[(nRows - 2) * nCols] + u[(nRows - 1) * nCols] + u[(nRows - 1) * nCols] + u[(nRows - 1) * nCols + 1]) / 4.0f - u[(nRows - 1) * nCols];
		v[(nRows - 1) * nCols] *= 0.99f;
		// Bottom right corner [nRows - 1, nCols - 1]
		v[(nRows - 1) * nCols + (nCols - 1)] += (u[(nRows - 2) * nCols + (nCols - 1)] + u[(nRows - 1) * nCols + (nCols - 1)] + u[(nRows - 1) * nCols + (nCols - 2)] + u[(nRows - 1) * nCols + (nCols - 1)]) / 4.0f - u[(nRows - 1) * nCols + (nCols - 1)];
		v[(nRows - 1) * nCols + (nCols - 1)] *= 0.99f;
		// Top edge less corners, bottom edge less corners
		for (int j = 1; j < nCols - 1; j++)
		{
			v[j] += (u[j] + u[nCols + j] + u[j - 1] + u[j + 1]) / 4.0f - u[j];
			v[j] *= 0.99f;
			v[(nRows - 1) * nCols + j] += (u[(nRows - 2) * nCols + j] + u[(nRows - 1) * nCols + j] + u[(nRows - 1) * nCols + (j - 1)] + u[(nRows - 1) * nCols + (j + 1)]) / 4.0f - u[(nRows - 1) * nCols + j];
			v[(nRows - 1) * nCols + j] *= 0.99f;
		}
		// Left edge less corners, right edge less corners
		for (int i = 1; i < nRows - 1; i++)
		{
			v[i * nCols] += (u[(i - 1) * nCols] + u[(i + 1) * nCols] + u[i * nCols] + u[i * nCols + 1]) / 4.0f - u[i * nCols];
			v[i * nCols] *= 0.99f;
			v[i * nCols + (nCols - 1)] += (u[(i - 1) * nCols + (nCols - 1)] + u[(i + 1) * nCols + (nCols - 1)] + u[i * nCols + (nCols - 2)] + u[i * nCols + (nCols - 1)]) / 4.0f - u[i * nCols + (nCols - 1)];
			v[i * nCols + (nCols - 1)] *= 0.99f;
		}
		// END MIRROR BOUNDARY CONDITION

		// Calculate new velocities - internal cells
		for (int i = 1; i < nRows - 1; i++)
		{
			for (int j = 1; j < nCols - 1; j++)
			{
				v[i * nCols + j] += (u[(i - 1) * nCols + j] + u[(i + 1) * nCols + j] + u[i * nCols + (j - 1)] + u[i * nCols + (j + 1)]) / 4.0f - u[i * nCols + j];
				v[i * nCols + j] *= 0.99f; // dampen
			}
		}
		//for (int i = 0; i < nRows; i++)
		//{
		//	for (int j = 1; j < nCols; j++)
		//	{
		//		v[i * nCols + j] += getVelocityChange(i, j);
		//		v[i * nCols + j] *= 0.99f; // dampen
		//	}
		//}
		// Update heights based on new velocities
		for (int i = 0; i < (nRows * nCols); i++)
		{
			//u[i] += v[i];
			u[i] += v[i] * fElapsedTime;
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

	//float getVelocityChange(const int& i, const int& j)
	//{
	//	
	//	float eastHeight, westHeight, northHeight, southHeight; // heights of neighbors

	//	// Mirrored boundary conditions
	//	// TODO can improve performance drastically here by putting this all in line manually so we don't check indices for all internal ones
	//	if (i == 0) northHeight = u[j];
	//	else northHeight = u[(i - 1) * nCols + j];
	//	if (i == nRows - 1) southHeight = u[i * nCols + j];
	//	else southHeight = u[(i + 1) * nCols + j];
	//	if (j == 0) westHeight = u[i * nCols + j];
	//	else westHeight = u[i * nCols + (j - 1)];
	//	if (j == nCols - 1) eastHeight = u[i * nCols + j];
	//	else eastHeight = u[i * nCols + (j + 1)];

	//	float result = (northHeight + southHeight + westHeight + eastHeight) / 4.0f - u[i * nCols + j];

	//	return result;
	//}
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
		int nRows = ScreenWidth();
		int nCols = ScreenHeight();
		initialHeights = new float[nRows * nCols];
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				float dist = sqrt((i - nRows / 2) * (i - nRows / 2) + (j - nCols / 2) * (j - nCols / 2));
				initialHeights[i * nCols + j] = dist * 2;
			}
		}
		hField = new HeightField(nRows, nCols, initialHeights);
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		Sleep(10);
		int nRows = hField->getNRows();
		int nCols = hField->getNCols();

		for (int x = 0; x < nCols; x++)
		{
			for (int y = 0; y < nRows; y++)
			{
				int t = hField->u[y * nCols + x];
				Draw(x, y, olc::Pixel(t, 0, 0));
			}
		}
		hField->step(fElapsedTime);
		return true;
	}
};

int main()
{
	PixelFluidEngine demo;
	if (demo.Construct(100, 100, 10, 10))
		demo.Start();

	/*int nRows = 100;
	int nCols = 100;

	float* initialHeights = new float[nRows * nCols];
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			initialHeights[i * nCols + j] = (float)((i+1) * (255 / nRows));
		}
	}
	HeightField* hField = new HeightField(nRows, nCols, initialHeights);


	for (int k = 0; k < 10000; k++)
	{
		for (int x = 0; x < nCols; x++)
		{
			for (int y = 0; y < nRows; y++)
			{
				int t = hField->u[y * nCols + x];
				hField->step(0.01f);
			}
		}
	}*/

	return 0;
}