#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include <math.h>


// TODO: 
// Buttons to set perlin parameters to be used when we reset to perlin
// Toggle for dampening, and be able to set dampening term
// Basic display/debug infor showing perlin params, etc
// Be able to flag arbitrary cells as boundary cells and simulate arbitrary domains
// Toggle periodic boundary conditions vs mirror
// Click to interact with it (2D)
// Different render modes for 2D - gradients, component velocities, etc?
// Be able to render a 3D plane with pixelGameEngine
// Render heightmap with pixelGameEngine (3D)
// Switch between 2D and 3D views (menu system?)


void perlinNoise2D(int nRows, int nCols, float* fSeed, int nOctaves, float fBias, float* fOutput)
{
	/*
	Populates the array fOutput, which is a flattened matrix of size nRows x nCols, with 2D perlin noise
	fSeed should also be a flattened matrix of size nRows x nCols, populated with random noise
	
	Made referencing Javidx9's Perlin source under the GNU GPLv3 license:
	https://github.com/OneLoneCoder/videos/blob/master/OneLoneCoder_PerlinNoise.cpp
	*/
	for (int x = 0; x < nCols; x++)
	{
		for (int y = 0; y < nRows; y++)
		{
			float fNoise = 0.0f;
			float fScaleAcc = 0.0f;
			float fScale = 1.0f;

			for (int o = 0; o < nOctaves; o++)
			{
				int nPitch = nCols >> o;
				int nSampleX1 = (x / nPitch) * nPitch;
				int nSampleY1 = (y / nPitch) * nPitch;

				int nSampleX2 = (nSampleX1 + nPitch) % nCols;
				int nSampleY2 = (nSampleY1 + nPitch) % nCols;

				float fBlendX = (float)(x - nSampleX1) / (float)nPitch;
				float fBlendY = (float)(y - nSampleY1) / (float)nPitch;

				float fSampleT = (1.0f - fBlendX) * fSeed[nSampleY1 * nCols + nSampleX1] + fBlendX * fSeed[nSampleY1 * nCols + nSampleX2];
				float fSampleB = (1.0f - fBlendX) * fSeed[nSampleY2 * nCols + nSampleX1] + fBlendX * fSeed[nSampleY2 * nCols + nSampleX2];
			
				fScaleAcc += fScale;
				fNoise += (fBlendY * (fSampleB - fSampleT) + fSampleT) * fScale;
				fScale += fScale / fBias;
			}

			// Scale down by dividing by accumulated scaling factors
			fOutput[y * nCols + x] = fNoise / fScaleAcc;
		}
	}
}

class HeightField
	/*
	Heightfield approximation of a fluid surface
	See Sigraph 2007 course notes by Bridson and Muller-Fischer, Chapter 8:
	https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf
	*/
{
public:
	HeightField(int rows, int cols)
	{
		nRows = rows;
		nCols = cols;
		u = new float[nRows * nCols];
		v = new float[nRows * nCols];
		for (int i = 0; i < (nRows * nCols); i++)
		{
			u[i] = 1.0f;
			v[i] = 0;
		}
	}

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

	void setHeights(float* heights)
	{
		for (int i = 0; i < (nRows * nCols); i++)
		{
			u[i] = heights[i];
		}
	}

	void setVelocities(float* velocities)
	{
		for (int i = 0; i < (nRows * nCols); i++)
		{
			v[i] = velocities[i];
		}
	}

	void step(const float& fElapsedTime)
	{
		// BEGIN MIRROR BOUDNARY CONDITION
		// Calculate new velocities - boundary cells
		// Top left corner [0,0]
		v[0] += (u[0] + u[(0 + 1) * nCols] + u[0] + u[1]) / 4.0f - u[0];
		//v[0] *= 0.99f;
		// Top right corner [0,nCols - 1]
		v[nCols - 1] += (u[nCols - 1] + u[1 * nCols + (nCols - 1)] + u[nCols - 2] + u[nCols - 1]) / 4.0f - u[nCols - 1];
		//v[nCols - 1] *= 0.99f;
		// Bottom left corner [nRows - 1, 0]
		v[(nRows - 1) * nCols] += (u[(nRows - 2) * nCols] + u[(nRows - 1) * nCols] + u[(nRows - 1) * nCols] + u[(nRows - 1) * nCols + 1]) / 4.0f - u[(nRows - 1) * nCols];
		//v[(nRows - 1) * nCols] *= 0.99f;
		// Bottom right corner [nRows - 1, nCols - 1]
		v[(nRows - 1) * nCols + (nCols - 1)] += (u[(nRows - 2) * nCols + (nCols - 1)] + u[(nRows - 1) * nCols + (nCols - 1)] + u[(nRows - 1) * nCols + (nCols - 2)] + u[(nRows - 1) * nCols + (nCols - 1)]) / 4.0f - u[(nRows - 1) * nCols + (nCols - 1)];
		//v[(nRows - 1) * nCols + (nCols - 1)] *= 0.99f;
		// Top edge less corners, bottom edge less corners
		for (int j = 1; j < nCols - 1; j++)
		{
			v[j] += (u[j] + u[nCols + j] + u[j - 1] + u[j + 1]) / 4.0f - u[j];
			v[j] *= 0.99f;
			v[(nRows - 1) * nCols + j] += (u[(nRows - 2) * nCols + j] + u[(nRows - 1) * nCols + j] + u[(nRows - 1) * nCols + (j - 1)] + u[(nRows - 1) * nCols + (j + 1)]) / 4.0f - u[(nRows - 1) * nCols + j];
			//v[(nRows - 1) * nCols + j] *= 0.99f;
		}
		// Left edge less corners, right edge less corners
		for (int i = 1; i < nRows - 1; i++)
		{
			v[i * nCols] += (u[(i - 1) * nCols] + u[(i + 1) * nCols] + u[i * nCols] + u[i * nCols + 1]) / 4.0f - u[i * nCols];
			//v[i * nCols] *= 0.99f;
			v[i * nCols + (nCols - 1)] += (u[(i - 1) * nCols + (nCols - 1)] + u[(i + 1) * nCols + (nCols - 1)] + u[i * nCols + (nCols - 2)] + u[i * nCols + (nCols - 1)]) / 4.0f - u[i * nCols + (nCols - 1)];
			//v[i * nCols + (nCols - 1)] *= 0.99f;
		}
		// END MIRROR BOUNDARY CONDITION

		// Calculate new velocities - internal cells
		for (int i = 1; i < nRows - 1; i++)
		{
			for (int j = 1; j < nCols - 1; j++)
			{
				v[i * nCols + j] += (u[(i - 1) * nCols + j] + u[(i + 1) * nCols + j] + u[i * nCols + (j - 1)] + u[i * nCols + (j + 1)]) / 4.0f - u[i * nCols + j];
				//v[i * nCols + j] *= 0.99f; // dampen
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
			// u[i] += v[i]; // This gives more extreme and visually interesting results but is unphysical.
			u[i] += v[i] * fElapsedTime;
		}
	}

	float getHeight(const int& x, const int& y)
	{
		return u[y * nCols + x];
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
	float* u = nullptr; // height matrix, flattened

	// Removed for now because slower by doing checks on internal cells;
	// However leaving in because in the future we may want to do arbitrary domain;
	// in which case we will have to check if neighbors are in or out of bounds! 
	// (which we can do by having another matrix of flags, or bitpacked flags in a smaller array if we really want)
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

private:
	HeightField* hField = nullptr;
	float* initialHeights = nullptr;
	
	// perlin noise parameters
	float* fNoiseSeed = nullptr;
	int nOctaveCount = 8;
	float fScalingBias = 2.0f;

	bool OnUserCreate() override
	{
		int nRows = ScreenWidth();
		int nCols = ScreenHeight();
		fNoiseSeed = new float[nRows * nCols];
		for (int i = 0; i < nRows * nCols; i++) fNoiseSeed[i] = (float)rand() / (float)RAND_MAX;

		// Populate initial conditions of PDE with perlin noise
		initialHeights = new float[nRows * nCols];
		perlinNoise2D(nRows, nCols, fNoiseSeed, nOctaveCount, fScalingBias, initialHeights);

		hField = new HeightField(nRows, nCols, initialHeights);
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		int nRows = hField->getNRows();
		int nCols = hField->getNCols();

		if (GetKey(olc::Key::SPACE).bPressed) hField->setHeights(initialHeights); // reinitialize to perlin noise

		for (int x = 0; x < nCols; x++)
		{
			for (int y = 0; y < nRows; y++)
			{
				int t = std::clamp(hField->getHeight(x,y) * 255.0f, 0.0f, 255.0f);
				Draw(x, y, olc::Pixel(t/3, t/2, t));
			}
		}
		hField->step(fElapsedTime);
		return true;
	}
};

int main()
{
	PixelFluidEngine demo;
	if (demo.Construct(256, 256, 3, 3))
		demo.Start();

	return 0;
}