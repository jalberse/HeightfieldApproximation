/*
	What is this?
	~~~~~~~~~~~~~

	This software is a tool for visualizing fluid simulation algorithms for computer graphics applications.
	This is built on top of the OneLoneCoder PixelGameEngine: https://github.com/OneLoneCoder/olcPixelGameEngine/wiki


	LICENSE (GNU GPL-3.0-or-later)
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	Copyright 2020 John Alberse


	This software is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This software is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this software.  If not, see <https://www.gnu.org/licenses/>.


	Author
	~~~~~~

	@author John Alberse
	Contact: alberse.john@gmail.com
*/

#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include <math.h>


// TODO: 
// Move heightfield to its own file
// Use perlin noise and thresholding to generate terrain
// Set "height" of water by changing threshold on terrain
// Buttons to change parameters/clean up UI generally. On own branch. Keep hotkeys for "expert" (my) use
// Why is, on octave 1, regens always increasing?
// Toggle periodic boundary conditions vs mirror (?)
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

	HeightField(int rows, int cols, float* heights)
	{
		nRows = rows;
		nCols = cols;
		u = new float[nRows*nCols];
		v = new float[nRows*nCols];
		bDomain = new bool[nRows * nCols];
		for (int i = 0; i < (nRows * nCols); i++)
		{
			u[i] = heights[i];
			v[i] = 0;
			bDomain[i] = true;
		}
	}

	HeightField(int rows, int cols, float* heights, bool* domain)
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

	void setHeights(float* heights)
	{
		for (int i = 0; i < (nRows * nCols); i++) u[i] = heights[i];
	}

	void step(const float& fElapsedTime, const float& fDamp = 1.0f)
	{
		// Calculate new velocities
		for (int i = 1; i < nRows - 1; i++)
		{
			for (int j = 1; j < nCols - 1; j++)
			{
				v[i * nCols + j] += getVelocityChange(j, i);
				v[i * nCols + j] *= fDamp; // dampen
			}
		}
		// Update heights based on new velocities
		for (int i = 0; i < (nRows * nCols); i++)
		{
			// u[i] += v[i]; // This gives more extreme and visually interesting results but is unphysical.
			u[i] += v[i] * fElapsedTime;
		}
	}

	void setHeight(const int& x, const int& y, const float& fHeight)
	{
		u[y * nCols + x] = fHeight;
	}

	void setDomain(bool* domain)
	{
		for (int i = 0; i < nRows * nCols; i++) bDomain = domain;
	}

	void setDomainCell(const int& x, const int& y, const bool& b)
	{
		if (x >= 0 && x < nCols && y >= 0 && y < nRows) bDomain[y * nCols + x] = b;
	}

	void zeroVelocities()
	{
		for (int i = 0; i < (nRows * nCols); i++) v[i] = 0;
	}

	void clearDomain()
	{
		for (int i = 0; i < (nRows * nCols); i++) bDomain[i] = true;
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

	bool isInDomain(const int& x, const int& y)
	{
		if (x < 0 || x >= nCols || y < 0 || y >= nRows) return false;
		return bDomain[y * nCols + x];
	}

private:
	int nRows;
	int nCols;
	float* v = nullptr; // velocity matrix, flattened
	float* u = nullptr; // height matrix, flattened
	bool* bDomain = nullptr; // bool representing fluid domain; 1 if can flow, 0 if not 

	float getVelocityChange(const int& x, const int& y)
	{
		
		float eastHeight, westHeight, northHeight, southHeight; // heights of neighbors

		// Mirrored boundary conditions
		if (y == 0 || !bDomain[(y - 1) * nCols + x]) northHeight = u[x];
		else northHeight = u[(y - 1) * nCols + x]; 
		if (y == nRows - 1 || !bDomain[(y+1) * nCols + x]) southHeight = u[y * nCols + x];
		else southHeight = u[(y + 1) * nCols + x];
		if (x == 0 || !bDomain[y * nCols + (x - 1)]) westHeight = u[y * nCols + x];
		else westHeight = u[y * nCols + (x - 1)];
		if (x == nCols - 1 || !bDomain[y * nCols + (x + 1)]) eastHeight = u[y * nCols + x];
		else eastHeight = u[y * nCols + (x + 1)];

		float result = (northHeight + southHeight + westHeight + eastHeight) / 4.0f - u[y * nCols + x];

		return result;
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

private:
	HeightField* hField = nullptr;
	float* initialHeights = nullptr;
	float nRows = 256;
	float nCols = 256;
	int nMouseX;
	int nMouseY;
	bool paused = false;
	float fDampMax = 1.0f;
	float fDamp = 1.0f;
	float fDampMin = 0.98f;
	float fDampStep = 0.001;
	bool bDomainModMode = false; // true if removing terrain (adding to fluid domain), false if removing terrain
	
	// perlin noise parameters
	float* fNoiseSeed = nullptr;
	int nOctaveMax = 8;
	int nOctave = 8;
	float fScalingBiasMin = 0.2f;
	float fScalingBias = 2.0f;
	float fScalingBiasStep = 0.2f;
	

	bool OnUserCreate() override
	{
		fNoiseSeed = new float[nRows * nCols];
		for (int i = 0; i < nRows * nCols; i++) fNoiseSeed[i] = (float)rand() / (float)RAND_MAX;

		// Populate initial conditions of PDE with perlin noise
		initialHeights = new float[nRows * nCols];
		perlinNoise2D(nRows, nCols, fNoiseSeed, nOctave, fScalingBias, initialHeights);

		hField = new HeightField(nRows, nCols, initialHeights);
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		// Reset heights to perlin noise with new random seed, and reset terrain
		if (GetKey(olc::Key::R).bReleased)
		{
			for (int i = 0; i < nRows * nCols; i++) fNoiseSeed[i] = (float)rand() / (float)RAND_MAX;
			perlinNoise2D(nRows, nCols, fNoiseSeed, nOctave, fScalingBias, initialHeights);
			hField->setHeights(initialHeights); 
			hField->zeroVelocities();
			hField->clearDomain();
		}

		// Set various perlin noise params
		if (GetKey(olc::Key::P).bReleased) (nOctave == nOctaveMax) ? nOctave = 1 : nOctave++;
		if (GetKey(olc::Key::O).bReleased) (nOctave == 1) ? nOctave = nOctaveMax : nOctave--;
		if (GetKey(olc::Key::L).bReleased) fScalingBias += fScalingBiasStep;
		if (GetKey(olc::Key::K).bReleased) if (fScalingBias >= fScalingBiasMin + fScalingBiasStep) fScalingBias -= fScalingBiasStep;
		if (GetKey(olc::Key::M).bReleased) if (fDamp <= fDampMax - fDampStep) fDamp += fDampStep;
		if (GetKey(olc::Key::N).bReleased) if (fDamp >= fDampMin + fDampStep) fDamp -= fDampStep;
		if (GetKey(olc::Key::T).bReleased) bDomainModMode = !bDomainModMode;
		if (GetKey(olc::Key::SPACE).bReleased) paused = !paused;

		nMouseX = GetMouseX();
		nMouseY = GetMouseY();
		if (GetMouse(0).bHeld)
		{
			if (nMouseX >= 0 && nMouseX < nCols && nMouseY >= 0 && nMouseY < nRows)
			{
				hField->setHeight(nMouseX, nMouseY, 1.0f);
			}
		}
		if (GetMouse(1).bHeld)
		{
			if (nMouseX >= 0 && nMouseX < nCols && nMouseY >= 0 && nMouseY < nRows)
			{
				hField->setDomainCell(nMouseX, nMouseY, bDomainModMode);
			}
		}

		Clear(olc::Pixel(0,0,0));

		for (int x = 0; x < nCols; x++)
		{
			for (int y = 0; y < nRows; y++)
			{
				int t = std::clamp(hField->getHeight(x,y) * 255.0f, 0.0f, 255.0f);
				Draw(x, y, olc::Pixel(t/3, t/2, t));
				if (!hField->isInDomain(x, y)) Draw(x, y, olc::Pixel(0, 0, 0));
				
			}
		}

		DrawString(260, 10, "SPACE to pause");
		DrawString(260, 30, "Click to perturb");
		if (!bDomainModMode) DrawString(260, 50, "RClick to add terrain");
		else DrawString(260, 50, "RClick to remove terrain");
		DrawString(260, 70, "(T to toggle terrain mode)");
		DrawString(260, 90, "R to reset");
		DrawString(260, 110, "Octave [O,P]: " + std::to_string(nOctave));
		DrawString(260, 130, "Scaling Bias [K,L]: " + std::to_string(fScalingBias));
		DrawString(260, 150, "Dampening [N,M]: " + std::to_string(fDamp));

		if (!paused) hField->step(fElapsedTime, fDamp);
		
		return true;
	}
};

int main()
{
	PixelFluidEngine demo;
	if (demo.Construct(512, 256, 3, 3))
		demo.Start();

	return 0;
}