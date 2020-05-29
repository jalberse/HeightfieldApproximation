/*
	What is this?
	~~~~~~~~~~~~~

	This software is a tool for visualizing fluid simulation algorithms for computer graphics applications.
	This is built on top of the OneLoneCoder PixelGameEngine: https://github.com/OneLoneCoder/olcPixelGameEngine/wiki

	LICENSE (MIT)
	~~~~~~~~~~~~~
	
	Copyright 2020 John Alberse

	Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
	to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
	and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
	FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
	WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

	Author
	~~~~~~

	@author John Alberse
	Contact: alberse.john@gmail.com
*/

#define OLC_PGE_APPLICATION

#include "HeightField.h"
#include "PerlinNoise.h"

#include "olcPixelGameEngine.h"

#include <math.h>

// TODO: 
// GPU support
// A button to make it start/stop "raining" would be fun
// Buttons to change parameters/clean up UI generally. On own branch. Keep hotkeys for "expert" (my) use
// Toggle periodic boundary conditions vs mirror (?)
// Different render modes for 2D - gradients, component velocities, etc?
// Be able to render a 3D plane with pixelGameEngine
// Render heightmap with pixelGameEngine (3D)
// Switch between 2D and 3D views (menu system?)
// More complex heightfield model which accounts for grid spacing h and wave speed c
// Once in 3D, add rudimentary interactions with particles - if particle collides with column, 
//	add height to column and delete the particle

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
	float* fTerrain = nullptr;
	int nRows = 256;
	int nCols = 256;
	int nMouseX;
	int nMouseY;
	bool paused = true;
	float fDampMax = 1.0f;
	float fDamp = 0.999f;
	float fDampMin = 0.98f;
	float fDampStep = 0.001;
	
	// perlin noise parameters - water
	float* fNoiseSeed = nullptr;
	int nOctaveMax = 8;
	int nOctave = 8;
	float fScalingBiasMin = 0.2f;
	float fScalingBias = 2.0f;
	float fScalingBiasStep = 0.2f;

	// perlin noise parameters - terrain
	int nOctaveTerrain = 4;
	float fScalingBiasTerrain = 2.0f;
	
	// Used to theshold terrain to generate domain of water
	// If terrain is above the fluid level, we take it out of the domain
	// If the terrain is below the fluid level, it is in the domain
	float fFluidLevel = 0.5f; 
	float fFluidLevelStep = 0.05f;
	float fFluidLevelMin = 0.1f;
	float fFluidLevelMax = 1.0f;

	bool OnUserCreate() override
	{
		fNoiseSeed = new float[nRows * nCols];
		initialHeights = new float[nRows * nCols];

		// Populate initial conditions of PDE with perlin noise
		for (int i = 0; i < nRows * nCols; i++) fNoiseSeed[i] = (float)rand() / (float)RAND_MAX;
		perlinNoise2D(nRows, nCols, fNoiseSeed, nOctave, fScalingBias, initialHeights);
		
		hField = new HeightField(nRows, nCols, initialHeights);

		// Generate some initial terrain and use it to set domain of heightfield
		fTerrain = new float[nRows * nCols];
		generateTerrain();
		applyTerrain();

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		handleInput(fElapsedTime);

		Clear(olc::Pixel(0,0,0));
		for (int x = 0; x < nCols; x++)
		{
			for (int y = 0; y < nRows; y++)
			{
				float height = hField->getHeight(x, y) * 255.0f;
				int t = height <= 0.0f ? 0.0f : height <= 255.0f ? height : 255.0f;
				Draw(x, y, olc::Pixel(t/3, t/2, t));
				if (!hField->isInDomain(x, y)) Draw(x, y, olc::Pixel(0, 0, 0));
				
			}
		}
		drawUI();
		
		if (!paused) hField->step(fElapsedTime);
		
		return true;
	}

	void handleInput(float fElapsedTime)
	{
		// Reset heights to perlin noise with new random seed
		if (GetKey(olc::Key::R).bReleased)
		{
			if (GetKey(olc::Key::SHIFT).bHeld) {
				generateTerrain(); applyTerrain();
			}
			generateFluidSurface(); applyFluidSurface();
		}
		if (GetKey(olc::Key::X).bReleased && fFluidLevel + fFluidLevelStep <= fFluidLevelMax)
		{
			fFluidLevel += fFluidLevelStep;
			applyTerrain();
			generateFluidSurface(); applyFluidSurface();

		}
		if (GetKey(olc::Key::Z).bReleased && fFluidLevel - fFluidLevelStep >= fFluidLevelMin)
		{
			fFluidLevel -= fFluidLevelStep;
			applyTerrain();
			generateFluidSurface(); applyFluidSurface();
		}
		// Set various perlin noise params
		if (GetKey(olc::Key::P).bReleased) {
			if (GetKey(olc::Key::SHIFT).bHeld) (nOctaveTerrain == nOctaveMax) ? nOctaveTerrain = 1 : nOctaveTerrain++;
			else (nOctave == nOctaveMax) ? nOctave = 1 : nOctave++;
		}
		if (GetKey(olc::Key::O).bReleased)
		{
			if (GetKey(olc::Key::SHIFT).bHeld) (nOctaveTerrain == 1) ? nOctaveTerrain = nOctaveMax : nOctaveTerrain--;
			else (nOctave == 1) ? nOctave = nOctaveMax : nOctave--;
		}
		if (GetKey(olc::Key::L).bReleased)
		{
			if (GetKey(olc::Key::SHIFT).bHeld) fScalingBiasTerrain += fScalingBiasStep;
			else fScalingBias += fScalingBiasStep;
		}
		if (GetKey(olc::Key::K).bReleased)
		{
			if (GetKey(olc::Key::SHIFT).bHeld) {
				if (fScalingBiasTerrain >= fScalingBiasMin + fScalingBiasStep) fScalingBiasTerrain -= fScalingBiasStep;
			}
			else if (fScalingBias >= fScalingBiasMin + fScalingBiasStep) fScalingBias -= fScalingBiasStep;
		}
		if (GetKey(olc::Key::M).bReleased && fDamp <= fDampMax - fDampStep) {
			fDamp += fDampStep;
		}
		if (GetKey(olc::Key::N).bReleased && fDamp >= fDampMin + fDampStep) {
			fDamp -= fDampStep;
		}
		if (GetKey(olc::Key::SPACE).bReleased) paused = !paused;
		if (GetKey(olc::Key::C).bReleased) hField->clearDomain();
		if (GetKey(olc::Key::S).bHeld && paused) hField->step(fElapsedTime);

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
				hField->setDomainCell(nMouseX, nMouseY, GetKey(olc::Key::SHIFT).bHeld);
			}
		}
	}

	void drawUI()
	{
		DrawString(260, 50, "Dampening      : " + std::to_string(fDamp));
		DrawString(260, 70, "Fluid level    : " + std::to_string(fFluidLevel));
		DrawString(260, 90, "Fluid Perlin Parameters");
		DrawString(260, 110, "  Octave       : " + std::to_string(nOctave));
		DrawString(260, 130, "  Scaling Bias : " + std::to_string(fScalingBias));
		DrawString(260, 150, "Terrain Perlin Parameters");
		DrawString(260, 170, "  Octave       : " + std::to_string(nOctaveTerrain));
		DrawString(260, 190, "  Scaling Bias : " + std::to_string(fScalingBiasTerrain));
	}

	void generateTerrain()
	{
		for (int i = 0; i < nRows * nCols; i++) fNoiseSeed[i] = (float)rand() / (float)RAND_MAX;
		perlinNoise2D(nRows, nCols, fNoiseSeed, nOctaveTerrain, fScalingBiasTerrain, fTerrain);
		
	}

	void applyTerrain()
	{
		for (int i = 0; i < nRows * nCols; i++) hField->setDomainCell(i % nCols, i / nCols, fTerrain[i] < fFluidLevel);
	}

	void generateFluidSurface()
	{
		for (int i = 0; i < nRows * nCols; i++) fNoiseSeed[i] = (float)rand() / (float)RAND_MAX;
		perlinNoise2D(nRows, nCols, fNoiseSeed, nOctave, fScalingBias, initialHeights);
	}

	void applyFluidSurface()
	{
		hField->setHeights(initialHeights);
		hField->zeroVelocities();
	}
};

int main()
{
	PixelFluidEngine demo;
	if (demo.Construct(1920 / 4, 1080 / 4, 3, 3))
		demo.Start();

	return 0;
}