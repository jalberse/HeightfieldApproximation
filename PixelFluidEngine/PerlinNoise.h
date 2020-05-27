#pragma once

void perlinNoise2D(int nRows, int nCols, float* fSeed, int nOctaves, float fBias, float* fOutput)
{
	/*
	Populates the array fOutput, which is a flattened matrix of size nRows x nCols, with 2D perlin noise [0,1]
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