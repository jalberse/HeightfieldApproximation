#pragma once

class HeightField
{
	/*
	Heightfield approximation of a fluid surface
	See Sigraph 2007 course notes by Bridson and Muller-Fischer, Chapter 8:
	https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf
	*/
public:
	HeightField(int rows, int cols);
	HeightField(int rows, int cols, float* heights);
	HeightField(int rows, int cols, float* heights, bool* domain);

	void step(const float& fElapsedTime, const float& fDamp = 1.0f);

	void setHeights(float* heights);
	void setHeight(const int& x, const int& y, const float& fHeight);
	void setDomain(bool* domain);
	void setDomainCell(const int& x, const int& y, const bool& b);
	void zeroVelocities();
	void clearDomain();
	float getHeight(const int& x, const int& y);
	int getNCols() { return nCols; }
	int getNRows() { return nRows; }
	bool isInDomain(const int& x, const int& y);

private:
	int nRows;
	int nCols;
	float* v = nullptr; // velocity matrix, flattened
	float* u = nullptr; // height matrix, flattened
	bool* bDomain = nullptr; // bool representing fluid domain; 1 if can flow, 0 if not 

	float getVelocityChange(const int& x, const int& y);
};