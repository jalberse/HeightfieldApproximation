#include "HeightField.h"

#include "thrust/host_vector.h"
#include "thrust/device_vector.h"
#include "thrust/iterator/constant_iterator.h"

HeightField::HeightField(int rows, int cols, float* heights) : 
	h_z(heights, heights + (rows * cols)), 
	h_dz(rows * cols, 0.0f), 
	h_bDomain(rows * cols, true),
	d_z(rows * cols),
	d_dz(rows * cols),
	d_ddz(rows * cols),
	d_bDomain(rows * cols),
	fDamp(0.999f)
{
	nRows = rows;
	nCols = cols;
	for (int y = 0; y < nRows; y++)
	{
		for (int x = 0; x < nCols; x++)
		{
			if (y == 0 || x == 0 || y == nRows - 1 || x == nCols - 1) h_bDomain[y * nCols + x] = false;
			else h_bDomain[y * nCols + x] = true;
			h_dz[y * nCols + x] = 0;
		}
	}
}

struct get_velocity_change {
	template <typename Tuple>
	__device__
	float operator()(Tuple t) const {
		/*
			Tuple consists of:
			center cell height
			north cell height
			south cell height
			west cell height
			east cell height
			north cell domain flag
			south cell domain flag
			west cell domain flag
			east cell domain flag
		*/
		
		float center = thrust::get<0>(t);
		float north = thrust::get<1>(t);
		float south = thrust::get<2>(t);
		float west = thrust::get<3>(t);
		float east = thrust::get<4>(t);

		// Mirror boundary condition - if NSWE cells not flagged as part of domain, treat height as equal to central cell
		if (!thrust::get<5>(t)) north = center;
		if (!thrust::get<6>(t)) south = center;
		if (!thrust::get<7>(t)) west = center;
		if (!thrust::get<8>(t)) east = center;
		
		thrust::get<9>(t) = (north + south + west + east) / 4.0f - center;
	}
};

struct update_velocity_with_dampening {
	float fDamp;
	update_velocity_with_dampening(float damp) : fDamp(damp) {}

	__device__
	float operator()(const float& fVel, const float& fAcc) const {
		return (fVel + fAcc) * fDamp;
	}
};

struct update_height_with_velocity {
	float fElapsedTime;
	update_height_with_velocity(float fTime) : fElapsedTime(fTime) {}

	__device__
	float operator()(const float& fHeight, const float& fVel) const {
		return fHeight + (fVel * fElapsedTime);
	}
};

void HeightField::step(const float& fElapsedTime)
{
	// Copy data onto device
	thrust::copy(h_z.begin(), h_z.end(), d_z.begin());
	thrust::copy(h_dz.begin(), h_dz.end(), d_dz.begin());
	thrust::copy(h_bDomain.begin(), h_bDomain.end(), d_bDomain.begin());
	
	auto start = thrust::make_zip_iterator(thrust::make_tuple(
		&d_z[nCols + 1],                         // Center
		&d_z[1],                                 // North
		&d_z[nCols + 1 + nCols],                 // South
		&d_z[nCols],                             // West
		&d_z[nCols + 2],                         // East
		&d_bDomain[1],                                 // North
		&d_bDomain[nCols + 1 + nCols],                 // South
		&d_bDomain[nCols],                             // West
		&d_bDomain[nCols + 2],						   // East
		d_ddz.begin() + nCols + 1));                   // Output
	auto finish = thrust::make_zip_iterator(thrust::make_tuple(
		&d_z[(nRows - 2) * nCols + (nCols - 1)], // Center - end(), so one-past-last element
		&d_z[(nRows - 3) * nCols + (nCols - 1)], // North
		&d_z[(nRows - 1) * nCols + (nCols - 1)], // South
		&d_z[(nRows - 2) * nCols + (nCols - 2)], // West
		&d_z[(nRows - 2) * nCols + (nCols)],     // East
		&d_bDomain[(nRows - 3) * nCols + (nCols - 1)], // North
		&d_bDomain[(nRows - 1) * nCols + (nCols - 1)], // South
		&d_bDomain[(nRows - 2) * nCols + (nCols - 2)], // West
		&d_bDomain[(nRows - 2) * nCols + (nCols)],     // East
		d_ddz.begin() + (nRows - 2) * nCols + (nCols - 1)));  // Output 

	// Apply transformations
	thrust::for_each(start, finish, get_velocity_change()); // This is our bottleneck
	thrust::transform(d_dz.begin(), d_dz.end(), d_ddz.begin(), d_dz.begin(), update_velocity_with_dampening(fDamp));
	thrust::transform(d_z.begin(), d_z.end(), d_dz.begin(), d_z.begin(), update_height_with_velocity(fElapsedTime));

	// Copy data back to host
	thrust::copy(d_z.begin(), d_z.end(), h_z.begin());
	thrust::copy(d_dz.begin(), d_dz.end(), h_dz.begin());
}

void HeightField::setHeights(float* heights)
{
	thrust::copy(heights, heights + (nRows * nCols), h_z.begin());
}

void HeightField::setHeight(const int& x, const int& y, const float& fHeight)
{
	h_z[y * nCols + x] = fHeight;
}

void HeightField::setDomain(bool* domain)
{
	thrust::copy(domain, domain + (nRows * nCols), h_bDomain.begin());
	// Ensure border is not part of domain
	for (int y = 0; y < nRows; y++)
	{
		for (int x = 0; x < nCols; x++)
		{
			if (y == 0 || x == 0 || y == nRows - 1 || x == nCols - 1) h_bDomain[y * nCols + x] = false;
		}
	}
}

void HeightField::setDomainCell(const int& x, const int& y, const bool& b)
{
	// Can't edit border cells
	if (x > 0 && x < nCols - 1 && y > 0 && y < nRows - 1) h_bDomain[y * nCols + x] = b;
}

void HeightField::zeroVelocities()
{
	for (int i = 0; i < (nRows * nCols); i++) h_dz[i] = 0;
}

void HeightField::clearDomain()
{
	for (int y = 0; y < nRows; y++)
	{
		for (int x = 0; x < nCols; x++)
		{
			if (y == 0 || x == 0 || y == nRows - 1 || x == nCols - 1) h_bDomain[y * nCols + x] = false;
			else h_bDomain[y * nCols + x] = true;
		}
	}
}

float HeightField::getHeight(const int& x, const int& y)
{
	return h_z[y * nCols + x];
}

bool HeightField::isInDomain(const int& x, const int& y)
{
	if (x < 0 || x >= nCols || y < 0 || y >= nRows) return false;
	return h_bDomain[y * nCols + x];
}