/*
 * SCC_GridFunctionNdSkin.h
 *
 * Classes : SCC::GridFunction1dSkin
 *           SCC::GridFunction2dSkin
 *           SCC::GridFunction3dSkin
 *
 *
 * These classes provide a GridFunctionNd wrapper (or skin) for existing data.
 *
 * It is assumed that the data associated with an X-panels x Y-panels x Z-panels
 * grid is a single double* array of size (X-panels + 1)(Y-panels + 1)(Z-panels + 1)
 * stored by ROWS (C convention).
 *
 * There is no internal validation on the assumption about the size and storage
 * format associated with the data pointer used to initialize an instance
 * of the GridFunctionNdSkin class.
 *
 * The methods are primarily provided as a utility for utilizing
 * methods whose input and output arugments are SCC::GridFunctionNd
 * instances.
 *
 * Typical usage
 *
 * G = a grid object with data stored in a double* objectData
 *
 * (1) GridFunctionNdSkin GS(G.objectData, xPanels, xMin, xMax, ....)
 *
 * (2) Call member function requiring GridFunctionNd arguments and specify GS instead of G
 *
 * Transformations of the data associated with GS will will be performed on the
 * data of the associated grid object G.
 *
 * No data copying is required.
 *
 *
 *  Created on: Sep 30, 2016
 *      Author: anderson
 */

/*
#############################################################################
#
# Copyright 2015-16 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/

#include "SCC_GridFunction1d.h"
#include "SCC_GridFunction2d.h"
#include "SCC_GridFunction3d.h"


#ifndef _SCC_GridFunctionNdSkin_
#define _SCC_GridFunctionNdSkin_


namespace SCC
{

class GridFunction1dSkin : public GridFunction1d
{
public:

	 GridFunction1dSkin() : GridFunction1d()
	 {
		 this->dataPtr = nullptr;
	 };


	 GridFunction1dSkin(const GridFunction1dSkin& G)
	 {
		 initialize(G.dataPtr, G.xPanels, G.xMin, G.xMax);
	 };


	 GridFunction1dSkin(double* data1d, long xPanels, double xMin, double xMax)
	 {
		 initialize(data1d, xPanels, xMin, xMax);
	 }

	 void initialize(double* data1d, long xPanels, double xMin, double xMax)
	 {
	    // Set internal data pointer to null ptr then call base class initialize

		this->dataPtr = nullptr;
		GridFunction1d::initialize();

		// Set internal data pointer to supplied data pointer

		this->dataPtr = data1d;

		// Pack structure informaiton

		// DoubleVector1d structure information

		this->index1Size = xPanels+1;

    	// GridFunction1d structure information

		this->xPanels = xPanels;

		this->xMin = xMin;
		this->xMax = xMax;

		this->hx  = (xMax-xMin)/(double)(xPanels);

	 }


	 virtual ~GridFunction1dSkin()
	 {
		// Don't delete the data associated with the instance

		 this->dataPtr = nullptr;
	 }
};


class GridFunction2dSkin : public GridFunction2d
{
public:

	 GridFunction2dSkin() : GridFunction2d()
	 {
		 this->dataPtr = nullptr;
	 };


	 GridFunction2dSkin(const GridFunction2dSkin& G)
	 {
		 initialize(G.dataPtr, G.xPanels, G.xMin, G.xMax, G.yPanels, G.yMin, G.yMax);
	 };


	 GridFunction2dSkin(double* data2d, long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax)
	 {
		 initialize(data2d, xPanels, xMin, xMax, yPanels, yMin, yMax);
	 }


	 void initialize(double* data2d, long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax)
	 {
	    // Set internal data pointer to null ptr then call base class initialize

		this->dataPtr = nullptr;
		GridFunction2d::initialize();

		// Set internal data pointer to supplied data pointer

		this->dataPtr = data2d;

		// Pack structure informaiton

		// DoubleVector3d structure information

		this->index1Size = xPanels+1;
    	this->index2Size = yPanels+1;

    	// GridFunction3d structure information

		this->xPanels = xPanels;
		this->yPanels = yPanels;

		this->xMin = xMin;
		this->xMax = xMax;
		this->yMin = yMin;
		this->yMax = yMax;

		this->hx  = (xMax-xMin)/(double)(xPanels);
		this->hy  = (yMax-yMin)/(double)(yPanels);
	 }


	 virtual ~GridFunction2dSkin()
	 {
		// Don't delete the data associated with the instance

		 this->dataPtr = nullptr;
	 }
};



class GridFunction3dSkin : public GridFunction3d
{
public:

	 GridFunction3dSkin() : GridFunction3d()
	 {
		 this->dataPtr = nullptr;
	 };


	 GridFunction3dSkin(const GridFunction3dSkin& G)
	 {
		 initialize(G.dataPtr, G.xPanels, G.xMin, G.xMax, G.yPanels, G.yMin, G.yMax, G.zPanels, G.zMin, G.zMax);
	 };


	 GridFunction3dSkin(double* data3d, long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax,
                        long zPanels, double zMin, double zMax)
	 {
		 initialize(data3d, xPanels, xMin, xMax, yPanels, yMin, yMax, zPanels, zMin, zMax);
	 }


	 void initialize(double* data3d, long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax,
                     long zPanels, double zMin, double zMax)
	 {
	    // Set internal data pointer to null ptr then call base class initialize

		this->dataPtr = nullptr;
		GridFunction3d::initialize();

		// Set internal data pointer to supplied data pointer

		this->dataPtr = data3d;

		// Pack structure informaiton

		// DoubleVector3d structure information

		this->index1Size = xPanels+1;
    	this->index2Size = yPanels+1;
    	this->index3Size = zPanels+1;

    	// GridFunction3d structure information

		this->xPanels = xPanels;
		this->yPanels = yPanels;
		this->zPanels = zPanels;

		this->xMin = xMin;
		this->xMax = xMax;
		this->yMin = yMin;
		this->yMax = yMax;
		this->zMin = zMin;
		this->zMax = zMax;

		this->hx  = (xMax-xMin)/(double)(xPanels);
		this->hy  = (yMax-yMin)/(double)(yPanels);
		this->hz  = (zMax-zMin)/(double)(zPanels);
	 }


	 virtual ~GridFunction3dSkin()
	 {
		// Don't delete the data associated with the instance

		 this->dataPtr = nullptr;
	 }
};
}




#endif /*_SCC_GridFunctionNdSkin_*/
