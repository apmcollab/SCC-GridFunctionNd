/*
 * SCC_GridFunction3d.h
 *
 *  Created on: Jun 27, 2015
 *      Author: anderson
 *
 * Decisions: Extending DoubleVector3D so that move semantics can be incorporated into
 * the underlying vector operations on grid values without having to duplicate all
 * member functions utilizing move semantics.
 *
 * Providing a Values member function so that data values can be accessed as if this
 * grid function implementation contained an instance of DoubleVector1D Values.
 *
 * Revised: Nov. 26, 2015 
 *          Jan. 26, 2016
 *
 *  Release : 18.09.04
 */

/*
#############################################################################
#
# Copyright  2015 Chris Anderson
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


#ifndef SCC_GRID_FUNCTION_3D_
#define SCC_GRID_FUNCTION_3D_


#include "SCC_GridFunction1d.h"
#include "SCC_GridFunction2d.h"

#include <cmath>
#include <functional>
#include <iostream>

#ifdef _MSC_VER
#include "iso646.h"          // So "and" is equivalenced to &&
typedef unsigned int uint;   // Define uint to be unsigned int
#undef min
#undef max
#endif

#include "../DoubleVectorNd/SCC_DoubleVector3d.h"


namespace SCC
{
class GridFunction3d : public DoubleVector3d
{

public :

GridFunction3d() : DoubleVector3d()
{
    this->xPanels = 0;
    this->yPanels = 0;
    this->zPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;
    this->hz      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;
    this->zMin = 0.0;
    this->zMax = 1.0;


    this->XYperiodicityFlag  = false;
    this->XYZperiodicityFlag = false;
}

GridFunction3d(const GridFunction3d& G) : DoubleVector3d(G)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;
    this->zPanels = G.zPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;
    this->hz     = G.hz;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
    this->zMin = G.zMin;
    this->zMax = G.zMax;

    this->XYperiodicityFlag  = G.XYperiodicityFlag;
    this->XYZperiodicityFlag = G.XYZperiodicityFlag;
}


GridFunction3d(DoubleVector3d&& G) : DoubleVector3d((DoubleVector3d&&)G)
{
    this->xPanels = 0;
    this->yPanels = 0;
    this->zPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;
    this->hz      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;
    this->zMin = 0.0;
    this->zMax = 1.0;


    this->XYperiodicityFlag  = false;
    this->XYZperiodicityFlag = false;
}


GridFunction3d(long xPanels, double hx, long yPanels, double hy, long zPanels, double hz)
: DoubleVector3d(xPanels+1,yPanels+1,zPanels+1)
{
    this->xPanels = xPanels;
    this->yPanels = yPanels;
    this->zPanels = zPanels;

    this->hx     = hx;
    this->hy     = hy;
    this->hz     = hz;

    this->xMin = -(xPanels*hx)/2.0;
    this->xMax =  (xPanels*hx)/2.0;
    this->yMin = -(yPanels*hy)/2.0;
    this->yMax =  (yPanels*hy)/2.0;
    this->zMin = -(zPanels*hz)/2.0;
    this->zMax =  (zPanels*hz)/2.0;

	this->setToValue(0.0);

    this->XYperiodicityFlag  = false;
    this->XYZperiodicityFlag = false;
}


GridFunction3d(long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax,
               long zPanels, double zMin, double zMax) : DoubleVector3d(xPanels+1,yPanels+1,zPanels+1)
{
    this->xPanels = xPanels;
    this->yPanels = yPanels;
    this->zPanels = zPanels;

    this->xMin = xMin;
    this->xMax = xMax;
    this->yMin = yMin;
    this->yMax = yMax;
    this->zMin = zMin;
    this->zMax = zMax;

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);
    this->hz     = (zMax-zMin)/(double)(zPanels);

	this->setToValue(0.0);

    this->XYperiodicityFlag  = false;
    this->XYZperiodicityFlag = false;
}


GridFunction3d(GridFunction3d&& G) : DoubleVector3d((DoubleVector3d&&)G)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;
    this->zPanels = G.zPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;
    this->hz     = G.hz;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
    this->zMin = G.zMin;
    this->zMax = G.zMax;

    this->XYperiodicityFlag  = G.XYperiodicityFlag;
    this->XYZperiodicityFlag = G.XYZperiodicityFlag;
}

virtual ~GridFunction3d(){}

GridFunction3d& operator=(const GridFunction3d& G)
{
//
// Only propagate domain information across assignment if the
// left hand side is a null instance.
//
if(this->dataPtr == nullptr)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;
    this->zPanels = G.zPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;
    this->hz     = G.hz;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
    this->zMin = G.zMin;
    this->zMax = G.zMax;


    this->XYperiodicityFlag  = G.XYperiodicityFlag;
    this->XYZperiodicityFlag = G.XYZperiodicityFlag;
}

// Propagate values

DoubleVector3d::operator=(G);

return *this;
}

GridFunction3d& operator=(GridFunction3d&& G)
{
//
// Only propagate domain information across assignment if the
// left hand side is a null instance.
//
if(this->dataPtr == nullptr)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;
    this->zPanels = G.zPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;
    this->hz     = G.hz;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
    this->zMin = G.zMin;
    this->zMax = G.zMax;
}

// Propagate values

DoubleVector3d::operator=((DoubleVector3d&&)G);
return *this;
}

GridFunction3d& operator=(DoubleVector3d&& G)
{
//
// Since underlying operations that result in DoubleVector3d&&G strip
// the domain information an assignment to a null instance generates
// an error.

if(this->dataPtr == nullptr)
{
	assert(domainCheck());
}

DoubleVector3d::operator=((DoubleVector3d&&)G);
return *this;
}


GridFunction3d& operator=(DoubleVector3d& G)
{
//
// Since underlying operations that result in DoubleVector3d&&G strip
// the domain information an assignment to a null instance generates
// an error.

if(this->dataPtr == nullptr)
{
	assert(domainCheck());
}

DoubleVector3d::operator=((DoubleVector3d&)G);
return *this;
}

//####  Incremental operators  ###

void operator*=(double alpha)
{SCC::DoubleVector3d::operator*=(alpha);}

void operator/=(double alpha)
{SCC::DoubleVector3d::operator/=(alpha);}


void operator+=(const GridFunction3d& G)
{SCC::DoubleVector3d::operator+=(G);}

void operator-=(const GridFunction3d& G)
{SCC::DoubleVector3d::operator-=(G);}

void operator*=(const GridFunction3d& G)
{SCC::DoubleVector3d::operator*=(G);}

void operator/=(const GridFunction3d& G)
{SCC::DoubleVector3d::operator/=(G);}


void operator+=(const DoubleVector3d& G)
{SCC::DoubleVector3d::operator+=(G);}

void operator-=(const DoubleVector3d& G)
{SCC::DoubleVector3d::operator-=(G);}

void operator*=(const DoubleVector3d& G)
{SCC::DoubleVector3d::operator*=(G);}

void operator/=(const DoubleVector3d& G)
{SCC::DoubleVector3d::operator/=(G);}


void operator*=(const  std::function<double(double,double,double)>& F)
{
    double xPos; double yPos; double zPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    Values(i,j,k) *= F(xPos,yPos,zPos);
    }}}
}

void operator/=(const std::function<double(double,double,double)>& F)
{
    double xPos; double yPos; double zPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    Values(i,j,k) /= F(xPos,yPos,zPos);
    }}}
}


void operator+=(const  std::function<double(double,double,double)>& F)
{
    double xPos; double yPos; double zPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    Values(i,j,k) += F(xPos,yPos,zPos);
    }}}
}

void operator-=(const  std::function<double(double,double,double)>& F)
{
    double xPos; double yPos; double zPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    Values(i,j,k) -= F(xPos,yPos,zPos);
    }}}
}

// ######################################################################

// Initialization

void initialize()
{
    DoubleVector3d::initialize();
    this->xPanels = 0;
    this->yPanels = 0;
    this->zPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;
    this->hz      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;
    this->zMin = 0.0;
    this->zMax = 1.0;

    this->XYperiodicityFlag  = false;
    this->XYZperiodicityFlag = false;
}

void initialize(const GridFunction3d& G)
{
	DoubleVector3d::initialize(G);
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;
    this->zPanels = G.zPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;
    this->hz     = G.hz;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;
    this->zMin = G.zMin;
    this->zMax = G.zMax;

    this->XYperiodicityFlag  = G.XYperiodicityFlag;
    this->XYZperiodicityFlag = G.XYZperiodicityFlag;
}

void initialize(long xPanels, double hx, long yPanels, double hy, long zPanels, double hz)
{
    DoubleVector3d::initialize(xPanels+1,yPanels+1,zPanels+1);
    this->xPanels = xPanels;
    this->yPanels = yPanels;
    this->zPanels = zPanels;

    this->hx     = hx;
    this->hy     = hy;
    this->hz     = hz;

    this->xMin = -(xPanels*hx)/2.0;
    this->xMax =  (xPanels*hx)/2.0;
    this->yMin = -(yPanels*hy)/2.0;
    this->yMax =  (yPanels*hy)/2.0;
    this->zMin = -(zPanels*hz)/2.0;
    this->zMax =  (zPanels*hz)/2.0;

	this->setToValue(0.0);

	this->XYperiodicityFlag  = false;
    this->XYZperiodicityFlag = false;
}


void initialize(long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax,
               long zPanels, double zMin, double zMax)
{
    DoubleVector3d::initialize(xPanels+1,yPanels+1,zPanels+1);
    this->xPanels = xPanels;
    this->yPanels = yPanels;
    this->zPanels = zPanels;

    this->xMin = xMin;
    this->xMax = xMax;
    this->yMin = yMin;
    this->yMax = yMax;
    this->zMin = zMin;
    this->zMax = zMax;

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);
    this->hz     = (zMax-zMin)/(double)(zPanels);

	this->setToValue(0.0);

	this->XYperiodicityFlag  = false;
    this->XYZperiodicityFlag = false;
}

GridFunction3d* newDuplicate() const
{
    GridFunction3d* Mptr = new GridFunction3d(*this);
    return Mptr;
}

bool isNull() const
{
if(dataPtr == nullptr) return true;
return false;
}

void setPeriodicity(bool val = true)
{
	if(val)
	{
		XYZperiodicityFlag = true;
	    enforcePeriodicity();
	}
}

// Clears all types of periodicity

void clearPeriodicity()
{
	XYZperiodicityFlag = false;
	XYperiodicityFlag = false;
}

void setXYperiodicity(bool val = true)
{
	if(val)
	{
		XYperiodicityFlag = true;
	    enforceXYperiodicity();
	}
}

void clearXYZperiodicity()
{
	XYZperiodicityFlag = false;
}

void setXYZperiodicity(bool val = true)
{
	if(val)
	{
		XYZperiodicityFlag = true;
	    enforcePeriodicity();
	}
}

void clearXYperiodicity()
{
	XYperiodicityFlag = false;
}

/*!  Returns the size of the number of independent values
     associated with the grid function, i.e. the dimension
     corresponding to the vector of independent function
     values.

     This dimension depends upon the specified boundary values.
*/

virtual long getDimension() const
{
	long dimension;

    if(not (XYperiodicityFlag || XYZperiodicityFlag))
    {
    	dimension = DoubleVector3d::getDimension();
    }
    else
    {
    	if(XYperiodicityFlag)  {dimension = (index1Size-1)*(index2Size-1)*(index3Size);}
    	if(XYZperiodicityFlag) {dimension = (index1Size-1)*(index2Size-1)*(index3Size-1);}
    }

    return dimension;
}

void createProductFunction(const GridFunction1d& funX, const GridFunction1d& funY, const GridFunction1d& funZ)
{
    initialize(funX.xPanels, funX.xMin, funX.xMax,
               funY.xPanels, funY.xMin, funY.xMax,
               funZ.xPanels, funZ.xMin, funZ.xMax);

	long i; long j; long k;

	double fX; double fY; double fZ;

	for(i = 0; i <= xPanels; i++)
	{
	fX = funX.Values(i);
	for(j = 0; j <= yPanels; j++)
	{
	fY = funY.Values(j);
	for(k = 0; k <= zPanels; k++)
	{
	fZ = funZ.Values(k);

	Values(i,j,k) = fX*fY*fZ;
	}}}

	if(XYperiodicityFlag){enforceXYperiodicity();}
	if(XYZperiodicityFlag){enforcePeriodicity();}

}

void createProductFunction(const GridFunction2d& funXY,  const GridFunction1d& funZ)
{
    initialize(funXY.xPanels, funXY.xMin, funXY.xMax,
               funXY.yPanels, funXY.yMin, funXY.yMax,
               funZ.xPanels,   funZ.xMin,  funZ.xMax);


	long i; long j; long k;

	double fXY; double fZ;

	for(i = 0; i <= xPanels; i++)
	{
    for(j = 0; j <= yPanels; j++)
	{
	fXY = funXY(i,j);
	for(k = 0; k <= zPanels; k++)
	{
	fZ = funZ(k);

	Values(i,j,k) = fXY*fZ;
	}}}

	if(XYperiodicityFlag){enforceXYperiodicity();}
	if(XYZperiodicityFlag){enforcePeriodicity();}
}

void createProductFunction(const GridFunction1d& funX,const  GridFunction2d& funYZ)
{

    initialize(funX.xPanels,   funX.xMin,  funX.xMax,
               funYZ.xPanels, funYZ.xMin, funYZ.xMax,
               funYZ.yPanels, funYZ.yMin, funYZ.yMax);


	long i; long j; long k;
	double fX; double fYZ;

	for(i = 0; i <= xPanels; i++)
	{
	fX = funX(i);
	for(j = 0; j <= yPanels; j++)
	{
	for(k = 0; k <= zPanels; k++)
	{
	fYZ = funYZ(j,k);
	Values(i,j,k) = fX*fYZ;
	}}}

    if(XYperiodicityFlag){enforceXYperiodicity();}
	if(XYZperiodicityFlag){enforcePeriodicity();}
}


void specify(std::function<double(double,double,double)> F)
{
    double xPos; double yPos; double zPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    for(long k = 0; k <= zPanels; k++)
    {
    zPos = zMin + k*hz;
    Values(i,j,k) = F(xPos,yPos,zPos);
    }}}

    if(XYperiodicityFlag){enforceXYperiodicity();}
	if(XYZperiodicityFlag){enforcePeriodicity();}
}

void squareValues()
{
    transformValues([](double x){return x*x;});
}

void zeroNegativePart()
{
    transformValues([](double x){if(x < 0.0){return 0.0;} return x;});
}


//
// Important: the dot product for this class is a mesh
// scaled inner product based upon the Trapezoidal rule.
//
// If the function is periodic, or satisfies homogeneous boundary
// conditions, then the inner product computed is identical
// to the standard inner product of all independent function
// values times a product of the mesh sizes.
//
using DoubleVector3d::dot;

double dot(const GridFunction3d& V) const
{
	return scaledDot(V);
}
//
// Dot product of all values scaled with mesh widths
//           (Trapezoidal Method)
//

double scaledDot(const GridFunction3d& V) const
{
    double dotVal = 0.0;

    long i; long j; long k;
    double hxhyhz = hx*hy*hz;

//  Interior points weighted by 1

    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz;
    }}}
//
//  Interior face points (weighted by 1/2)
//
    i = 0;
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}

    i = xPanels;
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}

    //

    j = 0;
    for(i = 1; i < xPanels; i++)
    {
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}

    j = yPanels;
    for(i = 1; i <  xPanels; i++)
    {
    for(k = 1; k <  zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}

    //
    k = 0;
    for(i = 1; i <  xPanels; i++)
    {
    for(j = 1; j <  yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}

    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.5;
    }}
//
//  Interior corner points (weighted by 1/4)
//
    i = 0;
    j = 0;
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = 0;
    j = yPanels;
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = xPanels;
    j = 0;
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = xPanels;
    j = yPanels;
    for(k = 1; k < zPanels; k++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }


    i = 0;
    k = 0;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = xPanels;
    k = 0;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = 0;
    k = zPanels;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    i = xPanels;
    k = zPanels;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }
//
    j = 0;
    k = 0;
    for(i = 1; i < xPanels; i++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    j = yPanels;
    k = 0;
    for(i = 1; i < xPanels; i++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    j = 0;
    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }

    j = yPanels;
    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    dotVal += Values(i,j,k)*V.Values(i,j,k)*hxhyhz*0.25;
    }
//
//  Corner points (weighted by 1/8)
//
    dotVal += Values(0,0,0)*V.Values(0,0,0)*hxhyhz*.125;
    dotVal += Values(xPanels,0,0)*V.Values(xPanels,0,0)*hxhyhz*.125;

    dotVal += Values(0,yPanels,0)*V.Values(0,yPanels,0)*hxhyhz*.125;
    dotVal += Values(xPanels,yPanels,0)*V.Values(xPanels,yPanels,0)*hxhyhz*.125;

    dotVal += Values(0,0,zPanels)*V.Values(0,0,zPanels)*hxhyhz*.125;
    dotVal += Values(xPanels,0,zPanels)*V.Values(xPanels,0,zPanels)*hxhyhz*.125;

    dotVal += Values(0,yPanels,zPanels)*V.Values(0,yPanels,zPanels)*hxhyhz*.125;
    dotVal += Values(xPanels,yPanels,zPanels)*V.Values(xPanels,yPanels,zPanels)*hxhyhz*.125;

    return dotVal;
}



//  Trapezoidal method

double integralTrapezoidal(std::function<double(double)> F) const
{
    double intVal = 0.0;
    long i; long j; long k;

    double hxhyhz = hx*hy*hz;
//
//  Interior points weighted by 1
//
    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz;
    }}}
//
//  Interior face points (weighted by 1/2)
//
    i = 0;
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}

    i = xPanels;
    for(j = 1; j < yPanels; j++)
    {
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}

    //

    j = 0;
    for(i = 1; i < xPanels; i++)
    {
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}

    j = yPanels;
    for(i = 1; i < xPanels; i++)
    {
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}

    //
    k = 0;
    for(i = 1; i <  xPanels; i++)
    {
    for(j = 1; j <  yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}

    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.5;
    }}
//
//  Interior corner points (weighted by 1/4)
//
    i = 0;
    j = 0;
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = 0;
    j = yPanels;
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = xPanels;
    j = 0;
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = xPanels;
    j = yPanels;
    for(k = 1; k < zPanels; k++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    //

    i = 0;
    k = 0;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = xPanels;
    k = 0;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = 0;
    k = zPanels;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    i = xPanels;
    k = zPanels;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }
//
    j = 0;
    k = 0;
    for(i = 1; i < xPanels; i++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    j = yPanels;
    k = 0;
    for(i = 1; i < xPanels; i++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    j = 0;
    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }

    j = yPanels;
    k = zPanels;
    for(i = 1; i < xPanels; i++)
    {
    intVal += F(Values(i,j,k))*hxhyhz*0.25;
    }
//
//  Corner points (weighted by 1/8)
//
    intVal += F(Values(0,0,0))*hxhyhz*.125;
    intVal += F(Values(xPanels,0,0))*hxhyhz*.125;

    intVal += F(Values(0,yPanels,0))*hxhyhz*.125;
    intVal += F(Values(xPanels,yPanels,0))*hxhyhz*.125;

    intVal += F(Values(0,0,zPanels))*hxhyhz*.125;
    intVal += F(Values(xPanels,0,zPanels))*hxhyhz*.125;

    intVal += F(Values(0,yPanels,zPanels))*hxhyhz*.125;
    intVal += F(Values(xPanels,yPanels,zPanels))*hxhyhz*.125;

    return intVal;
}

double norm1() const
{
  return integralTrapezoidal([](double x){return std::abs(x);});
}

// Trapezoidal Method Integral Approximation

double integralTrapezoidal() const
{
    return integralTrapezoidal([](double x){return x;});
}

// 2-Norm^2 computed using the trapezoidal rule approximation

double norm2() const
{
    return std::sqrt(std::abs(integralTrapezoidal([](double x){return x*x;})));
}

double nrm2() const
{
    return std::sqrt(std::abs(integralTrapezoidal([](double x){return x*x;})));
}


// 2-Norm^2 computed using the trapezoidal rule approximation

double norm2squared() const
{
    return integralTrapezoidal([](double x){return x*x;});
}

//
// Riemann Sum
//

double  integralRiemann() const
{
//  All points weighted by hxhy

    double intVal = 0.0;
    double hxhy = hx*hy;

    long dCount  = (xPanels+1)*(yPanels+1)*(zPanels+1);

    for(long k =0; k < dCount; k++)
    {
    intVal += dataPtr[k];
    }

    intVal *= hxhy;

    return intVal;
}

//  Average computed using Riemann sum approximation of the integral

double getRiemannAverage() const
{
    double avgValue;
    avgValue = integralRiemann();
    return avgValue/((xMax-xMin)*(yMax-yMin)*(zMax-zMin));
}

//  Average computed using Trapezoidal approximation of the integral

double getTrapezoidalAverage() const
{
    double avgValue;
    avgValue = integralTrapezoidal();
    return avgValue/((xMax-xMin)*(yMax-yMin)*(zMax-zMin));
}


double min() const
{
//  Compute function minimum

    long dCount      = (xPanels+1)*(yPanels+1)*(zPanels+1);
    double minValue  = dataPtr[0];
    for(long k = 1; k < dCount; k++)
    {
    minValue = ( minValue < dataPtr[k] ) ? minValue : dataPtr[k];
    }
    return minValue;

}

double max() const
{
//  Compute function maximum

    long dCount  = (xPanels+1)*(yPanels+1)*(zPanels+1);

    double maxValue  = dataPtr[0];
    for(long k = 1; k < dCount; k++)
    {
    maxValue = ( maxValue > dataPtr[k] ) ? maxValue : dataPtr[k];
    }
    return maxValue;
}

GridFunction2d getConstantZslice(long zIndex) const //(x-y function)
{
    GridFunction2d R(xPanels,xMin,xMax,yPanels,yMin,yMax);
    for(long i = 0; i <= xPanels; i++)
    {
    for(long j = 0; j <= yPanels; j++)
    {
    R.Values(i,j) = Values(i,j,zIndex);
    }}
    return R;
}

void getConstantZslice(long zIndex,GridFunction2d& R) const //(x-y function)
{
    for(long i = 0; i <= xPanels; i++)
    {
    for(long j = 0; j <= yPanels; j++)
    {
    R.Values(i,j) = Values(i,j,zIndex);
    }}
}

void setConstantZslice(long zIndex,const GridFunction2d& R) //(x-y function)
{
    for(long i = 0; i <= xPanels; i++)
    {
    for(long j = 0; j <= yPanels; j++)
    {
    Values(i,j,zIndex) = R.Values(i,j);
    }}
}

GridFunction2d getConstantYslice(long yIndex) const //(x-z function)
{
	GridFunction2d R(xPanels,xMin,xMax,zPanels,zMin,zMax);
    for(long i = 0; i <= xPanels; i++)
    {
    for(long k = 0; k <= zPanels; k++)
    {
    R.Values(i,k) = Values(i,yIndex,k);
    }}
    return R;
}

void getConstantYslice(long yIndex,GridFunction2d& R) const //(x-z function)
{
    for(long i = 0; i <= xPanels; i++)
    {
    for(long k = 0; k <= zPanels; k++)
    {
    R.Values(i,k) = Values(i,yIndex,k);
    }}
}

void setConstantYslice(long yIndex,const GridFunction2d& R)  //(x-z function)
{
    for(long i = 0; i <= xPanels; i++)
    {
    for(long k = 0; k <= zPanels; k++)
    {
    Values(i,yIndex,k) = R.Values(i,k);
    }}
}


GridFunction2d getConstantXslice(long xIndex) const //(y-z function)
{
	GridFunction2d R(yPanels,yMin,yMax,zPanels,zMin,zMax);
    for(long j = 0; j <= yPanels; j++)
    {
    for(long k = 0; k <= zPanels; k++)
    {
    R.Values(j,k) = Values(xIndex,j,k);
    }}
    return R;
}

void getConstantXslice(long xIndex, GridFunction2d& R) const //(y-z function)
{
    for(long j = 0; j <= yPanels; j++)
    {
    for(long k = 0; k <= zPanels; k++)
    {
    R.Values(j,k) = Values(xIndex,j,k);
    }}
}

void setConstantXslice(long xIndex, const GridFunction2d& R) //(y-z function)
{
    for(long j = 0; j <= yPanels; j++)
    {
    for(long k = 0; k <= zPanels; k++)
    {
    Values(xIndex,j,k) = R.Values(j,k);
    }}
}


GridFunction1d getConstantYZslice(long yIndex, long zIndex) const  // (x function)
{
	GridFunction1d R(xPanels,xMin,xMax);
	for(long i = 0; i <= xPanels; i++)
	{
	R.Values(i) = Values(i,yIndex,zIndex);
	}
	return R;
}

void getConstantYZslice(long yIndex, long zIndex, GridFunction1d& R) const  // (x function)
{
	for(long i = 0; i <= xPanels; i++)
	{
	R.Values(i) = Values(i,yIndex,zIndex);
	}
}

void setConstantYZslice(long yIndex, long zIndex,const GridFunction1d& R)   // (x function)
{
	for(long i = 0; i <= xPanels; i++)
	{
	Values(i,yIndex,zIndex) = R.Values(i);
	}
}


GridFunction1d getConstantXZslice(long xIndex, long zIndex) const  //( y function)
{
	GridFunction1d R(yPanels,yMin,yMax);
	for(long j = 0; j <= yPanels; j++)
	{
	R.Values(j) = Values(xIndex,j,zIndex);
	}
	return R;
}

void  getConstantXZslice(long xIndex, long zIndex,GridFunction1d& R) const  //( y function)
{
	for(long j = 0; j <= yPanels; j++)
	{
	R.Values(j) = Values(xIndex,j,zIndex);
	}
}

void  setConstantXZslice(long xIndex, long zIndex,const GridFunction1d& R)    //( y function)
{
	for(long j = 0; j <= yPanels; j++)
	{
	Values(xIndex,j,zIndex) = R.Values(j);
	}
}

GridFunction1d getConstantXYslice(long xIndex, long yIndex) const  //( z function)
{
	GridFunction1d R(zPanels,zMin,zMax);
	for(long k = 0; k <= zPanels; k++)
	{
	R.Values(k) = Values(xIndex,yIndex,k);
	}
	return R;
}

void getConstantXYslice(long xIndex, long yIndex,GridFunction1d& R) const  //( z function)
{
	for(long k = 0; k <= zPanels; k++)
	{
	R.Values(k) = Values(xIndex,yIndex,k);
	}
}

void setConstantXYslice(long xIndex, long yIndex,const GridFunction1d& R)    //( z function)
{
	for(long k = 0; k <= zPanels; k++)
	{
	Values(xIndex,yIndex,k) = R.Values(k);
	}
}


virtual bool isXperiodic()   const {return (XYperiodicityFlag  || XYZperiodicityFlag);}
virtual bool isYperiodic()   const {return (XYperiodicityFlag  || XYZperiodicityFlag);}
virtual bool isZperiodic()   const {return (XYZperiodicityFlag); }

virtual bool isXYperiodic()  const {return (XYperiodicityFlag  || XYZperiodicityFlag);}
virtual bool isXYZperiodic() const {return (XYZperiodicityFlag);}


double getHx()   const  {return            hx;}
double getXmin() const  {return          xMin;}
double getXmax() const  {return          xMax;}
long   getXpanelCount() const {return xPanels;}

double getHy()   const  {return            hy;}
double getYmin() const  {return          yMin;}
double getYmax() const  {return          yMax;}
long   getYpanelCount() const {return yPanels;}

double getHz()   const  {return            hz;}
double getZmin() const  {return          zMin;}
double getZmax() const  {return          zMax;}
long   getZpanelCount() const {return zPanels;}


DoubleVector3d getValues() const
{
    return DoubleVector3d(*this);
}

DoubleVector3d* getValuesPointer()
{
    return (DoubleVector3d*)this;
}

const DoubleVector3d* getValuesPointer() const
{
    return (DoubleVector3d*)this;
}

/// Set boundary values to specified value

void setBoundaryValues(double value)
{
    long i; long j; long k;

	i = 0;
	for(j = 0; j <= yPanels; j++)
	{
	for(k = 0; k <= zPanels; k++)
	{
		Values(i,j,k) = value;
	}}

    i = xPanels;
    for(j = 0; j <= yPanels; j++)
    {
    for(k = 0; k <= zPanels; k++)
    {
     Values(i,j,k) = value;
    }}


    j = 0;
    for(i = 0; i <= xPanels; i++)
    {
    for(k = 0; k <= zPanels; k++)
    {
    Values(i,j,k) = value;
    }}

    j = yPanels;
    for(i = 0; i <= xPanels; i++)
    {
    for(k = 0; k <= zPanels; k++)
    {
    Values(i,j,k) = value;
    }}

//

    k = 0;
    for(i = 0; i <=  xPanels; i++)
    {
    for(j = 0; j <=  yPanels; j++)
    {
    Values(i,j,k) = value;
    }}

    k = zPanels;
    for(i = 0; i <= xPanels; i++)
    {
    for(j = 0; j <= yPanels; j++)
    {
    Values(i,j,k) = value;
    }}
}

//
// Enforces periodicity in all three x,y,z coordinate directions
//
void enforcePeriodicity()
{
	XYZperiodicityFlag = true;
    if(this->isNull()){return;}

    long i; long j; long k;

    i = xPanels;
    for(j = 0; j <= yPanels; j++)
    {
    for(k = 0; k <= zPanels; k++)
    {
     Values(i,j,k) = Values(0,j,k);
    }}

    j = yPanels;
    for(i = 0; i <= xPanels; i++)
    {
    for(k = 0; k <= zPanels; k++)
    {
    Values(i,j,k) = Values(i,0,k);
    }}

    k = zPanels;
    for(i = 0; i <= xPanels; i++)
    {
    for(j = 0; j <= yPanels; j++)
    {
    Values(i,j,k) = Values(i,j,0);
    }}
}

// Enforces periodicity in the x-y coordinate directions

void enforceXYperiodicity()
{
	XYperiodicityFlag = true;
	if(this->isNull()){return;}

    long i; long j; long k;

    i = xPanels;
    for(j = 0; j <= yPanels; j++)
    {
    for(k = 0; k <= zPanels; k++)
    {
     Values(i,j,k) = Values(0,j,k);
    }}

    j = yPanels;
    for(i = 0; i <= xPanels; i++)
    {
    for(k = 0; k <= zPanels; k++)
    {
    Values(i,j,k) = Values(i,0,k);
    }}
}

//  Returns true if the input grid function is structurally identical,
//  e.g. a coincident computational domain and mesh size
//

bool isCoincident(const GridFunction3d& V)
{
	double domainDiffTol    = 1.0e-13;

    long discretizationDiff = std::abs(xPanels - V.xPanels)
                            + std::abs(yPanels - V.yPanels)
                            + std::abs(zPanels - V.zPanels);

    double domainDiff = (std::abs(xMin - V.xMin) + std::abs(xMax - V.xMax))/std::abs(xMax-xMin)
                      + (std::abs(yMin - V.yMin) + std::abs(yMax - V.yMax))/std::abs(yMax-yMin)
                      + (std::abs(zMin - V.zMin) + std::abs(zMax - V.zMax))/std::abs(zMax-zMin);

    if((discretizationDiff > 0) || (domainDiff > domainDiffTol)) return false;
    return true;
}

/*! BLAS axpby : this  <- alpha*v + beta*this  */

   void axpby(double alpha, const DoubleVector3d& G, double beta)
   {
   SCC::DoubleVector3d::axpby(alpha,G,beta);
   }

/*! BLAS axpy : this  <- alpha*v + this  */

   void axpy(double alpha, const DoubleVector3d& G)
   {
   SCC::DoubleVector3d::axpy(alpha,G);
   }

//  Grid Geometry

    double xMin; double xMax;   // Computational Region is [xMin,xMax]x
    double yMin; double yMax;   //                         [yMin,yMax]x
    double zMin; double zMax;   //                         [zMin,zMax]

    long   xPanels;             // Number of panels
    long   yPanels;
    long   zPanels;

    double hx;                  // mesh width
    double hy;
    double hz;

    bool  XYperiodicityFlag;     // indicates periodicity in the XY direction
    bool XYZperiodicityFlag;     // indicates periodicity in the XYZ direction

//###################################################################
//          Values Access (Alternative to operator())
//###################################################################

#ifdef _DEBUG
    double&  Values(long i1, long i2, long i3)
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    assert(boundsCheck(i3, 0, index3Size-1,3));
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };

    const double&  Values(long i1, long i2, long i3) const
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    assert(boundsCheck(i3, 0, index3Size-1,3));
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };
#else
    /*!
    Returns a reference to the element with index (i1,i2,i3) - indexing
    starting at (0,0,0).
    */
    inline double&  Values(long i1, long i2, long i3)
    {
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };

    /*!
    Returns a reference to the element with index (i1,i2,i3) - indexing
    starting at (0,0,0).
    */
    inline const double&  Values(long i1, long i2, long i3) const
    {
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };
#endif

//###################################################################
//                      Bounds Checking
//###################################################################

#ifdef _DEBUG
        bool domainCheck() const
        {
        std::cerr << "XXX SCC::GridFunction3d Error  XXX" << std::endl;
        std::cerr << "Left side of assignment must be a non-null GridFunction3d instance." << std::endl;
        std::cerr << std::endl;
        return false;
        }
#else
        bool domainCheck() const {return true;}
#endif



};

}


#endif /* SCC_GRIDFUNCTIONSCC_GRIDFUNCTION1D_H_ */
