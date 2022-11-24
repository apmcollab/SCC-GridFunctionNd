/*
 * SCC_GridFunction2d.h
 *
 *  Created on: Jun 27, 2015
 *      Author: anderson
 *
 * Release : 18.09.04
 *
 * Decisions: Extending DoubleVector2D so that move semantics can be incorporated into
 * the underlying vector operations on grid values without having to duplicate all
 * member functions utilizing move semantics.
 *
 * Providing a Values member function so that data values can be accessed as if this
 * grid function implementation contained an instance of DoubleVector1D Values.
 *
 * Revised: Nov. 26, 2015 
 *          Jan. 26, 2016
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

#ifndef SCC_GRID_FUNCTION_2D_
#define SCC_GRID_FUNCTION_2D_

#include <functional>
#include <iostream>
#include <cmath>

#ifdef _MSC_VER
#include "iso646.h"          // So "and" is equivalenced to &&
typedef unsigned int uint;   // Define uint to be unsigned int
#undef min
#undef max
#endif

#include "../GridFunctionNd/SCC_GridFunction1d.h"
#include "../DoubleVectorNd/SCC_DoubleVector2d.h"

namespace SCC
{
class GridFunction2d : public DoubleVector2d
{

public :

GridFunction2d() : DoubleVector2d()
{
    this->xPanels = 0;
    this->yPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;

    this->XYperiodicityFlag = false;
    this->XperiodicityFlag  = false;
    this->YperiodicityFlag  = false;
}

GridFunction2d(const GridFunction2d& G) : DoubleVector2d(G)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;

    this->XYperiodicityFlag = G.XYperiodicityFlag;
    this->XperiodicityFlag  = G.XperiodicityFlag;
    this->YperiodicityFlag  = G.YperiodicityFlag;
}


GridFunction2d(DoubleVector2d&& G) : DoubleVector2d((DoubleVector2d&&)G)
{
    this->xPanels = 0;
    this->yPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;

    this->XYperiodicityFlag = false;
    this->XperiodicityFlag  = false;
    this->YperiodicityFlag  = false;
}


GridFunction2d(long xPanels, double hx, long yPanels, double hy)
: DoubleVector2d(xPanels+1,yPanels+1)
{
    this->xPanels = xPanels;
    this->yPanels = yPanels;

    this->hx     = hx;
    this->hy     = hy;

    this->xMin = -(xPanels*hx)/2.0;
    this->xMax =  (xPanels*hx)/2.0;
    this->yMin = -(yPanels*hy)/2.0;
    this->yMax =  (yPanels*hy)/2.0;

    this->XYperiodicityFlag = false;
    this->XperiodicityFlag  = false;
    this->YperiodicityFlag  = false;
	this->setToValue(0.0);
}


GridFunction2d(long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax)
: DoubleVector2d(xPanels+1,yPanels+1)
{
    this->xPanels = xPanels;
    this->yPanels = yPanels;

    this->xMin = xMin;
    this->xMax = xMax;
    this->yMin = yMin;
    this->yMax = yMax;

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);

    this->XYperiodicityFlag = false;
    this->XperiodicityFlag  = false;
    this->YperiodicityFlag  = false;
	this->setToValue(0.0);
}


GridFunction2d(GridFunction2d&& G) : DoubleVector2d((DoubleVector2d&&)G)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;

    this->XYperiodicityFlag = G.XYperiodicityFlag;
    this->XperiodicityFlag  = G.XperiodicityFlag;
    this->YperiodicityFlag  = G.YperiodicityFlag;
}

virtual ~GridFunction2d(){}

// Initialization

void initialize()
{
    DoubleVector2d::initialize();
    this->xPanels = 0;
    this->yPanels = 0;

    this->hx      = 0.0;
    this->hy      = 0.0;

    this->xMin = 0.0;
    this->xMax = 1.0;
    this->yMin = 0.0;
    this->yMax = 1.0;

    this->XYperiodicityFlag = false;
    this->XperiodicityFlag  = false;
    this->YperiodicityFlag  = false;
}

void initialize(const GridFunction2d& G)
{
	DoubleVector2d::initialize(G);
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;

    this->XYperiodicityFlag = G.XYperiodicityFlag;
    this->XperiodicityFlag  = G.XperiodicityFlag;
    this->YperiodicityFlag  = G.YperiodicityFlag;

}

void initialize(long xPanels, double hx, long yPanels, double hy)
{
	DoubleVector2d::initialize(xPanels+1,yPanels+1);

    this->xPanels = xPanels;
    this->yPanels = yPanels;

    this->hx     = hx;
    this->hy     = hy;

    this->xMin = -(xPanels*hx)/2.0;
    this->xMax =  (xPanels*hx)/2.0;
    this->yMin = -(yPanels*hy)/2.0;
    this->yMax =  (yPanels*hy)/2.0;

    this->XYperiodicityFlag = false;
    this->XperiodicityFlag  = false;
    this->YperiodicityFlag  = false;

	this->setToValue(0.0);
}

void initialize(long xPanels, double xMin, double xMax, long yPanels, double yMin, double yMax)
{
    DoubleVector2d::initialize(xPanels+1,yPanels+1);

    this->xPanels = xPanels;
    this->yPanels = yPanels;

    this->xMin = xMin;
    this->xMax = xMax;
    this->yMin = yMin;
    this->yMax = yMax;

    this->hx     = (xMax-xMin)/(double)(xPanels);
    this->hy     = (yMax-yMin)/(double)(yPanels);

    this->XYperiodicityFlag = false;
    this->XperiodicityFlag  = false;
    this->YperiodicityFlag  = false;

	this->setToValue(0.0);
}


GridFunction2d* newDuplicate() const
{
    GridFunction2d* Mptr = new GridFunction2d(*this);
    return Mptr;
}

bool isNull() const
{
if(dataPtr == nullptr) return true;
return false;
}

double getHx()   const  {return hx;}
double getXmin() const  {return   xMin;}
double getXmax() const  {return   xMax;}
long   getXpanelCount() const {return xPanels;}

double getHy()   const  {return  hy;}
double getYmin() const  {return yMin;}
double getYmax() const  {return yMax;}
long   getYpanelCount() const {return yPanels;}

virtual bool isXYperiodic() const {return XYperiodicityFlag;}
virtual bool isXperiodic()  const {return XperiodicityFlag;}
virtual bool isYperiodic()  const {return YperiodicityFlag;}

DoubleVector2d getValues() const
{
    return DoubleVector2d(*this);
}

DoubleVector2d* getValuesPointer()
{
    return (DoubleVector2d*)this;
}

const DoubleVector2d* getValuesPointer() const
{
    return (DoubleVector2d*)this;
}


void setPeriodicity(bool val = true)
{
	if(val)
	{
		XYperiodicityFlag = true;
	    enforcePeriodicity();
	}
}

void clearPeriodicity()
{
	XYperiodicityFlag = false;
}

void setXYperiodicity(bool val = true)
{
	if(val)
	{
		XYperiodicityFlag = true;
	    enforcePeriodicity();
	}
}

void clearXYperiodicity()
{
	XYperiodicityFlag = false;
}

void setXperiodicity(bool val = true)
{
	if(val)
	{
		XperiodicityFlag = true;
	    enforceXperiodicity();
	}
}

void clearXperiodicity()
{
	XperiodicityFlag = false;
}

void setYperiodicity(bool val = true)
{
	if(val)
	{
		YperiodicityFlag = true;
	    enforceYperiodicity();;
	}
}

void clearYperiodicity()
{
	YperiodicityFlag = false;
}


/*!  Returns the size of the number of independent values
     associated with the grid function, i.e. the dimension
     corresponding to the vector of independent function
     values.

     This dimension depends upon the specified boundary values.
*/

virtual long getDimension() const
{
	long dimension = DoubleVector2d::getDimension();

    if(not XYperiodicityFlag)
    {
    	dimension = DoubleVector2d::getDimension();
    }
    else if(XYperiodicityFlag)
    {
    	dimension = (index1Size-1)*(index2Size-1);
    }
    else if(XperiodicityFlag)
    {
    	dimension = (index1Size-1)*(index2Size);
    }
    else if(YperiodicityFlag)
    {
    	dimension = (index1Size)*(index2Size-1);
    }

    return dimension;
}


//###################################################################
//                  Values Access
//
//       An alternative to using operator()(...)
//###################################################################

#ifdef _DEBUG
    double& Values(long i1, long i2)
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    return *(dataPtr +  i2 + i1*index2Size);
    };

    const double&  Values(long i1, long i2) const
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    return *(dataPtr +   i2  + i1*index2Size);
    };
#else
    /*!
    Returns a reference to the element with index (i1,i2) - indexing
    starting at (0,0).
    */
    inline double&  Values(long i1, long i2)
    {
    return *(dataPtr +  i2 + i1*index2Size);
    };

    /*!
    Returns a reference to the element with index (i1,i2) - indexing
    starting at (0,0).
     */
    inline const double&  Values(long i1, long i2) const
    {
    return *(dataPtr +  i2  + i1*index2Size);
    };
#endif



GridFunction2d& operator=(const GridFunction2d& G)
{
//
// Only propagate domain information across assignment if the
// left hand side is a null instance.
//
if(this->dataPtr == nullptr)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;

    this->XYperiodicityFlag = G.XYperiodicityFlag;
}

// Propagate values

DoubleVector2d::operator=(G);

return *this;
}

GridFunction2d& operator=(GridFunction2d&& G)
{
//
// Only propagate domain information across assignment if the
// left hand side is a null instance.
//
if(this->dataPtr == nullptr)
{
    this->xPanels = G.xPanels;
    this->yPanels = G.yPanels;

    this->hx     = G.hx;
    this->hy     = G.hy;

    this->xMin = G.xMin;
    this->xMax = G.xMax;
    this->yMin = G.yMin;
    this->yMax = G.yMax;

    this->XYperiodicityFlag = G.XYperiodicityFlag;
}

// Propagate values

DoubleVector2d::operator=((DoubleVector2d&&)G);
return *this;
}

GridFunction2d& operator=(DoubleVector2d& G)
{
//
// Since underlying operations that result in DoubleVector2d&&G strip
// the domain information an assignment to a null instance generates
// an error.

if(this->dataPtr == nullptr)
{
	assert(domainCheck());
}

DoubleVector2d::operator=((DoubleVector2d&)G);
return *this;
}

GridFunction2d& operator=(DoubleVector2d&& G)
{
//
// Since underlying operations that result in DoubleVector2d&&G strip
// the domain information an assignment to a null instance generates
// an error.

if(this->dataPtr == nullptr)
{
	assert(domainCheck());
}

DoubleVector2d::operator=((DoubleVector2d&&)G);
return *this;
}

//####  Incremental operators  ###


void operator*=(double alpha)
{SCC::DoubleVector2d::operator*=(alpha);}

void operator/=(double alpha)
{SCC::DoubleVector2d::operator/=(alpha);}


void operator+=(const GridFunction2d& G)
{SCC::DoubleVector2d::operator+=(G);}

void operator-=(const GridFunction2d& G)
{SCC::DoubleVector2d::operator-=(G);}

void operator*=(const GridFunction2d& G)
{SCC::DoubleVector2d::operator*=(G);}

void operator/=(const GridFunction2d& G)
{SCC::DoubleVector2d::operator/=(G);}



void operator+=(const DoubleVector2d& G)
{SCC::DoubleVector2d::operator+=(G);}

void operator-=(const DoubleVector2d& G)
{SCC::DoubleVector2d::operator-=(G);}

void operator*=(const DoubleVector2d& G)
{SCC::DoubleVector2d::operator*=(G);}

void operator/=(const DoubleVector2d& G)
{SCC::DoubleVector2d::operator/=(G);}


void operator*=(const  std::function<double(double,double)>& F)
{
    double xPos; double yPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    Values(i,j) *= F(xPos,yPos);
    }}
}

void operator/=(const std::function<double(double,double)>& F)
{
    double xPos; double yPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    Values(i,j) /= F(xPos,yPos);
    }}
}


void operator+=(const  std::function<double(double,double)>& F)
{
    double xPos; double yPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    Values(i,j) += F(xPos,yPos);
    }}
}


void operator-=(const  std::function<double(double,double)>& F)
{
    double xPos; double yPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    Values(i,j) -= F(xPos,yPos);
    }}
}

//############################################################################

void createProductFunction(const GridFunction1d& funX, const GridFunction1d& funY)
{
    initialize(funX.xPanels, funX.xMin, funX.xMax,
               funY.xPanels, funY.xMin, funY.xMax);

	long i; long j;

	double fX; double fY;

	for(i = 0; i <= xPanels; i++)
	{
	fX = funX.Values(i);
	for(j = 0; j <= yPanels; j++)
	{
	fY = funY.Values(j);
	Values(i,j) = fX*fY;
	}}

	if(XYperiodicityFlag) {enforcePeriodicity();}
	if(XperiodicityFlag)  {enforceXperiodicity();}
	if(YperiodicityFlag)  {enforceYperiodicity();}
}


void specify(const std::function<double(double,double)>& F)
{
    double xPos; double yPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    for(long j = 0; j <= yPanels; j++)
    {
    yPos = yMin + j*hy;
    Values(i,j) = F(xPos,yPos);
    }}

	if(XYperiodicityFlag) {enforcePeriodicity();}
	if(XperiodicityFlag)  {enforceXperiodicity();}
	if(YperiodicityFlag)  {enforceYperiodicity();}
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
using DoubleVector2d::dot;
    
virtual double dot(const GridFunction2d& V) const
{
	return scaledDot(V);
}


virtual double scaledDot(const GridFunction2d& V) const
{
    double dotVal = 0.0;

    long i; long j;
    double hxhy = hx*hy;

//  Interior points weighted by 1
//
    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j)*V.Values(i,j)*hxhy;
    }}
//
//  Interior face points (weighted by 1/2)
//
    i = 0;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.5;
    }

    i = xPanels;
    for(j = 1; j < yPanels; j++)
    {
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.5;
    }

    //

    j = 0;
    for(i = 1; i < xPanels; i++)
    {
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.5;
    }

    j = yPanels;
    for(i = 1; i <  xPanels; i++)
    {
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.5;
    }

//
//  Corner points (weighted by 1/4)
//
    i = 0;
    j = 0;
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.25;


    i = 0;
    j = yPanels;
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.25;

    i = xPanels;
    j = 0;
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.25;

    i = xPanels;
    j = yPanels;
    dotVal += Values(i,j)*V.Values(i,j)*hxhy*0.25;

    return dotVal;
}

// Trapezoidal method approximation of the integral of F(Values(i,j))

double integralTrapezoidal(std::function<double(double)> F) const
{
    double intVal = 0.0;

    long i; long j;
    double hxhy = hx*hy;

//  Interior points weighted by 1

    for(i = 1; i < xPanels; i++)
    {
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j))*hxhy;
    }}

//  Interior face points (weighted by 1/2)

    i = 0;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j))*hxhy*0.5;
    }

    i = xPanels;
    for(j = 1; j < yPanels; j++)
    {
    intVal += F(Values(i,j))*hxhy*0.5;
    }

    //

    j = 0;
    for(i = 1; i < xPanels; i++)
    {
    intVal += F(Values(i,j))*hxhy*0.5;
    }

    j = yPanels;
    for(i = 1; i <  xPanels; i++)
    {
    intVal += F(Values(i,j))*hxhy*0.5;
    }

//
//  Corner points (weighted by 1/4)
//
    i = 0;
    j = 0;
    intVal += F(Values(i,j))*hxhy*0.25;

    i = 0;
    j = yPanels;
    intVal += F(Values(i,j))*hxhy*0.25;

    i = xPanels;
    j = 0;
    intVal += F(Values(i,j))*hxhy*0.25;

    i = xPanels;
    j = yPanels;
    intVal += F(Values(i,j))*hxhy*0.25;

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

double nrm2() const
{
    return norm2();
}

double norm2() const
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

    for(long k =0; k < (xPanels+1)*(yPanels+1); k++)
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
    return avgValue/((xMax-xMin)*(yMax-yMin));
}

//  Average computed using Trapezoidal approximation of the integral

double getTrapezoidalAverage() const
{
    double avgValue;
    avgValue = integralTrapezoidal();
    return avgValue/((xMax-xMin)*(yMax-yMin));
}

double min() const
{
//  Compute function minimum

    double minValue  = dataPtr[0];
    for(long k = 1; k < (xPanels+1)*(yPanels+1); k++)
    {
    minValue = ( minValue < dataPtr[k] ) ? minValue : dataPtr[k];
    }
    return minValue;
}

double max() const
{
//  Compute function maximum

    double maxValue  = dataPtr[0];
    for(long k = 1; k < (xPanels+1)*(yPanels+1); k++)
    {
    maxValue = ( maxValue > dataPtr[k] ) ? maxValue : dataPtr[k];
    }
    return maxValue;
}

GridFunction1d getConstantYslice(long yIndex) const  // (x function)
{
	GridFunction1d R(xPanels,xMin,xMax);
	for(long i = 0; i <= xPanels; i++)
	{
	R.Values(i) = Values(i,yIndex);
	}
	return R;
}

void getConstantYslice(long yIndex, GridFunction1d& R) const  // (x function)
{
	for(long i = 0; i <= xPanels; i++)
	{
	R.Values(i) = Values(i,yIndex);
	}
}

void setConstantYslice(long yIndex, const GridFunction1d& R)  // (x function)
{
	for(long i = 0; i <= xPanels; i++)
	{
	Values(i,yIndex) = R.Values(i);
	}
}

GridFunction1d getConstantXslice(long xIndex) const  //( y function)
{
	GridFunction1d R(yPanels,yMin,yMax);
	for(long j = 0; j <= yPanels; j++)
	{
	R.Values(j) = Values(xIndex,j);
	}
	return R;
}

void getConstantXslice(long xIndex, GridFunction1d& R) const  //( y function)
{
	for(long j = 0; j <= yPanels; j++)
	{
	R.Values(j) = Values(xIndex,j);
	}
}

void setConstantXslice(long xIndex, const GridFunction1d& R)    //( y function)
{
	for(long j = 0; j <= yPanels; j++)
	{
	Values(xIndex,j) = R.Values(j);
	}
}



void setBoundaryValues(double value)
{
    long i; long j;

	i = 0;
	for(j = 0; j <= yPanels; j++)
	{
		Values(i,j) = value;
	}

    i = xPanels;
    for(j = 0; j <= yPanels; j++)
    {
     Values(i,j) = value;
    }

    j = 0;
    for(i = 0; i <= xPanels; i++)
    {
    Values(i,j) = value;
    }

    j = yPanels;
    for(i = 0; i <= xPanels; i++)
    {
    Values(i,j) = value;
    }
}

// Sets the periodicity flag and copies values
// to enforce periodicity in both coordinate directions.
//
// The values along the lower left corner
// boundary edges (xMin,y) and (x,yMin) are
// propagated to the upper right boundary edges
// (xMax,y) and (x,yMax) to enforce periodicity.
//

void enforceXYperiodicity()
{
	enforcePeriodicity();
}

void enforcePeriodicity()
{
	XYperiodicityFlag = true;

	if(this->isNull()){return;}

    long i; long j;

    i = xPanels;
    for(j = 0; j <= yPanels; j++)
    {
     Values(i,j) = Values(0,j);
    }

    j = yPanels;
    for(i = 0; i <= xPanels; i++)
    {
    Values(i,j) = Values(i,0);
    }
}

void enforceXperiodicity()
{
	XperiodicityFlag = true;

	if(this->isNull()){return;}

    long i; long j;

    i = xPanels;
    for(j = 0; j <= yPanels; j++)
    {
     Values(i,j) = Values(0,j);
    }
}

void enforceYperiodicity()
{
	YperiodicityFlag = true;

	if(this->isNull()){return;}

    long i; long j;

    j = yPanels;
    for(i = 0; i <= xPanels; i++)
    {
    Values(i,j) = Values(i,0);
    }
}

//  Returns true if the input grid function is structurally identical,
//  e.g. a coincident computational domain and mesh size
//

bool isCoincident(const GridFunction2d& V)
{
	double domainDiffTol    = 1.0e-13;

    long discretizationDiff = std::abs(xPanels - V.xPanels)
                            + std::abs(yPanels - V.yPanels);

    double domainDiff = (std::abs(xMin - V.xMin) + std::abs(xMax - V.xMax))/std::abs(xMax-xMin)
                      + (std::abs(yMin - V.yMin) + std::abs(yMax - V.yMax))/std::abs(yMax-yMin);

    if((discretizationDiff > 0) || (domainDiff > domainDiffTol)) return false;
    return true;
}

/*! BLAS axpby : this  <- alpha*v + beta*this  */

   void axpby(double alpha, const DoubleVector2d& G, double beta)
   {
   SCC::DoubleVector2d::axpby(alpha,G,beta);
   }

/*! BLAS axpy : this  <- alpha*v + this  */

   void axpy(double alpha, const DoubleVector2d& G)
   {
   SCC::DoubleVector2d::axpy(alpha,G);
   }

//  Grid Geometry

double xMin; double xMax;   // Computational Region is [xMin,xMax]x
double yMin; double yMax;   //                         [yMin,yMax]

long   xPanels;             // Number of panels
long   yPanels;


double hx;                  // mesh width
double hy;

bool XYperiodicityFlag;     // indicates periodicity in the XY direction
bool XperiodicityFlag;      // indicates periodicity in the X direction
bool YperiodicityFlag;      // indicates periodicity in the Y direction

//###################################################################
//                     Error Checking
//###################################################################

protected:

#ifdef _DEBUG
        bool domainCheck() const
        {
        std::cerr << "XXX SCC::GridFunction2d Error  XXX" << std::endl;
        std::cerr << "Left side of assignment must be a non-null GridFunction2d instance." << std::endl;
        std::cerr << std::endl;
        return false;
        }
#else
        bool domainCheck() const {return true;}
#endif

};

}


#endif /* SCC_GRIDFUNCTIONSCC_GRIDFUNCTION1D_H_ */
