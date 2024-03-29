/*
 * SCC_GridFunction1d.h
 *
 *  Created on: Jun 27, 2015
 *      Author: anderson
 *
 *  Release : 18.09.04
 *
 *
 *
 * Decisions: Extending DoubleVector1D so that move semantics can be incorporated into
 * the underlying vector operations on grid values without having to duplicate all
 * member functions utilizing move semantics.
 *
 * Providing a Values member function so that data values can be accessed as if this
 * grid function implementation contained an instance of DoubleVector1D Values.
 *
 * Using lambda functions to specify operations on grid values
 * Use of trapezoidal method when creating the scaled dot product and integral norms
 *
 * Revised: Nov. 26, 2015 
 *          Jan. 26, 2016
 */
/*
#############################################################################
#
# Copyright  2015-16 Chris Anderson
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

#ifndef SCC_GRID_FUNCTION_1D_
#define SCC_GRID_FUNCTION_1D_

#include <functional>
#include <iostream>
#include <cmath>

#ifdef _MSC_VER
#include "iso646.h"          // So "and" is equivalenced to &&
typedef unsigned int uint;   // Define uint to be unsigned int
#undef min
#undef max
#endif

#include "../DoubleVectorNd/SCC_DoubleVector1d.h"

namespace SCC
{

class GridFunction1d : public DoubleVector1d
{

public :

GridFunction1d() : DoubleVector1d()
{
	this->xMin    = 0.0;
	this->xMax    = 1.0;
	this->xPanels = 0;
	this->hx      = 0.0;

	this->XperiodicityFlag         = false;
	this->dirichletDimensionFlag   = false;
}

GridFunction1d(const GridFunction1d& G) : DoubleVector1d(G)
{
	this->xMin    = G.xMin;
	this->xMax    = G.xMax;
	this->xPanels = G.xPanels;
	this->hx      = G.hx;

    this->XperiodicityFlag        = G.XperiodicityFlag;
    this->dirichletDimensionFlag  = G.dirichletDimensionFlag;
}

GridFunction1d(DoubleVector1d&& G) : DoubleVector1d((DoubleVector1d&&)G)
{
	this->xMin    = 0.0;
	this->xMax    = 1.0;
	this->xPanels = 0;
	this->hx      = 0.0;

	this->XperiodicityFlag         = false;
	this->dirichletDimensionFlag   = false;
}


GridFunction1d(long xPanels, double hx) : DoubleVector1d(xPanels+1)
{
    this->xPanels =  xPanels;
    this->hx      =  hx;
    this->xMin    = -(xPanels*hx)/2.0;
    this->xMax    =  (xPanels*hx)/2.0;
	this->setToValue(0.0);

	this->XperiodicityFlag        = false;
	this->dirichletDimensionFlag   = false;
}


GridFunction1d(long xPanels, double xMin, double xMax) : DoubleVector1d(xPanels+1)
{
	this->xMin    = xMin;
	this->xMax    = xMax;
	this->xPanels = xPanels;
	this->hx      = (xMax-xMin)/(double)xPanels;
	this->setToValue(0.0);

	this->XperiodicityFlag         = false;
	this->dirichletDimensionFlag   = false;
}


GridFunction1d(GridFunction1d&& G) : DoubleVector1d((DoubleVector1d&&)G)
{
	this->xMin    = G.xMin;
	this->xMax    = G.xMax;
	this->xPanels = G.xPanels;
	this->hx      = G.hx;

	this->XperiodicityFlag         = G.XperiodicityFlag;
	this->dirichletDimensionFlag   = G.dirichletDimensionFlag;
}

virtual ~GridFunction1d(){}

// Initialization

void initialize()
{
    DoubleVector1d::initialize();
	xMin    = 0.0;
	xMax    = 1.0;
	xPanels = 0;
	hx      = 0.0;

	XperiodicityFlag         = false;
    dirichletDimensionFlag   = false;
}

void initialize(const GridFunction1d& G)
{
	DoubleVector1d::initialize(G);
	xMin    = G.xMin;
	xMax    = G.xMax;
	xPanels = G.xPanels;
	hx      = G.hx;

	XperiodicityFlag         = G.XperiodicityFlag;
	dirichletDimensionFlag   = G.dirichletDimensionFlag;
}

void initialize(long xPanels, double hx)
{
	DoubleVector1d::initialize(xPanels+1);
	this->xMin    = 0.0;
	this->xMax    = xPanels*hx;
	this->xPanels = xPanels;
	this->hx      = hx;

	this->XperiodicityFlag         = false;
	this->dirichletDimensionFlag   = false;
}

void initialize(long xPanels, double xMin, double xMax)
{
	DoubleVector1d::initialize(xPanels+1);
	this->xMin    = xMin;
	this->xMax    = xMax;
	this->xPanels = xPanels;
	this->hx      = (xMax-xMin)/(double)xPanels;

	this->XperiodicityFlag         = false;
    this->dirichletDimensionFlag   = false;
}

GridFunction1d* newDuplicate() const
{
    GridFunction1d* Mptr = new GridFunction1d(*this);
    return Mptr;
}

bool isNull() const
{
if(dataPtr == nullptr) return true;
return false;
}

long   getXpanelCount() const {return xPanels;}
double getHx()          const  {return     hx;}
double getXmin()        const  {return   xMin;}
double getXmax()        const  {return   xMax;}

DoubleVector1d getValues() const
{
    return DoubleVector1d(*this);
}

DoubleVector1d* getValuesPointer()
{
    return (DoubleVector1d*)this;
}

const DoubleVector1d* getValuesPointer() const
{
    return (DoubleVector1d*)this;
}

virtual bool isXperiodic()  const {return XperiodicityFlag;}


void setPeriodicity(bool val = true)
{
	if(val)
	{
		XperiodicityFlag = true;
	    enforcePeriodicity();
	}
}

void clearPeriodicity()
{
	XperiodicityFlag = false;
}


void setDirichletDimensionFlag(bool val = true)
{
     dirichletDimensionFlag = val;
     if(val)
     {XperiodicityFlag = false;}
}

void clearDirichletDimensionFlag()
{
	 dirichletDimensionFlag = false;
}


void setXperiodicity(bool val = true)
{
	if(val)
	{
		XperiodicityFlag = true;
	    enforcePeriodicity();
	}
}

void clearXperiodicity()
{
	XperiodicityFlag = false;
}

/*!  Returns the size of the number of independent values
     associated with the grid function, i.e. the dimension
     corresponding to the vector of independent function
     values.

     This dimension depends upon the specified boundary values.
*/

virtual long getDimension() const
{
	long dimension = DoubleVector1d::getDimension();

    if(XperiodicityFlag)
    {
    	dimension -=1;
    }
    if(dirichletDimensionFlag)
    {
         dimension -= 2;
    }
    return dimension;
}



//###################################################################
//                  Values Access
//
//       An alternative to using operator()(...)
//###################################################################


#ifndef NDEBUG
    double&  Values(long i1)
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    return *(dataPtr +  i1);
    };

    const double&  Values(long i1) const
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    return *(dataPtr +  i1);
    };

#else
    inline double&  Values(long i1)
    {
    return *(dataPtr + i1);
    };

    inline const double&  Values(long i1) const
    {
    return *(dataPtr + i1);
    };
#endif

GridFunction1d& operator=(const GridFunction1d& G)
{
//
// Only propagate domain information across assignment if the
// left hand side is a null instance.
//
if(this->dataPtr == nullptr)
{
	xMin    = G.xMin;
	xMax    = G.xMax;
	xPanels = G.xPanels;
	hx      = G.hx;

    XperiodicityFlag         = G.XperiodicityFlag;
	dirichletDimensionFlag   = G.dirichletDimensionFlag;
}

// Propagate values

DoubleVector1d::operator=(G);

return *this;
}

GridFunction1d& operator=(GridFunction1d&& G)
{
//
// Only propagate domain information across assignment if the
// left hand side is a null instance.
//
if(this->dataPtr == nullptr)
{
	xMin    = G.xMin;
	xMax    = G.xMax;
	xPanels = G.xPanels;
	hx      = G.hx;

    XperiodicityFlag         = G.XperiodicityFlag;
	dirichletDimensionFlag   = G.dirichletDimensionFlag;
}

// Propagate values

DoubleVector1d::operator=((DoubleVector1d&&)G);
return *this;
}


GridFunction1d& operator=(DoubleVector1d& G)
{
//
// Since underlying operations that result in DoubleVector1d&G strip
// the domain information an assignment to a null instance generates
// an error.

if(this->dataPtr == nullptr)
{
	assert(domainCheck());
}

DoubleVector1d::operator=((DoubleVector1d&)G);
return *this;
}

GridFunction1d& operator=(DoubleVector1d&& G)
{
//
// Since underlying operations that result in DoubleVector1d&&G strip
// the domain information an assignment to a null instance generates
// an error.

if(this->dataPtr == nullptr)
{
	assert(domainCheck());
}

DoubleVector1d::operator=((DoubleVector1d&&)G);
return *this;
}

//####  Incremental operators  ###

void operator*=(double alpha)
{SCC::DoubleVector1d::operator*=(alpha);}

void operator/=(double alpha)
{SCC::DoubleVector1d::operator/=(alpha);}


void operator+=(const GridFunction1d& G)
{SCC::DoubleVector1d::operator+=(G);}

void operator-=(const GridFunction1d& G)
{SCC::DoubleVector1d::operator-=(G);}

void operator*=(const GridFunction1d& G)
{SCC::DoubleVector1d::operator*=(G);}

void operator/=(const GridFunction1d& G)
{SCC::DoubleVector1d::operator/=(G);}



void operator+=(const DoubleVector1d& G)
{SCC::DoubleVector1d::operator+=(G);}

void operator-=(const DoubleVector1d& G)
{SCC::DoubleVector1d::operator-=(G);}

void operator*=(const DoubleVector1d& G)
{SCC::DoubleVector1d::operator*=(G);}

void operator/=(const DoubleVector1d& G)
{SCC::DoubleVector1d::operator/=(G);}


void operator*=(const std::function<double(double)>& F)
{
    double xPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    Values(i) *= F(xPos);
    }
}

void operator/=(const std::function<double(double)>& F)
{
    double xPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    Values(i) /= F(xPos);
    }
}

void operator+=(const std::function<double(double)>& F)
{
    double xPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    Values(i) += F(xPos);
    }
}

void operator-=(const std::function<double(double)>& F)
{
    double xPos;
	for(long i = 0; i <= xPanels; i++)
    {
    xPos = xMin + i*hx;
    Values(i) -= F(xPos);
    }
}

//#########################################################


void specify(std::function<double(double)> F)
{
	double xPos;
	for(long k = 0; k <= xPanels; k++)
	{
	xPos = xMin + k*hx;
	Values(k) = F(xPos);
	}
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
// values times the mesh size.
//
using DoubleVector1d::dot;
    
virtual double dot(const GridFunction1d& V) const
{
	return scaledDot(V);
}

// Dot product scaled with mesh widths
//
// Trapezoidal approximation of the integral
//

virtual double scaledDot(const GridFunction1d& V) const
{
    double intValue  = 0.0;
    long  i;

    for(i = 1; i <= xPanels-1; i++)
    {
        intValue += Values(i)*V.Values(i)*hx;
    }
    intValue += Values(0)*V.Values(0)*(hx/2.0);
    intValue += Values(xPanels)*V.Values(xPanels)*(hx/2.0);

    return intValue;
}
//
// Trapezoidal method applied to F(Values(i))
// (Accepts lambda functions as arguments)
//

double integralTrapezoidal(std::function<double(double)> F) const
{
    double intValue  = 0.0;
    long  i;

    // Trapezoidal method
    for(i = 1; i <= xPanels-1; i++)
    {
        intValue += F(Values(i))*hx;
    }
    intValue += F(Values(0))*(hx/2.0);
    intValue += F(Values(xPanels))*(hx/2.0);

    return intValue;
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


// Riemann sum approximation of integral

double integralRiemann() const
{
    double intValue  = 0.0;
    long  i;

    for(i = 0; i <= xPanels; i++)
    {
        intValue += Values(i)*hx;
    }

    return intValue;
}

//  Average computed using Riemann sum approximation

double getRiemannAverage() const
{
    double avgValue;
    avgValue = integralRiemann();
    return avgValue/(xMax-xMin);
}

//  Average computed using Trapezoidal approximation of the integral

double getTrapezoidalAverage() const
{
    double avgValue;
    avgValue = integralTrapezoidal();
    return avgValue/(xMax-xMin);
}


double min() const
{
//  Compute function minimum

    double minValue  = dataPtr[0];
    for(long k = 1; k <= xPanels; k++)
    {
    minValue = ( minValue < dataPtr[k] ) ? minValue : dataPtr[k];
    }
    return minValue;
}

double max() const
{

//  Compute function maximum

    double maxValue  = dataPtr[0];
    for(long k = 1; k <= xPanels; k++)
    {
    maxValue = ( maxValue > dataPtr[k] ) ? maxValue : dataPtr[k];
    }
    return maxValue;
}

void setBoundaryValues(double value)
{
    Values(0)       = value;
    Values(xPanels) = value;
}

void enforcePeriodicity()
{
	XperiodicityFlag = true;
	if(this->isNull()){return;}

    Values(xPanels) = Values(0);
}



//  Returns true if the input grid function is structurally identical,
//  e.g. a coincident computational domain and mesh size
//

bool isCoincident(const GridFunction1d& V)
{
	double domainDiffTol    = 1.0e-13;

    long discretizationDiff = std::abs(xPanels - V.xPanels);

    double domainDiff = (std::abs(xMin - V.xMin) + std::abs(xMax - V.xMax))/std::abs(xMax-xMin);

    if((discretizationDiff > 0) || (domainDiff > domainDiffTol)) return false;
    return true;
}

/*! BLAS axpby : this  <- alpha*v + beta*this  */

   void axpby(double alpha, const DoubleVector1d& G, double beta)
   {
   SCC::DoubleVector1d::axpby(alpha,G,beta);
   }

/*! BLAS axpy : this  <- alpha*v + this  */

   void axpy(double alpha, const DoubleVector1d& G)
   {
   SCC::DoubleVector1d::axpy(alpha,G);
   }

//
//  Grid Geometry
//
    double xMin; double xMax;      // Computational Region is [xMin, xMax]
    long   xPanels;                // Number of panels
    double hx;                     // mesh width

    bool   XperiodicityFlag;
    bool   dirichletDimensionFlag; // Dimension of unknowns that excludes boundary values

protected :

//###################################################################
//                     Error Checking
//###################################################################

#ifdef _DEBUG
        bool domainCheck() const
        {
        std::cerr << "XXX SCC::GridFunction1d Error  XXX" << std::endl;
        std::cerr << "Left side of assignment must be a non-null GridFunction1d instance." << std::endl;
        std::cerr << std::endl;
        return false;
        }
#else
        bool domainCheck() const {return true;}
#endif


};

}


#endif /* SCC_GRIDFUNCTIONSCC_GRIDFUNCTION1D_H_ */
