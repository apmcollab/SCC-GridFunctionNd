
/*
 * SCC_GridFunction2dUtility.h
 *
 *  Created on: Jun 28, 2015
 *      Author: anderson
 *
 *
 * Utility routines - mostly routines for file input and output
 * of GridFunction2d instances.
 *
 * Modifications: Jun 28, 2015
 *               Sept. 4, 2018
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

#ifndef _SCC_GridFunction2dUtility_
#define _SCC_GridFunction2dUtility_
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <functional>
using namespace std;

// MS compilers generate warnings if fopen is used instead of fopen_s (a Microsoft specific language
// extension, so this macro implements the appropriate equivalent to fopen that MS wants when
// MS compilers are being used. In both versions, the
// macro returns a non-zero value if the open fails (e.g. a non-zero error code).
//
#ifndef _MSC_VER
#ifndef OPENFILE
#define OPENFILE(dataFile,filename,mode) ((dataFile = fopen(filename,  mode)) == NULL)
#endif
#else
#ifndef OPENFILE
#define OPENFILE(dataFile,fileName,mode) ((fopen_s(&dataFile,fileName, mode)) != 0)
#pragma warning(push)
#pragma warning(disable: 4996)
#endif
#endif


#include "../GridFunctionNd/SCC_GridFunction2d.h"
#include "../DoubleVectorNd/SCC_DoubleVector2d.h"

namespace SCC
{
class  GridFunction2dUtility
{
public :

// Adds the values of F to the input gFun

void addToValues(GridFunction2d& gFun, std::function<double(double,double)>& F)
{
    long i; long j;

    double xMin     = gFun.getXmin();
    double yMin     = gFun.getYmin();

    long xPanels = gFun.getXpanelCount();
    long yPanels = gFun.getYpanelCount();

    double  hx   = gFun.getHx();
    double  hy   = gFun.getHy();

    double x; double y;

    for(i = 0; i <= xPanels; i++)
    {
    x = xMin + i*hx;
    for(j = 0; j <= yPanels; j++)
    {
    y = yMin + j*hy;
    gFun.Values(i,j) += F(x,y);
    }}
}
//
// Gnuplot output data structure
//
// # xPanels xMin xMax yPanels yMin yMax
// x_0 y_0       f(x_0, y_0)
// x_0 y_1       f(x_0, y_1)
//          ***
// x_0 y_yPanels f(x_0, y_yPanels)
//
// x_1 y_0       f(x_1, y_0)
// x_1 y_1       f(x_1, y_1))
//
void outputToGNUplot(GridFunction2d& gF, const string& fileName, const string& formatString = "%20.15e")
{
    ostringstream s;
    s.str("");
    s << formatString << "  " << formatString << "  " << formatString << " \n";
//
//  Open and then write to a file
//
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    long xPanel   = gF.getXpanelCount();
    double xMin   = gF.getXmin();
    double xMax   = gF.getXmax();

    long yPanel   = gF.getYpanelCount();
    double yMin   = gF.getYmin();
    double yMax   = gF.getYmax();

    long i; long j;


    fprintf(dataFile,"# %ld %20.15e %20.15e  %ld  %20.15e %20.15e \n",
            xPanel,xMin,xMax,yPanel, yMin, yMax);

    double x;  double y;

    double hx = (xMax-xMin)/(double)(xPanel);
    double hy = (yMax-yMin)/(double)(yPanel);

    for(i = 0; i <= xPanel; i++)
    {
    x = xMin + double(i)*hx;
	for(j = 0; j <= yPanel; j++)
	{
    y = yMin + double(j)*hy;
    fprintf(dataFile,(s.str()).c_str(),x,y,gF(i,j));
    }
    fprintf(dataFile,"\n");
    }

    fclose(dataFile);
}


 void inputFromGNUplot(GridFunction2d& gF, const string& fileName)
{
//
//  Open and then read from file
//
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "r" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    size_t rValue = 0;

    char poundChar;

    long i;       long j;
    long xPanels; long yPanels;

    double xMin; double xMax;
    double yMin; double yMax;

    rValue = fscanf(dataFile,"%c", &poundChar)  != 1 ? 1 : rValue;   // remove leading #
    rValue = fscanf(dataFile,"%ld",&xPanels)    != 1 ? 1 : rValue;
    rValue = fscanf(dataFile,"%lf",&xMin)       != 1 ? 1 : rValue;
    rValue = fscanf(dataFile,"%lf",&xMax)       != 1 ? 1 : rValue;
    rValue = fscanf(dataFile,"%ld",&yPanels)    != 1 ? 1 : rValue;
    rValue = fscanf(dataFile,"%lf",&yMin)       != 1 ? 1 : rValue;
    rValue = fscanf(dataFile,"%lf",&yMax)       != 1 ? 1 : rValue;

    rValue = (xMax < xMin) ? 1 : rValue;
    rValue = (yMax < yMin) ? 1 : rValue;

    rValue = (xPanels <= 0) ? 1 : rValue;
    rValue = (yPanels <= 0) ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction2d could not be initialized from gnuplot data file " + fileName + " \n");
    }


    gF.initialize(xPanels,xMin,xMax,yPanels, yMin, yMax);

    double x; double y;
    for(i = 0;  i <= xPanels; i++)
    {
	for(j = 0; j  <= yPanels; j++)
	{
	rValue = fscanf(dataFile,"%lf %lf %lf",&x,&y,&gF(i,j)) != 3 ? 1 : rValue;
    }
    }

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction2d could not be initialized from gnuplot data file " + fileName + " \n");
    }

    fclose(dataFile);
}



void outputDataToVTKfile(const GridFunction2d& gridFun, const string& fileName, const string& dataLabel)
{
	FILE* dataFile;

    double a  = gridFun.getXmin();
    double c  = gridFun.getYmin();
	double e  = 0.0;

    double hx = gridFun.getHx();
    double hy = gridFun.getHy();
	double hz = hx;
	hz        = (hz < hy) ? hy : hz;

    long mPt = gridFun.getXpanelCount() + 1;
    long nPt = gridFun.getYpanelCount() + 1;
	long pPt = 2;

    long dataCount = mPt*nPt*pPt;

    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      printf( "The file %s could not be  opened\n",fileName.c_str());
      exit(1);
    }

    long i; long j; long k;
    double xPos; double yPos; double zPos;
    //
    // output the regular positions
    //
	fprintf(dataFile, "# vtk DataFile Version 2.0\n");
	fprintf(dataFile, "%s \n",dataLabel.c_str());
    fprintf(dataFile, "ASCII\n");

    fprintf(dataFile, "DATASET RECTILINEAR_GRID\n");
    fprintf(dataFile, "DIMENSIONS %ld %ld %ld \n",mPt,nPt,pPt);
    fprintf(dataFile, "X_COORDINATES %ld float \n",mPt);
    for(i = 0; i < mPt; i++)
    {
    xPos = i*hx + a;
    fprintf(dataFile, "%10.5e ",xPos);
    }
    fprintf(dataFile, "\n");
    fprintf(dataFile, "Y_COORDINATES %ld float \n",nPt);
    for(j = 0; j < nPt; j++)
    {
    yPos = j*hy + c;
    fprintf(dataFile, "%10.5e ",yPos);
    }
    fprintf(dataFile, "\n");

    fprintf(dataFile, "Z_COORDINATES %ld float \n",pPt);
    for(k = 0; k < pPt; k++)
    {
    zPos = k*hz + e;
    fprintf(dataFile, "%10.5e ",zPos);
    }
    fprintf(dataFile, "\n");

    fprintf(dataFile, "POINT_DATA %ld\n",dataCount);
    fprintf(dataFile, "SCALARS %s float\n",dataLabel.c_str());
    fprintf(dataFile, "LOOKUP_TABLE default\n");
    for(k = 0; k <  pPt; k++)
    {
    for(j = 0; j < nPt; j++)
    {
    for(i = 0; i < mPt; i++)
    {
    fprintf(dataFile, "%15.8e ",gridFun(i,j));
    }
    fprintf(dataFile, "\n");
    }}

    fclose(dataFile);

}


//
// Output data structure for ASCII and Binary output
//
// xPanels
// yPanels
// xMin
// xMax
// yMin
// yMax
// F(0,0)
// F(0,1)
// F(0,2)
//
//  ***
//


void outputToDataFile(const GridFunction2d& gF, const string& fileName, const string& formatString = "%20.15e")
{
//
//  Create format string
//
    ostringstream s;
    s.str("");
    s << formatString << " ";
//
//  Open and then write to a file
//
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    long xPanels = gF.getXpanelCount();
    long yPanels = gF.getYpanelCount();

    double xMin  = gF.getXmin();
    double yMin  = gF.getYmin();

    double xMax    = gF.getXmax();
    double yMax    = gF.getYmax();

    fprintf(dataFile,"%ld \n",    xPanels);
	fprintf(dataFile,"%ld \n",    yPanels);
    fprintf(dataFile,"%20.15e \n",xMin);
	fprintf(dataFile,"%20.15e \n",xMax);
	fprintf(dataFile,"%20.15e \n",yMin);
	fprintf(dataFile,"%20.15e \n",yMax);

    for(long i = 0; i <= xPanels; i++)
    {
    for(long j = 0; j <= yPanels; j++)
    {
    fprintf(dataFile, s.str().c_str(),gF(i,j));
    }
    fprintf(dataFile, "\n");
    }


    fclose(dataFile);
}

void inputFromDataFile(GridFunction2d& gF, FILE* dataFile, string fileName = "")
{
	size_t rValue = 0;

    long xPanels;
    long yPanels;

    double xMin;
    double yMin;

    double xMax;
    double yMax;

    rValue = fscanf(dataFile,"%ld", &xPanels) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%ld", &yPanels) != 1 ? 1 : rValue;

    rValue = fscanf(dataFile,"%lf",&xMin) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%lf",&xMax) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%lf",&yMin) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%lf",&yMax) != 1 ? 1 : rValue;

    rValue = (xPanels <= 0) ? 1 : rValue;
    rValue = (yPanels <= 0) ? 1 : rValue;
	rValue = (xMax < xMin) ? 1 : rValue;
    rValue = (yMax < yMin) ? 1 : rValue;


    gF.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax);

    for(long i = 0; i <= xPanels; i++)
    {
	for(long j = 0; j <= yPanels; j++)
    {
    rValue = fscanf(dataFile,"%lf",&gF(i,j)) != 1 ? 1 : rValue;
    }}

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction2d could not be initialized from file " + fileName + " \n");
    }

}

void inputFromDataFile(GridFunction2d& gF, const string& fileName)
{
//
//  Open input file
//
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "r" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    inputFromDataFile(gF, dataFile,fileName);

	fclose(dataFile);
}

void outputToBinaryDataFile(const GridFunction2d& gF, FILE* dataFile)
{
    long dataSize;

    long xPanels = gF.getXpanelCount();
    long yPanels = gF.getYpanelCount();

    std::int64_t xPanels64 = (std::int64_t) xPanels;
    std::int64_t yPanels64 = (std::int64_t) yPanels;

    double xMin  = gF.getXmin();
	double xMax  = gF.getXmax();
    double yMin  = gF.getYmin();
    double yMax  = gF.getYmax();

	//
	//  Write out the grid structure information. Using std:int64
	//  for integer values to avoid problems with machines with
	//  alternate storage sizes for int's and long's
	//

    fwrite(&xPanels64,  sizeof(std::int64_t), 1, dataFile);
	fwrite(&yPanels64,  sizeof(std::int64_t), 1, dataFile);
	fwrite(&xMin,  sizeof(double), 1, dataFile);
	fwrite(&xMax,  sizeof(double), 1, dataFile);
	fwrite(&yMin,  sizeof(double), 1, dataFile);
    fwrite(&yMax,  sizeof(double), 1, dataFile);

//
//  Write ot the function values
//
    dataSize = (xPanels+1)*(yPanels+1);
    fwrite(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}


void inputFromBinaryDataFile(GridFunction2d& gF, const string& fileName)
{

//  Open input file for binary read

    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "rb" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

	inputFromBinaryDataFile(gF,dataFile,fileName);

    fclose(dataFile);
}


void inputFromBinaryDataFile(GridFunction2d& gF, FILE* dataFile, string fileName = "")
{
	size_t rValue;
    size_t dataSize;

    long    xPanels;
    long    yPanels;
    double xMin; double yMin;
    double xMax; double yMax;

	std::int64_t xPanels64;
	std::int64_t yPanels64;


	rValue = 0;

	rValue = fread(&xPanels64,  sizeof(std::int64_t), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&yPanels64,  sizeof(std::int64_t), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&xMin,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&xMax,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&yMin,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&yMax,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;

	xPanels = (long)xPanels64;
	yPanels = (long)yPanels64;

	rValue = (xPanels <= 0) ? 1 : rValue;
    rValue = (yPanels <= 0) ? 1 : rValue;
    rValue = (xMax < xMin)  ? 1 : rValue;
    rValue = (yMax < yMin)  ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction2d could not be initialized from file " + fileName + " \n");
    }

    // Initialize instance and then read in the data

	gF.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax);
	dataSize = (xPanels+1)*(yPanels+1);

	rValue = fread(gF.getDataPointer(),  sizeof(double), dataSize, dataFile) != dataSize ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction2d could not be initialized from file " + fileName + " \n");
    }

}


void inputValuesFromBinaryDataFile(GridFunction2d& gF, FILE* dataFile, string fileName = "")
{
	size_t rValue = 0;
    size_t dataSize;

    long xPanels = gF.getXpanelCount();
    long yPanels = gF.getYpanelCount();

    rValue = (xPanels <= 0) ? 1 : rValue;
    rValue = (yPanels <= 0) ? 1 : rValue;

	dataSize = (xPanels+1)*(yPanels+1);
	rValue = fread(gF.getDataPointer(),  sizeof(double), dataSize, dataFile) != dataSize ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nValues from SCC::GridFunction2d could not be read from file " + fileName + " \n");
    }
}


void appendValuesToBinaryDataFile(const GridFunction2d& gF, FILE* dataFile)
{

	long dataSize;
    long xPanels = gF.getXpanelCount();
    long yPanels = gF.getYpanelCount();
//
//  Write ot the function values
//
    dataSize = (xPanels+1)*(yPanels+1);
    fwrite(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}

//
// This routine opens up a new file and write GridFunction2d structure and data
// and then closes the file
//
void outputToBinaryDataFile(const GridFunction2d& gF, const string& fileName)
{
//
//  Open and then write to a file (remember to use the b mode to specify binary!!!!)
//
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "wb" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    outputToBinaryDataFile(gF, dataFile);

    fclose(dataFile);
}

};

}

#undef OPENFILE
#ifdef _MSC_VER
#pragma warning(pop)
#endif
#endif




 
