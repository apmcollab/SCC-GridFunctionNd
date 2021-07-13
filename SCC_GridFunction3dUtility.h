
/*
 * SCC_GridFunction3dUtility.h
 *
 *  Created on: Jun 28, 2015
 *      Author: anderson
 *
 *
 * Utility routines - mostly routines for file input and output
 * of GridFunction3d instances.
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

#ifndef SCC_GRID_FUNCTION_3D_UTILITY_
#define SCC_GRID_FUNCTION_3D_UTILITY_

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <functional>
#include <cassert>
#include <cstdint>

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


#include "../DoubleVectorNd/SCC_DoubleVector3d.h"
#include "../GridFunctionNd/SCC_GridFunction3d.h"


namespace SCC
{
class  GridFunction3dUtility
{
public :


void outputDataToVTKfile(const GridFunction3d& gridFun, const std::string& fileName, const std::string& dataLabel)
{
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    double a  = gridFun.getXmin();
    double c  = gridFun.getYmin();
	double e  = gridFun.getZmin();

    double hx = gridFun.getHx();
    double hy = gridFun.getHy();
	double hz = gridFun.getHz();

    long mPt = gridFun.getXpanelCount() + 1;
    long nPt = gridFun.getYpanelCount() + 1;
	long pPt = gridFun.getZpanelCount() + 1;

    long dataCount = mPt*nPt*pPt;

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
    fprintf(dataFile, "%15.8e ",gridFun.Values(i,j,k));
    }
    fprintf(dataFile, "\n");
    }}

    fclose(dataFile);
}


// Legacy interface for scaling w.r.t. x coordinate only. This is a deprecated interface
// and the calling code should be updated.

void outputDataToVTKfile(const GridFunction3d& gridFun, const std::string& fileName, const std::string& dataLabel,
double xScalingFactor)
{
    std::string scalingCoord = "x";
    if(xScalingFactor > 0)
	{
	outputDataToVTKfile(gridFun,fileName,dataLabel,scalingCoord,xScalingFactor);
	}
	else
	{
	outputDataToVTKfile(gridFun,fileName,dataLabel);
	}

}
//
// This routine scales the coordinate values of the scaling coordinate specified (one of "x", "y" or "z") so that
// the coordinate extent is the scaling value times the largest size of the domain in the other coordinate
// directions. This utility is added so that rectangular regions in which the specified coordinate is thin
// with respect to the others can be viewed with the thin region expanded.
//
void outputDataToVTKfile(const GridFunction3d& gridFun, const std::string& fileName, const std::string& dataLabel,
std::string scalingCoord, double scalingValue)
{
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    double a  = gridFun.getXmin();  double b  = gridFun.getXmax();
    double c  = gridFun.getYmin();  double d  = gridFun.getYmax();
	double e  = gridFun.getZmin();  double f  = gridFun.getZmax();

    double hx = gridFun.getHx();
    double hy = gridFun.getHy();
	double hz = gridFun.getHz();

    double xScalingFactor = -1.0;
    double yScalingFactor = -1.0;
    double zScalingFactor = -1.0;

    if((scalingCoord == "x") || (scalingCoord == "X")  ){xScalingFactor = scalingValue;}
    if((scalingCoord == "y") || (scalingCoord == "Y")  ){yScalingFactor = scalingValue;}
    if((scalingCoord == "z") || (scalingCoord == "Z")  ){zScalingFactor = scalingValue;}

    double transverseSizeMax;

    if(xScalingFactor > 0.0)
    {
    transverseSizeMax = d-c;
    transverseSizeMax = (transverseSizeMax > (f-e)) ? transverseSizeMax : (f-e);
    xScalingFactor = (xScalingFactor*transverseSizeMax)/(b-a);
    }
    else {xScalingFactor = 1.0;}

    if(yScalingFactor > 0.0)
    {
    transverseSizeMax = b-a;
    transverseSizeMax = (transverseSizeMax > (f-e)) ? transverseSizeMax : (f-e);
    yScalingFactor = (yScalingFactor*transverseSizeMax)/(d-c);
    }
    else {yScalingFactor = 1.0;}


    if(zScalingFactor > 0.0)
    {
    transverseSizeMax = b-a;
    transverseSizeMax = (transverseSizeMax > (d-c)) ? transverseSizeMax : (d-c);
    zScalingFactor = (zScalingFactor*transverseSizeMax)/(f-e);
    }
    else {zScalingFactor = 1.0;}


    long mPt = gridFun.getXpanelCount() + 1;
    long nPt = gridFun.getYpanelCount() + 1;
	long pPt = gridFun.getZpanelCount() + 1;

    long dataCount = mPt*nPt*pPt;

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
    xPos = i*hx*xScalingFactor + a*xScalingFactor;
    fprintf(dataFile, "%10.5e ",xPos);
    }
    fprintf(dataFile, "\n");
    fprintf(dataFile, "Y_COORDINATES %ld float \n",nPt);
    for(j = 0; j < nPt; j++)
    {
    yPos = j*hy*yScalingFactor + c*yScalingFactor;
    fprintf(dataFile, "%10.5e ",yPos);
    }
    fprintf(dataFile, "\n");

    fprintf(dataFile, "Z_COORDINATES %ld float \n",pPt);
    for(k = 0; k < pPt; k++)
    {
    zPos = k*hz*zScalingFactor + e*zScalingFactor;
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
    fprintf(dataFile, "%15.8e ",gridFun.Values(i,j,k));
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
// zPanels
// xMin
// xMax
// yMin
// yMax
// zMin
// zMax
// F(0,0,0)
// F(0,0,1)
// F(0,0,2)
//
//  ***
//



void outputToDataFile(const GridFunction3d& gF, const std::string& fileName, const std::string& formatString = "%20.15e")
{
//
//  Create format string
//
    std::ostringstream s;
    s.str("");
    s << formatString << " ";

    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    long xPanels = gF.getXpanelCount();
    long yPanels = gF.getYpanelCount();
    long zPanels = gF.getZpanelCount();

    double xMin  = gF.getXmin();
    double yMin  = gF.getYmin();
    double zMin  = gF.getZmin();

    double xMax    = gF.getXmax();
    double yMax    = gF.getYmax();
    double zMax    = gF.getZmax();

    fprintf(dataFile,"%ld \n",    xPanels);
	fprintf(dataFile,"%ld \n",    yPanels);
	fprintf(dataFile,"%ld \n",    zPanels);
    fprintf(dataFile,"%20.15e \n",xMin);
	fprintf(dataFile,"%20.15e \n",xMax);
	fprintf(dataFile,"%20.15e \n",yMin);
	fprintf(dataFile,"%20.15e \n",yMax);
    fprintf(dataFile,"%20.15e \n",zMin);
	fprintf(dataFile,"%20.15e \n",zMax);

    for(long i = 0; i <= xPanels; i++)
    {
    for(long j = 0; j <= yPanels; j++)
    {
    for(long k = 0; k <= zPanels; k++)
    {
    fprintf(dataFile, s.str().c_str(),gF(i,j,k));
    }
    fprintf(dataFile, "\n");
    }}


    fclose(dataFile);
}

void inputFromDataFile(GridFunction3d& gF, FILE* dataFile, std::string fileName = "")
{
	bool errFlag = false;

    long xPanels; long yPanels; long zPanels;

    double xMin; double yMin; double zMin;
    double xMax; double yMax; double zMax;

    errFlag = (fscanf(dataFile,"%ld", &xPanels) != 1  ) || errFlag;
	errFlag = (fscanf(dataFile,"%ld", &yPanels) != 1  ) || errFlag;
	errFlag = (fscanf(dataFile,"%ld", &zPanels) != 1  ) || errFlag;
    errFlag = (fscanf(dataFile,"%lf", &xMin)    != 1  ) || errFlag;
	errFlag = (fscanf(dataFile,"%lf", &xMax)    != 1  ) || errFlag;
	errFlag = (fscanf(dataFile,"%lf", &yMin)    != 1  ) || errFlag;
	errFlag = (fscanf(dataFile,"%lf", &yMax)    != 1  ) || errFlag;
    errFlag = (fscanf(dataFile,"%lf", &zMin)    != 1  ) || errFlag;
	errFlag = (fscanf(dataFile,"%lf", &zMax)    != 1  ) || errFlag;

	errFlag = (xPanels <= 0 ) || errFlag;
    errFlag = (yPanels <= 0 ) || errFlag;
    errFlag = (zPanels <= 0 ) || errFlag;

	errFlag = (xMax < xMin ) || errFlag;
    errFlag = (yMax < yMin ) || errFlag;
    errFlag = (zMax < zMin ) || errFlag;

    gF.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);

    for(long i = 0; i <= xPanels; i++)
    {
    for(long j = 0; j <= yPanels; j++)
    {
	for(long k = 0; k <= zPanels; k++)
    {
    errFlag = (fscanf(dataFile,"%lf",&gF(i,j,k)) != 1  ) || errFlag;
    }
    }}

    if(errFlag)
    {
    throw std::runtime_error("\nSCC::GridFunction3d could not be initialized from file " + fileName + " \n");
    }

}

void inputFromDataFile(GridFunction3d& gF, const std::string& fileName)
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


void outputToBinaryDataFile(const GridFunction3d& gF, FILE* dataFile)
{
    long dataSize;

    long xPanels = gF.getXpanelCount();
    long yPanels = gF.getYpanelCount();
    long zPanels = gF.getZpanelCount();

    std::int64_t xPanels64 = (std::int64_t) xPanels;
    std::int64_t yPanels64 = (std::int64_t) yPanels;
    std::int64_t zPanels64 = (std::int64_t) zPanels;

    double xMin  = gF.getXmin();
	double xMax  = gF.getXmax();
    double yMin  = gF.getYmin();
    double yMax  = gF.getYmax();
	double zMin  = gF.getZmin();
	double zMax  = gF.getZmax();

	//
	//  Write out the grid structure information. Using std:int64
	//  for integer values to avoid problems with machines with
	//  alternate storage sizes for int's and long's
	//

    fwrite(&xPanels64,  sizeof(std::int64_t), 1, dataFile);
	fwrite(&yPanels64,  sizeof(std::int64_t), 1, dataFile);
	fwrite(&zPanels64,  sizeof(std::int64_t), 1, dataFile);

	fwrite(&xMin,  sizeof(double), 1, dataFile);
	fwrite(&xMax,  sizeof(double), 1, dataFile);
	fwrite(&yMin,  sizeof(double), 1, dataFile);
    fwrite(&yMax,  sizeof(double), 1, dataFile);
	fwrite(&zMin,  sizeof(double), 1, dataFile);
	fwrite(&zMax,  sizeof(double), 1, dataFile);
//
//  Write ot the function values
//
    dataSize = (xPanels+1)*(yPanels+1)*(zPanels+1);
    fwrite(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}


void inputFromBinaryDataFile(GridFunction3d& gF, const std::string& fileName)
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

// Legacy call not using throw/catch error for file errors

void inputFromBinaryDataFile(GridFunction3d& gF, const std::string& fileName, int& noFileFlag)
{
	//
	//  Open input file (remember to use the b mode to specify binary!!!!)
	//
	FILE* dataFile = 0;
	noFileFlag     = 0;
	if(OPENFILE(dataFile,fileName.c_str(), "rb" ))
	{
		  noFileFlag = 1;
	      return;
	}

	inputFromBinaryDataFile(gF,dataFile,fileName);

    fclose(dataFile);
}


void inputFromBinaryDataFile(GridFunction3d& gF, FILE* dataFile, std::string fileName = "")
{
    size_t dataSize;

    long xPanels; long yPanels; long zPanels;

    double xMin; double yMin; double zMin;
    double xMax; double yMax; double zMax;

	std::int64_t xPanels64;
	std::int64_t yPanels64;
	std::int64_t zPanels64;

	bool errFlag = false;

	errFlag = (fread(&xPanels64,  sizeof(std::int64_t), 1, dataFile) != 1 ) || errFlag;
	errFlag = (fread(&yPanels64,  sizeof(std::int64_t), 1, dataFile) != 1 ) || errFlag;
	errFlag = (fread(&zPanels64,  sizeof(std::int64_t), 1, dataFile) != 1 ) || errFlag;

	errFlag = (fread(&xMin,  sizeof(double), 1, dataFile) != 1 ) || errFlag;
	errFlag = (fread(&xMax,  sizeof(double), 1, dataFile) != 1 ) || errFlag;

	errFlag = (fread(&yMin,  sizeof(double), 1, dataFile) != 1 ) || errFlag;
	errFlag = (fread(&yMax,  sizeof(double), 1, dataFile) != 1 ) || errFlag;

	errFlag = (fread(&zMin,  sizeof(double), 1, dataFile) != 1 ) || errFlag;
	errFlag = (fread(&zMax,  sizeof(double), 1, dataFile) != 1 ) || errFlag;

	xPanels = (long)xPanels64;
	yPanels = (long)yPanels64;
	zPanels = (long)zPanels64;

    errFlag = (xMax < xMin) || errFlag;
    errFlag = (yMax < yMin) || errFlag;
    errFlag = (zMax < zMin) || errFlag;

	errFlag = (xPanels <= 0) || errFlag;
    errFlag = (yPanels <= 0) || errFlag;
    errFlag = (zPanels <= 0) || errFlag;

    if(errFlag)
    {
    throw std::runtime_error("\nSCC::GridFunction3d could not be initialized from file " + fileName + " \n");
    }

    // Initialize instance and then read in the data

	gF.initialize(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels,zMin,zMax);

	dataSize = (xPanels+1)*(yPanels+1)*(zPanels+1);

	errFlag = (fread(gF.getDataPointer(),  sizeof(double), dataSize, dataFile) != dataSize) || errFlag;

    if(errFlag)
    {
    throw std::runtime_error("\nSCC::GridFunction3d could not be initialized from file " + fileName + " \n");
    }
}


void inputValuesFromBinaryDataFile(GridFunction3d& gF, FILE* dataFile, std::string fileName = "")
{
	bool errFlag = false;
    size_t dataSize;

    long xPanels = gF.getXpanelCount();
    long yPanels = gF.getYpanelCount();
	long zPanels = gF.getZpanelCount();

    errFlag = (xPanels <= 0)  || errFlag;
    errFlag = (yPanels <= 0)  || errFlag;
    errFlag = (zPanels <= 0)  || errFlag;

	dataSize = (xPanels+1)*(yPanels+1)*(zPanels+1);

	errFlag = (fread(gF.getDataPointer(),  sizeof(double), dataSize, dataFile) != dataSize  ) || errFlag;

    if(errFlag)
    {
    throw std::runtime_error("\nValues from SCC::GridFunction3d could not be read from file " + fileName + " \n");
    }
}


void appendValuesToBinaryDataFile(const GridFunction3d& gF, FILE* dataFile)
{

	long dataSize;
    long xPanels = gF.getXpanelCount();
    long yPanels = gF.getYpanelCount();
	long zPanels = gF.getZpanelCount();
//
//  Write ot the function values
//
    dataSize = (xPanels+1)*(yPanels+1)*(zPanels+1);
    fwrite(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}

//
// This routine opens up a new file and write GridFunction3d structure and data
// and then closes the file
//
void outputToBinaryDataFile(const GridFunction3d& gF, const std::string& fileName)
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
 
