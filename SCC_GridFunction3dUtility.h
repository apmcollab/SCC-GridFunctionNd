
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
 *  Release : 18.09.03
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

#ifndef _SCC_GridFunction3dUtility_
#define _SCC_GridFunction3dUtility_
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <functional>
#include <cassert>
#include <cstdint>
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


#include "../DoubleVectorNd/SCC_DoubleVector3d.h"
#include "../GridFunctionNd/SCC_GridFunction3d.h"


namespace SCC
{
class  GridFunction3dUtility
{
public :


void outputDataToVTKfile(const GridFunction3d& gridFun, const string& fileName, const string& dataLabel)
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

void outputDataToVTKfile(const GridFunction3d& gridFun, const string& fileName, const string& dataLabel,
double xScalingFactor)
{
    string scalingCoord = "x";
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
// This scales the coordinate values of the scaling coordinate specified (one of "x", "y" or "z") so that
// the coordinate extent is the scaling value times the largest size of the domain in the other coordinate
// directions. This utility is added so that rectangular regions in which the specified coordinate is thin
// with respect to the others can be viewed with the thin region expanded.
//
void outputDataToVTKfile(const GridFunction3d& gridFun, const string& fileName, const string& dataLabel,
string scalingCoord, double scalingValue)
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
    yScalingFactor = (yScalingFactor*transverseSizeMax)/(f-e);
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


void outputToDataFile(const GridFunction3d& gF, const string& fileName, const string& formatString = "%20.15e")
{
//
//  Create format string
//
    ostringstream s;
    s.str("");
    s << formatString << " ";

    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    long xPt = gF.getXpanelCount() + 1;
    long yPt = gF.getYpanelCount() + 1;
	long zPt = gF.getZpanelCount() + 1;

    double xMin  = gF.getXmin();
    double yMin  = gF.getYmin();
	double zMin  = gF.getZmin();

    double xMax    = gF.getXmax();
    double yMax    = gF.getYmax();
    double zMax    = gF.getZmax();

    fprintf(dataFile,"%ld \n", xPt);
	fprintf(dataFile,"%ld \n", yPt);
	fprintf(dataFile,"%ld \n", zPt);

    fprintf(dataFile,"%20.15e \n",xMin);
	fprintf(dataFile,"%20.15e \n",xMax);
	fprintf(dataFile,"%20.15e \n",yMin);
	fprintf(dataFile,"%20.15e \n",yMax);
    fprintf(dataFile,"%20.15e \n",zMin);
	fprintf(dataFile,"%20.15e \n",zMax);

    for(long i = 0; i < xPt; i++)
    {
    for(long j = 0; j < yPt; j++)
    {
    for(long k = 0; k < zPt; k++)
    {
    fprintf(dataFile, s.str().c_str(),gF(i,j,k));
    }
    fprintf(dataFile, "\n");
    }}


    fclose(dataFile);
}

void inputFromDataFile(GridFunction3d& gF, FILE* dataFile, string fileName = "")
{
	size_t rValue = 0;

    long xPt;
    long yPt;
	long zPt;

    double xMin;
    double yMin;
	double zMin;

    double xMax;
    double yMax;
    double zMax;

    rValue = fscanf(dataFile,"%ld", &xPt) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%ld", &yPt) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%ld", &zPt) != 1 ? 1 : rValue;

    rValue = (xPt <= 0) ? 1 : rValue;
    rValue = (yPt <= 0) ? 1 : rValue;
    rValue = (zPt <= 0) ? 1 : rValue;

    rValue = fscanf(dataFile,"%lf",&xMin) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%lf",&xMax) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%lf",&yMin) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%lf",&yMax) != 1 ? 1 : rValue;
    rValue = fscanf(dataFile,"%lf",&zMin) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%lf",&zMax) != 1 ? 1 : rValue;

	rValue = (xMax < xMin) ? 1 : rValue;
    rValue = (yMax < yMin) ? 1 : rValue;
    rValue = (zMax < zMin) ? 1 : rValue;

    gF.initialize(xPt-1,xMin,xMax,yPt-1,yMin,yMax,zPt-1,zMin,zMax);

    for(long i = 0; i < xPt; i++)
    {
    for(long j = 0; j < yPt; j++)
    {
	for(long k = 0; k < zPt; k++)
    {
    rValue = fscanf(dataFile,"%lf",&gF(i,j,k)) != 1 ? 1 : rValue;
    }
    }}

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction3d could not be initialized from file " + fileName + " \n");
    }

}

void inputFromDataFile(GridFunction3d& gF, const string& fileName)
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

//
//                    !!! Note !!!
//
// Date: Sept. 2, 2018
//
// The structure of data output changed to be consistent with ASCII output created by
// outputToDataFile(...). Binary data files written using older versions of this
// class cannot be read by this version.
//

void outputToBinaryDataFile(const GridFunction3d& gF, FILE* dataFile)
{
    long dataSize;

    std::int64_t Xpt64 = gF.getXpanelCount() + 1;
    std::int64_t Ypt64 = gF.getYpanelCount() + 1;
	std::int64_t Zpt64 = gF.getZpanelCount() + 1;

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

    fwrite(&Xpt64,  sizeof(std::int64_t), 1, dataFile);
	fwrite(&Ypt64,  sizeof(std::int64_t), 1, dataFile);
	fwrite(&Zpt64,  sizeof(std::int64_t), 1, dataFile);

	fwrite(&xMin,  sizeof(double), 1, dataFile);
	fwrite(&xMax,  sizeof(double), 1, dataFile);
	fwrite(&yMin,  sizeof(double), 1, dataFile);
    fwrite(&yMax,  sizeof(double), 1, dataFile);
	fwrite(&zMin,  sizeof(double), 1, dataFile);
	fwrite(&zMax,  sizeof(double), 1, dataFile);
//
//  Write ot the function values
//
    dataSize = Xpt64*Ypt64*Zpt64;
    fwrite(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}


void inputFromBinaryDataFile(GridFunction3d& gF, const string& fileName)
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

void inputFromBinaryDataFile(GridFunction3d& gF, const string& fileName, int& noFileFlag)
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




void inputFromBinaryDataFile(GridFunction3d& gF, FILE* dataFile, string fileName = "")
{
	size_t rValue;
    long dataSize;

    long    xPt;    long yPt;    long zPt;
    double xMin; double yMin; double zMin;
    double xMax; double yMax; double zMax;

	std::int64_t Xpt64;
	std::int64_t Ypt64;
	std::int64_t Zpt64;

	rValue = 0;

	rValue = fread(&Xpt64,  sizeof(std::int64_t), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&Ypt64,  sizeof(std::int64_t), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&Zpt64,  sizeof(std::int64_t), 1, dataFile) != 1 ? 1 : rValue;

	rValue = fread(&xMin,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&xMax,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;

	rValue = fread(&yMin,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&yMax,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;

	rValue = fread(&zMin,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&zMax,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;


    rValue = (xMax < xMin) ? 1 : rValue;
    rValue = (yMax < yMin) ? 1 : rValue;
    rValue = (zMax < zMin) ? 1 : rValue;


	xPt = (long)Xpt64;
	yPt = (long)Ypt64;
	zPt = (long)Zpt64;

	rValue = (xPt <= 0) ? 1 : rValue;
    rValue = (yPt <= 0) ? 1 : rValue;
    rValue = (zPt <= 0) ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction3d could not be initialized from file " + fileName + " \n");
    }

    // Initialize instance and then read in the data

	gF.initialize(xPt-1,xMin,xMax,yPt-1,yMin,yMax,zPt-1,zMin,zMax);
	dataSize = xPt*yPt*zPt;

	rValue = fread(gF.getDataPointer(),  sizeof(double), dataSize, dataFile) != (uint)dataSize ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction3d could not be initialized from file " + fileName + " \n");
    }

}


void inputValuesFromBinaryDataFile(GridFunction3d& gF, FILE* dataFile, string fileName = "")
{
	size_t rValue = 0;

    long dataSize;

    long xPt = gF.getXpanelCount() + 1;
    long yPt = gF.getYpanelCount() + 1;
	long zPt = gF.getZpanelCount() + 1;

    rValue = (xPt <= 0) ? 1 : rValue;
    rValue = (yPt <= 0) ? 1 : rValue;
    rValue = (zPt <= 0) ? 1 : rValue;

	dataSize = xPt*yPt*zPt;
	rValue = fread(gF.getDataPointer(),  sizeof(double), dataSize, dataFile) != (uint)dataSize ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nValues from SCC::GridFunction3d could not be read from file " + fileName + " \n");
    }
}


void appendValuesToBinaryDataFile(const GridFunction3d& gF, FILE* dataFile)
{

	long dataSize;
    long xPt = gF.getXpanelCount() + 1;
    long yPt = gF.getYpanelCount() + 1;
	long zPt = gF.getZpanelCount() + 1;
//
//  Write ot the function values
//
    dataSize = xPt*yPt*zPt;
    fwrite(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}

//
// This routine opens up a new file and write GridFunction3d structure and data
// and then closes the file
//
void outputToBinaryDataFile(const GridFunction3d& gF, const string& fileName)
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

/*

void GridFunction3dUtility::getDomainData(XML_ParameterListArray& paramList,
                      long& xPanels, double& xMin, double& xMax,
                      long& yPanels, double& yMin, double& yMax,
                      long& zPanels, double& zMin, double& zMax)
{
	if(paramList.isParameterList("ComputationalDomain") == 0)
    {
    string mesg = "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
    mesg       += "Error initializing GridFunction3d class \n";
    mesg       += "ComputationalDomain parameter list was not found \n";
    mesg       += "in input XML_ParameterListArray \n";
    mesg       += "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
    throw runtime_error(mesg);
    }

    if(paramList.isParameterList("GridParameters") == 0)
    {
    string mesg = "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
    mesg       += "Error initializing GridFunction3d class \n";
    mesg       += "GridParameters parameter list was not found \n";
    mesg       += "in input XML_ParameterListArray \n";
    mesg       += "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
    throw runtime_error(mesg);
    }

    xPanels  = paramList.getParameterValue("xPanels","GridParameters");
    yPanels  = paramList.getParameterValue("yPanels","GridParameters");
    zPanels  = paramList.getParameterValue("zPanels","GridParameters");

    xMin  = paramList.getParameterValue("xMin","ComputationalDomain");
	yMin  = paramList.getParameterValue("yMin","ComputationalDomain");
	zMin  = paramList.getParameterValue("zMin","ComputationalDomain");
	xMax  = paramList.getParameterValue("xMax","ComputationalDomain");
	yMax  = paramList.getParameterValue("yMax","ComputationalDomain");
	zMax  = paramList.getParameterValue("zMax","ComputationalDomain");
}

*/

};
}

#undef OPENFILE
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif
 
