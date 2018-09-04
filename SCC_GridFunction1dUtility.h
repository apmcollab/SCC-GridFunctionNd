/*
 * SCC_GridFunction1dUtility.h
 *
 *  Created on: Jun 28, 2015
 *      Author: anderson
 *
 *
 * Utility routines - mostly routines for file input and output
 * of GridFunction1d instances.
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

#ifndef _SCC_GridFunction1dUtility_
#define _SCC_GridFunction1dUtility_


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>
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

#include "SCC_GridFunction1d.h"


namespace SCC
{

class GridFunction1dUtility
{
public:

void outputToGNUplot(const GridFunction1d& gF, const string& fileName,
const string& formatString = "%20.15e")
{
//
//  Open and then write to a file
//
    ostringstream s;

    s.str("");
    s << formatString << "  " << formatString << " \n";

    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    long xPanels   = gF.getXpanelCount();
    double  xMin   = gF.getXmin();
    double  xMax   = gF.getXmax();
    double  hx     = (xMax-xMin)/(double)xPanels;

    fprintf(dataFile,"# %ld %20.15e %20.15e  \n", xPanels,xMin,xMax);

    double x;
    for(long i = 0;  i <= xPanels; i++)
    {
    x = xMin + i*hx;
    fprintf(dataFile,(s.str()).c_str(),x, gF.Values(i));
    }

    fclose(dataFile);
}

void outputToMatlab(const GridFunction1d& gF,  const string& fileName, const string& formatString = "%20.15e")
{
//
//  Open and then write to a file
//
    ostringstream s;
    s.str("");
    s << formatString << "  " << formatString << " \n";

    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }
    long i;

    double xp;

    double hx     = gF.getHx();
    long mPanel   = gF.getXpanelCount();
    double xMin   = gF.getXmin();

	fprintf(dataFile,"%ld \n",mPanel+1);

	// print out x values

    for(i = 0;  i <= mPanel; i++)
    {
    xp = xMin + i*hx;
    fprintf(dataFile,formatString.c_str(),xp);
	fprintf(dataFile,"\n");
    }

	// print out y values

	for(i = 0;  i <= mPanel; i++)
    {
    fprintf(dataFile,formatString.c_str(),gF.Values(i));
	fprintf(dataFile,"\n");
    }

    fclose(dataFile);
}

void outputFunction(const GridFunction1d& fun, const string& fileName, const string& outputFormat,
const string& formatString = "%20.15e")
{
    ostringstream s;
    s.str("");

    if(outputFormat.compare("GNUPLOT")==0)
    {
    s << fileName << ".dat";
    outputToGNUplot(fun, (s.str()).c_str(),formatString);
    }
    else if(outputFormat.compare("MATLAB")==0)
    {
    s << fileName << ".dat";
    outputToMatlab(fun, (s.str()).c_str(),formatString);
    }
	else if(outputFormat.compare("MATLAB_PLAIN")==0)
    {
    s << fileName << ".dat";
    outputToMatlab(fun, (s.str()).c_str(),formatString);
    }
    else
    {
    s << fileName << ".dat";
    outputToGNUplot(fun, (s.str()).c_str(),formatString);
    }
}


void appendToGNUplot(const GridFunction1d& gF,const string& fileName, const string& formatString = "%20.15e")
{
//
//  Open and then write to a file
//
	ostringstream s;

    s.str("");
    s << formatString << "  " << formatString << " \n";


    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "a+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    long i;

    double* xp;

    double hx     = gF.getHx();
    long mPanel   = gF.getXpanelCount();
    double xMin   = gF.getXmin();

    xp     = new double[mPanel+1];

    for(i = 0; i <= mPanel; i++)
    {
    xp[i] = xMin + i*hx;
    }

    fprintf(dataFile,"\n");
    for(i = 0;  i <= mPanel; i++)
    {
    fprintf(dataFile,(s.str()).c_str(),xp[i], gF.Values(i));
    }

    delete [] xp;

    fclose(dataFile);
}

void inputFromGNUplot(GridFunction1d& gF, const string& fileName)
{
	//  Open input file

    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "r" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    size_t rValue = 0;

    char poundChar;
    long xPanels;
    double xMin; double xMax;

    rValue = fscanf(dataFile,"%c", &poundChar)  != 1 ? 1 : rValue;   // remove leading #
    rValue = fscanf(dataFile,"%ld",&xPanels)    != 1 ? 1 : rValue;
    rValue = fscanf(dataFile,"%lf",&xMin)       != 1 ? 1 : rValue;
    rValue = fscanf(dataFile,"%lf",&xMax)       != 1 ? 1 : rValue;

    rValue = (xMax < xMin)  ? 1 : rValue;
    rValue = (xPanels <= 0) ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction1d could not be initialized from gnuplot data file " + fileName + " \n");
    }


    gF.initialize(xPanels,xMin,xMax);

    double x;
    for(long i = 0;  i <= xPanels; i++)
    {
	rValue = fscanf(dataFile,"%lf %lf ",&x,&gF(i)) != 2 ? 1 : rValue;
    }

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction1d could not be initialized from gnuplot data file " + fileName + " \n");
    }

	fclose(dataFile);
}

//
// Outputs the data to a file in gnuplot format, and then creates a very basic Veusz plot file that
// is linked to that data.
//
void outputToVeusz(const GridFunction1d& gF, const string& fileName, const string& formatString = "%20.15e")
{
//	 Output data file in gnuplot format

	outputToGNUplot(gF, fileName,formatString);

	int lastindex    = fileName.find_last_of(".");
    string baseName  = fileName.substr(0, lastindex);
    string veuszName = baseName;
    veuszName.append(".vsz");

    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

	ostringstream s;
    s.str("");

	fprintf(dataFile,"%s","#Veusz saved document (version 1.23.1)\n");
	fprintf(dataFile,"%s","AddImportPath(u'./')\n");

	s << "ImportFile(u'./" << fileName << "', u'', ignoretext=True, linked=True, prefix=u'" << baseName << "')\n";
	fprintf(dataFile,"%s",(s.str()).c_str());

	fprintf(dataFile,"%s","Add('page', name='page1', autoadd=False)\n");
	fprintf(dataFile,"%s","To('page1')\n");
	fprintf(dataFile,"%s","Add('graph', name='graph1', autoadd=False)\n");
	fprintf(dataFile,"%s","To('graph1')\n");
	fprintf(dataFile,"%s","Add('axis', name='x', autoadd=False)\n");
	fprintf(dataFile,"%s","Add('axis', name='y', autoadd=False)\n");
	fprintf(dataFile,"%s","To('y')\n");
	fprintf(dataFile,"%s","Set('direction', 'vertical')\n");
	fprintf(dataFile,"%s","To('..')\n");
	fprintf(dataFile,"%s","Add('xy', name='xy1', autoadd=False)\n");
	fprintf(dataFile,"%s","To('xy1')\n");

	s.str("");
	s << "Set('xData', u'" << baseName << "1" << "')\n";
	fprintf(dataFile,"%s",(s.str()).c_str());
	s.str("");
	s << "Set('yData', u'" << baseName << "2" << "')\n";
	fprintf(dataFile,"%s",(s.str()).c_str());

	fprintf(dataFile,"%s","Set('markerSize', u'1.5pt')\n");
	fprintf(dataFile,"%s","To('..')\n");
    fprintf(dataFile,"%s","To('..')\n");
    fprintf(dataFile,"%s","To('..')\n");

    fclose(dataFile);
}




//
// Output data structure for ASCII and Binary output
//
// xPanels
// xMin
// xMax
// F(0)
// F(1)
// F(2)
//
//  ***
//


void outputToDataFile(const GridFunction1d& gF, const string& fileName, const string& formatString = "%20.15e")
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

    double xMin  = gF.getXmin();
    double xMax  = gF.getXmax();

    fprintf(dataFile,"%ld \n", xPanels);

    fprintf(dataFile,"%20.15e \n",xMin);
	fprintf(dataFile,"%20.15e \n",xMax);


    for(long i = 0; i <= xPanels; i++)
    {
    fprintf(dataFile, s.str().c_str(),gF(i));
    }


    fclose(dataFile);
}

void inputFromDataFile(GridFunction1d& gF, FILE* dataFile, string fileName = "")
{
	size_t rValue = 0;

    long xPanels;
    double xMin; double xMax;

    rValue = fscanf(dataFile,"%ld", &xPanels) != 1 ? 1 : rValue;

    rValue = (xPanels <= 0) ? 1 : rValue;

    rValue = fscanf(dataFile,"%lf",&xMin) != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%lf",&xMax) != 1 ? 1 : rValue;

	rValue = (xMax < xMin) ? 1 : rValue;


    gF.initialize(xPanels,xMin,xMax);

    for(long i = 0; i <= xPanels; i++)
    {
    rValue = fscanf(dataFile,"%lf",&gF(i)) != 1 ? 1 : rValue;
    }

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction1d could not be initialized from file " + fileName + " \n");
    }

}

void inputFromDataFile(GridFunction1d& gF, const string& fileName)
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

void outputToBinaryDataFile(const GridFunction1d& gF, FILE* dataFile)
{
    long dataSize;

    long         xPanels   = gF.getXpanelCount();
    std::int64_t xPanels64 = (std::int64_t)xPanels;

    double xMin  = gF.getXmin();
	double xMax  = gF.getXmax();

	//
	//  Write out the grid structure information. Using std:int64
	//  for integer values to avoid problems with machines with
	//  alternate storage sizes for int's and long's
	//
    fwrite(&xPanels64,  sizeof(std::int64_t), 1, dataFile);
	fwrite(&xMin,  sizeof(double), 1, dataFile);
	fwrite(&xMax,  sizeof(double), 1, dataFile);

//
//  Write the function values
//
    dataSize = xPanels+1;
    fwrite(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}


void inputFromBinaryDataFile(GridFunction1d& gF, const string& fileName)
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



void inputFromBinaryDataFile(GridFunction1d& gF, FILE* dataFile, string fileName = "")
{
	size_t rValue;
    long dataSize;

    long    xPanels;
    double xMin;
    double xMax;

	std::int64_t xPanels64;

	rValue = 0;

	rValue  = fread(&xPanels64,  sizeof(std::int64_t), 1, dataFile) != 1 ? 1 : rValue;
	xPanels = (long)xPanels64;

	rValue = fread(&xMin,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&xMax,  sizeof(double), 1, dataFile) != 1 ? 1 : rValue;

    rValue = (xMax < xMin) ? 1 : rValue;
	rValue = (xPanels <= 0) ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction1d could not be initialized from file " + fileName + " \n");
    }

    // Initialize instance and then read in the data

	gF.initialize(xPanels,xMin,xMax);

	dataSize = xPanels+1;

	rValue = fread(gF.getDataPointer(),  sizeof(double), dataSize, dataFile) != (uint)dataSize ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nSCC::GridFunction1d could not be initialized from file " + fileName + " \n");
    }

}


void inputValuesFromBinaryDataFile(GridFunction1d& gF, FILE* dataFile, string fileName = "")
{
	size_t rValue = 0;

    long dataSize;

    long xPanels = gF.getXpanelCount();

    rValue = (xPanels <= 0) ? 1 : rValue;

	dataSize = xPanels+1;
	rValue = fread(gF.getDataPointer(),  sizeof(double), dataSize, dataFile) != (uint)dataSize ? 1 : rValue;

    if(rValue == 1)
    {
    throw std::runtime_error("\nValues from SCC::GridFunction1d could not be read from file " + fileName + " \n");
    }
}


void appendValuesToBinaryDataFile(const GridFunction1d& gF, FILE* dataFile)
{
	long dataSize;
    long xPanels = gF.getXpanelCount();

//  Write ot the function values

    dataSize = xPanels+1;
    fwrite(gF.getDataPointer(),  sizeof(double), dataSize, dataFile);
}

//
// This routine opens up a new file and write GridFunction1d structure and data
// and then closes the file
//
void outputToBinaryDataFile(const GridFunction1d& gF, const string& fileName)
{
//
//  Open and then write to a file (remember to use the b mode to specify binary!!!!)
//
    FILE* dataFile;

    if(OPENFILE(dataFile,fileName.c_str(), "wb" ))
    {
      printf( "The file %s could not be  opened\n",fileName.c_str());
      return;
    }
    outputToBinaryDataFile(gF, dataFile);

    fclose(dataFile);
}

//
// Veusz plot used to determine veusz plot data structure for the above member function
//
/*
# Veusz saved document (version 1.23.1)
# Saved at 2015-06-29T22:42:43.183012

AddImportPath(u'./')
ImportFile(u'./Mollifier.dat', u'', ignoretext=True, linked=True, prefix=u'M')
Add('page', name='page1', autoadd=False)
To('page1')
Add('graph', name='graph1', autoadd=False)
To('graph1')
Add('axis', name='x', autoadd=False)
Add('axis', name='y', autoadd=False)
To('y')
Set('direction', 'vertical')
To('..')
Add('xy', name='xy1', autoadd=False)
To('xy1')
Set('xData', u'M1')
Set('yData', u'M2')
Set('markerSize', u'1.5pt')
To('..')
To('..')
To('..')
 */






};
} // namespace


#undef  OPENFILE
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif /* SCC_GRIDFUNCTION_SCC_GRIDFUNCTION1DUTILITY_H_ */
