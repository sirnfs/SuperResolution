// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently

#pragma once
#include <stdlib.h>
#include <iostream>
#include <crtdbg.h>
#include <vector>
#include <tchar.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <bitset>
#include <limits>
#include <time.h>
#include <fstream>
//ITK headers
#include "itkOffset.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVTKImageIO.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkVector.h"
#include "itkPointSet.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkPadImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkZeroFluxNeumannPadImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkGaussianImageSource.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
//MATLAB header
#include "engine.h"
// OpenCV Headers
//#include <opencv2/imgproc/imgproc.hpp>  
//#include <opencv2/imgproc/imgproc_c.h> 
//#include <opencv2/core/core.hpp>      
//#include <opencv2/core/operations.hpp>    
//#include <opencv2/highgui/highgui.hpp>  
//#include <opencv2/legacy/compat.hpp>

#define PI 3.1415926535897932384626433