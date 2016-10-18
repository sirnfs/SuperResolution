#ifndef __FRC__
#define __FRC__

#include "stdafx.h"

class FRC
{  
  public:
	//--------------Public variables------------------
	//Dimension of work area
	static const unsigned int Dimension = 2;  //This '3' needs to be changed if we want to work in other dimensions -- this needs to be set some other way through main
    //Input image types and pointers
    typedef unsigned short PixelType;
	typedef itk::Image<PixelType, Dimension> ImageType; 
	//Output image pointer
	ImageType::Pointer outputImage;
	//For MV iterator
    typedef unsigned int ArrayType; 
    typedef itk::Image<ArrayType, Dimension> Vector_Array; 
    typedef itk::ImageRegionIterator<Vector_Array> MVIterator;
	//For validity array
	typedef itk::Image<ArrayType, Dimension> Validity; 
	typedef itk::ImageRegionIterator<Validity> ValIterator;
		
	//-------------Public functions-------------------
	//Pass names of two images and the dimensions of the images
	FRC(int validity_thres);
    void perform_FRC(std::vector<MVIterator>& currentMV_iterator, ImageType::Pointer p_frame, ImageType::Pointer c_frame, Validity::Pointer valid1, char *output_filename);
	void frame_difference(char *in_filename, char *out_filename);
	~FRC();
	
  private:
	//--------------Private variables-----------------
	//Input image -- used for comparing with estimated FRC image
	ImageType::Pointer inputImage;
	ImageType::Pointer FRC_diff;
	int val_threshold;
	
				
	//-----------Private functions-------------------
	//Function used to find matching block in search region
	
};

#endif