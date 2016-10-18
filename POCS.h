#ifndef __POCS__
#define __POCS__

#include "stdafx.h"

class POCS
{
  public:	
	//Dimension of work area
	static const unsigned int Dimension = 2;  //This '3' needs to be changed if we want to work in other dimensions -- this needs to be set some other way through main
	//Input image types and pointers
    typedef unsigned short PixelType;
	typedef double PixelType_Double;
    typedef itk::Image<PixelType, Dimension> ImageType; 
	//For MV iterator
    typedef unsigned int ArrayType; 
    typedef itk::Image<ArrayType, Dimension> Vector_Array; 
    typedef itk::ImageRegionIterator<Vector_Array> MVIterator;
	//For validity -- array and iterator
	typedef itk::Image<ArrayType, Dimension> Validity; 
	typedef itk::ImageRegionIterator<Validity> ValIterator;

	//--------------Public functions------------------
    POCS(char *anchorframe_file, int interp_factor, double d_iter, double variance, int blur_sz, int blur_scale, int v_threshold, int m_intensity);
	void process_frame(ImageType::Pointer c_frame, std::vector<MVIterator>& currentMV_iterator);
	void process_frame_reverse(ImageType::Pointer c_frame, std::vector<MVIterator>& currentMV_iterator, Validity::Pointer valid1);
	void write_output(char *output_filename);
	~POCS();

	//--------------Public variables------------------
	
  private:
	//--------------Private variables-----------------
	//For diagnostic output purposes
	std::ofstream fout; 
	//Double Image Type for HR image during POCS
	typedef itk::Image<PixelType_Double, Dimension> ImageType_Double; 
    ImageType::Pointer HR_image;
	ImageType_Double::Pointer HR_image_double;
	//This is the interpolation factor that will be used for HR image
	int i_factor;
	//Error allowed in POCS
	double delta_iter;
	//The sum of the squares of all the kernel coefficients
	double h_squared;
	//Blur kernel
	int blur_size;
	typedef double blur_dtype;
	typedef itk::Image<blur_dtype, Dimension> blur_image;
	typedef itk::GaussianImageSource<blur_image> Gaussian_Generator;
	Gaussian_Generator::Pointer kernel;
	//Validity threshold -- the MVs whose values are greater than this threshold will not be used in the POCS reconstruction
	int validity_threshold;
	//Max intensity value -- for example, 8-bit max intensity is 255, 16-bit is 65535;
	int max_intensity;
	//Iterators
	//Create iterators to go through pixels in the current frame and in the HR frame, and through the MVs
    typedef itk::ImageRegionIterator<ImageType_Double> ImageRegionIterator_Double;
    typedef itk::ImageRegionIterator<ImageType> ImageRegionIterator;
	//For the kernel iterator
    typedef itk::ImageRegionIterator<blur_image> KernelRegionIterator;
	//Filter for interpolating image
	typedef itk::ResampleImageFilter<ImageType, ImageType> FilterType;
	//Interpolation function for filter
	typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
	//Interpolator function for filter -- using bspline interpolation
	typedef itk::BSplineInterpolateImageFunction<ImageType, double, double> BSpline_Interpolator;
	//Use identity transform when doing interpolation -- we don't want to do anything BUT interpolate
	typedef itk::IdentityTransform<double, Dimension> TransformType;
	//Reader for input file
    typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader;
	//Write for output file
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer;

	//--------------Private functions-----------------
	void calculate_gaussian_coefficients(double variance, int blur_scale);
	void increment_HRindex(itk::ImageRegionIterator<POCS::ImageType>& HR_index, int shift_amount, itk::Size<POCS::Dimension> shift_size, int& done, itk::Size<POCS::Dimension> boundary);
	void increment_HRindex2(itk::ImageRegionIterator<POCS::ImageType_Double>& HR_index, int shift_amount, itk::Size<POCS::Dimension> shift_size, int& done, itk::Size<POCS::Dimension> boundary);
	void increment_VALindex(itk::ImageRegionIterator<POCS::Validity>& HR_index, int shift_amount, itk::Size<POCS::Dimension> shift_size, int& done, itk::Size<POCS::Dimension> boundary);
	double min(double x, double y);
	double max(double x, double y);
};

#endif