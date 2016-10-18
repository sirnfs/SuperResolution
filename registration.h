#ifndef __Registration__
#define __Registration__

#include "stdafx.h"

//Dimension of work area -- global variable
static const unsigned int Dimension = 2;  //This '3' needs to be changed if we want to work in other dimensions -- this needs to be set some other way through main

class MVs_PL //this class will be used for motion vectors at each level of the pyramid (PL = pyramid level)
{
  friend class Registration; //So that class Registration can access the member variables below

  public:
	//--------------Public variables------------------

    //For motion vector array
    typedef unsigned int ArrayType; 
    typedef itk::Image<ArrayType, Dimension> Vector_Array; 
    std::vector<Vector_Array::Pointer> MV_array; //This is the motion vector array

	//MV Iterator
    typedef itk::ImageRegionIterator<Vector_Array> MVIterator;
    std::vector<MVIterator> MVs; 

    //MV Neighorhood Iterator
	typedef itk::ConstantBoundaryCondition<Vector_Array> BoundaryConditionType; //boundaries are assumed to be 0
    typedef itk::NeighborhoodIterator<Vector_Array, BoundaryConditionType> NeighborhoodIteratorType; //uses zero flux boundary condition with neighborhood iterator
	std::vector<NeighborhoodIteratorType> Reg_MVs;
    NeighborhoodIteratorType::RadiusType radius; //will be used to set the size of the radius for the neighborhood iterator

};

class Registration
{  
  public:
	//--------------Public variables------------------

	//Pixel and image type
    typedef unsigned short PixelType;
    typedef itk::Image<PixelType, Dimension> ImageType; 

	//Current frame -- this is saved for POCS so we dont have to keep interpolating the frames
	ImageType::Pointer current_frame; //this will hold the interpolated version of the current frame
    //Previous frame -- this if for the first frame in the argument list -- will be passed to other functions
	ImageType::Pointer previous_frame;

    //Create a vector of classes of type MVs_PL -- these will hold the motion vectors
	std::vector<MVs_PL*> MV_PLNum; 	

	//Create a vector of classes of type MVs_PL again -- these will hold the scaled motion vectors for de-aliasing algorithm
	std::vector<MVs_PL*> MV_alias;

	//For assigning a validity to each motion vector
	typedef itk::Image<unsigned int, Dimension> Validity;  
	Validity::Pointer valid;
	//Iterator for validity
	typedef itk::ImageRegionIterator<Validity> ValIterator;

	//-------------Public functions-------------------
	//Pass names of two images and the dimensions of the images
    Registration(char *image1_file, char *image2_file, int * bsize, int * ssize, int num_levels, int interp_factor, int f_num);
	void motion_estimation(); //perform block matching to find the best matching cube
	void scale_MVs(int scale_factor); //filling in MVs for a scaled image by using the MVs from smaller image 
	void volume_difference(); //take the absolute value of pixel difference between every pixel in two volumes
	void draw_MC_volume(); //draw the motion-compensated volume using motion vectors and pixels from previous frame
	void print_MVs(int increment); //print the list of motion vectors -- only used for debugging/testing purposes
	void print_MVs_alias(int increment); //print the list of motion vectors for the upscaled case -- only used for debugging/testing purposes
	void draw_MVs_MATLAB(int increment); //draw the MVs using Matlab
	void write_output(char *output_filename); //write/draw the output volume
	//template <class T> T min(const T& x, const T& y); //find minimum between two variables 
	//template <class T> T max(const T& x, const T& y); //find maximum between two variables
	~Registration();
	
  private:
	//--------------Private variables-----------------

	//Input images
	ImageType::Pointer inputImage1;
	ImageType::Pointer inputImage2;	

	//For diagnostic output purposes
	std::ofstream fout; 

    //# of pyramid levels, block and search size
	int levels;
	int *b_size; //this will hold the block size array -- MUST FREE THIS!
	int *s_size; //this will hold the search size array -- MUST FREE THIS!
	//Keeps track of which level in the image pyramid we're at
	int current_level;
	//This is the interpolation factor that will be used for getting sub-pixel motion vectors
	int i_factor;	
	//This is the frame number that we are currently working on -- this is used in the write_const_MVs() function, and is only for testing purposes -- can be removed later.
	int frame_num;
	//Neighborhood size for regularization
	static const int neighborhood_size = 9;//27; -- use 9 for 2-D and 27 for 3-D

	//For multi-resolution image pyramid
	typedef itk::MultiResolutionPyramidImageFilter<ImageType, ImageType> ImagePyramid;
	ImagePyramid::Pointer pyramid_pointer1; //for the first volume
	ImagePyramid::Pointer pyramid_pointer2; //for the second volume
	//For padding images to make sure we have an even number of blocks
	typedef itk::ZeroFluxNeumannPadImageFilter<ImageType, ImageType> PadFilter;
	//typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadFilter;
	PadFilter::Pointer pad_pointer1;
	PadFilter::Pointer pad_pointer2;
	
	//For DICOM Images -- probably not needed
	/*typedef itk::GDCMImageIO ImageIOType;
    ImageIOType::Pointer gdcmImageIO; */
	//Output image pointers
	ImageType::Pointer outputImage;
	
	//Write for output file
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer;
	//Iterators
	typedef itk::ImageRegionIterator<ImageType> ImageRegionIterator;	
	
	//The parameters below are used in the block matching step
	ImageType::RegionType::SizeType block_size; //region of size block size
	ImageRegionIterator::OffsetType block_offset; //block size
    ImageRegionIterator::OffsetType search_offset; //search size
    ImageRegionIterator::OffsetType begin_offset; //just all zeros -- where the image starts
    ImageRegionIterator::OffsetType end_offset; //marks the end of the image 
	ImageType::RegionType::IndexType end_boundary; //to stop from going outside of boundaries in image 2
	ImageType::RegionType::IndexType end_boundary_prev; //this is the end boundary from the previous level of pyramid
	ImageType::RegionType::IndexType block_boundary; //to stop from going outside of block for splitting

	//To hold and reference the current block that is being searched in image 1
	ImageType::RegionType block_region1; //region size used in block matching for image 1
	ImageType::RegionType block_region2; //region size used in block matching for image 2

	//These parameters are used in block matching step as well, specifically in the find_block_match() function
	ImageType::RegionType::IndexType current_index; //this subtracts the search size from the current index to get the initial start index -- checked by variable below
	ImageType::RegionType::IndexType start_index; //this is the start index used in block matching, but there is a check to make sure it doesnt become negative
	ImageType::RegionType::SizeType search_size; //sets the search size, but makes sure that we do not go outside the boundaries of the volume
	ImageType::RegionType search_region; //region size set by search size above
		
	//These are for splitting of the blocks/MVs
	ImageType::RegionType::IndexType block_start;
	ImageType::RegionType split_region;

	ImageType::RegionType block_region1_reg;
	ImageType::RegionType block_region2_reg;
	
	ImageType::RegionType MC_region;
	ImageType::RegionType out_region;
				
	//-----------Private functions-------------------
	//Function used to find matching block in search region
	void create_pyramid(); //create a pyramid of multi-resolution volumes
	void pad_volumes(); //pad volumes to make sure that an even number of cubes fit in the volume	
	void find_block_match(ImageRegionIterator &image2, int *& match_vec);	
	void find_block_match_pyramid(ImageRegionIterator &image2, int *& match_vec, std::vector<MVs_PL::MVIterator>& previous_index);
	void init_MVs(std::vector<MVs_PL::MVIterator>& current_index);
	void init_MVs_block(std::vector<MVs_PL::MVIterator>& current_index);
	void increment_index(itk::ImageRegionIterator<Registration::ImageType> &current_index, int increment_amount, int &done, Registration::ImageType::RegionType::IndexType boundary);
	void increment_neigh_index(itk::NeighborhoodIterator<MVs_PL::Vector_Array, MVs_PL::BoundaryConditionType> &current_index, int increment_amount, int &done, MVs_PL::Vector_Array::RegionType::IndexType boundary);
	void increment_block(std::vector<MVs_PL::MVIterator> &current_index, int increment_amount, int &done, Registration::ImageType::RegionType::IndexType boundary);
	void increment_MVs(std::vector<MVs_PL::MVIterator>& current_index, int increment_amount, int &done, Registration::ImageType::RegionType::IndexType boundary);
	void increment_all(std::vector<MVs_PL::NeighborhoodIteratorType> &current_index, int b_size, int &done, Registration::ImageType::RegionType::IndexType boundary);
	void increment_allMVs(std::vector<MVs_PL::MVIterator>& MVs);
	void split_MVs(int split_bsize, MVs_PL *& class_element);
	void write_const_MVs();
	void regularization(int lambda);
	int calc_data(MVs_PL::NeighborhoodIteratorType::IndexType frame2_index, int j, int MV[][neighborhood_size]);//, int x_MV, int y_MV, int z_MV);
    int calc_reg(int current_candidate, int MV[][neighborhood_size]);
	int calc_reg_boundary(int current_candidate, int MV[][neighborhood_size], std::vector<int>& valid_faces);
	void calculate_validity();
	void find_min(int * array_vals, int size, int &min_index);
	int min(int x, int y);
	int max(int x, int y);	

};

int compare(const void * a, const void * b);

//template<class T> T Registration::min(const T& x, const T& y)
//{
//  return !(y < x) ? x : y;
//}
//template<class T> T Registration::max(const T& x, const T& y)
//{
//  return (x < y) ? y : x; 
//}

#endif