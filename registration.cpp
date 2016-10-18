#include "registration.h"

Registration::Registration(char *image1_file, char *image2_file, int * bsize, int * ssize, int num_levels, int interp_factor, int f_num)
{ 

  //Open output file which will just be used for diagnostic purposes
  fout.open("motion_vectors.txt");

  //We do a check here to make sure that the block size of the current level is smaller or equal to the block size of the previous level
  /*for(int i = 1; i < num_levels; ++i)
  {
    if (bsize[i-1] < bsize[i])
	{
	  std::cerr << "The block size at the current level must be the same or smaller than the block size at the previous level" << std::endl;
	  exit(1);
	}
  }*/

  //Set # of pyramid levels, block and search sizes -- note that search size only defines the number of pixels on a single side, not both sides.
  levels = num_levels;
  b_size = new int [levels]; //create dynamic array for block sizes
  s_size = new int [levels]; //create dynamic array for search sizes
  //Copy the arrays from main
  for(int i = 0; i < num_levels; ++i)
  {
    b_size[i] = bsize[i];
    s_size[i] = ssize[i];
  }

  //Saved interpolation factor
  i_factor = interp_factor;

  //Save frame number
  frame_num = f_num;

  //Reader for input files
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader1;
  ReaderType::Pointer reader2;
  reader1 = ReaderType::New(); //reader for input image 1
  reader2 = ReaderType::New(); //reader for input image 2

  //For DICOM images -- probably not needed
  //gdcmImageIO = ImageIOType::New();
  //reader1->SetImageIO(gdcmImageIO);

  //Pass file names to readers
  reader1->SetFileName(image1_file); 
  reader2->SetFileName(image2_file); 

  //Try to update the readers or catch the failure
  try
  {
    reader1->Update();
	reader2->Update();
	inputImage1 = reader1->GetOutput();
	outputImage = reader1->GetOutput();
	//write_output("./3d_volumes/sam.dcm");
	inputImage2 = reader2->GetOutput();
  }
  catch(itk::ExceptionObject& err)
  {
	std::cerr << "Error in one of file readers" << std::endl;
    std::cerr << err << std::endl;
	exit(1);
  }

}

void Registration::motion_estimation()
{
  //Create multi-resolution image pyramid
  create_pyramid(); //create a three-level pyramid
  
  int i,j,k; //just so we don't have to keep reinitializing the variables.

  //Handling the first (lowest-resolution) level of the pyramid
  current_level = 0;
  pad_volumes(); //pad the volumes for the first level of the pyramid
       
  //Iterator for image 2 -- this will be the control image -- "current image"
  ImageRegionIterator image2(pad_pointer2->GetOutput(), pad_pointer2->GetOutput()->GetRequestedRegion());
   
  //Block size -- set 3D block size
  //Block offset for incrementing blocks in image 2 
  //Create search offset -- we will specify how much to search to the left and to the right of each block as a number
  //Begin offset -- this is just so we can prevent the search area from going outside the beginning of the image
  //End offset -- this is just so we can prevent the search area from going outside the beginning of the image
  //Offset for blocks in block matching search step
  ImageType::RegionType::IndexType block_index; 
  for(i = 0; i < Dimension; ++i)
  {
	block_size[i] = b_size[current_level];
    block_offset[i] = b_size[current_level];
	search_offset[i] = s_size[current_level];
	begin_offset[i] = 0;
	end_offset[i] = pad_pointer2->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
	end_boundary[i] = (pad_pointer2->GetOutput()->GetLargestPossibleRegion().GetSize()[i]/b_size[current_level])*b_size[current_level];
	block_index[i] = b_size[current_level];
  }

  //set sizes of region that will be used in find_block_match function
  block_region1.SetSize(block_size);
  block_region2.SetSize(block_size);
    
  //These are the positions of the block found in image 1.
  int *match_vec = new int[Dimension]; //Need to free this
  
  int done = 0;
  int done2 = 0;
  //Do the block matching here.  
  for (image2.GoToBegin(), init_MVs(MV_PLNum[0]->MVs); !done; increment_index(image2, b_size[current_level], done, end_boundary), increment_MVs(MV_PLNum[0]->MVs, b_size[current_level], done2, end_boundary))
  {
	//For the block in image 2, we want to search all positions inside the search window in image 1
	find_block_match(image2, match_vec);

	for(i = 0; i < Dimension; ++i)
	{
	  MV_PLNum[0]->MVs[i].Set(match_vec[i]); //MV_PLNum[0] is the first level of the hierarchy -- we should always have a first level, so the '0' is fine
	}	
  }

  //We want to do the splitting of MVs here.  Note that if this is the last level of the hierarchy, we want to have a MV for every pixel in the block -- do this at the end
  int split_bsize; //to hold variable for splitting
  int temp_lambda; //for lagrange multiplier in regularization
  int temp_bsize; //for storing the block size
  
  //START OF CODE FOR HANDLING NEXT LEVELS IN PYRAMID
  //Note that we need more MVs at the next higher resolution level of the hierarchy, so we will have to get these
  //by splitting the blocks at the previous level of the hierarchy.  For example, if we have 8x8 blocks at the
  //previous level of the hierarchy, and we have 4x4 blocks at the current level of the hierarchy, then we 
  //actually need motion vectors for 2x2 blocks at the previous level of the hierarchy.  This makes sense because
  //as we go up to the next higher resolution level of the hierarchy, we are increasing the resolution by a factor
  //of two.  This increase in resolution means that a 2x2 block at the previous level of the hierarchy is a 4x4
  //block at the current level of the hierarchy.

  //We first split the blocks at the previous level.  The smaller/split blocks will take on the MVs of the 
  //parent blocks for now.  Later, we will apply regularization each time before we do any block splitting.
  //In order to split the blocks, we have to know the block size at the current level and the block size at the
  //previous level

  //We want to do the splitting of MVs here.  Note that if this is the last level of the hierarchy, we want to have a MV for every pixel in the block -- do this at the end
  //Perform regularization and splitting -- we go all the way down to block sizes of 2x2, 3 iterations overall, and 4 iterations at each block size = 12 total iterations
  for(int k = 0; k < 3; ++k)
  {
    temp_bsize = b_size[current_level];
	temp_lambda = (temp_bsize*3) >> 1;
	while(temp_bsize > 1)
	{
	  for(int l = 0; l < 4; ++l)
	  {
	    regularization((l+1)*temp_lambda);
	  }
	  split_MVs(temp_bsize, MV_PLNum[current_level]);
	  temp_bsize = (temp_bsize >> 1);
	  temp_lambda = temp_lambda >> 1;
	}
  }
  split_bsize = (b_size[1] >> 1); //b_size of next level divided by 2. Split bsize is the bsize at the next level divided by 2
  //split_MVs(split_bsize, MV_PLNum[0]); //pass the motion vectors to be split
    
  //We are now ready to handle any additional levels of the image pyramid -- assuming that there is more than one level
  for(i = 1; i < levels; ++i) //This for loop will be used to do block matching at the remaining levels (if any)
  { 	     			
	//Set current level and create padded volumes 
	current_level = i;
	pad_volumes(); //pad the volumes for current level of pyramid
		  
    //Block size -- set 3D block size
	//Block offset for incrementing blocks in image 2 
    //Create search offset -- we will specify how much to search to the left and to the right of each block as a number
    //Begin offset -- this is just so we can prevent the search area from going outside the beginning of the image
    //End offset -- this is just so we can prevent the search area from going outside the beginning of the image
    for(j = 0; j < Dimension; ++j)
    {
	  block_size[j] = b_size[current_level];
      block_offset[j] = b_size[current_level];
	  search_offset[j] = s_size[current_level];
	  begin_offset[j] = 0;
	  end_offset[j] = pad_pointer2->GetOutput()->GetLargestPossibleRegion().GetSize()[j];
	  end_boundary_prev[j] = end_boundary[j];
	  end_boundary[j] = (pad_pointer2->GetOutput()->GetLargestPossibleRegion().GetSize()[j]/b_size[current_level])*b_size[current_level];
	  block_index[j] = b_size[current_level];
    }

	//set sizes of region
    block_region1.SetSize(block_size);
    block_region2.SetSize(block_size);

	//Iterator for image 2
    ImageRegionIterator image2(pad_pointer2->GetOutput(), pad_pointer2->GetOutput()->GetRequestedRegion());

	int done = 0;
    int done2 = 0;
	int done3 = 0;
    //Do the block matching here.  
    for (image2.GoToBegin(), init_MVs(MV_PLNum[current_level]->MVs), init_MVs(MV_PLNum[current_level-1]->MVs); !done; increment_index(image2, b_size[current_level], done, end_boundary), increment_MVs(MV_PLNum[current_level]->MVs, b_size[current_level], done2, end_boundary), increment_MVs(MV_PLNum[current_level-1]->MVs, split_bsize, done3, end_boundary_prev))
    {
      //For the block in image 2, we want to search all positions inside the search window in image 1 using the MVs from previous level
	  //We need to pass the MV position from the previous frame to the function below
	  find_block_match_pyramid(image2, match_vec, MV_PLNum[current_level-1]->MVs);

	  for(k = 0; k < Dimension; ++k)
	  {
	    MV_PLNum[current_level]->MVs[k].Set(match_vec[k]);
	  }	
    }
	
	//Perform regularization and splitting -- we go all the way down to block sizes of 2x2, 3 iterations overall, and 4 iterations at each block size = 12 total iterations
	for(int k = 0; k < 3; ++k)
    {
      temp_bsize = b_size[current_level];
	  temp_lambda = (temp_bsize*3) >> 1;
	  while(temp_bsize > 1)
	  {
	    for(int l = 0; l < 4; ++l)
	    {
	      regularization((l+1)*temp_lambda);
	    }
	    split_MVs(temp_bsize, MV_PLNum[current_level]);
	    temp_bsize = (temp_bsize >> 1);
	    temp_lambda = temp_lambda >> 1;
	  }
    }
    split_bsize = (b_size[min(current_level + 1, (levels-1))] >> 1);
	//split_MVs(split_bsize, MV_PLNum[current_level]); //pass the motion vectors to be split

  }

  //At the very end, we want a MV for every pixel, so we need to do one more splitting of the MVs with an increment of one.
  split_MVs(1, MV_PLNum[levels-1]); //split MVs
  //write_const_MVs();
  //if (frame_num == 15) {
  //print_MVs(1);
  //std::cout << "Hit enter to go to next set of frames" << std::endl;
  //getchar(); }
  //draw_MVs_MATLAB(16);

  //Creating an array of values which indicate when the MV is valid or invalid
  valid = Validity::New();
  valid->SetRegions(pad_pointer2->GetOutput()->GetRequestedRegion());
  valid->Allocate(); 
  valid->FillBuffer(0); //initialize all values to 0

  calculate_validity(); //this only needs to be done at the end of the motion estimation process once we have a MV for each pixel

  delete [] match_vec;
}

void Registration::scale_MVs(int scale_factor)
{

  //need to create a MV array that is scale factor * size of the highest level of the hierarchy
  
  ImageType::RegionType new_array; //this will hold the region variable of the new array
  ImageType::RegionType::SizeType new_size; //this will hold the size of the region
  for (int j = 0; j < Dimension; ++j)
    new_size[j] = scale_factor*MV_PLNum[levels-1]->MV_array[0]->GetLargestPossibleRegion().GetSize()[j];

  new_array.SetSize(new_size);

  //Add new element to class that will hold the arrays and then fill the class with the array
  MV_alias.push_back(new MVs_PL); //This will need to be freed!
  //Fill classes with arrays -- we create an array for each component of the MV -- x, y, and z dimensions for the 3-D case
  for(int j = 0; j < Dimension; ++j)
  {
    MV_alias[0]->MV_array.push_back(MVs_PL::Vector_Array::New());
	
	MV_alias[0]->MV_array[j]->SetRegions();
	MV_alias[0]->MV_array[j]->Allocate();
	MV_alias[0]->MV_array[j]->FillBuffer(0);

	//Create iterators -- three for each MV 
	MVs_PL::MVIterator temp_it(MV_alias[0]->MV_array[j], MV_alias[0]->MV_array[j]->GetLargestPossibleRegion());
	MV_alias[0]->MVs.push_back(temp_it);
  }

 



  int done = 0;
  for (int j = 0; j < Dimension; ++j)
    end_boundary[j] = MV_PLNum[levels-1]->MV_array[0]->GetLargestPossibleRegion().GetSize()[j];
   
  for(init_MVs(MV_PLNum[levels-1]->MVs); !done; increment_MVs(MV_PLNum[levels-1]->MVs, 1, done, end_boundary))
  {  
  }

}

void Registration::create_pyramid()
{
  
  if (i_factor > 1)
  {
    //The code below is used for interpolating the volume !!NEED TO CHECK IF INTERPOLATION FACTOR == 1, IF SO, SKIP INTERPOLATION
    //Filter for interpolating images
    typedef itk::ResampleImageFilter<ImageType, ImageType> FilterType;
    //Create filter which will be used for interpolation
    FilterType::Pointer filter1 = FilterType::New();
    FilterType::Pointer filter2 = FilterType::New();
    //Use identity transform when doing interpolation -- we don't want to do anything BUT interpolate
    typedef itk::IdentityTransform<double, Dimension> TransformType;
    //Set up identity transform for the interpolation
    TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
    filter1->SetTransform(transform);
    filter2->SetTransform(transform);
    //Interpolation function for filter
    typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
    //Create interpolation function -- linear interpolater
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    filter1->SetInterpolator(interpolator);
    filter2->SetInterpolator(interpolator);
    //Set up spacing to do the interpolation
    ImageType::SpacingType inputSpacing = inputImage1->GetSpacing();
    ImageType::SizeType inputSize = inputImage1->GetLargestPossibleRegion().GetSize();
    ImageType::SpacingType outputSpacing = inputSpacing;
    ImageType::SizeType outputSize = inputSize; 
    for (int i = 0; i < Dimension; ++i)
    {
      outputSpacing[i] /= i_factor; //changes output spacing to be smaller -- more pixels
      outputSize[i] = static_cast<int>(outputSize[i] * i_factor);
    }
    //More parameters for interpolation filter
    filter1->SetSize(outputSize);
    filter2->SetSize(outputSize);
    filter1->SetOutputSpacing(outputSpacing);
    filter2->SetOutputSpacing(outputSpacing);
    filter1->SetInput(inputImage1);
    filter2->SetInput(inputImage2);
    filter1->SetDefaultPixelValue(0); //value of pixels outside of the image
    filter2->SetDefaultPixelValue(0); //value of pixels outside of the image
    filter1->SetOutputOrigin(inputImage1->GetOrigin());
    filter2->SetOutputOrigin(inputImage2->GetOrigin());

    try //run the interpolation process
    {
      filter1->Update();	
	  filter2->Update();	
	  previous_frame = filter1->GetOutput(); //this is used for the first frame in the argument list
	  current_frame = filter2->GetOutput(); //this is used so we don't have to interpolate again for the POCS function
	  //Set output image to be the one of the levels from pyramid_pointer
      outputImage = filter1->GetOutput();
	  //write_output("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/interpolated_only11.dcm");
    }
    catch(itk::ExceptionObject& err)
    {
	  std::cerr << "Not able to create interpolated volumes" << std::endl;
      std::cerr << err << std::endl;
	  getchar();
	  exit(1);
    }

	//Create image pyramid
    pyramid_pointer1 = ImagePyramid::New();
    pyramid_pointer2 = ImagePyramid::New();
    pyramid_pointer1->SetInput(filter1->GetOutput()); //the input to the pyramid will be the interpolated volumes
    pyramid_pointer2->SetInput(filter2->GetOutput());
    pyramid_pointer1->SetNumberOfLevels(levels);
    pyramid_pointer2->SetNumberOfLevels(levels);

  }
  else //no interpolation
  {
	previous_frame = inputImage1; //this is used for the first frame in the argument list
	current_frame = inputImage2; //this is used so we don't have to interpolate again for the POCS function
	
	//Set output image to be the one of the levels from pyramid_pointer
    outputImage = inputImage1;
    
	//Create image pyramid
    pyramid_pointer1 = ImagePyramid::New();
    pyramid_pointer2 = ImagePyramid::New();
    pyramid_pointer1->SetInput(inputImage1); //the input to the pyramid will be the interpolated volumes
    pyramid_pointer2->SetInput(inputImage2);
    pyramid_pointer1->SetNumberOfLevels(levels);
    pyramid_pointer2->SetNumberOfLevels(levels);
  } 

  try //run the process to create a pyramid of images
  {
    pyramid_pointer1->Update();	
	pyramid_pointer2->Update();	
  }
  catch(itk::ExceptionObject& err)
  {
	std::cerr << "Error in creation of multi-resolution pyramid" << std::endl;
    std::cerr << err << std::endl;
	getchar();
	exit(1);
  }

  //Create output image
  /*outputImage = ImageType::New();
  outputImage->SetRegions(inputImage1->GetLargestPossibleRegion());
  outputImage->Allocate(); 
  outputImage->FillBuffer(0);

  //Set output image to be the one of the levels from pyramid_pointer
  outputImage = pyramid_pointer->GetOutput(2);
  */

}

void Registration::pad_volumes()
{
  //First, we need to determine how much to pad the volumes in each dimension
  unsigned long pad_amount[Dimension];
  unsigned long zero_pad[Dimension];

  for(int i = 0; i < Dimension; ++i)
  {
    pad_amount[i] = (b_size[current_level] - (pyramid_pointer1->GetOutput(current_level)->GetLargestPossibleRegion().GetSize()[i] - (((int)pyramid_pointer1->GetOutput(current_level)->GetLargestPossibleRegion().GetSize()[i]/b_size[current_level])*b_size[current_level]))) % b_size[current_level];
	if (current_level > 0) //statement below basically says that we need to make sure that the next level in the pyramid is at LEAST twice the size of the previous level.  The statement just figures out how much padding to do in order for this condition to be true.
	  pad_amount[i] = max(pyramid_pointer1->GetOutput(current_level)->GetLargestPossibleRegion().GetSize()[i] + pad_amount[i], pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[i]*2) - pyramid_pointer1->GetOutput(current_level)->GetLargestPossibleRegion().GetSize()[i];
	//std::cout << "pad amount " << i << " is " << pad_amount[i] << std::endl;
	zero_pad[i] = 0; //lower bound padding will always be zero -- we pad the end of the volume
  }

  pad_pointer1 = PadFilter::New();
  pad_pointer2 = PadFilter::New();
  pad_pointer1->SetInput(pyramid_pointer1->GetOutput(current_level));
  pad_pointer2->SetInput(pyramid_pointer2->GetOutput(current_level));
  pad_pointer1->SetPadLowerBound(zero_pad);
  pad_pointer2->SetPadLowerBound(zero_pad);
  pad_pointer1->SetPadUpperBound(pad_amount);
  pad_pointer2->SetPadUpperBound(pad_amount);

  try
  {
    pad_pointer1->Update();	
	pad_pointer2->Update();	
  }
  catch(itk::ExceptionObject& err)
  {
	std::cerr << "Error in padding of volumes" << std::endl;
    std::cerr << err << std::endl;
	getchar();
	exit(1);
  }

  //Add new element to class that will hold the arrays and then fill the class with the array
  MV_PLNum.push_back(new MVs_PL); //This will need to be freed!
  //Fill classes with arrays -- we create an array for each component of the MV -- x, y, and z dimensions for the 3-D case
  for(int j = 0; j < Dimension; ++j)
  {
    MV_PLNum[current_level]->MV_array.push_back(MVs_PL::Vector_Array::New());
	MV_PLNum[current_level]->MV_array[j]->SetRegions(pad_pointer2->GetOutput()->GetRequestedRegion());
	MV_PLNum[current_level]->MV_array[j]->Allocate();
	MV_PLNum[current_level]->MV_array[j]->FillBuffer(0);

	//Create iterators -- three for each MV 
	MVs_PL::MVIterator temp_it(MV_PLNum[current_level]->MV_array[j], MV_PLNum[current_level]->MV_array[j]->GetLargestPossibleRegion());
	MV_PLNum[current_level]->MVs.push_back(temp_it);
  }

}

void Registration::find_block_match(ImageRegionIterator &image2, int *& match_vec)
{  
  int i; //Don't have to keep initializing

  //Current index -- The beginning index of the search window is the current index minus the search offset -- how many pixels to search on one side of the block
  //Start index -- Set the begin index so that it does not go outside the image (negative values)
  //Set size of search iterator -- min functions make sure that the search size does not go outside the image (larger than width, height, depth)
  //We really have ssize*2 + b_size but then we subtract out b_size so that we don't go beyond the search region
  for(int i = 0; i < Dimension; ++i)
  {
	current_index[i] = image2.GetIndex().GetElement(i) - search_offset[i];  
    start_index[i] = max(begin_offset[i], current_index[i]);
    search_size[i] = min(end_offset[i] - block_offset[i] + 1, current_index[i] + s_size[current_level]*2) - start_index[i]; //instead of + 1, we may need 0 here depending on how IsAtEnd() works in loop below.
  }
  search_region.SetIndex(start_index);
  search_region.SetSize(search_size);

  //Create search region iterator -- to navigate through the whole search region
  ImageRegionIterator searchregion_it(pad_pointer1->GetOutput(), search_region);

  //Create image 2 block iterator -- a block iterator inside the larger search region
  block_region2.SetIndex(image2.GetIndex()); //this sets the starting index, the size was set in the block_matching() function
  ImageRegionIterator image2block_it(pad_pointer2->GetOutput(), block_region2);

  //Initialize SAD values
  int min_SAD = std::numeric_limits<int>::max(); //set the minimum SAD to the maximum possible integer value
  int min_distance = std::numeric_limits<int>::max(); //set the minimum distance to the maximum possible integer value
  int SAD; //to hold the current SAD value
  int distance; //hold the distance from the block in the search window to the center of the search window
  for(i = 0; i < Dimension; ++i)
  {
    match_vec[i] = 0;
  }

  for(searchregion_it.GoToBegin(); !searchregion_it.IsAtEnd(); ++searchregion_it)
  { 
	//Create image 1 block iterator
	block_region1.SetIndex(searchregion_it.GetIndex()); //this sets the starting index, the site was set in the block_matching() function
    ImageRegionIterator image1block_it(pad_pointer1->GetOutput(), block_region1); //block region size set in the block_matching() function

	//Initialize SAD to 0
	SAD = 0;
	
	//Compute SAD between blocks in image1 and image2
	for(image2block_it.GoToBegin(); !image2block_it.IsAtEnd(); ++image2block_it, ++image1block_it)
	{
	  SAD += abs(image2block_it.Get() - image1block_it.Get());
	}

    //Check if the current SAD is less than the current minimum SAD -- if so, update the min_SAD value.
	if (SAD < min_SAD)
	{
	  //Store the beginning x, y, z position of the block with the minimum SAD
	  image1block_it.GoToBegin(); //this will let us set the beginning index of block 1 
	  image2block_it.GoToBegin(); //this will let us set the beginning index of block 2

	  min_distance = 0; //reset min distance so we can use the += operator below
	  for(i = 0; i < Dimension; ++i)
	  {
	    match_vec[i] = image1block_it.GetIndex()[i]; //this is the position in image 1
		min_distance += abs(match_vec[i] - image2block_it.GetIndex()[i]); //update the distance for the minimum SAD value;
	  }
	  min_SAD = SAD;
	}
	else if (SAD == min_SAD)
	{
      //if the SAD is the same as the minimum SAD, then we change the minimum SAD if the block is closer to the block in frame 2 -- i.e., the distance from the block to the center of the search window is smaller
	  image1block_it.GoToBegin();
	  image2block_it.GoToBegin();
	  	
	  distance = 0;
	  for(i = 0; i < Dimension; ++i)
	    distance += abs(image1block_it.GetIndex()[i] - image2block_it.GetIndex()[i]);

	  if (distance < min_distance)
	  {
		min_distance = distance;
	    for(i = 0; i < Dimension; ++i)
	      match_vec[i] = image1block_it.GetIndex()[i]; //this is the position in image 1
	  }
	}
  } 
}

void Registration::find_block_match_pyramid(ImageRegionIterator &image2, int *& match_vec, std::vector<MVs_PL::MVIterator>& previous_index)
{  
  int i; //to prevent having to keep re-initializing
  //Current index -- The beginning index of the search window is given from the previous level's motion vector times 2, minus the search offset -- how many pixels to search on one side of the block
  //Start index -- Set the begin index so that it does not go outside the image (negative values)
  //Set size of search iterator -- min functions make sure that the search size does not go outside the image (larger than width, height, depth)
  //We really have ssize*2 + b_size but then we subtract out b_size so that we don't go beyond the search region
  for(int i = 0; i < Dimension; ++i)
  {
	current_index[i] = previous_index[i].Get()*2 - search_offset[i];  //index in previous level multipled by 2 for next level - search offset of next level
    start_index[i] = max(begin_offset[i], current_index[i]); //make sure that current index doesn't go outside of volume boundaries
    search_size[i] = min(end_offset[i] - block_offset[i] + 1, current_index[i] + s_size[current_level]*2) - start_index[i];  
  }
  search_region.SetIndex(start_index);
  search_region.SetSize(search_size);

  //Create search region iterator
  ImageRegionIterator searchregion_it(pad_pointer1->GetOutput(), search_region);

  //Create image 2 block iterator
  block_region2.SetIndex(image2.GetIndex()); //this sets the starting index, the size was set in the block_matching() function
  ImageRegionIterator image2block_it(pad_pointer2->GetOutput(), block_region2);

  //Initialize SAD values
  int min_SAD = std::numeric_limits<int>::max(); //set the minimum SAD to the maximum possible integer value
  int min_distance = std::numeric_limits<int>::max(); //set the minimum distance to the maximum possible integer value
  int SAD; //to hold the current SAD value
  int distance; //hold the distance from the block in the search window to the center of the search window
  for(i = 0; i < Dimension; ++i)
  {
    match_vec[i] = 0;
  }

  for(searchregion_it.GoToBegin(); !searchregion_it.IsAtEnd(); ++searchregion_it)
  { 
	//Create image 1 block iterator
	block_region1.SetIndex(searchregion_it.GetIndex()); //this sets the starting index, the site was set in the block_matching() function
    ImageRegionIterator image1block_it(pad_pointer1->GetOutput(), block_region1); //block region size set in the block_matching() function

	//Initialize SAD to 0
	SAD = 0;

	//Compute SAD between blocks in image1 and image2
	for(image2block_it.GoToBegin(); !image2block_it.IsAtEnd(); ++image2block_it, ++image1block_it)
	{
	  SAD += abs(image2block_it.Get() - image1block_it.Get());
	}
    //Check if the current SAD is less than the current minimum SAD -- if so, update the min_SAD value.
	if (SAD < min_SAD)
	{
	  //Store the beginning x, y, z position of the block with the minimum SAD
	  image1block_it.GoToBegin(); //this will let us set the beginning index of block 1 
	  image2block_it.GoToBegin(); //this will let us set the beginning index of block 2

	  min_distance = 0; //reset min distance so we can use the += operator below
	  for(i = 0; i < Dimension; ++i)
	  {
	    match_vec[i] = image1block_it.GetIndex()[i]; //this is the position in image 1
		min_distance += abs(image1block_it.GetIndex()[i] - image2block_it.GetIndex()[i]); //update the distance for the minimum SAD value;
	  }
	  min_SAD = SAD;
	}
	else if (SAD == min_SAD)
	{
      //if the SAD is the same as the minimum SAD, then we change the minimum SAD if the block is closer to the block in frame 2 -- i.e., the distance from the block to the center of the search window is smaller
	  image1block_it.GoToBegin();
	  image2block_it.GoToBegin();
	  	
	  distance = 0;
	  for(i = 0; i < Dimension; ++i)
	    distance += abs(image1block_it.GetIndex()[i] - image2block_it.GetIndex()[i]);

	  if (distance < min_distance)
	  {
		min_distance = distance;
	    for(i = 0; i < Dimension; ++i)
	      match_vec[i] = image1block_it.GetIndex()[i]; //this is the position in image 1
	  }
	}
  } 
}

void Registration::regularization(int lambda) //need to change how MVs outside boundary are selected
{
  int local_bsize = b_size[current_level]; //so we don't have to keep accessing the memory location
  int i,j; //Avoid multiple declarations in loop
  //Validity::RegionType val_region; //block region for validity metric
  
  MV_PLNum[current_level]->radius.Fill(local_bsize); //create a radius of b_size in all dimensions -- this is to set up a neighborhood of MVs
 
  //Create regularization iterators -- need to modify this to handle any dimension, not just 3-D
  for(i = 0; i < Dimension; ++i)
  {
    MVs_PL::NeighborhoodIteratorType reg_temp(MV_PLNum[current_level]->radius, MV_PLNum[current_level]->MV_array[i], MV_PLNum[current_level]->MV_array[i]->GetLargestPossibleRegion());
	MV_PLNum[current_level]->Reg_MVs.push_back(reg_temp);
  }

  //Set up parameters for figuring out which of the faces of the neighborhood iterator have gone out of bounds.  We will use the ones NOT out of bounds in regularization
  MVs_PL::NeighborhoodIteratorType::OffsetType internalIndex;
  MVs_PL::NeighborhoodIteratorType::OffsetType offset_index;
  std::vector<int> valid_faces;

  //Set up offsets that will be needed to get the MVs we want -- we will use a large neighborhood 26 neighbors
  //Probably should use overloading of an assignmnet operator to do this more efficiently
  typedef MVs_PL::NeighborhoodIteratorType::OffsetType offset_type;
  offset_type offset[neighborhood_size];

  if (Dimension == 3)
  {
    itk::Offset<Dimension>::OffsetValueType test[3];
    
    //This is the rear slice
    test[0] = -local_bsize;
    test[1] = -local_bsize;
    test[2] = local_bsize;
    offset[0].SetOffset(test);
    test[0] = 0;
    test[1] = -local_bsize;
    test[2] = local_bsize;
    offset[1].SetOffset(test);
    test[0] = local_bsize;
    test[1] = -local_bsize;
    test[2] = local_bsize;
    offset[2].SetOffset(test);
    test[0] = -local_bsize;
    test[1] = 0;
    test[2] = local_bsize;
    offset[3].SetOffset(test);
    test[0] = 0;
    test[1] = 0;
    test[2] = local_bsize;
    offset[4].SetOffset(test);
    test[0] = local_bsize; 
    test[1] = 0;
    test[2] = local_bsize;
    offset[5].SetOffset(test);
    test[0] = -local_bsize;
    test[1] = local_bsize;
    test[2] = local_bsize;
    offset[6].SetOffset(test);
    test[0] = 0;
    test[1] = local_bsize;
    test[2] = local_bsize;
    offset[7].SetOffset(test);
    test[0] = local_bsize;
    test[1] = local_bsize;
    test[2] = local_bsize;
    offset[8].SetOffset(test);
  
    //This is the center slice
    test[0] = -local_bsize;
    test[1] = -local_bsize;
    test[2] = 0;
    offset[9].SetOffset(test);
    test[0] = 0;
    test[1] = -local_bsize;
    test[2] = 0;
    offset[10].SetOffset(test);
    test[0] = local_bsize;
    test[1] = -local_bsize;
    test[2] = 0;
    offset[11].SetOffset(test);
    test[0] = -local_bsize;
    test[1] = 0;
    test[2] = 0;
    offset[12].SetOffset(test);
    test[0] = 0;
    test[1] = 0;
    test[2] = 0;
    offset[13].SetOffset(test); //this is the center block
    test[0] = local_bsize;
    test[1] = 0;
    test[2] = 0;
    offset[14].SetOffset(test);
    test[0] = -local_bsize;
    test[1] = local_bsize;
    test[2] = 0;
    offset[15].SetOffset(test);
    test[0] = 0;
    test[1] = local_bsize;
    test[2] = 0;
    offset[16].SetOffset(test);
    test[0] = local_bsize;
    test[1] = local_bsize;
    test[2] = 0;
    offset[17].SetOffset(test);
 
    //This is the front slice
    test[0] = -local_bsize;
    test[1] = -local_bsize;
    test[2] = -local_bsize;
    offset[18].SetOffset(test);
    test[0] = 0;
    test[1] = -local_bsize;
    test[2] = -local_bsize;
    offset[19].SetOffset(test);
    test[0] = local_bsize;
    test[1] = -local_bsize;
    test[2] = -local_bsize;
    offset[20].SetOffset(test);
    test[0] = -local_bsize;
    test[1] = 0;
    test[2] = -local_bsize;
    offset[21].SetOffset(test);
    test[0] = 0;
    test[1] = 0;
    test[2] = -local_bsize;
    offset[22].SetOffset(test);
    test[0] = local_bsize;
    test[1] = 0;
    test[2] = -local_bsize;
    offset[23].SetOffset(test);
    test[0] = -local_bsize;
    test[1] = local_bsize;
    test[2] = -local_bsize;
    offset[24].SetOffset(test);
    test[0] = 0;
    test[1] = local_bsize;
    test[2] = -local_bsize;
    offset[25].SetOffset(test);
    test[0] = local_bsize;
    test[1] = local_bsize;
    test[2] = -local_bsize;
    offset[26].SetOffset(test);
  }
  else if (Dimension == 2)
  {
    itk::Offset<Dimension>::OffsetValueType test[2];

    test[0] = -local_bsize;
    test[1] = -local_bsize;
    offset[0].SetOffset(test);
    test[0] = 0;
    test[1] = -local_bsize;
    offset[1].SetOffset(test);
    test[0] = local_bsize;
    test[1] = -local_bsize;
    offset[2].SetOffset(test);
    test[0] = -local_bsize;
    test[1] = 0;
    offset[3].SetOffset(test);
    test[0] = 0;
    test[1] = 0;
    offset[4].SetOffset(test);
    test[0] = local_bsize; 
    test[1] = 0;
    offset[5].SetOffset(test);
    test[0] = -local_bsize;
    test[1] = local_bsize;
    offset[6].SetOffset(test);
    test[0] = 0;
    test[1] = local_bsize;
    offset[7].SetOffset(test);
    test[0] = local_bsize;
    test[1] = local_bsize;
    offset[8].SetOffset(test);

  }

  //Initialize the block size for the increment functions and set up the terminating "done" condition
  itk::Offset<Dimension> shift_amount;
  int done = 0;
  for(i = 0; i < Dimension; ++i)
  {
    shift_amount[i] = local_bsize;
	block_size[i] = local_bsize;
  }
  
  //For storing the motion vector components and the indices of each iterator -- xindex, yindex, zindex, etc.
  int MVcomp[Dimension][neighborhood_size];
  int indices[Dimension];
  
  //Data term, regularization term, number of candidates to test, and Lagrange multiple, lambda.
  //Min index is the index of the selected candidate
  int data[neighborhood_size];
  int reg[neighborhood_size];
  int candidate[neighborhood_size]; //this will hold the data + lambda*regularization term
  //int lambda = 1;
  int min_index; //the candidate with the minimum overall value which was chosen from the MV array

  //for handling boundaries in neighborhood iterator that go outside the image -- not being used right now
  //itk::Offset<Dimension>::OffsetType offset_index;
  //itk::Offset<Dimension>::OffsetType the_offset; 

  for(i = 0; i < Dimension; ++i)
  {
	MV_PLNum[current_level]->Reg_MVs[i].GoToBegin(); //set the iterators to the beginning positions
  }  

  for( ; !done; increment_all(MV_PLNum[current_level]->Reg_MVs, local_bsize, done, end_boundary))  
  {

	 for(i = 0; i < Dimension; ++i)
	  indices[i] = MV_PLNum[current_level]->Reg_MVs[i].GetIndex()[i]; //get the x, y, or z index...	

	//This actually causes the border blocks to be skipped -- I added this because it was causing the index to 
	//go outside the bounds of the image.  For the border blocks, we need to use only MVs that are in the image 
	//in the data + lambda*reg code.
	if (!MV_PLNum[current_level]->Reg_MVs[0].InBounds()) //continue;
	{
	  valid_faces.clear(); //empty the vector
	  //if we get here, we now need to check which of the indices are out of bounds  //continue;
	  for(i = 0; i < neighborhood_size; ++i)
	  {
	    if(MV_PLNum[current_level]->Reg_MVs[0].IndexInBounds(MV_PLNum[current_level]->Reg_MVs[0].GetNeighborhoodIndex(offset[i]), internalIndex, offset_index))
		{
		   //we get here if the index given by offset[i] is in bounds
		   valid_faces.push_back(i); //we are pushing back which one of the neighbors are valid, which will be used for the offsets below
		}
	  }
	 	 
      for(j = 0; j < valid_faces.size(); ++j)
	  {
	    for(i = 0; i < Dimension; ++i)
	    {
	      MVcomp[i][j] = MV_PLNum[current_level]->Reg_MVs[i].GetPixel(offset[valid_faces[j]]) - (indices[i] + offset[valid_faces[j]].GetOffset()[i]);
	    }
	
	    data[j] = calc_data(MV_PLNum[current_level]->Reg_MVs[0].GetIndex(), j, MVcomp);
	    reg[j] = calc_reg_boundary(j, MVcomp, valid_faces);
	    candidate[j] = data[j] + lambda*reg[j];

	  }   

	  //Find the minimum candidate
	  find_min(candidate, valid_faces.size(), min_index);

	}
	else
	{
      for(j = 0; j < neighborhood_size; ++j)
	  {
	    for(i = 0; i < Dimension; ++i)
	    {
	      MVcomp[i][j] = MV_PLNum[current_level]->Reg_MVs[i].GetPixel(offset[j]) - (indices[i] + offset[j].GetOffset()[i]);
	    }
	
	    data[j] = calc_data(MV_PLNum[current_level]->Reg_MVs[0].GetIndex(), j, MVcomp);
	    reg[j] = calc_reg(j, MVcomp);
	    candidate[j] = data[j] + lambda*reg[j];

	  }   

	  //Find the minimum candidate
	  find_min(candidate, neighborhood_size, min_index);

	}

    //Change the motion vector for all of the iterators (x,y,z)
	for(i = 0; i < Dimension; ++i)
	{
      MV_PLNum[current_level]->Reg_MVs[i].SetPixel(offset[neighborhood_size/2], MVcomp[i][min_index] + indices[i]);
	}

  }
}

int Registration::calc_data(MVs_PL::NeighborhoodIteratorType::IndexType frame2_index, int j, int MV[][neighborhood_size])
{
  //Set up start index for the region that will be used in the iterator
  for(int i = 0; i < Dimension; ++i)
  {
    //start_index[i] = max(0, min(MV[i][j] + frame2_index.GetIndex()[i], pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[i] - block_size[i]));
    start_index[i] = MV[i][j] + frame2_index.GetIndex()[i];
	if ((start_index[i] < 0) || start_index[i] > (pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[i] - block_size[i]))
	  return (std::numeric_limits<int>::max()/2); //return a really large value if it goes outside of the bounds so that this MV will not be used
  }

  //Set the indices and sizes for the region that will be used in the iterator
  block_region1_reg.SetIndex(start_index);
  block_region2_reg.SetIndex(frame2_index);
  block_region1_reg.SetSize(block_size);
  block_region2_reg.SetSize(block_size);

  //Create iterator
  ImageRegionIterator frame1_region(pad_pointer1->GetOutput(), block_region1_reg);
  ImageRegionIterator frame2_region(pad_pointer2->GetOutput(), block_region2_reg);

  //Calculate the sum of absolute differences (SAD) value between the two blocks in the corresponding volumes
  int SAD = 0; //initialize SAD value
  for(frame1_region.GoToBegin(), frame2_region.GoToBegin(); !frame1_region.IsAtEnd(); ++frame1_region, ++frame2_region)
  {
    SAD += abs(frame2_region.Get() - frame1_region.Get());
  }

  return SAD;
}

int Registration::calc_reg(int current_candidate, int MV[][neighborhood_size])
{  
  //The type of regularization used here basically takes one of the neighboring motion vectors and uses it to replace the motion vector of the center block.
  //The current MV that is in the center block is subtracted from all of the other neighbors (including its own, which gives a value of zero).  The L1 norm 
  //is used in the subtraction
  int sum = 0;

  //We break this function up into two for loops because we are relacing the center motion vector with one of the neighbors.  Therefore, the center motion vector should
  //not be used in any calculations except for when we assume the center motion vector is the correct MV.
  for(int j = 0; j < (neighborhood_size/2); ++j) //skip 13, which is the center MV -- this MV is replaced so it shouldnt be counted
  {
    for(int i = 0; i < Dimension; ++i)
	{
	  sum += abs(MV[i][current_candidate] - MV[i][j]);  
	}

  }

  for(int j = (neighborhood_size/2) + 1; j < neighborhood_size; ++j)
  {
    for(int i = 0; i < Dimension; ++i)
	{
      sum += abs(MV[i][current_candidate] - MV[i][j]);
	}
  }
  
  return sum;
}

int Registration::calc_reg_boundary(int current_candidate, int MV[][neighborhood_size], std::vector<int>& valid_faces)
{  
  //The type of regularization used here basically takes one of the neighboring motion vectors and uses it to replace the motion vector of the center block.
  //The current MV that is in the center block is subtracted from all of the other neighbors (including its own, which gives a value of zero).  The L1 norm 
  //is used in the subtraction
  int sum = 0;
  
  //We need to make sure that the current candidate is not subtracted from the MV of the center block in the eight-connected neighbors
  
  for(int j = 0; j < valid_faces.size(); ++j)
  {
	if (valid_faces[j] == neighborhood_size/2 && valid_faces[current_candidate] != neighborhood_size/2) //we are replacing center MV, so we don't include it in the calculation
	  continue;

    for(int i = 0; i < Dimension; ++i)
	{
	  sum += abs(MV[i][current_candidate] - MV[i][j]);  
	}
  }

  return sum;
}

void Registration::calculate_validity()
{
  int i; //prevent having to reinitialize
  //Using the overlap method (see paper "BLOCK-OVERLAP-BASED VALIDITY METRIC FOR HYBRID DE-INTERLACING" by Michael Santoro)
  //Set up an iterator for the counter to sum overlapped pixels in frame 1
  Validity::Pointer count;
  count = Validity::New();
  count->SetRegions(pad_pointer1->GetOutput()->GetLargestPossibleRegion());
  count->Allocate(); 
  count->FillBuffer(0);
  //Create iterator for count
  ValIterator counter(count, pad_pointer1->GetOutput()->GetLargestPossibleRegion());

  //We need to set up an image iterators to calculate the absolute value of pixel difference between frame 2 and frame 1
  ImageRegionIterator image1(pad_pointer1->GetOutput(), pad_pointer1->GetOutput()->GetLargestPossibleRegion());
  ImageRegionIterator image2(pad_pointer2->GetOutput(), pad_pointer2->GetOutput()->GetLargestPossibleRegion());

  //Create iterator to go through validity array
  ValIterator v_it(valid, pad_pointer1->GetOutput()->GetLargestPossibleRegion());
  
  //Set up an iterator that goes through the MV positions in frame 2 -- these should already be set up from previous block matching functions
  //The first loop is just going to count up the pixel overlap (if any).  The second loop will be needed to calculate the confidence value based on the total amount of overlap.
  Validity::RegionType::IndexType index_fromMVs;
  for(init_MVs(MV_PLNum[current_level]->MVs); !MV_PLNum[current_level]->MVs[0].IsAtEnd(); increment_allMVs(MV_PLNum[current_level]->MVs))//++MV_PLNum[current_level]->MVs[0], ++MV_PLNum[current_level]->MVs[1], ++MV_PLNum[current_level]->MVs[2])
  {
	for(i = 0; i < Dimension; ++i)
	{
	  index_fromMVs[i] = max(0, min(MV_PLNum[current_level]->MVs[i].Get(), pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[i] - 1));
	}
    //index_fromMVs[0] = max(0, min(MV_PLNum[current_level]->MVs[0].Get(), pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[0] - 1)); //this is the x-index in frame 1
	//index_fromMVs[1] = max(0, min(MV_PLNum[current_level]->MVs[1].Get(), pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[1] - 1)); //this is the y-index in frame 1
	//index_fromMVs[2] = max(0, min(MV_PLNum[current_level]->MVs[2].Get(), pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1)); //this is the z-index in frame 1

	//Count up the number of overlapping pixels
    counter.SetIndex(index_fromMVs);
	counter.Value()++;
  }

  unsigned int confidence_val;
  for(v_it.GoToBegin(), image2.GoToBegin(), init_MVs(MV_PLNum[current_level]->MVs); !v_it.IsAtEnd(); increment_allMVs(MV_PLNum[current_level]->MVs), ++image2, ++v_it)
  {
	for(i = 0; i < Dimension; ++i)
	{
	  index_fromMVs[i] = max(0, min(MV_PLNum[current_level]->MVs[i].Get(), pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[i] - 1));
	}
	//index_fromMVs[0] = max(0, min(MV_PLNum[current_level]->MVs[0].Get(), pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[0] - 1)); //this is the x-index in frame 1
	//index_fromMVs[1] = max(0, min(MV_PLNum[current_level]->MVs[1].Get(), pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[1] - 1)); //this is the y-index in frame 1
	//index_fromMVs[2] = max(0, min(MV_PLNum[current_level]->MVs[2].Get(), pad_pointer1->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1)); //this is the z-index in frame 1
	
	//Counter index is the position in frame 1 that counted all of the overlapping pixels
	counter.SetIndex(index_fromMVs);
	//Image 1 index is the same as the counter index -- we need this to calculate the absolute pixel difference
	image1.SetIndex(index_fromMVs);

    //Image 2 index is the same as the validity index -- don't need to set anything
	
	//The confidence value is basically the overlap multiplied by the absolute pixel difference, i.e., confidence = overlap*(1 + abs(pixel_diff));  The '1' is there so that the overlap will be taken into account even if the 
	//absolute value of the pixel difference is zero.
	confidence_val = counter.Get()*(1 + abs(image2.Get() - image1.Get()));
	//Set the validity value
    v_it.Set(confidence_val);
  } 
}

void Registration::find_min(int * array_vals, int size, int &min_index)
{
   int min_value = array_vals[0]; // initialize 
   min_index = 0;
   
   for (int i = 1; i < size; ++i)
   {
     if(min_value > array_vals[i])
     {
	    min_value = array_vals[i];
	    min_index = i;
     }
   }
}

void Registration::init_MVs(std::vector<MVs_PL::MVIterator>& current_index)
{
  for(int i = 0; i < Dimension; ++i)
  {	
    current_index[i].GoToBegin();	
  }
}

void Registration::init_MVs_block(std::vector<MVs_PL::MVIterator>& current_index)
{
  for(int i = 0; i < Dimension; ++i)
  {	
    current_index[i].GoToBegin();		
  }
}

void Registration::increment_index(itk::ImageRegionIterator<Registration::ImageType> &current_index, int increment_amount, int &done, Registration::ImageType::RegionType::IndexType boundary)
{
  done = 0;
  ImageType::RegionType::IndexType temp_index;

  for(int i = 0; i < Dimension; ++i)
  {
    if ((current_index.GetIndex()[i] + increment_amount) < end_boundary[i])
	{

	  temp_index = current_index.GetIndex();
	  current_index.GoToBegin();
	  for(int j = i-1; j > -1; --j)
	    temp_index[j] = current_index.GetIndex()[j];
	  temp_index[i] += increment_amount;
	  current_index.SetIndex(temp_index);
	  return;
	 
	}
  } 

  done = 1; //if we get all the way through without returning, set done = 1;
}

void Registration::increment_neigh_index(MVs_PL::NeighborhoodIteratorType &current_index, int increment_amount, int &done, MVs_PL::Vector_Array::RegionType::IndexType boundary)
{
  done = 0;
  ImageType::RegionType::IndexType temp_index;

  for(int i = 0; i < Dimension; ++i)
  {
    if ((current_index.GetIndex()[i] + increment_amount) < end_boundary[i])
	{
      temp_index = current_index.GetIndex();
	  current_index.GoToBegin();
	  for(int j = i-1; j > -1; --j)
	    temp_index[j] = current_index.GetIndex()[j];
	  temp_index[i] += increment_amount;
	  current_index.SetLocation(temp_index);  //this is supposed to be inefficient -- implement another method.
	  return;
	 
	}
  } 

  done = 1; //if we get all the way through without returning, set done = 1;
}

void Registration::increment_block(std::vector<MVs_PL::MVIterator> &current_index, int increment_amount, int &done, Registration::ImageType::RegionType::IndexType boundary)
{
  done = 0;
  ImageType::RegionType::IndexType temp_index;

  for(int i = 0; i < Dimension; ++i)
  {
    if ((current_index[0].GetIndex()[i] + increment_amount) < boundary[i])
	{

	  temp_index = current_index[0].GetIndex();
	  current_index[0].GoToBegin();
	  for(int j = i-1; j > -1; --j)
	    temp_index[j] = current_index[0].GetIndex()[j];
	  temp_index[i] += increment_amount;
	  current_index[0].SetIndex(temp_index);
	  //Assign the rest of the current_index elements (corresponding to y, z, etc.) to the same as current_index[0]
      for(int i = 1; i < Dimension; ++i)
        current_index[i].SetIndex(current_index[i-1].GetIndex());
	  return;

	}
  }  

  done = 1; //if we get all the way through without returning, set done = 1;
}

void Registration::increment_MVs(std::vector<MVs_PL::MVIterator>& current_index, int increment_amount, int &done, Registration::ImageType::RegionType::IndexType boundary) //set the beginning index for all of the MV components -- x,y,z
{ 
  done = 0;
  ImageType::RegionType::IndexType temp_index;

  for(int i = 0; i < Dimension; ++i)
  {
    if ((current_index[0].GetIndex()[i] + increment_amount) < boundary[i])
	{
	  temp_index = current_index[0].GetIndex();
	  current_index[0].GoToBegin();
	  for(int j = i-1; j > -1; --j)
	    temp_index[j] = current_index[0].GetIndex()[j]; //set all previous elements to zero (if they exist)
	  temp_index[i] += increment_amount; //this will overwrite the element at position zero if there is one.
	  current_index[0].SetIndex(temp_index);
	  //Assign the rest of the current_index elements (corresponding to y, z, etc.) to the same as current_index[0]
      for(int i = 1; i < Dimension; ++i)
        current_index[i].SetIndex(current_index[i-1].GetIndex());
	  return;
	}
  }

  done = 1; //if we get all the way through without returning, set done = 1;

 //The idea is to increment the x direction as long as x doesn't go out of bounds
 //If the x direction goes out of bounds, we try to increment y.  If y goes out of bounds, we increment z

}

void Registration::increment_all(std::vector<MVs_PL::NeighborhoodIteratorType> &current_index, int b_size, int &done, Registration::ImageType::RegionType::IndexType boundary)
{
  for(int i = 0; i < Dimension; ++i)
  {
    increment_neigh_index(current_index[i], b_size, done, boundary);
  }
}

void Registration::increment_allMVs(std::vector<MVs_PL::MVIterator>& MVs)
{
  for(int i = 0; i < Dimension; ++i)
  {
     ++MVs[i];
  }
}

void Registration::volume_difference()
{
  //Create output image
  outputImage = ImageType::New();
  outputImage->SetRegions(inputImage1->GetRequestedRegion());
  outputImage->Allocate(); 

  //Iterator for image 1
  ImageRegionIterator image1(inputImage1, inputImage1->GetRequestedRegion());

  //Iterator for image 2
  ImageRegionIterator image2(inputImage2, inputImage2->GetRequestedRegion());

  //Output iterator
  ImageRegionIterator out(outputImage, inputImage1->GetRequestedRegion());

  for (image1.GoToBegin(); !image1.IsAtEnd(); ++image1, ++image2, ++out)
  {
    out.Set(abs(image1.Get() - image2.Get()));
  }

}

void Registration::draw_MC_volume()
{
  //Create output image
  outputImage = ImageType::New();
  outputImage->SetRegions(pad_pointer2->GetOutput()->GetLargestPossibleRegion());
  outputImage->Allocate(); 
  outputImage->FillBuffer(0);

  int i;

  //Set size of motion-compensated region that will be used to copy pixels from Image 1 equal to N-D block size
  //set sizes of region
  for(i = 0; i < Dimension; ++i)
  {
    block_size[i] = b_size[current_level]; 
  }
  MC_region.SetSize(block_size);
  out_region.SetSize(block_size);
        
  int done = 0; //initlize end value to zero
  ImageType::RegionType::IndexType temp_index;
 
  //Go through all the blocks in Image 2, find the corresponding blocks in Image 1, and copy the pixels from the block in Image 1 to the output image.
  for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, b_size[current_level], done, end_boundary)) //not sure if out++ gets us to the right position, and do we need to check out end bound?
  {
    for(i = 0; i < Dimension; ++i)
	{	
	  //Get start position in Image 1
	  temp_index[i] = MV_PLNum[current_level]->MVs[i].Get(); //get the x, y, and z position in Image 1
	}
	
	MC_region.SetIndex(temp_index); //this will set up the region in Image 1 that we will copy the pixels from
	out_region.SetIndex(MV_PLNum[current_level]->MVs[0].GetIndex()); //set up an output region -- we have to do this so that we can increment with ++ in the loop below (wrap-around)
	
	//Create motion-compensated iterator
    ImageRegionIterator MC_it(pad_pointer1->GetOutput(), MC_region); 
	//Create output region iterator
	ImageRegionIterator out(outputImage, out_region);

	for(MC_it.GoToBegin(), out.GoToBegin(); !MC_it.IsAtEnd(); ++MC_it, ++out)
	{
	  out.Set(MC_it.Get()); //set output pixel to the pixel from Image 1
	}
  } 
}

void Registration::split_MVs(int split_bsize, MVs_PL *& class_element)
{
  //Is this function necessary -- can't we just do this while we are making the MV assignments?
  int done = 0; //reset done variable for first loop
  int done2;
  
  int i;

  ImageType::RegionType::IndexType temp_index;

  //Create temporary space to hold the MV that we will computer
  int *temp_MV = new int[Dimension]; //Need to free this!

  //Create vector to hold the iterators -- probably need to clear the vector at end of loop
  std::vector<MVs_PL::MVIterator> split_iterators; 

  //The first loop goes through the MVs at the previous level of the hierarchy
  for(init_MVs(class_element->MVs); !done; increment_MVs(class_element->MVs, b_size[current_level], done, end_boundary))
  {
	//This loop goes inside one of the blocks in order to split/assign new motion vectors inside each block
	//Need to create new iterator here which goes through each region and assigns the MVs
	//Create the iterator region
	for(i = 0; i < Dimension; ++i)
	{
	  temp_index[i] = class_element->MVs[0].GetIndex()[i];// + split_index[i]; //this is the starting index for the region.  we can skip the first element because it already has a MV
	  temp_MV[i] = class_element->MVs[i].Get() - temp_index[i]; //calculate the MV for the block using position in upper left hand corner
	  block_boundary[i] = class_element->MVs[0].GetIndex()[i] + block_size[i]; //Set end boundaries for the iterator -- these will be needed for the increment_index() function
	}
	split_region.SetIndex(temp_index);
	split_region.SetSize(block_size);
	
	//Set up iterators
	for(i = 0; i < Dimension; ++i)
	{
	  MVs_PL::MVIterator temp(class_element->MV_array[i], split_region); 
	  split_iterators.push_back(temp);
	}

	done2 = 0; //reset done variable
	for(init_MVs_block(split_iterators); !done2; increment_block(split_iterators, split_bsize, done2, block_boundary))
	{
      for(i = 0; i < Dimension; ++i) 
	  {
	    split_iterators[i].Set(split_iterators[0].GetIndex()[i] + temp_MV[i]); //assign the MV of the children inside the block to the parent MV
	  }
	}
	
	//Clear split_iterators vector
	split_iterators.clear();
  }

  delete [] temp_MV;
}

void Registration::print_MVs(int increment)
{
  int done = 0; //reset done value
  int count = 0;
  //All of the MV arrays are the same size, so we just choose the first one
  for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
  {	
	for(int i = 0; i < Dimension; ++i)
	{	
	  //fout << "Index [" << i << "]: " << MV_PLNum[current_level]->MVs[i].Get() << " " << std::endl;
		fout << "MV [" << i << "]: " << (int)(MV_PLNum[current_level]->MVs[i].Get() - MV_PLNum[current_level]->MVs[0].GetIndex()[i]) << " " << std::endl;
	}
	count++;
  } 
  std::cout << "count is " << count << std::endl;
}

void Registration::write_const_MVs()
{
  int done = 0; //reset done value
  int increment = 1;

  for(int i = 0; i < Dimension; ++i)
  {
    end_boundary[i] = current_frame->GetRequestedRegion().GetSize()[i]; //this will get the interpolated but unpadded image -- what we want.
  }
 
  if (frame_num == 1) //left 1
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] + 1);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1]);
    } 

  }
  else if (frame_num == 2) //right 1
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] - 1);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1]);
    } 

  }
  else if (frame_num == 3) //up 1
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0]);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] + 1);
    } 

  }
  else if (frame_num == 4) //down 1
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0]);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] - 1);
    } 

  }
  else if (frame_num == 5) //right 2
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] - 2);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1]);
    } 

  }
  else if (frame_num == 6) //up 2
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0]);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] + 2);	
    } 

  }
  else if (frame_num == 7) //left 1, up 1
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] + 1);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] + 1);
    } 

  }
  else if (frame_num == 8) //left 1, up 2
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] + 1);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] + 2);
    } 

  }
  else if (frame_num == 9) //left 2, up 1
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] + 2);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] + 1);
    } 

  }
  else if (frame_num == 10) //left 3, up 1
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] + 3);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] + 1);
    } 

  }
  else if (frame_num == 11) //left 1, down 1
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] + 1);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] - 1);
    } 

  }
  else if (frame_num == 12) //left 2, down 1
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] + 2);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] - 1);
    } 

  }
  else if (frame_num == 13) //left 2, down 2
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] + 2);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] - 2);
    } 

  }
  else if (frame_num == 14) //right 1, down 1
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] - 1);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] - 1);
    } 

  }
  else if (frame_num == 15) //right 1, down 2
  {
    //All of the MV arrays are the same size, so we just choose the first one
    for(init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
    {	
	  MV_PLNum[current_level]->MVs[0].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[0] - 1);
	  MV_PLNum[current_level]->MVs[1].Set(MV_PLNum[current_level]->MVs[0].GetIndex()[1] - 2);
    } 

  }


}

void Registration::draw_MVs_MATLAB(int increment)
{
	 
  //Set up end boundaries for the iterator we will use below
  for(int i = 0; i < Dimension; ++i)
    end_boundary[i] = (pad_pointer2->GetOutput()->GetLargestPossibleRegion().GetSize()[i]/b_size[current_level])*b_size[current_level];

  //Probably the easiest way to do this is to create three vectors of type int to be passed to Matlab variables
  std::vector<double> u; //positions in frame 1
  std::vector<double> v;
  //std::vector<double> w;

  std::vector<double> x; //positions in frame 2
  std::vector<double> y;
  //std::vector<double> z;
     
  //Save the iterators values into the vectors
  int done = 0;
  for (init_MVs(MV_PLNum[current_level]->MVs); !done; increment_MVs(MV_PLNum[current_level]->MVs, increment, done, end_boundary))
  {
	
    x.push_back((int)MV_PLNum[current_level]->MVs[0].GetIndex()[0]);
	y.push_back((int)MV_PLNum[current_level]->MVs[0].GetIndex()[1]);
	//z.push_back(MV_PLNum[current_level]->MVs[0].GetIndex()[2]);

	u.push_back((int)MV_PLNum[current_level]->MVs[0].GetIndex()[0] - (int)MV_PLNum[current_level]->MVs[0].Get());
    v.push_back((int)MV_PLNum[current_level]->MVs[0].GetIndex()[1] - (int)MV_PLNum[current_level]->MVs[1].Get());
	//w.push_back(MV_PLNum[current_level]->MVs[2].Get());
  }

  Engine *ep; //to hold handle to Matlab engine
  //Initialize Matlab variables
  mxArray *x_var = NULL;
  mxArray *y_var = NULL;
  //mxArray *z_var = NULL;
  mxArray *u_var = NULL;
  mxArray *v_var = NULL;
  //mxArray *w_var = NULL;
  mxArray *result = NULL;

   //* Start the MATLAB engine locally by executing the string
   //* "matlab"
   //*
   //* To start the session on a remote host, use the name of
   //* the host as the string rather than \0
   //*
   //* For more complicated cases, use any string with whitespace,
   //* and that string will be executed literally to start MATLAB
   //*
  if (!(ep = engOpen("\0"))) 
  {
    fprintf(stderr, "\nCan't start MATLAB engine\n");
	getchar();
	exit(1);
 	//return EXIT_FAILURE;
  }

  //Create variables for our data
  //Start with the points in frame 2
  x_var = mxCreateDoubleMatrix(1, x.size(), mxREAL);  //I THINK THESE ARE USING THE WRONG SIZES
  y_var = mxCreateDoubleMatrix(1, y.size(), mxREAL);
  //z_var = mxCreateDoubleMatrix(1, z.size(), mxREAL);

  //Now create variables for the mapped points in frame 1
  u_var = mxCreateDoubleMatrix(1, u.size(), mxREAL);
  v_var = mxCreateDoubleMatrix(1, v.size(), mxREAL);
  //w_var = mxCreateDoubleMatrix(1, w.size(), mxREAL);
  
  //Copy from vectors to Matlab variables
  std::copy(x.begin(), x.end(), mxGetPr(x_var));
  std::copy(y.begin(), y.end(), mxGetPr(y_var));
  //std::copy(z.begin(), z.end(), mxGetPr(z_var));
  std::copy(u.begin(), u.end(), mxGetPr(u_var));
  std::copy(v.begin(), v.end(), mxGetPr(v_var));
  //std::copy(w.begin(), w.end(), mxGetPr(w_var));

  //Place the variables into the MATLAB workspace
  engPutVariable(ep, "x", x_var);
  engPutVariable(ep, "y", y_var);
  //engPutVariable(ep, "z", z_var);
  engPutVariable(ep, "u", u_var);
  engPutVariable(ep, "v", v_var);
  //engPutVariable(ep, "w", w_var);

  //Plot using 3D Quiver
  //engEvalString(ep, "quiver3(x,y,z,u,v,w,1);");
  engEvalString(ep, "quiver(x,y,u,v,2);");
  engEvalString(ep, "set(gca, 'YDir', 'reverse');");
  getchar();

  //We're done! Free memory, close MATLAB figure.
  mxDestroyArray(x_var);
  mxDestroyArray(y_var);
  //mxDestroyArray(z_var);
  mxDestroyArray(u_var);
  mxDestroyArray(v_var);
  //mxDestroyArray(w_var);
  engEvalString(ep, "close;");

  //We're done! Free memory, close MATLAB engine and exit.
  printf("Done!\n");
  engClose(ep);

}

void Registration::write_output(char *output_filename)
{
  //Create writer
  writer = WriterType::New();
  writer->SetFileName(output_filename);
  writer->SetInput(outputImage);

  //Writer settings for DICOM images
  //writer->UseInputMetaDataDictionaryOff(); //this makes sure we use the dictionary from gdcmImageIO and not the input
  //writer->SetImageIO(gdcmImageIO);

  //Try to update the writer or catch error.
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception in file writer " << std::endl;
    std::cerr << err << std::endl;
	getchar();
    exit(1);
  }
}

int compare(const void * a, const void * b)
{
  return (*(int*)a - *(int*)b);
}

int Registration::min(int x, int y)
{
  return !(y < x) ? x : y;
}

int Registration::max(int x, int y)
{
  return (x < y) ? y : x; 
}

Registration::~Registration()
{
  //std::cout << "destructor called" << std::endl;
  for(int i = 0; i < levels; ++i)
    delete MV_PLNum[i];

  MV_PLNum.clear(); //delete vector pointer

  delete [] b_size;
  delete [] s_size;
}
