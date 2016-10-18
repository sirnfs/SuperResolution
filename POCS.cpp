#include "POCS.h"

POCS::POCS(char *anchorframe_file, int interp_factor, double d_iter, double variance, int blur_sz, int blur_scale, int v_threshold, int m_intensity)
{
  //Open output file which will just be used for diagnostic purposes
  fout.open("diagnostic.txt");

  //Save interpolation factor as a variable in the class
  i_factor = interp_factor;

  //Save validity threshold as a variable in the class
  validity_threshold = v_threshold;

  //Set the maximum intensity as a variable in the class
  max_intensity = m_intensity;

  //Create reader to read input files
  reader = ReaderType::New(); //reader for input image  
  
  //Pass file names to readers
  reader->SetFileName(anchorframe_file); 
 
  //Try to update the readers or catch the failure
  try
  {
    reader->Update();
	HR_image = reader->GetOutput();
  }
  catch(itk::ExceptionObject& err)
  {
	std::cerr << "Error in reading anchor frame from file" << std::endl;
    std::cerr << err << std::endl;
	getchar();
	exit(1);
  }

  //The code below is used for interpolating the volume
  //Create filter which will be used for interpolation
  FilterType::Pointer filter = FilterType::New();
  //Set up identity transform for the interpolation
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  filter->SetTransform(transform);
  //Create interpolation function -- linear interpolater -- uses next two lines
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  filter->SetInterpolator(interpolator);
  //Create interpolation function -- bspline interpolator -- uses next three lines
  //BSpline_Interpolator::Pointer bs_interpolator = BSpline_Interpolator::New();
  //bs_interpolator->SetSplineOrder(3);
  //filter->SetInterpolator(bs_interpolator);
  //Set up spacing
  ImageType::SpacingType inputSpacing = HR_image->GetSpacing();
  ImageType::SizeType inputSize = HR_image->GetLargestPossibleRegion().GetSize();
  ImageType::SpacingType outputSpacing = inputSpacing;
  ImageType::SizeType outputSize = inputSize; 
  for (int i = 0; i < Dimension; ++i)
  {
    outputSpacing[i] /= i_factor;
    outputSize[i] = static_cast<int>(outputSize[i] * i_factor);
  }
  filter->SetSize(outputSize);
  filter->SetOutputSpacing(outputSpacing);
  filter->SetInput(HR_image); 
  filter->SetDefaultPixelValue(0); //value of pixels outside of the image
  filter->SetOutputOrigin(HR_image->GetOrigin());

  try
  {
    filter->Update();	
	HR_image = filter->GetOutput();
  }
  catch(itk::ExceptionObject& err)
  {
	std::cerr << "Not able to create interpolated anchor frame" << std::endl;
    std::cerr << err << std::endl;
	getchar();
	exit(1);
  }

  //Cast the HR image from unsigned short to type double
  typedef itk::CastImageFilter<ImageType, ImageType_Double> CastFilterType;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(HR_image);

  //typedef itk::RescaleIntensityImageFilter<ImageType, ImageType_Double> RescaleFilterType;
  //RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  //rescaleFilter->SetInput(HR_image);
  //rescaleFilter->SetOutputMinimum(0);
  //rescaleFilter->SetOutputMaximum(65535);

  try
  {
    //rescaleFilter->Update();
	castFilter->Update();
	//HR_image_double = rescaleFilter->GetOutput();
	HR_image_double = castFilter->GetOutput();
  }
  catch(itk::ExceptionObject& err)
  {
	std::cerr << "Not able to convert HR image from unsigned short to double" << std::endl;
    std::cerr << err << std::endl;
	getchar();
	exit(1);
  }

  //Setup POCS parameters
  delta_iter = d_iter;
  blur_size = blur_sz;

  //Fill blur kernel
  calculate_gaussian_coefficients(variance, blur_scale); 

}

double POCS::min(double x, double y)
{
  return !(y < x) ? x : y;
}

double POCS::max(double x, double y)
{
  return (x < y) ? y : x; 
}

void POCS::calculate_gaussian_coefficients(double variance, int blur_scale)
{
  //float pi = 3.14159f; //26535897932384626433832795f;
  //Create Guassian Kernel
  kernel = Gaussian_Generator::New();
  //Set standard deviation, kernel size, and mean
  double sigma[Dimension];
  unsigned long size[Dimension];
  double mean[Dimension];
  for(int i = 0; i < Dimension; ++i)
  {
    sigma[i] = sqrt(variance);
	size[i] = blur_size;
	mean[i] = (double)(blur_size/2); //center the kernel
  }
  kernel->SetSigma(sigma);
  kernel->SetSize(size);
  kernel->SetMean(mean);
  //kernel->SetScale(blur_scale); //Values will go between 0 and blur_scale

  try
  {
    kernel->Update();
  }
  catch(itk::ExceptionObject& err)
  {
	std::cerr << "Could not create Gaussian kernel" << std::endl;
    std::cerr << err << std::endl;
	getchar();
	exit(1);
  }

  double sum = 0;
  //Sum up all of the kernel values
  KernelRegionIterator kernel_it(kernel->GetOutput(), kernel->GetOutput()->GetLargestPossibleRegion());
  for(kernel_it.GoToBegin(); !kernel_it.IsAtEnd(); ++kernel_it)
    sum += kernel_it.Get();

  //Divide all kernel values by the sum
  for(kernel_it.GoToBegin(); !kernel_it.IsAtEnd(); ++kernel_it)
  {
	//std::cout << "Before Value is " << kernel_it.Get() << std::endl;
    kernel_it.Set(kernel_it.Get()/sum);
    //std::cout << "Value is " << kernel_it.Get() << std::endl;
  }

  //Calculate h_squared, which will be needed in POCS
  h_squared = 0;
  for(kernel_it.GoToBegin(); !kernel_it.IsAtEnd(); ++kernel_it)
    h_squared += pow(kernel_it.Get(), 2);

  //std::cout << "hsquared: " << h_squared << std::endl;
}    

void POCS::increment_HRindex(itk::ImageRegionIterator<POCS::ImageType>& HR_index, int shift_amount, itk::Size<POCS::Dimension> shift_size, int& done, itk::Size<POCS::Dimension> boundary)
{
  done = 0;
  //To hold temp index while updating position
  ImageType::RegionType::IndexType temp_index; 

  for(int i = 0; i < Dimension; ++i)
  {
	if ((HR_index.GetIndex()[i] + shift_amount) < boundary[i]) 
	{
	  temp_index = HR_index.GetIndex();
	  HR_index.GoToBegin();
	  //HR_index.SetIndex(HR_index.GetIndex() + shift_size);
	  for(int j = i-1; j > -1; --j)
	    temp_index[j] = HR_index.GetIndex()[j];
	  temp_index[i] += shift_amount;
	  HR_index.SetIndex(temp_index);
	  return;
	}
  } 

  done = 1; //if we get all the way through without returning, set done = 1;
}

void POCS::increment_HRindex2(itk::ImageRegionIterator<POCS::ImageType_Double>& HR_index, int shift_amount, itk::Size<POCS::Dimension> shift_size, int& done, itk::Size<POCS::Dimension> boundary)
{
  done = 0;
  //To hold temp index while updating position
  ImageType_Double::RegionType::IndexType temp_index; 

  for(int i = 0; i < Dimension; ++i)
  {
	if ((HR_index.GetIndex()[i] + shift_amount) < boundary[i]) 
	{
	  temp_index = HR_index.GetIndex();
	  HR_index.GoToBegin();
	  //HR_index.SetIndex(HR_index.GetIndex() + shift_size);
	  for(int j = i-1; j > -1; --j)
	    temp_index[j] = HR_index.GetIndex()[j];
	  temp_index[i] += shift_amount;
	  HR_index.SetIndex(temp_index);
	  return;
	}
  } 

  done = 1; //if we get all the way through without returning, set done = 1;
}

void POCS::increment_VALindex(itk::ImageRegionIterator<POCS::Validity>& HR_index, int shift_amount, itk::Size<POCS::Dimension> shift_size, int& done, itk::Size<POCS::Dimension> boundary)
{
  done = 0;
  //To hold temp index while updating position
  ImageType::RegionType::IndexType temp_index; 

  for(int i = 0; i < Dimension; ++i)
  {
	if ((HR_index.GetIndex()[i] + shift_amount) < boundary[i]) 
	{
	  temp_index = HR_index.GetIndex();
	  HR_index.GoToBegin();
	  //HR_index.SetIndex(HR_index.GetIndex() + shift_size);
	  for(int j = i-1; j > -1; --j)
	    temp_index[j] = HR_index.GetIndex()[j];
	  temp_index[i] += shift_amount;
	  HR_index.SetIndex(temp_index);
	  return;
	}
  } 

  done = 1; //if we get all the way through without returning, set done = 1;
}

void POCS::process_frame(ImageType::Pointer c_frame, std::vector<MVIterator>& currentMV_iterator)
{
  //Create iterators to go through pixels in the current frame and in the HR frame, and through the MVs
  //For the adjacent frame -- the frame other than the anchor
  ImageRegionIterator current_it(c_frame, c_frame->GetLargestPossibleRegion());
  //Iterator for the HR image...works with type double
  ImageRegionIterator_Double HR_it(HR_image_double, HR_image_double->GetLargestPossibleRegion());

  //Create a local region iterator that will be used in the image to multiple the HR estimate by the kernel
  ImageType_Double::RegionType local_region;

  //For the kernel iterator
  typedef itk::ImageRegionIterator<blur_image> KernelRegionIterator;
  KernelRegionIterator kernel_it(kernel->GetOutput(), kernel->GetOutput()->GetLargestPossibleRegion());

  //For POCS 
  double residual = 0; //Will hold error between the pixel of the HR estimate and the pixel of the current frame
  double residual_diff = 0; //difference between residual and the allowed error, delta_iter
  double result = 0; //Holds sum of the variable of type double to hold the value we use for POCS
  unsigned int LR_pixel; //Current LR pixel we are working with

  //We have to skip some lines below because the 3x3 blur kernel will go outside the bounds at the beginning of the volume
  int shift_param = blur_size >> 1;
  itk::Size<Dimension> shift_amount;
  itk::Size<Dimension> kernel_size;
  itk::Size<Dimension> interp_shift; //how much to shift because of interpolated pixels
  itk::Index<Dimension> zero_index;
  for(int i = 0; i < Dimension; ++i)
  {
    shift_amount[i] = shift_param;
	kernel_size[i] = blur_size;
	zero_index[i] = 0;
	interp_shift[i] = i_factor; //We will skip any of the interpolated pixels
  }
   //To hold temp index while updating position
  ImageType::RegionType::IndexType temp_index;

  int done = 0; //initialize done to zero
  int skip_count = 0; //To keep track of how many pixels were not included in the HR image by POCS
  int total_count = 0; //Total number of pixels in all volumes that were processed, but not necessarily included in the HR iamge

  //set the boundaries of the volume
  itk::Size<Dimension> boundary;
  for(int i = 0; i < Dimension; ++i)
  {
    boundary[i] = HR_image_double->GetLargestPossibleRegion().GetSize()[i];
  }
  
  HR_it.GoToBegin(); //set iterator to beginning
  HR_it.SetIndex(HR_it.GetIndex() + shift_amount);  //This code skips the initial lines at the start of the volume volume
  
  for(; !done; increment_HRindex2(HR_it, i_factor, interp_shift, done, boundary))
  {
	//fout << "index is " << HR_it.GetIndex() << std::endl;

	//Set parameters for local region iterator
	local_region.SetIndex(HR_it.GetIndex() - shift_amount);
    local_region.SetSize(kernel_size);
    ImageRegionIterator_Double local_it(HR_image_double, local_region); //local_it uses the double data type

	result = 0;
	//Multiply the pixels in the HR estimate by the Gaussian kernel -- local_it is of data type double
    for(kernel_it.GoToBegin(), local_it.GoToBegin(); !kernel_it.IsAtEnd(); ++kernel_it, ++local_it)
	{
      result += kernel_it.Get()*local_it.Get();
	}

	//Get the LR pixel that corresponds to the HR pixel (found from the MV)
	for(int i = 0; i < Dimension; ++i)
	{
	  currentMV_iterator[i].SetIndex(HR_it.GetIndex()); //this will set up the x, y, and z positions in the MV iterator
	  temp_index[i] = currentMV_iterator[i].Get();
	  if((temp_index[i] < 0) || (temp_index[i] >= boundary[i]))
	    goto skip_pixel;
	}
	current_it.SetIndex(temp_index); //set the index of the current frame to the pixel index obtained from motion positions -- currentMV_iterator has positions and not MVs
	LR_pixel  = current_it.Get(); //get the current LR pixel
	
	//Subtract the current LR pixel from the result which was obtained by multiplying the HR pixels by the Gaussian kernel
	residual = (double)LR_pixel - result; 
	int count = 0; //to keep track of the number of iterations

	while((residual > delta_iter) || (residual < (-delta_iter)))
    {
	  count++;
      if (count > 1000) //in case sets do not converge
	  {  
		skip_count++;
	    //std::cout << count << std::endl;
	    break;  
	  }

      if (residual > delta_iter)
      {
        residual_diff = residual - delta_iter;
      }
      else //if (residual < delta_iter)
      {
        residual_diff = residual + delta_iter; 
      }

	  //Multiply the residual difference by the Gaussian kernel
      for(local_it.GoToBegin(), kernel_it.GoToBegin(); !kernel_it.IsAtEnd(); ++kernel_it, ++local_it)
	  {
	    local_it.Set(local_it.Get() + (residual_diff*kernel_it.Get())/h_squared); 
	  }

	  //Are we done with POCS on this pixel, or do we need to do more iterations to reduce the error?
	  //Let's find out...
	  //Initialize result back to zero
	  result = 0;
	  //Multiply the new pixel values in the HR estimate by the Gaussian kernel
      for(kernel_it.GoToBegin(), local_it.GoToBegin(); !kernel_it.IsAtEnd(); ++kernel_it, ++local_it)
	  {
        result += kernel_it.Get()*local_it.Get();
	  }
	
	  residual = (double)LR_pixel - result; 

    }

	//HR_it.Set(max(0, min(HR_it.Get(), std::numeric_limits<double>::max())));  //clamp the values to fall between 0 and 65535

	skip_pixel:  result = 0; //the result = 0 is really only here so i can have a label attached to a statement
  }

  //Apply amplitude constraints
  ImageRegionIterator_Double amplitude_it(HR_image_double, HR_image_double->GetLargestPossibleRegion()); //amplitude_it uses the double data type
  for(amplitude_it.GoToBegin(); !amplitude_it.IsAtEnd(); ++amplitude_it)
  {
	 amplitude_it.Set(max(0, min(amplitude_it.Get(), max_intensity))); //clamp the pixel values to be between 0 and 65535.
  }

  std::cout << "Total count is " << total_count << std::endl;
  std::cout << "Skip count is " << skip_count << std::endl;

}

void POCS::process_frame_reverse(ImageType::Pointer c_frame, std::vector<MVIterator>& currentMV_iterator, Validity::Pointer valid1)
{
  //Create iterators to go through pixels in the current frame and in the HR frame, and through the MVs
  //For the adjacent frame -- the frame other than the anchor
  ImageRegionIterator current_it(c_frame, c_frame->GetLargestPossibleRegion());
  //Iterator for the HR image...works with type double
  ImageRegionIterator_Double HR_it(HR_image_double, HR_image_double->GetLargestPossibleRegion());

  //Create a local region iterator that will be used in the image to multiple the HR estimate by the kernel
  ImageType_Double::RegionType local_region;

  //For the kernel iterator
  typedef itk::ImageRegionIterator<blur_image> KernelRegionIterator;
  KernelRegionIterator kernel_it(kernel->GetOutput(), kernel->GetOutput()->GetLargestPossibleRegion());

  //Set up the validity iterator
  ValIterator val_it(valid1, valid1->GetLargestPossibleRegion());

  //For POCS 
  double residual = 0; //Will hold error between the pixel of the HR estimate and the pixel of the current frame
  double residual_diff = 0; //difference between residual and the allowed error, delta_iter
  double result = 0; //Holds sum of the variable of type double to hold the value we use for POCS
  unsigned int LR_pixel; //Current LR pixel we are working with

  //We have to skip some lines below because the 3x3 blur kernel will go outside the bounds at the beginning of the volume
  int shift_param = blur_size >> 1;
  itk::Size<Dimension> shift_amount;
  itk::Size<Dimension> kernel_size;
  itk::Size<Dimension> interp_shift; //how much to shift because of interpolated pixels
  itk::Index<Dimension> zero_index;
  for(int i = 0; i < Dimension; ++i)
  {
    shift_amount[i] = shift_param;
	kernel_size[i] = blur_size;
	zero_index[i] = 0;
	interp_shift[i] = i_factor; //We will skip any of the interpolated pixels
  }

  //To hold temp index while updating position
  ImageType::RegionType::IndexType temp_index;

  //Set the boundaries of the volume
  itk::Size<Dimension> boundary;
  for(int i = 0; i < Dimension; ++i)
  {
    boundary[i] = HR_image->GetLargestPossibleRegion().GetSize()[i];
  }
  
  int skip_count = 0; //To keep track of how many pixels were not included in the HR image by POCS
  int total_count = 0; //Total number of pixels in all volumes that were processed, but not necessarily included in the HR iamge

  int done = 0; //indicates when we have finished processing the frame
  int done2 = 0; //not used, but passed.

  //Below, we will only include the LR pixels in the POCS framework
  //for(current_it.GoToBegin(), val_it.GoToBegin(); !current_it.IsAtEnd(); ++current_it, ++val_it) //go through all LR pixels and project onto HR grid
  for(current_it.GoToBegin(), val_it.GoToBegin(); !done; increment_HRindex(current_it, i_factor, interp_shift, done, boundary), increment_VALindex(val_it, i_factor, interp_shift, done2, boundary))
  {	
	total_count++;
	if (val_it.Get() > validity_threshold)
	{
	  skip_count++;
	  continue;
	}

	//Get the HR pixel that corresponds to the LR pixel (found from the MV)
	for(int i = 0; i < Dimension; ++i)
	{
	  currentMV_iterator[i].SetIndex(current_it.GetIndex()); //this will set up the x, y, and z positions in the MV iterator
	  temp_index[i] = currentMV_iterator[i].Get(); //this gets the index in frame 1 that was stored for the MV
	  //Figure out if the temp_index is inside the boundary of the HR image
	  //The code below basically just sets zero_index to some known index value so that the code after the for loop can run.
	  //The problem is that the 'goto skip_pixel' code has to come after the "ImageRegionIterator_Double local_it(HR_image_double, local_region)" line, or the compile will give an error.
	  if(((temp_index[i] - shift_param) < 0) || ((temp_index[i] + shift_param) >= boundary[i]))
		  zero_index[i] = shift_param;
	  else
		  zero_index[i] = temp_index[i];
	}

	//Set parameters for local region iterator
	HR_it.SetIndex(zero_index);
	local_region.SetIndex(HR_it.GetIndex() - shift_amount);
    local_region.SetSize(kernel_size);
    ImageRegionIterator_Double local_it(HR_image_double, local_region); //local_it uses the double data type

	//The code below checks to make sure that we stay inside the volume
	for(int i = 0; i < Dimension; ++i)
	{
   	  if(((temp_index[i] - shift_param) < 0) || ((temp_index[i] + shift_param) >= boundary[i]))
	  {
		skip_count++;
	    goto skip_pixel;
	  }
	}

	result = 0;
	//Multiply the pixels in the HR estimate by the Gaussian kernel -- local_it is of data type double
    for(kernel_it.GoToBegin(), local_it.GoToBegin(); !kernel_it.IsAtEnd(); ++kernel_it, ++local_it)
	{
      result += kernel_it.Get()*local_it.Get();
	}
		
	//Subtract the current LR pixel from the result which was obtained by multiplying the HR pixels by the Gaussian kernel
	LR_pixel  = current_it.Get(); //get the current LR pixel
	residual = (double)LR_pixel - result; 
	int count = 0; //to keep track of the number of iterations

	while((residual > delta_iter) || (residual < (-delta_iter)))
    {
	  count++;
      if (count > 1000) //in case sets do not converge
	  {  
		skip_count++;
	    //std::cout << count << std::endl;
	    break;  
	  }

      if (residual > delta_iter)
      {
        residual_diff = residual - delta_iter;
      }
      else //if (residual < delta_iter) -- this comment is here just to show what the other case is
      {
        residual_diff = residual + delta_iter; 
      }

	  //Multiply the residual difference by the Gaussian kernel
      for(local_it.GoToBegin(), kernel_it.GoToBegin(); !kernel_it.IsAtEnd(); ++kernel_it, ++local_it)
	  {
	    local_it.Set(local_it.Get() + (residual_diff*kernel_it.Get())/h_squared); 
	    //local_it.Set(max(0, min((local_it.Get() + (residual_diff*kernel_it.Get())/h_squared), max_intensity))); //clamp the pixel values to be between 0 and 65535.
	  }

	  //Are we done with POCS on this pixel, or do we need to do more iterations to reduce the error?
	  //Let's find out...
	  //Initialize result back to zero
	  result = 0;
	  //Multiply the new pixel values in the HR estimate by the Gaussian kernel
      for(kernel_it.GoToBegin(), local_it.GoToBegin(); !kernel_it.IsAtEnd(); ++kernel_it, ++local_it)
	  {
        result += kernel_it.Get()*local_it.Get();
	  }
	
	  residual = (double)LR_pixel - result; 

    }
	
	skip_pixel:  result = 0; //the result = 0 is really only here so i can have a label attached to a statement -- it serves no other purpose
  }

  //Apply amplitude constraints
  ImageRegionIterator_Double amplitude_it(HR_image_double, HR_image_double->GetLargestPossibleRegion()); //amplitude_it uses the double data type
  for(amplitude_it.GoToBegin(); !amplitude_it.IsAtEnd(); ++amplitude_it)
  {
	 amplitude_it.Set(max(0, min(amplitude_it.Get(), max_intensity))); //clamp the pixel values to be between 0 and 65535.
  }

  std::cout << "Total count is " << total_count << std::endl;
  std::cout << "Skip count is " << skip_count << std::endl;
}

void POCS::write_output(char *output_filename)
{
  //At the end of all of this, we convert the HR_image from double back to unsigned short
  /*typedef itk::RescaleIntensityImageFilter<ImageType_Double, ImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(HR_image_double);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(24400);

  try
  {
    rescaleFilter->Update();
	HR_image = rescaleFilter->GetOutput();
  }
  catch(itk::ExceptionObject& err)
  {
	std::cerr << "Not able to rescale HR image from unsigned short to double" << std::endl;
    std::cerr << err << std::endl;
	getchar();
	exit(1);
  }*/

  //Cast the HR image from unsigned short to type double
  typedef itk::CastImageFilter<ImageType_Double, ImageType> CastFilterType;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(HR_image_double);

  try
  {
    castFilter->Update();
	HR_image = castFilter->GetOutput();
  }
  catch(itk::ExceptionObject& err)
  {
	std::cerr << "Not able to convert HR image from double to unsigned short" << std::endl;
    std::cerr << err << std::endl;
	getchar();
	exit(1);
  }

  //Create writer
  writer = WriterType::New();
  writer->SetFileName(output_filename);
  writer->SetInput(HR_image); //may need to cast this?

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

POCS::~POCS()
{

}