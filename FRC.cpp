#include "FRC.h"

FRC::FRC(int validity_thres)
{
  val_threshold = validity_thres;
}

void FRC::choose_MV(std::vector<MVIterator>& unaliased_iterator, ImageType::Pointer p_frameUA, ImageType::Pointer c_frameUA, std::vector<MVIterator>& aliased_iterator, ImageType::Pointer p_frameA, ImageType::Pointer c_frameA, int levels, int threshold, char *output_filename)
{
  int i, total_cost;
  int cost[Dimension]; //this holds the difference between the MVs of the aliased and unaliased frames

  //To hold temp index while updating position
  ImageType::RegionType::IndexType temp_aliased;
  ImageType::RegionType::IndexType temp_unaliased;

  //To hold indices which are the positions in the current and previous frame that will be used to get pixels for the middle frame
  ImageType::RegionType::IndexType frame1_index;
  ImageType::RegionType::IndexType frame2_index;

  //Create temp iterators that will be used to access the pixels in previous and current frames of the unaliased case
  typedef itk::ImageRegionIterator<ImageType> ImageRegionIterator;
  ImageRegionIterator frame1UA_it(p_frameUA, p_frameUA->GetLargestPossibleRegion());
  ImageRegionIterator frame2UA_it(c_frameUA, c_frameUA->GetLargestPossibleRegion());

  //Create temp iterators that will be used to access the pixels in previous and current frames of the aliased case
  ImageRegionIterator frame1A_it(p_frameA, p_frameA->GetLargestPossibleRegion());
  ImageRegionIterator frame2A_it(c_frameA, c_frameA->GetLargestPossibleRegion());

  //Set up end boundary to make sure that we don't go outside of image
  for(i = 0; i < Dimension; ++i)
    end_boundary[i] = p_frameA->GetLargestPossibleRegion().GetSize()[i]; //to stop from going outside of boundaries in image 2

  //Create output image for the middle frame
  outputImage = ImageType::New();
  outputImage->SetRegions(p_frameA->GetLargestPossibleRegion());
  outputImage->Allocate(); 
  outputImage->FillBuffer(0);
  //Iterator for output image = middle image
  ImageRegionIterator out_it(outputImage, outputImage->GetLargestPossibleRegion());

  //need to improve the line below to handle more than 2 dimensions
  for(unaliased_iterator[0].GoToBegin(), unaliased_iterator[1].GoToBegin(), aliased_iterator[0].GoToBegin(), aliased_iterator[1].GoToBegin(), out_it.GoToBegin(); !unaliased_iterator[0].IsAtEnd(); ++unaliased_iterator[0], ++unaliased_iterator[1], ++aliased_iterator[0], ++aliased_iterator[1], ++out_it)
  {
	total_cost = 0;
    for(i = 0; i < Dimension; ++i)
	{
	  temp_aliased[i] = aliased_iterator[i].Get(); //this gets the index in previous frame that was stored for the MV for the aliased case
	  temp_unaliased[i] = unaliased_iterator[i].Get(); //this gets the index in previous frame that was stored for the MV for the unaliased case
	  cost[i] = abs(temp_aliased[i] - temp_unaliased[i]); //cost of motion vector in each direction
	  total_cost += cost[i]; //accumulated cost
	}

	if (total_cost < threshold)// && total_cost != 0)
	{
      for(i = 0; i < Dimension; ++i)
	  { //TO DO!! NOTE that the dividing by 2 we are doing below is not the optimal way.  We can either use round or interpolate from neighboring pixels
        frame1_index[i] = max(0, min(end_boundary[i], (temp_aliased[i] + aliased_iterator[0].GetIndex()[i])/2)); //this is the position in frame 1 (previous frame) that will be used for middle frame
	    frame2_index[i] = max(0, min(end_boundary[i], (3*aliased_iterator[0].GetIndex()[i] - temp_aliased[i])/2)); //this is the position in frame 2 (current frame) that will be used for middle frame
	  }
	  frame1A_it.SetIndex(frame1_index); //set the indices so we can acess the pixel values
	  frame2A_it.SetIndex(frame2_index);

	  //Set output pixel to be the avearge of the two pixels in frame 1 and frame 2
	  out_it.Set((frame1A_it.Get() + frame2A_it.Get())/2);
	}
	else
	{
	  for(i = 0; i < Dimension; ++i)
	  {
        frame1_index[i] = max(0, min(end_boundary[i], (temp_unaliased[i] + unaliased_iterator[0].GetIndex()[i])/2))/2; //this is the position in frame 1 (previous frame) that will be used for middle frame
	    frame2_index[i] = max(0, min(end_boundary[i], (3*unaliased_iterator[0].GetIndex()[i] - temp_unaliased[i])/2))/2; //this is the position in frame 2 (current frame) that will be used for middle frame
	  }
	  frame1UA_it.SetIndex(frame1_index); //set the indices so we can acess the pixel values
	  frame2UA_it.SetIndex(frame2_index);

	  //Set output pixel to be the avearge of the two pixels in frame 1 and frame 2
	  out_it.Set((frame1UA_it.Get() + frame2UA_it.Get())/2);
	}	
  }

  //Scale output image to be between 0 and 1
  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(outputImage);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(65535);

  //Write output frame
  //Create writer for output file
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer1;
  writer1 = WriterType::New();
  writer1->SetFileName(output_filename);
  writer1->SetInput(rescaleFilter->GetOutput());
    
  //Try to update the writer or catch error.
  try
  {
    writer1->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception in file writer " << std::endl;
    std::cerr << err << std::endl;
	getchar();
    exit(1);
  }

}

void FRC::perform_FRC(std::vector<MVIterator>& currentMV_iterator, ImageType::Pointer p_frame, ImageType::Pointer c_frame, Validity::Pointer valid1, char *output_filename)
{  
  //Create iterators to go through pixels in the current frame
  typedef itk::ImageRegionIterator<ImageType> ImageRegionIterator;
  ImageRegionIterator current_it(c_frame, c_frame->GetLargestPossibleRegion());

  //Create temp iterators that will be used to access the pixels in previous and current frames
  ImageRegionIterator frame1_it(p_frame, p_frame->GetLargestPossibleRegion());
  ImageRegionIterator frame2_it(c_frame, c_frame->GetLargestPossibleRegion());

  //Create output image for the middle frame
  outputImage = ImageType::New();
  outputImage->SetRegions(c_frame->GetLargestPossibleRegion());
  outputImage->Allocate(); 
  outputImage->FillBuffer(0);
  //Iterator for output image = middle image
  ImageRegionIterator out_it(outputImage, outputImage->GetLargestPossibleRegion());

  //To hold temp index while updating position
  ImageType::RegionType::IndexType temp_index;
  
  //To hold indices which are the positions in the current and previous frame that will be used to get pixels for the middle frame
  ImageType::RegionType::IndexType frame1_index;
  ImageType::RegionType::IndexType frame2_index;
  
  //Set up the validity iterator
  ValIterator val_it(valid1, valid1->GetLargestPossibleRegion());

  for(current_it.GoToBegin(), val_it.GoToBegin(), out_it.GoToBegin(); !current_it.IsAtEnd(); ++current_it, ++out_it, ++val_it)
  {

    //We need to add code right here to check the validity value -- if the validity is small, we need to do averaging.
    if (val_it.Get() > val_threshold)
	{
	  continue;
	}

    for(int i = 0; i < Dimension; ++i)
	{
	  currentMV_iterator[i].SetIndex(current_it.GetIndex()); //this will set up the x, y, and z positions in the MV iterator
	  temp_index[i] = currentMV_iterator[i].Get(); //this gets the index in previous frame that was stored for the MV
	  frame1_index[i] = (temp_index[i] + current_it.GetIndex()[i])/2; //this is the position in frame 1 (previous frame) that will be used for middle frame
	  frame2_index[i] = (3*current_it.GetIndex()[i] - temp_index[i])/2; //this is the position in frame 2 (current frame) that will be used for middle frame
	}

	frame1_it.SetIndex(frame1_index); //set the indices so we can acess the pixel values
	frame2_it.SetIndex(frame2_index);

	//Set output pixel to be the avearge of the two pixels in frame 1 and frame 2
	out_it.Set((frame1_it.Get() + frame2_it.Get())/2);

  }
  
  //Write output frame
  //Create writer for output file
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer;
  writer = WriterType::New();
  writer->SetFileName(output_filename);
  writer->SetInput(outputImage);
    
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

void FRC::frame_difference(char *in_filename, char *out_filename)
{
  //This function will be used to subtract the actual middle from from the middle frame that we estimated using FRC

  //Open the input file
  //Reader for input file
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader;
  reader = ReaderType::New(); //reader for input image
 
  //Pass file names to readers
  reader->SetFileName(in_filename); 
 
  //Try to update the readers or catch the failure
  try
  {
    reader->Update();
	inputImage = reader->GetOutput();
  }
  catch(itk::ExceptionObject& err)
  {
	std::cerr << "Error in one of file readers" << std::endl;
    std::cerr << err << std::endl;
	exit(1);
  }

  //Create iterator for input file
  typedef itk::ImageRegionIterator<ImageType> ImageRegionIterator;
  ImageRegionIterator input_it(inputImage, inputImage->GetLargestPossibleRegion());

  //Create iterator for FRC estimated frame
  ImageRegionIterator estimated_it(outputImage, outputImage->GetLargestPossibleRegion());

  //Create difference frame
  FRC_diff = ImageType::New();
  FRC_diff->SetRegions(outputImage->GetLargestPossibleRegion());
  FRC_diff->Allocate(); 
  FRC_diff->FillBuffer(0);
  //Create iterator for the difference frame
  ImageRegionIterator diff_it(FRC_diff, FRC_diff->GetLargestPossibleRegion());

  for(input_it.GoToBegin(), estimated_it.GoToBegin(), diff_it.GoToBegin(); !input_it.IsAtEnd(); ++input_it, ++estimated_it, ++diff_it)
  {
	  diff_it.Set(abs(input_it.Get() - estimated_it.Get()));
  }

  //Write difference frame
  //Create writer for output file
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer;
  writer = WriterType::New();
  writer->SetFileName(out_filename);
  writer->SetInput(FRC_diff);
    
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

int FRC::min(int x, int y)
{
  return !(y < x) ? x : y;
}

int FRC::max(int x, int y)
{
  return (x < y) ? y : x; 
}