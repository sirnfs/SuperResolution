#include "stdafx.h"
#include "registration.h"
#include "FRC.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *[], int nrhs, const mxArray *prhs[])
{	

  char* vin1 = mxArrayToString(prhs[0]);
  char* vin2 = mxArrayToString(prhs[1]);
  char* orig = mxArrayToString(prhs[2]);
    
  //clock_t start, end; //Variables for keeping track of the time elapsed
  //const int levels = 4; //Number of levels in pyramid
  //int block_size[levels] = {32, 32, 16, 8}; //arranged from lowest resolution level to highest resolution level
  //int search_size[levels] = {12, 16, 24, 12}; //arranged from lowest resolution level to highest resolution level

  const int levels = 3; //Number of levels in pyramid
  int block_size[levels] = {8, 8, 8}; //arranged from lowest resolution level to highest resolution level
  int search_size[levels] = {10, 4, 4}; //arranged from lowest resolution level to highest resolution level

  int interp_factor = 1; //Factor to multiple dimensions of original volume by
  int validity_threshold = 30000;//30000;//12;////30;//80;//100;//2560; (pixel-wise SAD of 5 + 1 = 6 * (overlap of 2) = 12
  int max_intensity = 255;

  //Create classes for each pair of frames that we will compute the motion of
  //NOTE:  PERFORMANCE MODIFICATION NEEDED -- WE REPEAT THE INTERPOLATION OF THE ANCHOR FRAME MULTIPLE TIMES IN THE REGISTRATION CLASS - ONLY NEED TO INTERPOLATE ONCE
  Registration* pair = new Registration(vin1, vin2, block_size, search_size, levels, interp_factor, 1); 
  //start = clock(); //Start timer

  //Perform motion estimation on the frame pair
  pair->motion_estimation();

  //These next two lines are optional and were only used to test an idea to prevent aliasing in undersampled MRI images
  //Need to scale MVs for an image of a different size -- the original size is the unaliased image and the new size will be the size of the aliased image
  int scale_factor = 2;
  pair->scale_MVs(scale_factor);
  
  //Set up frame rate conversion class
  FRC* frame = new FRC(validity_threshold);
   
  //Estimate middle frame using frame rate conversione
  frame->perform_FRC(pair->MV_PLNum[levels-1]->MVs, pair->previous_frame, pair->current_frame, pair->valid, "./FRC_data/out.png");
  frame->frame_difference(orig, "./FRC_data/diff.png");
  
  //end = clock(); //Stop timer
  //double time_elapsed = double(end - start)/CLOCKS_PER_SEC; //time elapsed
  //std::cout << "Time elapsed is " << time_elapsed << std::endl; //output time elapsed

  //std::cout << "Done!" << std::endl;
  //getchar();

  //return 0;
}
