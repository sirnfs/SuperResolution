#include "stdafx.h"
#include "pair_motion.h"

Pair_Motion::Pair_Motion(int * search_sizes, int * block_sizes)
{  
  //These are the search sizes and block sizes for the entire hierarchy of images
  level4_search_size = search_sizes[3];
  level4_block_size = block_sizes[3];
  level3_search_size = search_sizes[2];
  level3_block_size = block_sizes[2];
  level2_search_size = search_sizes[1]; 
  level2_block_size = block_sizes[1];
  level1_search_size = search_sizes[0]; 
  level1_block_size = block_sizes[0]; 
  start_pos_level4 = (level4_search_size >> 1) - (level4_block_size >> 1);
  start_pos_level3 = (start_pos_level4 << 1); //(level3_search_size >> 1) - (level3_block_size >> 1);
  start_pos_level2 = (start_pos_level3 << 1);
  start_pos_level1 = (start_pos_level2 << 1);  
  //Shifting block inside search window -- needed to prevent exceeding width/height
  v_shift1 = level1_search_size - level1_block_size; 
  v_shift2 = level2_search_size - level2_block_size;
  v_shift3 = level3_search_size - level3_block_size;
  v_shift4 = level4_search_size - level4_block_size;

}
void Pair_Motion::pad_images(IplImage*& imageA, IplImage*& imageB, IplImage*& imageA_pad, IplImage*& imageB_pad, int start_pos)
{		
	//add_height = (block_size - (image1->height - (int)floor(double(image1->height)/block_size)*block_size)) % block_size + ((start_pos << 1) - block_size);
	//add_width = (block_size - (image1->width - (int)floor(double(image1->width)/block_size)*block_size)) % block_size + ((start_pos << 1) - block_size);
	
	//std::cout << "add Height is " << add_height << std::endl;
	//std::cout << "add Width is " << add_width << std::endl;
	imageA_pad = cvCreateImage(cvSize(imageA->width + add_width, imageA->height + add_height), IPL_DEPTH_8U, 1);
	imageB_pad = cvCreateImage(cvSize(imageB->width + add_width, imageB->height + add_height), IPL_DEPTH_8U, 1);
	cvZero(imageA_pad);
	cvZero(imageB_pad);

	for(int i = start_pos; i < (imageA->height+start_pos); i++)
	{
	  for(int j = start_pos; j < (imageA->width+start_pos); j++)
	  {
	    imageA_pad->imageData[i*imageA_pad->widthStep+j] = imageA->imageData[(i-start_pos)*imageA->widthStep+(j-start_pos)];
		//image1_pad->imageData[i*image1_pad->widthStep+j*3+1] = image1->imageData[(i-start_pos)*image1->widthStep+(j-start_pos)*3+1];
		//image1_pad->imageData[i*image1_pad->widthStep+j*3+2] = image1->imageData[(i-start_pos)*image1->widthStep+(j-start_pos)*3+2];
		imageB_pad->imageData[i*imageB_pad->widthStep+j] = imageB->imageData[(i-start_pos)*imageB->widthStep+(j-start_pos)];	
		//image2_pad->imageData[i*image2_pad->widthStep+j*3+1] = image2->imageData[(i-start_pos)*image2->widthStep+(j-start_pos)*3+1];	
		//image2_pad->imageData[i*image2_pad->widthStep+j*3+2] = image2->imageData[(i-start_pos)*image2->widthStep+(j-start_pos)*3+2];	
	  }
	}
}
void Pair_Motion::pad_image_color(IplImage*& imageA, IplImage*& imageA_pad, int start_pos)
{			
	
  imageA_pad = cvCreateImage(cvSize(imageA->width + add_width, imageA->height + add_height), IPL_DEPTH_8U, 3);
	
  cvZero(imageA_pad);	

  for(int i = start_pos; i < (imageA->height+start_pos); i++)
  {
    for(int j = start_pos; j < (imageA->width+start_pos); j++)
	{
	  imageA_pad->imageData[i*imageA_pad->widthStep+j*3] = imageA->imageData[(i-start_pos)*imageA->widthStep+(j-start_pos)*3];
	  imageA_pad->imageData[i*imageA_pad->widthStep+j*3+1] = imageA->imageData[(i-start_pos)*imageA->widthStep+(j-start_pos)*3+1];
	  imageA_pad->imageData[i*imageA_pad->widthStep+j*3+2] = imageA->imageData[(i-start_pos)*imageA->widthStep+(j-start_pos)*3+2];
			
	}
  }
}
void Pair_Motion::start_calculation_nointerp(IplImage*& image1, IplImage*& image2, int lambda_value)
{           
  // Load input images
  Computed_Data.frame1 = image1;//cvCloneImage(image1); 
  frame2 = image2;//cvCloneImage(image2); 
  level1a = Computed_Data.frame1; 
  level1b = frame2;  
   
  //Convert color images to grayscale if necessary (Level1a)
  if (level1a->nChannels == 3)
  {
	level1a_grey = cvCreateImage(cvGetSize(level1a), IPL_DEPTH_8U, 1);
	// Convert to grayscale image
	cvCvtColor(level1a, level1a_grey, CV_RGB2GRAY);		
  }	
  else
  {
    //Means we already have a grayscale image
	level1a_grey = cvCreateImage(cvGetSize(level1a), IPL_DEPTH_8U, 1);
    cvCopyImage(level1a, level1a_grey); 
  }

  //Convert color images to grayscale if necessary (Level1a)
  if (level1b->nChannels == 3)
  {
	level1b_grey = cvCreateImage(cvGetSize(level1b), IPL_DEPTH_8U, 1);
	// Convert to grayscale image
	cvCvtColor(level1b, level1b_grey, CV_RGB2GRAY);		
  }	
  else
  {
    //Means we already have a grayscale image
	level1b_grey = cvCreateImage(cvGetSize(level1b), IPL_DEPTH_8U, 1);
    cvCopyImage(level1b, level1b_grey); 
  }  
         	
  //Level 2 Downsampling
  level2a_grey = cvCreateImage(cvSize(level1a_grey->width/2,level1a_grey->height/2), IPL_DEPTH_8U, 1);
  level2b_grey = cvCreateImage(cvSize(level1a_grey->width/2,level1a_grey->height/2), IPL_DEPTH_8U, 1);	
  cvPyrDown(level1a_grey, level2a_grey);
  cvPyrDown(level1b_grey, level2b_grey);   
  
  //Level 3 Downsampling 
  level3a_grey = cvCreateImage(cvSize(level2a_grey->width/2,level2a_grey->height/2), IPL_DEPTH_8U, 1);
  level3b_grey = cvCreateImage(cvSize(level2b_grey->width/2,level2b_grey->height/2), IPL_DEPTH_8U, 1); 
  cvPyrDown(level2a_grey, level3a_grey);
  cvPyrDown(level2b_grey, level3b_grey);

  //Level 4 Downsampling
  level4a_grey = cvCreateImage(cvSize(level3a_grey->width/2,level3a_grey->height/2), IPL_DEPTH_8U, 1);
  level4b_grey = cvCreateImage(cvSize(level3b_grey->width/2,level3b_grey->height/2), IPL_DEPTH_8U, 1); 
  cvPyrDown(level3a_grey, level4a_grey);
  cvPyrDown(level3b_grey, level4b_grey);

  //Pad Images
  //Level 4
  add_height = (level4a_grey->height - (int)floor(double(level4a_grey->height)/level4_block_size)*level4_block_size) % level4_block_size + ((start_pos_level4 << 1));
  add_width =  (level4a_grey->width - (int)floor(double(level4a_grey->width)/level4_block_size)*level4_block_size) % level4_block_size + ((start_pos_level4 << 1));
  add_height4 = add_height;
  add_width4 = add_width;
  pad_images(level4a_grey, level4b_grey, level4a_pad, level4b_pad, start_pos_level4);
  //Don't need the original images anymore, so let's delete them
  cvReleaseImage(&level4a_grey);
  cvReleaseImage(&level4b_grey); 
  //Assign data to padded images
  data_level4a = (uchar *)level4a_pad->imageData;
  data_level4b = (uchar *)level4b_pad->imageData;	
  height4 = level4a_pad->height; 
  width4 = level4a_pad->width; 
  step4 = level4a_pad->widthStep;

  //Level 3
  add_height = (add_height << 1);
  add_width = (add_width << 1);
  add_height3 = add_height;
  add_width3 = add_width;
  pad_images(level3a_grey, level3b_grey, level3a_pad, level3b_pad, start_pos_level3);
  //Don't need the original images anymore, so let's delete them
  cvReleaseImage(&level3a_grey);
  cvReleaseImage(&level3b_grey);
  //Assign data to padded images
  data_level3a = (uchar *)level3a_pad->imageData;
  data_level3b = (uchar *)level3b_pad->imageData;	
  height3 = level3a_pad->height; 
  width3 = level3a_pad->width; 
  step3 = level3a_pad->widthStep;

  //Level 2
  add_height = (add_height << 1);
  add_width = (add_width << 1);
  add_height2 = add_height;
  add_width2 = add_width;
  pad_images(level2a_grey, level2b_grey, level2a_pad, level2b_pad, start_pos_level2);
  //Don't need the original images anymore, so let's delete them
  cvReleaseImage(&level2a_grey);
  cvReleaseImage(&level2b_grey);
  //Assign data to padded images
  data_level2a = (uchar *)level2a_pad->imageData;
  data_level2b = (uchar *)level2b_pad->imageData;	
  height2 = level2a_pad->height; 
  width2 = level2a_pad->width; 
  step2 = level2a_pad->widthStep;

  //Level 1
  add_height = (add_height << 1);
  add_width = (add_width << 1);
  add_height1 = add_height;
  add_width1 = add_width;
  pad_images(level1a_grey, level1b_grey, level1a_pad, level1b_pad, start_pos_level1);
  //The next line is for the motion compensated frame, data_color1 below
  if (level1a->nChannels == 3)
  {
    pad_image_color(level1b, level1b_pad_color, start_pos_level1);
    pad_image_color(level1a, level1a_pad_color, start_pos_level1);
	//cvSaveImage("./stereo/level1bpad.png", level1b_pad_color);
	//cvSaveImage("./stereo/level1apad.png", level1a_pad_color);

  }

  width1 = level1a_pad->width; 
  height1 = level1a_pad->height;
  step1 = level1a_pad->widthStep;
  data_level1a = (uchar *)level1a_pad->imageData;
  data_level1b = (uchar *)level1b_pad->imageData;  
  
  //For displaying the motion vectors of the two frames 
  if (level1a->nChannels == 3)
  {
    debug = cvCloneImage(level1a_pad_color); //was level1b
  //debug2 = cvCloneImage(level1b_pad);  
  }
  else
    debug = cvCloneImage(level1a_pad);

  //For the predicted image (assuming highest level here)
  if (level1a->nChannels == 3)    
  {
    frame_predict = cvCreateImage(cvSize(level1a_grey->width, level1a_grey->height), IPL_DEPTH_8U, 3);
    data_color1 = (uchar *)level1a_pad_color->imageData; 
    color_step = level1a_pad_color->widthStep; 
  }
  else    
  { 
    frame_predict = cvCreateImage(cvSize(level1a_grey->width, level1a_grey->height), IPL_DEPTH_8U, 1);
	data_color1 = (uchar *)level1a_pad->imageData; 
    color_step = level1a_pad->widthStep; 
 }

  data_predict = (uchar *)frame_predict->imageData;
  predict_step = frame_predict->widthStep; 

  //Delete level1 greyscale images that we don't need
  cvReleaseImage(&level1a_grey);
  cvReleaseImage(&level1b_grey);
  
  // This is where we allocate our dynamic arrays used to calculate motion positions and motion vectors
  Computed_Data.v_x1 = new int[height1*width1];//step1];
  Computed_Data.v_y1 = new int[height1*width1];//step1];	
  Computed_Data.reliability = new double[height1*width1];//step1];
  overlap1 = new int[height1*width1];   
  overlap2 = new int[height2*width2];  
  overlap3 = new int[height3*width3]; 
  overlap4 = new int[height4*width4]; 
  v_x2 = new int[height2*width2];
  v_y2 = new int[height2*width2]; 
  v_x3 = new int[height3*width3];
  v_y3 = new int[height3*width3];  
  v_x4 = new int[height4*width4];
  v_y4 = new int[height4*width4];  

  //Initialize all values to zero for Level 1
  for(int i = 0; i < height1*width1; i++)
  {
    Computed_Data.v_x1[i] = 0;
	Computed_Data.v_y1[i] = 0;
	Computed_Data.reliability[i] = 0;	
	overlap1[i] = 0;		
  }
  //Initialize all values to zero for level2
  for(int i = 0; i < height2*width2; i++)
  {
    v_x2[i] = 0;
	v_y2[i] = 0;	  
	overlap2[i] = 0;		
  }
  //Initialize all values to zero for level3
  for(int i = 0; i < height3*width3; i++)
  {   
    v_x3[i] = 0;	  
	v_y3[i] = 0;	  	
	overlap3[i] = 0;	   
  }

  //Initialize all values to zero for level4
  for(int i = 0; i < height4*width4; i++)
  {   
    v_x4[i] = 0;	  
	v_y4[i] = 0;	 		  	
	overlap4[i] = 0;		
  }

  calculate_motion_vectors_overlap(lambda_value);

}
void Pair_Motion::calculate_motion_vectors_overlap(int lambda_value)
{
  //Applies to all levels
  int temp_lambda;

  //--------------------Level 4----------------------------
  //Initializations for Level 4
  level = 4;
  width = width4;
  height = height4;
  add_height = add_height4;
  add_width = add_width4;
  step = step4;   
  search_size = level4_search_size;
  b_size = level4_block_size;
  start_pos = start_pos_level4;
  v_shift = v_shift4;
  data_a = data_level4a;
  data_b = data_level4b;
  v_x = v_x4;   
  v_y = v_y4;   
  Computed_Data.overlap = overlap4;
  //End of initializations

  //onelevlspiral_BM_minoverlap(); 
  onelevlspiral_BM();
  
  for(int k = 0; k < 3; k++) //change back to 3
  {
	b_size = level4_block_size;
	temp_lambda = (level4_block_size*3) >> 2;
	while (b_size > 1) 
	{ 
	  for (int l = 0; l < 4; l++) //change back to 4
	  {
		calculate_pixel_reliability();
		add_smoothness8_overlap((l+1)*temp_lambda);
		//add_smoothness8_old((l+1)*temp_lambda);
	  }  
	  setMVs_iter();
	  b_size = (b_size >> 1);
	  temp_lambda = temp_lambda >> 1;
	}	
  } 
  b_size = (level3_block_size >> 1); // set to half of block size needed at next level
  
  //-----------------End of Level 4------------------------
  
  //--------------------Level 3----------------------------
  //Initializations for Level 3
  level = 3;
  width_previous = width4;
  width = width3;
  height_previous = height4;
  height = height3;
  add_height = add_height3;
  add_width = add_width3;
  step = step3;   
  step_previous = step4;
  search_size = level3_search_size;
  b_size_previous = b_size; 
  b_size = level3_block_size;
  start_pos_prev = start_pos_level4;
  start_pos = start_pos_level3;
  v_shift = v_shift3;
  v_shift_previous = v_shift4;
  data_a = data_level3a;
  data_b = data_level3b;
  v_xprev = v_x4;
  v_yprev = v_y4;
  v_x = v_x3;   
  v_y = v_y3;   
  Computed_Data.overlap = overlap3;  
  //End of initializations

  nextlevlspiral_BM(); 
  //nextlevlspiral_BM_minoverlap(1);  
  //nextlevlspiral_BM_weightedoverlap(-1);
  
  for(int k = 0; k < 3; k++) //change back to 3
  {
	b_size = level3_block_size;
	temp_lambda = (level3_block_size*3) >> 2;
	while (b_size > 1) 
	{ 
	  for (int l = 0; l < 4; l++) //change back to 4
	  {
	    calculate_pixel_reliability();
		add_smoothness8_overlap((l+1)*temp_lambda);
		//add_smoothness8_old((l+1)*temp_lambda);
	  }  
	  setMVs_iter();
	  b_size = (b_size >> 1);
	  temp_lambda = temp_lambda >> 1;
	}	
  } 
  b_size = (level2_block_size >> 1); // set to half of block size needed at next level
  
  //-----------------End of Level 3------------------------

  //---------------------Level 2---------------------------
  //Initializations for Level 2 
  level = 2;
  width_previous = width3;
  width = width2;
  height_previous = height3;
  height = height2;
  add_height = add_height2;
  add_width = add_width2;
  step = step2;
  step_previous = step3;
  search_size = level2_search_size;
  b_size_previous = b_size; 
  b_size = level2_block_size;
  start_pos_prev = start_pos_level3;
  start_pos = start_pos_level2;
  v_shift = v_shift2;
  v_shift_previous = v_shift3;
  data_a = data_level2a;
  data_b = data_level2b;
  v_xprev = v_x3;
  v_yprev = v_y3;
  v_x = v_x2;  
  v_y = v_y2;   
  Computed_Data.overlap = overlap2;
  //End of initializations
  
  nextlevlspiral_BM(); 
  //nextlevlspiral_BM_minoverlap(1);
  //nextlevlspiral_BM_weightedoverlap(1);
   
  for(int k = 0; k < 3; k++) //change back to 3
  {
    b_size = level2_block_size;
	temp_lambda = (level2_block_size*3) >> 2;
	while (b_size > 1) 
	{	
	  for (int l = 0; l < 4; l++) //change back to 4
	  {
	    calculate_pixel_reliability();
		add_smoothness8_overlap((l+1)*temp_lambda);
		//add_smoothness8_old((l+1)*temp_lambda);
	  }  	  	 
      setMVs_iter();
	  b_size = (b_size >> 1);
	  temp_lambda = temp_lambda >> 1;
	}
  }  
  b_size = (level1_block_size >> 1); //set to half of the block size needed at next level  
    
  //-----------------End of Level 2------------------------  

  //---------------------Level 1---------------------------
  //Initializations for Level 1  
  level = 1;
  width_previous = width2;
  width = width1;
  height_previous = height2;
  height = height1;
  add_height = add_height1;
  add_width = add_width1;
  step = step1;
  step_previous = step2;
  search_size = level1_search_size;
  b_size_previous = b_size;
  b_size = level1_block_size;
  //start_pos no change at first
  v_shift = v_shift1;
  v_shift_previous = v_shift2;
  data_a = data_level1a;
  data_b = data_level1b;
  v_xprev = v_x2;
  v_yprev = v_y2;
  v_x = Computed_Data.v_x1;
  v_y = Computed_Data.v_y1;
  Computed_Data.overlap = overlap1;
  //End of initializations    
  start_pos_prev = start_pos_level2;
  start_pos = start_pos_level1;
  //End of initializations

  nextlevlspiral_BM();    
  //nextlevlspiral_BM_minoverlap(1);
  //nextlevlspiral_BM_weightedoverlap(-1);
    
  for(int k = 0; k < 3; k++) //change back to 3
  {
	b_size = level1_block_size;
	temp_lambda = (level1_block_size*3) >> 2;
	while(b_size > 1)
	{
	  for (int l = 0; l < 4; l++) //change back to 4
	  {
	    calculate_pixel_reliability();
		add_smoothness8_overlap((l+1)*temp_lambda);
		//add_smoothness8_old((l+1)*temp_lambda);	  

	  }
	  setMVs_iter();
      b_size = (b_size >> 1);
	  temp_lambda = temp_lambda >> 1;
	}
  }  

  //validity_for_SR();

   
  //-----------------End of Level 1------------------------  

}
void Pair_Motion::calculate_reliability()
{
  int mod_factor = b_size - 1;
  int div_factor = int(log10((double)b_size)/log10(2.0));

  //clear overlap and reliability
  for(int i = 0; i < height1; i++)
  {
    for (int j = 0; j < width1; j++)
    {     
	  Computed_Data.reliability[i*step1+j] = 0;	
	  Computed_Data.overlap[i*step1+j] = 0;		  
    }
  }
  calculate_block_overlap();
  double max_value = find_overlap_max(); //(4*b_size*b_size) / (b_size*b_size + 4);
  double tl_overlap, bl_overlap, tr_overlap, br_overlap, total_overlap;
  int vx_pos, vy_pos, top_left, top_right, bottom_left, bottom_right;
  for (int i = start_pos; i < height - (add_height - start_pos); i+=b_size)//+)  //goes through all vertical motion positions
  {
    for (int j = start_pos; j < width - (add_width - start_pos); j+=b_size)//+) //goes through all horizontal motion positions
	{
 	  vx_pos = v_x[i*step+j] + j; 
	  vy_pos = v_y[i*step+j] + i; 

	  calculate_singleMV_overlap(b_size, vy_pos, vx_pos, top_left, top_right, bottom_left, bottom_right); 
	  
	  if (top_left == 0) //this shouldnt happen though -- safety check
	    tl_overlap = 0;
	  else
		tl_overlap = (double)(top_left % Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+((vx_pos >> div_factor) << div_factor)]) / (Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+((vx_pos >> div_factor) << div_factor)]);

	  if(bottom_left == 0)
		bl_overlap = 0;
	  else
	    bl_overlap = (double) (bottom_left % Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor) + b_size)*step+((vx_pos >> div_factor) << div_factor)]) / (Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor) + b_size)*step+((vx_pos >> div_factor) << div_factor)]); 
		  
	  if(top_right == 0)
	    tr_overlap = 0;
	  else
		tr_overlap = (double) (top_right % Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]) / (Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]); 
		  
	  if(bottom_right == 0)
		br_overlap = 0;
	  else
		br_overlap = (double) (bottom_right % Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor)+b_size)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]) / (Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor)+b_size)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]);
	  
	  total_overlap = tl_overlap + bl_overlap + tr_overlap + br_overlap;	  
	 
	 // total_overlap = (double)(top_left % Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+((vx_pos >> div_factor) << div_factor)]) / (Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+((vx_pos >> div_factor) << div_factor)])
		//  + (double) (bottom_left % Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor) + b_size)*step+((vx_pos >> div_factor) << div_factor)]) / (Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor) + b_size)*step+((vx_pos >> div_factor) << div_factor)]) 
		 // + (double) (top_right % Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]) / (Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]) 
		 // + (double) (bottom_right % Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor)+b_size)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]) / (Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor)+b_size)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]);
	  //total_overlap = (double) Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+((vx_pos >> div_factor) << div_factor)] / (top_left+1) + (double) Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor) + b_size)*step+((vx_pos >> div_factor) << div_factor)] / (bottom_left+1) + (double) Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+(((vx_pos >> div_factor) << div_factor)+b_size)] / (top_right+1) + (double) Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor)+b_size)*step+(((vx_pos >> div_factor) << div_factor)+b_size)] / (bottom_right+1);
	  
	  Computed_Data.reliability[i*step+j] = (total_overlap / max_value); //this is not actually the reliability, but a penalty weighting.
	}
  }
  
  //add_smoothness8_weighted(b_size << 1); //uses b_size = 2;  this is only postfiltering the MV field with weights on neighbors. 

}
void Pair_Motion::calculate_pixel_reliability()
{  
  //clear overlap and reliability
  for(int i = 0; i < height; i++)
  {
    for (int j = 0; j < width; j++)
    {     
	  //Computed_Data.reliability[i*step+j] = 0;	
	  Computed_Data.overlap[i*step+j] = 0;		  
    }
  }
    
  //double mean;
  calculate_pixel_overlap();
  
  /*
  int x_pos, y_pos;
  double block_overlap;
  int b_size_sq = 4;
  double temp_SAD;  

  if (b_size == 1 && level == 1)
  {
    b_size = 2;
	calculate_pixel_overlap();
	b_size = 1;
  for (int i = start_pos; i < height - (add_height - start_pos); i+=b_size)//+)  //goes through all vertical motion positions
  {
    for (int j = start_pos; j < width - (add_width - start_pos); j+=b_size)//+) //goes through all horizontal motion positions
	{
 	  x_pos = v_x[i*step+j] + j; 
	  y_pos = v_y[i*step+j] + i;
	  block_overlap = 0;
	  b_size = 2;
	  for(int k = y_pos; k < y_pos + b_size; k++)
	  {
	    for(int l = x_pos; l < x_pos + b_size; l++)
		{
		  block_overlap += Computed_Data.overlap[k*step+l]; 
		}
	  }
	  	 	  
	  temp_SAD = (double)calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
	  	  
	  //if(calculate_SAD(i, j, x_pos, y_pos) > (b_size*b_size*10))
        //Computed_Data.reliability[i*step+j] = 1;
	  //else
	    //Computed_Data.reliability[i*step+j] = (b_size_sq / block_overlap); 
	   Computed_Data.reliability[i*step+j] = b_size_sq/block_overlap; 
	   b_size = 1;
	  //std::cout << Computed_Data.reliability[i*step+j] << std::endl;
	}
  }
	b_size = 1;
  }

  //std::cout << "variance is " << variance << std::endl;
  */
}
void Pair_Motion::validity_for_SR()
{
  b_size = 2;
  double b_size_sq = b_size*b_size;

  //clear overlap and reliability
  for(int i = 0; i < height; i++)
  {
    for (int j = 0; j < width; j++)
    {     
	  Computed_Data.overlap[i*step+j] = 0;		  
    }
  }
    
  calculate_pixel_overlap();
  int x_pos, y_pos;
  double block_overlap, temp_SAD;

  for (int i = start_pos; i < height - (add_height - start_pos); i+=b_size) //goes through all vertical motion positions
  {
    for (int j = start_pos; j < width - (add_width - start_pos); j+=b_size) //goes through all horizontal motion positions
	{
 	  x_pos = v_x[i*step+j] + j; 
	  y_pos = v_y[i*step+j] + i;
	  block_overlap = 0;
	  
	  for(int k = y_pos; k < y_pos + b_size; k++)
	  {
	    for(int l = x_pos; l < x_pos + b_size; l++)
		{
		  block_overlap += Computed_Data.overlap[k*step+l]; 
		}
	  }
	  	 	  
	  temp_SAD = (double)calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
	  /*if (temp_SAD > 25)//18) //8
	  {
		for(int y = i; y < i + b_size; y++)
		{
		  for(int z = j; z < j + b_size; z++)
		  {
		    Computed_Data.reliability[y*step+z] = 0;
		  }
		}
		
	  }*/
	  //else
	  //{
	    for(int y = i; y < i + b_size; y++)
		{
		  for(int z = j; z < j + b_size; z++)
		  {
		    Computed_Data.reliability[y*step+z] = b_size_sq/(block_overlap*(1+temp_SAD));
	   	  }
		}	    
	  //}
	}
  }
	
  b_size = 1;

}
void Pair_Motion::calculate_block_overlap()
{
  int x_mod, y_mod, overlapx_opp, overlapy_opp, x_div, y_div;  
  int mod_factor = b_size - 1;
  int div_factor = int(log10((double)b_size)/log10(2.0));
  //12/1/11 -- changed to add_height - start_pos // was start_pos + b_size
  for (int i = start_pos; i < height-(add_height - start_pos); i+=b_size)//+)  //goes through all vertical motion positions
  {
    for (int j = start_pos; j < width-(add_width - start_pos); j+=b_size)//+) //goes through all horizontal motion positions
    {
      x_mod = (v_x[i*step+j] + j) & mod_factor;  //note that the & onyly works because (v_x[i*step+j] + i) & 3 <=> (v_x[i*step+j] + i) % 4;
	  y_mod = (v_y[i*step+j] + i) & mod_factor;  //x % 2 == x & 1, x % 4 == x & 3, x % 8 == x & 7

	  overlapy_opp = y_mod; //b_size - y_mod //NOTE:  PRETTY SURE THIS ONLY WORKS FOR 2x2 blocks!!
	  overlapx_opp = x_mod; //b_size - x_mod

	  x_div = ((v_x[i*step+j] + j) >> div_factor) << div_factor; //the part in parentheses causes a truncation -- don't think I can simplify this statement further.
	  y_div = ((v_y[i*step+j] + i) >> div_factor) << div_factor;

	  //y_div = ((v_x[i*step+j] + i) >> div_factor);// << div_factor; //the part in parentheses causes a truncation -- don't think I can simplify this statement further.
	  //x_div = ((v_y[i*step+j] + j) >> div_factor);// << div_factor;

	  if(x_mod == 0 && y_mod == 0)
	  {
	    Computed_Data.overlap[y_div*step+x_div] += (b_size*b_size);	//Total overlap = b_size + b_size
	  }
	  else if(x_mod != 0 && y_mod != 0) //overlap in four different blocks
	  {
        Computed_Data.overlap[y_div*step+x_div] += (b_size - x_mod)*(b_size - y_mod); //was just x_mod + y_mod
		Computed_Data.overlap[y_div*step+(x_div+b_size)] += overlapx_opp*(b_size - y_mod);
		Computed_Data.overlap[(y_div+b_size)*step+x_div] += (b_size - x_mod)*overlapy_opp;
		Computed_Data.overlap[(y_div+b_size)*step+(x_div+b_size)] += overlapx_opp*overlapy_opp;
	  }
	  else if(x_mod != 0) //overlap in x direction only
	  {
	    Computed_Data.overlap[y_div*step+x_div] += (b_size - x_mod)*b_size;
		Computed_Data.overlap[y_div*step+(x_div+b_size)] += b_size*overlapx_opp;
	  }
	  else //overlapy in y direction only
	  {
	    Computed_Data.overlap[y_div*step+x_div] += (b_size - y_mod)*b_size;
		Computed_Data.overlap[(y_div+b_size)*step+x_div] += b_size*overlapy_opp;
	  }
	}
  }
}
void Pair_Motion::calculate_pixel_overlap()
{
  int x_pos, y_pos;
  //double temp_SAD;
  //for variance
  //double n = 0;
  //double M2 = 0;
  //double delta;
  //double variance;
  //mean = 0;

  for (int i = start_pos; i < height-(add_height - start_pos); i+=b_size)
  {
    for (int j = start_pos; j < width-(add_width - start_pos); j+=b_size)
    {
      x_pos = v_x[i*step+j] + j;
	  y_pos = v_y[i*step+j] + i;

	  /*temp_SAD = calculate_SAD(i, j, x_pos, y_pos);
	  
      //Let's figure out the variance of the SAD values
	  n = n + 1;
	  delta = temp_SAD - mean;
	  mean = mean + delta/n;
	  M2 = M2 + delta*(temp_SAD - mean); */
	  
	  for(int k = y_pos; k < y_pos + b_size; k++)
	  {
	    for(int l = x_pos; l < x_pos + b_size; l++)
		{
		  Computed_Data.overlap[k*step+l] += 1;  //adding up all of the pixel overlaps
		}
	  }
    }	
  }
  //variance = M2/n;
}
double Pair_Motion::find_overlap_max()
{
  int vx_pos, vy_pos, top_left, top_right, bottom_left, bottom_right;
  double tl_overlap, bl_overlap, tr_overlap, br_overlap, total_overlap;
  int div_factor = int(log10((double)b_size)/log10(2.0));
  double max_value = 0;
 
  for (int i = start_pos; i < height - (add_height - start_pos); i+=b_size)//+)  //goes through all vertical motion positions
  {
	for (int j = start_pos; j < width - (add_width - start_pos); j+=b_size)//+) //goes through all horizontal motion positions
	{
	  vx_pos = v_x[i*step+j] + j; 
	  vy_pos = v_y[i*step+j] + i; 

	  calculate_singleMV_overlap(b_size, vy_pos, vx_pos, top_left, top_right, bottom_left, bottom_right); 
	  //total_overlap = (double)((top_left + 1) % (Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+((vx_pos >> div_factor) << div_factor)] + 1)) / (Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+((vx_pos >> div_factor) << div_factor)] + 1)
		 // + (double) ((bottom_left + 1) % (Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor) + b_size)*step+((vx_pos >> div_factor) << div_factor)] + 1)) / (Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor) + b_size)*step+((vx_pos >> div_factor) << div_factor)]+1) 
		//  + (double) ((top_right + 1) % (Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+(((vx_pos >> div_factor) << div_factor)+b_size)] + 1)) / (Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]+1) 
		//  + (double) ((bottom_right + 1) % (Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor)+b_size)*step+(((vx_pos >> div_factor) << div_factor)+b_size)] + 1)) / (Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor)+b_size)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]+1);
	  //	  (double) Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+((vx_pos >> div_factor) << div_factor)] / (top_left+1) + (double) Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor) + b_size)*step+((vx_pos >> div_factor) << div_factor)] / (bottom_left+1) + (double) Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+(((vx_pos >> div_factor) << div_factor)+b_size)] / (top_right+1) + (double) Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor)+b_size)*step+(((vx_pos >> div_factor) << div_factor)+b_size)] / (bottom_right+1);
	  
	  calculate_singleMV_overlap(b_size, vy_pos, vx_pos, top_left, top_right, bottom_left, bottom_right); 
	  
	  if (top_left == 0) //this shouldnt happen though -- safety check
	    tl_overlap = 0;
	  else
		tl_overlap = (double)(top_left % Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+((vx_pos >> div_factor) << div_factor)]) / (Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+((vx_pos >> div_factor) << div_factor)]);

	  if(bottom_left == 0)
		bl_overlap = 0;
	  else
	    bl_overlap = (double) (bottom_left % Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor) + b_size)*step+((vx_pos >> div_factor) << div_factor)]) / (Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor) + b_size)*step+((vx_pos >> div_factor) << div_factor)]); 
		  
	  if(top_right == 0)
	    tr_overlap = 0;
	  else
		tr_overlap = (double) (top_right % Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]) / (Computed_Data.overlap[((vy_pos >> div_factor) << div_factor)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]); 
		  
	  if(bottom_right == 0)
		br_overlap = 0;
	  else
		br_overlap = (double) (bottom_right % Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor)+b_size)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]) / (Computed_Data.overlap[(((vy_pos >> div_factor) << div_factor)+b_size)*step+(((vx_pos >> div_factor) << div_factor)+b_size)]);
	  
	  total_overlap = tl_overlap + bl_overlap + tr_overlap + br_overlap;
	  
	  if (total_overlap > max_value)
	    max_value = total_overlap;
	}
  }		
  return max_value;
}
int Pair_Motion::get_overlap_volume(int x_pos, int y_pos)
{
  int volume = 0;
 
  for(int k = y_pos; k < y_pos + b_size; k++)
  {
    for(int l = x_pos; l < x_pos + b_size; l++)
	{
	  volume += Computed_Data.overlap[k*step+l]; 
	}
  }
  return volume;
}
void Pair_Motion::calculate_singleMV_overlap(int b_size, int pos_y, int pos_x, int& top_left, int& top_right, int& bottom_left, int& bottom_right)
{
  int x_mod, y_mod, overlapx_opp, overlapy_opp;
  int mod_factor = b_size - 1;
    
  x_mod = (pos_x) & mod_factor;  //note that the & onyly works because (v_x[i*step+j] + i) & 3 <=> (v_x[i*step+j] + i) % 4;
  y_mod = (pos_y) & mod_factor;  //x % 2 == x & 1, x % 4 == x & 3, x % 8 == x & 7

  overlapy_opp = y_mod; //b_size - y_mod  //NOTE:  PRETTY SURE THIS ONLY WORKS FOR 2x2 blocks!!
  overlapx_opp = x_mod; //b_size - x_mod
  
  if(x_mod == 0 && y_mod == 0)
  {
    top_left = (b_size*b_size);	//Total overlap = b_size + b_size
	top_right = 0;
	bottom_left = 0;
	bottom_right = 0;
  }
  else if(x_mod != 0 && y_mod != 0) //overlap in four different blocks
  {
    top_left = (b_size - x_mod)*(b_size - y_mod);
    top_right = overlapx_opp*(b_size - y_mod);
    bottom_left = (b_size - x_mod)*overlapy_opp;
    bottom_right = overlapx_opp*overlapy_opp;
  }    
  else if(x_mod != 0) //overlap in x direction only
  {
    top_left = (b_size - x_mod)*b_size;
	top_right = b_size*overlapx_opp;
	bottom_left = 0;
	bottom_right = 0;
  }
  else //overlapy in y direction only
  {
    top_left = (b_size - y_mod)*b_size;
	top_right = 0;
	bottom_left = b_size*overlapy_opp;
	bottom_right = 0;
  }
}
void Pair_Motion::onelevl_BM()
{  
  int SAD_value, min_SAD_value;
  int motion_pos_x[1], motion_pos_y[1];		
  int i,j,m,n,l,k; //x,y;

  //Reduce computations
  //int b_size_div = b_size >> 1;
  //int search_size_div = search_size >> 1;

  //Shifting block inside search window
  //int v_shift;
  //v_shift = search_size - b_size; //this is if we shift by blocksize/2
  
  //These first two loops are for each pixel in the image.  For the meantime, we start at search_size/2 and ignore edges
  for (i = start_pos; i <= height-(start_pos+b_size/*+v_shift*/); i+=b_size)//+)  //goes through all vertical pixels
  {
	for (j = start_pos; j <= width-(start_pos+b_size/*+v_shift*/); j+=b_size)//+)  //goes through all horizontal pizels
    {
	   l = (j - start_pos); //only correct because first level -- should be vshift/2
	   k = (i - start_pos);

	   SAD_value = calculate_SAD(i,j,l,k);	  
	   //Initializations of SAD value
	   min_SAD_value = SAD_value;
	   motion_pos_x[0] = l;
	   motion_pos_y[0] = k;

      for (n = 0; n <= v_shift; n++) //n controls the vertical shift
      {		
	    for (m = 0; m <= v_shift; m++) //m controls the horizontal shift
	    {         	      
		  SAD_value = calculate_SAD(i,j,l,k);		  	  
          if (SAD_value < min_SAD_value)
		  {		    		    
		    min_SAD_value = SAD_value;
			motion_pos_x[0] = l;
            motion_pos_y[0] = k;
		  }		 
		  l++;
		}
		l = (j - start_pos);
		k++;
	  }		
	  	
	  // Here we do our assignments for the block with the lowest SAD.
	  v_x[i*step+j] = motion_pos_x[0] - j; //store motion vector x in array
	  v_y[i*step+j] = motion_pos_y[0] - i; //store motion vector y in array
	  //v_x_old[i*step+j] = motion_pos_x[0] - i; //store motion vector x in array -- this is for testing original MV in smoothness fn
	  //v_y_old[i*step+j] = motion_pos_y[0] - j; //store motion vector y in array -- this is for testing original MV in smoothness fn
    }  
  }
}
void Pair_Motion::onelevlspiral_BM()
{  
  int SAD_value, min_SAD_value;
  int motion_pos_x, motion_pos_y;		
  int i,j,m,l,k,t; 

   //These first two loops are for each pixel in the image.  For the meantime, we start at search_size/2 and ignore edges
  for (i = start_pos; i < height-(add_height - start_pos); i+=b_size)  //goes through all vertical pixels
  {
	for (j = start_pos; j < width-(add_width - start_pos); j+=b_size)  //goes through all horizontal pizels
    {    
	  l = j;
	  k = i;

	  //This is the same block as the block of current frame
	  SAD_value = calculate_SAD(i,j,l,k);
	  
	  //Initializations of SAD value
	  min_SAD_value = SAD_value;
	  motion_pos_x = l;// + (b_size_div-1);
	  motion_pos_y = k;// + (b_size_div-1);
	    
	  //This first outer loop is used to do the spiral search
	  //We are repeating patterns of moving right,down, left, then up.
	  //At the very end, we go right once more to finish things off.
	  //The algorithm is basically:
	  //right, down, left(m+1), up(m+1).  And then right(m+1) at the very end.
	  for (m = 1; m < v_shift; m+=2)
	  { 
	    //the variable m will tell us how much to shift each time
	    //the variable t is a counter.  if we have to shift 5 times, we will
	    //shift one position at a time and calculate the SAD for each shift.
	    for(t = 0; t < m; t++) 
		{
	      l = l + 1; //m;
		  SAD_value = calculate_SAD(i,j,l,k);
          if (SAD_value < min_SAD_value)
		  {		    		    
		    min_SAD_value = SAD_value;
			motion_pos_x = l;// + (b_size_div-1);
            motion_pos_y = k;// + (b_size_div-1);
		  }		 
		}
	    
		for(t = 0; t < m; t++)
		{
	      k = k + 1; //m;
          SAD_value = calculate_SAD(i,j,l,k);
          if (SAD_value < min_SAD_value)
		  {		    
		    min_SAD_value = SAD_value;
			motion_pos_x = l;// + (b_size_div-1);
            motion_pos_y = k;// + (b_size_div-1);
		  }		  
		}
	   
		for(t = 0; t < m+1; t++)
		{
	      l = l - 1; //(m + 1);
		  SAD_value = calculate_SAD(i,j,l,k);
          if (SAD_value < min_SAD_value)
		  {		    		    
		    min_SAD_value = SAD_value;
			motion_pos_x = l;// + (b_size_div-1);
            motion_pos_y = k;// + (b_size_div-1);
		  }		  
		}

		for(t = 0; t < m+1; t++)
		{
	      k = k - 1; //(m + 1);
		  SAD_value = calculate_SAD(i,j,l,k);
          if (SAD_value < min_SAD_value)
		  {		    		    
		    min_SAD_value = SAD_value;
			motion_pos_x = l;// + (b_size_div-1);
            motion_pos_y = k;// + (b_size_div-1);
		  }		  
		}
	  }

	  //This is what we do at the end to move across the top row.
      for(t = 0; t < (m-1); t++)
	  {
	    l = l + 1; //m;
		SAD_value = calculate_SAD(i,j,l,k);
		if (SAD_value < min_SAD_value)
		{
		  min_SAD_value = SAD_value;
		  motion_pos_x = l;// + (b_size_div-1);
          motion_pos_y = k;// + (b_size_div-1);
		}		
	  }
	  	
	  // Here we do our assignments for the block with the lowest SAD.
	  v_x[i*step+j] = motion_pos_x - j; //store motion vector x in array
	  v_y[i*step+j] = motion_pos_y - i; //store motion vector y in array
	  
    }  
  }
}
void Pair_Motion::onelevlspiral_BM_minoverlap()
{  
  int SAD_value, min_SAD_value, min_SAD_overlap, SAD_overlap;
  int motion_pos_x, motion_pos_y;		
  int i,j,m,l,k,t,y,z; 

   //These first two loops are for each pixel in the image.  For the meantime, we start at search_size/2 and ignore edges
  for (i = start_pos; i < height-(add_height - start_pos); i+=b_size)  //goes through all vertical pixels
  {
	for (j = start_pos; j < width-(add_width - start_pos); j+=b_size)  //goes through all horizontal pizels
    {    
	  l = j;
	  k = i;

	  //This is the same block as the block of current frame
	  min_SAD_value = calculate_SAD(i,j,l,k);
	  min_SAD_overlap = min_SAD_value + get_overlap_volume(l, k);
	  motion_pos_x = l;
	  motion_pos_y = k;
	    
	  //This first outer loop is used to do the spiral search
	  //We are repeating patterns of moving right,down, left, then up.
	  //At the very end, we go right once more to finish things off.
	  //The algorithm is basically:
	  //right, down, left(m+1), up(m+1).  And then right(m+1) at the very end.
	  for (m = 1; m < v_shift; m+=2)
	  { 
	    //the variable m will tell us how much to shift each time
	    //the variable t is a counter.  if we have to shift 5 times, we will
	    //shift one position at a time and calculate the SAD for each shift.
	    for(t = 0; t < m; t++) 
		{
	      l = l + 1; //m;
		  SAD_value = calculate_SAD(i,j,l,k);
          if (SAD_value < min_SAD_value)
		  {		    		    
		    min_SAD_value = SAD_value;
			min_SAD_overlap = SAD_value + get_overlap_volume(l, k);
			motion_pos_x = l;
            motion_pos_y = k;
		  }		 
		  else if (SAD_value == min_SAD_value)
		  {
		    SAD_overlap = SAD_value + get_overlap_volume(l, k);			

            if (SAD_overlap < min_SAD_overlap)
		    {
		      min_SAD_overlap = SAD_overlap;
			  motion_pos_x = l; 
              motion_pos_y = k; 
		    }
		  }
		}
	    
		for(t = 0; t < m; t++)
		{
	      k = k + 1; //m;
          SAD_value = calculate_SAD(i,j,l,k);
          if (SAD_value < min_SAD_value)
		  {		    		    
		    min_SAD_value = SAD_value;
			min_SAD_overlap = SAD_value + get_overlap_volume(l, k);
			motion_pos_x = l;
            motion_pos_y = k;
		  }		 
		  else if (SAD_value == min_SAD_value)
		  {
		    SAD_overlap = SAD_value + get_overlap_volume(l, k);			

            if (SAD_overlap < min_SAD_overlap)
		    {
		      min_SAD_overlap = SAD_overlap;
			  motion_pos_x = l; 
              motion_pos_y = k; 
		    }
		  }		  
		}
	   
		for(t = 0; t < m+1; t++)
		{
	      l = l - 1; //(m + 1);
		  SAD_value = calculate_SAD(i,j,l,k);
          if (SAD_value < min_SAD_value)
		  {		    		    
		    min_SAD_value = SAD_value;
			min_SAD_overlap = SAD_value + get_overlap_volume(l, k);
			motion_pos_x = l;
            motion_pos_y = k;
		  }		 
		  else if (SAD_value == min_SAD_value)
		  {
		    SAD_overlap = SAD_value + get_overlap_volume(l, k);			

            if (SAD_overlap < min_SAD_overlap)
		    {
		      min_SAD_overlap = SAD_overlap;
			  motion_pos_x = l; 
              motion_pos_y = k; 
		    }
		  }		  
		}

		for(t = 0; t < m+1; t++)
		{
	      k = k - 1; //(m + 1);
		  SAD_value = calculate_SAD(i,j,l,k);
          if (SAD_value < min_SAD_value)
		  {		    		    
		    min_SAD_value = SAD_value;
			min_SAD_overlap = SAD_value + get_overlap_volume(l, k);
			motion_pos_x = l;
            motion_pos_y = k;
		  }		 
		  else if (SAD_value == min_SAD_value)
		  {
		    SAD_overlap = SAD_value + get_overlap_volume(l, k);			

            if (SAD_overlap < min_SAD_overlap)
		    {
		      min_SAD_overlap = SAD_overlap;
			  motion_pos_x = l; 
              motion_pos_y = k; 
		    }
		  }		  
		}
	  }

	  //This is what we do at the end to move across the top row.
      for(t = 0; t < (m-1); t++)
	  {
	    l = l + 1; //m;
		SAD_value = calculate_SAD(i,j,l,k);
        if (SAD_value < min_SAD_value)
		{		    		    
		  min_SAD_value = SAD_value;
		  min_SAD_overlap = SAD_value + get_overlap_volume(l, k);
		  motion_pos_x = l;
          motion_pos_y = k;
		}		 
		else if (SAD_value == min_SAD_value)
		{
		  SAD_overlap = SAD_value + get_overlap_volume(l, k);			
		  
		  if (SAD_overlap < min_SAD_overlap)
		  {
		    min_SAD_overlap = SAD_overlap;
			motion_pos_x = l; 
            motion_pos_y = k; 
		  }
		}		
	  }
	  	
	  //Put in overlap for chosen block
	  for(y = motion_pos_y; y < motion_pos_y + b_size; y++)
	    for(z = motion_pos_x; z < motion_pos_x + b_size; z++)
		  Computed_Data.overlap[y*step+z] += 1;

	  // Here we do our assignments for the block with the lowest SAD.
	  v_x[i*step+j] = motion_pos_x - j; //store motion vector x in array
	  v_y[i*step+j] = motion_pos_y - i; //store motion vector y in array
	  
    }  
  }
}
void Pair_Motion::nextlevl_BM()
{  
  int SAD_value, min_SAD_value;
  int motion_pos_x[1], motion_pos_y[1];
 
  // These are all our initializations before running BM algorithm.
  SAD_value = 0; // result of SAD computation stored here
  		
  int i,j,m,n,l,k;
  int new_i, new_j; //these will be i and j scaled by level number
  int pos_x, pos_y;
        
  //These first two loops are for each MV calculated from previous.  
  for (i = start_pos_prev; i <= height_previous-(start_pos_prev+b_size_previous/*+v_shift_previous*/); i+=b_size_previous)//+)  //goes through all vertical MVs
  {
	for (j = start_pos_prev; j <= width_previous-(start_pos_prev+b_size_previous/*+v_shift_previous*/); j+=b_size_previous)//+)  //goes through all horizontal MVs
    {
	  new_i = (i << 1); //multiply by 2 for next level
	  new_j = (j << 1); //multiply by 2 for next level
	  
	  //position in frame 1
	  l = (new_j + (v_xprev[i*step_previous+j] << 1)); //multiply by 2 for next level -- frame 1 positions (change in search)
	  k = (new_i + (v_yprev[i*step_previous+j] << 1));
	  	       
	  l = l - (v_shift >> 1);
	  k = k - (v_shift >> 1);

	  if ((l < b_size) || (k < b_size)) //this shouldn't happen!  can remove?
	    continue;
	    
	  pos_x = l;
	  pos_y = k;
	  
	  //This is the same block as the block of current frame
	  SAD_value = calculate_SAD(new_i,new_j,l,k);
	  min_SAD_value = SAD_value;
	  motion_pos_x[0] = l;
      motion_pos_y[0] = k;
	  
      for (n = 0; n <= v_shift; n++) //m = m + b_size_div) //we could shift the block by 1 pixel, but we shift it by b_size/2 to reduce computation
      {		
	    for (m = 0; m <= v_shift; m++) //n = n + b_size_div)
	    { 
		  SAD_value = calculate_SAD(new_i,new_j,pos_x,pos_y);		  	  
          if (SAD_value < min_SAD_value)
		  {		    		    
		    min_SAD_value = SAD_value;
			motion_pos_x[0] = pos_x;
            motion_pos_y[0] = pos_y;
		  }		 
		  pos_x++;	
		}
		pos_x = l;
		pos_y++;
	  } 

	  // Here we do our assignments for the block with the lowest SAD.	  
	  v_x[new_i*step+new_j] = motion_pos_x[0] - new_j; //store motion vector x in array
	  v_y[new_i*step+new_j] = motion_pos_y[0] - new_i; //store motion vector y in array	  	 
	  //v_x_old[new_i*step+new_j] = motion_pos_x[0] - new_i; //store motion vector x in array -- this is for testing original MV in smoothness fn
	  //v_y_old[new_i*step+new_j] = motion_pos_y[0] - new_j; //store motion vector y in array -- this is for testing original MV in smoothness fn
    }  
  }
}
void Pair_Motion::nextlevlspiral_BM()
{  
  int SAD_value, min_SAD_value;
  int motion_pos_x, motion_pos_y; 
  		
  int i,j,m,l,k,t;
  int new_i, new_j; //these will be i and j scaled by level number
  int pos_x, pos_y; //these are positions of MVs calcualted from previous level
       
  //These first two loops are for each MV calculated from previous.  
  for (i = start_pos_prev; i < height_previous-((add_height >> 1) - start_pos_prev); i+=b_size_previous) //goes through all vertical MVs
  {
	for (j = start_pos_prev; j < width_previous-((add_width >> 1) - start_pos_prev); j+=b_size_previous) //goes through all horizontal MVs
    {
	  new_i = (i << 1); //multiply by 2 for next level
	  new_j = (j << 1); //multiply by 2 for next level

      //we want to form a block around motion vector position i,j.  
	  l = (new_j + (v_xprev[i*step_previous+j] << 1)); //multiply by 2 for next level
	  k = (new_i + (v_yprev[i*step_previous+j] << 1));
     
      pos_x = new_j;
      pos_y = new_i; 

	  //This is the same block as the block of current frame
	  SAD_value = calculate_SAD(pos_y,pos_x,l,k);
	  min_SAD_value = SAD_value;
	  motion_pos_x = l; 
      motion_pos_y = k; 
	  
	  //This first outer loop is used to do the spiral search
	  //We are repeating patterns of moving right,down, left, then up.
	  //At the very end, we go right once more to finish things off.
	  //The algorithm is basically:
	  //right, down, left(m+1), up(m+1).  And then right(m+1) at the very end.
	  for (m = 1; m < v_shift; m+=2)
	  { 
	    //the variable m will tell us how much to shift each time
	    //the variable t is a counter.  if we have to shift 5 times, we will
	    //shift one position at a time and calculate the SAD for each shift.
	    for(t = 0; t < m; t++) 
		{
	      l = l + 1; //m;
		  SAD_value = calculate_SAD(pos_y,pos_x,l,k);
          if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;
			motion_pos_x = l; 
            motion_pos_y = k; 
		  }
		}
	    
		for(t = 0; t < m; t++)
		{
	      k = k + 1; //m;
          SAD_value = calculate_SAD(pos_y,pos_x,l,k);
          if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;
			motion_pos_x = l; 
            motion_pos_y = k; 
		  }
		}
	   
		for(t = 0; t < m+1; t++)
		{
	      l = l - 1; //(m + 1);
          SAD_value = calculate_SAD(pos_y,pos_x,l,k);
          if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;
			motion_pos_x = l;
            motion_pos_y = k; 
		  }
		}

		for(t = 0; t < m+1; t++)
		{
	      k = k - 1; //(m + 1);
		  SAD_value = calculate_SAD(pos_y,pos_x,l,k);
          if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;
			motion_pos_x = l; 
            motion_pos_y = k;  
		  }
		}
	  }

	  //This is what we do at the end to move across the top row.
      for(t = 0; t < (m-1); t++)
	  {
	    l = l + 1; //m;
	    SAD_value = calculate_SAD(pos_y,pos_x,l,k);
        if (SAD_value < min_SAD_value)
		{
		  min_SAD_value = SAD_value;
		  motion_pos_x = l; 
          motion_pos_y = k; 
		}
	  }	  	  
	  
	  // Here we do our assignments for the block with the lowest SAD.	 
	  v_x[new_i*step+new_j] = motion_pos_x - new_j; //store motion vector x in array
	  v_y[new_i*step+new_j] = motion_pos_y - new_i; //store motion vector y in array	  
	  
    }  
  }
}
void Pair_Motion::nextlevlspiral_BM_minoverlap(int start_direction)
{  
  int SAD_value, min_SAD_value, min_SAD_overlap, SAD_overlap;
  int motion_pos_x, motion_pos_y; 
  //int overlap_volume = b_size*b_size;
  		
  //CvPoint pt1, pt2;

  int i,j,m,l,k,t,y,z;
  int new_i, new_j; //these will be i and j scaled by level number
  int pos_x, pos_y; //these are positions of MVs calcualted from previous level
       
  //These first two loops are for each MV calculated from previous.  
  for (i = start_pos_prev; i < height_previous-((add_height >> 1) - start_pos_prev/* + b_size_previous*/); i+=b_size_previous) //goes through all vertical MVs
  {
	for (j = start_pos_prev; j < width_previous-((add_width >> 1) - start_pos_prev/* + b_size_previous*/); j+=b_size_previous) //goes through all horizontal MVs
    {
	  new_i = (i << 1); //multiply by 2 for next level
	  new_j = (j << 1); //multiply by 2 for next level

      //we want to form a block around motion vector position i,j.  
	  l = (new_j + (v_xprev[i*step_previous+j] << 1)); //multiply by 2 for next level
	  k = (new_i + (v_yprev[i*step_previous+j] << 1));
     
      pos_x = new_j;
      pos_y = new_i; 

	  //This is the same block as the block of current frame
	  SAD_value = calculate_SAD(pos_y,pos_x,l,k);
	  min_SAD_value = SAD_value;

	  /*//Put in overlap for this block
	  for(y = k; y < k + b_size; y++)
	    for(z = l; z < l + b_size; z++)
		  Computed_Data.overlap[y*step+z] += 1; */

	  SAD_overlap = SAD_value + get_overlap_volume(l, k);// + overlap_volume;
	  min_SAD_overlap = SAD_overlap;
	  motion_pos_x = l; 
      motion_pos_y = k; 

	  /*//Take out overlap
	  for(y = k; y < k + b_size; y++)
	    for(z = l; z < l + b_size; z++)
		  Computed_Data.overlap[y*step+z] -= 1; */
	  
	  //This first outer loop is used to do the spiral search
	  //We are repeating patterns of moving right,down, left, then up.
	  //At the very end, we go right once more to finish things off.
	  //The algorithm is basically:
	  //right, down, left(m+1), up(m+1).  And then right(m+1) at the very end.
	  for (m = 1; m < v_shift; m+=2)
	  { 
	    //the variable m will tell us how much to shift each time
	    //the variable t is a counter.  if we have to shift 5 times, we will
	    //shift one position at a time and calculate the SAD for each shift.
	    for(t = 0; t < m; t++) 
		{
	      l = l + start_direction; //m;
		  SAD_value = calculate_SAD(pos_y,pos_x,l,k);
		  if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;						
			min_SAD_overlap = SAD_value + get_overlap_volume(l, k);
			motion_pos_x = l; 
            motion_pos_y = k; 		    
		  }
		  else if (SAD_value == min_SAD_value)
		  {  
			SAD_overlap = SAD_value + get_overlap_volume(l, k);			

            if (SAD_overlap < min_SAD_overlap)
		    {
		      min_SAD_overlap = SAD_overlap;
			  motion_pos_x = l; 
              motion_pos_y = k; 
		    }
		  }
		}
	    
		for(t = 0; t < m; t++)
		{
	      k = k + 1; //m;
          SAD_value = calculate_SAD(pos_y,pos_x,l,k);
		  if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;						
			min_SAD_overlap = SAD_value + get_overlap_volume(l, k);
			motion_pos_x = l; 
            motion_pos_y = k; 		    
		  }
		  else if (SAD_value == min_SAD_value)
		  {  
			SAD_overlap = SAD_value + get_overlap_volume(l, k);			

            if (SAD_overlap < min_SAD_overlap)
		    {
		      min_SAD_overlap = SAD_overlap;
			  motion_pos_x = l; 
              motion_pos_y = k; 
		    }
		  }
		}
	   
		for(t = 0; t < m+1; t++)
		{
	      l = l - start_direction; //(m + 1);
          SAD_value = calculate_SAD(pos_y,pos_x,l,k);
		  if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;						
			min_SAD_overlap = SAD_value + get_overlap_volume(l, k);
			motion_pos_x = l; 
            motion_pos_y = k; 		    
		  }
		  else if (SAD_value == min_SAD_value)
		  {  
			SAD_overlap = SAD_value + get_overlap_volume(l, k);			

            if (SAD_overlap < min_SAD_overlap)
		    {
		      min_SAD_overlap = SAD_overlap;
			  motion_pos_x = l; 
              motion_pos_y = k; 
		    }
		  }
		}

		for(t = 0; t < m+1; t++)
		{
	      k = k - 1; //(m + 1);
		  SAD_value = calculate_SAD(pos_y,pos_x,l,k);
		  if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;						
			min_SAD_overlap = SAD_value + get_overlap_volume(l, k);
			motion_pos_x = l; 
            motion_pos_y = k; 		    
		  }
		  else if (SAD_value == min_SAD_value)
		  {  
			SAD_overlap = SAD_value + get_overlap_volume(l, k);			

            if (SAD_overlap < min_SAD_overlap)
		    {
		      min_SAD_overlap = SAD_overlap;
			  motion_pos_x = l; 
              motion_pos_y = k; 
		    }
		  }
		}
	  }

	  //This is what we do at the end to move across the top row.
      for(t = 0; t < (m-1); t++)
	  {
	    l = l + start_direction; //m;
	    SAD_value = calculate_SAD(pos_y,pos_x,l,k);
		if (SAD_value < min_SAD_value)
		{
		  min_SAD_value = SAD_value;						
		  min_SAD_overlap = SAD_value + get_overlap_volume(l, k);
		  motion_pos_x = l; 
          motion_pos_y = k; 		    
		}
		else if (SAD_value == min_SAD_value)
		{  
		  SAD_overlap = SAD_value + get_overlap_volume(l, k);			

          if (SAD_overlap < min_SAD_overlap)
		  {
		    min_SAD_overlap = SAD_overlap;
			motion_pos_x = l; 
            motion_pos_y = k; 
		  }
		}
	  }

	  //Put in overlap for chosen block
	  for(y = motion_pos_y; y < motion_pos_y + b_size; y++)
	    for(z = motion_pos_x; z < motion_pos_x + b_size; z++)
		  Computed_Data.overlap[y*step+z] += 1; 

	  /*if(get_overlap_volume(motion_pos_x, motion_pos_y) > (b_size*b_size))
	  {
	    pt1.x = motion_pos_x;
		pt1.y = motion_pos_y;
		pt2.x = motion_pos_x + (b_size - 1);
		pt2.y = motion_pos_y + (b_size - 1);
		
		cvRectangle(debug, pt1, pt2, CV_RGB(94,13,91), 1, 8, 0);
		//cvShowImage("main", debug);
		//cvWaitKey(30);

		for(y = motion_pos_y; y < motion_pos_y + b_size; y++)
	      for(z = motion_pos_x; z < motion_pos_x + b_size; z++)
		    Computed_Data.overlap[y*step+z] -= 1; 
	  }*/
	  
	  // Here we do our assignments for the block with the lowest SAD.	 
	  v_x[new_i*step+new_j] = motion_pos_x - new_j; //store motion vector x in array
	  v_y[new_i*step+new_j] = motion_pos_y - new_i; //store motion vector y in array	  
	  
    }  
  }
}
void Pair_Motion::nextlevlspiral_BM_weightedoverlap(int start_direction)
{  
  double min_SAD_overlap, SAD_overlap;
  int motion_pos_x, motion_pos_y; 
  int bs_squared = b_size*b_size;
 
  int i,j,m,l,k,t,y,z;
  int new_i, new_j; //these will be i and j scaled by level number
  int pos_x, pos_y; //these are positions of MVs calcualted from previous level
       
  //These first two loops are for each MV calculated from previous.  
  for (i = start_pos_prev; i < height_previous-((add_height >> 1) - start_pos_prev/* + b_size_previous*/); i+=b_size_previous) //goes through all vertical MVs
  {
	for (j = start_pos_prev; j < width_previous-((add_width >> 1) - start_pos_prev/* + b_size_previous*/); j+=b_size_previous) //goes through all horizontal MVs
    {
	  new_i = (i << 1); //multiply by 2 for next level
	  new_j = (j << 1); //multiply by 2 for next level

      //we want to form a block around motion vector position i,j.  
	  l = (new_j + (v_xprev[i*step_previous+j] << 1)); //multiply by 2 for next level
	  k = (new_i + (v_yprev[i*step_previous+j] << 1));
     
      pos_x = new_j;
      pos_y = new_i; 

	  //This is the same block as the block of current frame
	  SAD_overlap = calculate_SAD(pos_y,pos_x,l,k)*(1 + (double)get_overlap_volume(l, k)/bs_squared);// + overlap_volume;
	  min_SAD_overlap = SAD_overlap;
	  motion_pos_x = l; 
      motion_pos_y = k; 	
	  
	  //This first outer loop is used to do the spiral search
	  //We are repeating patterns of moving right,down, left, then up.
	  //At the very end, we go right once more to finish things off.
	  //The algorithm is basically:
	  //right, down, left(m+1), up(m+1).  And then right(m+1) at the very end.
	  for (m = 1; m < v_shift; m+=2)
	  { 
	    //the variable m will tell us how much to shift each time
	    //the variable t is a counter.  if we have to shift 5 times, we will
	    //shift one position at a time and calculate the SAD for each shift.
	    for(t = 0; t < m; t++) 
		{
	      l = l + start_direction; //m;
		  SAD_overlap = calculate_SAD(pos_y,pos_x,l,k)*(1 + (double)get_overlap_volume(l, k)/bs_squared);
		  if (SAD_overlap < min_SAD_overlap)
		  {		    					
			min_SAD_overlap = SAD_overlap;
			motion_pos_x = l; 
            motion_pos_y = k; 		    
		  }		 
		}
	    
		for(t = 0; t < m; t++)
		{
	      k = k + 1; //m;
          SAD_overlap = calculate_SAD(pos_y,pos_x,l,k)*(1 + (double)get_overlap_volume(l, k)/bs_squared);
		  if (SAD_overlap < min_SAD_overlap)
		  {		    					
			min_SAD_overlap = SAD_overlap;
			motion_pos_x = l; 
            motion_pos_y = k; 		    
		  }
		}
	   
		for(t = 0; t < m+1; t++)
		{
	      l = l - start_direction; //(m + 1);
          SAD_overlap = calculate_SAD(pos_y,pos_x,l,k)*(1 + (double)get_overlap_volume(l, k)/bs_squared);
		  if (SAD_overlap < min_SAD_overlap)
		  {		    					
			min_SAD_overlap = SAD_overlap;
			motion_pos_x = l; 
            motion_pos_y = k; 		    
		  }
		}

		for(t = 0; t < m+1; t++)
		{
	      k = k - 1; //(m + 1);
		  SAD_overlap = calculate_SAD(pos_y,pos_x,l,k)*(1 + (double)get_overlap_volume(l, k)/bs_squared);
		  if (SAD_overlap < min_SAD_overlap)
		  {		    					
			min_SAD_overlap = SAD_overlap;
			motion_pos_x = l; 
            motion_pos_y = k; 		    
		  }
		}
	  }

	  //This is what we do at the end to move across the top row.
      for(t = 0; t < (m-1); t++)
	  {
	    l = l + start_direction; //m;
	    SAD_overlap = calculate_SAD(pos_y,pos_x,l,k)*(1 + (double)get_overlap_volume(l, k)/bs_squared);
		if (SAD_overlap < min_SAD_overlap)
		{		    					
		  min_SAD_overlap = SAD_overlap;
		  motion_pos_x = l; 
          motion_pos_y = k; 		    
		}
	  }

	  //Put in overlap for chosen block
	  for(y = motion_pos_y; y < motion_pos_y + b_size; y++)
	    for(z = motion_pos_x; z < motion_pos_x + b_size; z++)
		  Computed_Data.overlap[y*step+z] += 1; 

	  /*if(get_overlap_volume(motion_pos_x, motion_pos_y) > (b_size*b_size))
	  {
	    pt1.x = motion_pos_x;
		pt1.y = motion_pos_y;
		pt2.x = motion_pos_x + (b_size - 1);
		pt2.y = motion_pos_y + (b_size - 1);
		
		cvRectangle(debug, pt1, pt2, CV_RGB(94,13,91), 1, 8, 0);
		//cvShowImage("main", debug);
		//cvWaitKey(30);

		for(y = motion_pos_y; y < motion_pos_y + b_size; y++)
	      for(z = motion_pos_x; z < motion_pos_x + b_size; z++)
		    Computed_Data.overlap[y*step+z] -= 1; 
	  }*/
	  
	  // Here we do our assignments for the block with the lowest SAD.	 
	  v_x[new_i*step+new_j] = motion_pos_x - new_j; //store motion vector x in array
	  v_y[new_i*step+new_j] = motion_pos_y - new_i; //store motion vector y in array	  
	  
    }  
  }
}
void Pair_Motion::nextlevlspiral_BM_adapt()
{  
  int SAD_value, min_SAD_value;
  int minMV_x, minMV_y, minMV_SAD, minMV_cost, MV_cost;
  int motion_pos_x[1], motion_pos_y[1], v_xtemp, v_ytemp;
 		
  int i,j,m,l,k,t;
  int new_i, new_j; //these will be i and j scaled by level number
  int pos_x, pos_y; //these are positions of MVs calcualted from previous level    
    
  //All of these entires use priors to calculate best MV.  Need the start_pos+b_size_previous because we already did start_pos elements.
  for (i = start_pos_prev/*+b_size_previous*/; i < height_previous-((add_height >> 1) - start_pos_prev + b_size_previous); i+=b_size_previous)//+)  //goes through all vertical MVs
  {
	for (j = start_pos_prev/*+b_size_previous*/; j < width_previous-((add_width >> 1) - start_pos_prev + b_size_previous); j+=b_size_previous)//+)  //goes through all horizontal MVs
    {
	  new_i = (i << 1); //multiply by 2 for next level
	  new_j = (j << 1); //multiply by 2 for next level

      //we want to form a block around motion vector position i,j.  
	  l = (new_j + (v_xprev[i*step_previous+j] << 1));// - (b_size_div-1); //multiply by 2 for next level
	  k = (new_i + (v_yprev[i*step_previous+j] << 1));// - (b_size_div-1);

	  //if (l < v_shift)
		 // continue;
	  //if (k < v_shift)
		  //continue;
	       
      pos_x = new_j;
      pos_y = new_i;	 

	  //This is the same block as the block of current frame
	  SAD_value = calculate_SAD(pos_y,pos_x,l,k);
	  min_SAD_value = SAD_value;
	  motion_pos_x[0] = l; 
      motion_pos_y[0] = k; 

	  //Initializations for min. position search
	  minMV_x = 0;
	  minMV_y = 0;
	  minMV_SAD = SAD_value;
	  //Since the first MV is the 0,0 MV, no need to subtract the MVs below from 0,0.  
	  //minMV_cost = abs(v_x[new_i*step+(new_j-b_size)]) + abs(v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs(v_x[(new_i-b_size)*step+new_j]) + abs(v_y[new_i*step+(new_j-b_size)])  + abs(v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs(v_y[(new_i-b_size)*step+new_j]);
	  //new neighbor
	  //minMV_cost = abs(v_x[new_i*step+(new_j-b_size)]) + abs(v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs(v_x[(new_i-b_size)*step+new_j]) + abs(v_x[(new_i-b_size)*step+(new_j+b_size)]) + abs(v_y[new_i*step+(new_j-b_size)])  + abs(v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs(v_y[(new_i-b_size)*step+new_j]) + abs(v_y[(new_i-b_size)*step+(new_j+b_size)]);
	  //prev level
	  minMV_cost = (abs(v_xprev[i*step_previous+j] << 1) + abs(v_xprev[i*step_previous+(j-b_size_previous)] << 1) + abs(v_xprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1) + abs(v_xprev[(i-b_size_previous)*step_previous+j] << 1) + abs(v_xprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1) + abs(v_xprev[i*step_previous+(j+b_size_previous)] << 1) + abs(v_xprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1) + abs(v_xprev[(i+b_size_previous)*step_previous+j] << 1) + abs(v_xprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1) + abs(v_yprev[i*step_previous+j] << 1) + abs(v_yprev[i*step_previous+(j-b_size_previous)] << 1) + abs(v_yprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1) + abs(v_yprev[(i-b_size_previous)*step_previous+j] << 1) + abs(v_yprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1) + abs(v_yprev[i*step_previous+(j+b_size_previous)] << 1) + abs(v_yprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1) + abs(v_yprev[(i+b_size_previous)*step_previous+j] << 1) + abs(v_yprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1));
	  
	  //This first outer loop is used to do the spiral search
	  //We are repeating patterns of moving right,down, left, then up.
	  //At the very end, we go right once more to finish things off.
	  //The algorithm is basically:
	  //right, down, left(m+1), up(m+1).  And then right(m+1) at the very end.
	  for (m = 1; m < v_shift; m+=2)
	  { 
	    //the variable m will tell us how much to shift each time
	    //the variable t is a counter.  if we have to shift 5 times, we will
	    //shift one position at a time and calculate the SAD for each shift.
	    for(t = 0; t < m; t++) 
		{
	      l = l + 1; //m;
		  SAD_value = calculate_SAD(pos_y,pos_x,l,k);
		  //MV_cost = abs((l-new_i) - v_x[new_i*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+new_j]) + abs((k-new_j) - v_y[new_i*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+new_j]);
	      //new neighbor
		  //MV_cost = abs((l-new_i) - v_x[new_i*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+new_j]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j+b_size)]) + abs((k-new_j) - v_y[new_i*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+new_j]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j+b_size)]);
	      //prev level
		  MV_cost = (abs((l-new_j) - (v_xprev[i*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)));
	  	  //MV_cost = (abs((l-new_i) - (v_xprev[i*step_previous+j] << 1)) + abs((l-new_i) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_i) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((l-new_i) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_i) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_j) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_j) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)));
		  if (MV_cost < minMV_cost)
		  {
		    minMV_cost = MV_cost;
		    minMV_SAD = SAD_value;
		    minMV_x = (l-new_j);
		    minMV_y = (k-new_i);
		  }  
          if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;
			motion_pos_x[0] = l; 
            motion_pos_y[0] = k; 
		  }
		}
	    
		for(t = 0; t < m; t++)
		{
	      k = k + 1; //m;
          SAD_value = calculate_SAD(pos_y,pos_x,l,k);
		  //MV_cost = abs((l-new_i) - v_x[new_i*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+new_j]) + abs((k-new_j) - v_y[new_i*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+new_j]);
	      //new neighbor
		  //MV_cost = abs((l-new_i) - v_x[new_i*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+new_j]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j+b_size)]) + abs((k-new_j) - v_y[new_i*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+new_j]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j+b_size)]);
	      //prev level
		  MV_cost = (abs((l-new_j) - (v_xprev[i*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)));
	  	  //MV_cost = (abs((l-new_i) - (v_xprev[i*step_previous+j] << 1)) + abs((l-new_i) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_i) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((l-new_i) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_i) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_j) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_j) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)));
		  if (MV_cost < minMV_cost)
		  {
		    minMV_cost = MV_cost;
		    minMV_SAD = SAD_value;
		    minMV_x = (l-new_j);
		    minMV_y = (k-new_i);
		  }
          if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;
			motion_pos_x[0] = l; 
            motion_pos_y[0] = k; 
		  }
		}
	   
		for(t = 0; t < m+1; t++)
		{
	      l = l - 1; //(m + 1);
          SAD_value = calculate_SAD(pos_y,pos_x,l,k);
		  //MV_cost = abs((l-new_i) - v_x[new_i*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+new_j]) + abs((k-new_j) - v_y[new_i*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+new_j]);
	      //new neighbor
		  //MV_cost = abs((l-new_i) - v_x[new_i*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+new_j]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j+b_size)]) + abs((k-new_j) - v_y[new_i*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+new_j]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j+b_size)]);
	      //prev level
		  MV_cost = (abs((l-new_j) - (v_xprev[i*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)));
	  	  //MV_cost = (abs((l-new_i) - (v_xprev[i*step_previous+j] << 1)) + abs((l-new_i) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_i) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((l-new_i) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_i) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_j) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_j) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)));
		  if (MV_cost < minMV_cost)
		  {
		    minMV_cost = MV_cost;
		    minMV_SAD = SAD_value;
		    minMV_x = (l-new_j);
		    minMV_y = (k-new_i);
		  }
          if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;
			motion_pos_x[0] = l; 
            motion_pos_y[0] = k; 
		  }
		}

		for(t = 0; t < m+1; t++)
		{
	      k = k - 1; //(m + 1);
		  SAD_value = calculate_SAD(pos_y,pos_x,l,k);
		  //MV_cost = abs((l-new_i) - v_x[new_i*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+new_j]) + abs((k-new_j) - v_y[new_i*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+new_j]);
	      //new neighbor
		  //MV_cost = abs((l-new_i) - v_x[new_i*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+new_j]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j+b_size)]) + abs((k-new_j) - v_y[new_i*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+new_j]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j+b_size)]);
	      //prev level
		  MV_cost = (abs((l-new_j) - (v_xprev[i*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)));
	  	  //MV_cost = (abs((l-new_i) - (v_xprev[i*step_previous+j] << 1)) + abs((l-new_i) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_i) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((l-new_i) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_i) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_j) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_j) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)));
		  if (MV_cost < minMV_cost)
		  {
		    minMV_cost = MV_cost;
		    minMV_SAD = SAD_value;
		    minMV_x = (l-new_j);
		    minMV_y = (k-new_i);
		  }
          if (SAD_value < min_SAD_value)
		  {
		    min_SAD_value = SAD_value;
			motion_pos_x[0] = l; 
            motion_pos_y[0] = k; 
		  }
		}
	  }

	  //This is what we do at the end to move across the top row.
      for(t = 0; t < (m-1); t++)
	  {
	    l = l + 1; //m;
	    SAD_value = calculate_SAD(pos_y,pos_x,l,k);
		//MV_cost = abs((l-new_i) - v_x[new_i*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+new_j]) + abs((k-new_j) - v_y[new_i*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+new_j]);
	    //new neighbor
		//MV_cost = abs((l-new_i) - v_x[new_i*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((l-new_i) - v_x[(new_i-b_size)*step+new_j]) + abs((l-new_i) - v_x[(new_i-b_size)*step+(new_j+b_size)]) + abs((k-new_j) - v_y[new_i*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((k-new_j) - v_y[(new_i-b_size)*step+new_j]) + abs((k-new_j) - v_y[(new_i-b_size)*step+(new_j+b_size)]);
	    //prev level
		MV_cost = (abs((l-new_j) - (v_xprev[i*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((l-new_j) - (v_xprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((k-new_i) - (v_yprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)));
	  	//MV_cost = (abs((l-new_i) - (v_xprev[i*step_previous+j] << 1)) + abs((l-new_i) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((l-new_i) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((l-new_i) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((l-new_i) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((k-new_j) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((k-new_j) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((k-new_j) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)));
		if (MV_cost < minMV_cost)
		{
		  minMV_cost = MV_cost;
		  minMV_SAD = SAD_value;
		  minMV_x = (l-new_j);
		  minMV_y = (k-new_i);
		}
        if (SAD_value < min_SAD_value)
		{
		  min_SAD_value = SAD_value;
		  motion_pos_x[0] = l; //+ (b_size_div-1);
          motion_pos_y[0] = k; //+ (b_size_div-1);
		}
	  }
	  
	  v_xtemp = motion_pos_x[0] - new_j;
	  v_ytemp = motion_pos_y[0] - new_i;

	  //Calculate MV cost for the block that has the min SAD
	  //MV_cost = abs((v_xtemp) - v_x[new_i*step+(new_j-b_size)]) + abs((v_xtemp) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((v_xtemp) - v_x[(new_i-b_size)*step+new_j]) + abs((v_ytemp) - v_y[new_i*step+(new_j-b_size)]) + abs((v_ytemp) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((v_ytemp) - v_y[(new_i-b_size)*step+new_j]);
	  //new neighbor
	  //MV_cost = abs((v_xtemp) - v_x[new_i*step+(new_j-b_size)]) + abs((v_xtemp) - v_x[(new_i-b_size)*step+(new_j-b_size)]) + abs((v_xtemp) - v_x[(new_i-b_size)*step+new_j]) + abs((v_xtemp) - v_x[(new_i-b_size)*step+(new_j+b_size)]) + abs((v_ytemp) - v_y[new_i*step+(new_j-b_size)]) + abs((v_ytemp) - v_y[(new_i-b_size)*step+(new_j-b_size)]) + abs((v_ytemp) - v_y[(new_i-b_size)*step+new_j]) + abs((v_ytemp) - v_y[(new_i-b_size)*step+(new_j+b_size)]);
	  //prev level
	  MV_cost = (abs((v_xtemp) - (v_xprev[i*step_previous+j] << 1)) + abs((v_xtemp) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((v_xtemp) - (v_xprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((v_xtemp) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((v_xtemp) - (v_xprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((v_xtemp) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((v_xtemp) - (v_xprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((v_xtemp) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((v_xtemp) - (v_xprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((v_ytemp) - (v_yprev[i*step_previous+j] << 1)) + abs((v_ytemp) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((v_ytemp) - (v_yprev[(i-b_size_previous)*step_previous+(j-b_size_previous)] << 1)) + abs((v_ytemp) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((v_ytemp) - (v_yprev[(i-b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((v_ytemp) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((v_ytemp) - (v_yprev[(i+b_size_previous)*step_previous+(j+b_size_previous)] << 1)) + abs((v_ytemp) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((v_ytemp) - (v_yprev[(i+b_size_previous)*step_previous+(j-b_size_previous)] << 1)));
	  //MV_cost = (abs((v_xtemp) - (v_xprev[i*step_previous+j] << 1)) + abs((v_xtemp) - (v_xprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((v_xtemp) - (v_xprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((v_xtemp) - (v_xprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((v_xtemp) - (v_xprev[(i+b_size_previous)*step_previous+j] << 1)) + abs((v_ytemp) - (v_yprev[i*step_previous+j] << 1)) + abs((v_ytemp) - (v_yprev[i*step_previous+(j-b_size_previous)] << 1)) + abs((v_ytemp) - (v_yprev[(i-b_size_previous)*step_previous+j] << 1)) + abs((v_ytemp) - (v_yprev[i*step_previous+(j+b_size_previous)] << 1)) + abs((v_ytemp) - (v_yprev[(i+b_size_previous)*step_previous+j] << 1)));
	  //'fout << abs(minMV_SAD-min_SAD_value) << std::endl;
	  /*if (abs(minMV_SAD-min_SAD_value) < 20 && level == 1){
	    //{
		  //if ((minMV_x == v_xtemp) && (minMV_y == v_ytemp) && level == 1){
		  //fout << min_SAD_value << std::endl;
		  CvPoint p,q;
		  p.x = new_j - ((b_size/2) - 1);
		  p.y = new_i - ((b_size/2) - 1);
		  q.x = new_j + (b_size/2);
		  q.y = new_i + (b_size/2);
		  cvRectangle(debug, p, q, CV_RGB(117,47,156), 1, 8, 0);
	      //cvCircle(debug, p, 1, CV_RGB(117,47,156), -1);	  
	    }*/
	  	  
	  //total_count++; //for statistics
	  //Choose either the MV that gives min SAD or the MV that gives min cost.
	  if ((minMV_SAD + minMV_cost) < (min_SAD_value + MV_cost)) //equal sign gives preference on MV with lowest cost.
	  {	    
		 //if (!((minMV_SAD == min_SAD_value) && (minMV_cost == MV_cost))) {
		 //fout << "chosen!" << std::endl;
		 //fout << "minMV_SAD is " << minMV_SAD << " and minMV_cost is " << minMV_cost << std::endl;
		 //fout << "min SAD is " << min_SAD_value << " and min SAD cost is " << MV_cost << std::endl;
		 //}
		 //Here we do our assignment for the block with the lowest cost. This is for the MV with min cost
	     v_x[new_i*step+new_j] = minMV_x; //store motion vector x in array
	     v_y[new_i*step+new_j] = minMV_y; //store motion vector y in array
		 //v_x_old[new_i*step+new_j] = minMV_x; //store motion vector x in array -- this is for testing original MV in smoothness fn
	     //v_y_old[new_i*step+new_j] = minMV_y; //store motion vector y in array -- this is for testing original MV in smoothness fn
	  }
	  else
	  {			
	    // Here we do our assignments for the block with the lowest SAD.  This if for MV with min SAD.
	    v_x[new_i*step+new_j] = v_xtemp; //store motion vector x in array
	    v_y[new_i*step+new_j] = v_ytemp; //store motion vector y in array
		//v_x_old[new_i*step+new_j] = v_xtemp; //store motion vector x in array -- this is for testing original MV in smoothness fn
	    //v_y_old[new_i*step+new_j] = v_ytemp; //store motion vector y in array -- this is for testing original MV in smoothness fn
	  }  	 
    }  
  }  
}
void Pair_Motion::setMVs_iter()
{
  int b_size_div = (b_size >> 1);
  //int median_x, median_y; //average of surrounding MVs
  
  //This function is just taking the motion vector of a large block and assigning it to 4 smaller blocks
  //inside this large block.  We incorporate median for determining the motion vectors
 
  //everything in the middle that has surrounding blocks to use for averaging MVs
  //12/1/11 -- Added add_height and add_width to for loops -- wasn't there before.
  for (int i = start_pos; i < height-(add_height - start_pos)/*v_shift*/; i+=b_size)  //goes through all vertical motion positions
  {
    for (int j = start_pos; j < width-(add_width - start_pos)/*v_shift*/; j+=b_size) //goes through all horizontal motion positions
    {
      //int median_array_x[5] = {v_x[i*step_size+j],v_x[(i-block_size)*step_size+j],v_x[(i+block_size)*step_size+j],v_x[i*step_size+(j+block_size)],v_x[i*step_size+(j-block_size)]}; 
	  //int median_array_y[5] = {v_y[i*step_size+j],v_y[(i-block_size)*step_size+j],v_y[(i+block_size)*step_size+j],v_y[i*step_size+(j+block_size)],v_y[i*step_size+(j-block_size)]}; 
	  //median_x = opt_med5(median_array_x); 
	  //median_y = opt_med5(median_array_y); 

	  //median_y = (v_y[i*step_size+j] + v_y[(i-block_size)*step_size+j] + v_y[(i+block_size)*step_size+j] + v_y[i*step_size+(j-block_size)] + v_y[i*step_size+(j+block_size)])/5;
	  //median_x = (v_x[i*step_size+j] + v_x[(i-block_size)*step_size+j] + v_x[(i+block_size)*step_size+j] + v_x[i*step_size+(j-block_size)] + v_x[i*step_size+(j+block_size)])/5;
	  
	  v_x[(i+b_size_div)*step+j] = v_x[i*step+j]; //median_x; //take care of top right 1/4 block
	  v_x[i*step+(j+b_size_div)] = v_x[i*step+j]; //median_x; //take care of bottom left 1/4 block
	  v_x[(i+b_size_div)*step+(j+b_size_div)] = v_x[i*step+j]; //median_x; //take care of bottom right block

	  //old values to prevent propagation...
	  
	  //v_x_old[i*step+j] = v_x[i*step+j];
	  //v_x_old[(i+b_size_div)*step+j] = v_x[i*step+j]; //median_x; //take care of top right 1/4 block
	  //v_x_old[i*step+(j+b_size_div)] = v_x[i*step+j]; //median_x; //take care of bottom left 1/4 block
	  //v_x_old[(i+b_size_div)*step+(j+b_size_div)] = v_x[i*step+j]; //median_x; //take care of bottom right block
	  

	  v_y[(i+b_size_div)*step+j] = v_y[i*step+j]; //median_y; //take care of top right 1/4 block
	  v_y[i*step+(j+b_size_div)] = v_y[i*step+j]; //median_y; //take care of bottom left 1/4 block
	  v_y[(i+b_size_div)*step+(j+b_size_div)] = v_y[i*step+j]; //median_y; //take care of bottom right block

	  //old values to prevent propagation...
	  
	  //v_y_old[i*step+j] = v_y[i*step+j];
	  //v_y_old[(i+b_size_div)*step+j] = v_y[i*step+j]; //median_y; //take care of top right 1/4 block
	  //v_y_old[i*step+(j+b_size_div)] = v_y[i*step+j]; //median_y; //take care of bottom left 1/4 block
	  //v_y_old[(i+b_size_div)*step+(j+b_size_div)] = v_y[i*step+j]; //median_y; //take care of bottom right block
	}
  }
  /*
  //leftmost column
  int i = search_size_div-1;
  for (int j = search_size_div-1; j < height-(search_size_div-1); j+=block_size)  //goes through all vertical motion positions
  {
    v_x[(i+block_size_div)*step_size+j] = v_x[i*step_size+j]; //take care of top right 1/4 block
	v_x[i*step_size+(j+block_size_div)] = v_x[i*step_size+j]; //take care of bottom left 1/4 block
	v_x[(i+block_size_div)*step_size+(j+block_size_div)] = v_x[i*step_size+j]; //take care of bottom right block

	v_y[(i+block_size_div)*step_size+j] = v_y[i*step_size+j]; //take care of top right 1/4 block
	v_y[i*step_size+(j+block_size_div)] = v_y[i*step_size+j]; //take care of bottom left 1/4 block
	v_y[(i+block_size_div)*step_size+(j+block_size_div)] = v_y[i*step_size+j]; //take care of bottom right block
  }

  //top row
  int j = search_size_div-1;
  for (int i = search_size_div-1; i < width-(search_size_div-1); i+=block_size) //goes through all horizontal motion positions
  {
    v_x[(i+block_size_div)*step_size+j] = v_x[i*step_size+j]; //take care of top right 1/4 block
	v_x[i*step_size+(j+block_size_div)] = v_x[i*step_size+j]; //take care of bottom left 1/4 block
	v_x[(i+block_size_div)*step_size+(j+block_size_div)] = v_x[i*step_size+j]; //take care of bottom right block

	v_y[(i+block_size_div)*step_size+j] = v_y[i*step_size+j]; //take care of top right 1/4 block
	v_y[i*step_size+(j+block_size_div)] = v_y[i*step_size+j]; //take care of bottom left 1/4 block
	v_y[(i+block_size_div)*step_size+(j+block_size_div)] = v_y[i*step_size+j]; //take care of bottom right block
  }

  //rightmost column
  i = width-((search_size_div-1) + block_size);
  for (int j = search_size_div-1; j < height-(search_size_div-1); j+=block_size)  //goes through all vertical motion positions
  {
    v_x[(i+block_size_div)*step_size+j] = v_x[i*step_size+j]; //take care of top right 1/4 block
	v_x[i*step_size+(j+block_size_div)] = v_x[i*step_size+j]; //take care of bottom left 1/4 block
	v_x[(i+block_size_div)*step_size+(j+block_size_div)] = v_x[i*step_size+j]; //take care of bottom right block

	v_y[(i+block_size_div)*step_size+j] = v_y[i*step_size+j]; //take care of top right 1/4 block
	v_y[i*step_size+(j+block_size_div)] = v_y[i*step_size+j]; //take care of bottom left 1/4 block
	v_y[(i+block_size_div)*step_size+(j+block_size_div)] = v_y[i*step_size+j]; //take care of bottom right block
  }

  //bottom row
  j = height-((search_size_div-1) + block_size);
  for (int i = search_size_div-1; i < width-(search_size_div-1); i+=block_size)  //goes through all vertical motion positions
  {
    v_x[(i+block_size_div)*step_size+j] = v_x[i*step_size+j]; //take care of top right 1/4 block
	v_x[i*step_size+(j+block_size_div)] = v_x[i*step_size+j]; //take care of bottom left 1/4 block
	v_x[(i+block_size_div)*step_size+(j+block_size_div)] = v_x[i*step_size+j]; //take care of bottom right block

	v_y[(i+block_size_div)*step_size+j] = v_y[i*step_size+j]; //take care of top right 1/4 block
	v_y[i*step_size+(j+block_size_div)] = v_y[i*step_size+j]; //take care of bottom left 1/4 block
	v_y[(i+block_size_div)*step_size+(j+block_size_div)] = v_y[i*step_size+j]; //take care of bottom right block
  } */
}
void Pair_Motion::copyto_oldMVs()
{  
  
  //This function is just copying the current motion vectors into the old motion vector arrays so that we can make prevent bad
  //motion vectors from propagating -- we always check to see if the original (old) MVs in the neighborhood are a better choice
  //in the smoothness function.  

  for (int i = start_pos; i <= height-(start_pos + b_size)/*v_shift*/; i+=b_size)  //goes through all vertical motion positions
  {
    for (int j = start_pos; j <= width-(start_pos + b_size)/*v_shift*/; j+=b_size) //goes through all horizontal motion positions
    {
	  //old values to prevent propagation...
	  //v_x_old[i*step+j] = v_x[i*step+j];	  
	  //v_y_old[i*step+j] = v_y[i*step+j];  

	}
  }  
}
int Pair_Motion::calculate_SAD(int y, int x, int l, int k)
{
  int SAD_value = 0;
  int reset_l_value = l;
  
  for (int i = y; i < y + b_size; ++i) //this is the height of the block
  {
    for (int j = x; j < x + b_size; ++j) //this is the width of the block
	{
	  SAD_value += abs(data_b[i*step+j] - data_a[k*step+l]); //subtract current frame pixel from previous frame pixel
	  ++l; //increment horizontal position of pixel in previous block
	}			
	l = reset_l_value;
	++k; //increment vertical position of pixel in previous block
  }
  
  return SAD_value;
}
double Pair_Motion::calculate_SAD2(int y, int x, int l, int k)
{
  double SAD_value = 0;
  x = x - start_pos;
  y = y - start_pos;
  l = l - start_pos;
  k = k - start_pos;
  int reset_l_value = l; 
  
  for (int i = y; i < y + b_size; i++) //this is the height of the block
  {
    for (int j = x; j < x + b_size; j++) //this is the width of the block
	{
	  SAD_value = SAD_value + ((uchar)level1b_grey->imageData[i*level1b_grey->widthStep+j] - (uchar)level1a_grey->imageData[k*level1a_grey->widthStep+l])*((uchar)level1b_grey->imageData[i*level1b_grey->widthStep+j] - (uchar)level1a_grey->imageData[k*level1a_grey->widthStep+l]);
		  //(image2_data[i*step_inputs+j] - image1_data[k*step_inputs+l])*(image2_data[i*step_inputs+j] - image1_data[k*step_inputs+l]); //subtract current frame pixel from previous frame pixel
	  l++; //increment horizontal position of pixel in previous block
	}			
	l = reset_l_value;
	k++; //increment vertical position of pixel in previous block
  }  

  return SAD_value;
 
}
int Pair_Motion::calculate_local_SAD(int y, int x, int l, int k)
{
  int SAD_value = 0;
  int reset_l_value = l;

  for (int i = y; i < y + b_size; i++) //this is the height of the block
  {
    for (int j = x; j < x + b_size; j++) //this is the width of the block
	{
	  SAD_value = SAD_value + abs(data_b[i*step+j] - data_b[k*step+l]); //subtract current frame pixel from previous frame pixel
	  l++; //increment horizontal position of pixel in previous block
	}			
	l = reset_l_value;
	k++; //increment vertical position of pixel in previous block
  }
  
  return SAD_value;
}
void Pair_Motion::refine_MVs()
{
  int min_index, i, j, k, l, n;
  double lambda = 8; //4th iteration, bsize = 2;
  double candidate[12];  
  
  //The variables below are for removing duplicate MVs
  int duplicate_list[12] = {0};
  int hash[12] = {0};
  int MV_x[12];
  int MV_y[12];
    
  b_size = 2;
  int bs_squared = b_size;
  int total_overlap;
  int overlap_amount = b_size*b_size;
  int overlap_skip = overlap_amount + b_size*(b_size/2);  
  	  
  for (i = start_pos + b_size; i < height-(add_height - start_pos + b_size); i+=b_size) //goes through all vertical pixels
  {
	for (j = start_pos + b_size; j < width-(add_width - start_pos + b_size); j+=b_size) //goes through all horizontal pixels
    {	   

	  total_overlap = get_overlap_volume(v_x[i*step+j] + j, v_y[i*step+j] + i);
	  if (total_overlap < (overlap_skip))
	    continue;
	  else
	  {

	  //Initalize duplicate list to zero -- no duplicates
	  for(k = 0; k < 12; k++)
	  {
       duplicate_list[k] = 0;  
	   hash[k] = 0;
	  }
	  
	  MV_x[0] = v_x[i*step+j];
	  MV_y[0] = v_y[i*step+j];
	  MV_x[1] = v_x[(i-b_size)*step+j];
	  MV_y[1] = v_y[(i-b_size)*step+j];
	  MV_x[2] = v_x[i*step+(j-b_size)];
	  MV_y[2] = v_y[i*step+(j-b_size)];
	  MV_x[3] = v_x[(i+b_size)*step+j];
	  MV_y[3] = v_y[(i+b_size)*step+j];
	  MV_x[4] = v_x[i*step+(j+b_size)];
	  MV_y[4] = v_y[i*step+(j+b_size)];
	  MV_x[5] = v_x[(i-b_size)*step+(j-b_size)];
	  MV_y[5] = v_y[(i-b_size)*step+(j-b_size)];
	  MV_x[6] = v_x[(i+b_size)*step+(j-b_size)];
	  MV_y[6] = v_y[(i+b_size)*step+(j-b_size)];
	  MV_x[7] = v_x[(i-b_size)*step+(j+b_size)];
	  MV_y[7] = v_y[(i-b_size)*step+(j+b_size)];	  
	  MV_x[8] = v_x[(i+b_size)*step+(j+b_size)];
	  MV_y[8] = v_y[(i+b_size)*step+(j+b_size)];
	  MV_x[9] = v_x[i*step+j] + 1;
	  MV_y[9] = v_y[i*step+j] + 1;
	  MV_x[10] = v_x[i*step+j] - 1;
	  MV_y[10] = v_y[i*step+j] - 1;
	  MV_x[11] = v_x[i*step+j] + 2;
	  MV_y[11] = v_y[i*step+j] + 2;
	  MV_x[12] = v_x[i*step+j] - 2;
	  MV_y[12] = v_y[i*step+j] - 2;

	  //Check for duplicates, a duplicate is assigned a value of '1' in the duplicate list
      for(k = 1; k < 12; k++)
      {
        for(l = 0; l < 12; l++)
	    {
	      if ((l == k) || (duplicate_list[k] == 1))
	        continue;
          else if ((MV_x[k] == MV_x[l]) && (MV_y[k] == MV_y[l]))
	      {
	        duplicate_list[l] = 1;
	      } 	    
	    }
      }	  
	 
	  //Subtract out the overlap for the one already assigned prior to smoothness function.
	  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	    for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
		  Computed_Data.overlap[k*step+l] -= 1;

	  n = 0;
	  for(k = 0; k < 12; k++)
	  {
		if (duplicate_list[k] == 0)
		{
		  hash[n] = k; // this tell us which position of MV_x and MV_y the candidate is stored at.
	      candidate[n++] = (1 + calculate_SAD(i, j, MV_x[k] + j, MV_y[k] + i))*(1 + ((double)get_overlap_volume(MV_x[k] + j, MV_y[k] + i)/bs_squared)) + lambda*cost_fn(MV_y[k], MV_x[k], i, j);
		}
	  }

	  min_array(candidate, n, min_index); 	
	  v_x[i*step+j] = MV_x[hash[min_index]];
	  v_y[i*step+j] = MV_y[hash[min_index]];

	  //Next, add the overlap for the minimum found from min_array() above.
	  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	    for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
		  Computed_Data.overlap[k*step+l] += 1;    
	  }
	}	
  }  
} 
void Pair_Motion::add_smoothness8_old(double lambda_value) 
{
  int min_index, i, j;//, temp_x, temp_y;
  //double min_value;
  //double lambda = 20000; //should probably depend on pixel variance.  high variance = more weight on lambda. 
  double lambda;//, largest_val;
  double candidate[18];  
  double SAD_temp[9];
  double cost_temp[9];
  //int choice_x[9];
  //int choice_y[9];
  //double xs_temp, numerator, denominator;
  //double SAD_max, SAD_min, cost_max, cost_min, denominator;
  //double K = ((b_size*b_size)*64);//10000);//12500);
  //int boundary_value;

  lambda = lambda_value;
	  
  for (i = start_pos + b_size; i < height-(add_height - start_pos + b_size); i+=b_size)//+)  //goes through all vertical pixels
  {
	for (j = start_pos + b_size; j < width-(add_width - start_pos + b_size); j+=b_size)//+)  //goes through all horizontal pixels
    {
	        
      //lambda = lambda_value;	
	  
	  /*spatial_test[0] = (calculate_local_SAD(i, j, i-b_size, j));
	  spatial_test[1] = (calculate_local_SAD(i, j, i, j-b_size));
	  spatial_test[2] = (calculate_local_SAD(i, j, i+b_size, j));
	  spatial_test[3] = (calculate_local_SAD(i, j, i, j+b_size));
	  spatial_test[4] = (calculate_local_SAD(i, j, i-b_size, j-b_size));
	  spatial_test[5] = (calculate_local_SAD(i, j, i+b_size, j-b_size));
	  spatial_test[6] = (calculate_local_SAD(i, j, i-b_size, j+b_size));
	  spatial_test[7] = (calculate_local_SAD(i, j, i+b_size, j+b_size));	
	  //spatial_test[8] = (calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j));
	  
	  largest_val = K;//max_array_val(spatial_test, 9);*/

	  //lambda_value = 55;

	  /*lambda_test[0] = lambda_value;// + (largest_val/(1+calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j)));
	  lambda_test[1] = lambda_value + (largest_val/(1+calculate_local_SAD(i, j, i-b_size, j)));
	  lambda_test[2] = lambda_value + (largest_val/(1+calculate_local_SAD(i, j, i, j-b_size)));
	  lambda_test[3] = lambda_value + (largest_val/(1+calculate_local_SAD(i, j, i+b_size, j)));
	  lambda_test[4] = lambda_value + (largest_val/(1+calculate_local_SAD(i, j, i, j+b_size)));
	  lambda_test[5] = lambda_value + (largest_val/(1+calculate_local_SAD(i, j, i-b_size, j-b_size)));
	  lambda_test[6] = lambda_value + (largest_val/(1+calculate_local_SAD(i, j, i+b_size, j-b_size)));
	  lambda_test[7] = lambda_value + (largest_val/(1+calculate_local_SAD(i, j, i-b_size, j+b_size)));
	  lambda_test[8] = lambda_value + (largest_val/(1+calculate_local_SAD(i, j, i+b_size, j+b_size)));*/

	  //if(level == 3 || level == 2)
	  //lambda = (lambda_value);// << 1);// << 1); //minus 3 worked well for the foreman sequence.
	  //else
	  //lambda = calculate_local_SAD(i, j, i-b_size, j) + calculate_local_SAD(i, j, i, j-b_size) + calculate_local_SAD(i, j, i+b_size, j) + calculate_local_SAD(i, j, i, j+b_size) + calculate_local_SAD(i, j, i-b_size, j-b_size) + calculate_local_SAD(i, j, i+b_size, j-b_size) + calculate_local_SAD(i, j, i-b_size, j+b_size) + calculate_local_SAD(i, j, i+b_size, j+b_size);
	  	  
	  SAD_temp[0] = calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
	  cost_temp[0] = cost_fn8(i, j, i, j);	
	  //choice_x[0] = v_x[i*step+j];
	  //choice_y[0] = v_y[i*step+j];
	  SAD_temp[1] = calculate_SAD(i, j, v_x[(i-b_size)*step+j] + j, v_y[(i-b_size)*step+j] + i); 
	  cost_temp[1] = cost_fn8a(i-b_size,j, i, j);
	  //choice_x[1] = v_x[(i-b_size)*step+j];
	  //choice_y[1] = v_y[(i-b_size)*step+j];
	  SAD_temp[2] = calculate_SAD(i, j, v_x[i*step+(j-b_size)] + j, v_y[i*step+(j-b_size)] + i); 
	  cost_temp[2] = cost_fn8b(i,j-b_size, i, j); 	
	  //choice_x[2] = v_x[i*step+(j-b_size)];
	  //choice_y[2] = v_y[i*step+(j-b_size)];
	  SAD_temp[3] = calculate_SAD(i, j, v_x[(i+b_size)*step+j] + j, v_y[(i+b_size)*step+j] + i); 
	  cost_temp[3] = cost_fn8c(i+b_size,j, i, j); 	 
	  //choice_x[3] = v_x[(i+b_size)*step+j];
	  //choice_y[3] = v_y[(i+b_size)*step+j];
	  SAD_temp[4] = calculate_SAD(i, j, v_x[i*step+(j+b_size)] + j, v_y[i*step+(j+b_size)] + i); 
	  cost_temp[4] = cost_fn8d(i,j+b_size, i, j); 	 
	  //choice_x[4] = v_x[i*step+(j+b_size)];
	  //choice_y[4] = v_y[i*step+(j+b_size)];
	  SAD_temp[5] = calculate_SAD(i, j, v_x[(i-b_size)*step+(j-b_size)] + j, v_y[(i-b_size)*step+(j-b_size)] + i); 
	  cost_temp[5] = cost_fn8e(i-b_size,j-b_size, i, j); 
	  //choice_x[5] = v_x[(i-b_size)*step+(j-b_size)];
	  //choice_y[5] = v_y[(i-b_size)*step+(j-b_size)];
	  SAD_temp[6] = calculate_SAD(i, j, v_x[(i+b_size)*step+(j-b_size)] + j, v_y[(i+b_size)*step+(j-b_size)] + i); 
	  cost_temp[6] = cost_fn8f(i+b_size,j-b_size, i, j); 	 
	  //choice_x[6] = v_x[(i+b_size)*step+(j-b_size)];
	  //choice_y[6] = v_y[(i+b_size)*step+(j-b_size)];
	  SAD_temp[7] = calculate_SAD(i, j, v_x[(i-b_size)*step+(j+b_size)] + j, v_y[(i-b_size)*step+(j+b_size)] + i); 
	  cost_temp[7] = cost_fn8g(i-b_size,j+b_size, i, j); 	  
	  //choice_x[7] = v_x[(i-b_size)*step+(j+b_size)];
	  //choice_y[7] = v_y[(i-b_size)*step+(j+b_size)];
	  SAD_temp[8] = calculate_SAD(i, j, v_x[(i+b_size)*step+(j+b_size)] + j, v_y[(i+b_size)*step+(j+b_size)] + i); 
	  cost_temp[8] = cost_fn8h(i+b_size,j+b_size, i, j); 
	  //choice_x[8] = v_x[(i+b_size)*step+(j+b_size)];
	  //choice_y[8] = v_y[(i+b_size)*step+(j+b_size)];

	  //These are using the old MVs -- to make sure that we haven't eliminated any good motion vectors by propagating.
	  /*SAD_temp[9] = calculate_SAD(i, j, v_x_old[(i-b_size)*step+j] + i, v_y_old[(i-b_size)*step+j] + j); 
	  cost_temp[9] = cost_fn8a_old(i-b_size,j, i, j);
	  SAD_temp[10] = calculate_SAD(i, j, v_x_old[i*step+(j-b_size)] + i, v_y_old[i*step+(j-b_size)] + j); 
	  cost_temp[10] = cost_fn8b_old(i,j-b_size, i, j); 	
	  SAD_temp[11] = calculate_SAD(i, j, v_x_old[(i-b_size)*step+(j-b_size)] + i, v_y_old[(i-b_size)*step+(j-b_size)] + j); 
	  cost_temp[11] = cost_fn8e_old(i-b_size,j-b_size, i, j); 
	  SAD_temp[12] = calculate_SAD(i, j, v_x_old[(i-b_size)*step+(j+b_size)] + i, v_y_old[(i-b_size)*step+(j+b_size)] + j); 
	  cost_temp[12] = cost_fn8g_old(i-b_size,j+b_size, i, j); 	*/
	  
	  //Median Filtering
	  //int med_array_x[9] = {SAD_temp[0], SAD_temp[1], SAD_temp[2], SAD_temp[3], SAD_temp[4], SAD_temp[5], SAD_temp[6], SAD_temp[7], SAD_temp[8]};
	  //{v_x[i*step+j],v_x[(i-b_size)*step+j],v_x[i*step+(j-b_size)],v_x[(i+b_size)*step+j],v_x[i*step+(j+b_size)],v_x[(i-b_size)*step+(j-b_size)],v_x[(i+b_size)*step+(j-b_size)],v_x[(i-b_size)*step+(j+b_size)],v_x[(i+b_size)*step+(j+b_size)]};
      //int med_array_y[9] = {v_y[i*step+j],v_y[(i-b_size)*step+j],v_y[i*step+(j-b_size)],v_y[(i+b_size)*step+j],v_y[i*step+(j+b_size)],v_y[(i-b_size)*step+(j-b_size)],v_y[(i+b_size)*step+(j-b_size)],v_y[(i-b_size)*step+(j+b_size)],v_y[(i+b_size)*step+(j+b_size)]};
      //int med_x = opt_med9(med_array_x);
      //int med_y = opt_med9(med_array_y);

	  //if (twomin_diff(SAD_temp,9) > 2)
	    //lambda = 1;//0;
	  //else
	    //lambda = lambda_value; // (b_size << 1);
	  //lambda = 1 / (1 + twomin_diff(SAD_temp, 9));
	  
	  //best linear line of fit 
	  /*xs_temp = SAD_temp[0]+SAD_temp[1]+SAD_temp[2]+SAD_temp[3]+SAD_temp[4]+SAD_temp[5]+SAD_temp[6]+SAD_temp[7]+SAD_temp[8];
	  numerator = 9*(SAD_temp[0]*cost_temp[0] + SAD_temp[1]*cost_temp[1] + SAD_temp[2]*cost_temp[2] + SAD_temp[3]*cost_temp[3] + SAD_temp[4]*cost_temp[4] + SAD_temp[5]*cost_temp[5] + SAD_temp[6]*cost_temp[6] + SAD_temp[7]*cost_temp[7] + SAD_temp[8]*cost_temp[8])-xs_temp*(cost_temp[0]+cost_temp[1]+cost_temp[2]+cost_temp[3]+cost_temp[4]+cost_temp[5]+cost_temp[6]+cost_temp[7]+cost_temp[8]);
	  denominator = 9*(SAD_temp[0]*SAD_temp[0]+SAD_temp[1]*SAD_temp[1]+SAD_temp[2]*SAD_temp[2]+SAD_temp[3]*SAD_temp[3]+SAD_temp[4]*SAD_temp[4]+SAD_temp[5]*SAD_temp[5]+SAD_temp[6]*SAD_temp[6]+SAD_temp[7]*SAD_temp[7]+SAD_temp[8]*SAD_temp[8]) - xs_temp*xs_temp;
	  
	  if (denominator == 0)
	    lambda = (b_size << 2);
	  else
	    lambda = max((b_size << 2), (-1*(numerator/denominator)));	  */

	  //max min try
	  /*SAD_max = max_array_val(SAD_temp, 9);
	  SAD_min = min_array_val(SAD_temp, 9);
	  cost_max = max_array_val(cost_temp, 9);
	  cost_min = min_array_val(cost_temp, 9);
	  denominator = cost_max - cost_min;

	  if (denominator == 0)
	    lambda = (b_size << 2);
	  else
	    lambda = (SAD_max - SAD_min)/(denominator);//max((b_size << 1), (SAD_max - SAD_min) / (cost_max - cost_min));*/

	  //lambda = lambda_value; //min((b_size << 3), SAD_temp[0] + SAD_temp[1] + SAD_temp[2] + SAD_temp[3] + SAD_temp[4] + SAD_temp[5] + SAD_temp[6] + SAD_temp[7] + SAD_temp[8]);
	  //lambda = (b_size << 3);
	  //fout << "ratio is " << -1*(numerator/denominator) << std::endl;
	  //fout << "lambda value is " << lambda << std::endl;

	  //int med_array_x[5] = {v_x[i*step+j],v_x[(i-b_size)*step+j],v_x[i*step+(j-b_size)],v_x[(i+b_size)*step+j],v_x[i*step+(j+b_size)]};
      //int med_array_y[5] = {v_y[i*step+j],v_y[(i-b_size)*step+j],v_y[i*step+(j-b_size)],v_y[(i+b_size)*step+j],v_y[i*step+(j+b_size)]};
      //int med_x = opt_med5(med_array_x);
      //int med_y = opt_med5(med_array_y);
      
	  //candidate[0] = calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j) + lambda*(abs(v_x[i*step+j] - med_x) + abs(v_y[i*step+j] - med_y));
	  //candidate[1] = calculate_SAD(i, j, v_x[(i-b_size)*step+j] + i, v_y[(i-b_size)*step+j] + j) + lambda*(abs(v_x[(i-b_size)*step+j] - med_x) + abs(v_y[(i-b_size)*step+j] - med_y));
	  //candidate[2] = calculate_SAD(i, j, v_x[i*step+(j-b_size)] + i, v_y[i*step+(j-b_size)] + j) + lambda*(abs(v_x[i*step+(j-b_size)] - med_x) + abs(v_y[i*step+(j-b_size)] - med_y));
	  //candidate[3] = calculate_SAD(i, j, v_x[(i+b_size)*step+j] + i, v_y[(i+b_size)*step+j] + j) + lambda*(abs(v_x[(i+b_size)*step+j] - med_x) + abs(v_y[(i+b_size)*step+j] - med_y));
	  //candidate[4] = calculate_SAD(i, j, v_x[i*step+(j+b_size)] + i, v_y[i*step+(j+b_size)] + j) + lambda*(abs(v_x[i*step+(j+b_size)] - med_x) + abs(v_y[i*step+(j+b_size)] - med_y));
	  //candidate[5] = calculate_SAD(i, j, v_x[(i-b_size)*step+(j-b_size)] + i, v_y[(i-b_size)*step+(j-b_size)] + j) + lambda*(abs(v_x[(i-b_size)*step+(j-b_size)] - med_x) + abs(v_y[(i-b_size)*step+(j-b_size)] - med_y));
	  //candidate[6] = calculate_SAD(i, j, v_x[(i+b_size)*step+(j-b_size)] + i, v_y[(i+b_size)*step+(j-b_size)] + j) + lambda*(abs(v_x[(i+b_size)*step+(j-b_size)] - med_x) + abs(v_y[(i+b_size)*step+(j-b_size)] - med_y));
	  //candidate[7] = calculate_SAD(i, j, v_x[(i-b_size)*step+(j+b_size)] + i, v_y[(i-b_size)*step+(j+b_size)] + j) + lambda*(abs(v_x[(i-b_size)*step+(j+b_size)] - med_x) + abs(v_y[(i-b_size)*step+(j+b_size)] - med_y));
	  //candidate[8] = calculate_SAD(i, j, v_x[(i+b_size)*step+(j+b_size)] + i, v_y[(i+b_size)*step+(j+b_size)] + j) + lambda*(abs(v_x[(i+b_size)*step+(j+b_size)] - med_x) + abs(v_y[(i+b_size)*step+(j+b_size)] - med_y));
	  	   
	  candidate[0] = SAD_temp[0] + lambda*cost_temp[0];
	  candidate[1] = SAD_temp[1] + lambda*cost_temp[1];
	  candidate[2] = SAD_temp[2] + lambda*cost_temp[2];
	  candidate[3] = SAD_temp[3] + lambda*cost_temp[3];
	  candidate[4] = SAD_temp[4] + lambda*cost_temp[4];
	  candidate[5] = SAD_temp[5] + lambda*cost_temp[5];
	  candidate[6] = SAD_temp[6] + lambda*cost_temp[6];
	  candidate[7] = SAD_temp[7] + lambda*cost_temp[7];
	  candidate[8] = SAD_temp[8] + lambda*cost_temp[8]; 
	  
	  /*
	  double local_min = candidate[0];
	  int local_min_index = 0;
	  for(int z = 1; z < 9; z++)
	  {
        if (candidate[z] < local_min)
		{
		  local_min = candidate[z];
		  local_min_index = z;
		}
	  }

	  for(int z = 0; z < 9; z++)
	  {
		if (z == local_min_index)
		  continue;
	    else if (candidate[z] == candidate[local_min_index])
		{
           if ((SAD_temp[z] == SAD_temp[local_min_index]) && (cost_temp[z] == cost_temp[local_min_index]))
		     total_count = total_count;
		   else
		     total_count++;
		}
	  } */



	  //Candidates to prevent elimination of MV when propagating
	  /*candidate[9] = SAD_temp[9] + lambda*cost_temp[9]; 
	  candidate[10] = SAD_temp[10] + lambda*cost_temp[10]; 
	  candidate[11] = SAD_temp[11] + lambda*cost_temp[11]; 
	  candidate[12] = SAD_temp[12] + lambda*cost_temp[12]; */

	  //fout << "Lambda value is " << lambda << std::endl;
	  /*if (i >= 373 && i <= 400 && j >= 166 && j <= 195 && b_size == 8){
	  fout << "For position [" << i << ":" << j << "] here are the SADs and SAD + Costs:" << std::endl;
	  fout << cost_temp[0] << " " << SAD_temp[0] << std::endl;
	  fout << cost_temp[1] << " " << SAD_temp[1] << std::endl;
	  fout << cost_temp[2] << " " << SAD_temp[2] << std::endl;
	  fout << cost_temp[3] << " " << SAD_temp[3] << std::endl;
	  fout << cost_temp[4] << " " << SAD_temp[4] << std::endl;
	  fout << cost_temp[5] << " " << SAD_temp[5] << std::endl;
	  fout << cost_temp[6] << " " << SAD_temp[6] << std::endl;
	  fout << cost_temp[7] << " " << SAD_temp[7] << std::endl;
	  fout << cost_temp[8] << " " << SAD_temp[8] << std::endl;	  
	  
	  //fout << "Candidate[0] is " << candidate[0] << std::endl;
	  //fout << "Costs" << std::endl;
	  //fout << cost_temp[0] << std::endl;
	  //fout << cost_temp[1] << std::endl;
	  //fout << cost_temp[2] << std::endl;
	  //fout << cost_temp[3] << std::endl;
	  //fout << cost_temp[4] << std::endl;
	  //fout << cost_temp[5] << std::endl;
	  //fout << cost_temp[6] << std::endl;
	  //fout << cost_temp[7] << std::endl;
	  //fout << cost_temp[8] << std::endl;

	  //fout << "Candidate[1] is " << candidate[1] << std::endl;
	  
	  fout << " " << std::endl;
	  fout << " " << std::endl; 
	  } */

	  /*if (level != 3)
	  {
		//std::cout << "level is " << level << std::endl;
	    temp_x = v_x[i*step+j];
	    temp_y = v_y[i*step+j];
	    v_x[i*step+j] = (v_xprev[((i>>1)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);
	    v_y[i*step+j] = (v_yprev[((i>>1)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);	
	    candidate[9] = calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j) + lambda*cost_fn8(i, j, i, j);

		v_x[i*step+j] = (v_xprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);
	    v_y[i*step+j] = (v_yprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);
	    candidate[10] = calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j) + lambda*cost_fn8(i, j, i, j);

		v_x[i*step+j] = (v_xprev[((i>>1)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);
	    v_y[i*step+j] = (v_yprev[((i>>1)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);	
	    candidate[11] = calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j) + lambda*cost_fn8(i, j, i, j);

		v_x[i*step+j] = (v_xprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);
	    v_y[i*step+j] = (v_yprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);	
	    candidate[12] = calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j) + lambda*cost_fn8(i, j, i, j);

		v_x[i*step+j] = (v_xprev[((i>>1)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);
	    v_y[i*step+j] = (v_yprev[((i>>1)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);	
	    candidate[13] = calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j) + lambda*cost_fn8(i, j, i, j);

		v_x[i*step+j] = (v_xprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);
	    v_y[i*step+j] = (v_yprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);	
	    candidate[14] = calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j) + lambda*cost_fn8(i, j, i, j);

		v_x[i*step+j] = (v_xprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);
	    v_y[i*step+j] = (v_yprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);	
	    candidate[15] = calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j) + lambda*cost_fn8(i, j, i, j);

		v_x[i*step+j] = (v_xprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);
	    v_y[i*step+j] = (v_yprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);		
	    candidate[16] = calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j) + lambda*cost_fn8(i, j, i, j);

		v_x[i*step+j] = (v_xprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);
	    v_y[i*step+j] = (v_yprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);		
	    candidate[17] = calculate_SAD(i, j, v_x[i*step+j] + i, v_y[i*step+j] + j) + lambda*cost_fn8(i, j, i, j);

		v_x[i*step+j] = temp_x;
		v_y[i*step+j] = temp_y;

	    min_array(candidate, 18, min_index); //find the minimum value and array index of this value
	  }
	  else	*/
	    //min_array(candidate, 9, min_index);
	  
	  /*if (level == 1 && b_size == 2) //was b_size = 4
	  {
        if((min_array(candidate, 9, min_index) - lambda*cost_temp[min_index]) > reliability_threshold) //find the minimum value and array index of this value
	    {		
		//for(int x=i; x < i+b_size; x++)
		//{
		  //for(int y = j; y < j+b_size; y++)
		  //{
		    //v_x[i*step+j] = 0; //there is no motion vector due to occlusion.
			//v_y[i*step+j] = 0;
		  //}
		//}
		//std::cout << "Skipped" << std::endl;
	    //continue;
          for(int y = i; y < i+b_size; y++)
		    for(int x = j; x < j+b_size; x++)
			 Computed_Data.reliability[y*step+x] = 0;
	    }
	    //else	
	      //min_array(candidate, 9, min_index);
	  }
	  else	*/
	    min_array(candidate, 9, min_index);
		//v_x[i*step+j] = choice_x[min_index];
	    //v_y[i*step+j] = choice_y[min_index];

	  //this is a bad way to do this -- need to make this faster by using some clever indexing
	  switch (min_index) {
        case 0:
			//fout << "case 0 chosen" << std::endl;
		  //if (level != 3)
		  //{
            //v_x[i*step+j] = temp_x;
	        //v_y[i*step+j] = temp_y;
		  //}
          break;
		case 1:
			//fout << "case 1 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x[(i-b_size)*step+j];
	      v_y[i*step+j] = v_y[(i-b_size)*step+j];
		  break;
		case 2:
			//fout << "case 2 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x[i*step+(j-b_size)];
	      v_y[i*step+j] = v_y[i*step+(j-b_size)];
		  break;
		case 3:
			//fout << "case 3 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x[(i+b_size)*step+j];
	      v_y[i*step+j] = v_y[(i+b_size)*step+j];
		  break;
		case 4:
			//fout << "case 4 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x[i*step+(j+b_size)];
	      v_y[i*step+j] = v_y[i*step+(j+b_size)];
		  break;
		case 5:
			//fout << "case 5 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x[(i-b_size)*step+(j-b_size)];
	      v_y[i*step+j] = v_y[(i-b_size)*step+(j-b_size)];
		  break;
		case 6:
			//fout << "case 6 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x[(i+b_size)*step+(j-b_size)];
	      v_y[i*step+j] = v_y[(i+b_size)*step+(j-b_size)];
		  break;
	    case 7:
			//fout << "case 7 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x[(i-b_size)*step+(j+b_size)];
	      v_y[i*step+j] = v_y[(i-b_size)*step+(j+b_size)];
		  break;
	    case 8:
			//fout << "case 8 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x[(i+b_size)*step+(j+b_size)];
	      v_y[i*step+j] = v_y[(i+b_size)*step+(j+b_size)];
		  break; 
        //To prevent propagation errros...
		/*  
	    case 9:
			//fout << "case 9 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x_old[(i-b_size)*step+j];
	      v_y[i*step+j] = v_y_old[(i-b_size)*step+j];
		  
		  break;
	    case 10:
			//fout << "case 10 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x_old[i*step+(j-b_size)];
	      v_y[i*step+j] = v_y_old[i*step+(j-b_size)];
		  break;
	    case 11:
			//fout << "case 11 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x_old[(i-b_size)*step+(j-b_size)];
	      v_y[i*step+j] = v_y_old[(i-b_size)*step+(j-b_size)];
		  break;
	    case 12:
			//fout << "case 12 chosen" << std::endl;
		  //v_x_old[i*step+j] = v_x[i*step+j];
		  //v_y_old[i*step+j] = v_y[i*step+j];
		  v_x[i*step+j] = v_x_old[(i-b_size)*step+(j+b_size)];
	      v_y[i*step+j] = v_y_old[(i-b_size)*step+(j+b_size)];
		  break;
		  */
        /*
		case 9:
		  v_x[i*step+j] = (v_xprev[((i>>1)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);
	      v_y[i*step+j] = (v_yprev[((i>>1)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);	
		  break;
		case 10:
		  v_x[i*step+j] = (v_xprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);
	      v_y[i*step+j] = (v_yprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);
		  break;
		case 11:
		  v_x[i*step+j] = (v_xprev[((i>>1)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);
	      v_y[i*step+j] = (v_yprev[((i>>1)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);
		  break;
		case 12:
		  v_x[i*step+j] = (v_xprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);
	      v_y[i*step+j] = (v_yprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+((j>>1)+b_size_previous)] << 1);
		  break;
		case 13:
		  v_x[i*step+j] = (v_xprev[((i>>1)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);
	      v_y[i*step+j] = (v_yprev[((i>>1)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);
		  break;
		case 14:
		  v_x[i*step+j] = (v_xprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);
	      v_y[i*step+j] = (v_yprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);	
		  break;
		case 15:
		  v_x[i*step+j] = (v_xprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);
	      v_y[i*step+j] = (v_yprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+(((j>>1)-b_size_previous)+b_size_previous)] << 1);	
		  break;
		case 16:
		  v_x[i*step+j] = (v_xprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);
	      v_y[i*step+j] = (v_yprev[(((i>>1)-b_size_previous)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);	
		  break;
	    case 17:
		  v_x[i*step+j] = (v_xprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);
	      v_y[i*step+j] = (v_yprev[(((i>>1)+b_size_previous)+b_size_previous)*step_previous+(((j>>1)+b_size_previous)+b_size_previous)] << 1);
		  break; */
	  }      
	}
  }
  
  //Top row
  i = start_pos;
  for (j = start_pos + b_size; j < width-(add_width - start_pos + b_size); j+=b_size)
  {    
    SAD_temp[0] = calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
	cost_temp[0] = cost_fn8top(i, j, i, j);		
		 
	SAD_temp[1] = calculate_SAD(i, j, v_x[i*step+(j-b_size)] + j, v_y[i*step+(j-b_size)] + i); 
	cost_temp[1] = cost_fn8top(i,j-b_size, i, j); 	
	
	SAD_temp[2] = calculate_SAD(i, j, v_x[(i+b_size)*step+j] + j, v_y[(i+b_size)*step+j] + i); 
	cost_temp[2] = cost_fn8top(i+b_size,j, i, j); 	 
	
	SAD_temp[3] = calculate_SAD(i, j, v_x[i*step+(j+b_size)] + j, v_y[i*step+(j+b_size)] + i); 
	cost_temp[3] = cost_fn8top(i,j+b_size, i, j); 	 	
	 
	SAD_temp[4] = calculate_SAD(i, j, v_x[(i+b_size)*step+(j-b_size)] + j, v_y[(i+b_size)*step+(j-b_size)] + i); 
	cost_temp[4] = cost_fn8top(i+b_size,j-b_size, i, j); 	    
	 
	SAD_temp[5] = calculate_SAD(i, j, v_x[(i+b_size)*step+(j+b_size)] + j, v_y[(i+b_size)*step+(j+b_size)] + i); 
	cost_temp[5] = cost_fn8top(i+b_size,j+b_size, i, j); 
	 
	candidate[0] = SAD_temp[0] + lambda*cost_temp[0];
	candidate[1] = SAD_temp[1] + lambda*cost_temp[1];
	candidate[2] = SAD_temp[2] + lambda*cost_temp[2];
	candidate[3] = SAD_temp[3] + lambda*cost_temp[3];
	candidate[4] = SAD_temp[4] + lambda*cost_temp[4];
	candidate[5] = SAD_temp[5] + lambda*cost_temp[5]; 
	 
	min_array(candidate, 6, min_index);
		
	switch (min_index) {
      case 0:		
        break;
	  case 1:
		v_x[i*step+j] = v_x[i*step+(j-b_size)];
	    v_y[i*step+j] = v_y[i*step+(j-b_size)];
		break;
	  case 2:
		v_x[i*step+j] = v_x[(i+b_size)*step+j];
	    v_y[i*step+j] = v_y[(i+b_size)*step+j];
		break;
	  case 3:
		v_x[i*step+j] = v_x[i*step+(j+b_size)];
	    v_y[i*step+j] = v_y[i*step+(j+b_size)];
		break;
	  case 4:
		v_x[i*step+j] = v_x[(i+b_size)*step+(j-b_size)];
	    v_y[i*step+j] = v_y[(i+b_size)*step+(j-b_size)];
		break;
	  case 5:
	    v_x[i*step+j] = v_x[(i+b_size)*step+(j+b_size)];
	    v_y[i*step+j] = v_y[(i+b_size)*step+(j+b_size)];
		break; 
    }
  }

  //Left column
  j = start_pos;
  for (i = start_pos + b_size; i < height-(add_height - start_pos + b_size); i+=b_size)
  {    
    SAD_temp[0] = calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
	cost_temp[0] = cost_fn8left(i, j, i, j);	
	 
	SAD_temp[1] = calculate_SAD(i, j, v_x[(i-b_size)*step+j] + j, v_y[(i-b_size)*step+j] + i); 
	cost_temp[1] = cost_fn8left(i-b_size,j, i, j);
		
	SAD_temp[2] = calculate_SAD(i, j, v_x[(i+b_size)*step+j] + j, v_y[(i+b_size)*step+j] + i); 
	cost_temp[3] = cost_fn8left(i+b_size,j, i, j); 	 
	 
	SAD_temp[3] = calculate_SAD(i, j, v_x[i*step+(j+b_size)] + j, v_y[i*step+(j+b_size)] + i); 
	cost_temp[3] = cost_fn8left(i,j+b_size, i, j); 	 

	SAD_temp[4] = calculate_SAD(i, j, v_x[(i-b_size)*step+(j+b_size)] + j, v_y[(i-b_size)*step+(j+b_size)] + i); 
	cost_temp[4] = cost_fn8left(i-b_size,j+b_size, i, j); 	  
	  
	SAD_temp[5] = calculate_SAD(i, j, v_x[(i+b_size)*step+(j+b_size)] + j, v_y[(i+b_size)*step+(j+b_size)] + i); 
	cost_temp[5] = cost_fn8left(i+b_size,j+b_size, i, j);  
	 
	candidate[0] = SAD_temp[0] + lambda*cost_temp[0];
	candidate[1] = SAD_temp[1] + lambda*cost_temp[1];
	candidate[2] = SAD_temp[2] + lambda*cost_temp[2];
	candidate[3] = SAD_temp[3] + lambda*cost_temp[3];
	candidate[4] = SAD_temp[4] + lambda*cost_temp[4];
	candidate[5] = SAD_temp[5] + lambda*cost_temp[5]; 
	 
	min_array(candidate, 6, min_index);
		
	switch (min_index) {
      case 0:		
        break;
	  case 1:
		v_x[i*step+j] = v_x[(i-b_size)*step+j];
	    v_y[i*step+j] = v_y[(i-b_size)*step+j];
		break;
	  case 2:
		v_x[i*step+j] = v_x[(i+b_size)*step+j];
	    v_y[i*step+j] = v_y[(i+b_size)*step+j];
		break;
	  case 3:
		v_x[i*step+j] = v_x[i*step+(j+b_size)];
	    v_y[i*step+j] = v_y[i*step+(j+b_size)];
		break;
	  case 4:
	    v_x[i*step+j] = v_x[(i-b_size)*step+(j+b_size)];
	    v_y[i*step+j] = v_y[(i-b_size)*step+(j+b_size)];
		break;
	  case 5:
		v_x[i*step+j] = v_x[(i+b_size)*step+(j+b_size)];
	    v_y[i*step+j] = v_y[(i+b_size)*step+(j+b_size)];
		break;
    }
  }

  //Right column  
  if (((width - add_width) % b_size) == 0)
    j = width - (add_width - start_pos + b_size); 
  else
	j = width - (add_width - start_pos + ((width - add_width) % b_size));

  for (i = start_pos + b_size; i < height-(add_height - start_pos + b_size); i+=b_size)
  {   

    SAD_temp[0] = calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
	cost_temp[0] = cost_fn8right(i, j, i, j);	
	 
	SAD_temp[1] = calculate_SAD(i, j, v_x[(i-b_size)*step+j] + j, v_y[(i-b_size)*step+j] + i); 
	cost_temp[1] = cost_fn8right(i-b_size,j, i, j);
	  
	SAD_temp[2] = calculate_SAD(i, j, v_x[i*step+(j-b_size)] + j, v_y[i*step+(j-b_size)] + i); 
	cost_temp[2] = cost_fn8right(i,j-b_size, i, j); 	
	
	SAD_temp[3] = calculate_SAD(i, j, v_x[(i+b_size)*step+j] + j, v_y[(i+b_size)*step+j] + i); 
	cost_temp[3] = cost_fn8right(i+b_size,j, i, j); 	 
	 
	SAD_temp[4] = calculate_SAD(i, j, v_x[(i-b_size)*step+(j-b_size)] + j, v_y[(i-b_size)*step+(j-b_size)] + i); 
	cost_temp[4] = cost_fn8right(i-b_size,j-b_size, i, j); 
	 
	SAD_temp[5] = calculate_SAD(i, j, v_x[(i+b_size)*step+(j-b_size)] + j, v_y[(i+b_size)*step+(j-b_size)] + i); 
	cost_temp[5] = cost_fn8right(i+b_size,j-b_size, i, j); 	 	
    
	candidate[0] = SAD_temp[0] + lambda*cost_temp[0];
	candidate[1] = SAD_temp[1] + lambda*cost_temp[1];
	candidate[2] = SAD_temp[2] + lambda*cost_temp[2];
	candidate[3] = SAD_temp[3] + lambda*cost_temp[3];
	candidate[4] = SAD_temp[4] + lambda*cost_temp[4];
	candidate[5] = SAD_temp[5] + lambda*cost_temp[5]; 
	 
	min_array(candidate, 6, min_index);
		
	switch (min_index) {
      case 0:		
        break;
	  case 1:
		v_x[i*step+j] = v_x[(i-b_size)*step+j];
	    v_y[i*step+j] = v_y[(i-b_size)*step+j];
		break;	 
	  case 2:
		v_x[i*step+j] = v_x[i*step+(j-b_size)];
	    v_y[i*step+j] = v_y[i*step+(j-b_size)];
		break;
	  case 3:
	    v_x[i*step+j] = v_x[(i+b_size)*step+j];
	    v_y[i*step+j] = v_y[(i+b_size)*step+j];
		break;
	  case 4:
		v_x[i*step+j] = v_x[(i-b_size)*step+(j-b_size)];
	    v_y[i*step+j] = v_y[(i-b_size)*step+(j-b_size)];
		break;
	  case 5:
		v_x[i*step+j] = v_x[(i+b_size)*step+(j-b_size)];
	    v_y[i*step+j] = v_y[(i+b_size)*step+(j-b_size)];
		break;
    }  
  }

  //Bottom row
  if (((height - add_height) % b_size) == 0) 
    i = height - (add_height - start_pos + b_size);
  else
	i = height - (add_height - start_pos + ((height - add_height) % b_size));  

  for (j = start_pos + b_size; j < width-(add_width - start_pos + b_size); j+=b_size)
  { 
    SAD_temp[0] = calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
	cost_temp[0] = cost_fn8bottom(i, j, i, j);	
	 
	SAD_temp[1] = calculate_SAD(i, j, v_x[(i-b_size)*step+j] + j, v_y[(i-b_size)*step+j] + i); 
	cost_temp[1] = cost_fn8bottom(i-b_size,j, i, j);
	  
	SAD_temp[2] = calculate_SAD(i, j, v_x[i*step+(j-b_size)] + j, v_y[i*step+(j-b_size)] + i); 
	cost_temp[2] = cost_fn8bottom(i,j-b_size, i, j); 	
	
	SAD_temp[3] = calculate_SAD(i, j, v_x[i*step+(j+b_size)] + j, v_y[i*step+(j+b_size)] + i); 
	cost_temp[3] = cost_fn8bottom(i,j+b_size, i, j); 	 
	 
	SAD_temp[4] = calculate_SAD(i, j, v_x[(i-b_size)*step+(j-b_size)] + j, v_y[(i-b_size)*step+(j-b_size)] + i); 
	cost_temp[4] = cost_fn8bottom(i-b_size,j-b_size, i, j); 	 
	  
	SAD_temp[5] = calculate_SAD(i, j, v_x[(i-b_size)*step+(j+b_size)] + j, v_y[(i-b_size)*step+(j+b_size)] + i); 
	cost_temp[5] = cost_fn8bottom(i-b_size,j+b_size, i, j); 	  

	candidate[0] = SAD_temp[0] + lambda*cost_temp[0];
	candidate[1] = SAD_temp[1] + lambda*cost_temp[1];
	candidate[2] = SAD_temp[2] + lambda*cost_temp[2];
	candidate[3] = SAD_temp[3] + lambda*cost_temp[3];
	candidate[4] = SAD_temp[4] + lambda*cost_temp[4];
	candidate[5] = SAD_temp[5] + lambda*cost_temp[5]; 
	 
	min_array(candidate, 6, min_index);
		
	switch (min_index) {
      case 0:		
        break;
	  case 1:
		v_x[i*step+j] = v_x[(i-b_size)*step+j];
	    v_y[i*step+j] = v_y[(i-b_size)*step+j];
		break;	 
	  case 2:
		v_x[i*step+j] = v_x[i*step+(j-b_size)];
	    v_y[i*step+j] = v_y[i*step+(j-b_size)];
		break;
	  case 3:
	    v_x[i*step+j] = v_x[i*step+(j+b_size)];
	    v_y[i*step+j] = v_y[i*step+(j+b_size)];
		break;
	  case 4:
		v_x[i*step+j] = v_x[(i-b_size)*step+(j-b_size)];
	    v_y[i*step+j] = v_y[(i-b_size)*step+(j-b_size)];
		break;
	  case 5:
		v_x[i*step+j] = v_x[(i-b_size)*step+(j+b_size)];
	    v_y[i*step+j] = v_y[(i-b_size)*step+(j+b_size)];
		break;
    }
  }
  
  //Left corner
  i = start_pos;
  j = start_pos;

  SAD_temp[0] = calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
  cost_temp[0] = cost_fn8leftcorn(i, j, i, j);	
	  
  SAD_temp[1] = calculate_SAD(i, j, v_x[(i+b_size)*step+j] + j, v_y[(i+b_size)*step+j] + i); 
  cost_temp[1] = cost_fn8leftcorn(i+b_size,j, i, j); 	
	
  SAD_temp[2] = calculate_SAD(i, j, v_x[i*step+(j+b_size)] + j, v_y[i*step+(j+b_size)] + i); 
  cost_temp[2] = cost_fn8leftcorn(i,j+b_size, i, j); 

  SAD_temp[3] = calculate_SAD(i, j, v_x[(i+b_size)*step+(j+b_size)] + j, v_y[(i+b_size)*step+(j+b_size)] + i); 
  cost_temp[3] = cost_fn8leftcorn(i+b_size,j+b_size, i, j);

  min_array(candidate, 4, min_index);
		
  switch (min_index) {
    case 0:		
      break;
	case 1:
	  v_x[i*step+j] = v_x[(i+b_size)*step+j];
	  v_y[i*step+j] = v_y[(i+b_size)*step+j];
	  break;	 
	case 2:
	  v_x[i*step+j] = v_x[i*step+(j+b_size)];
	  v_y[i*step+j] = v_y[i*step+(j+b_size)];
	  break;
	case 3:
	  v_x[i*step+j] = v_x[(i+b_size)*step+(j+b_size)];
	  v_y[i*step+j] = v_y[(i+b_size)*step+(j+b_size)];
	  break;
  }

  //Right corner
  i = start_pos;  

  if (((width - add_width) % b_size) == 0)
    j = width - (add_width - start_pos + b_size); 
  else
	j = width - (add_width - start_pos + ((width - add_width) % b_size));

  SAD_temp[0] = calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
  cost_temp[0] = cost_fn8rightcorn(i, j, i, j);	
	  
  SAD_temp[1] = calculate_SAD(i, j, v_x[(i+b_size)*step+j] + j, v_y[(i+b_size)*step+j] + i); 
  cost_temp[1] = cost_fn8rightcorn(i+b_size,j, i, j); 	
	
  SAD_temp[2] = calculate_SAD(i, j, v_x[i*step+(j-b_size)] + j, v_y[i*step+(j-b_size)] + i); 
  cost_temp[2] = cost_fn8rightcorn(i,j-b_size, i, j); 

  SAD_temp[3] = calculate_SAD(i, j, v_x[(i+b_size)*step+(j-b_size)] + j, v_y[(i+b_size)*step+(j-b_size)] + i); 
  cost_temp[3] = cost_fn8rightcorn(i+b_size,j-b_size, i, j);

  min_array(candidate, 4, min_index);
		
  switch (min_index) {
    case 0:		
      break;
	case 1:
	  v_x[i*step+j] = v_x[(i+b_size)*step+j];
	  v_y[i*step+j] = v_y[(i+b_size)*step+j];
	  break;	 
	case 2:
	  v_x[i*step+j] = v_x[i*step+(j-b_size)];
	  v_y[i*step+j] = v_y[i*step+(j-b_size)];
	  break;
	case 3:
	  v_x[i*step+j] = v_x[(i+b_size)*step+(j-b_size)];
	  v_y[i*step+j] = v_y[(i+b_size)*step+(j-b_size)];
	  break;
  }

  //Bottom left corner
  if (((height - add_height) % b_size) == 0) 
    i = height - (add_height - start_pos + b_size);
  else
	i = height - (add_height - start_pos + ((height - add_height) % b_size)); 

  j = start_pos;

  SAD_temp[0] = calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
  cost_temp[0] = cost_fn8bottomleftcorn(i, j, i, j);	
	  
  SAD_temp[1] = calculate_SAD(i, j, v_x[(i-b_size)*step+j] + j, v_y[(i-b_size)*step+j] + i); 
  cost_temp[1] = cost_fn8bottomleftcorn(i-b_size,j, i, j); 	
	
  SAD_temp[2] = calculate_SAD(i, j, v_x[i*step+(j+b_size)] + j, v_y[i*step+(j+b_size)] + i); 
  cost_temp[2] = cost_fn8bottomleftcorn(i,j+b_size, i, j); 

  SAD_temp[3] = calculate_SAD(i, j, v_x[(i-b_size)*step+(j+b_size)] + j, v_y[(i-b_size)*step+(j+b_size)] + i); 
  cost_temp[3] = cost_fn8bottomleftcorn(i-b_size,j+b_size, i, j);

  min_array(candidate, 4, min_index);
		
  switch (min_index) {
    case 0:		
      break;
	case 1:
	  v_x[i*step+j] = v_x[(i-b_size)*step+j];
	  v_y[i*step+j] = v_y[(i-b_size)*step+j];
	  break;	 
	case 2:
	  v_x[i*step+j] = v_x[i*step+(j+b_size)];
	  v_y[i*step+j] = v_y[i*step+(j+b_size)];
	  break;
	case 3:
	  v_x[i*step+j] = v_x[(i-b_size)*step+(j+b_size)];
	  v_y[i*step+j] = v_y[(i-b_size)*step+(j+b_size)];
	  break;
  }
  
  //Bottom right corner
  if (((height - add_height) % b_size) == 0) 
    i = height - (add_height - start_pos + b_size);
  else
	i = height - (add_height - start_pos + ((height - add_height) % b_size));

  if (((width - add_width) % b_size) == 0)
    j = width - (add_width - start_pos + b_size); 
  else
	j = width - (add_width - start_pos + ((width - add_width) % b_size));

  SAD_temp[0] = calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i);
  cost_temp[0] = cost_fn8bottomrightcorn(i, j, i, j);	
	  
  SAD_temp[1] = calculate_SAD(i, j, v_x[(i-b_size)*step+j] + j, v_y[(i-b_size)*step+j] + i); 
  cost_temp[1] = cost_fn8bottomrightcorn(i-b_size,j, i, j); 	
	
  SAD_temp[2] = calculate_SAD(i, j, v_x[i*step+(j-b_size)] + j, v_y[i*step+(j-b_size)] + i); 
  cost_temp[2] = cost_fn8bottomrightcorn(i,j-b_size, i, j); 

  SAD_temp[3] = calculate_SAD(i, j, v_x[(i-b_size)*step+(j-b_size)] + j, v_y[(i-b_size)*step+(j-b_size)] + i); 
  cost_temp[3] = cost_fn8bottomrightcorn(i-b_size,j-b_size, i, j);

  min_array(candidate, 4, min_index);
		
  switch (min_index) {
    case 0:		
      break;
	case 1:
	  v_x[i*step+j] = v_x[(i-b_size)*step+j];
	  v_y[i*step+j] = v_y[(i-b_size)*step+j];
	  break;	 
	case 2:
	  v_x[i*step+j] = v_x[i*step+(j-b_size)];
	  v_y[i*step+j] = v_y[i*step+(j-b_size)];
	  break;
	case 3:
	  v_x[i*step+j] = v_x[(i-b_size)*step+(j-b_size)];
	  v_y[i*step+j] = v_y[(i-b_size)*step+(j-b_size)];
	  break;
  }	
}
void Pair_Motion::add_smoothness8_overlap(double lambda_value) 
{
  int min_index, i, j, k, l, n;
  double lambda;
  double candidate[9];  
  
  //The variables below are for removing duplicate MVs
  int duplicate_list[9] = {0};
  int hash[9] = {0};
  int MV_x[9];
  int MV_y[9];
    
  int bs_squared = b_size;
    
  lambda = lambda_value;
	  
  for (i = start_pos + b_size; i < height-(add_height - start_pos + b_size); i+=b_size) //goes through all vertical pixels
  {
	for (j = start_pos + b_size; j < width-(add_width - start_pos + b_size); j+=b_size) //goes through all horizontal pixels
    {	   

	  //Initalize duplicate list to zero -- no duplicates
	  for(k = 0; k < 9; k++)
	  {
       duplicate_list[k] = 0;  
	   hash[k] = 0;
	  }

	  MV_x[0] = v_x[i*step+j];
	  MV_y[0] = v_y[i*step+j];
	  MV_x[1] = v_x[(i-b_size)*step+j];
	  MV_y[1] = v_y[(i-b_size)*step+j];
	  MV_x[2] = v_x[i*step+(j-b_size)];
	  MV_y[2] = v_y[i*step+(j-b_size)];
	  MV_x[3] = v_x[(i+b_size)*step+j];
	  MV_y[3] = v_y[(i+b_size)*step+j];
	  MV_x[4] = v_x[i*step+(j+b_size)];
	  MV_y[4] = v_y[i*step+(j+b_size)];
	  MV_x[5] = v_x[(i-b_size)*step+(j-b_size)];
	  MV_y[5] = v_y[(i-b_size)*step+(j-b_size)];
	  MV_x[6] = v_x[(i+b_size)*step+(j-b_size)];
	  MV_y[6] = v_y[(i+b_size)*step+(j-b_size)];
	  MV_x[7] = v_x[(i-b_size)*step+(j+b_size)];
	  MV_y[7] = v_y[(i-b_size)*step+(j+b_size)];	  
	  MV_x[8] = v_x[(i+b_size)*step+(j+b_size)];
	  MV_y[8] = v_y[(i+b_size)*step+(j+b_size)];

	  //Check for duplicates, a duplicate is assigned a value of '1' in the duplicate list
      for(k = 1; k < 9; k++)
      {
        for(l = 0; l < 9; l++)
	    {
	      if ((l == k) || (duplicate_list[k] == 1))
	        continue;
          else if ((MV_x[k] == MV_x[l]) && (MV_y[k] == MV_y[l]))
	      {
	        duplicate_list[l] = 1;
	      } 	    
	    }
      }	  
	 
	  //Subtract out the overlap for the one already assigned prior to smoothness function.
	  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	    for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
		  Computed_Data.overlap[k*step+l] -= 1;

	  n = 0;
	  for(k = 0; k < 9; k++)
	  {
		if (duplicate_list[k] == 0)
		{
		  hash[n] = k; // this tell us which position of MV_x and MV_y the candidate is stored at.
	      candidate[n++] = (1 + calculate_SAD(i, j, MV_x[k] + j, MV_y[k] + i))*(1 + ((double)get_overlap_volume(MV_x[k] + j, MV_y[k] + i)/bs_squared)) + lambda*cost_fn(MV_y[k], MV_x[k], i, j);
		}
	  }

	  min_array(candidate, n, min_index); 	
	  v_x[i*step+j] = MV_x[hash[min_index]];
	  v_y[i*step+j] = MV_y[hash[min_index]];

	  //Next, add the overlap for the minimum found from min_array() above.
	  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	    for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
		  Computed_Data.overlap[k*step+l] += 1;    	
	}	
  }  
  
  //Top row
  i = start_pos;
  for (j = start_pos + b_size; j < width-(add_width - start_pos + b_size); j+=b_size)
  {    

    //Initalize duplicate list to zero -- no duplicates
	for(k = 0; k < 6; k++)
	{
      duplicate_list[k] = 0;  
	  hash[k] = 0;
	}

	MV_x[0] = v_x[i*step+j];
    MV_y[0] = v_y[i*step+j];
	MV_x[1] = v_x[i*step+(j-b_size)];
	MV_y[1] = v_y[i*step+(j-b_size)];
	MV_x[2] = v_x[(i+b_size)*step+j];
	MV_y[2] = v_y[(i+b_size)*step+j];
	MV_x[3] = v_x[i*step+(j+b_size)];
	MV_y[3] = v_y[i*step+(j+b_size)];
	MV_x[4] = v_x[(i+b_size)*step+(j-b_size)];
	MV_y[4] = v_y[(i+b_size)*step+(j-b_size)];		
	MV_x[5] = v_x[(i+b_size)*step+(j+b_size)];
	MV_y[5] = v_y[(i+b_size)*step+(j+b_size)];	  

	//Check for duplicates, a duplicate is assigned a value of '1' in the duplicate list
    for(k = 1; k < 6; k++)
    {
      for(l = 0; l < 6; l++)
	  {
	    if ((l == k) || (duplicate_list[k] == 1))
	      continue;
        else if ((MV_x[k] == MV_x[l]) && (MV_y[k] == MV_y[l]))
	    {
	      duplicate_list[l] = 1;
	    } 	    
	  }
    }	  
	
	//Subtract out the overlap for the one already assigned prior to smoothness function.
	for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	  for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	    Computed_Data.overlap[k*step+l] -= 1;

	n = 0;
	for(k = 0; k < 6; k++)
	{
	  if (duplicate_list[k] == 0)
	  {
		hash[n] = k; // this tell us which position of MV_x and MV_y the candidate is stored at.
	    candidate[n++] = (1 + calculate_SAD(i, j, MV_x[k] + j, MV_y[k] + i))*(1 + ((double)get_overlap_volume(MV_x[k] + j, MV_y[k] + i)/bs_squared)) + lambda*cost_fn_top_overlap(MV_y[k], MV_x[k], i, j);
	  }
	}
	
	min_array(candidate, n, min_index); 
	v_x[i*step+j] = MV_x[hash[min_index]];
	v_y[i*step+j] = MV_y[hash[min_index]];

	//Next, add the overlap for the minimum found from min_array() above.
	for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	  for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	    Computed_Data.overlap[k*step+l] += 1; 				  	
  }

  //Left column
  j = start_pos;
  for (i = start_pos + b_size; i < height-(add_height - start_pos + b_size); i+=b_size)
  {        
    //Initalize duplicate list to zero -- no duplicates
	for(k = 0; k < 6; k++)
	{
      duplicate_list[k] = 0;  
	  hash[k] = 0;
	}

	MV_x[0] = v_x[i*step+j];
	MV_y[0] = v_y[i*step+j];
	MV_x[1] = v_x[(i-b_size)*step+j];
	MV_y[1] = v_y[(i-b_size)*step+j];
	MV_x[2] = v_x[(i+b_size)*step+j];
	MV_y[2] = v_y[(i+b_size)*step+j];
	MV_x[3] = v_x[i*step+(j+b_size)];
	MV_y[3] = v_y[i*step+(j+b_size)];	 
	MV_x[4] = v_x[(i-b_size)*step+(j+b_size)];
	MV_y[4] = v_y[(i-b_size)*step+(j+b_size)];	
	MV_x[5] = v_x[(i+b_size)*step+(j+b_size)];
	MV_y[5] = v_y[(i+b_size)*step+(j+b_size)];	 
	
	//Check for duplicates, a duplicate is assigned a value of '1' in the duplicate list
    for(k = 1; k < 6; k++)
    {
      for(l = 0; l < 6; l++)
	  {
	    if ((l == k) || (duplicate_list[k] == 1))
	      continue;
        else if ((MV_x[k] == MV_x[l]) && (MV_y[k] == MV_y[l]))
	    {
	      duplicate_list[l] = 1;
	    } 	    
	  }
    }	  

	//Subtract out the overlap for the one already assigned prior to smoothness function.
	for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	  for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	    Computed_Data.overlap[k*step+l] -= 1;
	 
	n = 0;
	for(k = 0; k < 6; k++)
	{
	  if (duplicate_list[k] == 0)
	  {
		hash[n] = k; // this tell us which position of MV_x and MV_y the candidate is stored at.
	    candidate[n++] = (1 + calculate_SAD(i, j, MV_x[k] + j, MV_y[k] + i))*(1 + ((double)get_overlap_volume(MV_x[k] + j, MV_y[k] + i)/bs_squared)) + lambda*cost_fn_left_overlap(MV_y[k], MV_x[k], i, j);
	  }
	}
	
	min_array(candidate, n, min_index); 
	v_x[i*step+j] = MV_x[hash[min_index]];
	v_y[i*step+j] = MV_y[hash[min_index]];

	//Next, add the overlap for the minimum found from min_array() above.
	for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	  for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	    Computed_Data.overlap[k*step+l] += 1;     
  }

  //Right column
  if (((width - add_width) % b_size) == 0)
    j = width - (add_width - start_pos + b_size); 
  else
	j = width - (add_width - start_pos + ((width - add_width) % b_size));

  for (i = start_pos + b_size; i < height-(add_height - start_pos + b_size); i+=b_size)
  {       
    //Initalize duplicate list to zero -- no duplicates
	for(k = 0; k < 6; k++)
	{
      duplicate_list[k] = 0;  
	  hash[k] = 0;
	}

	MV_x[0] = v_x[i*step+j];
	MV_y[0] = v_y[i*step+j];
	MV_x[1] = v_x[(i-b_size)*step+j];
	MV_y[1] = v_y[(i-b_size)*step+j];
	MV_x[2] = v_x[i*step+(j-b_size)];
	MV_y[2] = v_y[i*step+(j-b_size)];
	MV_x[3] = v_x[(i+b_size)*step+j];
	MV_y[3] = v_y[(i+b_size)*step+j];
	MV_x[4] = v_x[(i-b_size)*step+(j-b_size)];
	MV_y[4] = v_y[(i-b_size)*step+(j-b_size)];		 
	MV_x[5] = v_x[(i+b_size)*step+(j-b_size)];
	MV_y[5] = v_y[(i+b_size)*step+(j-b_size)];	

	//Check for duplicates, a duplicate is assigned a value of '1' in the duplicate list
    for(k = 1; k < 6; k++)
    {
      for(l = 0; l < 6; l++)
	  {
	    if ((l == k) || (duplicate_list[k] == 1))
	      continue;
        else if ((MV_x[k] == MV_x[l]) && (MV_y[k] == MV_y[l]))
	    {
	      duplicate_list[l] = 1;
	    } 	    
	  }
    }	  
	 
	//Subtract out the overlap for the one already assigned prior to smoothness function.
	for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	  for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	    Computed_Data.overlap[k*step+l] -= 1;

	n = 0;
	for(k = 0; k < 6; k++)
	{
	  if (duplicate_list[k] == 0)
	  {
		hash[n] = k; // this tell us which position of MV_x and MV_y the candidate is stored at.
	    candidate[n++] = (1 + calculate_SAD(i, j, MV_x[k] + j, MV_y[k] + i))*(1 + ((double)get_overlap_volume(MV_x[k] + j, MV_y[k] + i)/bs_squared)) + lambda*cost_fn_right_overlap(MV_y[k], MV_x[k], i, j);
	  }
	}
	
	min_array(candidate, n, min_index); 
	v_x[i*step+j] = MV_x[hash[min_index]];
	v_y[i*step+j] = MV_y[hash[min_index]];

	//Next, add the overlap for the minimum found from min_array() above.
	for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	  for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	    Computed_Data.overlap[k*step+l] += 1;      
  }

  //Bottom row   
  if (((height - add_height) % b_size) == 0) 
    i = height - (add_height - start_pos + b_size);
  else
	i = height - (add_height - start_pos + ((height - add_height) % b_size));

  for (j = start_pos + b_size; j < width-(add_width - start_pos + b_size); j+=b_size)
  { 
    //Initalize duplicate list to zero -- no duplicates
	for(k = 0; k < 6; k++)
	{
      duplicate_list[k] = 0;  
	  hash[k] = 0;
	}

	MV_x[0] = v_x[i*step+j];
	MV_y[0] = v_y[i*step+j];
	MV_x[1] = v_x[(i-b_size)*step+j];
	MV_y[1] = v_y[(i-b_size)*step+j];
	MV_x[2] = v_x[i*step+(j-b_size)];
	MV_y[2] = v_y[i*step+(j-b_size)];
	MV_x[3] = v_x[i*step+(j+b_size)];
	MV_y[3] = v_y[i*step+(j+b_size)];
	MV_x[4] = v_x[(i-b_size)*step+(j-b_size)];
	MV_y[4] = v_y[(i-b_size)*step+(j-b_size)];	
	MV_x[5] = v_x[(i-b_size)*step+(j+b_size)];
	MV_y[5] = v_y[(i-b_size)*step+(j+b_size)];	

	//Check for duplicates, a duplicate is assigned a value of '1' in the duplicate list
    for(k = 1; k < 6; k++)
    {
      for(l = 0; l < 6; l++)
	  {
	    if ((l == k) || (duplicate_list[k] == 1))
	      continue;
        else if ((MV_x[k] == MV_x[l]) && (MV_y[k] == MV_y[l]))
	    {
	      duplicate_list[l] = 1;
	    } 	    
	  }
    }	  
	
	//Subtract out the overlap for the one already assigned prior to smoothness.
	for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	  for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	    Computed_Data.overlap[k*step+l] -= 1;

	n = 0;
	for(k = 0; k < 6; k++)
	{
	  if (duplicate_list[k] == 0)
	  {
		hash[n] = k; // this tell us which position of MV_x and MV_y the candidate is stored at.
	    candidate[n++] = (1 + calculate_SAD(i, j, MV_x[k] + j, MV_y[k] + i))*(1 + ((double)get_overlap_volume(MV_x[k] + j, MV_y[k] + i)/bs_squared)) + lambda*cost_fn_bottom_overlap(MV_y[k], MV_x[k], i, j);
	  }
	}
	  	 	  
	min_array(candidate, n, min_index); 
	v_x[i*step+j] = MV_x[hash[min_index]];
	v_y[i*step+j] = MV_y[hash[min_index]];   

	//Next, add the overlap for the minimum found from min_array() above.
	for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	  for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	    Computed_Data.overlap[k*step+l] += 1;
	
  }  

  //CORNERS!

  //Left corner
  i = start_pos;
  j = start_pos;

  //Initalize duplicate list to zero -- no duplicates
  for(k = 0; k < 4; k++)
  {
    duplicate_list[k] = 0;  
	hash[k] = 0;
  }

  MV_x[0] = v_x[i*step+j];
  MV_y[0] = v_y[i*step+j];
  MV_x[1] = v_x[(i+b_size)*step+j];
  MV_y[1] = v_y[(i+b_size)*step+j];	
  MV_x[2] = v_x[i*step+(j+b_size)];
  MV_y[2] = v_y[i*step+(j+b_size)];	     
  MV_x[3] = v_x[(i+b_size)*step+(j+b_size)];
  MV_y[3] = v_y[(i+b_size)*step+(j+b_size)];	 

  //Check for duplicates, a duplicate is assigned a value of '1' in the duplicate list
  for(k = 1; k < 4; k++)
  {
    for(l = 0; l < 4; l++)
	{
	  if ((l == k) || (duplicate_list[k] == 1))
	    continue;
      else if ((MV_x[k] == MV_x[l]) && (MV_y[k] == MV_y[l]))
	  {
	    duplicate_list[l] = 1;
	  } 	    
	}
  }	  
	 
  //Subtract out the overlap for the one already assigned prior to smoothness.
  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	  Computed_Data.overlap[k*step+l] -= 1;

  n = 0;
  for(k = 0; k < 4; k++)
  {
	if (duplicate_list[k] == 0)
	{
	  hash[n] = k; // this tell us which position of MV_x and MV_y the candidate is stored at.
	  candidate[n++] = (1 + calculate_SAD(i, j, MV_x[k] + j, MV_y[k] + i))*(1 + ((double)get_overlap_volume(MV_x[k] + j, MV_y[k] + i)/bs_squared)) + lambda*cost_fn_leftcorn_overlap(MV_y[k], MV_x[k], i, j);
	}
  }
	  	 	  
  min_array(candidate, n, min_index); 
  v_x[i*step+j] = MV_x[hash[min_index]];
  v_y[i*step+j] = MV_y[hash[min_index]];   

  //Next, add the overlap for the minimum found from min_array() above.
  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	  Computed_Data.overlap[k*step+l] += 1; 
  
  
  //Right corner
  i = start_pos;

  if (((width - add_width) % b_size) == 0)
    j = width - (add_width - start_pos + b_size); 
  else
	j = width - (add_width - start_pos + ((width - add_width) % b_size));
  
  //Initalize duplicate list to zero -- no duplicates
  for(k = 0; k < 4; k++)
  {
    duplicate_list[k] = 0;  
	hash[k] = 0;
  }

  MV_x[0] = v_x[i*step+j];
  MV_y[0] = v_y[i*step+j]; 
  MV_x[1] = v_x[(i+b_size)*step+j];
  MV_y[1] = v_y[(i+b_size)*step+j]; 
  MV_x[2] = v_x[i*step+(j-b_size)];
  MV_y[2] = v_y[i*step+(j-b_size)];    
  MV_x[3] = v_x[(i+b_size)*step+(j-b_size)];
  MV_y[3] = v_y[(i+b_size)*step+(j-b_size)]; 

  //Check for duplicates, a duplicate is assigned a value of '1' in the duplicate list
  for(k = 1; k < 4; k++)
  {
    for(l = 0; l < 4; l++)
	{
	  if ((l == k) || (duplicate_list[k] == 1))
	    continue;
      else if ((MV_x[k] == MV_x[l]) && (MV_y[k] == MV_y[l]))
	  {
	    duplicate_list[l] = 1;
	  } 	    
	}
  }	  
	 
  //Subtract out the overlap for the one already assigned prior to smoothness.
  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	  Computed_Data.overlap[k*step+l] -= 1;

  n = 0;
  for(k = 0; k < 4; k++)
  {
	if (duplicate_list[k] == 0)
	{
	  hash[n] = k; // this tell us which position of MV_x and MV_y the candidate is stored at.
	  candidate[n++] = (1 + calculate_SAD(i, j, MV_x[k] + j, MV_y[k] + i))*(1 + ((double)get_overlap_volume(MV_x[k] + j, MV_y[k] + i)/bs_squared)) + lambda*cost_fn_rightcorn_overlap(MV_y[k], MV_x[k], i, j);
	}
  }
	  	 	  
  min_array(candidate, n, min_index); 
  v_x[i*step+j] = MV_x[hash[min_index]];
  v_y[i*step+j] = MV_y[hash[min_index]];   

  //Next, add the overlap for the minimum found from min_array() above.
  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	  Computed_Data.overlap[k*step+l] += 1;


  //Bottom left corner
  if (((height - add_height) % b_size) == 0) 
    i = height - (add_height - start_pos + b_size);
  else
	i = height - (add_height - start_pos + ((height - add_height) % b_size));

  j = start_pos;

  //Initalize duplicate list to zero -- no duplicates
  for(k = 0; k < 4; k++)
  {
    duplicate_list[k] = 0;  
	hash[k] = 0;
  }

  MV_x[0] = v_x[i*step+j];
  MV_y[0] = v_y[i*step+j];
  MV_x[1] = v_x[(i-b_size)*step+j];
  MV_y[1] = v_y[(i-b_size)*step+j]; 
  MV_x[2] = v_x[i*step+(j+b_size)];
  MV_y[2] = v_y[i*step+(j+b_size)];    
  MV_x[3] = v_x[(i-b_size)*step+(j+b_size)];
  MV_y[3] = v_y[(i-b_size)*step+(j+b_size)];  

  //Check for duplicates, a duplicate is assigned a value of '1' in the duplicate list
  for(k = 1; k < 4; k++)
  {
    for(l = 0; l < 4; l++)
	{
	  if ((l == k) || (duplicate_list[k] == 1))
	    continue;
      else if ((MV_x[k] == MV_x[l]) && (MV_y[k] == MV_y[l]))
	  {
	    duplicate_list[l] = 1;
	  } 	    
	}
  }	  

  //Subtract out the overlap for the one already assigned prior to smoothness.
  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	  Computed_Data.overlap[k*step+l] -= 1;

  n = 0;
  for(k = 0; k < 4; k++)
  {
	if (duplicate_list[k] == 0)
	{
	  hash[n] = k; // this tell us which position of MV_x and MV_y the candidate is stored at.
	  candidate[n++] = (1 + calculate_SAD(i, j, MV_x[k] + j, MV_y[k] + i))*(1 + ((double)get_overlap_volume(MV_x[k] + j, MV_y[k] + i)/bs_squared)) + lambda*cost_fn_bottomleftcorn_overlap(MV_y[k], MV_x[k], i, j);
	}
  }
	  	 	  
  min_array(candidate, n, min_index); 
  v_x[i*step+j] = MV_x[hash[min_index]];
  v_y[i*step+j] = MV_y[hash[min_index]];   

  //Next, add the overlap for the minimum found from min_array() above.
  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	  Computed_Data.overlap[k*step+l] += 1;
   

  //Bottom right corner
  if (((height - add_height) % b_size) == 0) 
    i = height - (add_height - start_pos + b_size);
  else
	i = height - (add_height - start_pos + ((height - add_height) % b_size));

  if (((width - add_width) % b_size) == 0)
    j = width - (add_width - start_pos + b_size); 
  else
	j = width - (add_width - start_pos + ((width - add_width) % b_size));

  //Initalize duplicate list to zero -- no duplicates
  for(k = 0; k < 4; k++)
  {
    duplicate_list[k] = 0;  
	hash[k] = 0;
  }

  MV_x[0] = v_x[i*step+j];
  MV_y[0] = v_y[i*step+j];  
  MV_x[1] = v_x[(i-b_size)*step+j];
  MV_y[1] = v_y[(i-b_size)*step+j];
  MV_x[2] = v_x[i*step+(j-b_size)];
  MV_y[2] = v_y[i*step+(j-b_size)]; 
  MV_x[3] = v_x[(i-b_size)*step+(j-b_size)];
  MV_y[3] = v_y[(i-b_size)*step+(j-b_size)];

  //Check for duplicates, a duplicate is assigned a value of '1' in the duplicate list
  for(k = 1; k < 4; k++)
  {
    for(l = 0; l < 4; l++)
	{
	  if ((l == k) || (duplicate_list[k] == 1))
	    continue;
      else if ((MV_x[k] == MV_x[l]) && (MV_y[k] == MV_y[l]))
	  {
	    duplicate_list[l] = 1;
	  } 	    
	}
  }	  
	 
  //Subtract out the overlap for the one already assigned prior to smoothness.
  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	  Computed_Data.overlap[k*step+l] -= 1;

  n = 0;
  for(k = 0; k < 4; k++)
  {
	if (duplicate_list[k] == 0)
	{
	  hash[n] = k; // this tell us which position of MV_x and MV_y the candidate is stored at.
	  candidate[n++] = (1 + calculate_SAD(i, j, MV_x[k] + j, MV_y[k] + i))*(1 + ((double)get_overlap_volume(MV_x[k] + j, MV_y[k] + i)/bs_squared)) + lambda*cost_fn_bottomrightcorn_overlap(MV_y[k], MV_x[k], i, j);
	}
  }
	  	 	  
  min_array(candidate, n, min_index); 
  v_x[i*step+j] = MV_x[hash[min_index]];
  v_y[i*step+j] = MV_y[hash[min_index]];   

  //Next, add the overlap for the minimum found from min_array() above.
  for(k = v_y[i*step+j] + i; k < v_y[i*step+j] + i + b_size; k++)
	for(l = v_x[i*step+j] + j; l < v_x[i*step+j] + j + b_size; l++)
	  Computed_Data.overlap[k*step+l] += 1;

}
double Pair_Motion::calc_derivX(int k, int l)
{
  return (data_a[k*step+(l+1)] - data_a[k*step+l]);
}
double Pair_Motion::calc_derivY(int k, int l)
{
  return (data_a[(k+1)*step+l] - data_a[k*step+l]);
}
double Pair_Motion::cost_fn(int v_yval, int v_xval, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_5, sum_6, sum_7, sum_8;
  int entry_x = v_xval;
  int entry_y = v_yval;

  //assumes 8-connected neighborhood  

  //Function 3 in paper
  sum_1 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn_dynamic_smoothness(int v_yval, int v_xval, int orig_y, int orig_x)
{

  double sum_1, sum_2, sum_3, sum_4;
  int entry_x = v_xval;
  int entry_y = v_yval;

  //assumes top left, left, top, top right blocks

  //Function 3 in paper
  sum_1 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);  
  sum_3 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]); 
  sum_4 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]); 

  return (sum_1 + sum_2 + sum_3 + sum_4);     

}
double Pair_Motion::cost_fn_top_overlap(int v_yval, int v_xval, int orig_y, int orig_x)
{
  double sum_2, sum_3, sum_4, sum_5, sum_6;

  int entry_x = v_xval;
  int entry_y = v_yval;

  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);  
  sum_5 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]); 
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  return (sum_2 + sum_3 + sum_4 + sum_5 + sum_6);   

}
double Pair_Motion::cost_fn_left_overlap(int v_yval, int v_xval, int orig_y, int orig_x)
{
  double sum_2, sum_3, sum_4, sum_5, sum_6;
 
  int entry_x = v_xval;
  int entry_y = v_yval;
  
  sum_2 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);  
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);  
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  return (sum_2 + sum_3 + sum_4 + sum_5 + sum_6);    
}
double Pair_Motion::cost_fn_right_overlap(int v_yval, int v_xval, int orig_y, int orig_x)
{
  double sum_2, sum_3, sum_4, sum_5, sum_6;
 
  int entry_x = v_xval;
  int entry_y = v_yval;
  
  sum_2 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_3 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_4 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);  
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  
  return (sum_2 + sum_3 + sum_4 + sum_5 + sum_6);   
}
double Pair_Motion::cost_fn_bottom_overlap(int v_yval, int v_xval, int orig_y, int orig_x)
{
  double sum_2, sum_3, sum_4, sum_5, sum_6;
 
  int entry_x = v_xval;
  int entry_y = v_yval;
  
  sum_2 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_3 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);  
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);  
  sum_6 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]); 

  return (sum_2 + sum_3 + sum_4 + sum_5 + sum_6);
}
double Pair_Motion::cost_fn_leftcorn_overlap(int v_yval, int v_xval, int orig_y, int orig_x)
{
  double sum_2, sum_3, sum_4;

  int entry_x = v_xval;
  int entry_y = v_yval;    

  sum_2 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_3 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);  
  sum_4 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  return (sum_2 + sum_3 + sum_4); 
}
double Pair_Motion::cost_fn_rightcorn_overlap(int v_yval, int v_xval, int orig_y, int orig_x)
{
  double sum_2, sum_3, sum_4;

  int entry_x = v_xval;
  int entry_y = v_yval;    
  
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);  
  sum_4 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);  

  return (sum_2 + sum_3 + sum_4);   
}
double Pair_Motion::cost_fn_bottomleftcorn_overlap(int v_yval, int v_xval, int orig_y, int orig_x)
{
  double sum_2, sum_3, sum_4;

  int entry_x = v_xval;
  int entry_y = v_yval;   
  
  sum_2 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);  
  sum_3 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);  
  sum_4 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);  

  return (sum_2 + sum_3 + sum_4);	
}
double Pair_Motion::cost_fn_bottomrightcorn_overlap(int v_yval, int v_xval, int orig_y, int orig_x)
{
  double sum_2, sum_3, sum_4;

  int entry_x = v_xval;
  int entry_y = v_yval;   
 
  sum_2 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_3 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]); 
  sum_4 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);  

  return (sum_2 + sum_3 + sum_4);     
}
double Pair_Motion::cost_fn8(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_5, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood
  
  //Function 1 in paper
  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/

  //Function 2 in paper
  /*int sum_array[8]; 
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[2] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[3] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[4] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[5] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[6] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[7] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 8));*/

  //Function 3 in paper
  sum_1 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8top(int pos_y, int pos_x, int orig_y, int orig_x)
{
  double sum_2, sum_3, sum_4, sum_5, sum_6;
  //int boundary_min = 10;
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //sum_1 = abs(entry_x - v_x[orig_y*step+orig_x]) + abs(entry_y - v_y[orig_y*step+orig_x]);
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);  
  sum_5 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]); 
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  return (sum_2 + sum_3 + sum_4 + sum_5 + sum_6);   
}
double Pair_Motion::cost_fn8bottom(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_2, sum_3, sum_4, sum_5, sum_6;
 
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //sum_1 = abs(entry_x - v_x[orig_y*step+orig_x]) + abs(entry_y - v_y[orig_y*step+orig_x]);
  sum_2 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_3 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);  
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);  
  sum_6 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]); 

  return (sum_2 + sum_3 + sum_4 + sum_5 + sum_6);     
}
double Pair_Motion::cost_fn8left(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_2, sum_3, sum_4, sum_5, sum_6;
 
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];
  
  //sum_1 = abs(entry_x - v_x[orig_y*step+orig_x]) + abs(entry_y - v_y[orig_y*step+orig_x]);
  sum_2 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);  
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);  
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  return (sum_2 + sum_3 + sum_4 + sum_5 + sum_6);     
}
double Pair_Motion::cost_fn8right(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_2, sum_3, sum_4, sum_5, sum_6;
 
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //sum_1 = abs(entry_x - v_x[orig_y*step+orig_x]) + abs(entry_y - v_y[orig_y*step+orig_x]);  
  sum_2 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_3 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_4 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);  
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  
  return (sum_2 + sum_3 + sum_4 + sum_5 + sum_6);     
}
double Pair_Motion::cost_fn8leftcorn(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_2, sum_3, sum_4;

  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];    

  //sum_1 = abs(entry_x - v_x[orig_y*step+orig_x]) + abs(entry_y - v_y[orig_y*step+orig_x]);  
  sum_2 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_3 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);  
  sum_4 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  return (sum_2 + sum_3 + sum_4);     
}
double Pair_Motion::cost_fn8rightcorn(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_2, sum_3, sum_4;

  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];    

  //sum_1 = abs(entry_x - v_x[orig_y*step+orig_x]) + abs(entry_y - v_y[orig_y*step+orig_x]);  
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);  
  sum_4 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);  

  return (sum_2 + sum_3 + sum_4);     
}
double Pair_Motion::cost_fn8bottomleftcorn(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_2, sum_3, sum_4;

  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];    

  //sum_1 = abs(entry_x - v_x[orig_y*step+orig_x]) + abs(entry_y - v_y[orig_y*step+orig_x]);
  sum_2 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);  
  sum_3 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);  
  sum_4 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);  

  return (sum_2 + sum_3 + sum_4);     
}
double Pair_Motion::cost_fn8bottomrightcorn(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_2, sum_3, sum_4;

  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];    

  //sum_1 = abs(entry_x - v_x[orig_y*step+orig_x]) + abs(entry_y - v_y[orig_y*step+orig_x]);
  sum_2 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_3 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]); 
  sum_4 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);  

  return (sum_2 + sum_3 + sum_4);     
}
double Pair_Motion::cost_fn8_weighted(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_5, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  double entry_x = v_x[pos_y*step+pos_x];
  double entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood
  
  //Function 1 in paper
  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/

  //Function 2 in paper
  /*int sum_array[8]; 
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[2] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[3] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[4] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[5] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[6] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[7] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 8));*/

  //Function 3 in paper
  sum_1 = (abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y-b_size)*step+orig_x];
  sum_2 = (abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x-b_size)];
  sum_3 = (abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y+b_size)*step+orig_x];
  sum_4 = (abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x+b_size)];
  sum_5 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x-b_size)];
  sum_6 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x-b_size)];
  sum_7 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x+b_size)];
  sum_8 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x+b_size)];

  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8a(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_2, sum_3, sum_4, sum_5, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood
   
  //sum_1 = 0; // min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[1] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[2] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[3] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[4] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_2 + sum_3 + sum_4 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8a_weighted(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_2, sum_3, sum_4, sum_5, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  double entry_x = v_x[pos_y*step+pos_x];
  double entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood
   
  //sum_1 = 0; // min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_2 = (abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x-b_size)];
  sum_3 = (abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y+b_size)*step+orig_x];
  sum_4 = (abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x+b_size)];
  sum_5 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x-b_size)];
  sum_6 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x-b_size)];
  sum_7 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x+b_size)];
  sum_8 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x+b_size)];

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[1] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[2] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[3] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[4] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_2 + sum_3 + sum_4 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8b(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_3, sum_4, sum_5, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood
   
  sum_1 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  //sum_2 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  //sum_2 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[1] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[2] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[3] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[4] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  //sum_2 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_3 + sum_4 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8b_weighted(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_3, sum_4, sum_5, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  double entry_x = v_x[pos_y*step+pos_x]; 
  double entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood
   
  sum_1 = (abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y-b_size)*step+orig_x];
  //sum_2 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_3 = (abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y+b_size)*step+orig_x];
  sum_4 = (abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x+b_size)];
  sum_5 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x-b_size)];
  sum_6 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x-b_size)];
  sum_7 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x+b_size)];
  sum_8 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x+b_size)];

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  //sum_2 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[1] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[2] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[3] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[4] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  //sum_2 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_3 + sum_4 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8c(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_4, sum_5, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood
   
  sum_1 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  //sum_3 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  //sum_3 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[2] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[3] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[4] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  //sum_3 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_2 + sum_4 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8c_weighted(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_4, sum_5, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  double entry_x = v_x[pos_y*step+pos_x];
  double entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood
   
  sum_1 = (abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y-b_size)*step+orig_x];
  sum_2 = (abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x-b_size)];
  //sum_3 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_4 = (abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x+b_size)];
  sum_5 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x-b_size)];
  sum_6 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x-b_size)];
  sum_7 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x+b_size)];
  sum_8 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x+b_size)];

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  //sum_3 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[2] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[3] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[4] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  //sum_3 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_2 + sum_4 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8d(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_5, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood   
  sum_1 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  //sum_4 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[2] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  //sum_4 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[3] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[4] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  //sum_4 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_2 + sum_3 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8d_weighted(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_5, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  double entry_x = v_x[pos_y*step+pos_x]; 
  double entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood   
  sum_1 = (abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y-b_size)*step+orig_x];
  sum_2 = (abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x-b_size)];
  sum_3 = (abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y+b_size)*step+orig_x];
  //sum_4 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_5 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x-b_size)];
  sum_6 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x-b_size)];
  sum_7 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x+b_size)];
  sum_8 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x+b_size)];

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[2] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  //sum_4 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[3] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[4] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  //sum_4 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_2 + sum_3 + sum_5 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8e(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood   
  sum_1 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  //sum_5 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[2] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[3] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  //sum_5 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[4] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  //sum_5 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8e_weighted(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_6, sum_7, sum_8;
  //int boundary_min = 10;
  double entry_x = v_x[pos_y*step+pos_x]; 
  double entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood   
  sum_1 = (abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y-b_size)*step+orig_x];
  sum_2 = (abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x-b_size)];
  sum_3 = (abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y+b_size)*step+orig_x];
  sum_4 = (abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x+b_size)];
  //sum_5 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_6 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x-b_size)];
  sum_7 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x+b_size)];
  sum_8 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x+b_size)];

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[2] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[3] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  //sum_5 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[4] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  //sum_5 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_6 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8f(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_5, sum_7, sum_8;
  //int boundary_min = 10;
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood  
  sum_1 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  //sum_6 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_7 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[2] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[3] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[4] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  //sum_6 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  //sum_6 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]); */
  
  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_5 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8f_weighted(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_5, sum_7, sum_8;
  //int boundary_min = 10;
  double entry_x = v_x[pos_y*step+pos_x]; 
  double entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood  
  sum_1 = (abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y-b_size)*step+orig_x];
  sum_2 = (abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x-b_size)];
  sum_3 = (abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y+b_size)*step+orig_x];
  sum_4 = (abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x+b_size)];
  sum_5 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x-b_size)];
  //sum_6 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_7 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x+b_size)];
  sum_8 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x+b_size)];

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[2] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[3] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[4] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  //sum_6 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[5] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  //sum_6 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]); */
  
  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_5 + sum_7 + sum_8);     
}
double Pair_Motion::cost_fn8g(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_5, sum_6, sum_8;
  //int boundary_min = 10;
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood   
  sum_1 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  //sum_7 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_8 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]);

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[2] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[3] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[4] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[5] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  //sum_7 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  //sum_7 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_5 + sum_6 + sum_8);     
}
double Pair_Motion::cost_fn8g_weighted(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_5, sum_6, sum_8;
  //int boundary_min = 10;
  double entry_x = v_x[pos_y*step+pos_x]; 
  double entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood   
  sum_1 = (abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y-b_size)*step+orig_x];
  sum_2 = (abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x-b_size)];
  sum_3 = (abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y+b_size)*step+orig_x];
  sum_4 = (abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x+b_size)];
  sum_5 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x-b_size)];
  sum_6 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x-b_size)];
  //sum_7 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_8 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x+b_size)];

  //function 2 in paper
  /*int sum_array[7];
  sum_array[0] = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_array[1] = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_array[2] = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_array[3] = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_array[4] = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_array[5] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  //sum_7 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_array[6] = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]))*spatial_test[7];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  //sum_7 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));
  sum_8 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_5 + sum_6 + sum_8);     
}
double Pair_Motion::cost_fn8h(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_5, sum_6, sum_7;
  //int boundary_min = 10;
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood  
  sum_1 = abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_2 = abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_4 = abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]);
  //sum_8 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));

  //function 2 in paper
  /*int sum_array[7];
  sum_1 = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_2 = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_3 = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_4 = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_5 = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_6 = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_7 = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_5 + sum_6 + sum_7);     
}
double Pair_Motion::cost_fn8h_weighted(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  double sum_1, sum_2, sum_3, sum_4, sum_5, sum_6, sum_7;
  //int boundary_min = 10;
  double entry_x = v_x[pos_y*step+pos_x]; 
  double entry_y = v_y[pos_y*step+pos_x];

  //assumes 8-connected neighborhood  
  sum_1 = (abs(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y-b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y-b_size)*step+orig_x];
  sum_2 = (abs(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x-b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x-b_size)];
  sum_3 = (abs(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + abs(entry_y - v_y[(orig_y+b_size)*step+orig_x]))*Computed_Data.reliability[(orig_y+b_size)*step+orig_x];
  sum_4 = (abs(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + abs(entry_y - v_y[orig_y*step+(orig_x+b_size)]))*Computed_Data.reliability[orig_y*step+(orig_x+b_size)];
  sum_5 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x-b_size)];
  sum_6 = (abs(entry_x - v_x[(orig_y+b_size)*step+(orig_x-b_size)]) + abs(entry_y - v_y[(orig_y+b_size)*step+(orig_x-b_size)]))*Computed_Data.reliability[(orig_y+b_size)*step+(orig_x-b_size)];
  sum_7 = (abs(entry_x - v_x[(orig_y-b_size)*step+(orig_x+b_size)]) + abs(entry_y - v_y[(orig_y-b_size)*step+(orig_x+b_size)]))*Computed_Data.reliability[(orig_y-b_size)*step+(orig_x+b_size)];
  //sum_8 = 0; //min(boundary_min, abs(entry_x - v_x[orig_x*step+orig_y]) + abs(entry_y - v_y[orig_x*step+orig_y]));

  //function 2 in paper
  /*int sum_array[7];
  sum_1 = (abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]))*spatial_test[0];
  sum_2 = (abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]))*spatial_test[1];
  sum_3 = (abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]))*spatial_test[2];
  sum_4 = (abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]))*spatial_test[3];
  sum_5 = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]))*spatial_test[4];
  sum_6 = (abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]))*spatial_test[5];
  sum_7 = (abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]))*spatial_test[6];
  return (min_array_val(sum_array, 7));*/

  /*sum_1 = abs(entry_x - v_x[(orig_x-b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x-b_size)*step+orig_y]);
  sum_2 = abs(entry_x - v_x[orig_x*step+(orig_y-b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y-b_size)]);
  sum_3 = abs(entry_x - v_x[(orig_x+b_size)*step+orig_y]) + abs(entry_y - v_y[(orig_x+b_size)*step+orig_y]);
  sum_4 = abs(entry_x - v_x[orig_x*step+(orig_y+b_size)]) + abs(entry_y - v_y[orig_x*step+(orig_y+b_size)]);
  sum_5 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y-b_size)]);
  sum_6 = abs(entry_x - v_x[(orig_x+b_size)*step+(orig_y-b_size)]) + abs(entry_y - v_y[(orig_x+b_size)*step+(orig_y-b_size)]);
  sum_7 = abs(entry_x - v_x[(orig_x-b_size)*step+(orig_y+b_size)]) + abs(entry_y - v_y[(orig_x-b_size)*step+(orig_y+b_size)]);*/
  
  return (sum_1 + sum_2 + sum_3 + sum_4 + sum_5 + sum_6 + sum_7);     
}
double Pair_Motion::vector_cost1(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block

  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];
    
  //Function 2 in paper
  double sum_array[4]; 
  sum_array[0] = (entry_x - v_x[(orig_y-b_size)*step+orig_x])*(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + (entry_y - v_y[(orig_y-b_size)*step+orig_x])*(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_array[1] = (entry_x - v_x[orig_y*step+(orig_x-b_size)])*(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + (entry_y - v_y[orig_y*step+(orig_x-b_size)])*(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_array[2] = (entry_x - v_x[(orig_y+b_size)*step+orig_x])*(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + (entry_y - v_y[(orig_y+b_size)*step+orig_x])*(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_array[3] = (entry_x - v_x[orig_y*step+(orig_x+b_size)])*(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + (entry_y - v_y[orig_y*step+(orig_x+b_size)])*(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  return (min_array_val(sum_array, 4));

}
double Pair_Motion::vector_cost2(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block
  
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];
   
  //function 2 in paper
  double sum_array[3];
  sum_array[0] = (entry_x - v_x[orig_y*step+(orig_x-b_size)])*(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + (entry_y - v_y[orig_y*step+(orig_x-b_size)])*(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_array[1] = (entry_x - v_x[(orig_y+b_size)*step+orig_x])*(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + (entry_y - v_y[(orig_y+b_size)*step+orig_x])*(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_array[2] = (entry_x - v_x[orig_y*step+(orig_x+b_size)])*(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + (entry_y - v_y[orig_y*step+(orig_x+b_size)])*(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  return (min_array_val(sum_array, 3)); 

}
double Pair_Motion::vector_cost3(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block
 
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];
    
  //function 2 in paper
  double sum_array[3];
  sum_array[0] = (entry_x - v_x[(orig_y-b_size)*step+orig_x])*(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + (entry_y - v_y[(orig_y-b_size)*step+orig_x])*(entry_y - v_y[(orig_y-b_size)*step+orig_x]); 
  sum_array[1] = (entry_x - v_x[(orig_y+b_size)*step+orig_x])*(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + (entry_y - v_y[(orig_y+b_size)*step+orig_x])*(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  sum_array[2] = (entry_x - v_x[orig_y*step+(orig_x+b_size)])*(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + (entry_y - v_y[orig_y*step+(orig_x+b_size)])*(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  return (min_array_val(sum_array, 3));

}
double Pair_Motion::vector_cost4(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block
  
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];
     
  //function 2 in paper
  double sum_array[3];
  sum_array[0] = (entry_x - v_x[(orig_y-b_size)*step+orig_x])*(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + (entry_y - v_y[(orig_y-b_size)*step+orig_x])*(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_array[1] = (entry_x - v_x[orig_y*step+(orig_x-b_size)])*(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + (entry_y - v_y[orig_y*step+(orig_x-b_size)])*(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_array[2] = (entry_x - v_x[orig_y*step+(orig_x+b_size)])*(entry_x - v_x[orig_y*step+(orig_x+b_size)]) + (entry_y - v_y[orig_y*step+(orig_x+b_size)])*(entry_y - v_y[orig_y*step+(orig_x+b_size)]);
  return (min_array_val(sum_array, 3)); 

}
double Pair_Motion::vector_cost5(int pos_y, int pos_x, int orig_y, int orig_x)
{ 
  // posx, pos y are for the neighboring block we will try
  // orig_x and orig_y are for the original center block
 
  int entry_x = v_x[pos_y*step+pos_x];
  int entry_y = v_y[pos_y*step+pos_x];
   
  //function 2 in paper
  double sum_array[3];
  sum_array[0] = (entry_x - v_x[(orig_y-b_size)*step+orig_x])*(entry_x - v_x[(orig_y-b_size)*step+orig_x]) + (entry_y - v_y[(orig_y-b_size)*step+orig_x])*(entry_y - v_y[(orig_y-b_size)*step+orig_x]);
  sum_array[1] = (entry_x - v_x[orig_y*step+(orig_x-b_size)])*(entry_x - v_x[orig_y*step+(orig_x-b_size)]) + (entry_y - v_y[orig_y*step+(orig_x-b_size)])*(entry_y - v_y[orig_y*step+(orig_x-b_size)]);
  sum_array[2] = (entry_x - v_x[(orig_y+b_size)*step+orig_x])*(entry_x - v_x[(orig_y+b_size)*step+orig_x]) + (entry_y - v_y[(orig_y+b_size)*step+orig_x])*(entry_y - v_y[(orig_y+b_size)*step+orig_x]);
  return (min_array_val(sum_array, 3));

}
int Pair_Motion::calculate_block_value(int y, int x)
{  
  int block_sum = 0;
  x = x - start_pos;
  y = y - start_pos;
      
  for (int i = y; i < y + b_size; i++) //this is the height of the block
  {
    for (int j = x; j < x + b_size; j++) //this is the width of the block
	{
	  block_sum = block_sum + (uchar)level1b_grey->imageData[i*level1b_grey->widthStep+j]; //add up pixels in block	  
	}		
  }  
  return block_sum;
}
double Pair_Motion::angle_similarity(int y, int x, int i, int j)
{  
  double dot_product = v_x[i*step+j]*v_x[y*step+x] + v_y[i*step+j]*v_y[y*step+x];
  double mag_product = sqrt(double(v_x[i*step+j]*v_x[i*step+j] + v_y[i*step+j]*v_y[i*step+j]))*sqrt(double(v_x[y*step+x]*v_x[y*step+x] + v_y[y*step+x]*v_y[y*step+x]));
  if (mag_product == 0)
    return 1.0;
  else
    return (dot_product / mag_product);
}
double Pair_Motion::min_array(double * array_vals, int size, int &min_index)
{
   double min_value = array_vals[0]; // initialize 
   min_index = 0;
   
   for (int i = 1; i < size; i++)
   {
     if(min_value > array_vals[i])
     {
	    min_value = array_vals[i];
	    min_index = i;
     }
   }

   /*for(int i = 0; i < size; i++)
   {
     if (array_vals[i] == min_value)
	   total_count++;
   }

   total_count = total_count - 1; //because at least one value should be the same.
   */

   return min_value;
}
int Pair_Motion::min_array_val(int * array_vals, int size)
{
   int min_value = array_vals[0]; // initialize 
   
   for (int i = 1; i < size; i++)
   {
     if(min_value > array_vals[i])
     {
	    min_value = array_vals[i];
     }
   }
   return min_value;
}
double Pair_Motion::min_array_val(double * array_vals, int size)
{
   double min_value = array_vals[0]; // initialize 
   
   for (int i = 1; i < size; i++)
   {
     if(min_value > array_vals[i])
     {
	    min_value = array_vals[i];
     }
   }
   return min_value;
}
double Pair_Motion::max_array_val(double * array_vals, int size)
{
   double max_value = array_vals[0]; // initialize 
   
   for (int i = 1; i < size; i++)
   {
     if(max_value < array_vals[i])
     {
	    max_value = array_vals[i];
     }
   }
   return max_value;
}
int Pair_Motion::opt_med5(int * p)
{
    PIX_SORT(p[0],p[1]) ; PIX_SORT(p[3],p[4]) ; PIX_SORT(p[0],p[3]) ;
    PIX_SORT(p[1],p[4]) ; PIX_SORT(p[1],p[2]) ; PIX_SORT(p[2],p[3]) ;
    PIX_SORT(p[1],p[2]) ; return(p[2]) ;
}
int Pair_Motion::opt_med9(int * p)
{
    PIX_SORT(p[1], p[2]) ; PIX_SORT(p[4], p[5]) ; PIX_SORT(p[7], p[8]) ;
    PIX_SORT(p[0], p[1]) ; PIX_SORT(p[3], p[4]) ; PIX_SORT(p[6], p[7]) ;
    PIX_SORT(p[1], p[2]) ; PIX_SORT(p[4], p[5]) ; PIX_SORT(p[7], p[8]) ;
    PIX_SORT(p[0], p[3]) ; PIX_SORT(p[5], p[8]) ; PIX_SORT(p[4], p[7]) ;
    PIX_SORT(p[3], p[6]) ; PIX_SORT(p[1], p[4]) ; PIX_SORT(p[2], p[5]) ;
    PIX_SORT(p[4], p[7]) ; PIX_SORT(p[4], p[2]) ; PIX_SORT(p[6], p[4]) ;
    PIX_SORT(p[4], p[2]) ; return(p[4]) ;
}
void Pair_Motion::PIX_SORT(int a, int b)
{ 
  if (a > b)
    PIX_SWAP(a,b); 
}
void Pair_Motion::PIX_SWAP(int a, int b)
{
  int temp=a;
  a=b;
  b=temp; 
}
void Pair_Motion::draw_motion_vectors()
{    	
  //reset_draw_image();

  //to show a specific level
  /*b_size = level2_block_size;
  start_pos = start_pos_level2;
  step = step2;
  width = width2;
  height = height2;
  v_x = v_x2;
  v_y = v_y2;*/
  //end of specific level code

  //b_size = 7;
  //int do_once = 0;

  // I found and slightly modified this cool code on the internet for drawing a flow field
  for (int i = start_pos; i < height - (add_height - start_pos)/*width-(start_pos + v_shift)*/; i+=(b_size*10))// << 1))  
  {
    for (int j = start_pos; j < width - (add_width - start_pos)/*height-(start_pos + v_shift)*/; j+=(b_size*10))// << 1))  
    {     
	  int line_thickness = 1;	
	  
	  // CV_RGB(red, green, blue) is the red, green, and blue components
	  // of the color you want, each out of 255.
	
	  CvScalar line_color;//, line_color2, line_color3, line_color4;
	  line_color = CV_RGB(0,0,255); //CV_RGB(255,255,255);
	 
	  CvPoint p,q;//,r,s;

	  p.x = j;// + (b_size_div-1);
	  p.y = i;// + (b_size_div-1); //I suppose I could put a non-integer value for the position (like
		                        //if j = 4 and bsize = 4, then j + (b_size_div-1) = 5, but center of 
	                            //block is really 5.5

	  q.x = v_x[i*step+j] + p.x;
	  q.y = v_y[i*step+j] + p.y;

	  if (q.x < 0 || q.y < 0 || q.x >= (width - (add_width - start_pos)) || q.y >= (height - (add_height - start_pos)))
	    continue;
	  
	  // Now we draw the main line of the arrow.
	  // "frame1" is the frame to draw on.
	  // "p" is the point where the line begins.
	  // "q" is the point where the line stops.
	  // "CV_AA" means antialiased drawing.
	  // "0" means no fractional bits in the center cooridinate or radius.

	  //if (v_x[i*step+j] == 0 && v_y[i*step+j] == 0)
	  //continue;
		
	  cvLine(debug, p, q, line_color, line_thickness, CV_AA, 0);
	  //cvCircle(debug, p, 1, CV_RGB(255,0,0), -1);	  


	  //Let's draw the 9 boxes for each of the candidates

	  //candidate 0 
	 // r.x = v_y[192*step+224] + 224;//192; //q.x; //224;
	 // r.y = v_x[192*step+224] + 192;//224;  //192;
	  //s.x = r.x + 3;
	 // s.y = r.y + 3;
	  //cvRectangle(debug, r, s, CV_RGB(9,56,191), 1, 8, 0);	  

	/*  //candidate 1
	  r.x = v_y[188*step+224] + 224; //q.x; //224;
	  r.y = v_x[188*step+224] + 188;  //192;
	  s.x = r.x + 3;
	  s.y = r.y + 3;
	  cvRectangle(debug, r, s, CV_RGB(9,56,191), 1, 8, 0);	

	  //candidate 2
	  r.x = v_y[192*step+220] + 220; //q.x; //224;
	  r.y = v_x[192*step+220] + 192;  //192;
	  s.x = r.x + 3;
	  s.y = r.y + 3;
	  cvRectangle(debug, r, s, CV_RGB(9,56,191), 1, 8, 0);	

	  //candidate 3
	  r.x = v_y[196*step+224] + 224; //q.x; //224;
	  r.y = v_x[196*step+224] + 196;  //192;
	  s.x = r.x + 3;
	  s.y = r.y + 3;
	  cvRectangle(debug, r, s, CV_RGB(9,56,191), 1, 8, 0);	

	  //candidate 4
	  r.x = v_y[192*step+228] + 228; //q.x; //224;
	  r.y = v_x[192*step+228] + 192;  //192;
	  s.x = r.x + 3;
	  s.y = r.y + 3;
	  cvRectangle(debug, r, s, CV_RGB(9,56,191), 1, 8, 0);	

	  //candidate 5
	  r.x = v_y[188*step+220] + 220; //q.x; //224;
	  r.y = v_x[188*step+220] + 188;  //192;
	  s.x = r.x + 3;
	  s.y = r.y + 3;
	  cvRectangle(debug, r, s, CV_RGB(9,56,191), 1, 8, 0);	*/

	  //candidate 6
	 // r.x = v_y[196*step+220] + 220; //q.x; //224;
	  //r.y = v_x[196*step+220] + 196;  //192;
	  //s.x = r.x + 3;
	  //s.y = r.y + 3;
	  //cvRectangle(debug, r, s, CV_RGB(0,255,0), 1, 8, 0);	

	  /*//candidate 7
	  r.x = v_y[188*step+228] + 228; //q.x; //224;
	  r.y = v_x[188*step+228] + 188;  //192;
	  s.x = r.x + 3;
	  s.y = r.y + 3;
	  cvRectangle(debug, r, s, CV_RGB(9,56,191), 1, 8, 0);	

	  //candidate 8
	  r.x = v_y[196*step+228] + 228; //q.x; //224;
	  r.y = v_x[196*step+228] + 196;  //192;
	  s.x = r.x + 3;
	  s.y = r.y + 3;
	  cvRectangle(debug, r, s, CV_RGB(9,56,191), 1, 8, 0);	*/

	  //r.x = 224; //256 //224
	  //r.y = 192;  //180	   //192

	  //s.x = 227;
	  //s.y = 195;

	  //cvRectangle(debug, r, s, CV_RGB(255,0,0), 1, 8, 0);	  	  

	  /*if(do_once == 0){
	  CvPoint r,s;
		  r.x = j - ((64/2) - 1);
		  r.y = i - ((64/2) - 1);
		  s.x = j + (64/2);
		  s.y = i + (64/2);
		  cvRectangle(debug, r, s, CV_RGB(117,47,156), 1, 8, 0);}
	  do_once = 1;*/

	}
  }	
  //IplImage *debug_zoom = cvCreateImage(cvSize((level2a->width*2), (level2a->height*2)), IPL_DEPTH_8U, 1);
  //cvResize(debug, debug_zoom, CV_INTER_LINEAR);
  
}
void Pair_Motion::drawblock_currentframe()
{    	
	   
	int line_thickness = 1;	
	  
	// CV_RGB(red, green, blue) is the red, green, and blue components
	// of the color you want, each out of 255.
	
	CvScalar line_color;//, line_color2, line_color3, line_color4;
	line_color = CV_RGB(117,47,156); //CV_RGB(255,255,255);
	 
	CvPoint r,s;

	//p.y = i;// + (b_size_div-1);
	//p.x = j;// + (b_size_div-1); //I suppose I could put a non-integer value for the position (like
		                    //if j = 4 and bsize = 4, then j + (b_size_div-1) = 5, but center of 
	                        //block is really 5.5

	//q.x = v_y[i*step+j] + p.x;
	//q.y = v_x[i*step+j] + p.y;
	  
	// Now we draw the main line of the arrow.
	// "frame1" is the frame to draw on.
	// "p" is the point where the line begins.
	// "q" is the point where the line stops.
	// "CV_AA" means antialiased drawing.
	// "0" means no fractional bits in the center cooridinate or radius.
		
	//cvLine(debug, p, q, line_color, line_thickness, CV_AA, 0);
	//cvCircle(debug, p, 1, CV_RGB(117,47,156), -1);	  

	//Let's draw the 9 boxes for each of the candidates

	//candidate 0 
	//r.x = v_y[192*step+224] + 192; //q.x; //224;
	//r.y = v_x[192*step+224] + 224;  //192;
	//s.x = r.x + 3;
	//s.y = r.y + 3;
	//cvRectangle(debug2, r, s, CV_RGB(9,56,191), 1, 8, 0);	
	
	r.x = 224; //256 //224
	r.y = 192;  //180	   //192

	s.x = 227;
	s.y = 195;
	cvRectangle(debug2, r, s, CV_RGB(255,0,0), 1, 8, 0);	  	 

  
}
void Pair_Motion::draw_test_boxes()
{
  //int count = 0;
  int t = 0;
  CvPoint p, q, r, s;
  CvFont font;
  double hScale=0.25;
  double vScale=0.25;
  int lineWidth=1;
  cvInitFont(&font,CV_FONT_HERSHEY_SIMPLEX|CV_FONT_ITALIC, hScale,vScale,0,lineWidth);

  //char num[5];
  b_size = 8;
   
  for (int i = start_pos; i < height - (add_height - start_pos); i+=b_size)//(b_size)//+8))// << 1))  
  {
    for (int j = start_pos; j < width - (add_width - start_pos); j+=b_size)//(b_size)//+8))// << 1))  
    {     
	  //int line_thickness = 1;		
	  
	  // CV_RGB(red, green, blue) is the red, green, and blue components
	  // of the color you want, each out of 255.
	
	  //CvScalar line_color;//, line_color2, line_color3, line_color4;
	  //line_color = CV_RGB(117,47,156); //CV_RGB(255,255,255);
	 
	 

	  p.x = j;// + (b_size_div-1);
	  p.y = i;// + (b_size_div-1); //I suppose I could put a non-integer value for the position (like
		                        //if j = 4 and bsize = 4, then j + (b_size_div-1) = 5, but center of 
	                            //block is really 5.5

	  q.x = v_x[i*step+j] + p.x; //position in frame 1
	  q.y = v_y[i*step+j] + p.y; //position in frame 1
	  
	  // Now we draw the main line of the arrow.
	  // "frame1" is the frame to draw on.
	  // "p" is the point where the line begins.
	  // "q" is the point where the line stops.
	  // "CV_AA" means antialiased drawing.
	  // "0" means no fractional bits in the center cooridinate or radius.

	  r.x = q.x + (b_size-1);
	  r.y = q.y + (b_size-1);	

	  //For placing letters
	  s.x = q.x + 1;
	  s.y = r.y - 1;

	  /*
	  cvPoint s, t, u, v, w, x, y, z, a, b, c, d, e, f, g, h;
	  int shift = 2;

	  s.x = v_x[(i-shift)*step+(j-shift)] + (j-shift);
	  s.y = v_y[(i-shift)*step+(j-shift)] + (i-shift);

	  t.x = s.x + 1;
	  t.y = s.y + 1;

	  u.x = v_x[(i-shift)*step+j] + j;
	  u.y = v_y[(i-shift)*step+j] + (i-shift);

	  v.x = u.x + 1;
	  v.y = u.y + 1;

	  w.x = v_x[(i-shift)*step+(j+shift)] + (j+shift);
	  w.y = v_y[(i-shift)*step+(j+shift)] + (i-shift);

	  x.x = w.x + 1;
	  x.y = w.y + 1;

	  y.x = v_x[i*step+(j-shift)] + (j-shift);
	  y.y = v_y[i*step+(j-shift)] + i;

	  z.x = y.x + 1;
	  z.y = y.y + 1;

	  a.x = v_x[i*step+(j+shift)] + (j+shift);
	  a.y = v_y[i*step+(j+shift)] + i;

	  b.x = a.x + 1;
	  b.y = a.y + 1;

	  c.x = v_x[(i+shift)*step+(j-shift)] + (j-shift);
	  c.y = v_y[(i+shift)*step+(j-shift)] + (i+shift);

	  d.x = c.x + 1;
	  d.y = c.y + 1;

	  e.x = v_x[(i+shift)*step+j] + j;
	  e.y = v_y[(i+shift)*step+j] + (i+shift);

	  f.x = e.x + 1;
	  f.y = e.y + 1;

	  g.x = v_x[(i+shift)*step+(j+shift)] + (j+shift);
	  g.y = v_y[(i+shift)*step+(j+shift)] + (i+shift);

	  h.x = g.x + 1;
	  h.y = g.y + 1; */

	 //if (Computed_Data.overlap[((q.y >> 1) << 1)*step+((q.x >> 1) << 1)] > 12 || Computed_Data.overlap[(((q.y >> 1) << 1)+2)*step+((q.x >> 1) << 1)] > 12 || Computed_Data.overlap[((q.y >> 1) << 1)*step+(((q.x >> 1) << 1)+2)] > 12 || Computed_Data.overlap[(((q.y >> 1) << 1) + 2)*step+(((q.x >> 1) << 1) + 2)] > 12)	
	   // cvRectangle(debug, q, r, CV_RGB(117,47,156), 1, 8, 0);
	  //else
	  //if (Computed_Data.reliability[i*step+j] < 0.5)
	    //cvRectangle(debug, q, r, CV_RGB(19,138,214), 1/*CV_FILLED*/, 8, 0); //unreliable
		
		/*if (Computed_Data.reliability[i*step+j] > 0.6)
		{
		  cvRectangle(debug, s, t, CV_RGB(255,0,0), 1, 8, 0); //unreliable
		  cvRectangle(debug, u, v, CV_RGB(255,0,0), 1, 8, 0); //unreliable
		  cvRectangle(debug, w, x, CV_RGB(255,0,0), 1/, 8, 0); //unreliable
		  cvRectangle(debug, y, x, CV_RGB(255,0,0), 1, 8, 0); //unreliable
		  cvRectangle(debug, q, r, CV_RGB(19,138,214), 1, 8, 0); //unreliable
		  cvRectangle(debug, a, b, CV_RGB(255,0,0), 1, 8, 0); //unreliable
		  cvRectangle(debug, c, d, CV_RGB(255,0,0), 1, 8, 0); //unreliable
		  cvRectangle(debug, e, f, CV_RGB(255,0,0), 1, 8, 0); //unreliable
		  cvRectangle(debug, g, h, CV_RGB(255,0,0), 1, 8, 0); //unreliable
		}*/
	  //b_size = 2;
	    //if ((Computed_Data.reliability[i*step+j] <= 0.5))// && (calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i) > 10))

	  /*
	  double overlap_center = 0;
	  double overlap_left = 0;
	  double overlap_right = 0;
	  double overlap_top = 0;
	  double overlap_bottom = 0;

	  for (int k = q.y; k < q.y + b_size; k++)
	    for(int l = q.x; l < q.x + b_size; l++)
		  overlap_center += Computed_Data.overlap[k*step+l]; //calculate total overlap for the whole block

	  for (int k = q.y; k < q.y + b_size; k++)
	    for(int l = q.x-b_size; l < q.x; l++)
		  overlap_left += Computed_Data.overlap[k*step+l]; //calculate total overlap for the whole block

	  for (int k = q.y; k < q.y + b_size; k++)
	    for(int l = q.x + b_size; l < q.x + 2*b_size; l++)
		  overlap_right += Computed_Data.overlap[k*step+l]; //calculate total overlap for the whole block

	  for (int k = q.y - b_size; k < q.y; k++)
	    for(int l = q.x; l < q.x + b_size; l++)
		  overlap_top += Computed_Data.overlap[k*step+l]; //calculate total overlap for the whole block

	  for (int k = q.y + b_size; k < q.y + 2*b_size; k++)
	    for(int l = q.x; l < q.x + b_size; l++)
		  overlap_bottom += Computed_Data.overlap[k*step+l]; //calculate total overlap for the whole block
	  

	  if ((overlap_center > (b_size*b_size + 1)) || (((overlap_left > 5) &&  (overlap_right > 5)) || ((overlap_top > 5) && (overlap_bottom > 5))))// || ((overlap_left + overlap_right + overlap_top + overlap_bottom) > 16))//4*b_size*b_size))
	  {
		 cvRectangle(debug, p, r, CV_RGB(19,138,214), 1, 8, 0); 
		 t++; 
	  }*/

     // if (mult_vec[i*step+j].dup == 1)//(Computed_Data.reliability[i*step+j] < 0.5)// && (calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i) > 0))
	  //{
        if (Computed_Data.reliability[i*step+j] < 0.01)
	    cvRectangle(debug, q, r, CV_RGB(0,100,0), 1, 8, 0); 
		//itoa(t, num, 10);
		//cvPutText(debug, num, s, &font, CV_RGB(0,100,0));
		//t++;
	  //}

	 // else if (

	 // if ((overlap_center <= (b_size*b_size + 4)) && ((overlap_left + overlap_right + overlap_top + overlap_bottom) >= 4*b_size*b_size) && ((overlap_left + overlap_right + overlap_top + overlap_bottom) <= 4*b_size*b_size + 16)) //reliable
	    //t++;//cvRectangle(debug, q, r, CV_RGB(19,138,214), 1, 8, 0); 
	  //else
		//cvRectangle(debug, q, r, CV_RGB(19,138,214), 1, 8, 0); 
	
		//b_size = 1;
	    //if (Computed_Data.reliability[i*step+j] < 0.8)
     //	{
		  //count++;
	      //if((Computed_Data.overlap[i*step+j] == 0) || (Computed_Data.overlap[i*step+j] > 1))
			//cvCircle(debug, q, 1, CV_RGB(19,138,214), 1, 8, 0);
		    //cvRectangle(debug, q, r, CV_RGB(19,138,214), 1, 8, 0); //unreliable
		//}

		//if(calculate_SAD(i, j, v_x[i*step+j] + j, v_y[i*step+j] + i) > 30) 
		  //cvRectangle(debug, q, r, CV_RGB(255,0,0), 1, 8, 0);
	      //cvCircle(debug, s, 1, CV_RGB(255,0,0), -1);	  
	}
  }  
  //std::cout << t << std::endl;
}
void Pair_Motion::frame_prediction()
{
  int i,j,k,l,m,n;

  //to show a specific level
  /*b_size = level2_block_size >> 2;
  start_pos = start_pos_level2;
  step = step2;
  width = width2;
  height = height2;
  v_x = v_x2;
  v_y = v_y2;*/
  //end of spec level code
  //b_size = b_size >> 1;

  if (frame_predict->nChannels == 3)
  {

  for (i = start_pos; i < height - (add_height - start_pos + b_size)/*(start_pos + v_shift)*/; i+=b_size)//goes through all vertical pixels in current frame
  {
	for (j = start_pos; j < width - (add_width - start_pos + b_size)/*(start_pos + v_shift)*/; j+=b_size) //goes through all horizontal pixels in current frame
    {
      l = (v_x[i*step+j] + j);// - (block_size_div-1); //these are frame1 positions
	  k = (v_y[i*step+j] + i);// - (block_size_div-1);

	  if (l < 0 || k < 0 || l >= (width - (add_width - start_pos)) || k >= (height - (add_height - start_pos)))
	    continue;
      
	  for(n = i/* - (block_size_div-1)*/; n < i + b_size; /* (block_size_div+1);*/ n++)  
	  {
	    for(m = j/* - (block_size_div-1)*/; m < j + b_size; /*(block_size_div+1);*/ m++)
		{
          //so basically new position m will get the value of position k in old frame
	  
 	      data_predict[(n-start_pos)*predict_step+(m-start_pos)*3] = data_color1[k*color_step+l*3];	//3 is # channels
	      data_predict[(n-start_pos)*predict_step+(m-start_pos)*3 + 1] = data_color1[k*color_step+l*3 + 1];
	      data_predict[(n-start_pos)*predict_step+(m-start_pos)*3 + 2] = data_color1[k*color_step+l*3 + 2];
		  l++;
		}
		l = (v_x[i*step+j] + j); // - (block_size_div-1); //these are frame1 positions
		k++;
	  }	
	}
  }  

  }
  else
  {

  for (i = start_pos; i < height - (add_height - start_pos + b_size); i+=b_size)//goes through all vertical pixels in current frame
  {
	for (j = start_pos; j < width - (add_width - start_pos + b_size); j+=b_size) //goes through all horizontal pixels in current frame
    {
      l = (v_x[i*step+j] + j);// - (block_size_div-1); //these are frame1 positions
	  k = (v_y[i*step+j] + i);// - (block_size_div-1);

	  if (l < 0 || k < 0 || l >= (width - (add_width - start_pos)) || k >= (height - (add_height - start_pos)))
	    continue;
      
	  for(n = i/* - (block_size_div-1)*/; n < i + b_size; /* (block_size_div+1);*/ n++)  
	  {
	    for(m = j/* - (block_size_div-1)*/; m < j + b_size; /*(block_size_div+1);*/ m++)
		{
          //so basically new position m will get the value of position k in old frame
	  
 	      data_predict[(n-start_pos)*predict_step+(m-start_pos)] = data_color1[k*color_step+l];	
		  l++;
		}
		l = (v_x[i*step+j] + j); // - (block_size_div-1); //these are frame1 positions
		k++;
	  }	
	}
  }  

  }
  //cvShowImage("Motion Compensated Frame", frame_predict);
}
void Pair_Motion::clean_up_nointerp()
{  
  delete [] Computed_Data.v_x1; 
  delete [] Computed_Data.v_y1; 
  delete [] Computed_Data.reliability; 
  delete [] overlap1;
  delete [] overlap2;
  delete [] overlap3;
  delete [] overlap4;
  delete [] v_x2;
  delete [] v_y2; 
  delete [] v_x3;
  delete [] v_y3; 
  delete [] v_x4;
  delete [] v_y4; 

  if(level1a->nChannels == 3)
  {
    cvReleaseImage(&level1a_pad_color); 
    cvReleaseImage(&level1b_pad_color);
  }

  cvReleaseImage(&level1a);
  cvReleaseImage(&level1b);

  cvReleaseImage(&debug);
  //cvReleaseImage(&debug2);
  cvReleaseImage(&frame_predict);

  cvReleaseImage(&level2a_grey);
  cvReleaseImage(&level2b_grey);

  cvReleaseImage(&level1a_pad);
  cvReleaseImage(&level1b_pad);

  cvReleaseImage(&level2a_pad);
  cvReleaseImage(&level2b_pad);

  cvReleaseImage(&level3a_pad);
  cvReleaseImage(&level3b_pad);

  cvReleaseImage(&level4a_pad);
  cvReleaseImage(&level4b_pad);
}
double Pair_Motion::min(double x, double y)
{
  return (x < y) ? x : y;
}
int Pair_Motion::min2(int x, int y)
{
  return (x < y) ? x : y;
}
double Pair_Motion::max(double x, double y)
{
  return (x > y) ? x : y;
}
double Pair_Motion::twomin_diff(int * array_vals, int size)
{
   int min_value1 = array_vals[0]; // initialize 
   int min_value2 = array_vals[0];
   
   for (int i = 1; i < size; i++)
   {
     if(min_value1 > array_vals[i])
     {
		min_value2 = min_value1;
	    min_value1 = array_vals[i];
     }
   }
   return double(min_value2 - min_value1);
}
Pair_Motion::~Pair_Motion()
{
  
}

