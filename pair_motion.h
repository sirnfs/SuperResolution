#ifndef __Pair_Motion__
#define __Pair_Motion__

#include "stdafx.h"

//template <class A_type> class complex_data 
//{
//  public:
//    A_type r;
//	A_type c;
//};

class MF_data { //motion vectors and frame data
  public:
    int * v_x1, * v_y1, *overlap;
	double *reliability;
	IplImage *frame1;
};

class Pair_Motion {
  public:       
    MF_data Computed_Data;	
	IplImage *frame_predict;
	IplImage *debug, *debug2;
	int width2, height2, step1, step2, level2_search_size, level2_block_size, v_shift2, start_pos_level2;
	int width1, height1, level1_search_size, level1_block_size, v_shift1, start_pos_level1; 
	int add_height, add_width, add_height1, add_width1, add_height2, add_width2, add_height3, add_width3, add_height4, add_width4; //this is how much we pad the borders of the images
	
	Pair_Motion(int * search_sizes, int * block_sizes); //initializes both square block sizes and non-square block sizes together.
	void start_calculation_nointerp(IplImage*& image1, IplImage*& image2, int lambda_value); 
	void draw_motion_vectors();	
	void drawblock_currentframe();
	void draw_test_boxes();	
	void frame_prediction();	
	void clean_up_nointerp();
    ~Pair_Motion();
  
  private:	    
	void pad_images(IplImage*& imageA, IplImage*& imageB, IplImage*& imageA_pad, IplImage*& imageB_pad, int start_pos);	
	void pad_image_color(IplImage*& imageA, IplImage*& imageA_pad, int start_pos);	
	void calculate_motion_vectors_overlap(int lambda_value);	
	void calculate_reliability();
	void calculate_pixel_reliability();
	void validity_for_SR();	
    void calculate_block_overlap();
	void calculate_pixel_overlap();
	double find_overlap_max();
	int get_overlap_volume(int x_pos, int y_pos);
	void calculate_singleMV_overlap(int b_size, int pos_y, int pos_x, int& top_left, int& top_right, int& bottom_left, int& bottom_right);
	void onelevl_BM();
	void onelevlspiral_BM();
	void onelevlspiral_BM_minoverlap();
	void nextlevl_BM();
	void nextlevlspiral_BM();
	void nextlevlspiral_BM_minoverlap(int start_direction);
	void nextlevlspiral_BM_weightedoverlap(int start_direction);
	void nextlevlspiral_BM_adapt();
	void setMVs_iter();
	void copyto_oldMVs();
	int calculate_SAD(int y, int x, int l, int k);
	int calculate_local_SAD(int y, int x, int l, int k);
	double calculate_SAD2(int y, int x, int l, int k);
	void refine_MVs();	
	void add_smoothness8_old(double lambda_value);
	void add_smoothness8_overlap(double lambda_value); 
	double calc_derivX(int k, int l);
    double calc_derivY(int k, int l);
	double min_array(double * array_vals, int size, int &min_index);
	int min_array_val(int * array_vals, int size);
	double min_array_val(double * array_vals, int size);
	double max_array_val(double * array_vals, int size);
	int opt_med5(int * p);
	int opt_med9(int * p);
    void PIX_SORT(int a, int b);
    void PIX_SWAP(int a, int b);
	double cost_fn(int v_yval, int v_xval, int orig_y, int orig_x);
	double cost_fn_dynamic_smoothness(int v_yval, int v_xval, int orig_y, int orig_x);
	double cost_fn_top_overlap(int v_yval, int v_xval, int orig_y, int orig_x);
	double cost_fn_left_overlap(int v_yval, int v_xval, int orig_y, int orig_x);
	double cost_fn_right_overlap(int v_yval, int v_xval, int orig_y, int orig_x);
	double cost_fn_bottom_overlap(int v_yval, int v_xval, int orig_y, int orig_x);
	double cost_fn_leftcorn_overlap(int v_yval, int v_xval, int orig_y, int orig_x);
	double cost_fn_rightcorn_overlap(int v_yval, int v_xval, int orig_y, int orig_x);
	double cost_fn_bottomleftcorn_overlap(int v_yval, int v_xval, int orig_y, int orig_x);
	double cost_fn_bottomrightcorn_overlap(int v_yval, int v_xval, int orig_y, int orig_x);	
	double cost_fn8(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8top(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8bottom(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8left(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8right(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8leftcorn(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8rightcorn(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8bottomleftcorn(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8bottomrightcorn(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8_weighted(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8a(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8a_weighted(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8b(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8b_weighted(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8c(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8c_weighted(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8d(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8d_weighted(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8e(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8e_weighted(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8f(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8f_weighted(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8g(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8g_weighted(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8h(int pos_y, int pos_x, int orig_y, int orig_x);
	double cost_fn8h_weighted(int pos_y, int pos_x, int orig_y, int orig_x);
	double vector_cost1(int pos_y, int pos_x, int orig_y, int orig_x);
	double vector_cost2(int pos_y, int pos_x, int orig_y, int orig_x);
	double vector_cost3(int pos_y, int pos_x, int orig_y, int orig_x);
	double vector_cost4(int pos_y, int pos_x, int orig_y, int orig_x);
	double vector_cost5(int pos_y, int pos_x, int orig_y, int orig_x);
	int calculate_block_value(int y, int x);
	double angle_similarity(int y, int x, int i, int j);
	double min(double x, double y);
	int min2(int x, int y);
	double max(double x, double y);
	double twomin_diff(int * array_vals, int size);

    uchar * data_level1a, * data_level1b; //Level 1 -- highest resolution
	int lambda1;
	double * reliability1; 
    int * overlap1; 
	IplImage *level1a, *level1b, *level1a_pad, *level1a_pad_color, *level1b_pad_color, *level1b_pad, *level1a_grey, *level1b_grey;
		
	IplImage *frame2, *frame2_interp, *level2a_pad, *level2b_pad, *level2a_grey, *level2b_grey;  //Level 2	 
	uchar * data_level2a, * data_level2b;
	int lambda2;
	int * v_x2, * v_y2;
	int * overlap2;
				
	int width3, height3, step3, level3_search_size, level3_block_size, v_shift3; //Level 3 
	uchar * data_level3a, * data_level3b;
	int lambda3, start_pos_level3;
	int * v_x3, * v_y3;
	int * overlap3;
	IplImage *level3a_grey, *level3b_grey, *level3a_pad, *level3b_pad;

	int width4, height4, step4, level4_search_size, level4_block_size, v_shift4; //Level 4 -- lowest resolution
	uchar * data_level4a, * data_level4b;
	int lambda4, start_pos_level4;
	int * v_x4, * v_y4;	
	int * overlap4;
	IplImage *level4a_grey, *level4b_grey, *level4a_pad, *level4b_pad;

	//These are going to hold the arguments for each level while it's been worked on
	int width, height, width_previous, height_previous, predict_step, color_step, step, step_previous, search_size, b_size, b_size_previous, v_shift, v_shift_previous, start_pos_prev; //Holds current level values
	uchar * data_a, * data_b, * data_predict, *data_color1;
	int lambda, start_pos;
	int * v_x, * v_y, * v_xprev, * v_yprev;
	std::ofstream fout; //for debugging		
	int level; //debugging

};

#endif