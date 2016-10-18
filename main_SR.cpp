#include "stdafx.h"
#include "registration.h"
#include "POCS.h"

int _tmain(int argc, _TCHAR* argv[])
{	
  clock_t start, end; //Variables for keeping track of the time elapsed
  const int levels = 4; //Number of levels in pyramid
  int block_size[levels] = {32, 32, 16, 8}; //arranged from lowest resolution level to highest resolution level
  int search_size[levels] = {12, 16, 24, 12}; //arranged from lowest resolution level to highest resolution level

  //2-D parameters for 3x Magnification
  //const int levels = 4; //Number of levels in pyramid
  //int block_size[levels] = {16, 32, 32, 64}; //arranged from lowest resolution level to highest resolution level
  //int search_size[levels] = {10, 5, 16, 5}; //arranged from lowest resolution level to highest resolution level

  //const int levels = 3; //Number of levels in pyramid
  //int block_size[levels] = {32, 32, 16}; //arranged from lowest resolution level to highest resolution level
  //int search_size[levels] = {10, 16, 24}; //arranged from lowest resolution level to highest resolution level

  int interp_factor = 4; //Factor to multiple dimensions of original volume by
  int validity_threshold = 30000;//12;////30;//80;//100;//2560; (pixel-wise SAD of 5 + 1 = 6 * (overlap of 2) = 12
  int max_intensity = 255;

  //Create classes for each pair of frames that we will compute the motion of
  //NOTE:  PERFORMANCE MODIFICATION NEEDED -- WE REPEAT THE INTERPOLATION OF THE ANCHOR FRAME MULTIPLE TIMES IN THE REGISTRATION CLASS - ONLY NEED TO INTERPOLATE ONCE
  std::vector<Registration*> pair_num;

  //pair_num.push_back(new Registration("./BlockMatchingVerification/frame1.bmp", "./BlockMatchingVerification/frame2.bmp", block_size, search_size, levels, interp_factor, 1)); 

  //pair_num.push_back(new Registration("./BlockMatchingVerification/frame1.bmp", "./BlockMatchingVerification/frame2.bmp", block_size, search_size, levels, interp_factor)); 
  
  //16 frames
  /*pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0015.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0014.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0000.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0013.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0001.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0012.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0002.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0011.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0003.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0010.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0004.bmp", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0009.bmp", block_size, search_size, levels, interp_factor, 1));
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0005.bmp", block_size, search_size, levels, interp_factor, 1));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0008.bmp", block_size, search_size, levels, interp_factor, 1));
  pair_num.push_back(new Registration("./SimulatedSRDataset/EIA_sequence/image0007.bmp", "./SimulatedSRDataset/EIA_sequence/image0006.bmp", block_size, search_size, levels, interp_factor, 1)); */
  
  //16 frames - my synthetic dataset
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img16.tif", block_size, search_size, levels, interp_factor, 15));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img15.tif", block_size, search_size, levels, interp_factor, 14));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img14.tif", block_size, search_size, levels, interp_factor, 13));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img13.tif", block_size, search_size, levels, interp_factor, 12));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img12.tif", block_size, search_size, levels, interp_factor, 11));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img11.tif", block_size, search_size, levels, interp_factor, 10));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img10.tif", block_size, search_size, levels, interp_factor, 9));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img9.tif", block_size, search_size, levels, interp_factor, 8));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img8.tif", block_size, search_size, levels, interp_factor, 7));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img7.tif", block_size, search_size, levels, interp_factor, 6));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img6.tif", block_size, search_size, levels, interp_factor, 5));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img4.tif", block_size, search_size, levels, interp_factor, 4));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img4.tif", block_size, search_size, levels, interp_factor, 3));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img3.tif", block_size, search_size, levels, interp_factor, 2));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img1.tif", "./SimulatedSRDataset/my_dataset/img2.tif", block_size, search_size, levels, interp_factor, 1));  

  /*pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img16.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 15));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img15.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 14));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img14.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 13));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img13.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 12));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img12.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 11));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img11.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 10));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img10.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 9));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img9.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 8));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img8.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 7));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img7.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 6));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img6.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 5));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img5.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 4));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img4.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 3));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img3.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 2));  
  pair_num.push_back(new Registration("./SimulatedSRDataset/my_dataset/img2.tif", "./SimulatedSRDataset/my_dataset/img1.tif", block_size, search_size, levels, interp_factor, 1));  */

  //30 frames
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image00.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image29.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image01.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image28.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image02.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image27.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image03.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image26.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image04.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image25.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image05.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image24.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image06.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image23.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image07.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image22.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image08.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image21.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image09.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image20.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image10.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image19.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image11.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image18.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image12.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image17.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image13.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image16.bmp", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./SimulatedSRDataset/image15.bmp", "./SimulatedSRDataset/image14.bmp", block_size, search_size, levels, interp_factor));  
  


  //pair_num.push_back(new Registration("./3d_volumes/vol3.dcm", "./3d_volumes/vol1.dcm", block_size, search_size, levels, interp_factor)); //skipping this for now
  //pair_num.push_back(new Registration("./3d_volumes/vol3.dcm", "./3d_volumes/vol6.dcm", block_size, search_size, levels, interp_factor));
  //pair_num.push_back(new Registration("./3d_volumes/vol3.dcm", "./3d_volumes/vol5.dcm", block_size, search_size, levels, interp_factor));  
  //pair_num.push_back(new Registration("./3d_volumes/vol3.dcm", "./3d_volumes/vol2.dcm", block_size, search_size, levels, interp_factor));
  //pair_num.push_back(new Registration("./3d_volumes/vol3.dcm", "./3d_volumes/vol4.dcm", block_size, search_size, levels, interp_factor));  

  //7 frames
  /*pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame7.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame7.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame4.tif", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame7.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", block_size, search_size, levels, interp_factor, 1));
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame7.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame5.tif", block_size, search_size, levels, interp_factor, 1));  
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame7.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame8.tif", block_size, search_size, levels, interp_factor, 1));
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame7.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame6.tif", block_size, search_size, levels, interp_factor, 1));  
  */

  //11 frames
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame14.tif", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame4.tif", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame13.tif", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame5.tif", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame12.tif", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame6.tif", block_size, search_size, levels, interp_factor, 1)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame11.tif", block_size, search_size, levels, interp_factor, 1));
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame7.tif", block_size, search_size, levels, interp_factor, 1));  
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", block_size, search_size, levels, interp_factor, 1));
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame8.tif", block_size, search_size, levels, interp_factor, 1));  
  
  //15 frames
  /*pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame17.tif", block_size, search_size, levels, interp_factor)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame3.tif", block_size, search_size, levels, interp_factor)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame16.tif", block_size, search_size, levels, interp_factor)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame4.tif", block_size, search_size, levels, interp_factor)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame15.tif", block_size, search_size, levels, interp_factor)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame5.tif", block_size, search_size, levels, interp_factor)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame14.tif", block_size, search_size, levels, interp_factor)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame6.tif", block_size, search_size, levels, interp_factor)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame13.tif", block_size, search_size, levels, interp_factor)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame7.tif", block_size, search_size, levels, interp_factor)); 
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame12.tif", block_size, search_size, levels, interp_factor));
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame8.tif", block_size, search_size, levels, interp_factor));  
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame11.tif", block_size, search_size, levels, interp_factor));
  pair_num.push_back(new Registration("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame10.tif", "./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", block_size, search_size, levels, interp_factor));  
*/



  //Create POCS class
  POCS HR_POCS("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/IM_0007-Frame9.tif", interp_factor, 8.0, 1.0/*0.40*/, 3.0, 1.0, validity_threshold, max_intensity); //we will create the initial HR estimate (interpolated version of anchor frame) in this class.
  //POCS(char *anchorframe_file, int interp_factor, float d_iter, float variance, int blur_sz, int blur_scale, int v_threshold) //this line is just here to show what the values above are
    
  start = clock(); //Start timer

  //Calculate the motion vectors for each pair, then pass the motion vectors, low-resolution volume, and validity array to the POCS class for fusion.
  int temp_size = pair_num.size(); //this is needed because the vector will change size when im deleting elements
  for(int i = temp_size - 1; i > -1; --i) //skipping the first frame for now -- source of artifacts
  {
	//Perform motion estimation on the frame pair
    pair_num[i]->motion_estimation();
	//std::cout << "Frame done" << std::endl;

	//Call POCS to merge pixels into HR frame -- has current working frame, motion vectors, and validity array of 0s and 1s -- invalid and valid
	HR_POCS.process_frame_reverse(pair_num[i]->current_frame, pair_num[i]->MV_PLNum[levels-1]->MVs, pair_num[i]->valid); //We use the MVs from the highest level of the pyramid
	//HR_POCS.process_frame(pair_num[i]->current_frame, pair_num[i]->MV_PLNum[levels-1]->MVs); //We use the MVs from the highest level of the pyramid

	//Erase the data for the current pair since we are done with it
	delete pair_num[i];
	pair_num.erase(pair_num.end() - 1);
  }

  end = clock(); //Stop timer
  double time_elapsed = double(end - start)/CLOCKS_PER_SEC; //time elapsed
  std::cout << "Time elapsed is " << time_elapsed << std::endl; //output time elapsed

  //Write the output DICOM file after POCS
  //HR_POCS.write_output("./3d_volumes/POCS_after_all.dcm");
  HR_POCS.write_output("./3d_volumes/Corazon_Caro/DICOM/HighTemporal_HighSpatial/SR_image11.dcm");

  //Draw Motion-Compensated Volume
  //pair1.draw_MC_volume();

  //Print MV array
  //pair1.print_MVs();

  //Draw 3-D MVs using Matlab
  //pair1.draw_MVs_MATLAB();

  //Calculate absolute value of pixel differences between volumes
  //pair1.volume_difference(); 

  //Write the output DICOM file
  //pair1.write_output("./3d_volumes/output.dcm"); 

  std::cout << "Done!" << std::endl;
  getchar();

  return 0;
}
