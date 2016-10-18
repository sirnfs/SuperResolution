#include "stdafx.h"
#include "registration.h"

int _tmain(int argc, _TCHAR* argv[])
{	   
  //2-D parameters for 3x Magnification
  const int levels = 4; //Number of levels in pyramid
  int block_size[levels] = {16, 16, 16, 16}; //arranged from lowest resolution level to highest resolution level
  int search_size[levels] = {10, 8, 8, 8}; //arranged from lowest resolution level to highest resolution level
   
  int interp_factor = 4; //Factor to multiple dimensions of original volume by

  Registration Set1("./SimulatedSRDataset/sam2.img", "./SimulatedSRDataset/sam2.img", block_size, search_size, levels, interp_factor);        

  //Write the output DICOM file 
  Set1.write_output("./SimulatedSRDataset/sam2.dcm");

  std::cout << "Done!" << std::endl;
  getchar();

  return 0;
}