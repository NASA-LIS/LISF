//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#include <iostream>
#include <vector>
#include <string>

#include "gvf_binary_reader_modis.h"
#include "cppfunmodis.h"
#include "ftn_drv.h"

using namespace std;   

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Program:	read_gvf_binary_modis.cc
	
	Date:		Nov 2013
	
	Authors:	Jayanthi Srikishen/Jonathan Case
	        	
	Purpose:	This program will read in gzipped GVF datafile for
 			ingest into the NASA/LIS
		
	            NOTE: We needed to set the swap byte flag to 1. 
	                  
	Output: 	1D array of dimensions nx * ny.
				
	Modification History:  This file was modified from the original
	source code from NSSL used to read gzipped MRMS precipitation files.
	       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


/************************************************/
/***  F U N C T I O N  P R O T O T Y P E (S)  ***/
/************************************************/


float* gvf_binary_reader_modis(char *vfname, int &nelements, int swap_flag);
                   

/********************************/
/***  M A I N  P R O G R A M  ***/
/********************************/


int FTN(readgvfmodis)(char * infile, int * nx, int * ny, char * flag, float data_1D[], int len)
{

    /*---------------------------------------*/
    /*** 0. Process command-line arguments ***/ 
    /*---------------------------------------*/

/*    cout << infile << endl; */
/*    cout << flag << endl;   */
    
/*    cout << input_file << endl; */
    int swap_flag = atoi(flag);
/*    cout << swap_flag << endl;  */

    /*---------------------------------------------------*/
    /*** 1. Declare several variables needed by reader ***/
    /*---------------------------------------------------*/

     int nelements; 

    float* input_data_1D = 0;

    /*-------------------------*/
    /*** 2. Read input field ***/
    /*-------------------------*/

    //Read file
    nelements= *nx * *ny;
/*    cout << "nelements / nx / ny = "<< nelements<<" / "<< *nx<<" / "<<*ny <<endl; */
    input_data_1D = gvf_binary_reader_modis(infile, nelements, swap_flag);

/*    cout<<"back from read"<<endl; */

    //Perform some error checking
    if(input_data_1D == 0)
    {
      cout<<"+++ERROR: Failed to read "<<infile<<" Exiting!"<<endl;
      return 0;
    }
/*    cout<<"DONE reading file"<<endl<<endl; */
    
    /*----------------------------------------------*/
    /*** 3. Reassign to data_1D[]                 ***/
    /*----------------------------------------------*/

/*    cout<<"Reassigning to data_1D"<<endl; */
    
/*    for (int i = 0; i < nelements-1; i++) { */
    for (int i = 0; i < nelements; i++) {
      data_1D[i] = input_data_1D[i];
    }

    /*------------------------*/
    /*** 4. Free-up Memory, ***/
    /*------------------------*/

/*    cout<<"Freeing up memory"<<endl; */

    //memory clean-up
    if(input_data_1D != 0) delete [] input_data_1D;

    return 0;
    
}//end main function


/**************************/
/*** F U N C T I O N S  ***/
/**************************/

//see gvf_binary_reader_modis.h
