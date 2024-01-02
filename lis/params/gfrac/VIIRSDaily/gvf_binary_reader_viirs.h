//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

#ifndef GVF_BINARY_READER_H
#define GVF_BINARY_READER_H

#include <iostream>
#include <fstream>
#include <zlib.h>
#include <vector>
#include <string>
#include <ctime>

using namespace std;


/*------------------------------------------------------------------

	Function:	byteswap (version 1)
		
	Purpose:	Perform byte swap operation on various data
	            types		
				
	Input:		data = some value
					
	Output:		value byte swapped
	
------------------------------------------------------------------*/
template < class Data_Type >
inline void byteswap( Data_Type &data )
{
  unsigned int num_bytes = sizeof( Data_Type );
  char *char_data = reinterpret_cast< char * >( &data );
  char *temp = new char[ num_bytes ];

  for( unsigned int i = 0; i < num_bytes; i++ )
  {
    temp[ i ] = char_data[ num_bytes - i - 1 ];
  }

  for( unsigned int i = 0; i < num_bytes; i++ )
  {
    char_data[ i ] = temp[ i ];
  }
  delete [] temp;
}
  


/*------------------------------------------------------------------

	Function:	byteswap (version 2)
		
	Purpose:	Perform byte swap operation on an array of
	            various data types		
				
	Input:		data = array of values
					
	Output:		array of byte swapped values
	
------------------------------------------------------------------*/
template < class Data_Type >
inline void byteswap( Data_Type *data_array, int num_elements )
{
  int num_bytes = sizeof( Data_Type );
  char *temp = new char[ num_bytes ];
  char *char_data;

  for( int i = 0; i < num_elements; i++ )
  {
    char_data = reinterpret_cast< char * >( &data_array[ i ] );

    for( int i = 0; i < num_bytes; i++ )
    {
      temp[ i ] = char_data[ num_bytes - i - 1 ];
    }

    for( int i = 0; i < num_bytes; i++ )
    {
      char_data[ i ] = temp[ i ];
    }
  }
  delete [] temp;
}



/*------------------------------------------------------------------

	Function: gvf_binary_reader_viirs
		
	Purpose:  Read a GVF binary file and return
              it's header info and data.		
				
	Input:    vfname = input file name and path
	
              swap_flag = flag (= 0 or 1) indicating if values
                          read from the file should be byte 
                          swapped
                
                
	Output:   
	            
	          binary_data = Variable data  in 1D array
	
------------------------------------------------------------------*/
float* gvf_binary_reader_viirs(char *vfname, int &nelements, int swap_flag)
{

    /*--------------------------*/
    /*** 0. Declare variables ***/ 
    /*--------------------------*/
    
    float *binary_data = 0;


    int temp;



    /*-------------------------*/
    /*** 1. Open binary file ***/ 
    /*-------------------------*/
        
    char open_mode[3];
    gzFile   fp_gzip;

    sprintf(open_mode,"%s","rb");
    open_mode[2] = '\0';

    if ( (fp_gzip = gzopen(vfname,open_mode) ) == (gzFile) NULL )
    {
      cout<<"+++ERROR: Could not open "<<vfname<<endl;
      return binary_data;
    }



    /*---------------------------------------*/
    /*** 2. Read binary header information ***/ 
    /*---------------------------------------*/
    



    /*-------------------------*/
    /*** 3. Read binary data ***/ 
    /*-------------------------*/
    
    int num = nelements;
    binary_data = new  float[num];
      
    //read data array
/*    cout<<"before gzread"<<endl; */
    gzread(fp_gzip,binary_data,num*sizeof(float));
/*    cout<<"back from gzread"<<endl; */
    if (swap_flag==1) byteswap(binary_data,num);
      
      
      
    /*------------------------------*/
    /*** 4. Close file and return ***/ 
    /*------------------------------*/
   
    //close file      
    gzclose( fp_gzip );

    return binary_data;

}//end gvf_binary_reader_viirs function

#endif
