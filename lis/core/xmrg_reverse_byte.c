//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

#ifndef XMRG_RESERVE_BYTE
#define XMRG_RESERVE_BYTE 1
//BOP
// !ROUTINE: reverse_byte_order
// \label{reverse_byte_order}
// !INTERFACE:
void reverse_byte_order(int *in_array,int arraysize)
// !DESCRIPTION:
// This routine reverses byte order for 4-byte integer. Real number has to be 
// casted into (int *) type when calling reverse\_byte\_order.
// This code is adapted from \url{http://www.nws.noaa.gov/oh/hrl/dmip/2/src/read_xmrg2.c}
//EOP
{

    unsigned int   i,k;
    signed char *p_data;   /*data pointer*/
    signed char *p_temp;   /*temporaty data pointer */
    int temp;

    p_data = (signed char *) in_array - 1;
    for ( k = 0 ; k < arraysize ; k++ )
    {
        temp = *( in_array + k );
        p_temp = ( signed char * ) ( &temp ) + 4;

        for( i = 0 ; i < 4 ; i++ )
        {
            *(++p_data) = *(--p_temp);
        }
    }
}

//BOP
// !ROUTINE: reverse_byte_order_short
// \label{reverse_byte_order_short}
// !INTERFACE:
void reverse_byte_order_short(short *in_array,int arraysize)
// !DESCRIPTION:
// This routine reverses byte order for 2-byte integer. 
// This code is adpated from \url{http://www.nws.noaa.gov/oh/hrl/dmip/2/src/read_xmrg2.c}
//EOP
{

    unsigned int   i,k;
    signed char *p_data;   /*data pointer*/
    signed char *p_temp;   /*temporaty data pointer */
    short temp;

    p_data = (signed char *) in_array - 1;
    for ( k = 0 ; k < arraysize ; k++ )
    {
        temp = *( in_array + k );
        p_temp = ( signed char * ) ( &temp ) + 2;

        for  ( i = 0 ; i < 2 ; i++ )
        {
            *(++p_data) = *(--p_temp);
        }
    }
}
#endif
