!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

  ! SY: File indentation performed for integration into LIS
  ! By Rex L. Baum, 1 April 2004
  ! 19Aug2009 RLB added bkgrof
  module input_vars
        logical:: ans,outp(8),rodoc,lskip,lany,llus,lps0,unsat0,bkgrof
        logical:: lasc,lpge0 ! Added 4/14-15/2010
        logical:: igcapf ! added 2/15/2012 RLB
        logical,allocatable:: unsat(:), igcap(:)  ! added igcap(:) 2/15/2012 RLB
        integer:: imax,row,col,nwf,tx,nmax
        integer:: flag,nper,spcg ! Added spcg 2/14/2012 RLB
        integer:: nzs,mmax
        integer:: nzon,nout
        integer,allocatable:: ksav(:), uijz(:) !added uijz, 12/7/2010 RLB
        real:: uww,zmin,t,dep,czmax,crizero,slomin,deepz
        real,allocatable:: ths(:),thr(:),alp(:),dif(:),c(:),phi(:)
        real,allocatable:: ks(:),uws(:),capt(:),cri(:),tsav(:)
        character (len=5):: flowdir, el_or_dep
        character (len=4):: deepwat

  ! SY: Begin listing all the variables common across multiple routines.
  ! SY: Variables are str1_str2 where str1=calling routine name, str2=variabla name in calling routine 
         integer:: trigrs_grd,trigrs_imx1,trigrs_nodata
         integer:: trigrs_ncol,trigrs_nrow,trigrs_sctr
         character (len=14):: trigrs_header(6)
         integer:: trigrs_mnd
  ! SY: End listing all the variables common across multiple routines.
  
  end module input_vars
