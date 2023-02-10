!=============================i==========================================
!  MODULE, TOOLSUBS P.W.LIU, 09/12/18
!  Contains subroutines that are necessry for downscaling program
!-----------------------------------------------------------------------
!  AVHRR_NDVI
!  NOTES : Modified P.W. Liu, 09/12/18
!=======================================================================

#include "LDT_misc.h"

MODULE TOOLSUBS
    USE FUNCTIONS
#if (defined USE_HDF5)
    USE HDF5
#endif
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE FindInVector(n,TF,pos)

    ! Inlet variables
    INTEGER,INTENT(IN):: n      ! Dimension of logical vector
    LOGICAL,INTENT(IN):: TF(n)  ! Logical vector (True or False)
    ! Outlet variables
    INTEGER npos                ! number of "true" conditions
    INTEGER,DIMENSION(:),ALLOCATABLE  :: pos! position of "true" conditions
    ! Internal variables
    INTEGER i                   ! counter
    INTEGER v(n)                ! vector of all positions
    !pos = 0                     ! Initialize pos
    FORALL(i=1:n)   v(i) = i    ! Enumerate all positions
    npos  = COUNT(TF)           ! Count the elements of TF that are .True.
!    PRINT*, 'npos',npos
    ALLOCATE(pos(npos))
    pos= pack(v, TF)    ! With Pack function, verify position of true conditions
!    PRINT *, 'pos', pos
!    PRINT *, 'TEST', pos(npos+1)
    END SUBROUTINE FindInVector

!    SUBROUTINE AVHRR_NDVI(filename,lon_ind,lat_ind,avg_NDVI)
    SUBROUTINE AVHRR_NDVI(ssid,lon_ind,lat_ind,avg_NDVI)
         CHARACTER (len=100)    :: filename
         INTEGER :: sd_id, ssid, DFACC, status, sfstart,sfselect,sfrdata, sfendacc, sfend, N_NDVI
         INTEGER :: start(2), edges(2), stride(2)
         INTEGER*2,DIMENSION(:,:),ALLOCATABLE :: data
!         INTEGER*8 ::sum_NDVI
         INTEGER,DIMENSION(:),ALLOCATABLE  :: lat_ind, lon_ind
         REAL*4 ::sum_NDVI, avg_NDVI

         start(1) = lon_ind(1)-1
         start(2) = lat_ind(1)-1
         edges(1) = size(lon_ind)
         edges(2) = size(lat_ind)
         stride(1) = 1
         stride(2) = 1
         ALLOCATE(data(size(lon_ind),size(lat_ind)))
         !DFACC=1 !Read Only Access
         !sd_id = sfstart(trim(filename),DFACC)
         !ssid = sfselect(sd_id, 0)

#if (defined USE_HDF5)
         status = sfrdata(ssid,start, stride,edges, data)
#endif
        !PRINT *, "Filename", filename
          sum_NDVI=SUM(REAL(data,4),MASK=data.NE.-9999)
          N_NDVI=COUNT(MASK=data.NE.-9999)
          avg_NDVI=REAL(sum_NDVI,4)/REAL(N_NDVI,4)*0.0001
!        PRINT*,data
!        PRINT*,'Sum count avg', sum_NDVI,N_NDVI,avg_NDVI
        !status = sfendacc(ssid)
        !status = sfend(sd_id)
    END SUBROUTINE AVHRR_NDVI
    
    SUBROUTINE AVHRR_NDVI_MATRIX(ssid,lon_ind,lat_ind, NDVI_MAT)
         INTEGER :: sd_id, ssid, DFACC, status, sfstart,sfselect,sfrdata, sfendacc, sfend, N_NDVI
         INTEGER :: start(2), edges(2), stride(2)
         INTEGER*2,DIMENSION(:,:),ALLOCATABLE :: data
         REAL*4,DIMENSION(:,:),ALLOCATABLE :: NDVI_MAT
         INTEGER,DIMENSION(:),ALLOCATABLE  :: lat_ind, lon_ind
         REAL*4 :: sum_NDVI, avg_NDVI
         
         start(1) = lon_ind(1)-1
         start(2) = lat_ind(1)-1
         edges(1) = size(lon_ind)
         edges(2) = size(lat_ind)
         stride(1) = 1
         stride(2) = 1
         ALLOCATE(data(size(lon_ind),size(lat_ind)))
         !DFACC=1 !Read Only Access
         !sd_id = sfstart(trim(filename),DFACC)
         !ssid = sfselect(sd_id, 0)
#if (defined USE_HDF5)
         status = sfrdata(ssid,start, stride,edges, data)
#endif
         NDVI_MAT=-9999
!        WHERE(data.NE.-9999)
        WHERE(data*0.0001.GT.0.0.AND.data*0.0001.LT.1.0)
        NDVI_MAT=data*0.0001
        END WHERE
        !PRINT *, "Filename", filename
!          sum_NDVI=SUM(data,MASK=data.NE.-9999)
!          N_NDVI=COUNT(MASK=data.NE.-9999)
!          avg_NDVI=REAL(sum_NDVI,4)/REAL(N_NDVI,4)*0.0001
        !PRINT*,'Sum count avg', sum_NDVI,N_NDVI,avg_NDVI
        !status = sfendacc(ssid)
        !status = sfend(sd_id)
    END SUBROUTINE AVHRR_NDVI_MATRIX

    SUBROUTINE GetSMAP(filename,dataset,smap_row,smap_col,lon_ind,lat_ind,sm_mat)
    CHARACTER (len=100)    :: filename, dataset
#if (defined USE_HDF5)
    INTEGER(HID_T) :: file_id, dataset_id
#endif
    INTEGER        :: hdferr, smap_row, smap_col
    INTEGER,DIMENSION(:),ALLOCATABLE :: lon_ind, lat_ind
    REAL*8,DIMENSION(:,:),ALLOCATABLE  :: sm_mat, data_out
#if (defined USE_HDF5)
    INTEGER(HSIZE_T),DIMENSION(2):: data_dims
#endif
#if (defined USE_HDF5)
       CALL h5open_f(hdferr) !Initialize hdf5
       CALL h5fopen_f (trim(filename),H5F_ACC_RDONLY_F,file_id,hdferr) !Open file
!        PRINT*, 'sd_id, hdferr', file_id, hdferr
       CALL h5dopen_f (file_id, trim(dataset),dataset_id, hdferr) !Open dataset 
!       PRINT*,'ds_id hdferr', dataset_id, hdferr

        ALLOCATE(data_out(smap_row,smap_col),sm_mat(size(lon_ind),size(lat_ind)))
        data_dims(1) = smap_row
        data_dims(2) = smap_col
        CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, hdferr) !Read dataset
        CALL h5dclose_f(dataset_id,hdferr)
        CALL h5fclose_f(file_id,hdferr)
        CALL h5close_f(hdferr)
        sm_mat = data_out(lon_ind(1):lon_ind(size(lon_ind)),lat_ind(1):lat_ind(size(lat_ind)))        
       ! PRINT*, 'SMAP Data', sm_mat
#endif
    END SUBROUTINE GetSMAP

    SUBROUTINE GetSMAP_L2(filename,dataset1,dataset2,dataset3,data_out,ind)
    CHARACTER (len=100)    :: filename, dataset1, dataset2, dataset3
#if (defined USE_HDF5)
    INTEGER(HID_T) :: file_id, dataset_id1, dataset_id2, dataset_id3, dspace_id
#endif
    INTEGER        :: hdferr, m
    REAL*8,DIMENSION(:),ALLOCATABLE  :: data_out
    INTEGER*4,DIMENSION(:),ALLOCATABLE  :: row_temp_ind, col_temp_ind, ind
#if (defined USE_HDF5)
    INTEGER(HSIZE_T),DIMENSION(1):: dims, maxdims
#endif
#if (defined USE_HDF5)
       CALL h5open_f(hdferr) !Initialize hdf5
       CALL h5fopen_f (trim(filename),H5F_ACC_RDONLY_F,file_id,hdferr) !Open file
!        PRINT*, 'sd_id, hdferr', file_id, hdferr
       CALL h5dopen_f (file_id, trim(dataset1),dataset_id1, hdferr) !Open dataset
       call h5dget_space_f(dataset_id1,dspace_id,hdferr)
       call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr) 
       CALL h5dopen_f (file_id, trim(dataset2),dataset_id2, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset3),dataset_id3, hdferr) !Open dataset 
!       PRINT*,'ds_id hdferr', dataset_id1, hdferr
!       PRINT*,'dims maxdims', dims, maxdims
        m=dims(1) 
        ALLOCATE(data_out(m),row_temp_ind(m),col_temp_ind(m), ind(m))
!        print*, 'm dims', m, dims
        CALL h5dread_f(dataset_id1, H5T_NATIVE_DOUBLE, data_out, dims, hdferr) !Read dataset
        CALL h5dread_f(dataset_id2, H5T_NATIVE_INTEGER, col_temp_ind, dims, hdferr) !Read dataset
        CALL h5dread_f(dataset_id3, H5T_NATIVE_INTEGER, row_temp_ind, dims, hdferr) !Read dataset
        ind=row_temp_ind*3856+col_temp_ind+1 !2D Matrix index to 1D vector index
!        print *, size(data_out), data_out(192571), data_out(205664)
!        print *, row_temp_ind(192571), row_temp_ind(205664)
!        print *, col_temp_ind(192571), col_temp_ind(205664)
        CALL h5dclose_f(dataset_id1,hdferr)
!        CALL h5dclose_f(dataset_id2,hdferr)
!        CALL h5dclose_f(dataset_id3,hdferr)
        CALL h5fclose_f(file_id,hdferr)
        CALL h5close_f(hdferr)
       ! PRINT*, 'SMAP Data', sm_mat
#endif
    END SUBROUTINE GetSMAP_L2

    SUBROUTINE GetSMAP_L1B(filename, data1_out, data2_out, data3_out, data4_out, data5_out, data6_out, data7_out, data8_out, &
                           data9_out, data10_out, data11_out, data12_out, data13_out, data14_out, data15_out, n, m)

      use LDT_logMod, only: LDT_logunit ! EMK

    CHARACTER (len=100)    :: filename, dataset1, dataset2, dataset3, dataset4, dataset5, dataset6, dataset7
    CHARACTER (len=100)    :: dataset8, dataset9, dataset10, dataset11, dataset12, dataset13, dataset14, dataset15
#if (defined USE_HDF5)
    INTEGER(HID_T) :: file_id, dataset_id1, dataset_id2, dataset_id3, dataset_id4, dataset_id5, dataset_id6, dataset_id7, dspace_id
    INTEGER(HID_T) :: dataset_id8, dataset_id9, dataset_id10, dataset_id11, dataset_id12, dataset_id13, dataset_id14, dataset_id15
#endif
    INTEGER        :: hdferr, m, n
    REAL*4,DIMENSION(:,:),ALLOCATABLE  :: data1_out, data2_out, data3_out, data4_out, data5_out, data6_out
    REAL*4,DIMENSION(:,:),ALLOCATABLE  :: data7_out, data8_out, data9_out, data10_out, data11_out, data15_out
    REAL*4,DIMENSION(:),ALLOCATABLE  :: data14_out
    INTEGER*4,DIMENSION(:,:),ALLOCATABLE  :: data12_out, data13_out
#if (defined USE_HDF5)
    INTEGER(HSIZE_T),DIMENSION(2):: dims, maxdims
#endif
       !DEFINED DATA TYPE
       dataset1="/Brightness_Temperature/tb_time_seconds/"
       dataset2="/Brightness_Temperature/tb_v_surface_corrected/"          !Dataset name in smap
       dataset3="/Brightness_Temperature/tb_h_surface_corrected/"          !Dataset name in smap
       dataset4="/Brightness_Temperature/tb_v/"                                !Dataset name in smap
       dataset5="/Brightness_Temperature/tb_h/"                                !Dataset name in smap
       dataset6="/Brightness_Temperature/surface_water_fraction_mb_v/"    !Dataset name in smap
       dataset7="/Brightness_Temperature/surface_water_fraction_mb_h/"    !Dataset name in smap
       dataset8="/Brightness_Temperature/nedt_v/"    !Dataset name in smap
       dataset9="/Brightness_Temperature/nedt_h/"    !Dataset name in smap

       dataset10="/Brightness_Temperature/tb_lat/"
       dataset11="/Brightness_Temperature/tb_lon/"
       dataset12="/Brightness_Temperature/tb_qual_flag_v/"
       dataset13="/Brightness_Temperature/tb_qual_flag_h/"
       dataset14="/Spacecraft_Data/sc_nadir_angle/"
       dataset15="/Brightness_Temperature/antenna_scan_angle/"

#if (defined USE_HDF5)
       CALL h5open_f(hdferr) !Initialize hdf5
       CALL h5fopen_f (trim(filename),H5F_ACC_RDONLY_F,file_id,hdferr) !Open file
       !        PRINT*, 'sd_id, hdferr', file_id, hdferr
       CALL h5dopen_f (file_id, trim(dataset1),dataset_id1, hdferr) !Open dataset

       call h5dget_space_f(dataset_id1,dspace_id,hdferr)
       call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
       CALL h5dopen_f (file_id, trim(dataset2),dataset_id2, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset3),dataset_id3, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset4),dataset_id4, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset5),dataset_id5, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset6),dataset_id6, hdferr) !Open dataset       
       CALL h5dopen_f (file_id, trim(dataset7),dataset_id7, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset8),dataset_id8, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset9),dataset_id9, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset10),dataset_id10, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset11),dataset_id11, hdferr) !Open dataset        
       CALL h5dopen_f (file_id, trim(dataset12),dataset_id12, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset13),dataset_id13, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset14),dataset_id14, hdferr) !Open dataset 
       CALL h5dopen_f (file_id, trim(dataset15),dataset_id15, hdferr) !Open dataset       

       !PRINT*,'ds_id hdferr', dataset_id1, hdferr
       !PRINT*,'dims maxdims', dims, maxdims
        n=dims(1)
        m=dims(2)
        ALLOCATE(data1_out(n,m),data2_out(n,m),data3_out(n,m),data4_out(n,m),data5_out(n,m),data6_out(n,m),data7_out(n,m),data8_out(n,m))
        ALLOCATE(data9_out(n,m),data10_out(n,m),data11_out(n,m),data12_out(n,m),data13_out(n,m),data14_out(m),data15_out(n,m))
  
!        print*, 'm dims', m, dims
        CALL h5dread_f(dataset_id1, H5T_IEEE_F32LE, data1_out, dims, hdferr) !Read dataset
        CALL h5dread_f(dataset_id2, H5T_IEEE_F32LE, data2_out, dims, hdferr) !Read dataset
        CALL h5dread_f(dataset_id3, H5T_IEEE_F32LE, data3_out, dims, hdferr) !Read datase
        CALL h5dread_f(dataset_id4, H5T_IEEE_F32LE, data4_out, dims, hdferr)
        CALL h5dread_f(dataset_id5, H5T_IEEE_F32LE, data5_out, dims, hdferr) !Read datase
        CALL h5dread_f(dataset_id6, H5T_IEEE_F32LE, data6_out, dims, hdferr) !Read datase
        CALL h5dread_f(dataset_id7, H5T_IEEE_F32LE, data7_out, dims, hdferr) !Read dataset
        CALL h5dread_f(dataset_id8, H5T_IEEE_F32LE, data8_out, dims, hdferr) !Read datase
        CALL h5dread_f(dataset_id9, H5T_IEEE_F32LE, data9_out, dims, hdferr)
        CALL h5dread_f(dataset_id10, H5T_IEEE_F32LE, data10_out, dims, hdferr) !Read datase
        CALL h5dread_f(dataset_id11, H5T_IEEE_F32LE, data11_out, dims, hdferr) !Read datase
        CALL h5dread_f(dataset_id12, H5T_NATIVE_INTEGER, data12_out, dims, hdferr)
        CALL h5dread_f(dataset_id13, H5T_NATIVE_INTEGER, data13_out, dims, hdferr)
        CALL h5dread_f(dataset_id14, H5T_IEEE_F32LE, data14_out, dims, hdferr) !Read datase
        CALL h5dread_f(dataset_id15, H5T_IEEE_F32LE, data15_out, dims, hdferr) !Read datase

        !CALL h5dread_f(dataset_id4, H5T_STD_U16LE, data4_out, dims, hdferr) !Read dataset
!        print *, size(data_out), data_out(192571), data_out(205664)
!        print *, row_temp_ind(192571), row_temp_ind(205664)
!        print *, col_temp_ind(192571), col_temp_ind(205664)
        CALL h5dclose_f(dataset_id1,hdferr)
        CALL h5dclose_f(dataset_id2,hdferr)
        CALL h5dclose_f(dataset_id3,hdferr)
        CALL h5dclose_f(dataset_id4,hdferr)
        CALL h5dclose_f(dataset_id5,hdferr)
        CALL h5dclose_f(dataset_id6,hdferr)
        CALL h5dclose_f(dataset_id7,hdferr)
        CALL h5dclose_f(dataset_id8,hdferr)
        CALL h5dclose_f(dataset_id9,hdferr)
        CALL h5dclose_f(dataset_id10,hdferr)
        CALL h5dclose_f(dataset_id11,hdferr)
        CALL h5dclose_f(dataset_id12,hdferr)
        CALL h5dclose_f(dataset_id13,hdferr)
        CALL h5dclose_f(dataset_id14,hdferr)
        CALL h5dclose_f(dataset_id15,hdferr)
        CALL h5fclose_f(file_id,hdferr)
        CALL h5close_f(hdferr)
       ! PRINT*, 'SMAP Data', sm_mat
#endif
    END SUBROUTINE GetSMAP_L1B

      ! Forked version of GetSMAP_L1B to SMAP files in operations, processing
      ! a subset of the fields.  Also with fault tolerance and some special
      ! logic for older NRT files missing the tb_v_surface_corrected field.
      ! Eric Kemp/SSAI.
      SUBROUTINE GetSMAP_L1B_NRT_subset(filename, tb_time_seconds, &
           tb_v_surface_corrected, &
           tb_lat, tb_lon, tb_qual_flag_v, tb_qual_flag_h, sc_nadir_angle, &
           antenna_scan_angle, n, m, ierr)

      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun ! EMK

      ! Arguments
      character(*), intent(in) :: filename
      real*4, allocatable, intent(out) :: tb_time_seconds(:,:)
      real*4, allocatable, intent(out) :: tb_v_surface_corrected(:,:)
      real*4, allocatable, intent(out) :: tb_lat(:,:), tb_lon(:,:)
      integer*4, allocatable, intent(out) :: tb_qual_flag_v(:,:)
      integer*4, allocatable, intent(out) :: tb_qual_flag_h(:,:)
      real*4, allocatable, intent(out) :: sc_nadir_angle(:)
      real*4, allocatable, intent(out) :: antenna_scan_angle(:,:)
      integer, intent(out) :: m, n
      integer, intent(out) :: ierr

#if (defined USE_HDF5)

      ! Locals
      character(100) :: dataset
      integer(HID_T) :: file_id, dataset_id, dspace_id
      integer(HSIZE_T) :: dims(2), maxdims(2)
      integer :: hdferr
      logical :: exists, ishdf5

      ierr = 0
      m = 0
      n = 0

      ! Make sure file exists
      inquire(file=trim(filename), exist=exists)
      if (.not. exists) then
         write(LDT_logunit,*)'[ERR] Cannot find file ', trim(filename)
         ierr = 1
         return
      end if

      ! Initialize HDF5
      call h5open_f(hdferr)
      if (hdferr == -1) then
         write(LDT_logunit,*)'[ERR] Cannot initialize HDF5 Fortran interface!'
         call h5close_f(hdferr)
         ierr = 1
         return
      end if

      ! Make sure the file is HDF5
      call h5fis_hdf5_f(trim(filename), ishdf5, hdferr)
      if (hdferr == -1) then
         write(LDT_logunit,*)'[ERR] Problem checking if ', trim(filename), &
              ' is HDF5'
         call h5close_f(hdferr)
         ierr = 1
         return
      end if
      if (.not. ishdf5) then
         write(LDT_logunit,*)'[ERR] File ', trim(filename), ' is not HDF5!'
         call h5close_f(hdferr)
         ierr = 1
         return
      end if

      ! Open the file
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
      if (hdferr == -1) then
         write(LDT_logunit,*)'[ERR] Cannot open ', trim(filename)
         call h5close_f(hdferr)
         ierr = 1
         return
      end if

      ! Find Tb_lat, plus dimensions
      dataset = "/Brightness_Temperature/tb_lat/"
      call get_dataset_id(file_id, dataset, dataset_id, ierr)
      if (ierr == 1) return
      call h5dget_space_f(dataset_id, dspace_id, hdferr)
      if (hdferr == -1) then
         write(LDT_logunit,*)'[ERR] Cannot find dimensions for ', &
              trim(dataset)
         call h5dclose_f(dataset_id, hdferr)
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         ierr = 1
         return
      end if
      call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
      if (hdferr == -1) then
         write(LDT_logunit,*)'[ERR] Cannot get dimensions for ', &
              trim(dataset)
         call h5dclose_f(dataset_id, hdferr)
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         ierr = 1
         return
      end if

      ! We have the dimensions for the arrays, so let's allocate here.
      n = dims(1)
      m = dims(2)
      allocate(tb_v_surface_corrected(n,m)); tb_v_surface_corrected = 0
      allocate(tb_time_seconds(n,m)); tb_time_seconds = 0
      allocate(tb_lat(n,m)) ; tb_lat = 0
      allocate(tb_lon(n,m)) ; tb_lon = 0
      allocate(tb_qual_flag_v(n,m)) ; tb_qual_flag_v = 0
      allocate(tb_qual_flag_h(n,m)) ; tb_qual_flag_h = 0
      allocate(sc_nadir_angle(n)) ; sc_nadir_angle = 0
      allocate(antenna_scan_angle(n,m)) ; antenna_scan_angle = 0

      ! Get Tb_lat
      call read_dataset_real(file_id, dataset, dataset_id, tb_lat, &
           dims, ierr)
      if (ierr == 1) return

      ! Get tb_lon
      dataset = "/Brightness_Temperature/tb_lon/"
      call get_dataset_id(file_id, dataset, dataset_id, ierr)
      if (ierr == 1) return
      call read_dataset_real(file_id, dataset, dataset_id, tb_lon, &
           dims, ierr)
      if (ierr == 1) return

      ! Get tb_time_seconds
      dataset = "/Brightness_Temperature/tb_time_seconds/"
      call get_dataset_id(file_id, dataset, dataset_id, ierr)
      if (ierr == 1) return
      call read_dataset_real(file_id, dataset, dataset_id, tb_time_seconds, &
           dims, ierr)
      if (ierr == 1) return

      ! Get tb_qual_flag_v
      dataset = "/Brightness_Temperature/tb_qual_flag_v/"
      call get_dataset_id(file_id, dataset, dataset_id, ierr)
      if (ierr == 1) return
      call read_dataset_integer(file_id, dataset, dataset_id, tb_qual_flag_v, &
           dims, ierr)
      if (ierr == 1) return

      ! Get tb_qual_flag_h
      dataset = "/Brightness_Temperature/tb_qual_flag_h/"
      call get_dataset_id(file_id, dataset, dataset_id, ierr)
      if (ierr == 1) return
      call read_dataset_integer(file_id, dataset, dataset_id, tb_qual_flag_h, &
           dims, ierr)
      if (ierr == 1) return

      ! Get sc_nadir_angle
      dataset = "/Spacecraft_Data/sc_nadir_angle/"
      call get_dataset_id(file_id, dataset, dataset_id, ierr)
      if (ierr == 1) return
      call read_dataset_real_1d(file_id, dataset, dataset_id, sc_nadir_angle, &
           dims, ierr)
      if (ierr == 1) return

      ! Get antenna_scan_angle
      dataset = "/Brightness_Temperature/antenna_scan_angle/"
      call get_dataset_id(file_id, dataset, dataset_id, ierr)
      if (ierr == 1) return
      call read_dataset_real(file_id, dataset, dataset_id, &
           antenna_scan_angle, dims, ierr)
      if (ierr == 1) return

      ! We will try to get tb_v_surface_corrected, but this is missing
      ! in older NRT files.  If absent, we will substitute tb_v.
      ! NOTE:  If a *read* error occurs, the file will assumed to be
      ! corrupt all arrays will be nuked.
      ! Get tb_v_surface_corrected
      dataset = "/Brightness_Temperature/tb_v_surface_corrected/"
      call get_dataset_id(file_id, dataset, dataset_id, ierr, &
           handle_missing_nrt_field=.true.)
      if (ierr == 1) then
         ierr = 0
         dataset = "/Brightness_Temperature/tb_v/"
         write(LDT_logunit,*)'[WARN] Will try substituting ', trim(dataset)
         call get_dataset_id(file_id, dataset, dataset_id, ierr)
         if (ierr == 1) return
      end if
      call read_dataset_real(file_id, dataset, dataset_id, &
           tb_v_surface_corrected, &
           dims, ierr)
      if (ierr == 1) return

      ! Clean up
      call h5fclose_f(file_id, hdferr)
      call h5close_f(hdferr)
      ierr = 0

      return

    contains

      ! Internal subroutine
      subroutine get_dataset_id(file_id, dataset, dataset_id, ierr, &
           handle_missing_nrt_field)
          implicit none
          integer(HID_T), intent(in) :: file_id
          character(*), intent(in) :: dataset
          integer(HID_T), intent(out) :: dataset_id
          integer, intent(out) :: ierr
          logical, optional, intent(in) :: handle_missing_nrt_field
          logical :: link_exists
          integer :: hdferr
          call h5lexists_f(file_id, trim(dataset), link_exists, hdferr)
          if (hdferr == -1) then
             write(LDT_logunit,*)'[ERR] Problem finding ', trim(dataset)
             ierr = 1
             if (handle_missing_nrt_field) return
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             call freeall(ierr)
             return
          endif
          if (.not. link_exists) then
             write(LDT_logunit,*)'[ERR] Nonexistent dataset ', trim(dataset)
             ierr = 1
             if (handle_missing_nrt_field) return
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             call freeall(ierr)
             return
          endif
          call h5dopen_f(file_id, trim(dataset), dataset_id, hdferr)
          if (hdferr == -1) then
             write(LDT_logunit,*)'[ERR] Cannot open dataset ', trim(dataset)
             ierr = 1
             if (handle_missing_nrt_field) return
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             call freeall(ierr)
             return
          end if
        end subroutine get_dataset_id

        ! Internal subroutine
        subroutine read_dataset_real(file_id, dataset, dataset_id, buf, &
             dims, ierr)
          implicit none
          integer(HID_T), intent(in) :: file_id
          character(*), intent(in) :: dataset
          integer(HID_T), intent(in) :: dataset_id
          real*4, intent(inout) :: buf(:,:)
          integer(HSIZE_T), intent(in) :: dims(:)
          integer, intent(out) :: ierr
          integer :: hdferr
          call h5dread_f(dataset_id, H5T_IEEE_F32LE, buf, dims, hdferr)
          if (hdferr == -1) then
             write(LDT_logunit,*)'[ERR] Cannot read dataset ', trim(dataset)
             call h5dclose_f(dataset_id, hdferr)
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             call freeall(ierr)
             return
          endif
          call h5dclose_f(dataset_id, hdferr)
          if (hdferr == -1) then
             write(LDT_Logunit,*)'[ERR] Problem closing dataset ', &
                  trim(dataset)
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             call freeall(ierr)
             return
          end if
        end subroutine read_dataset_real

        ! Internal subroutine
        subroutine read_dataset_real_1d(file_id, dataset, dataset_id, buf, &
             dims, ierr)
          implicit none
          integer(HID_T), intent(in) :: file_id
          character(*), intent(in) :: dataset
          integer(HID_T), intent(in) :: dataset_id
          real*4, intent(inout) :: buf(:)
          integer(HSIZE_T), intent(in) :: dims(:)
          integer, intent(out) :: ierr
          integer :: hdferr
          call h5dread_f(dataset_id, H5T_IEEE_F32LE, buf, dims, hdferr)
          if (hdferr == -1) then
             write(LDT_logunit,*)'[ERR] Cannot read dataset ', trim(dataset)
             call h5dclose_f(dataset_id, hdferr)
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             call freeall(ierr)
             return
          endif
          call h5dclose_f(dataset_id, hdferr)
          if (hdferr == -1) then
             write(LDT_Logunit,*)'[ERR] Problem closing dataset ', &
                  trim(dataset)
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             call freeall(ierr)
             return
          end if
        end subroutine read_dataset_real_1d

        ! Internal subroutine
        subroutine read_dataset_integer(file_id, dataset, dataset_id, buf, &
             dims, ierr)
          implicit none
          integer(HID_T), intent(in) :: file_id
          character(*), intent(in) :: dataset
          integer(HID_T), intent(in) :: dataset_id
          integer*4, intent(inout) :: buf(:,:)
          integer(HSIZE_T), intent(in) :: dims(:)
          integer, intent(out) :: ierr
          integer :: hdferr
          call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, buf, dims, hdferr)
          if (hdferr == -1) then
             write(LDT_logunit,*)'[ERR] Cannot read dataset ', trim(dataset)
             call h5dclose_f(dataset_id, hdferr)
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             call freeall(ierr)
             return
          endif
          call h5dclose_f(dataset_id, hdferr)
          if (hdferr == -1) then
             write(LDT_Logunit,*)'[ERR] Problem closing dataset ', &
                  trim(dataset)
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             call freeall(ierr)
             return
          end if
        end subroutine read_dataset_integer

        ! Internal subroutine.  Warning -- deallocates memory in
        ! parent subroutine and resets two variables.  This is intended
        ! for gracefully handling errors returned from HDF5.
        subroutine freeall(ierr)
          implicit none
          integer, intent(out) :: ierr
          if (allocated(tb_v_surface_corrected)) &
               deallocate(tb_v_surface_corrected)
          if (allocated(tb_lat)) deallocate(tb_lat)
          if (allocated(tb_lon)) deallocate(tb_lon)
          if (allocated(tb_time_seconds)) deallocate(tb_time_seconds)
          if (allocated(tb_qual_flag_v)) deallocate(tb_qual_flag_v)
          if (allocated(tb_qual_flag_h)) deallocate(tb_qual_flag_h)
          if (allocated(sc_nadir_angle))deallocate(sc_nadir_angle)
          if (allocated(antenna_scan_angle)) deallocate(antenna_scan_angle)
          m = 0
          n = 0
          ierr = 1
          return
        end subroutine freeall
#else
        ! Dummy version if LDT was compiled w/o HDF5 support.
        write(LDT_logunit,*) &
             '[ERR] GetSMAP_L1B_NRT called without HDF5 support!'
        write(LDT_logunit,*) &
             '[ERR] Recompile LDT with HDF5 support and try again!'
        call LDE_endrun()
#endif
      end subroutine GetSMAP_L1B_NRT_Subset

    SUBROUTINE NEAREST_1d ( nd, xd, yd, ni, xi, yi )
     IMPLICIT NONE
     INTEGER ::  i, j, k, nd, ni
     REAL*8  :: d, d2, xd(nd), xi(ni)
     REAL*4  ::  yd(nd), yi(ni)

      DO i = 1, ni
         k = 1
         d = ABS( xi(i) - xd(k))
         DO j = 2, nd
            d2 = ABS( xi(i) - xd(j) )
            IF ( d2 < d ) THEN
               k = j
               d = d2
            END IF
         END DO

        yi(i) = yd(k)
     END DO
   END SUBROUTINE NEAREST_1d

   SUBROUTINE IDW_2D (nd, xd, yd, zd, p, xi, yi, zi)
!*****************************************************************************
!
!! IDW_2D evaluates a 2D Inverse Distance Interpolant.
!  Discussion: This code should be vectorized.
!  Modified: 24 Jan 2019 by P.-W. Liu
!  Reference:
!
!  Parameters:
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    Input, real ( kind = 8 ) XD(ND), YD(ND), the data points.
!    Input, real ( kind = 8 ) ZD(ND), the data values.
!    Input, real ( kind = 8 ) P, the power.
!    Input, real ( kind = 8 ) XI, YI the interpolation points.
!    Output, real ( kind = 8 ) ZI, the interpolated values.
!*****************************************************************************
     integer ( kind = 4 ) nd, j, z
     real ( kind = 4 ) p, s
     real ( kind = 8 ) w(nd)
     real ( kind = 4 ) xd(nd)
     real ( kind = 4 ) xi
     real ( kind = 4 ) yd(nd)
     real ( kind = 4 ) yi
     real ( kind = 4 ) zd(nd)
     real ( kind = 4 ) zi


    if ( p == 0.0D+00 ) then
      w(1:nd) = 1.0D+00 / real ( nd, kind = 8 )
    else

      z = -1
      do j = 1, nd
           w(j) = sqrt ( ( xi - xd(j) ) ** 2 + ( yi - yd(j) ) ** 2 )
        if ( w(j) == 0.0D+00 ) then
          z = j
          exit
        end if
      end do

      if ( z /= -1 ) then
        w(1:nd) = 0.0D+00
        w(z) = 1.0D+00
      else
        w(1:nd) = 1.0D+00 / w(1:nd) ** p
!        PRINT*,'w2',w
        s = sum ( w, MASK=zd.GT.-9999 )
!        PRINT*,'s',s
        WHERE(zd.LE.-9998)
        w = w*0
        END WHERE
        w(1:nd) = w(1:nd) / s
!        PRINT*,'w3',w, 'SUM_W',sum(w)

      end if
    end if

       IF (s.LE.0.00000001) THEN
           zi=-9999
       ELSE
           zi = dot_product ( w, zd )
       END IF
!  return
END SUBROUTINE IDW_2D


   SUBROUTINE LSQ_Linear(n, x, y, a, b, r, rmse)
   !a:slope, b: intercep, r: squared correlation coefficient
   IMPLICIT NONE
   INTEGER          ::  n            !!Number of data point
   REAL*4           ::  x(n), y(n)   !!Input data vectors of x and y
   REAL*8           ::  a, b, r, sumx, sumx2, sumxy, sumy, sumy2
   REAL*8 ::rmse
      x=REAL(x,8)
      y=REAL(y,8)
!      PRINT *,'n=',n,'x and Y', x, y
      sumx  = SUM(x)                              ! compute sum of x
      sumx2 = DOT_PRODUCT(x, x)                   ! compute sum of x**2
      sumxy = DOT_PRODUCT(x, y)                   ! compute sum of x * y
      sumy  = SUM(y)                              ! compute sum of y
      sumy2 = DOT_PRODUCT(y, y)                   ! compute sum of y**2
      
   a = (REAL(n,8) * sumxy  -  sumx * sumy) / (REAL(n,8) * sumx2 - sumx**2)     ! compute slope
   b = (sumy * sumx2  -  sumx * sumxy) / (REAL(n,8) * sumx2  -  sumx**2)       ! compute y-intercept
   r = ((sumxy - sumx * sumy / REAL(n,8)) / SQRT((sumx2 - sumx**2/REAL(n,8)) * (sumy2 - sumy**2/REAL(n,8))))**2 !Compute r square
   rmse =SQRT(DOT_PRODUCT((a*x+b-y),(a*x+b-y))/REAL(n-2,8))
   !PRINT*, 'RMSE', rmse
   !PRINT*, " Slope        a = ", a                        ! print results
   !PRINT*, " y-intercept  b = ", b
   !PRINT*, " Correlation  r = ", r
   END SUBROUTINE LSQ_Linear



END MODULE TOOLSUBS
