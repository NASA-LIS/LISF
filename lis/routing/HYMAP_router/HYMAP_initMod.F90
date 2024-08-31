!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module HYMAP_initMod
  
  implicit none
  
contains
  
  subroutine adjust_nextxy(nx,ny,imis,i2nextx,i2nexty,i2mask,zgx,zgy,zpx,zpy,zres)
    ! ================================================
    ! to adjust flow direction matrixes for smaller domain
    ! by Augusto GETIRANA
    ! on 13th Mar 2012
    ! at HSL/GSFC/NASA
    ! ================================================   
    
    implicit none       
  
    integer, intent(in)    :: nx                  ! number of grids in horizontal
    integer, intent(in)    :: ny                  ! number of grids in vertical
    integer, intent(in)    :: imis                ! integer undefined value
    integer, intent(inout) :: i2nextx(nx,ny)      ! point downstream horizontal
    integer, intent(inout) :: i2nexty(nx,ny)      ! point downstream horizontal
    integer, intent(in)    :: i2mask(nx,ny)       ! mask limiting modeled region (0: out; >=1: in)
    real*4,  intent(in)    :: zgx                 ! lon min of smaller domain
    real*4,  intent(in)    :: zgy                 ! lat min of smaller domain
    real*4,  intent(in)    :: zpx                 ! lon min of larger domain
    real*4,  intent(in)    :: zpy                 ! lat min of larger domain
    real*4,  intent(in)    :: zres                ! spatial resolution
    
    integer, parameter :: ibound = -9
    
    integer             ::  idx,idy
    
    idx=int((zgx-zpx)/zres)
    idy=int((zgy-zpy)/zres)
    
    where(i2nextx>0)i2nextx=i2nextx-idx
    where(i2nexty>0)i2nexty=i2nexty-idy
    
    i2nexty(1,:)=imis
    i2nexty(nx,:)=imis
    i2nexty(:,1)=imis
    i2nexty(:,ny)=imis
    
    i2nextx(1,:)=imis
    i2nextx(nx,:)=imis
    i2nextx(:,1)=imis
    i2nextx(:,ny)=imis
    
    where(i2nextx<=1.and.i2nextx/=imis.and.i2mask>0)
       i2nextx=ibound
       i2nexty=ibound
    endwhere
    
    where(i2nextx>=nx.and.i2mask>0)
       i2nextx=ibound
       i2nexty=ibound
    endwhere
    
    where(i2nexty<=1.and.i2nexty/=imis.and.i2mask>0)
       i2nextx=ibound
       i2nexty=ibound
    endwhere
    
    where(i2nexty>=ny.and.i2mask>0)
       i2nextx=ibound
       i2nexty=ibound
    endwhere
         
  end subroutine adjust_nextxy
  ! ================================================          
  ! ================================================            
  subroutine calc_1d_seq(nx,ny,imis,i2nextx,i2mask,i1seqx,i1seqy,nseqriv,nseqall)
    ! ================================================
    ! to calculate 1d sequnece file
    ! by Dai YAMAZAKI
    ! on 24th Mar 2008
    ! at IIS, UT
    ! ================================================
    use LIS_logMod,     only : LIS_logunit
    
    implicit none
    ! Input
    integer             ::  nx                  !! number of grids in horizontal
    integer             ::  ny                  !! number of grids in vertical
    integer             ::  imis                !! integer undefined value
    integer             ::  i2nextx(nx,ny)      !! point downstream horizontal
    integer             ::  i2mask(nx,ny)       !! mask limiting modeled region (0: out; >=1: in)
    ! Output
    integer             ::  i1seqx(nx*ny)       !! 1D sequence horizontal
    integer             ::  i1seqy(nx*ny)       !! 1D sequence vertical
    integer             ::  nseqriv             !! length of 1D sequnece for river
    integer             ::  nseqall             !! length of 1D sequnece for river and mouth
    ! Local
    integer             ::  ix
    integer             ::  iy
    integer             ::  iseq
    integer             ::  i
    ! ================================================
    iseq=0
    do iy=1, ny
       do ix=1, nx
          if( i2nextx(ix,iy).gt.0 .and.i2mask(ix,iy)>0 )then
             iseq=iseq+1
             i1seqx(iseq)=ix
             i1seqy(iseq)=iy
          endif
       end do
    end do  
    nseqriv=iseq   
    do iy=1, ny
       do ix=1, nx
          if( i2nextx(ix,iy).lt.0 .and. i2nextx(ix,iy).ne.imis .and.i2mask(ix,iy)>0 )then
             iseq=iseq+1
             i1seqx(iseq)=ix
             i1seqy(iseq)=iy
          endif
       end do
    end do
    nseqall=iseq    
    write(LIS_logunit,*)'[calc_1d_seq] number of cells',nseqall   
    return
  end subroutine calc_1d_seq
  ! ================================================          
  ! ================================================              
  subroutine set_fldstg(nx,ny,nz,i1seqx,i1seqy,nseqall,&
    r2fldhgt,r2areamat,r2rivlen,   &
    r2rivwth,r2rivstomax,          &
    r2fldstomax,r2fldgrd,r2rivare  )
    ! ================================================
    ! to   set floodplain staging
    ! by   Augusto Getirana after Dai YAMAZAKI
    ! on   27 Oct 2011
    ! at   NASA,GSFC
    ! Adapted for flow routing implementation in LIS 09 Nov 2011
    ! ================================================
    use LIS_logMod,     only : LIS_logunit
    
    implicit none
    ! index
    integer             ::  nx                  !! number of grids in horizontal
    integer             ::  ny                  !! number of grids in vertical
    integer             ::  nz                  !! number of stages in the sub-grid discretization
    integer             ::  i1seqx(nx*ny)       !! 1D sequence horizontal
    integer             ::  i1seqy(nx*ny)       !! 1D sequence vertical
    integer             ::  nseqall             !! length of 1D sequnece for river and mouth
    ! boundary
    real    ::  r2areamat(nx,ny)    !! area of the grid [m^2]
    real    ::  r2rivlen(nx,ny)     !! channel length [m]
    real    ::  r2rivwth(nx,ny)     !! river width [m]
    ! input
    real    ::  r2fldhgt(nx,ny,nz)  !! floodplain height
    real    ::  r2rivstomax(nx,ny)  !! maximum river storage [m3]
    ! output
    real    ::  r2fldstomax(nx,ny,nz)  !! maximum floodplain storage [m3]
    real    ::  r2fldgrd(nx,ny,nz)  !! floodplain gradient [-]
    ! local
    integer             ::  ix,iy, iseq, i
    real    ::  rstonow
    real    ::  rstopre
    real    ::  rhgtpre
    real    ::  rwthinc
    integer ::  i1,i2

    real    ::  r2rivare(nx,ny)
    real    ::  rrivfrc          !!fraction of river surface area [-]
    real    ::  r2hgtprv(nx,ny)
    ! ================================================
    do iseq=1, nseqall
       ix=i1seqx(iseq)
       iy=i1seqy(iseq)
       rrivfrc=r2rivare(ix,iy)/r2areamat(ix,iy)
       i1=max(1,int(rrivfrc*nz))
       i2=min(nz,i1+1)
       r2hgtprv(ix,iy)=r2fldhgt(ix,iy,i1)+(r2fldhgt(ix,iy,i2)-r2fldhgt(ix,iy,i1))*(rrivfrc*nz-i1)
       r2fldhgt(ix,iy,:)=r2fldhgt(ix,iy,:)-r2hgtprv(ix,iy)
       where(r2fldhgt(ix,iy,:)<0)r2fldhgt(ix,iy,:)=0.
    enddo
    r2fldstomax = 0.
    r2fldgrd    = 0.   
    do iseq=1, nseqall
        ix=i1seqx(iseq)
        iy=i1seqy(iseq)
        rstopre = r2rivstomax(ix,iy)
        rhgtpre = 0.
        rwthinc = r2areamat(ix,iy) / r2rivlen(ix,iy) * 0.1
        do i=1, nz
           rstonow = r2rivlen(ix,iy) * ( r2rivwth(ix,iy) + rwthinc*(real(i)-0.5) ) * (r2fldhgt(ix,iy,i)-rhgtpre)
           r2fldstomax(ix,iy,i) = rstopre + rstonow
           r2fldgrd(ix,iy,i) = (r2fldhgt(ix,iy,i)-rhgtpre) / rwthinc
           rstopre = r2fldstomax(ix,iy,i)
           rhgtpre = r2fldhgt(ix,iy,i)
        end do
     end do

   return
 end subroutine set_fldstg
 
end module HYMAP_initMod
