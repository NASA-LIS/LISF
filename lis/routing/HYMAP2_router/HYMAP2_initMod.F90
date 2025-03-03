!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module HYMAP2_initMod

  use LIS_logMod,     only : LIS_logunit
  implicit none
  
contains
  ! ================================================          
  ! ================================================              
  subroutine HYMAP2_grid2vector(nx,ny,nz,nseqall,imis,seqx,seqy,grid,vector)
    implicit none
    integer, intent(in)  :: nx             !number of grids in horizontal
    integer, intent(in)  :: ny             !number of grids in vertical
    integer, intent(in)  :: nz             !number of array layers
    integer, intent(in)  :: nseqall        !length of 1D sequnece for river and mouth
    integer, intent(in)  :: imis           !integer undefined value
    integer, intent(in)  :: seqx(nseqall)  !1D sequence horizontal
    integer, intent(in)  :: seqy(nseqall)  !1D sequence vertical
    real,    intent(in)  :: grid(nx,ny,nz)   !input matrix
    real,    intent(out) :: vector(nseqall,nz) 
    integer              :: i,ix,iy
    
    do i=1,nseqall
      ix=seqx(i)
      iy=seqy(i)
      vector(i,:)=grid(ix,iy,:)
    enddo
  
  end subroutine HYMAP2_grid2vector
  ! ================================================          
  ! ================================================              
  subroutine HYMAP2_vector2grid(nx,ny,nz,nseqall,imis,seqx,seqy,grid,vector)
    implicit none
    integer, intent(in)  :: nx             !number of grids in horizontal
    integer, intent(in)  :: ny             !number of grids in vertical
    integer, intent(in)  :: nz             !number of array layers
    integer, intent(in)  :: nseqall        !length of 1D sequnece for river and mouth
    integer, intent(in)  :: imis           !integer undefined value
    integer, intent(in)  :: seqx(nseqall)  !1D sequence horizontal
    integer, intent(in)  :: seqy(nseqall)  !1D sequence vertical
    real,    intent(out) :: grid(nx,ny,nz)   !input matrix
    real,    intent(in)   :: vector(nseqall,nz) 
    integer              :: i,ix,iy
    
    grid=real(imis)
    do i=1,nseqall
      ix=seqx(i)
      iy=seqy(i)
      grid(ix,iy,:)=vector(i,:)
    enddo
  
  end subroutine HYMAP2_vector2grid
  ! ================================================          
  ! ================================================              
  subroutine HYMAP2_get_vector_size(nx,ny,nxg,nyg,offx,offy,imis,nextx,mask,nseqall)
    implicit none
    integer, intent(in)  ::  nx                !number of grids in horizontal
    integer, intent(in)  ::  ny                !number of grids in vertical
    integer, intent(in)  ::  nxg                
    integer, intent(in)  ::  nyg                
    integer, intent(in)  ::  offx                
    integer, intent(in)  ::  offy               
    integer, intent(in)  ::  imis              !integer undefined value
    integer, intent(in)  ::  nextx(nxg,nyg)      !point downstream horizontal
    integer, intent(in)  ::  mask(nx,ny)       !mask limiting modeled region (0: out; >=1: in)
    integer, intent(out) ::  nseqall           !length of 1D sequnece for river and mouth
    integer              ::  ix,iy
    integer              ::  ix1,iy1
    
    nseqall=0
    do ix=1,nx
      do iy=1,ny
         ix1 = ix + offx - 1
         iy1 = iy + offy -1
         if(nextx(ix1,iy1)/=imis.and.mask(ix,iy)>0) then 
            nseqall=nseqall+1
         endif
      enddo
   enddo
   write(LIS_logunit,*)'[get_vector_size] number of cells',nseqall
   return
   
  end subroutine HYMAP2_get_vector_size
  ! ================================================   
!#if 0        
  ! ================================================              
  subroutine HYMAP2_get_seq(nx,ny,nxg,nyg,offx,offy,nseqall,&
       imis,nextx,nexty,mask,sindex,outlet,seqx,seqy,next)
    use LIS_coreMod
    implicit none
    integer, intent(in)   :: nx             !number of grids in horizontal
    integer, intent(in)   :: ny             !number of grids in vertical
    integer, intent(in)   :: nxg
    integer, intent(in)   :: nyg
    integer, intent(in)   :: offx
    integer, intent(in)   :: offy
    integer, intent(in)   :: nseqall        !length of 1D sequnece for river and mouth
    integer, intent(in)   :: imis           !integer undefined value
    integer, intent(in)   :: nextx(nxg,nyg)   !point downstream horizontal
    integer, intent(in)   :: nexty(nxg,nyg)   !point downstream vertical

    integer, intent(in)   :: mask(nxg,nyg)    !mask limiting modeled region (0: out; >=1: in)

    integer, intent(in)  :: sindex(nxg,nyg)    !2-D sequence index 
    integer, intent(out)  :: outlet(nseqall) !outlet flag: 0 - river; 1 - ocean
    integer, intent(out)  :: seqx(nseqall)   !1D sequence horizontal
    integer, intent(out)  :: seqy(nseqall)   !1D sequence vertical
    integer, intent(out)  :: next(nseqall)   !downstream grid cell
    integer               :: i,j,ix,iy,jx,jy,iloc(2)
    integer               :: ix1, iy1
    i=0

    do ix=1,nx
       do iy=1,ny
         ix1 = ix + offx - 1
         iy1 = iy + offy -1
         if(nextx(ix1,iy1)/=imis.and.mask(ix1,iy1)>0)then
            i=i+1
            seqx(i)=ix
            seqy(i)=iy
            if(nextx(ix1,iy1)>0)then
               outlet(i)=0
            else
               outlet(i)=1
            endif
         endif
      enddo
   enddo
    do i=1,nseqall
      if(outlet(i)==0)then
        ix=seqx(i)
        iy=seqy(i)
        ix1 = ix + offx - 1
        iy1 = iy + offy -1        
        jx=nextx(ix1,iy1)
        jy=nexty(ix1,iy1)
        next(i)=sindex(jx,jy) 
     else
        next(i)=imis
      endif
    enddo    
  end subroutine HYMAP2_get_seq
  ! ================================================              
  subroutine HYMAP2_get_sindex(nxg,nyg,nseqall,imis,nextx,nexty,mask,&
       sindex,outlet,next)
    use LIS_coreMod
    implicit none
    integer, intent(in)   :: nxg
    integer, intent(in)   :: nyg
    integer, intent(in)   :: nseqall
    integer, intent(in)   :: imis           !integer undefined value
    integer, intent(in)   :: nextx(nxg,nyg)   !point downstream horizontal
    integer, intent(in)   :: nexty(nxg,nyg)   !point downstream vertical
    integer, intent(in)   :: mask(nxg,nyg)    !mask limiting modeled region (0: out; >=1: in)

    integer, intent(out)  :: sindex(nxg,nyg)    !2-D sequence index 
    integer, intent(out)  :: outlet(nseqall) !outlet flag: 0 - river; 1 - ocean
    integer, intent(out)  :: next(nseqall)   !downstream grid cell

    integer               :: seqx(nseqall)
    integer               :: seqy(nseqall)
    
    integer               :: i,j,ix,iy,jx,jy,iloc(2)
    integer               :: ix1, iy1
    i=0

    sindex=imis
    do ix=1,nxg
       do iy=1,nyg
         if(nextx(ix,iy)/=imis.and.mask(ix,iy)>0)then
            i=i+1
            seqx(i) = ix
            seqy(i) = iy
            sindex(ix,iy)=i
            if(nextx(ix,iy)>0)then
               outlet(i)=0
            else
               outlet(i)=1
            endif
         endif
      enddo
   enddo
   do i=1,nseqall
      if(outlet(i)==0) then 
         ix = seqx(i)
         iy = seqy(i)
         jx = nextx(ix,iy)
         jy = nexty(ix,iy)
         next(i) = sindex(jx,jy)
      endif
   enddo

 end subroutine HYMAP2_get_sindex
  ! ================================================              
  ! ================================================              
  subroutine HYMAP2_set_fldstg(nz,nseqall,fldhgt,areamat,rivlen,&
                        rivwth,rivstomax,fldstomax,fldgrd,rivare)
    ! ================================================
    ! to   set floodplain staging
    ! by   Augusto Getirana after Dai YAMAZAKI
    ! on   27 Oct 2011
    ! at   NASA,GSFC
    ! Adapted for flow routing implementation in LIS 09 Nov 2011
    ! ================================================
    use LIS_logMod,     only : LIS_logunit

    implicit none
    integer, intent(in)   :: nz                  !number of stages in the sub-grid discretization
    integer, intent(in)  :: nseqall               !length of 1D sequnece for river and mouth
    real,    intent(in)  :: areamat(nseqall)      !area of the grid [m^2]
    real,    intent(in)  :: rivlen(nseqall)       !channel length [m]
    real,    intent(in)  :: rivwth(nseqall)       !river width [m]
    real,    intent(in) :: fldhgt(nseqall,nz)    !floodplain height
    real,    intent(in)  :: rivstomax(nseqall)    !maximum river storage [m3]
    real,    intent(out) :: fldstomax(nseqall,nz) !maximum floodplain storage [m3]
    real,    intent(out) :: fldgrd(nseqall,nz)    !floodplain gradient [-]
    real,    intent(out) :: rivare(nseqall)

    integer              :: iseq, i
    real                 :: stonow
    real                 :: stopre
    real                 :: hgtpre
    real                 :: wthinc
    integer              :: i1,i2
    real                 :: rivfrc                !fraction of river surface area [-]
    real                 :: hgtprv(nseqall)
    ! ================================================
!ag (24Feb2020)
!    do iseq=1,nseqall
!       !ix=seqx(iseq)
!       !iy=seqy(iseq)
!       rivfrc=rivare(iseq)/areamat(iseq)
!       i1=max(1,int(rivfrc*nz))
!       i2=min(nz,i1+1)
!       hgtprv(iseq)=fldhgt(iseq,i1)+(fldhgt(iseq,i2)-fldhgt(iseq,i1))*(rivfrc*nz-i1)
!       fldhgt(iseq,:)=fldhgt(iseq,:)-hgtprv(iseq)
!       where(fldhgt(iseq,:)<0)fldhgt(iseq,:)=0.
!    enddo

    fldstomax=0.
    fldgrd=0.
    do iseq=1,nseqall
        !ag (24Feb2020)
        rivare(iseq)=rivwth(iseq)*rivlen(iseq)
        stopre=rivstomax(iseq)
        hgtpre=0.
        !ag (12Feb2020)
        !wthinc=areamat(iseq)/rivlen(iseq)*0.1
        wthinc=(areamat(iseq)-rivare(iseq))/rivlen(iseq)/nz
        do i=1,nz
           stonow=rivlen(iseq)*(rivwth(iseq)+wthinc*(real(i)-0.5))*(fldhgt(iseq,i)-hgtpre)
           fldstomax(iseq,i)=stopre+stonow
           !           fldgrd(iseq,i)=(fldhgt(iseq,i)-hgtpre) / wthinc
           if(wthinc <=0.0) then
!              write(LIS_logunit,*) '[WARN] set_fldg withinc= 0'
              fldgrd(iseq,i) = 0.0
           else
              fldgrd(iseq,i) = (fldhgt(iseq,i)-hgtpre)/wthinc
           endif
           stopre=fldstomax(iseq,i)
           hgtpre=fldhgt(iseq,i)
        enddo
     enddo

 end subroutine HYMAP2_set_fldstg

end module HYMAP2_initMod
