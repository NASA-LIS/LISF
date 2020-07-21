!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
! 
! Copyright (c) 2015 United States Government as represented by the
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
         if(nextx(ix1,iy1)/=imis.and.mask(ix,iy)>0)nseqall=nseqall+1
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
    !real,    intent(in)   :: uparea(nx,ny)  !upstream area
    integer, intent(in)   :: mask(nxg,nyg)    !mask limiting modeled region (0: out; >=1: in)

    integer, intent(in)  :: sindex(nxg,nyg)    !2-D sequence index 
    integer, intent(out)  :: outlet(nseqall) !outlet flag: 0 - river; 1 - ocean
    integer, intent(out)  :: seqx(nseqall)   !1D sequence horizontal
    integer, intent(out)  :: seqy(nseqall)   !1D sequence vertical
    integer, intent(out)  :: next(nseqall)   !downstream grid cell

    !integer               :: tmp(nx,ny)   
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
!        if(jx.gt.nx.or.jy.gt.ny) then 
!           print*, 'H2:out of bounds ',ix,iy,sindex(jx,jy),LIS_localPet
!        endif
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

#if 0
  ! ================================================
  subroutine HYMAP2_get_seq_glb(nx,ny,nseqall,imis,nextx,nexty,mask,sindex,outlet,seqx,seqy,next)
    implicit none
    integer, intent(in)   :: nx             !number of grids in horizontal
    integer, intent(in)   :: ny             !number of grids in vertical
    integer, intent(in)   :: nseqall        !length of 1D sequnece for river and mouth
    integer, intent(in)   :: imis           !integer undefined value
    integer, intent(in)   :: nextx(nx,ny)   !point downstream horizontal
    integer, intent(in)   :: nexty(nx,ny)   !point downstream vertical
    !real,    intent(in)   :: uparea(nx,ny)  !upstream area
    integer, intent(in)   :: mask(nx,ny)    !mask limiting modeled region (0: out; >=1: in)

    integer, intent(out)  :: sindex(nx,ny)    !2-D sequence index
    integer, intent(out)  :: outlet(nseqall) !outlet flag: 0 - river; 1 - ocean
    integer, intent(out)  :: seqx(nseqall)   !1D sequence horizontal
    integer, intent(out)  :: seqy(nseqall)   !1D sequence vertical
    integer, intent(out)  :: next(nseqall)   !downstream grid cell

    !integer               :: tmp(nx,ny)
    integer               :: i,j,ix,iy,jx,jy,iloc(2)
    i=0
    sindex=imis
    do ix=1,nx
      do iy=1,ny
        if(nextx(ix,iy)/=imis.and.mask(ix,iy)>0)then
          i=i+1
          seqx(i)=ix
          seqy(i)=iy
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
      if(outlet(i)==0)then
        ix=seqx(i)
        iy=seqy(i)
        jx=nextx(ix,iy)
        jy=nexty(ix,iy)
        next(i)=sindex(jx,jy)
      else
        next(i)=imis
      endif
    enddo
  end subroutine HYMAP2_get_seq_glb
#endif
  ! ================================================          
  ! ================================================              
  !No longer used
  subroutine HYMAP2_get_seq0(nx,ny,nseqall,imis,uparea,nextx,nexty,mask,sindex,outlet,seqx,seqy,next)
    implicit none
    integer, intent(in)   :: nx             !number of grids in horizontal
    integer, intent(in)   :: ny             !number of grids in vertical
    integer, intent(in)   :: nseqall        !length of 1D sequnece for river and mouth
    integer, intent(in)   :: imis           !integer undefined value
    integer, intent(in)   :: nextx(nx,ny)   !point downstream horizontal
    integer, intent(in)   :: nexty(nx,ny)   !point downstream vertical
    real,    intent(in)   :: uparea(nx,ny)  !upstream area
    integer, intent(in)   :: mask(nx,ny)    !mask limiting modeled region (0: out; >=1: in)

    integer, intent(out)  :: sindex(nx,ny)    !2-D sequence index 
    integer, intent(out)  :: outlet(nseqall) !outlet flag: 0 - river; 1 - ocean
    integer, intent(out)  :: seqx(nseqall)   !1D sequence horizontal
    integer, intent(out)  :: seqy(nseqall)   !1D sequence vertical
    integer, intent(out)  :: next(nseqall)   !downstream grid cell

    integer               :: tmp(nx,ny)   
    integer               :: i,j,ix,iy,jx,jy,iloc(2)
    
    tmp=1
    i=0
    sindex=imis
    do while(count(tmp/=imis.and.mask/=imis.and.nextx/=imis.and.uparea/=imis)>0)
      i=i+1
      write(*,'(10i8)')i,count(tmp/=imis.and.mask/=imis.and.nextx/=imis.and.uparea/=imis),&
        count(tmp/=imis),count(mask/=imis),count(nextx/=imis),count(uparea/=imis)
      iloc(:)=minloc(uparea,tmp/=imis.and.mask/=imis.and.nextx/=imis.and.uparea/=imis)
      ix=iloc(1)
      iy=iloc(2)
      seqx(i)=ix
      seqy(i)=iy
      sindex(ix,iy)=i
      tmp(ix,iy)=imis
      if(nextx(ix,iy)>0)then
        outlet(i)=0
      else
        outlet(i)=1
      endif
    enddo
    
    do i=1,nseqall
      if(outlet(i)==0)then
        ix=seqx(i)
        iy=seqy(i)
        jx=nextx(ix,iy)
        jy=nexty(ix,iy)
        next(i)=sindex(jx,jy) 
      else
        next(i)=imis
      endif
    enddo
  
  end subroutine HYMAP2_get_seq0
  ! ================================================          
  ! ================================================              
  !No longer used
  subroutine HYMAP2_calc_1d_seq(nx,ny,imis,nextx,mask,seqx,seqy,nseqriv,nseqall)
    ! ================================================
    ! to calculate 1d sequnece file
    ! by Dai YAMAZAKI
    ! on 24th Mar 2008
    ! at IIS, UT
    ! ================================================
    use LIS_logMod,     only : LIS_logunit
    
    implicit none
    ! Input
    integer             ::  nx                  !number of grids in horizontal
    integer             ::  ny                  !number of grids in vertical
    integer             ::  imis                !integer undefined value
    integer             ::  nextx(nx,ny)      !point downstream horizontal
    integer             ::  mask(nx,ny)       !mask limiting modeled region (0: out; >=1: in)
    ! Output
    integer             ::  seqx(nx*ny)       !1D sequence horizontal
    integer             ::  seqy(nx*ny)       !1D sequence vertical
    integer             ::  nseqriv             !length of 1D sequnece for river
    integer             ::  nseqall             !length of 1D sequnece for river and mouth
    ! Local
    integer             ::  ix
    integer             ::  iy
    integer             ::  iseq
    integer             ::  i
    ! ================================================
    iseq=0
    do iy=1, ny
       do ix=1, nx		
          if( nextx(ix,iy).gt.0 .and.mask(ix,iy)>0 )then
             iseq=iseq+1
             seqx(iseq)=ix
             seqy(iseq)=iy
          endif
       end do
    end do  
    nseqriv=iseq   
    do iy=1, ny
       do ix=1, nx
          if( nextx(ix,iy).lt.0 .and. nextx(ix,iy).ne.imis .and.mask(ix,iy)>0 )then
             iseq=iseq+1
             seqx(iseq)=ix
             seqy(iseq)=iy
          endif
       end do
    end do
    nseqall=iseq    
    write(LIS_logunit,*)'[calc_1d_seq] number of cells',nseqall   
    return
  end subroutine HYMAP2_calc_1d_seq
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
    !integer, intent(in)   :: nx                  !number of grids in horizontal
    !integer, intent(in)   :: ny                  !number of grids in vertical
    integer, intent(in)   :: nz                  !number of stages in the sub-grid discretization
    !integer, intent(in)   :: seqx(nx*ny)         !1D sequence horizontal
    !integer, intent(in)   :: seqy(nx*ny)         !1D sequence vertical
    integer, intent(in)  :: nseqall               !length of 1D sequnece for river and mouth
    real,    intent(in)  :: areamat(nseqall)      !area of the grid [m^2]
    real,    intent(in)  :: rivlen(nseqall)       !channel length [m]
    real,    intent(in)  :: rivwth(nseqall)       !river width [m]
    real,    intent(out) :: fldhgt(nseqall,nz)    !floodplain height
    real,    intent(in)  :: rivstomax(nseqall)    !maximum river storage [m3]
    real,    intent(out) :: fldstomax(nseqall,nz) !maximum floodplain storage [m3]
    real,    intent(out) :: fldgrd(nseqall,nz)    !floodplain gradient [-]

    integer              :: iseq, i
    real                 :: stonow
    real                 :: stopre
    real                 :: hgtpre
    real                 :: wthinc
    integer              :: i1,i2
    real                 :: rivare(nseqall)
    real                 :: rivfrc                !fraction of river surface area [-]
    real                 :: hgtprv(nseqall)
    ! ================================================
    do iseq=1,nseqall
       !ix=seqx(iseq)
       !iy=seqy(iseq)
       rivfrc=rivare(iseq)/areamat(iseq)
       i1=max(1,int(rivfrc*nz))
       i2=min(nz,i1+1)
       hgtprv(iseq)=fldhgt(iseq,i1)+(fldhgt(iseq,i2)-fldhgt(iseq,i1))*(rivfrc*nz-i1)
       fldhgt(iseq,:)=fldhgt(iseq,:)-hgtprv(iseq)
       where(fldhgt(iseq,:)<0)fldhgt(iseq,:)=0.
    enddo

    fldstomax=0.
    fldgrd=0.   
    do iseq=1,nseqall
        !ix=seqx(iseq)
        !iy=seqy(iseq)
        stopre=rivstomax(iseq)
        hgtpre=0.
        wthinc=areamat(iseq)/rivlen(iseq)*0.1
        do i=1,nz
           stonow=rivlen(iseq)*(rivwth(iseq)+wthinc*(real(i)-0.5))*(fldhgt(iseq,i)-hgtpre)
           fldstomax(iseq,i)=stopre+stonow
           fldgrd(iseq,i)=(fldhgt(iseq,i)-hgtpre) / wthinc
           stopre=fldstomax(iseq,i)
           hgtpre=fldhgt(iseq,i)
        enddo
     enddo

 end subroutine HYMAP2_set_fldstg
 
end module HYMAP2_initMod
