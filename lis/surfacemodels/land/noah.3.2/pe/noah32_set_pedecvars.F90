!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noah32_set_pedecvars
!  \label{noah32_set_pedecvars}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah32_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_soilsMod,  only : LIS_soils
  use LIS_logMod,       only : LIS_logunit,LIS_verify
  use noah32_lsmMod, only : noah32_struc
  use noah32_peMod,  only : noah32_pe_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Noah's model variables. 
! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: varField
  type(ESMF_Field)       :: feasField
  real, allocatable          :: vardata(:)
  integer, allocatable       :: mod_flag(:)
  real                   :: vardata_min, vardata_max
  integer                :: i,t
  integer                :: status

  n = 1

  call ESMF_StateGet(Feas_State, "Feasibility Flag", feasField, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(feasField,localDE=0,farrayPtr=mod_flag,rc=status)
  call LIS_verify(status)

  mod_flag = 0
  
  do i=1,noah32_pe_struc(n)%nparams
     
     if(noah32_pe_struc(n)%param_select(i).eq.1) then 
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."SMCMAX")) then 

           call ESMF_StateGet(DEC_State,"SMCMAX",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%smcmax = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."PSISAT")) then 

           call ESMF_StateGet(DEC_State,"PSISAT",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%psisat = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."DKSAT")) then 

           call ESMF_StateGet(DEC_State,"DKSAT",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%dksat = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."DWSAT")) then 

           call ESMF_StateGet(DEC_State,"DWSAT",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%dwsat = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."BEXP")) then 

           call ESMF_StateGet(DEC_State,"BEXP",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%bexp = vardata(t) 
              endif
           enddo
              
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."QUARTZ")) then 

           call ESMF_StateGet(DEC_State,"QUARTZ",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%quartz = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."RSMIN")) then 

           call ESMF_StateGet(DEC_State,"RSMIN",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%rsmin = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."RGL")) then 

           call ESMF_StateGet(DEC_State,"RGL",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%rgl = vardata(t) 
              endif
           enddo
           
        endif
        
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."HS")) then 

           call ESMF_StateGet(DEC_State,"HS",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%hs = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."Z0")) then 

           call ESMF_StateGet(DEC_State,"Z0",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%z0 = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."LAI")) then 

           call ESMF_StateGet(DEC_State,"LAI",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%lai = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."CFACTR")) then 

           call ESMF_StateGet(DEC_State,"CFACTR",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%cfactr = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."CMCMAX")) then 

           call ESMF_StateGet(DEC_State,"CMCMAX",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%cmcmax = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."SBETA")) then 

           call ESMF_StateGet(DEC_State,"SBETA",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%sbeta = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."RSMAX")) then 

           call ESMF_StateGet(DEC_State,"RSMAX",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%rsmax = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."TOPT")) then 

           call ESMF_StateGet(DEC_State,"TOPT",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%topt = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."REFDK")) then 

           call ESMF_StateGet(DEC_State,"REFDK",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%refdk = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."FXEXP")) then 

           call ESMF_StateGet(DEC_State,"FXEXP",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%fxexp = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."REFKDT")) then 

           call ESMF_StateGet(DEC_State,"REFKDT",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%refkdt = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."CZIL")) then 

           call ESMF_StateGet(DEC_State,"CZIL",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%czil = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."CSOIL")) then 

           call ESMF_StateGet(DEC_State,"CSOIL",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%csoil = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."FRZK")) then 

           call ESMF_StateGet(DEC_State,"FRZK",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%frzk = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."SNUP")) then 

           call ESMF_StateGet(DEC_State,"SNUP",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%snup = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."SMCREF")) then 

           call ESMF_StateGet(DEC_State,"SMCREF",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%smcref = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."SMCDRY")) then 

           call ESMF_StateGet(DEC_State,"SMCDRY",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%smcdry = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."SMCWLT")) then 

           call ESMF_StateGet(DEC_State,"SMCWLT",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%smcwlt = vardata(t) 
              endif
           enddo
           
        endif

         if((trim(noah32_pe_struc(n)%param_name(i)).eq."F1")) then 

           call ESMF_StateGet(DEC_State,"F1",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%f1 = vardata(t) 
              endif
           enddo
           
        endif

         if((trim(noah32_pe_struc(n)%param_name(i)).eq."SLOPE")) then 

           call ESMF_StateGet(DEC_State,"SLOPE",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%slope = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."EMISS")) then 

           call ESMF_StateGet(DEC_State,"EMISS",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%emiss = vardata(t) 
              endif
           enddo
           
        endif
        if((trim(noah32_pe_struc(n)%param_name(i)).eq."SIGMA_FLX")) then 

           call ESMF_StateGet(DEC_State,"SIGMA_FLX",varField,rc=status)
           call LIS_verify(status)
           
           call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
           call LIS_verify(status)
           call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
           call LIS_verify(status)

           do t=1,LIS_rc%ntiles(n)
              if(vardata(t).lt.vardata_min) then 
                 vardata(t) = vardata_min
                 mod_flag(t) = 1
              endif
              if(vardata(t).gt.vardata_max) then 
                 vardata(t) = vardata_max
                 mod_flag(t) = 1
              endif
           enddo

           do t=1,LIS_rc%ntiles(n)
              if(mod_flag(t).eq.0) then 
                 noah32_struc(n)%noah(t)%sigma_flx = vardata(t) 
              endif
           enddo
           
        endif
     endif
  enddo


end subroutine noah32_set_pedecvars



