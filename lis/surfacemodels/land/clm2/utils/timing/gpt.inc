!
! Fortran include file to match gpt.h entries.
! Note that the values must match
!
      integer usrsys
      integer wall
      integer pcl_l1dcache_miss
      integer pcl_l2cache_miss
      integer pcl_cycles
      integer pcl_elapsed_cycles
      integer pcl_fp_instr
      integer pcl_loadstore_instr
      integer pcl_instr
      integer pcl_stall

      parameter (usrsys              = 1)
      parameter (wall                = 2)
      parameter (pcl_l1dcache_miss   = 4)
      parameter (pcl_l2cache_miss    = 5)
      parameter (pcl_cycles          = 6)
      parameter (pcl_elapsed_cycles  = 7)
      parameter (pcl_fp_instr        = 8)
      parameter (pcl_loadstore_instr = 9)
      parameter (pcl_instr           = 10)
      parameter (pcl_stall           = 11)
