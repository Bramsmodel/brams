  integer, parameter :: maxsrcfiles   = 1500 
  integer, parameter :: nsrc=4  !number_sources
  character(len=20), parameter :: src_name(nsrc)= (/&
       'antro               ', &    
       'bburn               ', &
       'bioge               ', &
       'geoge               '/)

  integer, parameter :: antro = 01 ! anthropogenic sources
  integer, parameter :: bburn = 02 ! biomass burning sources 
  integer, parameter :: bioge = 03 ! biogenic sources 
  integer, parameter :: geoge = 04 ! geogenic/volc sources ! must be equal to "nsrc"

  integer, parameter :: max_ntimes_src = 2  !- number maximum of src files for linterp.
  
