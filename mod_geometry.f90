module mod_geometry
  use mod_precision
  implicit none
  
  !Definition of domain size
  integer, parameter :: icell = 192, jcell = 128, kcell = 128
  real(kind=dp), parameter :: pi = acos(-1.0_dp)
  !Length of the domain
  real(kind=dp), parameter :: lx=2.0_dp*pi,ly=lx/3.0_dp,lz=1.0_dp
  
  !Inverse of dimensions of a cell 
  real(kind=dp), parameter :: dxi = (icell*1.0_dp)/lx, dyi = (jcell*1.0_dp)/ly, dzi = (kcell*1.0_dp)/lz  

  real(kind=dp), parameter :: dx = lx/(icell*1.0_dp), dy = ly/(1.0_dp*jcell), dz = lz/(1.0_dp*kcell)
  !Ratio for ffts
  real(kind=dp), parameter :: lxx = (2.0_dp*pi)/lx , lyy = (2.0_dp*pi)/ly

  !Velocity fields data directory 
  character(len=36), parameter :: predata = "/tennekes/manu/data/dns/Re360/sim_2/"
  character(len=39) :: postdata = "/tennekes/manu/data/lse/Re360/yplus_75/"

  !Number of files and first file number
  integer, parameter :: first_file
  integer, parameter :: last_file
  integer, parameter :: total_files

  type correlate
     real(kind=dp), dimension(icell,jcell,kcell) :: u,v,w
     real(kind=dp), dimension(icell,jcell,3)    :: tmp
  end type correlate

  type plane
     complex(kind=dp), dimension(icell/2+1,jcell) :: u,v,w
     integer :: loc
  end type plane

end module mod_geometry
