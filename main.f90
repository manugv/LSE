!Calculate sum over all snapshots : <u_iu_j> for
!And calculate linear stochastic estimate coefficients

program two_point_correlations
   use mod_precision          !define precision
   use mod_geometry           !define all geometry
   use mod_fftw3
   use mod_routines
   use mod_output
   use mod_lse_coeff
   implicit none

   !**************************************************
   
   !define correlations data types
   type(correlate)  :: ru,rv,rw

   !variable defining the event plane
   type(plane)      :: event

   !other variables
   integer          :: k, variable, ierror
   character(len=7) :: filename

   !define velocity variables
   real(kind=dp), dimension(:,:,:), allocatable :: uc, vc, wc

   !define average variables 
   real(kind=dp), dimension(kcell)   :: um, vm
   real(kind=dp), dimension(kcell+1) :: wm

   !**************************************************
   
   !allocate velocity arrays at same point
   allocate(uc(icell,jcell,kcell))
   allocate(vc(icell,jcell,kcell))
   allocate(wc(icell,jcell,kcell))
   
   !set all correlations to zero
   ru%u = 0.0_dp
   ru%v = 0.0_dp
   ru%w = 0.0_dp
   rv%u = 0.0_dp
   rv%v = 0.0_dp
   rv%w = 0.0_dp
   rw%u = 0.0_dp
   rw%v = 0.0_dp
   rw%w = 0.0_dp

   !Number of files, start of the file and end of the file 
   first_file  = 35
   last_file   = 200  
   total_files = 166

   !Input event plane
   !on run
   ! write(*,*)"grid number"
   ! read(*,*)event%loc
   !define at beginning
   event%loc = 36
   write(*,*)"event loc",event%loc

   !initalize all FFTs for computing correlations
   call init_fft()
   
   !Read average velocities and pressure from average file
   call load_average(um,vm,wm)
   write(*,*)"Read average",maxval(um),maxval(vm),maxval(wm)

   !**************************************************
   
   read_files: do variable = first_file, last_file
      write(filename,'(i7.7)')variable*1000
      write(*,*)filename

      !Reading data files
      !subtract mean velocity
      !interpolate to cell centers in wall-normal direction
      call load_velocity(filename,um,vm,wm,uc,vc,wc)
      write(*,*)maxval(uc),maxval(vc),maxval(wc)

      !2-d FFT of the velocities of the event planes
      !'fft_ref_plane' function in mod_routines
      !FFTs are dumped into event%u,v,w planes 
      call fft_ref_plane(uc(:,:,event%loc),vc(:,:,event%loc),wc(:,:,event%loc),event%u,event%v,event%w)

      !Compute 2-d correlations for all planes in wall-normal direction
      !'fft2d' function in mod_routines computes all cross correlations 
      !with the event plane and dumped into ru%tmp, rv%tmp and rw%tmp
      !Further these are added to the correlation average ru,rv,rw
      
      shift_wall: do k=1,kcell
         !Compute FFTs of velocity and multiply by event plane FFTs
         !Then do inverse FFTs to calculate correlations ru%tmp, rv%tmp, rw%tmp
         call fft2d(uc(:,:,k),vc(:,:,k),wc(:,:,k),event%u,event%v,event%w,ru%tmp,rv%tmp,rw%tmp)

         !Compute correlation average
         ru%u(:,:,k) = ru%u(:,:,k) + ru%tmp(:,:,1)
         ru%v(:,:,k) = ru%v(:,:,k) + ru%tmp(:,:,2)
         ru%w(:,:,k) = ru%w(:,:,k) + ru%tmp(:,:,3)
         rv%u(:,:,k) = rv%u(:,:,k) + rv%tmp(:,:,1)
         rv%v(:,:,k) = rv%v(:,:,k) + rv%tmp(:,:,2)
         rv%w(:,:,k) = rv%w(:,:,k) + rv%tmp(:,:,3)
         rw%u(:,:,k) = rw%u(:,:,k) + rw%tmp(:,:,1)
         rw%v(:,:,k) = rw%v(:,:,k) + rw%tmp(:,:,2)
         rw%w(:,:,k) = rw%w(:,:,k) + rw%tmp(:,:,3)
      end do shift_wall
   end do read_files   !End of main loop
   
   !Normalize the averages
   ru%u = ru%u/real(icell*jcell*total_files)
   rv%u = rv%u/real(icell*jcell*total_files)
   rw%u = rw%u/real(icell*jcell*total_files)
   ru%v = ru%v/real(icell*jcell*total_files)
   rv%v = rv%v/real(icell*jcell*total_files)
   rw%v = rw%v/real(icell*jcell*total_files)
   ru%w = ru%w/real(icell*jcell*total_files)
   rv%w = rv%w/real(icell*jcell*total_files)
   rw%w = rw%w/real(icell*jcell*total_files)

   !Save all final output files (save_data subroutine in mod_output.f90)
   call save_data(ru%u,ru%v,ru%w,postdata//"final_u")
   call save_data(rv%u,rv%v,rv%w,postdata//"final_v")
   call save_data(rw%u,rw%v,rw%w,postdata//"final_w")

   !Deallocate FFTs and velocity variables
   deallocate(uc,vc,wc)
   call destroy_fft          !In file mod_fftw3.f90

   !FINALLY
   !CALCULATE LSE COEFFICIENTS (function 'lse' in in file mod_lse_coeff)
   call lse(ru%u,ru%v,ru%w,rv%u,rv%v,rv%w,rw%u,rw%v,rw%w,event%loc)

end program two_point_correlations
