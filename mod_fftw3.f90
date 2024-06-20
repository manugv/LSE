module mod_fftw3
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
  type(C_PTR) :: fplan,bplan

contains

  subroutine init_fft
    use mod_geometry ,only : icell,jcell
    implicit none 
    call fftw3_fplan(icell,jcell,fplan)
    call fftw3_bplan(icell,jcell,bplan)
  end subroutine init_fft

  subroutine fftw3_fplan(nx,ny,plan)
    complex(C_DOUBLE_COMPLEX), dimension(:,:),allocatable :: output
    real(C_DOUBLE), dimension(:,:), allocatable :: input
    integer(C_INT), intent(in) :: nx,ny
    type(C_PTR) :: plan
    allocate(input(nx,ny))
    allocate(output(nx/2+1,ny))
    plan = fftw_plan_dft_r2c_2d(ny,nx,input,output,FFTW_EXHAUSTIVE)
    deallocate(input,output)
  end subroutine fftw3_fplan

  subroutine fftw3_bplan(nx,ny,plan)
    complex(C_DOUBLE_COMPLEX), dimension(:,:),allocatable :: input
    real(C_DOUBLE), dimension(:,:), allocatable :: output
    integer(C_INT), intent(in) :: nx,ny
    type(C_PTR) :: plan
    allocate(output(nx,ny))
    allocate(input(nx/2+1,ny))
    plan = fftw_plan_dft_c2r_2d(ny,nx,input,output,FFTW_EXHAUSTIVE)
    deallocate(input,output)
  end subroutine fftw3_bplan
  
  subroutine destroy_fft
    call fftw_destroy_plan(fplan)
    call fftw_destroy_plan(bplan)
    return
  end subroutine destroy_fft
  
  subroutine fft_forward(input,output)
    complex(C_DOUBLE_COMPLEX), dimension(:,:),intent(out) :: output
    real(C_DOUBLE), dimension(:,:) :: input
    output = 0.0d0
    call fftw_execute_dft_r2c(fplan,input,output)
    return
  end subroutine fft_forward
  
  subroutine fft_backward(nx,ny,input,output)
    integer,intent(in) :: nx,ny
    complex(C_DOUBLE_COMPLEX), dimension(:,:) :: input
    real(C_DOUBLE), dimension(:,:),intent(out) :: output
    output = 0.0d0
    call fftw_execute_dft_c2r(bplan,input,output)
    output = output/real(nx*ny)
    return
  end subroutine fft_backward
  
end module mod_fftw3
