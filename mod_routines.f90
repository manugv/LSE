module mod_routines
  use mod_precision
  use mod_geometry
  use mod_fftw3
  implicit none
contains
  
  subroutine fft_ref_plane(pu,pv,pw,event_u,event_v,event_w)
    real(kind=dp), dimension(:,:), intent(in) :: pu,pv,pw
    complex(kind=dp), dimension(:,:), intent(out) :: event_u,event_v,event_w
    integer :: i,j
    real(kind=dp), dimension(icell,jcell) :: tmp_u,tmp_v,tmp_w
    call fft_forward(pu,event_u)
    call fft_forward(pv,event_v)
    call fft_forward(pw,event_w)
    event_u = conjg(event_u)
    event_v = conjg(event_v)
    event_w = conjg(event_w)
    return
  end subroutine fft_ref_plane
  
  subroutine fft2d(pu, pv, pw,event_u, event_v, event_w, u_tmp, v_tmp, w_tmp)
    real(kind=dp), dimension(:,:), intent(in)     :: pu,pv,pw
    complex(kind=dp), dimension(:,:), intent(in)  :: event_u,event_v,event_w
    real(kind=dp), dimension(:,:,:), intent(out)  :: u_tmp,v_tmp,w_tmp
    complex(kind=dp), dimension(icell/2+1,jcell) :: ctmp_u,ctmp_v,ctmp_w,temp
    integer :: i,j
    call fft_forward(pu,ctmp_u)
    call fft_forward(pv,ctmp_v)
    call fft_forward(pw,ctmp_w)
    temp = ctmp_u*event_u
    call fft_backward(icell,jcell,temp,u_tmp(:,:,1))
    temp = ctmp_u*event_v
    call fft_backward(icell,jcell,temp,u_tmp(:,:,2))
    temp = ctmp_u*event_w
    call fft_backward(icell,jcell,temp,u_tmp(:,:,3))
    temp = ctmp_v*event_u
    call fft_backward(icell,jcell,temp,v_tmp(:,:,1))
    temp = ctmp_v*event_v
    call fft_backward(icell,jcell,temp,v_tmp(:,:,2))
    temp = ctmp_v*event_w
    call fft_backward(icell,jcell,temp,v_tmp(:,:,3))
    temp = ctmp_w*event_u
    call fft_backward(icell,jcell,temp,w_tmp(:,:,1))
    temp = ctmp_w*event_v
    call fft_backward(icell,jcell,temp,w_tmp(:,:,2))
    temp = ctmp_w*event_w
    call fft_backward(icell,jcell,temp,w_tmp(:,:,3))
    return
  end subroutine fft2d

end module mod_routines
