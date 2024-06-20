MODULE mod_output
  use mod_precision
  use mod_geometry
  implicit none
contains
  

  subroutine load_average(um,vm,wm)
    implicit none
    real(kind=dp), dimension(:), intent(out) :: um,vm
    real(kind=dp), dimension(:), intent(out) :: wm
    open(10, FILE=predata//'stats/average.bin', FORM='unformatted',ACCESS='stream')
    read(10)um,vm,wm
    close(10)
  end subroutine load_average
  
  
  subroutine load_velocity(filename,um,vm,wm,u,v,w)
    use mod_compact
    use mod_lapack
    implicit none
    real(kind=dp), dimension(:), intent(in) :: um,vm
    real(kind=dp), dimension(:), intent(in) :: wm
    real(kind=dp), dimension(:,:,:), intent(out)  :: u,v,w
    real(kind=dp), dimension(icell,jcell,kcell+1) :: wn
    character(len=*), intent(in) :: filename
    real(kind=8)       :: time,itr
    integer            :: i,j,k, iol
    real(kind=dp),dimension(kcell) :: dia
    real(kind=dp),dimension(kcell-1) :: sub,super

    !Initialize tridiagonal matrix
    call interpolate_node_cell_lhs(kcell,sub,dia,super)

    inquire(iolength=iol) time
    open(10, FILE=predata//filename, FORM='unformatted', &
         ACCESS='DIRECT', RECL=(2*(icell*jcell*kcell)+(icell*jcell*(kcell+1))+2)*iol )
    read(10,rec=1)(((u(i,j,k),i=1,icell),j=1,jcell),k=1,kcell), &
         (((v(i,j,k),i=1,icell),j=1,jcell),k=1,kcell), &
         (((wn(i,j,k),i=1,icell),j=1,jcell),k=1,kcell+1), &
         time,itr
    close(10)
    do k=1,kcell
       u(:,:,k) = u(:,:,k) - um(k)
       v(:,:,k) = v(:,:,k) - vm(k)
    enddo
    do k=1,kcell+1
       wn(:,:,k) = wn(:,:,k) - wm(k)
    enddo
    !Interpolate
    do j=1,jcell
       do i=1,icell
          call interpolate_node_cell_rhs(wn(i,j,:),w(i,j,:))
          call solve_tri(sub,dia,super,w(i,j,:))
       enddo
    enddo
    
    return
  end subroutine load_velocity
  

  SUBROUTINE save_data(a,b,c,filename)
    real(kind=dp), dimension(:,:,:),intent(in) :: a,b,c
    character(len=*), intent(in) :: filename
    open(10,file=filename,access='stream',form="unformatted")
    write(10)a,b,c
    close(10)
    return
  END SUBROUTINE save_data


END MODULE mod_output
