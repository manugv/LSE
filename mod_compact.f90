module mod_compact
   use mod_precision
   use mod_geometry ,only : dzi
   implicit none
contains

   subroutine interpolate_node_cell_lhs(n,sub,dia,super)
      integer, intent(in) :: n
      real(kind=dp), dimension(n), intent(out) :: dia
      real(kind=dp), dimension(n-1), intent(out) :: sub,super
      integer :: i
      dia = 0._dp
      sub = 0._dp
      super = 0._dp
      !Fourth order
      dia(1) = 1._dp
      super(1) = 0._dp
      !Sixth order
      do i=2,n-1
         sub(i-1) = 3._dp/10._dp
         dia(i)   = 1._dp
         super(i) = 3._dp/10._dp
      enddo
      !Fourth order
      i = n
      sub(i-1) = 0._dp
      dia(i)   = 1._dp
      return
   end subroutine interpolate_node_cell_lhs


   subroutine interpolate_node_cell_rhs(f,fi)
      real(kind=dp), dimension(:), intent(in) :: f
      real(kind=dp), dimension(:), intent(out) :: fi
      integer :: i,n
      real(kind=dp) :: c1,c2,c3,c4
      n = size(fi)
      fi = 0._dp
      fi(1) = (5._dp/16._dp)*f(1)+(15._dp/16._dp)*f(2)-(5._dp/16._dp)*f(3)+(1._dp/16._dp)*f(4)
      do i=2,n-1
         fi(i) = f(i+2)*(1._dp/20._dp)+f(i+1)*(3._dp/4._dp)+f(i)*(3._dp/4._dp)+f(i-1)*(1._dp/20._dp)
      enddo
      !Fourth order
      i = n
      fi(i)=  (5._dp/16._dp)*f(i+1)+(15._dp/16._dp)*f(i)-(5._dp/16._dp)*f(i-1)+(1._dp/16._dp)*f(i-2)
      return
   end subroutine interpolate_node_cell_rhs


   subroutine solve_tri(m_sub,m_dia,m_super,x)
      real(kind=dp),dimension(:),intent(in)    :: m_sub,m_dia,m_super
      real(kind=dp),dimension(:),intent(inout) :: x
      real(kind=dp),dimension(:),allocatable   :: sub,dia,super
      integer :: info,n
      n = size(m_dia)
      allocate(sub(n-1))
      allocate(dia(n))
      allocate(super(n-1))
      sub   = 0.0d0
      dia   = 0.0d0
      super = 0.0d0
      sub   = m_sub
      dia   = m_dia
      super = m_super
      call dgtsv(n,1,sub,dia,super,x,n,info)
      if (info < 0) then 
         write(*,*)"linear system solver failed"
      elseif(info > 0) then
         write(*,*)"linear system solver weird behaviour"
      endif
      deallocate(sub,super,dia)
      return
   end subroutine solve_tri

end module mod_compact
