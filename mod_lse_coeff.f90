module mod_lse_coeff
  use mod_precision
  use mod_output
  use mod_geometry
  implicit none
contains

  subroutine lse(ruu,rvu,rwu,ruv,rvv,rwv,ruw,rvw,rww,yplus)
    real(kind=dp), dimension(:,:,:), intent(in) :: ruu,rvu,rwu
    real(kind=dp), dimension(:,:,:), intent(in) :: ruv,rvv,rwv
    real(kind=dp), dimension(:,:,:), intent(in) :: ruw,rvw,rww
    integer, intent(in)  :: yplus
    real(kind=dp), dimension(:,:,:), allocatable :: l11,l12,l13
    real(kind=dp), dimension(:,:,:), allocatable :: l21,l22,l23
    real(kind=dp), dimension(:,:,:), allocatable :: l31,l32,l33
    integer       :: i, j, k
    real(kind=dp) :: detA,a1,a2,a3,b1,b2,b3,c1,c2,c3,sigma_u_2,sigma_v_2,sigma_w_2
    real(kind=dp) :: temp1,temp2,temp3,temp4,temp5,temp6

    allocate(l11(icell,jcell,kcell))
    allocate(l12(icell,jcell,kcell))
    allocate(l13(icell,jcell,kcell))
    allocate(l21(icell,jcell,kcell))
    allocate(l22(icell,jcell,kcell))
    allocate(l23(icell,jcell,kcell))
    allocate(l31(icell,jcell,kcell))
    allocate(l32(icell,jcell,kcell))
    allocate(l33(icell,jcell,kcell))
    
    !Read from files(do loop if many files exist)
    !call read_data(ruu,rvu,rwu,datadir//"final_u")
    write(*,*)maxval(ruu),maxval(rvu),maxval(rwu)
    !call read_data(ruv,rvv,rwv,datadir//"final_v")
    write(*,*)maxval(ruv),maxval(rvv),maxval(rwv)
    !call read_data(ruw,rvw,rww,datadir//"final_w")
    write(*,*)maxval(ruw),maxval(rvw),maxval(rww)
    
    !Find Lii: Linear stochastic estimate coefficients
    !Solving linear equation by Carmer's method
    
    temp1 = ruu(1,1,yplus)*rvv(1,1,yplus)*rww(1,1,yplus)
    write(*,*)ruu(1,1,yplus),rvv(1,1,yplus),rww(1,1,yplus)
    
    temp2 = rvu(1,1,yplus)*rwv(1,1,yplus)*ruw(1,1,yplus)
    write(*,*)rvu(1,1,yplus),rwv(1,1,yplus),ruw(1,1,yplus)
    
    temp3 = rwu(1,1,yplus)*ruv(1,1,yplus)*rvw(1,1,yplus)
    write(*,*)rwu(1,1,yplus),ruv(1,1,yplus),rvw(1,1,yplus)
    
    temp4 = ruw(1,1,yplus)*rvv(1,1,yplus)*rwu(1,1,yplus)
    write(*,*)ruw(1,1,yplus),rvv(1,1,yplus),rwu(1,1,yplus)
    
    temp5 = rvw(1,1,yplus)*rwv(1,1,yplus)*ruu(1,1,yplus)
    write(*,*)rvw(1,1,yplus),rwv(1,1,yplus),ruu(1,1,yplus)
    
    temp6 = rww(1,1,yplus)*ruv(1,1,yplus)*rvu(1,1,yplus)
    write(*,*)rww(1,1,yplus),ruv(1,1,yplus),rvu(1,1,yplus)
    
    detA = temp1 + temp2 + temp3 - temp4 - temp5 - temp6
    
    a1 = (rvv(1,1,yplus)*rww(1,1,yplus) - rvw(1,1,yplus)*rwv(1,1,yplus)) / detA
    a2 = (rwu(1,1,yplus)*rvw(1,1,yplus) - rvu(1,1,yplus)*rww(1,1,yplus)) / detA
    a3 = (rvu(1,1,yplus)*rwv(1,1,yplus) - rwu(1,1,yplus)*rvv(1,1,yplus)) / detA
    
    b1 = (ruw(1,1,yplus)*rwv(1,1,yplus) - ruv(1,1,yplus)*rww(1,1,yplus)) / detA
    b2 = (ruu(1,1,yplus)*rww(1,1,yplus) - ruw(1,1,yplus)*rwu(1,1,yplus)) / detA
    b3 = (ruv(1,1,yplus)*rwu(1,1,yplus) - ruu(1,1,yplus)*rwv(1,1,yplus)) / detA
    
    c1 = (ruv(1,1,yplus)*rvw(1,1,yplus) - ruw(1,1,yplus)*rvv(1,1,yplus)) / detA
    c2 = (ruw(1,1,yplus)*rvu(1,1,yplus) - ruu(1,1,yplus)*rvw(1,1,yplus)) / detA
    c3 = (ruu(1,1,yplus)*rvv(1,1,yplus) - ruv(1,1,yplus)*rvu(1,1,yplus)) / detA
    
    l11 = ruu*a1 + rvu*a2 + rwu*a3
    l12 = ruu*b1 + rvu*b2 + rwu*b3
    l13 = ruu*c1 + rvu*c2 + rwu*c3
    l21 = ruv*a1 + rvv*a2 + rwv*a3
    l22 = ruv*b1 + rvv*b2 + rwv*b3
    l23 = ruv*c1 + rvv*c2 + rwv*c3
    l31 = ruw*a1 + rvw*a2 + rww*a3
    l32 = ruw*b1 + rvw*b2 + rww*b3
    l33 = ruw*c1 + rvw*c2 + rww*c3
    
    !Write Linear stochastic estimate coefficients
    call save_data(l11,l12,l13,postdata//"lsecoeff_1j")
    call save_data(l21,l22,l23,postdata//"lsecoeff_2j")
    call save_data(l31,l32,l33,postdata//"lsecoeff_3j")
    
  end subroutine lse
  
end module mod_lse_coeff
