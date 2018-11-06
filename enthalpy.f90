subroutine enthalpy(T,a,H,type,rho)
  use constants_mod

  implicit none

  character(len=1),intent(in) :: type
  double precision,intent(in) :: T,a,rho
  double precision,intent(out) :: H
  double precision :: vol,c1,c2,c3,natom

  ! Guhathakurta & Draine 1989
  if (type .eq. 'C') then
     vol = fourthirdpi*(1d-4*a)**3 ! in cm3
     natom = rho*vol/(12d0*m_p)
     H = 4.15d-22 * T**3.3
     H = H/(1d0 + 6.51d-3*T + 1.5d-6*T**2 + 8.3d-7*T**2.3)
     H = H*natom*(1d0-2d0/natom)
  else if (type .eq. 'S') then
     vol = fourthirdpi*(1d-4*a)**3 ! in cm3
     natom = rho*vol/(172d0/7d0*m_p) ! assume Fe:Mg:Si:O = 1:1:1:4
     if (T .le. 50.0d0) then
        H = 1.4d3*T**3/3d0
     else if (T .le. 150d0) then
        c1 = 1.4d3*50d0**3/3d0
        H = c1 + 2.2d4*(T**2.3 - 50d0**2.3)/2.3d0
     else if (T .le. 500d0) then
        c1 = 1.4d3*50d0**3/3d0
        c2 = 2.2d4*(150d0**2.3 - 50d0**2.3)/2.3d0
        H = c1 + c2 + 4.8d5*(T**1.68 - 150d0**1.68)/1.68d0
     else
        c1 = 1.4d3*50d0**3/3d0
        c2 = 2.2d4*(150d0**2.3 - 50d0**2.3)/2.3d0
        c3 = 4.8d5*(500d0**1.68 - 150d0**1.68)/1.68d0
        H = c1 + c2 + c3 + 3.41d7*(T - 500d0)
     end if
     H = H*vol*(1d0-2d0/natom)
  else
     H = 0.0d0
  end if

end subroutine enthalpy

subroutine calcenthalpy(nTemp,T,H,a,rho,type)

  implicit none

  integer,intent(in) :: nTemp
  character(len=1),intent(in) :: type
  double precision,intent(in) :: T(nTemp),a,rho
  double precision,intent(out) :: H(nTemp)
  integer :: i

  do i=1,nTemp
     call enthalpy(T(i),a,H(i),type,rho)
  end do

end subroutine calcenthalpy
