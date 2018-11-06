subroutine grainemis(nEnthalpy,Tmid,P,a_g,N_a,nWav,lambda,Q_a,Jgrain)
  use constants_mod

  implicit none

  integer,intent(in) :: nEnthalpy,nWav
  double precision,intent(in) :: Tmid(nEnthalpy),P(nEnthalpy),a_g,N_a
  double precision,intent(in) :: lambda(nWav),Q_a(nWav)
  double precision,intent(out) :: Jgrain(nWav)
  double precision :: emis,Jbb
  integer :: i,j

  do i=1,nWav
     Jgrain(i) = 0.0d0
     do j=1,nEnthalpy
        emis = pi*(1d-4*a_g)**2*P(j)*Q_a(i)*Jbb(lambda(i),Tmid(j))
        Jgrain(i) = Jgrain(i) + emis
     end do
     Jgrain(i) = Jgrain(i)*N_a
     if (Jgrain(i) .lt. 1.0d-99) Jgrain(i) = 0.0d0
  end do

end subroutine grainemis
