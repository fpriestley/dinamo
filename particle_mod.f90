module particle_mod

  implicit none

  integer,parameter :: nEner = 1000
  double precision,parameter :: agra(0:4) = (/ -8.1070d0, 1.0596d0, -0.27838d0, 0.11741d0, -0.010731d0 /)
  double precision,parameter :: asil(0:4) = (/ -8.1245d0, 1.1102d0, -0.31900d0, 0.12908d0, -0.011757d0 /)
  double precision :: Nelec,Telec,N_H,T_atom
  character(len=2) :: atomtype
  double precision :: Eelec(nEner),Pelec(nEner),EzetaE(nEner)
  double precision :: E_H(nEner),P_H(nEner),EzetaE_H(nEner)
  double precision :: heatinfo(5)

end module particle_mod
