module constants_mod
  
  implicit none

  logical,parameter :: lgBenchmark = .false. ! Use benchmark data from Camps+2015
  logical,parameter :: lgContcool = .false. ! Treat dust cooling as continuous (from level i to i-1 only)
  logical,parameter :: lgElec = .true. ! Include electron collisional heating
  logical,parameter :: lgHyd = .true. ! Include heating by hydrogen nuclei
  logical,parameter :: lgEquibT = .false. ! Calculate the equilibrium temperature and emission
  logical,parameter :: lgOldRe = .false. ! Use electron range from Dwek 1987
  double precision,parameter :: pi = 3.14159265
  double precision,parameter :: twopi = 2d0*pi
  double precision,parameter :: fourpi = 4d0*pi
  double precision,parameter :: fourpisq = fourpi*pi
  double precision,parameter :: fourthirdpi = fourpi/3d0
  double precision,parameter :: c_s = 2.99792458d10 ! speed of light cm s-1
  double precision,parameter :: h_p = 6.62606957d-27 ! planck constant erg s
  double precision,parameter :: hc = h_p*c_s
  double precision,parameter :: k_b = 1.3806488d-16 ! boltzmann constant erg K-1
  double precision,parameter :: eV = 1.6021766d-12 ! electron volt erg
  double precision,parameter :: m_e = 9.1093836d-28 ! electron mass g
  double precision,parameter :: m_p = 1.6726216d-24 ! proton mass g
  double precision,parameter :: pc = 3.086d18

end module constants_mod
