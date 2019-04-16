! DINAMO v1.02
program dinamo
  use constants_mod
  use particle_mod

  implicit none

  character(len=50) :: runlabel
  ! start/end times
  real :: starttime,endtime
  ! grain species variables
  integer :: nSpecies
  character(len=1),allocatable :: graintype(:)
  character(len=100),allocatable :: ofile(:),tfile(:),sfile(:),gfile(:)
  ! enthalpy bin variables
  integer,parameter :: nEnthalpy = 500
  double precision :: P(nEnthalpy),Tmid(nEnthalpy)
  ! radiation field properties
  integer :: nWav
  double precision,allocatable :: lambda(:),Jrad(:),Jgrain(:),Jtot(:)
  double precision :: JMat,Jbb,Jpl,Jsync,Jcrab
  character(len=20) :: radtype
  double precision :: radc1,radc2,T_cmb
  ! grain radiative properties
  integer :: nQsize
  double precision,allocatable :: Qdata(:,:),Qsize(:)
  character(len=100) :: nkfile
  double precision,allocatable :: nrad(:),krad(:)
  ! grain properties
  integer :: nSizes
  double precision :: rho_g
  double precision :: amin,amax,nmrn,N_g,Tsub
  double precision,allocatable :: agrid(:),f_a(:),a_g(:),N_a(:)
  double precision,allocatable :: Q_a(:)
  integer :: nTemp
  double precision,allocatable :: T_H(:),H_T(:)
  double precision,allocatable :: H_a(:)
  integer :: i,j,s

  write(*,"('Beginning dust stuff')") 
  call cpu_time(starttime)

  if (lgBenchmark) then
     runlabel = 'benchmark'
     nSpecies = 3
     allocate(ofile(nSpecies))
     allocate(tfile(nSpecies))
     allocate(sfile(nSpecies))
     allocate(graintype(nSpecies))
     ofile(1) = 'input/suvSil_121_1201.dat'
     tfile(1) = 'input/Silicate_Calorimetry_1000.dat'
     sfile(1) = 'input/ZDA_BARE_GR_S_SzDist_Sil.dat'
     ofile(2) = 'input/Gra_121_1201.dat'
     tfile(2) = 'input/Graphitic_Calorimetry_1000.dat'
     sfile(2) = 'input/ZDA_BARE_GR_S_SzDist_Gra.dat'
     ofile(3) = 'input/PAH_28_1201_neu.dat'
     tfile(3) = 'input/Graphitic_Calorimetry_1000.dat'
     sfile(3) = 'input/ZDA_BARE_GR_S_SzDist_PAH.dat'
     graintype = (/ 'S', 'C', 'C' /)
     Nelec = 0.0d0
     N_H = 0.0d0
     Telec = 1.0d0
     T_atom = 1.0d0
     atomtype = 'H'
     radtype = 'mathis'
     radc1 = 1.0d0
     radc2 = 1.0d0
     T_cmb = 2.73d0
     Tsub = 2.5d3
  else
     call readinput(runlabel,nSpecies,Nelec,N_H,Telec,T_atom,atomtype,radtype,radc1,radc2,T_cmb)
     allocate(gfile(nSpecies))
     allocate(graintype(nSpecies))
     call readgrainfiles(nSpecies,gfile)
  end if

  if (lgElec) call createparticle(nEner,Telec,Eelec,Pelec)
  if (lgHyd) call createparticle(nEner,T_atom,E_H,P_H)

  call openoutput(nSpecies,runlabel)
  call writeparams(12,radtype,radc1,radc2,Nelec,N_H,Telec,T_atom,atomtype,nSpecies,T_cmb)

  do s=1,nSpecies

     write(*,"('Species ',I1,' of ',I1)") s,nSpecies

     if (lgBenchmark) then
        call readQparams(ofile(s),nQsize,nWav)
        allocate(lambda(nWav))
        allocate(Qdata(nQsize,nWav))
        allocate(Qsize(nQsize))
        call readQ(ofile(s),nQsize,nWav,Qdata,Qsize,lambda,rho_g)
     else
        call readprop(gfile(s),nkfile,rho_g,graintype(s),nSizes,amin,amax,nmrn,N_g,Tsub)
        call countlines('input/'//trim(nkfile),nWav)
        allocate(lambda(nWav))
        allocate(nrad(nWav))
        allocate(krad(nWav))
        call readnk(nkfile,nWav,lambda,nrad,krad)
        call writegrainparams(12,nkfile,rho_g,graintype(s),nSizes,amin,amax,nmrn,N_g,Tsub)
     end if

     allocate(Jrad(nWav))
     allocate(Jgrain(nWav))

     if (radtype .eq. 'mathis') then
        do i=1,nWav
           Jrad(i) = radc1*JMat(lambda(i))
        end do
     else if (radtype .eq. 'blackbody') then
        do i=1,nWav
           Jrad(i) = radc1*Jbb(lambda(i),radc2)
        end do
     else if (radtype .eq. 'pl') then
        do i=1,nWav
           Jrad(i) = Jpl(lambda(i),radc1,radc2)
        end do
     else if (radtype .eq. 'sync') then
        do i=1,nWav
           Jrad(i) = Jsync(lambda(i),radc1,radc2)
        end do
     else if (radtype .eq. 'crab') then
        do i=1,nWav
           Jrad(i) = Jcrab(lambda(i),radc1,radc2)
        end do
     else
        Jrad = 0.0d0
     end if

     do i=1,nWav
        Jrad(i) = Jrad(i) + Jbb(lambda(i),T_cmb)
     end do

     if (s .eq. 1) then
        allocate(Jtot(nWav))
        Jtot = 0.0d0
     end if

     if (lgBenchmark) then
        call readnTemp(tfile(s),nTemp)
        allocate(T_H(nTemp))
        allocate(H_T(nTemp))
        call readHT(tfile(s),nTemp,T_H,H_T)
        call readnSizes(sfile(s),nSizes)
        allocate(agrid(0:nSizes))
        allocate(f_a(0:nSizes))
        call readSizes(sfile(s),nSizes,agrid,f_a)
        allocate(a_g(nSizes))
        allocate(N_a(nSizes))
        call createsizegrid(nSizes,agrid,f_a,a_g,N_a)
     else
        allocate(agrid(0:nSizes))
        allocate(f_a(0:nSizes))
        call mrndist(nSizes,amin,amax,nmrn,agrid,f_a,N_g)
        allocate(a_g(nSizes))
        allocate(N_a(nSizes))
        call createsizegrid(nSizes,agrid,f_a,a_g,N_a)
        nTemp = 1000
        allocate(T_H(nTemp))
        do i=0,nTemp-1
           T_H(i+1) = Tsub**(dble(i)/dble(nTemp-1))
        end do
     end if

     allocate(Q_a(nWav))
     allocate(H_a(nTemp))

     write(10,'(I1,I3)') s,nSizes
     write(11,'(I1,I3)') s,nSizes

     do i=1,nSizes
        write(*,"('Size ',I3,' of ',I3)") i,nSizes
        ! Solve temperature distribution
        if (lgBenchmark) then
           call createQ_a(nWav,a_g(i),Q_a,nQsize,Qsize,Qdata)
           call createH_a(nTemp,H_T,H_a,a_g(i),rho_g)
        else
           call calcmie(nWav,lambda,nrad,krad,a_g(i),Q_a)
           call calcenthalpy(nTemp,T_H,H_a,a_g(i),rho_g,graintype(s))
        end if
        call solveTdist(nEnthalpy,P,nTemp,H_a,T_H,Tmid,nWav,lambda,Q_a,Jrad,a_g(i),rho_g,graintype(s))
        call writeTdist(10,nEnthalpy,a_g(i),Tmid,P)
        ! Calculate emission
        call grainemis(nEnthalpy,Tmid,P,a_g(i),N_a(i),nWav,lambda,Q_a,Jgrain)
        call writeemis(11,nWav,a_g(i),Jgrain)
        do j=1,nWav
           Jtot(j) = Jtot(j) + Jgrain(j)
        end do
     end do

     if (s .eq. nSpecies) then
        write(11,"('Wavelength / um')")
        write(11,'(10000ES10.3)') (lambda(j),j=1,nWav)
        write(11,"('Total dust emission / erg s-1 cm-3 sr-1 cm-1')")
        write(11,'(10000ES10.3)') (Jtot(j),j=1,nWav)
        deallocate(Jtot)
     end if

     if (lgBenchmark) then
        deallocate(lambda)
        deallocate(Qdata)
        deallocate(Qsize)
        deallocate(Jrad)
        deallocate(Jgrain)
        deallocate(T_H)
        deallocate(H_T)
        deallocate(agrid)
        deallocate(f_a)
        deallocate(a_g)
        deallocate(N_a)
        deallocate(Q_a)
        deallocate(H_a)
     else
        deallocate(lambda)
        deallocate(Jrad)
        deallocate(Jgrain)
        deallocate(nrad)
        deallocate(krad)
        deallocate(agrid)
        deallocate(f_a)
        deallocate(a_g)
        deallocate(N_a)
        deallocate(T_H)
        deallocate(Q_a)
        deallocate(H_a)
     end if

     write(*,*)

  end do

  call closeoutput()

  if (lgBenchmark) then
     deallocate(ofile)
     deallocate(tfile)
     deallocate(sfile)
     deallocate(graintype)
  else
     deallocate(gfile)
     deallocate(graintype)
  end if

  call cpu_time(endtime)
  write(*,"('Finished dust stuff')") 
  write(*,"('Time: ',F8.2,' seconds')") endtime-starttime

end program dinamo
