module precision
  implicit none
  save
  integer, parameter :: ip = selected_int_kind(15)
  integer, parameter :: dp = selected_real_kind(15)
end module precision

module param
  use precision
  implicit none
  integer(ip) :: npartdim, natom, nstep, istep
  real(dp) :: tempK, dt, boxl(3), alat, mass
  real(dp) :: avtemp, ke, kb, epsilon, sigma, scale
  real(dp),dimension(3,4) :: rcell = reshape( (/ &
       0.0D+00, 0.0D+00, 0.0D+00, &
       0.5D+00, 0.5D+00, 0.0D+00, &
       0.0D+00, 0.5D+00, 0.5D+00, &
       0.5D+00, 0.0D+00, 0.5D+00 /), (/ 3, 4 /) )
  character(len=2) :: pot
end module param

module potential
  use precision
  use param, only : natom, boxl
  implicit none
  real(dp) :: r2, r6, d2, d
  real(dp), parameter :: de = 0.176d0, a = 1.4d0, re = 1d0
  real(dp) :: exparre
  
contains
  subroutine lennard_jones(coord, force, pener)
    ! Lennard Jones Potential
    ! V = 4 * epsilon * [ (sigma/r)**12 - (sigma/r)**6 ]
    !   = 4 * epsilon * (sigma/r)**6 * [ (sigma/r)**2 - 1 ]
    !   = 4 * r**(-6) * [ r**(-2) - 1 ] for epsilon=sigma=1
    ! F_i = 48 * epsilon * (sigma/r)**6 * (1/r**2) * [ ( sigma/r)** 6 - 0.5 ] * i where i = x,y,z
    !     = 48 * r**(-8) * [ r**(-6) - 0.5 ] * i  for epsilon=sigma=1    implicit none
    implicit none
    real(dp), dimension(:,:), intent(in) :: coord
    real(dp), dimension(:,:), intent(out) :: force
    real(dp), intent(out) :: pener
    real(dp) :: r(3)
    integer(ip) :: i, j
    
    force = 0d0 ; pener = 0d0
    !$omp parallel do default(shared) private(r,r2,r6) reduction(+:pener)
    do i = 1, natom
       do j = 1, natom
          if ( i == j ) cycle
          r(:) = coord(i,:) - coord(j,:)
          r(:) = r(:) - nint(r(:)/boxl(:)) * boxl(:)
          r2 = 1.d0 / dot_product(r,r)
          r6 = r2 ** 3

          force(i,:) = force(i,:) + dvdr_lj(r2, r6) * r(:)
          pener = pener + pot_lj(r2, r6) / 2.d0
       end do
    end do
    !$omp end parallel do
  end subroutine lennard_jones

  subroutine morse(coord, force, pener)
    ! Morse Potential
    ! V = D * [ 1 - exp(-a*(r - re)) ]^2
    ! F_i = 2*D * [ 1 - exp(-a*(r - re)) ] * a exp(-a*(r-re)) * i / r  
    implicit none
    real(dp), dimension(:,:), intent(in) :: coord
    real(dp), dimension(:,:), intent(out) :: force
    real(dp), intent(out) :: pener
    real(dp) :: r(3)
    integer(ip) :: i, j

    force = 0d0 ; pener = 0d0
    !$omp parallel do default(shared) private(r,d,d2,exparre) reduction(+:pener)
    do i = 1, natom
       do j = 1, natom
          if ( i == j ) cycle
          r(:) = coord(i,:) - coord(j,:)
          r(:) = r(:) - nint(r(:)/boxl(:)) * boxl(:)
          d2 = dot_product(r,r)
          d = sqrt(d2)
          exparre = exp( -a * (d - re ))
          
          force(i,:) = force(i,:) + dvdr_mp(exparre) * r(:)
          pener = pener + pot_mp(exparre) / 2.d0
       end do
    end do
    !$omp end parallel do
  end subroutine morse

  function pot_lj(r2, r6)
    implicit none
    real(dp), intent(in) :: r2, r6
    real(dp) :: pot_lj
    pot_lj = 4d0 * r6 * ( r6 - 1.d0 )
  end function pot_lj
  function pot_mp(exparre)
    implicit none
    real(dp), intent(in) :: exparre
    real(dp) :: pot_mp
    pot_mp = de * ( 1d0 - exparre )**2
  end function pot_mp

  function dvdr_lj(r2,r6)
    implicit none
    real(dp), intent(in) :: r2, r6
    real(dp) :: dvdr_lj
    dvdr_lj = 48d0 * r2 * r6 * ( r6 - 0.5d0 ) 
  end function dvdr_lj
  function dvdr_mp(exparre)
    implicit none
    real(dp), intent(in) :: exparre
    real(dp) :: dvdr_mp
    dvdr_mp = 2d0 * de * a * (1d0 - exparre) * exparre
  end function dvdr_mp
end module potential

program md

  ! Molecular Dynamics code for equilibration of Liquid Argon
  ! Author: Alex Pacheco
  ! Date  : Jan 30, 2014

  ! This program simulates the equilibration of Liquid Argon 
  ! starting from a FCC crystal structure using Lennard-Jones
  ! potential and velocity verlet algorithm
  ! See md-orig.f90 for more details about this program

  ! Disclaimer: 
  ! This is code can be used as an introduction to molecular dynamics. There are lot more
  ! concepts in MD that are not covered here. 

  ! Parameters:
  ! npartdim : number of unit cells, uniform in all directions. change to nonuniform if you desire
  ! natom : number of atoms
  ! nstep : nummber of simulation time steps
  ! tempK : equilibration temperature
  ! dt : simulation time steps
  ! boxl : length of simulation box in all directions
  ! alat : lattice constant for fcc crystal
  ! kb : boltzmann constant, set to 1 for simplicity
  ! mass : mass of Ar atom, set to 1 for simplicity
  ! epsilon, sigma : LJ parameters, set to 1 for simplicity 
  ! rcell : FCC unit cell 
  ! coord, coord_t0 : nuclear positions for each step, current and initial
  ! vel, vel_t0 : nuclear velocities for each step
  ! acc, acc_t0 : nuclear acceleration for each step
  ! force, pener : force and potential energy at current step
  ! avtemp : average temperature at current time step
  ! scale : scaling factor to set current temperature to desired temperature
  ! gasdev : Returns a normally distributed deviate with zero mean and unit variance from Numerical recipes

  use precision
  use param
  implicit none
  integer(ip) :: n, i, j, k, l
  real(dp), dimension(:,:), allocatable :: coord_t0, vel_t0, acc_t0
  real(dp), dimension(:,:), allocatable :: coord, vel, acc, force
  real(dp) :: pener
  ! For timing analysis
  real(dp) :: init_time, start_time, end_time
  integer, dimension(8) :: value

  interface
     subroutine initialize(coord_t0, vel_t0, acc_t0)
       use precision
       implicit none
       real(dp), dimension(:,:), intent(out) :: coord_t0, vel_t0, acc_t0
     end subroutine initialize
     subroutine verlet(coord, coord_t0, vel, vel_t0, acc, acc_t0, force, pener)
       use precision
       implicit none
       real(dp), dimension(:,:), intent(in) :: coord_t0, vel_t0, acc_t0
       real(dp), dimension(:,:), intent(out) :: coord, vel, acc, force
       real(dp), intent(out) :: pener
     end subroutine verlet
     subroutine linearmom(vel)
       use precision
       implicit none
       real(dp), dimension(:,:), intent(in) :: vel
     end subroutine linearmom
     subroutine get_temp(vel)
       use precision
       implicit none
       real(dp), dimension(:,:), intent(in) :: vel
     end subroutine get_temp
  end interface

  ! Namelist for reading input
  namelist/moldyn/natom,npartdim,nstep,tempK,dt,pot

  ! Format statements
1000 format(i8)
1100 format(a,2x,3f12.6)
1200 format(a,2x,e15.8)
1300 format(a,2x,i8,2x,e15.8,1x,e15.8)
  
  alat = 2d0 ** (2d0/3d0)
  kb = 1.d0
  mass = 1.d0
  epsilon = 1.d0
  sigma = 1.d0
  ! Set default values
  natom = 0
  npartdim = 0
  nstep = 1
  tempK = 10d0
  dt = 0.001d0
  pot="lj"

  read(5,nml=moldyn)
  if( natom > 0 ) then
     npartdim = ( natom / 4.d0 ) ** ( 1.d0 / 3.d0 ) 
  else if ( npartdim > 0 ) then
     natom = 4.d0 * npartdim ** 3
  else 
     print *, 'Invalid input: no atom in simulation box'
     stop
  end if
  write(6,'(a)') 'Input Parameters:'
  write(6,nml=moldyn)

  boxl = npartdim * alat

  ! Allocate arrays
  allocate(coord(natom,3), coord_t0(natom,3))
  allocate(vel(natom,3), vel_t0(natom,3))
  allocate(acc(natom,3), acc_t0(natom,3))
  allocate(force(natom,3))

  !=================================================
  ! Initialize coordinates and random velocities
  !=================================================

  call date_and_time(VALUES=value)
  init_time = real(value(5)*3600,dp) + real(value(6)*60,dp) + real(value(7),dp) + real(value(8),dp)/1000d0
  call initialize(coord_t0, vel_t0, acc_t0)

  open(unit=1,file='atom.xyz',status='unknown')
  write(1,1000) natom
  write(1,*)
  do i = 1, natom
     write(1,1100) 'Ar', coord_t0(i,1), coord_t0(i,2), coord_t0(i,3)
  end do
  
  ! Set Linear Momentum to zero
  call linearmom(vel_t0)

  ! get current temperature
  call get_temp(vel_t0)
  print 1200, 'Initial Average Temperature: ', avtemp

  ! scale initial velocity to desired temperature
  scale = sqrt( tempK / avtemp )
  vel_t0 = vel_t0 * scale
  call get_temp(vel_t0)
  print 1200, 'Initial Scaled Average Temperature: ', avtemp

  call date_and_time(VALUES=value)
  start_time = real(value(5)*3600,dp) + real(value(6)*60,dp) + real(value(7),dp) + real(value(8),dp)/1000d0

  !=================================================
  ! MD Simulation
  !=================================================

  do istep = 1, nstep

     ! Get new atom positions from Velocity Verlet Algorithm
     call verlet(coord, coord_t0, vel, vel_t0, acc, acc_t0, force, pener)

     ! Set Linear Momentum to zero
     call linearmom(vel)
     
     ! compute average temperature
     call get_temp(vel)
     print 1300, 'Average Temperature: ' , istep, avtemp, pener

     scale = sqrt ( tempk / avtemp )
     ! Reset for next time step
     coord_t0 = coord
     acc_t0 = acc
     ! scale velocity to desired temperature
     vel_t0 = vel * scale

     ! Write current coordinates to xyz file for visualization
     write(1,1000) natom
     write(1,*)
     do i = 1, natom
        write(1,1100) 'Ar', coord_t0(i,1), coord_t0(i,2), coord_t0(i,3)
     end do
  end do
  close(1)
  deallocate(coord_t0,vel_t0,acc_t0,coord,vel,acc,force)

  call date_and_time(VALUES=value)
  end_time = real(value(5)*3600,dp) + real(value(6)*60,dp) + real(value(7),dp) + real(value(8),dp)/1000d0
  write(6,'(a,f12.3,/,a,f12.3)') 'Init Time: ', start_time - init_time, &
       'Sim  Time: ', end_time - start_time

end program md

subroutine initialize(coord_t0, vel_t0, acc_t0)
  use precision
  use param, only : natom, npartdim, alat, rcell
  implicit none
  real(dp), dimension(:,:), intent(out) :: coord_t0, vel_t0, acc_t0
  integer(ip) :: n, i, j, k, l

  ! Set initial coordinates, velocity and acceleration to zero
  coord_t0 = 0d0 ; vel_t0 = 0d0 ; acc_t0 = 0d0
  
  ! Create FCC unit cell
  rcell = rcell * alat
  
  ! Create a FCC crystal structure
  n = 1
  do i = 1, npartdim
     do j = 1, npartdim
        do k = 1, npartdim
           do l = 1, 4
              coord_t0(n,1) = alat * real(i - 1, dp) + rcell(1,l)
              coord_t0(n,2) = alat * real(j - 1, dp) + rcell(2,l)
              coord_t0(n,3) = alat * real(k - 1, dp) + rcell(3,l)
              n = n + 1
           end do
        end do
     end do
  end do
  
  ! Assign initial random velocities
  do i = 1, natom
     do j = 1, 3
        vel_t0(i,j) = gasdev()
     end do
  end do

contains
  function gasdev()
    use precision
    implicit none
    real(dp) :: gasdev
    real(dp) :: v1, v2, fac, rsq
    real(dp), save :: gset
    logical, save :: available = .false.
    
    if (available) then
       gasdev = gset
       available = .false.
    else
       do
          call random_number(v1)
          call random_number(v2)
          v1 = 2.d0 * v1 - 1.d0
          v2 = 2.d0 * v2 - 1.d0
          rsq = v1**2 + v2**2
          if ( rsq > 0.d0 .and. rsq < 1.d0 ) exit
       end do
       fac = sqrt(-2.d0 * log(rsq) / rsq)
       gasdev = v1 * fac
       gset = v2 * fac
       available = .true.
    end if
  end function gasdev

end subroutine initialize
  
subroutine verlet(coord, coord_t0, vel, vel_t0, acc, acc_t0, force, pener)
  use precision
  use potential, only : lennard_jones, morse
  use param, only : natom, mass, dt, boxl, pot
  implicit none
  real(dp), dimension(:,:), intent(in) :: coord_t0, vel_t0, acc_t0
  real(dp), dimension(:,:), intent(out) :: coord, vel, acc, force
  real(dp), intent(out) :: pener
  integer(ip) :: i
  
  ! Set coordinates, velocity, acceleration and force at next time step to zero
  coord = 0d0 ; vel = 0d0 ; acc = 0d0! ; force = 0d0 
!  pener = 0d0
  
  ! Get new atom positions from Velocity Verlet Algorithm  
  !$omp parallel do default(shared)
  do i = 1, natom
     coord(i,:) = coord_t0(i,:) + vel_t0(i,:) * dt + 0.5d0 * acc_t0(i,:) * dt ** 2
  end do
  !$omp end parallel do

  !$omp parallel do default(shared)
  do i = 1, natom
     ! Apply PBC to coordinates
     where ( coord(i,:) > boxl(:) )
        coord(i,:) = coord(i,:) - boxl(:)
     elsewhere ( coord(i,:) < 0d0 )
        coord(i,:) = coord(i,:) + boxl(:)
     end where
  end do
  !$omp end parallel do

  ! Get force and potential at new atom positions
  select case(pot)
  case('mp')
     call morse(coord, force, pener)
  case default
     call lennard_jones(coord, force, pener)
  end select

  ! Calculate Acceleration and Velocity  at current time step
  !$omp parallel do default(shared)
  do i = 1, natom
     acc(i,:) = force(i,:) / mass
     vel(i,:) = vel_t0(i,:) + 0.5d0 * ( acc(i,:) + acc_t0(i,:) ) * dt
  end do
  !$omp end parallel do
  
end subroutine verlet

subroutine linearmom(vel)
  use precision
  use param, only : natom
  implicit none
  real(dp), dimension(:,:), intent(inout) :: vel
  integer(ip) :: i, j
  real(dp) :: vcm(3)

  ! First get center of mass velocity
  vcm = 0d0
  !$omp parallel do default(shared) reduction(+:vcm)
  do i = 1, 3
     do j = 1, natom
        vcm(i) = vcm(i) + vel(j,i)
     end do
  end do
  !$omp end parallel do
  vcm = vcm / real(natom,dp)
  
  ! Now remove center of mass velocity from all atoms
  !$omp parallel do default(shared)
  do i = 1, natom
     vel(i,:) = vel(i,:) - vcm(:)
  end do
  !$omp end parallel do

end subroutine linearmom

subroutine get_temp(vel)
  use precision
  use param, only : natom, avtemp, mass, kb
  implicit none
  real(dp), dimension(:,:), intent(in) :: vel
  integer(ip) :: i
  real(dp) :: ke

  ke = 0d0
  !$omp parallel do default(shared) reduction(+:ke)
  do i = 1, natom
     ke = ke + dot_product(vel(i,:),vel(i,:))
  end do
  !$omp end parallel do
  avtemp = mass * ke / ( 3d0 * kb * real( natom - 1, dp))
  
end subroutine get_temp

