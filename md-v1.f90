module precision
  implicit none
  save
  integer, parameter :: ip = selected_int_kind(15)
  integer, parameter :: dp = selected_real_kind(15)
end module precision

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
  implicit none
  integer(ip), parameter :: npartdim = 10 
  integer(ip), parameter :: natom = 4.d0 * npartdim ** 3
  integer(ip), parameter :: nstep = 1000
  real(dp), parameter :: tempK = 10, dt = 1d-3
  integer(ip) :: istep
  real(dp) :: boxl(3), alat
  integer(ip) :: n, i, j, k, l

  real(dp) :: coord_t0(natom,3), coord(natom,3)
  real(dp) :: vel_t0(natom,3), vel(natom,3)
  real(dp) :: acc_t0(natom,3), acc(natom,3)
  real(dp) :: force(natom,3), pener, mass

  real(dp) :: vcm(3), r(3), rr, r2, r6, f(3)
  real(dp) :: avtemp, ke, kb, epsilon, sigma, scale
  real(dp) :: gasdev
  ! For timing analysis
  real(dp) :: init_time, start_time, end_time
  integer, dimension(8) :: value

  ! Format statements
1000 format(i8)
1100 format(a,2x,3f12.6)
1200 format(a,2x,e15.8)
1300 format(a,2x,i8,2x,e15.8,1x,e15.8)
  
  alat = 2d0 ** (2d0/3d0)
  boxl = npartdim * alat
  kb = 1.d0
  mass = 1.d0
  epsilon = 1.d0
  sigma = 1.d0

  !=================================================
  ! Initialize coordinates and random velocities
  !=================================================

  call date_and_time(VALUES=value)
  init_time = real(value(5)*3600,dp) + real(value(6)*60,dp) + real(value(7),dp) + real(value(8),dp)/1000d0
  ! Set initial coordinates, velocity and acceleration to zero
  coord_t0 = 0d0 ; vel_t0 = 0d0 ; acc_t0 = 0d0
  call initialize(natom, npartdim, alat, coord_t0, vel_t0)

  open(unit=1,file='atom.xyz',status='unknown')
  write(1,1000) natom
  write(1,*)
  do i = 1, natom
     write(1,1100) 'Ar', coord_t0(i,1), coord_t0(i,2), coord_t0(i,3)
  end do
  
  ! Set Linear Momentum to zero
  call linearmom(natom, vel_t0)

  ! get current temperature
  call get_temp(natom, vel_t0, avtemp, mass, kb)
  print 1200, 'Initial Average Temperature: ', avtemp

  ! scale initial velocity to desired temperature
  scale = sqrt( tempK / avtemp )
  vel_t0 = vel_t0 * scale
  call get_temp(natom, vel_t0, avtemp, mass, kb)
  print 1200, 'Initial Scaled Average Temperature: ', avtemp

  call date_and_time(VALUES=value)
  start_time = real(value(5)*3600,dp) + real(value(6)*60,dp) + real(value(7),dp) + real(value(8),dp)/1000d0

  !=================================================
  ! MD Simulation
  !=================================================

  do istep = 1, nstep

     ! Get new atom positions from Velocity Verlet Algorithm
     call verlet(coord, coord_t0, vel, vel_t0, acc, acc_t0, force, pener, natom, mass, dt, boxl)

     ! Set Linear Momentum to zero
     call linearmom(natom, vel)
     
     ! compute average temperature
     call get_temp(natom, vel, avtemp, mass, kb)
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

  call date_and_time(VALUES=value)
  end_time = real(value(5)*3600,dp) + real(value(6)*60,dp) + real(value(7),dp) + real(value(8),dp)/1000d0
  write(6,'(a,f12.3,/,a,f12.3)') 'Init Time: ', start_time - init_time, &
       'Sim  Time: ', end_time - start_time

end program md

subroutine initialize(natom, npartdim, alat, coord_t0, vel_t0)
  use precision
  implicit none
  integer(ip) :: natom, npartdim
  real(dp), dimension(natom,3) :: coord_t0, vel_t0
  integer(ip) :: n, i, j, k, l
  real(dp) :: alat, rcell(3,4)
  
  ! Create FCC unit cell
  rcell = 0d0
  rcell(1,2) = 0.5d0 * alat
  rcell(2,2) = 0.5d0 * alat
  rcell(2,3) = 0.5d0 * alat
  rcell(3,3) = 0.5d0 * alat
  rcell(1,4) = 0.5d0 * alat
  rcell(3,4) = 0.5d0 * alat
  
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
  
subroutine verlet(coord, coord_t0, vel, vel_t0, acc, acc_t0, force, pener, natom, mass, dt, boxl)
  use precision
  implicit none
  real(dp), dimension(natom,3) :: coord, vel, acc, force
  real(dp), dimension(natom,3) :: coord_t0, vel_t0, acc_t0
  real(dp) :: pener
  integer(ip) :: natom, i, j, k
  real(dp) :: mass, dt, boxl(3)
  real(dp) :: r(3), f(3)
  real(dp) :: rr, r2, r6

  ! Set coordinates, velocity, acceleration and force at next time step to zero
  coord = 0d0 ; vel = 0d0 ; acc = 0d0 ; force = 0d0 
  pener = 0d0
  
  ! Get new atom positions from Velocity Verlet Algorithm  
  do i = 1, natom
     coord(i,:) = coord_t0(i,:) + vel_t0(i,:) * dt + 0.5d0 * acc_t0(i,:) * dt ** 2
     ! Apply PBC to coordinates
     do j = 1, 3
        if ( coord(i,j) > boxl(j) ) then
           coord(i,j) = coord(i,j) - boxl(j)
        else if ( coord(i,j) < 0d0 ) then
           coord(i,j) = coord(i,j) + boxl(j)
        endif
     end do
  end do
  
  ! Get force at new atom positions
  ! Using Lennard Jones Potential
  ! Hint: you might want to also seperate the potential and force calculation into a separate subroutine
  ! this will be useful if you want to use other potentials
  
  do i = 1, natom - 1
     do j = i + 1, natom
        r(:) = coord(i,:) - coord(j,:)
        ! minimum image criterion
        r = r - nint( r / boxl ) * boxl
        ! Hint: Use dot_product
        rr = r(1) ** 2 + r(2) ** 2 + r(3) ** 2
        r2 = 1.d0 / rr
        r6 = r2 ** 3
        ! Lennard Jones Potential
        ! V = 4 * epsilon * [ (sigma/r)**12 - (sigma/r)**6 ]
        !   = 4 * epsilon * (sigma/r)**6 * [ (sigma/r)**2 - 1 ]
        !   = 4 * r**(-6) * [ r**(-2) - 1 ] for epsilon=sigma=1
        ! F_i = 48 * epsilon * (sigma/r)**6 * (1/r**2) * [ ( sigma/r)** 6 - 0.5 ] * i where i = x,y,z
        !     = 48 * r**(-8) * [ r**(-6) - 0.5 ] * i  for epsilon=sigma=1
        pener = pener + 4d0 * r6 * ( r6 - 1.d0 )
        f = 48d0 * r2 * r6 * ( r6 - 0.5d0 ) * r
        force(i,:) = force(i,:) + f(:)
        force(j,:) = force(j,:) - f(:)
     end do
  end do
  
  ! Calculate Acceleration and Velocity  at current time step
  acc = force / mass
  vel = vel_t0 + 0.5d0 * ( acc + acc_t0 ) * dt
  
end subroutine verlet

subroutine linearmom(natom, vel)
  use precision
  implicit none
  integer(ip) :: natom
  real(dp) :: vel(natom,3)
  integer(ip) :: i
  real(dp) :: vcm(3)

  ! First get center of mass velocity
  vcm = 0d0
  do i = 1, 3
     vcm(i) = sum(vel(:,i))
  end do
  vcm = vcm / real(natom,dp)
  
  ! Now remove center of mass velocity from all atoms
  do i = 1, natom
     vel(i,:) = vel(i,:) - vcm(:)
  end do

end subroutine linearmom

subroutine get_temp(natom, vel, avtemp, mass, kb)
  use precision
  implicit none
  integer(ip) :: natom
  real(dp) :: vel(natom,3)
  real(dp) :: avtemp
  real(dp) :: mass, kb
  integer(ip) :: i, j
  real(dp) :: ke

  ke = 0d0
  do i = 1, natom
     ke = ke + dot_product(vel(i,:),vel(i,:))
  end do
  avtemp = mass * ke / ( 3d0 * kb * real( natom - 1, dp))
  
end subroutine get_temp

