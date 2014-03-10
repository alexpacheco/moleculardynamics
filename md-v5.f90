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

module dynamic_data
  use precision
  implicit none

  type dynamics
     real(dp) :: x, y, z
  end type dynamics

  interface dot_product
     module procedure dprod
  end interface dot_product
  interface operator (*)
     module procedure tmulr, rmult, tmulr_3d, rmult_3d
  end interface operator (*)
  interface operator (/)
     module procedure tdivr, tdivr_3d
  end interface operator (/)
  interface operator (+)
     module procedure taddt, taddt_3d
  end interface operator (+)
  interface operator (-)
     module procedure tsubt, tsubt_3d
  end interface operator (-)
  interface assignment (=)
     module procedure assigntr, assigntr_3d
  end interface assignment (=)

contains
  function dprod(a,b) result(c)
    implicit none
    type(dynamics), intent(in) :: a, b
    real(dp) :: c
    c = a%x * b%x + a%y * b%y + a%z * b%z
  end function dprod
  function tmulr(a,b) result(c)
    implicit none
    type(dynamics), intent(in) :: a
    real(dp), intent(in) :: b
    type(dynamics) :: c
    c%x = a%x * b ; c%y = a%y * b ; c%z = a%z * b
  end function tmulr
  function rmult(b,a) result(c)
    implicit none
    type(dynamics), intent(in) :: a
    real(dp), intent(in) :: b
    type(dynamics) :: c
    c%x = b * a%x ; c%y = b * a%y ; c%z = b * a%z
  end function rmult
  function tmulr_3d(a,b) result(c)
    implicit none
    type(dynamics), dimension(:), intent(in) :: a
    real(dp), intent(in) :: b
    type(dynamics), dimension(size(a)) :: c
    c%x = a%x * b ; c%y = a%y * b ; c%z = a%z * b
  end function tmulr_3d
  function rmult_3d(b,a) result(c)
    implicit none
    type(dynamics), dimension(:), intent(in) :: a
    real(dp), intent(in) :: b
    type(dynamics), dimension(size(a)) :: c
    c%x = b * a%x ; c%y = b * a%y ; c%z = b * a%z
  end function rmult_3d
  function tdivr(a,b) result(c)
    implicit none
    type(dynamics), intent(in) :: a
    real(dp), intent(in) :: b
    type(dynamics) :: c
    c%x = a%x / b ; c%y = a%y / b ; c%z = a%z / b
  end function tdivr
  function tdivr_3d(a,b) result(c)
    implicit none
    type(dynamics), dimension(:), intent(in) :: a
    real(dp), intent(in) :: b
    type(dynamics), dimension(size(a)) :: c
    c%x = a%x / b ; c%y = a%y / b ; c%z = a%z / b
  end function tdivr_3d
  function taddt(a,b) result(c)
    implicit none
    type(dynamics), intent(in) :: a, b
    type(dynamics) :: c
    c%x = a%x + b%x ; c%y = a%y + b%y ; c%z = a%z + b%z
  end function taddt
  function taddt_3d(a,b) result(c)
    implicit none
    type(dynamics), dimension(:), intent(in) :: a, b
    type(dynamics),dimension(size(a)) :: c
    c%x = a%x + b%x ; c%y = a%y + b%y ; c%z = a%z + b%z
  end function taddt_3d
  function tsubt(a,b) result(c)
    implicit none
    type(dynamics), intent(in) :: a, b
    type(dynamics) :: c
    c%x = a%x - b%x ; c%y = a%y - b%y ; c%z = a%z - b%z
  end function tsubt
  function tsubt_3d(a,b) result(c)
    implicit none
    type(dynamics), dimension(:), intent(in) :: a, b
    type(dynamics),dimension(size(a)) :: c
    c%x = a%x - b%x ; c%y = a%y - b%y ; c%z = a%z - b%z
  end function tsubt_3d
  subroutine assigntr(a,b)
    implicit none
    type(dynamics), intent(out) :: a
    real(dp), intent(in) :: b
    a%x = b ; a%y = b ; a%z = b
  end subroutine assigntr
  subroutine assigntr_3d(a,b)
    implicit none
    type(dynamics), dimension(:), intent(out) :: a
    real(dp), intent(in) :: b
    a%x = b ; a%y = b ; a%z = b
  end subroutine assigntr_3d

end module dynamic_data

module potential
  use dynamic_data
  implicit none
  real(dp) :: r2, r6, d2, d
  real(dp), parameter :: de = 0.176d0, a = 1.4d0, re = 1d0
  real(dp) :: exparre
  
contains
  subroutine lennard_jones(r,f,p)
    ! Lennard Jones Potential
    ! V = 4 * epsilon * [ (sigma/r)**12 - (sigma/r)**6 ]
    !   = 4 * epsilon * (sigma/r)**6 * [ (sigma/r)**2 - 1 ]
    !   = 4 * r**(-6) * [ r**(-2) - 1 ] for epsilon=sigma=1
    ! F_i = 48 * epsilon * (sigma/r)**6 * (1/r**2) * [ ( sigma/r)** 6 - 0.5 ] * i where i = x,y,z
    !     = 48 * r**(-8) * [ r**(-6) - 0.5 ] * i  for epsilon=sigma=1    implicit none
    implicit none
    type(dynamics), intent(in) :: r
    type(dynamics), intent(out) :: f
    real(dp), intent(out) :: p

    r2 = 1.d0 / dot_product(r,r)
    r6 = r2 ** 3

    f = dvdr_lj(r2, r6) * r
    p = pot_lj(r2, r6)
  end subroutine lennard_jones

  subroutine morse(r,f,p)
    ! Morse Potential
    ! V = D * [ 1 - exp(-a*(r - re)) ]^2
    ! F_i = 2*D * [ 1 - exp(-a*(r - re)) ] * a exp(-a*(r-re)) * i / r  
    implicit none
    type(dynamics), intent(in) :: r
    type(dynamics), intent(out) :: f
    real(dp), intent(out) :: p

    d2 = dot_product(r,r)
    d = sqrt(d2)
    exparre = exp( -a * (d - re ))
    
    f = dvdr_mp(exparre) * r
    p = pot_mp(exparre)
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
  use dynamic_data
  implicit none
  integer(ip) :: n, i, j, k, l
  type(dynamics), dimension(:), allocatable :: coord_t0, vel_t0, acc_t0
  type(dynamics), dimension(:), allocatable :: coord, vel, acc, force
  real(dp) :: pener
  ! For timing analysis
  real(dp) :: init_time, start_time, end_time
  integer, dimension(8) :: value

  interface
     subroutine initialize(coord_t0, vel_t0, acc_t0)
       use dynamic_data
       implicit none
       type(dynamics), dimension(:), intent(out) :: coord_t0, vel_t0, acc_t0
     end subroutine initialize
     subroutine verlet(coord, coord_t0, vel, vel_t0, acc, acc_t0, force, pener)
       use dynamic_data
       implicit none
       type(dynamics), dimension(:), intent(in) :: coord_t0, vel_t0, acc_t0
       type(dynamics), dimension(:), intent(out) :: coord, vel, acc, force
       real(dp), intent(out) :: pener
     end subroutine verlet
     subroutine linearmom(vel)
       use dynamic_data
       implicit none
       type(dynamics), dimension(:), intent(in) :: vel
     end subroutine linearmom
     subroutine get_temp(vel)
       use dynamic_data
       implicit none
       type(dynamics), dimension(:), intent(in) :: vel
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
  allocate(coord(natom), coord_t0(natom))
  allocate(vel(natom), vel_t0(natom))
  allocate(acc(natom), acc_t0(natom))
  allocate(force(natom))

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
     write(1,1100) 'Ar', coord_t0(i)%x, coord_t0(i)%y, coord_t0(i)%z
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
        write(1,1100) 'Ar', coord_t0(i)%x, coord_t0(i)%y, coord_t0(i)%z
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
  use dynamic_data
  use param, only : natom, npartdim, alat, rcell
  implicit none
  type(dynamics), dimension(:), intent(out) :: coord_t0, vel_t0, acc_t0
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
              coord_t0(n)%x = alat * real(i - 1, dp) + rcell(1,l)
              coord_t0(n)%y = alat * real(j - 1, dp) + rcell(2,l)
              coord_t0(n)%z = alat * real(k - 1, dp) + rcell(3,l)
              n = n + 1
           end do
        end do
     end do
  end do
  
  ! Assign initial random velocities
  do i = 1, natom
     vel_t0(i) = dynamics( gasdev(), gasdev(), gasdev() )
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
  use dynamic_data
  use potential
  use param, only : natom, mass, dt, boxl, pot
  implicit none
  type(dynamics), dimension(:), intent(in) :: coord_t0, vel_t0, acc_t0
  type(dynamics), dimension(:), intent(out) :: coord, vel, acc, force
  real(dp), intent(out) :: pener
  integer(ip) :: i, j, k
  real(dp) :: epot
  type(dynamics) :: r, f

  ! Set coordinates, velocity, acceleration and force at next time step to zero
  coord = 0d0 ; vel = 0d0 ; acc = 0d0 ; force = 0d0 
  pener = 0d0
  r = 0d0 
  
  ! Get new atom positions from Velocity Verlet Algorithm  
  coord = coord_t0 + vel_t0 * dt + 0.5d0 * acc_t0 * dt ** 2
  do i = 1, natom
     ! Apply PBC to coordinates
     coord(i)%x = pbc(coord(i)%x, boxl(1))
     coord(i)%y = pbc(coord(i)%y, boxl(2))
     coord(i)%z = pbc(coord(i)%z, boxl(3))
  end do
  
  ! Get force at new atom positions
  do i = 1, natom - 1
     do j = i + 1, natom
        r = coord(i) - coord(j)
        ! minimum image criterion
        r = nearint(r, boxl)
        select case(pot)
        case('mp')
           call morse( r, f, epot )
        case default
           call lennard_jones( r, f, epot )
        end select
        pener = pener + epot
        force(i) = force(i) + f
        force(j) = force(j) - f
     end do
  end do
  
  ! Calculate Acceleration and Velocity  at current time step
  acc = force / mass
  vel = vel_t0 + 0.5d0 * ( acc + acc_t0 ) * dt

contains
  function pbc(a,b)
    real(dp), intent(in) :: a, b
    real(dp) :: pbc

    if ( a > b ) then
       pbc = a - b
    else if ( a < 0d0 ) then
       pbc = a + b
    else
       pbc = a
    end if
  end function pbc
  function nearint(a,b) result(c)
    type(dynamics), intent(in) :: a
    real(dp), intent(in) :: b(3)
    type(dynamics) :: c
    c%x = a%x - nint(a%x/b(1)) * b(1)
    c%y = a%y - nint(a%y/b(2)) * b(2)
    c%z = a%z - nint(a%z/b(3)) * b(3)
  end function nearint
end subroutine verlet

subroutine linearmom(vel)
  use dynamic_data
  use param, only : natom
  implicit none
  type(dynamics), dimension(:), intent(inout) :: vel
  integer(ip) :: i
  type(dynamics) :: vcm

  ! First get center of mass velocity
  vcm = 0d0
  do i = 1, natom
     vcm = vcm + vel(i)
  end do
  vcm = vcm / real(natom,dp)
  
  ! Now remove center of mass velocity from all atoms
  do i = 1, natom
     vel(i) = vel(i) - vcm
  end do

end subroutine linearmom

subroutine get_temp(vel)
  use dynamic_data
  use param, only : natom, avtemp, mass, kb
  implicit none
  type(dynamics), dimension(:), intent(in) :: vel
  integer(ip) :: i
  real(dp) :: ke

  ke = 0d0
  do i = 1, natom
     ke = ke + dot_product(vel(i),vel(i))
  end do
  avtemp = mass * ke / ( 3d0 * kb * real( natom - 1, dp))
  
end subroutine get_temp

