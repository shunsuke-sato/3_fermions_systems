module global_variables
  implicit none 
! math parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! Physical constants
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: fs = 1d0/0.024189d0
  real(8),parameter :: bohr = 0.52917721067d0

! Finite difference parameters
  real(8),parameter :: lc2 = -1d0/12d0, lc1 = 4d0/3d0, lc0 = -5d0/2d0
  real(8),parameter :: gc2 = -1d0/12d0, gc1 = 2d0/3d0

  integer :: nx, nt
  real(8) :: dx, dt
  real(8) :: Tprop

! material parameters
  real(8) :: lattice_constant
  real(8) :: bvc_lattice_constant


! laser parameters
  real(8) :: E0, omega, Tpulse. phi_CEP


! grids
  real(8), allocatable :: xn(:)

! wavefunction
  complex(8), allocatable :: zpsi(:,:,:)

! potentials
  real(8), allocatable :: vpot_1d(:), wpot_1d(:)
  real(8), allocatable :: vpot(:,:,:), wpot(:,:,:), tot_pot(:,:,:)


end module global
!-------------------------------------------------------
program main
  use global_variables
  implicit none





end program main
!-------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  integer :: ix
  
  call read_input_parameters

  call set_grids

 


  allocate(zpsi(0:nx-1, 0:nx-1, 0:nx-1))

  allocate(vpot(0:nx-1, 0:nx-1, 0:nx-1))
  allocate(wpot(0:nx-1, 0:nx-1, 0:nx-1))
  allocate(tot_pot(0:nx-1, 0:nx-1, 0:nx-1))


  call set_potentials


end subroutine initialize
!-------------------------------------------------------
subroutine read_input_parameters
  use global_variables
  implicit none
  real(8) :: Tprop_fs
  real(8) :: E0_MVm, omega_ev, Tpulse_fs. phi_CEP_2pi

  read(*,*) lattice_constant, nx
  read(*,*) Tprop_fs, dt
  read(*,*) E0_MVm, omega_ev, Tpulse_fs. phi_CEP_2pi

  Tprop = Tprop_fs*fs
  E0 = E0_MVm*1d6*ev/(bohr*1d10)
  omega = omega_ev*ev
  Tpulse = Tpulse_fs*fs
  phi_CEP = phi_CEP_2pi*2d0*pi


  bvc_lattice_constant = lattice_constant*3


end subroutine read_input_parameters
!-------------------------------------------------------
subroutine set_grids
  use global_variables
  implicit none
  integer :: ix

  allocate(xn(0:nx-1))
  dx = bvc_lattice_constant/nx


  do ix = 0, nx-1
    xn(ix) = ix*dx
  end do

end subroutine set_grids
!-------------------------------------------------------
subroutine set_potentials
  use global_variables
  implicit none
  integer :: ix1, ix2, ix3
  real(8) :: x1, x2, x3
  real(8),parameter :: v0 = 1d0
  real(8),parameter :: w0 = 1d0

  allocate(vpot_1d(0:nx-1))
  allocate(wpot_1d(0:nx-1))

  do ix1 = 0, nx-1
    x1 = xn(ix1)
    vpot_1d(ix1) = v0*( (cos(pi*x1/lattice_constant))**2 &
        + 0.25d0*sin(4d0*pi*x1/lattice_constant) )
  end do


  do ix1 = 0, nx-1
    x1 = xn(ix1)
    wpot_1d(ix1) = w0*cos(pi*x1/bvc_lattice_constant)**16
  end do    

  do ix1 = 1, nx-1
    do ix2 = 0, nx-1
      do ix3 = 0, nx-1
        vpot(ix1, ix2, ix3) = vpot_1d(ix1) + vpot_1d(ix2) + vpot_1d(ix3)
        wpot(ix1, ix2, ix3) = wpot_1d(abs(ix1-ix2)) &
            + wpot_1d(abs(ix2-ix3)) + wpot_1d(abs(ix3-ix1))
      end do
    end do
  end do

  tot_pot = vpot + wpot


end subroutine set_potentials
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
