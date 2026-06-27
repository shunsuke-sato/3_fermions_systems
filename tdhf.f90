module global_variables
  implicit none
! Build:
!   gfortran -std=f95 -O2 tdhf.f90 -o tdhf
! Run:
!   ./tdhf < input_tdse
! Output:
!   ground_state_hf.log : imaginary-time HF iteration, energy, residual, norms/overlaps
!   current_tdhf.out    : t, A(t), j(t), particle_number, energy, max_norm_error, max_overlap
! math parameters
  real(8),parameter :: pi = 3.141592653589793238462643383279502884197d0
  complex(8),parameter :: zi = (0d0, 1d0)

! Physical constants (atomic units)
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: fs = 1d0/0.024189d0
  real(8),parameter :: bohr = 0.52917721067d0

! Finite difference parameters (4th-order central stencils)
  real(8),parameter :: lc2 = -1d0/12d0, lc1 = 4d0/3d0, lc0 = -5d0/2d0
  real(8),parameter :: gc2 = -1d0/12d0, gc1 = 2d0/3d0

! Numerical parameters kept near the top for easy changes.
  integer,parameter :: nocc = 3
  integer,parameter :: hf_max_iter = 4000
  integer,parameter :: hf_output_stride = 20
  integer,parameter :: output_stride = 1
  real(8),parameter :: hf_energy_tol = 1d-11
  real(8),parameter :: hf_residual_tol = 1d-11

  integer :: nx, nt
  real(8) :: dx, dt
  real(8) :: Tprop

! Material parameters
  real(8) :: lattice_constant
  real(8) :: bvc_lattice_constant

! laser parameters
  real(8) :: E0, omega, Tpulse, phi_CEP

! grids
  real(8), allocatable :: xn(:)

! occupied single-particle orbitals
  complex(8), allocatable :: zphi(:,:)

! potentials
  real(8), allocatable :: vpot_1d(:), wpot_1d(:)

end module global_variables
!-------------------------------------------------------
program main
  use global_variables
  implicit none
  real(8) :: e0_hf, res_hf, norm_err, overlap_err

  call initialize
  call initialize_orbitals
  call ground_state_hf(e0_hf, res_hf)

  call orbital_diagnostics(norm_err, overlap_err)
  write(*,'(a,1pe16.8)') 'HF ground-state energy     = ', e0_hf
  write(*,'(a,1pe16.8)') 'HF residual norm           = ', res_hf
  write(*,'(a,1pe16.8)') 'Particle number            = ', particle_number()
  write(*,'(a,1pe16.8)') 'Max orbital norm error     = ', norm_err
  write(*,'(a,1pe16.8)') 'Max orbital overlap        = ', overlap_err

  call propagate_tdhf
  call finalize

contains
!-------------------------------------------------------
subroutine initialize
  implicit none

  call read_input_parameters
  call set_grids

  allocate(zphi(0:nx-1, nocc))

  call set_potentials
  call check_time_step

end subroutine initialize
!-------------------------------------------------------
subroutine read_input_parameters
  implicit none
  real(8) :: Tprop_fs
  real(8) :: E0_MVm, omega_ev, Tpulse_fs, phi_CEP_2pi

  read(*,*)lattice_constant, nx
  read(*,*)Tprop_fs, dt
  read(*,*)E0_MVm, omega_ev, Tpulse_fs, phi_CEP_2pi

  write(*,*)'lattice_constant = ', lattice_constant
  write(*,*)'nx = ', nx
  write(*,*)'Tprop_fs = ', Tprop_fs
  write(*,*)'dt = ', dt
  write(*,*)'E0_MVm = ', E0_MVm
  write(*,*)'omega_ev = ', omega_ev
  write(*,*)'Tpulse_fs = ', Tpulse_fs
  write(*,*)'phi_CEP_2pi = ', phi_CEP_2pi

  Tprop = Tprop_fs*fs
  nt = max(1, nint(Tprop/dt))+1
  dt = Tprop/dble(nt)
  write(*,*)'dt (refined) = ', dt
  write(*,*)'nt = ', nt

  E0 = E0_MVm*1d-6*ev/(bohr*1d-10)
  omega = omega_ev*ev
  Tpulse = Tpulse_fs*fs
  phi_CEP = phi_CEP_2pi*2d0*pi

  bvc_lattice_constant = lattice_constant*3d0

end subroutine read_input_parameters
!-------------------------------------------------------
subroutine set_grids
  implicit none
  integer :: ix

  allocate(xn(0:nx-1))
  dx = bvc_lattice_constant/dble(nx)

  do ix = 0, nx-1
    xn(ix) = ix*dx
  end do

end subroutine set_grids
!-------------------------------------------------------
subroutine set_potentials
  implicit none
  integer :: ix1
  real(8) :: x1
  real(8),parameter :: v0 = 0.11813d0
  real(8),parameter :: w0 = 0.01d0

  allocate(vpot_1d(0:nx-1))
  allocate(wpot_1d(0:nx-1))

  do ix1 = 0, nx-1
    x1 = xn(ix1)
    vpot_1d(ix1) = v0*(cos(2d0*pi*x1/lattice_constant) &
        + 0.5d0*sin(4d0*pi*x1/lattice_constant))
  end do

! Pair potential tabulated by the minimum Born-von Karman distance, as in tdse.f90.
  do ix1 = 0, nx-1
    x1 = dble(min(ix1, nx-ix1))*dx
    wpot_1d(ix1) = w0*cos(pi*x1/bvc_lattice_constant)**16
  end do

end subroutine set_potentials
!-------------------------------------------------------
subroutine check_time_step
  implicit none
  real(8) :: hmax_est

! One-particle fourth-order kinetic estimate.  The nonlinear HF field can add
! to this bound when w0 is changed from the tdse.f90 default.
  hmax_est = 8d0/(3d0*dx**2) + maxval(abs(vpot_1d)) &
      + dble(nocc)*maxval(abs(wpot_1d))
  if (dt*hmax_est > 2.5d0) then
    write(*,'(a,1pe12.4,a)') 'Warning: RK4 time step may be unstable; dt*Hmax ~= ', &
        dt*hmax_est, '.'
  end if

end subroutine check_time_step
!-------------------------------------------------------
integer function periodic_distance_index(i, j)
  implicit none
  integer,intent(in) :: i, j
  integer :: d

  d = abs(i-j)
  periodic_distance_index = min(d, nx-d)

end function periodic_distance_index
!-------------------------------------------------------
integer function ipbc(i)
  implicit none
  integer,intent(in) :: i

  ipbc = modulo(i, nx)

end function ipbc
!-------------------------------------------------------
subroutine initialize_orbitals
  implicit none
  integer :: ix
  real(8) :: k(nocc)

! Start from the same three periodic plane waves used to build the TDSE
! Slater determinant.  The following imaginary-time HF stage relaxes them
! toward a static self-consistent determinant.
  k(1) = -2d0*pi/bvc_lattice_constant
  k(2) = 0d0
  k(3) =  2d0*pi/bvc_lattice_constant

  do ix = 0, nx-1
    zphi(ix,1) = exp(zi*k(1)*xn(ix))/sqrt(bvc_lattice_constant)
    zphi(ix,2) = exp(zi*k(2)*xn(ix))/sqrt(bvc_lattice_constant)
    zphi(ix,3) = exp(zi*k(3)*xn(ix))/sqrt(bvc_lattice_constant)
  end do

  call orthonormalize_orbitals(zphi)

end subroutine initialize_orbitals
!-------------------------------------------------------
complex(8) function orbital_inner(a, b)
  implicit none
  complex(8),intent(in) :: a(0:nx-1), b(0:nx-1)

  orbital_inner = sum(conjg(a)*b)*dx

end function orbital_inner
!-------------------------------------------------------
subroutine orthonormalize_orbitals(phi)
  implicit none
  complex(8),intent(inout) :: phi(0:nx-1,nocc)
  integer :: i, j
  complex(8) :: zs
  real(8) :: nrm

! Modified Gram-Schmidt is sufficient for nocc=3 and is applied after each
! real-time RK4 step to remove nonlinear/time-discretization drift.
  do i = 1, nocc
    do j = 1, i-1
      zs = orbital_inner(phi(:,j), phi(:,i))
      phi(:,i) = phi(:,i) - zs*phi(:,j)
    end do
    nrm = sqrt(max(0d0, real(orbital_inner(phi(:,i), phi(:,i)))))
    if (nrm <= 0d0) stop 'Cannot normalize a zero orbital.'
    phi(:,i) = phi(:,i)/nrm
  end do

end subroutine orthonormalize_orbitals
!-------------------------------------------------------
subroutine orbital_diagnostics(norm_err, overlap_err)
  implicit none
  real(8),intent(out) :: norm_err, overlap_err
  integer :: i, j
  complex(8) :: zs

  norm_err = 0d0
  overlap_err = 0d0
  do i = 1, nocc
    zs = orbital_inner(zphi(:,i), zphi(:,i))
    norm_err = max(norm_err, abs(real(zs)-1d0))
    do j = 1, i-1
      zs = orbital_inner(zphi(:,j), zphi(:,i))
      overlap_err = max(overlap_err, abs(zs))
    end do
  end do

end subroutine orbital_diagnostics
!-------------------------------------------------------
real(8) function particle_number()
  implicit none
  integer :: i

  particle_number = 0d0
  do i = 1, nocc
    particle_number = particle_number + real(orbital_inner(zphi(:,i), zphi(:,i)))
  end do

end function particle_number
!-------------------------------------------------------
subroutine density_from_orbitals(phi, rho)
  implicit none
  complex(8),intent(in) :: phi(0:nx-1,nocc)
  real(8),intent(out) :: rho(0:nx-1)
  integer :: ix, i

  rho = 0d0
  do i = 1, nocc
    do ix = 0, nx-1
      rho(ix) = rho(ix) + abs(phi(ix,i))**2
    end do
  end do

end subroutine density_from_orbitals
!-------------------------------------------------------
subroutine apply_onebody(phi, hphi, avec)
  implicit none
  complex(8),intent(in) :: phi(0:nx-1,nocc)
  complex(8),intent(out) :: hphi(0:nx-1,nocc)
  real(8),intent(in) :: avec
  integer :: ix, i
  integer :: ixp1, ixp2, ixm1, ixm2
  complex(8) :: zc_0, zc_p1, zc_p2, zc_m1, zc_m2

  zc_0  = -0.5d0*lc0/dx**2 + 0.5d0*avec**2
  zc_p1 = -0.5d0*lc1/dx**2 - zi*gc1*avec/dx
  zc_p2 = -0.5d0*lc2/dx**2 - zi*gc2*avec/dx
  zc_m1 = -0.5d0*lc1/dx**2 + zi*gc1*avec/dx
  zc_m2 = -0.5d0*lc2/dx**2 + zi*gc2*avec/dx

!$omp parallel do private(ix, i, ixp1, ixp2, ixm1, ixm2)
  do i = 1, nocc
    do ix = 0, nx-1
      ixp1 = ipbc(ix+1)
      ixp2 = ipbc(ix+2)
      ixm1 = ipbc(ix-1)
      ixm2 = ipbc(ix-2)

      hphi(ix,i) = (zc_0 + vpot_1d(ix))*phi(ix,i) &
          + zc_p1*phi(ixp1,i) + zc_m1*phi(ixm1,i) &
          + zc_p2*phi(ixp2,i) + zc_m2*phi(ixm2,i)
    end do
  end do

end subroutine apply_onebody
!-------------------------------------------------------
subroutine apply_hf_hamiltonian(phi, hphi, avec)
  implicit none
  complex(8),intent(in) :: phi(0:nx-1,nocc)
  complex(8),intent(out) :: hphi(0:nx-1,nocc)
  real(8),intent(in) :: avec
  real(8),allocatable :: rho(:), vh(:)
  complex(8),allocatable :: fock(:,:)
  integer :: ix, jx, i, j, id
  complex(8) :: zsum

  allocate(rho(0:nx-1))
  allocate(vh(0:nx-1))
  allocate(fock(0:nx-1,nocc))

  call apply_onebody(phi, hphi, avec)
  call density_from_orbitals(phi, rho)

! Hartree field: V_H(x) = integral dx' w(x-x') rho(x').
  vh = 0d0
!$omp parallel do private(ix, jx, id) reduction(+:vh)
  do ix = 0, nx-1
    do jx = 0, nx-1
      id = periodic_distance_index(ix, jx)
      vh(ix) = vh(ix) + wpot_1d(id)*rho(jx)*dx
    end do
  end do

! Fock field: (K_j phi_i)(x) = phi_j(x) integral dx' w(x-x')
!                              phi_j*(x') phi_i(x').
  fock = (0d0,0d0)
!$omp parallel do private(i, j, ix, jx, id, zsum)
  do i = 1, nocc
    do ix = 0, nx-1
      do j = 1, nocc
        zsum = (0d0,0d0)
        do jx = 0, nx-1
          id = periodic_distance_index(ix, jx)
          zsum = zsum + wpot_1d(id)*conjg(phi(jx,j))*phi(jx,i)*dx
        end do
        fock(ix,i) = fock(ix,i) + phi(ix,j)*zsum
      end do
    end do
  end do

  do i = 1, nocc
    do ix = 0, nx-1
      hphi(ix,i) = hphi(ix,i) + vh(ix)*phi(ix,i) - fock(ix,i)
    end do
  end do

  deallocate(rho, vh, fock)

end subroutine apply_hf_hamiltonian
!-------------------------------------------------------
subroutine ground_state_hf(energy, residual_norm)
  implicit none
  real(8),intent(out) :: energy, residual_norm
  complex(8),allocatable :: hphi(:,:), phi_old(:,:)
  real(8) :: energy_old, dtau, norm_err, overlap_err
  integer :: iter

  allocate(hphi(0:nx-1,nocc))
  allocate(phi_old(0:nx-1,nocc))

! A compact imaginary-time HF preparation.  This is not a production
! diagonalization/SCF solver, but it gives a self-consistent determinant for
! the same lattice Hamiltonian without introducing external dependencies.
  dtau = min(0.05d0, 0.15d0/(8d0/(3d0*dx**2) + maxval(abs(vpot_1d)) &
      + dble(nocc)*maxval(abs(wpot_1d)) + 1d-12))

  call orthonormalize_orbitals(zphi)
  call apply_hf_hamiltonian(zphi, hphi, 0d0)
  energy = total_hf_energy(zphi, 0d0)
  residual_norm = hf_residual_norm(zphi, hphi)
  energy_old = energy

  open(20,file='ground_state_hf.log',status='replace')
  write(20,'(a)') '# iter energy residual_norm particle_number max_norm_error max_overlap'
  call orbital_diagnostics(norm_err, overlap_err)
  write(20,'(i8,5(1x,1pe20.12))') 0, energy, residual_norm, &
      particle_number(), norm_err, overlap_err

  do iter = 1, hf_max_iter
    phi_old = zphi
    zphi = zphi - dtau*hphi
    call orthonormalize_orbitals(zphi)

    call apply_hf_hamiltonian(zphi, hphi, 0d0)
    energy = total_hf_energy(zphi, 0d0)
    residual_norm = hf_residual_norm(zphi, hphi)

    if (mod(iter, hf_output_stride) == 0 .or. iter == 1) then
      call orbital_diagnostics(norm_err, overlap_err)
      write(20,'(i8,5(1x,1pe20.12))') iter, energy, residual_norm, &
          particle_number(), norm_err, overlap_err
    end if

    if (abs(energy-energy_old) < hf_energy_tol .and. &
        residual_norm < hf_residual_tol) exit
    if (energy > energy_old + 1d-8) then
! Back off once if the explicit imaginary-time step overshoots.
      zphi = phi_old
      dtau = 0.5d0*dtau
      call apply_hf_hamiltonian(zphi, hphi, 0d0)
      energy = total_hf_energy(zphi, 0d0)
      residual_norm = hf_residual_norm(zphi, hphi)
    end if
    energy_old = energy
  end do

  if (mod(iter, hf_output_stride) /= 0) then
    call orbital_diagnostics(norm_err, overlap_err)
    write(20,'(i8,5(1x,1pe20.12))') min(iter,hf_max_iter), energy, residual_norm, &
        particle_number(), norm_err, overlap_err
  end if
  close(20)

  deallocate(hphi, phi_old)

end subroutine ground_state_hf
!-------------------------------------------------------
real(8) function hf_residual_norm(phi, hphi)
  implicit none
  complex(8),intent(in) :: phi(0:nx-1,nocc), hphi(0:nx-1,nocc)
  complex(8) :: lambda(nocc,nocc)
  complex(8),allocatable :: res(:,:)
  integer :: i, j

! Stationary HF equations are satisfied up to rotations within the occupied
! subspace: h phi_i = sum_j phi_j lambda_ji.
  allocate(res(0:nx-1,nocc))
  do i = 1, nocc
    do j = 1, nocc
      lambda(j,i) = orbital_inner(phi(:,j), hphi(:,i))
    end do
  end do

  res = hphi
  do i = 1, nocc
    do j = 1, nocc
      res(:,i) = res(:,i) - phi(:,j)*lambda(j,i)
    end do
  end do

  hf_residual_norm = sqrt(max(0d0, real(sum(conjg(res)*res)*dx)))
  deallocate(res)

end function hf_residual_norm
!-------------------------------------------------------
real(8) function total_hf_energy(phi, avec)
  implicit none
  complex(8),intent(in) :: phi(0:nx-1,nocc)
  real(8),intent(in) :: avec
  complex(8),allocatable :: h0phi(:,:)
  real(8),allocatable :: rho(:), vh(:)
  real(8) :: e_one, e_hartree, e_exchange
  integer :: ix, jx, i, j, id
  complex(8) :: ztmp

  allocate(h0phi(0:nx-1,nocc))
  allocate(rho(0:nx-1))
  allocate(vh(0:nx-1))

  call apply_onebody(phi, h0phi, avec)
  e_one = 0d0
  do i = 1, nocc
    e_one = e_one + real(orbital_inner(phi(:,i), h0phi(:,i)))
  end do

  call density_from_orbitals(phi, rho)
  vh = 0d0
  do ix = 0, nx-1
    do jx = 0, nx-1
      id = periodic_distance_index(ix, jx)
      vh(ix) = vh(ix) + wpot_1d(id)*rho(jx)*dx
    end do
  end do
  e_hartree = 0.5d0*sum(rho*vh)*dx

  e_exchange = 0d0
  do i = 1, nocc
    do j = 1, nocc
      do ix = 0, nx-1
        ztmp = (0d0,0d0)
        do jx = 0, nx-1
          id = periodic_distance_index(ix, jx)
          ztmp = ztmp + wpot_1d(id)*conjg(phi(jx,j))*phi(jx,i)*dx
        end do
        e_exchange = e_exchange + real(conjg(phi(ix,i))*phi(ix,j)*ztmp)*dx
      end do
    end do
  end do
  e_exchange = -0.5d0*e_exchange

  total_hf_energy = e_one + e_hartree + e_exchange

  deallocate(h0phi, rho, vh)

end function total_hf_energy
!-------------------------------------------------------
subroutine propagate_tdhf
  implicit none
  integer :: it
  real(8) :: t, avec, norm_err, overlap_err
  complex(8),allocatable :: hphi(:,:)
! physics
  real(8),allocatable :: current_t(:), energy_t(:), pnum_t(:)
  real(8),allocatable :: normerr_t(:), overlap_t(:)

  allocate(current_t(0:nt))
  allocate(energy_t(0:nt))
  allocate(pnum_t(0:nt))
  allocate(normerr_t(0:nt))
  allocate(overlap_t(0:nt))

  allocate(hphi(0:nx-1,nocc))

  do it = 0, nt
    write(*,'(a,i8)')'it = ', it
    t = dble(it)*dt
    avec = vector_potential(t)
    call apply_hf_hamiltonian(zphi, hphi, avec)
    energy_t(it)  = total_hf_energy(zphi, avec)
    current_t(it) = total_current(zphi, avec)
    pnum_t(it)    = particle_number()
    call orbital_diagnostics(norm_err, overlap_err)
    normerr_t(it) = norm_err
    overlap_t(it) = overlap_err

    if (it < nt) call rk4_step(t, dt)
  end do

  open(30,file='current_tdhf.out',status='replace')
  do it = 0, nt
    t = dble(it)*dt
    write(30,"(999e26.16e3)")t, vector_potential(t), current_t(it), &
        pnum_t(it), energy_t(it), normerr_t(it), overlap_t(it)
  end do
  close(30)

  deallocate(hphi)

end subroutine propagate_tdhf
!-------------------------------------------------------
subroutine rk4_step(t, h)
  implicit none
  real(8),intent(in) :: t, h
  complex(8),allocatable :: y0(:,:), yt(:,:), k1(:,:), k2(:,:)
  complex(8),allocatable :: k3(:,:), k4(:,:)

  allocate(y0(0:nx-1,nocc))
  allocate(yt(0:nx-1,nocc))
  allocate(k1(0:nx-1,nocc))
  allocate(k2(0:nx-1,nocc))
  allocate(k3(0:nx-1,nocc))
  allocate(k4(0:nx-1,nocc))

  y0 = zphi
  call tdhf_rhs(y0, k1, t)
  yt = y0 + 0.5d0*h*k1
  call orthonormalize_orbitals(yt)
  call tdhf_rhs(yt, k2, t+0.5d0*h)
  yt = y0 + 0.5d0*h*k2
  call orthonormalize_orbitals(yt)
  call tdhf_rhs(yt, k3, t+0.5d0*h)
  yt = y0 + h*k3
  call orthonormalize_orbitals(yt)
  call tdhf_rhs(yt, k4, t+h)

  zphi = y0 + h*(k1 + 2d0*k2 + 2d0*k3 + k4)/6d0
  call orthonormalize_orbitals(zphi)

  deallocate(y0, yt, k1, k2, k3, k4)

end subroutine rk4_step
!-------------------------------------------------------
subroutine tdhf_rhs(phi, rhs, t)
  implicit none
  complex(8),intent(in) :: phi(0:nx-1,nocc)
  complex(8),intent(out) :: rhs(0:nx-1,nocc)
  real(8),intent(in) :: t

  call apply_hf_hamiltonian(phi, rhs, vector_potential(t))
  rhs = -zi*rhs

end subroutine tdhf_rhs
!-------------------------------------------------------
real(8) function vector_potential(t)
  implicit none
  real(8),intent(in) :: t
  real(8) :: env

  if (t < 0d0 .or. t > Tpulse) then
    vector_potential = 0d0
  else
    env = sin(pi*t/Tpulse)**4
! A(t) is chosen so that the field is approximately E(t)=-dA/dt for a
! slowly varying envelope; the exact A(t) is what enters the Hamiltonian.
    vector_potential = -(E0/omega)*env*sin(omega*(t-0.5d0*tpulse) + phi_CEP)
  end if

end function vector_potential
!-------------------------------------------------------
real(8) function total_current(phi, avec)
  implicit none
  complex(8),intent(in) :: phi(0:nx-1,nocc)
  real(8),intent(in) :: avec
  integer :: ix, i
  integer :: ixp1, ixp2, ixm1, ixm2
  real(8) :: curr_tmp
  complex(8) :: zc_p1, zc_p2, zc_m1, zc_m2

  zc_p1 = -zi*gc1/dx
  zc_p2 = -zi*gc2/dx
  zc_m1 =  zi*gc1/dx
  zc_m2 =  zi*gc2/dx

  curr_tmp = 0d0

! Sum of the single-particle current expectation values:
! <phi_i| -i d/dx + A(t) |phi_i>.
!$omp parallel do private(i, ix, ixp1, ixp2, ixm1, ixm2) reduction(+:curr_tmp)
  do i = 1, nocc
    do ix = 0, nx-1
      ixp1 = ipbc(ix+1)
      ixp2 = ipbc(ix+2)
      ixm1 = ipbc(ix-1)
      ixm2 = ipbc(ix-2)

      curr_tmp = curr_tmp + real(conjg(phi(ix,i))*(zc_p1*phi(ixp1,i) &
          + zc_m1*phi(ixm1,i) + zc_p2*phi(ixp2,i) &
          + zc_m2*phi(ixm2,i))) + avec*abs(phi(ix,i))**2
    end do
  end do

  total_current = curr_tmp*dx

end function total_current
!-------------------------------------------------------
subroutine finalize
  implicit none

  if (allocated(xn)) deallocate(xn)
  if (allocated(zphi)) deallocate(zphi)
  if (allocated(vpot_1d)) deallocate(vpot_1d)
  if (allocated(wpot_1d)) deallocate(wpot_1d)

end subroutine finalize
!-------------------------------------------------------
end program main
