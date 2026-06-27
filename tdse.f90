module global_variables
  implicit none
! Build:
!   gfortran -std=f95 -O2 tdse.f90 -o tdse
! Run:
!   ./tdse < input_tdse
! Output:
!   ground_state.log : CG iteration, energy, residual, norm, antisymmetry check
!   current.dat      : t, A(t), j(t), norm, energy, antisymmetry check
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
  integer,parameter :: cg_max_iter = 600
  integer,parameter :: output_stride = 1
  real(8),parameter :: cg_energy_tol = 1d-11
  real(8),parameter :: cg_residual_tol = 1d-9

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

! wavefunction
  complex(8), allocatable :: zpsi(:,:,:)

! potentials
  real(8), allocatable :: vpot_1d(:), wpot_1d(:)
  real(8), allocatable :: vpot(:,:,:), wpot(:,:,:), tot_pot(:,:,:)

end module global_variables
!-------------------------------------------------------
program main
  use global_variables
  implicit none
  real(8) :: e0_gs, res_gs, asym_gs

  call initialize
  call initialize_antisymmetric_state
  call ground_state_cg(e0_gs, res_gs)

  asym_gs = antisymmetry_error(zpsi)
  write(*,'(a,1pe16.8)') 'Ground-state energy        = ', e0_gs
  write(*,'(a,1pe16.8)') 'Ground-state residual norm = ', res_gs
  write(*,'(a,1pe16.8)') 'Ground-state norm          = ', wavefunction_norm(zpsi)
  write(*,'(a,1pe16.8)') 'Antisymmetry error         = ', asym_gs

  call propagate_tdse
  call finalize

contains
!-------------------------------------------------------
subroutine initialize
  implicit none

  call read_input_parameters
  call set_grids

  allocate(zpsi(0:nx-1, 0:nx-1, 0:nx-1))
  allocate(vpot(0:nx-1, 0:nx-1, 0:nx-1))
  allocate(wpot(0:nx-1, 0:nx-1, 0:nx-1))
  allocate(tot_pot(0:nx-1, 0:nx-1, 0:nx-1))

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
  integer :: ix1, ix2, ix3
  integer :: id12, id23, id31
  real(8) :: x1
  real(8),parameter :: v0 = 0.11813d0
  real(8),parameter :: w0 = 1d0*0d0

  allocate(vpot_1d(0:nx-1))
  allocate(wpot_1d(0:nx-1))

  do ix1 = 0, nx-1
    x1 = xn(ix1)
    vpot_1d(ix1) = v0*(cos(2d0*pi*x1/lattice_constant) &
        + 0.5d0*sin(4d0*pi*x1/lattice_constant))
  end do

! Pair potential tabulated by the minimum Born-von Karman distance.
  do ix1 = 0, nx-1
    x1 = dble(min(ix1, nx-ix1))*dx
    wpot_1d(ix1) = w0*cos(pi*x1/bvc_lattice_constant)**16
  end do

  do ix1 = 0, nx-1
    do ix2 = 0, nx-1
      do ix3 = 0, nx-1
        id12 = periodic_distance_index(ix1, ix2)
        id23 = periodic_distance_index(ix2, ix3)
        id31 = periodic_distance_index(ix3, ix1)
        vpot(ix1, ix2, ix3) = vpot_1d(ix1) + vpot_1d(ix2) + vpot_1d(ix3)
        wpot(ix1, ix2, ix3) = wpot_1d(id12) + wpot_1d(id23) + wpot_1d(id31)
      end do
    end do
  end do

  tot_pot = vpot + wpot

end subroutine set_potentials
!-------------------------------------------------------
subroutine check_time_step
  implicit none
  real(8) :: hmax_est

! For the fourth-order Laplacian, the largest one-particle kinetic eigenvalue
! is 8/(3 dx**2); three particles give 8/dx**2 before adding potentials.
  hmax_est = 8d0/dx**2 + maxval(abs(tot_pot))
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
subroutine initialize_antisymmetric_state
  implicit none
  integer :: ix1, ix2, ix3
  real(8) :: k(3)
  complex(8) :: a(3), b(3), c(3)

  k(1) = -2d0*pi/bvc_lattice_constant
  k(2) = 0d0
  k(3) =  2d0*pi/bvc_lattice_constant

  do ix1 = 0, nx-1
    do ix2 = 0, nx-1
      do ix3 = 0, nx-1
        a = exp(zi*k*xn(ix1))
        b = exp(zi*k*xn(ix2))
        c = exp(zi*k*xn(ix3))
        zpsi(ix1,ix2,ix3) = determinant3(a, b, c)
      end do
    end do
  end do

  call antisymmetrize(zpsi)
  call normalize(zpsi)

end subroutine initialize_antisymmetric_state
!-------------------------------------------------------
complex(8) function determinant3(a, b, c)
  implicit none
  complex(8),intent(in) :: a(3), b(3), c(3)

  determinant3 = a(1)*(b(2)*c(3)-b(3)*c(2)) &
      - a(2)*(b(1)*c(3)-b(3)*c(1)) &
      + a(3)*(b(1)*c(2)-b(2)*c(1))

end function determinant3
!-------------------------------------------------------
subroutine antisymmetrize(psi)
  implicit none
  complex(8),intent(inout) :: psi(0:nx-1,0:nx-1,0:nx-1)
  integer :: i, j, k
  complex(8) :: a

  do i = 0, nx-1
    psi(i,i,:) = (0d0,0d0)
    psi(i,:,i) = (0d0,0d0)
    psi(:,i,i) = (0d0,0d0)
  end do

!$omp parallel do private(i,j,k,a)
  do i = 0, nx-3
    do j = i+1, nx-2
      do k = j+1, nx-1
        a = (psi(i,j,k) + psi(j,k,i) + psi(k,i,j) &
           - psi(j,i,k) - psi(i,k,j) - psi(k,j,i))/6d0
        psi(i,j,k) =  a
        psi(j,k,i) =  a
        psi(k,i,j) =  a
        psi(j,i,k) = -a
        psi(i,k,j) = -a
        psi(k,j,i) = -a
      end do
    end do
  end do

end subroutine antisymmetrize
!-------------------------------------------------------
subroutine normalize(psi)
  implicit none
  complex(8),intent(inout) :: psi(0:nx-1,0:nx-1,0:nx-1)
  real(8) :: nrm

  nrm = wavefunction_norm(psi)
  if (nrm <= 0d0) stop 'Cannot normalize a zero wavefunction.'
  psi = psi/nrm

end subroutine normalize
!-------------------------------------------------------
real(8) function wavefunction_norm(psi)
  implicit none
  complex(8),intent(in) :: psi(0:nx-1,0:nx-1,0:nx-1)

  wavefunction_norm = sqrt(max(0d0, real(inner_product(psi, psi))))

end function wavefunction_norm
!-------------------------------------------------------
complex(8) function inner_product(a, b)
  implicit none
  complex(8),intent(in) :: a(0:nx-1,0:nx-1,0:nx-1)
  complex(8),intent(in) :: b(0:nx-1,0:nx-1,0:nx-1)

  inner_product = sum(conjg(a)*b)*dx**3

end function inner_product
!-------------------------------------------------------
subroutine apply_hamiltonian(psi, hpsi, avec)
  implicit none
  complex(8),intent(in) :: psi(0:nx-1,0:nx-1,0:nx-1)
  complex(8),intent(out) :: hpsi(0:nx-1,0:nx-1,0:nx-1)
  real(8),intent(in) :: avec
  integer :: ix1, ix2, ix3
  integer :: ix1p1, ix2p1, ix3p1
  integer :: ix1p2, ix2p2, ix3p2
  integer :: ix1m1, ix2m1, ix3m1
  integer :: ix1m2, ix2m2, ix3m2
  complex(8) :: zc_0, zc_p1, zc_p2, zc_m1, zc_m2


  zc_0  = 3d0*(-0.5d0*lc0/dx**2 + 0.5d0*avec**2)
  zc_p1 = -0.5d0*lc1/dx**2 - zi*gc1*avec/dx
  zc_p2 = -0.5d0*lc2/dx**2 - zi*gc2*avec/dx
  zc_m1 = -0.5d0*lc1/dx**2 + zi*gc1*avec/dx
  zc_m2 = -0.5d0*lc2/dx**2 + zi*gc2*avec/dx

!$omp parallel do private(ix1, ix2, ix3, ix1p1, ix2p1, ix3p1, &
!$omp    ix1p2, ix2p2, ix3p2, ix1m1, ix2m1, ix3m1, &
!$omp    ix1m2, ix2m2, ix3m2)
  do ix1 = 0, nx-1

    ix1p1 = ipbc(ix1+1)
    ix1p2 = ipbc(ix1+2)
    ix1m1 = ipbc(ix1-1)
    ix1m2 = ipbc(ix1-2)

    do ix2 = 0, nx-1

      ix2p1 = ipbc(ix2+1)
      ix2p2 = ipbc(ix2+2)
      ix2m1 = ipbc(ix2-1)
      ix2m2 = ipbc(ix2-2)

      do ix3 = 0, nx-1
        
        ix3p1 = ipbc(ix3+1)
        ix3p2 = ipbc(ix3+2)
        ix3m1 = ipbc(ix3-1)
        ix3m2 = ipbc(ix3-2)

        
        hpsi(ix1,ix2,ix3) = (zc_0+tot_pot(ix1, ix2, ix3))*psi(ix1,ix2,ix3) &
            + zc_p1*(psi(ix1p1,ix2,ix3) + psi(ix1,ix2p1,ix3) &
            + psi(ix1,ix2,ix3p1)) &
            + zc_m1*(psi(ix1m1,ix2,ix3) + psi(ix1,ix2m1,ix3) &
            + psi(ix1,ix2,ix3m1)) &
            + zc_p2*(psi(ix1p2,ix2,ix3) + psi(ix1,ix2p2,ix3) &
            + psi(ix1,ix2,ix3p2)) &
            + zc_m2*(psi(ix1m2,ix2,ix3) + psi(ix1,ix2m2,ix3) &
            + psi(ix1,ix2,ix3m2))
        
      end do
    end do
  end do
! The finite-difference Hamiltonian is permutation symmetric; the projection
! removes round-off drift and exactly zeros Pauli-forbidden diagonal components.
!  call antisymmetrize(hpsi)

end subroutine apply_hamiltonian
!-------------------------------------------------------
!> The algorithm of Jiang et al. Phys. Rev. B 68, 165337 (2003) is impleemnted with a sign fix of Eq. (21).
subroutine ground_state_cg(energy, residual_norm)
  implicit none
  real(8),intent(out) :: energy, residual_norm
  complex(8),allocatable :: hpsi(:,:,:)
  complex(8),allocatable :: xi(:,:,:), phi_t(:,:,:), phi_old(:,:,:)
  real(8) :: lambda, xixi, xixi_old, gamma, theta
  integer :: iter
  real(8) :: ss, aa, bb
  complex(8) :: zs

  allocate(hpsi(0:nx-1,0:nx-1,0:nx-1))
  allocate(xi(0:nx-1,0:nx-1,0:nx-1))
  allocate(phi_t(0:nx-1,0:nx-1,0:nx-1))
  allocate(phi_old(0:nx-1,0:nx-1,0:nx-1))

  call antisymmetrize(zpsi)
  call normalize(zpsi)
  call apply_hamiltonian(zpsi, hpsi, 0d0)
  lambda = real(inner_product(zpsi, hpsi))

  xi = lambda*zpsi - hpsi 
  call antisymmetrize(xi)
  xixi = real(inner_product(xi, xi))
  xixi_old = xixi

  residual_norm = xixi
  energy = lambda

  open(20,file='ground_state.log',status='replace')
  write(20,'(a)') '# iter energy residual_norm norm antisymmetry_error'
  write(20,'(i8,4(1x,1pe20.12))') 0, energy, residual_norm, &
      wavefunction_norm(zpsi), antisymmetry_error(zpsi)

  do iter = 0, cg_max_iter

    if(iter == 0)then
      gamma = 0d0
      phi_old = 0d0
    else
      gamma = xixi/xixi_old
      xixi_old = xixi
    end if

    phi_t = xi + gamma*phi_old
    phi_old = phi_t
    call antisymmetrize(phi_t)

    zs = inner_product(zpsi, phi_t)
    phi_t = phi_t - zs*zpsi
    call antisymmetrize(phi_t)
    
    ss = real(inner_product(phi_t, phi_t))
    phi_t = phi_t/sqrt(ss)


    bb = 2d0*real(inner_product(phi_t, hpsi))
    call apply_hamiltonian(phi_t, hpsi, 0d0)
    aa = real(inner_product(phi_t, hpsi)) - lambda
    aa = -aa ! fix: there is a sign error in Eq. (21) of Jiang et al. Phys. Rev. B 68, 165337 (2003)

!    theta = 0.5d0*atan2(bb, aa)
    if(aa /= 0d0)then
      theta = 0.5d0*atan(bb/aa)
    else
      if(bb > 0d0)then
        theta = 0.25d0*pi
      else
        theta = -0.25d0*pi
      end if
    end if
    zpsi = cos(theta)*zpsi + sin(theta)*phi_t
    call antisymmetrize(zpsi)
    call normalize(zpsi)

    if(iter == cg_max_iter)exit

! calc xi
    call apply_hamiltonian(zpsi, hpsi, 0d0)
    lambda = real(inner_product(zpsi, hpsi))
    xi = lambda*zpsi - hpsi
    call antisymmetrize(xi)
    xixi = real(inner_product(xi, xi))

    residual_norm = xixi
    energy = lambda

    write(20,'(i8,4(1x,1pe20.12))') iter+1, energy, residual_norm, &
        wavefunction_norm(zpsi), antisymmetry_error(zpsi)

  end do


  close(20)

end subroutine ground_state_cg
!-------------------------------------------------------
subroutine rayleigh_ritz_update(psi, p, hpsi, hp, energy)
  implicit none
  complex(8),intent(inout) :: psi(0:nx-1,0:nx-1,0:nx-1)
  complex(8),intent(in) :: p(0:nx-1,0:nx-1,0:nx-1)
  complex(8),intent(in) :: hpsi(0:nx-1,0:nx-1,0:nx-1)
  complex(8),intent(in) :: hp(0:nx-1,0:nx-1,0:nx-1)
  real(8),intent(out) :: energy
  real(8) :: a, c, delta, lambda
  complex(8) :: b, ratio

  a = real(inner_product(psi, hpsi))
  c = real(inner_product(p, hp))
  b = inner_product(psi, hp)
  delta = sqrt((a-c)**2 + 4d0*abs(b)**2)
  lambda = 0.5d0*(a + c - delta)

  if (abs(b) > 1d-300) then
    ratio = -(a-lambda)/b
    psi = psi + ratio*p
  else if (c < a) then
    psi = p
  end if

  call antisymmetrize(psi)
  energy = lambda

end subroutine rayleigh_ritz_update
!-------------------------------------------------------
subroutine propagate_tdse
  implicit none
  integer :: it
  real(8) :: t, avec
  complex(8),allocatable :: hpsi(:,:,:)
! physics
  real(8),allocatable :: current_t(:), energy_t(:), norm_t(:)

  allocate(current_t(0:nt))
  allocate(energy_t(0:nt))
  allocate(norm_t(0:nt))

  allocate(hpsi(0:nx-1,0:nx-1,0:nx-1))


  do it = 0, nt
    write(*,'(a,i8)')'it = ', it
    t = dble(it)*dt
    avec = vector_potential(t)
    call apply_hamiltonian(zpsi, hpsi, avec)
    energy_t(it)  = real(inner_product(zpsi, hpsi))
    current_t(it) = total_current(zpsi, avec)
    norm_t(it)    = wavefunction_norm(zpsi)

    if (it < nt) call rk4_step(t, dt)
  end do

  open(30,file='current.out',status='replace')
  do it = 0, nt
    t = dble(it)*dt
    write(30,"(999e26.16e3)")t, vector_potential(t), current_t(it), &
        norm_t(it), energy_t(it)
  end do
  close(30)

  deallocate(hpsi)

end subroutine propagate_tdse
!-------------------------------------------------------
subroutine rk4_step(t, h)
  implicit none
  real(8),intent(in) :: t, h
  complex(8),allocatable :: y0(:,:,:), yt(:,:,:), k1(:,:,:), k2(:,:,:)
  complex(8),allocatable :: k3(:,:,:), k4(:,:,:)

  allocate(y0(0:nx-1,0:nx-1,0:nx-1))
  allocate(yt(0:nx-1,0:nx-1,0:nx-1))
  allocate(k1(0:nx-1,0:nx-1,0:nx-1))
  allocate(k2(0:nx-1,0:nx-1,0:nx-1))
  allocate(k3(0:nx-1,0:nx-1,0:nx-1))
  allocate(k4(0:nx-1,0:nx-1,0:nx-1))

  y0 = zpsi
  call tdse_rhs(y0, k1, t)
  yt = y0 + 0.5d0*h*k1
  call antisymmetrize(yt)
  call tdse_rhs(yt, k2, t+0.5d0*h)
  yt = y0 + 0.5d0*h*k2
  call antisymmetrize(yt)
  call tdse_rhs(yt, k3, t+0.5d0*h)
  yt = y0 + h*k3
  call antisymmetrize(yt)
  call tdse_rhs(yt, k4, t+h)

  zpsi = y0 + h*(k1 + 2d0*k2 + 2d0*k3 + k4)/6d0
  call antisymmetrize(zpsi)

  deallocate(y0, yt, k1, k2, k3, k4)

end subroutine rk4_step
!-------------------------------------------------------
subroutine tdse_rhs(psi, rhs, t)
  implicit none
  complex(8),intent(in) :: psi(0:nx-1,0:nx-1,0:nx-1)
  complex(8),intent(out) :: rhs(0:nx-1,0:nx-1,0:nx-1)
  real(8),intent(in) :: t

  call apply_hamiltonian(psi, rhs, vector_potential(t))
  rhs = -zi*rhs

end subroutine tdse_rhs
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
real(8) function total_current(psi, avec)
  implicit none
  complex(8),intent(in) :: psi(0:nx-1,0:nx-1,0:nx-1)
  real(8),intent(in) :: avec
  integer :: ix1, ix2, ix3
  integer :: ix1p1, ix2p1, ix3p1
  integer :: ix1p2, ix2p2, ix3p2
  integer :: ix1m1, ix2m1, ix3m1
  integer :: ix1m2, ix2m2, ix3m2

  real(8) :: curr_tmp
  complex(8) :: zc_p1, zc_p2, zc_m1, zc_m2

  zc_p1 = -zi*gc1/dx
  zc_p2 = -zi*gc2/dx
  zc_m1 =  zi*gc1/dx
  zc_m2 =  zi*gc2/dx

  curr_tmp = 0d0

!$omp parallel do private(ix1, ix2, ix3, ix1p1, ix2p1, ix3p1, &
!$omp   ix1p2, ix2p2, ix3p2, ix1m1, ix2m1, ix3m1, &
!$omp   ix1m2, ix2m2, ix3m2) reduction(+:curr_tmp)
  do ix1 = 0, nx-1
    ix1p1 = ipbc(ix1+1)
    ix1p2 = ipbc(ix1+2)
    ix1m1 = ipbc(ix1-1)
    ix1m2 = ipbc(ix1-2)

    do ix2 = 0, nx-1
      ix2p1 = ipbc(ix2+1)
      ix2p2 = ipbc(ix2+2)
      ix2m1 = ipbc(ix2-1)
      ix2m2 = ipbc(ix2-2)

      do ix3 = 0, nx-1
        ix3p1 = ipbc(ix3+1)
        ix3p2 = ipbc(ix3+2)
        ix3m1 = ipbc(ix3-1)
        ix3m2 = ipbc(ix3-2)

        curr_tmp = curr_tmp + real(conjg(psi(ix1,ix2,ix3))*(zc_p1*(psi(ix1p1,ix2,ix3) &
            + psi(ix1,ix2p1,ix3) + psi(ix1,ix2,ix3p1)) &
            + zc_m1*(psi(ix1m1,ix2,ix3) + psi(ix1,ix2m1,ix3) &
            + psi(ix1,ix2,ix3m1)) &
            + zc_p2*(psi(ix1p2,ix2,ix3) + psi(ix1,ix2p2,ix3) &
            + psi(ix1,ix2,ix3p2)) &
            + zc_m2*(psi(ix1m2,ix2,ix3) + psi(ix1,ix2m2,ix3) &
            + psi(ix1,ix2,ix3m2)))) &
            + 3d0*avec*abs(psi(ix1,ix2,ix3))**2
      end do
    end do
  end do

  curr_tmp = curr_tmp*dx**3

  total_current = curr_tmp

end function total_current
!-------------------------------------------------------
real(8) function antisymmetry_error(psi)
  implicit none
  complex(8),intent(in) :: psi(0:nx-1,0:nx-1,0:nx-1)
  integer :: i, j, k
  real(8) :: err

  err = 0d0
  do i = 0, nx-1
    do j = 0, nx-1
      do k = 0, nx-1
        err = max(err, abs(psi(i,j,k) + psi(j,i,k)))
        err = max(err, abs(psi(i,j,k) + psi(i,k,j)))
        err = max(err, abs(psi(i,j,k) + psi(k,j,i)))
      end do
    end do
  end do
  antisymmetry_error = err

end function antisymmetry_error
!-------------------------------------------------------
subroutine finalize
  implicit none

  if (allocated(xn)) deallocate(xn)
  if (allocated(zpsi)) deallocate(zpsi)
  if (allocated(vpot_1d)) deallocate(vpot_1d)
  if (allocated(wpot_1d)) deallocate(wpot_1d)
  if (allocated(vpot)) deallocate(vpot)
  if (allocated(wpot)) deallocate(wpot)
  if (allocated(tot_pot)) deallocate(tot_pot)

end subroutine finalize
!-------------------------------------------------------
end program main
