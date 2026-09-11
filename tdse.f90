module global_variables
  implicit none
! Build:
!   gfortran -std=f95 -O2 tdse.f90 -o tdse
! Run:
!   ./tdse < input_tdse
! Output:
!   ground_state.log : CG iteration, energy, residual, norm, antisymmetry check
!   eigenstates.log  : convergence data for all computed eigenstates
!   eigenstate_currents.out : field-free particle velocity and charge current
!   current.out      : total single-run current, norm, and energy versus time
!   state_populations.out : raw projections and norm-corrected populations
!   state_manifold_populations.out : normalized degenerate-manifold populations
! Optional triplet mode also writes current_plus/minus/zero.out and the
! zero-field-subtracted even response in current_second_order.out.
! A single-run current contains all response orders; it is not, by itself, the
! shift current.  An even-in-field component requires separate +E0 and -E0 runs.
! math parameters
  real(8),parameter :: pi = 3.141592653589793238462643383279502884197d0
  complex(8),parameter :: zi = (0d0, 1d0)

! Physical constants (atomic units)
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: fs = 1d0/0.024189d0
  real(8),parameter :: bohr = 0.52917721067d0
! Atomic unit of electric field in V/m (CODATA value used for input conversion).
  real(8),parameter :: electric_field_au_Vpm = 5.14220674763d11

! Finite difference parameters (4th-order central stencils)
  real(8),parameter :: lc2 = -1d0/12d0, lc1 = 4d0/3d0, lc0 = -5d0/2d0
  real(8),parameter :: gc2 = -1d0/12d0, gc1 = 2d0/3d0

! Numerical parameters kept near the top for easy changes.
  integer,parameter :: output_stride = 1
  real(8),parameter :: cg_energy_tol = 1d-11
  real(8),parameter :: cg_residual_tol = 1d-9
! States separated by less than one microhartree are grouped.  This is far
! below the optical scale but above typical numerical splittings inside a
! converged degenerate multiplet.
  real(8),parameter :: degeneracy_energy_tol = 1d-6
  integer :: cg_max_iter = 2000
  integer :: num_eigenstates = 5
  logical :: require_converged_ground_state = .true.
  logical :: run_second_order_triplet = .false.
  logical :: run_field_scaling_check = .false.

  integer :: nx, nt
  real(8) :: dx, dt
  real(8) :: Tprop



! Material parameters
  real(8) :: lattice_constant
  real(8) :: bvc_lattice_constant
! Pair-interaction strength.  The default w0=0 is the intentional
! noninteracting baseline; an optional input record enables interacting studies.
  real(8) :: w0 = 0d0

! laser parameters
  real(8) :: E0, omega, Tpulse, phi_CEP

! grids
  real(8), allocatable :: xn(:)

! wavefunction
  complex(8), allocatable :: zpsi(:,:,:)
! eigenvectors(ix1,ix2,ix3,state_index), with state_index=1 the ground state
  real(8), allocatable :: eigenvalues(:), eigenstate_residuals(:)
  integer, allocatable :: eigenstate_iterations(:)
  logical, allocatable :: eigenstate_converged(:)
  complex(8), allocatable :: eigenvectors(:,:,:,:)

! potentials
  real(8), allocatable :: vpot_1d(:), wpot_1d(:)
  real(8), allocatable :: vpot(:,:,:), wpot(:,:,:), tot_pot(:,:,:)
! Field-free one-body eigenbasis used only to construct high-quality
! antisymmetric many-body initial guesses.
  real(8), allocatable :: onebody_values(:),onebody_vectors(:,:)

end module global_variables
!-------------------------------------------------------
program main
  use global_variables
  implicit none
  real(8) :: e0_gs, res_gs, asym_gs

  call initialize
  call compute_lowest_fermion_states(num_eigenstates, eigenvalues, &
      eigenvectors, eigenstate_residuals)
  call output_energy_gaps(num_eigenstates,eigenvalues)
  call output_eigenstate_currents(num_eigenstates, eigenvalues, eigenvectors)
  zpsi = eigenvectors(:,:,:,1)
  e0_gs = eigenvalues(1)
  res_gs = eigenstate_residuals(1)

  asym_gs = antisymmetry_error(zpsi)
  write(*,'(a,1pe16.8)') 'Ground-state energy        = ', e0_gs
  write(*,'(a,1pe16.8)') 'Ground-state residual norm = ', res_gs
  write(*,'(a,1pe16.8)') 'Ground-state norm          = ', wavefunction_norm(zpsi)
  write(*,'(a,1pe16.8)') 'Antisymmetry error         = ', asym_gs

  call check_current_operator(zpsi)

  if (require_converged_ground_state .and. .not.eigenstate_converged(1)) then
    write(*,'(a)') 'ERROR: ground state did not meet the strict residual target.'
    write(*,'(a)') 'Time propagation disabled by require_converged_ground_state.'
    call finalize
    stop 2
  end if

  if (run_second_order_triplet) then
    call propagate_second_order
  else
    call propagate_single_run
  end if
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
  allocate(eigenvalues(num_eigenstates))
  allocate(eigenstate_residuals(num_eigenstates))
  allocate(eigenstate_iterations(num_eigenstates))
  allocate(eigenstate_converged(num_eigenstates))
  allocate(eigenvectors(0:nx-1,0:nx-1,0:nx-1,num_eigenstates))

  call set_potentials
  call compute_onebody_basis
  call check_time_step

end subroutine initialize
!-------------------------------------------------------
subroutine read_input_parameters
  implicit none
  integer :: ios
  real(8) :: Tprop_fs
  real(8) :: E0_MVm, omega_ev, Tpulse_fs, phi_CEP_2pi

  read(*,*)lattice_constant, nx
  read(*,*)Tprop_fs, dt
  read(*,*)E0_MVm, omega_ev, Tpulse_fs, phi_CEP_2pi

! Optional trailing records preserve compatibility with the original input:
!   record 4: w0 [Hartree] (default 0, the noninteracting baseline)
!   record 5: run_triplet, run_scaling, require_converged, nstates, max_iter
  read(*,*,iostat=ios) w0
  if (ios /= 0) w0 = 0d0
  if (ios == 0) then
    read(*,*,iostat=ios) run_second_order_triplet,run_field_scaling_check, &
        require_converged_ground_state,num_eigenstates,cg_max_iter
  end if

! Input units are intentionally mixed for backward compatibility:
! lattice_constant and dt are in atomic units (bohr and atomic time), while
! Tprop_fs/Tpulse_fs are fs, E0_MVm is MV/m, and omega_ev is eV.
  if (lattice_constant <= 0d0) stop 'lattice_constant must be > 0 bohr.'
  if (nx < 5) stop 'nx must be >= 5 for the fourth-order finite differences.'
  if (dt <= 0d0) stop 'dt must be > 0 atomic units of time.'
  if (Tprop_fs <= 0d0) stop 'Tprop_fs must be > 0 fs.'
  if (Tpulse_fs <= 0d0) stop 'Tpulse_fs must be > 0 fs.'
  if (omega_ev <= 0d0) stop 'omega_ev must be > 0 eV.'
  if (num_eigenstates < 5) stop 'num_eigenstates must be at least 5.'
  if (cg_max_iter < 1) stop 'cg_max_iter must be positive.'
  if (run_field_scaling_check) run_second_order_triplet = .true.

  write(*,'(a)') 'Input-unit convention: lattice_constant [bohr], nx [grid points]'
  write(*,'(a)') '  Tprop_fs/Tpulse_fs [fs], dt [a.u. time], E0_MVm [MV/m],'
  write(*,'(a)') '  omega_ev [eV], phi_CEP_2pi [cycles].'
  write(*,*)'lattice_constant [bohr] = ', lattice_constant
  write(*,*)'nx = ', nx
  write(*,*)'Tprop_fs = ', Tprop_fs
  write(*,*)'requested dt [a.u. time] = ', dt
  write(*,*)'phi_CEP_2pi = ', phi_CEP_2pi
  Tprop = Tprop_fs*fs
  nt = max(1, nint(Tprop/dt))
  dt = Tprop/dble(nt)
  write(*,*)'refined dt [a.u. time] = ', dt
  write(*,*)'number of time steps nt = ', nt
  write(*,*)'Tprop and nt*dt [a.u. time] = ', Tprop, dble(nt)*dt

  E0 = E0_MVm*1d6/electric_field_au_Vpm
  omega = omega_ev*ev
  Tpulse = Tpulse_fs*fs
  phi_CEP = phi_CEP_2pi*2d0*pi

  bvc_lattice_constant = lattice_constant*3d0
  write(*,'(a)') 'Input field amplitude:'
  write(*,'(a,1pe20.12,a)') '  E0 = ',E0_MVm,' MV/m'
  write(*,'(a,1pe20.12,a)') '     = ',E0_MVm*1d6,' V/m'
  write(*,'(a,1pe20.12,a)') '     = ',E0,' a.u.'
  write(*,'(a,1pe20.12,a)') 'omega = ',omega_ev,' eV'
  write(*,'(a,1pe20.12,a)') '      = ',omega,' a.u.'
  write(*,'(a,1pe20.12,a)') 'Tpulse = ',Tpulse_fs,' fs'
  write(*,'(a,1pe20.12,a)') '       = ',Tpulse,' a.u.'
  write(*,'(a,1pe16.8)') 'ring length L=3a [bohr]  = ', bvc_lattice_constant
  write(*,'(a,1pe16.8)') 'pair strength w0 [a.u.]  = ', w0
  if (abs(w0) <= tiny(1d0)) then
    write(*,'(a)') 'Interaction: w0=0 (noninteracting baseline).'
  end if
  if (mod(nx,3) /= 0) then
    write(*,'(a)') 'WARNING:'
    write(*,'(a)') 'nx is not divisible by 3.'
    write(*,'(a)') 'The three-unit-cell translation symmetry is not represented exactly'// &
        ' on the discrete grid.'
    write(*,'(a)') 'Use nx = 3*m for production shift-current calculations.'
  end if
  write(*,'(a,i8)') 'num_eigenstates = ',num_eigenstates
  write(*,'(a,i8)') 'cg_max_iter = ',cg_max_iter
  write(*,'(a,1pe12.4)') 'Target residual norm = ',cg_residual_tol
  write(*,'(a,l1)') 'require_converged_ground_state = ', &
      require_converged_ground_state
  write(*,'(a,l1)') 'run_second_order_triplet = ',run_second_order_triplet
  write(*,'(a,l1)') 'run_field_scaling_check = ',run_field_scaling_check

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
subroutine compute_onebody_basis
  implicit none
  real(8),allocatable :: hmat(:,:)
  real(8) :: c0,c1,c2
  integer :: i

  allocate(hmat(0:nx-1,0:nx-1))
  allocate(onebody_values(0:nx-1),onebody_vectors(0:nx-1,0:nx-1))
  hmat=0d0
  c0=-0.5d0*lc0/dx**2
  c1=-0.5d0*lc1/dx**2
  c2=-0.5d0*lc2/dx**2
  do i=0,nx-1
    hmat(i,i)=hmat(i,i)+c0+vpot_1d(i)
    hmat(i,ipbc(i+1))=hmat(i,ipbc(i+1))+c1
    hmat(i,ipbc(i-1))=hmat(i,ipbc(i-1))+c1
    hmat(i,ipbc(i+2))=hmat(i,ipbc(i+2))+c2
    hmat(i,ipbc(i-2))=hmat(i,ipbc(i-2))+c2
  end do
  call jacobi_real_symmetric(hmat,onebody_values,onebody_vectors)
  call sort_onebody_eigenpairs(onebody_values,onebody_vectors)
  deallocate(hmat)

end subroutine compute_onebody_basis
!-------------------------------------------------------
subroutine jacobi_real_symmetric(matrix,values,vectors)
  implicit none
  real(8),intent(inout) :: matrix(0:nx-1,0:nx-1)
  real(8),intent(out) :: values(0:nx-1),vectors(0:nx-1,0:nx-1)
  integer :: sweep,p,q,k
  real(8) :: app,aqq,apq,tau,t,c,s,akp,akq,vkp,vkq,max_offdiag

  vectors=0d0
  do p=0,nx-1
    vectors(p,p)=1d0
  end do
  do sweep=1,100
    max_offdiag=0d0
    do p=0,nx-2
      do q=p+1,nx-1
        apq=matrix(p,q)
        max_offdiag=max(max_offdiag,abs(apq))
        if (abs(apq) <= 1d-15) cycle
        app=matrix(p,p)
        aqq=matrix(q,q)
        tau=(aqq-app)/(2d0*apq)
        if (tau >= 0d0) then
          t=1d0/(tau+sqrt(1d0+tau*tau))
        else
          t=-1d0/(-tau+sqrt(1d0+tau*tau))
        end if
        c=1d0/sqrt(1d0+t*t)
        s=t*c
        matrix(p,p)=app-t*apq
        matrix(q,q)=aqq+t*apq
        matrix(p,q)=0d0
        matrix(q,p)=0d0
        do k=0,nx-1
          if (k /= p .and. k /= q) then
            akp=matrix(k,p)
            akq=matrix(k,q)
            matrix(k,p)=c*akp-s*akq
            matrix(p,k)=matrix(k,p)
            matrix(k,q)=s*akp+c*akq
            matrix(q,k)=matrix(k,q)
          end if
          vkp=vectors(k,p)
          vkq=vectors(k,q)
          vectors(k,p)=c*vkp-s*vkq
          vectors(k,q)=s*vkp+c*vkq
        end do
      end do
    end do
    if (max_offdiag < 1d-14) exit
  end do
  do p=0,nx-1
    values(p)=matrix(p,p)
  end do
  if (max_offdiag >= 1d-12) then
    write(*,'(a,1pe16.8)') 'Warning: one-body Jacobi off-diagonal norm = ', &
        max_offdiag
  end if

end subroutine jacobi_real_symmetric
!-------------------------------------------------------
subroutine sort_onebody_eigenpairs(values,vectors)
  implicit none
  real(8),intent(inout) :: values(0:nx-1),vectors(0:nx-1,0:nx-1)
  integer :: i,j
  real(8) :: value_work
  real(8),allocatable :: vector_work(:)

  allocate(vector_work(0:nx-1))
  do i=0,nx-2
    do j=i+1,nx-1
      if (values(j) < values(i)) then
        value_work=values(i)
        values(i)=values(j)
        values(j)=value_work
        vector_work=vectors(:,i)
        vectors(:,i)=vectors(:,j)
        vectors(:,j)=vector_work
      end if
    end do
  end do
  deallocate(vector_work)

end subroutine sort_onebody_eigenpairs
!-------------------------------------------------------
subroutine check_time_step
  implicit none
  real(8) :: hmax_est, avec_max_est, pmax_one_est

! This is a conservative spectral-scale estimate, not a rigorous RK4 bound.
! The fourth-order Laplacian contributes at most 8/dx**2 for three particles.
! For H(A)=H(0)+A*sum_i p_i+3*A**2/2, use |A|<=|E0/omega| and
! a triangle-inequality estimate of the fourth-order derivative spectrum.
  avec_max_est = abs(E0/omega)
  pmax_one_est = 2d0*(abs(gc1)+abs(gc2))/dx
  hmax_est = 8d0/dx**2 + maxval(abs(tot_pot)) &
      + 3d0*avec_max_est*pmax_one_est + 1.5d0*avec_max_est**2
  write(*,'(a,1pe12.4)') 'RK4 field-inclusive dt*H scale estimate = ',dt*hmax_est
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
subroutine compute_lowest_fermion_states(nstates, eigenvalues_out, &
    eigenvectors_out, residuals_out)
  implicit none
  integer,intent(in) :: nstates
  real(8),intent(out) :: eigenvalues_out(nstates), residuals_out(nstates)
  complex(8),intent(out) :: eigenvectors_out(0:nx-1,0:nx-1,0:nx-1,nstates)
  complex(8),allocatable :: psi(:,:,:)
  integer :: istate

  if (nx < 4 .and. nstates >= 4) then
    stop 'At least nx=4 is required for four antisymmetric states.'
  end if

  allocate(psi(0:nx-1,0:nx-1,0:nx-1))
  eigenvectors_out = (0d0,0d0)

  open(20,file='ground_state.log',status='replace')
  write(20,'(a)') '# iter energy residual_norm norm antisymmetry_error'
  open(21,file='eigenstates.log',status='replace')
  write(21,'(a)') '# state iter energy residual_norm norm antisymmetry_error'// &
      ' best_residual residual_increased'

  do istate = 1, nstates
    call initialize_eigenstate_guess(istate, psi)
    call orthogonalize_against_states(psi, eigenvectors_out, istate-1)
    call normalize(psi)
    call minimize_projected_rayleigh(istate, psi, eigenvectors_out, &
        istate-1, eigenvalues_out(istate), residuals_out(istate), &
        eigenstate_iterations(istate),eigenstate_converged(istate))
    eigenvectors_out(:,:,:,istate) = psi
  end do

  close(20)
  close(21)
  deallocate(psi)

  call sort_eigenpairs(nstates, eigenvalues_out, eigenvectors_out, residuals_out)
  call verify_eigenstates(nstates, eigenvalues_out, eigenvectors_out, residuals_out)

end subroutine compute_lowest_fermion_states
!-------------------------------------------------------
subroutine initialize_eigenstate_guess(istate, psi)
  implicit none
  integer,intent(in) :: istate
  complex(8),intent(out) :: psi(0:nx-1,0:nx-1,0:nx-1)
  integer :: ix1, ix2, ix3
  integer :: modes(3)
  complex(8) :: a(3), b(3), c(3)

! Select a different low-energy one-body configuration for every state.
! At w0=0 these determinants are exact discrete many-body eigenstates; for
! w0/=0 they remain physically motivated, linearly independent CG seeds.
  call select_slater_configuration(istate,modes)

  do ix1 = 0, nx-1
    do ix2 = 0, nx-1
      do ix3 = 0, nx-1
        a = cmplx(onebody_vectors(ix1,modes),0d0,kind(zi))
        b = cmplx(onebody_vectors(ix2,modes),0d0,kind(zi))
        c = cmplx(onebody_vectors(ix3,modes),0d0,kind(zi))
        psi(ix1,ix2,ix3) = determinant3(a,b,c)
      end do
    end do
  end do

  call antisymmetrize(psi)

end subroutine initialize_eigenstate_guess
!-------------------------------------------------------
subroutine select_slater_configuration(rank,modes)
  implicit none
  integer,intent(in) :: rank
  integer,intent(out) :: modes(3)
  integer :: i,j,k,nconfig,index,best_index,selection
  integer,allocatable :: configurations(:,:)
  real(8),allocatable :: configuration_energy(:)
  logical,allocatable :: selected(:)
  real(8) :: best_energy

  nconfig=nx*(nx-1)*(nx-2)/6
  if (rank > nconfig) stop 'Too many eigenstates for antisymmetric grid space.'
  allocate(configurations(3,nconfig),configuration_energy(nconfig))
  allocate(selected(nconfig))
  index=0
  do i=0,nx-3
    do j=i+1,nx-2
      do k=j+1,nx-1
        index=index+1
        configurations(:,index)=(/i,j,k/)
        configuration_energy(index)=onebody_values(i)+onebody_values(j)+ &
            onebody_values(k)
      end do
    end do
  end do
  selected=.false.
  best_index=1
  do selection=1,rank
    best_energy=huge(1d0)
    do index=1,nconfig
      if (.not.selected(index) .and. configuration_energy(index)<best_energy) then
        best_energy=configuration_energy(index)
        best_index=index
      end if
    end do
    selected(best_index)=.true.
  end do
  modes=configurations(:,best_index)
  deallocate(configurations,configuration_energy,selected)

end subroutine select_slater_configuration
!-------------------------------------------------------
subroutine orthogonalize_against_states(psi, states, nstates_done)
  implicit none
  complex(8),intent(inout) :: psi(0:nx-1,0:nx-1,0:nx-1)
  complex(8),intent(in) :: states(0:nx-1,0:nx-1,0:nx-1,*)
  integer,intent(in) :: nstates_done
  integer :: ipass, jstate
  complex(8) :: overlap

! Two-pass modified Gram-Schmidt is robust also for nearly degenerate states.
  do ipass = 1, 2
    do jstate = 1, nstates_done
      overlap = inner_product(states(:,:,:,jstate),psi)
      psi = psi - overlap*states(:,:,:,jstate)
    end do
  end do
  call antisymmetrize(psi)

end subroutine orthogonalize_against_states
!-------------------------------------------------------
subroutine minimize_projected_rayleigh(istate, psi, states, nstates_done, &
    energy, residual_norm,iterations_used,converged)
  implicit none
  integer,intent(in) :: istate, nstates_done
  complex(8),intent(inout) :: psi(0:nx-1,0:nx-1,0:nx-1)
  complex(8),intent(in) :: states(0:nx-1,0:nx-1,0:nx-1,*)
  real(8),intent(out) :: energy, residual_norm
  integer,intent(out) :: iterations_used
  logical,intent(out) :: converged
  complex(8),allocatable :: hpsi(:,:,:), residual(:,:,:), xi(:,:,:)
  complex(8),allocatable :: direction(:,:,:), direction_old(:,:,:), hp(:,:,:)
  complex(8) :: overlap
  real(8) :: xixi, xixi_old, gamma, direction_norm, energy_old
  real(8) :: best_residual,previous_residual
  integer :: iter, stagnation_count, increase_count
  logical :: residual_increased

  allocate(hpsi(0:nx-1,0:nx-1,0:nx-1))
  allocate(residual(0:nx-1,0:nx-1,0:nx-1))
  allocate(xi(0:nx-1,0:nx-1,0:nx-1))
  allocate(direction(0:nx-1,0:nx-1,0:nx-1))
  allocate(direction_old(0:nx-1,0:nx-1,0:nx-1))
  allocate(hp(0:nx-1,0:nx-1,0:nx-1))

  direction_old = (0d0,0d0)
  xixi_old = 1d0
  converged = .false.
  iterations_used = 0
  stagnation_count = 0
  increase_count = 0

  call apply_hamiltonian(psi,hpsi,0d0)
  energy = real(inner_product(psi,hpsi))
  residual = hpsi - energy*psi
  residual_norm = wavefunction_norm(residual)
  best_residual = residual_norm
  previous_residual = residual_norm
  write(21,'(2i8,5(1x,1pe20.12),1x,l1)') istate,0,energy,residual_norm, &
      wavefunction_norm(psi),antisymmetry_error(psi),best_residual,.false.
  if (istate == 1) write(20,'(i8,4(1x,1pe20.12))') 0,energy, &
      residual_norm,wavefunction_norm(psi),antisymmetry_error(psi)

  do iter = 1, cg_max_iter
    if (residual_norm < cg_residual_tol) then
      converged = .true.
      exit
    end if

! Jacobi preconditioning uses the field-free Hamiltonian diagonal.  The
! positive floor avoids amplifying components near a shifted diagonal zero.
    call apply_diagonal_preconditioner(residual,xi,energy)
    xi = -xi
    call orthogonalize_against_states(xi,states,nstates_done)
    overlap = inner_product(psi,xi)
    xi = xi - overlap*psi
    call antisymmetrize(xi)
    xixi = -real(inner_product(residual,xi))
    if (xixi <= 1d-28) exit

! Periodic restart limits loss of conjugacy in nearly degenerate subspaces.
    if (iter == 1 .or. modulo(iter-1,25) == 0) then
      gamma = 0d0
    else
      gamma = xixi/xixi_old
    end if
    direction = xi + gamma*direction_old
    call orthogonalize_against_states(direction,states,nstates_done)
    overlap = inner_product(psi,direction)
    direction = direction - overlap*psi
    call antisymmetrize(direction)
    direction_norm = wavefunction_norm(direction)
    if (direction_norm <= 1d-14) then
      direction = xi
      direction_norm = wavefunction_norm(direction)
      if (direction_norm <= 1d-14) exit
    end if
    direction_old = direction
    direction = direction/direction_norm

    call apply_hamiltonian(direction,hp,0d0)
    energy_old = energy
    call rayleigh_ritz_update(psi,direction,hpsi,hp,energy)
    call orthogonalize_against_states(psi,states,nstates_done)
    call normalize(psi)

    call apply_hamiltonian(psi,hpsi,0d0)
    energy = real(inner_product(psi,hpsi))
    residual = hpsi - energy*psi
    residual_norm = wavefunction_norm(residual)
    residual_increased = residual_norm > previous_residual*(1d0+1d-12)
    if (residual_increased) increase_count = increase_count+1
    if (residual_norm < best_residual*(1d0-1d-4)) then
      best_residual = residual_norm
      stagnation_count = 0
    else
      stagnation_count = stagnation_count+1
    end if
    write(21,'(2i8,5(1x,1pe20.12),1x,l1)') istate,iter,energy, &
        residual_norm,wavefunction_norm(psi),antisymmetry_error(psi), &
        best_residual,residual_increased
    if (istate == 1) write(20,'(i8,4(1x,1pe20.12))') iter,energy, &
        residual_norm,wavefunction_norm(psi),antisymmetry_error(psi)

    xixi_old = xixi
    previous_residual = residual_norm
    iterations_used = iter
    if (residual_norm < cg_residual_tol .and. &
        abs(energy-energy_old) < cg_energy_tol) then
      converged = .true.
      exit
    end if
  end do

  iterations_used = min(iterations_used,cg_max_iter)
  write(*,'(a,i4)') 'Eigenstate solver state   : ',istate
  write(*,'(a,1pe16.8)') '  Target residual norm   : ',cg_residual_tol
  write(*,'(a,1pe16.8)') '  Final residual norm    : ',residual_norm
  write(*,'(a,i8)') '  Iterations             : ',iterations_used
  if (converged) then
    write(*,'(a)') '  Status                 : CONVERGED'
  else
    write(*,'(a)') '  Status                 : NOT CONVERGED'
  end if
  write(*,'(a,i8)') '  Residual increases     : ',increase_count
  if (stagnation_count >= 50) then
    write(*,'(a,i8,a)') '  Diagnostic             : residual stagnated for ', &
        stagnation_count,' final iterations.'
  end if

  if (.not.converged .and. residual_norm >= cg_residual_tol) then
    write(*,'(a,i4,a,i8,a,1pe16.8,a,1pe16.8)') &
        'Warning: eigenstate ',istate,' did not converge in ', &
        min(iter,cg_max_iter), &
        ' iterations; energy=',energy,', residual=',residual_norm
  end if

  deallocate(hpsi,residual,xi,direction,direction_old,hp)

end subroutine minimize_projected_rayleigh
!-------------------------------------------------------
subroutine apply_diagonal_preconditioner(residual,preconditioned,energy)
  implicit none
  complex(8),intent(in) :: residual(0:nx-1,0:nx-1,0:nx-1)
  complex(8),intent(out) :: preconditioned(0:nx-1,0:nx-1,0:nx-1)
  real(8),intent(in) :: energy
  integer :: ix1,ix2,ix3
  real(8) :: diagonal,denominator,denominator_floor

  denominator_floor=0.1d0/dx**2
  do ix1=0,nx-1
    do ix2=0,nx-1
      do ix3=0,nx-1
        diagonal=3d0*(-0.5d0*lc0/dx**2)+tot_pot(ix1,ix2,ix3)
        denominator=max(abs(diagonal-energy),denominator_floor)
        preconditioned(ix1,ix2,ix3)=residual(ix1,ix2,ix3)/denominator
      end do
    end do
  end do

end subroutine apply_diagonal_preconditioner
!-------------------------------------------------------
subroutine sort_eigenpairs(nstates, values, states, residuals)
  implicit none
  integer,intent(in) :: nstates
  real(8),intent(inout) :: values(nstates), residuals(nstates)
  complex(8),intent(inout) :: states(0:nx-1,0:nx-1,0:nx-1,nstates)
  complex(8),allocatable :: work(:,:,:)
  real(8) :: value_work, residual_work
  integer :: iteration_work
  logical :: converged_work
  integer :: i, j

  allocate(work(0:nx-1,0:nx-1,0:nx-1))
  do i = 1, nstates-1
    do j = i+1, nstates
      if (values(j) < values(i)) then
        value_work = values(i)
        values(i) = values(j)
        values(j) = value_work
        residual_work = residuals(i)
        residuals(i) = residuals(j)
        residuals(j) = residual_work
        iteration_work = eigenstate_iterations(i)
        eigenstate_iterations(i) = eigenstate_iterations(j)
        eigenstate_iterations(j) = iteration_work
        converged_work = eigenstate_converged(i)
        eigenstate_converged(i) = eigenstate_converged(j)
        eigenstate_converged(j) = converged_work
        work = states(:,:,:,i)
        states(:,:,:,i) = states(:,:,:,j)
        states(:,:,:,j) = work
      end if
    end do
  end do
  deallocate(work)

end subroutine sort_eigenpairs
!-------------------------------------------------------
subroutine verify_eigenstates(nstates, values, states, residuals)
  implicit none
  integer,intent(in) :: nstates
  real(8),intent(inout) :: values(nstates), residuals(nstates)
  complex(8),intent(in) :: states(0:nx-1,0:nx-1,0:nx-1,nstates)
  complex(8),allocatable :: hpsi(:,:,:), residual(:,:,:)
  complex(8) :: overlap, target
  real(8) :: max_orth_error
  integer :: i, j

  allocate(hpsi(0:nx-1,0:nx-1,0:nx-1))
  allocate(residual(0:nx-1,0:nx-1,0:nx-1))
  max_orth_error = 0d0
  do i = 1, nstates
    do j = 1, nstates
      target = (0d0,0d0)
      if (i == j) target = (1d0,0d0)
      overlap = inner_product(states(:,:,:,i),states(:,:,:,j))
      max_orth_error = max(max_orth_error,abs(overlap-target))
    end do
    call apply_hamiltonian(states(:,:,:,i),hpsi,0d0)
    values(i) = real(inner_product(states(:,:,:,i),hpsi))
    residual = hpsi - values(i)*states(:,:,:,i)
    residuals(i) = wavefunction_norm(residual)
    eigenstate_converged(i) = residuals(i) < cg_residual_tol
  end do

  write(*,'(a)') 'Lowest antisymmetric eigenstates:'
  do i = 1, nstates
    write(*,'(a,i2,3(a,1pe16.8),a,i8,a,l1)') ' state=',i-1,' energy=',values(i), &
        ' residual=',residuals(i),' antisymmetry_error=', &
        antisymmetry_error(states(:,:,:,i)),' iterations=', &
        eigenstate_iterations(i),' converged=',eigenstate_converged(i)
  end do
  write(*,'(a,l2)') ' Energies in ascending order = ', &
      all(values(2:nstates) >= values(1:nstates-1))
  write(*,'(a,1pe16.8)') ' Maximum orthonormality error = ',max_orth_error

  deallocate(hpsi,residual)

end subroutine verify_eigenstates
!-------------------------------------------------------
subroutine output_energy_gaps(nstates,values)
  implicit none
  integer,intent(in) :: nstates
  real(8),intent(in) :: values(nstates)
  integer :: i
  real(8) :: gap,gap_ev,detuning_ev

  open(23,file='eigenstate_energies.out',status='replace')
  write(23,'(a)') '# state E_n(Ha) DeltaE(Ha) DeltaE(eV) detuning(eV)'
  write(*,'(a)') 'Eigenstate energies and laser detunings:'
  write(*,'(a)') ' state       E_n(Ha)        DeltaE(Ha)'// &
      '      DeltaE(eV)     detuning(eV)'
  do i=1,nstates
    gap=values(i)-values(1)
    gap_ev=gap/ev
    detuning_ev=(gap-omega)/ev
    write(23,'(i8,4(1x,1pe20.12))') i-1,values(i),gap,gap_ev,detuning_ev
    write(*,'(i6,4(1x,1pe16.8))') i-1,values(i),gap,gap_ev,detuning_ev
  end do
  close(23)

end subroutine output_energy_gaps
!-------------------------------------------------------
subroutine output_eigenstate_currents(nstates, values, states)
  implicit none
  integer,intent(in) :: nstates
  real(8),intent(in) :: values(nstates)
  complex(8),intent(in) :: states(0:nx-1,0:nx-1,0:nx-1,nstates)
  real(8) :: particle_velocity, charge_current
  integer :: istate

! The eigenstates belong to the field-free Hamiltonian, so A=0 here.
  open(22,file='eigenstate_currents.out',status='replace')
  write(22,'(a)') '# Field-free normalized expectations; atomic units.'
  write(22,'(a)') '# columns: state energy <sum_i(p_i+A)> -<sum_i(p_i+A)>/L'
  write(*,'(a)') 'Field-free current expectation values:'
  do istate = 1, nstates
    particle_velocity = total_particle_velocity(states(:,:,:,istate),0d0)
    charge_current = charge_current_density(states(:,:,:,istate),0d0)
    write(22,'(i8,3(1x,1pe20.12))') istate,values(istate), &
        particle_velocity,charge_current
    write(*,'(a,i2,3(a,1pe16.8))') ' state=',istate, &
        ' energy=',values(istate),' particle velocity=',particle_velocity, &
        ' charge current density=',charge_current
  end do
  close(22)

end subroutine output_eigenstate_currents
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
subroutine propagate_single_run
  implicit none
  real(8),allocatable :: current_t(:)

  allocate(current_t(0:nt))
  call propagate_trajectory(E0,'current.out',.true.,current_t)
  if (abs(E0) <= tiny(1d0)) then
    write(*,'(a)') 'Zero-field numerical-current diagnostic:'
    write(*,'(a,1pe16.8)') '  max_t |J0(t)| = ',maxval(abs(current_t))
    write(*,'(a,1pe16.8)') '  RMS[J0]       = ', &
        sqrt(sum(current_t*current_t)/dble(nt+1))
  end if
  deallocate(current_t)

end subroutine propagate_single_run
!-------------------------------------------------------
subroutine propagate_second_order
  implicit none
  real(8),allocatable :: jplus(:),jminus(:),jzero(:)
  real(8),allocatable :: jplus_half(:),jminus_half(:),jzero_half(:)
  real(8),allocatable :: response(:),response_half(:)

  allocate(jplus(0:nt),jminus(0:nt),jzero(0:nt),response(0:nt))
! Every trajectory is reset to the same already-computed ground-state vector.
  call propagate_trajectory(E0,'current_plus.out',.true.,jplus)
  call propagate_trajectory(-E0,'current_minus.out',.false.,jminus)
  call propagate_trajectory(0d0,'current_zero.out',.false.,jzero)
  call write_second_order_file('current_second_order.out',E0,jplus,jminus, &
      jzero,response)
  call report_zero_field_floor(jzero,response,E0)

  if (run_field_scaling_check) then
    allocate(jplus_half(0:nt),jminus_half(0:nt),jzero_half(0:nt))
    allocate(response_half(0:nt))
    call propagate_trajectory(0.5d0*E0,'current_plus_half.out',.false., &
        jplus_half)
    call propagate_trajectory(-0.5d0*E0,'current_minus_half.out',.false., &
        jminus_half)
    call propagate_trajectory(0d0,'current_zero_half.out',.false.,jzero_half)
    call write_second_order_file('current_second_order_half.out',0.5d0*E0, &
        jplus_half,jminus_half,jzero_half,response_half)
    call write_scaling_comparison(response,response_half,E0)
    deallocate(jplus_half,jminus_half,jzero_half,response_half)
  end if
  deallocate(jplus,jminus,jzero,response)

end subroutine propagate_second_order
!-------------------------------------------------------
subroutine propagate_trajectory(field_amplitude,current_file,write_pop,current_t)
  implicit none
  real(8),intent(in) :: field_amplitude
  character(*),intent(in) :: current_file
  logical,intent(in) :: write_pop
  real(8),intent(out) :: current_t(0:nt)
  integer :: it,istate,imanifold,nmanifolds
  integer,allocatable :: manifold(:)
  real(8) :: t,avec,norm_squared,energy,particle_velocity
  real(8) :: raw_sum,normalized_sum
  real(8),allocatable :: raw_population(:),normalized_population(:)
  real(8),allocatable :: manifold_population(:)
  complex(8),allocatable :: hpsi(:,:,:),amplitude(:)

  allocate(hpsi(0:nx-1,0:nx-1,0:nx-1))
  allocate(amplitude(num_eigenstates),raw_population(num_eigenstates))
  allocate(normalized_population(num_eigenstates),manifold(num_eigenstates))
  call assign_energy_manifolds(manifold,nmanifolds)
  allocate(manifold_population(nmanifolds))
  zpsi = eigenvectors(:,:,:,1)

  open(30,file=current_file,status='replace')
  write(30,'(a)') '# Physical charge current density for one trajectory.'
  write(30,'(a)') '# Atomic units; q=-1; columns: t A particle_velocity'// &
      ' charge_current_density norm_squared normalized_energy'
  if (write_pop) then
    open(31,file='state_populations.out',status='replace')
    write(31,'(a)') '# columns: t norm, then for each state: Re(overlap)'// &
        ' Im(overlap) raw_overlap_squared normalized_population'
    write(31,'(a)') '# normalized_population=|<n|psi>|^2/<psi|psi>.'
    write(31,'(a)') '# final columns: raw_population_sum'// &
        ' normalized_population_sum.'
    open(32,file='state_manifold_populations.out',status='replace')
    write(32,'(a,1pe12.4)') '# Energy grouping tolerance [Hartree] = ', &
        degeneracy_energy_tol
    write(32,'(a)') '# columns: t(a.u.) norm manifold_0 manifold_1 ...'// &
        ' (normalized populations)'
  end if

  do it = 0,nt
    if (mod(it,max(1,nt/100)) == 0) write(*,'(a,a,a,i8)') &
        'trajectory ',trim(current_file),' it = ',it
    t = dble(it)*dt
    avec = vector_potential_for_field(t,field_amplitude)
    norm_squared = real(inner_product(zpsi,zpsi))
    if (norm_squared <= 0d0) stop 'TDSE state has non-positive norm squared.'
    call apply_hamiltonian(zpsi,hpsi,avec)
    energy = real(inner_product(zpsi,hpsi))/norm_squared
    particle_velocity = total_particle_velocity(zpsi,avec)
    current_t(it) = charge_current_density(zpsi,avec)
    write(30,"(999e26.16e3)") t,avec,particle_velocity,current_t(it), &
        norm_squared,energy
    if (write_pop .and. modulo(it,output_stride) == 0) then
      do istate=1,num_eigenstates
        amplitude(istate)=inner_product(eigenvectors(:,:,:,istate),zpsi)
        raw_population(istate)=abs(amplitude(istate))**2
      end do
      normalized_population=raw_population/norm_squared
      raw_sum=sum(raw_population)
      normalized_sum=sum(normalized_population)
      write(31,"(999e26.16e3)") t,norm_squared, &
          (real(amplitude(istate)),aimag(amplitude(istate)), &
          raw_population(istate),normalized_population(istate), &
          istate=1,num_eigenstates),raw_sum,normalized_sum
      manifold_population=0d0
      do istate=1,num_eigenstates
        imanifold=manifold(istate)
        manifold_population(imanifold)=manifold_population(imanifold)+ &
            normalized_population(istate)
      end do
      write(32,"(999e26.16e3)") t,norm_squared,manifold_population
    end if
    if (it < nt) call rk4_step(t,dt,field_amplitude)
  end do
  close(30)
  if (write_pop) then
    close(31)
    close(32)
  end if
  deallocate(hpsi,amplitude,raw_population,normalized_population)
  deallocate(manifold,manifold_population)

end subroutine propagate_trajectory
!-------------------------------------------------------
subroutine assign_energy_manifolds(manifold,nmanifolds)
  implicit none
  integer,intent(out) :: manifold(num_eigenstates),nmanifolds
  integer :: i

  nmanifolds=1
  manifold(1)=1
  do i=2,num_eigenstates
    if (abs(eigenvalues(i)-eigenvalues(i-1)) >= degeneracy_energy_tol) &
        nmanifolds=nmanifolds+1
    manifold(i)=nmanifolds
  end do
  write(*,'(a,i8,a,i8,a)') 'Grouped ',num_eigenstates,' states into ', &
      nmanifolds,' manifolds.'
  do i=1,num_eigenstates
    write(*,'(a,i4,a,i4)') '  state ',i-1,' -> manifold ',manifold(i)-1
  end do

end subroutine assign_energy_manifolds
!-------------------------------------------------------
subroutine write_second_order_file(filename,field_amplitude,jplus,jminus, &
    jzero,response)
  implicit none
  character(*),intent(in) :: filename
  real(8),intent(in) :: field_amplitude
  real(8),intent(in) :: jplus(0:nt),jminus(0:nt),jzero(0:nt)
  real(8),intent(out) :: response(0:nt)
  integer :: it
  real(8) :: t,jeven,jodd,scaled

  open(33,file=filename,status='replace')
  write(33,'(a)') '# Physical charge current density, atomic units.'
  write(33,'(a)') '# columns: t J_plus J_minus J_zero J_even'// &
      ' J_even_induced J_even_induced_over_E0_squared J_odd'
  do it=0,nt
    t=dble(it)*dt
    jeven=0.5d0*(jplus(it)+jminus(it))
    jodd=0.5d0*(jplus(it)-jminus(it))
    response(it)=jeven-jzero(it)
    scaled=0d0
    if (abs(field_amplitude) > 0d0) scaled=response(it)/field_amplitude**2
    write(33,"(999e26.16e3)") t,jplus(it),jminus(it),jzero(it), &
        jeven,response(it),scaled,jodd
  end do
  close(33)

end subroutine write_second_order_file
!-------------------------------------------------------
subroutine report_zero_field_floor(jzero,response,field_amplitude)
  implicit none
  real(8),intent(in) :: jzero(0:nt),response(0:nt),field_amplitude
  real(8) :: max_zero,rms_zero,max_signal,ratio

  max_zero=maxval(abs(jzero))
  rms_zero=sqrt(sum(jzero*jzero)/dble(nt+1))
  max_signal=maxval(abs(response))
  ratio=huge(1d0)
  if (max_zero > 0d0) ratio=max_signal/max_zero
  write(*,'(a)') 'Zero-field numerical-current diagnostic:'
  write(*,'(a,1pe16.8)') '  max_t |J0(t)| = ',max_zero
  write(*,'(a,1pe16.8)') '  RMS[J0]       = ',rms_zero
  write(*,'(a,1pe16.8)') '  max |J_even_induced| = ',max_signal
  write(*,'(a,1pe16.8)') '  signal / max|J0|     = ',ratio
  write(*,'(a,1pe16.8)') '  field amplitude used = ',field_amplitude
  write(*,'(a)') '  This ratio is a numerical-floor diagnostic, not an error bar.'

end subroutine report_zero_field_floor
!-------------------------------------------------------
subroutine write_scaling_comparison(response,response_half,field_amplitude)
  implicit none
  real(8),intent(in) :: response(0:nt),response_half(0:nt),field_amplitude
  integer :: it
  real(8) :: scaled_full,scaled_half,relative_difference,denominator
  real(8) :: maximum_relative_difference

  open(34,file='second_order_scaling.out',status='replace')
  write(34,'(a)') '# Compare J_even_induced/E^2 at E0 and E0/2.'
  write(34,'(a)') '# columns: t scaled_E0 scaled_E0_over_2 relative_difference'
  maximum_relative_difference=0d0
  do it=0,nt
    scaled_full=0d0
    scaled_half=0d0
    if (abs(field_amplitude) > 0d0) then
      scaled_full=response(it)/field_amplitude**2
      scaled_half=response_half(it)/(0.5d0*field_amplitude)**2
    end if
    denominator=max(abs(scaled_full),abs(scaled_half),tiny(1d0))
    relative_difference=abs(scaled_full-scaled_half)/denominator
    maximum_relative_difference=max(maximum_relative_difference, &
        relative_difference)
    write(34,"(999e26.16e3)") dble(it)*dt,scaled_full,scaled_half, &
        relative_difference
  end do
  close(34)
  write(*,'(a,1pe16.8)') 'Maximum pointwise second-order scaling difference = ', &
      maximum_relative_difference

end subroutine write_scaling_comparison
!-------------------------------------------------------
subroutine rk4_step(t, h, field_amplitude)
  implicit none
  real(8),intent(in) :: t, h, field_amplitude
  complex(8),allocatable :: y0(:,:,:), yt(:,:,:), k1(:,:,:), k2(:,:,:)
  complex(8),allocatable :: k3(:,:,:), k4(:,:,:)

  allocate(y0(0:nx-1,0:nx-1,0:nx-1))
  allocate(yt(0:nx-1,0:nx-1,0:nx-1))
  allocate(k1(0:nx-1,0:nx-1,0:nx-1))
  allocate(k2(0:nx-1,0:nx-1,0:nx-1))
  allocate(k3(0:nx-1,0:nx-1,0:nx-1))
  allocate(k4(0:nx-1,0:nx-1,0:nx-1))

  y0 = zpsi
  call tdse_rhs(y0, k1, t, field_amplitude)
  yt = y0 + 0.5d0*h*k1
  call antisymmetrize(yt)
  call tdse_rhs(yt, k2, t+0.5d0*h, field_amplitude)
  yt = y0 + 0.5d0*h*k2
  call antisymmetrize(yt)
  call tdse_rhs(yt, k3, t+0.5d0*h, field_amplitude)
  yt = y0 + h*k3
  call antisymmetrize(yt)
  call tdse_rhs(yt, k4, t+h, field_amplitude)

  zpsi = y0 + h*(k1 + 2d0*k2 + 2d0*k3 + k4)/6d0
  call antisymmetrize(zpsi)

  deallocate(y0, yt, k1, k2, k3, k4)

end subroutine rk4_step
!-------------------------------------------------------
subroutine tdse_rhs(psi, rhs, t, field_amplitude)
  implicit none
  complex(8),intent(in) :: psi(0:nx-1,0:nx-1,0:nx-1)
  complex(8),intent(out) :: rhs(0:nx-1,0:nx-1,0:nx-1)
  real(8),intent(in) :: t, field_amplitude

  call apply_hamiltonian(psi, rhs, vector_potential_for_field(t,field_amplitude))
  rhs = -zi*rhs

end subroutine tdse_rhs
!-------------------------------------------------------
real(8) function vector_potential(t)
  implicit none
  real(8),intent(in) :: t
  vector_potential=vector_potential_for_field(t,E0)

end function vector_potential
!-------------------------------------------------------
real(8) function vector_potential_for_field(t,field_amplitude)
  implicit none
  real(8),intent(in) :: t,field_amplitude
  real(8) :: env

  if (t < 0d0 .or. t > Tpulse) then
    vector_potential_for_field = 0d0
  else
    env = sin(pi*t/Tpulse)**4
! A(t) is chosen so that the field is approximately E(t)=-dA/dt for a
! slowly varying envelope; the exact A(t) is what enters the Hamiltonian.
    vector_potential_for_field = -(field_amplitude/omega)*env* &
        sin(omega*(t-0.5d0*tpulse) + phi_CEP)
  end if

end function vector_potential_for_field
!-------------------------------------------------------
real(8) function particle_velocity_numerator(psi, avec)
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

  particle_velocity_numerator = curr_tmp

end function particle_velocity_numerator
!-------------------------------------------------------
real(8) function total_particle_velocity(psi, avec)
  implicit none
  complex(8),intent(in) :: psi(0:nx-1,0:nx-1,0:nx-1)
  real(8),intent(in) :: avec
  real(8) :: norm_squared

! Normalized expectation of the total mechanical velocity,
! <sum_i(p_i+A)>.  This is the raw particle-current quantity, not charge
! current density and not an unnormalized matrix element.
  norm_squared = real(inner_product(psi,psi))
  if (norm_squared <= 0d0) stop 'Cannot evaluate current for zero norm.'
  total_particle_velocity = particle_velocity_numerator(psi,avec)/norm_squared

end function total_particle_velocity
!-------------------------------------------------------
real(8) function charge_current_density(psi, avec)
  implicit none
  complex(8),intent(in) :: psi(0:nx-1,0:nx-1,0:nx-1)
  real(8),intent(in) :: avec

! Electron charge q=-1 in atomic units; divide the total ring current by L.
  charge_current_density = -total_particle_velocity(psi,avec) &
      /bvc_lattice_constant

end function charge_current_density
!-------------------------------------------------------
real(8) function normalized_energy(psi, avec)
  implicit none
  complex(8),intent(in) :: psi(0:nx-1,0:nx-1,0:nx-1)
  real(8),intent(in) :: avec
  complex(8),allocatable :: hpsi(:,:,:)
  real(8) :: norm_squared

  norm_squared = real(inner_product(psi,psi))
  if (norm_squared <= 0d0) stop 'Cannot evaluate energy for zero norm.'
  allocate(hpsi(0:nx-1,0:nx-1,0:nx-1))
  call apply_hamiltonian(psi,hpsi,avec)
  normalized_energy = real(inner_product(psi,hpsi))/norm_squared
  deallocate(hpsi)

end function normalized_energy
!-------------------------------------------------------
subroutine check_current_operator(psi)
  implicit none
  complex(8),intent(in) :: psi(0:nx-1,0:nx-1,0:nx-1)
  complex(8),allocatable :: scaled_psi(:,:,:),trial_psi(:,:,:)
  real(8),parameter :: avec_test = 1d-2
  real(8),parameter :: delta_avec = 1d-5
  real(8),parameter :: test_scale = 1.234d0
  real(8) :: derivative_fd, particle_velocity, abs_error, tolerance
  real(8) :: energy_scale_error, current_scale_error
  real(8) :: trial_derivative_fd,trial_momentum,trial_error,boost_k
  integer :: ix1,ix2,ix3

! For a fixed state, central differentiation is exact for the quadratic A
! dependence up to floating-point cancellation.  This tests that the current
! operator uses precisely the same fourth-order stencil as H(A).
  derivative_fd = (normalized_energy(psi,avec_test+delta_avec) &
      - normalized_energy(psi,avec_test-delta_avec))/(2d0*delta_avec)
  particle_velocity = total_particle_velocity(psi,avec_test)
  abs_error = abs(derivative_fd-particle_velocity)
  tolerance = 1d-9*max(1d0,abs(derivative_fd),abs(particle_velocity))
  write(*,'(a)') 'Ground-state current-operator test at finite A:'
  write(*,'(a,1pe16.8)') '  d<H>/dA                   = ',derivative_fd
  write(*,'(a,1pe16.8)') '  <sum_i(p_i+A)>            = ',particle_velocity
  write(*,'(a,1pe16.8)') '  absolute error            = ',abs_error
  if (abs_error > tolerance) then
    write(*,'(a,1pe16.8)') 'Warning: current-operator test tolerance = ',tolerance
  end if

! Deliberately rescale the fixed state to confirm that reported observables,
! unlike the diagnostic norm and raw projections, do not inherit norm drift.
  allocate(scaled_psi(0:nx-1,0:nx-1,0:nx-1))
  scaled_psi = test_scale*psi
  energy_scale_error = abs(normalized_energy(scaled_psi,avec_test) &
      - normalized_energy(psi,avec_test))
  current_scale_error = abs(total_particle_velocity(scaled_psi,avec_test) &
      - total_particle_velocity(psi,avec_test))

! A ring-compatible center-of-mass phase preserves both PBC and fermionic
! antisymmetry while producing a reproducible nonzero canonical momentum.
  allocate(trial_psi(0:nx-1,0:nx-1,0:nx-1))
  boost_k=2d0*pi/bvc_lattice_constant
  do ix1=0,nx-1
    do ix2=0,nx-1
      do ix3=0,nx-1
        trial_psi(ix1,ix2,ix3)=psi(ix1,ix2,ix3)* &
            exp(zi*boost_k*(xn(ix1)+xn(ix2)+xn(ix3)))
      end do
    end do
  end do
  call antisymmetrize(trial_psi)
  call normalize(trial_psi)
  trial_derivative_fd=(normalized_energy(trial_psi,delta_avec)- &
      normalized_energy(trial_psi,-delta_avec))/(2d0*delta_avec)
  trial_momentum=total_particle_velocity(trial_psi,0d0)
  trial_error=abs(trial_derivative_fd-trial_momentum)
  tolerance=1d-9*max(1d0,abs(trial_derivative_fd),abs(trial_momentum))
  write(*,'(a)') 'Nonzero-momentum antisymmetric trial-state test at A=0:'
  write(*,'(a,1pe16.8)') '  d<H>/dA                   = ',trial_derivative_fd
  write(*,'(a,1pe16.8)') '  <sum_i p_i>               = ',trial_momentum
  write(*,'(a,1pe16.8)') '  absolute error            = ',trial_error
  write(*,'(a,1pe16.8)') '  antisymmetry error        = ', &
      antisymmetry_error(trial_psi)
  if (abs(trial_momentum) < 1d-6) then
    write(*,'(a)') 'Warning: trial-state momentum unexpectedly small.'
  end if
  if (trial_error > tolerance) then
    write(*,'(a,1pe16.8)') 'Warning: trial current-test tolerance = ',tolerance
  end if
  write(*,'(a,1pe16.8)') 'Norm-scaling test energy error    = ',energy_scale_error
  write(*,'(a,1pe16.8)') 'Norm-scaling test current error   = ',current_scale_error
  deallocate(scaled_psi,trial_psi)

end subroutine check_current_operator
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
  if (allocated(eigenvalues)) deallocate(eigenvalues)
  if (allocated(eigenstate_residuals)) deallocate(eigenstate_residuals)
  if (allocated(eigenstate_iterations)) deallocate(eigenstate_iterations)
  if (allocated(eigenstate_converged)) deallocate(eigenstate_converged)
  if (allocated(eigenvectors)) deallocate(eigenvectors)
  if (allocated(vpot_1d)) deallocate(vpot_1d)
  if (allocated(wpot_1d)) deallocate(wpot_1d)
  if (allocated(onebody_values)) deallocate(onebody_values)
  if (allocated(onebody_vectors)) deallocate(onebody_vectors)
  if (allocated(vpot)) deallocate(vpot)
  if (allocated(wpot)) deallocate(wpot)
  if (allocated(tot_pot)) deallocate(tot_pot)

end subroutine finalize
!-------------------------------------------------------
end program main
