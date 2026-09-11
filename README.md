# Three spinless fermions on a one-dimensional ring

Build and run:

```text
gfortran -std=f95 -O2 tdse.f90 -o tdse
./tdse < input_tdse
```

The input uses three required records and up to two optional records:

```text
lattice_constant(bohr)  nx
Tprop(fs)                dt(atomic time)
E0(MV/m) omega(eV) Tpulse(fs) CEP/(2*pi)
w0(Hartree)
run_triplet run_scaling require_converged_ground_state num_eigenstates cg_max_iter
```

`w0` defaults to `0.0`, the intentional noninteracting baseline. The last
record defaults to `.false. .false. .true. 5 2000`. At least five states are
retained for the ground state and four-state first-excited manifold. The strict
residual target remains `1e-9`; by default, propagation is disabled if the
ground state does not meet it.

A single-run `current.out` contains all response orders and is not itself called
a shift current. Triplet mode reuses one ground state for `+E0`, `-E0`, and zero
field and writes `current_second_order.out` with

```text
J_even_induced = (J(+E0) + J(-E0) - 2 J(0))/2.
```

Optional scaling mode repeats the diagnostic at `E0/2` and writes
`second_order_scaling.out`.
