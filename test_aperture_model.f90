! test_aperture_model.f90
!
! Test of the Aperture fluid model (fluid_model_aperture.f90) on a binary
! written by aperture_to_binary.py.  At exact cell centers the bilinear
! weights are 0/1, so aperture_vals must return the stored values.  Checks
! (each a FAIL if violated; tolerances are for the single-precision outputs):
!   1. n_e and theta_e at the cell equal the stored n_e, theta_e (relative
!      1e-6; cells with n_e > 0 that come back empty are counted, not failed:
!      the per-sample sigma cut of an APL8 binary, or a non-timelike u);
!   2. u and b equal the stored Kerr-Schild u^mu, b^mu of that cell, taken to
!      Boyer-Lindquist with uks2ubl and then renormalized and projected as
!      aperture_vals does (relative 1e-5): this tests the indexing and the
!      weights, not the transform;
!   3. the transform itself: the raw (unnormalized) stored vectors have the
!      same u.u, b.b and b.u in the Kerr-Schild metric (kerr.f90,
!      ksmetric_cov) as their uks2ubl images in the Boyer-Lindquist metric
!      (relative 1e-8 of the vectors' norms), which only holds if uks2ubl is
!      the coordinate change of that metric;
!   4. the returned u.u = -1 and b.u = 0 (by construction, as before).
! It also prints the stored vectors' own |u.u + 1| (nonzero for coarse-
! grained files, whose block means are not unit vectors) and the ranges.
!
! Usage: ./test_aperture_model [fluid_binary [camera_spin]]
!   default binary: Cylco-Synchrotron/aperture_fluid.bin (the shared scratch
!   file of the last trace).  Without camera_spin the binary's own spin is
!   used; with it, initialize_aperture_model is called with that spin, so a
!   mismatch stops the program (exit 1) as it stops grtrans.
!
! Build (from grtrans/, after make grtrans):
!   gfortran -O2 -std=legacy -fopenmp test_aperture_model.f90 \
!       fluid_model_aperture.o jnu_table.o cyclotron_emis.o kerr.o \
!       class_four_vector.o interpolate.o interpolate_aux.o math.o \
!       phys_constants.o -o test_aperture_model
!   (fluid_model_aperture uses jnu_table's set_hybrid_dispatch, which needs
!   cyclotron_emis and OpenMP; interpolate needs interpolate_aux's dhunt)

program test_aperture_model
  use fluid_model_aperture
  use class_four_vector
  use kerr, only: kerr_metric, uks2ubl
  implicit none

  integer, parameter :: NP = 20000
  type(four_vector), allocatable :: x0(:), u(:), b(:), uks(:), bks(:), ubl(:), bbl(:)
  real, allocatable :: rho(:), p(:), bmag(:)
  real(kind=8), allocatable :: rr(:), tt(:), pp(:)
  integer, allocatable :: kk(:)
  real :: a
  integer :: i, ir, ith, nvalid, nempty, nfail
  double precision :: udotu, bdotu, udotmax, bnorm, brelmax, rawmax
  double precision :: nerr, terr, uerr, berr, inv_err, sc, uu, ub, d
  character(len=512) :: path, arg

  path = '/home/aaron/Research/Cylco-Synchrotron/aperture_fluid.bin'
  if (command_argument_count() >= 1) call get_command_argument(1, path)
  if (command_argument_count() >= 2) then
     call get_command_argument(2, arg)
     read(arg, *) a
     call initialize_aperture_model(dble(a), trim(path))
  else
     call read_aperture_data(trim(path))
     a = real(abins(1)%asim)
  endif
  ! the first (here the only) binary read (fluid_model_aperture, abins)
  associate (nr => abins(1)%nr, nth => abins(1)%nth, r_arr => abins(1)%r_arr, &
       th_arr => abins(1)%th_arr, n_e => abins(1)%n_e, theta_e => abins(1)%theta_e, &
       u0 => abins(1)%u0, u1 => abins(1)%u1, u2 => abins(1)%u2, u3 => abins(1)%u3, &
       b0 => abins(1)%b0, b1 => abins(1)%b1, b2 => abins(1)%b2, b3 => abins(1)%b3)

  allocate(x0(NP), u(NP), b(NP), rho(NP), p(NP), bmag(NP), kk(NP))
  allocate(uks(NP), bks(NP), ubl(NP), bbl(NP), rr(NP), tt(NP), pp(NP))

  ! exact cell centers, spread over the grid
  do i = 1, NP
     ir = 1 + mod(i*37, nr); ith = 1 + mod(i*53, nth)
     kk(i) = (ith-1)*nr + ir
     rr(i) = r_arr(ir); tt(i) = th_arr(ith); pp(i) = 0d0
     x0(i)%data = (/0d0, rr(i), tt(i), 0d0/)
     uks(i)%data = (/u0(kk(i)), u1(kk(i)), u2(kk(i)), u3(kk(i))/)
     bks(i)%data = (/b0(kk(i)), b1(kk(i)), b2(kk(i)), b3(kk(i))/)
  enddo

  call aperture_vals(x0, a, rho, p, b, u, bmag)

  ! the expected values: the stored vectors through uks2ubl, in both metrics
  ubl = uks2ubl(uks, rr, dble(a))
  bbl = uks2ubl(bks, rr, dble(a))
  call assign_metric(ubl, transpose(kerr_metric(rr, tt, dble(a))))
  call assign_metric(bbl, transpose(kerr_metric(rr, tt, dble(a))))
  call assign_metric(uks, transpose(kerr_metric(rr, tt, pp, dble(a))))
  call assign_metric(bks, transpose(kerr_metric(rr, tt, pp, dble(a))))

  udotmax = 0d0; brelmax = 0d0; rawmax = 0d0
  nerr = 0d0; terr = 0d0; uerr = 0d0; berr = 0d0; inv_err = 0d0
  nvalid = 0; nempty = 0
  do i = 1, NP
     ! 3. invariance of the scalar products under uks2ubl (every cell with a
     !    timelike stored u, whatever its n_e)
     if (uks(i)*uks(i) < 0d0) then
        sc = abs(uks(i)%data(1))**2 + sum(abs(uks(i)%data(2:4)))**2 &
             + sum(abs(bks(i)%data))**2 + 1d-300
        inv_err = max(inv_err, abs(uks(i)*uks(i) - ubl(i)*ubl(i))/sc, &
             abs(bks(i)*bks(i) - bbl(i)*bbl(i))/sc, abs(bks(i)*uks(i) - bbl(i)*ubl(i))/sc)
     endif
     if (.not. (n_e(kk(i)) > 0d0)) cycle
     if (rho(i) <= 0.0) then
        nempty = nempty + 1
        cycle
     endif
     nvalid = nvalid + 1
     ! 1. scalars
     nerr = max(nerr, abs(rho(i) - n_e(kk(i)))/n_e(kk(i)))
     if (theta_e(kk(i)) > 0d0) terr = max(terr, abs(p(i) - theta_e(kk(i)))/theta_e(kk(i)))
     ! 2. vectors: normalize/project the expected ones as aperture_vals does
     uu = ubl(i)*ubl(i)
     rawmax = max(rawmax, abs(uks(i)*uks(i) + 1d0))
     ubl(i)%data = ubl(i)%data/sqrt(-uu)
     ub = bbl(i)*ubl(i)
     bbl(i)%data = bbl(i)%data + ub*ubl(i)%data
     d = maxval(abs(u(i)%data - ubl(i)%data))/maxval(abs(ubl(i)%data))
     uerr = max(uerr, d)
     bnorm = maxval(abs(bbl(i)%data))
     if (bnorm > 0d0) berr = max(berr, maxval(abs(b(i)%data - bbl(i)%data))/bnorm)
     ! 4. invariants of the returned vectors
     udotu = u(i)*u(i)
     bdotu = b(i)*u(i)
     udotmax = max(udotmax, abs(udotu + 1d0))
     bnorm = sqrt(max(b(i)*b(i), 0d0))
     if (bnorm > 0d0) brelmax = max(brelmax, abs(bdotu)/bnorm)
  enddo

  print *, '--- Aperture fluid model test (exact cell centers) ---'
  print '(A,A)',       ' binary              : ', trim(path)
  print '(A,F10.6)',   ' spin                : ', a
  print '(A,I0)',      ' sample points       : ', NP
  print '(A,I0)',      ' cells with n_e > 0  : ', nvalid + nempty
  print '(A,I0)',      '   returned empty    : ', nempty
  print '(A,ES12.4)',  ' 1. max rel err n_e    : ', nerr
  print '(A,ES12.4)',  '    max rel err theta_e: ', terr
  print '(A,ES12.4)',  ' 2. max rel err u      : ', uerr
  print '(A,ES12.4)',  '    max rel err b      : ', berr
  print '(A,ES12.4)',  ' 3. KS vs BL invariants: ', inv_err
  print '(A,ES12.4)',  ' 4. max |u.u + 1|      : ', udotmax
  print '(A,ES12.4)',  '    max |b.u| / |b|    : ', brelmax
  print '(A,ES12.4)',  ' stored |u.u + 1| (KS) : ', rawmax
  print '(A,2ES12.4)', ' n_e     range         : ', minval(rho),  maxval(rho)
  print '(A,2ES12.4)', ' theta_e range         : ', minval(p),    maxval(p)
  print '(A,2ES12.4)', ' |b|     range         : ', minval(bmag), maxval(bmag)

  nfail = 0
  if (nvalid == 0) then
     print *, 'FAIL: no cell with n_e > 0 was checked'; nfail = nfail + 1
  endif
  if (nerr > 1d-6 .or. terr > 1d-6) then
     print *, 'FAIL: n_e or theta_e at a cell center differs from the stored value'
     nfail = nfail + 1
  endif
  if (uerr > 1d-5 .or. berr > 1d-5) then
     print *, 'FAIL: u or b at a cell center differs from the stored vector'
     nfail = nfail + 1
  endif
  if (inv_err > 1d-8) then
     print *, 'FAIL: uks2ubl does not preserve the scalar products (KS -> BL)'
     nfail = nfail + 1
  endif
  if (udotmax > 1d-5 .or. brelmax > 1d-5) then
     print *, 'FAIL: returned u.u /= -1 or b.u /= 0'
     nfail = nfail + 1
  endif
  if (nfail == 0) then
     print *, 'PASS: cell-center values, KS->BL transform and invariants'
  else
     stop 1
  endif

  end associate
  call del_aperture_data()
end program test_aperture_model
