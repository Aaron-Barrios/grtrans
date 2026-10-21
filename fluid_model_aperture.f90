! fluid_model_aperture.f90
!
! Stage 1 (data layer) for ray tracing Aperture4 GR-PIC output with grtrans.
! Reads the stream binary produced by aperture_to_binary.py (no HDF5 in Fortran)
! and bilinearly interpolates the reconstructed emitter variables onto ray
! points.  Mirrors fluid_model_harm.f90 but the grid is Kerr-Schild (r, theta),
! which shares the SAME r and theta as Boyer-Lindquist, so x0%data(2)/data(3)
! index the grid directly.
!
! Conventions (all verified against reconstruct_fluid.cpp and the Aperture4
! source):
!   * grid: axis0 = theta (polar angle, [0, pi]), axis1 = r (log-spaced radius).
!     flat index k = ith*nr + ir  (theta slowest, radius fastest).
!   * u^mu and b^mu are CONTRAVARIANT, in Kerr-Schild spherical coordinates.
!     They are transformed to Boyer-Lindquist (uks2ubl) before use, because
!     grtrans geodesics/metrics live in BL.
!   * rho  slot carries n_e  (electron rest-frame density n_e_rest, code units).
!   * p    slot carries theta_e (electron rest-frame temperature P_e/n_e_rest).
!   * bmag = |b| = sqrt(b^mu b_mu), code units: the binary's bmag array is
!     skipped when reading; aperture_vals takes |b| from the interpolated b^mu
!     (see below).
!   * Interpolated u^mu is renormalized to u.u = -1 and b^mu projected
!     orthogonal to it.  A weighted average of unit timelike vectors pointing
!     in different directions is not a unit vector (|u.u| > 1), so bilinear
!     interpolation alone gives a wrong u^t and hence a wrong Doppler factor:
!     averaged over 2x2 cells, 11-13% of the points in the late dumps had
!     |u.u + 1| > 0.1 (docs/history.tex, "Interpolated four-velocities were
!     not normalized").  Since 2026-10-10 |b| is taken from the
!     interpolated, projected b^mu (Aperture interpolates field components),
!     and theta_e as <n_e theta_e>/<n_e> (see aperture_vals); the binary's
!     bmag array is skipped.
!   * Binary magic 'APL7' or 'APL8' (aperture_to_binary.py's docstring has the
!     layout).  APL7: any sigma cut was applied per cell by the writer
!     (--sigma-cut-mode cell: cells above it have n_e = theta_e = 0).  APL8
!     (--sigma-cut-mode point): the cells are written uncut and one float64
!     after the sigma array holds the cut, which aperture_vals applies at each
!     ray sample (point_sigcut).  The pipeline's store tags mark the two as
!     _sc<cut> (cell) and _sp<cut> (point).
!   * The camera spin (grtrans's &geodata spin) must equal the spin the fluid
!     was reconstructed with (the header's asim, i.e. this binary's own spin):
!     a mismatch stops grtrans (initialize_aperture_model; it only warned
!     before 2026-10-15), and raytrace_q6.py refuses it before grtrans runs
!     (check_camera_spin).
!   * Several binaries in one run (2026-10-15): &harm fdfile naming a '.lst'
!     file (one binary path per line) loads every binary listed into abins;
!     with nmdot equal to their number, grtrans_driver traces binary m along
!     the same ray paths into the cameras of model parameter m (one geodesic
!     per pixel for all of them; IMPLEMENTATION.md, "Several binaries in one
!     run").  The binaries may have different grids, lambda_u, mass ratios
!     and sigma cuts, but must share the spin (the camera's) and the
!     emissivity dispatch.  A single binary behaves as before.

      module fluid_model_aperture

      use class_four_Vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, uks2ubl
      use jnu_table, only: set_hybrid_dispatch
      implicit none

      ! One fluid binary: its grid, header values and the fields aperture_vals
      ! reads.  The writer's |b| and validity arrays are skipped when reading
      ! (|b| is recomputed from the interpolated b^mu; invalid cells are zero
      ! in every field, see aperture_vals), which saves 2 of the 16 arrays.
      type aperture_bin
         character(len=1024) :: path = ''
         integer :: nr = 0, nth = 0
         real(kind=8) :: asim = 0d0
         ! cm^-3 per code density unit, from the binary's header (physical_scale.py
         ! computes it from the run's first dump; fluid.f90 converts with it)
         real(kind=8) :: lambda_u = 0d0
         ! code particle mass in units of m_e, from the header: 1 for the pair-
         ! plasma reading, N = m_p/m_e when relabeled (relabel.py); fluid.f90's
         ! |b| -> B carries sqrt(mass_ratio)
         real(kind=8) :: mass_ratio = 0d0
         ! the emissivity dispatch of the header (jnu_table.set_hybrid_dispatch)
         real(kind=8) :: th_hyb = 0d0, nth_max = 0d0, th_syn = 0d0
         ! sigma cut applied at each ray sample (APL8 binaries, written with
         ! pipeline_config.TRACE_SIGMA_CUT_MODE = "point"): a sample whose sigma,
         ! interpolated over the corners that hold fluid (sigma_arr > 0), exceeds
         ! it is left empty (no n_e, Theta_e, u, b: it neither emits nor absorbs),
         ! as ipole's geodesic sigma cut.  0: none here (an APL7
         ! binary was cut per cell by the writer, or not at all)
         real(kind=8) :: point_sigcut = 0d0
         real(kind=8), dimension(:), allocatable :: r_arr, th_arr
         real(kind=8), dimension(:), allocatable :: n_e, theta_e
         real(kind=8), dimension(:), allocatable :: u0, u1, u2, u3
         real(kind=8), dimension(:), allocatable :: b0, b1, b2, b3
         ! the pipeline's plasma beta_p per cell (relabel.beta_p; 0 where no
         ! field), interpolated and emission-weighted for extra slot 20
         real(kind=8), dimension(:), allocatable :: beta_arr
         ! the pipeline's magnetization sigma = b^2/(2 rho_0) per cell (pair fluid;
         ! 0 where no fluid), interpolated and emission-weighted for extra slot 21
         real(kind=8), dimension(:), allocatable :: sigma_arr
      end type aperture_bin

      ! The loaded binaries.  One, as before, when &harm fdfile names a binary;
      ! several when fdfile ends in '.lst' (a text file, one binary path per
      ! line; blank lines and lines starting with '#' are skipped): grtrans
      ! then traces each binary along the same ray paths, binary m in the
      ! cameras of model parameter m (nmdot must equal the number of binaries;
      ! grtrans_driver selects it with select_aperture_binary).
      type (aperture_bin), dimension(:), allocatable, target :: abins
      integer :: naperture = 0
      ! the binary aperture_vals and convert_fluidvars_aperture use, per
      ! OpenMP thread (each thread traces its own rays); 1 unless several
      ! binaries are loaded
      integer :: aperture_active = 1
!$omp threadprivate(aperture_active)

      interface aperture_vals
        module procedure aperture_vals
      end interface

      interface initialize_aperture_model
        module procedure initialize_aperture_model
      end interface

      interface del_aperture_data
        module procedure del_aperture_data
      end interface

      contains

        subroutine init_aperture_data(bn, nr_in, nth_in)
        type (aperture_bin), intent(inout) :: bn
        integer, intent(in) :: nr_in, nth_in
        bn%nr = nr_in; bn%nth = nth_in
        allocate(bn%r_arr(nr_in), bn%th_arr(nth_in))
        allocate(bn%n_e(nr_in*nth_in), bn%theta_e(nr_in*nth_in))
        allocate(bn%u0(nr_in*nth_in), bn%u1(nr_in*nth_in), bn%u2(nr_in*nth_in), bn%u3(nr_in*nth_in))
        allocate(bn%b0(nr_in*nth_in), bn%b1(nr_in*nth_in), bn%b2(nr_in*nth_in), bn%b3(nr_in*nth_in))
        allocate(bn%beta_arr(nr_in*nth_in), bn%sigma_arr(nr_in*nth_in))
        end subroutine init_aperture_data

        ! Read one binary into bn (stream access: consecutive reads continue
        ! where the previous one stopped).
        subroutine read_one_binary(dfile, bn)
        character(len=*), intent(in) :: dfile
        type (aperture_bin), intent(inout) :: bn
        integer :: magic, nrr, nthh, ios
        real(kind=8), dimension(:), allocatable :: skip
        open(unit=8, file=dfile, form='unformatted', access='stream', &
             status='old', action='read', iostat=ios)
        if (ios /= 0) then
           write(6,*) 'ERROR: cannot open aperture binary ', trim(dfile)
           stop 1
        endif
        read(8) magic
        if (magic == int(z'41504C46') .or. magic == int(z'41504C32') &
             .or. magic == int(z'41504C33') .or. magic == int(z'41504C34') &
             .or. magic == int(z'41504C35') .or. magic == int(z'41504C36')) then
           write(6,*) 'ERROR: ', trim(dfile), ' has an old layout (APLF..APL6,'
           write(6,*) '  without lambda_u, the mass ratio, beta_p, sigma or the'
           write(6,*) '  emissivity dispatch); rewrite it with the current'
           write(6,*) '  aperture_to_binary.py'
           close(8)
           stop 1
        endif
        if (magic /= int(z'41504C37') .and. magic /= int(z'41504C38')) then
           write(6,*) 'ERROR: bad magic in aperture binary ', trim(dfile)
           close(8)
           stop 1
        endif
        bn%path = dfile
        read(8) nrr, nthh, bn%asim, bn%lambda_u, bn%mass_ratio, &
             bn%th_hyb, bn%nth_max, bn%th_syn
        if (.not. (bn%lambda_u > 0d0) .or. .not. (bn%mass_ratio > 0d0)) then
           write(6,*) 'ERROR: lambda_u or mass_ratio <= 0 in aperture binary ', trim(dfile)
           close(8)
           stop 1
        endif
        call init_aperture_data(bn, nrr, nthh)
        allocate(skip(nrr*nthh))
        ! layout: r, th, n_e, theta_e, |b| (skipped), u0..u3, b0..b3,
        ! valid (skipped), beta_p, sigma
        read(8) bn%r_arr, bn%th_arr, bn%n_e, bn%theta_e
        read(8) skip
        read(8) bn%u0, bn%u1, bn%u2, bn%u3, bn%b0, bn%b1, bn%b2, bn%b3
        read(8) skip
        read(8) bn%beta_arr, bn%sigma_arr
        deallocate(skip)
        bn%point_sigcut = 0d0
        if (magic == int(z'41504C38')) read(8) bn%point_sigcut
        close(8)
        write(6,*) 'read aperture data: nr=', bn%nr, ' nth=', bn%nth
        write(6,*) '  r in [', bn%r_arr(1), ',', bn%r_arr(bn%nr), ']'
        write(6,*) '  th in [', bn%th_arr(1), ',', bn%th_arr(bn%nth), ']'
        write(6,*) '  lambda_u = ', bn%lambda_u, ' cm^-3 per code unit'
        write(6,*) '  code particle mass = ', bn%mass_ratio, ' m_e'
        if (bn%point_sigcut > 0d0) write(6,*) '  sigma cut at each ray sample: ', bn%point_sigcut
        end subroutine read_one_binary

        ! Read the binary dfile, or every binary listed in dfile when its name
        ! ends in '.lst' (see abins).  All binaries must carry the same
        ! emissivity dispatch (it is one setting of jnu_table, shared by the
        ! threads); grids, lambda_u, mass ratio and sigma cut may differ.
        subroutine read_aperture_data(dfile)
        character(len=*), intent(in) :: dfile
        character(len=1024), dimension(:), allocatable :: paths
        integer :: k, n
        call del_aperture_data()
        n = len_trim(dfile)
        if (n >= 4) then
           if (dfile(n-3:n) == '.lst') then
              call read_binary_list(dfile, paths)
           endif
        endif
        if (.not. allocated(paths)) then
           allocate(paths(1)); paths(1) = dfile
        endif
        naperture = size(paths)
        allocate(abins(naperture))
        do k = 1, naperture
           if (naperture > 1) write(6,'(a,i0,a,i0,2a)') ' aperture binary ', k, ' of ', &
                naperture, ': ', trim(paths(k))
           call read_one_binary(trim(paths(k)), abins(k))
           if (k > 1) then
              if (abins(k)%th_hyb /= abins(1)%th_hyb .or. abins(k)%nth_max /= abins(1)%nth_max &
                   .or. abins(k)%th_syn /= abins(1)%th_syn) then
                 write(6,*) 'ERROR: aperture binary ', trim(paths(k)), ' has another', &
                      ' emissivity dispatch (THETA_HYBRID, NTHETA_MAX, THETA_SYNCH)', &
                      ' than ', trim(paths(1)), '; trace them in separate runs'
                 stop 1
              endif
           endif
        enddo
        deallocate(paths)
        aperture_active = 1
        ! the emissivity dispatch is the pipeline's (pipeline_config.py)
        call set_hybrid_dispatch(abins(1)%th_hyb, abins(1)%nth_max, abins(1)%th_syn)
        end subroutine read_aperture_data

        ! The binary paths in the list file lfile, one per line (leading and
        ! trailing blanks removed; blank lines and '#' comments skipped).
        subroutine read_binary_list(lfile, paths)
        character(len=*), intent(in) :: lfile
        character(len=1024), dimension(:), allocatable, intent(out) :: paths
        character(len=1024) :: line
        integer :: u, ios, n, pass
        open(newunit=u, file=lfile, form='formatted', status='old', action='read', iostat=ios)
        if (ios /= 0) then
           write(6,*) 'ERROR: cannot open aperture binary list ', trim(lfile)
           stop 1
        endif
        do pass = 1, 2
           n = 0
           rewind(u)
           do
              read(u, '(a)', iostat=ios) line
              if (ios /= 0) exit
              line = adjustl(line)
              if (len_trim(line) == 0) cycle
              if (line(1:1) == '#') cycle
              if (len_trim(line) >= len(line)) then
                 write(6,*) 'ERROR: path too long in aperture binary list ', trim(lfile)
                 stop 1
              endif
              n = n + 1
              if (pass == 2) paths(n) = line
           enddo
           if (pass == 1) then
              if (n == 0) then
                 write(6,*) 'ERROR: no binary in aperture binary list ', trim(lfile)
                 stop 1
              endif
              allocate(paths(n))
           endif
        enddo
        close(u)
        end subroutine read_binary_list

        ! Make binary k the one aperture_vals and convert_fluidvars_aperture use
        ! on the calling thread.
        subroutine select_aperture_binary(k)
        integer, intent(in) :: k
        if (k < 1 .or. k > naperture) then
           write(6,*) 'ERROR: select_aperture_binary: no binary ', k, ' (', naperture, ' loaded)'
           stop 1
        endif
        aperture_active = k
        end subroutine select_aperture_binary

        ! cm^-3 per code density unit and code particle mass [m_e] of the
        ! calling thread's binary (convert_fluidvars_aperture, fluid.f90)
        real(kind=8) function aperture_lambda_u()
        aperture_lambda_u = abins(aperture_active)%lambda_u
        end function aperture_lambda_u

        real(kind=8) function aperture_mass_ratio()
        aperture_mass_ratio = abins(aperture_active)%mass_ratio
        end function aperture_mass_ratio

        ! Bisection bracket: return idx in [1, n-1] with xx(idx) <= x < xx(idx+1)
        subroutine bracket(xx, n, x, idx)
        integer, intent(in) :: n
        real(kind=8), intent(in) :: xx(n), x
        integer, intent(out) :: idx
        integer :: lo, hi, mid
        lo = 1; hi = n
        do while (hi - lo > 1)
           mid = (lo + hi) / 2
           if (xx(mid) <= x) then
              lo = mid
           else
              hi = mid
           endif
        enddo
        idx = max(1, min(lo, n-1))
        end subroutine bracket

        ! Interpolation as Aperture4 does it (utils/interpolation.hpp, interp_t<1>
        ! in coord_policy_gr_ks_sph.hpp): first order (cloud-in-cell) in the
        ! code's own coordinates x1 = ln r and x2 = theta, with each value at
        ! the center of the cell it describes (aperture_to_binary.py writes the
        ! block centers, not grid.h5's first-native-cell labels), and applied
        ! to the raw components, never to derived quantities:
        !   * n_e and the electron pressure n_e theta_e are interpolated, and
        !     theta_e = <n_e theta_e>/<n_e> (a cell without electrons adds
        !     nothing to either, instead of pulling theta_e toward 0);
        !   * u^mu and b^mu component by component (then u renormalized and b
        !     projected orthogonal to it), and |b| from the interpolated b^mu,
        !     not an interpolated |b|;
        !   * no validity gate and no renormalization of the weights, as in
        !     Aperture's block means: invalid cells (no timelike current) are
        !     zero in every quantity; electron-free cells and sigma-masked
        !     cells (cell mode) count as zero only in n_e and n_e theta_e and
        !     keep a velocity and field (the pair fluid's, or the electrons'),
        !     which enter neighboring points.
        ! beta_p and sigma (diagnostics) are interpolated as given.
        subroutine aperture_vals(x0, a, rho, p, b, u, bmag, beta, sigma)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real, dimension(size(x0)), intent(out) :: rho, p, bmag
        ! the pipeline's beta_p and sigma, interpolated (optional: grtrans's
        ! driver asks for them, the stand-alone test does not)
        real, dimension(size(x0)), intent(out), optional :: beta, sigma
        type (four_Vector), intent(out), dimension(size(x0)) :: u, b
        integer :: npts, i, lr, ur, lth, uth, k00, k01, k10, k11
        real(kind=8) :: rq, thq, rd, td, w00, w01, w10, w11
        real(kind=8) :: ne_i, pe_i, sw, ss, bs
        real(kind=8) :: rmin, rmax, thmin, thmax
        type (four_Vector), dimension(size(x0)) :: u_ks, b_ks
        real(kind=8), dimension(size(x0)) :: rr, tth
        logical, dimension(size(x0)) :: filled
        real(kind=8) :: uu, ub, bb

        ! the calling thread's binary (select_aperture_binary)
        associate (bn => abins(aperture_active))
        associate (nr => bn%nr, nth => bn%nth, r_arr => bn%r_arr, th_arr => bn%th_arr, &
             n_e => bn%n_e, theta_e => bn%theta_e, u0 => bn%u0, u1 => bn%u1, &
             u2 => bn%u2, u3 => bn%u3, b0 => bn%b0, b1 => bn%b1, b2 => bn%b2, &
             b3 => bn%b3, beta_arr => bn%beta_arr, sigma_arr => bn%sigma_arr, &
             point_sigcut => bn%point_sigcut)
        npts = size(x0)
        rmin = r_arr(1); rmax = r_arr(nr)
        thmin = th_arr(1); thmax = th_arr(nth)
        rr = x0%data(2)    ! BL r == KS r
        tth = x0%data(3)   ! BL theta == KS theta

        ! placeholders for every point (outside the cell centers, or no electrons)
        filled = .false.
        do i = 1, npts
           rho(i) = 0.; p(i) = 0.; bmag(i) = 0.
           if (present(beta)) beta(i) = 0.
           if (present(sigma)) sigma(i) = 0.
           u_ks(i)%data = (/1d0, 0d0, 0d0, 0d0/)
           b_ks(i)%data = 0d0
        enddo

        do i = 1, npts
           rq = rr(i); thq = tth(i)
           if (rq < rmin .or. rq > rmax .or. thq < thmin .or. thq > thmax) then
              cycle
           endif
           call bracket(r_arr, nr, rq, lr);  ur = lr + 1
           call bracket(th_arr, nth, thq, lth); uth = lth + 1
           ! linear in x1 = ln r and x2 = theta (Aperture's cell coordinates)
           rd = log(rq / r_arr(lr)) / log(r_arr(ur) / r_arr(lr))
           td = (thq - th_arr(lth)) / (th_arr(uth) - th_arr(lth))
           ! flat corners (theta slowest): k = ith*nr + ir
           k00 = (lth-1)*nr + lr
           k01 = (uth-1)*nr + lr
           k10 = (lth-1)*nr + ur
           k11 = (uth-1)*nr + ur
           w00 = (1d0-rd)*(1d0-td)
           w01 = (1d0-rd)*td
           w10 = rd*(1d0-td)
           w11 = rd*td

           ne_i = w00*n_e(k00) + w01*n_e(k01) + w10*n_e(k10) + w11*n_e(k11)
           ! beta_p and sigma (extra slots 20, 21) from the corners that hold
           ! fluid (sigma_arr > 0), renormalized, as the per-sample sigma cut
           ! below tests it: an empty cell is sparse sampling, not a zero
           ! (muchado.md L1-R18; until 2026-10-19 empty corners counted as 0)
           sw = 0d0; ss = 0d0; bs = 0d0
           if (sigma_arr(k00) > 0d0) then
              sw = sw + w00; ss = ss + w00*sigma_arr(k00); bs = bs + w00*beta_arr(k00)
           endif
           if (sigma_arr(k01) > 0d0) then
              sw = sw + w01; ss = ss + w01*sigma_arr(k01); bs = bs + w01*beta_arr(k01)
           endif
           if (sigma_arr(k10) > 0d0) then
              sw = sw + w10; ss = ss + w10*sigma_arr(k10); bs = bs + w10*beta_arr(k10)
           endif
           if (sigma_arr(k11) > 0d0) then
              sw = sw + w11; ss = ss + w11*sigma_arr(k11); bs = bs + w11*beta_arr(k11)
           endif
           if (sw > 0d0) then
              if (present(beta)) beta(i) = real(bs / sw)
              if (present(sigma)) sigma(i) = real(ss / sw)
           endif
           if (.not. (ne_i > 0d0)) cycle
           if (point_sigcut > 0d0) then
              ! sigma at the sample from the corners that hold fluid (an empty
              ! cell, sigma_arr = 0, is sparse sampling, not vacuum)
              sw = 0d0; ss = 0d0
              if (sigma_arr(k00) > 0d0) then
                 sw = sw + w00; ss = ss + w00*sigma_arr(k00)
              endif
              if (sigma_arr(k01) > 0d0) then
                 sw = sw + w01; ss = ss + w01*sigma_arr(k01)
              endif
              if (sigma_arr(k10) > 0d0) then
                 sw = sw + w10; ss = ss + w10*sigma_arr(k10)
              endif
              if (sigma_arr(k11) > 0d0) then
                 sw = sw + w11; ss = ss + w11*sigma_arr(k11)
              endif
              if (sw > 0d0) then
                 if (ss > point_sigcut*sw) cycle
              endif
           endif
           pe_i = w00*n_e(k00)*theta_e(k00) + w01*n_e(k01)*theta_e(k01) &
                + w10*n_e(k10)*theta_e(k10) + w11*n_e(k11)*theta_e(k11)
           rho(i) = real(ne_i)
           p(i)   = real(max(pe_i, 0d0) / ne_i)

           u_ks(i)%data(1) = w00*u0(k00) + w01*u0(k01) + w10*u0(k10) + w11*u0(k11)
           u_ks(i)%data(2) = w00*u1(k00) + w01*u1(k01) + w10*u1(k10) + w11*u1(k11)
           u_ks(i)%data(3) = w00*u2(k00) + w01*u2(k01) + w10*u2(k10) + w11*u2(k11)
           u_ks(i)%data(4) = w00*u3(k00) + w01*u3(k01) + w10*u3(k10) + w11*u3(k11)

           b_ks(i)%data(1) = w00*b0(k00) + w01*b0(k01) + w10*b0(k10) + w11*b0(k11)
           b_ks(i)%data(2) = w00*b1(k00) + w01*b1(k01) + w10*b1(k10) + w11*b1(k11)
           b_ks(i)%data(3) = w00*b2(k00) + w01*b2(k01) + w10*b2(k10) + w11*b2(k11)
           b_ks(i)%data(4) = w00*b3(k00) + w01*b3(k01) + w10*b3(k10) + w11*b3(k11)
           filled(i) = .true.
        enddo
        end associate
        end associate

        ! Kerr-Schild (contravariant) -> Boyer-Lindquist (contravariant)
        u = uks2ubl(u_ks, rr, dble(a))
        b = uks2ubl(b_ks, rr, dble(a))
        call assign_metric(u, transpose(kerr_metric(rr, tth, dble(a))))
        call assign_metric(b, transpose(kerr_metric(rr, tth, dble(a))))

        ! Restore u.u = -1 and b.u = 0 at the interpolated points (see header),
        ! then |b| from the interpolated, projected field
        do i = 1, npts
           if (.not. filled(i)) cycle
           uu = u(i) * u(i)
           if (uu < 0d0) then
              u(i)%data = u(i)%data / sqrt(-uu)
              ub = b(i) * u(i)
              b(i)%data = b(i)%data + ub * u(i)%data
              bb = b(i) * b(i)
              bmag(i) = real(sqrt(max(bb, 0d0)))
           else
              rho(i) = 0.; p(i) = 0.      ! no timelike velocity: no emitter
           endif
        enddo
        end subroutine aperture_vals

        subroutine initialize_aperture_model(a, dfile)
        real(kind=8), intent(in) :: a
        character(len=*), intent(in) :: dfile
        integer :: k
        call read_aperture_data(dfile)
        ! the metric of the camera and of the fluid's u, b must be the same
        ! (for every binary of a list)
        do k = 1, naperture
           if (abs(abins(k)%asim - a) > 1d-4) then
              write(6,*) 'ERROR: aperture spin mismatch: the binary ', trim(abins(k)%path), &
                   ' has a =', abins(k)%asim, ', inputs.in (&geodata spin) a =', a
              stop 1
           endif
        enddo
        end subroutine initialize_aperture_model

        subroutine del_aperture_data()
        if (allocated(abins)) deallocate(abins)   ! frees the components too
        naperture = 0
        aperture_active = 1
        end subroutine del_aperture_data

      end module fluid_model_aperture
