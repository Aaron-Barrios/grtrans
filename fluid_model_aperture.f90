! fluid_model_aperture.f90
!
! The fluid model for ray tracing Aperture4 GR-PIC output with grtrans.
! Reads a stream binary of the reconstructed fluid (layout below; no HDF5 in
! Fortran) and bilinearly interpolates the reconstructed emitter variables onto ray
! points.  Mirrors fluid_model_harm.f90 but the grid is Kerr-Schild (r, theta),
! which shares the SAME r and theta as Boyer-Lindquist, so x0%data(2)/data(3)
! index the grid directly.
!
! Conventions (all verified against the Aperture4 source and the
! reconstruction that writes the binary):
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
!   * u^mu is the electrons' mean velocity: the electron current n_e u^mu is
!     interpolated and divided by the interpolated n_e (interpolating u^mu
!     itself would let a hot dense cell's electrons take the velocity of a
!     nearly empty neighbor).
!   * Interpolated u^mu is renormalized to u.u = -1 and b^mu projected
!     orthogonal to it.  A weighted average of unit timelike vectors pointing
!     in different directions is not a unit vector (|u.u| > 1), so bilinear
!     interpolation alone gives a wrong u^t and hence a wrong Doppler factor:
!     averaged over 2x2 cells, 11-13% of the points in the late dumps had
!     |u.u + 1| > 0.1.  |b| is taken from the
!     interpolated, projected b^mu (Aperture interpolates field components),
!     and theta_e as <n_e theta_e>/<n_e> (see aperture_vals); the binary's
!     bmag array is skipped.
!   * Binary layout (little-endian, stream access, no record markers):
!       int32   magic        'APL7' (0x41504C37) or 'APL8' (0x41504C38); older
!                            magics (APLF, APL2..APL6) are refused
!       int32   nr, nth      radial and polar cell counts
!       float64 spin         Kerr-Schild a of the fluid
!       float64 lambda_u     cm^-3 per code density unit (n_e and |b| -> CGS)
!       float64 mass_ratio   code particle mass / m_e (1 for a pair plasma;
!                            N = m_p/m_e when the particles are read as ions;
!                            |b| -> B carries sqrt(mass_ratio))
!       float64 theta_hybrid, ntheta_max, theta_synch   the emissivity dispatch
!                            (jnu_table.set_hybrid_dispatch)
!       float64 r_arr[nr]    radii of the cell CENTERS, log-spaced, increasing
!       float64 th_arr[nth]  polar angles of the cell centers, increasing
!       then 14 arrays of nth*nr float64 (k = ith*nr + ir): n_e, theta_e,
!       |b| (skipped), u0..u3, b0..b3 (electron frame, contravariant KS),
!       valid (1/0, skipped), beta_p (2P/b^2 of the fluid, 0 where none;
!       extra slot 20), sigma (b^2/(2 rho_0), 0 where no fluid; extra slot 21)
!       APL8 only: float64 sigma_cut, the last 8 bytes.
!     APL7: any sigma cut was applied per cell by the writer (cells above it
!     have n_e = theta_e = 0).  APL8: the cells are written uncut, and
!     aperture_vals applies sigma_cut at each ray sample (point_sigcut).
!   * The camera spin (grtrans's &geodata spin) must equal the spin the fluid
!     was reconstructed with (the header's asim, i.e. this binary's own spin):
!     a mismatch stops grtrans (initialize_aperture_model).
!   * Several binaries in one run: &harm fdfile naming a '.lst'
!     file (one binary path per line) loads every binary listed into abins;
!     with nmdot equal to their number, grtrans_driver traces binary m along
!     the same ray paths into the cameras of model parameter m (one geodesic
!     per pixel for all of them).  The binaries may have different grids, lambda_u, mass ratios
!     and sigma cuts, but must share the spin (the camera's) and the
!     emissivity dispatch.  A plain binary path traces that one binary.
!   * Slow light: when the lines of the '.lst' carry a time
!     after the path ("path t", that binary's simulation time in M, i.e.
!     r_g/c, in increasing order), the binaries are time slices of ONE fluid,
!     not separate models.  A line "tframe <t>" (required) sets the frame's
!     time; an optional line "tref <L>" sets the look-back that corresponds to
!     it, one value for every ray (default r_out = 1/uout of inputs.in, passed
!     in by pgrtrans with aperture_set_rout, so tframe is the fluid's time at
!     the hole's plane to within a few M; not each ray's first-sample radius,
!     which differs between rays by up to ~10% of r_out and would break the
!     common clock).  Each
!     ray sample is evaluated at the Kerr-Schild time t_emit = tframe + tref - L,
!     with L its Kerr-Schild look-back time from the ray's first sample
!     (ks_time), linearly between the two slices around t_emit
!     (aperture_vals_slow); outside the slices the nearest one holds, and
!     aperture_slow_report prints how many emitting samples that was.  The
!     look-back needs the geodesics' t and grtrans's t0 pass that puts all
!     rays on one observer clock: set nload > 1 in inputs.in (any value > 1;
!     the APERTURE model does not use it otherwise).  The slices must share
!     the grid, lambda_u, mass ratio, sigma cut and emissivity dispatch;
!     naperture is then 1 (one model) and nslice holds their number.
!     Example slices.lst:
!         tframe 1450
!         /path/to/fluid_01350.bin 1350
!         /path/to/fluid_01360.bin 1360
!         ...
!         /path/to/fluid_01550.bin 1550
!     With all slices the same binary the image equals the fast-light one.

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
         ! cm^-3 per code density unit, from the binary's header (fluid.f90
         ! converts with it)
         real(kind=8) :: lambda_u = 0d0
         ! code particle mass in units of m_e, from the header: 1 for the pair-
         ! plasma reading, N = m_p/m_e when the particles are read as ions; fluid.f90's
         ! |b| -> B carries sqrt(mass_ratio)
         real(kind=8) :: mass_ratio = 0d0
         ! the emissivity dispatch of the header (jnu_table.set_hybrid_dispatch)
         real(kind=8) :: th_hyb = 0d0, nth_max = 0d0, th_syn = 0d0
         ! sigma cut applied at each ray sample (APL8 binaries, written uncut
         ! with the cut in their last 8 bytes): a sample whose sigma,
         ! interpolated over the corners that hold fluid (sigma_arr > 0), exceeds
         ! it is left empty (no n_e, Theta_e, u, b: it neither emits nor absorbs),
         ! as ipole's geodesic sigma cut.  0: none here (an APL7
         ! binary was cut per cell by the writer, or not at all)
         real(kind=8) :: point_sigcut = 0d0
         real(kind=8), dimension(:), allocatable :: r_arr, th_arr
         ! ln r_arr(1) and the inverse mean steps in ln r and theta, for the
         ! O(1) first guess of grid_interval
         real(kind=8) :: lr1 = 0d0, dlr_inv = 0d0, dth_inv = 0d0
         real(kind=8), dimension(:), allocatable :: n_e, theta_e
         real(kind=8), dimension(:), allocatable :: u0, u1, u2, u3
         real(kind=8), dimension(:), allocatable :: b0, b1, b2, b3
         ! the plasma beta_p per cell from the binary (2P/b^2 of the fluid; 0 where
         ! no field), interpolated and emission-weighted for extra slot 20
         real(kind=8), dimension(:), allocatable :: beta_arr
         ! the magnetization sigma = b^2/(2 rho_0) per cell from the binary (0 where
         ! no fluid), interpolated and emission-weighted for extra slot 21
         real(kind=8), dimension(:), allocatable :: sigma_arr
         ! slow light: this slice's time [r_g/c] (the dump's)
         real(kind=8) :: t = 0d0
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
      ! slow light (see the header): the time slices in abins(1:nslice), the
      ! frame's time, the look-back that corresponds to it, and the samples
      ! with fluid evaluated outside the slices (clamped to the nearest)
      logical :: aperture_slow = .false.
      integer :: nslice = 0
      real(kind=8) :: aperture_tframe = 0d0, aperture_tref = 0d0
      logical :: have_tframe = .false., have_tref = .false.
      ! the camera's starting radius 1/uout (aperture_set_rout): tref's default
      real(kind=8) :: aperture_rout = 0d0
      integer(kind=8) :: nslow_samples = 0, nslow_early = 0, nslow_late = 0

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
           write(6,*) '  emissivity dispatch); rewrite it in the APL7/APL8 layout'
           write(6,*) '  (header of fluid_model_aperture.f90)'
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
        bn%lr1 = log(bn%r_arr(1))
        bn%dlr_inv = dble(nrr-1) / (log(bn%r_arr(nrr)) - bn%lr1)
        bn%dth_inv = dble(nthh-1) / (bn%th_arr(nthh) - bn%th_arr(1))
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
        real(kind=8), dimension(:), allocatable :: times
        logical :: timed
        integer :: k, n
        call del_aperture_data()
        timed = .false.
        n = len_trim(dfile)
        if (n >= 4) then
           if (dfile(n-3:n) == '.lst') then
              call read_binary_list(dfile, paths, times, timed)
           endif
        endif
        if (.not. allocated(paths)) then
           allocate(paths(1)); paths(1) = dfile
        endif
        n = size(paths)
        allocate(abins(n))
        do k = 1, n
           if (n > 1) write(6,'(a,i0,a,i0,2a)') ' aperture binary ', k, ' of ', &
                n, ': ', trim(paths(k))
           call read_one_binary(trim(paths(k)), abins(k))
           if (timed) abins(k)%t = times(k)
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
        if (timed) then
           ! slow light: one fluid in time slices
           if (n < 2) then
              write(6,*) 'ERROR: slow light needs at least two time slices in ', trim(dfile)
              stop 1
           endif
           if (.not. have_tframe) then
              write(6,*) 'ERROR: slow light: no "tframe <t>" line in ', trim(dfile)
              stop 1
           endif
           call check_slices(n)
           aperture_slow = .true.
           nslice = n
           naperture = 1
           if (have_tref) then
              write(6,'(a,i0,a,f10.3,a,f10.3,a,f10.3,a,f10.3)') ' slow light: ', n, ' time slices, t = ', &
                   abins(1)%t, ' to ', abins(n)%t, '; tframe ', aperture_tframe, ', tref ', aperture_tref
           else
              if (.not. (aperture_rout > 0d0)) then
                 write(6,*) 'ERROR: slow light: no "tref <L>" line and no r_out passed in (aperture_set_rout)'
                 stop 1
              endif
              write(6,'(a,i0,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a)') ' slow light: ', n, ' time slices, t = ', &
                   abins(1)%t, ' to ', abins(n)%t, '; tframe ', aperture_tframe, ', tref ', aperture_rout, &
                   ' (r_out)'
           endif
        else
           aperture_slow = .false.
           nslice = 0
           naperture = n
        endif
        if (allocated(times)) deallocate(times)
        aperture_active = 1
        ! the emissivity dispatch is the one the binary's header carries
        call set_hybrid_dispatch(abins(1)%th_hyb, abins(1)%nth_max, abins(1)%th_syn)
        end subroutine read_aperture_data

        ! The binary paths in the list file lfile, one per line (leading and
        ! trailing blanks removed; blank lines and '#' comments skipped).
        ! Slow light: "path t" lines (every line or none) and the lines
        ! "tframe <t>", "tref <L>" (see the header); timed is then true.
        subroutine read_binary_list(lfile, paths, times, timed)
        character(len=*), intent(in) :: lfile
        character(len=1024), dimension(:), allocatable, intent(out) :: paths
        real(kind=8), dimension(:), allocatable, intent(out) :: times
        logical, intent(out) :: timed
        character(len=1024) :: line
        integer :: u, ios, n, pass, sp, nt
        real(kind=8) :: tv
        open(newunit=u, file=lfile, form='formatted', status='old', action='read', iostat=ios)
        if (ios /= 0) then
           write(6,*) 'ERROR: cannot open aperture binary list ', trim(lfile)
           stop 1
        endif
        timed = .false.
        do pass = 1, 2
           n = 0
           nt = 0
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
              if (line(1:7) == 'tframe ') then
                 read(line(8:), *) aperture_tframe
                 have_tframe = .true.
                 cycle
              endif
              if (line(1:5) == 'tref ') then
                 read(line(6:), *) aperture_tref
                 have_tref = .true.
                 cycle
              endif
              n = n + 1
              sp = index(trim(line), ' ')
              if (sp > 0) then
                 nt = nt + 1
                 if (pass == 2) then
                    read(line(sp+1:), *, iostat=ios) tv
                    if (ios /= 0) then
                       write(6,*) 'ERROR: bad time in aperture binary list ', trim(lfile), ': ', trim(line)
                       stop 1
                    endif
                    paths(n) = line(1:sp-1)
                    times(n) = tv
                 endif
              else
                 if (pass == 2) paths(n) = line
              endif
           enddo
           if (pass == 1) then
              if (n == 0) then
                 write(6,*) 'ERROR: no binary in aperture binary list ', trim(lfile)
                 stop 1
              endif
              if (nt /= 0 .and. nt /= n) then
                 write(6,*) 'ERROR: aperture binary list ', trim(lfile), ': give a time on every line or none'
                 stop 1
              endif
              timed = (nt == n)
              allocate(paths(n), times(n))
              times = 0d0
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

        ! The interval idx in [1, n-1] with xx(idx) <= x < xx(idx+1) (idx = n-1
        ! at x = xx(n)), for xx(1) <= x <= xx(n) and an increasing xx: a first
        ! guess from the axis's mean step (f = (x - xx(1)) * step_inv; the axes
        ! are uniform in ln r and theta up to rounding), then moved to the
        ! interval that holds x.  The same interval as bisection, in O(1).
        subroutine grid_interval(xx, n, x, f, idx)
        integer, intent(in) :: n
        real(kind=8), intent(in) :: xx(n), x, f
        integer, intent(out) :: idx
        idx = max(1, min(int(f) + 1, n-1))
        do while (idx < n-1)
           if (xx(idx+1) > x) exit
           idx = idx + 1
        enddo
        do while (idx > 1)
           if (xx(idx) <= x) exit
           idx = idx - 1
        enddo
        end subroutine grid_interval

        ! Interpolation as Aperture4 does it (utils/interpolation.hpp, interp_t<1>
        ! in coord_policy_gr_ks_sph.hpp): first order (cloud-in-cell) in the
        ! code's own coordinates x1 = ln r and x2 = theta, with each value at
        ! the center of the cell it describes (the binary's r_arr, th_arr are the
        ! cell centers, not Aperture's grid.h5 first-native-cell labels), and applied
        ! to the raw components, never to derived quantities:
        !   * n_e and the electron pressure n_e theta_e are interpolated, and
        !     theta_e = <n_e theta_e>/<n_e> (a cell without electrons adds
        !     nothing to either, instead of pulling theta_e toward 0);
        !   * u^mu as the electron current n_e u^mu divided by the interpolated
        !     n_e (the electrons' mean velocity, as theta_e is their mean
        !     temperature), b^mu component by component (then u
        !     renormalized and b projected orthogonal to it), and |b| from the
        !     interpolated b^mu, not an interpolated |b|;
        !   * no validity gate and no renormalization of the weights, as in
        !     Aperture's block means: invalid cells (no timelike current) are
        !     zero in every quantity; electron-free cells and sigma-masked
        !     cells (cell mode) count as zero in n_e, n_e theta_e and n_e u^mu
        !     and keep a field, which enters neighboring points.
        ! beta_p and sigma (diagnostics) are averaged over the corners that hold
        ! fluid (sigma_arr > 0), with the weights renormalized over them.
        subroutine aperture_vals(x0, a, rho, p, b, u, bmag, beta, sigma)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real, dimension(size(x0)), intent(out) :: rho, p, bmag
        ! the binary's beta_p and sigma, interpolated (optional: grtrans's
        ! driver asks for them, the stand-alone test does not)
        real, dimension(size(x0)), intent(out), optional :: beta, sigma
        type (four_Vector), intent(out), dimension(size(x0)) :: u, b
        integer :: npts, i, lr, ur, lth, uth, k00, k01, k10, k11
        real(kind=8) :: rq, thq, rd, td, w00, w01, w10, w11
        real(kind=8) :: ne_i, pe_i, sw, ss, bs, wn00, wn01, wn10, wn11
        real(kind=8) :: rmin, rmax, thmin, thmax, rlo, rhi, thlo, thhi
        type (four_Vector), dimension(size(x0)) :: u_ks, b_ks
        real(kind=8), dimension(size(x0)) :: rr, tth
        logical, dimension(size(x0)) :: filled
        real(kind=8) :: uu, ub, bb
        real(kind=8), dimension(size(x0),10) :: metric

        if (aperture_slow) then
           call aperture_vals_slow(x0, a, rho, p, b, u, bmag, beta, sigma)
           return
        endif
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
        ! the grid's own extent: half a cell beyond the outermost cell centers
        ! (half a step in ln r and in theta; at the axis this is theta = 0 and pi)
        rlo = rmin * sqrt(r_arr(1) / r_arr(2)); rhi = rmax * sqrt(r_arr(nr) / r_arr(nr-1))
        thlo = thmin - 0.5d0 * (th_arr(2) - th_arr(1))
        thhi = thmax + 0.5d0 * (th_arr(nth) - th_arr(nth-1))
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
           ! beyond the grid there is no fluid; between the outermost cell
           ! centers and the grid's edge (the half cell at the axis and at the
           ! inner and outer radii) the edge cell's values hold, as inside a
           ! cell.
           ! Only the interpolation uses the clamped position; the metric
           ! below uses the true one.
           if (rq < rlo .or. rq > rhi .or. thq < thlo .or. thq > thhi) then
              cycle
           endif
           rq = min(max(rq, rmin), rmax)
           thq = min(max(thq, thmin), thmax)
           call grid_interval(r_arr, nr, rq, (log(rq) - bn%lr1)*bn%dlr_inv, lr);  ur = lr + 1
           call grid_interval(th_arr, nth, thq, (thq - thmin)*bn%dth_inv, lth); uth = lth + 1
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

           ! the electron current n_e u^mu over n_e (ne_i > 0 here)
           wn00 = w00*n_e(k00); wn01 = w01*n_e(k01); wn10 = w10*n_e(k10); wn11 = w11*n_e(k11)
           u_ks(i)%data(1) = (wn00*u0(k00) + wn01*u0(k01) + wn10*u0(k10) + wn11*u0(k11))/ne_i
           u_ks(i)%data(2) = (wn00*u1(k00) + wn01*u1(k01) + wn10*u1(k10) + wn11*u1(k11))/ne_i
           u_ks(i)%data(3) = (wn00*u2(k00) + wn01*u2(k01) + wn10*u2(k10) + wn11*u2(k11))/ne_i
           u_ks(i)%data(4) = (wn00*u3(k00) + wn01*u3(k01) + wn10*u3(k10) + wn11*u3(k11))/ne_i

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
        metric = kerr_metric(rr, tth, dble(a))
        call assign_metric(u, transpose(metric))
        call assign_metric(b, transpose(metric))

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

        ! Slow light: the n slices must be one fluid on one grid, in increasing time.
        subroutine check_slices(n)
        integer, intent(in) :: n
        integer :: k
        do k = 2, n
           if (abins(k)%nr /= abins(1)%nr .or. abins(k)%nth /= abins(1)%nth) then
              write(6,*) 'ERROR: slow light: slice ', trim(abins(k)%path), ' has another grid'
              stop 1
           endif
           if (any(abins(k)%r_arr /= abins(1)%r_arr) .or. any(abins(k)%th_arr /= abins(1)%th_arr) &
                .or. abins(k)%lambda_u /= abins(1)%lambda_u .or. abins(k)%mass_ratio /= abins(1)%mass_ratio &
                .or. abins(k)%point_sigcut /= abins(1)%point_sigcut) then
              write(6,*) 'ERROR: slow light: slice ', trim(abins(k)%path), &
                   ' differs from the first in grid, lambda_u, mass ratio or sigma cut'
              stop 1
           endif
           if (.not. (abins(k)%t > abins(k-1)%t)) then
              write(6,*) 'ERROR: slow light: slice times must increase (', abins(k-1)%t, abins(k)%t, ')'
              stop 1
           endif
        enddo
        end subroutine check_slices

        ! Kerr-Schild time of the Boyer-Lindquist event (r, t), up to a constant
        ! (M = 1): t + ln Delta + ln((r - r+)/(r - r-)) / sqrt(1 - a^2), the
        ! integral of 2r/Delta (kerr.f90's bl2ks_time has half that last
        ! coefficient; it is not used here)
        real(kind=8) function ks_time(r, t, a)
        real(kind=8), intent(in) :: r, t, a
        real(kind=8) :: s, rp, rm, rr
        s = sqrt(1d0 - a*a)
        rp = 1d0 + s; rm = 1d0 - s
        rr = max(r, rp*(1d0 + 1d-12))
        ks_time = t + log(rr*rr - 2d0*rr + a*a) + log((rr - rp)/(rr - rm)) / s
        end function ks_time

        ! aperture_vals for slow light: the same interpolation in space as
        ! aperture_vals (see there), done in the two time slices around each
        ! sample's emission time and combined linearly in time, for the raw
        ! quantities: n_e, n_e theta_e, n_e u^mu, b^mu, and beta_p and sigma
        ! over the corners that hold fluid.  The sigma cut tests the
        ! time-interpolated sigma.
        subroutine aperture_vals_slow(x0, a, rho, p, b, u, bmag, beta, sigma)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real, dimension(size(x0)), intent(out) :: rho, p, bmag
        real, dimension(size(x0)), intent(out), optional :: beta, sigma
        type (four_Vector), intent(out), dimension(size(x0)) :: u, b
        integer :: npts, i, lr, ur, lth, uth, k00, k01, k10, k11, ks, q, kk
        real(kind=8) :: rq, thq, rd, td, w00, w01, w10, w11, wt(2), tem, l0, tref
        real(kind=8) :: ne_i, pe_i, sw, ss, bs, sg, bt, sgw, cur(4), bb4(4)
        real(kind=8) :: rmin, rmax, thmin, thmax, rlo, rhi, thlo, thhi
        type (four_Vector), dimension(size(x0)) :: u_ks, b_ks
        real(kind=8), dimension(size(x0)) :: rr, tth
        logical, dimension(size(x0)) :: filled
        real(kind=8) :: uu, ub, bb
        real(kind=8), dimension(size(x0),10) :: metric
        integer(kind=8) :: nsmp, nearly, nlate
        integer :: clamp

        npts = size(x0)
        if (npts > 1) then
           if (all(x0%data(1) == 0d0)) then
              write(6,*) 'ERROR: slow light needs the geodesics'' t: set nload > 1 in inputs.in'
              stop 1
           endif
        endif
        associate (nr => abins(1)%nr, nth => abins(1)%nth, r_arr => abins(1)%r_arr, &
             th_arr => abins(1)%th_arr, point_sigcut => abins(1)%point_sigcut)
        rmin = r_arr(1); rmax = r_arr(nr)
        thmin = th_arr(1); thmax = th_arr(nth)
        rlo = rmin * sqrt(r_arr(1) / r_arr(2)); rhi = rmax * sqrt(r_arr(nr) / r_arr(nr-1))
        thlo = thmin - 0.5d0 * (th_arr(2) - th_arr(1))
        thhi = thmax + 0.5d0 * (th_arr(nth) - th_arr(nth-1))
        rr = x0%data(2)
        tth = x0%data(3)
        ! the look-back time is measured from the ray's first sample at t = 0
        ! (the t0 pass puts every ray's samples on one observer clock)
        l0 = ks_time(rr(1), 0d0, dble(a))
        tref = aperture_tref
        if (.not. have_tref) tref = aperture_rout
        nsmp = 0; nearly = 0; nlate = 0

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
           if (rq < rlo .or. rq > rhi .or. thq < thlo .or. thq > thhi) cycle
           rq = min(max(rq, rmin), rmax)
           thq = min(max(thq, thmin), thmax)
           call grid_interval(r_arr, nr, rq, (log(rq) - abins(1)%lr1)*abins(1)%dlr_inv, lr);  ur = lr + 1
           call grid_interval(th_arr, nth, thq, (thq - thmin)*abins(1)%dth_inv, lth); uth = lth + 1
           rd = log(rq / r_arr(lr)) / log(r_arr(ur) / r_arr(lr))
           td = (thq - th_arr(lth)) / (th_arr(uth) - th_arr(lth))
           k00 = (lth-1)*nr + lr
           k01 = (uth-1)*nr + lr
           k10 = (lth-1)*nr + ur
           k11 = (uth-1)*nr + ur
           w00 = (1d0-rd)*(1d0-td)
           w01 = (1d0-rd)*td
           w10 = rd*(1d0-td)
           w11 = rd*td
           ! the emission time and the two slices around it
           tem = aperture_tframe + tref - (l0 - ks_time(x0(i)%data(2), x0(i)%data(1), dble(a)))
           clamp = 0
           if (tem <= abins(1)%t) then
              ks = 1; wt = (/1d0, 0d0/)
              if (tem < abins(1)%t) clamp = -1
           else if (tem >= abins(nslice)%t) then
              ks = nslice - 1; wt = (/0d0, 1d0/)
              if (tem > abins(nslice)%t) clamp = 1
           else
              ks = 1
              do while (abins(ks+1)%t <= tem)
                 ks = ks + 1
              enddo
              wt(2) = (tem - abins(ks)%t) / (abins(ks+1)%t - abins(ks)%t)
              wt(1) = 1d0 - wt(2)
           endif
           ne_i = 0d0; pe_i = 0d0; cur = 0d0; bb4 = 0d0; sg = 0d0; bt = 0d0; sgw = 0d0
           do q = 1, 2
              if (wt(q) == 0d0) cycle
              kk = ks + q - 1
              associate (n_e => abins(kk)%n_e, theta_e => abins(kk)%theta_e, u0 => abins(kk)%u0, &
                   u1 => abins(kk)%u1, u2 => abins(kk)%u2, u3 => abins(kk)%u3, b0 => abins(kk)%b0, &
                   b1 => abins(kk)%b1, b2 => abins(kk)%b2, b3 => abins(kk)%b3, &
                   beta_arr => abins(kk)%beta_arr, sigma_arr => abins(kk)%sigma_arr)
              ne_i = ne_i + wt(q)*(w00*n_e(k00) + w01*n_e(k01) + w10*n_e(k10) + w11*n_e(k11))
              pe_i = pe_i + wt(q)*(w00*n_e(k00)*theta_e(k00) + w01*n_e(k01)*theta_e(k01) &
                   + w10*n_e(k10)*theta_e(k10) + w11*n_e(k11)*theta_e(k11))
              cur(1) = cur(1) + wt(q)*(w00*n_e(k00)*u0(k00) + w01*n_e(k01)*u0(k01) &
                   + w10*n_e(k10)*u0(k10) + w11*n_e(k11)*u0(k11))
              cur(2) = cur(2) + wt(q)*(w00*n_e(k00)*u1(k00) + w01*n_e(k01)*u1(k01) &
                   + w10*n_e(k10)*u1(k10) + w11*n_e(k11)*u1(k11))
              cur(3) = cur(3) + wt(q)*(w00*n_e(k00)*u2(k00) + w01*n_e(k01)*u2(k01) &
                   + w10*n_e(k10)*u2(k10) + w11*n_e(k11)*u2(k11))
              cur(4) = cur(4) + wt(q)*(w00*n_e(k00)*u3(k00) + w01*n_e(k01)*u3(k01) &
                   + w10*n_e(k10)*u3(k10) + w11*n_e(k11)*u3(k11))
              bb4(1) = bb4(1) + wt(q)*(w00*b0(k00) + w01*b0(k01) + w10*b0(k10) + w11*b0(k11))
              bb4(2) = bb4(2) + wt(q)*(w00*b1(k00) + w01*b1(k01) + w10*b1(k10) + w11*b1(k11))
              bb4(3) = bb4(3) + wt(q)*(w00*b2(k00) + w01*b2(k01) + w10*b2(k10) + w11*b2(k11))
              bb4(4) = bb4(4) + wt(q)*(w00*b3(k00) + w01*b3(k01) + w10*b3(k10) + w11*b3(k11))
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
                 sg = sg + wt(q)*ss/sw; bt = bt + wt(q)*bs/sw; sgw = sgw + wt(q)
              endif
              end associate
           enddo
           if (sgw > 0d0) then
              if (present(beta)) beta(i) = real(bt / sgw)
              if (present(sigma)) sigma(i) = real(sg / sgw)
           endif
           if (.not. (ne_i > 0d0)) cycle
           nsmp = nsmp + 1
           if (clamp < 0) nearly = nearly + 1
           if (clamp > 0) nlate = nlate + 1
           if (point_sigcut > 0d0 .and. sgw > 0d0) then
              if (sg / sgw > point_sigcut) cycle
           endif
           rho(i) = real(ne_i)
           p(i)   = real(max(pe_i, 0d0) / ne_i)
           u_ks(i)%data = cur / ne_i
           b_ks(i)%data = bb4
           filled(i) = .true.
        enddo
        end associate
!$omp atomic update
        nslow_samples = nslow_samples + nsmp
!$omp atomic update
        nslow_early = nslow_early + nearly
!$omp atomic update
        nslow_late = nslow_late + nlate

        u = uks2ubl(u_ks, rr, dble(a))
        b = uks2ubl(b_ks, rr, dble(a))
        metric = kerr_metric(rr, tth, dble(a))
        call assign_metric(u, transpose(metric))
        call assign_metric(b, transpose(metric))
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
              rho(i) = 0.; p(i) = 0.
           endif
        enddo
        end subroutine aperture_vals_slow

        ! Slow light: the camera's starting radius 1/uout, tref's default (pgrtrans
        ! calls this before loading the fluid)
        subroutine aperture_set_rout(rout)
        real(kind=8), intent(in) :: rout
        aperture_rout = rout
        end subroutine aperture_set_rout

        ! Slow light: how many samples with electrons were evaluated, and how
        ! many of them fell before the first slice or after the last (clamped)
        subroutine aperture_slow_report()
        if (.not. aperture_slow) return
        write(6,'(a,i0,a,i0,a,i0,a)') ' slow light: ', nslow_samples, ' samples with electrons; ', &
             nslow_early, ' of them before the first slice and ', nslow_late, &
             ' after the last were clamped to it'
        end subroutine aperture_slow_report

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
        aperture_slow = .false.; nslice = 0
        nslow_samples = 0; nslow_early = 0; nslow_late = 0
        have_tframe = .false.; have_tref = .false.
        end subroutine del_aperture_data

      end module fluid_model_aperture
