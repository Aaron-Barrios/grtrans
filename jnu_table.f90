! jnu_table.f90
!
! Stage 4: trilinear interpolation of the precomputed emissivity lookup table
! produced by make_jnu_table.py.
!
! The exact jnu integrator is too slow for per-ray-step evaluation, so we
! exploit exact linearity: for a thermal (Maxwell-Juettner) plasma,
!
!       j_nu(nu, B, n_e) = n_e * B * F(x, theta_B, theta_e),
!       x = nu / nu_c,   nu_c = e B / (2 pi m_e c),
!
! where F is a dimensionless (in the sense of being a pure function of the
! dimensionless x, theta_B, theta_e) reduced emissivity.  The table stores
!
!       logF(i,j,k) = log10( F(x_i, thetaB_j, thetaE_k) )
!
! over a log-spaced x grid, a linear theta_B grid in [0.02, pi/2], and a
! log-spaced theta_e grid.  F was computed at B=1 G, n_e=1 cm^-3, so F is
! numerically j_nu in cgs; the table file stores log10(F).
!
! Interpolation is trilinear in (log x, theta_B, log theta_e), with values
! clamped at the theta_B and theta_e edges and at x_min (saturating, not
! extrapolating).  Above x_max (10^5) the table is continued with the shape of
! the ultra-relativistic thermal-synchrotron fit of Leung, Gammie & Noble
! (2011), anchored to the table value at x_max:
!       F(x) = F(x_max) * f(X)/f(X_max),  f(X) = (X^1/2 + 2^(11/12) X^1/6)^2 exp(-X^1/3),
!       X = x / ((2/9) theta_e^2 sin theta_B),
! instead of the saturating edge value, which grossly overestimated j_nu
! there.  Below x_max nothing changes.  The theta_B
! symmetry j_I(theta_B) = j_I(pi - theta_B) is exploited so only [0, pi/2] is
! stored.
!
! Absorption is added on the fly from Kirchhoff's law (LTE):
!       alpha_nu = j_nu / B_nu(T_e),   T_e = theta_e m_e c^2 / k_B.
!
! File layout (Fortran stream binary, written by make_jnu_table.py):
!   i32 magic (0x4A4E5554 'JNUT'), i32 nx, i32 nthb, i32 nte
!   f64 x_axis[nx], f64 thB_axis[nthb], f64 te_axis[nte]
!   f64 logF[nx*nthb*nte]   numpy C-order flatten of logF(i_x, j_thB, k_te):
!                           offset (0-based) = (i*nthb + j)*nte + k
!                           i.e. x SLOWEST, theta_e FASTEST.  (2026-09-23: the
!                           indexing below was transposed relative to this --
!                           it assumed x fastest -- and has been corrected to
!                           match the bytes; verified against the jnu binary:
!                           0.46 dex median |dlog10 F| after, 34.8 dex before.)

      module jnu_table

      use phys_constants, only: c, e, m, k, h, pi
      use cyclotron_emis, only: cyclotron_emissivity

      implicit none

      integer :: nx, nthb, nte
      real(kind=8), dimension(:), allocatable :: x_axis, thB_axis, te_axis
      real(kind=8), dimension(:), allocatable :: logF
      ! log of x_axis: the x axis may be piecewise log-uniform (denser around
      ! the cyclotron harmonics), so its interval is found by bisection
      real(kind=8), dimension(:), allocatable :: lx_axis
      real(kind=8) :: lxmin, lxmax, ltemin, ltemax, dxinv, dteinv, dthbinv
      logical :: loaded = .false.

      integer, parameter :: JMAGIC = int(z'4A4E5554')

      interface init_jnu_table
        module procedure init_jnu_table
      end interface

      interface del_jnu_table
        module procedure del_jnu_table
      end interface

      interface jnu_table_emissivity
        module procedure jnu_table_emissivity
      end interface

      interface hybrid_emissivity
        module procedure hybrid_emissivity
      end interface

      ! set by the first clamped evaluation, so each edge warns once per run
      ! (shared by the threads; see jnu_table_emissivity)
      logical :: warned_hot = .false., warned_cold = .false.

      ! hybrid_emissivity's dispatch.  The defaults are the cyclotron/table
      ! split alone; the APERTURE model sets all three from the fluid binary's
      ! header (pipeline_config.THETA_HYBRID, NTHETA_MAX, THETA_SYNCH), so the
      ! Fortran holds no copy of the pipeline's choices:
      !   theta_e < th_hybrid and (nu/nu_B) theta_e <= ntheta_max: cyclotron
      !   theta_e >= th_synch: thermal synchrotron fit (Leung et al. 2011)
      !   otherwise: the exact table
      ! th_synch <= 0 means "beyond the table": the fit only above the table's
      ! largest theta_e (te_axis(nte)), where the table would clamp
      real(kind=8) :: th_hybrid = 5d-3, ntheta_max = 0.1d0, th_synch = huge(1d0)

      ! The lookup interpolates log F linearly in (log x, theta_B, log theta_e)
      ! below theta_e = th_cubic and with a tensor-product Catmull-Rom cubic over
      ! the 4 x 4 x 4 surrounding entries from th_cubic up.  log F curves
      ! downward between the entries, so the linear lookup is low: by 2.0-2.4%
      ! in the emission-weighted sum of relabeled 230 GHz gas (Theta_e >= 5),
      ! two thirds of it from the theta_e axis; the cubic brings that to
      ! +0.04-0.13% (docs/history.tex, Stage 79).  Below theta_e ~ 1 the cyclotron
      ! harmonics make log F too rough for a cubic on this grid (it overshoots
      ! at the fundamental), so the linear lookup is kept there.  An edge
      ! interval of an axis is linear in that axis, and a neighborhood holding
      ! a floor entry (log F <= LF_FLOOR, i.e. no emission) falls back to the
      ! linear lookup.  th_cubic <= 0: the cubic everywhere; huge: never.
      real(kind=8) :: th_cubic = 1d0
      real(kind=8), parameter :: LF_FLOOR = -290d0

      contains

        subroutine init_jnu_table(fname)
        character(len=*), intent(in) :: fname
        character(len=300) :: table_path
        integer :: magic, ios, stat
        real(kind=8) :: nuc1
! The grtrans pixel loop is OpenMP-parallel and every thread calls
! initialize_emissivity (hence this) once per ray.  Guard the one-time load
! with a critical section so only the first thread reads/allocates the table
! and the rest wait and see `loaded == .true.`.
!$omp critical (jnu_table_load)
        if (.not. loaded) then
           ! Allow override via JNU_TABLE_PATH, default to fname argument.
           call get_environment_variable('JNU_TABLE_PATH', table_path, status=stat)
           if (stat /= 0 .or. len_trim(table_path) == 0) table_path = fname
           open(unit=8, file=trim(table_path), form='unformatted', access='stream', &
                status='old', action='read', iostat=ios)
           if (ios /= 0) then
              write(6,*) 'ERROR: cannot open jnu table ', trim(table_path)
              stop 1
           endif
           read(8) magic, nx, nthb, nte
           if (magic /= JMAGIC) then
              write(6,*) 'ERROR: bad magic in jnu table ', trim(table_path)
              close(8); stop 1
           endif
           allocate(x_axis(nx), thB_axis(nthb), te_axis(nte))
           allocate(logF(nx*nthb*nte))
           read(8) x_axis, thB_axis, te_axis, logF
           close(8)

           allocate(lx_axis(nx))
           lx_axis = log(x_axis)
           lxmin = log(x_axis(1)); lxmax = log(x_axis(nx))
           ltemin = log(te_axis(1)); ltemax = log(te_axis(nte))
           dxinv = 1d0/(lxmax-lxmin)
           dteinv = 1d0/(ltemax-ltemin)
           dthbinv = 1d0/(thB_axis(nthb)-thB_axis(1))

           nuc1 = e/(2d0*pi*m*c)   ! nu_c at B = 1 Gauss
           loaded = .true.
           write(6,*) 'loaded jnu table: ', trim(table_path)
           write(6,*) '  nx,nthb,nte =', nx, nthb, nte
           write(6,*) '  x    in [', x_axis(1), ',', x_axis(nx), ']  (nu/nu_c)'
           write(6,*) '  thB  in [', thB_axis(1), ',', thB_axis(nthb), '] rad'
           write(6,*) '  te   in [', te_axis(1), ',', te_axis(nte), ']'
           write(6,*) '  nu_c(B=1G) =', nuc1, ' Hz'
        endif
!$omp end critical (jnu_table_load)
        end subroutine init_jnu_table

        subroutine del_jnu_table()
        ! Keep the table resident for the whole run: it is looked up per ray
        ! step and reloading per ray would dominate the cost.  Intentionally
        ! a no-op; the process exit releases the memory.
        end subroutine del_jnu_table

        subroutine jnu_table_emissivity(nu, n_e, B, theta_B, theta_e, j, alpha)
        real(kind=8), dimension(:), intent(in) :: nu, n_e, B, theta_B, theta_e
        real(kind=8), dimension(:), intent(out) :: j, alpha
        real(kind=8) :: nuc, x, lx, te, thb
        real(kind=8) :: fx, ft, fb, lf
        real(kind=8) :: w000, w100, w010, w110, w001, w101, w011, w111
        real(kind=8) :: TeK, bnu, xs, fext
        real(kind=8) :: wx(4), wb(4), wt(4), lfc, v
        logical :: warn, flo
        integer :: i, ix, ith, ite, k000, k100, k010, k110, k001, k101, k011, k111
        integer :: a1, b1, c1

        do i = 1, size(nu)
           nuc = e*B(i)/(2d0*pi*m*c)
           x = nu(i)/nuc
           ! map theta_B into [0, pi/2] using j_I symmetry
           thb = theta_B(i)
           if (thb < 0d0) thb = -thb
           thb = min(thb, pi - thb)
           te = max(theta_e(i), 1d-300)

           ! the table's theta_e range is finite and the lookup clamps at its
           ! edges, so emitting gas outside it would silently get the edge
           ! value: say so, once per run for each edge.  The flags are shared by
           ! all threads: they are read atomically and set in a critical section
           ! (one warning per edge, not one per thread).  Through
           ! hybrid_emissivity the hot edge is reached only when THETA_SYNCH is
           ! set above the table's largest theta_e (otherwise that gas goes to
           ! the synchrotron fit)
           if (n_e(i) > 0d0 .and. theta_e(i) > te_axis(nte)) then
!$omp atomic read
              warn = warned_hot
              if (.not. warn) then
!$omp critical (jnu_table_warn)
                 if (.not. warned_hot) then
                    write(6,*) 'WARNING jnu_table: theta_e =', theta_e(i), ' above the table', &
                         ' maximum', te_axis(nte), '; j_nu is clamped to the edge (lower', &
                         ' THETA_SYNCH to at most the table maximum, or extend the table:', &
                         ' make_jnu_table.py --extend)'
!$omp atomic write
                    warned_hot = .true.
                 endif
!$omp end critical (jnu_table_warn)
              endif
           endif
           if (n_e(i) > 0d0 .and. theta_e(i) > 0d0 .and. theta_e(i) < te_axis(1)) then
!$omp atomic read
              warn = warned_cold
              if (.not. warn) then
!$omp critical (jnu_table_warn)
                 if (.not. warned_cold) then
                    write(6,*) 'WARNING jnu_table: theta_e =', theta_e(i), ' below the table', &
                         ' minimum', te_axis(1), '; j_nu is clamped to the edge (extend the', &
                         ' table: make_jnu_table.py --extend)'
!$omp atomic write
                    warned_cold = .true.
                 endif
!$omp end critical (jnu_table_warn)
              endif
           endif
           lx = log(max(x, tiny(1d0)))
           fb = clamp01((thb - thB_axis(1))*dthbinv)
           ft = clamp01((log(te) - ltemin)*dteinv)

           ! x: the interval by bisection (the axis may be piecewise uniform)
           call x_interval(min(max(lx, lxmin), lxmax), ix, fx)
           ith = min(int(fb*(nthb-1))+1, nthb-1)
           fb = fb*(nthb-1) - (ith-1)
           ite = min(int(ft*(nte-1))+1, nte-1)
           ft = ft*(nte-1) - (ite-1)

           ! on-disk layout (0-based): offset = (i*nthb + j)*nte + k,
           ! i = x, j = theta_B, k = theta_e  (x slowest, theta_e fastest)
           k000 = ((ix-1)*nthb + (ith-1))*nte + ite
           k100 = k000 + nthb*nte
           k010 = k000 + nte
           k110 = k010 + nthb*nte
           k001 = k000 + 1
           k101 = k100 + 1
           k011 = k010 + 1
           k111 = k011 + nthb*nte

           w000 = (1d0-fx)*(1d0-fb)*(1d0-ft)
           w100 = fx*(1d0-fb)*(1d0-ft)
           w010 = (1d0-fx)*fb*(1d0-ft)
           w110 = fx*fb*(1d0-ft)
           w001 = (1d0-fx)*(1d0-fb)*ft
           w101 = fx*(1d0-fb)*ft
           w011 = (1d0-fx)*fb*ft
           w111 = fx*fb*ft

           lf = w000*logF(k000) + w100*logF(k100) + w010*logF(k010) &
              + w110*logF(k110) + w001*logF(k001) + w101*logF(k101) &
              + w011*logF(k011) + w111*logF(k111)

           ! hot gas: the tensor-product cubic over entries ix-1..ix+2 (etc.)
           if (te >= th_cubic) then
              call cr_weights(fx, ix, nx, wx)
              ! a cubic in x only where the three spacings around the interval
              ! are equal (inside one log-uniform run); otherwise linear in x
              if (ix > 1 .and. ix < nx-1) then
                 if (abs(lx_axis(ix)-lx_axis(ix-1) - (lx_axis(ix+1)-lx_axis(ix))) > &
                      1d-9*(lx_axis(ix+1)-lx_axis(ix)) .or. &
                      abs(lx_axis(ix+2)-lx_axis(ix+1) - (lx_axis(ix+1)-lx_axis(ix))) > &
                      1d-9*(lx_axis(ix+1)-lx_axis(ix))) wx = (/0d0, 1d0-fx, fx, 0d0/)
              endif
              call cr_weights(fb, ith, nthb, wb)
              call cr_weights(ft, ite, nte, wt)
              lfc = 0d0; flo = .false.
              do a1 = 1, 4
                 if (wx(a1) == 0d0) cycle
                 do b1 = 1, 4
                    if (wb(b1) == 0d0) cycle
                    do c1 = 1, 4
                       if (wt(c1) == 0d0) cycle
                       v = logF(((ix+a1-3)*nthb + (ith+b1-3))*nte + ite+c1-2)
                       if (v <= LF_FLOOR) flo = .true.
                       lfc = lfc + wx(a1)*wb(b1)*wt(c1)*v
                    enddo
                 enddo
              enddo
              if (.not. flo) lf = lfc
           endif

           j(i) = n_e(i)*B(i)*10d0**lf
           ! beyond the table in x: continue with the synchrotron fall-off
           ! (only for finite x: B = 0 gives x = Inf and j = 0 already)
           if (lx > lxmax .and. j(i) > 0d0 .and. x < 1d300) then
              xs = (2d0/9d0)*te*te*max(sin(thb), 1d-6)
              fext = lsynch(min(x/xs, 1d30)) - lsynch(x_axis(nx)/xs)
              j(i) = j(i)*exp(max(fext, -700d0))
           endif

           ! Kirchhoff absorption: alpha = j / B_nu(T_e), T_e = theta_e m_e c^2/k
           TeK = theta_e(i)*m*c*c/k
           bnu = planck(TeK, nu(i))
           if (bnu > 0d0) then
              alpha(i) = j(i)/bnu
           else
              alpha(i) = 0d0
           endif
        enddo
        end subroutine jnu_table_emissivity

        ! the interval i (entries i, i+1) of lx_axis holding lv, and the
        ! fraction f within it; lv is clamped to the axis by the caller
        subroutine x_interval(lv, i, f)
        real(kind=8), intent(in) :: lv
        integer, intent(out) :: i
        real(kind=8), intent(out) :: f
        integer :: lo, hi, mid
        lo = 1; hi = nx
        do while (hi - lo > 1)
           mid = (lo + hi)/2
           if (lx_axis(mid) <= lv) then
              lo = mid
           else
              hi = mid
           endif
        enddo
        i = lo
        f = (lv - lx_axis(i))/(lx_axis(i+1) - lx_axis(i))
        end subroutine x_interval

        ! Catmull-Rom weights over the entries i-1, i, i+1, i+2 of an axis of n
        ! entries for the fraction f within interval i (entries i, i+1); in the
        ! first and last interval, linear weights on i and i+1 (the others 0)
        subroutine cr_weights(f, i, n, w)
        real(kind=8), intent(in) :: f
        integer, intent(in) :: i, n
        real(kind=8), intent(out) :: w(4)
        if (i <= 1 .or. i >= n-1) then
           w = (/0d0, 1d0-f, f, 0d0/)
        else
           w(1) = 0.5d0*(-f + 2d0*f*f - f**3)
           w(2) = 0.5d0*(2d0 - 5d0*f*f + 3d0*f**3)
           w(3) = 0.5d0*(f + 4d0*f*f - 3d0*f**3)
           w(4) = 0.5d0*(-f*f + f**3)
        endif
        end subroutine cr_weights

        ! log of the Leung et al. (2011) thermal-synchrotron spectral shape
        real(kind=8) function lsynch(xx)
        real(kind=8), intent(in) :: xx
        lsynch = 2d0*log(sqrt(xx) + 2d0**(11d0/12d0)*xx**(1d0/6d0)) - xx**(1d0/3d0)
        end function lsynch

        real(kind=8) function clamp01(x)
        real(kind=8), intent(in) :: x
        clamp01 = min(max(x, 0d0), 1d0)
        end function clamp01

        real(kind=8) function planck(T, nu)
        real(kind=8), intent(in) :: T, nu
        real(kind=8) :: x
        x = h*nu/(k*T)
        if (x < 1d-6) then
           planck = 2d0*nu*nu*k*T/(c**2)
        else
           planck = 2d0*h*nu/ (c**2) * nu*nu &
                    / (exp(x)-1d0)
        endif
        end function planck

        ! Hybrid dispatch.  The analytic cyclotron formula (cyclotron_emis) is
        ! the leading-order low-temperature expansion of the harmonic sum, valid
        ! only when BOTH theta_e << 1 and n theta_e << 1, with n = nu/nu_B the
        ! harmonic the point is observed at (nu_B = e B/(2 pi m_e c), the
        ! cyclotron frequency, called nu_c in this file).  hybrid_emissivity
        ! uses it only where
        !     theta_e < th_hybrid (5e-3)   and   (nu/nu_B) * theta_e <= ntheta_max (0.1)
        ! (and not where n_e <= 0 or theta_e <= 0, which get j = alpha = 0
        ! from every branch); the thermal-synchrotron fit
        ! (synch_emissivity) where theta_e >= th_synch, or, for th_synch <= 0,
        ! beyond the table's largest theta_e; and the exact Maxwell-Juettner
        ! table (theta_e from ~7e-5 up) everywhere else.  Each branch is
        ! evaluated only on the points that use it, through pack/unpack.
        subroutine set_hybrid_dispatch(theta_hybrid, ntheta, theta_synch)
        real(kind=8), intent(in) :: theta_hybrid, ntheta, theta_synch
        th_hybrid = theta_hybrid; ntheta_max = ntheta; th_synch = theta_synch
        ! the fit is the ultrarelativistic form (Leung et al. 2011): below
        ! theta_e ~ 1 its shape, not only K2, departs from the exact sum
        if (th_synch > 0d0 .and. th_synch < 1d0) then
           write(6,*) 'WARNING: THETA_SYNCH =', th_synch, ' < 1: the thermal-synchrotron', &
                ' fit is ultrarelativistic and inaccurate for theta_e of order 1 or below'
        endif
        if (th_synch > 0d0) then
           write(6,*) 'emissivity dispatch: cyclotron below theta_e =', th_hybrid, &
                ' (n theta_e <=', ntheta_max, '), thermal synchrotron from', th_synch
        else
           write(6,*) 'emissivity dispatch: cyclotron below theta_e =', th_hybrid, &
                ' (n theta_e <=', ntheta_max, '), thermal synchrotron beyond the table'
        endif
        end subroutine set_hybrid_dispatch

        ! The thermal-synchrotron emissivity of Leung, Gammie & Noble (2011),
        ! absolutely normalized (the shape lsynch is the table's x > x_max
        ! continuation), with Kirchhoff absorption.  For theta_e >= th_synch,
        ! where it agrees with jnu_calculator to 0.2-2% (history.tex, Stage 24):
        !   j = sqrt(2) pi e^2 n_e nu_s /(3 c K2(1/theta_e)) exp(lsynch(nu/nu_s)),
        !   nu_s = (2/9) nu_c theta_e^2 sin(theta_B).
        ! K2(1/theta_e) is computed (besselk2 below, to ~1e-14 relative), not
        ! approximated by 2 theta_e^2 - 1/2,
        ! which is exact at the default (fit only above the table's largest
        ! theta_e, 4.6e4) but 7.7% low at theta_e = 1 and zero (100% low) at 0.5
        ! (review item C8, fixed 2026-10-15; history.tex, Stage 46).  Where K2 underflows (theta_e < ~1/700) the
        ! point gets j = alpha = 0.
        subroutine synch_emissivity(nu, n_e, B, theta_B, theta_e, j, alpha)
        real(kind=8), dimension(:), intent(in) :: nu, n_e, B, theta_B, theta_e
        real(kind=8), dimension(:), intent(out) :: j, alpha
        real(kind=8) :: nuc, nus, k2, bnu
        integer :: i
        do i = 1, size(nu)
           j(i) = 0d0; alpha(i) = 0d0
           if (.not. (n_e(i) > 0d0 .and. B(i) > 0d0 .and. theta_e(i) > 0d0)) cycle
           nuc = e*B(i)/(2d0*pi*m*c)
           nus = (2d0/9d0)*nuc*theta_e(i)**2*max(abs(sin(theta_B(i))), 1d-6)
           k2 = besselk2(1d0/theta_e(i))
           if (.not. (k2 > 0d0 .and. k2 < huge(1d0))) cycle
           j(i) = sqrt(2d0)*pi*e*e*n_e(i)*nus/(3d0*c*k2) &
                * exp(max(lsynch(min(nu(i)/nus, 1d30)), -700d0))
           bnu = planck(theta_e(i)*m*c*c/k, nu(i))
           if (bnu > 0d0) alpha(i) = j(i)/bnu
        enddo
        end subroutine synch_emissivity

        ! K_2(x) = int_0^inf exp(-x cosh t) cosh(2t) dt by the trapezoid rule in
        ! t, which converges exponentially for this integrand: step
        ! min(0.1, 0.5/sqrt(x)) (the peak at t = 0 has width ~1/sqrt(x)), summed
        ! until the terms fall e^-40 below the first.  Agrees with
        ! scipy.special.kn(2, x) to < 2e-14 for x = 2e-5 to 500; underflows to 0
        ! above x ~ 745.  (Kept here rather than bessel.f90's besselk so that
        ! jnu_table.o needs no new Makefile dependency.)
        real(kind=8) function besselk2(x)
        real(kind=8), intent(in) :: x
        real(kind=8) :: hh, t, a
        hh = min(0.1d0, 0.5d0/sqrt(x))
        besselk2 = 0.5d0*exp(-x)
        t = 0d0
        do
           t = t + hh
           a = x*cosh(t)
           besselk2 = besselk2 + exp(-a)*cosh(2d0*t)
           if (a - 2d0*t > x + 40d0) exit
        enddo
        besselk2 = besselk2*hh
        end function besselk2

        subroutine hybrid_emissivity(nu, n_e, B, theta_B, theta_e, j, alpha)
        real(kind=8), dimension(:), intent(in) :: nu, n_e, B, theta_B, theta_e
        real(kind=8), dimension(:), intent(out) :: j, alpha
        real(kind=8), dimension(size(nu)) :: ntheta
        real(kind=8), dimension(:), allocatable :: jc, ac
        logical, dimension(size(nu)) :: cyc, syn, tab
        integer :: nc
        j = 0d0; alpha = 0d0
        where (B > 0d0)
           ntheta = nu/(e*B/(2d0*pi*m*c))*theta_e
        elsewhere
           ntheta = huge(1d0)
        endwhere
        ! Points with n_e <= 0 (electron-free or sigma-cut cells) or theta_e <= 0
        ! (the writer floors a negative theta_e to 0) go to no branch and keep
        ! j = alpha = 0.  They used to go to the table, which clamped
        ! theta_e = 0 to its smallest theta_e and gave n_e > 0 points j > 0
        ! with alpha = 0 (emission without absorption, log F <= -60 there;
        ! docs/muchado.md, L8-F5); the cyclotron sum's harmonic cap
        ! 2/(e^2 theta_e) is infinite at theta_e = 0.  A NaN n_e or theta_e
        ! still reaches the table and gives a NaN j, which the driver reports
        cyc = theta_e < th_hybrid .and. ntheta <= ntheta_max .and. n_e > 0d0 .and. theta_e > 0d0
        if (th_synch > 0d0) then
           syn = (.not. cyc) .and. theta_e >= th_synch .and. n_e > 0d0
        else
           syn = (.not. cyc) .and. theta_e > te_axis(nte) .and. n_e > 0d0
        endif
        tab = .not. (cyc .or. syn .or. n_e <= 0d0 .or. theta_e <= 0d0)
        ! each branch is evaluated only where it is used (the table also warns
        ! about clamping, which gas sent to the synchrotron fit must not trigger)
        nc = count(tab)
        if (nc > 0) then
           allocate(jc(nc), ac(nc))
           call jnu_table_emissivity(pack(nu, tab), pack(n_e, tab), pack(B, tab), &
                pack(theta_B, tab), pack(theta_e, tab), jc, ac)
           j = unpack(jc, tab, j); alpha = unpack(ac, tab, alpha)
           deallocate(jc, ac)
        endif
        nc = count(syn)
        if (nc > 0) then
           allocate(jc(nc), ac(nc))
           call synch_emissivity(pack(nu, syn), pack(n_e, syn), pack(B, syn), &
                pack(theta_B, syn), pack(theta_e, syn), jc, ac)
           j = unpack(jc, syn, j); alpha = unpack(ac, syn, alpha)
           deallocate(jc, ac)
        endif
        ! evaluate the cyclotron formula only where it is used: its harmonic
        ! sum runs to ~nu/nu_c, which at high frequency in cold gas is ~10^6
        ! terms per point that the dispatch would then discard
        nc = count(cyc)
        if (nc > 0) then
           allocate(jc(nc), ac(nc))
           call cyclotron_emissivity(pack(nu, cyc), pack(n_e, cyc), pack(B, cyc), &
                pack(theta_B, cyc), pack(theta_e, cyc), jc, ac)
           j = unpack(jc, cyc, j); alpha = unpack(ac, cyc, alpha)
           deallocate(jc, ac)
        endif
        end subroutine hybrid_emissivity

      end module jnu_table
