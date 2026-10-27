! cyclotron_emis.f90
!
! Non-relativistic (cyclotron) thermal emissivity, valid for theta_e << 1.
! This is the exact low-temperature limit of the Maxwell-Juettner synchrotron
! emissivity, from the low-temperature expansion of the Schott formula
! (docs/CS_Derivations_Sum26.tex, "Non-relativistic (cyclotron) limit"; ported from the
! jnu_calculator/symph2_jnu_plots.ipynb `jnu_cyclotron` function).
!
!       j_nu(nu) = sum_{n>=1} j_n^(int) * phi_n(nu)
!
! with the harmonic weight
!       j_n^(int) = (2 pi n_e e^2 nu_c^2 / c)
!                   * (n+1) n^(2n+1) / (2n+1)! * 2^n n! * theta_e^n
! and a Gaussian line phi_n with the mean and variance of the thermal line:
!       nu_n       = n nu_c (1 - k theta_e),           k = n + 3/2,
!       Delta_nu_n = n nu_c sqrt(theta_e cos^2(theta_B) + k theta_e^2),
!       phi_n(nu)  = exp[-(nu - nu_n)^2 / (2 Delta_nu_n^2)]
!                    / (sqrt(2 pi) Delta_nu_n) .
! The cos^2 term is the Doppler width.  The others are the relativistic mass
! shift, nu = n nu_c / gamma: the electrons that emit harmonic n have kinetic
! energies gamma - 1 distributed as epsilon^(n+1/2) exp(-epsilon/theta_e), a
! Gamma distribution of mean k theta_e and variance k theta_e^2, so the line
! sits below n nu_c and has a width of its own, which dominates within
! |cos theta_B| <~ sqrt(theta_e) of 90 degrees.  With the Doppler width alone the
! line there collapsed to a spike at n nu_c (with no overlap with the exact line
! at 90 degrees); the moment-matched line overlaps it by 83-92% at every angle
! (theta_e = 1e-3 to 5e-3, harmonics 1-3; history.tex, Stage 81).
!
! Direction: j_n^(int) is the angle-averaged weight (total power / 4 pi).  The
! emission toward theta_B carries the single-particle radiation pattern of
! harmonic n, which is leading order (not a relativistic correction):
!       D_n(theta_B) = sin^(2n-2)(theta_B) (1 + cos^2 theta_B) / <same>,
!       <same> = 4^n n (n+1) ((n-1)!)^2 / (2n+1)!   (average over cos theta_B),
! so that D_1 = 3/2 along the field and 3/4 across it, and every n >= 2
! vanishes along the field (CS_Derivations_Sum26.tex, "Doppler line profile"
! and "Convergence and validity").  j_n = j_n^(int) D_n phi_n.  theta_B also
! sets the Doppler line width.  In the product j_n^(int) D_n the factorials
! (2n+1)! cancel:
!       j_n^(int) D_n = (2 pi n_e e^2 nu_c^2 / c) n^(2n+2) theta_e^n / (2^n n!)
!                       * sin^(2n-2)(theta_B) (1 + cos^2 theta_B),
! so the sum needs log(n!) only, one log_gamma per point and a recurrence in n
! (until 2026-10-27, five log_gamma calls per harmonic; the same values to
! ~1e-14).
!
! Absorption is added from Kirchhoff's law (valid for a thermal, isotropic
! Maxwellian distribution -- which is our case):
!       alpha_nu = j_nu / B_nu(T_e),   T_e = theta_e m_e c^2 / k_B.
! The Rayleigh-Jeans limit of B_nu is used, which overestimates B_nu (and so
! underestimates alpha_nu) by about x/2, x = h nu/(k_B T_e) = 4.8e-11 nu[Hz]/T_e[K]:
! ~1e-9 at 1.6 MHz for theta_e >~ 1e-5, but at 230 GHz ~2e-5 at the largest
! theta_e this branch takes there (n theta_e <= 0.1, for nu_B ~ 1e8 Hz), and
! ~1% only for theta_e ~ 1e-7.
!
! Units: cgs.  j_nu [erg s^-1 cm^-3 Hz^-1 sr^-1], alpha_nu [cm^-1].

      module cyclotron_emis

      use phys_constants, only: c, e, m, k, pi

      implicit none

      ! small |cos(theta_B)| floor in the Doppler width (the mass-shift width
      ! keeps the line finite at 90 degrees in any case)
      real(kind=8), parameter :: COSTH_FLOOR = 1d-4

      interface cyclotron_emissivity
        module procedure cyclotron_emissivity
      end interface

      contains

        ! log(n!) = log_gamma(n+1)
        real(kind=8) function log_factorial(n)
        integer, intent(in) :: n
        log_factorial = log_gamma(dble(n) + 1d0)
        end function log_factorial

        subroutine cyclotron_emissivity(nu, n_e, B, theta_B, theta_e, j, alpha)
        real(kind=8), dimension(:), intent(in) :: nu, n_e, B, theta_B, theta_e
        real(kind=8), dimension(:), intent(out) :: j, alpha
        integer :: i, n, nmax, nlo
        real(kind=8) :: nuc, costh, sig0, lnW, jn, dnu, g, TeK, bnu, sin2, lnD
        real(kind=8) :: kk, nun, sigm, lfn, lnc
        do i = 1, size(nu)
           nuc = e*B(i)/(2d0*pi*m*c)
           if (nuc <= 0d0 .or. n_e(i) <= 0d0 .or. theta_e(i) <= 0d0 .or. &
               nu(i) <= 0d0) then
              j(i) = 0d0; alpha(i) = 0d0
              cycle
           endif
           costh = cos(theta_B(i))
           sin2 = max(0d0, 1d0 - costh*costh)
           sig0 = sqrt(theta_e(i))*max(abs(costh), COSTH_FLOOR)

           ! only harmonics near nu/nu_c contribute (Gaussian ~ 0 far away)
           nmax = nint(nu(i)/nuc) + 20
           ! Harmonic cap: the low-temperature weight below behaves like
           ! (n e theta_e/2)^n, whose successive ratio -> (e^2/2) n theta_e, so
           ! it stops decreasing at n* = 2/(e^2 theta_e) and grows without bound
           ! (overflowing to Inf) beyond it.  The true emission there is
           ! negligible, so stop the sum at n* (history.tex, "The cyclotron branch
           ! breaks down at high harmonics").  hybrid_emissivity only sends
           ! points with n theta_e <= 0.1 here, well below n*; the cap is a
           ! safety net for any other caller.
           nmax = min(nmax, max(1, int(2d0/(exp(2d0)*theta_e(i)))))
           if (nmax < 1) nmax = 1
           ! harmonics below nu/nu_c/(1 + 40 sigm) sit more than 40 line widths
           ! below nu (the mass shift only lowers them further): their Gaussian is
           ! exp(-800), exactly 0 in double precision, so skipping them changes
           ! nothing and avoids summing ~nu/nu_c terms.  sigm: the relative width
           ! of the widest line summed (n = nmax)
           sigm = sqrt(sig0*sig0 + (dble(nmax)+1.5d0)*theta_e(i)*theta_e(i))
           nlo = max(1, int(min(nu(i)/nuc/(1d0 + 40d0*sigm), dble(nmax))))

           j(i) = 0d0
           ! log(n!), carried from n = nlo by lfn = lfn + log(n)
           lfn = log_factorial(nlo - 1)
           lnc = log(1d0 + costh*costh)
           do n = nlo, nmax
              lfn = lfn + log(dble(n))
              ! directional factor D_n (see header); along the field only
              ! the fundamental survives
              if (n > 1 .and. sin2 <= 0d0) cycle
              ! log(j_n^(int) D_n) without the prefactor (see header)
              lnW = (2d0*n+2d0)*log(dble(n)) + dble(n)*(log(theta_e(i)) - log(2d0)) - lfn
              lnD = lnc
              if (n > 1) lnD = lnD + dble(n-1)*log(sin2)
              jn = 2d0*pi*n_e(i)*e*e*nuc*nuc/c * exp(lnW + lnD)
              kk = dble(n) + 1.5d0
              nun = dble(n)*nuc*(1d0 - kk*theta_e(i))
              dnu = dble(n)*nuc*sqrt(sig0*sig0 + kk*theta_e(i)*theta_e(i))
              g = exp(-0.5d0*((nu(i)-nun)/dnu)**2) &
                  / (sqrt(2d0*pi)*dnu)
              j(i) = j(i) + jn*g
           enddo

           ! Kirchhoff absorption (Rayleigh-Jeans Planck, h nu << k T_e)
           TeK = theta_e(i)*m*c*c/k
           bnu = 2d0*nu(i)*nu(i)*k*TeK/c/c
           if (bnu > 0d0) then
              alpha(i) = j(i)/bnu
           else
              alpha(i) = 0d0
           endif
        enddo
        end subroutine cyclotron_emissivity

      end module cyclotron_emis
