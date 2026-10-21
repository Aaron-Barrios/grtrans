
       module grtrans

       use omp_lib

       use fluid_model, only: initialize_fluid_model, del_fluid_model, &
        get_fluid_vars, fluid, source_params, convert_fluid_vars, &
        initialize_source_params,del_source_params, APERTURE
       use class_rad_trans, only: initialize_rad_trans, del_rad_trans, &
        rad_trans 
       use geodesics, only: initialize_geodesic, del_geodesic, &
         geokerr_args, geo
       use emissivity
       use ray_trace, only: save_raytrace_camera_pixel, ray_set
!       use interpolate, only: get_weight, locate
       use phys_constants, only: GC => G, EC => E, CC => C, C2, MSUN, &
            ME => m, KB => k, PI => pi
       use kerr, only: comoving_ortho, calc_polar_psi, calc_kb_ang, calcg, lnrf_frame, &
            comoving_ortho_debug, kerr_metric, calc_nullp
       use class_four_vector, only: four_vector, assign_metric
       use grtrans_inputs, only: refine_on, refine_rtol, refine_atol_frac, &
            refine_tau_cut, refine_max_depth, refine_max_points
       use chandra_tab24, only: interp_chandra_tab24
       use radtrans_integrate, only: integrate, del_radtrans_integrate_data, &
            init_radtrans_integrate_data, intensity
       use fluid_model_aperture, only: naperture, select_aperture_binary

       implicit none

! grtrans object is made up of fluid, geodesic, emissivity, and rad trans objects

!        integer, parameter :: IS_LINEAR_STOKES=1
       integer :: EXTRA_QUANTS=0, WRITE_GEO=0
!        integer, dimension(3) :: stats
!        real(kind=8), parameter :: MAX_TAU = 10d0
       real(kind=8) :: LBH, MBH
!        real(kind=8) :: ortol=1.d-6, oatol=1.d-8, hmax
!        real(kind=8), dimension(:), allocatable :: ss
       type (fluid) :: f
       type (rad_trans) :: r
       type (geo) :: g
       type (emis) :: e
!         type (emis_params) :: eparams
!         type (source_params) :: sparams
!$omp threadprivate(f,r,g,e,EXTRA_QUANTS,WRITE_GEO,LBH,MBH)
! Adaptive ray refinement bookkeeping (refine_ray), indexed by ray number i and
! shared between threads (each ray belongs to one thread): points added, summed
! over frequencies and cameras, and a flag with bit 1 = an interval was still
! unconverged after refine_max_depth bisections, bit 2 = refine_max_points ran
! out.  Allocated by init_refine_stats only when refine_on.
       integer, dimension(:), allocatable :: refine_added, refine_flag
! Rays whose emissivity had a NaN somewhere (shared, incremented atomically by
! grtrans_driver; reset and reported by grtrans_main, pgrtrans.f90).  Such a
! ray's fac is NaN, so its pixel and extras are stored as 0, except extra
! slot 2 (docs/muchado.md, L8-G4).
       integer :: nan_rays = 0

!       interface grtrans_driver_jac_form
!         module procedure grtrans_driver_jac_form
!       end interface

!       interface grtrans_driver_rhs_form
!         module procedure grtrans_driver_rhs_form
!       end interface

       interface grtrans_driver
         module procedure grtrans_driver
       end interface

       contains

         subroutine grtrans_driver(gargs,gunit,c,i,l, &
          iname,ename,fname,sparams,eparams,nfreq,nparams,freqs,nup,extra,debug)
! Driver routine for GRTrans. JAD 2/27/2011
         implicit none
         type (geokerr_args), intent(in) :: gargs
         type (ray_set), intent(inout), dimension(:) :: c
         type (emis_params), intent(in) :: eparams
         type (source_params), dimension(:), intent(in) :: sparams
         type (source_params) :: sp
         type (emis_params) :: ep
!         real(kind=8), intent(in) :: mbh1,mdot
         real(kind=8), intent(in), dimension(:) :: freqs
         integer, intent(in) :: gunit,i,l,nup, nfreq, nparams, extra, debug
         character(len=100), intent(in) :: iname,fname,ename
         real(kind=8) :: fac
         real(kind=8), dimension(:), allocatable :: s2xi,c2xi,s2psi,c2psi, &
          rshift,ang,nu,cosne,tau,tau_temp,dummy,vrl,vtl,vpl,cosne2, &
         aat,aar,aath,aaph,kht,khr,khth,khph,bht,bhr,bhth,bhph,dlp,subw
         real(kind=8), dimension(:,:), allocatable :: aahat
         real(kind=8), dimension(:,:), allocatable :: tau_arr
         integer :: k,npow,status,j,nf,nitems,m,kk,taudex,ii,nvals,nnan
         ! the one-pass extra outputs (slots 7-14, 20, 21)
         integer :: np_
         real(kind=8) :: acc(0:10), wgt, b2one, q20, q21
         integer, dimension(:), allocatable :: inds
         ! several APERTURE binaries in one run (fdfile = a '.lst' file,
         ! fluid_model_aperture): binary m fills the cameras of parameter m
         logical :: multibin
         npow=3; j=0; status=0; k=0
         ep=eparams
         if(extra==1) EXTRA_QUANTS=1
! Need to make this a decision of load geodesics from file or generate new ones, and pass 
!gunit as last argument if loading.
!         write(6,*) 'grtrans_driver: ',i,debug
         call initialize_geodesic(g,gargs,i,status)
         if(debug==1) then
            WRITE_GEO=1
         else
            WRITE_GEO=0
         endif
         if (WRITE_GEO==1) then
            nitems=18
! Dump geodesic information for debugging
            open(unit=9,file='geodebug.out',form='formatted')
              ! header with number of points, number of items, mdot, mbh, frequency
            write(9,*) g%npts, nitems, g%gk%alpha(1), g%gk%beta(1),g%gk%su,g%gk%sm, &
                          sparams%mdot,sparams%mbh, fac, g%gk%a, g%gk%mu0
! geodesic properties
            write(9,*) g%lambda
            write(9,*) g%x%data(1)
            write(9,*) g%x%data(2)
            write(9,*) g%x%data(3)
            write(9,*) g%x%data(4)
            write(9,*) g%k%data(1)
            write(9,*) g%k%data(2)
            write(9,*) g%k%data(3)
            write(9,*) g%k%data(4)
            write(9,*) g%tprarr
            write(9,*) g%tpmarr
         endif
         if(status==1.and.g%npts.gt.0) then
!            hmax=(g%lambda(1)-g%lambda(2))/3.
!            if(hmax.gt.10.) hmax=10.
!            hmax=0.1
!            write(6,*) 'hmax: ',hmax
!         write(6,*) 'after init', g%npts,c(1)%nvals
            allocate(nu(g%npts))
            allocate(s2xi(g%npts)); allocate(c2xi(g%npts))
            s2xi=0d0; c2xi=0d0
            allocate(rshift(g%npts)); allocate(ang(g%npts))
            allocate(cosne(g%npts)); allocate(tau(g%npts))
            call initialize_fluid_model(f,fname,gargs%a,g%npts)
            call select_emissivity_values(e,ename)
!         write(6,*) 'evals: ',e%nk,e%neq
            multibin = (f%model==APERTURE .and. naperture>1)
            if (multibin) call select_aperture_binary(1)
            ! Get fluid variables:
            call get_fluid_vars(g%x,g%k,g%gk%a,f)
            !write(6,*) 'fluid model: '!,f%rho
            call comoving_ortho(g%x%data(2),g%x%data(3),g%gk%a, &
                 g%gk%alpha(1),g%gk%beta(1),g%gk%mu0,f%u,f%b,g%k,s2xi,c2xi, &
                 ang,rshift,cosne)
            !write(6,*) 'comoving ortho grtrans driver: '!,s2xi,c2xi,'\n'
            ! an emissivity with one Stokes component cannot fill nvals = 4
            ! (e.g. JNUTABLE: e%j would be read as four columns in integrate)
            if (e%neq==1 .and. c(1)%nvals==4) then
               write(6,*) 'ERROR: emissivity ', trim(ename), ' gives Stokes I only;', &
                    ' set nvals=1 in inputs.in'
               stop 1
            endif
            call initialize_rad_trans(r,iname,g%npts,c(1)%nvals,extra)
            call init_radtrans_integrate_data(r%iflag,r%neq,g%npts,r%npts)
            do m=1,nparams
! Several APERTURE binaries: everything after the geodesic depends on the
! binary, so the fluid at the samples and the comoving frame are redone for
! binary m (m=1 was done above) along the same geodesic g; the rest of this
! loop is then exactly a single-binary trace of binary m.
               if (multibin .and. m>1) then
                  call select_aperture_binary(m)
                  call get_fluid_vars(g%x,g%k,g%gk%a,f)
                  s2xi=0d0; c2xi=0d0
                  call comoving_ortho(g%x%data(2),g%x%data(3),g%gk%a, &
                       g%gk%alpha(1),g%gk%beta(1),g%gk%mu0,f%u,f%b,g%k,s2xi,c2xi, &
                       ang,rshift,cosne)
               endif
               sp=sparams(m)
               !write(6,*) 'grtrans_driver sp: ',sp%nfac,sp%gminval
               call initialize_source_params(sp,g%npts)
               call initialize_emis_params(ep,g%npts)
               call initialize_emissivity(e,g%npts,f%nfreq,rshift,ang,cosne,f%nrelbin,f%bingammamin,f%bingammamax)!,emisargs)
               !write(6,*) 'f%nrelbin',f%nrelbin
               !               write(6,*) 'm: ',m,nparams,sparams(m)%mdot,size(sparams)
!               allocate(sparams(m)%gmin(npts))
               MBH=sp%mbh
               LBH=GC*MBH*MSUN/C2
!               write(6,*) 'grtrans_driver convert: ',sp%nfac,sp%gminval
               call convert_model(sp)
!               write(6,*) 'gmin mu: ',sp%gmin
!               write(6,*) 'gmin mu: ',sp%mu
               ep%gmin=sp%gmin
               ep%mu=sp%mu
              !write(6,*) 'after convert'
!              write(6,*) 'emis model',FREQS
!              write(6,*) 'emis model'
               call emis_model(e,ep)
! Compute emissivity:
!               write(6,*) 'emis',NFREQ
               do k=1,NFREQ
                  nu=FREQS(k)/rshift
                  call calc_emissivity(nu,e,ep)
                  nnan=count(isnan(e%j(:,1)))
                  if(nnan > 0) then
                     write(6,*) 'NaN in emissivity at i: ',i
!$omp atomic update
                     nan_rays=nan_rays+1
                  endif

!                  if(any(isnan(e%K))) then
                     !write(6,*) 'NaN in opacity at i: ',i
!                     write(6,*) 'NaN in emissivity where: ',maxloc(isnan(e%j))
!                     write(6,*) 'NaN in emissivity where: ',maxloc(isnan(e%K))
                     !write(6,*) 'NaN in emissivity fac LBH: ',fac,LBH
                     !write(6,*) 'fluid properties: ',e%ncgs,e%tcgs,e%bcgs
!                  endif
   
!               write(6,*) 'emis i: ',i
!               write(6,*) 'emis j: ',minval(e%j(:,1))
!!               write(6,*) 'emis n: ',minval(e%ncgsnth)
!               write(6,*) 'emis b: ',minval(e%bcgs)
!               write(6,*) 'emis ang: ',e%incang
! don't integrate if there's no emission:
!                  write(6,*) 'ej: ',sum(e%j(:,1))
!               write(6,*) 'ortho: ',s2xi
!                  write(6,*) 'tiny: ',epsilon(sum(e%j(:,1))),g%npts*epsilon(sum(e%j(:,1)))
!                  if (sum(e%j(:,1)).gt.(g%npts*tiny(sum(e%j(:,1))))) then
!                  fac=sum(e%j(:,1))/g%npts
                  if (g%npts.gt.1) then
                     fac=sum(e%j(:,1)*(g%lambda(1)-g%lambda(1:g%npts)))
                  else
                     fac=sum(e%j(:,1)/g%npts)
                  endif
                  !write(6,*) 'i fac: ',i,fac,1d0/(1d0/fac)
!                  if(g%npts.gt.1) then
!                     fac=tsum(g%lambda(1)-g%lambda(1:g%npts),e%j(:,1))
!                  else
!                     fac=1d0
!                  endif
!                  if (1d0/(1d0/(fac/1d10)).gt.0d0) then
                  if (fac/1d16.gt.0d0) then
!                     write(6,*) 'emis i: ',i
!                     write(6,*) 'emis j: ',maxval(e%j(:,1))
!                     write(6,*) 'emis n: ',minval(e%ncgs)
!                     write(6,*) 'emis b: ',minval(e%bcgs)
!                     write(6,*) 'emis te: ',i,e%tcgs(maxloc(e%ncgs))
                     if(e%neq.eq.4) call rotate_emis(e,s2xi,c2xi)
! Call integration routine:
                     !write(6,*) 'integrate'
!                     if(any(isnan(e%j))) then
                        !write(6,*) 'NaN in emissivity after rotate at i: ',i
                        !write(6,*) 'fluid properties: ',e%ncgs,e%tcgs,e%bcgs
!                     endif
!                     if(any(isnan(e%K))) then
                        !write(6,*) 'NaN in opacity after rotate at i: ',i
                        !write(6,*) 'fluid properties: ',e%ncgs,e%tcgs,e%bcgs
!                     endif

                     if(g%npts.ne.1) then
                        call invariant_emis(e,rshift)
!                        fac=sum(e%j(:,1))/g%npts
! scale emission close to 1 so that LSODA integrates accurately and put anu in cgs units
                        e%j=e%j/fac; e%K=e%K*LBH
!                        write(6,*) 'i: ',i
                        call calc_opt_depth(g%lambda(1:g%npts),e,tau,1)
                        tau = -tau
!                        call init_radtrans_integrate_data(r%iflag,r%neq,g%npts,r%npts)
                        !write(6,*) 'maxtau: ',maxval(e%K(:,1)),maxval(abs(sqrt(e%K(:,5)**2.+e%K(:,6)**2.))),maxval(abs(tau))
                        call integrate(g%lambda(1:g%npts),e%j,e%K,tau,r%npts)
                        r%I(1:r%neq,1:r%npts) = intensity(1:r%neq,1:r%npts)
!                        call del_radtrans_integrate_data()
!                        write(6,*) 'integrate: ',g%npts,r%neq,r%npts,r%I(1,r%npts)
! added factor of LBH for cgs intensity units ergs / cm^2 / s / Hz / ster
                        r%I=r%I*fac*LBH
! Adaptive ray refinement (off by default; &refine in inputs.in): resample the
! ray where it is under-resolved and replace the explin intensity.
                        if (refine_on .and. r%iflag==4 .and. r%neq==1) &
                             call refine_ray(i,freqs(k),sparams(m),eparams,fac, &
                             fname,ename,rshift)
!                        write(6,*) 'grtrans driver extra quants: ',EXTRA_QUANTS,r%npts
                        if (EXTRA_QUANTS==1) then
!                           if(r%npts.gt.1) then
                           allocate(tau_arr(g%npts,6))
                           allocate(tau_temp(g%npts))
                           allocate(dummy(g%npts)); allocate(inds(6)); allocate(dlp(r%npts-1))
                           inds=(/1,2,3,4,5,7/)
!                           write(6,*) 'tau_arr', size(g%lambda), size(tau_temp)
! Slots 2-6 are the optical depths of the polarized coefficients (K columns
! 2,3,4,5,7).  With nvals=1, e%K has a single column, so those optical depths
! are set to zero rather than read past the end of e%K.
                           do ii=1,6
                              if (inds(ii) > size(e%K,2)) then
                                 tau_arr(:,ii)=0d0
                                 cycle
                              endif
                              call calc_opt_depth(g%lambda(1)-g%lambda(1:g%npts),e,tau_temp,inds(ii))
!                              write(6,*) 'extra after opt_depth ',ii, &
!                                   size(tau_arr,1), size(tau_arr,2), size(tau_temp)
                              tau_arr(:,ii)=tau_temp
                           enddo
                           ! now calculate photosphere location, values there:
                           ! tau is now the running optical depth from the camera
                           ! (calc_opt_depth); the ray is optically thin, and the
                           ! whole ray counts, when its total tau is below 1.
                           ! (Was minval(tau), which is always 0 for a running
                           ! integral.)
                           taudex=minloc(abs(tau-1d0),1)
                           if(maxval(tau).lt.1d0) taudex=g%npts
                           if(taudex.eq.1) taudex=2
!                           write(6,*) 'taudex: ',taudex,tau(1),tau(g%npts),tau(taudex)
!                           write(6,*) 'ej: ',e%j(:,1), fac, r%I, tau
                           r%tau(1:6)=tau_arr(taudex,:)
                           ! slot 1 = TOTAL optical depth of the ray in Stokes I
                           ! (turnover_frequency.py reads it); slots 2-6 keep the
                           ! stock values at the photosphere taudex.
                           r%tau(1)=tau_arr(g%npts,1)
!                           write(6,*) 'rtau',size(r%tau(1:6)),size(tau_arr(taudex,:))
                           ! now do emissivity-weighted ray averages of other quantities
                           dummy=e%j(:,1)*exp(-tau(1:g%npts))
                           if(r%neq==4) then
                              dlp=sqrt(r%I(1,r%npts:2:-1)**2.+r%I(2,r%npts:2:-1)**2.)-&
                                sqrt(r%I(1,r%npts-1:1:-1)**2.+r%I(2,r%npts-1:1:-1)**2.)
                           else
                              dlp(:)=0.
                           endif
!                           write(6,*) 'dlp: ',minval(dlp),maxval(dlp)
!                           write(6,*) 'rtau', size(r%tau),size(r%tau(7:7))
!                           write(6,*) 'tsum: ',tsum(g%lambda(1:taudex),dummy(1:taudex))
!                           write(6,*) 'tsum: ',tsum(g%lambda(1:taudex),dummy(1:taudex)*g%x(1:taudex)%data(2))
! Slots 7-14 and 20-21 in ONE pass over the ray (2026-10-17, docs/history.tex,
! Stage 62).  Each slot is <q> = int(w q dlambda) / int(w dlambda) with the
! weight w = j exp(-tau); the trapezoid rule makes either integral
! sum_i c_i y_i with c_i half the sum of the two neighboring steps of lambda,
! so one loop accumulates all the numerators and the denominator.  Until then
! each slot called tsum (math.f90) on its own temporary arrays, about a dozen
! passes per ray, which cost several times the transfer itself.  The ratio is
! the same quantity (the sign of the decreasing lambda cancels; only the
! summation order differs), and the per-point expressions are those of the
! old code, single-precision fluid fields included.  A ray whose weight
! integrates to 0 (e.g. exp(-tau) underflowing wherever j > 0) still gets 0 in
! all of them rather than 0/0 (docs/muchado.md, L8-G8).
! Slots 20-21 (APERTURE): the PIPELINE's beta_p (f%rho2) and sigma (f%Be),
! carried in the binary (relabel.beta_p; b^2/(2 rho_0) of the pair fluid);
! grtrans's own 4 n_e theta_e/b^2 and b_e^2/(4 n_e_rest) are not those once
! the temperature is relabeled.  For other models the stock beta and
! sigma, masked as before (b^2 > 0, rho > 0, floored at 0) because bilinear
! interpolation can undershoot below zero at the torus edge.
                           allocate(subw(g%npts))
                           np_=g%npts
                           subw=0d0          ! here: the trapezoid weights c_i
                           if (np_ > 1) then
                              subw(1:np_-1)=subw(1:np_-1)+0.5d0*(g%lambda(2:np_)-g%lambda(1:np_-1))
                              subw(2:np_)=subw(2:np_)+0.5d0*(g%lambda(2:np_)-g%lambda(1:np_-1))
                           endif
                           acc=0d0
                           do ii=1,np_
                              wgt=subw(ii)*dummy(ii)
                              acc(0)=acc(0)+wgt
                              acc(1)=acc(1)+wgt*g%x(ii)%data(2)
                              acc(2)=acc(2)+wgt*g%x(ii)%data(3)
                              acc(3)=acc(3)+wgt*g%x(ii)%data(4)
                              acc(4)=acc(4)+wgt*e%ncgs(ii)
                              acc(5)=acc(5)+wgt*e%tcgs(ii)
                              acc(6)=acc(6)+wgt*e%bcgs(ii)
! slot 13 is \beta (the stock 2p/b^2)
                              acc(7)=acc(7)+subw(ii)*(dummy(ii)*f%p(ii)*2./f%bmag(ii)**2.)
! slot 14: fraction of emission produced above/below the midplane (e.g. counter-jet or forward jet images)
                              acc(8)=acc(8)+wgt*sign(1d0,cos(g%x(ii)%data(3)))
                              b2one=max(f%bmag(ii)**2., 0d0)
                              if (f%model==APERTURE) then
                                 q20=max(dble(f%rho2(ii)),0d0)
                                 q21=max(dble(f%Be(ii)),0d0)
                              else
                                 q20=0d0
                                 if (b2one > 0d0) q20=4d0*max(f%rho(ii),0d0)*max(f%p(ii),0d0)/b2one
                                 q21=0d0
                                 if (f%rho(ii) > 0d0) q21=b2one/(4d0*f%rho(ii))
                              endif
                              acc(9)=acc(9)+wgt*q20
                              acc(10)=acc(10)+wgt*q21
                           enddo
                           if (abs(acc(0)) > 0d0) then
                              r%tau(7:14)=acc(1:8)/acc(0)
                              r%tau(20)=acc(9)/acc(0)
                              r%tau(21)=acc(10)/acc(0)
                           else
                              r%tau(7:14)=0d0
                              r%tau(20)=0d0
                              r%tau(21)=0d0
                           endif
                           deallocate(subw)
! slots 15-19: weighted r,th,\tau_FR,\tau_FC for LP
! Slots 15-19 are weighted by the change in linear polarization, which is zero
! (or has no elements) in total intensity: output zero instead of 0/0.
                           if (sum(dlp) /= 0d0) then
                              r%tau(15)=sum(dlp*g%x(1:r%npts-1)%data(2))/sum(dlp)
                              r%tau(16)=sum(dlp*g%x(1:r%npts-1)%data(3))/sum(dlp)
                              r%tau(17)=sum(dlp*tau_arr(1:r%npts-1,5))/sum(dlp)
                              r%tau(18)=sum(dlp*tau_arr(1:r%npts-1,6))/sum(dlp)
                              r%tau(19)=sum(dlp*sign(1d0,cos(g%x(1:r%npts-1)%data(3))))/sum(dlp)
                           else
                              r%tau(15:19)=0d0
                           endif
                           deallocate(tau_arr); deallocate(tau_temp); deallocate(dummy)
                           deallocate(inds); deallocate(dlp)
!                           else
                        endif
                     else
!             write(6,*) 'rshift :', size(rshift)
                        call invariant_emis(e,rshift,npow)
!             write(6,*) 'compute intensity'
                        call grtrans_compute_intensity()
                        if(EXTRA_QUANTS==1) then
!                           write(6,*) 'grtrans driver extra quants npts 1'
                           allocate(s2psi(g%npts)); s2psi=0d0
                           allocate(c2psi(g%npts)); c2psi=0d0
                           allocate(cosne2(g%npts)); cosne2=0d0
                           call transpol(rshift,s2psi,c2psi,cosne2)
                           r%tau(1)=s2xi(1)
                           r%tau(2)=s2psi(1)
                           r%tau(3)=cosne2(1)
                           r%tau(4)=ang(1)
                           r%tau(5)=cosne(1)
                           r%tau(6)=c2xi(1)
                           r%tau(7)=c2psi(1)
                           deallocate(s2psi); deallocate(cosne2); deallocate(c2psi)
                        endif
                     endif
                  else
                     r%I(:,:)=0d0; r%npts=1
                     if(EXTRA_QUANTS==1) then
                        r%tau(:)=0d0
                        ! a ray skipped because of a NaN in j (fac is NaN) is
                        ! flagged in slot 2 with its number of NaN samples, so it
                        ! can be told from a ray without emission (slot 2 is
                        ! otherwise 0 with one Stokes component; L8-G4)
                        if(nnan > 0 .and. e%neq==1) r%tau(2)=dble(nnan)
                     endif
!                     write(6,*) 'zero: ',status, size(r%I,1), size(r%I,2), r%I(:,r%npts)
                  endif
!           write(6,*) 'after intensity', i, r%npts, r%I(:,r%npts), sum(e%j(:,1))!, e%j(:,1), e%K(:,1)
!                  if(g%npts.gt.1.and.mod(i,500).eq.0.and.mod(k,nfreq).eq.0.and.mod(m,nparams).eq.0) &
!                       write(6,*) 'save ', i, r%I(:,r%npts),l,k,m
                  if(EXTRA_QUANTS==0) then
! this is normal situation
                     call save_raytrace_camera_pixel(c((l-1)*nfreq*nparams+(m-1)*nfreq+k), &
                       (/real(gargs%alpha(i)),real(gargs%beta(i))/),real(r%I(:,r%npts)),i)
                  else
! this is where we want to save additional info
!                     write(6,*) 'save: ',size(r%I(:,r%npts)),size(r%tau),&
!                          size((/real(r%I(:,r%npts)),real(r%tau)/))
!                     write(6,*) 'cam: ',size(c((l-1)*nfreq*nparams+(m-1)*nfreq+k)%pixvals)
!                     write(6,*) 'save cam extra: ',r%tau,s2xi(1),c2xi(1),status,g%npts,r%npts,EXTRA_QUANTS
!                     write(6,*) 'save cam extra: ',r%tau,real(r%tau)
!                     write(6,*) 'save cam extra: ',(/real(r%I(:,r%npts)),real(r%tau)/)
!                     write(6,*) 'save cam extra: ',sum(e%j(:,1))
                     call save_raytrace_camera_pixel(c((l-1)*nfreq*nparams+(m-1)*nfreq+k), &
                       (/real(gargs%alpha(i)),real(gargs%beta(i))/),(/real(r%I(:,r%npts)),&
                       real(r%tau)/),i)
!                     write(6,*) 'extra'
                  endif
                  !write(6,*) 'after save'
                  if (WRITE_GEO==1) then
! fluid properties
                     write(9,*) f%rho
                     write(9,*) f%p
                     write(9,*) f%bmag
                     write(9,*) f%u%data(1)
                     write(9,*) f%u%data(2)
                     write(9,*) f%u%data(3)
                     write(9,*) f%u%data(4)
                     write(9,*) f%b%data(1)
                     write(9,*) f%b%data(2)
                     write(9,*) f%b%data(3)
                     write(9,*) f%b%data(4)
! emission properties
                     write(9,*) e%ncgs
                     write(9,*) e%tcgs
                     write(9,*) e%bcgs
                     write(9,*) e%rshift
                     write(9,*) e%incang
                     write(9,*) e%j(:,1)
                     write(9,*) e%K(:,1)
! Faraday rotation, conversion coefs
                     write(9,*) e%j(:,2)
                     write(9,*) e%j(:,4)
                     write(9,*) e%K(:,2)
                     write(9,*) e%K(:,4)
                     write(9,*) e%K(:,5)
                     write(9,*) e%K(:,7)
! JAD 10/28/2019 removed /LBH/fac part here not sure ??
                     write(9,*) r%I(1,:)/fac/LBH
                     write(9,*) r%I(2,:)/fac/LBH
                     write(9,*) r%I(3,:)/fac/LBH
                     write(9,*) r%I(4,:)/fac/LBH
                     write(9,*) tau
                     write(9,*) e%j(:,3)
                     write(9,*) e%K(:,3)
                     write(9,*) e%K(:,6)
                     write(9,*) s2xi
                     write(9,*) c2xi
! calculate quantities with other methods:
                     write(9,*) calc_kb_ang(g%k,f%b,f%u,g%x%data(2),g%x%data(3),g%gk%a)
!                     call calc_polar_psi(g%x%data(2),cos(g%x%data(3)),g%gk%q2,g%gk%a, &
!                          g%gk%alpha,g%gk%beta,rshift,g%gk%mu0,g%k, &
!                          c2xi,s2xi,cosne)
!                     write(9,*) s2xi
!                     write(9,*) c2xi
                     allocate(vrl(g%npts)); allocate(vtl(g%npts)); allocate(vpl(g%npts))
                     call lnrf_frame(f%u%data(2)/f%u%data(1),f%u%data(3)/f%u%data(1), &
                          f%u%data(4)/f%u%data(1),g%x%data(2),g%gk%a,g%x%data(3), &
                          vrl,vtl,vpl)
                     write(9,*) calcg(1d0/g%x%data(2),cos(g%x%data(3)),g%gk%q2, &
                          g%gk%l,g%gk%a, &
                          g%tpmarr,g%tprarr,g%gk%su,g%gk%sm,vrl,vtl,vpl)
                     deallocate(vrl); deallocate(vtl); deallocate(vpl)
                     ! WRITE OTHER DEBUG QUANTITIES FROM COMOVING ORTHO
                     allocate(aahat(g%npts,3)); allocate(aat(g%npts))
                     allocate(bht(g%npts)); allocate(kht(g%npts))
                     allocate(bhr(g%npts)); allocate(khr(g%npts)); allocate(aar(g%npts))
                     allocate(bhth(g%npts)); allocate(khth(g%npts)); allocate(aath(g%npts))
                     allocate(bhph(g%npts)); allocate(khph(g%npts)); allocate(aaph(g%npts))
                     call comoving_ortho_debug(g%x%data(2),g%x%data(3),g%gk%a, &
                          g%gk%alpha(1),g%gk%beta(1),g%gk%mu0,f%u,f%b,g%k,s2xi,c2xi, &
                          ang,rshift,cosne,aahat,aat,aar,aath,aaph,kht,khr,khth,khph, &
                          bht,bhr,bhth,bhph)
                     write(9,*) s2xi
                     write(9,*) c2xi
                     write(9,*) aahat
                     write(9,*) aat
                     write(9,*) aar
                     write(9,*) aath
                     write(9,*) aaph
                     write(9,*) kht
                     write(9,*) khr
                     write(9,*) khth
                     write(9,*) khph
                     write(9,*) bht
                     write(9,*) bhr
                     write(9,*) bhth
                     write(9,*) bhph
                     write(9,*) fac
                     close(unit=9)
                     deallocate(aahat); deallocate(aat)
                     deallocate(bht); deallocate(kht)
                     deallocate(bhr); deallocate(khr); deallocate(aar)
                     deallocate(bhth); deallocate(khth); deallocate(aath)
                     deallocate(bhph); deallocate(khph); deallocate(aaph)
                  endif
!                  call del_rad_trans(r)
           !           write(6,*) 'after rad'
               enddo
!         write(6,*) 'before emis'
               call del_source_params(sp)
               call del_emis_params(ep)
               call del_emissivity(e)
        !         write(6,*) 'after emis'
!               deallocate(sparams(m)%gmin)
            enddo

            call del_rad_trans(r)
            call del_radtrans_integrate_data()
            call del_fluid_model(f)
            if (multibin) call select_aperture_binary(1)
            
     !         write(6,*) 'after fluid'
            deallocate(nu); deallocate(s2xi); deallocate(c2xi)
            deallocate(rshift); deallocate(ang); deallocate(cosne); deallocate(tau)
!         write(6,*) 'after xi'
!         if(g%npts==1.and.e%neq==4) then
!            deallocate(jout)
!            deallocate(jtemp)
!         endif
!         write(6,*) 'after polarpsi'
         else
! Geodesic failure. Save intensity as zero
!         write(6,*) 'status m1'
            do m=1,NPARAMS
               do k=1,NFREQ
                  call save_raytrace_camera_pixel(c((l-1)*nfreq*nparams+(m-1)*nfreq+k),(/real(gargs%alpha(i)), &
                       real(gargs%beta(i))/),(/(0.,j=1,c((l-1)*nfreq*nparams+(m-1)*nfreq+k)%nvals+c((l-1)* &
                       nfreq*nparams+(m-1)*nfreq+k)%nextra)/),i)
               enddo
            enddo
         endif
!         write(6,*) 'del geo: ',allocated(g%gk%muf)
         call del_geodesic(g)
  !        write(6,*) 'del geo: ',allocated(g%gk%muf)
         end subroutine grtrans_driver

         subroutine refine_ray(iray,fr,spm,epar,fac,fname,ename,rshift)
! Adaptive ray refinement (RT_background.tex, "How to resolve the rays").  The ray has
! been integrated with explin on the nn(3) geokerr points; e%j, e%K hold the
! invariant j/fac and K*LBH there (the optical depth is recomputed here,
! running_tau).  Rays with fewer than four points are left alone.  Intervals are bisected where the sampling is too coarse:
!   triggers (first pass, original intervals): Delta tau > 0.5; a cyclotron
!     harmonic crossing, floor(nu/(g nu_c)) different at the two ends, counted
!     only where harmonics are separate (n <= 2/(e^2 Theta_e), the cap of
!     cyclotron_emis.f90); or j or alpha changing by more than 10x;
!   only intervals with camera-side tau < refine_tau_cut, whose explin
!     contribution is at least refine_atol_frac of the ray total;
!   accuracy target: a bisected interval is accepted when its contribution
!     changes by less than refine_rtol, otherwise both halves are bisected
!     again, up to refine_max_depth levels and refine_max_points points.
! New points take their position from a cubic in the affine parameter through
! the neighboring geokerr points, their photon momentum from calc_nullp at
! that position (with the radial and polar signs of the nearer end point, so
! an interval containing a turning point is split at the right place), and
! then pass through the same fluid, tetrad and emissivity chain as the ray.
! The intensity is then the explin sum over the refined points; the extra
! quantities (EXTRA_QUANTS) stay on the original points.
         integer, intent(in) :: iray
         real(kind=8), intent(in) :: fr, fac
         type (source_params), intent(in) :: spm
         type (emis_params), intent(in) :: epar
         character(len=100), intent(in) :: fname, ename
         real(kind=8), dimension(:), intent(in) :: rshift
         real(kind=8), dimension(:), allocatable :: lam, jv, kv, tv, cold, &
              lam2, jv2, kv2, lmid, jmid, kmid
         real(kind=8), dimension(:), allocatable :: hv
         integer, dimension(:), allocatable :: bix, bix2, imid
         logical, dimension(:), allocatable :: act, act2
         real(kind=8) :: total, cnew, dl
         integer :: n0, n, cap, lev, ii, jj, nnew, nadd, flag
         n0=g%npts
         ! the cubic of refine_eval_points needs four points (i0 = n0-3 >= 1)
         if (n0 < 4) return
         cap=max(refine_max_points,n0)
         allocate(lam(cap),jv(cap),kv(cap),tv(cap),cold(cap),bix(cap),act(cap))
         allocate(hv(n0))
         n=n0
         lam(1:n)=g%lambda(1:n); jv(1:n)=e%j(:,1); kv(1:n)=abs(e%K(:,1))
         do ii=1,n
            bix(ii)=min(ii,n0-1)
         enddo
         hv=harmonic_index(fr/rshift(1:n0),e%bcgs,e%tcgs,f%model==APERTURE)
         call running_tau()
         act=.false.
         do ii=1,n-1
            act(ii)=(tv(ii+1)-tv(ii) > 0.5d0) .or. &
                 (hv(ii) /= hv(ii+1) .and. max(hv(ii),hv(ii+1)) >= 0d0) .or. &
                 jump(jv(ii),jv(ii+1)) .or. jump(kv(ii),kv(ii+1))
         enddo
         flag=0; nadd=0
         do lev=1,refine_max_depth
            if (lev > 1) call running_tau()
            total=explin_sum()
            do ii=1,n-1
               if (act(ii)) then
                  cold(ii)=term(ii)
                  if (tv(ii) >= refine_tau_cut .or. &
                       abs(cold(ii)) < refine_atol_frac*abs(total)) act(ii)=.false.
               endif
            enddo
            nnew=count(act(1:n-1))
            if (nnew == 0) exit
            if (n+nnew > cap) then
! Point budget: keep the intervals nearest the camera (least attenuated).
               flag=ior(flag,2)
               nnew=cap-n
               jj=0
               do ii=1,n-1
                  if (act(ii)) then
                     jj=jj+1
                     if (jj > nnew) act(ii)=.false.
                  endif
               enddo
               if (nnew == 0) exit
            endif
            allocate(lmid(nnew),imid(nnew),jmid(nnew),kmid(nnew))
            jj=0
            do ii=1,n-1
               if (act(ii)) then
                  jj=jj+1; imid(jj)=ii; lmid(jj)=0.5d0*(lam(ii)+lam(ii+1))
               endif
            enddo
            call refine_eval_points(lmid,bix(imid),fr,spm,epar,fac,fname,ename,jmid,kmid)
! Merge the new points and test each bisected interval against its old
! contribution (same camera-side tau at its near end).
            allocate(lam2(cap),jv2(cap),kv2(cap),bix2(cap),act2(cap))
            act2=.false.
            jj=0
            do ii=1,n
               lam2(ii+jj)=lam(ii); jv2(ii+jj)=jv(ii); kv2(ii+jj)=kv(ii); bix2(ii+jj)=bix(ii)
               if (ii < n) then
                  if (act(ii)) then
                     jj=jj+1
                     lam2(ii+jj)=lmid(jj); jv2(ii+jj)=jmid(jj); kv2(ii+jj)=kmid(jj)
                     bix2(ii+jj)=bix(ii)
                     dl=0.5d0*abs(lmid(jj)-lam(ii))*(kv(ii)+kmid(jj))
                     cnew=exp(-tv(ii))*(local_term(lam(ii),lmid(jj),jv(ii),jmid(jj),kv(ii),kmid(jj)) &
                          +exp(-dl)*local_term(lmid(jj),lam(ii+1),jmid(jj),jv(ii+1),kmid(jj),kv(ii+1)))
                     if (abs(cnew-cold(ii)) > refine_rtol*abs(cnew) .and. &
                          abs(cnew) >= refine_atol_frac*abs(total)) then
                        if (lev < refine_max_depth) then
                           act2(ii+jj-1)=.true.; act2(ii+jj)=.true.
                        else
                           flag=ior(flag,1)
                        endif
                     endif
                  endif
               endif
            enddo
            n=n+nnew; nadd=nadd+nnew
            lam(1:n)=lam2(1:n); jv(1:n)=jv2(1:n); kv(1:n)=kv2(1:n); bix(1:n)=bix2(1:n)
            act(1:n)=act2(1:n)
            deallocate(lam2,jv2,kv2,bix2,act2,lmid,imid,jmid,kmid)
         enddo
         if (nadd > 0) then
            call running_tau()
            r%I(1,:)=explin_sum()*fac*LBH
         endif
         if (allocated(refine_added)) then
            if (iray >= 1 .and. iray <= size(refine_added)) then
               refine_added(iray)=refine_added(iray)+nadd
               refine_flag(iray)=ior(refine_flag(iray),flag)
            endif
         endif
         deallocate(lam,jv,kv,tv,cold,bix,act,hv)

         contains

           subroutine running_tau()
! Same trapezoid as calc_opt_depth, from the camera (point 1).
           integer :: q
           tv(1)=0d0
           do q=2,n
              tv(q)=tv(q-1)+0.5d0*abs(lam(q)-lam(q-1))*(kv(q)+kv(q-1))
           enddo
           end subroutine running_tau

           function local_term(la,lb,ja,jb,ka,kb) result(c)
! explin contribution of one interval without the camera-side attenuation
! (radtrans_integrate_explin).
           real(kind=8), intent(in) :: la,lb,ja,jb,ka,kb
           real(kind=8) :: c,d,ed,w0,w1,sa,sb
           d=0.5d0*abs(lb-la)*(ka+kb)
           if (d < 1d-4) then
              c=0.5d0*abs(lb-la)*(ja+jb*exp(-d))
           else
              ed=exp(-d); w0=1d0-ed; w1=(1d0-ed-d*ed)/d
              sa=0d0; sb=0d0
              if (ka > 0d0) sa=ja/ka
              if (kb > 0d0) sb=jb/kb
              c=sa*(w0-w1)+sb*w1
           endif
           end function local_term

           function term(q) result(c)
           integer, intent(in) :: q
           real(kind=8) :: c
           c=exp(-tv(q))*local_term(lam(q),lam(q+1),jv(q),jv(q+1),kv(q),kv(q+1))
           end function term

           function explin_sum() result(tot)
           real(kind=8) :: tot
           integer :: q
           tot=0d0
           do q=1,n-1
              if (tv(q) > 700d0) exit
              tot=tot+term(q)
           enddo
           end function explin_sum

           logical function jump(a,b)
           real(kind=8), intent(in) :: a,b
           jump=max(abs(a),abs(b)) > 10d0*min(abs(a),abs(b))
           end function jump

         end subroutine refine_ray

         function harmonic_index(nu,bcgs,tcgs,tdimless) result(h)
! floor(nu/nu_c) with nu_c = eB/(2 pi m c), or -1 where the harmonics overlap
! (n > 2/(e^2 Theta_e), the cap of cyclotron_emis.f90) or B = 0.
! tcgs is the temperature in kelvin, or already Theta_e when tdimless (the
! APERTURE model, convert_fluidvars_aperture).  Before 2026-10-15 it was
! always read as kelvin, so for APERTURE n* ~ 1.6e9/Theta_e and every change
! of floor(nu/nu_c) triggered a split (docs/muchado.md, M4-F1).
         real(kind=8), dimension(:), intent(in) :: nu,bcgs,tcgs
         logical, intent(in) :: tdimless
         real(kind=8), dimension(size(nu)) :: h
         real(kind=8) :: nh, nstar, th
         integer :: q
         do q=1,size(nu)
            h(q)=-1d0
            if (bcgs(q) > 0d0 .and. tcgs(q) > 0d0) then
               nh=nu(q)/(EC*bcgs(q)/(2d0*PI*ME*CC))
               if (tdimless) then
                  th=tcgs(q)
               else
                  th=KB*tcgs(q)/(ME*CC*CC)
               endif
               nstar=max(1d0,2d0/(exp(2d0)*th))
               if (nh <= nstar) h(q)=aint(nh)
            endif
         enddo
         end function harmonic_index

         subroutine refine_eval_points(lm,ib,fr,spm,epar,fac,fname,ename,jout,kout)
! Invariant j/fac and |K|*LBH at new points with affine parameters lm inside
! original intervals ib, through the ray's fluid, tetrad and emissivity chain
! on private fluid and emissivity objects (the ray's own f and e are kept).
         real(kind=8), dimension(:), intent(in) :: lm
         integer, dimension(:), intent(in) :: ib
         real(kind=8), intent(in) :: fr, fac
         type (source_params), intent(in) :: spm
         type (emis_params), intent(in) :: epar
         character(len=100), intent(in) :: fname, ename
         real(kind=8), dimension(:), intent(out) :: jout, kout
         type (fluid) :: fl
         type (emis) :: el
         type (source_params) :: sps
         type (emis_params) :: eps
         type (four_vector), dimension(size(lm)) :: xs, ks
         real(kind=8), dimension(size(lm)) :: s2x, c2x, angs, rsh, cosn, nus, sr, sm
         real(kind=8), dimension(size(lm),10) :: met
         real(kind=8), dimension(size(lm)) :: ncgs, bcgs, tcgs, ncgsnth
         real(kind=8), dimension(:,:), allocatable :: nnthcgs, fnu
         real(kind=8), dimension(:), allocatable :: freqarr
         real(kind=8) :: w(4)
         integer :: ns, q, i0, a, b, n0, it
         ns=size(lm); n0=g%npts
         do q=1,ns
! cubic (Lagrange) in lambda through four geokerr points around the interval
            i0=min(max(ib(q)-1,1),n0-3)
            do a=1,4
               w(a)=1d0
               do b=1,4
                  if (b /= a) w(a)=w(a)*(lm(q)-g%lambda(i0+b-1))/(g%lambda(i0+a-1)-g%lambda(i0+b-1))
               enddo
            enddo
            xs(q)%data=0d0
            do a=1,4
               xs(q)%data=xs(q)%data+w(a)*g%x(i0+a-1)%data
            enddo
            it=ib(q)
            if (abs(lm(q)-g%lambda(ib(q)+1)) < abs(lm(q)-g%lambda(ib(q)))) it=ib(q)+1
            sr(q)=g%gk%su(1)*merge(-1d0,1d0,mod(nint(g%tprarr(it)),2)==1)
            sm(q)=g%gk%sm(1)*merge(-1d0,1d0,mod(nint(g%tpmarr(it)),2)==1)
         enddo
         ks=calc_nullp(g%gk%q2(1),g%gk%l(1),g%gk%a,xs%data(2),cos(xs%data(3)),sr,sm,1,1)
         met=kerr_metric(xs%data(2),xs%data(3),g%gk%a)
         call assign_metric(xs,transpose(met))
         call assign_metric(ks,transpose(met))
! Start from copies of the ray's own objects so every scalar setting (nfreq,
! nrelbin, ...) matches, then reallocate them for the ns new points.
         fl=f
         call del_fluid_model(fl)
         call initialize_fluid_model(fl,fname,g%gk%a,ns)
         el=e
         call del_emissivity(el)
         call select_emissivity_values(el,ename)
         call get_fluid_vars(xs,ks,g%gk%a,fl)
         s2x=0d0; c2x=0d0
         call comoving_ortho(xs%data(2),xs%data(3),g%gk%a, &
              g%gk%alpha(1),g%gk%beta(1),g%gk%mu0,fl%u,fl%b,ks,s2x,c2x, &
              angs,rsh,cosn)
         sps=spm
         call initialize_source_params(sps,ns)
         eps=epar
         call initialize_emis_params(eps,ns)
         call initialize_emissivity(el,ns,fl%nfreq,rsh,angs,cosn,fl%nrelbin,fl%bingammamin,fl%bingammamax)
         allocate(nnthcgs(ns,fl%nrelbin))
         call convert_fluid_vars(fl,ncgs,ncgsnth,nnthcgs,bcgs,tcgs,fnu,freqarr,sps)
         call assign_emis_params(el,ncgs,ncgsnth,nnthcgs,bcgs,tcgs,fnu,freqarr,fl%nfreq)
         eps%gmin=sps%gmin
         eps%mu=sps%mu
         call emis_model(el,eps)
         nus=fr/rsh
         call calc_emissivity(nus,el,eps)
         call invariant_emis(el,rsh)
         jout=el%j(:,1)/fac
         kout=abs(el%K(:,1))*LBH
         call del_source_params(sps)
         call del_emis_params(eps)
         call del_emissivity(el)
         call del_fluid_model(fl)
         end subroutine refine_eval_points

         subroutine init_refine_stats(nrays)
         integer, intent(in) :: nrays
         if (allocated(refine_added)) deallocate(refine_added,refine_flag)
         allocate(refine_added(nrays),refine_flag(nrays))
         refine_added=0; refine_flag=0
         end subroutine init_refine_stats

         subroutine write_refine_stats(outname)
! Summary to stdout and per-ray counts to outname: ray index, points added,
! cap flag (1 = refine_max_depth reached unconverged, 2 = refine_max_points
! ran out, 3 = both), for every ray that was refined or capped.
         character(len=*), intent(in) :: outname
         integer :: q, u
         if (.not. allocated(refine_added)) return
         write(6,'(a,i8,a,i10,a,i8,a,i8,a,i8)') 'refine: rays refined ', count(refine_added > 0), &
              '  points added ', sum(int(refine_added,8)), '  max per ray ', maxval(refine_added), &
              '  depth cap ', count(iand(refine_flag,1) /= 0), '  point cap ', count(iand(refine_flag,2) /= 0)
         open(newunit=u,file=outname,form='formatted',status='replace')
         write(u,'(a)') '# ray  points_added  flag(1=max_depth,2=max_points)'
         do q=1,size(refine_added)
            if (refine_added(q) > 0 .or. refine_flag(q) /= 0) &
                 write(u,'(i8,1x,i8,1x,i2)') q, refine_added(q), refine_flag(q)
         enddo
         close(u)
         end subroutine write_refine_stats

         subroutine convert_model(sp)
!         type (fluid), intent(in) :: f
!         type (emis), intent(inout) :: e
         type (source_params), intent(inout) :: sp
         real(kind=8), dimension(e%npts) :: ncgs,bcgs,tcgs,ncgsnth
         real(kind=8), dimension(e%npts,f%nrelbin) :: nnthcgs !AC f%nrelbin or e%nrelbin?
         real(kind=8), dimension(:), allocatable :: freqarr
         real(kind=8), dimension(:,:), allocatable :: fnu
         call convert_fluid_vars(f,ncgs,ncgsnth,nnthcgs,bcgs,tcgs,fnu,freqarr,sp)
!         write(6,*) 'CONVERTED'
!         write(6,*) 'sz: ',size(ncgs), size(bcgs), size(tcgs)
!         write(6,*) 'alloc: ',allocated(e%ncgs),allocated(e%bcgs),allocated(e%tcgs)
         call assign_emis_params(e,ncgs,ncgsnth,nnthcgs,bcgs,tcgs,fnu,freqarr,f%nfreq)
         end subroutine convert_model

         subroutine transpol(rshift,s2psi,c2psi,cosne)
! Subroutine to do CPS80 polarization for single points
! For now just uses Chandrasekhar table for degree, but could change that
! JAD 4/20/2012
         real(kind=8), dimension(:), intent(in) :: rshift
         real, dimension(g%npts) :: interpI,interpdel
         real(kind=8), dimension(g%npts) :: polarpsi
         real(kind=8), dimension(g%npts), intent(out) :: s2psi,cosne,c2psi
!         write(6,*) 'rmu: ',g%x%data(2), g%x%data(3)
           call calc_polar_psi(g%x%data(2),cos(g%x%data(3)), &
                  g%gk%q2,g%gk%a, &
                  g%gk%alpha, &
                  g%gk%beta,rshift,g%gk%mu0,g%k, &
                  c2psi,s2psi,cosne)
!         write(6,*) 'polar psi: ',cosne,polarpsi
           call interp_chandra_tab24(real(cosne),interpI,interpdel)
           e%j(:,1)=e%j(:,1)*dble(interpI)
!           e%j(:,2)=e%j(:,1)*cos(polarpsi)*dble(interpdel)
!           e%j(:,3)=e%j(:,1)*sin(polarpsi)*dble(interpdel)
           e%j(:,2)=e%j(:,1)*c2psi*dble(interpdel)
           e%j(:,3)=e%j(:,1)*s2psi*dble(interpdel)
           e%j(:,4)=0.
         end subroutine transpol

         subroutine grtrans_compute_intensity()
! Subroutine to compute intensity at single point rather than integrate
! JAD 10/13/2011
!         write(6,*) 'grtrans intensity: ',e%neq,e%npts,r%neq,g%npts
!         write(6,*) 'grtrans intensity: ',size(r%I), size(e%j)
         r%I(1:e%neq,1:e%npts)=transpose(e%j(1:e%npts,1:e%neq))
         end subroutine grtrans_compute_intensity

       end module grtrans
