       ! Module for radiative transfer objects including intensities
       ! and methods for computing emissivities and for 
       ! integrating the radiative transfer equation
       module class_rad_trans
       use grtrans_inputs
       implicit none
       type rad_trans
         real(kind=8), dimension(:,:), allocatable :: I
         real(kind=8), dimension(:), allocatable :: tau
         real(kind=8) :: nu
         integer :: neq,iflag,onpts,npts,ename,idex
       end type rad_trans
        
       interface initialize_rad_trans
         module procedure initialize_rad_trans
       end interface

       interface del_rad_trans
         module procedure del_rad_trans
       end interface

       contains
         
         subroutine initialize_rad_trans(r,iname,npts,neq,extra)
         type (rad_trans), intent(out) :: r
         integer, intent(in) :: npts, neq,extra
         character(len=20), intent(in) :: iname
         integer :: nextra
         if(iname=='lsoda') then
            r%iflag=0
         else if (iname=='delo') then
            r%iflag=1
         else if (iname=='formal') then
            r%iflag=2
         else if (iname=='lsodasph') then
            r%iflag=3
         else if (iname=='explin') then
            ! exponential transfer step, source function linear in tau
            ! (radtrans_integrate_explin; total intensity only)
            r%iflag=4
         endif
!         write(6,*) 'rneq: ',neq
         r%neq=neq
! Determine which integration routine we're using:
         r%onpts=npts
! Allocate intensity array including intermediate points:
         allocate(r%I(r%neq,r%onpts))
         r%I(:,:)=0d0; r%npts=r%onpts
         if(extra.eq.1) then
! 21 = 19 stock slots + (20) emission-weighted plasma beta + (21) sigma.
! MUST match nextra in pgrtrans.f90 (which sizes the camera's pixvals), or the
! (/real(r%I), real(r%tau)/) constructor passed to save_raytrace_camera_pixel
! has the wrong length.  Always 21, including the npts==1 branch: that branch
! only fills slots 1-7 and the rest stay zeroed, which is self-consistent and
! avoids the pre-existing size mismatch when a ray has a single point.
            nextra=21
            allocate(r%tau(nextra)); r%tau(:)=0d0
         endif
!         write(6,*) 'r%tau, r%I init: ',minval(r%tau),maxval(r%tau),minval(r%I),maxval(r%I)
         r%ename=0
         return
         end subroutine initialize_rad_trans

         subroutine del_rad_trans(r)
         type (rad_trans), intent(inout) :: r
         deallocate(r%I)
         if(allocated(r%tau)) deallocate(r%tau)
         r%ename=0
         r%iflag=0
         end subroutine del_rad_trans

       end module class_rad_trans
