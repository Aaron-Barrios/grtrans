! fits.f90  (cfitsio STUB — Option A)
!
! The original module wrapped cfitsio for FITS file I/O.  We only need the
! binary (cflag=0) and in-memory (pgrtrans) output paths, so every cfitsio call
! is replaced with a no-op that reports success (status=0).  grtrans now builds
! and links without libcfitsio.  Attempting FITS output (cflag=1) is a silent
! no-op.

module fits

  implicit none

  interface create_fits
    module procedure create_fits
  end interface

  interface delete_fits
    module procedure delete_fits
  end interface

  interface open_fits
    module procedure open_fits
  end interface

  interface close_fits
    module procedure close_fits
  end interface

  interface write_fits_array
    module procedure write_fits_array_int
    module procedure write_fits_image_int
    module procedure write_fits_array_dp
    module procedure write_fits_image_dp
    module procedure write_fits_array_real
    module procedure write_fits_image_real
  end interface

  interface read_fits_array
    module procedure read_fits_array_int
    module procedure read_fits_image_int
    module procedure read_fits_array_dp
    module procedure read_fits_image_dp
    module procedure read_fits_array_real
    module procedure read_fits_image_real
  end interface

  interface new_fits_extension
    module procedure new_fits_extension
  end interface

  interface write_fits_extra_keyword
    module procedure write_fits_extra_keyword_int
    module procedure write_fits_extra_keyword_real
    module procedure write_fits_extra_keyword_dp
  end interface

  interface write_fits_header
    module procedure write_fits_header_ps
    module procedure write_fits_header_pr
  end interface

contains

  subroutine create_fits(unit, filename, status)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: filename
    integer, intent(out) :: status
    status = 0
  end subroutine create_fits

  subroutine delete_fits(filename, status)
    character(len=*), intent(in) :: filename
    integer, intent(inout) :: status
    status = 0
  end subroutine delete_fits

  subroutine open_fits(unit, filename, readwrite, status)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: filename
    integer, intent(in) :: readwrite
    integer, intent(out) :: status
    status = 0
  end subroutine open_fits

  subroutine close_fits(unit, status)
    integer, intent(in) :: unit
    integer, intent(inout) :: status
    status = 0
  end subroutine close_fits

  subroutine new_fits_extension(unit, status)
    integer, intent(in) :: unit
    integer, intent(inout) :: status
    status = 0
  end subroutine new_fits_extension

  subroutine write_fits_header_ps(unit, bitpix, naxis, naxes, status)
    integer, intent(in) :: unit, bitpix, naxis
    integer, intent(in), dimension(naxis) :: naxes
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_header_ps

  subroutine write_fits_header_pr(unit, simple, bitpix, naxis, naxes, &
       pcount, gcount, extend, status)
    integer, intent(in) :: unit, bitpix, naxis, pcount, gcount
    integer, intent(in), dimension(naxis) :: naxes
    logical, intent(in) :: simple, extend
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_header_pr

  subroutine write_fits_array_int(unit, group, fpixel, nelements, values, status)
    integer, intent(in) :: unit, group, fpixel, nelements
    integer, intent(in), dimension(:) :: values
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_array_int

  subroutine write_fits_image_int(unit, group, fpixel, nelements, values, status)
    integer, intent(in) :: unit, group, fpixel, nelements
    integer, intent(in), dimension(:,:) :: values
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_image_int

  subroutine write_fits_array_real(unit, group, fpixel, nelements, values, status)
    integer, intent(in) :: unit, group, fpixel, nelements
    real, intent(in), dimension(:) :: values
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_array_real

  subroutine write_fits_image_real(unit, group, fpixel, nelements, values, status)
    integer, intent(in) :: unit, group, fpixel, nelements
    real, intent(in), dimension(:,:) :: values
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_image_real

  subroutine write_fits_array_dp(unit, group, fpixel, nelements, values, status)
    integer, intent(in) :: unit, group, fpixel, nelements
    real(kind=8), intent(in), dimension(:) :: values
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_array_dp

  subroutine write_fits_image_dp(unit, group, fpixel, nelements, values, status)
    integer, intent(in) :: unit, group, fpixel, nelements
    real(kind=8), intent(in), dimension(:,:) :: values
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_image_dp

  subroutine read_fits_array_int(unit, group, fpixel, nelements, nullval, values, anyf, status)
    integer, intent(in) :: unit, group, fpixel, nelements, nullval
    integer, intent(out), dimension(:) :: values
    integer, intent(out) :: anyf
    integer, intent(inout) :: status
    values = 0; anyf = 0; status = 0
  end subroutine read_fits_array_int

  subroutine read_fits_image_int(unit, group, fpixel, nelements, nullval, values, anyf, status)
    integer, intent(in) :: unit, group, fpixel, nelements, nullval
    integer, intent(out), dimension(:,:) :: values
    integer, intent(out) :: anyf
    integer, intent(inout) :: status
    values = 0; anyf = 0; status = 0
  end subroutine read_fits_image_int

  subroutine read_fits_array_real(unit, group, fpixel, nelements, nullval, values, anyf, status)
    integer, intent(in) :: unit, group, fpixel, nelements, nullval
    real, intent(out), dimension(:) :: values
    integer, intent(out) :: anyf
    integer, intent(inout) :: status
    values = 0.; anyf = 0; status = 0
  end subroutine read_fits_array_real

  subroutine read_fits_image_real(unit, group, fpixel, nelements, nullval, values, anyf, status)
    integer, intent(in) :: unit, group, fpixel, nelements, nullval
    real, intent(out), dimension(:,:) :: values
    integer, intent(out) :: anyf
    integer, intent(inout) :: status
    values = 0.; anyf = 0; status = 0
  end subroutine read_fits_image_real

  subroutine read_fits_array_dp(unit, group, fpixel, nelements, nullval, values, anyf, status)
    integer, intent(in) :: unit, group, fpixel, nelements, nullval
    real(kind=8), intent(out), dimension(:) :: values
    integer, intent(out) :: anyf
    integer, intent(inout) :: status
    values = 0d0; anyf = 0; status = 0
  end subroutine read_fits_array_dp

  subroutine read_fits_image_dp(unit, group, fpixel, nelements, nullval, values, anyf, status)
    integer, intent(in) :: unit, group, fpixel, nelements, nullval
    real(kind=8), intent(out), dimension(:,:) :: values
    integer, intent(out) :: anyf
    integer, intent(inout) :: status
    values = 0d0; anyf = 0; status = 0
  end subroutine read_fits_image_dp

  subroutine write_fits_extra_keyword_int(unit, kname, kval, decimal, kdesc, status)
    integer, intent(in) :: unit, kval, decimal
    character(len=20), intent(in) :: kname, kdesc
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_extra_keyword_int

  subroutine write_fits_extra_keyword_real(unit, kname, kval, decimal, kdesc, status)
    integer, intent(in) :: unit, decimal
    real, intent(in) :: kval
    character(len=20), intent(in) :: kname, kdesc
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_extra_keyword_real

  subroutine write_fits_extra_keyword_dp(unit, kname, kval, decimal, kdesc, status)
    integer, intent(in) :: unit, decimal
    real(kind=8), intent(in) :: kval
    character(len=20), intent(in) :: kname, kdesc
    integer, intent(inout) :: status
    status = 0
  end subroutine write_fits_extra_keyword_dp

  subroutine aftpkyj(unit, kname, kval, kdesc, status)
    integer, intent(in) :: unit, kval
    character(len=20), intent(in) :: kname, kdesc
    integer, intent(inout) :: status
    status = 0
  end subroutine aftpkyj

  subroutine aftpkyd(unit, kname, kval, decimal, kdesc, status)
    integer, intent(in) :: unit, decimal
    real(kind=8), intent(in) :: kval
    character(len=20), intent(in) :: kname, kdesc
    integer, intent(inout) :: status
    status = 0
  end subroutine aftpkyd

  subroutine aftpkys(unit, kname, kval, kdesc, status)
    integer, intent(in) :: unit
    character(len=20), intent(in) :: kname, kval, kdesc
    integer, intent(inout) :: status
    status = 0
  end subroutine aftpkys

end module fits
