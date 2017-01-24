module FatesGlobals
  ! NOTE(bja, 201608) This is a temporary hack module to store global
  ! data used inside fates. It's main use it to explicitly call out
  ! global data that needs to be dealt with, but doesn't have an
  ! immediately obvious home.

  use FatesConstantsMod         , only : r8 => fates_r8
!  use EDTypesMod                , only : cp_nclmax, cp_nlevcan, numpft_ed
   
  implicit none

  public :: FatesGlobalsInit
  public :: fates_log
  public :: fates_global_verbose
  public :: SetFatesTime
  public :: set_fates_global_elements


  ! for setting number of patches per gridcell and number of cohorts per patch
  ! for I/O and converting to a vector

  integer, parameter :: maxPatchesPerSite  = 10   ! maximum number of patches to live on a site
  integer, parameter :: maxCohortsPerPatch = 160  ! maximum number of cohorts to live on a patch


  ! Variables mostly used for dimensioning host land model (HLM) array spaces

  integer, protected :: maxElementsPerPatch       ! maxElementsPerPatch is the value that is ultimately
                                                  ! used to set the size of the largest arrays necessary
                                                  ! in things like restart files (probably hosted by the 
                                                  ! HLM). The size of these arrays are not a parameter
                                                  ! because it is simply the maximum of several different
                                                  ! dimensions. It is possible that this would be the
                                                  ! maximum number of cohorts per patch, but
                                                  ! but it could be other things.

  integer, protected :: maxElementsPerSite        ! This is the max number of individual items one can store per 
                                                  ! each grid cell and effects the striding in the ED restart 
                                                  ! data as some fields are arrays where each array is
                                                  ! associated with one cohort

  integer, protected :: maxCohortsPerSite         ! Maximum number of cohorts that can exist in a given
                                                  ! site. Its possible this is not used.


  integer, parameter :: cp_nclmax = 2       ! Maximum number of canopy layers

  integer, parameter :: cp_nlevcan = 40     ! number of leaf layers in canopy layer

  integer , parameter :: numpft_ed = 2      ! number of PFTs used in ED. 


  ! -------------------------------------------------------------------------------------
  ! Timing Variables
  ! It is assumed that all of the sites on a given machine will be synchronous.
  ! It is also assumed that the HLM will control time.
  ! -------------------------------------------------------------------------------------
  integer, protected  :: current_year    ! Current year
  integer, protected  :: current_month   ! month of year
  integer, protected  :: current_day     ! day of month
  integer, protected  :: current_tod     ! time of day (seconds past 0Z)
  integer, protected  :: current_date    ! time of day (seconds past 0Z)
  integer, protected  :: reference_date  ! YYYYMMDD
  real(r8), protected :: model_day       ! elapsed days between current date and reference
  integer, protected  :: day_of_year     ! The integer day of the year
  integer, protected  :: days_per_year   ! The HLM controls time, some HLMs may include a leap
  real(r8), protected :: freq_day        ! fraction of year for daily time-step (1/days_per_year)
                                         ! this is a frequency

  integer, private :: fates_log_
  logical, private :: fates_global_verbose_

contains

  subroutine set_fates_global_elements()
     implicit none

     maxElementsPerPatch = max(maxCohortsPerPatch, &
           numpft_ed * cp_nclmax * cp_nlevcan)
     
     maxCohortsPerSite = maxPatchesPerSite * maxCohortsPerPatch
     
     maxElementsPerSite = maxPatchesPerSite * maxElementsPerPatch

  end subroutine set_fates_global_elements

  ! =====================================================================================

  subroutine FatesGlobalsInit(log_unit,global_verbose)

    implicit none

    integer, intent(in) :: log_unit
    logical, intent(in) :: global_verbose

    fates_log_ = log_unit
    fates_global_verbose_ = global_verbose

  end subroutine FatesGlobalsInit

  ! =====================================================================================

  integer function fates_log()
    fates_log = fates_log_
  end function fates_log

  logical function fates_global_verbose()
    fates_global_verbose = fates_global_verbose_
  end function fates_global_verbose

  subroutine fates_endrun(msg) 

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    ! This subroutine was derived from CLM's
    ! endrun_vanilla() in abortutils.F90
    !
    use shr_sys_mod , only: shr_sys_abort
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in), optional :: msg    ! string to be printed
    !-----------------------------------------------------------------------

    if (present (msg)) then
       write(fates_log(),*)'ENDRUN:', msg
    else
       write(fates_log(),*)'ENDRUN: called without a message string'
    end if

    call shr_sys_abort()

  end subroutine fates_endrun

  ! =====================================================================================

  subroutine SetFatesTime(current_year_in, current_month_in, &
                          current_day_in, current_tod_in, &
                          current_date_in, reference_date_in, &
                          model_day_in, day_of_year_in, &
                          days_per_year_in, freq_day_in)

     ! This subroutine should be called directly from the HLM
     
     integer,  intent(in) :: current_year_in
     integer,  intent(in) :: current_month_in
     integer,  intent(in) :: current_day_in
     integer,  intent(in) :: current_tod_in
     integer,  intent(in) :: current_date_in
     integer,  intent(in) :: reference_date_in
     real(r8), intent(in) :: model_day_in
     integer,  intent(in) :: day_of_year_in
     integer,  intent(in) :: days_per_year_in
     real(r8), intent(in) :: freq_day_in

     current_year   = current_year_in
     current_month  = current_month_in
     current_day    = current_day_in
     current_tod    = current_tod_in
     current_date   = current_date_in
     reference_date = reference_date_in
     model_day      = model_day_in
     day_of_year    = day_of_year_in
     days_per_year  = days_per_year_in
     freq_day       = freq_day_in

  end subroutine SetFatesTime

end module FatesGlobals
