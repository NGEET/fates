module FatesDispersalMod

   use shr_log_mod           , only : errMsg => shr_log_errMsg
   use FatesGlobals          , only : endrun => fates_endrun
   use FatesGlobals          , only : fates_log
   use FatesConstantsMod     , only : r8 => fates_r8
   use FatesConstantsMod     , only : pi_const
   use FatesInterfaceTypesMod, only : numpft

   implicit none
   
   ! Neighbor node
   type, public :: neighbor_type
      
      ! Grid cell neighbor
      type(neighbor_type), pointer :: next_neighbor => null() 
   
      integer               :: gindex           ! grid cell index
      real(r8)              :: gc_dist          ! distance between source and neighbor
      real(r8), allocatable :: density_prob(:)  ! probability density from source per pft
      
   end type neighbor_type

   ! Neighborhood linked list
   type, public :: neighborhood_type

      ! Linked list of neighbors for a given source grid cell
      type(neighbor_type), pointer :: first_neighbor => null()
      type(neighbor_type), pointer :: last_neighbor => null()
     
      integer                :: neighbor_count   ! total neighbors near source
      integer, allocatable  :: neighbor_indices(:) ! list of gridcell indices

   end type neighborhood_type
   
   ! Dispersal type
   type, public :: dispersal_type
   
      real(r8), allocatable :: outgoing_local(:,:)    ! local buffer array of outgoing seeds, local gridcell x pft
      real(r8), allocatable :: outgoing_global(:,:)   ! global accumulation buffer array of outgoing seeds, global gridcell x pft
      real(r8), allocatable :: incoming_global(:,:)   ! local buffer array used to calculate incoming seeds based on nearest neighbors
      integer,  allocatable :: ncells_array(:)        ! local array with the number of gridcells per process for each rank index
      integer,  allocatable :: begg_array(:)          ! local array with the starting index of each gridcell for each rank index
           
      contains
      
         procedure :: init
      
   end type dispersal_type

   type(neighborhood_type), public, pointer :: lneighbors(:)

   public :: ProbabilityDensity
   public :: IsItDispersalTime
   
   integer :: dispersal_date = 0       ! Last candence date in which there was a dispersal
   logical :: dispersal_flag = .false. ! Have seeds been disperesed globally
  
   character(len=*), parameter, private :: sourcefile = __FILE__

contains

   ! ====================================================================================

   subroutine init(this, numprocs, numgc_global, numgc_local, numpft)
      
      ! Use
      use EDPftvarcon           , only : EDPftvarcon_inst 
      use FatesConstantsMod     , only : fates_check_param_set, fates_unset_int

      ! Arguments
      class(dispersal_type), intent(inout) :: this

      integer, intent(in) ::  numprocs      ! number of processors (across all nodes)
      integer, intent(in) ::  numgc_global  ! number of gridcells across all processors
      integer, intent(in) ::  numgc_local   ! number of gridcells on this processor
      integer, intent(in) ::  numpft        ! number of FATES pfts
      
      allocate(this%outgoing_local(numpft,numgc_local))
      allocate(this%outgoing_global(numpft,numgc_global))
      allocate(this%incoming_global(numpft,numgc_global))
      allocate(this%begg_array(numprocs))
      allocate(this%ncells_array(numprocs))
   
      this%outgoing_local(:,:) = 0._r8
      this%outgoing_global(:,:) = 0._r8
      this%incoming_global(:,:) = 0._r8
      this%ncells_array(:) = fates_unset_int
      this%begg_array(:) = fates_unset_int
      
      ! Set the dispersal date to the current date.  Dispersal will start at the end of
      ! current initial date
      dispersal_date = GetCadenceDate()           
      
   end subroutine init

   ! ====================================================================================

   subroutine ProbabilityDensity(pd, ipft, dist)

      ! Main subroutine that calls difference routines based on case select mode

      ! Use
      use FatesInterfaceTypesMod, only : fates_dispersal_kernel_exponential, &
                                         fates_dispersal_kernel_exppower, &
                                         fates_dispersal_kernel_logsech, &
                                         fates_dispersal_kernel_mode

      ! Arguments
      real(r8), intent(out) :: pd   ! Probability density
      integer, intent(in)   :: ipft ! pft index
      real(r8), intent(in)  :: dist ! distance

      ! Select the function to use based on the kernel mode
      ! Note that hlm_seeddisp_cadence is checked prior to this call being made
      select case(fates_dispersal_kernel_mode)

      case (fates_dispersal_kernel_exponential)
         pd = PD_exponential(dist,ipft)
      case (fates_dispersal_kernel_exppower)
         pd = PD_exppower(dist,ipft)
      case (fates_dispersal_kernel_logsech)
         pd = PD_logsech(dist,ipft)
      case default
         write(fates_log(),*) 'ERROR: An undefined dispersal kernel was specified: ', fates_dispersal_kernel_mode
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select

   end subroutine ProbabilityDensity

   ! ====================================================================================

   real(r8) function PD_exponential(dist, ipft)
      
      use EDPftvarcon           , only : EDPftvarcon_inst
   
      ! Arguments
      real(r8), intent(in) :: dist
      integer, intent(in)  :: ipft
      
      ! Assuming simple exponential decay.  In the future perhaps this could be an interface
      ! for different weight calculations (and could be held only in fates)
      
      PD_exponential = exp(-EDPftvarcon_inst%seed_dispersal_pdf_scale(ipft)*dist)

   end function PD_exponential

   ! ====================================================================================

   real(r8) function PD_exppower(dist, ipft)
      
      use EDPftvarcon           , only : EDPftvarcon_inst
   
      ! Arguments
      real(r8), intent(in) :: dist
      integer, intent(in)  :: ipft
      
      associate(&
         param_a => EDPftvarcon_inst%seed_dispersal_pdf_scale(ipft), &
         param_b => EDPftvarcon_inst%seed_dispersal_pdf_shape(ipft))      
      
      PD_exppower = (param_b / (2*pi_const*gamma(2/param_b))) * &
                    exp(-(dist**param_b)/(param_a**param_b))
                    
      end associate

   end function PD_exppower

   ! ====================================================================================

   real(r8) function PD_logsech(dist, ipft)
      
      use EDPftvarcon           , only : EDPftvarcon_inst
   
      ! Arguments
      real(r8), intent(in) :: dist
      integer, intent(in)  :: ipft

      associate(&
         param_a => EDPftvarcon_inst%seed_dispersal_pdf_scale(ipft), &
         param_b => EDPftvarcon_inst%seed_dispersal_pdf_shape(ipft))      
      
      PD_logsech = (1/(pi_const**2 * param_b * dist**2)) / &
                    ((dist/param_a)**(1/param_b) + &
                    (dist/param_a)**(-1/param_b))
                    
      end associate

   end function PD_logsech

   ! ====================================================================================
   
   logical function IsItDispersalTime(setdispersedflag)
   
   ! Determine if seeds should be dispersed across gridcells.  This eventually could be
   ! driven by plant reproduction dynamics.  For now this is based strictly on a calendar.
   ! This function attempts to wrap up all the logic for the dispersal code, which is
   ! takes place at multiple points in the code, into a single function.
   ! The logic for the function is as follows:
   !     - If new date and seeds not globally dispersed, pass local seed to global bufffer 
   !       (see wrap_update_hlmfates_dyn)
   !     - If new date and seeds not globally dispersed, globally disperse seeds (call WrapSeedGlobal) and
   !       set dispersed flag to true.  Note that since this must happen outside of threaded region,
   !       this comes after fates dynamic_driv procedure.
   !     - If new date or not and seeds have been globally dispersed, call wrap_seed_disperse
   !       to pass dispersed seeds to fates.  Set dispersed flag to false.  Given that this must 
   !       happen after WrapSeedGlobal, but can be threaded this takes place at the top of the 
   !       dynamics_driv call.
   
   ! Arguments
   logical, optional :: setdispersedflag ! Has the global dispersal been completed?

   ! Local
   logical :: setflag
   
   ! The default return value is false
   IsItDispersalTime = .false.
   
   ! Check if set dispersal flag is provided.  This should be provided during a check
   ! when the flag should be set to true after the global dispersal
   setflag = .false.
   if (present(setdispersedflag)) then
      setflag = setdispersedflag
   end if
      
   ! If dispersal flag is true, regardless of the date, pass dispersed seeds to fates and reset flag
   ! If dispersal flag is false, check if it is time to disperse
   ! If it's time to disperse, check to see if the dispersal flag should be set true and last
   ! dispersal date updated
   if (dispersal_flag) then
      IsItDispersalTime = .true.
      dispersal_flag = .false.
   else
      if (GetCadenceDate() .ne. dispersal_date) then
         IsItDispersalTime = .true.
         if (setflag) then
            dispersal_flag = .true.
            dispersal_date = GetCadenceDate()
         end if
      end if
   end if

   end function IsItDispersalTime
   
   ! ====================================================================================

   integer function GetCadenceDate()
   
   use FatesInterfaceTypesMod, only : hlm_current_day, &
                                      hlm_current_month, & 
                                      hlm_current_year, &
                                      hlm_current_date, &
                                      hlm_seeddisp_cadence, &
                                      fates_dispersal_cadence_daily, &
                                      fates_dispersal_cadence_monthly, &
                                      fates_dispersal_cadence_yearly

   ! Select the date type to check against based on the dispersal candence
   select case(hlm_seeddisp_cadence)
   case (fates_dispersal_cadence_daily)
      GetCadenceDate = hlm_current_day
   case (fates_dispersal_cadence_monthly)
      GetCadenceDate = hlm_current_month
   case (fates_dispersal_cadence_yearly)
      GetCadenceDate = hlm_current_year
   case default
      write(fates_log(),*) 'ERROR: An undefined dispersal cadence was specified: ', hlm_seeddisp_cadence
      call endrun(msg=errMsg(sourcefile, __LINE__))
   end select
                                      
   end function GetCadenceDate

      
end module FatesDispersalMod
