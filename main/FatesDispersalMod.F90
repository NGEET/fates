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
     
      integer  :: neighbor_count   ! total neighbors near source

   end type neighborhood_type
   
   ! Dispersal type
   type, public :: dispersal_type
   
      real(r8), allocatable :: outgoing_local(:,:)    ! local gridcell array of outgoing seeds, gridcell x pft
      real(r8), allocatable :: outgoing_global(:,:)   ! global accumulation array of outgoing seeds, gridcell x pft
      real(r8), allocatable :: incoming_global(:,:)   ! 
      
      contains
      
         procedure :: init
      
   end type dispersal_type

   type(neighborhood_type), public, pointer :: lneighbors(:)

   public :: ProbabilityDensity
   public :: IsItDispersalTime

   integer :: prev_dispersal_date = 0
   
   character(len=*), parameter, private :: sourcefile = __FILE__

contains

   ! ====================================================================================

   subroutine init(this, numprocs, numpft)
      
      ! Use
      use EDPftvarcon           , only : EDPftvarcon_inst 
      use FatesConstantsMod     , only : fates_check_param_set
      use FatesInterfaceTypesMod, only : fates_dispersal_kernel_mode
      use FatesInterfaceTypesMod, only : fates_dispersal_kernel_none
      
      ! Arguments
      class(dispersal_type), intent(inout) :: this
      integer, intent(in) ::  numprocs
      integer, intent(in) ::  numpft
      
      ! Check if seed dispersal mode is 'turned on' by checking the parameter values
      ! This assumes we consistency in the parameter file across all pfts, i.e. either
      ! all 'on' or all 'off'
      if (fates_dispersal_kernel_mode .eq. fates_dispersal_kernel_none) return 
      
      allocate(this%outgoing_local(numprocs,numpft))
      allocate(this%outgoing_global(numprocs,numpft))
      allocate(this%incoming_global(numprocs,numpft))
   
      this%outgoing_local(:,:) = 0._r8
      this%outgoing_global(:,:) = 0._r8
      this%incoming_global(:,:) = 0._r8
      
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
      ! Note that fates_dispersal_kernel_none mode is checked prior to this call being made
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

      write(fates_log(),*) 'ipft,dist,pd: ', ipft, dist, pd

   end subroutine ProbabilityDensity

   ! ====================================================================================

   real(r8) function PD_exponential(dist, ipft)
      
      use EDPftvarcon           , only : EDPftvarcon_inst
   
      ! Arguments
      real(r8), intent(in) :: dist
      integer, intent(in)  :: ipft
      
      ! Assuming simple exponential decay.  In the future perhaps this could be an interface
      ! for different weight calculations (and could be held only in fates)
      
      PD_exponential = exp(-EDPftvarcon_inst%seed_dispersal_param_A(ipft)*dist)
      write(fates_log(),*) 'ipft,dist,PD_exp: ', ipft, dist, PD_exponential

   end function PD_exponential

   ! ====================================================================================

   real(r8) function PD_exppower(dist, ipft)
      
      use EDPftvarcon           , only : EDPftvarcon_inst
   
      ! Arguments
      real(r8), intent(in) :: dist
      integer, intent(in)  :: ipft
      
      associate(&
         param_a => EDPftvarcon_inst%seed_dispersal_param_A(ipft), &
         param_b => EDPftvarcon_inst%seed_dispersal_param_B(ipft))      
      
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
         param_a => EDPftvarcon_inst%seed_dispersal_param_A(ipft), &
         param_b => EDPftvarcon_inst%seed_dispersal_param_B(ipft))      
      
      PD_logsech = (1/(pi_const**2 * param_b * dist**2)) / &
                    ((dist/param_a)**(1/param_b) + &
                    (dist/param_a)**(-1/param_b))
                    
      end associate

   end function PD_logsech

   ! ====================================================================================
   
   logical function IsItDispersalTime()
   
   ! Determine if seeds should be dispersed across gridcells.  This eventually could be
   ! driven by plant reproduction dynamics.  For now this is based strictly on a calendar
   
   use FatesInterfaceTypesMod,      only : hlm_current_day, &
                                      hlm_current_month, & 
                                      hlm_current_year, &
                                      hlm_current_date, &
                                      fates_dispersal_cadence, &
                                      fates_dispersal_cadence_daily, &
                                      fates_dispersal_cadence_monthly, &
                                      fates_dispersal_cadence_yearly

   ! LOCAL
   integer :: check_date = 0
   
   ! initialize the return value as false
   IsItDispersalTime = .false.
   
   ! Select the date type to check against based on the dispersal candence
   select case(fates_dispersal_cadence)
   case (fates_dispersal_cadence_daily)
      check_date = hlm_current_day
   case (fates_dispersal_cadence_monthly)
      check_date = hlm_current_month
   case (fates_dispersal_cadence_yearly)
      check_date = hlm_current_year
   case default
      write(fates_log(),*) 'ERROR: An undefined dispersal cadence was specified: ', fates_dispersal_cadence
      call endrun(msg=errMsg(sourcefile, __LINE__))
   end select
   
   ! Determine if it is a new day, month, or year.  If true update the previous date to the current value   
   if (check_date .ne. prev_dispersal_date) then
      IsItDispersalTime = .true.
      prev_dispersal_date = check_date
   end if
                                                                            
   end function IsItDispersalTime

      
end module FatesDispersalMod
