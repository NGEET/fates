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

   character(len=*), parameter, private :: sourcefile = __FILE__

contains

   ! ====================================================================================

   subroutine init(this, numprocs, numpft)
      
      class(dispersal_type), intent(inout) :: this
      integer, intent(in) ::  numprocs
      integer, intent(in) ::  numpft
      
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
      use FatesInterfaceTypesMod, only : hlm_dispersal_kernel_exponential, &
                                         hlm_dispersal_kernel_exppower, &
                                         hlm_dispersal_kernel_logsech, &
                                         hlm_dispersal_kernel_mode

      ! Arguments
      real(r8), intent(out) :: pd   ! Probability density
      integer, intent(in)   :: ipft ! pft index
      real(r8), intent(in)  :: dist ! distance

      ! Local - temp
      ! real(r8) :: param_a = 1._r8
      ! real(r8) :: param_b = 1._r8

      hlm_dispersal_kernel_mode = 1

      select case(hlm_dispersal_kernel_mode)

      case (hlm_dispersal_kernel_exponential)
         pd = PD_exponential(dist)
      case (hlm_dispersal_kernel_exppower)
         pd = PD_exppower(dist)
      case (hlm_dispersal_kernel_logsech)
         pd = PD_logsech(dist)
      case default
         write(fates_log(),*) 'ERROR: An undefined dispersal kernel was specified: ', hlm_dispersal_kernel_mode
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select

   end subroutine ProbabilityDensity

   ! ====================================================================================

   real(r8) function PD_exponential(dist)
      
      use EDPftvarcon           , only : EDPftvarcon_inst
   
      ! Arguments
      real(r8), intent(in) :: dist
      ! real(r8), intent(in) :: param_a
      
      ! Assuming simple exponential decay.  In the future perhaps this could be an interface
      ! for different weight calculations (and could be held only in fates)
      
      PD_exponential = exp(-EDPftvarcon_inst%seed_dispersal_param_A(ipft)*dist)

   end function PD_exponential

   ! ====================================================================================

   real(r8) function PD_exppower(dist,param_a,param_b)
      
      use EDPftvarcon           , only : EDPftvarcon_inst
   
      ! Arguments
      real(r8), intent(in) :: dist
      ! real(r8), intent(in) :: param_a
      ! real(r8), intent(in) :: param_b
      
      associate(&
         param_a => EDPftvarcon_inst%seed_dispersal_param_A, &
         param_b => EDPftvarcon_inst%seed_dispersal_param_B)      
      
      PD_exppower = (param_b / (2*pi_const*gamma(2/param_b))) * &
                    exp(-(dist**param_b)/(param_a**param_b))
                    
      end associate

   end function PD_exppower

   ! ====================================================================================

   real(r8) function PD_logsech(dist)
      
      use EDPftvarcon           , only : EDPftvarcon_inst
   
      ! Arguments
      real(r8), intent(in) :: dist
      ! real(r8), intent(in) :: param_a
      ! real(r8), intent(in) :: param_b

      associate(&
         param_a => EDPftvarcon_inst%seed_dispersal_param_A, &
         param_b => EDPftvarcon_inst%seed_dispersal_param_B)      
      
      PD_logsech = (1/(pi_const**2 * param_b * dist**2)) / &
                    ((dist/param_a)**(1/param_b) + &
                    (dist/param_a)**(-1/param_b))
                    
      end associate

   end function PD_logsech

   ! ====================================================================================

      
end module FatesDispersalMod
