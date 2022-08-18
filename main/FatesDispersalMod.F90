module FatesDispersalMod

   use shr_log_mod           , only : errMsg => shr_log_errMsg
   use FatesGlobals          , only : endrun => fates_endrun
   use FatesGlobals          , only : fates_log
   use FatesConstantsMod     , only : r8 => fates_r8
   use FatesConstantsMod     , only : pi_const

   implicit none
   
   ! Neighbor node
   type, public :: neighbor_type
      
      ! Grid cell neighbor
      type(neighbor_type), pointer :: next_neighbor => null() 
   
      integer  :: gindex       ! grid cell index
      real(r8) :: gc_dist      ! distance between source and neighbor
      real(r8) :: density_prob ! probability density from source
      
      contains
      
        procedure :: ProbabilityDensity
  
   end type neighbor_type

   ! Neighborhood linked list
   type, public :: neighborhood_type

      ! Linked list of neighbors for a given source grid cell
      type(neighbor_type), pointer :: first_neighbor => null()
      type(neighbor_type), pointer :: last_neighbor => null()
     
      integer  :: neighbor_count   ! total neighbors near source

   end type neighborhood_type

   type(neighborhood_type), public, pointer :: lneighbors(:)


contains

   ! ====================================================================================

   subroutine ProbabilityDensity(pd, dist)

      ! Main subroutine that calls difference routines based on case select mode

      ! Use
      use FatesInterfaceTypesMod, only : hlm_dispersal_kernel_exponential, &
                                         hlm_dispersal_kernel_exppower, &
                                         hlm_dispersal_kernel_logsech, &
                                         hlm_dispersal_kernel_mode

      ! Arguments
      real(r8), intent(out) :: pd   ! Probability density
      ! integer, intent(in)   :: ipft ! pft index - future arg
      real(r8), intent(in)  :: dist ! distance

      ! Local - temp
      param_a = 1._r8
      param_b = 1._r8
      hlm_dispersal_kernel_mode = 1
      hlm_dispersal_kernel_exponential = 1
      hlm_dispersal_kernel_exppower = 2
      hlm_dispersal_kernel_logsech = 3

      select case(hlm_dispersal_kernel_mode)

      case (hlm_dispersal_kernel_exponential)
         pd = PD_exponential(dist,param_a)
      case (hlm_dispersal_kernel_exppower)
         pd = PD_exppower(dist,param_a, param_b)
      case (hlm_dispersal_kernel_logsech)
         pd = PD_logsech(dist,param_a, param_b)
      case default
         write(fates_log(),*) 'ERROR: An undefined dispersal kernel was specified: ', hlm_dispersal_kernel_mode
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select

   end subroutine ProbabilityDensity

   ! ====================================================================================

   real(r8) function PD_exponential(dist,param_a)
      
      ! Arguments
      real(r8), intent(in) :: dist
      real(r8), intent(in) :: param_a
      
      ! Assuming simple exponential decay.  In the future perhaps this could be an interface
      ! for different weight calculations (and could be held only in fates)
      
      PD_exponential = exp(-param_a*dist)

   end function PD_exponential

   ! ====================================================================================

   real(r8) function PD_exppower(dist,param_a,param_b)
      
      ! Arguments
      real(r8), intent(in) :: dist
      real(r8), intent(in) :: param_a
      real(r8), intent(in) :: param_b
      
      ! Assuming simple exponential decay.  In the future perhaps this could be an interface
      ! for different weight calculations (and could be held only in fates)
      
      PD_exppower = (param_b / (2*pi_const*gamma(2/param_b))) * &
                    exp(-(dist**param_b)/(param_a**param_b))

   end function PD_exppower

   ! ====================================================================================

   real(r8) function PD_logsech(dist,param_a,param_b)
      
      ! Arguments
      real(r8), intent(in) :: dist
      real(r8), intent(in) :: param_a
      real(r8), intent(in) :: param_b
      
      ! Assuming simple exponential decay.  In the future perhaps this could be an interface
      ! for different weight calculations (and could be held only in fates)
      
      PD_logsech = (1/(pi_const**2 * param_b * dist**2)) / &
                    ((dist/param_a)**(1/param_b) + &
                    (dist/param_a)**(-1/param_b))

   end function PD_logsech

   ! ====================================================================================

      
end module FatesDispersalMod