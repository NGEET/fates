module FatesParameterDerivedMod

  ! -------------------------------------------------------------------------------------
  ! This module contains all procedures types and settings for any quantities that are
  ! statically derived from static model parameters.  These are unchanging quantities
  ! and are based off of simple relationships from parameters that the user can
  ! vary.  This should be called once, and early in the model initialization call
  ! sequence immediately after FATES parameters are read in.
  !
  ! -------------------------------------------------------------------------------------

  use FatesConstantsMod,     only : r8 => fates_r8
  use FatesConstantsMod,     only : umolC_to_kgC
  use FatesConstantsMod,     only : g_per_kg
  
  type param_derived_type

     real(r8), allocatable :: vcmax25top(:) ! canopy top: maximum rate of carboxylation
                                            ! at 25C (umol CO2/m**2/s)
     real(r8), allocatable :: jmax25top(:)  ! canopy top: maximum electron transport 
                                            ! rate at 25C (umol electrons/m**2/s)
     real(r8), allocatable :: tpu25top(:)   ! canopy top: triose phosphate utilization
                                            ! rate at 25C (umol CO2/m**2/s)
     real(r8), allocatable :: kp25top(:)    ! canopy top: initial slope of CO2 response
                                            ! curve (C4 plants) at 25C
   contains
     
     procedure :: Init
     procedure :: InitAllocate
     
  end type param_derived_type
  
  type(param_derived_type) :: param_derived
  
contains
  
  subroutine InitAllocate(this,numpft)
    
    class(param_derived_type), intent(inout) :: this
    integer, intent(in)                      :: numpft
    
    allocate(this%vcmax25top(numpft))
    allocate(this%jmax25top(numpft))
    allocate(this%tpu25top(numpft))
    allocate(this%kp25top(numpft))
    
    return
  end subroutine InitAllocate

  ! =====================================================================================
  
  subroutine Init(this,numpft)

    use EDPftvarcon, only: EDPftvarcon_inst

    class(param_derived_type), intent(inout) :: this
    integer, intent(in)                      :: numpft
    
    ! local variables
    integer  :: ft                 ! pft index
    real(r8) :: lnc                ! leaf N concentration (gN leaf/m^2)
    
    associate( vcmax25top => EDPftvarcon_inst%vcmax25top ) 
    
      call this%InitAllocate(numpft)
      
      do ft = 1,numpft
         
         ! Parameters derived from vcmax25top. 
         ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
         ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of 
         ! Experimental Botany 44:907-920.  Here use a factor "1.67", from 
         ! Medlyn et al (2002) Plant, Cell and Environment 25:1167-1179
         
         ! RF - copied this from the CLM trunk code, but where did it come from, 
         ! and how can we make these consistant? 
         ! jmax25top(ft) =  &
         ! (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrzc),11._r8),35._r8)) * vcmax25top(ft)
         
         this%jmax25top(ft) = 1.67_r8   * vcmax25top(ft)
         this%tpu25top(ft)  = 0.167_r8  * vcmax25top(ft)
         this%kp25top(ft)   = 20000._r8 * vcmax25top(ft)
         
         
         
      end do !ft 

    end associate
    return
  end subroutine Init
  

end module FatesParameterDerivedMod
