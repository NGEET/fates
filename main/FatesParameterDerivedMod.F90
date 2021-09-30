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
  use FatesInterfaceTypesMod,     only : nleafage
  use FatesInterfaceTypesMod,     only : ncrowndamage
  use FatesGlobals     ,     only : fates_log
  
  implicit none
  private

  type, public :: param_derived_type

     real(r8), allocatable :: jmax25top(:,:)  ! canopy top: maximum electron transport 
                                              ! rate at 25C (umol electrons/m**2/s)
     real(r8), allocatable :: tpu25top(:,:)   ! canopy top: triose phosphate utilization
                                              ! rate at 25C (umol CO2/m**2/s)
     real(r8), allocatable :: kp25top(:,:)    ! canopy top: initial slope of CO2 response
                                              ! curve (C4 plants) at 25C

     real(r8), allocatable :: branch_frac(:)  ! fraction of aboveground biomass in branches (as
                                              ! oppose to stems) - for use in damage allometries

     real(r8), allocatable :: damage_transitions(:,:,:) ! matrix of transition probabilities between
                                                      ! damage classes - one per PFT
     
   contains
     
     procedure :: Init
     procedure :: InitDamageTransitions
     procedure :: InitAllocate
     procedure :: InitAllocateDamageTransitions
     
  end type param_derived_type
  
  type(param_derived_type), public :: param_derived
  
contains 

  ! ===================================================================================
  subroutine InitAllocate(this,numpft)
    
    class(param_derived_type), intent(inout) :: this
    integer, intent(in)                      :: numpft
    
    allocate(this%jmax25top(numpft,nleafage))
    allocate(this%tpu25top(numpft,nleafage))
    allocate(this%kp25top(numpft,nleafage))

    allocate(this%branch_frac(numpft))
    
    
    return
  end subroutine InitAllocate

  ! =====================================================================================

 ! ===================================================================================
  subroutine InitAllocateDamageTransitions(this,ncrowndamage, numpft)
    
    class(param_derived_type), intent(inout) :: this
    integer, intent(in)                      :: ncrowndamage
    integer, intent(in)                      :: numpft

    allocate(this%damage_transitions(ncrowndamage,ncrowndamage, numpft))
    
    return
  end subroutine InitAllocateDamageTransitions

  ! =====================================================================================
 
  subroutine Init(this,numpft)

    use EDPftvarcon, only: EDPftvarcon_inst
    use SFParamsMod, only: SF_val_CWD_frac
    use FatesLitterMod, only : ncwd
    
    class(param_derived_type), intent(inout) :: this
    integer, intent(in)                      :: numpft
    
    ! local variables
    integer  :: ft                 ! pft index
    integer  :: iage               ! leaf age class index
    integer  :: c                  ! cwd index

    associate( vcmax25top => EDPftvarcon_inst%vcmax25top ) 
    
      call this%InitAllocate(numpft)
      
      do ft = 1,numpft
         
         do iage = 1, nleafage

            ! Parameters derived from vcmax25top. 
            ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
            ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of 
            ! Experimental Botany 44:907-920.  Here use a factor "1.67", from 
            ! Medlyn et al (2002) Plant, Cell and Environment 25:1167-1179
            
            ! RF - copied this from the CLM trunk code, but where did it come from, 
            ! and how can we make these consistant? 
            ! jmax25top(ft) =  &
            ! (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrzc),11._r8),35._r8)) * vcmax25top(ft)
            
            this%jmax25top(ft,iage) = 1.67_r8   * vcmax25top(ft,iage)
            this%tpu25top(ft,iage)  = 0.167_r8  * vcmax25top(ft,iage)
            this%kp25top(ft,iage)   = 20000._r8 * vcmax25top(ft,iage)
         
         end do

         ! Allocate fraction of biomass in branches
         this%branch_frac(ft) = sum(SF_val_CWD_frac(1:3))
         
      end do !ft

    end associate
    return
  end subroutine Init

!=========================================================================
  
  subroutine InitDamageTransitions(this, ncrowndamage, numpft)

    use EDPftvarcon, only: EDPftvarcon_inst


    class(param_derived_type), intent(inout) :: this
    integer, intent(in)                      :: ncrowndamage
    integer, intent(in)                      :: numpft

    ! local variables
    integer  :: ft                ! pft index
    integer  :: i                 ! crowndamage index
    real(r8) :: damage_frac       ! damage fraction 
   
   
    call this%InitAllocateDamageTransitions(ncrowndamage, numpft)
    
     do ft = 1, numpft

       damage_frac = EDPftvarcon_inst%damage_frac(ft)

       do i = 1, ncrowndamage

          ! zero the column
          this%damage_transitions(i,:,ft) = 0._r8
          ! 1 - damage rate stay the same
          this%damage_transitions(i,i,ft) = 1.0_r8 - damage_frac

          if(i < ncrowndamage) then
             ! fraction damaged get evenly split between higher damage classes
             this%damage_transitions(i,i+1:ncrowndamage,ft) = damage_frac/(ncrowndamage - i)
          end if
          ! Make sure it sums to one - they have to go somewhere
          this%damage_transitions(i, :, ft) = this%damage_transitions(i, :, ft)/SUM(this%damage_transitions(i, :, ft))
       end do

        write(fates_log(),'(a/,5(F12.6,1x))') 'annual transition matrix : ', this%damage_transitions(:,:,ft)
     end do



    
    return
  end subroutine InitDamageTransitions

end module FatesParameterDerivedMod
