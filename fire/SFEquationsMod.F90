module SFEquationsMod
  ! ============================================================================
  ! Helper methods for the FATES SPITFIRE model
  ! Most equations come from:
  !  Thonicke et al. 2010, Biogeosciences 7:1991-2011
  !  Rothermel et al. 1972, Research Paper INT 115
  !  Albini et al. 1976, Research Report INT 30
  ! ============================================================================
  
  use FatesConstantsMod, only : r8 => fates_r8
  
  implicit none
  private
  
  public :: OptimumPackingRatio
  
  contains 
  
    real(r8) function OptimumPackingRatio(SAV)
    !
    !  DESCRIPTION:
    !  Calculates optimum packing ratio [unitless]
    !  Equation A6 in Thonicke et al. 2010
    !
    
    ! ARGUMENTS:
    real(r8), intent(in) :: SAV ! surface area to volume ratio of fuel [/cm]
    
    ! CONSTANTS:
    real(r8), parameter :: a = 0.200395_r8
    real(r8), parameter :: b = -0.8189_r8
          
    if (SAV < nearzero) then
      beta_op = 0.0_r8 
    else  
      beta_op = a*(SAV**b)
    end if
          
    end function OptimumPackingRatio
    
    !-------------------------------------------------------------------------------------
  
    real(r8) function MaximumReactionVelocity(SAV)
      !
      !  DESCRIPTION:
      !  Calculates maximum reaction velocity in /min
      !  From Equation 36 in Rothermel 1972; Fig. 12
      !
      
      ! ARGUMENTS:
      real(r8), intent(in) :: SAV ! fuel surface area to volume ratio [/cm]
      
      MaximumReactionVelocity = 1.0_r8/(0.0591_r8 + 2.926_r8*(SAV**(-1.5_r8)))
    
    end function MaximumReactionVelocity
   
    !---------------------------------------------------------------------------------------
  
    real(r8) function OptimumReactionVelocity(max_reaction_vel, SAV, beta_ratio)
      !
      !  DESCRIPTION:
      !  Calculates optimum reaction velocity in /min
      !  From Equation 38 in Rothermel 1972; Fig. 11
      !
      
      ! ARGUMENTS:
      real(r8), intent(in) :: max_reaction_vel ! maximum reaction velocity [/min]
      real(r8), intent(in) :: SAV              ! fuel surface area to volume ratio [/cm]
      real(r8), intent(in) :: beta_ratio       ! ratio of packing ratio to optimum packing ratio [0-1]
      
      ! LOCALS:
      real (r8) :: a, a_beta ! intermediate variables
      
      ! Equations in Table A1 Thonicke et al. 2010
      a = 8.9033_r8*(SAV**(-0.7913_r8))
      a_beta = exp(a*(1.0_r8 - beta_ratio)) 
      
      OptimumReactionVelocity = max_reaction_vel*(beta_ratio**a)*a_beta
      
    end function OptimumReactionVelocity
   
    !---------------------------------------------------------------------------------------
  
    real(r8) function MoistureCoefficient(moisture, MEF)
      !
      !  DESCRIPTION:
      !  Calculates the moisture dampening coefficient for reaction intensity
      !  Based on Equation in table A1 Thonicke et al. 2010. 
      !
      
      ! ARGUMENTS:
      real(r8), intent(in) :: moisture ! fuel moisture [m3/m3]
      real(r8), intent(in) :: MEF      ! fuel moisture of extinction [m3/m3]
      
      ! LOCALS:
      real(r8) :: mw_weight ! relative fuel moisture/fuel moisture of extinction
      
      ! average values for litter pools (dead leaves, twigs, small and large branches), plus grass
      mw_weight = moisture/MEF
      
      ! moist_damp is unitless
      MoistureCoefficient = max(0.0_r8, (1.0_r8 - (2.59_r8*mw_weight) +                   &
        (5.11_r8*(mw_weight**2.0_r8)) - (3.52_r8*(mw_weight**3.0_r8))))
        
      if (MoistureCoefficient > 1.0_r8) MoistureCoefficient = 1.0_r8
      
    end function MoistureCoefficient
   
    !---------------------------------------------------------------------------------------

    real(r8) function ReactionIntensity(fuel_loading, SAV,  beta_ratio, moisture, MEF)
      !
      !  DESCRIPTION:
      !  Calculates reaction intensity in kJ/m2/min
      !
      
      ! USES
      use SFParamsMod, only : SF_val_fuel_energy, SF_val_miner_damp
      
      ! ARGUMENTS:
      real(r8), intent(in) :: fuel_loading ! net fuel loading [kg/m2]
      real(r8), intent(in) :: SAV          ! fuel surface area to volume ratio [/cm]
      real(r8), intent(in) :: beta_ratio   ! ratio of packing ratio to optimum packing ratio [0-1]
      real(r8), intent(in) :: moisture     ! fuel moisture [m3/m3]
      real(r8), intent(in) :: MEF          ! fuel moisture of extinction [m3/m3]
      
      ! LOCALS:
      real(r8) :: max_reaction_vel ! maximum reaction velocity
      real(r8) :: opt_reaction_vel ! optimum reaction velocity
      real(r8) :: moist_coeff      ! moisture dampening coefficient [0-1]
      
      ! calculate maximum reaction velocity [/min]
      max_reaction_vel = MaximumReactionVelocity(SAV)
      
      ! calculate optimum reacion velocity [/min]
      opt_reaction_vel = OptimumReactionVelocity(max_reaction_vel, SAV, beta_ratio)
      
      ! calculate moisture dampening coefficient [0-1]
      moist_coeff = MoistureCoefficient(moisture, MEF)
      
    ReactionIntensity = opt_reaction_vel*fuel_loading*SF_val_fuel_energy*                &
      moist_coeff*SF_val_miner_damp

  end function ReactionIntensity
   
end module SFEquationsMod
