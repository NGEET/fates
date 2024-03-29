module FatesBstressMod
   
   !-------------------------------------------------------------------------------------
   ! Description: calculate the stress impact on transpiration from salinity and sulphide in soils
   !              note that water stress is calculated in EDBtranMod or HYDRO
   ! ------------------------------------------------------------------------------------
   !
   use FatesConstantsMod , only : tfrz => t_water_freeze_k_1atm 
   use FatesConstantsMod , only : itrue,ifalse
   use EDParamsMod,        only : maxpft
   use EDTypesMod        , only : ed_site_type
   use FatesPatchMod,      only : fates_patch_type
   use FatesCohortMod    , only : fates_cohort_type
   use shr_kind_mod      , only : r8 => shr_kind_r8
   use FatesInterfaceTypesMod , only : bc_in_type, &
                                  bc_out_type, &
                                  numpft
   use FatesInterfaceTypesMod , only : hlm_use_planthydro
   use FatesGlobals      , only : fates_log
   use EDBtranMod        , only : check_layer_water
   use FatesAllometryMod , only : set_root_fraction

   implicit none
   private
   
   public :: btran_sal_stress_fates

contains
 ! =====================================================================================

  subroutine btran_sal_stress_fates( nsites, sites, bc_in)

      
      ! ---------------------------------------------------------------------------------
      ! Calculate the transpiration wetness function (BTRAN) and the root uptake
      ! distribution (ROOTR).
      ! Boundary conditions in: bc_in(s)%salinity_sl(j)        salinity concontration[ppm]
      ! Boundary conditions in: bc_in(s)%sulphide_sl(j)        sulphide concontration[ppm]
      ! Output cpatch%bstress_sal_ft(ft)
      ! Output cpatch%bstress_sul_ft(ft)
      ! ---------------------------------------------------------------------------------
      
      ! Arguments
      
      integer,intent(in)                      :: nsites
      type(ed_site_type),intent(inout),target :: sites(nsites)
      type(bc_in_type),intent(in)             :: bc_in(nsites)
     
      !
      ! !LOCAL VARIABLES:
      type(fates_patch_type),pointer             :: cpatch ! Current Patch Pointer
      type(fates_cohort_type),pointer            :: ccohort ! Current cohort pointer
      integer  :: s                 ! site
      integer  :: j                 ! soil layer
      integer  :: ft                ! plant functional type index
      real(r8) :: salinity_node     ! salinity in the soil water [ppt]
      real(r8) :: rresis            ! salinity limitation to transpiration independent

      !------------------------------------------------------------------------------
        
        do s = 1,nsites

           cpatch => sites(s)%oldest_patch
           do while (associated(cpatch))                 
              
              do ft = 1,numpft
                 cpatch%bstress_sal_ft(ft) = 0.0_r8

                 call set_root_fraction(sites(s)%rootfrac_scr, ft, &
                      sites(s)%zi_soil, &
                      bc_in(s)%max_rooting_depth_index_col )

                 do j = 1,bc_in(s)%nlevsoil
                    
                    ! Calculations are only relevant where liquid water exists
                    ! see clm_fates%wrap_btran for calculation with CLM/ELM
                    
                    if ( check_layer_water(bc_in(s)%h2o_liqvol_sl(j),bc_in(s)%tempk_sl(j)) )  then
                       
                       salinity_node =  bc_in(s)%salinity_sl(j)
                       
                       rresis  = min( 1.244_r8/(1+exp((0.186_r8-salinity_node)/(-0.132_r8))), 1._r8)
                       
                       cpatch%bstress_sal_ft(ft) = cpatch%bstress_sal_ft(ft)+sites(s)%rootfrac_scr(j)*rresis

                    end if
                    
                 end do !j
                 
              end do !PFT              
                 
              cpatch => cpatch%younger
	      
           end do
          
        end do
           
        return
    end subroutine btran_sal_stress_fates
          
  ! ====================================================================================

end module FatesBstressMod
