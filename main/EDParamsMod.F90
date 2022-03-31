module EDParamsMod

   !
   ! module that deals with reading the ED parameter file
   !

   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : nearzero
   use FatesParametersInterface, only : param_string_length
   use FatesGlobals        , only : fates_log
   use FatesGlobals        , only : endrun => fates_endrun

   ! CIME Globals
   use shr_log_mod         , only : errMsg => shr_log_errMsg

   implicit none
   private
   save

   !
   ! this is what the user can use for the actual values
   !
   
   real(r8),protected, public :: fates_mortality_disturbance_fraction ! the fraction of canopy mortality that results in disturbance
   real(r8),protected, public :: ED_val_comp_excln
   real(r8),protected, public :: ED_val_init_litter
   real(r8),protected, public :: ED_val_nignitions
   real(r8),protected, public :: ED_val_understorey_death
   real(r8),protected, public :: ED_val_cwd_fcel
   real(r8),protected, public :: ED_val_cwd_flig
   real(r8),protected, public :: ED_val_base_mr_20
   real(r8),protected, public :: ED_val_phen_drought_threshold
   real(r8),protected, public :: ED_val_phen_doff_time
   real(r8),protected, public :: ED_val_phen_a
   real(r8),protected, public :: ED_val_phen_b
   real(r8),protected, public :: ED_val_phen_c
   real(r8),protected, public :: ED_val_phen_chiltemp
   real(r8),protected, public :: ED_val_phen_mindayson
   real(r8),protected, public :: ED_val_phen_ncolddayslim
   real(r8),protected, public :: ED_val_phen_coldtemp
   real(r8),protected, public :: ED_val_cohort_size_fusion_tol
   real(r8),protected, public :: ED_val_cohort_age_fusion_tol
   real(r8),protected, public :: ED_val_patch_fusion_tol
   real(r8),protected, public :: ED_val_canopy_closure_thresh ! site-level canopy closure point where trees take on forest (narrow) versus savannah (wide) crown allometry
   integer,protected, public  :: stomatal_model  !switch for choosing between stomatal conductance models, 1 for Ball-Berry, 2 for Medlyn

   
   logical,protected, public :: active_crown_fire        ! flag, 1=active crown fire 0=no active crown fire
   character(len=param_string_length),parameter :: fates_name_active_crown_fire = "fates_fire_active_crown_fire"

   real(r8), protected, public :: cg_strikes             ! fraction of cloud to ground lightning strikes (0-1)
   character(len=param_string_length),parameter :: fates_name_cg_strikes="fates_fire_cg_strikes"
   
   real(r8),protected,public  :: q10_mr     ! Q10 for respiration rate (for soil fragmenation and plant respiration)    (unitless)
   real(r8),protected,public  :: q10_froz   ! Q10 for frozen-soil respiration rates (for soil fragmentation)            (unitless)

   ! two special parameters whose size is defined in the parameter file
   real(r8),protected,allocatable,public :: ED_val_history_sizeclass_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_ageclass_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_height_bin_edges(:)
   real(r8),protected,allocatable,public :: ED_val_history_coageclass_bin_edges(:)
   
   character(len=param_string_length),parameter,public :: ED_name_mort_disturb_frac = "fates_mort_disturb_frac"
   character(len=param_string_length),parameter,public :: ED_name_comp_excln = "fates_comp_excln"
   character(len=param_string_length),parameter,public :: ED_name_init_litter = "fates_init_litter"
   character(len=param_string_length),parameter,public :: ED_name_nignitions = "fates_fire_nignitions"
   character(len=param_string_length),parameter,public :: ED_name_understorey_death = "fates_mort_understorey_death"
   character(len=param_string_length),parameter,public :: ED_name_cwd_fcel= "fates_cwd_fcel"   
   character(len=param_string_length),parameter,public :: ED_name_cwd_flig= "fates_cwd_flig"   
   character(len=param_string_length),parameter,public :: ED_name_base_mr_20= "fates_base_mr_20"   
   character(len=param_string_length),parameter,public :: ED_name_phen_drought_threshold= "fates_phen_drought_threshold"   
   character(len=param_string_length),parameter,public :: ED_name_phen_doff_time= "fates_phen_doff_time"   
   character(len=param_string_length),parameter,public :: ED_name_phen_a= "fates_phen_a"   
   character(len=param_string_length),parameter,public :: ED_name_phen_b= "fates_phen_b"   
   character(len=param_string_length),parameter,public :: ED_name_phen_c= "fates_phen_c"   
   character(len=param_string_length),parameter,public :: ED_name_phen_chiltemp= "fates_phen_chiltemp"   
   character(len=param_string_length),parameter,public :: ED_name_phen_mindayson= "fates_phen_mindayson"   
   character(len=param_string_length),parameter,public :: ED_name_phen_ncolddayslim= "fates_phen_ncolddayslim"   
   character(len=param_string_length),parameter,public :: ED_name_phen_coldtemp= "fates_phen_coldtemp"   
   character(len=param_string_length),parameter,public :: ED_name_cohort_size_fusion_tol= "fates_cohort_size_fusion_tol"
   character(len=param_string_length),parameter,public :: ED_name_cohort_age_fusion_tol = "fates_cohort_age_fusion_tol"
   character(len=param_string_length),parameter,public :: ED_name_patch_fusion_tol= "fates_patch_fusion_tol"
   character(len=param_string_length),parameter,public :: ED_name_canopy_closure_thresh= "fates_canopy_closure_thresh"      
   character(len=param_string_length),parameter,public :: ED_name_stomatal_model= "fates_leaf_stomatal_model"

   ! Resistance to active crown fire
  

   character(len=param_string_length),parameter :: fates_name_q10_mr="fates_q10_mr"
   character(len=param_string_length),parameter :: fates_name_q10_froz="fates_q10_froz"


   ! non-scalar parameter names
   character(len=param_string_length),parameter,public :: ED_name_history_sizeclass_bin_edges= "fates_history_sizeclass_bin_edges"      
   character(len=param_string_length),parameter,public :: ED_name_history_ageclass_bin_edges= "fates_history_ageclass_bin_edges"      
   character(len=param_string_length),parameter,public :: ED_name_history_height_bin_edges= "fates_history_height_bin_edges"
   character(len=param_string_length),parameter,public :: ED_name_history_coageclass_bin_edges = "fates_history_coageclass_bin_edges"

   ! Hydraulics Control Parameters (ONLY RELEVANT WHEN USE_FATES_HYDR = TRUE)
   ! ----------------------------------------------------------------------------------------------
   real(r8),protected,public :: hydr_kmax_rsurf1         !  maximum conducitivity for unit root surface 
                                                  !  soil to root direction (kg water/m2 root area/Mpa/s)
   character(len=param_string_length),parameter,public :: hydr_name_kmax_rsurf1 = "fates_hydr_kmax_rsurf1"  
   
   real(r8),protected,public :: hydr_kmax_rsurf2         !  maximum conducitivity for unit root surface 
                                                  !  root to soil direciton (kg water/m2 root area/Mpa/s)
   character(len=param_string_length),parameter,public :: hydr_name_kmax_rsurf2 = "fates_hydr_kmax_rsurf2" 

   real(r8),protected,public :: hydr_psi0          !  sapwood water potential at saturation (MPa)
   character(len=param_string_length),parameter,public :: hydr_name_psi0 = "fates_hydr_psi0"

   real(r8),protected,public :: hydr_psicap        !  sapwood water potential at which capillary reserves exhausted (MPa)
   character(len=param_string_length),parameter,public :: hydr_name_psicap = "fates_hydr_psicap"

   !Soil BGC parameters, mostly used for testing FATES when not coupled to the dynamics bgc hlm
   ! ----------------------------------------------------------------------------------------------
   real(r8),protected,public :: bgc_soil_salinity ! site-level soil salinity for FATES when not coupled to dynamic soil BGC of salinity
   character(len=param_string_length),parameter,public :: bgc_name_soil_salinity= "fates_soil_salinity"      
   
   ! Logging Control Parameters (ONLY RELEVANT WHEN USE_FATES_LOGGING = TRUE)
   ! ----------------------------------------------------------------------------------------------

   real(r8),protected,public :: logging_dbhmin              ! Minimum dbh at which logging is applied (cm)
                                                            ! Typically associated with harvesting
   character(len=param_string_length),parameter,public :: logging_name_dbhmin = "fates_logging_dbhmin"

   real(r8),protected,public :: logging_dbhmax              ! Maximum dbh at which logging is applied (cm)
                                                            ! Typically associated with fire suppression
                                                            ! (THIS PARAMETER IS NOT USED YET)
   character(len=param_string_length),parameter,public :: logging_name_dbhmax = "fates_logging_dbhmax"


   real(r8),protected,public :: logging_collateral_frac     ! Ratio of collateral mortality to direct logging mortality
   character(len=param_string_length),parameter,public :: logging_name_collateral_frac = "fates_logging_collateral_frac"

   real(r8),protected,public :: logging_coll_under_frac ! Fraction of understory plants that die when logging disturbance
                                                 ! is generated
   character(len=param_string_length),parameter,public :: logging_name_coll_under_frac = "fates_logging_coll_under_frac"
   
   real(r8),protected,public :: logging_direct_frac         ! Fraction of stems logged per event
   character(len=param_string_length),parameter,public :: logging_name_direct_frac = "fates_logging_direct_frac"

   real(r8),protected,public :: logging_mechanical_frac         ! Fraction of stems logged per event
   character(len=param_string_length),parameter,public :: logging_name_mechanical_frac = "fates_logging_mechanical_frac"

   real(r8),protected,public :: logging_event_code          ! Code that options how logging events are structured 
   character(len=param_string_length),parameter,public :: logging_name_event_code = "fates_logging_event_code"
   
   real(r8),protected,public :: logging_dbhmax_infra        ! "Tree diameter, above which infrastructure from logging does not impact damage or mortality.
   character(len=param_string_length),parameter,public :: logging_name_dbhmax_infra = "fates_logging_dbhmax_infra"
   
   real(r8),protected,public :: logging_export_frac        ! "fraction of trunk product being shipped offsite, the 
                                                    ! leftovers will be left onsite as large CWD
   character(len=param_string_length),parameter,public :: logging_name_export_frac ="fates_logging_export_frac"   

   public :: FatesParamsInit
   public :: FatesRegisterParams
   public :: FatesReceiveParams
   public :: FatesReportParams
  
contains

  !-----------------------------------------------------------------------
  subroutine FatesParamsInit()
    ! Initialize all parameters to nan to ensure that we get valid
    ! values back from the host.
    
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)

    implicit none

    fates_mortality_disturbance_fraction  = nan
    ED_val_comp_excln                     = nan
    ED_val_init_litter                    = nan
    ED_val_nignitions                     = nan
    ED_val_understorey_death              = nan
    ED_val_cwd_fcel                       = nan
    ED_val_cwd_flig                       = nan
    ED_val_base_mr_20                     = nan
    ED_val_phen_drought_threshold         = nan
    ED_val_phen_doff_time                 = nan
    ED_val_phen_a                         = nan
    ED_val_phen_b                         = nan
    ED_val_phen_c                         = nan
    ED_val_phen_chiltemp                  = nan
    ED_val_phen_mindayson                 = nan
    ED_val_phen_ncolddayslim              = nan
    ED_val_phen_coldtemp                  = nan
    ED_val_cohort_size_fusion_tol         = nan
    ED_val_cohort_age_fusion_tol          = nan
    ED_val_patch_fusion_tol               = nan
    ED_val_canopy_closure_thresh          = nan
    stomatal_model                        = -9
    hydr_kmax_rsurf1                      = nan
    hydr_kmax_rsurf2                      = nan

    hydr_psi0                             = nan
    hydr_psicap                           = nan
    
    bgc_soil_salinity                     = nan

    logging_dbhmin                        = nan
    logging_dbhmax                        = nan
    logging_collateral_frac               = nan
    logging_direct_frac                   = nan
    logging_mechanical_frac               = nan
    logging_event_code                    = nan
    logging_dbhmax_infra                  = nan
    logging_export_frac                   = nan
    q10_mr                                = nan
    q10_froz                              = nan

  end subroutine FatesParamsInit

  !-----------------------------------------------------------------------
  subroutine FatesRegisterParams(fates_params)
    ! Register the parameters we want the host to provide, and
    ! indicate whether they are fates parameters or host parameters
    ! that need to be synced with host values.

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar, dimension_shape_1d
    use FatesParametersInterface, only : dimension_name_history_size_bins, dimension_name_history_age_bins
    use FatesParametersInterface, only : dimension_name_history_height_bins
    use FatesParametersInterface, only : dimension_name_history_coage_bins
    use FatesParametersInterface, only : dimension_shape_scalar


    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names_scalar(1) = (/dimension_name_scalar/)
    character(len=param_string_length), parameter :: dim_names_sizeclass(1) = (/dimension_name_history_size_bins/)
    character(len=param_string_length), parameter :: dim_names_ageclass(1) = (/dimension_name_history_age_bins/)
    character(len=param_string_length), parameter :: dim_names_height(1) = (/dimension_name_history_height_bins/)
    character(len=param_string_length), parameter :: dim_names_coageclass(1) = (/dimension_name_history_coage_bins/)

       
    call FatesParamsInit()

    call fates_params%RegisterParameter(name=ED_name_mort_disturb_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_comp_excln, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_init_litter, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_nignitions, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_understorey_death, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_cwd_fcel, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_cwd_flig, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_base_mr_20, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_drought_threshold, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_doff_time, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_a, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_b, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_c, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_chiltemp, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_mindayson, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_ncolddayslim, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_phen_coldtemp, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_cohort_size_fusion_tol, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_cohort_age_fusion_tol, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_patch_fusion_tol, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_canopy_closure_thresh, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=ED_name_stomatal_model, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
	 
    call fates_params%RegisterParameter(name=hydr_name_kmax_rsurf1, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=hydr_name_kmax_rsurf2, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=hydr_name_psi0, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=hydr_name_psicap, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=bgc_name_soil_salinity, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar) 

    call fates_params%RegisterParameter(name=logging_name_dbhmin, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_dbhmax, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_collateral_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_coll_under_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_direct_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_mechanical_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_event_code, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_dbhmax_infra, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=logging_name_export_frac, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=fates_name_q10_mr, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=fates_name_q10_froz, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    ! non-scalar parameters
    call fates_params%RegisterParameter(name=ED_name_history_sizeclass_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_sizeclass)

    call fates_params%RegisterParameter(name=ED_name_history_ageclass_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_ageclass)

    call fates_params%RegisterParameter(name=ED_name_history_height_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_height)

    call fates_params%RegisterParameter(name=fates_name_active_crown_fire, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=fates_name_cg_strikes, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)
    
    call fates_params%RegisterParameter(name=ED_name_history_coageclass_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_coageclass)


  end subroutine FatesRegisterParams

  
  !-----------------------------------------------------------------------
  subroutine FatesReceiveParams(fates_params)
    
    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    real(r8) :: tmpreal ! local real variable for changing type on read
    
    call fates_params%RetreiveParameter(name=ED_name_mort_disturb_frac, &
          data=fates_mortality_disturbance_fraction)

    call fates_params%RetreiveParameter(name=ED_name_comp_excln, &
         data=ED_val_comp_excln)

    call fates_params%RetreiveParameter(name=ED_name_init_litter, &
         data=ED_val_init_litter)

    call fates_params%RetreiveParameter(name=ED_name_nignitions, &
         data=ED_val_nignitions)

    call fates_params%RetreiveParameter(name=ED_name_understorey_death, &
         data=ED_val_understorey_death)

    call fates_params%RetreiveParameter(name=ED_name_cwd_fcel, &
         data=ED_val_cwd_fcel)

    call fates_params%RetreiveParameter(name=ED_name_cwd_flig, &
         data=ED_val_cwd_flig)

    call fates_params%RetreiveParameter(name=ED_name_base_mr_20, &
         data=ED_val_base_mr_20)

    call fates_params%RetreiveParameter(name=ED_name_phen_drought_threshold, &
         data=ED_val_phen_drought_threshold)

    call fates_params%RetreiveParameter(name=ED_name_phen_doff_time, &
         data=ED_val_phen_doff_time)

    call fates_params%RetreiveParameter(name=ED_name_phen_a, &
         data=ED_val_phen_a)

    call fates_params%RetreiveParameter(name=ED_name_phen_b, &
         data=ED_val_phen_b)

    call fates_params%RetreiveParameter(name=ED_name_phen_c, &
         data=ED_val_phen_c)

    call fates_params%RetreiveParameter(name=ED_name_phen_chiltemp, &
         data=ED_val_phen_chiltemp)

    call fates_params%RetreiveParameter(name=ED_name_phen_mindayson, &
         data=ED_val_phen_mindayson)

    call fates_params%RetreiveParameter(name=ED_name_phen_ncolddayslim, &
         data=ED_val_phen_ncolddayslim)

    call fates_params%RetreiveParameter(name=ED_name_phen_coldtemp, &
         data=ED_val_phen_coldtemp)

    call fates_params%RetreiveParameter(name=ED_name_cohort_size_fusion_tol, &
         data=ED_val_cohort_size_fusion_tol)

    call fates_params%RetreiveParameter(name=ED_name_cohort_age_fusion_tol, &
         data=ED_val_cohort_age_fusion_tol)

    call fates_params%RetreiveParameter(name=ED_name_patch_fusion_tol, &
         data=ED_val_patch_fusion_tol)
    
    call fates_params%RetreiveParameter(name=ED_name_canopy_closure_thresh, &
         data=ED_val_canopy_closure_thresh)

    call fates_params%RetreiveParameter(name=ED_name_stomatal_model, &
         data=tmpreal)
    stomatal_model = nint(tmpreal)
    
    call fates_params%RetreiveParameter(name=hydr_name_kmax_rsurf1, &
          data=hydr_kmax_rsurf1)

    call fates_params%RetreiveParameter(name=hydr_name_kmax_rsurf2, &
          data=hydr_kmax_rsurf2)	 
    
    call fates_params%RetreiveParameter(name=hydr_name_psi0, &
          data=hydr_psi0)

    call fates_params%RetreiveParameter(name=hydr_name_psicap, &
          data=hydr_psicap)
	  
    call fates_params%RetreiveParameter(name=bgc_name_soil_salinity, &
          data=bgc_soil_salinity)	  

    call fates_params%RetreiveParameter(name=logging_name_dbhmin, &
          data=logging_dbhmin)

    call fates_params%RetreiveParameter(name=logging_name_dbhmax, &
          data=logging_dbhmax)

    call fates_params%RetreiveParameter(name=logging_name_collateral_frac, &
          data=logging_collateral_frac)
    
    call fates_params%RetreiveParameter(name=logging_name_coll_under_frac, &
          data=logging_coll_under_frac)

    call fates_params%RetreiveParameter(name=logging_name_direct_frac, &
          data=logging_direct_frac)

    call fates_params%RetreiveParameter(name=logging_name_mechanical_frac, &
          data=logging_mechanical_frac)
    
    call fates_params%RetreiveParameter(name=logging_name_event_code, &
          data=logging_event_code)

    call fates_params%RetreiveParameter(name=logging_name_dbhmax_infra, &
          data=logging_dbhmax_infra)

    call fates_params%RetreiveParameter(name=logging_name_export_frac, &
          data=logging_export_frac)

    call fates_params%RetreiveParameter(name=fates_name_q10_mr, &
          data=q10_mr)
    
    call fates_params%RetreiveParameter(name=fates_name_q10_froz, &
          data=q10_froz)

    call fates_params%RetreiveParameter(name=fates_name_active_crown_fire, & 
          data=tmpreal)
    active_crown_fire = (abs(tmpreal-1.0_r8)<nearzero)

    call fates_params%RetreiveParameter(name=fates_name_cg_strikes, &
          data=cg_strikes)

    ! parameters that are arrays of size defined within the params file and thus need allocating as well
    call fates_params%RetreiveParameterAllocate(name=ED_name_history_sizeclass_bin_edges, &
          data=ED_val_history_sizeclass_bin_edges)

    call fates_params%RetreiveParameterAllocate(name=ED_name_history_ageclass_bin_edges, &
          data=ED_val_history_ageclass_bin_edges)

    call fates_params%RetreiveParameterAllocate(name=ED_name_history_height_bin_edges, &
          data=ED_val_history_height_bin_edges)

    call fates_params%RetreiveParameterAllocate(name=ED_name_history_coageclass_bin_edges, &
         data=ED_val_history_coageclass_bin_edges)


  end subroutine FatesReceiveParams
  
  ! =====================================================================================

  subroutine FatesReportParams(is_master)

     logical,intent(in) :: is_master

     character(len=32),parameter :: fmt0 = '(a,(F12.4))'
     logical, parameter :: debug_report = .false.
     
     if(debug_report .and. is_master) then
        
        write(fates_log(),*) '-----------  FATES Scalar Parameters -----------------'
        write(fates_log(),fmt0) 'fates_mortality_disturbance_fraction = ',fates_mortality_disturbance_fraction
        write(fates_log(),fmt0) 'ED_val_comp_excln = ',ED_val_comp_excln
        write(fates_log(),fmt0) 'ED_val_init_litter = ',ED_val_init_litter
        write(fates_log(),fmt0) 'ED_val_nignitions = ',ED_val_nignitions
        write(fates_log(),fmt0) 'ED_val_understorey_death = ',ED_val_understorey_death
        write(fates_log(),fmt0) 'ED_val_cwd_fcel = ',ED_val_cwd_fcel
        write(fates_log(),fmt0) 'ED_val_cwd_flig = ',ED_val_cwd_flig
        write(fates_log(),fmt0) 'ED_val_base_mr_20 = ', ED_val_base_mr_20
        write(fates_log(),fmt0) 'ED_val_phen_drought_threshold = ',ED_val_phen_drought_threshold
        write(fates_log(),fmt0) 'ED_val_phen_doff_time = ',ED_val_phen_doff_time
        write(fates_log(),fmt0) 'ED_val_phen_a = ',ED_val_phen_a
        write(fates_log(),fmt0) 'ED_val_phen_b = ',ED_val_phen_b
        write(fates_log(),fmt0) 'ED_val_phen_c = ',ED_val_phen_c
        write(fates_log(),fmt0) 'ED_val_phen_chiltemp = ',ED_val_phen_chiltemp
        write(fates_log(),fmt0) 'ED_val_phen_mindayson = ',ED_val_phen_mindayson
        write(fates_log(),fmt0) 'ED_val_phen_ncolddayslim = ',ED_val_phen_ncolddayslim
        write(fates_log(),fmt0) 'ED_val_phen_coldtemp = ',ED_val_phen_coldtemp
        write(fates_log(),fmt0) 'ED_val_cohort_size_fusion_tol = ',ED_val_cohort_size_fusion_tol
        write(fates_log(),fmt0) 'ED_val_cohort_age_fusion_tol = ',ED_val_cohort_age_fusion_tol
        write(fates_log(),fmt0) 'ED_val_patch_fusion_tol = ',ED_val_patch_fusion_tol
        write(fates_log(),fmt0) 'ED_val_canopy_closure_thresh = ',ED_val_canopy_closure_thresh
        write(fates_log(),fmt0) 'stomatal_model = ',stomatal_model
        write(fates_log(),fmt0) 'hydr_kmax_rsurf1 = ',hydr_kmax_rsurf1
        write(fates_log(),fmt0) 'hydr_kmax_rsurf2 = ',hydr_kmax_rsurf2  
        write(fates_log(),fmt0) 'hydr_psi0 = ',hydr_psi0
        write(fates_log(),fmt0) 'hydr_psicap = ',hydr_psicap
        write(fates_log(),fmt0) 'bgc_soil_salinity = ', bgc_soil_salinity
        write(fates_log(),fmt0) 'logging_dbhmin = ',logging_dbhmin
        write(fates_log(),fmt0) 'logging_dbhmax = ',logging_dbhmax
        write(fates_log(),fmt0) 'logging_collateral_frac = ',logging_collateral_frac
        write(fates_log(),fmt0) 'logging_coll_under_frac = ',logging_coll_under_frac
        write(fates_log(),fmt0) 'logging_direct_frac = ',logging_direct_frac
        write(fates_log(),fmt0) 'logging_mechanical_frac = ',logging_mechanical_frac
        write(fates_log(),fmt0) 'logging_event_code = ',logging_event_code
        write(fates_log(),fmt0) 'logging_dbhmax_infra = ',logging_dbhmax_infra
        write(fates_log(),fmt0) 'q10_mr = ',q10_mr
        write(fates_log(),fmt0) 'q10_froz = ',q10_froz
        write(fates_log(),fmt0) 'cg_strikes = ',cg_strikes
        write(fates_log(),'(a,L)') 'active_crown_fire = ',active_crown_fire
        write(fates_log(),*) '------------------------------------------------------'

     end if

  end subroutine FatesReportParams

  
end module EDParamsMod
