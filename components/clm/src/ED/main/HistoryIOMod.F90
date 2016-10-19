Module HistoryIOMod

  
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : fates_avg_flag_length, fates_short_string_length, fates_long_string_length
  use FatesGlobals    , only : fates_log
  use EDTypesMod      , only : cp_hio_ignore_val
  use pftconMod       , only : pftcon

  implicit none

  ! These variables hold the index of the history output structure so we don't
  ! have to constantly do name lookup when we want to populate the dataset
  ! These indices are set during "define_history_vars()" call to "set_history_var()"
  ! during the initialize phase.  Definitions are not provide, for an explanation of
  ! the variable go to its registry.  (IH_ signifies "index history")
  
  ! Indices to 1D Patch variables

  integer, private :: ih_trimming_pa
  integer, private :: ih_area_plant_pa
  integer, private :: ih_area_treespread_pa
  integer, private :: ih_canopy_spread_pa
  integer, private :: ih_nesterov_fire_danger_pa
  integer, private :: ih_spitfire_ROS_pa
  integer, private :: ih_effect_wspeed_pa
  integer, private :: ih_TFC_ROS_pa
  integer, private :: ih_fire_intensity_pa
  integer, private :: ih_fire_area_pa
  integer, private :: ih_scorch_height_pa
  integer, private :: ih_fire_fuel_bulkd_pa
  integer, private :: ih_fire_fuel_eff_moist_pa
  integer, private :: ih_fire_fuel_sav_pa
  integer, private :: ih_fire_fuel_mef_pa
  integer, private :: ih_sum_fuel_pa
  integer, private :: ih_litter_in_pa
  integer, private :: ih_litter_out_pa

  integer, private :: ih_efpot_pa        ! NA
  integer, private :: ih_rb_pa           ! NA

  integer, private :: ih_daily_temp
  integer, private :: ih_daily_rh
  integer, private :: ih_daily_prec
  integer, private :: ih_seed_bank_si
  integer, private :: ih_seeds_in_pa
  integer, private :: ih_seed_decay_pa
  integer, private :: ih_seed_germination_pa
  integer, private :: ih_bstore_pa
  integer, private :: ih_bdead_pa
  integer, private :: ih_balive_pa
  integer, private :: ih_bleaf_pa
  integer, private :: ih_btotal_pa
  integer, private :: ih_npp_pa
  integer, private :: ih_gpp_pa
  integer, private :: ih_aresp_pa
  integer, private :: ih_maint_resp_pa
  integer, private :: ih_growth_resp_pa
  
  ! Indices to (patch x pft) variables   (using nlevgrnd as surrogate)

  integer, private :: ih_biomass_pa_pft
  integer, private :: ih_leafbiomass_pa_pft
  integer, private :: ih_storebiomass_pa_pft
  integer, private :: ih_nindivs_pa_pft

  ! Indices to (site) variables


  integer, private :: ih_nep_si
  integer, private :: ih_nep_timeintegrated_si
  integer, private :: ih_npp_timeintegrated_si
  integer, private :: ih_hr_timeintegrated_si
  integer, private :: ih_nbp_si
  integer, private :: ih_npp_si
  integer, private :: ih_fire_c_to_atm_si
  integer, private :: ih_ed_to_bgc_this_edts_si
  integer, private :: ih_ed_to_bgc_last_edts_si
  integer, private :: ih_totecosysc_si
  integer, private :: ih_totecosysc_old_si
  integer, private :: ih_totedc_si
  integer, private :: ih_totedc_old_si
  integer, private :: ih_totbgcc_si
  integer, private :: ih_totbgcc_old_si
  integer, private :: ih_biomass_stock_si
  integer, private :: ih_litter_stock_si
  integer, private :: ih_cwd_stock_si
  integer, private :: ih_cbal_err_fates_si
  integer, private :: ih_cbal_err_bgc_si
  integer, private :: ih_cbal_err_tot_si
  integer, private :: ih_npatches_si
  integer, private :: ih_ncohorts_si
  
  ! Indices to (site x scpf) variables
  integer, private :: ih_nplant_si_scpf
  integer, private :: ih_gpp_si_scpf
  integer, private :: ih_npp_totl_si_scpf
  integer, private :: ih_npp_leaf_si_scpf
  integer, private :: ih_npp_seed_si_scpf
  integer, private :: ih_npp_fnrt_si_scpf
  integer, private :: ih_npp_bgsw_si_scpf
  integer, private :: ih_npp_bgdw_si_scpf
  integer, private :: ih_npp_agsw_si_scpf
  integer, private :: ih_npp_agdw_si_scpf
  integer, private :: ih_npp_stor_si_scpf
  integer, private :: ih_litt_leaf_si_scpf
  integer, private :: ih_litt_fnrt_si_scpf
  integer, private :: ih_litt_sawd_si_scpf
  integer, private :: ih_litt_ddwd_si_scpf
  integer, private :: ih_r_leaf_si_scpf
  integer, private :: ih_r_stem_si_scpf
  integer, private :: ih_r_root_si_scpf
  integer, private :: ih_r_stor_si_scpf

  integer, private :: ih_ddbh_si_scpf
  integer, private :: ih_ba_si_scpf
  integer, private :: ih_m1_si_scpf
  integer, private :: ih_m2_si_scpf
  integer, private :: ih_m3_si_scpf
  integer, private :: ih_m4_si_scpf
  integer, private :: ih_m5_si_scpf


  ! The number of variable dim/kind types we have defined (static)
  integer, parameter                :: n_iovar_dk = 6


  ! This structure is not allocated by thread, but the upper and lower boundaries
  ! of the dimension for each thread is saved in the clump_ entry
  type iovar_dim_type
     character(fates_short_string_length) :: name  ! This should match the name of the dimension
     integer :: lb                       ! lower bound
     integer :: ub                       ! upper bound
     integer,allocatable :: clump_lb(:)  ! lower bound of thread's portion of HIO array
     integer,allocatable :: clump_ub(:)  ! upper bound of thread's portion of HIO array
  end type iovar_dim_type
  

  
  ! This structure is allocated by thread, and must be calculated after the FATES
  ! sites are allocated, and their mapping to the HLM is identified.  This structure
  ! is not combined with iovar_bounds, because that one is multi-instanced.  This
  ! structure is used more during the update phase, wherease _bounds is used
  ! more for things like flushing
  type iovar_map_type
     integer, allocatable :: site_index(:)   ! maps site indexes to the HIO site position
     integer, allocatable :: patch1_index(:) ! maps site index to the HIO patch 1st position
  end type iovar_map_type


   
  ! This structure is not multi-threaded
  type iovar_dimkind_type
     character(fates_short_string_length) :: name  ! String labelling this IO type
     integer              :: ndims       ! number of dimensions in this IO type
     integer, allocatable :: dimsize(:)  ! The size of each dimension
     logical              :: active
     type(iovar_dim_type), pointer :: dim1_ptr
     type(iovar_dim_type), pointer :: dim2_ptr
  end type iovar_dimkind_type


  
  ! This type is instanteated in the HLM-FATES interface (clmfates_interfaceMod.F90)
  type iovar_def_type
     character(len=fates_short_string_length) :: vname
     character(len=fates_short_string_length) :: units
     character(len=fates_long_string_length) :: long
     character(len=fates_short_string_length) :: use_default  ! States whether a variable should be turned
                                         ! on the output files by default (active/inactive)
                                         ! It is a good idea to set inactive for very large
                                         ! or infrequently used output datasets
     character(len=fates_short_string_length) :: vtype
     character(len=fates_avg_flag_length) :: avgflag
     integer              :: upfreq  ! Update frequency (this is for checks and flushing)
                                     ! 1 = dynamics "dyn" (daily)
                                     ! 2 = production "prod" (prob model tstep)
     real(r8)             :: flushval
     type(iovar_dimkind_type),pointer :: iovar_dk_ptr
     ! Pointers (only one of these is allocated per variable)
     real(r8), pointer     :: r81d(:)
     real(r8), pointer     :: r82d(:,:)
     real(r8), pointer     :: r83d(:,:,:)
     integer,  pointer     :: int1d(:)
     integer,  pointer     :: int2d(:,:)
     integer,  pointer     :: int3d(:,:,:)
  end type iovar_def_type


  type, public :: fates_hio_interface_type
     
     ! Instance of the list of history output varialbes
     type(iovar_def_type), pointer :: hvars(:)
     integer                       :: n_hvars
     
     ! Instanteat one registry of the different dimension/kinds (dk)
     ! All output variables will have a pointer to one of these dk's
     type(iovar_dimkind_type), pointer :: iovar_dk(:)
     
     ! This is a structure that explains where FATES patch boundaries
     ! on each thread point to in the host IO array, this structure
     ! is allocated by number of threads
     type(iovar_dim_type) :: iopa_dim
     
     ! This is a structure that explains where FATES patch boundaries
     ! on each thread point to in the host IO array, this structure
     ! is allocated by number of threads
     type(iovar_dim_type) :: iosi_dim
     
     ! This is a structure that contains the boundaries for the
     ! ground level (includes rock) dimension
     type(iovar_dim_type) :: iogrnd_dim

     ! This is a structure that contains the boundaries for the
     ! number of size-class x pft dimension
     type(iovar_dim_type) :: ioscpf_dim


     type(iovar_map_type), pointer :: iovar_map(:)
     
   contains
     
     procedure, public :: update_history_dyn
     procedure, public :: update_history_prod
     procedure, public :: update_history_cbal
     procedure, public :: define_history_vars
     procedure, public :: set_history_var
     procedure, public :: init_iovar_dk_maps
     procedure, public :: iotype_index
     procedure, public :: set_dim_ptrs
     procedure, public :: get_hvar_bounds
     procedure, public :: dim_init
     procedure, public :: set_dim_thread_bounds
     procedure, private :: flush_hvars

  end type fates_hio_interface_type
   


contains

   ! ===================================================================================
  
   subroutine update_history_cbal(this,nc,nsites,sites)

     use EDtypesMod          , only : ed_site_type
     
     ! Arguments
     class(fates_hio_interface_type)                 :: this
     integer                 , intent(in)            :: nc   ! clump index
     integer                 , intent(in)            :: nsites
     type(ed_site_type)      , intent(inout), target :: sites(nsites)

     ! Locals
     integer  :: s        ! The local site index
     integer  :: io_si     ! The site index of the IO array
     
     
     associate( hio_nep_si            => this%hvars(ih_nep_si)%r81d, &
                 hio_nbp_si            => this%hvars(ih_nbp_si)%r81d, &
                 hio_fire_c_to_atm_si  => this%hvars(ih_fire_c_to_atm_si)%r81d, &
                 hio_totecosysc_si     => this%hvars(ih_totecosysc_si)%r81d, &
                 hio_cbal_err_fates_si => this%hvars(ih_cbal_err_fates_si)%r81d, &
                 hio_cbal_err_bgc_si   => this%hvars(ih_cbal_err_bgc_si)%r81d, &
                 hio_cbal_err_tot_si   => this%hvars(ih_cbal_err_tot_si)%r81d, &
                 hio_biomass_stock_si  => this%hvars(ih_biomass_stock_si)%r81d, &
                 hio_litter_stock_si   => this%hvars(ih_litter_stock_si)%r81d, &
                 hio_cwd_stock_si      => this%hvars(ih_cwd_stock_si)%r81d )

        ! ---------------------------------------------------------------------------------
        ! Flush arrays to values defined by %flushval (see registry entry in
        ! subroutine define_history_vars()
        ! ---------------------------------------------------------------------------------
        call this%flush_hvars(nc,upfreq_in=3)
        
        
        do s = 1,nsites
         
           io_si  = this%iovar_map(nc)%site_index(s)

           hio_nep_si(io_si) = sites(s)%nep
           hio_nbp_si(io_si) = sites(s)%nbp
           hio_fire_c_to_atm_si(io_si) = sites(s)%fire_c_to_atm
           hio_totecosysc_si(io_si) = sites(s)%totecosysc
           hio_cbal_err_fates_si(io_si) = sites(s)%cbal_err_fates
           hio_cbal_err_bgc_si(io_si) = sites(s)%cbal_err_bgc
           hio_cbal_err_tot_si(io_si) = sites(s)%cbal_err_tot
           hio_biomass_stock_si(io_si) = sites(s)%biomass_stock
           hio_litter_stock_si(io_si) = sites(s)%ed_litter_stock
           hio_cwd_stock_si(io_si) = sites(s)%cwd_stock

        end do

      end associate

   end subroutine update_history_cbal
   

  ! ====================================================================================
  
  subroutine update_history_dyn(this,nc,nsites,sites)
    
    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after Ecosystem Dynamics have been processed.
    ! ---------------------------------------------------------------------------------
    
    use EDtypesMod          , only : ed_site_type,   &
                                     ed_cohort_type, &
                                     ed_patch_type,  &
                                     AREA,           &
                                     sclass_ed,      &
                                     nlevsclass_ed
    use EDParamsMod      , only : ED_val_ag_biomass

    ! Arguments
    class(fates_hio_interface_type)                 :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    
    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa      ! The local "I"ndex of "PA"tches 
    integer  :: io_pa    ! The patch index of the IO array
    integer  :: io_pa1   ! The first patch index in the IO array for each site
    integer  :: io_soipa 
    integer  :: lb1,ub1,lb2,ub2  ! IO array bounds for the calling thread
    integer  :: ivar             ! index of IO variable object vector
    integer  :: ft               ! functional type index
    integer  :: scpf             ! index of the size-class x pft bin
    integer  :: sc               ! index of the size-class bin

    real(r8) :: n_density   ! individual of cohort per m2.
    real(r8) :: n_perm2     ! individuals per m2 for the whole column
    real(r8) :: patch_scaling_scalar ! ratio of canopy to patch area for counteracting patch scaling
    real(r8) :: dbh         ! diameter ("at breast height")

    type(iovar_def_type),pointer :: hvar
    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort

    real(r8), parameter :: daysecs = 86400.0_r8 ! What modeler doesn't recognize 86400?
    real(r8), parameter :: yeardays = 365.0_r8  ! ALM/CLM do not use leap-years
    
    associate( hio_npatches_si         => this%hvars(ih_npatches_si)%r81d, &
               hio_ncohorts_si         => this%hvars(ih_ncohorts_si)%r81d, &
               hio_trimming_pa         => this%hvars(ih_trimming_pa)%r81d, &
               hio_area_plant_pa       => this%hvars(ih_area_plant_pa)%r81d, &
               hio_area_treespread_pa  => this%hvars(ih_area_treespread_pa)%r81d, & 
               hio_canopy_spread_pa    => this%hvars(ih_canopy_spread_pa)%r81d, &
               hio_biomass_pa_pft      => this%hvars(ih_biomass_pa_pft)%r82d, &
               hio_leafbiomass_pa_pft  => this%hvars(ih_leafbiomass_pa_pft)%r82d, &
               hio_storebiomass_pa_pft => this%hvars(ih_storebiomass_pa_pft)%r82d, &
               hio_nindivs_pa_pft      => this%hvars(ih_nindivs_pa_pft)%r82d, &
               hio_nesterov_fire_danger_pa => this%hvars(ih_nesterov_fire_danger_pa)%r81d, &
               hio_spitfire_ros_pa     => this%hvars(ih_spitfire_ROS_pa)%r81d, &
               hio_tfc_ros_pa          => this%hvars(ih_TFC_ROS_pa)%r81d, &
               hio_effect_wspeed_pa    => this%hvars(ih_effect_wspeed_pa)%r81d, &
               hio_fire_intensity_pa   => this%hvars(ih_fire_intensity_pa)%r81d, &
               hio_fire_area_pa        => this%hvars(ih_fire_area_pa)%r81d, &
               hio_scorch_height_pa    => this%hvars(ih_scorch_height_pa)%r81d, &
               hio_fire_fuel_bulkd_pa  => this%hvars(ih_fire_fuel_bulkd_pa)%r81d, &
               hio_fire_fuel_eff_moist_pa => this%hvars(ih_fire_fuel_eff_moist_pa)%r81d, &
               hio_fire_fuel_sav_pa    => this%hvars(ih_fire_fuel_sav_pa)%r81d, &
               hio_fire_fuel_mef_pa    => this%hvars(ih_fire_fuel_mef_pa)%r81d, &
               hio_sum_fuel_pa         => this%hvars(ih_sum_fuel_pa)%r81d,  &
               hio_litter_in_pa        => this%hvars(ih_litter_in_pa)%r81d, &
               hio_litter_out_pa       => this%hvars(ih_litter_out_pa)%r81d, &
               hio_seed_bank_si        => this%hvars(ih_seed_bank_si)%r81d, &
               hio_seeds_in_pa         => this%hvars(ih_seeds_in_pa)%r81d, &
               hio_seed_decay_pa       => this%hvars(ih_seed_decay_pa)%r81d, &
               hio_seed_germination_pa => this%hvars(ih_seed_germination_pa)%r81d, &
               hio_bstore_pa           => this%hvars(ih_bstore_pa)%r81d, &
               hio_bdead_pa            => this%hvars(ih_bdead_pa)%r81d, &
               hio_balive_pa           => this%hvars(ih_balive_pa)%r81d, &
               hio_bleaf_pa            => this%hvars(ih_bleaf_pa)%r81d, &
               hio_btotal_pa           => this%hvars(ih_btotal_pa)%r81d, &
               hio_gpp_si_scpf         => this%hvars(ih_gpp_si_scpf)%r82d, &
               hio_npp_totl_si_scpf    => this%hvars(ih_npp_totl_si_scpf)%r82d, &
               hio_npp_leaf_si_scpf    => this%hvars(ih_npp_leaf_si_scpf)%r82d, &
               hio_npp_seed_si_scpf    => this%hvars(ih_npp_seed_si_scpf)%r82d, &
               hio_npp_fnrt_si_scpf    => this%hvars(ih_npp_fnrt_si_scpf)%r82d, &
               hio_npp_bgsw_si_scpf    => this%hvars(ih_npp_bgsw_si_scpf)%r82d, &
               hio_npp_bgdw_si_scpf    => this%hvars(ih_npp_bgdw_si_scpf)%r82d, &
               hio_npp_agsw_si_scpf    => this%hvars(ih_npp_agsw_si_scpf)%r82d, &
               hio_npp_agdw_si_scpf    => this%hvars(ih_npp_agdw_si_scpf)%r82d, &
               hio_npp_stor_si_scpf    => this%hvars(ih_npp_stor_si_scpf)%r82d, &
               hio_ddbh_si_scpf        => this%hvars(ih_ddbh_si_scpf)%r82d, &
               hio_ba_si_scpf          => this%hvars(ih_ba_si_scpf)%r82d, &
               hio_nplant_si_scpf      => this%hvars(ih_nplant_si_scpf)%r82d, &
               hio_m1_si_scpf          => this%hvars(ih_m1_si_scpf)%r82d, &
               hio_m2_si_scpf          => this%hvars(ih_m2_si_scpf)%r82d, &
               hio_m3_si_scpf          => this%hvars(ih_m3_si_scpf)%r82d, &
               hio_m4_si_scpf          => this%hvars(ih_m4_si_scpf)%r82d, &
               hio_m5_si_scpf          => this%hvars(ih_m5_si_scpf)%r82d )
               
      ! ---------------------------------------------------------------------------------
      ! Flush arrays to values defined by %flushval (see registry entry in
      ! subroutine define_history_vars()
      ! ---------------------------------------------------------------------------------
      call this%flush_hvars(nc,upfreq_in=1)

      ! ---------------------------------------------------------------------------------
      ! Loop through the FATES scale hierarchy and fill the history IO arrays
      ! ---------------------------------------------------------------------------------
      
      do s = 1,nsites
         
         io_si  = this%iovar_map(nc)%site_index(s)
         io_pa1 = this%iovar_map(nc)%patch1_index(s)
         io_soipa = io_pa1-1
         
         ! Set trimming on the soil patch to 1.0
         hio_trimming_pa(io_soipa) = 1.0_r8

         ! The seed bank is a site level variable
         hio_seed_bank_si(io_si) = sum(sites(s)%seed_bank) * 1.e3_r8

         ipa = 0
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            
            io_pa = io_pa1 + ipa

            ! Increment the number of patches per site
            hio_npatches_si(io_si) = hio_npatches_si(io_si) + 1._r8
            
            ccohort => cpatch%shortest
            do while(associated(ccohort))
               
               ft = ccohort%pft
               
               ! Increment the number of cohorts per site
               hio_ncohorts_si(io_si) = hio_ncohorts_si(io_si) + 1._r8
               
               if ((cpatch%area .gt. 0._r8) .and. (cpatch%total_canopy_area .gt. 0._r8)) then
                  
                  ! for quantities that are at the CLM patch level, because of the way 
                  ! that CLM patches are weighted for radiative purposes this # density needs 
                  ! to be over either ED patch canopy area or ED patch total area, whichever is less
                  n_density = ccohort%n/min(cpatch%area,cpatch%total_canopy_area) 
                  
                  ! for quantities that are natively at column level, calculate plant 
                  ! density using whole area
                  n_perm2   = ccohort%n/AREA   
                  
               else
                  n_density = 0.0_r8
                  n_perm2   = 0.0_r8
               endif
               
               if(associated(cpatch%tallest))then
                  hio_trimming_pa(io_pa) = cpatch%tallest%canopy_trim
               else
                  hio_trimming_pa(io_pa) = 0.0_r8
               endif
               
               hio_area_plant_pa(io_pa) = 1._r8
               
               if (min(cpatch%total_canopy_area,cpatch%area)>0.0_r8) then
                  hio_area_treespread_pa(io_pa) = cpatch%total_tree_area  &
                       / min(cpatch%total_canopy_area,cpatch%area)
               else
                  hio_area_treespread_pa(io_pa) = 0.0_r8
               end if
               
               ! Update biomass components
               hio_bleaf_pa(io_pa)  = hio_bleaf_pa(io_pa)  + n_density * ccohort%bl       * 1.e3_r8
               hio_bstore_pa(io_pa) = hio_bstore_pa(io_pa) + n_density * ccohort%bstore   * 1.e3_r8
               hio_btotal_pa(io_pa) = hio_btotal_pa(io_pa) + n_density * ccohort%b        * 1.e3_r8
               hio_bdead_pa(io_pa)  = hio_bdead_pa(io_pa)  + n_density * ccohort%bdead    * 1.e3_r8
               hio_balive_pa(io_pa) = hio_balive_pa(io_pa) + n_density * ccohort%balive   * 1.e3_r8
               
               ! Update PFT partitioned biomass components
               hio_biomass_pa_pft(io_pa,ft) = hio_biomass_pa_pft(io_pa,ft) + &
                    n_density * ccohort%b * 1.e3_r8
               
               hio_leafbiomass_pa_pft(io_pa,ft) = hio_leafbiomass_pa_pft(io_pa,ft) + &
                    n_density * ccohort%bl       * 1.e3_r8
             
               hio_storebiomass_pa_pft(io_pa,ft) = hio_storebiomass_pa_pft(io_pa,ft) + &
                    n_density * ccohort%bstore   * 1.e3_r8
               
               hio_nindivs_pa_pft(io_pa,ft) = hio_nindivs_pa_pft(io_pa,ft) + &
                    ccohort%n

               ! Site by Size-Class x PFT (SCPF) 
               ! ------------------------------------------------------------------------

               dbh = ccohort%dbh !-0.5*(1./365.25)*ccohort%ddbhdt
               sc  = count(dbh-sclass_ed.ge.0.0)
               scpf = (ft-1)*nlevsclass_ed+sc

               ! Flux Variables (cohorts must had experienced a day before any of these values
               ! have any meaning, otherwise they are just inialization values
               if( .not.(ccohort%isnew) ) then

                  hio_gpp_si_scpf(io_si,scpf)      = hio_gpp_si_scpf(io_si,scpf)      + &
                                                     n_perm2*ccohort%gpp_acc_hold ! [kgC/m2/yr]
                  hio_npp_totl_si_scpf(io_si,scpf) = hio_npp_totl_si_scpf(io_si,scpf) + &
                                                     ccohort%npp_acc_hold*n_perm2
                  hio_npp_leaf_si_scpf(io_si,scpf) = hio_npp_leaf_si_scpf(io_si,scpf) + &
                                                     ccohort%npp_leaf*n_perm2
                  hio_npp_fnrt_si_scpf(io_si,scpf) = hio_npp_fnrt_si_scpf(io_si,scpf) + &
                                                     ccohort%npp_froot*n_perm2
                  hio_npp_bgsw_si_scpf(io_si,scpf) = hio_npp_bgsw_si_scpf(io_si,scpf) + &
                                                     ccohort%npp_bsw*(1._r8-ED_val_ag_biomass)*n_perm2
                  hio_npp_agsw_si_scpf(io_si,scpf) = hio_npp_agsw_si_scpf(io_si,scpf) + &
                                                     ccohort%npp_bsw*ED_val_ag_biomass*n_perm2
                  hio_npp_bgdw_si_scpf(io_si,scpf) = hio_npp_bgdw_si_scpf(io_si,scpf) + &
                                                     ccohort%npp_bdead*(1._r8-ED_val_ag_biomass)*n_perm2
                  hio_npp_agdw_si_scpf(io_si,scpf) = hio_npp_agdw_si_scpf(io_si,scpf) + &
                                                     ccohort%npp_bdead*ED_val_ag_biomass*n_perm2
                  hio_npp_seed_si_scpf(io_si,scpf) = hio_npp_seed_si_scpf(io_si,scpf) + &
                                                     ccohort%npp_bseed*n_perm2
                  hio_npp_stor_si_scpf(io_si,scpf) = hio_npp_stor_si_scpf(io_si,scpf) + &
                                                     ccohort%npp_store*n_perm2

                  if( abs(ccohort%npp_acc_hold-(ccohort%npp_leaf+ccohort%npp_froot+ &
                        ccohort%npp_bsw+ccohort%npp_bdead+ &
                        ccohort%npp_bseed+ccohort%npp_store))>1.e-9) then
                     write(fates_log(),*) 'NPP Partitions are not balancing'
                     write(fates_log(),*) 'Fractional Error: ', &
                          abs(ccohort%npp_acc_hold-(ccohort%npp_leaf+ccohort%npp_froot+ &
                           ccohort%npp_bsw+ccohort%npp_bdead+ &
                           ccohort%npp_bseed+ccohort%npp_store))/ccohort%npp_acc_hold
                     write(fates_log(),*) 'Terms: ',ccohort%npp_acc_hold,ccohort%npp_leaf,ccohort%npp_froot, &
                           ccohort%npp_bsw,ccohort%npp_bdead, &
                           ccohort%npp_bseed,ccohort%npp_store
                     write(fates_log(),*) ' NPP components during FATES-HLM linking does not balance '
                     stop ! we need termination control for FATES!!!
                     ! call endrun(msg=errMsg(__FILE__, __LINE__))
                  end if
                  
                  ! Woody State Variables (basal area and number density and mortality)
                  if (pftcon%woody(ft) == 1) then

                     hio_m1_si_scpf(io_si,scpf) = hio_m1_si_scpf(io_si,scpf) + ccohort%bmort*n_perm2*AREA
                     hio_m2_si_scpf(io_si,scpf) = hio_m2_si_scpf(io_si,scpf) + ccohort%hmort*n_perm2*AREA
                     hio_m3_si_scpf(io_si,scpf) = hio_m3_si_scpf(io_si,scpf) + ccohort%cmort*n_perm2*AREA
                     hio_m4_si_scpf(io_si,scpf) = hio_m4_si_scpf(io_si,scpf) + ccohort%imort*n_perm2*AREA
                     hio_m5_si_scpf(io_si,scpf) = hio_m5_si_scpf(io_si,scpf) + ccohort%fmort*n_perm2*AREA

                     ! basal area  [m2/ha]
                     hio_ba_si_scpf(io_si,scpf) = hio_ba_si_scpf(io_si,scpf) + &
                           0.25_r8*3.14159_r8*((dbh/100.0_r8)**2.0_r8)*n_perm2*AREA
                     
                     ! number density [/ha]
                     hio_nplant_si_scpf(io_si,scpf) = hio_nplant_si_scpf(io_si,scpf) + AREA*n_perm2
                     
                     ! Growth Incrments must have NaN check and woody check
                     if(ccohort%ddbhdt == ccohort%ddbhdt) then
                        hio_ddbh_si_scpf(io_si,scpf) = hio_ddbh_si_scpf(io_si,scpf) + &
                              ccohort%ddbhdt*n_perm2*AREA
                     else
                        hio_ddbh_si_scpf(io_si,scpf) = -999.9
                     end if
                  end if

               end if
               
               ccohort => ccohort%taller
            enddo ! cohort loop
            
            ! Patch specific variables that are already calculated
            ! These things are all duplicated. Should they all be converted to LL or array structures RF? 
            ! define scalar to counteract the patch albedo scaling logic for conserved quantities
            
            if (cpatch%area .gt. 0._r8 .and. cpatch%total_canopy_area .gt.0 ) then
               patch_scaling_scalar  = min(1._r8, cpatch%area / cpatch%total_canopy_area)
            else
               patch_scaling_scalar = 0._r8
            endif
            
            ! Update Fire Variables
            hio_nesterov_fire_danger_pa(io_pa) = sites(s)%acc_NI
            hio_spitfire_ros_pa(io_pa)         = cpatch%ROS_front 
            hio_effect_wspeed_pa(io_pa)        = cpatch%effect_wspeed
            hio_tfc_ros_pa(io_pa)              = cpatch%TFC_ROS
            hio_fire_intensity_pa(io_pa)       = cpatch%FI
            hio_fire_area_pa(io_pa)            = cpatch%frac_burnt
            hio_scorch_height_pa(io_pa)        = cpatch%SH
            hio_fire_fuel_bulkd_pa(io_pa)      = cpatch%fuel_bulkd
            hio_fire_fuel_eff_moist_pa(io_pa)  = cpatch%fuel_eff_moist
            hio_fire_fuel_sav_pa(io_pa)        = cpatch%fuel_sav
            hio_fire_fuel_mef_pa(io_pa)        = cpatch%fuel_mef
            hio_sum_fuel_pa(io_pa)             = cpatch%sum_fuel * 1.e3_r8 * patch_scaling_scalar
            
            ! Update Litter Flux Variables
            hio_litter_in_pa(io_pa)            = (sum(cpatch%CWD_AG_in) +sum(cpatch%leaf_litter_in)) &
                 * 1.e3_r8 * 365.0_r8 * daysecs * patch_scaling_scalar
            hio_litter_out_pa(io_pa)           = (sum(cpatch%CWD_AG_out)+sum(cpatch%leaf_litter_out)) &
                 * 1.e3_r8 * 365.0_r8 * daysecs * patch_scaling_scalar
            
            hio_seeds_in_pa(io_pa)             = sum(cpatch%seeds_in) * &
                 1.e3_r8 * 365.0_r8 * daysecs * patch_scaling_scalar
            hio_seed_decay_pa(io_pa)           = sum(cpatch%seed_decay) &
                 * 1.e3_r8 * 365.0_r8 * daysecs * patch_scaling_scalar
            hio_seed_germination_pa(io_pa)     = sum(cpatch%seed_germination) &
                 * 1.e3_r8 * 365.0_r8 * daysecs * patch_scaling_scalar

            
            hio_canopy_spread_pa(io_pa)        = cpatch%spread(1) 
            
            
            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop
       
      enddo ! site loop
      
    end associate

    return
  end subroutine update_history_dyn
 
 ! ======================================================================================

 subroutine update_history_prod(this,nc,nsites,sites,dt_tstep)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after rapid timescale productivity calculations (gpp and respiration).
    ! ---------------------------------------------------------------------------------
    
    use EDtypesMod          , only : ed_site_type,   &
                                     ed_cohort_type, &
                                     ed_patch_type,  &
                                     AREA
    ! Arguments
    class(fates_hio_interface_type)                 :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    real(r8)                , intent(in)            :: dt_tstep
    
    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa      ! The local "I"ndex of "PA"tches 
    integer  :: io_pa    ! The patch index of the IO array
    integer  :: io_pa1   ! The first patch index in the IO array for each site
    integer  :: io_soipa 
    integer  :: lb1,ub1,lb2,ub2  ! IO array bounds for the calling thread
    integer  :: ivar             ! index of IO variable object vector

    real(r8) :: n_density   ! individual of cohort per m2.
    real(r8) :: n_perm2     ! individuals per m2 for the whole column

    type(iovar_def_type),pointer :: hvar
    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort

    real(r8), parameter :: daysecs = 86400.0_r8 ! What modeler doesn't recognize 86400?
    real(r8), parameter :: yeardays = 365.0_r8  ! Should this be 365.25?
    
    associate( hio_gpp_pa         => this%hvars(ih_gpp_pa)%r81d, &
               hio_npp_pa         => this%hvars(ih_npp_pa)%r81d, &
               hio_aresp_pa       => this%hvars(ih_aresp_pa)%r81d, &
               hio_maint_resp_pa  => this%hvars(ih_maint_resp_pa)%r81d, &
               hio_growth_resp_pa => this%hvars(ih_growth_resp_pa)%r81d, &
               hio_npp_si         => this%hvars(ih_npp_si)%r81d )

      ! Flush the relevant history variables 
      call this%flush_hvars(nc,upfreq_in=2)

      do s = 1,nsites
         
         io_si  = this%iovar_map(nc)%site_index(s)
         io_pa1 = this%iovar_map(nc)%patch1_index(s)
         io_soipa = io_pa1-1
         
         ipa = 0
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            
            io_pa = io_pa1 + ipa

            ccohort => cpatch%shortest
            do while(associated(ccohort))
               
               ! TODO: we need a standardized logical function on this (used lots, RGK)
               if ((cpatch%area .gt. 0._r8) .and. (cpatch%total_canopy_area .gt. 0._r8)) then
                  n_density = ccohort%n/min(cpatch%area,cpatch%total_canopy_area) 
                  n_perm2   = ccohort%n/AREA   
               else
                  n_density = 0.0_r8
                  n_perm2   = 0.0_r8
               endif
               
               if ( .not. ccohort%isnew ) then
                  
                  ! scale up cohort fluxes to their patches
                  hio_npp_pa(io_pa) = hio_npp_pa(io_pa) + &
                        ccohort%npp_tstep * 1.e3_r8 * n_density / dt_tstep
                  hio_gpp_pa(io_pa) = hio_gpp_pa(io_pa) + &
                        ccohort%gpp_tstep * 1.e3_r8 * n_density / dt_tstep
                  hio_aresp_pa(io_pa) = hio_aresp_pa(io_pa) + &
                        ccohort%resp_tstep * 1.e3_r8 * n_density / dt_tstep
                  hio_growth_resp_pa(io_pa) = hio_growth_resp_pa(io_pa) + &
                        ccohort%resp_g * 1.e3_r8 * n_density / dt_tstep
                  hio_maint_resp_pa(io_pa) = hio_maint_resp_pa(io_pa) + &
                        ccohort%resp_m * 1.e3_r8 * n_density / dt_tstep
                  
                  ! map ed cohort-level npp fluxes to clm column fluxes
                  hio_npp_si(io_si) = hio_npp_si(io_si) + ccohort%npp_tstep * n_perm2 * 1.e3_r8 /dt_tstep

               endif

               ccohort => ccohort%taller
            enddo ! cohort loop
            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop
         
      enddo ! site loop

    end associate
 
  end subroutine update_history_prod

 ! ======================================================================================

 subroutine flush_hvars(this,nc,upfreq_in)
 
   class(fates_hio_interface_type)        :: this
   integer,intent(in)                     :: nc
   integer,intent(in)                     :: upfreq_in

   integer                      :: ivar
   type(iovar_def_type),pointer :: hvar
   integer                      :: lb1,ub1,lb2,ub2


   do ivar=1,ubound(this%hvars,1)
      hvar => this%hvars(ivar)
      if (hvar%upfreq==upfreq_in) then ! Only flush variables with update on dynamics step
         call this%get_hvar_bounds(hvar,nc,lb1,ub1,lb2,ub2)
         select case(trim(hvar%iovar_dk_ptr%name))
         case('PA_R8') 
            hvar%r81d(lb1:ub1) = hvar%flushval
         case('SI_R8') 
            hvar%r81d(lb1:ub1) = hvar%flushval
         case('PA_GRND_R8') 
            hvar%r82d(lb1:ub1,lb2:ub2) = hvar%flushval
         case('PA_SCPF_R8') 
            hvar%r82d(lb1:ub1,lb2:ub2) = hvar%flushval
         case('SI_GRND_R8') 
            hvar%r82d(lb1:ub1,lb2:ub2) = hvar%flushval
         case('SI_SCPF_R8') 
            hvar%r82d(lb1:ub1,lb2:ub2) = hvar%flushval
         case('PA_INT')
            hvar%int1d(lb1:ub1) = nint(hvar%flushval)
         case default
            write(fates_log(),*) 'iotyp undefined while flushing history variables'
            stop
            !end_run
         end select
      end if
   end do
   
 end subroutine flush_hvars

 ! ====================================================================================
  
  subroutine define_history_vars(this,callstep,nvar)
    
    ! ---------------------------------------------------------------------------------
    ! 
    !                    REGISTRY OF HISTORY OUTPUT VARIABLES
    !
    ! This subroutine is called in two contexts, either in count mode or inialize mode
    ! In count mode, we just walk through the list of registerred variables, compare
    ! if the variable of interest list the current host model and add it to the count
    ! if true.  This count is used just to allocate the variable space.  After this
    ! has been done, we go through the list a second time populating a memory structure.
    ! This phase is the "initialize" phase.  These two phases are differntiated by the
    ! string "callstep", which should be either "count" or "initialize".
    !
    ! Note 1 there are different ways you can flush or initialize the output fields.
    ! If you flush to a native type, (such as zero), the entire slab which covers
    ! indices which may not be relevant to FATES, are flushed to this value.  So
    ! in that case, lakes and crops that are not controlled by FATES will zero'd
    ! and when values are scaled up to the land-grid, the zero's for non FATES will
    ! be included.  This is good and correct if nothing is there.  
    !
    ! But, what if crops exist in the host model and occupy a fraction of the land-surface
    ! shared with natural vegetation? In that case, you want to flush your arrays
    ! with a value that the HLM treats as "do not average"
    ! 
    ! If your HLM makes use of, and you want, INTEGER OUTPUT, pass the flushval as
    ! a real.  The applied flush value will use the NINT() intrinsic function
    ! ---------------------------------------------------------------------------------
    
    class(fates_hio_interface_type)        :: this
    character(len=*),intent(in)            :: callstep  ! are we 'count'ing or 'initializ'ing?
    integer,optional,intent(out)           :: nvar
    integer                                :: ivar
    
    if(.not. (trim(callstep).eq.'count' .or. trim(callstep).eq.'initialize') ) then
       write(fates_log(),*) 'defining history variables in FATES requires callstep count or initialize'
       ! end_run('MESSAGE')
    end if
    
    ivar=0
    
    ! Site level counting variables
    call this%set_history_var(vname='ED_NPATCHES',units='none',                &
         long='Total number of ED patches per site', use_default='active',      &
         avgflag='A',vtype='SI_R8',hlms='CLM:ALM',flushval=1.0_r8, upfreq=1,    &
         ivar=ivar, callstep=callstep,index = ih_npatches_si)

    call this%set_history_var(vname='ED_NCOHORTS',units='none',                &
         long='Total number of ED cohorts per site', use_default='active',      &
         avgflag='A',vtype='SI_R8',hlms='CLM:ALM',flushval=1.0_r8, upfreq=1,    &
         ivar=ivar, callstep=callstep,index = ih_ncohorts_si)
    
    ! Patch variables
    call this%set_history_var(vname='TRIMMING',units='none',                   &
         long='Degree to which canopy expansion is limited by leaf economics',  & 
         use_default='active', &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',flushval=1.0_r8, upfreq=1,    &
         ivar=ivar, callstep=callstep,index = ih_trimming_pa)
    
    call this%set_history_var(vname='AREA_PLANT',units='m2',                   &
         long='area occupied by all plants', use_default='active',              &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, callstep=callstep,index = ih_area_plant_pa)
    
    call this%set_history_var(vname='AREA_TREES',units='m2',                   &
         long='area occupied by woody plants', use_default='active',            &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, callstep=callstep,index = ih_area_treespread_pa)

    call this%set_history_var(vname='CANOPY_SPREAD',units='0-1',               &
         long='Scaling factor between tree basal area and canopy area',         &
         use_default='active',                                                  &
         avgflag='A',vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, callstep=callstep,index = ih_canopy_spread_pa)

    call this%set_history_var(vname='PFTbiomass',units='gC/m2',                   &
         long='total PFT level biomass', use_default='active',                     &
         avgflag='A', vtype='PA_GRND_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1, &
         ivar=ivar, callstep=callstep, index = ih_biomass_pa_pft )

    call this%set_history_var(vname='PFTleafbiomass', units='gC/m2',              &
         long='total PFT level leaf biomass', use_default='active',                &
         avgflag='A', vtype='PA_GRND_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1, &
         ivar=ivar, callstep=callstep, index = ih_leafbiomass_pa_pft )

    call this%set_history_var(vname='PFTstorebiomass',  units='gC/m2',            &
         long='total PFT level stored biomass', use_default='active',              &
         avgflag='A', vtype='PA_GRND_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1, &
         ivar=ivar, callstep=callstep, index = ih_storebiomass_pa_pft )
    
    call this%set_history_var(vname='PFTnindivs',  units='indiv / m2',            &
         long='total PFT level number of individuals', use_default='active',       &
         avgflag='A', vtype='PA_GRND_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1, &
         ivar=ivar, callstep=callstep, index = ih_nindivs_pa_pft )

    ! Fire Variables

    call this%set_history_var(vname='FIRE_NESTEROV_INDEX', units='none',       &
         long='nesterov_fire_danger index', use_default='active',               &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_nesterov_fire_danger_pa)

    call this%set_history_var(vname='FIRE_ROS', units='m/min',                 &
         long='fire rate of spread m/min', use_default='active',                &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_spitfire_ROS_pa)

    call this%set_history_var(vname='EFFECT_WSPEED', units='none',             &
         long ='effective windspeed for fire spread', use_default='active',     &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_effect_wspeed_pa )

    call this%set_history_var(vname='FIRE_TFC_ROS', units='none',              &
         long ='total fuel consumed', use_default='active',                     &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_TFC_ROS_pa )

    call this%set_history_var(vname='FIRE_INTENSITY', units='kJ/m/s',          &
         long='spitfire fire intensity: kJ/m/s', use_default='active',          &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_fire_intensity_pa )

    call this%set_history_var(vname='FIRE_AREA', units='fraction',             &
         long='spitfire fire area:m2', use_default='active',                    &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_fire_area_pa )

    call this%set_history_var(vname='SCORCH_HEIGHT', units='m',                &
         long='spitfire fire area:m2', use_default='active',                    &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_scorch_height_pa )

    call this%set_history_var(vname='fire_fuel_mef', units='m',                &
         long='spitfire fuel moisture',  use_default='active',                  &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_fire_fuel_mef_pa )

    call this%set_history_var(vname='fire_fuel_bulkd', units='m',              &
         long='spitfire fuel bulk density',  use_default='active',              &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_fire_fuel_bulkd_pa )

    call this%set_history_var(vname='FIRE_FUEL_EFF_MOIST', units='m',          &
         long='spitfire fuel moisture', use_default='active',                   &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_fire_fuel_eff_moist_pa )

    call this%set_history_var(vname='fire_fuel_sav', units='m',                &
         long='spitfire fuel surface/volume ',  use_default='active',           &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_fire_fuel_sav_pa )

    call this%set_history_var(vname='SUM_FUEL', units='gC m-2',                &
         long='total ground fuel related to ros (omits 1000hr fuels)',          & 
         use_default='active',                                                  & 
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_sum_fuel_pa )

    ! Litter Variables

    call this%set_history_var(vname='LITTER_IN', units='gC m-2 s-1',           &
         long='Litter flux in leaves',  use_default='active',                   &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_litter_in_pa )

    call this%set_history_var(vname='LITTER_OUT', units='gC m-2 s-1',          &
         long='Litter flux out leaves',  use_default='active',                  & 
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_litter_out_pa )

    call this%set_history_var(vname='SEED_BANK', units='gC m-2',               &
         long='Total Seed Mass of all PFTs',  use_default='active',             &
         avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_seed_bank_si )

    call this%set_history_var(vname='SEEDS_IN', units='gC m-2 s-1',            &
         long='Seed Production Rate',  use_default='active',                    &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_seeds_in_pa )

    call this%set_history_var(vname='SEED_GERMINATION', units='gC m-2 s-1',    &
         long='Seed mass converted into new cohorts',   use_default='active',   &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_seed_germination_pa )

    call this%set_history_var(vname='SEED_DECAY', units='gC m-2 s-1',          &
         long='Seed mass decay', use_default='active',                          &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_seed_decay_pa )
    
    call this%set_history_var(vname='ED_bstore', units='gC m-2',                  &
         long='Storage biomass', use_default='active',                          &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_bstore_pa )

    call this%set_history_var(vname='ED_bdead', units='gC m-2',                   &
         long='Dead (structural) biomass (live trees, not CWD)',                &
         use_default='active',                                                  &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_bdead_pa )

    call this%set_history_var(vname='ED_balive', units='gC m-2',                  &
         long='Live biomass', use_default='active',                             &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_balive_pa )

    call this%set_history_var(vname='ED_bleaf', units='gC m-2',                   &
         long='Leaf biomass',  use_default='active',                            &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_bleaf_pa )

    call this%set_history_var(vname='ED_biomass', units='gC m-2',                  &
         long='Total biomass',  use_default='active',                           &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=1,   &
         ivar=ivar,callstep=callstep, index = ih_btotal_pa )

    
    ! Ecosystem Carbon Fluxes (updated rapidly, upfreq=2)

    call this%set_history_var(vname='NPP_column', units='gC/m^2/s',                &
         long='net primary production on the site',  use_default='active',      &
         avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=2,   &
         ivar=ivar,callstep=callstep, index = ih_npp_si )

    call this%set_history_var(vname='GPP', units='gC/m^2/s',                   &
         long='gross primary production',  use_default='active',                &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=2,   &
         ivar=ivar,callstep=callstep, index = ih_gpp_pa )

    call this%set_history_var(vname='NPP', units='gC/m^2/s',                   &
         long='net primary production', use_default='active',                   &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=2,   &
         ivar=ivar,callstep=callstep, index = ih_npp_pa )

    call this%set_history_var(vname='AR', units='gC/m^2/s',                 &
         long='autotrophic respiration', use_default='active',                  &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=2,   &
         ivar=ivar,callstep=callstep, index = ih_aresp_pa )

    call this%set_history_var(vname='GROWTH_RESP', units='gC/m^2/s',           &
         long='growth respiration', use_default='active',                       &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=2,   &
         ivar=ivar,callstep=callstep, index = ih_growth_resp_pa )

    call this%set_history_var(vname='MAINT_RESP', units='gC/m^2/s',            &
         long='maintenance respiration', use_default='active',                  &
         avgflag='A', vtype='PA_R8',hlms='CLM:ALM',flushval=0.0_r8, upfreq=2,   &
         ivar=ivar,callstep=callstep, index = ih_maint_resp_pa )


    ! Carbon Flux (grid dimension x scpf) (THESE ARE DEFAULT INACTIVE!!!
    !                                     (BECAUSE THEY TAKE UP SPACE!!!
    ! ===================================================================================

    call this%set_history_var(vname='GPP_SCPF',units='kgC/m2/yr',            &
          long='gross primary production', use_default='inactive',           &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_gpp_si_scpf )

    call this%set_history_var(vname='NPP_SCPF',units='kgC/m2/yr',            &
          long='total net primary production', use_default='inactive',       &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_npp_totl_si_scpf )


    call this%set_history_var(vname='NPP_LEAF_SCPF',units='kgC/m2/yr',       &
          long='NPP flux into leaves', use_default='inactive',               &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_npp_leaf_si_scpf )

    call this%set_history_var(vname='NPP_SEED_SCPF',units='kgC/m2/yr',       &
          long='NPP flux into seeds', use_default='inactive',                &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_npp_seed_si_scpf )

    call this%set_history_var(vname='NPP_FNRT_SCPF',units='kgC/m2/yr',       &
          long='NPP flux into fine roots', use_default='inactive',           &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_npp_fnrt_si_scpf )

    call this%set_history_var(vname='NPP_BGSW_SCPF',units='kgC/m2/yr',       &
          long='NPP flux into below-ground sapwood', use_default='inactive', &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_npp_bgsw_si_scpf )

    call this%set_history_var(vname='NPP_BGDW_SCPF',units='kgC/m2/yr',       &
          long='NPP flux into below-ground deadwood', use_default='inactive',&
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_npp_bgdw_si_scpf )

    call this%set_history_var(vname='NPP_AGSW_SCPF',units='kgC/m2/yr',       &
          long='NPP flux into above-ground sapwood', use_default='inactive', &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_npp_agsw_si_scpf )

    call this%set_history_var(vname = 'NPP_AGDW_SCPF', units='kgC/m2/yr',    &
          long='NPP flux into above-ground deadwood', use_default='inactive',&
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_npp_agdw_si_scpf )

    call this%set_history_var(vname = 'NPP_STOR_SCPF', units='kgC/m2/yr',    &
          long='NPP flux into storage', use_default='inactive',              &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_npp_stor_si_scpf )

    call this%set_history_var(vname='DDBH_SCPF', units = 'cm/yr/ha',         &
          long='diameter growth increment and pft/size',use_default='inactive',&
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_ddbh_si_scpf )

    call this%set_history_var(vname='BA_SCPF',units = 'm2/ha',               &
          long='basal area by patch and pft/size', use_default='inactive',   &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_ba_si_scpf )

    call this%set_history_var(vname='NPLANT_SCPF',units = 'N/ha',         &
          long='stem number density by patch and pft/size', use_default='inactive', &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_nplant_si_scpf )

    call this%set_history_var(vname='M1_SCPF',units = 'N/ha/yr',          &
          long='background mortality count by patch and pft/size', use_default='inactive',&
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_m1_si_scpf )
    
    call this%set_history_var(vname='M2_SCPF',units = 'N/ha/yr',          &
          long='hydraulic mortality count by patch and pft/size',use_default='inactive',&
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_m2_si_scpf )

    call this%set_history_var(vname='M3_SCPF',units = 'N/ha/yr',          &
          long='carbon starvation mortality count by patch and pft/size',use_default='inactive', &
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_m3_si_scpf )

    call this%set_history_var(vname='M4_SCPF',units = 'N/ha/yr',          &
          long='impact mortality count by patch and pft/size',use_default='inactive',&
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_m4_si_scpf )

    call this%set_history_var(vname='M5_SCPF',units = 'N/ha/yr',          &
          long='fire mortality count by patch and pft/size',use_default='inactive',&
          avgflag='A', vtype='SI_SCPF_R8',hlms='CLM:ALM',flushval=0.0_r8,    &
          upfreq=1, ivar=ivar,callstep=callstep, index = ih_m5_si_scpf )


    ! CARBON BALANCE VARIABLES THAT DEPEND ON HLM BGC INPUTS

    call this%set_history_var(vname='NEP', units='gC/m^2/s', &
          long='net ecosystem production', use_default='active', &
          avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=cp_hio_ignore_val,    &
          upfreq=3, ivar=ivar,callstep=callstep, index = ih_nep_si )

    call this%set_history_var(vname='Fire_Closs', units='gC/m^2/s', &
          long='ED/SPitfire Carbon loss to atmosphere', use_default='active', &
          avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=cp_hio_ignore_val,    &
          upfreq=3, ivar=ivar,callstep=callstep, index = ih_fire_c_to_atm_si )
   
    call this%set_history_var(vname='NBP', units='gC/m^2/s', &
          long='net biosphere production', use_default='active', &
          avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=cp_hio_ignore_val,    &
          upfreq=3, ivar=ivar,callstep=callstep, index = ih_nbp_si )
   
    call this%set_history_var(vname='TOTECOSYSC', units='gC/m^2',  &
         long='total ecosystem carbon', use_default='active', &
         avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=cp_hio_ignore_val,    &
         upfreq=3, ivar=ivar,callstep=callstep, index = ih_totecosysc_si )
    
    call this%set_history_var(vname='CBALANCE_ERROR_ED', units='gC/m^2/s',  &
         long='total carbon balance error on ED side', use_default='active', &
         avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=cp_hio_ignore_val,    &
         upfreq=3, ivar=ivar,callstep=callstep, index = ih_cbal_err_fates_si )

    call this%set_history_var(vname='CBALANCE_ERROR_BGC', units='gC/m^2/s',  &
         long='total carbon balance error on HLMs BGC side', use_default='active', &
         avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=cp_hio_ignore_val,    &
         upfreq=3, ivar=ivar,callstep=callstep, index = ih_cbal_err_bgc_si )
    
    call this%set_history_var(vname='CBALANCE_ERROR_TOTAL', units='gC/m^2/s', &
          long='total carbon balance error total', use_default='active', &
          avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=cp_hio_ignore_val,    &
          upfreq=3, ivar=ivar,callstep=callstep, index = ih_cbal_err_tot_si )
    
    call this%set_history_var(vname='BIOMASS_STOCK_COL', units='gC/m^2',  &
          long='total ED biomass carbon at the column level', use_default='active', &
          avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=cp_hio_ignore_val,    &
          upfreq=3, ivar=ivar,callstep=callstep, index = ih_biomass_stock_si )
    
    call this%set_history_var(vname='ED_LITTER_STOCK_COL', units='gC/m^2', &
          long='total ED litter carbon at the column level', use_default='active', &
          avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=cp_hio_ignore_val,    &
          upfreq=3, ivar=ivar,callstep=callstep, index = ih_litter_stock_si )
    
    call this%set_history_var(vname='CWD_STOCK_COL', units='gC/m^2', &
          long='total CWD carbon at the column level', use_default='active', &
          avgflag='A', vtype='SI_R8',hlms='CLM:ALM',flushval=cp_hio_ignore_val,    &
          upfreq=3, ivar=ivar,callstep=callstep, index = ih_cwd_stock_si )
   

    ! Must be last thing before return
    if(present(nvar)) nvar = ivar
    
    return
    
  end subroutine define_history_vars
  
  ! =====================================================================================
   
  subroutine set_history_var(this,vname,units,long,use_default,avgflag,vtype,hlms, &
                             flushval,upfreq,ivar,callstep,index)


    use FatesUtilsMod, only : check_hlm_list
    use EDTypesMod, only    : cp_hlm_name

    ! arguments
    class(fates_hio_interface_type) :: this
    character(len=*),intent(in)  :: vname
    character(len=*),intent(in)  :: units
    character(len=*),intent(in)  :: long
    character(len=*),intent(in)  :: use_default
    character(len=*),intent(in)  :: avgflag
    character(len=*),intent(in)  :: vtype
    character(len=*),intent(in)  :: hlms
    real(r8),intent(in)          :: flushval ! IF THE TYPE IS AN INT WE WILL round with NINT
    integer,intent(in)           :: upfreq
    character(len=*),intent(in)  :: callstep
    integer, intent(inout)       :: ivar
    integer, intent(inout)       :: index  ! This is the index for the variable of
                                           ! interest that is associated with an
                                           ! explict name (for fast reference during update)
                                           ! A zero is passed back when the variable is
                                           ! not used

    ! locals
    type(iovar_def_type),pointer :: hvar
    integer :: ub1,lb1,ub2,lb2    ! Bounds for allocating the var
    integer :: ityp
    
    if( check_hlm_list(trim(hlms),trim(cp_hlm_name)) ) then
       
       ivar  = ivar+1
       index = ivar    
       
       if(trim(callstep).eq.'initialize')then
          
          hvar => this%hvars(ivar)
          hvar%vname = vname
          hvar%units = units
          hvar%long  = long
          hvar%use_default = use_default
          hvar%vtype = vtype
          hvar%avgflag = avgflag
          hvar%flushval = flushval
          hvar%upfreq = upfreq
          ityp=this%iotype_index(trim(vtype))
          hvar%iovar_dk_ptr => this%iovar_dk(ityp)
          this%iovar_dk(ityp)%active = .true.
          
          nullify(hvar%r81d)
          nullify(hvar%r82d)
          nullify(hvar%r83d)
          nullify(hvar%int1d)
          nullify(hvar%int2d)
          nullify(hvar%int3d)
          
          call this%get_hvar_bounds(hvar,0,lb1,ub1,lb2,ub2)
          
          ! currently, all array spaces are flushed each time
          ! the update is called. The flush here on the initialization
          ! may be redundant, but will prevent issues in the future
          ! if we have host models where not all threads are updating
          ! the HIO array spaces. (RGK:09-2016)

          select case(trim(vtype))
          case('PA_R8')
             allocate(hvar%r81d(lb1:ub1));hvar%r81d(:)=flushval
          case('SI_R8')
             allocate(hvar%r81d(lb1:ub1));hvar%r81d(:)=flushval
          case('PA_GRND_R8')
             allocate(hvar%r82d(lb1:ub1,lb2:ub2));hvar%r82d(:,:)=flushval
          case('PA_SCPF_R8')
             allocate(hvar%r82d(lb1:ub1,lb2:ub2));hvar%r82d(:,:)=flushval
          case('SI_GRND_R8')
             allocate(hvar%r82d(lb1:ub1,lb2:ub2));hvar%r82d(:,:)=flushval
          case('SI_SCPF_R8')
             allocate(hvar%r82d(lb1:ub1,lb2:ub2));hvar%r82d(:,:)=flushval
          case default
             write(fates_log(),*) 'Incompatible vtype passed to set_history_var'
             write(fates_log(),*) 'vtype = ',trim(vtype),' ?'
             stop
             ! end_run
          end select
          
       end if
    else
       
       index = 0
    end if
    
    return
  end subroutine set_history_var
  
  ! =====================================================================================

  subroutine get_hvar_bounds(this,hvar,thread,lb1,ub1,lb2,ub2)

     class(fates_hio_interface_type) :: this
     type(iovar_def_type),target,intent(in) :: hvar
     integer,intent(in)              :: thread
     integer,intent(out)             :: lb1
     integer,intent(out)             :: ub1
     integer,intent(out)             :: lb2
     integer,intent(out)             :: ub2

     ! local
     integer :: ndims

     lb1 = 0
     ub1 = 0
     lb2 = 0
     ub2 = 0

     ndims = hvar%iovar_dk_ptr%ndims

     ! The thread = 0 case is the boundaries for the whole proc/node
     if (thread==0) then
        lb1 = hvar%iovar_dk_ptr%dim1_ptr%lb
        ub1 = hvar%iovar_dk_ptr%dim1_ptr%ub
        if(ndims>1)then
           lb2 = hvar%iovar_dk_ptr%dim2_ptr%lb
           ub2 = hvar%iovar_dk_ptr%dim2_ptr%ub
        end if
     else
        lb1 = hvar%iovar_dk_ptr%dim1_ptr%clump_lb(thread)
        ub1 = hvar%iovar_dk_ptr%dim1_ptr%clump_ub(thread)
        if(ndims>1)then
           lb2 = hvar%iovar_dk_ptr%dim2_ptr%clump_lb(thread)
           ub2 = hvar%iovar_dk_ptr%dim2_ptr%clump_ub(thread)
        end if
     end if
     
     return
  end subroutine get_hvar_bounds


  ! ====================================================================================
  
  subroutine init_iovar_dk_maps(this)
    
    ! ----------------------------------------------------------------------------------
    ! This subroutine simply initializes the structures that define the different
    ! array and type formats for different IO variables
    !
    ! PA_R8   : 1D patch scale 8-byte reals
    ! SI_R8   : 1D site scale 8-byte reals
    !
    ! The allocation on the structures is not dynamic and should only add up to the
    ! number of entries listed here.
    !
    ! note (RGK) %active is not used yet. Was intended as a check on the HLM->FATES
    ! control parameter passing to ensure all active dimension types received all
    ! dimensioning specifications from the host, but we currently arent using those
    ! passing functions..
    ! ----------------------------------------------------------------------------------
    
    ! Arguments
    class(fates_hio_interface_type) :: this
       
    ! Locals
    integer            :: ityp
    integer, parameter :: unset_int = -999
    
    allocate(this%iovar_dk(n_iovar_dk))

    ! 1d Patch
    ityp = 1
    this%iovar_dk(ityp)%name  = 'PA_R8'
    this%iovar_dk(ityp)%ndims = 1
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.  
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)

    ! 1d Site
    ityp = 2
    this%iovar_dk(ityp)%name  = 'SI_R8'
    this%iovar_dk(ityp)%ndims = 1
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)

    ! patch x ground
    ityp = 3
    this%iovar_dk(ityp)%name = 'PA_GRND_R8'
    this%iovar_dk(ityp)%ndims = 2
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)

    ! patch x size-class/pft
    ityp = 4
    this%iovar_dk(ityp)%name = 'PA_SCPF_R8'
    this%iovar_dk(ityp)%ndims = 2
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)

    ! site x ground
    ityp = 5
    this%iovar_dk(ityp)%name = 'SI_GRND_R8'
    this%iovar_dk(ityp)%ndims = 2
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)

    ! site x size-class/pft
    ityp = 6
    this%iovar_dk(ityp)%name = 'SI_SCPF_R8'
    this%iovar_dk(ityp)%ndims = 2
    allocate(this%iovar_dk(ityp)%dimsize(this%iovar_dk(ityp)%ndims))
    this%iovar_dk(ityp)%dimsize(:) = unset_int
    this%iovar_dk(ityp)%active = .false.
    nullify(this%iovar_dk(ityp)%dim1_ptr)
    nullify(this%iovar_dk(ityp)%dim2_ptr)


   
    
    
    
    return
  end subroutine init_iovar_dk_maps
  
  ! ===================================================================================
  
  subroutine set_dim_ptrs(this,dk_name,idim,dim_target)
    
    ! arguments
    class(fates_hio_interface_type) :: this
    character(len=*),intent(in)     :: dk_name
    integer,intent(in)              :: idim  ! dimension index
    type(iovar_dim_type),target     :: dim_target
    
    
    ! local
    integer                         :: ityp
    
    ityp = this%iotype_index(trim(dk_name))
    
    ! First check to see if the dimension is allocated
    if(this%iovar_dk(ityp)%ndims<idim)then
       write(fates_log(),*)'Trying to define dimension size to a dim-type structure'
       write(fates_log(),*)'but the dimension index does not exist'
       write(fates_log(),*)'type: ',dk_name,' ndims: ',this%iovar_dk(ityp)%ndims,' input dim:',idim
       stop
       !end_run
    end if
    
    if(idim==1) then
       this%iovar_dk(ityp)%dim1_ptr => dim_target
    elseif(idim==2) then
       this%iovar_dk(ityp)%dim2_ptr => dim_target
    end if

    ! With the map, we can set the dimension size
    this%iovar_dk(ityp)%dimsize(idim) = dim_target%ub - dim_target%lb + 1

    
    return
 end subroutine set_dim_ptrs
  
  ! ====================================================================================
  
  function iotype_index(this,iotype_name) result(ityp)
    
    ! argument
    class(fates_hio_interface_type) :: this
    character(len=*),intent(in)     :: iotype_name

    ! local
    integer :: ityp
    
    do ityp=1,n_iovar_dk
       if(trim(iotype_name).eq.trim(this%iovar_dk(ityp)%name))then
          return
       end if
    end do
    write(fates_log(),*) 'An IOTYPE THAT DOESNT EXIST WAS SPECIFIED'
    !end_run
    
  end function iotype_index
   
  ! =====================================================================================

  subroutine dim_init(this,iovar_dim,dim_name,nthreads,lb_in,ub_in)

    ! arguments
    class(fates_hio_interface_type) :: this
    type(iovar_dim_type),target     :: iovar_dim
    character(len=*),intent(in)     :: dim_name
    integer,intent(in)              :: nthreads
    integer,intent(in)              :: lb_in
    integer,intent(in)              :: ub_in

    allocate(iovar_dim%clump_lb(nthreads))
    allocate(iovar_dim%clump_ub(nthreads))
    
    iovar_dim%name = trim(dim_name)
    iovar_dim%lb = lb_in
    iovar_dim%ub = ub_in

    return
  end subroutine dim_init

  ! =====================================================================================

  subroutine set_dim_thread_bounds(this,iovar_dim,nc,lb_in,ub_in)

    class(fates_hio_interface_type) :: this
    type(iovar_dim_type),target     :: iovar_dim
    integer,intent(in)              :: nc    ! Thread index
    integer,intent(in)              :: lb_in
    integer,intent(in)              :: ub_in
    
    iovar_dim%clump_lb(nc) = lb_in
    iovar_dim%clump_ub(nc) = ub_in

    return
  end subroutine set_dim_thread_bounds

   ! ====================================================================================
   ! DEPRECATED, TRANSITIONAL OR FUTURE CODE SECTION
   ! ====================================================================================

   !subroutine set_fates_hio_str(tag,iotype_name,iostr_val)

!       ! Arguments
!       character(len=*),intent(in)           :: tag
!       character(len=*), optional,intent(in) :: iotype_name
!       integer, optional, intent(in)         :: iostr_val

!       ! local variables
!       logical              :: all_set
!       integer,  parameter  :: unset_int = -999
!       real(r8), parameter  :: unset_double = -999.9
!       integer              :: ityp, idim

!       select case (trim(tag))
!       case('flush_to_unset')
!          write(*,*) ''
!          write(*,*) 'Flushing FATES IO types prior to transfer from host'
!          do ityp=1,ubound(iovar_str,1)
!             iovar_str(ityp)%dimsize = unset_int
!             iovar_str(ityp)%active  = .false.
!          end do

!       case('check_allset')
!          do ityp=1,ubound(iovar_str,1)
!             write(*,*) 'Checking to see if ',iovar_str(ityp)%name,' IO communicators were sent to FATES'
!             if(iovar_str(ityp)%active)then
!                if(iovar_str(ityp)%offset .eq. unset_int) then
!                   write(*,*) 'FATES offset information of IO type:',iovar_str(ityp)%name
!                   write(*,*) 'was never set'
!                   ! end_run('MESSAGE')
!                end if
!                do idim=1,iovar_str(ityp)%ndims
!                   if(iovar_str(ityp)%dimsize(idim) .eq. unset_int) then
!                      write(*,*) 'FATES dimension information of IO type:',iovar_str(ityp)%name
!                      write(*,*) 'was never set'
!                      ! end_run('MESSAGE')
!                   end if
!                end do
!             end if
!          end do
!          write(*,*) 'Checked. All history IO specifications properly sent to FATES.'
!       case default

!          ! Must have two arguments if this is not a check or flush
!          if(present(iostr_val) .and. present(iotype_name))then
!
!             ! Tag in this case is dimsize or offset
!             select case (trim(tag))
!
!             case('offset')
!                ityp=iotype_index(trim(iotype_name))
!                iovar_str(ityp)%offset = iostr_val
!                write(*,*) 'Transfering offset for IOTYPE',iotype_name,' to FATES'

!             case('dimsize1')
!                ityp=iotype_index(trim(iotype_name))
!                iovar_str(ityp)%dimsize(1) = iostr_val
!                write(*,*) 'Transfering 1st dimension size for IOTYPE',iotype_name,' to FATES'

!             case('dimsize2')
!                ityp=iotype_index(trim(iotype_name))
!                if(ubound(iovar_str(ityp)%dimsize,1)==1)then
!                   write(fates_log(),*) 'Transfering second dimensional bound to unallocated space'
!                   write(fates_log(),*) 'type:',iotype_name
!                   ! end_run
!                end if
!                iovar_str(ityp)%dimsize(2) = iostr_val
!                write(*,*) 'Transfering 2nd dimension size for IOTYPE',iotype_name,' to FATES'

!             case('dimsize3')
!                ityp=iotype_index(trim(iotype_name))
!                if(ubound(iovar_str(ityp)%dimsize,1)<3)then
!                   write(fates_log(),*) 'Transfering third dimensional bound to unallocated space'
!                   write(fates_log(),*) 'type:',iotype_name
!                   ! end_run
!                end if
!                iovar_str(ityp)%dimsize(3) = iostr_val
!                write(*,*) 'Transfering 3rd dimension size for IOTYPE',iotype_name,' to FATES'

!             case default
!                write(*,*) 'IO parameter not recognized:',trim(tag)
!                ! end_run
!             end select
!          else
!             write(*,*) 'no value was provided for the tag'
!          end if
!
!       end select
!       return
!     end subroutine set_fates_hio_str



end module HistoryIOMod
