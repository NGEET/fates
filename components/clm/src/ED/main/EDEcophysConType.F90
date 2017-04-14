module EDEcophysConType

  !----------------------------------------------------
  ! ED ecophysiological constants 
  !----------------------------------------------------
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod , only      : errMsg => shr_log_errMsg
  use FatesGlobals, only      : endrun => fates_endrun
  use FatesGlobals, only      : fates_log

  use FatesHydraulicsMemMod , only : n_porous_media
  use FatesHydraulicsMemMod , only : porous_media
  use FatesHydraulicsMemMod , only : npool_tot
  use FatesHydraulicsMemMod , only : npool_leaf
  use FatesHydraulicsMemMod , only : npool_stem
  use FatesHydraulicsMemMod , only : npool_aroot
  use FatesHydraulicsMemMod , only : npool_troot

  use EDTypesMod, only : use_fates_plant_hydro
  
  !
  implicit none
  save
  private

  character(len=*), parameter, private :: sourcefile = &
        __FILE__

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: EDecophysconInit
  !
  ! !PUBLIC TYPES:
  type, public :: EDecophyscon_type
     real(r8), pointer :: max_dbh            (:) ! maximum dbh at which height growth ceases... 
     real(r8), pointer :: freezetol          (:) ! minimum temperature tolerance... 
     real(r8), pointer :: wood_density       (:) ! wood density  g cm^-3  ...  
     real(r8), pointer :: alpha_stem         (:) ! live stem turnover rate. y-1 
     real(r8), pointer :: hgt_min            (:) ! sapling height m 
     real(r8), pointer :: cushion            (:) ! labile carbon storage target as multiple of leaf pool. 
     real(r8), pointer :: leaf_stor_priority (:) ! leaf turnover vs labile carbon use prioritisation. ! (1=lose leaves, 0=use store). 
     real(r8), pointer :: leafwatermax       (:) ! amount of water allowed on leaf   surfaces
     real(r8), pointer :: rootresist         (:)
     real(r8), pointer :: soilbeta           (:)
     real(r8), pointer :: crown              (:) ! fraction of the height of the plant that is occupied by crown. For fire model. 
     real(r8), pointer :: bark_scaler        (:) ! scaler from dbh to bark thickness. For fire model. 
     real(r8), pointer :: crown_kill         (:) ! scaler on fire death. For fire model. 
     real(r8), pointer :: initd              (:) ! initial seedling density 
     real(r8), pointer :: sd_mort            (:) ! rate of death of seeds produced from reproduction. 
     real(r8), pointer :: seed_rain          (:) ! seeds that come from outside the gridbox.  
     real(r8), pointer :: BB_slope           (:) ! ball berry slope parameter
     real(r8), pointer :: root_long          (:) ! root longevity (yrs)
     real(r8), pointer :: clone_alloc        (:) ! fraction of carbon balance allocated to clonal reproduction. 
     real(r8), pointer :: seed_alloc         (:) ! fraction of carbon balance allocated to seeds. 
     real(r8), pointer :: sapwood_ratio      (:) ! amount of sapwood per unit leaf carbon and m height 



     ! pft parameters for plant hydraulics (PFT)
     real(r8), pointer :: wd                 (:)   ! wood density (distinct from wood_density for testing)    [g m-3]
     real(r8), pointer :: lma                (:)   ! leaf mass per area                                       [g m-2]   
                                                   ! ~ 90 for tropical angiosperms, cf Patino et al. 2012  
                                                   ! (existing param 'slatop' is biased high)
     real(r8), pointer :: n                  (:)   ! leaf nitrogen                                            [mg g-1]
     real(r8), pointer :: p                  (:)   ! leaf phosphorus                                          [mg g-1]
     real(r8), pointer :: ldmc               (:)   ! leaf dry matter content                                  [g g-1]
     real(r8), pointer :: lmv                (:)   ! leaf mass per volume                                     [g m-3]
     real(r8), pointer :: psi0               (:)   ! sapwood water potential at saturation                    [MPa]
     real(r8), pointer :: psicap             (:)   ! sapwood water potential at rwcft                         [MPa]
                                                   ! BOC...rhoc, rint_petiole, rint_jansenchoat, ccontent and (maybe) 
                                                   ! rs2 and rfrac_stem should really be global constants, not pft parameters
     real(r8), pointer :: rhoc               (:)   ! dry matter (or cell wall) density of wood                [g cm-3]  Siau 1984
     real(r8), pointer :: rint_petiole       (:)   ! radius of xylem conduits in petioles                     [um]
     real(r8), pointer :: rint_jansenchoat   (:)   ! average radius of xylem conduits where ks mmts were made [um]      
                                                   ! taken from choat & jansen XFT database for tropical angiosperms only
     real(r8), pointer :: Amaxh              (:)   ! light-saturated photosynthesis rate                      [umol m-2 s-1]
     real(r8), pointer :: rs2                (:)   ! mean absorbing fine root radius                          [m]       ~ 0.001 m?
     real(r8), pointer :: srl                (:)   ! specific root length                                     [m kg-1]  
                                                   ! ~ 15000 for tropical angiosperms, cf Metcalfe et al. 2008 Plant Soil Fig. 2b; 
                                                   ! rootdens = 500 kg m-3 is biased high by an order of magnitude 
                                                   ! (cf Comas et al. 2002 Oecologia); SPA rootdens implies a SRL of only 637 m kg-1.
     real(r8), pointer :: ccontent           (:)   ! carbon content (fraction of dry mass)                    [-] 
                                                   ! ~ 0.47 for tropical angiosperms, cf Thomas & Martin (2012) Forests
     real(r8), pointer :: latosa             (:)   ! leaf to sapwood area ratio                               [m2  m-2] 
                                                   ! ~ 8000 for tropical angiosperms, cf Patino et al. 2012
     real(r8), pointer :: rfrac_stem         (:)   ! fraction of total tree resistance (under well-watered conditions) 
                                                   ! from troot to canopy (i.e., aboveground) [-] ~ 0.625 for tropical angiosperms, 
                                                   ! cf BC re-analysis of Fisher et al. 2006
     real(r8), pointer :: rootshoot          (:)   ! root:shoot ratio (belowground-to-aboveground biomass)    [-]       
                                                   ! ~ 0.20 for tropical forests (see Houghton et al. 2001 Table 3,
                                                   !  Cairns et al. 1997 Table 2, Jackson et al. 1996 Table 3)
     real(r8), pointer :: avuln_gs           (:)   ! stomata PLC: vulnerability curve shape parameter         [-]
     real(r8), pointer :: p50_gs             (:)   ! stomata PLC: water potential at 50% loss of conductivity [Pa]    
     
     ! pft parameters for plant hydraulics (PFT x tissue type (leaf, stem, troot, aroot))
     real(r8), pointer :: kmax_node          (:,:) ! xylem PLC: maximum xylem hydraulic conductivity          [kg m-1 s-1 Pa-1]
     real(r8), pointer :: avuln_node         (:,:) ! xylem PLC: vulnerability curve shape parameter           [-]
     real(r8), pointer :: p50_node           (:,:) ! xylem PLC: water potential at 50% loss of conductivity   [Pa]    
     real(r8), pointer :: thetas_node        (:,:) ! P-V curve: saturated volumetric water content for node   [m3 m-3]
     real(r8), pointer :: epsil_node         (:,:) ! P-V curve: bulk elastic modulus                          [MPa]
     real(r8), pointer :: pinot_node         (:,:) ! P-V curve: osmotic potential at full turgor              [MPa]
     real(r8), pointer :: pitlp_node         (:,:) ! P-V curve: osmotic potential at turgor loss              [MPa]
     real(r8), pointer :: resid_node         (:,:) ! P-V curve: residual fraction                             [-]
     real(r8), pointer :: rwctlp_node        (:,:) ! P-V curve: total relative water content at turgor loss   [g or m3 H2O / g or m3 H2O, sat]
     real(r8), pointer :: fcap_node          (:,:) ! P-V curve: fraction of (1-resid_node) that is capillary in source [-]
     real(r8), pointer :: rwcft_node         (:,:) ! P-V curve: total RWC @ which elastic drainage begins     [-]
     real(r8), pointer :: rwccap_node        (:,:) ! P-V curve: total RWC @ which capillary reserves exhausted
     real(r8), pointer :: slp_node           (:,:) ! P-V curve: slope of capillary region of curve (sapwood only)
     real(r8), pointer :: intercept_node     (:,:) ! P-V curve: intercept of capillary region of curve (sapwood only)
     real(r8), pointer :: corrInt_node       (:,:) ! P-V curve: correction for nonzero psi0
     


  end type EDecophyscon_type

  type(EDecophyscon_type), public :: EDecophyscon ! ED ecophysiological constants structure
  !------------------------------------------------------------------------

  



contains

  !------------------------------------------------------------------------
  subroutine EDecophysconInit(EDpftvarcon_inst, numpft)
    !
    ! !USES:
    use EDPftvarcon, only : EDPftvarcon_type
    !
    ! !ARGUMENTS:
    type(EDpftVarCon_type) , intent(in) :: EDpftvarcon_inst
    integer                , intent(in) :: numpft
    !
    ! !LOCAL VARIABLES:
    integer :: m, ib, n, k
    !------------------------------------------------------------------------

    allocate( EDecophyscon%max_dbh            (0:numpft)); EDecophyscon%max_dbh            (:) = nan
    allocate( EDecophyscon%freezetol          (0:numpft)); EDecophyscon%freezetol          (:) = nan
    allocate( EDecophyscon%wood_density       (0:numpft)); EDecophyscon%wood_density       (:) = nan           
    allocate( EDecophyscon%alpha_stem         (0:numpft)); EDecophyscon%alpha_stem         (:) = nan             
    allocate( EDecophyscon%hgt_min            (0:numpft)); EDecophyscon%hgt_min            (:) = nan                
    allocate( EDecophyscon%cushion            (0:numpft)); EDecophyscon%cushion            (:) = nan                
    allocate( EDecophyscon%leaf_stor_priority (0:numpft)); EDecophyscon%leaf_stor_priority (:) = nan     
    allocate( EDecophyscon%leafwatermax       (0:numpft)); EDecophyscon%leafwatermax       (:) = nan           
    allocate( EDecophyscon%rootresist         (0:numpft)); EDecophyscon%rootresist         (:) = nan             
    allocate( EDecophyscon%soilbeta           (0:numpft)); EDecophyscon%soilbeta           (:) = nan               
    allocate( EDecophyscon%crown              (0:numpft)); EDecophyscon%crown              (:) = nan                  
    allocate( EDecophyscon%bark_scaler        (0:numpft)); EDecophyscon%bark_scaler        (:) = nan            
    allocate( EDecophyscon%crown_kill         (0:numpft)); EDecophyscon%crown_kill         (:) = nan             
    allocate( EDecophyscon%initd              (0:numpft)); EDecophyscon%initd              (:) = nan                  
    allocate( EDecophyscon%sd_mort            (0:numpft)); EDecophyscon%sd_mort            (:) = nan                
    allocate( EDecophyscon%seed_rain          (0:numpft)); EDecophyscon%seed_rain          (:) = nan              
    allocate( EDecophyscon%BB_slope           (0:numpft)); EDecophyscon%BB_slope           (:) = nan               
    allocate( EDecophyscon%root_long          (0:numpft)); EDecophyscon%root_long          (:) = nan                
    allocate( EDecophyscon%seed_alloc         (0:numpft)); EDecophyscon%seed_alloc         (:) = nan               
    allocate( EDecophyscon%clone_alloc        (0:numpft)); EDecophyscon%clone_alloc        (:) = nan              
    allocate( EDecophyscon%sapwood_ratio      (0:numpft)); EDecophyscon%sapwood_ratio      (:) = nan            

    do m = 0,numpft
       EDecophyscon%max_dbh(m)               = EDPftvarcon_inst%max_dbh(m)
       EDecophyscon%freezetol(m)             = EDPftvarcon_inst%freezetol(m)
       EDecophyscon%wood_density(m)          = EDPftvarcon_inst%wood_density(m)
       EDecophyscon%alpha_stem(m)            = EDPftvarcon_inst%alpha_stem(m)
       EDecophyscon%hgt_min(m)               = EDPftvarcon_inst%hgt_min(m)
       EDecophyscon%cushion(m)               = EDPftvarcon_inst%cushion(m)
       EDecophyscon%leaf_stor_priority(m)    = EDPftvarcon_inst%leaf_stor_priority(m)
       EDecophyscon%leafwatermax(m)          = EDPftvarcon_inst%leafwatermax(m)
       EDecophyscon%rootresist(m)            = EDPftvarcon_inst%rootresist(m)
       EDecophyscon%soilbeta(m)              = EDPftvarcon_inst%soilbeta(m)
       EDecophyscon%crown(m)                 = EDPftvarcon_inst%crown(m)
       EDecophyscon%bark_scaler(m)           = EDPftvarcon_inst%bark_scaler(m)
       EDecophyscon%crown_kill(m)            = EDPftvarcon_inst%crown_kill(m)
       EDecophyscon%initd(m)                 = EDPftvarcon_inst%initd(m)
       EDecophyscon%sd_mort(m)               = EDPftvarcon_inst%sd_mort(m)
       EDecophyscon%seed_rain(m)             = EDPftvarcon_inst%seed_rain(m)
       EDecophyscon%bb_slope(m)              = EDPftvarcon_inst%bb_slope(m)
       EDecophyscon%root_long(m)             = EDPftvarcon_inst%root_long(m)
       EDecophyscon%seed_alloc(m)            = EDPftvarcon_inst%seed_alloc(m)
       EDecophyscon%clone_alloc(m)           = EDPftvarcon_inst%clone_alloc(m)
       EDecophyscon%sapwood_ratio(m)         = EDPftvarcon_inst%sapwood_ratio(m)
    end do


    if (use_fates_plant_hydro) then
       allocate( EDecophyscon%wd               (0:numpft) );                  EDecophyscon%wd               (:) = nan
       allocate( EDecophyscon%lma              (0:numpft) );                  EDecophyscon%lma              (:) = nan
       allocate( EDecophyscon%n                (0:numpft) );                  EDecophyscon%n                (:) = nan
       allocate( EDecophyscon%p                (0:numpft) );                  EDecophyscon%p                (:) = nan
       allocate( EDecophyscon%ldmc             (0:numpft) );                  EDecophyscon%ldmc             (:) = nan
       allocate( EDecophyscon%lmv              (0:numpft) );                  EDecophyscon%lmv              (:) = nan
       allocate( EDecophyscon%psi0             (0:numpft) );                  EDecophyscon%psi0             (:) = nan
       allocate( EDecophyscon%psicap           (0:numpft) );                  EDecophyscon%psicap           (:) = nan
       allocate( EDecophyscon%rhoc             (0:numpft) );                  EDecophyscon%rhoc             (:) = nan
       allocate( EDecophyscon%rint_petiole     (0:numpft) );                  EDecophyscon%rint_petiole     (:) = nan
       allocate( EDecophyscon%rint_jansenchoat (0:numpft) );                  EDecophyscon%rint_jansenchoat (:) = nan
       allocate( EDecophyscon%Amaxh            (0:numpft) );                  EDecophyscon%Amaxh            (:) = nan
       allocate( EDecophyscon%rs2              (0:numpft) );                  EDecophyscon%rs2              (:) = nan
       allocate( EDecophyscon%srl              (0:numpft) );                  EDecophyscon%srl              (:) = nan
       allocate( EDecophyscon%ccontent         (0:numpft) );                  EDecophyscon%ccontent         (:) = nan
       allocate( EDecophyscon%latosa           (0:numpft) );                  EDecophyscon%latosa           (:) = nan
       allocate( EDecophyscon%rfrac_stem       (0:numpft) );                  EDecophyscon%rfrac_stem       (:) = nan
       allocate( EDecophyscon%rootshoot        (0:numpft) );                  EDecophyscon%rootshoot        (:) = nan
       allocate( EDecophyscon%avuln_gs         (0:numpft) );                  EDecophyscon%avuln_gs         (:) = nan
       allocate( EDecophyscon%p50_gs           (0:numpft) );                  EDecophyscon%p50_gs           (:) = nan
       
       allocate( EDecophyscon%kmax_node        (0:numpft,1:n_porous_media) ); EDecophyscon%kmax_node        (:,:) = nan
       allocate( EDecophyscon%avuln_node       (0:numpft,1:n_porous_media) ); EDecophyscon%avuln_node       (:,:) = nan
       allocate( EDecophyscon%p50_node         (0:numpft,1:n_porous_media) ); EDecophyscon%p50_node         (:,:) = nan
       allocate( EDecophyscon%thetas_node      (0:numpft,1:n_porous_media) ); EDecophyscon%thetas_node      (:,:) = nan
       allocate( EDecophyscon%epsil_node       (0:numpft,1:n_porous_media) ); EDecophyscon%epsil_node       (:,:) = nan
       allocate( EDecophyscon%pinot_node       (0:numpft,1:n_porous_media) ); EDecophyscon%pinot_node       (:,:) = nan
       allocate( EDecophyscon%pitlp_node       (0:numpft,1:n_porous_media) ); EDecophyscon%pitlp_node       (:,:) = nan
       allocate( EDecophyscon%resid_node       (0:numpft,1:n_porous_media) ); EDecophyscon%resid_node       (:,:) = nan
       allocate( EDecophyscon%rwctlp_node      (0:numpft,1:n_porous_media) ); EDecophyscon%rwctlp_node      (:,:) = nan
       allocate( EDecophyscon%fcap_node        (0:numpft,1:n_porous_media) ); EDecophyscon%fcap_node        (:,:) = nan
       allocate( EDecophyscon%rwcft_node       (0:numpft,1:n_porous_media) ); EDecophyscon%rwcft_node       (:,:) = nan
       allocate( EDecophyscon%rwccap_node      (0:numpft,1:n_porous_media) ); EDecophyscon%rwccap_node      (:,:) = nan
       allocate( EDecophyscon%slp_node         (0:numpft,1:n_porous_media) ); EDecophyscon%slp_node         (:,:) = nan
       allocate( EDecophyscon%intercept_node   (0:numpft,1:n_porous_media) ); EDecophyscon%intercept_node   (:,:) = nan
       allocate( EDecophyscon%corrInt_node     (0:numpft,1:n_porous_media) ); EDecophyscon%corrInt_node     (:,:) = nan

       ! ------------------------------------------------------------------------------------------------
       ! Until the hydraulics parameter are added to the parameter file, they need a location to be set.
       ! This happens here until further notice.
       ! ------------------------------------------------------------------------------------------------
       call SetHydraulicsTestingParams(EDecophyscon)

    end if

  end subroutine EDecophysconInit

  subroutine SetHydraulicsTestingParams(EDEcophyscon)
     
     ! Arguments
     type(EDecophyscon_type), intent(inout) :: EDEcophyscon
     
     write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
     call endrun(msg=errMsg(sourcefile, __LINE__))
     
     
  end subroutine SetHydraulicsTestingParams

end module EDEcophysConType
