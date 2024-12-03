Module TwoStreamMLPEMod

  ! This module holds the routines to calculate two-tream
  ! radiation scattering of vegetation in
  ! "M"ultiple "L"ayers with "P"arellel "E"lements "MLPE"
  !
  ! In summary,
  ! there may numerous canopy layers. In each canopy layer,
  ! plant media for different functional types are grouped
  ! so that they inhabit their own exclusive footprint
  ! within the layer. Within these exclusive functional
  ! columns, there are further sub-layer discretizations,
  ! which are organized by top-down integrated vegetation
  ! area index.
  !
  ! Note that there is a separate allocation and call
  ! sequence for each broad band.  In other words, the
  ! two_stream_type is instantiated for each broad band.
  !
  ! 
  !
  ! Assumptions: band index 1 = visible (vis)
  !                         2 = near infrared (nir)
  !                         3 = thermal (not used at the moment)
  !

  use shr_log_mod   , only: errMsg => shr_log_errMsg
  use shr_sys_mod   , only: shr_sys_abort
  use FatesConstantsMod, only : r8 => fates_r8
  use shr_infnan_mod, only : shr_infnan_isnan
  
  implicit none
  private

  real(r8),parameter :: nearzero = 1.e-20_r8
  logical, parameter :: debug=.true.
  real(r8), parameter :: unset_r8 = 1.e-36_r8
  real(r8), parameter :: unset_int = -999
  integer, parameter :: twostr_vis = 1         ! Named index of visible shortwave radiation
  integer, parameter :: twostr_nir = 2         ! Named index for near infrared shortwave radiation

  integer, parameter :: max_bands = 2          ! maximum number of bands (for scratch space)
  
  ! Allowable error, as a fraction of total incident for total canopy
  ! radiation balance checks

  real(r8), public, parameter :: rel_err_thresh = 1.e-6_r8
  real(r8), public, parameter :: area_err_thresh = rel_err_thresh*0.1_r8
  
  ! These are the codes for how the upper boundary is specified, normalized or absolute
  integer,public, parameter :: normalized_upper_boundary = 1
  integer,public, parameter :: absolute_upper_boundary   = 2
  
  integer :: log_unit ! fortran output unit for logging
  
  ! These are parameter constants, ie things that are specific to the plant material
  ! and radiation band.  Not all of these need to be used. 2-stream ultimately wants
  ! optical depth, scattering coefficient and backscatter fractions for diffuse and
  ! direct light. So there are various ways to get to these parameters, depending
  ! on the host model's available parameters.  The rho,tau,xl and clumping parameters
  ! are standard elm/clm parameters, and provided as a convenience.


  ! Snow optical parameter constants for visible (index=1) and NIR (index=2)

  real(r8), parameter :: betad_snow(1:2) = (/0.5, 0.5/)   ! Diffuse backscatter fraction    (CLM50 Tech Man)
  real(r8), parameter :: betab_snow(1:2) = (/0.5, 0.5/)   ! Beam backscatter fraction       (CLM50 Tech Man)  
  real(r8), parameter :: om_snow(1:2)    = (/0.8, 0.4/)   ! Scattering coefficient for snow (CLM50 Tech Man)
  !real(r8), parameter :: om_snow(1:2)    = (/0.85, 0.75/) ! Tarboton 95

  ! Cap the maximum optical depth. After 30 or so, its
  ! so close to zero, if the values get too large, then
  ! it will blow up the exponents and cause math problems

  real(r8), parameter :: kb_max = 30._r8
  

  ! For air, use a nominal values to prevent div0s
  ! the key is that vai = 0
  
  real(r8), parameter :: k_air = 0.5_r8  
  real(r8), parameter :: om_air  = 0.5_r8
  real(r8), parameter :: beta_air = 0.5_r8
  integer, public, parameter :: air_ft = 0 

  type, public :: rad_params_type

     ! From the parameter file
     real(r8), allocatable :: rhol(:,:)         ! leaf material reflectance:   (band x pft)
     real(r8), allocatable :: rhos(:,:)         ! stem material reflectance:   (band x pft)
     real(r8), allocatable :: taul(:,:)         ! leaf material transmittance: (band x pft)
     real(r8), allocatable :: taus(:,:)         ! stem material transmittance: (band x pft)
     real(r8), allocatable :: xl(:)             ! leaf/stem orientation (pft)
     real(r8), allocatable :: clumping_index(:) ! clumping index 0-1, when
                                                ! leaves stick together (pft)

     ! Derived parameters
     real(r8), allocatable :: phi1(:)       ! intermediate term for kd and kb
     real(r8), allocatable :: phi2(:)       ! intermediate term for kd and kb
     real(r8), allocatable :: avmu(:)       ! average "av" inverse optical depth "mu" per unit leaf and stem area
     real(r8), allocatable :: kd_leaf(:)    ! Mean optical depth per unit area leaves in diffuse
     real(r8), allocatable :: kd_stem(:)    ! Mean optical depth per unit area stems in diffuse
     real(r8), allocatable :: om_leaf(:,:)  ! Leaf scattering coefficient (band x pft)
     real(r8), allocatable :: om_stem(:,:)  ! Stem scattering coefficient (band x pft)
  end type rad_params_type

  type(rad_params_type),public :: rad_params


  ! Information describing the scattering elements
  ! that is based on "g"eometry, and independent of wavelength

  type scelg_type
     integer  :: pft      ! pft index
     real(r8) :: area     ! m2 col/m2 ground
     real(r8) :: lai      ! m2 of leaf area / m2 col
     real(r8) :: sai      ! m2 of stem area / m2 col
     real(r8) :: Kb       ! Optical depth of beam radiation
     real(r8) :: Kb_leaf  ! Optical depth of just leaves in beam radiation
     real(r8) :: Kd       ! Optical depth of diffuse radiation
     real(r8) :: area_squeeze ! This is the ratio of the element area to the
     ! the area of its constituents. Ideally this
     ! should be 1.0, but if the host model does not
     ! do a good job of filling up a canopy with 100% space,
     ! and instead is fractionally more than 100%, we must
     ! squeeze the area of 1 or more elements to get an exact
     ! space usage.
  end type scelg_type


  ! Information describing the scattering elemnets that
  ! is dependent on wavelengths, ie "b"ands (this is allocated for each broad band)

  type scelb_type

     ! Terms used in the final solution, also used for decomposing solution
     real(r8) :: Au      ! Compound intercept term
     real(r8) :: Ad      ! Compound intercept term
     real(r8) :: B1      ! Compound term w/ lambdas (operates on e^{av})
     real(r8) :: B2      ! Compound term w/ lambdas (operates on e^{-av})
     real(r8) :: lambda1_diff ! Compount term w/ B for diffuse forcing
     real(r8) :: lambda2_diff ! Compound term w/ B for diffuse forcing
     real(r8) :: lambda1_beam ! Compount term w/ B for beam forcing
     real(r8) :: lambda2_beam ! Compound term w/ B for beam forcing
     
     real(r8) :: a       ! Complex term operating on veg area index
     real(r8) :: om      ! scattering coefficient for media as a whole
     real(r8) :: betad   ! backscatter fraction of diffuse radiation for media as a whole
     real(r8) :: betab   ! backscatter fraction of beam radiation for media as a whole
     real(r8) :: Rbeam0  ! Normalized downwelling beam radiation at
                         ! top of the element (relative to downwelling atmospheric beam) [-]

  end type scelb_type


  type band_type

     type(scelb_type), pointer :: scelb(:,:)       ! array of scattering coefficients (layer, column)
                                                   ! can be sparse, will only solve indices up to
     integer                   :: ib               ! band index, should be consistent with rad_params
     real(r8)                  :: Rbeam_atm        ! Downwelling beam radiation from atmosphere [W/m2 ground]
     real(r8)                  :: Rdiff_atm        ! Downwelling diffuse radiation from atmosphere [W/m2 ground]
     real(r8)                  :: albedo_grnd_diff ! Ground albedo diffuse
     real(r8)                  :: albedo_grnd_beam ! Ground albedo direct 

  end type band_type


  ! This type contains the pre-processed scattering coefficients
  ! and routines.  This is the parent type that holds almost everything
  ! in the two-stream solver.
  ! The scelg structure describes the scattering elements, these are values
  ! that need to be defined by the ecosystem model, somewhat of
  ! an input to the solver. Since this is a Perfect Plasticity Approximation
  ! enabled system, we partition the scattering media into "columns" and "layers"
  ! Layers are canopy layers, think understory, mid-story and upper canopy. Columns
  ! are divisions of horizontal space, ie literal columns of space. The current
  ! implementation limits this space to media that has uniform scattering coefficients.
  ! So there could not be different PFTs in the same column, because they would undoubtedly
  ! have different joint scattering coefficients at different height levels in
  ! the column.  Therefore, every column is connected with a PFT.


  type, public :: twostream_type

     type(scelg_type), pointer :: scelg(:,:)    ! array of scattering elements (layer, column)
                                                ! can be sparse, will only solve indices up to
                                                ! n_lyr,n_col(n_lyr). This is for band (wavelength)
                                                ! independent information

     type(band_type), pointer  :: band(:)       ! Holds scattering coefficients for each band
                                                ! vis,nir,etc (nothing that emits though, no thermal)

     integer                   :: n_bands       ! number of bands (allocation size of band(:))
     integer                   :: n_lyr         ! number of (vertical) scattering element layers
     integer, allocatable      :: n_col(:)      ! number of (horizontal) scattering element columns per layer
     integer                   :: n_scel        ! total number of scattering elements
     logical                   :: force_prep    ! Some coefficients are only updated
                                                ! when the canopy composition changes, ie
                                                ! changes in leaf, stem or snow structure.
                                                ! If so, this sets to true, signalling that diffuse
                                                ! scattering coefficients should be updated.
                                                ! Otherwise, we only updated zenith dependent
                                                ! parameters on short sub-daily timesteps
     real(r8)                  :: frac_snow     ! Current mean snow-fraction of the canopy
     real(r8)                  :: frac_snow_old ! Previous mean snow-fraction of the canopy
     real(r8)                  :: cosz          ! Current cosine of the zenith angle
     
   contains

     procedure :: ZenithPrep     ! Update coefficients as zenith changes
     procedure :: CanopyPrep     ! Update coefficients as canopy changes
     procedure :: Solve          ! Perform the scattering solution
     procedure :: Dump           ! Dump out (print out) the site of interest
     procedure :: GetNSCel
     procedure :: AllocInitTwoStream
     procedure :: DeallocTwoStream
     
     procedure :: GetRdUp
     procedure :: GetRdDn
     procedure :: GetRb
     procedure :: GetAbsRad
     

  end type twostream_type

  public :: RadParamPrep
  public :: AllocateRadParams
  public :: TwoStreamLogInit
  
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  
contains

  subroutine TwoStreamLogInit(log_unit_in)
    integer,intent(in) :: log_unit_in
    
    log_unit = log_unit_in
    
  end subroutine TwoStreamLogInit

  subroutine endrun(msg) 

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    ! This subroutine was derived from CLM's
    ! endrun_vanilla() in abortutils.F90
    !
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: msg    ! string to be printed
    !-----------------------------------------------------------------------

    write(log_unit,*)'ENDRUN:', msg
    call shr_sys_abort()

  end subroutine endrun

  
  ! ===============================================================================================
  
  subroutine AllocInitTwoStream(this,band_indices,ncan,ncol)

    class(twostream_type) :: this
    integer               :: band_indices(:)
    integer               :: ncan
    integer               :: ncol

    integer :: nbands
    integer :: ib

    nbands = ubound(band_indices,1)

    allocate(this%n_col(ncan))
    allocate(this%scelg(ncan,ncol))
    allocate(this%band(nbands))

    this%n_col(1:ncan)    = unset_int
    this%n_bands          = nbands
    this%n_lyr            = ncan
    this%frac_snow        = unset_r8
    this%frac_snow_old    = unset_r8

    do ib = 1,nbands

       allocate(this%band(ib)%scelb(ncan,ncol))
       this%band(ib)%albedo_grnd_diff = unset_r8
       this%band(ib)%albedo_grnd_beam = unset_r8
       this%band(ib)%ib = band_indices(ib)

    end do



    return
  end subroutine AllocInitTwoStream

  ! ===============================================================================================

  subroutine DeallocTwoStream(this)

    class(twostream_type) :: this

    integer               :: nbands
    integer               :: ib

    nbands = ubound(this%band,1)

    deallocate(this%scelg)
    deallocate(this%n_col)
    do ib = 1,nbands
       deallocate(this%band(ib)%scelb)
    end do
    deallocate(this%band)

    return
  end subroutine DeallocTwoStream

  ! ===============================================================================================

  subroutine AllocateRadParams(n_pft,n_bands)

    integer,intent(in) :: n_pft
    integer,intent(in) :: n_bands

    ! Include the zeroth pft index for air
    
    allocate(rad_params%rhol(n_bands,n_pft))
    allocate(rad_params%rhos(n_bands,n_pft))
    allocate(rad_params%taul(n_bands,n_pft))
    allocate(rad_params%taus(n_bands,n_pft))
    allocate(rad_params%xl(n_pft))
    allocate(rad_params%clumping_index(n_pft))

    allocate(rad_params%phi1(n_pft))
    allocate(rad_params%phi2(n_pft))
    allocate(rad_params%avmu(n_pft))
    allocate(rad_params%kd_leaf(n_pft))
    allocate(rad_params%kd_stem(n_pft))
    allocate(rad_params%om_leaf(n_bands,n_pft))
    allocate(rad_params%om_stem(n_bands,n_pft))
        
  end subroutine AllocateRadParams

  ! ================================================================================================

  function GetRdDn(this,ican,icol,ib,vai) result(r_diff_dn)

    class(twostream_type) :: this
    real(r8),intent(in)   :: vai
    integer,intent(in)    :: ican
    integer,intent(in)    :: icol
    integer,intent(in)    :: ib
    real(r8)              :: r_diff_dn

    ! Rdn = Ad e−(Kbv) + Re + λ1 B2 e^(av) + λ2 B1 e^(−av)

    associate(scelb => this%band(ib)%scelb(ican,icol), &
         scelg => this%scelg(ican,icol) )

      r_diff_dn = this%band(ib)%Rbeam_atm*( &
           scelb%Ad*exp(-scelg%Kb*vai) + &
           scelb%B2*scelb%lambda1_beam*exp(scelb%a*vai) + &
           scelb%B1*scelb%lambda2_beam*exp(-scelb%a*vai)) + &
           this%band(ib)%Rdiff_atm*( & 
           scelb%B2*scelb%lambda1_diff*exp(scelb%a*vai) + &
           scelb%B1*scelb%lambda2_diff*exp(-scelb%a*vai))

      if(debug)then
         ! if(isnan(r_diff_dn))then  !RGK: NVHPC HAS A BUG IN THIS INTRINSIC (01-2024)
         ! if(r_diff_dn /= r_diff_dn) then
         if(shr_infnan_isnan(r_diff_dn)) then
            write(log_unit,*)"GETRDN"
            write(log_unit,*)scelg%Kb
            write(log_unit,*)scelb%a
            write(log_unit,*)vai
            write(log_unit,*)scelb%Ad
            write(log_unit,*)scelb%B1,scelb%B2
            write(log_unit,*)scelb%lambda1_beam,scelb%lambda2_beam
            write(log_unit,*)scelb%lambda1_diff,scelb%lambda2_diff
            write(log_unit,*)this%band(ib)%Rbeam_atm
            write(log_unit,*)this%band(ib)%Rdiff_atm
            write(log_unit,*)exp(-scelg%Kb*vai)
            write(log_unit,*)exp(scelb%a*vai)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
      end if
      
    end associate
  end function GetRdDn

  function GetRdUp(this,ican,icol,ib,vai) result(r_diff_up)

    class(twostream_type) :: this
    real(r8),intent(in)   :: vai
    integer,intent(in)    :: ican
    integer,intent(in)    :: icol
    integer,intent(in)    :: ib
    real(r8)              :: r_diff_up

    ! Rup = Au e−(Kbv) + Re + λ1 B1 e^(av) + λ2 B2 e^(−av)

    associate(scelb => this%band(ib)%scelb(ican,icol), &
         scelg => this%scelg(ican,icol) )

      r_diff_up = this%band(ib)%Rbeam_atm*( &
           scelb%Au*exp(-scelg%Kb*vai) + &
           scelb%B1*scelb%lambda1_beam*exp(scelb%a*vai) + &
           scelb%B2*scelb%lambda2_beam*exp(-scelb%a*vai)) + &
           this%band(ib)%Rdiff_atm*( & 
           scelb%B1*scelb%lambda1_diff*exp(scelb%a*vai) + &
           scelb%B2*scelb%lambda2_diff*exp(-scelb%a*vai)) 
           
    end associate
  end function GetRdUp

  function GetRb(this,ican,icol,ib,vai) result(r_beam_dn)

    class(twostream_type) :: this
    real(r8),intent(in)   :: vai
    integer,intent(in)    :: ican
    integer,intent(in)    :: icol
    integer,intent(in)    :: ib
    real(r8)              :: r_beam_dn

    r_beam_dn = this%band(ib)%Rbeam_atm * &
         this%band(ib)%scelb(ican,icol)%Rbeam0*exp(-this%scelg(ican,icol)%Kb*vai)

  end function GetRb

  subroutine GetAbsRad(this,ican,icol,ib,vai_top,vai_bot, &
       Rb_abs,Rd_abs,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem,R_abs_snow,leaf_sun_frac)

    ! This routine is used to help decompose radiation scattering
    ! and return the amount of absorbed radiation.  The canopy layer and column
    ! index identify the element of interest. The other arguments are the upper and
    ! lower bounds within the element over which to evaluate absorbed radiation.
    ! The assumption is that the vegetation area index is zero at the top of the
    ! element, and increases going downwards.  As with all assumptions in this
    ! module, the scattering parameters are uniform within the element itself,
    ! which includes an assumption of the leaf/stem proportionality.
    ! ---------------------------------------------------------------------------
    ! Solution for radiative intensity of diffuse up and down at tai=v
    ! Rup = Au e−(Kbv) + Re + λ1 B1 e^(av) + λ2 B2 e^(−av)
    ! Rdn = Ad e−(Kbv) + Re + λ1 B2 e^(av) + λ2 B1 e^(−av)
    ! ---------------------------------------------------------------------------

    ! Arguments
    class(twostream_type) :: this
    integer,intent(in)    :: ican
    integer,intent(in)    :: icol
    integer,  intent(in)  :: ib            ! broad band index
    real(r8), intent(in)  :: vai_top       ! veg area index (from the top of element) to start
    real(r8), intent(in)  :: vai_bot       ! veg area index (from the top of element) to finish
    real(r8), intent(out) :: Rb_abs        ! total absorbed beam radiation [W/m2 ground]
    real(r8), intent(out) :: Rd_abs        ! total absorbed diffuse radiation [W/m2 ground]
    real(r8), intent(out) :: Rb_abs_leaf   ! Absorbed beam radiation from leaves [W/m2 ground]
    real(r8), intent(out) :: Rd_abs_leaf   ! Absorbed diff radiation from leaves [W/m2 ground]
    real(r8), intent(out) :: R_abs_stem    ! Absorbed beam+diff radiation stems  [W/m2 ground]
    real(r8), intent(out) :: R_abs_snow    ! Absorbed beam+diff radiation snow   [W/m2 ground]
    real(r8), intent(out) :: leaf_sun_frac ! Fraction of leaves in the interval exposed
                                           ! to sunlight

    real(r8)              :: dvai,dlai     ! Amount of VAI and LAI in this interval [m2/m2]
    real(r8)              :: Rd_net        ! Difference in diffuse radiation at upper and lower boundaries [W/m2]
    real(r8)              :: Rb_net        ! Difference in beam radiation at upper and lower boundaries [W/m2]
    real(r8)              :: vai_max       ! total integrated (leaf+stem) area index of the current element
    real(r8)              :: frac_abs_snow ! fraction of radiation absorbed by snow
    real(r8)              :: diff_wt_leaf  ! diffuse absorption weighting for leaves
    real(r8)              :: diff_wt_stem  ! diffuse absorption weighting for stems
    real(r8)              :: beam_wt_leaf  ! beam absorption weighting for leaves
    real(r8)              :: beam_wt_stem  ! beam absorption weighting for stems
    real(r8)              :: lai_bot,lai_top
    
    associate(scelb => this%band(ib)%scelb(ican,icol), &
         scelg => this%scelg(ican,icol), &
         ft => this%scelg(ican,icol)%pft )

      ! If this is air, trivial solutions
      if(ft==air_ft) then
         Rb_abs        = 0._r8
         Rd_abs        = 0._r8
         Rb_abs_leaf   = 0._r8
         Rd_abs_leaf   = 0._r8
         R_abs_stem    = 0._r8
         R_abs_snow    = 0._r8
         leaf_sun_frac = 0._r8
         return
      end if
      
      ! The total vegetation area index of the element
      vai_max = scelg%lai +  scelg%sai

      dvai = vai_bot - vai_top

      lai_top = vai_top*scelg%lai/( scelg%lai+ scelg%sai)
      lai_bot = vai_bot*scelg%lai/( scelg%lai+ scelg%sai)
      dlai    = dvai *  scelg%lai/( scelg%lai+ scelg%sai)

      
      if(dlai>nearzero)then
         leaf_sun_frac = max(0.001_r8,min(0.999_r8,scelb%Rbeam0/(dlai*scelg%Kb_leaf/rad_params%clumping_index(ft)) &
              *(exp(-scelg%Kb_leaf*lai_top) - exp(-scelg%Kb_leaf*lai_bot))))
      else
         leaf_sun_frac = 0001._r8
      end if
         
      !leaf_sun_frac = max(0.001_r8,min(0.999_r8,scelb%Rbeam0/(dvai*scelg%Kb/rad_params%clumping_index(ft)) &
      !     *(exp(-scelg%Kb*vai_top) - exp(-scelg%Kb*vai_bot))))

      
      if(debug) then
         if(leaf_sun_frac>1.0_r8 .or. leaf_sun_frac<0._r8) then
            write(log_unit,*)"impossible leaf sun fraction"
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
      end if

      ! We have to disentangle the absorption between leaves and stems, we give them both
      ! a weighting fraction of total absorption of  area*K*(1-om)

      frac_abs_snow = this%frac_snow*(1._r8-om_snow(ib)) / (1._r8-scelb%om)

      diff_wt_leaf =  scelg%lai*(1._r8-rad_params%om_leaf(ib,ft))*rad_params%Kd_leaf(ft)
      diff_wt_stem =  scelg%sai*(1._r8-rad_params%om_stem(ib,ft))*rad_params%Kd_stem(ft)

      beam_wt_leaf =  scelg%lai*(1._r8-rad_params%om_leaf(ib,ft))*scelg%Kb_leaf
      beam_wt_stem =  scelg%sai*(1._r8-rad_params%om_stem(ib,ft))*1._r8

      ! Mean element transmission coefficients adding snow scattering

      if(debug) then
         if( (vai_bot-vai_max)>rel_err_thresh)then
            write(log_unit,*)"During decomposition of the 2-stream radiation solution"
            write(log_unit,*)"A vegetation area index (VAI) was requested in GetAbsRad()"
            write(log_unit,*)"that is larger than the total integrated VAI of the "
            write(log_unit,*)"computation element of interest."
            write(log_unit,*)"vai_max: ",vai_max
            write(log_unit,*)"vai_bot: ",vai_bot
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         if( (vai_bot-vai_top)<-rel_err_thresh ) then
            write(log_unit,*)"During decomposition of the 2-stream radiation solution"
            write(log_unit,*)"the vegetation area index at the lower position was set"
            write(log_unit,*)"as greater than the value at the upper position."
            write(log_unit,*)"vai_max: ",vai_max
            write(log_unit,*)"vai_bot: ",vai_bot
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
      end if

      ! Amount of absorbed radiation is retrieved by doing an energy
      ! balance on this boundaries over the depth of interest (ie net)
      ! Result is Watts / m2 of the element's area footprint NOT
      ! per m2 of tissue (at least not in this step)

      Rb_net = this%GetRb(ican,icol,ib,vai_top)-this%GetRb(ican,icol,ib,vai_bot)

      Rd_net = (this%GetRdDn(ican,icol,ib,vai_top) - this%GetRdDn(ican,icol,ib,vai_bot)) + &
           (this%GetRdUp(ican,icol,ib,vai_bot) - this%GetRdUp(ican,icol,ib,vai_top))

      ! The net beam radiation includes that which is absorbed, but also,
      ! that which is re-scattered, the re-scattered acts as a source
      ! to the net diffuse balance and adds to the absorbed, and a sink
      ! on the beam absorbed term.

      Rb_abs = Rb_net * (1._r8-this%band(ib)%scelb(ican,icol)%om)
      Rd_abs = Rd_net +  Rb_net * this%band(ib)%scelb(ican,icol)%om


      Rb_abs_leaf = (1._r8-frac_abs_snow)*Rb_abs * beam_wt_leaf / (beam_wt_leaf+beam_wt_stem)
      Rd_abs_leaf = (1._r8-frac_abs_snow)*Rd_abs * diff_wt_leaf / (diff_wt_leaf+diff_wt_stem)

      R_abs_snow = (Rb_abs+Rd_abs)*frac_abs_snow

      R_abs_stem = (1._r8-frac_abs_snow)* &
           (Rb_abs*beam_wt_stem / (beam_wt_leaf+beam_wt_stem) + &
           Rd_abs*diff_wt_stem / (diff_wt_leaf+diff_wt_stem))




    end associate
    return
  end subroutine GetAbsRad

  ! ================================================================================================
  
  subroutine Dump(this,ib,lat,lon)

    ! Dump out everything we know about these two-stream elements

    class(twostream_type) :: this
    integer,intent(in)    :: ib
    real(r8),optional,intent(in)   :: lat
    real(r8),optional,intent(in)   :: lon
    integer  :: ican
    integer  :: icol

    write(log_unit,*) 'Dumping Two-stream elements for band ', ib
    write(log_unit,*)
    write(log_unit,*) 'rbeam atm: ',this%band(ib)%Rbeam_atm
    write(log_unit,*) 'rdiff_atm: ',this%band(ib)%Rdiff_atm
    write(log_unit,*) 'alb grnd diff: ',this%band(ib)%albedo_grnd_diff
    write(log_unit,*) 'alb grnd beam: ',this%band(ib)%albedo_grnd_beam
    write(log_unit,*) 'cosz: ',this%cosz
    write(log_unit,*) 'snow fraction: ',this%frac_snow
    if(present(lat)) write(log_unit,*) 'lat: ',lat
    if(present(lon)) write(log_unit,*) 'lon: ',lon

    do_can: do ican = 1,this%n_lyr
       do_col: do icol = 1,this%n_col(ican)
          associate(scelg => this%scelg(ican,icol), &
                    scelb => this%band(ib)%scelb(ican,icol))
            write(log_unit,*) '--',ican,icol,'--'
            write(log_unit,*) 'pft:',scelg%pft
            write(log_unit,*) 'area: ',scelg%area
            write(log_unit,*) 'lai,sai: ',scelg%lai,scelg%sai
            write(log_unit,*) 'Kb: ',scelg%Kb
            write(log_unit,*) 'Kb leaf: ',scelg%Kb_leaf
            write(log_unit,*) 'Kd: ',scelg%Kd
            write(log_unit,*) 'Rb0: ',scelb%Rbeam0
            write(log_unit,*) 'om: ',scelb%om
            write(log_unit,*) 'betad: ',scelb%betad
            write(log_unit,*) 'betab:',scelb%betab
            write(log_unit,*) 'a: ',scelb%a
            this%band(ib)%Rbeam_atm = 1.0_r8
            this%band(ib)%Rdiff_atm = 1.0_r8
            write(log_unit,*)'RDiff Down @ bottom: ',this%GetRdDn(ican,icol,ib,scelg%lai+scelg%sai)
            write(log_unit,*)'RDiff Up @ bottom: ',this%GetRdUp(ican,icol,ib,scelg%lai+scelg%sai)
            write(log_unit,*)'Rbeam @ bottom: ',this%GetRb(ican,icol,ib,scelg%lai+scelg%sai)
          end associate
       end do do_col
    end do do_can

  end subroutine Dump


  ! ================================================================================================

  subroutine RadParamPrep()

    integer  :: ft
    integer  :: nbands
    integer  :: numpft
    integer  :: ib

    numpft = ubound(rad_params%om_leaf,2)
    nbands = ubound(rad_params%om_leaf,1)
    
    do ft = 1,numpft

       ! The non-band specific parameters here will be re-derived for each
       ! band, which is inefficient, however this is an incredibly cheap
       ! routine to begin with, its only called during initialization, so
       ! just let it go, dont worry about it.

       if(rad_params%xl(ft)<-0.4_r8 .or. rad_params%xl(ft)>0.6_r8) then
          write(log_unit,*)"Leaf orientation factors (xl) should be between -0.4 and 0.6"
          write(log_unit,*)"ft: ",ft,"xl: ",rad_params%xl(ft)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       ! There is a singularity of leaf orientation is exactly 0
       ! phi1 = 0.5
       ! phi2 = 0.0
       ! avmu = 1/0  (1 - 0.5/0 * ln(0.5/0.5) ) but the limit approaches 1
       ! a value of 0.0001 does not break numerics and generates an avmu of nearly 1

       if( abs(rad_params%xl(ft)) <0.0001) rad_params%xl(ft)=0.0001_r8
       
       ! There must be protections on xl to prevent div0 and other weirdness
       rad_params%phi1(ft) = 0.5_r8 - 0.633_r8*rad_params%xl(ft) - 0.330_r8*rad_params%xl(ft)*rad_params%xl(ft)
       rad_params%phi2(ft) = 0.877_r8 * (1._r8 - 2._r8*rad_params%phi1(ft)) !0 = horiz leaves, 1 - vert leaves.

       ! Eq. 3.4 CLM50 Tech Man
       rad_params%avmu(ft) = (1._r8/rad_params%phi2(ft))* &
            (1._r8-(rad_params%phi1(ft)/rad_params%phi2(ft))* &
            log((rad_params%phi2(ft)+rad_params%phi1(ft))/rad_params%phi1(ft)))
       
       do ib = 1, nbands
          rad_params%Kd_leaf(ft) = rad_params%clumping_index(ft)/rad_params%avmu(ft)
          rad_params%Kd_stem(ft) = 1._r8 
          
          rad_params%om_leaf(ib,ft) = rad_params%rhol(ib,ft) + rad_params%taul(ib,ft)
          rad_params%om_stem(ib,ft) = rad_params%rhos(ib,ft) + rad_params%taus(ib,ft)

          if( rad_params%om_leaf(ib,ft) > 0.99_r8 ) then
             write(log_unit,*) "In: TwoStreamMLPEMod.F90:RadParamPrep()"
             write(log_unit,*) "An extremely high leaf scattering coefficient was generated:"
             write(log_unit,*) "om = tau + rho"
             write(log_unit,*) "band = ",ib
             write(log_unit,*) "pft = ",ft
             write(log_unit,*) "om_leaf = ",rad_params%om_leaf(ib,ft)
             write(log_unit,*) "rhol = ",rad_params%rhol(ib,ft)
             write(log_unit,*) "taul = ",rad_params%taul(ib,ft)
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
           if( rad_params%om_stem(ib,ft) > 0.99_r8 ) then
             write(log_unit,*) "In: TwoStreamMLPEMod.F90:RadParamPrep()"
             write(log_unit,*) "An extremely high stem scattering coefficient was generated:"
             write(log_unit,*) "om = tau + rho"
             write(log_unit,*) "band = ",ib
             write(log_unit,*) "pft = ",ft
             write(log_unit,*) "om_stem = ",rad_params%om_stem(ib,ft)
             write(log_unit,*) "rhos = ",rad_params%rhos(ib,ft)
             write(log_unit,*) "taus = ",rad_params%taus(ib,ft)
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          
       end do
       
    end do

    return
  end subroutine RadParamPrep

  ! ================================================================================================


  ! ================================================================================================

  subroutine CanopyPrep(this,frac_snow)

    ! Pre-process things that change with canopy-geometry or snow cover
    ! We try to only run this when necessary. For instance we only
    ! run this when the canopy vegetation composition changes, or
    ! when the amount of snow-cover changes.  

    class(twostream_type) :: this

    real(r8)              :: frac_snow   ! The fraction (in terms of vegetation area index)
    ! of vegetation covered with snow

    ! But we check if the snow conditions
    ! change during the high frequency calls
    ! as well.
    integer :: ib     ! The band of interest
    integer :: ican  ! scattering element canopy layer index (top down)
    integer :: icol  ! scattering element column
    real(r8) :: rho  ! element mean material reflectance
    real(r8) :: tau  ! element mean material transmittance
    real(r8) :: vai  ! vegetation area index lai+sai
    real(r8) :: om_veg     ! scattering coefficient for vegetation (no snow)
    real(r8) :: betad_veg  ! diffuse backscatter for vegetation (no snow)
    real(r8) :: betad_om   ! multiplication of diffuse backscatter and reflectance
    real(r8) :: area_check ! Checks to make sure each layer has 100% coverage
    real(r8) :: a2         ! The "a" term squared
    
    this%frac_snow = frac_snow

    if(.not.this%force_prep) then
       if(abs(this%frac_snow-this%frac_snow_old)<nearzero) then
          this%frac_snow_old = this%frac_snow
          return
       end if
    end if

    this%frac_snow_old = this%frac_snow

    do_can: do ican = 1,this%n_lyr

       area_check = 0._r8
       do_col: do icol = 1,this%n_col(ican)

          associate(lai => this%scelg(ican,icol)%lai, &
               sai => this%scelg(ican,icol)%sai, &
               ft  => this%scelg(ican,icol)%pft, &
               scelg => this%scelg(ican,icol))

            vai = lai + sai

            ! Mean element transmission coefficients w/o snow effects

            if(ft==air_ft) then
               scelg%Kd = k_air
            else
               if(debug)then
                  if(vai<nearzero)then
                     write(log_unit,*)"zero vai in non-air element?"
                     write(log_unit,*) lai,sai,ican,icol,this%n_lyr,this%n_col(ican)
                     write(log_unit,*)"TwoStreamMLPEMod.F90:CanopyPrep"
                     call endrun(msg=errMsg(sourcefile, __LINE__))
                  end if
               end if
               scelg%Kd =  (lai * rad_params%Kd_leaf(ft) + &
                    sai * rad_params%Kd_stem(ft))/vai
            end if

            area_check = area_check + scelg%area

            do_bands: do ib = 1, this%n_bands

               associate(scelb => this%band(ib)%scelb(ican,icol))

                 if (ft==air_ft) then

                    scelb%om = om_air
                    scelb%betad = beta_air

                 else

                    ! Material reflectance (weighted average of leaf stem and snow)

                    ! Eq. 3.11 and 3.12 ClM5.0 Tech Man
                    om_veg  =  (lai*rad_params%om_leaf(ib,ft) + &
                         sai*rad_params%om_stem(ib,ft))/vai

                    ! Eq. 3.5 ClM5.0 Tech Man
                    scelb%om = this%frac_snow*om_snow(ib) + (1._r8-this%frac_snow)*om_veg

                    ! Diffuse backscatter, taken from G. Bonan's code

                    rho = (lai * rad_params%rhol(ib,ft) + &
                         sai * rad_params%rhos(ib,ft))/vai
                    tau = (lai * rad_params%taul(ib,ft) + &
                         sai * rad_params%taus(ib,ft))/vai

                    ! Eq 3.13 from CLM5.0 Tech Man
                    betad_veg  = 0.5_r8 / scelb%om * &
                         ( scelb%om + (rho-tau) * ((1._r8+rad_params%xl(ft))/2._r8)**2._r8 )

                    ! Eq. 3.6 from CLM5.0 Tech Man
                    betad_om = betad_veg*om_veg*(1._r8-this%frac_snow) + &
                         om_snow(ib)*betad_snow(ib)*this%frac_snow

                    scelb%betad = betad_om / scelb%om

                    if(debug)then
                       !if(isnan(scelb%betad))then !RGK: NVHPC HAS A BUG IN THIS INTRINSIC (01-2024)
                       !if(scelb%betad /= scelb%betad) then
                       if(shr_infnan_isnan(scelb%betad))then
                          write(log_unit,*)"nans in canopy prep"
                          write(log_unit,*) ib,ican,icol,ft
                          write(log_unit,*) scelb%betad,scelb%om,lai,sai
                          write(log_unit,*) this%frac_snow,om_snow(ib),vai,om_veg
                          write(log_unit,*)"TwoStreamMLPEMod.F90:CanopyPrep"
                          call endrun(msg=errMsg(sourcefile, __LINE__))
                       end if
                    end if                    
                    
                 end if

                 a2 = scelg%Kd*scelg%Kd*(1._r8-scelb%om)*(1._r8-scelb%om+2._r8*scelb%om*scelb%betad)
                 if(a2<0._r8) then
                    write(log_unit,*)'a^2 is less than zero'
                    call endrun(msg=errMsg(sourcefile, __LINE__))
                 end if
                 
                 ! We also have to avoid singularities, see Ad and Au below,
                 ! where a^2-Kb^2 is in the denominator
                 
                 scelb%a  = sqrt(a2)
                 
               end associate
            end do do_bands
          end associate
       end do do_col

       ! RE-ENABLE THIS CHECK WHEN FATES IS BETTER AT CONSERVING AREA!!
       if(.false.)then
          !if( abs(area_check-1._r8) > 10._r8*area_err_thresh  )then
          write(log_unit,*)"Only a partial canopy was specified"
          write(log_unit,*)"Scattering elements must constitute 100% of the ground cover."
          write(log_unit,*)"for open spaces, create an air element with the respective area."
          write(log_unit,*)"total area (out of 1): ",area_check,ican
          write(log_unit,*)"layer: ",ican," of: ",this%n_lyr
          do icol = 1,this%n_col(ican)
             write(log_unit,*)this%scelg(ican,icol)%area,this%scelg(ican,icol)%pft
          end do
          write(log_unit,*)"TwoStreamMLPEMod.F90:CanopyPrep"
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
            
    end do do_can

    return
  end subroutine CanopyPrep

  ! ================================================================================================

  subroutine ZenithPrep(this,cosz_in)

    ! Pre-process things that change with the zenith angle
    ! i.e. the beam optical properties

    ! Important !!!!
    ! This should always be called after CanopyPrep() has been
    ! called.  This routine relies on the results of that routine
    ! notably the scattering coefficient "om".

    class(twostream_type) :: this
    
    real(r8),intent(in)   :: cosz_in ! Un-protected cosine of the zenith angle

    real(r8) :: cosz                          ! the near-zero protected cosz
    integer :: ican                           ! scattering element canopy layer index (top down)
    integer :: icol                           ! scattering element column
    integer :: ib                             ! band index, matches indexing of rad_params
    integer :: ib2                            ! band inner loop index while testing for singularity
    real(r8) :: asu                           ! single scattering albedo
    real(r8) :: gdir
    real(r8) :: tmp0,tmp1,tmp2
    real(r8) :: betab_veg                     ! beam backscatter for vegetation (no snow)
    real(r8) :: betab_om                      ! multiplication of beam backscatter and reflectance
    real(r8) :: om_veg                        ! scattering coefficient for vegetation (no snow)
    real(r8) :: Kb_sing(max_bands)            ! the KB_leaf that would generate a singularity
                                              ! with the scelb%a parameter
    real(r8) :: Kb_stem                       ! actual optical depth of stem with not planar geometry effects
                                              ! usually the base value
    real(r8), parameter :: Kb_stem_base = 1.0_r8
    real(r8), parameter :: sing_tol = 0.01_r8 ! allowable difference between
                                              ! the Kb_leaf that creates
                                              ! a singularity and the actual
    logical :: is_sing   ! use this to control if we are actively trying to remove a singularity
    integer :: iter_sing ! iterator check to ensure we don't try to fix a singularity indefinitely
    real(r8) :: Kb_eff   ! When testing for singularity, this is either the stem or stem and leaf optical depth
    
    if( (cosz_in-1.0) > nearzero ) then
       write(log_unit,*)"The cosine of the zenith angle cannot exceed 1"
       write(log_unit,*)"cosz: ",cosz
       write(log_unit,*)"TwoStreamMLPEMod.F90:ZenithPrep"
       call endrun(msg=errMsg(sourcefile, __LINE__))
    elseif(cosz_in<0._r8)then
       write(log_unit,*)"The cosine of the zenith angle should not be less than zero"
       write(log_unit,*)"It can be exactly zero, but not less than"
       write(log_unit,*)"cosz: ",cosz
       write(log_unit,*)"TwoStreamMLPEMod.F90:ZenithPrep"
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
       
    cosz = max(0.001,cosz_in)

    this%cosz = cosz
    
    do_ican: do ican = 1,this%n_lyr
       do_ical: do icol = 1,this%n_col(ican)

          associate(ft => this%scelg(ican,icol)%pft, &
               scelg => this%scelg(ican,icol))

            if(ft==air_ft)then
               ! Simple provisions for a ghost element (air)
               scelg%Kb_leaf = k_air
               scelg%Kb = k_air
            else
               gdir = rad_params%phi1(ft) + rad_params%phi2(ft) * cosz

               Kb_stem = Kb_stem_base
               
               !how much direct light penetrates a singleunit of lai?
               scelg%Kb_leaf = min(kb_max,rad_params%clumping_index(ft) * gdir / cosz)

               ! To avoid singularities, we need to make sure that Kb =/ a
               ! for any of the bands...
               ! If they are too similar, it will create a very large
               ! term in the linear solution and generate solution errors
               ! Lets identify the Kb_leaf that gives a singularity.
               ! We don't need to include the min() function
               ! a will never be that large.
               !
               ! kb = a = (lai*kb_leaf + sai*kb_stem)/(lai+sai)
               ! (a*(lai+sai) - sai*kb_stem)/lai = Kb_sing
               ! or.. adjust stem Kb?
               ! (a*(lai+sai) - lai*kb_leaf)/sai = kb_stem_sing

               if(scelg%lai>nearzero) then
                  Kb_eff = scelg%Kb_leaf
               else
                  Kb_eff = Kb_stem
               end if
                                  
               ! Assume there is a singularity so that we test for it
               is_sing = .true.
               iter_sing = 0

               ! Compute the singularity for all bands
               do ib = 1,this%n_bands
                  Kb_sing(ib) = this%band(ib)%scelb(ican,icol)%a
                  if (scelg%lai>nearzero) then
                     Kb_sing(ib) = (Kb_sing(ib) * (scelg%lai+scelg%sai) - scelg%sai*Kb_stem)/scelg%lai
                  end if
               end do

               do_test_sing: do while(is_sing)
                  ! Now that we have commited to testing it, assume the solution works
                  is_sing = .false.
                  iter_sing = iter_sing
                  if(iter_sing==10)then
                     write(log_unit,*) 'error trying to remove singularity',iter_sing,scelg%Kb_leaf,Kb_stem,Kb_sing(:)
                     call endrun(msg=errMsg(sourcefile, __LINE__))
                  end if
                  ! Test to see if there is a singularity and make corrections if needed
                  if (any((abs(Kb_sing(:) - Kb_eff)) < sing_tol)) then
                     Kb_eff = Kb_eff + sing_tol
                     is_sing = .true.
                  end if
               end do do_test_sing
               
               if(scelg%lai>nearzero) then
                  scelg%Kb_leaf = Kb_eff
               else
                  Kb_stem = Kb_eff
               end if
               
               
               ! RGK: My sense is that snow should be adding optical depth
               !      but we don't have any precedent for that in the FATES
               !      code or old CLM. Re-view this.
               !!scelbp%Kb = this%frac_snow*k_snow + scelbp%Kb

               scelg%Kb = min(kb_max,(scelg%lai*scelg%Kb_leaf + scelg%sai*Kb_stem)/(scelg%lai+scelg%sai))
               
               ! Component terms for asu (single scatering albedo)
               tmp0 = gdir +  rad_params%phi2(ft) * cosz
               tmp1 =  rad_params%phi1(ft) * cosz
               tmp2 = 1._r8 - tmp1/tmp0 * log((tmp1+tmp0)/tmp1)

            end if

            do_ib: do ib = 1,this%n_bands

               associate( scelb => this%band(ib)%scelb(ican,icol) )

                 if(ft==air_ft)then

                    ! Simple provisions for a ghost element (air)
                    scelb%betab = beta_air

                 else

                    ! betab - upscatter parameter for direct beam radiation, from G. Bonan
                    ! Eq. 3.16 CLM50 Tech Man
                    ! asu is the single scattering albedo per om_veg (material reflectance)

                    asu = 0.5_r8 * gdir / tmp0 * tmp2

                    betab_veg = (1._r8 + rad_params%avmu(ft)*scelg%Kb) / (rad_params%avmu(ft)*scelg%Kb) * asu

                    om_veg  =  (scelg%lai*rad_params%om_leaf(ib,ft) + &
                         scelg%sai*rad_params%om_stem(ib,ft))/(scelg%lai+scelg%sai)

                    ! Eq. 3.7 CLM50 Tech Man
                    betab_om = betab_veg*om_veg*(1._r8-this%frac_snow) + &
                         om_snow(ib)*betab_snow(ib)*this%frac_snow

                    scelb%betab = betab_om / scelb%om

                    if(debug)then
                       if( .not.(scelb%betab==scelb%betab))then
                          write(log_unit,*)"Beam backscatter fraction is NaN"
                          write(log_unit,*) betab_om,scelb%om,om_veg,this%frac_snow,betab_veg,asu,rad_params%avmu(ft),scelg%Kb
                          call endrun(msg=errMsg(sourcefile, __LINE__))
                       end if
                    end if

                 end if

               end associate
            end do do_ib
          end associate
       end do do_ical
    end do do_ican

    return
  end subroutine ZenithPrep

  ! ================================================================================================

  subroutine GetNSCel(this)

    ! Simply return the total number
    ! of scattering elements from the
    ! multi-layer scattering element array

    class(twostream_type) :: this
    integer :: ican

    this%n_scel = 0
    do ican = 1,this%n_lyr
       this%n_scel = this%n_scel + this%n_col(ican)
    end do
    return
  end subroutine GetNSCel

  ! ===============================================================

  subroutine Solve(this, ib, &
       upper_boundary_type, & 
       Rbeam_atm, & 
       Rdiff_atm, &
       taulamb,   &
       omega,     &
       ipiv,      &
       albedo_beam, & 
       albedo_diff, &
       consv_err,   &
       frac_abs_can_beam, &
       frac_abs_can_diff, &
       frac_beam_grnd_beam, &
       frac_diff_grnd_beam, &
       frac_diff_grnd_diff)

    ! Find the scattering coefficients for two-stream radiation in the canopy.

    ! Note that these scattering coefficients are separated for scattering
    ! generated by a beam radiation boundary condition, and a diffuse radiation
    ! boundary conditions. Thus, we need not provide the magnitude of the forcing
    ! for this step. If the user provides values of 1 for the Rbeam_atm and Rdiff_atm
    ! boundary condition. It is assumed this is a normalized solution.  If values
    ! other than 1 are passed, we assume that it is not a normalized solution,
    ! and we update the data structure values this%band(ib)%Rbeam_atm and
    ! this%band(ib)%Rdiff_atm.  In a normalized solution, we will leave this
    ! as unset.
    ! In ELM and CLM, the land-model requests an albedo and other
    ! normalized output from from this algorithm for the NEXT STEP. This is
    ! due to the atmospheric model needing an albedo to calculate the downwelling
    ! radiation on the next step. THus, the asynchronous nature of things.  That is
    ! why we allow a normalized solution here. When actual absorption or flux values are
    ! desired, the scattering coefficients that were determined during the normalized
    ! solution are still valid when the magnitude of the downwelling beam and diffuse
    ! radiation boundary conditions to the vegetation canopy are known.

    
    class(twostream_type) :: this
    integer   :: ib                 ! Band of interest, matches indexing of rad_params
    integer   :: upper_boundary_type ! Is this a normalized(1) or absolute(2) solution?
    
    real(r8)  :: Rbeam_atm          ! Intensity of beam radiation at top of canopy [W/m2 ground]
    real(r8)  :: Rdiff_atm          ! Intensity of diffuse radiation at top of canopy [W/m2 ground]
                                    ! 
    real(r8)  :: taulamb(:)         ! both the coefficient vector and constant side of the linear equation
    real(r8)  :: omega(:,:)         ! the square matrix to be inverted
    integer   :: ipiv(:)            ! pivot indices for LAPACK (not optional output, we don't use)
    
    real(r8) :: albedo_beam    ! Mean albedo at canopy top generated from beam radiation [-]
    real(r8) :: albedo_diff    ! Mean albedo at canopy top generated from downwelling diffuse [-]

    real(r8) :: temp_err  ! Used to build the other error terms, a temp
    real(r8) :: consv_err ! radiation canopy balance conservation
                          ! error, fraction of incident
     
    real(r8) :: frac_abs_can_beam ! Fraction of incident beam radiation absorbed by the vegetation [-]
    real(r8) :: frac_abs_can_diff ! Fraction of incident diffuse radiation absorbed by the vegetation [-]
    real(r8) :: frac_beam_grnd_beam  ! fraction of beam radiation at ground resulting from of beam at canopy top [-]
    real(r8) :: frac_diff_grnd_beam  ! fraction of down diffuse radiation at ground resulting from beam at canopy top
    real(r8) :: frac_diff_grnd_diff  ! fraction of down diffuse radiation at ground resulting from down diffuse at canopy top [-]

    ! These arrays are only used if we run in debug mode, and are
    ! looking to report the error on the linear solution e = TAU - OMEGA*LAMBDA
    real(r8),allocatable :: tau_temp(:)
    real(r8),allocatable :: omega_temp(:,:)
    
    ! Two stream solution arrays
    ! Each of these are given generic names, because
    ! they are assemblages of many terms. But generally
    ! they fit the linear algebra formulation:
    !
    ! TAU(:) = OMEGA(:,:) * LAMBDA(:)
    !
    ! Where, we invert to solve for the coefficients LAMBDA

    integer :: isol  ! Solution index loop (beam, beam+diff)
    integer :: ican  ! Loop index for canopy layers
    integer :: ibot  ! layer index for top side of layer divide
    integer :: itop  ! layer index for bottom side of layer divide
    integer :: icol  ! Loop index for canopy columns
    integer :: jcol  ! Another loop index for canopy columns
    integer :: ilem  ! Index for scattering elements
    integer :: k1,k2 ! Indices for the lambda terms in the OMEGA and LAMBDA array
    integer :: qp    ! Equation position index
    integer :: n_eq  ! Total number of equations

    integer :: ilem_off ! Offset, or total number of elements above layer of interest
    real(r8) :: b1,b2,nu_sqrd    ! intermediate terms, see documentation
    real(r8) :: Rbeam_top           ! Mean beam radiation at top of layer      [W/m2]
    real(r8) :: Rbeam_bot           ! Mean beam radiation at bottom of layer   [W/m2]
    real(r8) :: vai                 ! Vegetation area index [m2 vegetation / m2 ground]
    real(r8) :: rb_abs              ! beam absorbed over an element    [W/m2 ground]
    real(r8) :: rd_abs              ! diffuse absorbed over an element [W/m2 ground]
    real(r8) :: rd_abs_leaf         ! diffuse absorbed over leaves (dummy)
    real(r8) :: rb_abs_leaf         ! beam absorbed by leaves (dummy)
    real(r8) :: r_abs_stem          ! total absorbed by stems (dummy)
    real(r8) :: r_abs_snow          ! total absorbed by snow (dummy)
    real(r8) :: leaf_sun_frac       ! sunlit fraction of leaves (dummy)
   

    real(r8) :: beam_err,diff_err   ! error partitioned by beam and diffuse
    type(scelg_type),pointer :: scelgp   ! Pointer to the scelg data structure
    type(scelb_type),pointer :: scelbp   ! Pointer to the scelb data structure

    ! Parameters for solving via LAPACK  DGESV() and DGESVXX()
    integer :: info                                 ! Procedure diagnostic ouput

    ! Testing switch
    ! If true, then allow elements
    ! of different layers, but same row, to have priority
    ! flux into the other element, instead of a mix
    logical, parameter :: continuity_on = .true.    

    logical, parameter :: albedo_corr = .true.

    ! ------------------------------------------------------------------------------------
    ! Example system of equations for 2 parallel columns in each of two canopy
    ! layers.  Each line is one of the balanc equations. And the x's are
    ! the unknown coefficients used in those equations.  2 coefficients
    ! map to each element, and read left to right.
    ! EL1 is the element in top layer left column.
    ! EL2 is the element in the top layer, right column
    ! EL3 is the element in the bottom layer, left column
    ! EL4 is the element in the bottom layer, right column
    !
    !                                              EL1 EL2 EL3 EL4
    ! EQ: Idn balance with upper BC can1, col 1:   x x 
    ! EQ: Idn balance with upper BC can1, col 2:       x x
    ! EQ: Idn balance between upper & lower        x x x x x x
    ! EQ: Idn balance between upper & lower        x x x x     x x
    ! EQ: Iup balance between lower & upper        x x x x x x x x
    ! EQ: Iup balance between lower & upper        x x x x x x x x
    ! EQ: Iup/Idn balance with ground, 1st col:            x x
    ! EQ: Iup/Idn Balance with ground, 2nd lower col:          x x  
    !
    !     Note: The Iup balance between layers requires ALL
    !     terms, because light comes out of both
    !     upper canopy elements and reflects off soil
    !     AND, light upwells from both lower elements.
    !
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! Beam Scattering
    ! First do the direct beam stuff.  It is a trivial solution
    ! and is required as a boundary condition to the diffuse solver
    ! All parallel layers recieve downwelling form the
    ! atmosphere.
    ! Rbeam0 is the upper boundary condition provided by data or another
    ! model.
    ! Rbeam() is the incident beam radiation at the top of each layer
    ! upper canopy.
    ! --------------------------------------------------------------------------

    if((Rbeam_atm+Rdiff_atm)<nearzero)then
       write(log_unit,*)"No radiation"
       write(log_unit,*)"Two stream should not had been called"
       write(log_unit,*)Rbeam_atm,Rdiff_atm
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    
    Rbeam_top = 1.0_r8
    do ican = 1,this%n_lyr
       Rbeam_bot = 0._r8
       do icol = 1,this%n_col(ican)
          scelgp => this%scelg(ican,icol)
          scelbp => this%band(ib)%scelb(ican,icol)
          scelbp%Rbeam0 = Rbeam_top
          Rbeam_bot = Rbeam_bot + &
               Rbeam_top*scelgp%area*exp(-scelgp%Kb*(scelgp%lai+scelgp%sai))
       end do
       Rbeam_top = Rbeam_bot
    end do

    ! Calculate element-level intermediate terms to the solve
    ! These are dependent on leaf level scattering and beam scattering
    ! These values will be used to populate the matrix solve
    ! =====================================================================

    do ican = 1,this%n_lyr
       do icol = 1,this%n_col(ican)

          scelgp => this%scelg(ican,icol)
          scelbp => this%band(ib)%scelb(ican,icol)

          b2 = -(scelgp%Kd*(1._r8-scelbp%om)*(1._r8-2._r8*scelbp%betab)+scelgp%Kb) * &
               scelbp%om*scelgp%Kb*scelbp%Rbeam0
          
          b1 = -(scelgp%Kd*(1._r8-scelbp%om+2._r8*scelbp%om*scelbp%betad) + &
               (1._r8-2._r8*scelbp%betab)*scelgp%Kb) * &
               scelbp%om*scelgp%Kb*scelbp%Rbeam0
          
          nu_sqrd = (1._r8-scelbp%om)/(1._r8-scelbp%om+2._r8*scelbp%om*scelbp%betad)
          
          if(nu_sqrd<0._r8)then
             write(log_unit,*)'nu_sqrd is less than zero'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          
          ! B_1 term from documentation:
          scelbp%B1  = 0.5_r8*(1._r8+sqrt(nu_sqrd))
          
          ! B_2 term from documentation
          scelbp%B2 = 0.5_r8*(1._r8-sqrt(nu_sqrd))
          
          ! A_2 term from documentation
          scelbp%Ad    = -0.5_r8*(b1+b2)/(scelbp%a*scelbp%a-scelgp%Kb*scelgp%Kb)   ! aka half b2 minus b1

          ! A_1 term from documentation
          scelbp%Au    = -0.5_r8*(b1-b2)/(scelbp%a*scelbp%a-scelgp%Kb*scelgp%Kb)   ! aka half b1 plus b2
          
       end do
    end do

    ! =====================================================================
    ! Set up the linear systems solver
    !
    ! [TAU] = [OMEGA]*[LAMBDA]
    ! OMEGA(n_equations,n_coefficients)
    ! TAU(n_equations)
    ! LAMBDA (n_coefficients) (the solution)
    !
    ! Indexing Variables
    ! ilem : element position
    ! k1 : coefficient 1 position
    ! k2 : coefficient 2 position
    ! qp : equation position, this continues to increment
    ! =====================================================================

    n_eq = 2*this%n_scel
    
    ! TO-DO: MAKE THIS SCRATCH SPACE AT THE SITE LEVEL?
    !!allocate(OMEGA(2*this%n_scel,2*this%n_scel),stat=alloc_err)
    !!allocate(TAU(2*this%n_scel),stat=alloc_err)
    !!allocate(LAMBDA(2*this%n_scel),stat=alloc_err)

    ! We come up with two solutions:
    ! First: we run with now diffuse downwelling
    ! radiation, this allows us to calculate
    ! the canopy top albedo for beam radiation only
    ! which is useful for coupling with the atmosphere
    ! Second: we run with bot simultaneously, and
    ! use that solution to understand everything
    ! else, including the absorbed radiation

    do_isol: do isol = 1,2

       
       ! This is temporary (these need to be set
       ! because this routine makes a call to get normalized
       ! absorbtions to get total noramalized canopy absorbtion)
       ! We will set it back to unknown following that call
       
       if(isol==1)then
          this%band(ib)%Rbeam_atm = 1.0_r8
          this%band(ib)%Rdiff_atm = 0.0_r8
       else
          this%band(ib)%Rbeam_atm = 0.0_r8
          this%band(ib)%Rdiff_atm = 1.0_r8
       end if

       omega(1:n_eq,1:n_eq) = 0._r8
       taulamb(1:n_eq)      = 0._r8
       
       ! --------------------------------------------------------------------
       ! I. Flux equations with the atmospheric boundary
       ! These balance with all elements in the upper
       ! canopy, only.  The upper canopy is layer 1.
       ! --------------------------------------------------------------------

       qp = 0    ! e"Q"uation "P"osition
       do icol = 1,this%n_col(1)
          scelgp => this%scelg(1,icol)
          scelbp => this%band(ib)%scelb(1,icol)
          ilem = icol
          qp   = qp   + 1
          k1 = 2*(ilem-1)+1
          k2 = k1+1
          taulamb(qp)      =  this%band(ib)%Rdiff_atm - this%band(ib)%Rbeam_atm*scelbp%Ad
          omega(qp,k1) =  scelbp%B2
          omega(qp,k2) =  scelbp%B1
       end do


       if_understory: if(this%n_lyr>1) then


          ! -------------------------------------------------------------------
          ! II. Flux equations between canopy layers, DOWNWELLING
          ! We only perform flux balancing between layers
          ! if we have any understory, this is true if ican>1
          ! -------------------------------------------------------------------
          ! Refer to Equation X in technical document
          ! ------------------------------------------------------------

          ! This is the index offset for the layer above the
          ! current layer of interest. We start by evaluating
          ! Layer 2, so the offset refers to layer 1, and a
          ! value of 0

          ilem_off = 0
          do_dn_ican: do ican = 2,this%n_lyr

             itop = ican-1  ! Top layer of the balance
             ibot = ican    ! Bottom layer of the balance

             ! Downwelling, includes all members from top for
             ! each independant member below

             do jcol = 1,this%n_col(ibot)

                qp = qp + 1
                ilem = ilem_off + this%n_col(itop) + jcol
                k1 = 2*(ilem-1)+1
                k2 = k1 + 1

                ! Include the self terms for the current element
                ! This term is at v=0

                taulamb(qp) = this%band(ib)%Rbeam_atm*this%band(ib)%scelb(ibot,jcol)%Ad
                omega(qp,k1) = omega(qp,k1) - this%band(ib)%scelb(ibot,jcol)%B2
                omega(qp,k2) = omega(qp,k2) - this%band(ib)%scelb(ibot,jcol)%B1

                ! We need to include the terms from
                ! all elements above the current element of interest
                ! (this can be moved out of jcol loop for efficiency)
                do icol = 1,this%n_col(itop)

                   ilem = ilem_off + icol
                   k1 = 2*(ilem-1)+1
                   k2 = k1 + 1

                   scelgp => this%scelg(itop,icol)
                   scelbp => this%band(ib)%scelb(itop,icol)

                   vai = scelgp%lai + scelgp%sai

                   taulamb(qp) = taulamb(qp) - scelgp%area * this%band(ib)%Rbeam_atm*scelbp%Ad *exp(-scelgp%Kb*vai)
                   omega(qp,k1) = omega(qp,k1) + scelgp%area * scelbp%B2*exp(scelbp%a*vai)
                   omega(qp,k2) = omega(qp,k2) + scelgp%area * scelbp%B1*exp(-scelbp%a*vai)

                end do

             end do

             ilem_off = ilem_off + this%n_col(itop)

          end do do_dn_ican


          ! -------------------------------------------------------------------
          ! III. Flux equations between canopy layers, UPWELLING
          ! -------------------------------------------------------------------
          ! Refer to equation X in the technical documentation.
          ! Note the upwelling balance is performed on the upper layer,
          ! one equation for each element in the upper layer.
          ! Note that since we use "ghost elements" or air elements
          ! we don't have to factor in reflections from exposed ground.
          ! These effects will be mediated through the ghost elements
          ! -------------------------------------------------------------------

          ilem_off = 0

          do_up_ican: do ican = 2,this%n_lyr

             itop = ican-1
             ibot = ican

             do icol = 1,this%n_col(itop)

                qp = qp + 1

                ! Self terms (ie the upwelling evaluated at the bottom edge of each top element)
                ilem = ilem_off + icol
                k1   = 2*(ilem-1)+1
                k2   = k1 + 1
                scelgp => this%scelg(itop,icol)
                scelbp => this%band(ib)%scelb(itop,icol)

                vai = scelgp%lai + scelgp%sai
                taulamb(qp) = this%band(ib)%Rbeam_atm*scelbp%Au*exp(-scelgp%Kb*vai)
                omega(qp,k1) = omega(qp,k1) - scelbp%B1*exp(scelbp%a*vai)
                omega(qp,k2) = omega(qp,k2) - scelbp%B2*exp(-scelbp%a*vai)

                ! Terms for mean diffuse exiting lower elements (move out of this loop for efficiency)
                do jcol = 1,this%n_col(ibot)
                   ilem = ilem_off + this%n_col(itop) + jcol
                   k1 = 2*(ilem-1)+1
                   k2 = k1 + 1
                   scelgp => this%scelg(ibot,jcol)
                   scelbp => this%band(ib)%scelb(ibot,jcol)

                   taulamb(qp) = taulamb(qp) - this%band(ib)%Rbeam_atm*scelgp%area*scelbp%Au
                   omega(qp,k1) = omega(qp,k1) + scelgp%area*scelbp%B1
                   omega(qp,k2) = omega(qp,k2) + scelgp%area*scelbp%B2
                end do

             end do

             ilem_off = ilem_off + this%n_col(itop)
          end do do_up_ican


       end if if_understory


       ! Flux balance equations between the understory elements, and
       ! the ground below them
       ilem_off = 0
       do ican=1,this%n_lyr-1
          ilem_off = ilem_off + this%n_col(ican)
       end do

       do jcol = 1,this%n_col(this%n_lyr)

          ilem = ilem_off + jcol
          qp = qp + 1
          k1 = 2*(ilem-1)+1
          k2 = k1 + 1

          scelgp => this%scelg(this%n_lyr,jcol)
          scelbp => this%band(ib)%scelb(this%n_lyr,jcol)

          vai = scelgp%lai + scelgp%sai

          taulamb(qp) = this%band(ib)%Rbeam_atm*(scelbp%Au*exp(-scelgp%Kb*vai)  &
               - this%band(ib)%albedo_grnd_diff*scelbp%Ad*exp(-scelgp%Kb*vai) &
               - this%band(ib)%albedo_grnd_beam*scelbp%Rbeam0*exp(-scelgp%Kb*vai))

          omega(qp,k1) = omega(qp,k1) - scelbp%B1*exp(scelbp%a*vai)
          omega(qp,k2) = omega(qp,k2) - scelbp%B2*exp(-scelbp%a*vai)

          omega(qp,k1) = omega(qp,k1) + this%band(ib)%albedo_grnd_diff*scelbp%B2*exp(scelbp%a*vai)
          omega(qp,k2) = omega(qp,k2) + this%band(ib)%albedo_grnd_diff*scelbp%B1*exp(-scelbp%a*vai)

       end do

       ! dgesv will overwrite TAU with LAMBDA
       ! ie, left side of TAU = OMEGA*LAMBDA
       ! lets dave it temporarily

       if(debug)then
          allocate(tau_temp(n_eq),omega_temp(n_eq,n_eq))
          tau_temp(1:n_eq) = taulamb(1:n_eq)
          omega_temp(1:n_eq,1:n_eq) = omega(1:n_eq,1:n_eq) 
       end if

       ! the desired precision of the iterative algorithm is:
       ! RNRM < SQRT(N)*XNRM*ANRM*EPS
       ! eps is the machine precision of the real number type
       ! ANRM is the "infinity-norm", ie the infinity norm is the abs() maximum row sum
       ! XNRM is the abs() maximum value in that column
       
       ! Find the solution
       call dgesv(n_eq, 1, omega(1:n_eq,1:n_eq), n_eq, ipiv(1:n_eq),  taulamb(1:n_eq), n_eq, info)
       
       if(info.ne.0)then
          write(log_unit,*) 'Could not find a solution via dgesv'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       ! Perform a forward check on the solution error
       do ilem = 1,n_eq
          temp_err = tau_temp(ilem) - sum(taulamb(1:n_eq)*omega_temp(ilem,1:n_eq))
          if(abs(temp_err)>rel_err_thresh)then
             write(log_unit,*) 'Poor forward solution on two-stream solver'
             write(log_unit,*) 'isol (1=beam or 2=diff): ',isol
             write(log_unit,*) 'i (equation): ',ilem
             write(log_unit,*) 'band index (1=vis,2=nir): ',ib
             write(log_unit,*) 'error (tau(i) - omega(i,:)*lambda(:)) ',temp_err
             this%band(ib)%Rbeam_atm = 1._r8
             this%band(ib)%Rdiff_atm = 1._r8
             call this%Dump(ib)
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end do
       deallocate(tau_temp,omega_temp)
              

       ! Save the solution terms

       ilem_off = 0
       if(isol==1)then  !Beam
          do ican = 1,this%n_lyr
             do icol = 1,this%n_col(ican)
                ilem = ilem_off + icol
                k1 = 2*(ilem-1)+1
                k2 = k1 + 1
                scelgp => this%scelg(ican,icol)
                scelbp => this%band(ib)%scelb(ican,icol)
                scelbp%lambda1_beam = taulamb(k1)
                scelbp%lambda2_beam = taulamb(k2)
                ! The lambda diff terms will be
                ! multiplied by zero before we use them
                ! but, we dont want things like nan's
                ! or weird math, so we set them to zero too
                scelbp%lambda1_diff = 0._r8
                scelbp%lambda2_diff = 0._r8
             end do
             ilem_off = ilem_off + this%n_col(ican)
          end do
       else
          do ican = 1,this%n_lyr
             do icol = 1,this%n_col(ican)
                ilem = ilem_off + icol
                k1 = 2*(ilem-1)+1
                k2 = k1 + 1
                scelgp => this%scelg(ican,icol)
                scelbp => this%band(ib)%scelb(ican,icol)
                scelbp%lambda1_diff = taulamb(k1)
                scelbp%lambda2_diff = taulamb(k2)
             end do
             ilem_off = ilem_off + this%n_col(ican)
          end do
       end if
          
       ! Process the total canopy absorbed radiation in the
       ! two types of radiation, as well as the downwelling
       ! flux at the ground interface
       ! --------------------------------------------------------------------------------
       
       if_beam: if(isol==1)then

          ican = 1
          albedo_beam = 0._r8
          do icol = 1,this%n_col(ican)
             scelgp => this%scelg(ican,icol)
             scelbp => this%band(ib)%scelb(ican,icol)
             albedo_beam = albedo_beam + &
                  scelgp%area * this%GetRdUp(ican,icol,ib,0._r8)
          end do

          frac_diff_grnd_beam = 0._r8
          frac_beam_grnd_beam = 0._r8
          ican = this%n_lyr
          do icol = 1,this%n_col(ican)
             scelgp => this%scelg(ican,icol)
             scelbp => this%band(ib)%scelb(ican,icol)
             frac_diff_grnd_beam = frac_diff_grnd_beam + &
                  scelgp%area*this%GetRdDn(ican,icol,ib,scelgp%lai+scelgp%sai)
             frac_beam_grnd_beam = frac_beam_grnd_beam + &
                  scelgp%area*scelbp%Rbeam0*exp(-scelgp%Kb*(scelgp%lai+scelgp%sai))
          end do

          
          frac_abs_can_beam = 0._r8
          do ican = 1,this%n_lyr
             do icol = 1,this%n_col(ican)
                scelgp => this%scelg(ican,icol)
                scelbp => this%band(ib)%scelb(ican,icol)
                call this%GetAbsRad(ican,icol,ib, 0._r8,scelgp%lai+scelgp%sai, &
                     rb_abs,rd_abs,rd_abs_leaf,rb_abs_leaf,r_abs_stem,r_abs_snow,leaf_sun_frac)
                frac_abs_can_beam = frac_abs_can_beam + scelgp%area*(rb_abs+rd_abs)
             end do
          end do
          
       else  ! Diffuse
             
          albedo_diff = 0._r8
          do icol = 1,this%n_col(1)
             scelgp => this%scelg(1,icol)
             scelbp => this%band(ib)%scelb(1,icol)
             albedo_diff = albedo_diff + &
                  scelgp%area * this%GetRdUp(1,icol,ib,0._r8)
          end do

          frac_abs_can_diff = 0._r8
          
          do ican = 1,this%n_lyr
             do icol = 1,this%n_col(ican)
                scelgp => this%scelg(ican,icol)
                scelbp => this%band(ib)%scelb(ican,icol)
                call this%GetAbsRad(ican,icol,ib,0._r8,scelgp%lai+scelgp%sai, &
                     rb_abs,rd_abs,rd_abs_leaf,rb_abs_leaf,r_abs_stem,r_abs_snow,leaf_sun_frac)
                frac_abs_can_diff = frac_abs_can_diff + scelgp%area*rd_abs
             end do
          end do
          
          frac_diff_grnd_diff = 0._r8
          ican = this%n_lyr
          do icol = 1,this%n_col(ican)
             scelgp => this%scelg(ican,icol)
             scelbp => this%band(ib)%scelb(ican,icol)
             frac_diff_grnd_diff = frac_diff_grnd_diff + &
                  scelgp%area*this%GetRdDn(ican,icol,ib,scelgp%lai+scelgp%sai)
          end do
          
       end if if_beam
       
    end do do_isol

    
    ! Check the error balance
    ! ---------------------------------------------------------------------------------------------

    ! Source = upwelling + canopy absorbed + ground absorbed
    
    consv_err = ((Rbeam_atm + Rdiff_atm) - &
               (albedo_diff  + albedo_beam ) - & 
               (frac_abs_can_diff   + frac_abs_can_beam) - & 
               ((frac_diff_grnd_diff+frac_diff_grnd_beam)*(1._r8-this%band(ib)%albedo_grnd_diff)) - &
               (frac_beam_grnd_beam*(1._r8-this%band(ib)%albedo_grnd_beam)) ) / (Rbeam_atm + Rdiff_atm)

    ! This is an error magnitude, not a bias
    consv_err = abs(consv_err)
    
    beam_err = Rbeam_atm - (albedo_beam + frac_abs_can_beam + &
         frac_diff_grnd_beam*(1._r8-this%band(ib)%albedo_grnd_diff) + &
         frac_beam_grnd_beam*(1._r8-this%band(ib)%albedo_grnd_beam))

    diff_err = Rdiff_atm - (albedo_diff + frac_abs_can_diff + &
         frac_diff_grnd_diff*(1._r8-this%band(ib)%albedo_grnd_diff))
    
    if( consv_err > rel_err_thresh ) then
       write(log_unit,*)"Total canopy flux balance not closing in TwoStrteamMLPEMod:Solve"
       write(log_unit,*)"Relative Error, delta/(Rbeam_atm+Rdiff_atm) :",consv_err
       write(log_unit,*)"Max Error: ",rel_err_thresh
       write(log_unit,*)"ib: ",ib
       write(log_unit,*)"scattering coeff: ",(2*rad_params%om_leaf(ib,1)+0.5*rad_params%om_stem(ib,1))/2.5
       write(log_unit,*)"Breakdown:",this%n_lyr
       do ican = 1,this%n_lyr
          do icol = 1,this%n_col(ican)
             scelgp => this%scelg(ican,icol)
             scelbp => this%band(ib)%scelb(ican,icol)
             write(log_unit,*)"    ",ican,icol
             write(log_unit,*)"    ",scelgp%lai+scelgp%sai,scelgp%pft,scelgp%area
             write(log_unit,*)"    ",scelbp%om,scelgp%Kb,scelgp%Kd,scelbp%betab,scelbp%betad
             write(log_unit,*)"    ",scelbp%om*(1.0-scelbp%betad)
             write(log_unit,*)"    ",scelbp%lambda1_beam,scelbp%lambda2_beam
             write(log_unit,*)"    ",scelbp%lambda1_diff,scelbp%lambda2_diff
             write(log_unit,*)"AB TERMS: ",scelbp%Ad,scelbp%Au,scelbp%B1,scelbp%B2,scelbp%a
             write(log_unit,*)"LAMBDA TERMS: ",scelbp%lambda1_diff,scelbp%lambda2_diff,scelbp%lambda1_beam,scelbp%lambda2_beam
          end do
       end do
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! Re-cast the abledos so they are direct result of the components.
    ! CESM and E3SM have higher tolerances. We close to 1e-6 but they
    ! close to 1e-8, which is just very difficult when the canopies
    ! get complex
    if(albedo_corr)then

       albedo_beam = Rbeam_atm - (frac_abs_can_beam + &
            frac_diff_grnd_beam*(1._r8-this%band(ib)%albedo_grnd_diff) + &
            frac_beam_grnd_beam*(1._r8-this%band(ib)%albedo_grnd_beam))
       
       albedo_diff = Rdiff_atm - (frac_abs_can_diff + &
            frac_diff_grnd_diff*(1._r8-this%band(ib)%albedo_grnd_diff))

    end if
    
    ! Set the boundary conditions back to unknown for a normalized solution
    ! This prevents us from calling the absorption and flux query routines incorrectly.
    ! For non-normalized, set it to the actual input boundary conditions
    
    if(upper_boundary_type.eq.normalized_upper_boundary) then
       this%band(ib)%Rbeam_atm = unset_r8
       this%band(ib)%Rdiff_atm = unset_r8
    else
       this%band(ib)%Rbeam_atm = Rbeam_atm
       this%band(ib)%Rdiff_atm = Rdiff_atm
    end if

    
    return
  end subroutine Solve


end Module TwoStreamMLPEMod
