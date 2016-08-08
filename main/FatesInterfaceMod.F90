module FatesInterfaceMod

   ! ------------------------------------------------------------------------------------
   ! This is the FATES public API
   ! A host land model has defined and allocated a structure "fates" as
   ! defined by fates_interface_type
   !
   ! It is also likely/possible that this type is defined as a vector
   ! which is allocated by thread
   ! ------------------------------------------------------------------------------------

   ! ------------------------------------------------------------------------------------
   ! Used CLM Modules
   ! INTERF-TODO:  NO CLM MODULES SHOULD BE ACCESSIBLE BY THE FATES
   ! PUBLIC API!!!!
   ! ------------------------------------------------------------------------------------

   use EDtypesMod            , only : ed_site_type,      &
                                      numPatchesPerCol,  &
                                      ctrl_parms

   use shr_kind_mod          , only : r8 => shr_kind_r8  ! INTERF-TODO: REMOVE THIS
   use clm_varpar            , only : nlevgrnd
   

   ! ------------------------------------------------------------------------------------
   ! Notes on types
   ! For floating point arrays, it is sometimes the convention to define the arrays as
   ! POINTER instead of ALLOCATABLE.  This usually achieves the same result with subtle
   ! differences.  POINTER arrays can point to scalar values, discontinuous array slices
   ! or alias other variables, ALLOCATABLES cannnot.  According to S. Lionel 
   ! (Intel-Forum Post), ALLOCATABLES are better perfomance wise as long as they point 
   ! to contiguous memory spaces and do not alias other variables, the case here.
   ! Naming conventions:   _gl  means ground layer dimensions
   !                       _pa  means patch dimensions
   !                       _rb  means radiation band
   ! ------------------------------------------------------------------------------------
   
   type, public :: bc_in_type

      ! The actual number of FATES' ED patches
      integer :: npatches

      ! Downwelling direct beam radiation (patch,radiation-band) [W/m2]
      real(r8), allocatable :: solad_parb(:,:)  

      ! Downwelling diffuse (I-ndirect) radiation (patch,radiation-band) [W/m2]
      real(r8), allocatable :: solai_parb(:,:)

      ! Soil suction potential of layers in each site, negative, [mm]
      real(r8), allocatable :: smp_gl(:)

      ! Effective porosity = porosity - vol_ic, of layers in each site [-]
      real(r8), allocatable :: eff_porosity_gl(:)

      ! volumetric soil water at saturation (porosity)
      real(r8), allocatable :: watsat_gl(:)

      ! Temperature of ground layers [K]
      real(r8), allocatable :: tempk_gl(:)

      ! Liquid volume in ground layer
      real(r8), allocatable :: h2o_liqvol_gl(:)

      ! Site level filter for uptake response functions
      logical               :: filter_btran

      ! the index of the deepest model soil level where roots may be, due to permafrost or bedrock constraints
      integer  :: max_rooting_depth_index_col


   end type bc_in_type


   type, public :: bc_out_type

      ! Sunlit fraction of the canopy for this patch [0-1]
      real(r8),allocatable :: fsun_pa(:)

      ! Logical stating whether a ground layer can have water uptake by plants
      ! The only condition right now is that liquid water exists
      ! The name (suction) is used to indicate that soil suction should be calculated
      logical, allocatable :: active_suction_gl(:)

      ! Effective fraction of roots in each soil layer 
      real(r8), allocatable :: rootr_pagl(:,:)

      ! Integrated (vertically) transpiration wetness factor (0 to 1) 
      ! (diagnostic, should not be used by HLM)
      real(r8), allocatable :: btran_pa(:)

     ! litterfall fluxes of C from FATES patches to BGC columns
     real(r8), allocatable :: FATES_c_to_litr_lab_c_col(:)      !total labile    litter coming from ED. gC/m3/s
     real(r8), allocatable :: FATES_c_to_litr_cel_c_col(:)      !total cellulose litter coming from ED. gC/m3/s
     real(r8), allocatable :: FATES_c_to_litr_lig_c_col(:)      !total lignin    litter coming from ED. gC/m3/s


   end type bc_out_type


   type, public :: fates_interface_type
      
      ! This is the root of the ED/FATES hierarchy of instantaneous state variables
      ! ie the root of the linked lists. Each path list is currently associated with a 
      ! grid-cell, this is intended to be migrated to columns 

      integer                         :: nsites

      type(ed_site_type), allocatable :: sites(:)

      ! These are boundary conditions that the FATES models are required to be filled.  
      ! These values are filled by the driver or HLM.  Once filled, these have an 
      ! intent(in) status.  Each site has a derived type structure, which may include 
      ! a scalar for site level data, a patch vector, potentially cohort vectors (but 
      ! not yet atm) and other dimensions such as soil-depth or pft.  These vectors 
      ! are initialized by maximums, and the allocations are static in time to avoid
      ! having to allocate/de-allocate memory

      type(bc_in_type), allocatable   :: bc_in(:)

      ! These are the boundary conditions that the FATES model returns to its HLM or 
      ! driver. It has the same allocation strategy and similar vector types.
      
      type(bc_out_type), allocatable  :: bc_out(:)

   contains
      
      procedure, public :: zero_bcs

   end type fates_interface_type

  
   public :: set_fates_ctrlparms


contains

   ! ====================================================================================

   ! INTERF-TODO: THIS IS A PLACE-HOLDER ROUTINE, NOT CALLED YET...
   subroutine fates_clean(this)
      
      implicit none
      
      ! Input Arguments
      class(fates_interface_type), intent(inout) :: this
      
      ! Incrementally walk through linked list and deallocate
      
      
      
      ! Deallocate the site list
      deallocate (this%sites)
      
      return
   end subroutine fates_clean


   ! ====================================================================================
   

   subroutine allocate_bcin(bc_in)
      
      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_in_type), intent(inout) :: bc_in
      
      ! Allocate input boundaries
      
      ! Radiation
      allocate(bc_in%solad_parb(numPatchesPerCol,ctrl_parms%numSWBands))
      allocate(bc_in%solai_parb(numPatchesPerCol,ctrl_parms%numSWBands))
      
      ! Hydrology
      allocate(bc_in%smp_gl(ctrl_parms%numlevgrnd))
      allocate(bc_in%eff_porosity_gl(ctrl_parms%numlevgrnd))
      allocate(bc_in%watsat_gl(ctrl_parms%numlevgrnd))
      allocate(bc_in%tempk_gl(ctrl_parms%numlevgrnd))
      allocate(bc_in%h2o_liqvol_gl(ctrl_parms%numlevgrnd))
      
      return
   end subroutine allocate_bcin
   
   subroutine allocate_bcout(bc_out)

      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_out_type), intent(inout) :: bc_out
      
      
      ! Radiation
      allocate(bc_out%fsun_pa(numPatchesPerCol))
      
      ! Hydrology
      allocate(bc_out%active_suction_gl(ctrl_parms%numlevgrnd))
      allocate(bc_out%rootr_pagl(numPatchesPerCol,ctrl_parms%numlevgrnd))
      allocate(bc_out%btran_pa(numPatchesPerCol))

      ! biogeochemistry
      allocate(bc_out%FATES_c_to_litr_lab_c_col(ctrl_parms%numlevdecomp_full))        
      allocate(bc_out%FATES_c_to_litr_cel_c_col(ctrl_parms%numlevdecomp_full))
      allocate(bc_out%FATES_c_to_litr_lig_c_col(ctrl_parms%numlevdecomp_full))

      return
   end subroutine allocate_bcout

   ! ====================================================================================

   subroutine zero_bcs(this,s)

      implicit none
      class(fates_interface_type), intent(inout) :: this
      integer, intent(in) :: s

      ! Input boundaries

      this%bc_in(s)%solad_parb(:,:)     = 0.0_r8
      this%bc_in(s)%solai_parb(:,:)     = 0.0_r8
      this%bc_in(s)%smp_gl(:)           = 0.0_r8
      this%bc_in(s)%eff_porosity_gl(:)  = 0.0_r8
      this%bc_in(s)%watsat_gl(:)        = 0.0_r8
      this%bc_in(s)%tempk_gl(:)         = 0.0_r8
      this%bc_in(s)%h2o_liqvol_gl(:)    = 0.0_r8
      this%bc_in(s)%max_rooting_depth_index_col = 0
      
      ! Output boundaries
      this%bc_out(s)%active_suction_gl(:) = .false.
      this%bc_out(s)%fsun_pa(:)      = 0.0_r8
      this%bc_out(s)%rootr_pagl(:,:) = 0.0_r8
      this%bc_out(s)%btran_pa(:)     = 0.0_r8
      this%bc_out(s)%FATES_c_to_litr_lab_c_col(:) = 0.0_r8
      this%bc_out(s)%FATES_c_to_litr_cel_c_col(:) = 0.0_r8
      this%bc_out(s)%FATES_c_to_litr_lig_c_col(:) = 0.0_r8

      return
   end subroutine zero_bcs
   
   ! ==================================================================================== 

   subroutine set_fates_ctrlparms(tag,dimval)
      
      ! ---------------------------------------------------------------------------------
      ! INTERF-TODO:  NEED ALLOWANCES FOR REAL AND CHARACTER ARGS..
      ! Certain model control parameters and dimensions used by FATES are dictated by 
      ! the the driver or the host mode. To see which parameters should be filled here
      ! please also look at the ctrl_parms_type in FATESTYpeMod, in the section listing
      ! components dictated by the host model.
      !
      ! Some important points:
      ! 1. Calls to this function are likely from the clm_fates module in the HLM.
      ! 2. The calls should be preceeded by a flush function.
      ! 3. All values in ctrl_parm (FATESTypesMod.F90) that are classified as 
      !    'dictated by the HLM' must be listed in this subroutine
      ! 4. Should look like this:
      ! 
      ! call set_fates_ctrlparms('flush_to_unset')
      ! call set_fates_ctrlparms('num_sw_bbands',numrad)  ! or other variable
      ! ...
      ! call set_fates_ctrlparms('num_lev_ground',nlevgrnd)   ! or other variable
      ! call set_fates_ctrlparms('check_allset') 
      !
      ! RGK-2016
      ! ---------------------------------------------------------------------------------
      
      ! Arguments
      integer, optional, intent(in)  :: dimval
      character(len=*),intent(in)    :: tag
      
      ! local variables
      logical              :: all_set
      integer,  parameter  :: unset_int = -999
      real(r8), parameter  :: unset_double = -999.9
      
      select case (trim(tag))
      case('flush_to_unset')

         write(*,*) 'Flushing FATES control parameters prior to transfer from host'
         ctrl_parms%numSwBands = unset_int
         ctrl_parms%numlevgrnd = unset_int
         ctrl_parms%numlevdecomp_full = unset_int

      case('check_allset')
         
         if(ctrl_parms%numSWBands .eq. unset_int) then
            write(*,*) 'FATES dimension/parameter unset: num_sw_rad_bbands'
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if
         
         if(ctrl_parms%numlevgrnd .eq. unset_int) then
            write(*,*) 'FATES dimension/parameter unset: numlevground'
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         if(ctrl_parms%numlevdecomp_full .eq. unset_int) then
            write(*,*) 'FATES dimension/parameter unset: numlevdecomp_full'
            ! INTERF-TODO: FATES NEEDS INTERNAL end_run
            ! end_run('MESSAGE')
         end if

         write(*,*) 'Checked. All control parameters sent to FATES.'
         
      case default

         if(present(dimval))then
            select case (trim(tag))
               
            case('num_sw_bbands')
               
               ctrl_parms%numSwBands = dimval
               write(*,*) 'Transfering num_sw_bbands = ',dimval,' to FATES'
               
            case('num_lev_ground')
               
               ctrl_parms%numlevgrnd = dimval
               write(*,*) 'Transfering num_lev_ground = ',dimval,' to FATES'
               
            case('num_levdecomp_full')
               
               ctrl_parms%numlevdecomp_full = dimval
               write(*,*) 'Transfering num_levdecomp_full = ',dimval,' to FATES'
               
            case default
               write(*,*) 'tag not recognized:',trim(tag)
               ! end_run
            end select
         else
            write(*,*) 'no value was provided for the tag'
         end if
         
      end select
      
      
      return
   end subroutine set_fates_ctrlparms



end module FatesInterfaceMod
