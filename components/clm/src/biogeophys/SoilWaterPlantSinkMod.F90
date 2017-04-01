module SoilWaterPlantSinkMod

   use clm_varctl       , only : use_hydrstress
   use clm_varctl       , only : use_ed
   use decompMod        , only : bounds_type
   use shr_kind_mod          , only : r8 => shr_kind_r8
   use shr_log_mod           , only : errMsg => shr_log_errMsg
   use abortutils            , only : endrun
   use clm_varctl            , only : iulog
   use landunit_varcon       , only : istsoil,istcrop
   use column_varcon         , only : icol_road_perv
   implicit none

   character(len=*), parameter, private :: sourcefile = &
         __FILE__

contains
   
   subroutine Compute_EffecRootFrac_And_VertTranSink(bounds, num_hydrologyc, &
         filter_hydrologyc, soilstate_inst, canopystate_inst, waterflux_inst, clm_fates)
      
      ! ---------------------------------------------------------------------------------
      ! This is a wrapper for calculating the effective root fraction and soil
      ! water sink due to plant transpiration. 
      ! Calculate Soil Water Sink to Roots over different types
      ! of columns and for different process modules
      ! The super-set of all columns that should have a root water sink
      ! is filter_hydrologyc
      ! There are three groups of columns:
      ! 1) impervious roads, 2) non-natural vegetation and 3) natural vegetation
      ! There are several methods available.
      ! 1) the default version, 2) hydstress version and 3) fates boundary conditions
      !
      ! There are only two quantities that are the result of this routine, and its
      ! children:
      !   waterflux_inst%qflx_rootsoi_col(c,j)
      !   soilstate_inst%rootr_col(c,j)
      !
      !
      ! ---------------------------------------------------------------------------------

      use SoilStateType       , only : soilstate_type
      use WaterFluxType       , only : waterflux_type
      use CanopyStateType     , only : canopystate_type
      use CLMFatesInterfaceMod, only : hlm_fates_interface_type
      use ColumnType          , only : col 
      use LandunitType        , only : lun

      ! Arguments
      type(bounds_type)       , intent(in)    :: bounds               ! bounds
      integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
      integer                 , intent(in)    :: filter_hydrologyc(num_hydrologyc) ! column filter for soil points
      type(hlm_fates_interface_type), intent(inout) :: clm_fates
      type(soilstate_type)    , intent(inout) :: soilstate_inst
      type(waterflux_type)    , intent(inout) :: waterflux_inst
      type(canopystate_type)  , intent(in)    :: canopystate_inst
    

      ! Local Variables
      integer  :: filterc(bounds%endc-bounds%begc+1)           !column filter
      integer  :: num_filterc
      integer  :: num_filterc_tot
      integer  :: fc
      integer  :: c
      integer  :: l

      num_filterc_tot = 0
   

      ! 1) impervious roads
      num_filterc = 0
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if (col%itype(c) == icol_road_perv) then
            num_filterc = num_filterc + 1
            filterc(num_filterc) = c
         end if
      end do
      num_filterc_tot = num_filterc_tot+num_filterc
      if(use_hydrstress) then
         call Compute_EffecRootFrac_And_VertTranSink_HydStress_Roads(bounds, &
               num_filterc,filterc, soilstate_inst, waterflux_inst)
      else
         call Compute_EffecRootFrac_And_VertTranSink_Default(bounds, &
               num_filterc,filterc, soilstate_inst, waterflux_inst)
      end if



         
      ! 2) not-road or natural vegetation, everything else
      num_filterc = 0
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         l = col%landunit(c)
         if ( (col%itype(c) /= icol_road_perv) .and. (lun%itype(l) /= istsoil) ) then
            num_filterc = num_filterc + 1
            filterc(num_filterc) = c
         end if
      end do
      num_filterc_tot = num_filterc_tot+num_filterc
      if(use_hydrstress) then
         call Compute_EffecRootFrac_And_VertTranSink_HydStress(bounds, &
               num_filterc, filterc, waterflux_inst, soilstate_inst, canopystate_inst)
      else
         call Compute_EffecRootFrac_And_VertTranSink_Default(bounds, &
               num_filterc,filterc, soilstate_inst, waterflux_inst)
      end if
      


      
      ! 3) Natural vegetation
      num_filterc = 0
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         l = col%landunit(c)
         if ( (lun%itype(l) == istsoil) ) then
            num_filterc = num_filterc + 1
            filterc(num_filterc) = c
         end if
      end do

      num_filterc_tot = num_filterc_tot+num_filterc

      if( .not. use_ed ) then

         if(use_hydrstress) then
            call Compute_EffecRootFrac_And_VertTranSink_HydStress(bounds, &
               num_filterc, filterc, waterflux_inst, soilstate_inst, canopystate_inst)
         else
            call Compute_EffecRootFrac_And_VertTranSink_Default(bounds, &
                  num_filterc,filterc, soilstate_inst, waterflux_inst)
         end if

      else

         call Compute_EffecRootFrac_And_VertTranSink_Default(bounds, num_filterc, &
               filterc, soilstate_inst, waterflux_inst)

         call clm_fates%TransferRootSoilFlux(bounds, num_filterc, filterc, soilstate_inst, waterflux_inst)
                  
      end if

      if (num_hydrologyc /= num_filterc_tot) then
          write(iulog,*) 'The total number of columns flagged to root water uptake'
          write(iulog,*) 'did not match the total number calculated'
          write(iulog,*) 'This is likely a problem with the interpretation of column/lu filters.'
          call endrun(msg=errMsg(sourcefile, __LINE__))
      end if


      return
   end subroutine Compute_EffecRootFrac_And_VertTranSink

   ! ====================================================================================

   subroutine Compute_EffecRootFrac_And_VertTranSink_HydStress_Roads(bounds, &
         num_filterc,filterc, soilstate_inst, waterflux_inst)
      
      use SoilStateType    , only : soilstate_type
      use WaterFluxType    , only : waterflux_type
      use clm_varpar       , only : nlevsoi
      use clm_varpar       , only : max_patch_per_col
      use PatchType        , only : patch
      use ColumnType       , only : col

      ! Arguments
      type(bounds_type)       , intent(in)    :: bounds      
      integer                 , intent(in)    :: num_filterc
      integer                 , intent(in)    :: filterc(:) 
      type(soilstate_type)    , intent(inout) :: soilstate_inst
      type(waterflux_type)    , intent(inout) :: waterflux_inst

      ! Locals
      integer :: j
      integer :: c
      integer :: fc
      integer :: pi
      integer :: p
      real(r8) :: temp(bounds%begc:bounds%endc)                         ! accumulator for rootr weighting


      associate(& 
            qflx_rootsoi_col    => waterflux_inst%qflx_rootsoi_col    , & ! Output: [real(r8) (:,:) ]  vegetation/soil water exchange (mm H2O/s) (+ = to atm)
            qflx_tran_veg_patch => waterflux_inst%qflx_tran_veg_patch , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm) 
            qflx_tran_veg_col   => waterflux_inst%qflx_tran_veg_col   , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
            rootr_patch         => soilstate_inst%rootr_patch         , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer  
            rootr_col           => soilstate_inst%rootr_col             & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer  
            )

        ! First step is to calculate the column-level effective rooting
        ! fraction in each soil layer. This is done outside the usual
        ! PATCH-to-column averaging routines because it is not a simple
        ! weighted average of the PATCH level rootr arrays. Instead, the
        ! weighting depends on both the per-unit-area transpiration
        ! of the PATCH and the PATCHEs area relative to all PATCHES.
        
        temp(bounds%begc : bounds%endc) = 0._r8
    

        do j = 1, nlevsoi
           do fc = 1, num_filterc
              c = filterc(fc)
              rootr_col(c,j) = 0._r8
           end do
        end do
        
        do pi = 1,max_patch_per_col
           do j = 1,nlevsoi
              do fc = 1, num_filterc
                 c = filterc(fc)
                 if (pi <= col%npatches(c)) then
                    p = col%patchi(c) + pi - 1
                    if (patch%active(p)) then
                       rootr_col(c,j) = rootr_col(c,j) + rootr_patch(p,j) * &
                             qflx_tran_veg_patch(p) * patch%wtcol(p)
                    end if
                 end if
              end do
           end do
           do fc = 1, num_filterc
              c = filterc(fc)
              if (pi <= col%npatches(c)) then
                 p = col%patchi(c) + pi - 1
                 if (patch%active(p)) then
                    temp(c) = temp(c) + qflx_tran_veg_patch(p) * patch%wtcol(p)
                 end if
              end if
           end do
        end do

   
        do j = 1, nlevsoi
           do fc = 1, num_filterc
              c = filterc(fc)
              if (temp(c) /= 0._r8) then
                 rootr_col(c,j) = rootr_col(c,j)/temp(c)
              end if
              qflx_rootsoi_col(c,j) = rootr_col(c,j)*qflx_tran_veg_col(c)
           end do
        end do
      end associate
      return
   end subroutine Compute_EffecRootFrac_And_VertTranSink_HydStress_Roads
   
   ! ==================================================================================
   
   subroutine Compute_EffecRootFrac_And_VertTranSink_HydStress( bounds, &
           num_filterc, filterc, waterflux_inst, soilstate_inst, &
           canopystate_inst)

        !
        ! Generic routine to apply transpiration as a sink condition that
        ! is vertically distributed over the soil column. Plant hydraulic
        ! stress version
        ! (Previously Named "Compute_VertTranSink_PHS", moved and renamed
        !  rgk 02-03-2017)
        !
        !USES:
        use decompMod        , only : bounds_type
        use clm_varpar       , only : nlevsoi
        use clm_varpar       , only : max_patch_per_col
        use SoilStateType    , only : soilstate_type
        use WaterFluxType    , only : waterflux_type
        use CanopyStateType  , only : canopystate_type
        use PatchType        , only : patch
        use ColumnType       , only : col
        use clm_varctl       , only : iulog
        use PhotosynthesisMod, only : plc, params_inst
        use column_varcon    , only : icol_road_perv
        use shr_infnan_mod   , only : isnan => shr_infnan_isnan

        !
        ! !ARGUMENTS:
        type(bounds_type)    , intent(in)    :: bounds          ! bounds
        integer              , intent(in)    :: num_filterc     ! number of column soil points in column filter
        integer              , intent(in)    :: filterc(:)      ! column filter for soil points
        type(waterflux_type) , intent(inout) :: waterflux_inst
        type(soilstate_type) , intent(inout) :: soilstate_inst
        type(canopystate_type) , intent(in)  :: canopystate_inst

        !
        ! !LOCAL VARIABLES:
        real(r8) , pointer :: vegwp(:,:)  ! vegetation water matric potential (mm)
        integer  :: p,c,fc,j                                              ! do loop indices
        integer  :: pi                                                    ! patch index
        real(r8) :: temp(bounds%begc:bounds%endc)                         ! accumulator for rootr weighting
        !KO Eventually krmax will be 2D with the first index being patch%itype(p) and
        !KO  the second index corresponding to 1:nlevsoi
        real(r8) :: krmax(nlevsoi)        !
        real(r8) :: fs(nlevsoi)
        real(r8) :: rai(nlevsoi)          ! 
        real(r8) :: grav2                 !
        !-----------------------------------------------------------------------   
        
        associate(&
              qflx_rootsoi_col => waterflux_inst%qflx_rootsoi_col    , & ! Output: [real(r8) (:)]  col root and soil water exchange [mm H2O/s] [+ into root]
              rootr_col        => soilstate_inst%rootr_col           , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer  
              smp              => soilstate_inst%smp_l_col           , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]
              djk              => soilstate_inst%djk_l_col           , & ! Output: [real(r8) (:,:) ] col soil transpiration sink by layer
              bsw              => soilstate_inst%bsw_col             , & ! Input: [real(r8) (:,:) ]  Clapp and Hornberger "b"
              hk_l             => soilstate_inst%hk_l_col            , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity (mm/s)
              hksat            => soilstate_inst%hksat_col           , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
              sucsat           => soilstate_inst%sucsat_col          , & ! Input: [real(r8) (:,:) ]  minimum soil suction (mm)
              tsai             => canopystate_inst%tsai_patch        , & ! Input: [real(r8) (:)   ]  patch canopy one-sided stem area index, no burying by snow
              frac_veg_nosno   => canopystate_inst%frac_veg_nosno_patch , & ! Input:  [integer  (:)  ] fraction of vegetation not covered by snow (0 OR 1) [-]  
              rootfr           => soilstate_inst%rootfr_patch        , & ! Input: [real(r8) (:,:) ]  fraction of roots in each soil layer
              ivt              => patch%itype                        , & ! Input: [integer (:)    ]  patch vegetation type
              z                => col%z                              , & ! Input: [real(r8) (:,:) ]  layer node depth (m)
              vegwp            => canopystate_inst%vegwp_patch         & ! Input/Output: [real(r8) (:,:) ]  vegetation water matric potential (mm)
              )

          krmax(:) = 2.e-9_r8

          do fc = 1, num_filterc
             c = filterc(fc)
             
             do j = 1, nlevsoi
                grav2 = z(c,j) * 1000._r8
                temp(c) = 0._r8
                do pi = 1,max_patch_per_col
                   if (pi <= col%npatches(c)) then
                      p = col%patchi(c) + pi - 1
                      if (patch%active(p).and.frac_veg_nosno(p)>0) then 
                         if (patch%wtcol(p) > 0._r8) then
                            if (.not.isnan(smp(c,j))) then
                               rai(j) = tsai(p) * rootfr(p,j)
                               fs(j)=  min(1._r8,hk_l(c,j)/(hksat(c,j)* &
                                     plc(params_inst%psi_soil_ref(ivt(p)),p,c,j,1,bsw(c,j),sucsat(c,j))))
                               temp(c) = temp(c) + rai(j) * krmax(j) * fs(j) * &
                                     (smp(c,j) - vegwp(p,4) - grav2)* patch%wtcol(p)
                            endif
                         end if
                         !new jawn, zqz should be updated if hk formula changes
                         !zqz, also: is this actually right and good wrt PFTs sharing a column??
                      end if
                   end if
                end do
                qflx_rootsoi_col(c,j) = temp(c)
                djk(c,j) = temp(c)
             end do
             
          end do
          
          ! Back out the effective root density
          do fc = 1, num_filterc
             c = filterc(fc)
             if( sum(qflx_rootsoi_col(c,:))>0.0_r8 ) then
                do j = 1, nlevsoi
                   rootr_col(c,j) = qflx_rootsoi_col(c,j)/sum( qflx_rootsoi_col(c,:))
                end do
             else
                rootr_col(c,:) = 0.0_r8
             end if
          end do
        end associate
        return
     end subroutine Compute_EffecRootFrac_And_VertTranSink_HydStress
     
     ! ==================================================================================

     subroutine Compute_EffecRootFrac_And_VertTranSink_Default(bounds, num_filterc, &
           filterc, soilstate_inst, waterflux_inst)

    !
    ! Generic routine to apply transpiration as a sink condition that
    ! is vertically distributed over the soil column. Should be
    ! applicable to any Richards solver that is not coupled to plant
    ! hydraulics.
    !
    !USES:
    use decompMod        , only : bounds_type
    use shr_kind_mod     , only : r8 => shr_kind_r8
    use clm_varpar       , only : nlevsoi, max_patch_per_col
    use SoilStateType    , only : soilstate_type
    use WaterFluxType    , only : waterflux_type
    use PatchType        , only : patch
    use ColumnType       , only : col
    use clm_varctl       , only : use_hydrstress
    use column_varcon    , only : icol_road_perv
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds                          ! bounds
    integer              , intent(in)    :: num_filterc                     ! number of column soil points in column filter
    integer              , intent(in)    :: filterc(num_filterc)            ! column filter for soil points
    type(waterflux_type) , intent(inout) :: waterflux_inst
    type(soilstate_type) , intent(inout) :: soilstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,fc,j                                              ! do loop indices
    integer  :: pi                                                    ! patch index
    real(r8) :: temp(bounds%begc:bounds%endc)                         ! accumulator for rootr weighting
    associate(& 
          qflx_rootsoi_col    => waterflux_inst%qflx_rootsoi_col    , & ! Output: [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/s) (+ = to atm)
          qflx_tran_veg_patch => waterflux_inst%qflx_tran_veg_patch , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm) 
          qflx_tran_veg_col   => waterflux_inst%qflx_tran_veg_col   , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
          rootr_patch         => soilstate_inst%rootr_patch         , & ! Input: [real(r8) (:,:) ]  effective fraction of roots in each soil layer  
          rootr_col           => soilstate_inst%rootr_col             & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer  
          )
      
      ! First step is to calculate the column-level effective rooting
      ! fraction in each soil layer. This is done outside the usual
      ! PATCH-to-column averaging routines because it is not a simple
      ! weighted average of the PATCH level rootr arrays. Instead, the
      ! weighting depends on both the per-unit-area transpiration
      ! of the PATCH and the PATCHEs area relative to all PATCHES.
      
      temp(bounds%begc : bounds%endc) = 0._r8
      
      do j = 1, nlevsoi
         do fc = 1, num_filterc
            c = filterc(fc)
            rootr_col(c,j) = 0._r8
         end do
      end do
      
      do pi = 1,max_patch_per_col
         do j = 1,nlevsoi
            do fc = 1, num_filterc
               c = filterc(fc)
               if (pi <= col%npatches(c)) then
                  p = col%patchi(c) + pi - 1
                  if (patch%active(p)) then
                     rootr_col(c,j) = rootr_col(c,j) + rootr_patch(p,j) * &
                           qflx_tran_veg_patch(p) * patch%wtcol(p)
                  end if
               end if
            end do
         end do
         do fc = 1, num_filterc
            c = filterc(fc)
            if (pi <= col%npatches(c)) then
               p = col%patchi(c) + pi - 1
               if (patch%active(p)) then
                  temp(c) = temp(c) + qflx_tran_veg_patch(p) * patch%wtcol(p)
               end if
            end if
         end do
      end do
      
      do j = 1, nlevsoi
         do fc = 1, num_filterc
            c = filterc(fc)
            if (temp(c) /= 0._r8) then
               rootr_col(c,j) = rootr_col(c,j)/temp(c)
            end if
            qflx_rootsoi_col(c,j) = rootr_col(c,j)*qflx_tran_veg_col(c)

         end do
      end do
    end associate
    return
 end subroutine Compute_EffecRootFrac_And_VertTranSink_Default

end module SoilWaterPlantSinkMod

