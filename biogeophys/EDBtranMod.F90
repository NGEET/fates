module EDBtranMod
   
   !-------------------------------------------------------------------------------------
   ! Description:
   ! 
   ! ------------------------------------------------------------------------------------

  use pftconMod        , only : pftcon
  use EDTypesMod       , only : ed_patch_type, ed_cohort_type, numpft_ed
  use EDEcophysContype , only : EDecophyscon
  !
  implicit none
  private

  public :: btran_ed
  !
  type(ed_cohort_type), pointer  :: currentCohort ! current cohort
  type(ed_patch_type) , pointer  :: currentPatch ! current patch
  !------------------------------------------------------------------------------
   
contains 

  !------------------------------------------------------------------------------
  subroutine btran_ed( sites, nsites, bc_in, bc_out)


    use shr_kind_mod       , only : r8 => shr_kind_r8
    use EDtypesMod        , only  : ed_patch_type, ed_site_type
    use FatesInterfaceMod , only  : bc_in_type, &
                                    bc_out_type

    ! Arguments

    type(ed_site_type),intent(inout),target :: sites(nsites)
    integer,intent(in)                      :: nsites
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(out)           :: bc_out(nsites)

    !
    ! !LOCAL VARIABLES:
    integer :: iv !leaf layer
    integer :: s  !site
    integer :: c  !column
    integer :: j  !soil layer
    integer :: ft ! plant functional type index
    
    
    ! The root effective root fraction for different pfts and soil layers
    real(r8) :: rootr(numpft_ed,nlevgrnd)
    
    !----------------------------------------------------------------------

    ! Inputs to model from CLM. To be read in through an input file for the purposes of this program.      
    integer, parameter    :: nv = 5           ! Number of canopy layers
    real(r8) :: xksat                         ! maximum hydraulic conductivity of soil [mm/s]
    real(r8) :: s1                            ! HC intermediate
    real(r8) :: swp_mpa(nlevgrnd)             ! matrix potential - MPa
    real(r8) :: hk(nlevgrnd)                  ! hydraulic conductivity [mm h2o/s] 
    real(r8) :: rootxsecarea                  ! root X-sectional area (m2)
    real(r8) :: rootmass(nlevgrnd)            ! root mass in each layer (g) 
    real(r8) :: rootlength(nlevgrnd)          ! root length in each layer (m) 
    real(r8) :: soilr1(nlevgrnd)              ! soil-to-root resistance in each layer (MPa s m2 mmol-1)
    real(r8) :: soilr2(nlevgrnd)              ! internal root resistance in each layer (MPa s m2 mmol-1)
    real(r8) :: rs                            ! intermediate variable
    real(r8) :: soilr_z(nlevgrnd)             ! soil-to-xylem resistance in each layer (MPa s m2 mmol-1) 
    real(r8) :: lsoil(nlevgrnd)               ! hydraulic conductivity in each soil layer  

    real(r8) :: estevap(nlevgrnd)             ! potential suction from each soil layer (mmol m-2 s-1)
    real(r8) :: totestevap                    ! potential suction from each soil layer (mmol m-2 s-1)
    real(r8) :: fraction_uptake(nlevgrnd)     ! Uptake of water from each soil layer (-) 
    real(r8) :: maxevap(nlevgrnd)             ! potential suction from each soil layer (mmol m-2 s-1)
    real(r8) :: totmaxevap                    ! potential suction from each soil layer (mmol m-2 s-1)
    real(r8) :: fleaf                         ! fraction of leaves in each canopy layer

                                              ! Model parameters
    real(r8) :: head         = 0.009807_r8    ! head of pressure  (MPa/m)   
    real(r8) :: rootdens     = 0.5e6_r8       ! root density, g biomass m-3 root
    real(r8) :: pi           = shr_const_pi
    real(r8) :: vol_ice                       ! partial volume of ice lens in layer
    real(r8) :: eff_porosity                  ! effective porosity in layer
    real(r8) :: vol_liq                       ! partial volume of liquid water in layer
    real(r8) :: s_node                        ! vol_liq/eff_porosity
    real(r8) :: smp_node                      ! matrix potential

                                              ! To be read in from pft file ultimately.       
    real(r8) :: minlwp  =  -2.5_r8            ! minimum leaf water potential in MPa
    real(r8) :: rootrad  =  0.001_r8          ! root radius in metres

    ! Outputs to CLM_SPA
    real(r8) :: weighted_SWP                  ! weighted apparent soil water potential: MPa.
    real(r8) :: canopy_soil_resistance(nv)    ! Resistance experienced by each canopy layer: MPa s m2 mmol-1

    ! SPA Pointers from CLM type. 
    logical, parameter :: SPA_soil=.false.    ! Is the BTRAN model SPA or CLM? FIX(SPM,032414) ed - make this a namelist var

    real(r8) :: rresis_ft(numpft_ed,nlevgrnd) ! resistance to water uptake per pft and soil layer.
    real(r8) :: pftgs(numpft_ed)              ! pft weighted stomatal conductance s/m
    real(r8) :: temprootr                   
    !------------------------------------------------------------------------------


    do s = 1,nsites
       
       ifp = 0
       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))                 
          ifp=ifp+1
          
          do FT = 1,numpft_ed
             currentPatch%btran_ft(FT) = 0.0_r8
             do j = 1,nlevgrnd
                
                if(bc_in(s)%smp_sl(j) > -900.0_r8 ) then  ! The flag for frozen or no water is -999
                   
                   smp_node = max(smpsc(FT), bc_in(s)%smp_sl(j))
                   
                   rresis_ft = min( (bc_in(s)%eff_porosity_sl(j)/bc_in(s)%watsat(j))*               &
                         (smp_node - smpsc(FT)) / (smpso(FT) - smpsc(FT)), 1._r8)


                   currentPatch%rootr_ft(FT,j) = currentPatch%rootfr_ft(FT,j)*rresis_ft

                   ! root water uptake is not linearly proportional to root density,
                   ! to allow proper deep root funciton. Replace with equations from SPA/Newman. FIX(RF,032414)
                   ! currentPatch%rootr_ft(FT,j) = currentPatch%rootfr_ft(FT,j)**0.3*rresis_FT(FT,j)/ &
                   ! sum(currentPatch%rootfr_ft(FT,1:nlevgrnd)**0.3)
                   currentPatch%btran_ft(FT) = currentPatch%btran_ft(FT) + currentPatch%rootr_ft(FT,j)

                else
                   currentPatch%rootr_ft(FT,j) = 0._r8
                end if

            end do !j

            ! Normalize root resistances to get layer contribution to ET
            do j = 1,nlevgrnd    
               if (currentPatch%btran_ft(FT)  >  0.0_r8) then
                  currentPatch%rootr_ft(FT,j) = currentPatch%rootr_ft(FT,j)/currentPatch%btran_ft(FT)
               else
                  currentPatch%rootr_ft(FT,j) = 0._r8
               end if
            end do

         end do !PFT

         ! PFT-averaged point level root fraction for extraction purposese.
         ! This probably needs to be weighted by actual transpiration from each pft. FIX(RF,032414).
         pftgs(:) = 0._r8
         currentCohort => currentPatch%tallest
         do while(associated(currentCohort))
            pftgs(currentCohort%pft) = pftgs(currentCohort%pft) + currentCohort%gscan * currentCohort%n    
            currentCohort => currentCohort%shorter
         enddo


         ! Process the boundary output, this is necessary for calculating the soil-moisture
         ! sink term across the different layers in driver/host.  Photosynthesis will
         ! pass the host a total transpiration for the patch.  This needs rootr to be
         ! distributed over the soil layers.

         do j = 1,nlevgrnd
            bc_out(s)%rootr_pa(ifp,j) = 0._r8
            do FT = 1,numpft_ed
               if(sum(pftgs) > 0._r8)then !prevent problem with the first timestep - might fail
                  !bit-retart test as a result? FIX(RF,032414)  
                  bc_out(s)%rootr_pa(ifp,j) = bc_out(s)%rootr_pa(ifp,j) + &
                                              currentPatch%rootr_ft(FT,j) * pftgs(ft)/sum(pftgs)
               else
                  bc_out(s)%rootr_pa(ifp,j) = bc_out(s)%rootr_pa(ifp,j) + &
                                              currentPatch%rootr_ft(FT,j) * 1./numpft_ed
               end if
            enddo
         enddo
         


         cpatch => cpatch%younger
      end do
          
   end do
   



         
 do FT = 1,numpft_ed
             currentPatch%btran_ft(FT) = 0.0_r8
             do j = 1,nlevgrnd

               !Root resistance factors
               vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
               eff_porosity = watsat(c,j)-vol_ice
               vol_liq = min(eff_porosity, h2osoi_liq(c,j)/(dz(c,j)*denh2o))
               if (vol_liq  <=  0._r8 .or. t_soisno(c,j)  <=  tfrz-2._r8) then
                  currentPatch%rootr_ft(FT,j) = 0._r8
               else
                  s_node = max(vol_liq/eff_porosity,0.01_r8)
                  smp_node = max(smpsc(FT), -sucsat(c,j)*s_node**(-bsw(c,j)))
                  !FIX(RF,032414) for junipers
                  rresis_ft(FT,j) = min( (eff_porosity/watsat(c,j))* &
                       (smp_node - smpsc(FT)) / (smpso(FT) - smpsc(FT)), 1._r8)

                  currentPatch%rootr_ft(FT,j) = currentPatch%rootfr_ft(FT,j)*rresis_FT(FT,j)
                  ! root water uptake is not linearly proportional to root density,
                  ! to allow proper deep root funciton. Replace with equations from SPA/Newman. FIX(RF,032414)
                  ! currentPatch%rootr_ft(FT,j) = currentPatch%rootfr_ft(FT,j)**0.3*rresis_FT(FT,j)/ &
                  ! sum(currentPatch%rootfr_ft(FT,1:nlevgrnd)**0.3)
                  currentPatch%btran_ft(FT) = currentPatch%btran_ft(FT) + currentPatch%rootr_ft(FT,j)
               end if
            end do !j

            ! Normalize root resistances to get layer contribution to ET 
            do j = 1,nlevgrnd    
               if (currentPatch%btran_ft(FT)  >  0.0_r8) then
                  currentPatch%rootr_ft(FT,j) = currentPatch%rootr_ft(FT,j)/currentPatch%btran_ft(FT)
               else
                  currentPatch%rootr_ft(FT,j) = 0._r8
               end if
            end do

         end do !PFT

         ! PFT-averaged point level root fraction for extraction purposese.
         ! This probably needs to be weighted by actual transpiration from each pft. FIX(RF,032414).
         pftgs(:) = 0._r8
         currentCohort => currentPatch%tallest
         do while(associated(currentCohort))
            pftgs(currentCohort%pft) = pftgs(currentCohort%pft) + currentCohort%gscan * currentCohort%n    
            currentCohort => currentCohort%shorter
         enddo

         do j = 1,nlevgrnd       
            rootr(p,j) = 0._r8
            do FT = 1,numpft_ed
               if(sum(pftgs) > 0._r8)then !prevent problem with the first timestep - might fail
                  !bit-retart test as a result? FIX(RF,032414)  
                  rootr(p,j) = rootr(p,j) + currentPatch%rootr_ft(FT,j) * pftgs(ft)/sum(pftgs)
               else
                  rootr(p,j) = rootr(p,j) + currentPatch%rootr_ft(FT,j) * 1./numpft_ed
               end if
            enddo
         enddo
         


         !---------------------------------------------------------------------------------------
         ! SPA based recalculation of BTRAN and water uptake. 
         !---------------------------------------------------------------------------------------

         if (SPA_soil) then   ! normal case don't run this.
            rootr(p,:) = 0._r8
            do FT = 1,numpft_ed 

               ! Soil Physics  
               do j = 1,nlevgrnd
                  ! CLM water retention curve. Clapp and Hornberger equation.    
                  s1 = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
                  s1 = min(1.0_r8,s1)  
                  smp_node = -sucsat(c,j)*s1**(-bsw(c,j))
                  swp_mpa(j)  = smp_node *10.0_r8/1000000.0_r8  !convert from mm to Mpa 

                  ! CLM hydraulic conductivity curve. 
                  ! As opposed to the Richard's equation solution in SoilHydrology.Mod 
                  ! the conductivity here is defined in the middle of the layer in question, not at the edge... 
                  xksat   = 0.0070556_r8 * (10._r8**(-0.884_r8+0.0153_r8*sand(p)) ) 
                  hk(j)    =   xksat*s1**(2._r8*bsw(c,j)+2._r8)  !removed the ice from here to avoid 1st ts crashing        
               enddo

               ! Root resistance
               rootxsecarea=3.14159*rootrad**2
               do j = 1,nlevgrnd      
                  rootmass(j) =  EDecophyscon%soilbeta(FT) * currentPatch%rootfr_ft(FT,j)
                  rootlength(j) = rootmass(j)/(rootdens*rootxsecarea)   !m m-3 soil
                  Lsoil(j)     = hk(j)/1000/head !converts from mms-1 to ms-1 and then to m2 s-1 MPa-1    
                  if(Lsoil(j) < 1e-35_r8.or.currentPatch%rootfr_ft(ft,j) <= 0.0_r8)then   !prevent floating point error
                     soilr_z(j) = 1e35_r8
                     soilr2(j)  = 1e35_r8
                  else 
                     ! Soil-to-root water uptake from Newman (1969). 
                     rs = sqrt (1._r8  / (rootlength(j) * pi)) 
                     soilr1(j) = log(rs/rootrad) / (2.0_r8 * pi * rootlength(j) * Lsoil(j) * dz(c,j)) 
                     ! convert from MPa s m2 m-3 to MPa s m2 mmol-1     
                     soilr1(j) = soilr1(j) * 1E-6_r8 * 18_r8 * 0.001_r8  
                     ! second component of below ground resistance is related to root hydraulics
                     soilr2(j) = EDecophyscon%rootresist(FT)/(rootmass(j)*dz(c,j))
                     soilr_z(j) = soilr1(j)+soilr2(j)
                  end if
               enddo

               ! Aggregate soil layers
               totestevap=0._r8
               weighted_SWP=0._r8
               estevap=0._r8
               fraction_uptake=0._r8
               canopy_soil_resistance=0._r8  !Reset Counters
               totmaxevap = 0._r8

               ! Estimated max transpiration from LWP gradient / soil resistance
               do j = 1,nlevgrnd    
                  estevap(j) = (swp_mpa(j) - minlwp)/(soilr_z(j))
                  estevap(j) = max(0._r8,estevap(j))         ! no negative uptake 
                  maxevap(j) = (0.0_r8 - minlwp)/(soilr2(j)) 
               enddo
               totestevap = sum(estevap)
               totmaxevap = sum(maxevap)  

               ! Weighted soil water potential
               do j = 1,nlevgrnd 
                  if(totestevap > 0._r8)then
                     fraction_uptake(j) = estevap(j)/totestevap   !Fraction of total ET taken from this soil layer 
                  else
                     estevap(j) = 0._r8
                     fraction_uptake(j)=1._r8/nlevgrnd
                  end if
                  weighted_SWP = weighted_SWP + swp_mpa(j) * estevap(j)    
               enddo


               if(totestevap > 0._r8)then
                  weighted_swp = weighted_swp/totestevap 
                  ! weight SWP for the total evaporation 
               else   
                  write(iulog,*) 'empty soil', totestevap
                  ! error check
                  weighted_swp = minlwp
               end if

               ! Weighted soil-root resistance. Aggregate the conductances (1/soilR) for each soil layer
               do iv = 1,nv !leaf layers
                  fleaf = 1.0_r8/nv
                  do j = 1,nlevgrnd !root layers
                     ! Soil resistance for each canopy layer is related to leaf area
                     ! The conductance of the root system to the 
                     ! whole canopy is reduced by the fraction of leaves in this layer...
                     canopy_soil_resistance(iv) = canopy_soil_resistance(iv)+fleaf * 1.0_r8/(soilr_z(j))         
                  enddo
                  ! Turn aggregated conductance back into resistance. mmol MPa-1 s-1 m-2  to  MPa s m2 mmol-1
                  canopy_soil_resistance(iv) = 1./canopy_soil_resistance(iv)
               enddo

               currentPatch%btran_ft(FT) =  totestevap/totmaxevap         
               do j = 1,nlevgrnd       
                  if(sum(pftgs) > 0._r8)then !prevent problem with the first timestep - might fail
                     !bit-retart test as a result? FIX(RF,032414)   
                     rootr(p,j) = rootr(p,j) + fraction_uptake(j) * pftgs(ft)/sum(pftgs)
                  else
                     rootr(p,j) = rootr(p,j) + fraction_uptake(j) * 1./numpft_ed
                  end if
               enddo

            enddo !pft loop

         end if !
         !---------------------------------------------------------------------------------------
         ! end of SPA based recalculation of BTRAN and water uptake. 
         !---------------------------------------------------------------------------------------

         !weight patch level output BTRAN for the
         btran(p) = 0.0_r8
         do FT = 1,numpft_ed
            if(sum(pftgs) > 0._r8)then !prevent problem with the first timestep - might fail
               !bit-retart test as a result? FIX(RF,032414)   
               btran(p)   = btran(p) + currentPatch%btran_ft(FT)  * pftgs(ft)/sum(pftgs)
            else
               btran(p)   = btran(p) + currentPatch%btran_ft(FT) * 1./numpft_ed
            end if
         enddo

         temprootr = sum(rootr(p,:))
         if(temprootr /= 1.0_r8)then
            !write(iulog,*) 'error with rootr in canopy fluxes',sum(rootr(p,:))
            if(temprootr > 0._r8)then
               do j = 1,nlevgrnd
                  rootr(p,j) = rootr(p,j) / temprootr
               enddo
            end if
         end if

      else ! edpatch
         currentPatch%btran_ft(1:numpft_ed) = 1._r8
      end if ! edpatch

    end associate

  end subroutine btran_ed

end module EDBtranMod
