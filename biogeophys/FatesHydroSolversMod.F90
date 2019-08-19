module FatesHydroSolversMod



contains


  subroutine UpdatePlantKMax(ccohort_hydr,ccohort,csite_hydr,bc_in)

      ! ---------------------------------------------------------------------------------
      !
      ! This routine sets the maximum conductance of all compartments in the plant, from
      ! leaves, to stem, to transporting root, to the absorbing roots.
      ! These properties are dependent only on the materials (conductivity) and the
      ! geometry of the compartments.
      ! The units of all K_max values are [kg H2O s-1 MPa-1]
      !
      ! There are some different ways to represent overall conductance from node-to-node
      ! throughout the hydraulic system. Universally, all can make use of a system
      ! where we separate the hydraulic compartments of the nodes into the upper (closer
      ! to the sky) and lower (away from the sky) portions of the compartment. It is
      ! possible that due to things like xylem taper, the two portions may have different
      ! conductivity, and therefore differnet conductances.
      !
      ! Assumption 0.  This routine calculates maximum conductivity for 1 plant.
      ! Assumption 1.  The compartment volumes, heights and lengths have all been 
      !                determined, probably called just before this routine.
      !
      ! Steudle, E. Water uptake by roots: effects of water deficit. 
      ! J Exp Bot 51, 1531-1542, doi:DOI 10.1093/jexbot/51.350.1531 (2000).
      ! ---------------------------------------------------------------------------------

      ! Arguments

      type(ed_cohort_hydr_type),intent(inout),target :: ccohort_hydr
      type(ed_cohort_type),intent(in),target         :: ccohort
      type(ed_site_hydr_type),intent(in),target      :: csite_hydr
      type(bc_in_type),intent(in)                    :: bc_in



      ! Locals
      integer :: k    ! Compartment (node) index
      integer :: k_ag ! Compartment index for above-ground indexed array
      integer  :: pft          ! Plant Functional Type index
      real(r8) :: c_sap_dummy  ! Dummy variable (unused) with sapwood carbon [kg]
      real(r8) :: z_lower      ! distance between lower edge and mean petiole height [m]
      real(r8) :: z_upper      ! distance between upper edge and mean petiole height [m]
      real(r8) :: z_node       ! distance between compartment center and mph [m]
      real(r8) :: kmax_lower   ! Max conductance from compartment edge to mph [kg s-1 Mpa-1]
      real(r8) :: kmax_node    ! Max conductance from compartment edge to mph [kg s-1 Mpa-1]
      real(r8) :: kmax_upper   ! Max conductance from compartment edge to mph [kg s-1 Mpa-1]
      real(r8) :: a_sapwood    ! Mean cross section area of sapwood   [m2]
      real(r8) :: rmin_ag      ! Minimum total resistance of all above ground pathways
                               ! [kg-1 s MPa]
      real(r8) :: kmax_bg      ! Total maximum conductance of all below-ground pathways 
                               ! from the absorbing roots center nodes to the 
                               ! transporting root center node
      real(r8) :: rootfr       ! fraction of absorbing root in each soil layer
                               ! assumes propotion of absorbing root is equal
                               ! to proportion of total root
      real(r8) :: kmax_layer   ! max conductance between transporting root node
                               ! and absorbing root node in each layer [kg s-1 MPa-1]
      real(r8) :: surfarea_aroot_layer ! Surface area of absorbing roots in each
                                       ! soil layer [m2]

      real(r8),parameter :: taper_exponent = 1._r8/3._r8 ! Savage et al. (2010) xylem taper exponent [-]


      pft = ccohort%pft

      ! Get the cross-section of the plant's sapwood area [m2]
      call bsap_allom(ccohort%dbh,pft,ccohort%canopy_trim,a_sapwood,c_sap_dummy)


      ! Leaf Maximum Hydraulic Conductance
      ! The starting hypothesis is that there is no resistance inside the
      ! leaf, between the petiole and the center of storage.  To override
      ! this, make provisions by changing the kmax to a not-absurdly high 
      ! value.  It is assumed that the conductance in this default case,
      ! is regulated completely by the stem conductance from the stem's
      ! center of storage, to the petiole.

      ccohort_hydr%kmax_petiole_to_leaf = 1.e12_r8


      ! Stem Maximum Hydraulic Conductance
        
      do k=1, n_hypool_stem

         ! index for "above-ground" arrays, that contain stem and leaf
         ! in one vector
         k_ag = k+n_hypool_leaf
         
         ! Depth from the petiole to the lower, node and upper compartment edges

         z_lower = ccohort_hydr%z_node_ag(n_hypool_leaf) - ccohort_hydr%z_lower_ag(k_ag)
         z_node  = ccohort_hydr%z_node_ag(n_hypool_leaf) - ccohort_hydr%z_node_ag(k_ag)
         z_upper = ccohort_hydr%z_node_ag(n_hypool_leaf) - ccohort_hydr%z_upper_ag(k_ag)


         ! Then we calculate the maximum conductance from each the lower, node and upper 
         ! edges of the compartment to the petiole. The xylem taper factor requires
         ! that the kmax it is scaling is from the point of interest to the mean height
         ! of the petioles.  Then we can back out the conductance over just the path
         ! of the upper and lower compartments, but subtracting them as resistors in
         ! series.

         ! max conductance from upper edge to mean petiole height
         kmax_upper = EDPftvarcon_inst%hydr_kmax_node(pft,2) * &
                      xylemtaper(taper_exponent, z_upper) * &
                      a_sapwood / z_upper

         ! max conductance from node to mean petiole height
         kmax_node  = EDPftvarcon_inst%hydr_kmax_node(pft,2) * &
                      xylemtaper(taper_exponent, z_node) * &
                      a_sapwood / z_node

         ! max conductance from lower edge to mean petiole height
         kmax_lower = EDPftvarcon_inst%hydr_kmax_node(pft,2) * &
                      xylemtaper(taper_exponent, z_lower) * &
                      a_sapwood / z_lower

         ! Max conductance over the path of the upper side of the compartment
         ccohort_hydr%kmax_stem_upper(k_ag) = (1._r8/kmax_node - 1._r8/kmax_upper)**-1._r8

         ! Max conductance over the path on the loewr side of the compartment
         ccohort_hydr%kmax_stem_lower(k_ag) = (1._r8/kmax_lower - 1._r8/kmax_node)**-1._r8


       enddo

       ! Maximum conductance of the upper compartment in the transporting root
       ! that connects to the lowest stem (btw: z_lower_ag(n_hypool_ag) == 0)

       z_upper = ccohort_hydr%z_lower_ag(n_hypool_leaf)
       z_node  = ccohort_hydr%z_lower_ag(n_hypool_leaf)-ccohort_hydr%z_node_troot

       
       kmax_node = EDPftvarcon_inst%hydr_kmax_node(pft,2) * &
                   xylemtaper(taper_exponent, z_node) * &
                   a_sapwood / z_node

       kmax_upper = EDPftvarcon_inst%hydr_kmax_node(pft,2) * &
                    xylemtaper(taper_exponent, z_upper) * &
                    a_sapwood / z_upper
       
       ccohort_hydr%kmax_troot_upper = (1._r8/kmax_node - 1._r8/kmax_upper)**-1._r8


       ! The maximum conductance between the center node of the transporting root
       ! compartment, and the center node of the absorbing root compartment, is calculated
       ! as a residual.  Specifically, we look at the total resistance the plant has in
       ! the stem so far, by adding those resistances in series.
       ! Then we use a parameter to specify what fraction of the resistance
       ! should be below-ground between the transporting root node and the absorbing roots.
       ! After that total is calculated, we then convert to a conductance, and split the
       ! conductance in parallel between root layers, based on the root fraction.
       ! Note* The inverse of max conductance (KMax) is minimum resistance:
       
       
       rmin_ag = 1._r8/ccohort_hydr%kmax_petiole_to_leaf + &
                 sum(1._r8/ccohort_hydr%kmax_stem_upper(1:n_hypool_stem)) + &
                 sum(1._r8/ccohort_hydr%kmax_stem_lower(1:n_hypool_stem)) + &
                 1._r8/kmax_troot_upper

       ! Calculate the residual resistance below ground, as a resistor
       ! in series with the existing above ground
       ! Invert to find below-ground kmax
       kmax_bg = 1._r8/(rmin_ag * (1._r8/EDPftvarcon_inst%hydr_rfrac_stem(pft) - 1._r8))

       ! The max conductance of each layer is in parallel, therefore
       ! the kmax terms of each layer, should sum to kmax_bg
       do j=1,nlevsoi_hyd
          if(j == 1) then
               rootfr = zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j))
            else
               rootfr = zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j)) - &
                        zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j-1))
           end if

           kmax_layer = rootfr*kmax_bg

           ! Two transport pathways, in two compartments exist in each layer.
           ! These pathways are connected in serial.
           ! For simplicity, we simply split the resistance between the two.
           ! Mathematically, this results in simply doubling the conductance
           ! and applying to both paths.  Here are the two paths:
           ! 1) is the path between the transporting root's center node, to
           !    the boundary of the transporting root with the boundary of
           !    the absorbing root  (kmax_troot_lower)
           ! 2) is the path between the boundary of the absorbing root and
           !    transporting root, with the absorbing root's center node
           !    (kmax_aroot_upper)
           
           ccohort_hydr%kmax_troot_lower(j) = 2.0_r8 * kmax_layer
           ccohort_hydr%kmax_aroot_upper(j) = 2.0_r8 * kmax_layer

       end do

       ! Finally, we calculate maximum radial conductance from the root
       ! surface to its center node.  This transport is not a xylem transport
       ! like the calculations prior to this. This transport is through the
       ! exodermis, cortex, casparian strip and endodermis.  The actual conductance
       ! will possibly depend on the potential gradient (whether out-of the root,
       ! or in-to the root).  So we calculate the kmax's for both cases,
       ! and save them for the final conductance calculation.
       
       do j=1,nlevsoi_hyd
          
          ! Surface area of the absorbing roots for this cohort in this layer [m2]
          surfarea_aroot_layer = 2._r8 * pi_const *csite_hydr%rs1(j) * ccohort_hydr%l_aroot_layer(j)

          ! Convert from surface conductivity [kg H2O m-2 s-1 MPa-1] to [kg H2O s-1 MPa-1]
          ccohort_hydr%kmax_aroot_radial_in(j) = hydr_kmax_rsurf1 * surfarea_aroot_layer
          
          ccohort_hydr%kmax_aroot_radial_out(j) = hydr_kmax_rsurf2 * surfarea_aroot_layer

       end do


      return
    end subroutine UpdatePlantKMax



    subroutine UpdateHDiffCond1D(cohort_hydr,site_hydr,jpaths,ilayer,psi_node,flc_node,dflcdpsi_node, &
          hdiff_bound,k_bound,dhdiffdpsi0,dhdiffdpsi1,dkbounddpsi0,dkbounddpsi1)

        ! ------------------------------------------------------------------------------------------
        ! Calculate the hydraulic conductances across a list of paths.  The list is a 1D vector, and
        ! the list need not be across the whole path from stomata to the last rhizosphere shell, but
        ! it can only be 1d, which is part of a path through the plant and into 1 soil layer.
        ! ------------------------------------------------------------------------------------------

        ! !ARGUMENTS
        type(ed_cohort_hydr_type), intent(in),target :: cohort_hydr
        type(ed_site_hydr_type), intent(in),target   :: site_hydr
        integer           , intent(in)  :: jpaths(:)              ! The path indices that are to be calculated
        integer           , intent(in)  :: ilayer                 ! soil layer index of interest
        real(r8)          , intent(in)  :: flc_node(:)            ! fractional loss of conductivity at water storage nodes          [-]
        real(r8)          , intent(in)  :: dflcdpsi_node(:)       ! derivative of fractional loss of conductivity wrt psi           [MPa-1]
        real(r8)          , intent(out) :: hdiff_bound(:)         ! total water potential difference across lower boundary          [MPa]
        real(r8)          , intent(out) :: k_bound(:)             ! lower boundary hydraulic conductance of compartments            [kg s-1 MPa-1]
        real(r8)          , intent(out) :: dhdiffdpsi0(:)         ! derivative of total water potential difference wrt psi above    [-]
        real(r8)          , intent(out) :: dhdiffdpsi1(:)         ! derivative of total water potential difference wrt psi below    [-]
        real(r8)          , intent(out) :: dkbounddpsi0(:)        ! derivative of lower boundary conductance wrt psi above          [kg s-1 MPa-2]
        real(r8)          , intent(out) :: dkbounddpsi1(:)        ! derivative of lower boundary conductance wrt psi below          [kg s-1 MPa-2]


        ! !LOCAL VARIABLES:
        integer  :: inode_up                                 ! node index closest to atmosphere for the path of interest
        integer  :: inode_lo                                 ! node index further from atmosphere for path of interest
        integer  :: jpath                                    ! path index
        real(r8) :: k_bound_aroot_soil1                      ! radial conductance of absorbing roots                           [kg s-1 MPa-1]
        real(r8) :: k_bound_aroot_soil2                      ! conductance to root surface from innermost rhiz shell           [kg s-1 MPa-1]
        real(r8) :: k_lower                                  ! conductance node k to lower boundary                            [kg s-1 MPa-1]
        real(r8) :: k_upper                                  ! conductance node k+1 to upper boundary                          [kg s-1 MPa-1]
        !----------------------------------------------------------------------

        do j=1,size(jpaths,1)

            jpath=jpaths(j)

            inode_up = jpath
            inode_lo = jpath+1

            if(inode_up == 1) then

               ! Path is between the leaf node and first stem node

               istem = 1
               
               znode_up   = cohort_hydr%z_node_ag(inode_up)
               znode_lo   = cohort_hydr%z_node_ag(inode_lo)
               kmax_up    = cohort_hydr%kmax_petiole_to_leaf
               kmax_lo    = cohort_hydr%kmax_stem_upper(1)
               
            elseif(inode_up < n_hypool_ag) then

               ! Path is between stem compartments
               ! This condition is only possible if n_hypool_ag>2
               
               znode_up   = cohort_hydr%z_node_ag(inode_up)
               znode_lo   = cohort_hydr%z_node_ag(inode_lo)
               kmax_up    = cohort_hydr%kmax_stem_lower(inode_up-n_hypool_leaf)
               kmax_lo    = cohort_hydr%kmax_stem_upper(inode_lo-n_hypool_leaf)
               
            elseif(inode_up == n_hpool_ag) then
              
               ! Path is between lowest stem and transporting root
               
               znode_up = cohort_hydr%z_node_ag(n_hpool_ag)
               znode_lo = cohort_hydr%z_node_troot
               kmax_up  = cohort_hydr%kmax_stem_lower(n_hpool_ag)
               kmax_lo  = cohort_hydr%kmax_troot_upper

            elseif(inode_up == n_hpool_ag) then

               ! Path is between the transporting root 
               ! and the absorbing root nodes

               znode_up   = cohort_hydr%z_node_troot
               znode_lo   = bc_in%z_sisl(ilayer)
               kmax_up    = cohort_hydr%kmax_troot_lower(ilayer)
               kmax_lo    = cohort_hydr%kmax_aroot_upper(ilayer)

            else
               
               ! Path is between the absorbing root
               ! and the first rhizosphere shell nodes
               
               znode_up   = bc_in%z_sisl(ilayer)
               znode_lo   = bc_in%z_sisl(ilayer)

               ! Special case. Maximum conductance depends on the 
               ! potential gradient (same elevation, no geopotential
               ! required.
               
               if(cohort_hydr%psi_aroot(ilayer) < site_hydr%psisoi_liq_innershell(j)) then
                  kmax_up = cohort_hydr%kmax_aroot_radial_in(ilayer)
               else
                  kmax_up = cohort_hydr%kmax_aroot_radial_out(ilayer)
               end if

               kmax_lo = site_hydr%kmax_upper_shell(ilayer,1)

            else

               ! Path is between rhizosphere shells
               
               znode_up = bc_in%z_sisl(ilayer)
               znode_lo = bc_in%z_sisl(ilayer)

               ishell_up = inode_up - n_hypool_ag + 2 ! Remove total number of plant pools from index
               ishell_lo = ishell_up + 1

               kmax_up = site_hydr%kmax_outer_shell(ilayer,ishell_up)
               kmax_lo = site_hydr%kmax_inner_shell(ilayer,ishell_lo)
               
            end if
            
            
            ! This is the potential difference between the nodes (matric and geopotential)
            hdiff_bound(jpath) = mpa_per_pa*denh2o*grav_earth*(znode_up-znode_lo) + (psinode(inode_up)-psinode(inode_lo))


 

            if(do_kbound_upstream) then

               ! Examine direction of water flow; use the upstream node's k for the boundary k.
               ! (as suggested by Ethan Coon, LANL)
               
               if(hdiff_bound(jpath) < 0._r8) then
                  ! More potential in the lower node, use its fraction of conductivity loss
                  k_bound(jpath)       = flc_node(inode_lo) / &
                                         (1._r8/kmax_lo + 1._r8/kmax_up)
                       (1._r8/k_bound_aroot_soil1 + 1._r8/k_bound_aroot_soil2) * flc_node(k+1)  ! water moving towards atmosphere
                  dkdpsi0(jpath)  = 0._r8
                  dkdpsi1(jpath)  = kmax_bound(jpath) * dflcdpsi_node(inode_lo)

                else

                   k_path(jpath) = (

                    



                end if
            end if


            if(do_kbound_upstream) then

                ! absorbing root-1st rhizosphere shell boundary. 
                ! Comprised of two distinct conductance terms each with distinct water potentials

                if(k == (k_arootsoil)) then  

                    k_bound_aroot_soil1 =  kmax_bound_aroot_soil1 * flc_node(k)
                    k_bound_aroot_soil2 =  kmax_bound_aroot_soil2 * flc_node(k+1)

                    k_bound(k)          =  1._r8/(1._r8/k_bound_aroot_soil1 + 1._r8/k_bound_aroot_soil2)

                    dkbounddpsi0(k)     =  ((k_bound(k)/k_bound_aroot_soil1)**2._r8) * & 
                          kmax_bound_aroot_soil1*dflcdpsi_node(k)
                    dkbounddpsi1(k)     =  ((k_bound(k)/k_bound_aroot_soil2)**2._r8) * &
                          kmax_bound_aroot_soil2*dflcdpsi_node(k+1)
                else
                    ! examine direction of water flow; use the upstream node's k for the boundary k.
                    ! (as suggested by Ethan Coon, LANL)
                    if(hdiff_bound(k) < 0._r8) then
                        k_bound(k)       = kmax_bound(k) * flc_node(k+1)  ! water moving towards atmosphere
                        dkbounddpsi0(k)  = 0._r8
                        dkbounddpsi1(k)  = kmax_bound(k) * dflcdpsi_node(k+1)
                    else                                           
                        k_bound(k)       = kmax_bound(k) * flc_node(k)    ! water moving towards soil
                        dkbounddpsi0(k)  = kmax_bound(k) * dflcdpsi_node(k)
                        dkbounddpsi1(k)  = 0._r8
                    end if
                end if
            else
                k_lower                =  kmax_lower(k)   * flc_node(k)
                k_upper                =  kmax_upper(k+1) * flc_node(k+1)
                k_bound(k)             =  1._r8/(1._r8/k_lower + 1._r8/k_upper)
                dkbounddpsi0(k)        =  ((k_bound(k)/k_lower)**2._r8) * kmax_lower(k)  * dflcdpsi_node(k)
                dkbounddpsi1(k)        =  ((k_bound(k)/k_upper)**2._r8) * kmax_upper(k+1)* dflcdpsi_node(k+1)
            end if
            dhdiffdpsi0(k)  =  1.0
            dhdiffdpsi1(k)  = -1.0
        enddo
        k               = size(z_node)
        k_bound(k)      = 0._r8
        dkbounddpsi0(k) = 0._r8
        dkbounddpsi1(k) = 0._r8

      end subroutine HDiffK1D


    subroutine PlantKmax()







    end subroutine PlantKmax





end module FatesHydroSolversMod
