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

    ! ===================================================================================

    subroutine ImTaylorSolverTermsCond1D(cohort_hydr,site_hydr,ilayer,trisolve_terms)

        ! -------------------------------------------------------------------------------
        ! Calculate the hydraulic conductances across a list of paths.  The list is a 1D vector, and
        ! the list need not be across the whole path from stomata to the last rhizosphere shell, but
        ! it can only be 1d, which is part of a path through the plant and into 1 soil layer.
        ! -------------------------------------------------------------------------------

        ! !ARGUMENTS
        type(ed_cohort_hydr_type), intent(in),target :: cohort_hydr
        type(ed_site_hydr_type), intent(in),target   :: site_hydr
        integer           , intent(in)  :: ilayer              ! soil layer index of interest
        real(r8)          , intent(in)  :: psi_node(:)         ! matric potential of nodes [Mpa]
        real(r8)          , intent(in)  :: flc_node(:)         ! fractional loss of conductivity at water storage nodes          [-]
        real(r8), intent(out),optional  :: trisolve_terms(:,:) ! This contains the terms for a tri-diagonal matrix
                                                               ! which contains the constant term on the left side [col=1]
                                                               ! and for each node's solution, the terms for that node [col=3]
                                                               ! and its flanking nodes [i-1: col=2] and [i+1: col=3]

        ! Locals
        
        integer :: inode   ! node index "i"
        integer :: jpath   ! path index "j"
        integer :: ishell  ! rhizosphere shell index of the node
        integer :: i_dn    ! downstream node of current flow-path
        integer :: i_up    ! upstream node of current flow-path
        real(r8) :: kmax_up  ! maximum conductance of the upstream half of path [kg s-1 Mpa-1]
        real(r8) :: kmax_dn  ! maximum conductance of the downstream half of path [kg s-1 MPa-1]
        real(r8) :: th_node                         ! "theta" i.e. water content of node [m3 m-3]
        real(r8) :: z_node                          ! elevation of node [m]
        real(r8) :: psi_node(n_hypool_tot)          ! matric potential on node [Mpa]
        real(r8) :: ftc_node(n_hypool_tot)          ! frac total conductance on node [-]
        real(r8) :: h_node(n_hypool_tot)            ! total potential on node [Mpa]
        real(r8) :: dftc_dtheta_node(n_hypool_tot)  ! deriv FTC w.r.t. theta
        real(r8) :: dpsi_dtheta_node(n_hypool_tot)  ! deriv psi w.r.t. theta
        real(r8) :: k_eff(n_hypool_tot-1)           ! effective (used) conductance over path [kg s-1 MPa-1]
        real(r8) :: a_term(n_hypool_tot-1)          ! "A" term in the tri-diagonal implicit solve [-]
        real(r8) :: b_term(n_hypool_tot-1)          ! "B" term in the tri-diagonal implicit solve [-]

        ! -------------------------------------------------------------------------------
        ! Part 1.  Calculate node quantities:
        !          matric potential: psi_node
        !          fraction of total conductance: ftc_node
        !          total potential (matric + elevatio) h_node
        !          deriv. ftc  wrt  theta: dftc_dtheta_node
        !          deriv. psi  wrt  theta: dpsi_dtheta_node
        ! -------------------------------------------------------------------------------
        
     

        ! For leaf and stem pools

        do inode = 1,n_hypool_tot

           if (inode<=n_hypool_ag) then
              th_node = ccohort_hydr%th_ag(inode)
              z_node  = ccohort_hydr%z_node_ag(inode)
           elseif (inode==n_hypool_ag+1) then
              th_node = ccohort_hydr%th_troot(1)
              z_node  = ccohort_hydr%z_node_troot
           elseif (inode==n_hyppol_ag+2) then
              th_node = ccohort_hyd%th_aroot(ilayer)
              z_node  = bc_in(s)%z_sisl(ilayer)
           else
              ishell  = inode-(n_hypool_tot+2)
              th_node = site_hydr%h2osoi_liqvol_shell(ilayer,ishell)
              z_node  = bc_in(s)%z_sisl(ilayer)
           end if

           ! Get matric potential [Mpa]
           call psi_from_th(currentCohort%pft, porous_media(inode), ccohort_hydr%th_ag(inode), &
                psi_node(inode), site_hydr, bc_in)

           ! Get total potential [Mpa]
           h_node(inode) =  mpa_per_pa*denh2o*grav_earth*z_node(inode) + psi_node(inode)

           ! Get Fraction of Total Conductivity [-]
           call flc_from_psi(currentCohort%pft, porous_media(inode), ccohort_hydr%psi_ag(inode), &
                ftc_node(inode), site_hydr, bc_in) 

           ! deriv ftc wrt theta
           call dpsidth_from_th(currentCohort%pft, porous_media(inode), ccohort_hydr%th_ag(inode), & 
                dpsi_dtheta_node(inode), site_hydr, bc_in)
           
           call dflcdpsi_from_psi(currentCohort%pft, porous_media(inode), psi_node(inode), & 
                dftc_dpsi, site_hydr, bc_in)
           
           dftc_dtheta_node(inode) = dftc_psi * dpsi_dtheta_node(inode) 

        end do


        !--------------------------------------------------------------------------------
        ! Part 2.  Effective conductances over the path-length and Flux terms
        !          over the node-to-node paths
        !--------------------------------------------------------------------------------

        ! Path is between the leaf node and first stem node
        ! -------------------------------------------------------------------------------

        jpath    = 1
        i_dn     = 1
        i_up     = 2
        kmax_dn  = cohort_hydr%kmax_petiole_to_leaf
        kmax_up  = cohort_hydr%kmax_stem_upper(1)

        call GetImTaylorKAB(kmax_up,kmax_dn,       &
             ftc_node(i_up),ftc_node(i_dn),        & 
             h_node(i_up),h_node(i_dn),            & 
             dftc_dtheta(i_up), dftc_dtheta(i_dn), &
             dpsi_dtheta(i_up), dpsi_dtheta(i_dn), &
             k_eff(jpath),                         &
             A_term(jpath),                        & 
             B_term(jpath))
        

        ! Path is between stem nodes
        ! -------------------------------------------------------------------------------

        do jpath=2,n_hypool_ag-1

           i_dn = jpath
           i_up = jpath+1
           kmax_up    = cohort_hydr%kmax_stem_lower(inode_up-n_hypool_leaf)
           kmax_lo    = cohort_hydr%kmax_stem_upper(inode_lo-n_hypool_leaf)
           
           call GetImTaylorKAB(kmax_up,kmax_dn,       &
                ftc_node(i_up),ftc_node(i_dn),        & 
                h_node(i_up),h_node(i_dn),            & 
                dftc_dtheta(i_up), dftc_dtheta(i_dn), &
                dpsi_dtheta(i_up), dpsi_dtheta(i_dn), &
                k_eff(jpath),                         &
                A_term(jpath),                        & 
                B_term(jpath))
           
        end do

        
        ! Path is between lowest stem and transporting root

        jpath = n_hypool_ag
        i_dn  = jpath
        i_up  = jpath+1
        kmax_up  = cohort_hydr%kmax_stem_lower(n_hpool_ag)
        kmax_lo  = cohort_hydr%kmax_troot_upper

        call GetImTaylorKAB(kmax_up,kmax_dn,       &
             ftc_node(i_up),ftc_node(i_dn),        & 
             h_node(i_up),h_node(i_dn),            & 
             dftc_dtheta(i_up), dftc_dtheta(i_dn), &
             dpsi_dtheta(i_up), dpsi_dtheta(i_dn), &
             k_eff(jpath),                         &
             A_term(jpath),                        & 
             B_term(jpath))
        

        ! Path is between the absorbing root and the first
        ! rhizosphere

        jpath   = n_hypool_ag+1
        i_dn    = jpath
        i_up    = jpath+1
        kmax_up = cohort_hydr%kmax_troot_lower(ilayer)
        kmax_lo = cohort_hydr%kmax_aroot_upper(ilayer)

        call GetImTaylorKAB(kmax_up,kmax_dn,       &
             ftc_node(i_up),ftc_node(i_dn),        & 
             h_node(i_up),h_node(i_dn),            & 
             dftc_dtheta(i_up), dftc_dtheta(i_dn), &
             dpsi_dtheta(i_up), dpsi_dtheta(i_dn), &
             k_eff(jpath),                         &
             A_term(jpath),                        & 
             B_term(jpath))


        ! Path is between the absorbing root
        ! and the first rhizosphere shell nodes
        
        jpath = n_hypool_ag+2
        i_dn  = jpath
        i_up  = jpath+1

        ! Special case. Maximum conductance depends on the 
        ! potential gradient (same elevation, no geopotential
        ! required.
        if(cohort_hydr%psi_aroot(ilayer) < site_hydr%psisoi_liq_innershell(j)) then
           kmax_up = cohort_hydr%kmax_aroot_radial_in(ilayer)
        else
           kmax_up = cohort_hydr%kmax_aroot_radial_out(ilayer)
        end if
        kmax_lo = site_hydr%kmax_upper_shell(ilayer,1)
        
        call GetImTaylorKAB(kmax_up,kmax_dn,       &
             ftc_node(i_up),ftc_node(i_dn),        & 
             h_node(i_up),h_node(i_dn),            & 
             dftc_dtheta(i_up), dftc_dtheta(i_dn), &
             dpsi_dtheta(i_up), dpsi_dtheta(i_dn), &
             k_eff(jpath),                         &
             A_term(jpath),                        & 
             B_term(jpath))
        

        ! Path is between rhizosphere shells
        
        do jpath = n_hypool_ag+3,n_hpool_tot-1

           i_dn = jpath
           i_up = jpath+1
           ishell_dn = i_dn - (n_hypool_ag+2)
           ishell_up = i_up - (n_hypool_ag+2)
           kmax_up = site_hydr%kmax_outer_shell(ilayer,ishell_up)
           kmax_lo = site_hydr%kmax_inner_shell(ilayer,ishell_lo)

           call GetImTaylorKAB(kmax_up,kmax_dn,       &
             ftc_node(i_up),ftc_node(i_dn),        & 
             h_node(i_up),h_node(i_dn),            & 
             dftc_dtheta(i_up), dftc_dtheta(i_dn), &
             dpsi_dtheta(i_up), dpsi_dtheta(i_dn), &
             k_eff(jpath),                         &
             A_term(jpath),                        & 
             B_term(jpath))


        end do

        ! -------------------------------------------------------------------------------
        ! Part 3.
        ! Loop through nodes again, build matrix
        ! -------------------------------------------------------------------------------

        do inode = 1,n_hypool_tot
           a(inode)
           b(inode)
           c(inode)
           r(inode) 
        end do
          


        return
      end subroutine ImTaylorSolverTermsCond1D
        
      ! =================================================================================

      subroutine GetImTaylorKAB(kmax_up,kmax_dn, &
                                  ftc_up,ftc_dn, &
                                  h_up,h_dn, &
                                  dftc_dtheta_up, dftc_dtheta_dn, &
                                  dpsi_dtheta_up, dpsi_dtheta_dn, &
                                  k_eff,   &
                                  A_term,  & 
                                  B_term)

          ! -----------------------------------------------------------------------------
          ! This routine will return the effective conductance "K", as well
          ! as two terms needed to calculate the implicit solution (using taylor
          ! first order expansion).  The two terms are generically named A & B.
          ! Thus the name "KAB".  These quantities are specific not to the nodes
          ! themselves, but to the path between the nodes, defined as positive
          ! direction from upstream node to downstream node.
          ! -----------------------------------------------------------------------------

          real(r8),intent(in) :: kmax_up, kmax_dn  ! max conductance [kg s-1 Mpa-1]
          real(r8),intent(in) :: ftc_up, ftc_dn    ! frac total conductance [-]
          real(r8),intent(in) :: h_up, h_dn        ! total potential [Mpa]
          real(r8),intent(in) :: dftc_dtheta_up, dftc_dtheta_dn ! Derivative
                                                                ! of FTC wrt relative water content
                                                 
          real(r8),intent(in) :: dpsi_dtheta_up, dpsi_dtheta_dn ! Derivative of matric potential
                                                                ! wrt relative water content

          real(r8),intent(in) :: k_eff                ! effective conductance over path [kg s-1 Mpa-1]
          real(r8),intent(in) :: a_term               ! "A" term for path (See tech note)
          real(r8),intent(in) :: b_term               ! "B" term for path (See tech note)


          ! Calculate total effective conductance over path  [kg s-1 MPa-1]
          k_eff = 1._r8/(1._r8/(ftc_dn*kmax_dn)+1._r8/(ftc_up*kmax_up))

          ! Calculate difference in total potential over the path [MPa]
          h_diff  = h_up - h_dn
          
          ! "A" term, which operates on the down-stream node
          A_term = (k_eff**2.0_r8) * h_diff * (kmax_dn**-1.0_r8) * (ftc_dn**-2.0_r8) &
               * dftc_dtheta_dn - k_eff * dpsi_dtheta_dn
          
          ! "B" term, which operates on the up-stream node
          B_term = (k_eff**2.0_r8) * h_diff * (kmax_up**-1.0_r8) * (ftc_up**-2.0_r8) & 
               * dftc_dtheta_up + k_eff * dpsi_dtheta_up
          


          return
        end subroutine GetImTaylorTerms



end module FatesHydroSolversMod
