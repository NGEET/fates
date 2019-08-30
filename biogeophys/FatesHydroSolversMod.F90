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

      ccohort_hydr%kmax_petiole_to_leaf = 1.e8_r8


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

    subroutine ImTaylorSolverTermsCond1D(cohort_hydr,site_hydr,ilayer,dt_step,q_top,d_th_node,sapflow,rootuptake,wb_err)

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
        real(r8)          , intent(in)  :: flc_node(:)         ! fractional loss of conductivity at water storage nodes [-]
        real(r8)          , intent(in)  :: dt_step             ! time [seconds] over-which to calculate solution
        real(r8)          , intent(in)  :: q_top               ! transpiration flux rate at upper boundary [kg -s]
        real(r8),intent(out) :: d_th_node(n_hypool_tot)        ! change in theta over the timestep
        real(r8),intent(out) :: sapflow                        ! time integrated mass flux between transp-root and stem [kg]
        real(r8),intent(out) :: rootuptake                     ! time integrated mass flux between rhizosphere and aroot [kg]
        real(r8),intent(out) :: wb_err                         ! transpiration should match change in storage [kg]
        
        ! Locals
        
        integer :: inode   ! node index "i"
        integer :: jpath   ! path index "j"
        integer :: ishell  ! rhizosphere shell index of the node
        integer :: i_dn    ! downstream node of current flow-path
        integer :: i_up    ! upstream node of current flow-path
        integer :: iter    ! iteration count for sub-steps
        logical :: solution_found ! logical set to true if a solution was found within error tolerance
        real(r8) :: kmax_up  ! maximum conductance of the upstream half of path [kg s-1 Mpa-1]
        real(r8) :: kmax_dn  ! maximum conductance of the downstream half of path [kg s-1 MPa-1]
        real(r8) :: wb_step_err
        real(r8) :: wb_err                          ! sum of water balance error over substeps
        real(r8) :: leaf_water  ! kg of water in the leaf
        real(r8) :: stem_water  ! kg of water in the stem
        real(r8) :: root_water  ! kg of water in the transp and absorbing roots
        real(r8) :: th_node_init(n_hypool_tot)      ! "theta" i.e. water content of node [m3 m-3]
        real(r8) :: th_node(n_hypool_tot)
        real(r8) :: z_node(n_hypool_tot)            ! elevation of node [m]
        real(r8) :: v_node(n_hypool_tot)            ! volume of the node [m3]
        real(r8) :: psi_node(n_hypool_tot)          ! matric potential on node [Mpa]
        real(r8) :: ftc_node(n_hypool_tot)          ! frac total conductance on node [-]
        real(r8) :: h_node(n_hypool_tot)            ! total potential on node [Mpa]
        real(r8) :: dftc_dtheta_node(n_hypool_tot)  ! deriv FTC w.r.t. theta
        real(r8) :: dpsi_dtheta_node(n_hypool_tot)  ! deriv psi w.r.t. theta
        real(r8) :: k_eff(n_hypool_tot-1)           ! effective (used) conductance over path [kg s-1 MPa-1]
        real(r8) :: a_term(n_hypool_tot-1)          ! "A" term in the tri-diagonal implicit solve [-]
        real(r8) :: b_term(n_hypool_tot-1)          ! "B" term in the tri-diagonal implicit solve [-]
        real(r8) :: tris_a(n_hypool_tot)                   ! left of diagonal terms for tri-diagonal matrix solving delta theta
        real(r8) :: tris_b(n_hypool_tot)                   ! center diagonal terms for tri-diagonal matrix solving delta theta
        real(r8) :: tris_c(n_hypool_tot)       ! right of diaongal terms for tri-diagonal matrix solving delta theta
        real(r8) :: tris_r(n_hypool_tot) 

        integer, parameter  :: max_iter = 5
        real(r8), parameter :: max_wb_step_err = 1.e-6_r8 
        real(r8), parameter :: max_wb_err      = 1.e-4_r8  ! threshold for water balance error (stop model)   [mm h2o]

        ! -------------------------------------------------------------------------------
        ! Part 1.  Calculate node quantities:
        !          matric potential: psi_node
        !          fraction of total conductance: ftc_node
        !          total potential (matric + elevatio) h_node
        !          deriv. ftc  wrt  theta: dftc_dtheta_node
        !          deriv. psi  wrt  theta: dpsi_dtheta_node
        ! -------------------------------------------------------------------------------
        
     

        ! For all nodes leaf through rhizosphere
        ! Send node heights and compartment volumes to a node-based array
        
        do inode = 1,n_hypool_tot

           if (inode<=n_hypool_ag) then
              z_node(inode)  = ccohort_hydr%z_node_ag(inode)
              v_node(inode)  = ccohort_hydr%v_node_ag(inode)
              th_node_init(inode) = ccohort_hydr%th_ag(inode)
           elseif (inode==n_hypool_ag+1) then
              z_node(inode)  = ccohort_hydr%z_node_troot
              v_node(inode)  = ccohort_hydr%v_troot
              th_node_init(inode) = ccohort_hydr%th_troot
           elseif (inode==n_hyppol_ag+2) then
              z_node(inode)  = bc_in(s)%z_sisl(ilayer)
              v_node(inode)  = ccohort_hydr%v_aroot_layer(:)
              th_node_init(inode) = ccohort_hyd%th_aroot(ilayer)
           else
              ishell  = inode-(n_hypool_tot+2)
              z_node(inode)  = bc_in(s)%z_sisl(ilayer)
              v_node(inode)  = csite_hydr%v_shell(ilayer,ishell)
              th_node_init(inode) = site_hydr%h2osoi_liqvol_shell(ilayer,ishell)
           end if

        end do


        ! Outer iteration loop
        ! This cuts timestep in half and resolve the solution with smaller substeps
        ! This loop is cleared when the model has found a solution

        solution_found = .false.
        iter = 0
        do while( .not.solution_found ) 

           ! These are diagnostics that must be calculated.
           ! in this routine (uses differentials and actual fluxes)
           ! So we need to zero them, as they are incremented
           ! over the sub-steps
           
           sapflow        = 0._r8
           rootuptake     = 0._r8
           wb_err         = 0._r8

           ! Gracefully quit if too many iterations have been used
           if(iter>max_iter)then
              write(fates_log(),*) 'Could not find a stable solution for hydro 1D solve'
              write(fates_log(),*) ''
              leaf_water = sum(ccohort_hydr%th_ag(1:n_hypool_leaf)* &
                   ccohort_hydr%v_ag(1:n_hypool_leaf))*denh2o
              stem_water = sum(ccohort_hydr%th_ag(n_hypool_leaf+1:n_hypool_ag) * &
                   ccohort_hydr%v_ag(n_hypool_leaf+1:n_hypool_ag))*denh2o
              root_water = (ccohort_hydr%th_troot*ccohort_hydr%v_troot) + &
                    sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:))) * denh2o
              write(fates_log(),*) 'leaf water: ',leaf_water,' kg/plant'
              write(fates_log(),*) 'stem_water: ',stem_water,' kg/plant'
              write(fates_log(),*) 'root_water: ',root_water,' kg/plant'
              write(fates_log(),*) 'LWP: ',ccohort_hydr%psi_ag(1)
              write(fates_log(),*) 'dbh: ',ccohort%dbh
              write(fates_log(),*) 'pft: ',ccohort%pft
              write(fates_log(),*) 'tree lai: ',ccohort%treelai,' m2/m2 crown'
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if
           

           ! For each attempt, we want to reset theta with the initial value
           th_node(:) = th_node_init(:)


           ! Determine how many substeps, and how long they are

           nsteps = max(imult*iter,1)   ! Factor by which we divide through the timestep
                                          ! start with full step (ie dt_fac = 1)
                                          ! Then increase per the "imult" value.

           dt_substep = dt_step/real(nsteps,r8) ! This is the sub-stem length in seconds
           
           ! Walk through sub-steps
           do istep = 1,nsteps

              ! Total water mass in the plant at the beginning of this solve [kg h2o]
              w_tot_beg = sum(th_node(:)*v_node(:))*denh2o

              ! Calculate on-node quantities: potential, and derivatives
              do inode = 1,n_hypool_tot
              
                 ! Get matric potential [Mpa]
                 call psi_from_th(currentCohort%pft, porous_media(inode), th_node(inode), &
                      psi_node(inode), site_hydr, bc_in)
                 
                 ! Get total potential [Mpa]
                 h_node(inode) =  mpa_per_pa*denh2o*grav_earth*z_node(inode) + psi_node(inode)
              
                 ! Get Fraction of Total Conductivity [-]
                 call flc_from_psi(currentCohort%pft, porous_media(inode), psi_node(inode), &
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
              
              
              ! Path is between the transporting root 
              ! and the absorbing root for this layer
              
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
              
              tris_a(1) = 0._r8
              tris_b(1) = A_term(1) - denh20*vol_node(1)/dt_substep
              tris_c(1) = B_term(1)
              tris_r(1) = q_top - k_eff(1)*(h_node(2)-h_node(1))
              
              
              do inode = 2,n_hypool_tot-1
                 jpath = inode
                 tris_a(inode) = -A_term(jpath-1)
                 tris_b(inode) = A_term(jpath) - B_term(jpath-1) - denh2o*vol_node(inode)/dt_substep
                 tris_c(inode) = B_term(jpath)
                 tris_r(inode) = -k_eff(jpath)*(h_node(inode+1)-h_node(inode)) + &
                      k_eff(jpath-1)*(h_node(inode)-h_node(inode-1))
                 
              end do
              
              inode = n_hypool_tot
              jpath = n_hypool_tot
              tris_a(inode) = -A_term(jpath-1)
              tris_b(inode) = -B_term(jpath-1) - denh2o*vol_node(inode)/dt_substep
              tris_c(inode) = 0._r8
              tris_r(inode) = k_eff(jpath-1)*(h_node(inode)-h_node(inode-1))

              
              ! Calculate the change in theta
              
              call Hydraulics_Tridiagonal(tris_a, tris_b, tris_c, tris_r, dth_node)
              
              
              ! Catch super-saturated and sub-residual water contents

              ! Mass error (flux - change)
              ! Total water mass in the plant at the beginning of this solve [kg h2o]
              w_tot_end = sum((th_node(:)+dth_node(:))*v_node(:))*denh2o

              wb_step_err = (q_top*dt_substep) - (w_tot_beg-w_tot_end)

              if(abs(wb_step_err)>max_wb_step_err)then
                 solution_found = .false.
                 exit
              else
                 ! Note: this is somewhat of a default true. And the sub-steps
                 ! will keep going unless its changed and broken out of
                 ! the loop.
                 solution_found = .true.
              end if

              ! If we have not broken from the substep loop,
              ! that means this sub-step has been acceptable, and we may 
              ! go ahead and update the water content for the integrator

              th_node(:) = th_node(:) + dth_node(:)

              ! Accumulate the water balance error for diagnostic purposes
              wb_err = wb_err + wb_step_err

              ! -------------------------------------------------------------------------
              ! Diagnostics
              ! -------------------------------------------------------------------------

              ! Sapflow at the base of the tree is the flux rate
              ! between the transporting root node and the first stem node
              ! (note: a path j is between node i and i+1)
              ! [kg] = [kg/s] * [s]
              
              inode = n_hypool_ag
              sapflow = sapflow + dt_substep * & 
                   (k_eff(inode)*(h_node(inode+1)-h_node(inode)) + &  ! flux at (t) 
                   A_term(inode)*dth_node(inode)                 + &  ! dq at node i
                   B_term(inode)*dth_node(inode+1))                   ! dq at node i+1
              
              ! Root uptake is the integrated flux between the first rhizosphere
              ! shell and the absorbing root

              inode = h_hypool_ag+2
              rootuptake = rootuptake + dt_substep * & 
                   (k_eff(inode)*(h_node(inode+1)-h_node(inode)) + &  ! flux at (t) 
                   A_term(inode)*dth_node(inode)                 + &  ! dq at node i
                   B_term(inode)*dth_node(inode+1))                   ! dq at node i+1
              
              
           end do  ! do istep = 1,nsteps  (substep loop)
           
           iterh1=iterh1+1
           
        end do

        ! Save the number of times we refined our sub-step counts (iterh1)
        ccohort_hydr%iterh1 = real(iterh1)
        ! Save the number of sub-steps we ultimately used
        ccohort_hydr%iterh2 = real(nsteps)

        ! -----------------------------------------------------------
        ! To a final check on water balance error sumed over sub-steps
        ! ------------------------------------------------------------
        if ( abs(wb_err) > max_wb_err ) then

            write(fates_log(),*)'EDPlantHydraulics water balance error exceeds threshold of = ', max_wb_err
            write(fates_log(),*)'transpiration demand: ', dtime*qtop,' kg/step/plant'

            leaf_water = ccohort_hydr%th_ag(1)*ccohort_hydr%v_ag(1)*denh2o
            stem_water = sum(ccohort_hydr%th_ag(2:n_hypool_ag) * &
                  ccohort_hydr%v_ag(2:n_hypool_ag))*denh2o
            root_water = ( ccohort_hydr%th_troot*ccohort_hydr%v_troot + &
                  sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:))) * denh2o
            
            write(fates_log(),*) 'leaf water: ',leaf_water,' kg/plant'
            write(fates_log(),*) 'stem_water: ',stem_water,' kg/plant'
            write(fates_log(),*) 'root_water: ',root_water,' kg/plant'
            write(fates_log(),*) 'LWP: ',ccohort_hydr%psi_ag(1)
            write(fates_log(),*) 'dbh: ',ccohort%dbh
            write(fates_log(),*) 'pft: ',ccohort%pft
            write(fates_log(),*) 'tree lai: ',ccohort%treelai,' m2/m2 crown'
            call endrun(msg=errMsg(sourcefile, __LINE__))
        end if

        
        ! Adjust final water balance by adding back in the error term
        ! ------------------------------------------------------------
        
        w_tot_end_outer   = sum(th_node(:)*v_node(:))*denh2o                            ! kg
        dw_tot_outer      = w_tot_end_outer - w_tot_beg_outer                           ! kg/timestep
        we_tot_outer      = dw_tot_outer + (qtop_dt + dqtopdth_dthdt)                   ! kg/timestep
        we_area_outer     = we_tot_outer/(cCohort%c_area / cCohort%n)                   ! kg/m2 ground/individual
        if(abs(we_tot_outer*cCohort%n)/AREA>1.0e-7_r8) then
            if(debug) then
                write(fates_log(),*)'WARNING: plant hydraulics water balance error exceeds 1.0e-7 and is ajusted for error'
            endif
            !dump the error water to the bin with largest water storage
            max_l  = maxloc(th_node(:)*v_node(:),dim=1)
            th_node(max_l) = th_node(max_l)-  &
                  we_tot_outer/(v_node(max_l)*denh2o)
            th_node(max_l) = min (th_node(max_l),&
                  ths_node(max_l)-small_theta_num) 
            th_node(max_l) = max(th_node(max_l),&
                  thr_node(max_l)+small_theta_num)	  
            w_tot_end_outer   = sum(th_node(:)*v_node(:))*denh2o                            ! kg
            dw_tot_outer      = w_tot_end_outer - w_tot_beg_outer                           ! kg/timestep
            we_tot_outer      = dw_tot_outer + (qtop_dt + dqtopdth_dthdt)                   ! kg/timestep
            we_area_outer     = we_tot_outer/(cCohort%c_area / cCohort%n)                   ! kg/m2 ground/individual   
        end if
        


        ! If we have made it to this point, supposedly we have completed the whole time-step
        ! for this cohort x layer combination.  It is now safe to save the delta theta
        ! value and pass it back to the calling routine.  The value passed back is the
        ! change in theta over all sub-steps.
        
        dth_node(:) = th_node(:)-th_node_init(:)
        


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
