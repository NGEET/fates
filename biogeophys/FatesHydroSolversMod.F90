module FatesHydroSolversMod



contains




    subroutine HydraulicsMatrixSolvePHS(   )








    end subroutine HydraulicsMatrixSolvePHS

    subroutine HDiffK1D(cohort_hydr,site_hydr,inodes,psi_node,flc_node,dflcdpsi_node, &
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

            if(inode_up < n_hypool_ag+n_hypool_troot+1) then
                ! Path is between compartments within the plant
                ! (leaves,stems,transporting roots and absorbing root)

                znode_up = cohort_hydr%z_node(inode_up)
                znode_lo = cohort_hydr%z_node(inode_lo)
                psinode_up = psi_node(inode_up)
                psinode_lo = psi_node(inode_lo)
                kmax_up   = 
                kmax_lo   = 
                kmax_surf = 1.e20_r8  ! There are no surface conductances

            elseif(inode_up == n_hpool_ag+n_hypool_troot+1) then
                ! Path is between the absorbing root and the 1st
                ! rhizosphere shell compartment

                znode_up = bc_in(s)%z_sisl(ilayer)
                znode_lo = bc_in(s)%z_sisl(ilayer)
                psinode_up = psi_node(inode_up)
                psinode_lo = psi_node(inode_lo)

                kmax_surf = 
                kmax_up=kmax_bound
                kmax_lo=site_hydr%kmax_bound_shell(ilayer,inode_lo)

            else
                ! Path is between rhizosphere shells

                znode_up = bc_in(s)%z_sisl(ilayer)
                znode_lo = bc_in(s)%z_sisl(ilayer)
                psinode_up = psi_node(inode_up)
                psinode_lo = psi_node(inode_lo)

                kmax_up=site_hydr%kmax_bound_shell(ilayer,inode_up)
                kmax_lo=site_hydr%kmax_bound_shell(ilayer,inode_lo)

            end if

            hdiff_bound(jpath) = mpa_per_pa*denh2o*grav_earth*(znode_up-znode_lo) + (psinode_up-psinode_lo)


            ! examine direction of water flow; use the upstream node's k for the boundary k.
            ! (as suggested by Ethan Coon, LANL)
            if(do_kbound_upstream) then

                if(hdiff_bound(jpath) < 0._r8) then
                    ! More potential in the lower node, use its fraction of conductivity loss
                    k_bound(jpath)       = flc_node(inode_lo) / &
                          (1._r8/k_bound_aroot_soil1 + 1._r8/k_bound_aroot_soil2) * flc_node(k+1)  ! water moving towards atmosphere
                    dkdpsi0(jpath)  = 0._r8
                    dkdpsi1(jpath)  = kmax_bound(jpath) * dflcdpsi_node(inode_lo)

                else


                    



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


        ! ------------------------------------------------------------------------------
        ! Part II. Set maximum (size-dependent) hydraulic conductances
        ! ------------------------------------------------------------------------------
        
        ! first estimate cumulative (petiole to node k) conductances 
        ! without taper as well as the chi taper function
        
       do k=n_hypool_leaf,n_hypool_ag
          dz_node1_lowerk          = ccohort_hydr%z_node_ag(n_hypool_leaf) &
               - ccohort_hydr%z_lower_ag(k)
          if(k < n_hypool_ag) then
             dz_node1_nodekplus1   = ccohort_hydr%z_node_ag(n_hypool_leaf) &
                  - ccohort_hydr%z_node_ag(k+1)
          else
             dz_node1_nodekplus1   = ccohort_hydr%z_node_ag(n_hypool_leaf) &
                  - ccohort_hydr%z_node_troot(1)
          end if
          kmax_node1_nodekplus1(k) = EDPftvarcon_inst%hydr_kmax_node(ft,2) * a_sapwood / dz_node1_nodekplus1
          kmax_node1_lowerk(k)     = EDPftvarcon_inst%hydr_kmax_node(ft,2) * a_sapwood / dz_node1_lowerk
          chi_node1_nodekplus1(k)  = xylemtaper(taper_exponent, dz_node1_nodekplus1)
          chi_node1_lowerk(k)      = xylemtaper(taper_exponent, dz_node1_lowerk)
          if(.not.do_kbound_upstream) then
             if(crown_depth == 0._r8) then 
                write(fates_log(),*) 'do_kbound_upstream requires a nonzero canopy depth '
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
          end if
       enddo


       ! then calculate the conductances at node boundaries as the difference of cumulative conductances
       do k=n_hypool_leaf,n_hypool_ag
          if(k == n_hypool_leaf) then
             ccohort_hydr%kmax_bound(k)    = kmax_node1_nodekplus1(k)  * chi_node1_nodekplus1(k)
          else
             ccohort_hydr%kmax_bound(k)    = ( 1._r8/(kmax_node1_nodekplus1(k)  *chi_node1_nodekplus1(k)  ) - &
                  1._r8/(kmax_node1_nodekplus1(k-1)*chi_node1_nodekplus1(k-1))     ) ** (-1._r8)
          end if

       enddo

       ! finally, estimate the remaining tree conductance belowground as a residual

       kmax_treeag_tot = sum(1._r8/ccohort_hydr%kmax_bound(n_hypool_leaf:n_hypool_ag))** (-1._r8)

       kmax_tot        = EDPftvarcon_inst%hydr_rfrac_stem(ft)* kmax_treeag_tot

       ccohort_hydr%kmax_treebg_tot      = ( 1._r8/kmax_tot - 1._r8/kmax_treeag_tot ) ** (-1._r8)

       do j=1,nlevsoi_hyd
           if(j == 1) then
               rootfr = zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j))
           else
               rootfr = zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j)) - &
                     zeng2001_crootfr(roota, rootb, bc_in%zi_sisl(j-1))
           end if
           ccohort_hydr%kmax_treebg_layer(j) = rootfr*ccohort_hydr%kmax_treebg_tot
       end do




    end subroutine PlantKmax





end module FatesHydroSolversMod
