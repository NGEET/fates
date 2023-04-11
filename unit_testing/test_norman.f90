program test_norman
    !
    ! DESCRIPTION:
    !		Test the FATES radiation schemes
    !

    use IOMod
    use OrbitalValsMod, only : SHR_KIND_R8, get_orbital_vals, shr_orb_cosz,    &
                               shr_orb_decl

    implicit none

    ! Data dictionary: declare variable types, definitions, and units
    real(r8) :: rhol(num_pft,num_swb)   ! leaf reflectance (0-1)
    real(r8) :: rhos(num_pft,num_swb)   ! stem reflectance (0-1)
    real(r8) :: taul(num_pft,num_swb)   ! leaf transmittance (0-1)
    real(r8) :: taus(num_pft,num_swb)   ! stem transmittance (0-1)
    real(r8) :: xl(num_pft)             ! leaf/stem orientation index (-1-1)
    real(r8) :: clumping_index(num_pft) ! clumping index (0-1)
    real(r8) :: fcansno                 ! fraction of canopy covered by snow
    real(r8) :: canopy_area_profile(num_can,num_pft,nlevleaf) ! fraction of crown area per canopy area in each layer
    real(r8) :: elai_profile(num_can,num_pft,nlevleaf)       ! exposed leaf area in each canopy layer, pft, and leaf layer
    real(r8) :: esai_profile(num_can,num_pft,nlevleaf)       ! exposed stem area in each canopy layer, pft, and leaf layer
    real(r8) :: nrad_r(num_can,num_pft)   ! number of exposed leaf layers for each canopy layer and pft
    integer  :: nrad(num_can,num_pft)   ! number of exposed leaf layers for each canopy layer and pft
    real(r8) :: k_dir(num_pft,48)          ! direct beam extinction coefficient
    real(SHR_KIND_R8) :: eccen    ! orbital eccentricity
    real(SHR_KIND_R8) :: mvelp    ! moving vernal equinox long
    real(SHR_KIND_R8) :: obliqr   ! Earths obliquity in rad
    real(SHR_KIND_R8) :: lambm0   ! Mean long of perihelion at vernal equinox (radians)
    real(SHR_KIND_R8) :: mvelpp   ! moving vernal equinox long of perihelion plus pi (rad)
    real(SHR_KIND_R8) :: eccf     ! Earth-sun distance factor (ie. (1/r)**2)
    real(SHR_KIND_R8) :: calday   ! calendar day (including fraction)
    real(SHR_KIND_R8) :: declin   ! solar declination (radians)
    real(SHR_KIND_R8) :: cosz(48) ! cosine of solar zenith angle (radians)
    real(r8)          :: lat, lon
    integer           :: jday     ! julian day
    integer           :: year     ! year
    integer           :: i        ! looping index
    character(len=100) :: file_in, patch_file, out_file

    ! file_names
    file_in = "fates_params.nc"
    patch_file = "patch_data.nc"
    out_file = "nrad_out.nc"
    logf = open_file("log.txt")

    ! set julian day, year, and lat/lon here
    jday = 165
    year = 2000
    lat = 45.76_r8*pi_const/180.0_r8
    lon = 237.67_r8*pi_const/180.0_r8
    fcansno = 0.0_r8

    ! get patch and parameter values, as well as orbital parameters
    call get_parameters(file_in, rhol, rhos, taul, taus, xl, clumping_index)
    call read_patch_data(patch_file, canopy_area_profile, elai_profile,        &
        esai_profile, nrad_r)
    nrad = int(nrad_r)
    call get_orbital_vals(year, logf, eccen, mvelp, obliqr, lambm0, mvelpp)

    cosz(:) = 0.0
    k_dir(:,:) = 0.0
    
    ! for each half-hourly time step in the day
    do i = 0, 47
        calday = jday + i*0.02083333_SHR_KIND_R8
        call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, declin, eccf)
        cosz(i) = shr_orb_cosz(calday, lat, lon, declin)

        if (cosz(i) > 0.0_r8) then 
            ! call norman radiation scheme
            call PatchNormanRadiation(rhol, rhos, taul, taus, xl, clumping_index, &
                canopy_area_profile, elai_profile, esai_profile, fcansno,         &
                cosz(i), nrad, 2, k_dir(:,i))
        end if
    end do 

    call write_radiation_data(out_file, k_dir, cosz)

end program test_norman