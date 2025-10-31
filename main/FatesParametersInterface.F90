module FatesParametersInterface

  ! This module simply holds the instantiation
  ! of the generalized FATES parameter file that
  ! is provide by the JSON parser. ie pstruct
  ! This module is intentionaly low dependency
  ! Note also that JSONParameterUtilsMod
  ! only uses shr libraries

  use JSONParameterUtilsMod, only: params_type

  implicit none
  private
  
  type(params_type) :: pstruct

  public :: pstruct
  public :: Transp2dInt
  public :: Transp2dReal

contains


  subroutine Transp2dInt(i_2d_in,i_2d_out)

    ! The FATES JSON parameter files have a legacy from the netcdf 
    ! of having the pfts as the column (inner) indices. At least
    ! that is how it is presented in the text files.
    ! [ [pft1, pft2, pft3],[pf1, pft2, pft3] ]
    !
    ! In fortran, this would have the pft index as the second index
    ! because it goes row x column.
    !
    ! However, we the data arrays in fates use
    ! the PFT as the first dimension... so:
    
    integer :: i_2d_in(:,:)
    integer :: i_2d_transp(:,:)
    integer :: i,j

    do i=1,size(i_2d_in,dim=1)
       do j = 1,size(i_2d_in,dim=2)
          i_2d_out(j,i) = i_2d_in(i,j)
       end do
    end do
    
  end subroutine Transp2dInt
  
  subroutine Transp2dReal(r_2d_in,r_2d_out)

    ! The FATES JSON parameter files have a legacy from the netcdf 
    ! of having the pfts as the column (inner) indices. At least
    ! that is how it is presented in the text files.
    ! [ [pft1, pft2, pft3],[pf1, pft2, pft3] ]
    !
    ! In fortran, this would have the pft index as the second index
    ! because it goes row x column.
    !
    ! However, we the data arrays in fates use
    ! the PFT as the first dimension... so:
    
    integer :: r_2d_in(:,:)
    integer :: r_2d_transp(:,:)
    integer :: i,j
    
    do i=1,size(r_2d_in,dim=1)
       do j = 1,size(r_2d_in,dim=2)
          r_2d_out(j,i) = r_2d_in(i,j)
       end do
    end do

  end subroutine Transp2dReal
  
end module FatesParametersInterface

