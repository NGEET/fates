module FatesLitterMod
   
   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : i4 => fates_int
   use FatesConstantsMod, only : nearzero
   use FatesConstantsMod, only : calloc_abs_error
   use FatesGlobals     , only : endrun => fates_endrun
   use FatesGlobals     , only : fates_log 
   use shr_log_mod      , only : errMsg => shr_log_errMsg

   implicit none
   private
   
   
   ! Section 1: Public Definitions

   ! -------------------------------------------------------------------------------------
   ! Organ types
   ! These are public indices used to map the organs
   ! in each hypothesis to organs that acknowledged in the calling model
   ! -------------------------------------------------------------------------------------
   
   integer, public, parameter :: num_litt_types = 5
   integer, public, parameter :: all_litt     = 0    ! index for all litter pools
   integer, public, parameter :: leaf_litt    = 1    ! index for leaf litter
   integer, public, parameter :: fnrt_litt    = 2    ! index for fine root litter
   integer, public, parameter :: acwd_litt    = 3    ! index for above-ground coarse woody debris
   integer, public, parameter :: bcwd_litt    = 4    ! index for below-ground coarse woody debris
   integer, public, parameter :: repr_litt    = 5    ! index for reproductive tissues litter
   
   
   ! Section 2: Generic Types
   
   ! Each litter variable (state and fluxes) is dimensioned by (PFT x spatial coordinate)
   type litt_vartype
      
      real(r8),allocatable :: val(:,:)        ! Instantaneous state variable              [kg]
      real(r8),allocatable :: val0(:,:)       ! State variable at the beginning of step   [kg]
      real(r8),allocatable :: fluxin(:,:)     ! Input fluxes generated from turnover in   
                                              ! living plants [kg]
      real(r8),allocatable :: burned(:,:)     ! Output flux from burning
      real(r8),allocatable :: fraged(:,:)     ! Output flux from fragmentation
      real(r8),allocatable :: germed(:,:)     ! Output loss from germination
      real(r8),allocatable :: herbed(:,:)     ! Output flux from herbivory
      
   end type litt_vartype

   
   type litt_type
      
      type(litt_vartype),allocatable :: variables(:)    ! The state variables and fluxes
      
   contains
      
      ! These are extendable procedures that have specialized
      ! content in each of the different hypotheses
      
      procedure :: InitLitt   => InitLittBase
      procedure :: DailyLitt  => DailyLittBase
      procedure :: FastLitt   => FastLittBase
      
      ! These are generic functions that should work on all hypotheses
      
      procedure, non_overridable :: InitAllocate
      
   end type litt_type

   type litt_global_type


      


   end type litt_global_type



   ! Part 3: Public extended types

   type, public, extends(litt_vartype) :: c_litt_vartype
      
   contains
      
      procedure :: InitLitt   => InitLittCarbon
      procedure :: DailyLitt  => DailyLittCarbon
      procedure :: FastLitt   => FastLittCarbon
      
   end type c_litt_vartype
   

   character(len=*), parameter, private :: sourcefile = __FILE__

   
contains




end module FatesLitterMod
