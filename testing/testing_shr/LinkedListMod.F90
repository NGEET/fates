module LinkedListMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesCohortMod,    only : fates_cohort_type
  use FatesPatchMod,     only : fates_patch_type

  implicit none 
  private
  
  public sort_cohorts
  
  contains 
  
  subroutine sort_cohorts(patch)
    !
    ! DESCRIPTION: sort cohorts in patch's linked list
    ! uses insertion sort to build a new list
    !
  
    ! ARGUMENTS:
    type(fates_patch_type), intent(inout) :: patch ! patch
    
    ! LOCALS:
    type(fates_cohort_type), pointer :: currentCohort
    type(fates_cohort_type), pointer :: nextCohort
    type(fates_cohort_type), pointer :: sorted_head
    type(fates_cohort_type), pointer :: currentSorted

    ! initialize
    sorted_head => null()
    currentCohort => patch%shortest
    
    do while (associated(currentCohort))
        
      ! store the next cohort to sort
      nextCohort => currentCohort%taller
      
      ! disconnect from list
      currentCohort%taller => null()
      currentCohort%shorter => null()
      
      ! insert into the sorted part of list 
      if (.not. associated(sorted_head)) then
      
        ! sorted is null, insert at head of list
        sorted_head => currentCohort
        
      else if (sorted_head%height >= currentCohort%height) then
      
        ! insert at head of list
        currentCohort%taller => sorted_head
        sorted_head%shorter => currentCohort
        sorted_head => currentCohort
              
      else
        
        ! traverse sorted list to find where to place
        currentSorted => sorted_head
        do while (associated(currentSorted%taller))
          if (currentCohort%height < currentSorted%taller%height) exit
          currentSorted => currentSorted%taller
        end do
        
        ! insert cohort after currentSorted
        currentCohort%taller => currentSorted%taller
        
        ! set shorter cohort if it is not inserted at the end
        if (associated(currentSorted%taller)) then
          currentSorted%taller%shorter => currentCohort
        end if
        
        ! set taller to currentCohort
        currentSorted%taller => currentCohort
        currentCohort%shorter => currentSorted
      end if
      
      currentCohort => nextCohort
    end do

    ! update patch pointers for shortest and tallest
    patch%shortest => sorted_head
    currentCohort => sorted_head
    do while (associated(currentCohort%taller))
      currentCohort => currentCohort%taller
    end do
    patch%tallest => currentCohort
    
  end subroutine sort_cohorts
  
end module LinkedListMod


