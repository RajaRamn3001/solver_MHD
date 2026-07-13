module array_processes
contains
  SUBROUTINE QUICKSORT(ARRAY, N)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                            :: N
    INTEGER, INTENT(INOUT)                         :: ARRAY(N)
    
    ! Stack for storing partition boundaries
    INTEGER, PARAMETER :: MAXSTACK = 100
    INTEGER :: STACK(MAXSTACK, 2)
    INTEGER :: TOP, LOW, HIGH, PIVOT
    INTEGER :: I, J, TEMP
    
    ! Initialize stack
    TOP = 1
    STACK(1, 1) = 1
    STACK(1, 2) = N
    
    ! Main sorting loop
    DO WHILE (TOP > 0)
        ! Pop from stack
        LOW = STACK(TOP, 1)
        HIGH = STACK(TOP, 2)
        TOP = TOP - 1
        
        ! Only sort if more than one element
        IF (LOW < HIGH) THEN
            ! Partition the array
            CALL PARTITION(ARRAY, LOW, HIGH, PIVOT)
            
            ! Push left partition onto stack
            IF (PIVOT - 1 > LOW) THEN
                TOP = TOP + 1
                STACK(TOP, 1) = LOW
                STACK(TOP, 2) = PIVOT - 1
            END IF
            
            ! Push right partition onto stack
            IF (PIVOT + 1 < HIGH) THEN
                TOP = TOP + 1
                STACK(TOP, 1) = PIVOT + 1
                STACK(TOP, 2) = HIGH
            END IF
        END IF
    END DO
 END SUBROUTINE QUICKSORT

 SUBROUTINE PARTITION(ARRAY, LOW, HIGH, PIVOT_INDEX)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT)                             :: ARRAY(*)
    INTEGER, INTENT(IN)                                :: LOW, HIGH
    INTEGER, INTENT(OUT)                               :: PIVOT_INDEX
    
    INTEGER :: PIVOT_VALUE, I, J, TEMP
    
    ! Choose last element as pivot
    PIVOT_VALUE = ARRAY(HIGH)
    I = LOW - 1
    
    ! Partition loop
    DO J = LOW, HIGH - 1
        IF (ARRAY(J) <= PIVOT_VALUE) THEN
            I = I + 1
            ! Swap elements
            TEMP = ARRAY(I)
            ARRAY(I) = ARRAY(J)
            ARRAY(J) = TEMP
        END IF
    END DO
    
    ! Place pivot in correct position
    TEMP = ARRAY(I + 1)
    ARRAY(I + 1) = ARRAY(HIGH)
    ARRAY(HIGH) = TEMP
    
    PIVOT_INDEX = I + 1
 END SUBROUTINE PARTITION

  subroutine append_to_array(array, new_element)
    integer, allocatable, intent(inout) :: array(:)
    integer, intent(in) :: new_element
    integer, allocatable :: temp_array(:)
    integer :: n
    
    if (.not. allocated(array)) then
        allocate(array(1))
        array(1) = new_element
    else
        n = size(array)
        allocate(temp_array(n+1))
        temp_array(1:n) = array(1:n)
        temp_array(n+1) = new_element
        call move_alloc(temp_array, array)
    endif
  end subroutine
end module array_processes
