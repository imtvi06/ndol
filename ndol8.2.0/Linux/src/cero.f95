PROGRAM CommonArrayExample
    IMPLICIT NONE
    INTEGER :: i, arr
    COMMON arr(10)  ! Declare a COMMON integer array of size 10

    ! Initialize the array to zero
    arr = 0

    ! Print the array to verify initialization
    PRINT *, arr

END PROGRAM CommonArrayExample
