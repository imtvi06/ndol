program ReadWords
    implicit none
    character(len=100) :: inputString
    character(len=100) :: word
    integer :: i, j

    ! Prompt the user for input
    print *, 'Enter a string (max 100 characters):'
    read *, inputString

    ! Initialize the index for the word extraction
    j = 1

    ! Loop through the input string
    do i = 1, len_trim(inputString)
        ! Check if the character is not a blank
        if (inputString(i:i) /= ' ') then
            word(j:j) = inputString(i:i)
            j = j + 1
        else
            ! If we encounter a blank, print the current word if j > 1
            if (j > 1) then
                print *, 'Word found:', trim(word(1:j-1))
                j = 1  ! Reset for the next word
            end if
        end if
    end do

    ! Print the last word if the string does not end with a blank
    if (j > 1) then
        print *, 'Word found:', trim(word(1:j-1))
    end if

end program ReadWords
