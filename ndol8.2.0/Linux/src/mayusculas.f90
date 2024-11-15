    character(len=5) :: keyword1
    keyword1 = 'AaAaA'
    print *, keyword1
    call upcase(keyword1)
    print *, keyword1
    stop
    end
