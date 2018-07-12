program main
    use mydata
    implicit none

    logical, parameter :: master = .false.

    call init_world(master)
    call create_pools(master)

    call processdata()
    
    call mpi_finalize()

    contains
    subroutine processdata()
    integer(kind=i4b) :: recvbuf(2)
    integer(kind=i4b) :: sendbuf
    character(len=4) :: ftype

    !variables for managing the mpi part
    type(mpi_status) :: my_status
    type(mpi_request) :: my_request

    infty_loop: do
        !wait for handshake
        !write(*,'(A,I2,1x,A,I2,1x,A)') 'worker(',my_master,'): my_rank=',my_rank,'WAITING'
        call mpi_recv(ftype, 4, mpi_character, my_master, &
                      tag_ftype, my_pool, my_status)
        !write(*,'(A,I2,1x,A,I2,1x,A,A)') 'worker(',my_master,'): my_rank=',my_rank,'ftype=',ftype

        if (ftype .eq. 'func') then
            !initial blocking receive from its own master
            call mpi_recv(recvbuf, 2, mpi_integer, my_master, &
                          tag_indata, my_pool, my_status)
            !write(*,'(A,I2,1x,A,I2,1x,A,2(I0,",",1x))') 'worker(',my_master,'): my_rank=',my_rank,'recvbuf=',recvbuf

            !process data, emulate long process
            sendbuf = sum(recvbuf)
            call sleep(1)

            !tell the master that you're done 
            call mpi_send(1, 1, mpi_integer, my_master, &
                           tag_oudata, my_pool)

            !now send the results
            call mpi_send(sendbuf, 1, mpi_integer, my_master, &
                          tag_indata, my_pool)

            !write(*,'(A,I2,1x,A,I2,1x,A,I0)') 'worker(',my_master,'): my_rank=',my_rank,'sendbuf=',sendbuf
        else
            exit infty_loop
        endif
    enddo infty_loop
    end subroutine processdata
end program main
