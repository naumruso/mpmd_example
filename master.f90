program main
    use mydata
    implicit none

    logical, parameter :: master = .true.

    integer(kind=i4b), allocatable :: x(:,:)
    integer(kind=i4b) :: n = 2
    integer(kind=i4b) :: m = 100
    integer(kind=i4b) :: i
    integer(kind=i4b), parameter :: nxsum = 5
    integer(kind=i4b) :: xsum(nxsum)

    character(len=100) :: strformat

    call init_world(master)
    call create_pools(master)
    write(strformat,'(A,I0,A)') '(A,I2,1x,A,',size(my_workers),'(I0,",",1x))'
    write(*,trim(strformat)) 'master(', my_rank, ') my_workers=', my_workers

    !load test data
    allocate(x(n,m))
    do i=1, m, 1
        x(:,i) = i
    enddo

    !test subroutine
    call masterfunc(x, n, m, xsum, nxsum)
    if (my_rank .eq. main_proc) then
        write(strformat, '(A)') '(A,I2,1x,A,I0,1x,A,I0)'
        write(*,trim(strformat)) 'mfunc:my_rank=',my_rank,'xsum=',xsum(1),'sum(x)=',sum(x)
    endif

    !shutdown all workers
    call exitworkers()

    !quit yourself
    call mpi_finalize()

    contains
    subroutine masterfunc(x, n, m, xsum, nxsum)
    integer(kind=i4b), intent(in) :: n, m
    integer(kind=i4b), intent(in) :: x(n, m)
    integer(kind=i4b), intent(in) :: nxsum
    integer(kind=i4b), intent(inout) :: xsum(nxsum)

    integer(kind=i4b) :: p
    integer(kind=i4b) :: irank
    integer(kind=i4b), allocatable :: idx(:)
    integer(kind=i4b), allocatable :: sizes(:)
    integer(kind=i4b), allocatable :: xrecv(:,:)
    integer(kind=i4b), allocatable :: xres(:)

    type(mpi_status) :: my_status

    !contains the column indices of the main array
    !their indices start from zero because process
    !ranks go from zero on nmasters - 1
    allocate(idx(0:nmasters - 1))
    allocate(sizes(0:nmasters - 1))

    !temp. array that keeps results from each master
    allocate(xres(nxsum))

    !set how many columns we need to send to the masters
    !the last process might have to work over larger number
    !of columns than the rest of the workers
    p = m/nmasters

    !index of the start column for each master
    do i=0, nmasters - 1, 1
        idx(i) = i*p + 1
        sizes(i) = p
    enddo

    !check that the last process gets everything
    !to the end of the input array
    if (sum(sizes) .lt. m) sizes(nmasters - 1) = m - (nmasters - 1)*p

    !send or receive the initial data
    if (my_rank .eq. main_proc) then
        !Start the initial send from the main process to the other masters
        !
        !I do blocking sends, but it is also
        !okay to do non-blocking sends. All
        !results should be completely the same, but
        !we'd need to use mpi_wait
        do irank=1, nmasters - 1, 1
            call mpi_send(x(1,idx(irank)), sizes(irank)*n, mpi_integer, &
                           irank, tag_ftype, master_pool)
        enddo

        !the master process doesn't copy the data into a
        !temporary array, it sends the whole array to its workers
        !xrecv = x(1:n,idx(my_rank):sizes(my_rank))
        !write(*,*) 'my_rank=',my_rank,'xrecv=',xrecv
        call workerfunc(x(1,idx(my_rank)), n, sizes(my_rank), xres, nxsum)
        !write(*,*) 'wfunc:my_rank=',my_rank,'xres=',xres
        !call workerfunc(xrecv, n, sizes(my_rank), xres, nxsum)
    else
        !allocate the receive buffer for the input data
        !for more complicated setups this can be done only once
        allocate(xrecv(n,sizes(my_rank)))

        !blocking receives
        call mpi_recv(xrecv, n*sizes(my_rank), mpi_integer, main_proc, &
                      tag_ftype, master_pool, my_status)

        !call the subroutine that does the work
        !individually with each worker
        !write(*,*) 'my_rank=',my_rank,'xrecv=',xrecv
        call workerfunc(xrecv, n, sizes(my_rank), xres, nxsum)
        !write(*,*) 'wfunc:my_rank=',my_rank,'xres=',xres
    endif
    !write(*,*) 'my_rank=',my_rank,'xres=',xres

    !synchronize all masters
    call mpi_barrier(master_pool)

    !start post-processing of the results
    !Here we calculate the sum of the xres vectors
    !element by element using mpi_reduce.
    if (my_rank .eq. main_proc) then
        call mpi_reduce(xres, xsum, nxsum, mpi_integer, MPI_SUM, main_proc, master_pool)
    else
        call mpi_reduce(xres, xsum, nxsum, mpi_integer, MPI_SUM, main_proc, master_pool)
    endif
    end subroutine masterfunc


    subroutine workerfunc(x, n, m, xsum, nxsum)
    integer(kind=i4b), intent(in) :: n, m
    integer(kind=i4b), intent(in) :: x(n, m)
    integer(kind=i4b), intent(in) :: nxsum
    integer(kind=i4b), intent(out) :: xsum(nxsum)

    !array for keeping results from each worker
    integer(kind=i4b), allocatable :: xres(:)

    !variables for managing the mpi part
    integer(kind=i4b) :: itask, irank, irequest
    integer(kind=i4b) :: nmintasks, ntasks
    type(mpi_status) :: my_status
    type(mpi_request), allocatable :: my_requests(:)
    integer(kind=i4b), allocatable :: my_tasks(:)
    integer(kind=i4b) :: recvbuf

    !same as ntasks = m, this is the total
    !number of tasks that need to be processed
    !by this master
    ntasks = size(x, dim=2)
    allocate(xres(ntasks)) !array for the results

    !Below is an implementation of a basic task manager
    !Step 1: Load the workers with initial tasks
    !Step 2: Cycle over the tasks and the free workers
    !Step 3: Quit when all tasks have been executed
    !write(*,*) 'my_rank=',my_rank,'xrecv=',x,'ntasks=',ntasks

    !We cannot have more than nmintasks workers
    !doing calculations at any time
    nmintasks = min(ntasks, size(my_workers))
    allocate(my_requests(nmintasks))
    allocate(my_tasks(nmintasks))

    !Load the workers with the minimal number of tasks
    do itask=1, nmintasks, 1
        !on reply this will tell us where in xres to put the result
        my_tasks(itask) = itask

        !the rank of the worker with index itask
        irank = my_workers(itask)

        !initiate the handshake
        call mpi_send('func', 4, mpi_character, irank, tag_ftype, my_pool)
        !write(*,'(A,I0,1x,A,I0)') 'my_rank=',my_rank,'handshake to ', irank

        !send initial data to the worker
        call mpi_send(x(1,itask), 2, mpi_integer, irank, tag_indata, my_pool)
        !write(*,'(A,I0,1x,A,I0)') 'my_rank=',my_rank,'data sent to ', irank

        !request results using non-blocking receive.
        !Finished requests are processed below with mpi_waitany.
        call mpi_irecv(recvbuf, 1, mpi_integer, irank, tag_oudata,&
                      my_pool, my_requests(itask))
        !write(*,'(A,I0,1x,A,I0)') 'my_rank=',my_rank,'data req. to ', irank
    enddo
    !write(*,*) 'my_itask=',itask


    !Loop over each task as long as there is something to calculate
    !If all tasks have finished,i.e., aren't resupplied, then this
    !call returns MPI_UNDEFINED in itask
    call mpi_waitany(nmintasks, my_requests, irequest, my_status)
    !write(*,*) 'recvbuf=',recvbuf
    !write(*,*) 'mpi_waitany', irequest, MPI_UNDEFINED
    do while (irequest .ne. MPI_UNDEFINED)
        !The results are being transfered using another
        !blocking send to emulate more complicated mpi setup
        call mpi_recv(xres(my_tasks(irequest)), 1, mpi_integer, &
            my_workers(irequest), tag_indata, my_pool, my_status)
        !write(*,*) 'my_rank=',my_rank,'save to:', irequest, my_tasks(irequest), xres(my_tasks(irequest))

        !check if more work needs to be uploaded to the idling worker
        if (itask .le. ntasks) then
            irank = my_workers(irequest)
            my_tasks(irequest) = itask

            !initiate the handshake
            call mpi_send('func', 4, mpi_character, irank, tag_ftype, my_pool)

            !send data for processing
            call mpi_send(x(1,itask), 2, mpi_integer, irank, tag_indata, my_pool)

            !this is a second handshake
            !when we get a response it will mean
            !that we can pick up the results.
            call mpi_irecv(recvbuf, 1, mpi_integer, irank, tag_oudata,&
                           my_pool, my_requests(irequest))

            !up the increment
            itask = itask + 1
        endif

        !loop again over workers
        call mpi_waitany(nmintasks, my_requests, irequest, my_status)

        !save results that have been sent by the worker
        !xres(my_tasks(irequest)) = recvbuf
        !write(*,*) 'recvbuf=',recvbuf,'irequest=',irequest,'MPI_UNDEFINED=',MPI_UNDEFINED
    enddo
    !write(*,*) 'my_rank=',my_rank,'xres=',xres

    !post process results
    xsum = sum(xres)
    end subroutine workerfunc


    subroutine exitworkers()
    integer(kind=i4b) :: i, irank

    !tell all workers that we're done
    do i=1, size(my_workers)
        irank = my_workers(i)

        !tell worker to quit. Job's done.
        call mpi_send('exit', 4, mpi_character, irank, tag_ftype, my_pool)
    enddo
    end subroutine exitworkers
end program main