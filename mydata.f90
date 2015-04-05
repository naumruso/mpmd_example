module mydata
    use mpi_f08
    use iso_fortran_env
    implicit none

    integer, parameter :: sp  = REAL32  !< single precision real (float)
    integer, parameter :: dp  = REAL64  !< double precision real (double)
    integer, parameter :: i4b = INT32   !< 4 byte integer (int)
    
    type(mpi_comm) :: my_world !global communicator
    type(mpi_comm) :: my_pool !master <-> workers communicator
    type(mpi_comm) :: master_pool !masters communicator
    integer(kind=i4b) :: my_rank
    integer(kind=i4b) :: my_master 
    integer(kind=i4b), allocatable :: my_workers(:) !contains the ranks of the workers for the given master
    integer(kind=i4b) :: maxprocs
    integer(kind=i4b) :: nmasters
    integer(kind=i4b) :: nworkers
    integer(kind=i4b), parameter :: main_proc = 0 !this is the main process

    integer(kind=i4b), parameter :: TAG_ANY = 10
    integer(kind=i4b), parameter :: TAG_FTYPE = 100
    integer(kind=i4b), parameter :: TAG_INDATA= 101
    integer(kind=i4b), parameter :: TAG_OUDATA = 102

    contains
    subroutine init_world(ierror)
    integer(kind=i4b), intent(out) :: ierror

    integer(kind=i4b), parameter :: strlen = 1024
    integer(kind=i4b) :: valuelen
    character(len=strlen) :: value
    character(len=strlen) :: field
    integer(kind=i4b) :: num_apps
    integer(kind=i4b) :: nproc(2)
    logical :: flag

    ierror = 0

    !Start the MPI library
    call mpi_init()
    my_world = mpi_comm_world

    !get total number of processors
    field = 'maxprocs'
    call MPI_INFO_GET_VALUELEN(MPI_INFO_ENV, trim(field), valuelen, flag)
    call MPI_INFO_GET(MPI_INFO_ENV, trim(field), valuelen, VALUE(1:valuelen), FLAG)
    read(value(1:valuelen),*) maxprocs
    !write(*,*) 'maxprocs=',value(1:valuelen),flag

    !number of executables launched 
    field = 'ompi_num_apps'
    call MPI_INFO_GET_VALUELEN(MPI_INFO_ENV, trim(field), valuelen, flag)
    call MPI_INFO_GET(MPI_INFO_ENV, trim(field), valuelen, VALUE(1:valuelen), FLAG)
    !write(*,*) value(1:valuelen)
    read(value(1:valuelen),*) num_apps
    if (num_apps .ne. 2) then
        ierror = -1
        return
    endif

    !number of masters and workers
    field = 'ompi_np'
    call MPI_INFO_GET_VALUELEN(MPI_INFO_ENV, trim(field), valuelen, flag)
    call MPI_INFO_GET(MPI_INFO_ENV, trim(field), valuelen, VALUE(1:valuelen), FLAG)
    read(value(1:valuelen), *)  nproc
    nmasters = nproc(1)
    nworkers = nproc(2) !the total number of workers
    end subroutine init_world


    subroutine create_pools(master)
    logical, intent(in) :: master

    integer(kind=i4b) :: color
    integer(kind=i4b) :: irank, iworker
    integer(kind=i4b) :: new_rank

    type(mpi_status) :: my_status

    !get the rank in the global communicator
    call mpi_comm_rank(my_world, my_rank)

    !Now create the new master<->worker communicators
    color = mod(my_rank, nmasters)
    call mpi_comm_split(my_world, color, 0, my_pool)

    !Now masters have different ranks in my_pool, we need
    !to pass this information to the workers, that's all 
    !they need to know.
    !
    !The code below works when we start everything like this
    !mpirun -n NMASTERS master.out : -n MWORKERS worker.out
    !In this case the ranks look like this
    !|0,1,2,...,NMASTERS-1,NMASTERS,...,NMASTERS+NWORKERS-1|
    !|----MASTERS---------|--------WORKERS-----------------|
    !which means that we can use the function mod to separate
    !the processes into required groups
    if(master) then
        !find how many workers we have for this master
        iworker = 0
        do irank=nmasters, maxprocs - 1, 1
            color = mod(irank, nmasters)
            if (color .eq. my_rank) iworker = iworker + 1
        enddo
        !write(*,*) iworker
        allocate(my_workers(iworker)) !we keep their ranks here

        !Now load the array with their ranks in my_world comm.
        iworker = 1
        do irank=nmasters, maxprocs - 1, 1
            color = mod(irank, nmasters)
            if (color .eq. my_rank) then
                my_workers(iworker) = irank
                iworker = iworker + 1
            endif
        enddo
        !write(*,*) my_rank, my_workers

        !The rank of this master in my_pool
        call mpi_comm_rank(my_pool, new_rank)

        !Finally use sends to pass this information
        do iworker=1, size(my_workers), 1
            irank = my_workers(iworker)
            call mpi_send(new_rank, 1, mpi_integer, irank, &
                          tag_any, my_world)

        enddo

        !reload my_workers with the ranks of the workers
        !in the communicator my_pool. This is how we're going
        !to be communicating with them.
        iworker = 1
        do irank=0, size(my_workers), 1
            !write(*,*) iworker, irank, new_rank, my_rank
            if (irank .ne. new_rank) then
                my_workers(iworker) = irank
                iworker = iworker + 1
            endif
        enddo
        !write(*,*) my_rank, new_rank, my_workers
        !my_rank = new_rank
    else
        !find the rank of the master in my_world
        my_master = mod(my_rank, nmasters)

        !get its new rank in my_pool
        call mpi_recv(new_rank, 1, mpi_integer, my_master,&
                      tag_any, my_world, my_status)
        my_master = new_rank
    endif

    !Create the communicator for the masters
    !We might not need this communicator in the future. It
    !breaks the simple hierarchy that we're having now.
    if (master) then
        call mpi_comm_split(my_world, main_proc, 0, master_pool)

        call mpi_comm_rank(my_world, irank)
        if (main_proc .eq. irank) then
            call mpi_comm_rank(master_pool, new_rank)

            do irank=1, nmasters - 1, 1
                call mpi_send(new_rank, 1, mpi_integer, irank, &
                              tag_any, my_world)
            enddo

            my_master = new_rank
        else
            call mpi_recv(new_rank, 1, mpi_integer, main_proc,&
                          tag_any, my_world, my_status)
            my_master = new_rank
        endif
    else
        !For the workers the value of master_pool is MPI_COMM_NULL (opt out)
        call mpi_comm_split(my_world, MPI_UNDEFINED, 0, master_pool)
    endif

    call mpi_comm_rank(my_pool, new_rank)
    !write(*,'(A,I2,1x,A,I2)') 'oldrank=',my_rank,'new_rank=',new_rank
    end subroutine create_pools
end module mydata
