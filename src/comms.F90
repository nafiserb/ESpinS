!######################################################################
! This routine is part of
! ESpinS - Esfahan Spin Simulation 
! (c) 2016-2019 Dr. Nafise Rezaei and Dr. Mojtaba Alaei
! Physics Department, Isfahan University of Technology, Isfahan, Iran
!
! This program is free software: you can redistribute it and/or modify it 
! under the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or 
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
! or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
! for more details.
!
! You should have received a copy of the GNU General Public License along 
! with this program. If not, see http://www.gnu.org/licenses. 
!######################################################################


module mc_comms

  use mc_constants, only : dp
  use mc_io, only: io_error

  implicit none


  private

#ifdef MPI
  include 'mpif.h'
#endif

  logical, public, save      :: on_root
  integer, public, save      :: num_nodes,my_node_id
  integer, public, parameter :: root_id = 0

  integer,         parameter :: mpi_send_tag = 77 !abitrary

  public :: comms_setup
  public :: comms_end
!  public :: comms_abort     ! [GP]: do not use, use io_error instead
  public :: comms_bcast      ! send data from the root node
  public :: comms_send       ! send data from one node to another
  public :: comms_recv       ! accept data from one node to another
  public :: comms_reduce     ! reduce data onto root node (n.b. not allreduce);
                             ! note that on all other nodes, the data is lost
!  public :: comms_allreduce  ! reduce data onto all nodes
  public :: comms_barrier    ! puts a barrier so that the code goes on only when all nodes reach the barrier
  public :: comms_gatherv    ! gets chunks of an array from all nodes and gathers them on the root node
  public :: comms_scatterv    ! sends chunks of an array to all nodes scattering them from the root node

  public :: comms_array_split
  public :: comms_init_random

  interface comms_bcast
     module procedure comms_bcast_int
     module procedure comms_bcast_logical
     module procedure comms_bcast_real
!     module procedure comms_bcast_cmplx
     module procedure comms_bcast_char
  end interface comms_bcast

  interface comms_send
     module procedure comms_send_int
!     module procedure comms_send_logical
     module procedure comms_send_real
!     module procedure comms_send_cmplx
!     module procedure comms_send_char
  end interface comms_send

  interface comms_recv
     module procedure comms_recv_int
!     module procedure comms_recv_logical
     module procedure comms_recv_real
!     module procedure comms_recv_cmplx
!     module procedure comms_recv_char
  end interface comms_recv
  
  interface comms_reduce
     module procedure comms_reduce_int 
     module procedure comms_reduce_real
!     module procedure comms_reduce_cmplx
  end interface comms_reduce
!
!  interface comms_allreduce
!     module procedure comms_allreduce_int    ! to be done
!     module procedure comms_allreduce_real
!     module procedure comms_allreduce_cmplx
!  end interface comms_allreduce

  interface comms_gatherv
     module procedure comms_gatherv_int    ! to be done
     module procedure comms_gatherv_real
!     module procedure comms_gatherv_cmplx
  end interface comms_gatherv

  interface comms_scatterv
     module procedure comms_scatterv_int  
     module procedure comms_scatterv_real
!     module procedure comms_scatterv_cmplx
  end interface comms_scatterv

contains

  !==================================================================!
  subroutine comms_setup
    !==================================================================!
    !                                                                  !
    !==================================================================!
 
    implicit none

#ifdef MPI
    integer :: ierr

    call mpi_init(ierr)
    if (ierr.ne.0) call io_error('MPI initialisation error')
    call mpi_comm_rank(mpi_comm_world, my_node_id, ierr)
    call mpi_comm_size(mpi_comm_world, num_nodes, ierr)
#else
    num_nodes = 1
    my_node_id = 0
#endif

    on_root = .false.
    if(my_node_id==root_id) on_root = .true.
    
  end subroutine comms_setup

  !==================================================================!
  !> Given an array of size numpoints, we want to split on num_nodes nodes. This function returns
  !> two arrays: count and displs.
  !> The i-th element of the count array gives the number of elements
  !> that must be calculated by the process with id (i-1).
  !> The i-th element of the displs array gives the displacement of the array calculated locally on
  !> the process with id (i-1) with respect to the global array.
  !>
  !> \note These values are those to be passed to the functions MPI_Scatterv, MPI_Gatherv and MPI_Alltoallv.
  !>
  !> \note one can use the following do loop to run over the needed elements, if the full array is stored
  !> on all nodes:
  !> do i=displs(my_node_id)+1,displs(my_node_id)+counts(my_node_id)
  !> 
  !> \param numpoints Number of elements of the array to be scattered
  !> \param counts    Array (of size num_nodes) with the number of elements of the array on each node
  !> \param displs    Array (of size num_nodes) with the displacement relative to the global array

  !==================================================================!
  subroutine comms_array_split(numpoints,counts,displs)
    !==================================================================!
    !                                                                  !
    !==================================================================!  
    integer, intent(in) :: numpoints
    integer, dimension(0:num_nodes-1), intent(out) :: counts
    integer, dimension(0:num_nodes-1), intent(out) :: displs

    integer :: ratio, remainder, i

    ratio = numpoints / num_nodes
    remainder = MOD(numpoints, num_nodes)

    do i=0,num_nodes-1
       if (i < remainder) then
          counts(i) = ratio+1
          displs(i) = i*(ratio+1)
       else
          counts(i) = ratio
          displs(i) = remainder*(ratio+1) + (i-remainder)*ratio
       end if
    end do

  end subroutine comms_array_split

  !==================================================================!
  subroutine comms_end
    !==================================================================!
    !                                                                  !
    !==================================================================!  
 
    implicit none

#ifdef MPI
    integer :: ierr

    call mpi_finalize(ierr)
#endif
    
  end subroutine comms_end

!==================================================================!
  subroutine comms_barrier
    !==================================================================!
    !                                                                  !
    !==================================================================!  
 
    implicit none

#ifdef MPI
    integer :: ierr

    call mpi_barrier(mpi_comm_world, ierr)
#endif
    
  end subroutine comms_barrier

!==================================================================!
!  subroutine comms_abort
!
!    implicit none
!
!    integer :: ierr
!
!#ifdef MPI
!    call MPI_abort(MPI_comm_world,1,ierr)
!#else
!    STOP
!#endif
!
!  end subroutine comms_abort

    !============================Bcast=================================!
    !                                                                  !
    !==================================================================!  

  subroutine comms_bcast_int(array,size)

    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array,size,MPI_integer,root_id,mpi_comm_world,error)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_bcast_int')
    end if
#endif

    return

  end subroutine comms_bcast_int

  subroutine comms_bcast_real(array,size)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer,       intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array,size,MPI_double_precision,root_id,mpi_comm_world,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_bcast_real')
    end if
#endif

    return

  end subroutine comms_bcast_real

  subroutine comms_bcast_logical(array,size)

    implicit none

    logical, intent(inout) :: array
    integer, intent(in)    :: size

#ifdef MPI
    integer :: error

    call MPI_bcast(array,size,MPI_logical,root_id,mpi_comm_world,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_bcast_logical')
    end if
#endif

    return

  end subroutine comms_bcast_logical

  subroutine comms_bcast_char(array,size)

    implicit none

    character(len=*), intent(inout) :: array
    integer,          intent(in)    :: size


#ifdef MPI
    integer :: error

    call MPI_bcast(array,size,MPI_character,root_id,mpi_comm_world,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_bcast_char')
    end if
#endif

    return

  end subroutine comms_bcast_char

!  subroutine comms_bcast_cmplx(array,size)
!
!    implicit none
!
!    complex(kind=dp), intent(inout) :: array
!    integer, intent(in)    :: size
!
!
!#ifdef MPI
!    integer :: error
!
!    call MPI_bcast(array,size,MPI_double_complex,root_id,mpi_comm_world,error)
!
!    if(error.ne.MPI_success) then
!       call io_error('Error in comms_bcast_cmplx')
!    end if
!#endif
!
!    return
!
!  end subroutine comms_bcast_cmplx

    !============================Send==================================!
    !                                                                  !
    !==================================================================!  

!  subroutine comms_send_logical(array,size,to)
!
!    implicit none
!
!    logical, intent(inout) :: array
!    integer, intent(in)    :: size
!    integer, intent(in)    :: to
!
!#ifdef MPI
!    integer :: error
!
!    call MPI_send(array,size,MPI_logical,to, &
!         mpi_send_tag,mpi_comm_world,error)
!
!    if(error.ne.MPI_success) then
!       call io_error('Error in comms_send_logical')
!    end if
!#endif
!
!    return
!
!  end subroutine comms_send_logical


  subroutine comms_send_int(array,size,to)

    implicit none

    integer(8), intent(inout) :: array
    integer,    intent(in)    :: size
    integer,    intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array,size,MPI_Integer8,to, &
         mpi_send_tag,mpi_comm_world,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_send_int')
    end if
#endif

    return

  end subroutine comms_send_int


!  subroutine comms_send_char(array,size,to)
!
!    implicit none
!
!    character(len=*), intent(inout) :: array
!    integer, intent(in)    :: size
!    integer, intent(in)    :: to
!
!#ifdef MPI
!    integer :: error
!
!    call MPI_send(array,size,MPI_character,to, &
!         mpi_send_tag,mpi_comm_world,error)
!
!    if(error.ne.MPI_success) then
!       call io_error('Error in comms_send_char')
!    end if
!#endif
!
!    return
!
!  end subroutine comms_send_char


  subroutine comms_send_real(array,size,to)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer,       intent(in)    :: size
    integer,       intent(in)    :: to

#ifdef MPI
    integer :: error

    call MPI_send(array,size,MPI_double_precision,to, &
         mpi_send_tag,mpi_comm_world,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_send_real')
    end if
#endif

    return

  end subroutine comms_send_real


!  subroutine comms_send_cmplx(array,size,to)
!
!    implicit none
!
!    complex(kind=dp), intent(inout) :: array
!    integer, intent(in)    :: size
!    integer, intent(in)    :: to
!
!
!#ifdef MPI
!    integer :: error
!
!    call MPI_send(array,size,MPI_double_complex,to, &
!         mpi_send_tag,mpi_comm_world,error)
!
!    if(error.ne.MPI_success) then
!       call io_error('Error in comms_send_cmplx')
!    end if
!#endif
!
!    return
!
!  end subroutine comms_send_cmplx


    !============================Recive================================!
    !                                                                  !
    !==================================================================!  
!
!  subroutine comms_recv_logical(array,size,from)
!
!    implicit none
!
!    logical, intent(inout) :: array
!    integer, intent(in)    :: size
!    integer, intent(in)    :: from
!
!#ifdef MPI
!    integer :: error
!    integer :: status(MPI_status_size)
!
!    call MPI_recv(array,size,MPI_logical,from, &
!         mpi_send_tag,mpi_comm_world,status,error)
!
!    if(error.ne.MPI_success) then
!       call io_error('Error in comms_recv_logical')
!    end if
!#endif
!
!    return
!
!  end subroutine comms_recv_logical


  subroutine comms_recv_int(array,size,from)

    implicit none

    integer(8), intent(inout) :: array
    integer,    intent(in)    :: size
    integer,    intent(in)    :: from

#ifdef MPI
    integer :: error
    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_Integer8,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_recv_int')
    end if
#endif

    return

  end subroutine comms_recv_int


!  subroutine comms_recv_char(array,size,from)
!
!    implicit none
!
!    character(len=*), intent(inout) :: array
!    integer, intent(in)    :: size
!    integer, intent(in)    :: from
!
!#ifdef MPI
!    integer :: error
!    integer :: status(MPI_status_size)
!
!    call MPI_recv(array,size,MPI_character,from, &
!         mpi_send_tag,mpi_comm_world,status,error)
!
!    if(error.ne.MPI_success) then
!       call io_error('Error in comms_recv_char')
!    end if
!#endif
!
!    return
!
!  end subroutine comms_recv_char


  subroutine comms_recv_real(array,size,from)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer,       intent(in)    :: size
    integer,       intent(in)    :: from

#ifdef MPI
    integer :: error
    integer :: status(MPI_status_size)

    call MPI_recv(array,size,MPI_double_precision,from, &
         mpi_send_tag,mpi_comm_world,status,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_recv_real')
    end if
#endif

    return

  end subroutine comms_recv_real


!  subroutine comms_recv_cmplx(array,size,from)
!
!    implicit none
!
!    complex(kind=dp), intent(inout) :: array
!    integer, intent(in)    :: size
!    integer, intent(in)    :: from
!
!#ifdef MPI
!    integer :: error
!
!    integer :: status(MPI_status_size)
!
!    call MPI_recv(array,size,MPI_double_complex,from, &
!         mpi_send_tag,mpi_comm_world,status,error)
!
!    if(error.ne.MPI_success) then
!       call io_error('Error in comms_recv_cmplx')
!    end if
!
!#endif
!
!    return
!
!  end subroutine comms_recv_cmplx


    !============================Error=================================!
    !                                                                  !
    !==================================================================!  
!  subroutine comms_error
!
!    implicit none
!       
!#ifdef MPI
!    integer :: error
!
!    call MPI_abort(MPI_comm_world,1,error)
!
!#endif
!
!  end subroutine comms_error

  
    !============================Reduce================================!
    !                                                                  !
    !==================================================================!  
  ! COMMS_REDUCE (collect data on the root node)

  subroutine comms_reduce_int(array,size,op)

    implicit none

    integer, intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error,ierr

    integer, allocatable :: array_red(:)

    allocate(array_red(size),stat=ierr)
    if (ierr/=0) then
       call io_error('failure to allocate array_red in comms_reduce_int')
    end if

    select case(op)

    case ('SUM')
       call MPI_reduce(array,array_red,size,MPI_integer,MPI_sum,root_id,mpi_comm_world,error)
    case ('PRD')
       call MPI_reduce(array,array_red,size,MPI_integer,MPI_prod,root_id,mpi_comm_world,error)
    case default
       call io_error('Unknown operation in comms_reduce_int')

    end select

    call my_icopy(size,array_red,1,array,1)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_reduce_int')
    end if

    if (allocated(array_red)) deallocate(array_red)
#endif

    return

  end subroutine comms_reduce_int


  subroutine comms_reduce_real(array,size,op)

    implicit none

    real(kind=dp), intent(inout) :: array
    integer, intent(in)    :: size
    character(len=*), intent(in) :: op

#ifdef MPI
    integer :: error,ierr

    real(kind=dp), allocatable :: array_red(:)

    allocate(array_red(size),stat=ierr)
    if (ierr/=0) then
       call io_error('failure to allocate array_red in comms_reduce_real')
    end if

    select case(op)

    case ('SUM')
       call MPI_reduce(array,array_red,size,MPI_double_precision,MPI_sum,root_id,mpi_comm_world,error)
    case ('PRD')
       call MPI_reduce(array,array_red,size,MPI_double_precision,MPI_prod,root_id,mpi_comm_world,error)
    case ('MIN')
       call MPI_reduce(array,array_red,size,MPI_double_precision,MPI_MIN,root_id,mpi_comm_world,error)
    case ('MAX')
       call MPI_reduce(array,array_red,size,MPI_double_precision,MPI_max,root_id,mpi_comm_world,error)
    case default
       call io_error('Unknown operation in comms_reduce_real')

    end select

    call dcopy(size,array_red,1,array,1)

    if(error.ne.MPI_success) then
       call io_error('Error in comms_reduce_real')
    end if

    if (allocated(array_red)) deallocate(array_red)
#endif

    return

  end subroutine comms_reduce_real


!  subroutine comms_reduce_cmplx(array,size,op)
!
!    implicit none
!
!    complex(kind=dp), intent(inout) :: array
!    integer, intent(in)    :: size
!    character(len=*), intent(in) :: op
!
!#ifdef MPI
!    integer :: error,ierr
!
!    complex(kind=dp), allocatable :: array_red(:)
!
!    allocate(array_red(size),stat=ierr)
!    if (ierr/=0) then
!       call io_error('failure to allocate array_red in comms_reduce_cmplx')
!    end if
!
!    select case(op)
!
!    case ('SUM')
!       call MPI_reduce(array,array_red,size,MPI_double_complex,MPI_sum,root_id,mpi_comm_world,error)
!    case ('PRD')
!       call MPI_reduce(array,array_red,size,MPI_double_complex,MPI_prod,root_id,mpi_comm_world,error)
!    case default
!       call io_error('Unknown operation in comms_reduce_cmplx')
!
!    end select
!
!    call zcopy(size,array_red,1,array,1)
!
!    if(error.ne.MPI_success) then
!       call io_error('Error in comms_reduce_cmplx')
!    end if
!
!    if (allocated(array_red)) deallocate(array_red)
!#endif
!
!    return
!
!  end subroutine comms_reduce_cmplx
!
!  subroutine comms_allreduce_real(array,size,op)
!
!    implicit none
!
!    real(kind=dp), intent(inout) :: array
!    integer, intent(in)    :: size
!    character(len=*), intent(in) :: op
!
!#ifdef MPI
!    integer :: error,ierr
!
!    real(kind=dp), allocatable :: array_red(:)
!
!    allocate(array_red(size),stat=ierr)
!    if (ierr/=0) then
!       call io_error('failure to allocate array_red in comms_allreduce_real')
!    end if
!
!    select case(op)
!
!    case ('SUM')
!       call MPI_allreduce(array,array_red,size,MPI_double_precision,MPI_sum,mpi_comm_world,error)
!    case ('PRD')
!       call MPI_allreduce(array,array_red,size,MPI_double_precision,MPI_prod,mpi_comm_world,error)
!    case ('MIN')
!       call MPI_allreduce(array,array_red,size,MPI_double_precision,MPI_MIN,mpi_comm_world,error)
!    case ('MAX')
!       call MPI_allreduce(array,array_red,size,MPI_double_precision,MPI_max,mpi_comm_world,error)
!    case default
!       call io_error('Unknown operation in comms_allreduce_real')
!
!    end select
!
!    call dcopy(size,array_red,1,array,1)
!
!    if(error.ne.MPI_success) then
!       call io_error('Error in comms_allreduce_real')
!    end if
!
!    if (allocated(array_red)) deallocate(array_red)
!#endif
!
!    return
!
!  end subroutine comms_allreduce_real
!
!  subroutine comms_allreduce_cmplx(array,size,op)
!
!    implicit none
!
!    complex(kind=dp), intent(inout) :: array
!    integer, intent(in)    :: size
!    character(len=*), intent(in) :: op
!
!#ifdef MPI
!    integer :: error,ierr
!
!    complex(kind=dp), allocatable :: array_red(:)
!
!    allocate(array_red(size),stat=ierr)
!    if (ierr/=0) then
!       call io_error('failure to allocate array_red in comms_allreduce_cmplx')
!    end if
!
!    select case(op)
!
!    case ('SUM')
!       call MPI_allreduce(array,array_red,size,MPI_double_complex,MPI_sum,mpi_comm_world,error)
!    case ('PRD')
!       call MPI_allreduce(array,array_red,size,MPI_double_complex,MPI_prod,mpi_comm_world,error)
!    case default
!       call io_error('Unknown operation in comms_allreduce_cmplx')
!
!    end select
!
!    call zcopy(size,array_red,1,array,1)
!
!    if(error.ne.MPI_success) then
!       call io_error('Error in comms_allreduce_cmplx')
!    end if
!
!    if (allocated(array_red)) deallocate(array_red)
!#endif
!
!    return
!
!  end subroutine comms_allreduce_cmplx

  !============================Gather================================!
  !                                                                  !
  !==================================================================!

  ! Array: local array for sending data; localcount elements will be sent
  !        to the root node
  ! rootglobalarray: array on the root node to which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split
  subroutine comms_gatherv_int(array,localcount,rootglobalarray,counts,displs)

    implicit none

    integer(8),                 intent(inout) :: array
    integer,                       intent(in) :: localcount
    integer(8),                 intent(inout) :: rootglobalarray
    integer, dimension(num_nodes), intent(in) :: counts
    integer, dimension(num_nodes), intent(in) :: displs

#ifdef MPI
    integer :: error

    call MPI_gatherv(array,localcount,MPI_Integer8,rootglobalarray,counts,&
         displs,MPI_Integer8,root_id,mpi_comm_world,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_gatherv_int')
    end if

#else
    call my_icopy(localcount,array,1,rootglobalarray,1)
#endif

    return

  end subroutine comms_gatherv_int


  ! Array: local array for sending data; localcount elements will be sent
  !        to the root node
  ! rootglobalarray: array on the root node to which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split
  subroutine comms_gatherv_real(array,localcount,rootglobalarray,counts,displs)

    implicit none

    real(kind=dp),              intent(inout) :: array
    integer,                       intent(in) :: localcount
    real(kind=dp),              intent(inout) :: rootglobalarray
    integer, dimension(num_nodes), intent(in) :: counts
    integer, dimension(num_nodes), intent(in) :: displs

#ifdef MPI
    integer :: error

    call MPI_gatherv(array,localcount,MPI_double_precision,rootglobalarray,counts,&
         displs,MPI_double_precision,root_id,mpi_comm_world,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_gatherv_real')
    end if

#else
    call dcopy(localcount,array,1,rootglobalarray,1)
#endif

    return

  end subroutine comms_gatherv_real

  !============================Scatter===============================!
  !                                                                  !
  !==================================================================!
  ! Array: local array for getting data; localcount elements will be fetched
  !        from the root node
  ! rootglobalarray: array on the root node from which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split
  subroutine comms_scatterv_int(array,localcount,rootglobalarray,counts,displs)

    implicit none

    integer,                    intent(inout) :: array
    integer,                    intent(in)    :: localcount
    integer,                    intent(inout) :: rootglobalarray
    integer, dimension(num_nodes), intent(in) :: counts
    integer, dimension(num_nodes), intent(in) :: displs

#ifdef MPI
    integer :: error

    call MPI_scatterv(rootglobalarray,counts,displs,MPI_Integer,&
         Array,localcount,MPI_Integer,root_id,mpi_comm_world,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_scatterv_real')
    end if

#else
    call my_icopy(localcount,rootglobalarray,1,array,1)
#endif

    return

  end subroutine comms_scatterv_int



  ! Array: local array for getting data; localcount elements will be fetched
  !        from the root node
  ! rootglobalarray: array on the root node from which data will be sent
  ! counts, displs : how data should be partitioned, see MPI documentation or
  !                  function comms_array_split
  subroutine comms_scatterv_real(array,localcount,rootglobalarray,counts,displs)

    implicit none

    real(kind=dp),              intent(inout) :: array
    integer,                       intent(in) :: localcount
    real(kind=dp),              intent(inout) :: rootglobalarray
    integer, dimension(num_nodes), intent(in) :: counts
    integer, dimension(num_nodes), intent(in) :: displs

#ifdef MPI
    integer :: error

!    call MPI_scatterv(array,localcount,MPI_double_precision,rootglobalarray,counts,&
!         displs,MPI_double_precision,root_id,mpi_comm_world,error)
    call MPI_scatterv(rootglobalarray,counts,displs,MPI_double_precision,&
         array,localcount,MPI_double_precision,root_id,mpi_comm_world,error)

    if (error.ne.MPI_success) then
       call io_error('Error in comms_scatterv_real')
    end if

#else
    call dcopy(localcount,rootglobalarray,1,array,1)
#endif

    return

  end subroutine comms_scatterv_real

  !==================================================================!
  subroutine comms_init_random(lseed,seeds,state)
    !==================================================================!
    !                                                                  !
    !                                                                  !
    !                                                                  !
    !===================================================================  

    use stdtypes
    use mtprng,        only : mtprng_init,mtprng_state

    implicit none
    logical,             intent(in)    :: lseed
    integer,             intent(inout) :: seeds(num_nodes)
    type(mtprng_state),  intent(inout) :: state

    if (.not.lseed) then
      open(89,FILE='/dev/urandom',ACCESS='stream',FORM='UNFORMATTED')
      read(89) seeds(my_node_id+1)
      close(89)
    endif

    call mtprng_init(seeds(my_node_id+1),state)

    ! gather random seeds on to root_node
    if(.not.lseed) call comms_reduce(seeds(1), num_nodes, 'SUM')

    return

  end subroutine comms_init_random

end module mc_comms

subroutine my_ICOPY(N,ZX,INCX,ZY,INCY)
  !     .. Scalar Arguments ..
  integer INCX,INCY,N
  !     ..
  !     .. Array Arguments ..
  integer ZX(*),ZY(*)
  !     ..
  !
  !  Purpose
  !  =======
  !
  !     copies a vector, x, to a vector, y.
  !     jack dongarra, linpack, 4/11/78.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  !
  !     .. Local Scalars ..
  integer I,IX,IY
  !     ..
  if (N.le.0) return
  if (INCX.eq.1 .and. INCY.eq.1) GO TO 20
  !
  !        code for unequal increments or equal increments
  !          not equal to 1
  !
  IX = 1
  IY = 1
  if (INCX.lt.0) IX = (-N+1)*INCX + 1
  if (INCY.lt.0) IY = (-N+1)*INCY + 1
  do I = 1,N
     ZY(IY) = ZX(IX)
     IX = IX + INCX
     IY = IY + INCY
  end do
  return
  !
  !        code for both increments equal to 1
  !
20 do I = 1,N
     ZY(I) = ZX(I)
  end do
  return
end subroutine my_ICOPY


