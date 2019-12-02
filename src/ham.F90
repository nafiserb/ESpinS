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

module mc_ham

  use mc_constants, only : dp
  use mc_parameters
  use mc_neighbors, only : dnn,multi,natlist,nsplist,ncell

  implicit none

  private



  public :: ham_get
  public :: ham_write
  public :: ham_dealloc

  real(kind=dp), allocatable :: coef_j(:,:,:)

contains 
  !==================================================================!
  subroutine ham_get()
    !==================================================================!
    !                                                                  !
    !                                                                  ! 
    !                                                                  !
    !===================================================================  
    use mc_io,        only : io_error,io_stopwatch
    use mc_neighbors

    implicit none

    ! Variables that are private
    
    integer       :: nsp,nat,ndnn,nnx,nsp1,nat1,ierr         
    


    call io_stopwatch('ham: get',1)

    allocate(coef_j(num_species,num_species,maxval(shells(1,:))), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating coef_j in ham_get')
    coef_j=0.0_dp

    call neighbors_get
    call neighbors_write
    call neighbors_check

    do nsp=1,num_species
       do nat=1,atoms_species_num(nsp)
          do ndnn = 1, shells(1,nsp)
             do nnx=1 , multi(nsp,nat,ndnn)
                nat1 =  natlist(nsp,nat,ndnn,nnx)
                nsp1 =  nsplist(nsp,nat,ndnn,nnx)
                coef_j(nsp,nsp1,ndnn) = coef_j(nsp,nsp1,ndnn) + & 
                                       mag_moments(nsp,nat)*mag_moments(nsp1,nat1)/&
                                    abs(mag_moments(nsp,nat)*mag_moments(nsp1,nat1))
             end do !nnx
          enddo
       enddo
    enddo
    ! We have double counting for each pair
    coef_j = coef_j*0.5_dp


    call io_stopwatch('ham: get',2)
    
    return
    
  end subroutine ham_get

  !==================================================================!
  subroutine ham_write()
    !==================================================================!
    !                                                                  !
    !                 Writes mcarlo.ham file                           !
    !                                                                  ! 
    !==================================================================! 
    use mc_io,        only : io_file_unit,seedname,io_date,io_stopwatch

    implicit none

    integer           :: nsp,nat,nsp1,ndnn,nnx,i
    integer           :: ham_unit
    logical           :: pair(num_species,num_species,maxval(shells(1,:)))
    character (len=9) :: cdate,ctime

    call io_stopwatch('ham: write',1)

    ham_unit=io_file_unit()
    open(unit=ham_unit,file=trim(seedname)//'.ham',form='formatted')

    ! Date and time
    call io_date(cdate,ctime)
    write(ham_unit,'(4(a),/)') 'File written on ',cdate,' at ',ctime

    write(ham_unit,*)
    write(ham_unit,'(36x,a6)') '------'
    write(ham_unit,'(36x,a6)') 'SYSTEM'
    write(ham_unit,'(36x,a6)') '------'
    write(ham_unit,*)
    if (lenconfac.eq.1.0_dp) then
       write(ham_unit,'(30x,a21)') 'Lattice Vectors (Ang)'
    else
       write(ham_unit,'(28x,a22)') 'Lattice Vectors (Bohr)'
    endif
    write(ham_unit,101) 'a_1',(real_lattice(1,I)*lenconfac, i=1,3)
    write(ham_unit,101) 'a_2',(real_lattice(2,I)*lenconfac, i=1,3)
    write(ham_unit,101) 'a_3',(real_lattice(3,I)*lenconfac, i=1,3)
    write(ham_unit,*)
    write(ham_unit,'(20x,a17,3x,f11.5)',advance='no') &
         'Unit Cell Volume:',cell_volume*lenconfac**3
    if (lenconfac.eq.1.0_dp) then
       write(ham_unit,'(2x,a7)') '(Ang^3)'
    else
       write(ham_unit,'(2x,a8)') '(Bohr^3)'
    endif
    write(ham_unit,*)
   ! Atoms
    if(num_atoms>0) then
       write(ham_unit,'(1x,a)') '*-----------------------------------------------------------------------------------*'
       if (lenconfac.eq.1.0_dp) then
          write(ham_unit,'(1x,a)') '| Atom  Type   Fractional Coordinate             Cartesian Coordinate(Ang)   Moment |'
       else
          write(ham_unit,'(1x,a)') '| Atom  Type   Fractional Coordinate             Cartesian Coordinate(Bohr)  Moment |'
       endif
       write(ham_unit,'(1x,a)') '+-----------------------------------------------------------------------------------+'
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             write(ham_unit,'(1x,a1,1x,a3,2x,i3,3F10.5,3x,a1,3F10.5,2x,a1,F6.3,1x,a1)') '|',atoms_symbol(nsp),nsp,&
                   atoms_pos_frac(:,nsp,nat),'|',atoms_pos_cart(:,nsp,nat)*lenconfac,'|',mag_moments(nsp,nat),'|'
          end do
       end do
       write(ham_unit,'(1x,a)') '*-----------------------------------------------------------------------------------*'
    else
       write(ham_unit,'(25x,a)') 'No atom positions specified'
    end if

    write(ham_unit,*) ' '



    ! Neighbours
    if(num_atoms>0) then
    write(ham_unit,'(1x,a)') '*-----------------------------------------------------------------------------------*'
    write(ham_unit,'(1x,a)') '|                       Distance to Nearest-Neighbour Shells                        |'
    write(ham_unit,'(1x,a)') '|                       ------------------------------------                        |'
    if (lenconfac.eq.1.0_dp) then
       write(ham_unit,'(1x,a)') '|              Shell               Distance (Ang)            Multiplicity           |'
       write(ham_unit,'(1x,a)') '|              -----               --------------            ------------           |'
    else
       write(ham_unit,'(1x,a)') '|              Shell               Distance (Bohr)           Multiplicity           |'
       write(ham_unit,'(1x,a)') '|              -----               ---------------           ------------           |'
    endif

    do nsp=1,num_species
       do nat=1,atoms_species_num(nsp)
          do ndnn = 1, shells(1,nsp)
             write(ham_unit,'(1x,a,14x,i3,17x,f10.6,19x,i4,16x,a)') '|',ndnn,dnn(nsp,nat,ndnn)*lenconfac,&
                   multi(nsp,nat,ndnn),'|'
          enddo
       enddo
    enddo
    write(ham_unit,'(1x,a)') '*-----------------------------------------------------------------------------------*'
    endif

    write(ham_unit,'(1x,a)') '*-----------------------------------------------------------------------------------*'
    write(ham_unit,'(1x,a)') '|   Neighbors      Atom1     Type1     Atom2     Type2          coefficient         |'
    write(ham_unit,'(1x,a)') '+-----------------------------------------------------------------------------------+'

    pair = .False.
    do nsp=1,num_species
       do nat=1,atoms_species_num(nsp)
          do ndnn = 1, shells(1,nsp)
             do nnx=1 , multi(nsp,nat,ndnn)
                nsp1 = nsplist(nsp,nat,ndnn,nnx)
                pair(nsp,nsp1,ndnn) = .True.
             enddo
          enddo
       enddo
    enddo

    do nsp=1,num_species
       do ndnn = 1, shells(1,nsp)
          do nsp1=1,num_species
             if( pair(nsp,nsp1,ndnn)) &
                   write(ham_unit,'(1x,a1,6x,a2,i1,10x,a3,6x,i3,8x,a3,6x,i3,14x,F6.2,12x,a1)') '|','J_',ndnn,&
                         atoms_symbol(nsp),nsp,atoms_symbol(nsp1),nsp1,coef_j(nsp,nsp1,ndnn),'|'
          enddo
       enddo
    enddo


    write(ham_unit,'(1x,a)') '*-----------------------------------------------------------------------------------*'

    close(ham_unit)

    call io_stopwatch('ham: write',2)

101 format(20x,a3,2x,3F11.6)


    return

  end subroutine ham_write

  !==================================================================!
  subroutine ham_dealloc()
    !==================================================================!
    !                                                                  !
    !                                                                  ! 
    !                                                                  !
    !===================================================================  
    use mc_io,   only : io_error
    implicit none
    integer :: ierr


    ! Deallocate integer arrays that are public

    deallocate(coef_j, stat=ierr )
    if (ierr/=0) call io_error('Error in deallocating coef_j in ham_dealloc')

    return

  end subroutine ham_dealloc


end module mc_ham



