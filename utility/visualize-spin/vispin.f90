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

program visualize_spin

  implicit none

  integer, parameter            :: dp=selected_real_kind(15,300)
  real(dp), parameter           :: pi=3.141592653589793238462643383279_dp
  real(dp)                      :: t1,t2
  integer                       :: natoms_cell,num_tems
  integer                       :: num_supercell,supercell_size(3),ssize(3)
  real(dp)                      :: lattice(3,3),pos(3)
  integer                       :: i,j,n,n1,n2,n3,num_cell,num_atom 
  character(len=3), allocatable :: atoms_label(:)
  real(kind=dp), allocatable    :: atoms_pos(:,:),spin(:,:,:),tems(:)     
  integer                       :: sconfigunit,xsfunit
  integer, allocatable          :: visunit(:)                      
  character(len=30)             :: filename,filename1,dummy
  integer                       :: iostat

  call cpu_time(t1)
  call get_command_argument(1,filename)
  sconfigunit = io_file_unit() 
  open(unit=sconfigunit,file=adjustl(trim(filename))//"_sconfig.dat",action='read',iostat=iostat)
  if (iostat /= 0) then
     write(0,'(a)') "Error: Counldn't find '"// trim(filename)//"_sconfig.dat' file."
     stop 
  end if

  xsfunit = io_file_unit() 
  open(unit=xsfunit,file=adjustl(trim(filename))//".xsf",action='read',iostat=iostat)
  if (iostat /= 0) then
     write(0,'(a)') "Error: Counldn't find '" // trim(filename)//".xsf' file."
     stop 
  end if


  
  ! read xsf file        
  read(xsfunit,*) dummy ! 
  read(xsfunit,*) dummy ! 
  read(xsfunit,*) lattice(1,1), lattice(1,2),lattice(1,3)   ! A1                
  read(xsfunit,*) lattice(2,1), lattice(2,2),lattice(2,3)   ! A2                
  read(xsfunit,*) lattice(3,1), lattice(3,2),lattice(3,3)   ! A3                
  do i=1,5
     read(xsfunit,*) dummy ! 
  enddo
  read(xsfunit,*) natoms_cell, n ! number of atoms in unit cell and number pof types                

  allocate(atoms_label(natoms_cell))
  allocate(atoms_pos(natoms_cell,3))
  do i=1,natoms_cell
     read(xsfunit,*) atoms_label(i),atoms_pos(i,1),atoms_pos(i,2),atoms_pos(i,3)
  enddo
  close(xsfunit)

  read(sconfigunit,*) dummy ! 
  read(sconfigunit,*) supercell_size(1),supercell_size(2),supercell_size(3)
  read(sconfigunit,*) num_supercell,num_tems
  num_supercell=num_supercell/natoms_cell 
  if(supercell_size(1)*supercell_size(2)*supercell_size(3) .ne. num_supercell) then
    write(0,'(a)') "ERROR: inconsistancy between supercell_size and num_total_atoms in "//trim(filename)//"_sconfig.dat file!"
    stop 
  endif

  allocate(spin(num_tems,num_supercell,3))
  allocate(tems(num_tems))
  do i=1,num_tems
     read(sconfigunit,*) tems(i) 
     do j=1,num_supercell
        read(sconfigunit,*) spin(i,j,1),spin(i,j,2),spin(i,j,3)
     enddo
  enddo

  write(*,*) 'Please eneter the supercell size for visualization (for example: 2 2 2):'
  read (*,*) ssize(1),ssize(2),ssize(3)
  if ((ssize(1) .gt. supercell_size(1)) .or. (ssize(2) .gt. supercell_size(2)) &
      .or. (ssize(3) .gt. supercell_size(3))) &
     write(0,'(a)') "ERROR: supercell size for visualization should be less or equal to " &
                     // "supercell_size in "//trim(filename)//"_sconfig.dat file!"


  allocate(visunit(num_tems))
  do n=1,num_tems
     n1 = n/100
     n2 = (n-n1*100)/10
     n3 = n-n1*100-n2*10
     visunit(n) = io_file_unit()
     filename1 = trim(filename)//'_vispin'&
                //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat'
     open(visunit(n),file=filename1,status='unknown',form='formatted',action='write')
     write(visunit(n),'(a)') '# Running vispin. Written by N Rezaei. (c) 2019.'
     write(visunit(n),'(a,1f12.6)') '# T = ', tems(n)                                  
     write(visunit(n),'(a)') 'ATOMS'
  enddo

  do j=1,num_tems
     do n3=0,ssize(3)-1
        do n2=0,ssize(2)-1
           do n1=0,ssize(1)-1 
              pos = 0.0_dp
              num_cell = n3*supercell_size(2)*supercell_size(1)+n2*supercell_size(1)+n1
              do i=1,natoms_cell
                    num_atom=num_cell*natoms_cell+i
                    pos(:)=atoms_pos(i,:)+n1*lattice(1,:)+n2*lattice(2,:)+n3*lattice(3,:)
                    write(visunit(j),'(1a,6E16.8)') atoms_label(i), pos(1),pos(2),pos(3),&
                                                 spin(j,num_atom,1),spin(j,num_atom,2),spin(j,num_atom,3)  
              enddo
           enddo
        enddo
     enddo
  enddo
 

  deallocate(atoms_label)
  deallocate(atoms_pos)
  deallocate(spin)
  deallocate(tems) 
  close(sconfigunit)
  do n=1,num_tems
     close(visunit(n))
  enddo
  deallocate(visunit)

  call cpu_time(t2)
  print*,'time :',t2-t1,'(s)'

contains

  !==================================================================!
  function io_file_unit()
  !==================================================================!
  !                                                                  !
  ! Returns an unused unit number                                    !
  ! (so we can open a file on that unit                              !
  !                                                                  !
  !=================================================================== 
  implicit none

  integer :: io_file_unit,unit
  logical :: file_open

  unit = 9
  file_open = .true.
  do while ( file_open )
     unit = unit + 1
     inquire( unit, OPENED = file_open )
  end do

  io_file_unit = unit


  return
end function io_file_unit
 
end program visualize_spin
