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

module mc_dij      

  use mc_constants,  only : dp
  use mc_parameters, only : num_total_atoms,num_species,dij_num_nbors
  use mc_io,         only : io_error

  implicit none

  private


  real(kind=dp), allocatable, public, save :: dij_matrix(:,:,:)
  integer      , allocatable, public, save :: dij_nbors_matrix(:,:)
  integer      , allocatable, public, save :: dij_num_nbors_matrix(:)

  public :: dij_get
  public :: dij_setup
  public :: dij_dealloc


  ! Private Data
  real(kind=dp) , allocatable :: dij_tot_matrix(:,:,:)
  integer       , allocatable :: dij_tot_nbors_matrix(:,:)
  integer                     :: ierr

contains 

  !============================================!
  subroutine dij_setup()
    !============================================!

    implicit none

    allocate(dij_num_nbors_matrix(num_total_atoms), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating dij_num_nbors_matrix in dij_setup')

    return
  end subroutine dij_setup

  !============================================!
  subroutine dij_dealloc()
    !============================================!

    implicit none
    
    ! Deallocate arrays that are public 

    if ( allocated( dij_matrix ) ) then
       deallocate( dij_matrix, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating dij_matrix in param_dealloc')
    end if
    if ( allocated( dij_nbors_matrix ) ) then
       deallocate( dij_nbors_matrix, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating dij_nbors_matrix in param_dealloc')
    end if
    if ( allocated( dij_num_nbors_matrix ) ) then
       deallocate( dij_num_nbors_matrix, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating dij_num_nbors_matrix in param_dealloc')
    end if


    return
  end subroutine dij_dealloc


  !==================================================================!
  subroutine dij_get()
    !==================================================================!
    !                                                                  !
    !                                                                  ! 
    !                                                                  !
    !===================================================================  
    use mc_constants,  only : eps4,eps5,k_B
    use mc_io,         only : io_error,io_stopwatch,seedname,stdout
    use mc_parameters, only : real_lattice,atoms_pos_frac,recip_lattice,&
                              atoms_species_num,supercell_size,num_supercell,num_atoms,&
                              dij_1atom_site,dij_2atom_site,dij,dij_vectors_cart,&
                              coordinate,paramconfac,lenconfac,lboundary_open
    use mc_utility,    only : utility_frac_to_cart,utility_cart_to_frac
    use mc_jij,        only : mag_moments_matrix

    implicit none

    ! Variables that are private

    integer              :: l,m,n,nnx                             
    integer              :: isupercell,ibors,jbors,i,j         
    integer              :: nsp,nat          
    integer              :: lmn(3,num_supercell)
    integer              :: num(num_total_atoms)
    integer              :: num_first_atom,num_second_atom          
    integer              :: num_first_supercell,num_second_supercell          
    real(kind=dp)        :: pos_first_atom(3),pos_second_atom(3)
    real(kind=dp)        :: pos(3),pos_frac(3),param_tmp
    logical              :: atom_found,found

    call io_stopwatch('dij: get',1)

    allocate(dij_tot_matrix(3,num_total_atoms,dij_num_nbors), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating dij_tot_matrix in dij_get')
    allocate(dij_tot_nbors_matrix(num_total_atoms,dij_num_nbors), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating dij_tot_matrix in dij_get')


    ! Convert parameter to Kelvin 
    dij = dij/k_B

    ! Marking supercell 
    nnx = 0
    do n=0,supercell_size(3)-1
       do m=0,supercell_size(2)-1  
          do l=0,supercell_size(1)-1
             nnx = nnx+1
             lmn(1,nnx) = l
             lmn(2,nnx) = m
             lmn(3,nnx) = n
          enddo
       enddo
    enddo

    dij_tot_matrix = 0.0_dp
    dij_tot_nbors_matrix = 0
    num = 0

    do isupercell=1,num_supercell  ! Supercell
     do ibors=1,dij_num_nbors
        pos_first_atom(:)  = dij_1atom_site(:,ibors)+lmn(:,isupercell)
        pos_second_atom(:) = dij_2atom_site(:,ibors)+lmn(:,isupercell)

       !! Now we should specify the numeber of first atom
       if (abs(pos_first_atom(1)-nint(pos_first_atom(1))) .lt. eps5)  pos_first_atom(1) = nint(pos_first_atom(1))
       if (abs(pos_first_atom(2)-nint(pos_first_atom(2))) .lt. eps5)  pos_first_atom(2) = nint(pos_first_atom(2))
       if (abs(pos_first_atom(3)-nint(pos_first_atom(3))) .lt. eps5)  pos_first_atom(3) = nint(pos_first_atom(3))

       l = floor(pos_first_atom(1))
       m = floor(pos_first_atom(2))
       n = floor(pos_first_atom(3))

       !! Position of first atom respect to own supercell 
       pos(1) = pos_first_atom(1)-l
       pos(2) = pos_first_atom(2)-m
       pos(3) = pos_first_atom(3)-n

       !! Number of suprecell of first atom
       l = mod(l+10*supercell_size(1),supercell_size(1))
       m = mod(m+10*supercell_size(2),supercell_size(2))
       n = mod(n+10*supercell_size(3),supercell_size(3))
       num_first_supercell = n*supercell_size(2)*supercell_size(1)+m*supercell_size(1)+l+1
       !! Number of first atom 
       nnx = 0
       atom_found = .false.
       ok: do nsp=1,num_species
              do nat=1,atoms_species_num(nsp)
                 nnx = nnx+1
                 if (abs(pos(1)-atoms_pos_frac(1,nsp,nat)) .le. eps4 .and. &
                     abs(pos(2)-atoms_pos_frac(2,nsp,nat)) .le. eps4 .and. &
                     abs(pos(3)-atoms_pos_frac(3,nsp,nat)) .le. eps4 ) then
                    num_first_atom = (num_first_supercell-1)*num_atoms+nnx
                    atom_found = .true.
                    exit ok
                 endif
              enddo
       enddo ok
       if (.not.atom_found) then
          write(stdout,'(1x,a)') 'The Dij_parameters block of file '//trim(seedname)//'.mcin '&
                    //'contained Unrecognised Atomic Position '
          if (coordinate .eq. 'frac') then
             write(stdout,'(1x,a)') ' in Fractional Coordinate:'
             write(stdout,'(1x,a3,3F12.9)') 'f1=',dij_1atom_site(:,ibors)
          else
             call utility_frac_to_cart (dij_1atom_site(:,ibors),pos(:),real_lattice)
             if (lenconfac .eq. 1.0_dp) then
                write(stdout,'(1x,a)') ' in Cartesian Coordinate (Ang):'
             else
                write(stdout,'(1x,a)') ' in Cartesian Coordinate (Bohr):'
             endif
             write(stdout,'(1x,a3,3F12.9)') 'c1=',pos(:)*lenconfac
          endif
          call io_error('Unrecognised Atomic Position found in Dij_parameters block')
       endif

       num(num_first_atom) = num(num_first_atom) + 1
       !! Now we should specify the numeber of second atom
       if (abs(pos_second_atom(1)-nint(pos_second_atom(1))) .lt. eps5)  pos_second_atom(1) = nint(pos_second_atom(1))
       if (abs(pos_second_atom(2)-nint(pos_second_atom(2))) .lt. eps5)  pos_second_atom(2) = nint(pos_second_atom(2))
       if (abs(pos_second_atom(3)-nint(pos_second_atom(3))) .lt. eps5)  pos_second_atom(3) = nint(pos_second_atom(3))

       l = floor(pos_second_atom(1))
       m = floor(pos_second_atom(2))
       n = floor(pos_second_atom(3))

       if (lboundary_open .and. &
           (l.lt.0 .or. l.eq.supercell_size(1) .or. &
            m.lt.0 .or. m.eq.supercell_size(2)  .or. &
            n.lt.0 .or. n.eq.supercell_size(3))) cycle

       !! Position of second atom respect to own supercell 
       pos(1) = pos_second_atom(1)-l
       pos(2) = pos_second_atom(2)-m
       pos(3) = pos_second_atom(3)-n

       !! Number of suprecell of second atom
       l = mod(l+10*supercell_size(1),supercell_size(1))
       m = mod(m+10*supercell_size(2),supercell_size(2))
       n = mod(n+10*supercell_size(3),supercell_size(3))
       num_second_supercell = n*supercell_size(2)*supercell_size(1)+m*supercell_size(1)+l+1
       !! Number of second atom 
       nnx = 0
       atom_found = .false.
       ok1: do nsp=1,num_species
              do nat=1,atoms_species_num(nsp)
                 nnx = nnx+1
                 if (abs(pos(1)-atoms_pos_frac(1,nsp,nat)) .le. eps4 .and. &
                     abs(pos(2)-atoms_pos_frac(2,nsp,nat)) .le. eps4 .and. &
                     abs(pos(3)-atoms_pos_frac(3,nsp,nat)) .le. eps4 ) then
                    num_second_atom = (num_second_supercell-1)*num_atoms+nnx
                    dij_tot_matrix(:,num_first_atom,num(num_first_atom)) = & 
                         dij_vectors_cart(:,ibors)*dij(ibors) & 
                          /(abs(mag_moments_matrix(num_first_atom)*mag_moments_matrix(num_second_atom)))
                    dij_tot_nbors_matrix(num_first_atom,num(num_first_atom)) = num_second_atom 
                    atom_found = .true.
                    exit ok1
                 endif
              enddo
       enddo ok1
       if (.not.atom_found) then
          write(stdout,'(1x,a)') 'The Dij_parameters block of file '//trim(seedname)//'.mcin '&
                    //'contained Unrecognised Atomic Position '
          if (coordinate .eq. 'frac') then
             write(stdout,'(1x,a)') ' in Fractional Coordinate:'
             write(stdout,'(1x,a3,3F12.9)') 'f2=',dij_2atom_site(:,ibors)
          else
             call utility_frac_to_cart (dij_2atom_site(:,ibors),pos(:),real_lattice)
             if (lenconfac .eq. 1.0_dp) then
                write(stdout,'(1x,a)') ' in Cartesian Coordinate (Ang):'
             else
                write(stdout,'(1x,a)') ' in Cartesian Coordinate (Bohr):'
             endif
             write(stdout,'(1x,a3,3F12.9)') 'c2=',pos(:)*lenconfac
          endif
          call io_error('Unrecognised Atomic Position found in Dij_parameters block')
       endif
     enddo ! ibors
    enddo  ! isupercell

    call dij_neighbors_get()

    ! check the matrix is symmetric
    if (.not.lboundary_open) then
       do i=1,num_total_atoms
          do ibors=1,dij_num_nbors_matrix(i)
             found = .False.
             j = dij_tot_nbors_matrix(i,ibors)
             do jbors=1,dij_num_nbors_matrix(j)
                if (dij_nbors_matrix(j,jbors) .eq. i) then
                   found = .True.
                   exit
                endif
             enddo
             if (.not.found) then
                write(stdout,'(2x,a)') 'Dij_parameters block is not symmetric.'
                write(stdout,'(2x,a,i4,a,i4,a,5x,a,i4,a,i4,a)') 'Asymmetric pairs are (',i,',',j,')','(',j,',',i,')'
                call io_error('Found asymmetric pairs in Dij_parameters block')
             endif
             if (any(abs(dij_matrix(:,i,ibors)+dij_matrix(:,j,jbors)) > eps4)) then
                param_tmp = sqrt(dot_product(dij_matrix(:,i,ibors),dij_matrix(:,i,ibors)))
                pos(:) = dij_matrix(:,i,ibors)/param_tmp
                if (paramconfac.eq.1.0_dp) then
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'Dij ',param_tmp*k_B*paramconfac,&
                                                          ' (eV) between atoms ',i,' and ',j
                else
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'Dij is ',param_tmp*k_B*paramconfac,&
                                                          ' (Ryd) between atoms ',i,' and ',j
                endif
                if (coordinate .eq. 'frac') then
                   call utility_cart_to_frac (pos(:),pos_frac(:),recip_lattice)
                   write(stdout,'(1x,a,3F12.8)') ' with Dij vecotrs in fractional coordinate:',pos_frac(:)
                else
                   if (lenconfac .eq. 1.0_dp) then
                      write(stdout,'(1x,a,3F12.6)') ' with Dij vecotrs in Cartesian coordinate (Ang) is:',pos(:)*lenconfac
                   else
                      write(stdout,'(1x,a,3F12.6)') ' with Dij vecotrs in Cartesian coordinate (Bohr) is:',pos(:)*lenconfac
                   endif
                endif
                write(stdout,'(1x,a)') ' ' 
                param_tmp = sqrt(dot_product(dij_matrix(:,j,jbors),dij_matrix(:,j,jbors)))
                pos(:) = dij_matrix(:,j,jbors)/param_tmp
                if (paramconfac.eq.1.0_dp) then
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'but Dij ',param_tmp*k_B*paramconfac,&
                                                          ' (eV) between atoms ',j,' and ',i
                else
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'but Dij is ',param_tmp*k_B*paramconfac,&
                                                          ' (Ryd) between atoms ',j,' and ',i
                endif
                if (coordinate .eq. 'frac') then
                   call utility_cart_to_frac (pos(:),pos_frac(:),recip_lattice)
                   write(stdout,'(1x,a,3F12.8)') ' with Dij vecotrs in fractional coordinate:',pos_frac(:)
                else
                   if (lenconfac .eq. 1.0_dp) then
                      write(stdout,'(1x,a,3F12.6)') ' with Dij vecotrs in Cartesian coordinate (Ang) is:',pos(:)*lenconfac
                   else
                      write(stdout,'(1x,a,3F12.6)') ' with Dij vecotrs in Cartesian coordinate (Bohr) is:',pos(:)*lenconfac
                   endif
                endif
                call io_error('Found asymmetric pairs in Dij/Dij_vectors block')
             endif
          enddo
       enddo
    endif

    deallocate( dij_tot_matrix, stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating dij_tot_matrix in dij_get')
    deallocate( dij_tot_nbors_matrix, stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating dij_tot_matrix in dij_get')

    call io_stopwatch('dij: get',2)
    
    return
    
  contains

    !==================================================================!
    subroutine dij_neighbors_get()
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  
      use mc_constants,   only : zero

      implicit none
      integer       :: nbors_max 
      integer       :: nat1,nat2,counter 
      real(kind=dp) :: d_temp(3)

      dij_num_nbors_matrix = 0
      do nat1=1,num_total_atoms
         counter = 0
         do nat2=1,dij_num_nbors
            d_temp(:) = dij_tot_matrix(:,nat1,nat2)
            if(any(d_temp .ne. zero)) then
               counter = counter+1
            endif 
         enddo
         dij_num_nbors_matrix(nat1) = counter
      enddo

      ! Maximunm number of neighbours
      nbors_max = maxval(dij_num_nbors_matrix)
      allocate(dij_nbors_matrix(num_total_atoms,nbors_max), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating dij_nbors_matrix in dij_neighbors_get')
      allocate(dij_matrix(3,num_total_atoms,nbors_max), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating dij_matrix in dij_neighbors_get')
      dij_matrix = 0.0_dp
      dij_nbors_matrix = 0

      do nat1=1,num_total_atoms
         counter = 0
         do nat2=1,dij_num_nbors
            d_temp(:) = dij_tot_matrix(:,nat1,nat2)
            if(any(d_temp .ne. zero)) then
              counter = counter+1
              dij_nbors_matrix(nat1,counter) = dij_tot_nbors_matrix(nat1,nat2)
              dij_matrix(:,nat1,counter) = dij_tot_matrix(:,nat1,nat2)
            endif
         enddo
      enddo

      return

    end subroutine dij_neighbors_get 

  end subroutine dij_get

end module mc_dij      



