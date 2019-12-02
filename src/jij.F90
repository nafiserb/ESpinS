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

module mc_jij      

  use mc_constants,  only : dp
  use mc_io,         only : io_error
  use mc_parameters, only : num_total_atoms,num_species,jij_num_nbors,&
                            jij_shell_max,lspincorr,spin_glass

  implicit none

  private

  real(kind=dp), allocatable, public, save :: mag_moments_matrix(:)
  real(kind=dp), allocatable, public, save :: jij_matrix(:,:)
  integer,       allocatable, public, save :: jij_nbors_matrix(:,:)
  integer,       allocatable, public, save :: jij_num_nbors_matrix(:)
  integer,       allocatable, public, save :: shell_nbors_matrix(:,:,:)
  integer,       allocatable, public, save :: shell_num_nbors_matrix(:,:)
  integer,       allocatable, public, save :: shell_num_nbors_kind(:,:,:)
  integer,       allocatable, public, save :: kind_atoms(:)

  public :: jij_get
  public :: jij_setup
  public :: jij_dealloc

  ! Private Data
  integer      , allocatable :: shell_matrix(:,:)
  real(kind=dp), allocatable :: jij_tot_matrix(:,:)
  real(kind=dp), allocatable :: sigma_tot_matrix(:,:)
  real(kind=dp), allocatable :: sigma_matrix(:,:)
  integer      , allocatable :: jij_tot_nbors_matrix(:,:)
  integer                    :: ierr


contains 

  !============================================!
  subroutine jij_setup()
    !============================================!

    implicit none

    allocate(mag_moments_matrix(num_total_atoms), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating mag_moments_matrix in jij_setup')
    allocate(jij_num_nbors_matrix(num_total_atoms), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating jij_num_nbors_matrix in jij_setup')
    if (lspincorr) then
      allocate(shell_num_nbors_matrix(jij_shell_max,num_total_atoms), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating shell_num_nbors_matrix in jij_setup')
      allocate(shell_num_nbors_kind(jij_shell_max,num_species,num_species), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating shell_num_nbors_kind in jij_setup')
      allocate(kind_atoms(num_total_atoms), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating kind_atoms in jij_setup')
    endif

    return
  end subroutine jij_setup

  !============================================!
  subroutine jij_dealloc()
    !============================================!

    implicit none

    ! Deallocate arrays that are public 

    if ( allocated( jij_matrix ) ) then
       deallocate( jij_matrix, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating jij_matrix in jij_dealloc')
    end if
    if ( allocated( mag_moments_matrix ) ) then
       deallocate( mag_moments_matrix, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating mag_moments_matrix in jij_dealloc')
    end if
    if ( allocated( jij_nbors_matrix ) ) then
       deallocate( jij_nbors_matrix, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating jij_nbors_matrix in jij_dealloc')
    end if
    if ( allocated( jij_num_nbors_matrix ) ) then
       deallocate( jij_num_nbors_matrix, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating jij_num_nbors_matrix in jij_dealloc')
    end if
    if ( allocated ( kind_atoms ) ) then
       deallocate (kind_atoms, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating kind_atoms in jij_dealloc')
    end if
    if ( allocated( shell_nbors_matrix ) ) then
       deallocate( shell_nbors_matrix, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating shell_nbors_matrix in jij_dealloc')
    end if
    if ( allocated( shell_num_nbors_matrix ) ) then
       deallocate( shell_num_nbors_matrix, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating shell_num_nbors_matrix in jij_dealloc')
    end if
    if ( allocated( shell_num_nbors_kind ) ) then
       deallocate( shell_num_nbors_kind, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating shell_num_nbors_kind in jij_dealloc')
    end if

    return
  end subroutine jij_dealloc

  !==================================================================!
  subroutine jij_get()
    !==================================================================!
    !                                                                  !
    !                                                                  ! 
    !                                                                  !
    !===================================================================  
    use mc_constants,  only : eps4,eps5,k_B
    use mc_io,         only : io_stopwatch,stdout,seedname
    use mc_parameters, only : real_lattice,atoms_pos_frac,atoms_species_num,&
                              supercell_size,num_atoms,num_supercell,&
                              jij_1atom_site,jij_2atom_site,jij,jij_shell,&
                              spin_glass_sigma,mag_moments,coordinate,&
                              paramconfac,lenconfac,lboundary_open
    use mc_utility,    only : utility_frac_to_cart
    implicit none

    ! Variables that are private

    integer              :: l,m,n,nnx,nsp,nat                        
    integer              :: isupercell,ibors,jbors,i,j
    integer              :: lmn(3,num_supercell)
    integer              :: num(num_total_atoms)
    integer              :: num_first_atom,num_second_atom
    integer              :: num_first_supercell,num_second_supercell          
    real(kind=dp)        :: pos_first_atom(3),pos_second_atom(3)
    real(kind=dp)        :: pos(3)
    logical              :: atom_found,found

    call io_stopwatch('jij: get',1)

    allocate(jij_tot_matrix(num_total_atoms,jij_num_nbors), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating jij_tot_matrix in jij_get')
    allocate(jij_tot_nbors_matrix(num_total_atoms,jij_num_nbors), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating jij_tot_nbors_matrix in jij_get')
    allocate(shell_matrix(num_total_atoms,jij_num_nbors), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating shell_matrix in jij_get')
    allocate(sigma_tot_matrix(num_total_atoms,jij_num_nbors), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating sigma_tot_matrix in jij_get')

    ! Convert parameter(eV) to Kelvin 
    jij            = jij/k_B
    spin_glass_sigma = spin_glass_sigma/k_B

    ! Magnetic mpoments matrix
    nnx = 0
    do isupercell=1,num_supercell
       do nsp=1,num_species   
          do nat=1, atoms_species_num(nsp)
             nnx = nnx+1
             mag_moments_matrix(nnx) = mag_moments(nsp,nat)
             if(lspincorr) kind_atoms(nnx) = nsp
          enddo
       enddo
    enddo

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
                    
    jij_tot_matrix = 0.0_dp
    jij_tot_nbors_matrix = 0
    num = 0
    shell_matrix = 0
    sigma_tot_matrix = 0.0_dp
    do isupercell=1,num_supercell  ! Supercell
       do ibors=1,jij_num_nbors
          pos_first_atom(:) = jij_1atom_site(:,ibors)+lmn(:,isupercell)
          pos_second_atom(:) = jij_2atom_site(:,ibors)+lmn(:,isupercell)

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
             write(stdout,'(1x,a)') 'The Jij_parameters block of file '//trim(seedname)//'.mcin '&
                      //'contained Unrecognised Atomic Position '
             if (coordinate .eq. 'frac') then
                write(stdout,'(1x,a)') ' in Fractional Coordinate:'
                write(stdout,'(1x,a3,3F12.9)') 'f1=',jij_1atom_site(:,ibors)
             else
                call utility_frac_to_cart (jij_1atom_site(:,ibors),pos(:),real_lattice)
                if (lenconfac .eq. 1.0_dp) then
                   write(stdout,'(1x,a)') ' in Cartesian Coordinate (Ang):'
                else
                   write(stdout,'(1x,a)') ' in Cartesian Coordinate (Bohr):'
                endif
                write(stdout,'(1x,a3,3F12.9)') 'c1=',pos(:)*lenconfac
             endif
             call io_error('Unrecognised Atomic Position found in Jij_parameters block')
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
                       jij_tot_matrix(num_first_atom,num(num_first_atom)) = & 
                            jij(ibors)/(abs(mag_moments_matrix(num_first_atom)*mag_moments_matrix(num_second_atom)))
                       jij_tot_nbors_matrix(num_first_atom,num(num_first_atom)) = num_second_atom
                       shell_matrix(num_first_atom,num(num_first_atom)) = jij_shell(ibors)
                       sigma_tot_matrix(num_first_atom,num(num_first_atom)) = spin_glass_sigma(ibors)
                       atom_found = .true.
                       exit ok1
                    endif
                 enddo
          enddo ok1
          if (.not.atom_found) then
             write(stdout,'(1x,a)') 'The Jij_parameters block of file '//trim(seedname)//'.mcin '&
                       //'contained Unrecognised Atomic Position '
             if (coordinate .eq. 'frac') then
                write(stdout,'(1x,a)') ' in Fractional Coordinate:'
                write(stdout,'(1x,a3,3F12.9)') 'f2=',jij_2atom_site(:,ibors)
             else
                call utility_frac_to_cart (jij_2atom_site(:,ibors),pos(:),real_lattice)
                if (lenconfac .eq. 1.0_dp) then
                   write(stdout,'(1x,a)') ' in Cartesian Coordinate (Ang):'
                else
                   write(stdout,'(1x,a)') ' in Cartesian Coordinate (Bohr):'
                endif
                write(stdout,'(1x,a3,3F12.9)') 'c2=',pos(:)*lenconfac
             endif
             call io_error('Unrecognised Atomic Position found in Jij_parameters block')
          endif
       enddo ! ibors=1,jij_num_nbors
    enddo   ! isupercell=1,num_supercells

    call jij_neighbors_get()
    if (lspincorr) call shell_neighbors_get()

    ! check the matrix is symmetric
    if (.not.lboundary_open) then
       do i=1,num_total_atoms
          do ibors=1,jij_num_nbors_matrix(i)
             found = .False.
             j = jij_nbors_matrix(i,ibors)
             do jbors=1,jij_num_nbors_matrix(j)
                if (jij_nbors_matrix(j,jbors) .eq. i) then
                   found = .True.
                   exit
                endif
             enddo
             if (.not.found) then
                write(stdout,'(2x,a)') 'Jij_parameters block is not symmetric.'
                write(stdout,'(2x,a,i4,a,i4,a,5x,a,i4,a,i4,a)') 'Asymmetric pairs are (',i,',',j,')','(',j,',',i,')'
                call io_error('Found asymmetric pairs in Jij_parameters block')
             endif
             if (abs(jij_matrix(i,ibors)-jij_matrix(j,jbors)) > eps4) then
                if (paramconfac.eq.1.0_dp) then
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'Jij is ',jij_matrix(i,ibors)*k_B*paramconfac,&
                                                          ' (eV) between atoms ',i,' and ',j
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'but it is ',jij_matrix(j,jbors)*k_B*paramconfac,&
                                                          ' (eV) between atoms ',j,' and ',i
                else
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'Jij is ',jij_matrix(i,ibors)*k_B*paramconfac,&
                                                          ' (Ryd) between atoms ',i,' and ',j
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'but it is ',jij_matrix(j,jbors)*k_B*paramconfac,&
                                                          ' (Ryd) between atoms ',j,' and ',i
                endif
                call io_error('Found asymmetric pairs in Jij_parameters block')
             endif
             if (spin_glass .and. &
                 abs(sigma_matrix(i,ibors)-sigma_matrix(j,jbors)) > eps4) then
                if (paramconfac.eq.1.0_dp) then
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'Sigma is ',sigma_matrix(i,ibors)*k_B*paramconfac,&
                                                          ' (eV) between atoms ',i,' and ',j
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'but it is ',sigma_matrix(j,jbors)*k_B*paramconfac,&
                                                          ' (eV) between atoms ',j,' and ',i
                else
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'Sigma is ',sigma_matrix(i,ibors)*k_B*paramconfac,&
                                                          ' (Ryd) between atoms ',i,' and ',j
                   write(stdout,'(2x,a,F12.8,a,i4,a,i4)') 'but it is ',sigma_matrix(j,jbors)*k_B*paramconfac,&
                                                          ' (Ryd) between atoms ',j,' and ',i
                endif
                call io_error('Found asymmetric pairs in Jij_parameters block')
             endif
          enddo
       enddo
    endif

    if (spin_glass) call jij_spin_glass()

    deallocate( jij_tot_matrix, stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating jij_tot_matrix in jij_get')
    deallocate( jij_tot_nbors_matrix, stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating jij_tot_nbors_matrix in jij_get')
    deallocate( shell_matrix, stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating shell_matrix in jij_get')
    deallocate( sigma_tot_matrix, stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating sigma_tot_matrix in jij_get')
    if ( allocated( sigma_matrix ) ) then
       deallocate( sigma_matrix, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating sigma_matrix in jij_get')
    endif
    call io_stopwatch('jij: get',2)
     
    return
     
    contains

    !==================================================================!
    subroutine jij_neighbors_get()
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  
      use mc_constants,   only : zero
     
      implicit none
      integer :: nbors_max 
      integer :: nat1,nat2,counter 
     
      jij_num_nbors_matrix = 0
      do nat1=1,num_total_atoms
         counter = 0
         do nat2=1,jij_num_nbors
            if(jij_tot_matrix(nat1,nat2) .ne. zero) then
              counter = counter+1
            endif 
         enddo
         jij_num_nbors_matrix(nat1) = counter
      enddo
     
      ! Maximunm number of neighbours
      nbors_max = maxval(jij_num_nbors_matrix)
      allocate(jij_nbors_matrix(num_total_atoms,nbors_max), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating jij_nbors_matrix in jij_neighbors_get')
      allocate(jij_matrix(num_total_atoms,nbors_max), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating jij_matrix in jij_neighbors_get')
      if (spin_glass) then
         allocate(sigma_matrix(num_total_atoms,nbors_max), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating sigma_matrix in jij_neighbors_get')
         sigma_matrix = 0.0_dp
      endif
      jij_matrix = 0.0_dp
      jij_nbors_matrix = 0
     
      do nat1=1,num_total_atoms
         counter = 0
         do nat2=1,jij_num_nbors
            if( jij_tot_matrix(nat1,nat2) .ne. zero ) then
              counter = counter+1
              jij_nbors_matrix(nat1,counter) = jij_tot_nbors_matrix(nat1,nat2)
              jij_matrix(nat1,counter) = jij_tot_matrix(nat1,nat2)
              if(spin_glass) sigma_matrix(nat1,counter) = sigma_tot_matrix(nat1,nat2)
            endif 
         enddo
      enddo
     
     
      return
     
    end subroutine jij_neighbors_get 

    !==================================================================!
    subroutine shell_neighbors_get()
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !==================================================================!  
      use mc_constants,   only : zero

      implicit none
      integer :: nat1,nat2,nsp1,nsp2,counter(jij_shell_max),ishell
      integer :: nbors_max,counter_kind(jij_shell_max,num_species) 

      shell_num_nbors_matrix = 0
      shell_num_nbors_kind = 0

      do nat1=1,num_total_atoms
         counter = 0
         counter_kind = 0
         nsp1 = kind_atoms(nat1)
         do nat2=1,size(jij) !jij_num_nbors   
            nsp2 = kind_atoms(jij_tot_nbors_matrix(nat1,nat2))
            do ishell=1,jij_shell_max
               if ((shell_matrix(nat1,nat2) .ne. 0) .and. &
                   (shell_matrix(nat1,nat2) .eq. ishell)) then
                  counter(ishell) = counter(ishell)+1
                  counter_kind(ishell,nsp2) = counter_kind(ishell,nsp2)+1
               endif 
            enddo
         enddo
         shell_num_nbors_matrix(:,nat1) = counter(:)
         shell_num_nbors_kind(:,:,nsp1) = counter_kind(:,:)
      enddo

      ! Maximunm number of neighbours
      nbors_max = maxval(shell_num_nbors_matrix)
      allocate(shell_nbors_matrix(jij_shell_max,num_total_atoms,nbors_max), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating shell_nbors_matrix in shell_get')
      shell_nbors_matrix = 0

      do nat1=1,num_total_atoms
         counter = 0
         do nat2=1,size(jij) !jij_num_nbors
            do ishell=1,jij_shell_max
               if ((shell_matrix(nat1,nat2) .ne. 0) .and. &
                   (shell_matrix(nat1,nat2) .eq. ishell)) then
                  counter(ishell) = counter(ishell)+1
                  shell_nbors_matrix(ishell,nat1,counter(ishell)) = jij_tot_nbors_matrix(nat1,nat2)
               endif 
            enddo
         enddo
      enddo
  
      return

    end subroutine shell_neighbors_get 

    !==================================================================!
    subroutine jij_spin_glass()   
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  
      use mc_parameters,   only : spin_glass_seed,lspin_glass_seed
      use mtprng,          only : mtprng_init,mtprng_state,mtprng_normal
      use stdtypes

      implicit none
      integer            :: i,ibors,j,jbors
      logical            :: found
      real(kind=dp)      :: rand
      type(mtprng_state) :: state1

      if (.not.lspin_glass_seed) then
        open(89,FILE='/dev/urandom',ACCESS='stream',FORM='UNFORMATTED')
        read(89) spin_glass_seed
        close(89)
      endif

      call mtprng_init(spin_glass_seed,state1)

      do i=1,num_total_atoms
         do ibors=1,jij_num_nbors_matrix(i)
            found = .False.
            j = jij_nbors_matrix(i,ibors)
            if ( i.lt.j) then
               rand = mtprng_normal(state1)
               jij_matrix(i,ibors) = rand*sigma_matrix(i,ibors) + jij_matrix(i,ibors) ! 
               do jbors=1,jij_num_nbors_matrix(j)
                  if (jij_nbors_matrix(j,jbors) .eq. i) then
                     found = .True.
                     jij_matrix(j,jbors) = jij_matrix(i,ibors)
                     exit
                  endif
               enddo

               if (.not.found) then
                  write(stdout,'(2x,a)') 'Jij_parameters block is not symmetric.'
                  write(stdout,'(2x,a,i4,a,i4,a,5x,a,i4,a,i4,a)') 'Asymmetric pairs are (',i,',',j,')','(',j,',',i,')'
                  call io_error('Found asymmetric pairs in Jij_parameters block')
               endif
            endif
         enddo
      enddo

      return

    end subroutine jij_spin_glass   

  end subroutine jij_get

end module mc_jij      


