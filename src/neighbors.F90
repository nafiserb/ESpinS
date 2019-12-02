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

module mc_neighbors

  use mc_constants, only : dp
  use mc_parameters

  implicit none

  private



  public :: neighbors_get
  public :: neighbors_check
  public :: neighbors_write
  public :: mcin_write
  public :: inp2mcin_write
  public :: neighbors_dealloc


  integer, parameter                :: nsupcell=5
  integer                           :: lmn(3,(2*nsupcell+1)**3)
  integer, public,allocatable       :: natlist(:,:,:,:), nsplist(:,:,:,:),ncell(:,:,:,:,:)
  integer, public,allocatable       :: multi(:,:,:)
  real(kind=dp), public,allocatable :: dnn(:,:,:)

contains 
  !==================================================================!
  subroutine neighbors_get()
    !==================================================================!
    !                                                                  !
    !                                                                  ! 
    !                                                                  !
    !===================================================================  
    use mc_io,        only : io_error,io_stopwatch
    implicit none

    ! Variables that are private

    integer       :: nlist,l,m,n,nsp,nat,nsp1,nat1,shells_max
!    integer :: nlist,l,m,n,ndnn,nsp,nat,ndnntot,nsp1,nat1
    integer       :: counter,loop,nnx,ierr
    real(kind=dp) :: vkpp(3),vkpp2(3)
    real(kind=dp) :: dist, dnn0,dnn1
    real(kind=dp), parameter :: eta=99999999.0_dp    ! eta = very large 



    call io_stopwatch('neighbors: get',1)

    allocate(multi(num_species,maxval(atoms_species_num),maxval(shells)), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating multi in neighbors_get')
    multi=0

    allocate(dnn(num_species,maxval(atoms_species_num),maxval(shells)), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating dnn in neighbors_get')
    dnn=0.0_dp

    ! Sort the cell neighbours so we loop in order of distance from the home shell
    call neighbors_supercell_sort


    ! find the distance between a atoms and its nearest-neighbour shells

    dnn0 = 0.0_dp  
    dnn1 = eta  
    do nsp=1,num_species
       do nat=1,atoms_species_num(nsp)
          dnn0 = 0.0_dp  
          dnn1 = eta 
          shells_max = maxval(shells(:,nsp))  
          do nlist = 1, shells_max
             do loop=1,(2*nsupcell+1)**3
                l=lmn(1,loop);m=lmn(2,loop);n=lmn(3,loop)
                !
                vkpp2=matmul(lmn(:,loop),real_lattice)
                do nsp1=1,num_species
                   do nat1=1,atoms_species_num(nsp1)
                      vkpp=atoms_pos_cart(:,nsp1,nat1)+vkpp2
                      dist= sqrt((atoms_pos_cart(1,nsp,nat)-vkpp(1))**2 &
                                  + (atoms_pos_cart(2,nsp,nat)-vkpp(2))**2 + (atoms_pos_cart(3,nsp,nat)-vkpp(3))**2 )
                      if ( (dist.gt.neighbors_tol) .and. (dist.gt.dnn0 + neighbors_tol) ) then
                         if(dist.lt.dnn1-neighbors_tol) then
                           dnn1=dist  ! found a closer shell
                           counter=0
                         end if
                         if(dist.gt.(dnn1-neighbors_tol) .and. dist.lt.(dnn1+neighbors_tol)) then
                           counter=counter+1 ! count the multiplicity of the shell
                         end if
                      end if
                   enddo
                enddo
             enddo
             dnn(nsp,nat,nlist) = dnn1  
             multi(nsp,nat,nlist)=counter
             dnn0 = dnn1  
             dnn1 = eta  
          enddo
       enddo
    enddo

    allocate(ncell(3,num_species,maxval(atoms_species_num),maxval(shells),maxval(multi)), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating ncell in neighbors_get')
    allocate(nsplist(num_species,maxval(atoms_species_num),maxval(shells),maxval(multi)), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating nsplist in neighbors_get')
    allocate(natlist(num_species,maxval(atoms_species_num),maxval(shells),maxval(multi)), stat=ierr )
    if (ierr/=0) call io_error('Error in allocating natlist in neighbors_get')
    ncell=0;nsplist=0;natlist=0


    do nsp=1,num_species
       do nat=1,atoms_species_num(nsp)
          shells_max = maxval(shells(:,nsp))
          ok: do nlist = 1, shells_max
                 nnx=0
                 dnn1=dnn(nsp,nat,nlist)
                 do loop=1,(2*nsupcell+1)**3
                    l=lmn(1,loop);m=lmn(2,loop);n=lmn(3,loop)
                    !
                    vkpp2=matmul(lmn(:,loop),real_lattice)
                    do nsp1=1,num_species
                       do nat1=1,atoms_species_num(nsp1)
                          vkpp=atoms_pos_cart(:,nsp1,nat1)+vkpp2
                          dist= sqrt((atoms_pos_cart(1,nsp,nat)-vkpp(1))**2 &
                                  + (atoms_pos_cart(2,nsp,nat)-vkpp(2))**2 + (atoms_pos_cart(3,nsp,nat)-vkpp(3))**2 )
!                          if ( (dist.ge.dnn1*(1-neighbors_tol)) .and. (dist.le.dnn1*(1+neighbors_tol)) ) then
                          if ( (dist.ge.dnn1-neighbors_tol) .and. (dist.le.dnn1+neighbors_tol) ) then
                             nnx = nnx + 1
                             natlist(nsp,nat,nlist,nnx) = nat1
                             nsplist(nsp,nat,nlist,nnx) = nsp1
                             ncell(1,nsp,nat,nlist,nnx) = l
                             ncell(2,nsp,nat,nlist,nnx) = m
                             ncell(3,nsp,nat,nlist,nnx) = n
                          endif
                          !if we have the right number of neighbours we can exit
                          if(nnx==multi(nsp,nat,nlist)) cycle ok
                       enddo
                    enddo
                 enddo
          end do ok
       enddo
    enddo


    call io_stopwatch('neighbors: get',2)
    
    return
    
  end subroutine neighbors_get


  !==================================================================!


  !==================================================================!
  subroutine neighbors_check()
    !==================================================================!
    !                                                                  !
    !                                                                  !
    !==================================================================!  
    use mc_io,        only : seedname,stdout,io_error,io_stopwatch
    use mc_utility,   only : utility_sort

    implicit none

    integer              :: nsp,nat,nat1,ndnn,ierr,shells_max
    integer              :: loop,loop2,counter,num_dis
    integer, allocatable :: num_atoms_temp(:),num_atoms(:,:),freq(:)
    call io_stopwatch('neighbors: check',1)
    
    do nsp=1, num_species
       do nat=1, atoms_species_num(nsp)
          do nat1=nat+1, atoms_species_num(nsp)
             shells_max = maxval(shells(:,nsp))
             do ndnn=1,shells_max
                if (multi(nsp,nat,ndnn) .ne. multi(nsp,nat1,ndnn) .or. &
                    abs(dnn(nsp,nat,ndnn)-dnn(nsp,nat1,ndnn))>neighbors_tol) then
                    ! Now I count the number of different positions of nsp type.
                    allocate(num_atoms_temp(atoms_species_num(nsp)), stat=ierr )
                    if (ierr/=0) call io_error('Error in allocating num_atoms_temp in neighbors_check')
                    allocate(num_atoms(atoms_species_num(nsp),atoms_species_num(nsp)), stat=ierr )
                    if (ierr/=0) call io_error('Error in allocating num_atoms in neighbors_check')
                    allocate(freq(atoms_species_num(nsp)), stat=ierr )
                    if (ierr/=0) call io_error('Error in allocating freq in neighbors_check')
                    num_atoms_temp=0;num_atoms=0;freq=0;num_dis=0
                    do loop=1,atoms_species_num(nsp)
                       counter=1
                       num_atoms_temp(counter) = loop
                       do loop2=loop+1,atoms_species_num(nsp)
                          if(multi(nsp,loop,ndnn) .eq. multi(nsp,loop2,ndnn) .and. &
                             abs(dnn(nsp,loop,ndnn)-dnn(nsp,loop2,ndnn))<neighbors_tol) then
                             counter = counter+1
                             num_atoms_temp(counter) = loop2
                             freq(loop2) = -1
                          endif
                        enddo
                        if (freq(loop) .ne. -1) then
                           freq(loop) = counter
                           num_atoms(loop,:)=num_atoms_temp(:)
                           num_dis = num_dis+1
                        endif
                     enddo
                     write(stdout,'(1x,a,i3,a,a,a)')  'The Wyckoff positions of ',nsp, ' type (',trim(atoms_label(nsp)),&
                           ') are not same.'
                     write(stdout,'(1x,a,i3,a)')  'They have',num_dis,' different Wyckoff positions.'
                     do loop=1,atoms_species_num(nsp)
                        if (freq(loop) .ne. -1) then
                           write(stdout,'(1x,a,i3,a)')  'The number of atoms of ',nsp,&
                                 ' type that have same Wyckoff positions are:'
                           write(stdout,*) (num_atoms(loop,loop2),loop2=1,freq(loop))
                        endif
                     enddo 
                     write(stdout,'(a)') ' '
                     deallocate(num_atoms_temp, stat=ierr )
                     if (ierr/=0) call io_error('Error in deallocating num_atoms_temp in neighbors_check')
                     deallocate(num_atoms, stat=ierr )
                     if (ierr/=0) call io_error('Error in deallocating num_atoms in neighbors_check')
                     deallocate(freq, stat=ierr )
                     if (ierr/=0) call io_error('Error in deallocating freq in neighbors_check')
                     call io_error ('Error: inconsistency in atoms type.' &
                                     //'See '//trim(seedname)//'.mcout file for more details')
                endif
             enddo
          enddo
       enddo
    enddo

    call io_stopwatch('neighbors: check',2)
    
    return
    
  end subroutine neighbors_check

  !==================================================================!
  subroutine neighbors_write()
    !==================================================================!
    !                                                                  !
    !                 Writes mcarlo.neigh file                         !
    !                                                                  !
    !==================================================================!  
    use mc_io,        only : io_file_unit,seedname,io_date,io_stopwatch
    use mc_utility,   only : utility_cart_to_frac

    implicit none

    integer           :: nsp,nat,nsp1,nat1,ndnn,vkpp(3),nnx,i,shells_max
    integer           :: neigh_unit
    character (len=9) :: cdate,ctime
    real(kind=dp)     :: neigh_cart(3),neigh_frac(3)

    call io_stopwatch('neighbors: write',1)

    neigh_unit=io_file_unit()
    open(unit=neigh_unit,file=trim(seedname)//'.neigh',form='formatted')

    ! Date and time
    call io_date(cdate,ctime)
    write(neigh_unit,'(4(a),/)') 'File written on ',cdate,' at ',ctime

    write(neigh_unit,*)
    write(neigh_unit,'(36x,a6)') '------'
    write(neigh_unit,'(36x,a6)') 'SYSTEM'
    write(neigh_unit,'(36x,a6)') '------'
    write(neigh_unit,*)
    if (lenconfac.eq.1.0_dp) then
       write(neigh_unit,'(30x,a21)') 'Lattice Vectors (Ang)'
    else
       write(neigh_unit,'(28x,a22)') 'Lattice Vectors (Bohr)'
    endif
    write(neigh_unit,101) 'a_1',(real_lattice(1,I)*lenconfac, i=1,3)
    write(neigh_unit,101) 'a_2',(real_lattice(2,I)*lenconfac, i=1,3)
    write(neigh_unit,101) 'a_3',(real_lattice(3,I)*lenconfac, i=1,3)
    write(neigh_unit,*)
    write(neigh_unit,'(20x,a17,3x,f11.5)',advance='no') &
         'Unit Cell Volume:',cell_volume*lenconfac**3
    if (lenconfac.eq.1.0_dp) then
       write(neigh_unit,'(2x,a7)') '(Ang^3)'
    else
       write(neigh_unit,'(2x,a8)') '(Bohr^3)'
    endif
    write(neigh_unit,*)
    if (lenconfac.eq.1.0_dp) then
       write(neigh_unit,'(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    else
       write(neigh_unit,'(22x,a34)') 'Reciprocal-Space Vectors (Bohr^-1)'
    endif
    write(neigh_unit,101) 'b_1',(recip_lattice(1,I)/lenconfac, i=1,3)
    write(neigh_unit,101) 'b_2',(recip_lattice(2,I)/lenconfac, i=1,3)
    write(neigh_unit,101) 'b_3',(recip_lattice(3,I)/lenconfac, i=1,3)
    write(neigh_unit,*)   ' '

   ! Atoms
    if(num_atoms>0) then
       write(neigh_unit,'(1x,a)') '*----------------------------------------------------------------------------*'
       if (lenconfac.eq.1.0_dp) then
          write(neigh_unit,'(1x,a)') '| Atom Type    Fractional Coordinate          Cartesian Coordinate (Ang)     |'
       else
          write(neigh_unit,'(1x,a)') '| Atom Type    Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
       endif
       write(neigh_unit,'(1x,a)') '+----------------------------------------------------------------------------+'
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             write(neigh_unit,'(1x,a1,1x,a3,1x,i3,3F10.5,3x,a1,1x,3F10.5,3x,a1)') '|',atoms_symbol(nsp),nsp,&
                   atoms_pos_frac(:,nsp,nat),'|',atoms_pos_cart(:,nsp,nat)*lenconfac,'|'
          end do
       end do
       write(neigh_unit,'(1x,a)') '*----------------------------------------------------------------------------*'
    else
       write(neigh_unit,'(25x,a)') 'No atom positions specified'
    end if
    write(neigh_unit,*) ' '

    ! Neighbours
    if(num_atoms>0) then
    write(neigh_unit,'(1x,a)') '+----------------------------------------------------------------------------+'
    write(neigh_unit,'(1x,a)') '|                    Distance to Nearest-Neighbour Shells                    |'
    write(neigh_unit,'(1x,a)') '|                    ------------------------------------                    |'
    if (lenconfac.eq.1.0_dp) then
       write(neigh_unit,'(1x,a)') '|           Shell               Distance (Ang)            Multiplicity       |'
       write(neigh_unit,'(1x,a)') '|           -----               --------------            ------------       |'
    else
       write(neigh_unit,'(1x,a)') '|           Shell               Distance (Bohr)           Multiplicity       |'
       write(neigh_unit,'(1x,a)') '|           -----               ---------------           ------------       |'
    endif

    do nsp=1,num_species
       do nat=1,atoms_species_num(nsp)
          shells_max = maxval(shells(:,nsp))
          do ndnn = 1, shells_max
             write(neigh_unit,'(1x,a,11x,i3,17x,f10.6,19x,i4,12x,a)') '|',ndnn,dnn(nsp,nat,ndnn)*lenconfac,&
                   multi(nsp,nat,ndnn),'|'
          enddo
       enddo
    enddo
    write(neigh_unit,'(1x,a)') '+----------------------------------------------------------------------------+'


       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             write(neigh_unit,'(1x,a)') '*----------------------------------------------------------------------------*'
             if (lenconfac.eq.1.0_dp) then
                write(neigh_unit,'(1x,a)') '| Atom Type    Fractional Coordinate          Cartesian Coordinate (Ang)     |'
             else 
                write(neigh_unit,'(1x,a)') '| Atom Type    Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
             endif
             write(neigh_unit,'(1x,a)') '+----------------------------------------------------------------------------+'
             write(neigh_unit,'(1x,a1,1x,a3,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',atoms_symbol(nsp),nsp,&
                   atoms_pos_frac(:,nsp,nat),'|',atoms_pos_cart(:,nsp,nat)*lenconfac,'|'
             shells_max = maxval(shells(:,nsp))
             do ndnn = 1, shells_max
                !write(neigh_unit,'(1x,a)') '+----------------------------------------------------------------------------+'
                write(neigh_unit,'(1x,a)') '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                if (lenconfac.eq.1.0_dp) then
                   write(neigh_unit,'(1x,a,14x,i3,4x,2a,F10.6,16x,a)') '|',ndnn,' Neighbors & ','Distance(Ang) = '&
                        ,dnn(nsp,nat,ndnn)*lenconfac,'|'
                else
                   write(neigh_unit,'(1x,a,14x,i3,4x,2a,F10.6,15x,a)') '|',ndnn,' Neighbors & ','Distance(Bohr) = '&
                         ,dnn(nsp,nat,ndnn)*lenconfac,'|'
                endif
                write(neigh_unit,'(1x,a)') '|                ---------------------------------------------               |'

                do nnx=1 , multi(nsp,nat,ndnn)
                   nat1 =  natlist(nsp,nat,ndnn,nnx) 
                   nsp1 =  nsplist(nsp,nat,ndnn,nnx)
                   vkpp(1)  =  ncell(1,nsp,nat,ndnn,nnx)
                   vkpp(2)  =  ncell(2,nsp,nat,ndnn,nnx)
                   vkpp(3)  =  ncell(3,nsp,nat,ndnn,nnx)
                   neigh_cart(:)=matmul(vkpp(:),real_lattice)+atoms_pos_cart(:,nsp1,nat1)
                   call   utility_cart_to_frac(neigh_cart(:),neigh_frac(:),recip_lattice) 
                   write(neigh_unit,'(1x,a1,1x,a3,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',atoms_symbol(nsp1),nsp1,neigh_frac(:),&
                        '|',neigh_cart(:)*lenconfac,'|'
                end do !nnx
             enddo !ndnn
          end do ! nat
       enddo ! nsp
       write(neigh_unit,'(1x,a)') '*----------------------------------------------------------------------------*'
    else  !!!
       write(neigh_unit,'(25x,a)') 'No atom positions specified'
    end if
    write(neigh_unit,*) ' '
    close(neigh_unit)

    call xsf_write()
    call io_stopwatch('neighbors: write',2)

101 format(20x,a3,2x,3F11.6)


    return

  end subroutine neighbors_write

  !==================================================================!
  subroutine inp2mcin_write()
    !==================================================================!
    !                                                                  !
    ! Writes mcarlo.nat file                                           !
    !                                                                  ! 
    !                                                                  !
    !                                                                  !
    !===================================================================  
    use mc_io,        only : io_file_unit,seedname,io_stopwatch
    use mc_utility,   only : utility_cart_to_frac

    implicit none

    integer           :: nsp,nat,ndnn,nnx,i,nsp1
    integer           :: inp2mcin_unit
    integer           :: nat_pair(num_species,num_species,maxval(shells))
    logical           :: pair(num_species,num_species,maxval(shells))

    call io_stopwatch('inp2mcin: write',1)

    inp2mcin_unit=io_file_unit()
    open(unit=inp2mcin_unit,file=trim(seedname)//'.inp2.mcin',form='formatted')

    write(inp2mcin_unit,'(2x,a)') 'Begin Unit_Cell_Cart'
    if (lenconfac.ne.1) then
       write(inp2mcin_unit,'(3x,a4)') 'Bohr'
    endif
    write(inp2mcin_unit,'(3x,3f13.8)') (real_lattice(1,I)*lenconfac, i=1,3)
    write(inp2mcin_unit,'(3x,3f13.8)') (real_lattice(2,I)*lenconfac, i=1,3)
    write(inp2mcin_unit,'(3x,3f13.8)') (real_lattice(3,I)*lenconfac, i=1,3)
    write(inp2mcin_unit,'(2x,a)') 'End Unit_Cell_Cart'
    write(inp2mcin_unit,'(2x,a)') ''

    if (coordinate .eq. 'frac') then
       write(inp2mcin_unit,'(2x,a)') 'Begin Atoms_Frac'
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             write(inp2mcin_unit,'(3x,a3,3x,3f12.7,3x,f5.2)') atoms_symbol(nsp),&
                   (atoms_pos_frac(i,nsp,nat),i=1,3),mag_moments(nsp,nat)
          end do
       end do
       write(inp2mcin_unit,'(2x,a)') 'End Atoms_Frac'
    else
       write(inp2mcin_unit,'(2x,a)') 'Begin Atoms_Cart'
       if (lenconfac.ne.1) then
          write(inp2mcin_unit,'(3x,a4)') 'Bohr'
       endif
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             write(inp2mcin_unit,'(3x,a3,3x,3f12.7,3x,f5.2)') atoms_symbol(nsp),&
                   (atoms_pos_cart(i,nsp,nat)*lenconfac,i=1,3),mag_moments(nsp,nat)
          end do
       end do
       write(inp2mcin_unit,'(2x,a)') 'End Atoms_Cart'
    endif
    write(inp2mcin_unit,'(2x,a)') ''


    pair = .false.
    nat_pair = 0
    do nsp=1,num_species
       do nat=1,atoms_species_num(nsp)
          do ndnn = 1, shells(1,nsp)
             do nnx=1 , multi(nsp,nat,ndnn)
                nsp1 = nsplist(nsp,nat,ndnn,nnx)
                pair(nsp,nsp1,ndnn) = .true.
                nat_pair(nsp,nsp1,ndnn) = nat
             enddo
          enddo
       enddo
    enddo

    if (lenconfac  .ne.1.0_dp) then
       write(inp2mcin_unit,'(2x,a)')  'Length_unit     = Bohr'
    !else
    !   write(inp2mcin_unit,'(2x,a)')  '!! Length_unit     = Bohr'
    endif
    if (paramconfac.ne.1.0_dp) then
       write(inp2mcin_unit,'(2x,a)')  'Parameter_unit  = Ryd'
    !else
    !   write(inp2mcin_unit,'(2x,a)')  '!! Parameter_unit  = Ryd'
    endif
    if (coordinate .eq. 'cart') then
       write(inp2mcin_unit,'(2x,a/)')  'Coordinate      = Cart'
    !else
    !   write(inp2mcin_unit,'(2x,a/)')  '!! Coordinate      = Cart'
    endif

    write(inp2mcin_unit,'(2x,a)')    '!! Order_parameter  = .True.'
    write(inp2mcin_unit,'(2x,a)')    '!! Sfactor          = .True.'
    write(inp2mcin_unit,'(2x,a)')    '!! Staggered_m      = .True.'
    write(inp2mcin_unit,'(2x,a)')    '!! Binning_error    = .True.'     
    write(inp2mcin_unit,'(2x,a)')    '!! Spin_correlation = .True.' 
    write(inp2mcin_unit,'(2x,a/)')   '!! Energy_write     = .True.'     

    write(inp2mcin_unit,'(2x,a)')    '## Hamiltonian'
    write(inp2mcin_unit,'(2x,a)')    '!! Boundary         =  Open '
    if(have_biquad) &
       write(inp2mcin_unit,'(2x,a)')    'Ham_bij          = .True.'
    if(have_dm) &
    write(inp2mcin_unit,'(2x,a)')    'Ham_dij          = .True.'
    write(inp2mcin_unit,'(2x,a)')    '!! Ham_singleion    = .True.'
    write(inp2mcin_unit,'(2x,a)')    '!! Ham_field        = .True.'
    if (spin_glass) then
       write(inp2mcin_unit,'(2x,a/)')'Spin_glass       = .True.' 
    else
       write(inp2mcin_unit,'(2x,a/)')    '!! Spin_glass       = .True.  !Add the sigma parameters as sig=.. in Parameters_Jij Block' 
    endif

    write(inp2mcin_unit,'(2x,a)') 'Begin Parameters_Jij'
    if (lenconfac.eq.1.0_dp) then
       if (paramconfac.ne.1.0_dp) then
          write(inp2mcin_unit,'(3x,a)') 'Ryd'
       endif
    else
       if (paramconfac.eq.1.0_dp) then
          write(inp2mcin_unit,'(3x,a)') 'Bohr'
       else
          write(inp2mcin_unit,'(3x,a)') 'Ryd, Bohr'
       endif
    endif

    do nsp=1,num_species
       do ndnn = 1, shells(1,nsp)
          do nsp1=1,num_species
             nat=nat_pair(nsp,nsp1,ndnn)
             if( pair(nsp,nsp1,ndnn)) then
                if (spin_glass) then
                   write(inp2mcin_unit,'(3x,a3,i3,a4,i3,a4,i2,a22,a4,F11.8)') &
                   't1=',nsp,':t2=',nsp1,':sh=',ndnn,':Jij= ??????:sig=?????','!:d=',dnn(nsp,nat,ndnn)*lenconfac
                else
                   write(inp2mcin_unit,'(3x,a3,i3,a4,i3,a4,i2,a23,a4,F11.8)') &
                   't1=',nsp,':t2=',nsp1,':sh=',ndnn,':Jij= ??????!:sig=?????','!:d=',dnn(nsp,nat,ndnn)*lenconfac
                endif
             endif
          enddo
       enddo
    enddo
    write(inp2mcin_unit,'(2x,a)') 'End Parameters_Jij'

    if (have_biquad) then
       pair = .false.
       nat_pair = 0
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             do ndnn = 1, shells(2,nsp)
                do nnx=1 , multi(nsp,nat,ndnn)
                   nsp1 = nsplist(nsp,nat,ndnn,nnx)
                   pair(nsp,nsp1,ndnn) = .true.
                   nat_pair(nsp,nsp1,ndnn) = nat
                enddo
             enddo
          enddo
       enddo

       write(inp2mcin_unit,'(2x,a)') ''
       write(inp2mcin_unit,'(2x,a)') 'Begin Parameters_Bij'
       if (lenconfac.eq.1) then
          if (paramconfac.ne.1) then
             write(inp2mcin_unit,'(3x,a)') 'Ryd'
          endif
       else
          if (paramconfac.eq.1) then
             write(inp2mcin_unit,'(3x,a)') 'Bohr'
          else
             write(inp2mcin_unit,'(3x,a)') 'Ryd,Bohr'
          endif
       endif
       do nsp=1,num_species
          do ndnn = 1, shells(2,nsp)
             do nsp1=1,num_species
                nat=nat_pair(nsp,nsp1,ndnn)
                if( pair(nsp,nsp1,ndnn)) &
                   write(inp2mcin_unit,'(3x,a3,i3,a4,i3,a4,i2,a12,a4,F11.8)') &
                  't1=',nsp,':t2=',nsp1,':sh=',ndnn,':Bij= ??????','!:d=',dnn(nsp,nat,ndnn)*lenconfac
             enddo
          enddo
       enddo
       write(inp2mcin_unit,'(2x,a)') 'End Parameters_Bij'
    endif

    if (have_dm) then
       pair = .false.
       nat_pair = 0
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             do ndnn = 1, shells(3,nsp)
                do nnx=1 , multi(nsp,nat,ndnn)
                   nsp1 = nsplist(nsp,nat,ndnn,nnx)
                   pair(nsp,nsp1,ndnn) = .true.
                   nat_pair(nsp,nsp1,ndnn) = nat
                enddo
             enddo
          enddo
       enddo

       write(inp2mcin_unit,'(2x,a)') ''
       write(inp2mcin_unit,'(2x,a)') 'Begin Parameters_Dij'
       if (lenconfac.eq.1) then
          if (paramconfac.ne.1) then
             write(inp2mcin_unit,'(3x,a)') 'Ryd'
          endif
       else
          if (paramconfac.eq.1) then
             write(inp2mcin_unit,'(3x,a)') 'Bohr'
          else
             write(inp2mcin_unit,'(3x,a)') 'Ryd,Bohr'
          endif
       endif
       do nsp=1,num_species
          do ndnn = 1, shells(3,nsp)
             do nsp1=1,num_species
                nat=nat_pair(nsp,nsp1,ndnn)
                if( pair(nsp,nsp1,ndnn)) &
                   write(inp2mcin_unit,'(3x,a3,i3,a4,i3,a4,i2,a12,a4,F11.8)') &
                  't1=',nsp,':t2=',nsp1,':sh=',ndnn,':Dij= ??????','!:d=',dnn(nsp,nat,ndnn)*lenconfac
             enddo
          enddo
       enddo
       write(inp2mcin_unit,'(2x,a)') 'End Parameters_Dij'
    endif

    close(inp2mcin_unit)

    call io_stopwatch('inp2mcin: write',2)

    return

  end subroutine inp2mcin_write




  !==================================================================!
  subroutine mcin_write()
    !==================================================================!
    !                                                                  !
    ! Writes mcarlo.nat file                                           !
    !                                                                  ! 
    !                                                                  !
    !                                                                  !
    !===================================================================  
    use mc_io,        only : stdout,io_file_unit,seedname,io_stopwatch,&
                             io_error
    use mc_utility,   only : utility_cart_to_frac,utility_symmetric,&
                             utility_asymmetric_pairs

    implicit none

    integer           :: nsp,nat,nsp1,nat1,ndnn,vkpp(3),nnx
    integer           :: i,ierr,loop,num_asym,mcin_unit
    real(kind=dp)     :: neigh_cart(3),neigh_frac(3)
    integer,allocatable:: asym_matrix(:,:)


    call io_stopwatch('mcin: write',1)

    ! check the parameters_jij/bij/diij matrix is symmetric
    do ndnn=1,maxval(shells(1,:))
       call utility_symmetric(num_species,parameters_jij(:,:,ndnn),num_asym)
       if ( num_asym > 0 ) then
          allocate(asym_matrix(num_species,num_asym), stat=ierr )
          if (ierr/=0) call io_error('Error in allocating asym_matrix in mcin_write')
          call utility_asymmetric_pairs(num_species,num_asym,parameters_jij(:,:,ndnn),asym_matrix)
          write(stdout,'(2x,a,i8,a)') 'parameters_jij have',num_asym,' asymmetric pairs'
          do nsp=1,num_species
             do loop=1,num_asym
                if (asym_matrix(nsp,loop) > 0 ) then
                   nsp1=asym_matrix(nsp,loop)
                   write(stdout,'(2x,a,i3,a,i3,a,5x,a,i3,a,i3,a)') '(',nsp,',',nsp1,')','(',nsp1,',',nsp,')'
                endif
             enddo
          enddo
          deallocate(asym_matrix, stat=ierr )
          if (ierr/=0) call io_error('Error in deallocating asym_matrix in mcin_write')
          call io_error('Error: parameters_jij matrix is not symmetric. '&
                        //'See '//trim(seedname)//'.mcout for more details.')
       endif
    enddo

    if (spin_glass) then
       do ndnn=1,maxval(shells(1,:))
          call utility_symmetric(num_species,parameters_sigma(:,:,ndnn),num_asym)
          if ( num_asym > 0 ) then
             allocate(asym_matrix(num_species,num_asym), stat=ierr )
             if (ierr/=0) call io_error('Error in allocating asym_matrix in mcin_write')
             call utility_asymmetric_pairs(num_species,num_asym,parameters_sigma(:,:,ndnn),asym_matrix)
             write(stdout,'(2x,a,i8,a)') 'parameters_sigma have',num_asym,' asymmetric pairs'
             do nsp=1,num_species
                do loop=1,num_asym
                   if (asym_matrix(nsp,loop) > 0 ) then
                      nsp1=asym_matrix(nsp,loop)
                      write(stdout,'(2x,a,i3,a,i3,a,5x,a,i3,a,i3,a)') '(',nsp,',',nsp1,')','(',nsp1,',',nsp,')'
                   endif
                enddo
             enddo
             deallocate(asym_matrix, stat=ierr )
             if (ierr/=0) call io_error('Error in deallocating asym_matrix in mcin_write')
             call io_error('Error: parameters_sigma matrix is not symmetric. '&
                           //'See '//trim(seedname)//'.mcout for more details.')
          endif
       enddo
    endif

    if (have_biquad) then
       do ndnn=1,maxval(shells(2,:))
          call utility_symmetric(num_species,parameters_bij(:,:,ndnn),num_asym)
          if ( num_asym > 0 ) then
             allocate(asym_matrix(num_species,num_asym), stat=ierr )
             if (ierr/=0) call io_error('Error in allocating asym_matrix in mcin_write')
             call utility_asymmetric_pairs(num_species,num_asym,parameters_bij(:,:,ndnn),asym_matrix)
             write(stdout,'(2x,a,i8,a)') 'parameters_bij have',num_asym,'asymmetric pairs'
             do nsp=1,num_species
                do loop=1,num_asym
                   if (asym_matrix(nsp,loop) > 0 ) then
                      nsp1=asym_matrix(nsp,loop)
                      write(stdout,'(2x,a,i3,a,i3,a,5x,a,i3,a,i3,a)') '(',nsp,',',nsp1,')','(',nsp1,',',nsp,')'
                   endif
                enddo
             enddo
             deallocate(asym_matrix, stat=ierr )
             if (ierr/=0) call io_error('Error in deallocating asym_matrix in mcin_write')
             call io_error('Error: parameters_bij matrix is not symmetric. '&
                          //'See '//trim(seedname)//'.mcout for more details.')
          endif
       enddo
    endif

    if (have_dm) then
       do ndnn=1,maxval(shells(3,:))
          call utility_symmetric(num_species,parameters_dij(:,:,ndnn),num_asym)
          if ( num_asym > 0 ) then
             allocate(asym_matrix(num_species,num_asym), stat=ierr )
             if (ierr/=0) call io_error('Error in allocating asym_matrix in mcin_write')
             call utility_asymmetric_pairs(num_species,num_asym,parameters_dij(:,:,ndnn),asym_matrix)
             write(stdout,'(2x,a,i8,a)') 'parameters_dij have',num_asym,'asymmetric pairs'
             do nsp=1,num_species
                do loop=1,num_asym
                   if (asym_matrix(nsp,loop) > 0 ) then
                      nsp1=asym_matrix(nsp,loop)
                      write(stdout,'(2x,a,i3,a,i3,a,5x,a,i3,a,i3,a)') '(',nsp,',',nsp1,')','(',nsp1,',',nsp,')'
                   endif
                enddo
             enddo
             deallocate(asym_matrix, stat=ierr )
             if (ierr/=0) call io_error('Error in deallocating asym_matrix in mcin_write')
             call io_error('Error: parameters_dij matrix is not symmetric. '&
                           //'See '//trim(seedname)//'.mcout for more details.')
          endif
       enddo
    endif

    mcin_unit=io_file_unit()
    open(unit=mcin_unit,file=trim(seedname)//'.mcin',form='formatted')
    write(mcin_unit,'(1x,a)') ''
    write(mcin_unit,'(1x,a)') 'Begin Unit_Cell_Cart'
    if (lenconfac.ne.1) then
       write(mcin_unit,'(2x,a4)') 'Bohr'
    endif
    write(mcin_unit,'(2x,3f13.8)') (real_lattice(1,I)*lenconfac, i=1,3)
    write(mcin_unit,'(2x,3f13.8)') (real_lattice(2,I)*lenconfac, i=1,3)
    write(mcin_unit,'(2x,3f13.8)') (real_lattice(3,I)*lenconfac, i=1,3)
    write(mcin_unit,'(1x,a)') 'End Unit_Cell_Cart'
    write(mcin_unit,'(1x,a)') ''

    if (coordinate .eq. 'frac') then
       write(mcin_unit,'(1x,a)') 'Begin Atoms_Frac'
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             write(mcin_unit,'(2x,a3,3x,3f12.7,3x,f5.2)') atoms_symbol(nsp),&
                   (atoms_pos_frac(i,nsp,nat),i=1,3),mag_moments(nsp,nat)
          end do
       end do
       write(mcin_unit,'(1x,a)') 'End Atoms_Frac'
    else
       write(mcin_unit,'(1x,a)') 'Begin Atoms_Cart'
       if (lenconfac.ne.1) then
          write(mcin_unit,'(2x,a4)') 'Bohr'
       endif
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             write(mcin_unit,'(2x,a3,3x,3f12.7,3x,f5.2)') atoms_symbol(nsp),&
                   (atoms_pos_cart(i,nsp,nat)*lenconfac,i=1,3),mag_moments(nsp,nat)
          end do
       end do
       write(mcin_unit,'(1x,a)') 'End Atoms_Cart'
    endif
    write(mcin_unit,'(1x,a)') ''

    write(mcin_unit,'(1x,a)')     'tem_start          =   5'     
    write(mcin_unit,'(1x,a)')     'tem_end            =   5'     
    write(mcin_unit,'(1x,a)')     'tems_num           =   1'     
    write(mcin_unit,'(1x,a)')     '!! tems_mode          = man'     
    write(mcin_unit,'(1x,a)')     '!! tems               = 5.00 10.00 15.00 20.00'     
    write(mcin_unit,'(1x,a)')     ''

    write(mcin_unit,'(1x,a)') '!! Pt                 = .True.'     
    write(mcin_unit,'(1x,a)') '!! Pt_steps_swap      = 10'     
    write(mcin_unit,'(1x,a)') ''


    write(mcin_unit,'(1x,a,i12)') 'steps_warmup      = ',steps_warmup
    write(mcin_unit,'(1x,a,i12)') 'steps_mc          = ',steps_mc
    write(mcin_unit,'(1x,a,i12)') 'steps_measure     = ',steps_measure
    write(mcin_unit,'(1x,a)') ''

    write(mcin_unit,'(1x,a,a12)') 'initial_sconfig   = ',trim(initial_sconfig)
    write(mcin_unit,'(1x,a,a12)') 'mcarlo_mode       = ',trim(mcarlo_mode)
    write(mcin_unit,'(1x,a)') ''

    write(mcin_unit,'(1x,a,3i6)') 'supercell_size    = ',(supercell_size(i),i=1,3)
    write(mcin_unit,'(1x,a)') ''

    if (lspincorr) &
        write(mcin_unit,'(1x,a/)')'Spin_correlation  = .True.'                                                          
    if (energy_write) &
       write(mcin_unit,'(1x,a/)')  'Energy_write      = .True.'     

    if (lbinerror) then
       write(mcin_unit,'(1x,a)')  'Binning_error     = .True.'     
       write(mcin_unit,'(1x,a)')  'Binning_level     =  1-10 '     
       write(mcin_unit,'(1x,a)')   ''
    endif

    if (lstaggered) then
       write(mcin_unit,'(1x,a)')       'Staggered_m       = .True.'     
       write(mcin_unit,'(1x,a,100i2)') 'Staggered_m_coeff = ',(staggered_coeff(i),i=1,num_atoms)
       write(mcin_unit,'(1x,a)') ''
    endif

    if (lenconfac  .ne.1.0_dp) then
       write(mcin_unit,'(1x,a)')  'Length_unit       = Bohr'
       write(mcin_unit,'(1x,a)') ''
    !else
    !   write(mcin_unit,'(1x,a)')  '!! Length_unit         = Bohr'
    endif
    if (paramconfac.ne.1.0_dp) then
       write(mcin_unit,'(1x,a)')  'Parameter_unit    = Ryd'
       write(mcin_unit,'(1x,a)') ''
    !else
    !   write(mcin_unit,'(1x,a)')  '!! Parameter_unit      = Ryd'
    endif
    if (coordinate .ne. 'frac') then
       write(mcin_unit,'(1x,a)') 'Coordinate        = Cart'
       write(mcin_unit,'(1x,a)') ''
    !else
    !   write(mcin_unit,'(1x,a)')  '!! Coordinate          = Cart'
    endif

    if (lsfactor) then
       write(mcin_unit,'(1x,a)') 'Sfactor               = .True.'     
       write(mcin_unit,'(1x,a)') 'Sfactor_corner        = 0.000000      0.000000     0.000000'     
       write(mcin_unit,'(1x,a)') 'Sfactor_q1            = 1.000000      0.000000     0.000000'     
       write(mcin_unit,'(1x,a)') 'Sfactor_q2            = 0.000000      0.000000     1.000000'     
       write(mcin_unit,'(1x,a)') 'Sfactor_2Dqmesh       = 50     50     '     
       write(mcin_unit,'(1x,a)') 'Sfactor_steps_measure = 10  '     
       write(mcin_unit,'(1x,a)') ''
    endif

    if (lorderparam) then
       write(mcin_unit,'(1x,a)') 'Order_parameter   = .True.'     
       if (coordinate .eq. 'cart') then
          write(mcin_unit,'(1x,a)') 'Begin Order_Parameter_Axes_Cart'
           if (lenconfac.ne.1) &
              write(mcin_unit,'(3x,a)') 'Bohr'
       else
          write(mcin_unit,'(1x,a)') 'Begin Order_Parameter_Axes_Frac'
       endif
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             write(mcin_unit,'(1x,a,3x,a3,3x,a)')'!!   ', atoms_symbol(nsp),' axes       axes        axes'
          enddo
       enddo
       if (coordinate .eq. 'cart') then
          write(mcin_unit,'(1x,a)') 'End Order_Parameter_Axes_Cart'
       else
          write(mcin_unit,'(1x,a)') 'End Order_Parameter_Axes_Frac'
       endif
       write(mcin_unit,'(1x,a)') ''
    endif

    write(mcin_unit,'(2x,a)')    '## Hamiltonian'
    if (spin_glass) then
       write(mcin_unit,'(1x,a)')                     'Spin_glass        = .True.'                                                              
    !else
    !   write(mcin_unit,'(1x,a)')  '!! Spin_glass       = .True.  !Add the sigma parameters as sig=.. in Jij_parameters Block' 
    endif
    if (lboundary_open) &
       write(mcin_unit,'(1x,a)')                     'Boundary          = Open  '                                                              
    if (have_biquad)    write(mcin_unit,'(1x,a)')    'Ham_bij           = .True.'
    if (have_dm)        write(mcin_unit,'(1x,a)')    'Ham_dij           = .True.'
    if (have_singleion) write(mcin_unit,'(1x,a)')    'Ham_singleion     = .True.'
    if (have_field)     write(mcin_unit,'(1x,a)')    'Ham_field         = .True.'

    if (have_singleion) then
       if (coordinate .eq. 'cart') then
          write(mcin_unit,'(1x,a)') 'Begin SingleIon_Axes_Cart'
          if (lenconfac.eq.1) then
             if (paramconfac.ne.1)  write(mcin_unit,'(3x,a)') 'Ryd'
          else
             if (paramconfac.eq.1) then
                write(mcin_unit,'(3x,a)') 'Bohr'
             else
                write(mcin_unit,'(3x,a)') 'Ryd, Bohr'
             endif
          endif
       else
          write(mcin_unit,'(1x,a)') 'Begin SingleIon_Axes_Frac'
          if (paramconfac.ne.1)  write(mcin_unit,'(3x,a)') 'Ryd'
       endif
       i=0
       do nsp=1,num_species
           do nat=1,atoms_species_num(nsp)
              i=i+1
              !write(mcin_unit,'(6x,a3,2x,a,f10.5)') atoms_symbol(nsp),' 0.000000      0.000000  ' &
                    !//'   0.000000  ',singleion_parameters(i)
              write(mcin_unit,'(6x,a3,2x,a)') atoms_symbol(nsp),' axes        axes        axes        Delta'        
           enddo
       enddo
       if (coordinate .eq. 'cart') then
          write(mcin_unit,'(1x,a)') 'End SingleIon_Axes_Cart'
       else
          write(mcin_unit,'(1x,a)') 'End SingleIon_Axes_Frac'
       endif
       write(mcin_unit,'(1x,a)') ''
    endif

    if (have_field) then
       if (coordinate .eq. 'frac') then
          write(mcin_unit,'(1x,a)') 'Begin Field_Axes_Frac'
       else
          write(mcin_unit,'(1x,a)') 'Begin Field_Axes_Cart'
          if (lenconfac .ne. 1.0_dp) write(mcin_unit,'(3x,a)') 'Bohr'
       endif
       i=0
       do nsp=1,num_species
           do nat=1,atoms_species_num(nsp)
              i=i+1
              !write(mcin_unit,'(6x,a3,2x,a,f10.5)') atoms_symbol(nsp),' 0.000000      0.000000   '&
              !      //'   0.000000  ',field_parameters(i)
              write(mcin_unit,'(6x,a3,2x,a)') atoms_symbol(nsp),' axes        axes        axes        field'
           enddo
       enddo
       if (coordinate .eq. 'frac') then
          write(mcin_unit,'(1x,a)') 'End Field_Axes_Frac'
       else
          write(mcin_unit,'(1x,a)') 'End Field_Axes_Cart'
       endif
       write(mcin_unit,'(1x,a)') ''
    endif

    write(mcin_unit,'(1x,a)') 'Begin Jij_parameters'
    if (coordinate .eq. 'cart') then    
       if (lenconfac.eq.1) then
          if (paramconfac.ne.1)  write(mcin_unit,'(3x,a)') 'Ryd'
       else
          if (paramconfac.eq.1) then
             write(mcin_unit,'(2x,a)') 'Bohr'
          else
             write(mcin_unit,'(2x,a)') 'Ryd, Bohr'
          endif
       endif
    else
       if (paramconfac.ne.1)  write(mcin_unit,'(2x,a)') 'Ryd'
    endif
    do nsp=1,num_species
       do nat=1,atoms_species_num(nsp)
          do ndnn = 1, shells(1,nsp)
             do nnx=1 , multi(nsp,nat,ndnn)
                nat1 =  natlist(nsp,nat,ndnn,nnx) 
                nsp1 =  nsplist(nsp,nat,ndnn,nnx)
                if (parameters_jij(nsp,nsp1,ndnn).ne.0.0_dp) then
                   vkpp(1)  =  ncell(1,nsp,nat,ndnn,nnx)
                   vkpp(2)  =  ncell(2,nsp,nat,ndnn,nnx)
                   vkpp(3)  =  ncell(3,nsp,nat,ndnn,nnx)
                   neigh_cart(:)=matmul(vkpp(:),real_lattice)+atoms_pos_cart(:,nsp1,nat1)
                   if (coordinate .eq. 'frac')then    
                      call  utility_cart_to_frac(neigh_cart(:),neigh_frac(:),recip_lattice)
                      if (spin_glass) then
                         if (lspincorr) then
                            write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a5,F12.8,a4,i3,a4,i3,a4,i3)') 'f1=',&
                                  atoms_pos_frac(1,nsp,nat),',',atoms_pos_frac(2,nsp,nat),',',atoms_pos_frac(3,nsp,nat),&
                                  ':','f2=',neigh_frac(1),',',neigh_frac(2),',',neigh_frac(3),':','jij=',&
                                  parameters_jij(nsp,nsp1,ndnn)*paramconfac,':sig=',parameters_sigma(nsp,nsp1,ndnn)*paramconfac,&
                                  ':sh=',ndnn,'!t1=',nsp,':t2=',nsp1
                         else
                            write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a5,F12.8,a5,i3,a4,i3,a4,i3)') 'f1=',&
                                  atoms_pos_frac(1,nsp,nat),',',atoms_pos_frac(2,nsp,nat),',',atoms_pos_frac(3,nsp,nat),&
                                  ':','f2=',neigh_frac(1),',',neigh_frac(2),',',neigh_frac(3),':','jij=',&
                                  parameters_jij(nsp,nsp1,ndnn)*paramconfac,':sig=',parameters_sigma(nsp,nsp1,ndnn)*paramconfac,&
                                  '!:sh=',ndnn,'!t1=',nsp,':t2=',nsp1
                         endif
                      else
                         if (lspincorr) then
                            write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a4,i3,a5,i3,a4,i3)') 'f1=',&
                                  atoms_pos_frac(1,nsp,nat),',',atoms_pos_frac(2,nsp,nat),',',atoms_pos_frac(3,nsp,nat),&
                                  ':','f2=',neigh_frac(1),',',neigh_frac(2),',',neigh_frac(3),':','jij=',&
                                  parameters_jij(nsp,nsp1,ndnn)*paramconfac,':sh=',ndnn,'!:t1=',nsp,':t2=',nsp1
                         else
                            write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a5,i3,a5,i3,a4,i3)') 'f1=',&
                                  atoms_pos_frac(1,nsp,nat),',',atoms_pos_frac(2,nsp,nat),',',atoms_pos_frac(3,nsp,nat),&
                                  ':','f2=',neigh_frac(1),',',neigh_frac(2),',',neigh_frac(3),':','jij=',&
                                  parameters_jij(nsp,nsp1,ndnn)*paramconfac,'!:sh=',ndnn,'!:t1=',nsp,':t2=',nsp1
                         endif
                      endif
                   else
                      if (spin_glass) then
                         if (lspincorr) then
                            write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a5,F12.8,a4,i3,a4,i3,a4,i3)') 'c1=',&
                              atoms_pos_cart(1,nsp,nat)*lenconfac,',',atoms_pos_cart(2,nsp,nat)*lenconfac,',',&
                              atoms_pos_cart(3,nsp,nat)*lenconfac,':','c2=',neigh_cart(1)*lenconfac,',',&
                              neigh_cart(2)*lenconfac,',',neigh_cart(3)*lenconfac,':','jij=',&
                              parameters_jij(nsp,nsp1,ndnn)*paramconfac,':sig=',&
                              parameters_sigma(nsp,nsp1,ndnn)*paramconfac,':sh=',ndnn,'!t1=',nsp,':t2=',nsp1
                         else
                            write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a5,F12.8,a5,i3,a4,i3,a4,i3)') 'c1=',&
                              atoms_pos_cart(1,nsp,nat)*lenconfac,',',atoms_pos_cart(2,nsp,nat)*lenconfac,',',&
                              atoms_pos_cart(3,nsp,nat)*lenconfac,':','c2=',neigh_cart(1)*lenconfac,',',&
                              neigh_cart(2)*lenconfac,',',neigh_cart(3)*lenconfac,':','jij=',&
                              parameters_jij(nsp,nsp1,ndnn)*paramconfac,':sig=',&
                              parameters_sigma(nsp,nsp1,ndnn)*paramconfac,'!:sh=',ndnn,'!t1=',nsp,':t2=',nsp1
                         endif
                      else
                         if (lspincorr) then
                            write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a4,i3,a5,i3,a4,i3)') 'c1=',&
                              atoms_pos_cart(1,nsp,nat)*lenconfac,',',atoms_pos_cart(2,nsp,nat)*lenconfac,',',&
                              atoms_pos_cart(3,nsp,nat)*lenconfac,':','c2=',neigh_cart(1)*lenconfac,',',&
                              neigh_cart(2)*lenconfac,',',neigh_cart(3)*lenconfac,':','jij=',&
                              parameters_jij(nsp,nsp1,ndnn)*paramconfac,':sh=',ndnn,'!:t1=',nsp,':t2=',nsp1
                         else
                            write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a5,i3,a5,i3,a4,i3)') 'c1=',&
                              atoms_pos_cart(1,nsp,nat)*lenconfac,',',atoms_pos_cart(2,nsp,nat)*lenconfac,',',&
                              atoms_pos_cart(3,nsp,nat)*lenconfac,':','c2=',neigh_cart(1)*lenconfac,',',&
                              neigh_cart(2)*lenconfac,',',neigh_cart(3)*lenconfac,':','jij=',&
                              parameters_jij(nsp,nsp1,ndnn)*paramconfac,'!:sh=',ndnn,'!:t1=',nsp,':t2=',nsp1
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    enddo

    write(mcin_unit,'(1x,a)') 'End Jij_parameters'

    if (have_biquad) then
       write(mcin_unit,'(1x,a)') ''
       write(mcin_unit,'(1x,a)') 'Begin Bij_parameters'
       if (coordinate .eq. 'cart')then
          if (lenconfac.eq.1) then
             if (paramconfac.ne.1)  write(mcin_unit,'(2x,a)') 'Ryd'
          else
             if (paramconfac.eq.1) then
                write(mcin_unit,'(2x,a)') 'Bohr'
             else
                write(mcin_unit,'(2x,a)') 'Ryd, Bohr'
             endif
          endif
       else
          if (paramconfac.ne.1)  write(mcin_unit,'(2x,a)') 'Ryd'
       endif

       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             do ndnn = 1, shells(2,nsp)
                do nnx=1 , multi(nsp,nat,ndnn)
                   nat1 =  natlist(nsp,nat,ndnn,nnx)
                   nsp1 =  nsplist(nsp,nat,ndnn,nnx)
                   if (parameters_bij(nsp,nsp1,ndnn).ne.0.0_dp) then
                      vkpp(1)  =  ncell(1,nsp,nat,ndnn,nnx)
                      vkpp(2)  =  ncell(2,nsp,nat,ndnn,nnx)
                      vkpp(3)  =  ncell(3,nsp,nat,ndnn,nnx)
                      neigh_cart(:)=matmul(vkpp(:),real_lattice)+atoms_pos_cart(:,nsp1,nat1)
                      if (coordinate .eq. 'frac')then
                         call   utility_cart_to_frac(neigh_cart(:),neigh_frac(:),recip_lattice)
                         write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a5,i3,a4,i3)') 'f1=',&
                               atoms_pos_frac(1,nsp,nat),',',atoms_pos_frac(2,nsp,nat),',',atoms_pos_frac(3,nsp,nat),&
                               ':','f2=',neigh_frac(1),',',neigh_frac(2),',',neigh_frac(3),':','bij=',&
                               parameters_bij(nsp,nsp1,ndnn)*paramconfac,' !t1=',nsp,':t2=',nsp1
                      else
                         write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a5,i3,a4,i3)') 'c1=',&
                               atoms_pos_cart(1,nsp,nat)*lenconfac,',',atoms_pos_cart(2,nsp,nat)*lenconfac,',',&
                               atoms_pos_cart(3,nsp,nat)*lenconfac,':','c2=',neigh_cart(1)*lenconfac,',',&
                               neigh_cart(2)*lenconfac,',',neigh_cart(3)*lenconfac,':','bij=',&
                               parameters_bij(nsp,nsp1,ndnn)*paramconfac,' !t1=',nsp,':t2=',nsp1
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo

       write(mcin_unit,'(1x,a)') 'End Bij_parameters'

    endif

    if (have_dm) then
       write(mcin_unit,'(1x,a)') ''
       write(mcin_unit,'(1x,a)') 'Begin Dij_parameters'
       if (coordinate .eq. 'cart')then
          if (lenconfac.eq.1) then
             if (paramconfac.ne.1)  write(mcin_unit,'(2x,a)') 'Ryd'
          else
             if (paramconfac.eq.1) then
                write(mcin_unit,'(2x,a)') 'Bohr'
             else
                write(mcin_unit,'(2x,a)') 'Ryd, Bohr'
             endif
          endif
       else
          if (paramconfac.ne.1)  write(mcin_unit,'(2x,a)') 'Ryd'
       endif

       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             do ndnn = 1, shells(3,nsp)
                do nnx=1 , multi(nsp,nat,ndnn)
                   nat1 =  natlist(nsp,nat,ndnn,nnx)
                   nsp1 =  nsplist(nsp,nat,ndnn,nnx)
                   if (parameters_dij(nsp,nsp1,ndnn).ne.0.0_dp) then
                      vkpp(1)  =  ncell(1,nsp,nat,ndnn,nnx)
                      vkpp(2)  =  ncell(2,nsp,nat,ndnn,nnx)
                      vkpp(3)  =  ncell(3,nsp,nat,ndnn,nnx)
                      neigh_cart(:)=matmul(vkpp(:),real_lattice)+atoms_pos_cart(:,nsp1,nat1)
                      if (coordinate .eq. 'frac') then
                         call   utility_cart_to_frac(neigh_cart(:),neigh_frac(:),recip_lattice)
                         write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a5,i3,a4,i3)') 'f1=',&
                               atoms_pos_frac(1,nsp,nat),',',atoms_pos_frac(2,nsp,nat),',',atoms_pos_frac(3,nsp,nat),&
                               ':','f2=',neigh_frac(1),',',neigh_frac(2),',',neigh_frac(3),':','dij=',&
                               parameters_dij(nsp,nsp1,ndnn)*paramconfac,' !t1=',nsp,':t2=',nsp1
                      else
                         write(mcin_unit,'(2x,a3,3(F12.6,a),a3,3(F12.6,a),a4,F12.8,a5,i3,a4,i3)') 'c1=',&
                               atoms_pos_cart(1,nsp,nat)*lenconfac,',',atoms_pos_cart(2,nsp,nat)*lenconfac,',',&
                               atoms_pos_cart(3,nsp,nat)*lenconfac,':','c2=',neigh_cart(1)*lenconfac,',',&
                               neigh_cart(2)*lenconfac,',',neigh_cart(3)*lenconfac,':','dij=',&
                               parameters_dij(nsp,nsp1,ndnn)*paramconfac,' !t1=',nsp,':t2=',nsp1
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo

       write(mcin_unit,'(1x,a)') 'End Dij_parameters'

       write(mcin_unit,'(1x,a)') ''
       if (coordinate .eq. 'frac') then
          write(mcin_unit,'(1x,a)') 'Begin Dij_Vectors_frac'
       else
          write(mcin_unit,'(1x,a)') 'Begin Dij_Vectors_Cart'
          if (lenconfac .ne. 1.0_dp)  write(mcin_unit,'(3x,a4)') 'Bohr'
       endif
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             do ndnn = 1, shells(3,nsp)
                do nnx=1 , multi(nsp,nat,ndnn)
                   write(mcin_unit,'(1x,a)') '  0.000000      0.000000     0.000000  '
                enddo
             enddo
          enddo
       enddo
       if (coordinate .eq. 'frac') then
          write(mcin_unit,'(1x,a)') 'End Dij_Vectors_frac'
       else
          write(mcin_unit,'(1x,a)') 'End Dij_Vectors_Cart'
       endif



    endif


    close(mcin_unit)

    call io_stopwatch('mcin: write',2)

    return

  end subroutine mcin_write

  !==================================================================!
  subroutine xsf_write()
    !==================================================================!
    !                                                                  !
    ! Writes mcarlo.nat file                                           !
    !                                                                  ! 
    !                                                                  !
    !                                                                  !
    !===================================================================  
    use mc_io,        only : io_file_unit,seedname,io_stopwatch

    implicit none

    integer           :: nsp,nat,i
    integer           :: xsf_unit

    call io_stopwatch('xsf: write',1)

    xsf_unit=io_file_unit()
    open(unit=xsf_unit,file=trim(seedname)//'.xsf',form='formatted')

    write(xsf_unit,'("CRYSTAL")')
    write(xsf_unit,'("PRIMVEC")')
    write(xsf_unit,'(3f12.7)') real_lattice(1,1),real_lattice(1,2),real_lattice(1,3)
    write(xsf_unit,'(3f12.7)') real_lattice(2,1),real_lattice(2,2),real_lattice(2,3)
    write(xsf_unit,'(3f12.7)') real_lattice(3,1),real_lattice(3,2),real_lattice(3,3)
    write(xsf_unit,'("CONVVEC")')
    write(xsf_unit,'(3f12.7)') real_lattice(1,1),real_lattice(1,2),real_lattice(1,3)
    write(xsf_unit,'(3f12.7)') real_lattice(2,1),real_lattice(2,2),real_lattice(2,3)
    write(xsf_unit,'(3f12.7)') real_lattice(3,1),real_lattice(3,2),real_lattice(3,3)
    write(xsf_unit,'("PRIMCOORD")')
    write(xsf_unit,'(i6,"  1")') num_atoms
    do nsp=1,num_species
       do nat=1,atoms_species_num(nsp)
          write(xsf_unit,'(a2,3x,3f12.7)') atoms_symbol(nsp),(atoms_pos_cart(i,nsp,nat),i=1,3)
       end do
    end do

    close(xsf_unit)

    call io_stopwatch('xsf: write',2)

    return

  end subroutine xsf_write

  !==================================================================!
  subroutine neighbors_dealloc()
    !==================================================================!
    !                                                                  !
    !                                                                  ! 
    !                                                                  !
    !===================================================================  
    use mc_io,   only : io_error
    implicit none
    integer :: ierr


    ! Deallocate integer arrays that are public

    deallocate(ncell, stat=ierr )
    if (ierr/=0) call io_error('Error in deallocating ncell in neighbors_dealloc')
    deallocate(nsplist, stat=ierr )
    if (ierr/=0) call io_error('Error in deallocating nsplist in neighbors_dealloc')
    deallocate(natlist, stat=ierr )
    if (ierr/=0) call io_error('Error in deallocating natlist in neighbors_dealloc')

    return

  end subroutine neighbors_dealloc

  !==================================================================!
  subroutine neighbors_supercell_sort
    !==================================================================!
    !                                                                  !
    ! We look for atoms  neighbours in a large supercell of            !
    ! unit cells. Done sequentially this is very slow.                 !
    ! Here we order the cells by the distance from the origin          !
    ! Doing the search in this order gives a dramatic speed up         !
    !                                                                  !
    !==================================================================!  
    use mc_io,   only : io_stopwatch
    implicit none

    integer       :: counter,l,m,n,loop
    integer       :: lmn_cp( 3,(2*nsupcell+1)**3),indx(1)
    real(kind=dp) :: pos(3)
    real(kind=dp) :: dist((2*nsupcell+1)**3)
    real(kind=dp) :: dist_cp((2*nsupcell+1)**3)

    call io_stopwatch('neighbors: supercell_sort',1)

    counter=1
    lmn(:,counter)=0
    dist(counter)=0.0_dp
    do l = -nsupcell,  nsupcell 
       do m = -nsupcell,   nsupcell
          do n = -nsupcell,  nsupcell
             if(l==0 .and. m==0 .and. n==0) cycle
             counter=counter+1
             lmn(1,counter)=l;lmn(2,counter)=m;lmn(3,counter)=n
             pos=matmul(lmn(:,counter),real_lattice)
             dist(counter)=sqrt(dot_product(pos,pos))
          end do
       end do
    end do

    do loop=(2*nsupcell+1)**3,1,-1
       indx=internal_maxloc(dist)
       dist_cp(loop)=dist(indx(1))
       lmn_cp(:,loop)=lmn(:,indx(1))
       dist(indx(1))=-1.0_dp
    end do

    lmn=lmn_cp
    dist=dist_cp

    call io_stopwatch('neighbors: supercell_sort',2)

  end subroutine neighbors_supercell_sort


    !=========================================================================!
     function internal_maxloc(dist)
    !=========================================================================!
    !                                                                         !
    !  A predictable maxloc.                                                  !
    !                                                                         !
    !=========================================================================!

    use mc_constants, only : eps8   
    implicit none


    real(kind=dp), intent(in)  :: dist((2*nsupcell+1)**3)
    integer                    :: internal_maxloc
    integer                    :: guess(1),loop,counter
    integer                    :: list((2*nsupcell+1)**3)
    
    list=0
    counter=1
    
    guess=maxloc(dist)
    list(1)=guess(1)
    ! look for any degenerate values
    do loop=1,(2*nsupcell+1)**3
       if (loop==guess(1)) cycle
       if ( abs(dist(loop)-dist(guess(1))) < eps8 ) then
          counter=counter+1
          list(counter)=loop
       endif
    end do
    ! and always return the lowest index
    internal_maxloc=minval(list(1:counter))


    end function internal_maxloc



  end module mc_neighbors



