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

module mc_common

!==============================================================================
! This contains the common variables  to set up a Monte Carlo or parallel tempering simulations
!==============================================================================

  ! Should we remove this 'use mc_comms' and invoke in individual routines 
  ! when needed?
  !
  use mc_comms
  use mc_constants, only  : dp

  implicit none

  ! This 'save' statement could probably be ommited, since this module 
  ! is USEd by the main program 'wannier_parint'
  !
  save
 
  ! Default accessibility is PUBLIC
  !

  contains
  
  !===========================================================!
  subroutine mcint_param_dist
  !===========================================================!
  !                                                           !
  ! distribute the parameters across processors               !
  ! NOTE: we only send the ones mc or pt, not all in inputfile!
  !                                                           !
  !===========================================================!

    use mc_io,         only : io_error
    use mc_parameters

    integer :: ierr,max_sites
    call comms_bcast(real_lattice(1,1),9)
    call comms_bcast(recip_lattice(1,1),9)
    call comms_bcast(num_total_atoms,1)
    call comms_bcast(num_atoms,1)
    call comms_bcast(num_species,1) 
    call comms_bcast(num_supercell,1) 
    call comms_bcast(supercell_size(1),3) 
    call comms_bcast(energy_write,1) 
    call comms_bcast(energy_num_print,1) 
    call comms_bcast(initial_sconfig,len(initial_sconfig))
    call comms_bcast(mcarlo_mode,len(mcarlo_mode))
    call comms_bcast(steps_warmup,1)
    call comms_bcast(steps_measure,1)
    call comms_bcast(steps_mc,1)
    call comms_bcast(tilt_angles_max,1)
    call comms_bcast(pt,1)
    call comms_bcast(tems_num,1)
    if (pt) then 
       call comms_bcast(pt_steps_swap,1)
       call comms_bcast(pt_print_swap,1)
    endif
    call comms_bcast(have_biquad,1)
    call comms_bcast(have_dm,1)
    call comms_bcast(have_singleion,1)
    call comms_bcast(have_field,1)
    call comms_bcast(lorderparam,1)
    call comms_bcast(lstaggered,1)
    call comms_bcast(lbinerror,1)
    call comms_bcast(num_binning_level,1)
    call comms_bcast(lspincorr,1)
    call comms_bcast(jij_shell_max,1)
    call comms_bcast(lseed,1) 
    call comms_bcast(lsfactor,1)
    call comms_bcast(lsfactor_polar,1)
    call comms_bcast(sfactor_steps_measure,1)
    call comms_bcast(sfactor_2dqmesh(1),2)
    call comms_bcast(sfactor_nqpts,1)
    call comms_bcast(sfactor_corner(1),3)
    call comms_bcast(sfactor_q1(1),3)
    call comms_bcast(sfactor_q2(1),3)
    call comms_bcast(sfactor_polar(1),3)

    ! These variables are different from the ones above in that they are 
    ! allocatable, and in param_read they were allocated on the root node only
    !
    if (.not.on_root) then
      allocate(binning_level(num_binning_level),stat=ierr)
      if (ierr/=0) call io_error('Error allocating binning_level in mcint_param_dist')
      allocate(seeds(num_nodes),stat=ierr)
      if (ierr/=0) call io_error('Error allocating seeds in mcint_param_dist')
      allocate(atoms_species_num(num_species),stat=ierr)
      if (ierr/=0) call io_error('Error allocating atoms_species_num in mcint_param_dist')
      if(lorderparam) then
        allocate(orderparam_axes_cart(3,num_atoms),stat=ierr)
        if (ierr/=0)&
            call io_error('Error allocating orderparam_axes_cart in mcint_param_dist')
      endif
      if (lstaggered) then
         allocate(staggered_coeff(num_atoms),stat=ierr)
         if (ierr/=0)&
            call io_error('Error allocating staggered_coeff in mcint_param_dist')
      end if

      allocate(tems(tems_num),stat=ierr)
      if (ierr/=0) call io_error('error allocating tems in mcint_data_dist')
    endif ! .not.on_root

    call comms_bcast(tems(1),tems_num)
    if (lbinerror) call comms_bcast(binning_level(1),num_binning_level)
    call comms_bcast(seeds(1),num_nodes) 
    call comms_bcast(atoms_species_num(1),num_species) 
    max_sites = maxval(atoms_species_num)
    if(lorderparam)  call comms_bcast(orderparam_axes_cart(1,1),3*num_atoms)
    if(lstaggered)   call comms_bcast(staggered_coeff(1),num_atoms)         

    if (.not.on_root) then
      allocate(atoms_pos_frac(3,num_species,max_sites),stat=ierr)
      if (ierr/=0) call io_error('Error allocating atoms_pos_frac in param_get_atoms_moments')
    endif
    call comms_bcast(atoms_pos_frac(1,1,1),3*num_species*max_sites) 


  end subroutine mcint_param_dist


  !===========================================================!
  subroutine mcint_data_dist
  !===========================================================!
  !                                                           !
  ! Distribute the (jij-bij-dij-sion-field) matrix           !
  !                                                           !
  !===========================================================!

    use mc_io,         only : io_error
    use mc_parameters, only : have_biquad,have_dm,num_total_atoms,&
                              have_singleion,singleion_parameters,&
                              singleion_axes_cart,num_atoms,jij_shell_max,&
                              have_field,field_axes_cart,&
                              field_parameters,num_species,lspincorr
    use mc_jij,        only : jij_matrix,mag_moments_matrix,&
                              jij_num_nbors_matrix,jij_nbors_matrix,&
                              shell_num_nbors_matrix,shell_nbors_matrix,&
                              shell_num_nbors_kind,kind_atoms

    use mc_bij,        only : bij_matrix,bij_num_nbors_matrix,bij_nbors_matrix
    use mc_dij,        only : dij_matrix,dij_num_nbors_matrix,dij_nbors_matrix

    implicit none

    integer :: ierr

    ! allocate on all nodes
    if (.not.on_root) then

      if (.not.allocated(jij_matrix)) then
        allocate(jij_num_nbors_matrix(num_total_atoms),stat=ierr)
        if (ierr/=0)&
           call io_error('error allocating jij_num_nbors_matrix in mcint_data_dist')    
      endif

      if (lspincorr .and. .not.allocated(shell_num_nbors_matrix)) then
        allocate(shell_num_nbors_matrix(jij_shell_max,num_total_atoms),stat=ierr)
        if (ierr/=0)&
           call io_error('error allocating jij_num_nbors_matrix in mcint_data_dist')    
      endif

      if (have_biquad .and. .not.allocated(bij_matrix)) then
        allocate(bij_num_nbors_matrix(num_total_atoms),stat=ierr)
        if (ierr/=0)&
           call io_error('error allocating bij_num_nbors_matrix in mcint_data_dist')    
      endif

      if (have_dm .and. .not.allocated(dij_matrix)) then
        allocate(dij_num_nbors_matrix(num_total_atoms),stat=ierr)
        if (ierr/=0)&
           call io_error('error allocating dij_num_nbors_matrix in mcint_data_dist')    
      endif
 
    endif  ! .not. on_root

    call comms_bcast(jij_num_nbors_matrix(1),num_total_atoms)
    if (lspincorr)    call comms_bcast(shell_num_nbors_matrix(1,1),jij_shell_max*num_total_atoms)
    if (have_biquad)  call comms_bcast(bij_num_nbors_matrix(1),num_total_atoms)
    if (have_dm)      call comms_bcast(dij_num_nbors_matrix(1),num_total_atoms)

    if (.not.on_root) then
       if (.not.allocated(jij_matrix)) then
          allocate(jij_matrix(num_total_atoms,num_total_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating jij_matrix in mcint_data_dist')    
          allocate(jij_nbors_matrix(num_total_atoms,num_total_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating jij_nbors_matrix in mcint_data_dist')    
          allocate(mag_moments_matrix(num_total_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating mag_moments_matrix in mcint_data_dist')    
       endif

       if (lspincorr .and. .not.allocated(shell_num_nbors_kind)) then
          allocate(kind_atoms(num_total_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating kind_atoms in mcint_data_dist')    
          allocate(shell_num_nbors_kind(jij_shell_max,num_species,num_species),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating shell_num_nbors_kind in mcint_data_dist')    
          allocate(shell_nbors_matrix(jij_shell_max,num_total_atoms,maxval(shell_num_nbors_matrix)),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating jij_nbors_matrix in mcint_data_dist')    
       endif

       if (have_biquad .and. .not.allocated(bij_matrix)) then
          allocate(bij_matrix(num_total_atoms,num_total_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating bij_matrix in mcint_data_dist')    
          allocate(bij_nbors_matrix(num_total_atoms,num_total_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating bij_nbors_matrix in mcint_data_dist')    
       endif

       if (have_dm .and. .not.allocated(dij_matrix)) then
          allocate(dij_matrix(3,num_total_atoms,num_total_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating dij_matrix in mcint_data_dist')    
          allocate(dij_nbors_matrix(num_total_atoms,num_total_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating dij_nbors_matrix in mcint_data_dist')    
       endif
  
       if (have_singleion .and. .not.allocated(singleion_parameters))then
          allocate(singleion_parameters(num_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating singleion_parameters in mcint_data_dist')
          allocate(singleion_axes_cart(3,num_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating singleion_axes_cart in mcint_data_dist')
       endif

       if (have_field .and. .not.allocated(field_parameters))then
          allocate(field_parameters(num_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating field_parameters in mcint_data_dist')
          allocate(field_axes_cart(3,num_atoms),stat=ierr)
          if (ierr/=0)&
             call io_error('error allocating field_axes_cart in mcint_data_dist')
       endif
 
    endif ! .not. on_root


    call comms_bcast(jij_matrix(1,1),num_total_atoms*maxval(jij_num_nbors_matrix))
    call comms_bcast(jij_nbors_matrix(1,1),num_total_atoms*maxval(jij_num_nbors_matrix))
    call comms_bcast(mag_moments_matrix(1),num_total_atoms)

    if (lspincorr) then
       call comms_bcast(kind_atoms(1),num_total_atoms)
       call comms_bcast(shell_nbors_matrix(1,1,1),jij_shell_max*num_total_atoms&
                        *maxval(shell_num_nbors_matrix))
       call comms_bcast(shell_num_nbors_kind(1,1,1),jij_shell_max*num_species*num_species)
    endif

    if (have_biquad) then
       call comms_bcast(bij_matrix(1,1),num_total_atoms*maxval(bij_num_nbors_matrix))
       call comms_bcast(bij_nbors_matrix(1,1),num_total_atoms*maxval(bij_num_nbors_matrix))
    endif

    if (have_dm) then
       call comms_bcast(dij_matrix(1,1,1),3*num_total_atoms*maxval(dij_num_nbors_matrix))
       call comms_bcast(dij_nbors_matrix(1,1),num_total_atoms*maxval(dij_num_nbors_matrix))
    endif   

    if (have_singleion) then
       call comms_bcast(singleion_parameters(1),num_atoms)
       call comms_bcast(singleion_axes_cart(1,1),3*num_atoms)
    endif

    if (have_field) then
       call comms_bcast(field_parameters(1),num_atoms)
       call comms_bcast(field_axes_cart(1,1),3*num_atoms)
    endif

  end subroutine mcint_data_dist

!=======================================================================


end module mc_common
