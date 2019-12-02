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

module mc_pt

  use mc_constants,  only : dp,zero,pi,twopi
  use mc_parameters, only : &
       steps_warmup,steps_mc,steps_measure,&
       num_supercell,num_atoms,num_total_atoms,num_species,&
       atoms_symbol,initial_sconfig,tems,tems_num,&
       pt_steps_swap,pt_print_swap,state,energy_write,&
       energy_num_print,lorderparam,lstaggered,lsfactor,&
       sfactor_steps_measure,sfactor_nqpts,lbinerror,supercell_size,&
       num_binning_level,binning_level,lspincorr,jij_shell_max
  use mc_io,         only : io_error,stdout,io_stopwatch,&
                            io_file_unit,seedname,io_date 
  use mc_common
  use mc_comms
  use mc_get_quant
!  use stdtypes
  use mtprng,        only : mtprng_rand_real2

  implicit none

  private 
  public :: pt_main

  ! Constants to identify the nine components of a tensor when it is stored in packed form
  integer, parameter :: E  = 1
  integer, parameter :: E2 = 2
  integer, parameter :: E4 = 3
  integer, parameter :: M  = 4
  integer, parameter :: M2 = 5
  integer, parameter :: M4 = 6
  integer, parameter :: OP = 7
  integer, parameter :: OP2 = 8
  integer, parameter :: OP4 = 9
  integer            :: SM
  integer            :: SM2
  integer            :: SM4


contains 

  !============================================!
  subroutine pt_main()
    !============================================!

    implicit none

    character(len=9)                  :: cdate, ctime
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes-1) :: counts
    integer, dimension(0:num_nodes-1) :: displs
    integer                           :: ierr,mc_unit,orderparam_unit
    integer                           :: binerror_unit,pm_unit,sconfig_unit
    integer                           :: spincorr_unit,staggered_unit
    integer                           :: num_freq,pt_num_swap
    integer                           :: ndim_data,ndim_sum,ndim_pm
    integer,       allocatable        :: energy_unit(:)
    integer,       allocatable        :: map(:)
    integer(8),    allocatable        :: local_rejected(:),global_rejected(:)
    integer(8),    allocatable        :: num_swaps_tot(:,:)
    real(kind=dp)                     :: d_m
    real(kind=dp), allocatable        :: local_data(:,:),global_data(:,:)
    real(kind=dp), allocatable        :: local_spin_matrix(:,:,:),global_spin_matrix(:,:,:)
    real(kind=dp), allocatable        :: local_energy(:,:),global_energy(:,:)
    real(kind=dp), allocatable        :: atoms_supercell_pos_cart(:,:)
    real(kind=dp), allocatable        :: local_sfactor_matrix(:,:,:),global_sfactor_matrix(:,:,:)
    real(kind=dp), allocatable        :: local_spincorr(:,:,:,:),global_spincorr(:,:,:,:)
    real(kind=dp), allocatable        :: local_spincorr_abs(:,:,:,:),global_spincorr_abs(:,:,:,:)
    real(kind=dp), allocatable        :: local_spincorr_tot_abs(:,:,:,:),global_spincorr_tot_abs(:,:,:,:)
    real(kind=dp), allocatable        :: local_sum(:,:,:),local_binning(:,:,:),global_binning(:,:,:)
    real(kind=dp), allocatable        :: local_E(:,:),global_E(:,:)
    real(kind=dp), allocatable        :: m_array(:),local_pm(:,:,:),global_pm(:,:,:)
    real(kind=dp), allocatable        :: tmp_spin_matrix(:,:)

    ! Initial energy=local_energy(1:);;one configuration energy=energy_oc=local_energy(2,:);;final direct energy=local_energy(3,:)
    ! Energy_tot_ave=local_data(1,:);;energy2_tot_ave=local_data(2,:);;energy4_tot_ave=local_data(3,:)
    ! m_tot_ave=local_data(4,:);;m2_tot_ave=local_data(5,:);;m4_tot_ave=local_data(6,:)
    ! If we have onle order parameter:
    ! orderparam_ave=local_data(7,:);;orderparam2_ave=local_data(8,:);;orderparam4_ave=local_data(9,:)   
    ! If we have onle staggered magnetization:
    ! staggeredm_ave=local_data(7,:);;staggeredm2_ave=local_data(8,:);;staggeredm4_ave=local_data(9,:)   
    ! If we have both order parameter and staggered magnetization:
    ! orderparam_ave=local_data(7,:);;orderparam2_ave=local_data(8,:);;orderparam4_ave=local_data(9,:)   
    ! staggeredm_ave=local_data(10,:);;staggeredm2_ave=local_data(11,:);;staggeredm4_ave=local_data(12,:)  

    if(on_root) call io_stopwatch('pt_main',1)
    call comms_array_split(tems_num,counts,displs)

    call pt_setup()

    !! After setup spin matrix , I can overwrite the prefix_sconfig.dat if it's exist!!
    if(on_root)then   
       sconfig_unit = io_file_unit()
       open(sconfig_unit,file=trim(seedname)//'_sconfig.dat',status='unknown',form='formatted',action='write')

       call io_date(cdate,ctime)
       write(sconfig_unit,*) 'Created on '//cdate//' at '//ctime ! Date and time
       write(sconfig_unit,'(3I6)') supercell_size(1),supercell_size(2),supercell_size(3)
       write(sconfig_unit,'(2I6)') num_total_atoms,tems_num

       write(stdout,*) 
       write(stdout,'(1x,a)') '*-----------------------------------------------------------------------------------*'
       write(stdout,'(1x,a)') '|                            Parallel Tempering (Pt module)                         |'
       write(stdout,'(1x,a)') '*-----------------------------------------------------------------------------------*'

       mc_unit = io_file_unit()
       open(mc_unit,file=trim(seedname)//'_mc.dat',status='unknown',form='formatted',action='write')
       call io_date(cdate,ctime)
       write(mc_unit,*) '#written on '//cdate//' at '//ctime ! Date and time
       write(mc_unit,'(1x,a)') '#----------------------------------------------- MONTECARLO ----------'&
                                        //'-----------------------------------------#'
       write(mc_unit,'(1x,a)') '# Temp      Magnetization     Energy_ave          C_M              Sus'&
                                        //'              U_E              U_M       #'
       write(mc_unit,'(1x,a)') '#---------------------------------------------------------------------'&
                                        //'-----------------------------------------#'

       if(lspincorr) then
         spincorr_unit = io_file_unit()
         open(spincorr_unit,file=trim(seedname)//'_spincorr.dat',status='unknown',form='formatted',action='write')
         call io_date(cdate,ctime)
         write(spincorr_unit,*) '#written on '//cdate//' at '//ctime ! Date and time
         write(spincorr_unit,'(1x,a)')    '#---------------------------------------------------- MONTECARLO '&
                                        //'----------------------------------------------------#'
         write(spincorr_unit,'(1x,a)') '#  Temp    Shell   Atom1   Type1   Atom2   Type2   '&
                                        //'<\sum_{ij} Si.Sj/N>   <\sum_{ij}|Si.Sj|/N>   <|\sum_{ij}Si.Sj|/N> #'
         write(spincorr_unit,'(1x,a)') '#----------------------------------------------------------------'&
                                        //'----------------------------------------------------#'
       endif

       pm_unit = io_file_unit()
       open(pm_unit,file=trim(seedname)//'_pm.dat',status='unknown',form='formatted',action='write')
       call io_date(cdate,ctime)
       write(pm_unit,*) '#written on '//cdate//' at '//ctime ! Date and time
       if (lbinerror) then
         binerror_unit = io_file_unit()
         open(binerror_unit,file=trim(seedname)//'_binerror.dat',status='unknown',form='formatted',action='write')
         call io_date(cdate,ctime)
       endif
       if (lorderparam.and.lstaggered) then
          write(pm_unit,'(1x,a)') '#-------------------------------------- MONTECARLO ------'&
                                   //'------------------------------------#'
          write(pm_unit,'(1x,a)') '#       M or OP                 P(M)                     '&
                                   //'P(OP)               P(Staggered_m)  #'
          write(pm_unit,'(1x,a)') '#--------------------------------------------------------'&
                                   //'------------------------------------#'
          if (lbinerror) then
             write(binerror_unit,'(1x,a)') '#------------------------------------- MONTECARLO ---------------------'&
                                            //'--------------------#'
             write(binerror_unit,'(1x,a)') '#  Binning_level     Error(E)          Error(M)        Error(OP)       '&
                                            //'Error(Staggered_m)  #'
             write(binerror_unit,'(1x,a)') '#----------------------------------------------------------------------'&
                                            //'--------------------#'
          endif
       elseif (lorderparam) then
          write(pm_unit,'(1x,a)') '#---------------------------- MONTECARLO -----------------------------#'
          write(pm_unit,'(1x,a)') '#       M or OP                 P(M)                     P(OP)        #'
          write(pm_unit,'(1x,a)') '#---------------------------------------------------------------------#'
          if (lbinerror) then
             write(binerror_unit,'(1x,a)') '#--------------------------- MONTECARLO -----------------------------#'
             write(binerror_unit,'(1x,a)') '#  Binning_level     Error(E)          Error(M)         Error(OP)    #'
             write(binerror_unit,'(1x,a)') '#--------------------------------------------------------------------#'
          endif
       elseif (lstaggered) then
          write(pm_unit,'(1x,a)') '#--------------------------- MONTECARLO-----------------------------#'
          write(pm_unit,'(1x,a)') '# M or Staggered                P(M)                P(Staggered_m)  #'
          write(pm_unit,'(1x,a)') '#-------------------------------------------------------------------#'
          if (lbinerror) then
             write(binerror_unit,'(1x,a)') '#----------------------------- MONTECARLO -------------------------------#'
             write(binerror_unit,'(1x,a)') '#  Binning_level     Error(E)          Error(M)      Error(Staggered_m)  #'
             write(binerror_unit,'(1x,a)') '#------------------------------------------------------------------------#'
          endif
       else
          write(pm_unit,'(1x,a)') '#----------------- MONTECARLO -------------------#'
          write(pm_unit,'(1x,a)') '#        M                      P(M)             #'
          write(pm_unit,'(1x,a)') '#------------------------------------------------#'
          if (lbinerror) then
             write(binerror_unit,'(1x,a)') '#------------------- MONTECARLO -----------------#'
             write(binerror_unit,'(1x,a)') '#  Block_length      Error(E)          Error(M)  #'
             write(binerror_unit,'(1x,a)') '#------------------------------------------------#'
          endif
       endif

       if (lorderparam) then
         orderparam_unit = io_file_unit()
         open(orderparam_unit,file=trim(seedname)//'_op.dat',status='unknown',form='formatted',action='write')
         call io_date(cdate,ctime)
         write(orderparam_unit,*) '#written on '//cdate//' at '//ctime ! Date and time
         write(orderparam_unit,'(1x,a)') '#------------------------------- MONTECARLO -----------------------------------#'
         write(orderparam_unit,'(1x,a)') '#    Temp              OP                 Sus_OP                 U_OP          #'
         write(orderparam_unit,'(1x,a)') '#------------------------------------------------------------------------------#'
       endif ! lorderparam

       if (lstaggered) then
         staggered_unit = io_file_unit()
         open(staggered_unit,file=trim(seedname)//'_staggered.dat',status='unknown',form='formatted',action='write')
         call io_date(cdate,ctime)
         write(staggered_unit,*) '#written on '//cdate//' at '//ctime ! Date and time
         write(staggered_unit,'(1x,a)') '#------------------------------- MONTECARLO --------------------------------#'
         write(staggered_unit,'(1x,a)') '#    Temp         Staggered_m          Sus_Staggered_m       U_Staggered_m  #'
         write(staggered_unit,'(1x,a)') '#---------------------------------------------------------------------------#'
       endif ! lstaggered

       if (energy_write) then
          pt_num_swap = ceiling(real(steps_mc,dp)/real(pt_steps_swap,dp))
          call open_E_unit(pt_num_swap*pt_steps_swap,energy_unit)
       endif

    end if ! on_root

    call pt_mcarlo()

    ! Before ending, I deallocate memory

    call pt_dealloc()

    if (on_root)                        close(mc_unit)
    if (on_root)                        close(pm_unit)
    if (on_root)                        close(sconfig_unit)
    if (on_root .and. lbinerror)        close(binerror_unit)
    if (on_root .and. lspincorr)        close(spincorr_unit)
    if (on_root .and. lorderparam)      close(orderparam_unit)
    if (on_root .and. lstaggered)       close(staggered_unit)


    return

  contains

    !============================================!
    subroutine pt_setup() 
      !============================================!
 
      use mc_jij, only    : mag_moments_matrix
      use stdtypes
      use mtprng  

      implicit none

      integer                    :: loop,loop_counts
      integer                    :: loop_nodes,localidx,globalidx
      character(len=50)          :: filename 
      real(kind=dp)              :: r1,r2
      integer                    :: ierr,ifreq
      real(kind=dp)              :: m_min,m_max
!      integer,allocatable        :: seed_pt(:)


      if (lorderparam.and.lstaggered) then
         ndim_data = 12; ndim_pm = 3; ndim_sum = 6
         SM = 10; SM2 = 11; SM4 = 12
      elseif (lorderparam.or.lstaggered) then
         if (lstaggered) SM = 7; SM2 = 8; SM4 = 9
         ndim_data = 9 ; ndim_pm = 2; ndim_sum = 5
      else
         ndim_data = 6 ; ndim_pm = 1; ndim_sum = 4
      endif
  
      allocate(local_spin_matrix(3,num_total_atoms,counts(my_node_id)), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating local_spin_matrix in pt_setup')
      allocate(tmp_spin_matrix(3,num_total_atoms), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating tmp_spin_matrix in pt_setup')
      allocate(local_rejected(counts(my_node_id)), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating local_rejected in pt_setup')
      allocate(local_energy(3,counts(my_node_id)), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating local_energy in pt_setup')

      allocate(local_data(ndim_data,counts(my_node_id)), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating local_data in pt_setup')

      if (lbinerror) then
         allocate(local_sum(ndim_sum,num_binning_level,counts(my_node_id)), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating local_sum in pt_setup')
         allocate(local_binning(ndim_sum,num_binning_level,counts(my_node_id)), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating local_binning in pt_setup')
      endif

      if (lspincorr) then
        allocate(local_spincorr(jij_shell_max,num_species,num_species,counts(my_node_id)), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating local_spincorr in pt_setup')
        allocate(local_spincorr_abs(jij_shell_max,num_species,num_species,counts(my_node_id)), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating local_spincorr_abs in pt_setup')
        allocate(local_spincorr_tot_abs(jij_shell_max,num_species,num_species,counts(my_node_id)), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating local_spincorr_tot_abs in pt_setup')
      endif

      if (energy_write)then
        ! Since I have a loop over the pt_steps_swap, 
        ! I need to change energy_num_print according to the pt_steps_swap.
        energy_num_print = (energy_num_print/pt_steps_swap)*pt_steps_swap
        allocate(local_E(0:energy_num_print-1,0:counts(my_node_id)-1), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating local_E in pt_setup')
      endif

      if (lsfactor) then
        allocate(atoms_supercell_pos_cart(3,num_total_atoms), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating atoms_supercell_pos_cart in pt_setup')
        allocate(local_sfactor_matrix(3,sfactor_nqpts,counts(my_node_id)), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating local_sfactor_matrix in pt_setup')
        local_sfactor_matrix=0.0_dp
        call mcarlo_get_atoms_pos(atoms_supercell_pos_cart)
      endif


      allocate(local_spin_matrix(3,num_total_atoms,counts(my_node_id)), stat=ierr )

      if (on_root) then
        allocate(global_rejected(tems_num), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating global_rejected in pt_setup')
        allocate(global_spin_matrix(3,num_total_atoms,tems_num), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating global_spin_matrix in pt_setup')
        allocate(global_data(ndim_data,tems_num), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating global_data in pt_setup')
        allocate(global_energy(3,tems_num), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating global_energy in pt_setup')
        allocate(map(tems_num), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating map in pt_setup')
        allocate(num_swaps_tot(2,tems_num), stat=ierr )
        if (ierr/=0) call io_error('Error in allocating num_swaps_tot in pt_setup')
        if (lspincorr) then
          allocate(global_spincorr(jij_shell_max,num_species,num_species,tems_num), stat=ierr )
          if (ierr/=0) call io_error('Error in allocating global_spincorr in pt_setup')
          allocate(global_spincorr_abs(jij_shell_max,num_species,num_species,tems_num), stat=ierr )
          if (ierr/=0) call io_error('Error in allocating global_spincorr_abs in pt_setup')
          allocate(global_spincorr_tot_abs(jij_shell_max,num_species,num_species,tems_num), stat=ierr )
          if (ierr/=0) call io_error('Error in allocating global_spincorr_tot_abs in pt_setup')
        endif
        if (lbinerror) then
          allocate(global_binning(ndim_sum,num_binning_level,tems_num), stat=ierr )
          if (ierr/=0) call io_error('Error in allocating global_binning in pt_setup')
        endif
        if (lsfactor) then
          allocate(global_sfactor_matrix(3,sfactor_nqpts,tems_num), stat=ierr )
          if (ierr/=0) call io_error('Error in allocating global_sfactor_matrix in pt_setup')
          global_sfactor_matrix=0.0_dp
        endif
        if (energy_write) then
           allocate(energy_unit(tems_num), stat=ierr )
           if (ierr/=0) call io_error('Error in allocating energy_unit in pt_setup')
           allocate(global_E(0:energy_num_print-1,0:tems_num-1), stat=ierr )
           if (ierr/=0) call io_error('Error in allocating global_E in pt_main')
        endif
      else
         ! In principle, this should not be needed, because we use global_rejected,
         ! global_spin_matrix, global_energy, map, golbal_data, global_pm only on the root node. However, since all
         ! processors call comms_gatherv a few lines below, and one argument
         ! is global_rejected(1), some compilers complain.
         allocate(global_rejected(1),stat=ierr)
         if (ierr/=0) call io_error('Error in allocating global_rejected in pt_setup')
         allocate(global_spin_matrix(1,1,1),stat=ierr)
         if (ierr/=0) call io_error('Error in allocating global_spin_matrix in pt_setup')
         allocate(global_energy(1,1),stat=ierr)
         if (ierr/=0) call io_error('Error in allocating global_energy in pt_setup')
         allocate(map(1),stat=ierr)
         if (ierr/=0) call io_error('Error in allocating map in pt_setup')      
         allocate(num_swaps_tot(1,1),stat=ierr)
         if (ierr/=0) call io_error('Error in allocating num_swaps_tot in pt_setup')      
         allocate(global_data(1,1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_data in pt_setup')
         if (lspincorr) then 
           allocate(global_spincorr(1,1,1,1),stat=ierr)
           if (ierr/=0) call io_error('Error in allocating global_spincorr in pt_setup')
           allocate(global_spincorr_abs(1,1,1,1),stat=ierr)
           if (ierr/=0) call io_error('Error in allocating global_spincorr_abs in pt_setup')
           allocate(global_spincorr_tot_abs(1,1,1,1),stat=ierr)
           if (ierr/=0) call io_error('Error in allocating global_spincorr_tot_abs in pt_setup')
         endif
         if (lbinerror) then 
           allocate(global_binning(1,1,1),stat=ierr)
           if (ierr/=0) call io_error('Error in allocating global_binning pt_setup')      
         endif
         if (lsfactor) then
           allocate(global_sfactor_matrix(1,1,1), stat=ierr )
           if (ierr/=0) call io_error('Error in allocating global_sfactor_matrix in pt_setup')
           global_sfactor_matrix=0.0_dp
         endif
         if (energy_write) then
            allocate(energy_unit(1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating energy_unit in pt_setup')
           allocate(global_E(1,1), stat=ierr )
           if (ierr/=0) call io_error('Error in allocating global_E in pt_setup')
         endif
      end if ! on_root

      ! grid for distribution function
      m_min = zero; d_m = 0.01_dp
      m_max = zero
      do loop=1,num_atoms
         m_max = m_max+abs(mag_moments_matrix(loop))
      enddo
      m_min = -m_max
      num_freq = nint((m_max-m_min)/d_m)+1
      if (num_freq==1) num_freq = 2
      d_m = (m_max-m_min)/(num_freq-1)
      allocate(m_array(num_freq),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating m_array in pt_setup subroutine')
      allocate(local_pm(ndim_pm,num_freq,counts(my_node_id)),stat=ierr)
      if (ierr/=0 ) call io_error('Error in allocating pm in pt_setup subroutine')
      if (on_root) then
         allocate(global_pm(ndim_pm,num_freq,tems_num), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_pm in pt_setup')
      else
         allocate(global_pm(1,1,1),stat=ierr)
         if (ierr/=0) call io_error('Error in allocating global_pm in pt_setup')
      endif

      do ifreq=1,num_freq
         m_array(ifreq) = m_min+real(ifreq-1,dp)*d_m
      enddo
!!      allocate(seed_pt(num_nodes),stat=ierr)
!!      if (ierr/=0) call io_error('Error in allocating seed in pt_setup')
!!
!!      call random_seed()
!      if(.not.lseed) then
!        open(89,FILE='/dev/urandom',ACCESS='stream',FORM='UNFORMATTED')
!        read(89) seed(my_node_id+1)
!        close(89)
!      endif
!      call mtprng_init(seed(my_node_id+1),state)

      ! Initiale configuration

      if (index(initial_sconfig,'ferro')>0 ) then
         local_spin_matrix(1,:,:)=zero
         local_spin_matrix(2,:,:)=zero
         do loop_counts=1,counts(my_node_id)
            local_spin_matrix(3,:,loop_counts)=1.0_dp*mag_moments_matrix(:)
         enddo
      elseif (index(initial_sconfig,'rand')>0 ) then
         do loop=1,num_total_atoms

            r1=mtprng_rand_real2(state)
!            call random_number(r1)

            r2=mtprng_rand_real2(state)
!            call random_number(r2)
            local_spin_matrix(1,loop,:)=mag_moments_matrix(loop)*sin(r1*pi)*cos(r2*twopi)
            local_spin_matrix(2,loop,:)=mag_moments_matrix(loop)*sin(r1*pi)*sin(r2*twopi)
            local_spin_matrix(3,loop,:)=mag_moments_matrix(loop)*cos(r1*pi)
          enddo
      elseif (index(initial_sconfig,'file')>0 ) then

        if(on_root) then
 
          filename = trim(seedname)//'_sconfig.dat'
          !! Output is the normalized(spin_matrix)
          call mcarlo_read_sconfig(num_total_atoms,tems,global_spin_matrix,filename,tems_num)
          do loop=1,tems_num
             do loop_counts=1,3
                global_spin_matrix(loop_counts,:,loop) = global_spin_matrix(loop_counts,:,loop)*mag_moments_matrix(:)
             enddo
          enddo

          do loop_nodes=1,num_nodes-1
            do localidx = 1, counts(loop_nodes)
               globalidx = displs(loop_nodes)+ localidx
               local_spin_matrix(:,:,localidx) = global_spin_matrix(:,:,globalidx)
            enddo
            call comms_send(local_spin_matrix(1,1,1),3*num_total_atoms*counts(loop_nodes),loop_nodes)
          enddo

          do localidx = 1,counts(0)
             local_spin_matrix(:,:,localidx) = global_spin_matrix(:,:,localidx)
          enddo

        endif ! on_root

        if(.not. on_root) then
          call comms_recv(local_spin_matrix(1,1,1),3*num_total_atoms*counts(my_node_id),root_id)
        end if

      end if
 
      if (on_root) then 
         do loop=1,tems_num
            map(loop) = loop
         enddo
      endif

      return
 
    end subroutine pt_setup

    !============================================!
    subroutine pt_mcarlo()
      !============================================!
 
      implicit none
 
      integer(8)                 :: num_steps_tot,num_steps
      integer(8)                 :: num_sfsteps_tot(counts(my_node_id))
      integer                    :: swap_start
      integer                    :: pt_num_swap_warmup,pt_num_swap      
      integer                    :: loop_at,loop_mc,loop_msure
      integer                    :: loop,loop_sw,i
      integer                    :: localidx,globalidx,ifreq
      integer                    :: l,m_l,excc,nsp1,nsp2
      integer                    :: sconfig_warmup_unit
      logical                    :: sconfig_found 
      real(kind=dp)              :: local2_oc,beta
      real(kind=dp)              :: specific_heat,susceptibility
      real(kind=dp)              :: U_E,U_M,Delta_E,Delta_M,normalize
      real(kind=dp)              :: sus_orderparam,U_orderparam
      real(kind=dp)              :: sus_staggered,U_staggered
 
      ! Here I determine the number of swaps.
      !! First I allow to warm up the system as simple monte carlo.
      !! I use half of steps_warmup
      !! for the simple monte carlo and the half of it for the PT.
      steps_warmup=ceiling(real(steps_warmup,dp)/2.0_dp)
      pt_num_swap_warmup = ceiling(real(steps_warmup,dp)/real(pt_steps_swap,dp))
      if (on_root .and.(pt_steps_swap*pt_num_swap_warmup .ne. steps_warmup)) then
         write(stdout,'(1x,a)')     '***** The number of Warm-up steps is not diviable to pt_steps_swap : '
         write(stdout,'(1x,a,I9)')  '***** For consistancy the number of Warm-up steps increases to     : '&
                                    ,2*pt_num_swap_warmup*pt_steps_swap
      endif 
      steps_warmup = pt_num_swap_warmup*pt_steps_swap 
 
      pt_num_swap = ceiling(real(steps_mc,dp)/real(pt_steps_swap,dp))
      if (on_root .and. ((pt_num_swap*pt_steps_swap) .ne. steps_mc)) then 
         write(stdout,'(1x,a)')     '***** The number of Monte Carlo steps is not diviable to pt_steps_swap : '
         write(stdout,'(1x,a,I9)')  '***** For consistancy the number of Monte Carlo steps increases to     : '&
                                    ,pt_num_swap*pt_steps_swap 
      endif
      steps_mc = pt_num_swap*pt_steps_swap 

 
      local_rejected = 0; local_data = zero
      local_energy = zero; local_pm = zero
      if (lbinerror) then
         local_binning = zero; local_sum = zero
      endif
      if (lspincorr) then
         local_spincorr = zero; local_spincorr_abs = zero ; local_spincorr_tot_abs = zero
      endif

      if (on_root) call io_stopwatch('pt_main: warmup1',1)
      if (steps_warmup>0) then !
         do localidx = 1, counts(my_node_id)
            globalidx = displs(my_node_id)+ localidx
            beta = 1.0_dp/tems(globalidx)
            !! Initial energy 
            call energy_get(local_spin_matrix(:,:,localidx),local_energy(2,localidx))
            local_energy(1,localidx)=local_energy(2,localidx) 
            !! Thermalization
            num_steps_tot = 0
            do loop_mc=1, steps_warmup   ! monte carlo steps
               do loop_msure=1, steps_measure
                  do loop_at=1,num_total_atoms
                     call metropolis(loop_at,beta,local_spin_matrix(:,:,localidx),&
                                local_energy(2,localidx),local_rejected(localidx))
                     num_steps_tot = num_steps_tot+1
                  enddo ! atom
               enddo    ! measure
            enddo ! mc
         enddo ! localidx
       
         ! I calculate the final Energy at this step by direct method and save at local_energy(3,:)
         do localidx=1,counts(my_node_id)
            call energy_get(local_spin_matrix(:,:,localidx),local_energy(3,localidx))
         enddo
         ! gather energy & spin_matrix
         call comms_barrier       
         call comms_gatherv(local_energy(1,1),3*counts(my_node_id),global_energy(1,1), &
                            3*counts,3*displs) 
         call comms_gatherv(local_spin_matrix(1,1,1),3*num_total_atoms*counts(my_node_id)&
                            ,global_spin_matrix(1,1,1),3*num_total_atoms*counts,3*num_total_atoms*displs) 
         call comms_gatherv(local_rejected(1),counts(my_node_id),global_rejected(1), &
                 counts,displs)                                              
         if(on_root) then
           write(stdout,'(/1x,a)')'*------------------------------------ WARMUP1 --------------------------------------*'
           write(stdout,'(1x,a)') '|     Temp      Total Steps    Rejected Steps    Accepted Steps   Acc/Tot Ratio     |'
           write(stdout,'(1x,a)') '+-----------------------------------------------------------------------------------+'
       
           do loop=1,tems_num
              write(stdout,124) tems(loop),num_steps_tot,global_rejected(loop),&
                    num_steps_tot-global_rejected(loop),&
                    (num_steps_tot-global_rejected(loop))/real(num_steps_tot,dp),'  <-- WARMUP1'
           enddo
       
           write(stdout,'(/1x,a)')'*------------------------------------ WARMUP1 --------------------------------------*'
           write(stdout,'(1x,a)') '|   Temp             E_init             E_final(by adding)       E_final(direct)    |'
           write(stdout,'(1x,a)') '+-----------------------------------------------------------------------------------+'
           do loop=1,tems_num
              write(stdout,127) tems(loop),(global_energy(i,loop),i=1,3),'  <-- ENERGY1'
           enddo
       
           write(stdout,'(1x,a)')'*----------------------------------- END of WARMUP1 --------------------------------*'
           write(stdout,*) ' '

           inquire(file=trim(seedname)//'_sconfig_warmup.dat',exist=sconfig_found)
           if (sconfig_found) then
              sconfig_warmup_unit=io_file_unit()
              open(unit=sconfig_warmup_unit,file=trim(seedname)//'_sconfig_warmup.dat',status='old',position='append')
              close(sconfig_warmup_unit,status='delete')
           end if
           
           sconfig_warmup_unit=io_file_unit()
           open(sconfig_warmup_unit,file=trim(seedname)//'_sconfig_warmup.dat',status='unknown',form='formatted',action='write')
           call io_date(cdate,ctime)
           write(sconfig_warmup_unit,*) 'Created on '//cdate//' at '//ctime ! Date and time
           write(sconfig_warmup_unit,'(2I6)') num_total_atoms,tems_num
           
           do loop=1,tems_num                                         
              write(sconfig_warmup_unit,'(F12.6)') tems(loop)
              do l=1,num_total_atoms
                 write(sconfig_warmup_unit,'(3F12.6)') (global_spin_matrix(i,l,loop),i=1,3)
              enddo
           enddo
           close(sconfig_warmup_unit)
         endif
       
         call comms_barrier
         ! scatter spin_matrix,energy
!         call comms_scatterv(local_spin_matrix(1,1,1),3*num_total_atoms*counts(my_node_id),&
!                            global_spin_matrix(1,1,1),3*num_total_atoms*counts,3*num_total_atoms*displs) 
!       
!         call comms_scatterv(local_energy(1,1),3*counts(my_node_id),&
!                               global_energy(1,1),3*counts,3*displs) 
       
         if(on_root) call io_stopwatch('pt_main: warmup1',2)
       
         if(on_root) call io_stopwatch('pt_main: warmup2',1)
       
         ! After serial thermalization, I have thermalization with swap configuration
         if(on_root) then
           write(stdout,'(/1x,a)')'*---------------------------------- SWAP in WARMUP2 --------------------------------*'
           write(stdout,'(1x,a)') '|     Num_Steps                     Configuration Numbers                           |'
           write(stdout,'(1x,a)') '+-----------------------------------------------------------------------------------+'
           if (on_root) write(stdout,*) 'Initial configuration = ',map
         endif
       
         local_energy(1,:) = local_energy(2,:)
         local_rejected = 0; num_steps_tot = 0; num_swaps_tot = 0
         do loop_sw=1, pt_num_swap_warmup   ! swap configurations
            do localidx = 1, counts(my_node_id)
               globalidx = displs(my_node_id)+ localidx
               beta = 1.0_dp/tems(globalidx)
               num_steps_tot = 0
               do loop_mc=1, pt_steps_swap   ! monte carlo steps
                  num_steps = (loop_sw-1)*pt_steps_swap+loop_mc
                  if (on_root .and. ((mod(num_steps,pt_print_swap).eq. 0) .or.&
                      (num_steps.eq.1) .or. (num_steps.eq.steps_warmup))) write(stdout,*) num_steps,' | ',map
                  do loop_msure=1, steps_measure
                     do loop_at=1,num_total_atoms
                        call metropolis(loop_at,beta,local_spin_matrix(:,:,localidx),& 
                                       local_energy(2,localidx),local_rejected(localidx))
                        num_steps_tot = num_steps_tot+1
                     enddo ! atom
                  enddo    ! measure
               enddo !mc
            enddo ! localidx
            ! Gather energy & spin_matrix
            call comms_barrier       
            call comms_gatherv(local_energy(1,1),3*counts(my_node_id),global_energy(1,1), &
                               3*counts,3*displs) 
            call comms_gatherv(local_spin_matrix(1,1,1),3*num_total_atoms*counts(my_node_id)&
                               ,global_spin_matrix(1,1,1),3*num_total_atoms*counts,3*num_total_atoms*displs) 
            ! SWAP configurations
            swap_start = mod(loop_sw+1,2)+1
            if (on_root) call pt_swap(swap_start,global_spin_matrix,global_energy(2,:),map,num_swaps_tot)
       
            call comms_barrier
            ! Scatter spin_matrix,energy
            call comms_scatterv(local_energy(1,1),3*counts(my_node_id),&
                            global_energy(1,1),3*counts,3*displs) 
            call comms_scatterv(local_spin_matrix(1,1,1),3*num_total_atoms*counts(my_node_id),&
                            global_spin_matrix(1,1,1),3*num_total_atoms*counts,3*num_total_atoms*displs) 
         enddo  !! loop_sw
         if (on_root .and. (num_steps.eq.steps_warmup)) write(stdout,*) 'Final configuration = ',map
         if (on_root) write(stdout,'(1x,85a)') '*------------------------------------------------&
                                                &-----------------------------------*'
       
         ! I calculate the final Energy at this step by direct method and save at local_energy(3,:)
         do localidx=1,counts(my_node_id)
            call energy_get(local_spin_matrix(:,:,localidx),local_energy(3,localidx))
         enddo
         ! Gather energy & spin_matrix
         call comms_barrier       
         call comms_gatherv(local_energy(1,1),3*counts(my_node_id),global_energy(1,1), &
                            3*counts,3*displs) 
         call comms_gatherv(local_spin_matrix(1,1,1),3*num_total_atoms*counts(my_node_id)&
                            ,global_spin_matrix(1,1,1),3*num_total_atoms*counts,3*num_total_atoms*displs) 
         call comms_gatherv(local_rejected(1),counts(my_node_id),global_rejected(1), &
                 counts,displs)                                              
         if (on_root) then
           write(stdout,'(/1x,85a)')'*----------------------------------- WARMUP2 ---------------------------------------*'
           write(stdout,'(1x,85a)') '|  Temp_i     Temp_f         Tot-Swaps       Rej-Swaps      Acc-Swaps  Acc/Tot Ratio|'
           write(stdout,'(1x,85a)') '+-----------------------------------------------------------------------------------+'
           do loop=1,tems_num
              write(stdout,128) tems(loop),tems(map(loop)),num_swaps_tot(1,loop)&
                                ,num_swaps_tot(1,loop)-num_swaps_tot(2,loop),num_swaps_tot(2,loop),&
                                 num_swaps_tot(2,loop)/real(num_swaps_tot(1,loop),dp)
           enddo
       
           write(stdout,'(/1x,a)')'*------------------------------------ WARMUP2 --------------------------------------*'
           write(stdout,'(1x,a)') '|     Temp      Total Steps    Rejected Steps    Accepted Steps   Acc/Tot Ratio     |'
           write(stdout,'(1x,a)') '+-----------------------------------------------------------------------------------+'
       
           do loop=1,tems_num
              write(stdout,124) tems(loop),pt_num_swap_warmup*num_steps_tot,global_rejected(loop),&
                 pt_num_swap_warmup*num_steps_tot-global_rejected(loop),&
                 (pt_num_swap_warmup*num_steps_tot-global_rejected(loop))/real(pt_num_swap_warmup*num_steps_tot,dp),'  <-- WARMUP2'
           enddo
       
           write(stdout,'(/1x,a)')'*------------------------------------ WARMUP2 --------------------------------------*'
           write(stdout,'(1x,a)') '|   Temp             E_init             E_final(by adding)       E_final(direct)    |'
           write(stdout,'(1x,a)') '+-----------------------------------------------------------------------------------+'
           do loop=1,tems_num
              write(stdout,127) tems(loop),(global_energy(i,loop),i=1,3),'  <-- ENERGY2'
           enddo
       
           write(stdout,'(1x,a)')'*---------------------------------- END of WARMUP2 ---------------------------------*'
         endif
       
!         call comms_scatterv(local_spin_matrix(1,1,1),3*num_total_atoms*counts(my_node_id),&
!                            global_spin_matrix(1,1,1),3*num_total_atoms*counts,3*num_total_atoms*displs) 
!    
!         call comms_scatterv(local_energy(1,1),3*counts(my_node_id),&
!                               global_energy(1,1),3*counts,3*displs) 
    
         if (on_root) call io_stopwatch('pt_main: warmup2',2)
      else
         if (on_root) then
            write(stdout,'(/1x,85a)') '+------------------------------------ WARMUP --------------------------------------+'
            write(stdout,'(1x,85a)') '|                                                                                   |'                                                             
            write(stdout,'(1x,85a)') '+-----------------------------------------------------------------------------------+'
            write(stdout,'(1x,85a)') ' No Warm-up step is done! '                                                              
            write(stdout,'(1x,a)')   '*--------------------------------- END of WARMUP -----------------------------------*'
         endif
         ! I calculate the Initial  Energy if no warm-up is done!
         do localidx=1,counts(my_node_id)
            call energy_get(local_spin_matrix(:,:,localidx),local_energy(2,localidx))
         enddo
      endif
   
      if (on_root) call io_stopwatch('pt_main: main',1)
   
      !! Monte Carlo
      if (on_root) then
        write(stdout,'(/1x,85a)')'*------------------------------ SWAP in Monte-Carlo --------------------------------*'
        write(stdout,'(1x,85a)') '|     Num_Steps                     Configuration Numbers                           |'
        write(stdout,'(1x,85a)') '+-----------------------------------------------------------------------------------+'
        if (on_root) write(stdout,*) 'Initial configuration = ',map
      endif
   
      local_energy(1,:) = local_energy(2,:)
      local_rejected = 0
      num_steps_tot = 0; num_swaps_tot = 0; num_sfsteps_tot = 0

      !! I need excc only when energy_write is true 
      !excc = (energy_num_print/pt_steps_swap)*pt_steps_swap
 
      do loop_sw=1, pt_num_swap   ! swap configurations
         do localidx = 1, counts(my_node_id)
            globalidx = displs(my_node_id)+localidx
            beta = 1.0_dp/tems(globalidx)
            num_steps_tot = 0
            do loop_mc=1,pt_steps_swap   ! monte carlo steps
               num_steps = (loop_sw-1)*pt_steps_swap+loop_mc
               if (on_root .and. ((mod(num_steps,pt_print_swap).eq. 0) .or.&
                   (num_steps.eq.1) .or. (num_steps.eq. steps_mc))) write(stdout,*) num_steps,' | ',map
               do loop_msure=1,steps_measure
                  do loop_at=1,num_total_atoms
                     call metropolis(loop_at,beta,local_spin_matrix(:,:,localidx),& 
                                    local_energy(2,localidx),local_rejected(localidx))
                     num_steps_tot = num_steps_tot+1
                  enddo ! atom
               enddo    ! measure
               local_data(E,localidx)  = local_data(E,localidx)+local_energy(2,localidx)
               local_data(E2,localidx) = local_data(E2,localidx)+local_energy(2,localidx)**2
               local_data(E4,localidx) = local_data(E4,localidx)+local_energy(2,localidx)**4
      
               call magnetization_get(local_spin_matrix(:,:,localidx),local2_oc)
               local_data(M,localidx)  = local_data(M,localidx)+sqrt(local2_oc)
               local_data(M2,localidx) = local_data(M2,localidx)+local2_oc
               local_data(M4,localidx) = local_data(M4,localidx)+local2_oc**2

               if (lspincorr) call spincorr_get(local_spin_matrix(:,:,localidx),&
                                 local_spincorr(:,:,:,localidx),local_spincorr_abs(:,:,:,localidx),&
                                 local_spincorr_tot_abs(:,:,:,localidx))

               if (lbinerror) then
                  local_sum(1,:,localidx) = local_sum(1,:,localidx)+&
                                            (local_energy(2,localidx)/real(num_total_atoms,dp))
                  local_sum(2,:,localidx) = local_sum(2,:,localidx)+sqrt(local2_oc)         
               endif

               if (lorderparam) then
                  call orderparam_get(local_spin_matrix(:,:,localidx),local2_oc)
                  local_data(OP,localidx)  = local_data(OP,localidx)+sqrt(local2_oc)
                  local_data(OP2,localidx) = local_data(OP2,localidx)+local2_oc
                  local_data(OP4,localidx) = local_data(OP4,localidx)+local2_oc**2
                  if (lbinerror) &
                     local_sum(3,:,localidx) = local_sum(3,:,localidx)+sqrt(local2_oc)         
               endif

               if (lstaggered) then
                  call staggered_get(local_spin_matrix(:,:,localidx),local2_oc)
                  local_data(SM,localidx)  = local_data(SM,localidx)+sqrt(local2_oc)
                  local_data(SM2,localidx) = local_data(SM2,localidx)+local2_oc
                  local_data(SM4,localidx) = local_data(SM4,localidx)+local2_oc**2
                  ! ndim_sum = 4  when lorderparam = True
                  ! ndim_sum = 3  when lorderparam = False 
                  if (lbinerror) &
                     local_sum(ndim_sum,:,localidx) = local_sum(ndim_sum,:,localidx)+sqrt(local2_oc)         
               endif

               if (lorderparam.and.lstaggered) then
                  call Pm_get(local_spin_matrix(:,:,localidx),num_freq,local_pm(1,:,localidx),&
                               pdf_order=local_pm(2,:,localidx),pdf_staggered=local_pm(3,:,localidx))
               elseif (lorderparam) then
                  call Pm_get(local_spin_matrix(:,:,localidx),num_freq,local_pm(1,:,localidx),&
                               pdf_order=local_pm(2,:,localidx))
               elseif (lstaggered) then
                  call Pm_get(local_spin_matrix(:,:,localidx),num_freq,local_pm(1,:,localidx),&
                               pdf_staggered=local_pm(2,:,localidx))
               else
                  call Pm_get(local_spin_matrix(:,:,localidx),num_freq,local_pm(1,:,localidx))
               endif

               if (lbinerror) then
                  do l=1,num_binning_level
                     if (mod(num_steps,(2**(binning_level(l)))).eq. 0) then
                        local_sum(:,l,localidx)     = local_sum(:,l,localidx)/real(2**l,dp)
                        local_binning(:,l,localidx) = local_binning(:,l,localidx)+local_sum(:,l,localidx)**2
                        local_sum(:,l,localidx)     = zero
                     endif
                  enddo
               endif

               if (lsfactor .and. (mod(num_steps,sfactor_steps_measure) .eq. 0)) then
                  call mcarlo_structure_factor(local_spin_matrix(:,:,localidx),&
                             atoms_supercell_pos_cart,local_sfactor_matrix(:,:,localidx))
                  num_sfsteps_tot(localidx) = num_sfsteps_tot(localidx)+1
               endif

               if (energy_write) then
                  !local_E(mod(num_steps,excc),localidx-1) = local_energy(2,localidx)
                  local_E(mod(num_steps,energy_num_print),localidx-1) = local_energy(2,localidx)

               endif

            enddo !mc
         enddo ! localidx

         if (energy_write) then
            !if (mod(pt_steps_swap*loop_sw,excc).eq. 0) then
            if (mod(pt_steps_swap*loop_sw,energy_num_print).eq. 0) then
               call comms_gatherv(local_E(0,0),energy_num_print*counts(my_node_id),global_E(0,0), &
                       energy_num_print*counts,energy_num_print*displs)
               if (on_root) &
               !call write_E(pt_steps_swap*loop_sw,tems_num,excc,global_E(0:excc-1,0:tems_num-1),energy_unit)
               call write_E(pt_steps_swap*loop_sw,tems_num,energy_num_print,&
                            global_E(0:energy_num_print-1,0:tems_num-1),energy_unit)
            endif
            if ((loop_sw .eq. pt_num_swap) .and. & 
                !( steps_mc .gt. ((steps_mc/excc)*excc)))then
                ( steps_mc .gt. ((steps_mc/energy_num_print)*energy_num_print)))then
               !excc =  steps_mc-(steps_mc/excc)*excc  
               excc =  steps_mc-(steps_mc/energy_num_print)*energy_num_print
               local_E(0,:) = local_E(excc,:)
               call comms_gatherv(local_E(0,0),energy_num_print*counts(my_node_id),global_E(0,0), &
                       energy_num_print*counts,energy_num_print*displs)
               if (on_root) &
               call write_E(steps_mc,tems_num,excc,global_E(0:excc-1,0:tems_num-1),energy_unit)
            endif
         endif

         ! Gather energy & spin_matrix
         call comms_barrier       
         call comms_gatherv(local_energy(1,1),3*counts(my_node_id),global_energy(1,1), &
                            3*counts,3*displs) 
         call comms_gatherv(local_spin_matrix(1,1,1),3*num_total_atoms*counts(my_node_id)&
                            ,global_spin_matrix(1,1,1),3*num_total_atoms*counts,3*num_total_atoms*displs) 
         call comms_gatherv(local_rejected(1),counts(my_node_id),global_rejected(1), &
                 counts,displs)                                              
         !! SWAP Configurations
         swap_start = mod(loop_sw+1,2)+1
         if (on_root) call pt_swap(swap_start,global_spin_matrix,global_energy(2,:),map,num_swaps_tot)
 
         call comms_barrier
         ! Scatter spin_matrix,energy
         call comms_scatterv(local_spin_matrix(1,1,1),3*num_total_atoms*counts(my_node_id),&
                         global_spin_matrix(1,1,1),3*num_total_atoms*counts,3*num_total_atoms*displs) 
  
         call comms_scatterv(local_energy(1,1),3*counts(my_node_id),&
                         global_energy(1,1),3*counts,3*displs) 
      enddo ! loop_sw
      if (on_root .and. (num_steps.eq. steps_mc)) write(stdout,*) 'Final configuration = ',map
      if (on_root) write(stdout,'(1x,85a)') '*-----------------------------------'&
                                          //'------------------------------------------------*'


      if (lbinerror) then
         !! if (mod(steps_mc,2^l) .ne. 0) i.e. the last binning box is incomplete
         !! Therfore, the local_sum is non zero and I should remove them from the average <E> or <M> ,...
         do l=1,num_binning_level
            local_sum(1,l,:) = (local_data(E,:)/real(num_total_atoms,dp))-local_sum(1,l,:)
            local_sum(2,l,:) = local_data(M,:)-local_sum(2,l,:)
            if (lorderparam.and.lstaggered) then
               local_sum(3,l,:) = local_data(OP,:)-local_sum(3,l,:)
               local_sum(4,l,:) = local_data(SM,:)-local_sum(4,l,:)
            elseif (lorderparam.or.lstaggered) then
               ! if lorderparam= T and lstaggered= F, OP=7
               ! if lorderparam= F and lstaggered= T, SM=OP=7
               local_sum(3,l,:) = local_data(OP,:)-local_sum(3,l,:)
            endif
         
            m_l = (steps_mc)/(2**(binning_level(l)))
            local_sum(:,l,:)     = local_sum(:,l,:)/real((2**(binning_level(l))),dp)
            local_sum(:,l,:)     = local_sum(:,l,:)/real(m_l,dp)
            local_binning(:,l,:) = local_binning(:,l,:)/real(m_l,dp)
            local_binning(:,l,:) = (local_binning(:,l,:)-local_sum(:,l,:)**2)/real((m_l-1),dp)
            local_binning(:,l,:) = local_binning(:,l,:)/real(num_total_atoms,dp)
         enddo
      endif !autocorr


      ! I calculate the final Energy by direct method and save at local_energy(3,:)
      do localidx=1,counts(my_node_id)
         call energy_get(local_spin_matrix(:,:,localidx),local_energy(3,localidx))
      enddo
   
      ! Now I gather all data in monte carlo step on 
      ! root node(local_data,local_energy,local_rejected,local_sfactor_matrix)
      !NOTE****** I changed the sfactor_matrix dimension **************
      !*******from (0:product()-1,:) to (product(),:)  ************
      call comms_gatherv(local_rejected(1),counts(my_node_id),global_rejected(1), &
                         counts,displs)                                              
      call comms_gatherv(local_energy(1,1),3*counts(my_node_id),global_energy(1,1), &
                         3*counts,3*displs)                                            
      call comms_gatherv(local_spin_matrix(1,1,1),3*num_total_atoms*counts(my_node_id)&
                         ,global_spin_matrix(1,1,1),3*num_total_atoms*counts,3*num_total_atoms*displs) 
      call comms_gatherv(local_data(1,1),ndim_data*counts(my_node_id)&
                         ,global_data(1,1),ndim_data*counts,ndim_data*displs)
      call comms_gatherv(local_pm(1,1,1),ndim_pm*num_freq*counts(my_node_id)&
                         ,global_pm(1,1,1),ndim_pm*num_freq*counts,ndim_pm*num_freq*displs)
       if (lbinerror) &
          call comms_gatherv(local_binning(1,1,1),ndim_sum*num_binning_level*counts(my_node_id),global_binning(1,1,1), &
                             ndim_sum*num_binning_level*counts,ndim_sum*num_binning_level*displs)                                            
       if (lsfactor) then
          call comms_gatherv(local_sfactor_matrix(1,1,1),3*sfactor_nqpts*counts(my_node_id),& 
                             global_sfactor_matrix(1,1,1),3*sfactor_nqpts*counts,&
                             3*sfactor_nqpts*displs) 
       endif                                 
       if (lspincorr) then
          call comms_gatherv(local_spincorr(1,1,1,1),jij_shell_max*num_species*num_species*counts(my_node_id)&
                             ,global_spincorr(1,1,1,1),jij_shell_max*num_species*num_species*counts,&
                             jij_shell_max*num_species*num_species*displs) 
          call comms_gatherv(local_spincorr_abs(1,1,1,1),jij_shell_max*num_species*num_species*counts(my_node_id)&
                             ,global_spincorr_abs(1,1,1,1),jij_shell_max*num_species*num_species*counts,&
                             jij_shell_max*num_species*num_species*displs) 
          call comms_gatherv(local_spincorr_tot_abs(1,1,1,1),jij_shell_max*num_species*num_species*counts(my_node_id)&
                             ,global_spincorr_tot_abs(1,1,1,1),jij_shell_max*num_species*num_species*counts,&
                             jij_shell_max*num_species*num_species*displs) 
       endif                                 

       call comms_barrier

       if (on_root) then                                                    

          write(stdout,*)'' 
          write(stdout,'(/1x,a)')'*----------------------------------- MONTECARLO ------------------------------------*'
          write(stdout,'(1x,a)') '|  Temp_i     Temp_f         Tot-Swaps       Rej-Swaps      Acc-Swaps  Acc/Tot Ratio|'
          write(stdout,'(1x,a)') '+-----------------------------------------------------------------------------------+'
          do loop=1,tems_num
             write(stdout,128) tems(loop),tems(map(loop)),num_swaps_tot(1,loop)&
                               ,num_swaps_tot(1,loop)-num_swaps_tot(2,loop),num_swaps_tot(2,loop),&
                                num_swaps_tot(2,loop)/real(num_swaps_tot(1,loop),dp)
          enddo
          
          write(stdout,'(/1x,a)')'*----------------------------------- MONTECARLO ------------------------------------*'
          write(stdout,'(1x,a)') '|     Temp      Total Steps    Rejected Steps    Accepted Steps   Acc/Tot Ratio     |'
          write(stdout,'(1x,a)') '+-----------------------------------------------------------------------------------+'
          normalize   = 1.0_dp/real(steps_mc,dp)
          global_data =  global_data*normalize
          global_pm   =  global_pm*normalize/real(num_supercell,dp)
          if (lspincorr) then
            global_spincorr         = global_spincorr*normalize/real(num_total_atoms,dp)
            global_spincorr_abs     = global_spincorr_abs*normalize/real(num_total_atoms,dp)
            global_spincorr_tot_abs = global_spincorr_tot_abs*normalize/real(num_total_atoms,dp)
          endif
          
          do loop=1,tems_num                                         
             beta = 1.0_dp/tems(loop)
             specific_heat = beta*beta*(global_data(E2,loop)-global_data(E,loop)**2)/real(num_total_atoms,dp)
             susceptibility = beta*(global_data(M2,loop)-global_data(M,loop)**2)*num_total_atoms
             U_E = 1.0_dp-(1.0_dp/3.0_dp)*(global_data(E4,loop)/global_data(E2,loop)**2)
             U_M = 1.0_dp-(1.0_dp/3.0_dp)*(global_data(M4,loop)/global_data(M2,loop)**2)
             Delta_E = sqrt((global_data(E2,loop)-global_data(E,loop)**2)*normalize)/real(num_total_atoms,dp)
             Delta_M = sqrt((global_data(M2,loop)-global_data(M,loop)**2)*normalize)
          
             write(stdout,124) tems(loop),pt_num_swap*num_steps_tot,global_rejected(loop),&
                   pt_num_swap*num_steps_tot-global_rejected(loop),&                   
                   (pt_num_swap*num_steps_tot-global_rejected(loop))/real(pt_num_swap*num_steps_tot,dp),'  <-- MCARLO'
             write(stdout,'(4x,a,E15.7,a)')  ' Initial Energy                           = ',global_energy(1,loop),'  <-- ENERGY3'
             write(stdout,'(4x,a,E15.7,a)')  ' Final Energy by adding differences       = ',global_energy(2,loop),'  <-- ENERGY3'
             write(stdout,'(4x,a,E15.7,a)')  ' Final Energy by direct calculations      = ',global_energy(3,loop),'  <-- ENERGY3'
             write(stdout,'(4x,a,f15.9)')    '                             Delta E      = +/- ',Delta_E
             write(stdout,'(4x,a,f15.9)')    '                             Delta M      = +/- ',Delta_M
             write(stdout,'(1x,a85)') repeat('-',85)                           

             write(mc_unit,125) tems(loop),global_data(M,loop),&
                   global_data(E,loop)/real(num_total_atoms,dp),specific_heat,susceptibility,U_E,U_M
          
             if(lspincorr)then
               do l=1,jij_shell_max
                  do nsp1=1,num_species
                     do nsp2=1,num_species
                        write(spincorr_unit,'(f9.4,1x,I5,6x,a3,4x,i3,6x,a3,4x,i3,8x,f12.9,10x,f12.9,12x,f12.9)') &
                              tems(loop),l,atoms_symbol(nsp1),nsp1,atoms_symbol(nsp2),nsp2,&
                              global_spincorr(l,nsp1,nsp2,loop),global_spincorr_abs(l,nsp1,nsp2,loop),&
                              global_spincorr_tot_abs(l,nsp1,nsp2,loop)
                     enddo
                  enddo
               enddo               
             endif
          
             if (lorderparam) then
                sus_orderparam = beta*(global_data(OP2,loop)-global_data(OP,loop)**2)*num_total_atoms
                U_orderparam   = 1.0_dp-(1.0_dp/3.0_dp)*(global_data(OP4,loop)/global_data(OP2,loop)**2)
                write(orderparam_unit,126) tems(loop),global_data(OP,loop),sus_orderparam,U_orderparam
             endif

             if (lstaggered) then
                sus_staggered = beta*(global_data(SM2,loop)-global_data(SM,loop)**2)*num_total_atoms
                U_staggered   = 1.0_dp-(1.0_dp/3.0_dp)*(global_data(SM4,loop)/global_data(SM2,loop)**2)
                write(staggered_unit,126) tems(loop),global_data(SM,loop),sus_staggered,U_staggered
             endif
         
             if(lbinerror)then
               write(binerror_unit,'(1x,a,2x,f6.2)') '#T=',tems(loop)
               do l=1,num_binning_level 
                    if (lorderparam.and.lstaggered) then
                       write(binerror_unit,'(2x,I8,8x,E15.8,3x,E15.8,3x,E15.8,3x,E15.8)') binning_level(l),&
                       !     (sqrt(global_binning(i,l,loop)/real(num_total_atoms,dp)),i=1,4)
                            (sqrt(global_binning(i,l,loop)),i=1,4)
                    elseif (lorderparam.or.lstaggered) then
                       write(binerror_unit,'(2x,I8,8x,E15.8,3x,E15.8,3x,E15.8)') binning_level(l),&
                       !     (sqrt(global_binning(i,l,loop)/real(num_total_atoms,dp)),i=1,3)
                            (sqrt(global_binning(i,l,loop)),i=1,3)
                    else
                       write(binerror_unit,'(2x,I8,8x,E15.8,3x,E15.8)') binning_level(l), &
                       !     (sqrt(global_binning(i,l,loop)/real(num_total_atoms,dp)),i=1,2)
                            (sqrt(global_binning(i,l,loop)),i=1,2)
                    endif
               enddo
             endif
          
             write(pm_unit,'(1x,a,2x,f9.4)') '#T=',tems(loop)
          
             do ifreq=1,num_freq
                if (lorderparam.and.lstaggered) then
                   write(pm_unit,'(2x,E15.8,10x,E15.8,10x,E15.8,10x,E15.8)') m_array(ifreq),(global_pm(i,ifreq,loop),i=1,3)
                elseif (lorderparam.or.lstaggered) then
                   write(pm_unit,'(2x,E15.8,10x,E15.8,10x,E15.8)') m_array(ifreq),(global_pm(i,ifreq,loop),i=1,2)
                else
                   write(pm_unit,'(2x,E15.8,10x,E15.8)') m_array(ifreq),global_pm(1,ifreq,loop)
                endif
             enddo
          
             write(sconfig_unit,'(F12.6)') tems(loop)
             do l=1,num_total_atoms
                write(sconfig_unit,'(3F12.6)') (global_spin_matrix(i,l,loop),i=1,3)
             enddo
          
          enddo ! loop                                                            
         
          if (lorderparam.and.lstaggered) then
             write(stdout,'(4x,a)')  ' Order Parameters and Staggered Magnetization caculating ..done   '
             write(stdout,'(1x,a85)') repeat('-',85)
          elseif (lorderparam) then
             write(stdout,'(4x,a)')  ' Order Parameters caculating ..done   '
             write(stdout,'(1x,a85)') repeat('-',85)
          elseif (lstaggered) then
             write(stdout,'(4x,a)')  ' Staggered Magnetization caculating ..done   '
             write(stdout,'(1x,a85)') repeat('-',85)
          endif
 
          if (lsfactor) then
             global_sfactor_matrix = global_sfactor_matrix/real(num_total_atoms,dp)
             global_sfactor_matrix = global_sfactor_matrix/real(num_sfsteps_tot(1),dp)
             call mcarlo_structure_factor_write_pt(global_sfactor_matrix)
             write(stdout,'(4x,a,I14,17x,a)')  ' Total used steps at S Factor calculation = ',num_sfsteps_tot(1),'   <-- SFAC'
             write(stdout,'(1x,a85)') repeat('-',85)
          endif
 
       endif  ! on_root        
 
      call comms_barrier       
      if (on_root) call io_stopwatch('pt_main: main',2)

      if (on_root) then
         write(stdout,'(3X,A)') "Transport properties calculated."
         if(energy_write) call close_E_unit(energy_unit)
         write(stdout,*)
         write(stdout,'(1x,a)') '*-----------------------------------------------------------------------------------*'
         write(stdout,'(1x,a)') '|                             End of Parallel Tempering                             |'
         write(stdout,'(1x,a)') '*-----------------------------------------------------------------------------------*'
      end if
 
      if (on_root) call io_stopwatch('pt_main',2)


!124   format(2x,f9.4,1x,i14,1x,i14,5x,i14,7x,f14.8,1x,a)
124   format(3x,f9.4,1x,i14,1x,i14,5x,i14,7x,f10.8,4x,a)
125   format(f9.4,2x,E15.8,2x,E15.8,2x,E15.8,2x,E15.8,2x,E15.8,2x,E15.8)
126   format(3x,f9.4,5x,E15.8,6x,E15.8,8x,E15.8,5x)
127   format(2x,f9.4,8x,E15.8,8x,E15.8,8x,E15.8,2x,a)
128   format(2x,f9.4,2x,f9.4,1x,i14,1x,i14,1x,i14,3x,f14.8,2x)

      return

    end subroutine pt_mcarlo

    !==================================================================!
    subroutine pt_dealloc()
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  

      implicit none

      ! Deallocate real arrays that are private

      deallocate(global_spin_matrix, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating global_spin_matrix in pt_dealloc')
      deallocate(local_spin_matrix, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating spin_matrix in pt_dealloc')
      deallocate(tmp_spin_matrix, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating tmp_spin_matrix in pt_dealloc')

      deallocate(global_energy, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating global_energy in pt_dealloc')
      deallocate(local_energy, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating local_energy in pt_dealloc')

      deallocate(global_rejected, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating global_rejected in pt_dealloc')
      deallocate(local_rejected, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating local_rejected in pt_dealloc')

      deallocate(global_data, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating global_data in pt_dealloc')
      deallocate(local_data, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating local_data in pt_dealloc')

      deallocate(global_pm, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating global_pm in pt_dealloc')
      deallocate(local_pm, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating local_pm in pt_dealloc')
      deallocate(m_array, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating m_array in pt_dealloc')

      if (allocated( map ) ) then
         deallocate( map, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating map in pt_dealloc')
      end if
      if (allocated( num_swaps_tot )) then
         deallocate( num_swaps_tot, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating num_swaps_tot in pt_dealloc')
      end if

      if (lspincorr) then
        deallocate(global_spincorr, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating global_spincorr in pt_dealloc')
        deallocate(local_spincorr, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating local_spincorr in pt_dealloc')
        deallocate(global_spincorr_abs, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating global_spincorr_abs in pt_dealloc')
        deallocate(local_spincorr_abs, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating local_spincorr_abs in pt_dealloc')
        deallocate(global_spincorr_tot_abs, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating global_spincorr_tot_abs in pt_dealloc')
        deallocate(local_spincorr_tot_abs, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating local_spincorr_tot_abs in pt_dealloc')
      endif

      if (lbinerror) then
        deallocate(global_binning, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating global_binning in pt_dealloc')
        deallocate(local_binning, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating local_binning in pt_dealloc')
        deallocate(local_sum, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating local_sum in pt_dealloc')
      endif

      if (lsfactor) then
         deallocate( global_sfactor_matrix, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating global_sfactor_matrix in pt_dealloc')
         deallocate( local_sfactor_matrix, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating local_sfactor_matrix in pt_dealloc')
         deallocate( atoms_supercell_pos_cart, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating atoms_supercell_pos_cart in pt_dealloc')
      end if

      if (energy_write) then
         deallocate( energy_unit, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating energy_unit in pt_dealloc')
         deallocate( global_E, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating global_E in pt_dealloc')
         deallocate( local_E, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating local_E in pt_dealloc')
      end if

      return

    end subroutine pt_dealloc

  end subroutine pt_main

    !==================================================================!
    subroutine pt_swap(start,spin_matrix,energy,map1,num_swap)    
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  
   
      implicit none

      integer                      :: loop_swap,i
      logical                      :: swap
      real(kind=dp),intent(inout)  :: spin_matrix(3,num_total_atoms,tems_num)     
      real(kind=dp),intent(inout)  :: energy(tems_num)     
      integer(8),intent(inout)     :: num_swap(2,tems_num)     
      integer,intent(inout)        :: map1(tems_num)     
      integer,intent(in)           :: start     
      real(kind=dp)                :: pt_random,pt_fac     
      real(kind=dp)                :: tmp_energy     
      real(kind=dp)                :: spin_tmp(3,num_total_atoms)

      if(on_root) then
        do loop_swap=start,tems_num-1,2
           swap=.False.
           pt_fac=(1.0_dp/tems(loop_swap)-1.0_dp/tems(loop_swap+1))*&
                   (energy(loop_swap)-energy(loop_swap+1))
           num_swap(1,loop_swap)=num_swap(1,loop_swap)+1
           num_swap(1,loop_swap+1)=num_swap(1,loop_swap+1)+1
           if(pt_fac .gt. zero) then
             swap=.True.
           else
             pt_random=mtprng_rand_real2(state)
!            call random_number(pt_random)
             if(pt_random  .lt. exp(pt_fac)) swap=.True.
           endif      
 
           if(swap) then
             spin_tmp(:,:)=spin_matrix(:,:,loop_swap)
             spin_matrix(:,:,loop_swap)=spin_matrix(:,:,loop_swap+1)
             spin_matrix(:,:,loop_swap+1)=spin_tmp(:,:)
             tmp_energy                   = energy(loop_swap)
             energy(loop_swap)   = energy(loop_swap+1)
             energy(loop_swap+1) = tmp_energy
             i=map1(loop_swap)
             map1(loop_swap)=map1(loop_swap+1)
             map1(loop_swap+1)=i
             num_swap(2,loop_swap)=num_swap(2,loop_swap)+1
             num_swap(2,loop_swap+1)=num_swap(2,loop_swap+1)+1
           endif
        enddo ! loop_swap
      endif
 
      return

    end subroutine pt_swap   

 

end module mc_pt
