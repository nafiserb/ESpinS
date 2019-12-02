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

module mc_mcarlo 

  use mc_constants,  only : dp,zero,pi,twopi
  use mc_parameters, only : seeds,&
       steps_warmup,steps_mc,steps_measure,&
       num_supercell,num_atoms,num_total_atoms,num_species,&
       atoms_symbol,initial_sconfig,tems,tems_num,state,&
       energy_write,energy_num_print,lorderparam,lstaggered,&
       lsfactor,sfactor_steps_measure,sfactor_nqpts,supercell_size,&
       lbinerror,num_binning_level,binning_level,lspincorr,jij_shell_max
  use mc_io,         only : io_error,stdout,io_stopwatch,&
                            io_file_unit,seedname  
  use mc_common
  use mc_comms
  use mc_get_quant
!  use stdtypes
  use mtprng,        only : mtprng_rand_real2


  implicit none
  private 
  public :: mcarlo_main

  ! Constants to identify the nine components of a tensor when it is stored in packed form
  integer, parameter :: E   = 1
  integer, parameter :: E2  = 2
  integer, parameter :: E4  = 3
  integer, parameter :: M   = 4
  integer, parameter :: M2  = 5
  integer, parameter :: M4  = 6
  integer, parameter :: OP  = 7
  integer, parameter :: OP2 = 8 
  integer, parameter :: OP4 = 9 
  integer            :: SM  
  integer            :: SM2  
  integer            :: SM4  


contains 

  !============================================!
  subroutine mcarlo_main()
    !============================================!

    use mc_io,               only : io_date
    implicit none

    character(len=9)           :: cdate, ctime
    integer                    :: ierr,mc_unit,orderparam_unit
    integer                    :: binerror_unit,pm_unit,sconfig_unit
    integer                    :: spincorr_unit,staggered_unit
    integer                    :: num_freq,ndim_data,ndim_sum,ndim_pm
    integer,       allocatable :: energy_unit(:)
    integer(8),    allocatable :: global_rejected(:,:)
    real(kind=dp)              :: d_m                   
    real(kind=dp), allocatable :: global_tem(:)
    real(kind=dp), allocatable :: local_data(:),global_data(:,:)
    real(kind=dp), allocatable :: local_spin_matrix(:,:),global_spin_matrix(:,:,:)
    real(kind=dp), allocatable :: local_spin_matrix_file(:,:,:)
    real(kind=dp), allocatable :: global_energy(:,:)
    real(kind=dp), allocatable :: atoms_supercell_pos_cart(:,:)
    real(kind=dp), allocatable :: local_sfactor_matrix(:,:),global_sfactor_matrix(:,:,:)
    real(kind=dp), allocatable :: local_spincorr(:,:,:),global_spincorr(:,:,:,:)
    real(kind=dp), allocatable :: local_spincorr_abs(:,:,:),global_spincorr_abs(:,:,:,:)
    real(kind=dp), allocatable :: local_spincorr_tot_abs(:,:,:),global_spincorr_tot_abs(:,:,:,:)
    real(kind=dp), allocatable :: local_sum(:,:),local_binning(:,:),global_binning(:,:,:)
    real(kind=dp), allocatable :: local_E(:),global_E(:,:)
    real(kind=dp), allocatable :: m_array(:),local_pm(:,:),global_pm(:,:,:)

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

    if (on_root .and. mod(tems_num,num_nodes) .ne. 0) & 
        call io_error('number of temperature should be divisable and larger than number of nodes ')

    ! allocating the variables
    call mcarlo_setup()

    if (on_root) then
       sconfig_unit = io_file_unit()
       open(sconfig_unit,file=trim(seedname)//'_sconfig.dat',status='unknown',form='formatted',action='write')

       call io_date(cdate,ctime)
       write(sconfig_unit,*) 'Created on '//cdate//' at '//ctime ! Date and time
       write(sconfig_unit,'(3I6)') supercell_size(1),supercell_size(2),supercell_size(3)
       write(sconfig_unit,'(2I6)') num_total_atoms,tems_num

       !write(stdout,*)
       write(stdout,'(/1x,a)')'*--------------------------------- MONTECARLO --------------------------------------*'
       write(stdout,'(1x,a)') '|     Temp      Total Steps    Rejected Steps    Accepted Steps   Acc/Tot Ratio     |'
       write(stdout,'(1x,a)') '+-----------------------------------------------------------------------------------+'

       mc_unit = io_file_unit()
       open(mc_unit,file=trim(seedname)//'_mc.dat',status='unknown',form='formatted',action='write')
       call io_date(cdate,ctime)
       write(mc_unit,*) '#written on '//cdate//' at '//ctime ! Date and time
       write(mc_unit,'(1x,a)') '#------------------------------------------------- MONTECARLO --------'&
                                      //'-----------------------------------------#'
       write(mc_unit,'(1x,a)') '# Temp      Magnetization     Energy_ave          C_M              Sus'&
                                      //'              U_E              U_M       #'
       write(mc_unit,'(1x,a)') '#-----------------------------------------------------------------------'&
                                   //'---------------------------------------#'
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

       if(lorderparam) then
         orderparam_unit = io_file_unit()
         open(orderparam_unit,file=trim(seedname)//'_op.dat',status='unknown',form='formatted',action='write')
         call io_date(cdate,ctime)
         write(orderparam_unit,*) '#written on '//cdate//' at '//ctime ! Date and time
         write(orderparam_unit,'(1x,a)') '#------------------------------- MONTECARLO -----------------------------------#'
         write(orderparam_unit,'(1x,a)') '#    Temp              OP                 Sus_OP                 U_OP          #'
         write(orderparam_unit,'(1x,a)') '#------------------------------------------------------------------------------#'
       endif ! lorderparam

       if(lstaggered) then
         staggered_unit = io_file_unit()
         open(staggered_unit,file=trim(seedname)//'_staggered.dat',status='unknown',form='formatted',action='write')
         call io_date(cdate,ctime)
         write(staggered_unit,*) '#written on '//cdate//' at '//ctime ! Date and time
         write(staggered_unit,'(1x,a)') '#------------------------------- MONTECARLO --------------------------------#'
         write(staggered_unit,'(1x,a)') '#    Temp         Staggered_m          Sus_Staggered_m       U_Staggered_m  #'
         write(staggered_unit,'(1x,a)') '#---------------------------------------------------------------------------#'
       endif ! lstaggered 


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
          write(pm_unit,'(1x,a)') '#------------------------------------- MONTECARLO ---------'&
                                   //'----------------------------------#'
          write(pm_unit,'(1x,a)') '#       M or OP                 P(M)                     P('&
                                   //'OP)               P(Staggered_m)  #'
          write(pm_unit,'(1x,a)') '#----------------------------------------------------------'&
                                   //'----------------------------------#'
          if (lbinerror) then
             write(binerror_unit,'(1x,a)') '#------------------------------------ MONTECARLO ------'&
                                            //'------------------------------------#'
             write(binerror_unit,'(1x,a)') '#  Binning_level     Error(E)          Error(M)        '&
                                            //'Error(OP)       Error(Staggered_m)  #'
             write(binerror_unit,'(1x,a)') '#------------------------------------------------------'&
                                           //' ------------------------------------#'
          endif
       elseif (lorderparam) then
          write(pm_unit,'(1x,a)') '#-------------------------- MONTECARLO -------------------------------#'
          write(pm_unit,'(1x,a)') '#       M or OP                 P(M)                     P(OP)        #'
          write(pm_unit,'(1x,a)') '#---------------------------------------------------------------------#'
          if (lbinerror) then
             write(binerror_unit,'(1x,a)') '#------------------------- MONTECARLO -------------------------------#'
             write(binerror_unit,'(1x,a)') '#  Binning_level     Error(E)          Error(M)         Error(OP)    #'
             write(binerror_unit,'(1x,a)') '#--------------------------------------------------------------------#'
          endif
       elseif (lstaggered) then
          write(pm_unit,'(1x,a)') '#-------------------------- MONTECARLO -----------------------------#'
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

       if (energy_write) then
         call open_E_unit(steps_mc,energy_unit)
       endif  
    end if ! on_root

    call mcarlo()

    ! Before ending, I deallocate memory
    call mcarlo_dealloc()

    if (on_root)                   close(mc_unit)
    if (on_root)                   close(pm_unit)
    if (on_root)                   close(sconfig_unit)
    if (on_root .and. lbinerror)   close(binerror_unit)
    if (on_root .and. lspincorr)   close(spincorr_unit)
    if (on_root .and. lorderparam) close(orderparam_unit)
    if (on_root .and. lstaggered)  close(staggered_unit)
       
       
    return

   contains

    !============================================!
    subroutine mcarlo_setup() 
      !============================================!

 
      use mc_jij, only    : mag_moments_matrix

      implicit none

      character(len=50)          :: filename 
      integer                    :: loop,loop_nodes
      integer                    :: ierr,ifreq
      integer                    :: localidx,globalidx,globalidx_root
      real(kind=dp)              :: r1,r2
      real(kind=dp)              :: m_min,m_max
      real(kind=dp), allocatable :: global_spin_matrix_file(:,:,:)

      if (lorderparam.and.lstaggered) then
         ndim_data = 12; ndim_pm = 3; ndim_sum = 4
         SM = 10; SM2 = 11; SM4 = 12
      elseif (lorderparam.or.lstaggered) then
         if (lstaggered) SM = 7; SM2 = 8; SM4 = 9
         ndim_data = 9 ; ndim_pm = 2; ndim_sum = 3
      else
         ndim_data = 6 ; ndim_pm = 1; ndim_sum = 2
      endif

      ! Initiale configuration
      allocate(local_spin_matrix(3,num_total_atoms), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating spin_matrix in mcarlo_setup')
      if (index(initial_sconfig,'ferro')>0 ) then
         local_spin_matrix(1,:)=zero
         local_spin_matrix(2,:)=zero
         local_spin_matrix(3,:)=1.0_dp*mag_moments_matrix(:)
      elseif (index(initial_sconfig,'rand')>0 ) then
         do loop=1,num_total_atoms

            r1=mtprng_rand_real2(state)
!            call random_number(r1)

            r2=mtprng_rand_real2(state)
!            call random_number(r2)
            local_spin_matrix(1,loop)=mag_moments_matrix(loop)*sin(r1*pi)*cos(r2*twopi)
            local_spin_matrix(2,loop)=mag_moments_matrix(loop)*sin(r1*pi)*sin(r2*twopi)
            local_spin_matrix(3,loop)=mag_moments_matrix(loop)*cos(r1*pi)
          enddo
      elseif (index(initial_sconfig,'file')>0 ) then
         allocate(local_spin_matrix_file(3,num_total_atoms,int(tems_num/num_nodes)), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating spin_matrix in mcarlo_setup')
         ! Root node read the total spin matrix
         if (on_root) then
            allocate(global_spin_matrix_file(3,num_total_atoms,tems_num), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_spin_matrix_file in mcarlo_setup')
            filename = trim(seedname)//'_sconfig.dat'
            !! Output is the normalized(spin_matrix)
            call mcarlo_read_sconfig(num_total_atoms,tems,global_spin_matrix_file,filename,tems_num)
         else
            allocate(global_spin_matrix_file(1,1,1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_spin_matrix_file in mcarlo_setup')
         endif ! root

         ! Send the spin matrix each node
         do localidx=1, int(tems_num/num_nodes)
            globalidx = (localidx-1)*num_nodes+my_node_id+1
            if (on_root) then
               local_spin_matrix_file(:,:,localidx) = global_spin_matrix_file(:,:,globalidx)
               do loop_nodes=1,num_nodes-1
                  globalidx_root = globalidx-root_id+loop_nodes
                  call comms_send(global_spin_matrix_file(1,1,globalidx_root),3*num_total_atoms,loop_nodes)
               enddo
            else
              call comms_recv(local_spin_matrix_file(1,1,localidx),3*num_total_atoms*1,root_id)
            endif
              local_spin_matrix_file(1,:,localidx) = local_spin_matrix_file(1,:,localidx)*mag_moments_matrix(:)
              local_spin_matrix_file(2,:,localidx) = local_spin_matrix_file(2,:,localidx)*mag_moments_matrix(:)
              local_spin_matrix_file(3,:,localidx) = local_spin_matrix_file(3,:,localidx)*mag_moments_matrix(:)
         enddo
         if (allocated(global_spin_matrix_file)) then
            deallocate(global_spin_matrix_file, stat=ierr )
            if (ierr/=0) call io_error('Error in deallocating global_spin_matrix_file in mcarlo_setup')
         endif
      end if

      allocate(local_data(ndim_data), stat=ierr )
      if (ierr/=0) call io_error('Error in allocating local_data in mcarlo_setup')
      if (lbinerror) then
         allocate(local_sum(ndim_sum,num_binning_level), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating local_sum in mcarlo_setup')
         allocate(local_binning(ndim_sum,num_binning_level), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating local_binning in mcarlo_setup')
      endif 

      if (lspincorr) then
         allocate(local_spincorr(jij_shell_max,num_species,num_species), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating local_spincorr in mcarlo_setup')
         allocate(local_spincorr_abs(jij_shell_max,num_species,num_species), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating local_spincorr_abs in mcarlo_setup')
         allocate(local_spincorr_tot_abs(jij_shell_max,num_species,num_species), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating local_spincorr_tot_abs in mcarlo_setup')
      endif

      if (lsfactor) then
         allocate(local_sfactor_matrix(3,sfactor_nqpts), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating local_sfactor_matrix in mcarlo_setup')
         local_sfactor_matrix=0.0_dp
         allocate(atoms_supercell_pos_cart(3,num_total_atoms), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating atoms_supercell_pos_cart in mcarlo_setup')
         atoms_supercell_pos_cart=0.0_dp
         call mcarlo_get_atoms_pos(atoms_supercell_pos_cart)
      endif

      if (energy_write) then
         allocate(local_E(0:energy_num_print-1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating local_E in mcarlo_setup')
      endif
 
      if (on_root) then
         allocate(global_spin_matrix(3,num_total_atoms,0:num_nodes-1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_spin_matrix in mcarlo_setup')
         allocate(global_tem(0:num_nodes-1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_tem in mcarlo_setup')
         allocate(global_rejected(2,0:num_nodes-1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_rejected in mcarlo_setup')
         allocate(global_energy(3,0:num_nodes-1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_energy in mcarlo_setup')
         if (lspincorr) then
            allocate(global_spincorr(jij_shell_max,num_species,num_species,0:num_nodes-1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_spincorr in mcarlo_setup')
            allocate(global_spincorr_abs(jij_shell_max,num_species,num_species,0:num_nodes-1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_spincorr_abs in mcarlo_setup')
            allocate(global_spincorr_tot_abs(jij_shell_max,num_species,num_species,0:num_nodes-1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_spincorr_tot_abs in mcarlo_setup')
         endif
     
         if (lbinerror) then
            allocate(global_binning(ndim_sum,num_binning_level,0:num_nodes-1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_binning mcarlo_setup')
         endif
     
         allocate(global_data(ndim_data,0:num_nodes-1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_data mcarlo_setup')
     
         if (lsfactor) then
            allocate(global_sfactor_matrix(3,sfactor_nqpts,0:num_nodes-1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_sfactor_matrix in mcarlo_setup')
            global_sfactor_matrix=0.0_dp
         endif
         if (energy_write) then
            allocate(energy_unit(tems_num), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating energy_unit in mcarlo_setup')
            allocate(global_E(0:energy_num_print-1,0:num_nodes-1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_E in mcarlo_setup')
         endif
      else
         allocate(global_spin_matrix(1,1,1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_spin_matrix in mcarlo_setup')
         allocate(global_tem(1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_tem in mcarlo_setup')
         allocate(global_rejected(1,1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_rejected in mcarlo_setup')
         allocate(global_energy(1,1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_energy in mcarlo_setup')
         allocate(global_data(1,1), stat=ierr )
         if (ierr/=0) call io_error('Error in allocating global_data mcarlo_setup')
         if (lspincorr) then
            allocate(global_spincorr(1,1,1,1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_spincorr in mcarlo_setup')
            allocate(global_spincorr_abs(1,1,1,1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_spincorr_abs in mcarlo_setup')
            allocate(global_spincorr_tot_abs(1,1,1,1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_spincorr_tot_abs in mcarlo_setup')
         endif
         if (lbinerror) then
            allocate(global_binning(1,1,1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_binning mcarlo_setup')
         endif
         if (lsfactor) then
            allocate(global_sfactor_matrix(1,1,1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_sfactor_matrix in mcarlo_setup')
            global_sfactor_matrix = 0.0_dp
         endif
         if (energy_write) then
            allocate(energy_unit(1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating energy_unit in mcarlo_setup')
            allocate(global_E(1,1), stat=ierr )
            if (ierr/=0) call io_error('Error in allocating global_E in mcarlo_setup')
         endif

      endif

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
      if (ierr/=0) call io_error('Error in allocating m_array in mcarlo subroutine')
      allocate(local_pm(ndim_pm,num_freq),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating local_pm in mcarlo subroutine')
      if(on_root)then
        allocate(global_pm(ndim_pm,num_freq,0:num_nodes-1),stat=ierr)
        if (ierr/=0) call io_error('Error in allocating global_pm in mcarlo subroutine')
      else
        allocate(global_pm(1,1,1),stat=ierr)
        if (ierr/=0) call io_error('Error in allocating global_pm in mcarlo subroutine')
      endif

      do ifreq=1,num_freq
         m_array(ifreq) = m_min+real(ifreq-1,dp)*d_m
      enddo

      return
 
    end subroutine mcarlo_setup

    !==================================================================!
    subroutine mcarlo()
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  

      implicit none
      integer(8)           :: steps_mc_tot,num_sfsteps_tot
      integer(8)           :: steps_warmup_tot
      integer              :: loop_nodes       
      integer              :: loop_at,loop_mc,loop_msure
      integer              :: i,nsp1,nsp2
      integer              :: m_l,l,excc
      integer              :: localidx,globalidx
      integer              :: ifreq
      real(kind=dp)        :: local2_oc,beta
      real(kind=dp)        :: specific_heat,susceptibility
      real(kind=dp)        :: U_E,U_M,Delta_E,Delta_M,normalize
      real(kind=dp)        :: sus_orderparam,U_orderparam
      real(kind=dp)        :: sus_staggered,U_staggered
      integer(8)           :: local_rejected(2)
      real(kind=dp)        :: local_energy(3),local_tem
      integer              :: num_at     

      if(on_root) call io_stopwatch('mcarlo_main',1)

      do localidx=1, int(tems_num/num_nodes)
         local_data     = zero; local_energy = zero; local_pm = zero
         local_rejected = 0
         if (lsfactor) local_sfactor_matrix = zero
         if (lbinerror) then
            local_sum = zero; local_binning = zero
         endif
         if (lspincorr) then
            local_spincorr = zero; local_spincorr_abs = zero ; local_spincorr_tot_abs = zero
         endif
         globalidx      = (localidx-1)*num_nodes+my_node_id+1
         beta = 1.0_dp/tems(globalidx)
         local_tem = tems(globalidx)
         if (index(initial_sconfig,'file')>0 ) then
            local_spin_matrix = local_spin_matrix_file(:,:,localidx)
         endif
         !! Initial energy 
         call energy_get(local_spin_matrix,local_energy(2))
         local_energy(1) = local_energy(2) 
         !! Thermalization
         steps_warmup_tot=0
         do loop_mc=1, steps_warmup   ! monte carlo steps
            do loop_msure=1, steps_measure
               do loop_at=1,num_total_atoms
                  !num_at = num_total_atoms*mtprng_rand_real2(state)+1
                  num_at = loop_at
                  call metropolis(num_at,beta,local_spin_matrix,&
                             local_energy(2),local_rejected(1))
                  steps_warmup_tot=steps_warmup_tot+1
               enddo ! atom
            enddo    ! measure
         enddo ! mc

         !! Monte Carlo
         steps_mc_tot=0
         num_sfsteps_tot=0
         do loop_mc=1, steps_mc   ! monte carlo steps
            do loop_msure=1, steps_measure
               do loop_at=1,num_total_atoms
                  !num_at = num_total_atoms*mtprng_rand_real2(state)+1
                  num_at = loop_at
                  call metropolis(num_at,beta,local_spin_matrix,&
                                 local_energy(2),local_rejected(2))
                  steps_mc_tot = steps_mc_tot+1
               enddo ! atom
            enddo    ! measure
            local_data(E)  = local_data(E)+local_energy(2)
            local_data(E2) = local_data(E2)+local_energy(2)**2
            local_data(E4) = local_data(E4)+local_energy(2)**4

            call magnetization_get(local_spin_matrix,local2_oc)
            local_data(M)  = local_data(M)+sqrt(local2_oc)
            local_data(M2) = local_data(M2)+local2_oc
            local_data(M4) = local_data(M4)+local2_oc**2

            if (lspincorr) call spincorr_get(local_spin_matrix,local_spincorr,&
                                             local_spincorr_abs,local_spincorr_tot_abs)

            if (lbinerror) then
              local_sum(1,:) = local_sum(1,:)+(local_energy(2)/real(num_total_atoms,dp))
              local_sum(2,:) = local_sum(2,:)+sqrt(local2_oc)            
            endif

            if (lorderparam) then
               call orderparam_get(local_spin_matrix,local2_oc)
               local_data(OP)  = local_data(OP)+sqrt(local2_oc)
               local_data(OP2) = local_data(OP2)+local2_oc
               local_data(OP4) = local_data(OP4)+local2_oc**2
               if (lbinerror)  local_sum(3,:) = local_sum(3,:)+sqrt(local2_oc)         
            endif
 
            if (lstaggered) then
               call staggered_get(local_spin_matrix,local2_oc)
               local_data(SM)  = local_data(SM)+sqrt(local2_oc)
               local_data(SM2) = local_data(SM2)+local2_oc
               local_data(SM4) = local_data(SM4)+local2_oc**2
               ! ndim_sum = 4 If lorderparam = True
               ! ndim_sum = 3  when lorderparam = False 
               if(lbinerror)  local_sum(ndim_sum,:) = local_sum(ndim_sum,:)+sqrt(local2_oc)         
            endif

            if (lorderparam.and.lstaggered) then
               call Pm_get(local_spin_matrix,num_freq,local_pm(1,:),pdf_order=local_pm(2,:),&
                           pdf_staggered=local_pm(3,:))
            elseif (lorderparam) then 
               call Pm_get(local_spin_matrix,num_freq,local_pm(1,:),pdf_order=local_pm(2,:))
            elseif (lstaggered) then
               call Pm_get(local_spin_matrix,num_freq,local_pm(1,:),pdf_staggered=local_pm(2,:))
            else 
               call Pm_get(local_spin_matrix,num_freq,local_pm(1,:))
            endif
 
            if (lbinerror) then
               do l=1,num_binning_level
                  if (mod(loop_mc,(2**(binning_level(l)))).eq. 0) then
                     local_sum(:,l)     = local_sum(:,l)/real(2**(binning_level(l)),dp)
                     local_binning(:,l) = local_binning(:,l)+local_sum(:,l)**2
                     local_sum(:,l)     = zero
                  endif
               enddo
            endif

            if (lsfactor .and. mod(loop_mc,sfactor_steps_measure) .eq. 0.0) then
               call mcarlo_structure_factor(local_spin_matrix,&
                              atoms_supercell_pos_cart,local_sfactor_matrix)
               num_sfsteps_tot=num_sfsteps_tot+1
            endif

            if (energy_write) then
               call comms_barrier()
               local_E(mod(loop_mc,energy_num_print)) = local_energy(2)
               if (mod(loop_mc,energy_num_print).eq. 0) then
                  if (.not. on_root) &
                     call comms_send(local_E(0),energy_num_print,root_id)
                  if (on_root) then
                     do loop_nodes=1,num_nodes-1
                        call comms_recv(global_E(0,loop_nodes),energy_num_print,loop_nodes)
                     enddo
                     global_E(:,0) = local_E
                     call write_E(loop_mc,num_nodes,energy_num_print,global_E,&
                                  energy_unit((localidx-1)*num_nodes+1:localidx*num_nodes))
                     global_E = zero
                  endif
               endif
               ! This is for a last numbers less than energy_num_print
               if ((loop_mc .eq. steps_mc) .and. &
                   (steps_mc .gt. ((steps_mc/energy_num_print)*energy_num_print))) then
                  excc = steps_mc-(steps_mc/energy_num_print)*energy_num_print  
                  local_E(0) = local_energy(2)
                  if (.not. on_root) &
                     call comms_send(local_E(0),excc,root_id)
                  if (on_root) then
                     do loop_nodes=1,num_nodes-1
                        call comms_recv(global_E(0,loop_nodes),excc,loop_nodes)
                     enddo
                     global_E(0:excc-1,0)=local_E(0:excc-1)
                     call write_E(steps_mc,num_nodes,excc,global_E(0:excc-1,0:num_nodes-1),&
                                  energy_unit((localidx-1)*num_nodes+1:localidx*num_nodes))
                     global_E = zero
                  endif
               endif
            endif

         enddo !mc

         if (lbinerror) then
            !! If (mod(steps_mc,2^l) .ne. 0) i.e. the last binning box is incomplete
            !! the local_sum is non zero and I should remove them from the average <E> or <M> ,...
            local_sum(1,:) = (local_data(E)/real(num_total_atoms,dp))-local_sum(1,:)
            local_sum(2,:) = local_data(M)-local_sum(2,:)
            if (lorderparam.and.lstaggered) then
               local_sum(3,:) = local_data(OP)-local_sum(3,:)
               local_sum(4,:) = local_data(SM)-local_sum(4,:)
            elseif (lorderparam.or.lstaggered) then
               ! if lorderparam= T and lstaggered= F, OP2=8
               ! if lorderparam= F and lstaggered= T, SM2=OP2=8
               local_sum(3,:) = local_data(OP)-local_sum(3,:)
            endif

            do l=1,num_binning_level
               m_l = (steps_mc)/(2**(binning_level(l)))
               local_sum(:,l) = local_sum(:,l)/real((2**(binning_level(l))),dp)
               local_sum(:,l) = local_sum(:,l)/real(m_l,dp)
               local_binning(:,l) = local_binning(:,l)/real(m_l,dp)
               local_binning(:,l) = (local_binning(:,l)-local_sum(:,l)**2)/real((m_l-1),dp)
            enddo
         endif

         ! I calculate the final Energy by direct method and save at local_energy(3,:)
         call energy_get(local_spin_matrix,local_energy(3))
 
         ! Now I gather all data in monte carlo step on 
         ! root node(local_data,local_energy,local_rejected,local_sfactor_matrix)
         !NOTE****** I changed the sfactor_matrix dimension **************
         !*******from (0:product()-1,:) to (product(),:)  ************
!         call comms_barrier()
         if (.not.on_root) then
            call comms_send(local_rejected(1),2,root_id)
            call comms_send(local_energy(1),3,root_id)
            call comms_send(local_tem,1,root_id)
            call comms_send(local_spin_matrix(1,1),3*num_total_atoms,root_id)
            call comms_send(local_data(1),ndim_data,root_id)
            call comms_send(local_pm(1,1),ndim_pm*num_freq,root_id)

            if (lspincorr) then
               call comms_send(local_spincorr(1,1,1),jij_shell_max*num_species*num_species,root_id)
               call comms_send(local_spincorr_abs(1,1,1),jij_shell_max*num_species*num_species,root_id)
               call comms_send(local_spincorr_tot_abs(1,1,1),jij_shell_max*num_species*num_species,root_id)
            endif
            if (lsfactor) &
               call comms_send(local_sfactor_matrix(1,1),3*sfactor_nqpts,root_id)
            if (lbinerror) &
               call comms_send(local_binning(1,1),ndim_sum*num_binning_level,root_id)
         endif

         if (on_root) then
            do loop_nodes=1,num_nodes-1
               call comms_recv(global_rejected(1,loop_nodes),2,loop_nodes)
               call comms_recv(global_energy(1,loop_nodes),3,loop_nodes)
               call comms_recv(global_tem(loop_nodes),1,loop_nodes)
               call comms_recv(global_spin_matrix(1,1,loop_nodes),3*num_total_atoms,loop_nodes)
               call comms_recv(global_data(1,loop_nodes),ndim_data,loop_nodes)
               call comms_recv(global_pm(1,1,loop_nodes),ndim_pm*num_freq,loop_nodes)
               if (lspincorr) then
                  call comms_recv(global_spincorr(1,1,1,loop_nodes),jij_shell_max*num_species*num_species,loop_nodes)
                  call comms_recv(global_spincorr_abs(1,1,1,loop_nodes),jij_shell_max*num_species*num_species,loop_nodes)
                  call comms_recv(global_spincorr_tot_abs(1,1,1,loop_nodes),jij_shell_max*num_species*num_species,loop_nodes)
               endif

               if (lsfactor) call comms_recv(global_sfactor_matrix(1,1,loop_nodes),3*sfactor_nqpts,loop_nodes)
               if (lbinerror) call comms_recv(global_binning(1,1,loop_nodes),ndim_sum*num_binning_level,loop_nodes)
            enddo
            ! Until now I save all data of all node except root node on global variables,
            ! now I save the root data on global variables. 
            global_data(:,0)          = local_data
            global_pm(:,:,0)          = local_pm
            global_rejected(:,0)      = local_rejected
            global_energy(:,0)        = local_energy
            global_tem(0)            = local_tem
            global_spin_matrix(:,:,0) = local_spin_matrix
            if (lsfactor)  global_sfactor_matrix(:,:,0) = local_sfactor_matrix(:,:)
            if (lbinerror) global_binning(:,:,0) = local_binning

            normalize   = 1.0_dp/real(steps_mc,dp)
            global_data = global_data*normalize
            global_pm   = global_pm*normalize/real(num_supercell,dp)
            if (lspincorr) then
               global_spincorr(:,:,:,0) = local_spincorr
               global_spincorr_abs(:,:,:,0) = local_spincorr_abs
               global_spincorr_tot_abs(:,:,:,0) = local_spincorr_tot_abs
               global_spincorr = global_spincorr*normalize/real(num_total_atoms,dp)
               global_spincorr_abs = global_spincorr_abs*normalize/real(num_total_atoms,dp)
               global_spincorr_tot_abs = global_spincorr_tot_abs*normalize/real(num_total_atoms,dp)
            endif

            ! write the results 
            do loop_nodes=0,num_nodes-1                                         
               beta = 1.0_dp/global_tem(loop_nodes)
               specific_heat = beta*beta*(global_data(E2,loop_nodes)-global_data(E,loop_nodes)**2)/real(num_total_atoms,dp)
               susceptibility = beta*(global_data(M2,loop_nodes)-global_data(M,loop_nodes)**2)*num_total_atoms
               U_E = 1.0_dp-(1.0_dp/3.0_dp)*(global_data(E4,loop_nodes)/global_data(E2,loop_nodes)**2)
               U_M = 1.0_dp-(1.0_dp/3.0_dp)*(global_data(M4,loop_nodes)/global_data(M2,loop_nodes)**2)
               Delta_E = sqrt((global_data(E2,loop_nodes)-global_data(E,loop_nodes)**2)*normalize)/real(num_total_atoms,dp)
               Delta_M = sqrt((global_data(M2,loop_nodes)-global_data(M,loop_nodes)**2)*normalize)

               if (steps_warmup>0) then 
                  write(stdout,124) global_tem(loop_nodes),steps_warmup_tot,global_rejected(1,loop_nodes),&
                        steps_warmup_tot-global_rejected(1,loop_nodes),&                   
                        (steps_warmup_tot-global_rejected(1,loop_nodes))/real(steps_warmup_tot,dp),'  <-- WARMUP'
               else
                  write(stdout,'(1x,85a)') ' No Warm-up step is done! ' 
               endif  

               write(stdout,124) global_tem(loop_nodes),steps_mc_tot,global_rejected(2,loop_nodes),&
                     steps_mc_tot-global_rejected(2,loop_nodes),&                   
                     (steps_mc_tot-global_rejected(2,loop_nodes))/real(steps_mc_tot,dp),'  <-- MCARLO'
               write(stdout,'(4x,a,E15.7,a)')  ' Initial Energy                           = ',&
                                                 global_energy(1,loop_nodes),'  <-- ENERGY'
               write(stdout,'(4x,a,E15.7,a)')  ' Final Energy by adding differences       = ',&
                                                 global_energy(2,loop_nodes),'  <-- ENERGY'
               write(stdout,'(4x,a,E15.7,a)')  ' Final Energy by direct calculations      = ',&
                                                 global_energy(3,loop_nodes),'  <-- ENERGY'
               write(stdout,'(4x,a,f15.9)')    '                             Delta E      = +/- ',Delta_E
               write(stdout,'(4x,a,f15.9)')    '                             Delta M      = +/- ',Delta_M
               write(stdout,'(1x,a85)') repeat('-',85)                           

               write(mc_unit,125) global_tem(loop_nodes),global_data(M,loop_nodes),&
                     global_data(E,loop_nodes)/real(num_total_atoms,dp),specific_heat,susceptibility,U_E,U_M

               if (lspincorr) then
                  do l=1,jij_shell_max
                     do nsp1=1,num_species
                        do nsp2=1,num_species
                           write(spincorr_unit,'(f9.4,1x,I5,6x,a3,4x,i3,6x,a3,4x,i3,8x,f12.9,10x,f12.9,12x,f12.9)')  &
                                 global_tem(loop_nodes),l,atoms_symbol(nsp1),nsp1,atoms_symbol(nsp2),nsp2,&
                                 global_spincorr(l,nsp1,nsp2,loop_nodes),global_spincorr_abs(l,nsp1,nsp2,loop_nodes),&
                                 global_spincorr_tot_abs(l,nsp1,nsp2,loop_nodes)
                        enddo
                     enddo
                  enddo
               endif

               if (lorderparam) then
                  sus_orderparam = beta*(global_data(OP2,loop_nodes)-global_data(OP,loop_nodes)**2)*num_total_atoms
                  U_orderparam = 1.0_dp-(1.0_dp/3.0_dp)*(global_data(OP4,loop_nodes)/global_data(OP2,loop_nodes)**2)
                  write(orderparam_unit,126) global_tem(loop_nodes),global_data(OP,loop_nodes),sus_orderparam,U_orderparam
               endif
 
               if (lstaggered) then
                  sus_staggered = beta*(global_data(SM2,loop_nodes)-global_data(SM,loop_nodes)**2)*num_total_atoms
                  U_staggered   = 1.0_dp-(1.0_dp/3.0_dp)*(global_data(SM4,loop_nodes)/global_data(SM2,loop_nodes)**2)
                  write(staggered_unit,126) global_tem(loop_nodes),global_data(SM,loop_nodes),sus_staggered,U_staggered
               endif
 
              if (lbinerror) then
                 write(binerror_unit,'(1x,a,2x,f9.4)') '#T=',global_tem(loop_nodes)
                 do l=1,num_binning_level
                    if (lorderparam.and.lstaggered) then
                       write(binerror_unit,'(2x,I3,8x,E15.8,3x,E15.8,3x,E15.8,3x,E15.8)') binning_level(l),&
                            (sqrt(global_binning(i,l,loop_nodes)/real(num_total_atoms,dp)),i=1,4)
                    elseif (lorderparam.or.lstaggered) then
                       write(binerror_unit,'(2x,I8,8x,E15.8,3x,E15.8,3x,E15.8)') binning_level(l),&
                            (sqrt(global_binning(i,l,loop_nodes)/real(num_total_atoms,dp)),i=1,3)
                    else
                       write(binerror_unit,'(2x,I8,8x,E15.8,3x,E15.8)') binning_level(l), &
                            (sqrt(global_binning(i,l,loop_nodes)/real(num_total_atoms,dp)),i=1,2)
                    endif
                 enddo
              endif

              write(pm_unit,'(1x,a,2x,f9.4)') '#T=',global_tem(loop_nodes)
     
              do ifreq=1,num_freq
                 if (lorderparam.and.lstaggered) then
                    write(pm_unit,'(2x,E15.8,10x,E15.8,10x,E15.8,10x,E15.8)') m_array(ifreq),&
                          (global_pm(i,ifreq,loop_nodes),i=1,3)
                 elseif (lorderparam.or.lstaggered) then
                    write(pm_unit,'(2x,E15.8,10x,E15.8,10x,E15.8)') m_array(ifreq),(global_pm(i,ifreq,loop_nodes),i=1,2)
                 else
                    write(pm_unit,'(2x,E15.8,10x,E15.8)') m_array(ifreq),global_pm(1,ifreq,loop_nodes)
                 endif
              enddo

              write(sconfig_unit,'(F12.6)') global_tem(loop_nodes)
              do l=1,num_total_atoms
                 write(sconfig_unit,'(3F12.6)') (global_spin_matrix(i,l,loop_nodes),i=1,3)
              enddo

            enddo ! loop_nodes                                                            

          if (lsfactor) then
             global_sfactor_matrix = global_sfactor_matrix/real(num_total_atoms,dp)
             global_sfactor_matrix = global_sfactor_matrix/real(num_sfsteps_tot,dp)
             call mcarlo_structure_factor_write_mc(global_sfactor_matrix,localidx)
             if (localidx .eq. int(tems_num/num_nodes)) then
                write(stdout,'(4x,a,I14,17x,a)')  ' Total used steps at S Factor calculation = ',num_sfsteps_tot,'   <-- SFAC'
                write(stdout,'(1x,a85)') repeat('-',85)
             endif
          endif

         endif  ! on_root                                                              
 
       enddo ! localidx

       
       if (on_root) then
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
          write(stdout,'(3X,A)') "Transport properties calculated."

          if(energy_write) call close_E_unit(energy_unit)
       endif

       if(on_root) call io_stopwatch('mcarlo_main',2)

124   format(3x,f9.4,1x,i14,1x,i14,5x,i14,7x,f10.8,1x,a)
125   format(f9.4,2x,E15.8,2x,E15.8,2x,E15.8,2x,E15.8,2x,E15.8,2x,E15.8)
126   format(3x,f9.4,5x,E15.8,6x,E15.8,7x,E15.8,5x)

      return

    end subroutine mcarlo

    !==================================================================!
    subroutine mcarlo_dealloc()
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  

      implicit none

      ! Deallocate real arrays that are private

      deallocate(global_spin_matrix, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating global_spin_matrix in mcarlo_dealloc')
      deallocate(local_spin_matrix, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating local_spin_matrix in mcarlo_dealloc')
 
      if(allocated(local_spin_matrix_file)) then
       deallocate(local_spin_matrix_file, stat=ierr )
       if (ierr/=0) call io_error('Error in deallocating local_spin_matrix_file in mcarlo_dealloc')
      endif

      deallocate(global_tem, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating spin_matrix in mcarlo_dealloc')

      deallocate(global_energy, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating local_energy in mcarlo_dealloc')

      deallocate(global_rejected, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating local_rejected in mcarlo_dealloc')

      deallocate(global_data, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating global_data in mcarlo_dealloc')
      deallocate(local_data, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating local_data in mcarlo_dealloc')

      deallocate(global_pm, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating global_pm in mcarlo_dealloc')
      deallocate(local_pm, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating local_pm in mcarlo_dealloc')

      deallocate(m_array, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating m_array in mcarlo_dealloc')

      if (lspincorr) then
        deallocate(global_spincorr, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating global_spincorr in mcarlo_dealloc')
        deallocate(local_spincorr, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating local_spincorr in mcarlo_dealloc')
        deallocate(global_spincorr_abs, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating global_spincorr_abs in mcarlo_dealloc')
        deallocate(local_spincorr_abs, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating local_spincorr_abs in mcarlo_dealloc')
        deallocate(global_spincorr_tot_abs, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating global_spincorr_tot_abs in mcarlo_dealloc')
        deallocate(local_spincorr_tot_abs, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating local_spincorr_tot_abs in mcarlo_dealloc')
      endif

      if (lbinerror) then
        deallocate(global_binning, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating global_binning in mcarlo_dealloc')
        deallocate(local_binning, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating local_binning in mcarlo_dealloc')
        deallocate(local_sum, stat=ierr )
        if (ierr/=0) call io_error('Error in deallocating local_sum in mcarlo_dealloc')
      endif

      if (lsfactor) then
         deallocate( global_sfactor_matrix, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating global_sfactor_matrix in mcarlo_dealloc')
         deallocate( local_sfactor_matrix, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating local_sfactor_matrix in mcarlo_dealloc')
         deallocate( atoms_supercell_pos_cart, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating atoms_supercell_pos_cart in mcarlo_dealloc')
      end if

      if (energy_write) then
         deallocate( energy_unit, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating energy_unit in mcarlo_dealloc')
         deallocate( global_E, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating global_E in mcarlo_dealloc')
         deallocate( local_E, stat=ierr  )
         if (ierr/=0) call io_error('Error in deallocating local_E in mcarlo_dealloc')
      end if

      return

    end subroutine mcarlo_dealloc

  end subroutine mcarlo_main


end module mc_mcarlo
