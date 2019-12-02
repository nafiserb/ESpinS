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

module mc_parameters

  use mc_constants, only : dp
  use mc_io,        only : stdout,maxlen
  use mc_comms,     only : num_nodes
  use stdtypes
  use mtprng,       only : mtprng_state


  implicit none

  private

  !Input
  character(len=20),             public, save :: length_unit
  character(len=20),             public, save :: param_unit
  character(len=20),             public, save :: coordinate
  ! command
  logical,                       public, save :: inputfile_setup1
  logical,                       public, save :: inputfile_setup2
  logical,                       public, save :: hamfile_setup

 ! Atom sites
  character(len=maxlen), allocatable, public, save :: atoms_label(:)
  real(kind=dp),                 public, save :: real_lattice(3,3)
  real(kind=dp),    allocatable, public, save :: atoms_pos_frac(:,:,:)
  real(kind=dp),    allocatable, public, save :: atoms_pos_cart(:,:,:)
  integer,          allocatable, public, save :: atoms_species_num(:)
  character(len=3), allocatable, public, save :: atoms_symbol(:)
  integer,                       public, save :: num_species
  integer,                       public, save :: num_atoms

  real(kind=dp),                 public, save :: neighbors_tol
  integer,       allocatable,    public, save :: shells(:,:)           !shells_jij/bij/dij  save in this variable
  real(kind=dp), allocatable,    public, save :: parameters_sigma(:,:,:)  ! if Spin_glass is true.
  real(kind=dp), allocatable,    public, save :: parameters_jij(:,:,:)
  real(kind=dp), allocatable,    public, save :: parameters_bij(:,:,:)
  real(kind=dp), allocatable,    public, save :: parameters_dij(:,:,:)

  integer,                       public, save :: num_total_atoms
  integer,                       public, save :: supercell_size(3)
  integer,                       public, save :: num_supercell
  integer,                       public, save :: steps_mc
  integer,                       public, save :: steps_measure
  integer,                       public, save :: steps_warmup  
                         
  integer,                       public, save :: tems_num 
  character(len=20)                           :: tems_mode
  real(kind=dp)                               :: tem_start
  real(kind=dp)                               :: tem_end 
  real(kind=dp)                               :: tem_steps
  real(kind=dp)                               :: tem_inv_steps
  real(kind=dp), allocatable,    public, save :: tems(:)        

  !Module PARALLEL-TEMPERING
  logical,                       public, save :: pt
  integer,                       public, save :: pt_steps_swap
  integer,                       public, save :: pt_print_swap 

  character(len=20),             public, save :: initial_sconfig    
  character(len=20),             public, save :: mcarlo_mode    
  real(kind=dp),                 public, save :: tilt_angles_max 

  ! Boundary 
  character(len=20)                           :: boundary
  logical,                       public, save :: lboundary_open

  ! Random numbers
  logical,                       public, save :: lseed
  integer, allocatable,          public, save :: seeds(:) 
  type(mtprng_state),            public, save :: state

  ! SPin Glass
  logical,                       public, save :: spin_glass
  integer,                       public, save :: spin_glass_seed 
  logical,                       public, save :: lspin_glass_seed
  real(kind=dp), allocatable,    public, save :: spin_glass_sigma(:)

  ! Atom Neighbors
  !Heisenberg term
  integer,                        public, save :: jij_num_nbors
  real(kind=dp), allocatable,     public, save :: jij_1atom_site(:,:)
  real(kind=dp), allocatable,     public, save :: jij_2atom_site(:,:)
  real(kind=dp), allocatable,     public, save :: jij(:)
  real(kind=dp), allocatable,     public, save :: mag_moments(:,:)
  integer      , allocatable,     public, save :: jij_shell(:)
  integer                   ,     public, save :: jij_shell_max

  ! Biquadratic term
  logical,                        public, save :: have_biquad  
  integer,                        public, save :: bij_num_nbors
  real(kind=dp), allocatable,     public, save :: bij_1atom_site(:,:)
  real(kind=dp), allocatable,     public, save :: bij_2atom_site(:,:)
  real(kind=dp), allocatable,     public, save :: bij(:)

  ! DM term
  logical,                        public, save :: have_dm      
  integer,                        public, save :: dij_num_nbors
  real(kind=dp), allocatable,     public, save :: dij_1atom_site(:,:)
  real(kind=dp), allocatable,     public, save :: dij_2atom_site(:,:)
  real(kind=dp), allocatable,     public, save :: dij(:)
  real(kind=dp), allocatable,     public, save :: dij_vectors_cart(:,:)
  real(kind=dp), allocatable,     public, save :: dij_vectors_frac(:,:)

  !Single-Ion term
  logical,                        public, save :: have_singleion  
  real(kind=dp),    allocatable,  public, save :: singleion_parameters(:)
  !real(kind=dp),    allocatable,  public, save :: parameters_singleion(:)
  real(kind=dp),    allocatable,  public, save :: singleion_axes_cart(:,:)
  character(len=3), allocatable,  public, save :: singleion_atoms_symbol(:)

  !H field         
  logical,                        public, save :: have_field
  real(kind=dp),    allocatable,  public, save :: field_parameters(:)
  real(kind=dp),    allocatable,  public, save :: field_axes_cart(:,:)
  character(len=3), allocatable,  public, save :: field_atoms_symbol(:)

  ! Order Parameter
  logical,                        public, save :: lorderparam
  real(kind=dp),    allocatable,  public, save :: orderparam_axes_cart(:,:)
  character(len=3), allocatable,  public, save :: orderparam_atoms_symbol(:)

  ! staggered magnetization
  logical,                        public, save :: lstaggered
  integer,          allocatable,  public, save :: staggered_coeff(:)

  ! Error from binnning Analysis 
  logical,                        public, save :: lbinerror
  integer,                        public, save :: num_binning_level
  integer,          allocatable,  public, save :: binning_level(:)

  ! Spin_correlations
  logical,                        public, save :: lspincorr
                                  
  ! write ENERGY                  
  logical,                        public, save :: energy_write
  integer,                        public, save :: energy_num_print
                                  
  ! Structure Factor             
  logical,                        public, save :: lsfactor
  logical,                        public, save :: lsfactor_polar
  real(kind=dp),                  public, save :: sfactor_corner(3)
  real(kind=dp),                  public, save :: sfactor_q1(3)
  real(kind=dp),                  public, save :: sfactor_q2(3)
  integer,                        public, save :: sfactor_2dqmesh(2)
  integer,                        public, save :: sfactor_nqpts     
  integer,                        public, save :: sfactor_steps_measure
  real(kind=dp),                  public, save :: sfactor_polar(3)

  !parameters dervied from input
  real(kind=dp),                  public, save :: recip_lattice(3,3)
  real(kind=dp),                  public, save :: cell_volume
  real(kind=dp),                  public, save :: lenconfac
  real(kind=dp),                  public, save :: paramconfac

  ! Private data
  integer                            :: num_lines
  character(len=maxlen), allocatable :: in_data(:)
  integer,               allocatable :: convert_matrix_atoms(:,:)
  logical                            :: inputfile_setup,file_setup         

  public :: param_read
  public :: param_write
  public :: param_write_seed
  public :: param_first_dealloc
  public :: param_second_dealloc
  public :: param_write_header
  public :: param_memory_estimate

contains

  !==================================================================!
  subroutine param_read ( )
    !==================================================================!
    !                                                                  !
    ! Read parameters and calculate derived values                     !
    !                                                                  !
    !===================================================================  
    use mc_constants, only : bohr,hart,eps4,eps8,k_B,bohr_magn,zero
    use mc_utility,   only : utility_recip_lattice,utility_frac_to_cart,&
                             utility_cart_to_frac,utility_sort,&
                             utility_cross
    use mc_io,        only : io_error,seedname,input1_file_flag,&
                             input2_file_flag,ham_file_flag
    implicit none

    !local variables
    logical                    :: found,found2,lunits,sconfig_found
    integer                    :: i_tmp,loop,ierr
    integer, allocatable       :: type_tmp(:,:),shells_tmp(:)
    real(kind=dp), allocatable :: param_tmp(:),sigma_tmp(:)
    real(kind=dp)              :: real_lattice_tmp(3,3),vec_tmp(3)
    real(kind=dp)              :: r_tmp


    call param_in_file

   !%%%%%%%%%%%%%%%%
   !System variables
   !%%%%%%%%%%%%%%%%

    param_unit      = 'ev'          !
    paramconfac     = 1.0_dp
    call param_get_keyword('parameter_unit',found,c_value=param_unit)
    if (param_unit.ne.'ev' .and. param_unit.ne.'electronvolt' .and. index(param_unit,'ryd').ne.1 ) &
         call io_error('Error: Value of parameter_unit not recognised in param_read')
    if (index(param_unit,'ryd').eq.1) paramconfac = 2.0_dp/hart

    coordinate   = 'frac'          !
    call param_get_keyword('coordinate',found,c_value=coordinate)
    if (index(coordinate,'frac').ne.1 .and. index(coordinate,'cart').ne.1) &
         call io_error('Error: Value of coordinate not recognised in param_read')

    length_unit     = 'ang'         !
    lenconfac       = 1.0_dp
    call param_get_keyword('length_unit',found,c_value=length_unit)
    if (index(length_unit,'ang').ne.1 .and. length_unit.ne.'bohr') &
         call io_error('Error: Value of length_unit not recognised in param_read')
    if (length_unit.eq.'bohr') lenconfac = 1.0_dp/bohr

    lboundary_open  = .False.
    boundary        = 'periodic'
    call param_get_keyword('boundary',found,c_value=boundary)
    if ((index(boundary,'peri').ne.1) .and. (index(boundary,'open').ne.1)) & 
       call io_error('Error: boundary not recognised in param_read')
    if (index(boundary,'open').eq.1) lboundary_open = .true.

    inputfile_setup1   = .false.            ! set to true to write .neigh &  inp2.mcin files and exit
    call param_get_keyword('inputfile_setup1',found,l_value=inputfile_setup1)
    ! We allow this keyword to be overriden by a command line arg -inp1
    if (input1_file_flag) inputfile_setup1 = .true.

    inputfile_setup2   = .false.            ! set to true to write .neigh &  .mcin files and exit
    call param_get_keyword('inputfile_setup2',found,l_value=inputfile_setup2)
    ! We allow this keyword to be overriden by a command line arg -inp2
    if (input2_file_flag) inputfile_setup2 = .true.

    hamfile_setup = .false.            ! set to true to write .ham file and exit
    call param_get_keyword('hamfile_setup',found,l_value=hamfile_setup)
    ! We allow this keyword to be overriden by a command line arg -ham
    if (ham_file_flag) hamfile_setup = .true.

    inputfile_setup    = .false.
    if (inputfile_setup1 .or. inputfile_setup2) inputfile_setup = .true.

    file_setup    = .false.
    if (hamfile_setup .or. inputfile_setup1 .or. inputfile_setup2) file_setup = .true.

    tems_mode       = 'lin'
    call param_get_keyword('tems_mode',found,c_value=tems_mode)
    if ((index(tems_mode,'man').ne.1) .and. (index(tems_mode,'lin').ne.1) &
         .and. (index(tems_mode,'inv').ne.1) .and. (index(tems_mode,'log').ne.1)) & 
       call io_error('Error: tems_mode not recognised')

    tem_start = 5.0_dp
    call param_get_keyword('tem_start',found,r_value=tem_start)

    tem_end   = 5.0_dp
    call param_get_keyword('tem_end',found,r_value=tem_end)

    tems_num   =   1
    call param_get_keyword('tems_num',found,i_value=tems_num)
    if (tems_num < 1) &
         call io_error('Error: tems_num must be greater than zero')

    allocate(tems(tems_num),stat=ierr)
    if (ierr/=0) call io_error('Error allocating tems in param_read')

    call param_get_keyword_vector('tems',found,tems_num,r_value=tems)

    if (index(tems_mode,'lin').eq.1 .or. index(tems_mode,'inv').eq.1 .or. index(tems_mode,'log').eq.1) then
       if (tem_start.le.zero) &
         call io_error('Error: tem_start must be greater than zero')
 
       if (tem_end.le.zero) &
         call io_error('Error: tem_end must be greater than zero')

       tems(1)=tem_start
       if (tems_num>1) then 
          if (index(tems_mode,'lin').eq.1)then
             tem_steps = real((tem_end-tem_start),dp)/real((tems_num-1),dp)
             do loop=1,tems_num
                tems(loop) = tem_start+(loop-1)*tem_steps
             enddo
          elseif (index(tems_mode,'inv').eq.1)then
             tem_inv_steps = real((1.0_dp/tem_end-1.0_dp/tem_start),dp)/real((tems_num-1),dp)
             do loop=1,tems_num
                tems(loop) = 1.0_dp/tem_start+(loop-1)*tem_inv_steps
                tems(loop) = 1.0_dp/tems(loop)
             enddo
          elseif (index(tems_mode,'log').eq.1) then
             do loop=1,tems_num
                tems(loop) = tem_start*exp(log(tem_end/tem_start)*(1.0_dp*(loop-1)/(tems_num-1)))
             enddo
          endif
       endif
    else!if(index(tems_mode,'man')>0) then

       if (.not.found) &
          call io_error('Error: tems_mode set to manual but no tems found')
       if (any(tems.le.zero))&
         call io_error('Error: tems must be greater than zero')

    endif

    steps_warmup  = 100000
    call param_get_keyword('steps_warmup',found,i_value=steps_warmup)
    if (steps_warmup < 0) steps_warmup  = 0

    steps_mc      = 200000
    call param_get_keyword('steps_mc',found,i_value=steps_mc)
    if (steps_mc < 1) call io_error('Error: steps_mc must be greater than zero')

    steps_measure = 2  
    call param_get_keyword('steps_measure',found,i_value=steps_measure)
    if (steps_measure < 0) call io_error('Error: steps_measure must be greater than zero')

    supercell_size    = 4
    call param_get_keyword_vector('supercell_size',found,3,i_value=supercell_size)
    if (any(supercell_size < 1)) then
           call io_error('Error: supercell_size must be greater than zero')
    end if
    num_supercell = supercell_size(1)*supercell_size(2)*supercell_size(3)

    !! Initial configuration
    initial_sconfig   = 'ferro'
    call param_get_keyword('initial_sconfig',found,c_value=initial_sconfig)
    !if (found) then
    !   if ((index(initial_sconfig,'ferro').ne.1) .and. (index(initial_sconfig,'rand').ne.1) &
    !        .and. (index(initial_sconfig,'file').ne.1) ) then
    !      call io_error('Error: Initial_sconfig not recognised')
    !   elseif (index(initial_sconfig,'file').eq.1) then
    !      inquire(file=trim(seedname)//'_sconfig.dat',exist=sconfig_found)
    !      if (.not. sconfig_found) &
    !           call io_error('Error: Initial_sconfig requested from file '& 
    !                          //'but '//trim(seedname)//'_sconfig.dat file not found')
    !   endif
    !endif
    if ((index(initial_sconfig,'ferro').ne.1) .and. (index(initial_sconfig,'rand').ne.1) &
            .and. (index(initial_sconfig,'file').ne.1) ) &
          call io_error('Error: Initial_sconfig not recognised')
    if (index(initial_sconfig,'file').eq.1) then
       inquire(file=trim(seedname)//'_sconfig.dat',exist=sconfig_found)
       if (.not. sconfig_found) &
          call io_error('Error: Initial_sconfig requested from file '& 
                        //'but '//trim(seedname)//'_sconfig.dat file not found')
    endif

    mcarlo_mode       = 'random'
    call param_get_keyword('mcarlo_mode',found,c_value=mcarlo_mode)
    if ((index(mcarlo_mode,'const').ne.1) .and. (index(mcarlo_mode,'rand').ne.1) ) &
       call io_error('Error: mcarlo_mode not recognised')

    tilt_angles_max   = 0.125
    call param_get_keyword('tilt_angles_max',found,r_value=tilt_angles_max)
    if ((index(mcarlo_mode,'const') .eq. 1) .and. tilt_angles_max < 0) &
        call io_error('Error: tilt_angles_max must be positive')

   ! lseed             = .false.
   ! allocate( seeds(num_nodes) ,stat=ierr)
   ! if (ierr/=0) call io_error('Error allocating seeds in param_read')
   ! seeds = 0
   ! call param_get_keyword_vector('seeds',found,num_nodes,i_value=seeds)
   ! if (found) lseed = .true.

    lseed             = .false.
    allocate( seeds(num_nodes) ,stat=ierr)
    if (ierr/=0) call io_error('Error allocating seeds in param_read')
    seeds = 0
    call param_get_range_vector('seeds',found,i_tmp,lcount=.true.)
    if (found .and. i_tmp.ne.num_nodes) &
       call io_error('Error: number of seeds should be equal to number of nodes')
    if (found) then
       call param_get_range_vector('seeds',found,num_nodes,.false.,seeds)
       lseed = .true.
    end if



    !%%%%%%%%%%%%%%%%
    !  Other Stuff
    !%%%%%%%%%%%%%%%%

    ! Lattice
    call param_get_block_length('unit_cell_cart',found,i_tmp,lunits)
    if (.not. found) call io_error('Error: Did not find the cell information in the input file')
    if (lunits) i_tmp = i_tmp-1
    if (found .and. (i_tmp .ne. 3)) call io_error('Error: Not correct the cell information in the input file')
    call param_get_keyword_block('unit_cell_cart',found,3,3,r_value=real_lattice_tmp)
    real_lattice = transpose(real_lattice_tmp)
    call utility_recip_lattice(real_lattice,recip_lattice,cell_volume)

    ! Atoms and magnetic moments
    lunits            = .false.
    call param_get_block_length('atoms_frac',found,i_tmp)
    call param_get_block_length('atoms_cart',found2,i_tmp,lunits)
    if (found.and.found2) call io_error('Error: Cannot specify both atoms_frac and atoms_cart')
    if (.not.found  .and. .not.found2)  call io_error('Error: No atoms_frac or atoms_cart block found')
    if (found .or. found2) then
       if (lunits) i_tmp = i_tmp-1         
       num_atoms = i_tmp
    end if
    if (num_atoms <1) call io_error('Error: No atom positions specified')
    allocate(convert_matrix_atoms(num_atoms,num_atoms),stat=ierr)
    if (ierr/=0) call io_error('Error allocating convert_matrix_atoms in param_read')
    ! read the postion of atoms
    call param_get_atoms_moments(lunits)
    num_total_atoms = num_supercell*num_atoms

    ! Neighbors
    neighbors_tol     = 0.001_dp
    call param_get_keyword('neighbors_tol',found,r_value=neighbors_tol)
    if (neighbors_tol < 0.0_dp) call io_error('Error: neighbors_tol must be positive')

    !!! Parameters for main program

    !Staggered_magnetization
    lstaggered = .false.
    call param_get_keyword('staggered_m',found,l_value=lstaggered)
    !if (found .and. (num_atoms.eq.1)) then
    if (lstaggered .and. (num_atoms.eq.1)) then
       lstaggered = .false.
       write(stdout,'(1x,a)')&
           'WARNING: Calculation of staggered magnetization for one atom is'
       write(stdout,'(9x,a)')&
           ' meaningless staggered_m was changed to .False.'
    endif
    allocate(staggered_coeff(num_atoms),stat=ierr)
    if (ierr/=0) call io_error('Error allocating staggered_coeff in param_read')
    staggered_coeff = 1
    call param_get_keyword_vector('staggered_m_coeff',found,num_atoms,i_value=staggered_coeff)
    !if (found .and. any(abs(staggered_coeff).ne.1)) &
    if (lstaggered .and.found .and. any(abs(staggered_coeff).ne.1)) &
       call io_error('Error: staggered_m_coeff must be +1 or -1')
    if (lstaggered .and.found .and. all(staggered_coeff.eq.1)) then
       write(stdout,'(1x,a)')&
            'WARNING: All coefficent of staggered magnetization are +1, this means staggered'
       write(stdout,'(1x,a)')&
            ' magnetization is equal to magnetization and it is a extra calculation'
    endif
    staggered_coeff = matmul(staggered_coeff,convert_matrix_atoms)

    ! Order parameter
    lorderparam = .false.
    lunits      = .false.
    call param_get_keyword('order_parameter',found,l_value=lorderparam)
    call param_get_block_length('order_parameter_axes_frac',found,i_tmp)
    call param_get_block_length('order_parameter_axes_cart',found2,i_tmp,lunits)
    if (found .and. found2) &
       call io_error('Error: Cannot specify both order_parameter_axes_frac and order_parameter_axes_cart')
    if (lorderparam .and. .not.file_setup .and. .not.found .and. .not.found2) &
       call io_error('Error: No order_parameter_axes_frac or order_parameter_axes_cart block found')

    if (found .or. found2) then
      if (lunits) i_tmp = i_tmp-1
      if (i_tmp .ne. num_atoms) &
         call io_error('Error: Wrong number of lines in order_parameter_axes_frac or order_parameter_axes_cart block ') 
      allocate(orderparam_axes_cart(3,num_atoms),stat=ierr)
      if (ierr/=0) call io_error('Error allocating order_parameter_axes_cart in param_read')
      allocate(orderparam_atoms_symbol(num_atoms),stat=ierr)
      if (ierr/=0) call io_error('Error allocating order_parameter_atoms_symbol in param_read')
      if (found) then
         call param_get_axes('order_parameter_axes_frac',found,num_atoms,&
                             c_value=orderparam_atoms_symbol,r_value=&
                             orderparam_axes_cart,lsort=.true.,lnorm=.true.)
      elseif (found2) then
         call param_get_axes('order_parameter_axes_cart',found,num_atoms,&
                              c_value=orderparam_atoms_symbol,r_value=&
                              orderparam_axes_cart,lsort=.true.,lnorm=.true.)
      end if
    end if

    ! Error from binning analysis
    lbinerror         = .false.
    call param_get_keyword('binning_error',found,l_value=lbinerror)

    num_binning_level=5
    call param_get_range_vector('binning_level',found,num_binning_level,lcount=.true.)
    if (num_binning_level<1) call io_error('Error: problem reading binning_level')
    allocate(binning_level(num_binning_level),stat=ierr)
    if (ierr/=0) call io_error('Error allocating binning_level in param_read')
    if (found) then
       call param_get_range_vector('binning_level',found,num_binning_level,.false.,binning_level)
    else
       do i_tmp=1,num_binning_level
          binning_level(i_tmp)=i_tmp
       enddo
    end if
    if (lbinerror .and. (any(binning_level<=0) .or. any(int(steps_mc/(2**binning_level)).le.1))) &
         call io_error('Error: binning_level contain non-valid numbers. '&
                        //'It must be greater than zero and as ' &
                       //'(steps_mc/2^(binning_level))-1>0.')

    !num_binning_level=0
    !call param_get_range_vector('binning_level',found,num_binning_level,lcount=.true.)
    !if (lbinerror .and. .not.found) &
    !   call io_error('Error: binning_error is true but no binning_level found')
    !if (found) then
    !   if (num_binning_level<1) call io_error('Error: problem reading binning_level')
    !   allocate(binning_level(num_binning_level),stat=ierr)
    !   if (ierr/=0) call io_error('Error allocating binning_level in param_read')
    !   call param_get_range_vector('binning_level',found,num_binning_level,.false.,binning_level)
    !   if (any(binning_level<=0) .or. any(int(steps_mc/(2**binning_level)).le.1)) &
    !        call io_error('Error: binning_level contain non-valid numbers. '&
    !                       //'It must be greater than zero and as ' &
    !                       //'(steps_mc/2^(binning_level))-1>0.')
    !end if

    ! Spin_correlations
    lspincorr         = .false.
    call param_get_keyword('spin_correlation',found,l_value=lspincorr)

    ! Write energy     
    energy_write      = .false.
    call param_get_keyword('energy_write',found,l_value=energy_write)
    energy_num_print =   1000          ! frequency to write
    call param_get_keyword('energy_num_print',found,i_value=energy_num_print)
    if (energy_num_print < 1) call io_error('Error: energy_num_print must be greater than zero')

    ! Structure Factor
    lsfactor          = .false.
    call param_get_keyword('sfactor',found,l_value=lsfactor)

    lsfactor_polar    = .false.
    vec_tmp(1) = 1.0_dp
    vec_tmp(2) = 0.0_dp
    vec_tmp(3) = 0.0_dp
    call param_get_keyword_vector('sfactor_polar',found,3,r_value=vec_tmp)
    if (found) lsfactor_polar = .true.
    ! convert fractional coordinate to cartesian
    sfactor_polar = matmul(vec_tmp,recip_lattice)
    
    vec_tmp = 0.0_dp
    call param_get_keyword_vector('sfactor_corner',found,3,r_value=vec_tmp)
    ! convert fractional coordinate to cartesian
    sfactor_corner = matmul(vec_tmp,recip_lattice)

    vec_tmp(1) = 1.0_dp
    vec_tmp(2) = 0.0_dp
    vec_tmp(3) = 0.0_dp
    call param_get_keyword_vector('sfactor_q1',found,3,r_value=vec_tmp)
    sfactor_q1 = matmul(vec_tmp,recip_lattice)

    vec_tmp(1) = 0.0_dp
    vec_tmp(2) = 0.0_dp
    vec_tmp(3) = 1.0_dp
    call param_get_keyword_vector('sfactor_q2',found,3,r_value=vec_tmp)
    sfactor_q2 = matmul(vec_tmp,recip_lattice)

    ! Check if sfactor_q1 and sfactor_q2 not to be parallel
    call utility_cross(sfactor_q1,sfactor_q2,vec_tmp)
    ! Area (modulus sfactor_q1 x sfactor_q2 = vec_tmp)
    r_tmp = sqrt(vec_tmp(1)**2 + vec_tmp(2)**2 + vec_tmp(3)**2)
    if (r_tmp < eps8) call io_error( &
      'Error: Vectors sfactor_q1 and sfactor_q2 ' &
      //'not linearly independent')

    sfactor_steps_measure = 100
    call param_get_keyword('sfactor_steps_measure',found,i_value=sfactor_steps_measure)
    !if (lsfactor .and. (sfactor_steps_measure<=0 .or. sfactor_steps_measure>steps_mc)) & 
    if (lsfactor .and. (sfactor_steps_measure<1 .or. sfactor_steps_measure>steps_mc)) & 
       call io_error('Error: sfactor_steps_measure must '&
                      //'be greater than zero and less than steps_mc.')

    sfactor_2dqmesh(1:2) = 50
    call param_get_keyword_vector('sfactor_2dqmesh',found,2,i_value=sfactor_2dqmesh)
    if (lsfactor .and. any(sfactor_2dqmesh<1)) then
           call io_error('Error: sfactor_2dqmesh must be greater than zero')
    end if
    sfactor_nqpts = (sfactor_2dqmesh(1)+1)*(sfactor_2dqmesh(2)+1)
 
    !Parallel Tempering
    pt                  = .false.
    call param_get_keyword('pt',found,l_value=pt)

    if (pt .and. tems_num<=1) &
       call io_error('Error: tems_num must be greater than one in parallel tempering algorithm')

    if (pt) call utility_sort(tems_num,r_array=tems)

    pt_steps_swap = 1
    call param_get_keyword('pt_steps_swap',found,i_value=pt_steps_swap)
    !if (pt_steps_swap.le.0) call io_error('Error: pt_steps_swap must be greater than zero')
    if (pt .and. pt_steps_swap<1) call io_error('Error: pt_steps_swap must be greater than zero')
    if (pt .and. energy_write .and. (pt_steps_swap .gt. energy_num_print)) then
       energy_num_print = pt_steps_swap
    endif

    pt_print_swap =   10000          ! frequency to write
    call param_get_keyword('pt_print_swap',found,i_value=pt_print_swap)
    if (pt .and. pt_print_swap<1) call io_error('Error: pt_print_swap must be greater than zero')

    ! Spin_glass
    spin_glass    = .false.
    call param_get_keyword('spin_glass',found,l_value=spin_glass)
    lspin_glass_seed = .false.
    call param_get_keyword('spin_glass_seed',found,i_value=spin_glass_seed)
    if (found) lspin_glass_seed = .true.


    ! shells(1,:)= # of jij shells,shells(2,:)= # of bij shells,shells(3,:)= # of dij shells
    allocate(shells(3,num_species),stat=ierr) 
    if (ierr/=0) call io_error('Error allocating shells in param_read')
    shells = 0

    ! this parameters is for -inp1 command
    shells(1,:) = 1    
    call param_get_keyword('shells_jij',found,i_value=i_tmp)
    if (found) shells(1,:) = i_tmp
    if (shells(1,1) < 1) call io_error('Error: shells_jij must be greater than zero')

    !Parameters_Jij
    lunits        = .false.
    call param_get_block_length('parameters_jij',found,i_tmp,lunits)
    if (.not.found  .and. inputfile_setup2) &
       call io_error('param_read: no parameters_jij block found')
    if (found) then
       if (lunits) i_tmp = i_tmp-1
       allocate(param_tmp(i_tmp),stat=ierr)
       if (ierr/=0) call io_error('Error first allocating param_tmp in param_read')
       allocate(type_tmp(2,i_tmp),stat=ierr)
       if (ierr/=0) call io_error('Error first allocating type_tmp in param_read')
       allocate(shells_tmp(i_tmp),stat=ierr)
       if (ierr/=0) call io_error('Error first allocating shells_tmp in param_read')
       allocate(sigma_tmp(i_tmp),stat=ierr)
       if (ierr/=0) call io_error('Error first allocating sigma_tmp in param_read')
       if (spin_glass) then
          call param_get_inp2jbd('parameters_jij',found,i_tmp,param_tmp,type_tmp,shells_tmp,sigma_tmp)
          if (any(shells_tmp .le. 0)) &
              call io_error('param_read: sh in parameters_jij must be greater than zero')
          allocate(parameters_sigma(num_species,num_species,maxval(shells_tmp(:))),stat=ierr)
          if (ierr/=0) call io_error('Error allocating parameters_sigma in param_read')
          parameters_sigma = 0.0_dp
       else
          call param_get_inp2jbd('parameters_jij',found,i_tmp,param_tmp,type_tmp,shells_tmp)
          if (any(shells_tmp .le. 0)) &
              call io_error('param_read: sh in parameters_jij must be greater than zero')
       endif
       allocate(parameters_jij(num_species,num_species,maxval(shells_tmp(:))),stat=ierr)
       if (ierr/=0) call io_error('Error allocating parameters_jij in param_read')
       parameters_jij = 0.0_dp;shells(1,:) = 0
       do loop=1,i_tmp
          parameters_jij(type_tmp(1,loop),type_tmp(2,loop),shells_tmp(loop)) = param_tmp(loop)
          if (spin_glass) parameters_sigma(type_tmp(1,loop),type_tmp(2,loop),shells_tmp(loop)) = sigma_tmp(loop)
          if (shells(1,type_tmp(1,loop)) .lt. shells_tmp(loop)) shells(1,type_tmp(1,loop)) = shells_tmp(loop)
       enddo
       deallocate(param_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error first deallocating param_tmp in param_read')
       deallocate(type_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error first deallocating type_tmp in param_read')
       deallocate(shells_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error first deallocating shells_tmp in param_read')
       deallocate(sigma_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error first deallocating sigma_tmp in param_read')

    endif

    ! Biguadratic
    have_biquad       = .false.
    call param_get_keyword('ham_bij',found,l_value=have_biquad)

    if (inputfile_setup1 .and. have_biquad) shells(2,:) = 1
    call param_get_keyword('shells_bij',found,i_value=i_tmp)
    if (found) shells(2,:) = i_tmp
    if (inputfile_setup1 .and. found .and. shells(2,1) < 1)    have_biquad = .false.

    !Parameters_Bij
    lunits        = .false.
    call param_get_block_length('parameters_bij',found,i_tmp,lunits)
    if (have_biquad .and. inputfile_setup2 .and. .not.found) &
       call io_error('param_read: ham_bij is true but no parameters_bij block found')
    if (found) then
       if (lunits) i_tmp = i_tmp-1
       allocate(param_tmp(i_tmp),stat=ierr)
       if (ierr/=0) call io_error('Error second allocating param_tmp in param_read')
       allocate(type_tmp(2,i_tmp),stat=ierr)
       if (ierr/=0) call io_error('Error second allocating type_tmp in param_read')
       allocate(shells_tmp(i_tmp),stat=ierr)
       if (ierr/=0) call io_error('Error second allocating shells_tmp in param_read')
       call param_get_inp2jbd('parameters_bij',found,i_tmp,param_tmp,type_tmp,shells_tmp)
       if (any(shells_tmp .le. 0)) &
              call io_error('param_read: sh in parameters_bij must be greater than zero')
       allocate(parameters_bij(num_species,num_species,maxval(shells_tmp(:))),stat=ierr)
       if (ierr/=0) call io_error('Error allocating parameters_bij in param_read')
       parameters_bij = 0.0_dp;shells(2,:) = 0
       do loop=1,i_tmp
          parameters_bij(type_tmp(1,loop),type_tmp(2,loop),shells_tmp(loop)) = param_tmp(loop)
          if (shells(2,type_tmp(1,loop)) .lt. shells_tmp(loop)) shells(2,type_tmp(1,loop))=shells_tmp(loop)
       enddo
       !if (inputfile_setup .and. maxval(shells(2,:)) <= 0)  have_biquad = .false.
       deallocate(param_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error second deallocating param_tmp in param_read')
       deallocate(type_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error second deallocating type_tmp in param_read')
       deallocate(shells_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error second deallocating shells_tmp in param_read')
    endif

    ! DM
    have_dm           = .false.
    call param_get_keyword('ham_dij',found,l_value=have_dm)

    if (inputfile_setup1 .and. have_dm) shells(3,:) = 1
    call param_get_keyword('shells_dij',found,i_value=i_tmp)
    if (found) shells(3,:) = i_tmp
    if (inputfile_setup1 .and. found .and. shells(3,1) < 1) have_dm = .false.

    !Parameters_Dij
    lunits        = .false.
    call param_get_block_length('parameters_dij',found,i_tmp,lunits)
    if (have_dm .and. inputfile_setup2 .and. .not.found) &
       call io_error('param_read: ham_dij is true but no parameters_dij block found')
    if (found) then
       if (lunits) i_tmp = i_tmp-1
       allocate(param_tmp(i_tmp),stat=ierr)
       if (ierr/=0) call io_error('Error third allocating param_tmp in param_read')
       allocate(type_tmp(2,i_tmp),stat=ierr)
       if (ierr/=0) call io_error('Error third allocating type_tmp in param_read')
       allocate(shells_tmp(i_tmp),stat=ierr)
       if (ierr/=0) call io_error('Error third allocating shells_tmp in param_read')
       call param_get_inp2jbd('parameters_dij',found,i_tmp,param_tmp,type_tmp,shells_tmp)
       if (any(shells_tmp .le. 0)) &
              call io_error('param_read: sh in parameters_dij must be greater than zero')
       allocate(parameters_dij(num_species,num_species,maxval(shells_tmp(:))),stat=ierr)
       if (ierr/=0) call io_error('Error allocating parameters_dij in param_read')
       parameters_dij = 0.0_dp;shells(3,:) = 0
       do loop=1,i_tmp
          parameters_dij(type_tmp(1,loop),type_tmp(2,loop),shells_tmp(loop)) = param_tmp(loop)
          if (shells(3,type_tmp(1,loop)) .lt. shells_tmp(loop)) shells(3,type_tmp(1,loop))=shells_tmp(loop)
       enddo
       deallocate(param_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error third deallocating param_tmp in param_read')
       deallocate(type_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error third deallocating type_tmp in param_read')
       deallocate(shells_tmp,stat=ierr)
       if (ierr/=0) call io_error('Error third deallocating shells_tmp in param_read')
    endif

    !Exchange Jij
    jij_num_nbors = 0
    lunits        = .false.  
    call param_get_block_length('jij_parameters',found,i_tmp,lunits)
    if (.not.found  .and. .not.file_setup) &
       call io_error('param_read: no jij_parameters block found')
    if (found) then
       if (lunits) i_tmp = i_tmp-1                
       jij_num_nbors = i_tmp

       allocate(jij_1atom_site(3,jij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating jij_1atom_site in param_read')
       allocate(jij_2atom_site(3,jij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating jij_2atom_site in param_read')
       allocate(jij(jij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating jij in param_read')
       allocate(jij_shell(jij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating jij_shell in param_read')
       jij_shell = 0
       allocate(spin_glass_sigma(jij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating spin_glass_sigma in param_read')
       spin_glass_sigma = 0.0_dp
       ! read jij in eV. jij_1atom_site in fractional coordinate
       if (lspincorr .and. spin_glass) then
          call param_get_parameters('jij_parameters',found,jij_num_nbors,jij_1atom_site,&
                                    jij_2atom_site,jij,shell=jij_shell,sigma=spin_glass_sigma)
          jij_shell_max = maxval(jij_shell)
       elseif (lspincorr) then
          call param_get_parameters('jij_parameters',found,jij_num_nbors,jij_1atom_site,&
                                    jij_2atom_site,jij,shell=jij_shell)
          jij_shell_max = maxval(jij_shell)
       elseif (spin_glass) then
          call param_get_parameters('jij_parameters',found,jij_num_nbors,jij_1atom_site,&
                                    jij_2atom_site,jij,sigma=spin_glass_sigma)

       else
          call param_get_parameters('jij_parameters',found,jij_num_nbors,jij_1atom_site,&
                                    jij_2atom_site,jij)
        
       endif
       if (lspincorr .and. any(jij_shell.le.0) )&
             call io_error('param_read: sh in jij_parameters block must be greater than zero')
       if (spin_glass .and. any(spin_glass_sigma.lt. zero) )&
             call io_error('param_read: sig in jij_parameters block must be positive')
    endif

    !Biquad Bij
    bij_num_nbors = 0
    lunits        = .false.  
    call param_get_block_length('bij_parameters',found,i_tmp,lunits)
    if (have_biquad .and. .not.file_setup .and. .not.found) &
       call io_error('param_read: ham_bij is true but no bij_parameters block found')
    if (found) then
       if (lunits) i_tmp = i_tmp-1
       bij_num_nbors = i_tmp

       allocate(bij_1atom_site(3,bij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating bij_1atom_site in param_read')
       allocate(bij_2atom_site(3,bij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating bij_2atom_site in param_read')
       allocate(bij(bij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating bij in param_read')
       ! the unit of bij is eV
       call param_get_parameters('bij_parameters',found,bij_num_nbors,bij_1atom_site,&
                                  bij_2atom_site,bij)
    endif

    !DM Dij
    dij_num_nbors = 0
    lunits        = .false.  
    call param_get_block_length('dij_parameters',found,i_tmp,lunits)
    if (have_dm .and. .not.file_setup .and. .not.found) &
       call io_error('param_read: ham_dij is true but no dij_parameters block found')
    if (found) then
       if (lunits) i_tmp = i_tmp-1
       dij_num_nbors = i_tmp

       allocate(dij_1atom_site(3,dij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating dij_1atom_site in param_read')
       allocate(dij_2atom_site(3,dij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating dij_2atom_site in param_read')
       allocate(dij(dij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating dij in param_read')
       allocate (dij_vectors_cart(3,dij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating dij_vectors_cart in param_read')
       allocate (dij_vectors_frac(3,dij_num_nbors),stat=ierr)
       if (ierr/=0) call io_error('Error allocating dij_vectors_frac in param_read')

       ! DM parameters 
       ! the unit of dij is eV
       call param_get_parameters('dij_parameters',found,dij_num_nbors,dij_1atom_site,&
                                 dij_2atom_site,dij)
    endif
    
    ! Dij vectors
    lunits        = .false.  
    call param_get_block_length('dij_vectors_frac',found,i_tmp)
    call param_get_block_length('dij_vectors_cart',found2,i_tmp,lunits)
    if (found .and. found2) &
       call io_error('Error: Cannot specify both dij_vectors_frac and dij_vectors_cart')
    if (have_dm .and. .not.file_setup .and. .not.found .and. .not.found2) &
       call io_error('Error: no dij_vectors_frac or dij_vectors_cart block found')             

    if (found .or. found2) then
       if(lunits) i_tmp = i_tmp-1
       if(i_tmp .ne. dij_num_nbors) &
           call io_error('Error: Wrong number of lines in dij_vectors_cart or dij_vectors_frac block')

       if (found) then
          call param_get_keyword_block('dij_vectors_frac',found,dij_num_nbors,3,r_value=dij_vectors_frac)
          if (have_dm .and. .not.file_setup) then
             do loop=1,dij_num_nbors
                call utility_frac_to_cart (dij_vectors_frac(:,loop),dij_vectors_cart(:,loop),real_lattice)
                ! Normalize Vectors
                if (all(abs(dij_vectors_cart(:,loop)) .lt. eps4)) &
                     call io_error('Error: dij_vectors_frac cannot be zero')
                dij_vectors_cart(:,loop)=dij_vectors_cart(:,loop)& 
                                        /sqrt(dot_product(dij_vectors_cart(:,loop),dij_vectors_cart(:,loop)))
             end do
          endif
       elseif (found2) then
          call param_get_keyword_block('dij_vectors_cart',found,dij_num_nbors,3,r_value=dij_vectors_cart)
          if (have_dm .and. .not.file_setup) then
             do loop=1,dij_num_nbors
                call utility_cart_to_frac (dij_vectors_cart(:,loop),dij_vectors_frac(:,loop),recip_lattice)
                ! Normalize Vectors
                if ( all(abs(dij_vectors_cart(:,loop)) .lt. eps4)) &
                     call io_error('Error: dij_vectors_cart cannot be zero')
                dij_vectors_cart(:,loop)=dij_vectors_cart(:,loop) & 
                                        /sqrt(dot_product(dij_vectors_cart(:,loop),dij_vectors_cart(:,loop)))
             end do
          endif
       endif
    endif

!    if(found .or. found2) then
!      if(lunits) i_tmp = i_tmp-1
!      if(i_tmp .ne. dij_num_nbors) &
!           call io_error('Error: Wrong number of lines in dij_vectors_cart or dij_vectors_frac block')
!    endif
!
!    call param_get_keyword_block('dij_vectors_frac',found,dij_num_nbors,3,r_value=dij_vectors_frac)
!    call param_get_keyword_block('dij_vectors_cart',found2,dij_num_nbors,3,r_value=dij_vectors_cart)
!    if (have_dm .and. .not.file_setup) then
!       if (found) then
!          do loop=1,dij_num_nbors
!             call utility_frac_to_cart (dij_vectors_frac(:,loop),dij_vectors_cart(:,loop),real_lattice)
!             ! Normalize Vectors
!             if (all(abs(dij_vectors_cart(:,loop)) .lt. eps4)) &
!                  call io_error('Error: dij_vectors_frac cannot be zero')
!             dij_vectors_cart(:,loop)=dij_vectors_cart(:,loop)& 
!                                     /sqrt(dot_product(dij_vectors_cart(:,loop),dij_vectors_cart(:,loop)))
!          end do
!       elseif (found2) then
!          do loop=1,dij_num_nbors
!             call utility_cart_to_frac (dij_vectors_cart(:,loop),dij_vectors_frac(:,loop),recip_lattice)
!             ! Normalize Vectors
!             if ( all(abs(dij_vectors_cart(:,loop)) .lt. eps4)) &
!                  call io_error('Error: dij_vectors_cart cannot be zero')
!             dij_vectors_cart(:,loop)=dij_vectors_cart(:,loop) & 
!                                     /sqrt(dot_product(dij_vectors_cart(:,loop),dij_vectors_cart(:,loop)))
!          end do
!       endif
!    endif

    ! Single-ion
    have_singleion    = .false.
    call param_get_keyword('ham_singleion',found,l_value=have_singleion)

    lunits = .false.
    call param_get_block_length('singleion_axes_frac',found,i_tmp,lunits)
    call param_get_block_length('singleion_axes_cart',found2,i_tmp,lunits)
    if (found .and. found2) &
       call io_error('Error: Cannot specify both singleion_axes_frac and singleion_axes_cart')
    if (have_singleion .and. .not.file_setup .and. .not.found .and. .not.found2) &
       call io_error('Error: ham_singeleion is true but no singleion_axes_frac or singleion_axes_cart block found')

    if (found .or. found2) then
       if (lunits) i_tmp = i_tmp-1
       if (i_tmp .ne. num_atoms)  &
            call io_error('Error: Wrong number of lines in singleion_axes_frac or singleion_axes_cart block ')
       allocate(singleion_parameters(num_atoms),stat=ierr)
       if (ierr/=0) call io_error('Error allocating singleion_parameters in param_read')
       allocate(singleion_axes_cart(3,num_atoms),stat=ierr)
       if (ierr/=0) call io_error('Error allocating singleion_axes_cart in param_read')
       allocate(singleion_atoms_symbol(num_atoms),stat=ierr)
       if (ierr/=0) call io_error('Error allocating singleion_atoms_symbol in param_read')
       if (found) then
          call param_get_axes('singleion_axes_frac',found,num_atoms,&
                              c_value=singleion_atoms_symbol,r_value=singleion_axes_cart,&
                              param_value=singleion_parameters,lsort=.true.,lnorm=.true.)
       elseif (found2) then
          call param_get_axes('singleion_axes_cart',found,num_atoms,&
                              c_value=singleion_atoms_symbol,r_value=singleion_axes_cart,&
                              param_value=singleion_parameters,lsort=.true.,lnorm=.true.)
       end if
       !! eV to Kelvin
       singleion_parameters = singleion_parameters/k_B  
    endif

    ! H-field    
    have_field    = .false.
    call param_get_keyword('ham_field',found,l_value=have_field)

    lunits = .false.
    call param_get_block_length('field_axes_frac',found,i_tmp)
    call param_get_block_length('field_axes_cart',found2,i_tmp,lunits)
    if (found.and.found2) &
       call io_error('Error: Cannot specify both field_axes_frac and field_axes_cart')
    if (have_field .and. .not.file_setup .and. .not.found .and. .not.found2) &
       call io_error('Error: singeleion is true but no field_axes_frac or field_axes_cart found')

    if (found .or. found2) then
       if (lunits) i_tmp = i_tmp-1
       if (i_tmp .ne. num_atoms)  &
            call io_error('Error: Wrong number of lines in field_axes_frac or field_axes_cart block ')
       allocate(field_parameters(num_atoms),stat=ierr)
       if (ierr/=0) call io_error('Error allocating field_parameters in param_read')
       allocate(field_axes_cart(3,num_atoms),stat=ierr)
       if (ierr/=0) call io_error('Error allocating field_axes_cart in param_read')
       allocate(field_atoms_symbol(num_atoms),stat=ierr)
       if (ierr/=0) call io_error('Error allocating field_atoms_symbol in param_read')
       
       if (found) then
          call param_get_axes('field_axes_frac',found,num_atoms,&
                              c_value=field_atoms_symbol,r_value=field_axes_cart,&
                              param_value=field_parameters,lsort=.true.,lnorm=.true.)
       elseif (found2) then
          call param_get_axes('field_axes_cart',found,num_atoms,&
                              c_value=field_atoms_symbol,r_value=field_axes_cart,&
                              param_value=field_parameters,lsort=.true.,lnorm=.true.)
       end if
       ! Tesla to Kelvin
       field_parameters(:) = field_parameters(:)*bohr_magn/k_B
    endif

   ! check to see that there are no unrecognised keywords

    if (any(len_trim(in_data(:))>0 )) then
       if (inputfile_setup1) then
       write(stdout,'(1x,a)') 'The following section of file '//trim(seedname)//'.inp1.mcin '&
                               //'contained unrecognised keywords'
       write(*,'(1x,a)') 'The following section of file '//trim(seedname)//'.inp1.mcin '&
                              //'contained unrecognised keywords'
       elseif (inputfile_setup2) then
       write(stdout,'(1x,a)') 'The following section of file '//trim(seedname)//'.inp2.mcin '&
                               //'contained unrecognised keywords'
       write(*,'(1x,a)') 'The following section of file '//trim(seedname)//'.inp2.mcin '&
                              //'contained unrecognised keywords'
       elseif (hamfile_setup) then
       write(stdout,'(1x,a)') 'The following section of file '//trim(seedname)//'.ham.mcin '&
                               //'contained unrecognised keywords'
       write(*,'(1x,a)') 'The following section of file '//trim(seedname)//'.ham.mcin '&
                              //'contained unrecognised keywords'
       else
       write(stdout,'(1x,a)') 'The following section of file '//trim(seedname)//'.mcin '&
                               //'contained unrecognised keywords'
       write(*,'(1x,a)') 'The following section of file '//trim(seedname)//'.mcin '&
                              //'contained unrecognised keywords'
       endif
       write(stdout,*)
       do loop=1,num_lines
          if (len_trim(in_data(loop))>0) then
             write(stdout,'(1x,a)') trim(in_data(loop))
             write(*,'(1x,a)') trim(in_data(loop))
          end if
       end do
       write(stdout,*)
       call io_error('Unrecognised keyword(s) in input file')
    end if

   ! For aesthetic purposes, convert some things to uppercase
    call param_uppercase()

    deallocate(in_data,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating in_data in param_read')

    deallocate(convert_matrix_atoms,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating convert_matrix_atoms in param_read')


    return
  end subroutine param_read


  !===================================================================!
  subroutine param_uppercase
    !===================================================================!
    !                                                                   !
    ! Convert a few things to uppercase to look nice in the output      !
    !                                                                   !
    !===================================================================!  

    implicit none

    integer :: nsp,ic

    ! Atom labels (eg, si --> Si)
    do nsp=1,num_species
       ic = ichar(atoms_label(nsp)(1:1))                           
       if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
            atoms_label(nsp)(1:1) = char(ic+ichar('Z')-ichar('z'))
       ic = ichar(atoms_symbol(nsp)(1:1))
       if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
            atoms_symbol(nsp)(1:1) = char(ic+ichar('Z')-ichar('z'))
    enddo

    return

  end subroutine param_uppercase


  !===================================================================!
  subroutine param_write
    !===================================================================!
    !                                                                  !
    ! write parameters to stdout                                       !
    !                                                                  !
    !===================================================================  

    use mc_utility,     only : utility_cart_to_frac,utility_frac_to_cart
    use mc_constants,   only : k_B,bohr_magn

    implicit none

    integer       :: i,loop,nsp,nat
    real(kind=dp) :: vec_tmp1(3),vec_tmp2(3)

    ! System
    write(stdout,*)
    write(stdout,'(40x,a6)') '------'
    write(stdout,'(40x,a6)') 'SYSTEM'
    write(stdout,'(40x,a6)') '------'
    write(stdout,*)
    if (lenconfac.eq.1.0_dp) then
       write(stdout,'(33x,a21)') 'Lattice Vectors (Ang)'
    else
       write(stdout,'(32x,a22)') 'Lattice Vectors (Bohr)'
    endif
    write(stdout,101) 'a_1',(real_lattice(1,I)*lenconfac, i=1,3)
    write(stdout,101) 'a_2',(real_lattice(2,I)*lenconfac, i=1,3)
    write(stdout,101) 'a_3',(real_lattice(3,I)*lenconfac, i=1,3)
    write(stdout,*)
    write(stdout,'(24x,a17,3x,f11.5)',advance='no') &
         'Unit Cell Volume:',cell_volume*lenconfac**3
    if (lenconfac.eq.1.0_dp) then
       write(stdout,'(2x,a7)') '(Ang^3)'
    else
       write(stdout,'(2x,a8)') '(Bohr^3)'
    endif
    write(stdout,*)
    if (lenconfac.eq.1.0_dp) then
       write(stdout,'(27x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    else
       write(stdout,'(27x,a34)') 'Reciprocal-Space Vectors (Bohr^-1)'
    endif
    write(stdout,101) 'b_1',(recip_lattice(1,I)/lenconfac, i=1,3)
    write(stdout,101) 'b_2',(recip_lattice(2,I)/lenconfac, i=1,3)
    write(stdout,101) 'b_3',(recip_lattice(3,I)/lenconfac, i=1,3)
    write(stdout,*)   ' '


   ! Atoms
   write(stdout,'(1x,a)') '*-----------------------------------------------------------------------------------*'
   if (lenconfac.eq.1.0_dp) then
      write(stdout,'(1x,a)') '| Atom  Type   Fractional Coordinate             Cartesian Coordinate(Ang)   Moment |'
   else
      write(stdout,'(1x,a)') '| Atom  Type   Fractional Coordinate             Cartesian Coordinate(Bohr)  Moment |'
   endif
   write(stdout,'(1x,a)') '+-----------------------------------------------------------------------------------+'
   do nsp=1,num_species
      do nat=1,atoms_species_num(nsp)
         write(stdout,'(1x,a1,1x,a3,2x,i3,3F10.5,3x,a1,3F10.5,2x,a1,F6.3,1x,a1)') '|',atoms_symbol(nsp),nsp,&
               atoms_pos_frac(:,nsp,nat),'|',atoms_pos_cart(:,nsp,nat)*lenconfac,'|',mag_moments(nsp,nat),'|'
      end do
   end do
   write(stdout,'(1x,a)') '*-----------------------------------------------------------------------------------*'
    
   if (.not.file_setup) then
     ! Super_cell
     write(stdout,'(36x,a)') '--------------'
     write(stdout,'(36x,a)') 'SUPERCELL SIZE'
     write(stdout,'(36x,a)') '--------------'
     write(stdout,'(29x,a,i3,1x,a1,i3,1x,a1,i3,6x)') 'Supercell size =',&
           supercell_size(1),'x',supercell_size(2),'x',supercell_size(3)

     write(stdout,*)
     write(stdout,'(1x,a85)') '*---------------------------------- MONTECARLO -------------------------------------*'
     write(stdout,'(1x,a48,10x,L9,17x,a1)')      '|  Using Parallel Tempering algorithm         :  ',pt,'|'
     write(stdout,'(1x,a48,10x,a9,17x,a1)')      '|  Initial spins configuration                : ',trim(initial_sconfig),'|'
     write(stdout,'(1x,a48,10x,a9,17x,a1)')      '|  Monte Carlo Mode                           : ',trim(mcarlo_mode),'|'
     if (index(mcarlo_mode,'const') > 0) then
        write(stdout,'(1x,a48,10x,f9.3,17x,a1)')  '|  Tilting angles maximum                     : ',tilt_angles_max,'|'
     endif
     write(stdout,'(1x,a48,10x,I9,17x,a1)')      '|  Number of Warm-up steps                    : ',steps_warmup,'|'
     write(stdout,'(1x,a48,10x,I9,17x,a1)')      '|  Number of Monte Carlo sampling steps       : ',steps_mc,'|'
     write(stdout,'(1x,a48,10x,I9,17x,a1)')      '|  Number of steps to take a measurement      : ',steps_measure,'|'
     if (pt) &
     write(stdout,'(1x,a48,10x,i9,17x,a1)')      '|  Number of MC steps to swap configurations  : ',pt_steps_swap,'|'
     write(stdout,'(1x,a48,10x,L9,17x,a1)')      '|  Order parameter (calculate/write)          :  ',lorderparam,'|'
     write(stdout,'(1x,a48,10x,L9,17x,a1)')      '|  Staggered Magnetization (calculate/write)  :  ',lstaggered,'|'
     write(stdout,'(1x,a48,10x,L9,17x,a1)')      '|  Structure Factor (calculate/write)         :  ',lsfactor,'|'
     write(stdout,'(1x,a48,10x,L9,17x,a1)')      '|  Energy               (write)               :  ',energy_write,'|'
     write(stdout,'(1x,a48,10x,L9,17x,a1)')      '|  Spin correlation (calculate/write)         :  ',lspincorr,'|'
     write(stdout,'(1x,a48,10x,L9,17x,a1)')      '|  Auto correlation (calculate/write)         :  ',lbinerror,'|'

     write(stdout,'(1x,a85)') '*--------------------------------- Temperatures ------------------------------------*'
     write(stdout,'(1x,a48,10x,a9,17x,a1)')    '|  Temperature mode                           :  ',trim(tems_mode),'|'
     write(stdout,'(1x,a48,10x,i9,17x,a1)')    '|  Number of temperatures                     : ',tems_num,'|'
     if (.not.pt) then
       write(stdout,'(1x,a48,10x,f9.3,17x,a1)')  '|  Starting temperature                       : ',tems(1),'|'
       write(stdout,'(1x,a48,10x,f9.3,17x,a1)')  '|  Ending temperature                         : ',tems(tems_num),'|'
     else
       write(stdout,'(1x,a48,10x,f9.3,17x,a1)')  '|  Lowest temperature                         : ',tems(1),'|'
       write(stdout,'(1x,a48,10x,f9.3,17x,a1)')  '|  Highest temperature                        : ',tems(tems_num),'|'
     endif
     write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
     write(stdout,'(1x,a85)') '*---------------------------------- Hamiltonian ------------------------------------*'
     write(stdout,'(1x,a85)') '|  Exchange                                   :                   T                 |'
     write(stdout,'(1x,a48,10x,L9,17x,a1)') '|  Using random exchange parameters           :  ',spin_glass,'|'
     write(stdout,'(1x,a48,10x,L9,17x,a1)') '|  Biquadratic                                :  ',have_biquad,'|'
     write(stdout,'(1x,a48,10x,L9,17x,a1)') '|  Dzyaloshinskii-Moriya                      :  ',have_dm,'|'
     write(stdout,'(1x,a48,10x,L9,17x,a1)') '|  Single-ion                                 :  ',have_singleion,'|'
     write(stdout,'(1x,a48,10x,L9,17x,a1)') '|  Magnetic field                             :  ',have_field,'|'
     write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
     write(stdout,*) ' '

     if (jij_num_nbors>0) then
       write(stdout,'(1x,a85)') '*------------------------------------ Exchange -------------------------------------*'
       if (coordinate .eq. 'frac') then
          write(stdout,'(1x,a85)') '|   In Fractional Coordinate:                                                       |'
       else
          if (lenconfac .eq. 1.0_dp) then
             write(stdout,'(1x,a85)') '|   In Cartesian Coordinate(Ang):                                                       |'
          else
             write(stdout,'(1x,a85)') '|   In Cartesian Coordinate(Bohr):                                                       |'
          endif
       end if
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       if (paramconfac.eq.1.0_dp) then
          write(stdout,'(1x,a85)') '|              ATOM 1                             ATOM 2                  J(eV)     |'
       else
          write(stdout,'(1x,a85)') '|              ATOM 1                             ATOM 2                  J(Ryd)    |'
       endif
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       do loop=1,jij_num_nbors              
          if (coordinate .eq. 'frac') then
             write(stdout,102) '|',jij_1atom_site(:,loop),'|',jij_2atom_site(:,loop),'|',jij(loop)*paramconfac,'|'
          else
             call utility_frac_to_cart (jij_1atom_site(:,loop),vec_tmp1,real_lattice)
             call utility_frac_to_cart (jij_2atom_site(:,loop),vec_tmp2,real_lattice)
             write(stdout,102) '|',vec_tmp1*lenconfac,'|',vec_tmp2,'|',jij(loop)*paramconfac,'|'
          endif
       end do
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       write(stdout,*) ' '
       
       if (spin_glass) then
          write(stdout,'(1x,a85)') '*-------------------------------------- Sigma --------------------------------------*'
          if (coordinate .eq. 'frac') then
             write(stdout,'(1x,a85)') '|   In Fractional Coordinate:                                                       |'
          else
             if (lenconfac .eq. 1.0_dp) then
                write(stdout,'(1x,a85)') '|   In Cartesian Coordinate(Ang):                                                       |'
             else
                write(stdout,'(1x,a85)') '|   In Cartesian Coordinate(Bohr):                                                       |'
             endif
          end if
          write(stdout,'(1x,a85)')    '*-----------------------------------------------------------------------------------*'
          if (paramconfac.eq.1.0_dp) then
             write(stdout,'(1x,a85)') '|              ATOM 1                             ATOM 2               Sigma(eV)    |'
          else
             write(stdout,'(1x,a85)') '|              ATOM 1                             ATOM 2               Sigma(Ryd)   |'
          endif
          write(stdout,'(1x,a85)')    '*-----------------------------------------------------------------------------------*'
          do loop=1,jij_num_nbors
             if (coordinate .eq. 'frac') then
                write(stdout,102) '|',jij_1atom_site(:,loop),'|',jij_2atom_site(:,loop),'|',&
                                   spin_glass_sigma(loop)*paramconfac,'|'
             else
                call utility_frac_to_cart (jij_1atom_site(:,loop),vec_tmp1,real_lattice)
                call utility_frac_to_cart (jij_2atom_site(:,loop),vec_tmp2,real_lattice)
                write(stdout,102) '|',vec_tmp1*lenconfac,'|',vec_tmp2,'|',&
                                     spin_glass_sigma(loop)*paramconfac,'|'
             endif
          end do
          write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
        elseif (lspincorr) then
          write(stdout,'(1x,a85)') '*------------------------------------- Shell ---------------------------------------*'
          if (coordinate .eq. 'frac') then
             write(stdout,'(1x,a85)') '|   In Fractional Coordinate:                                                       |'
          else
             if (lenconfac .eq. 1.0_dp) then
                write(stdout,'(1x,a85)') '|   In Cartesian Coordinate(Ang):                                                       |'
             else
                write(stdout,'(1x,a85)') '|   In Cartesian Coordinate(Bohr):                                                       |'
             endif
          end if
          write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
          write(stdout,'(1x,a85)') '|              ATOM 1                             ATOM 2                   Shell    |'
          write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
          do loop=1,jij_num_nbors
             if (coordinate .eq. 'frac') then
                write(stdout,103) '|',jij_1atom_site(:,loop),'|',jij_2atom_site(:,loop),'|',jij_shell(loop),'|'
             else
                call utility_frac_to_cart (jij_1atom_site(:,loop),vec_tmp1,real_lattice)
                call utility_frac_to_cart (jij_2atom_site(:,loop),vec_tmp2,real_lattice)
                write(stdout,103) '|',vec_tmp1*lenconfac,'|',vec_tmp2,'|',jij_shell(loop),'|'
             endif
          end do
          write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       endif
     end if
 
     if (have_biquad) then
       write(stdout,'(1x,a85)') '*------------------------------------ Biquadratic ----------------------------------*'
       if (coordinate .eq. 'frac') then
          write(stdout,'(1x,a85)') '|   In Fractional Coordinate:                                                       |'
       else
          if (lenconfac .eq. 1.0_dp) then
             write(stdout,'(1x,a85)') '|   In Cartesian Coordinate(Ang):                                                       |'
          else
             write(stdout,'(1x,a85)') '|   In Cartesian Coordinate(Bohr):                                                       |'
          endif
       end if
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       if (paramconfac.eq.1.0_dp) then
          write(stdout,'(1x,a85)') '|              ATOM 1                             ATOM 2                  B(eV)     |'
       else
          write(stdout,'(1x,a85)') '|              ATOM 1                             ATOM 2                  B(Ryd)    |'
       endif
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       do loop=1,bij_num_nbors              
          if (coordinate .eq. 'frac') then
             write(stdout,102) '|',bij_1atom_site(:,loop),'|',bij_2atom_site(:,loop),'|',bij(loop)*paramconfac,'|'
          else
             call utility_frac_to_cart (bij_1atom_site(:,loop),vec_tmp1,real_lattice)
             call utility_frac_to_cart (bij_2atom_site(:,loop),vec_tmp2,real_lattice)
             write(stdout,102) '|',vec_tmp1*lenconfac,'|',vec_tmp2,'|',bij(loop)*paramconfac,'|'
          endif
       end do
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       write(stdout,*) ' '
     end if
 
     if (have_dm) then
       write(stdout,'(1x,a85)') '*------------------------------- Dzyaloshinskii-Moriya -----------------------------*'
       if (coordinate .eq. 'frac') then
          write(stdout,'(1x,a85)') '|   In Fractional Coordinate:                                                       |'
       else
          if (lenconfac .eq. 1.0_dp) then
             write(stdout,'(1x,a85)') '|   In Cartesian Coordinate(Ang):                                                       |'
          else
             write(stdout,'(1x,a85)') '|   In Cartesian Coordinate(Bohr):                                                       |'
          endif
       end if
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       if (paramconfac.eq.1.0_dp) then
          write(stdout,'(1x,a85)') '|              ATOM 1                             ATOM 2                  D(eV)     |'
       else
          write(stdout,'(1x,a85)') '|              ATOM 1                             ATOM 2                  D(Ryd)    |'
       endif
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       do loop=1,dij_num_nbors              
          if (coordinate .eq. 'frac') then
             write(stdout,102) '|',dij_1atom_site(:,loop),'|',dij_2atom_site(:,loop),'|',dij(loop)*paramconfac,'|'
          else
             call utility_frac_to_cart (dij_1atom_site(:,loop),vec_tmp1,real_lattice)
             call utility_frac_to_cart (dij_2atom_site(:,loop),vec_tmp2,real_lattice)
             write(stdout,102) '|',vec_tmp1*lenconfac,'|',vec_tmp2,'|',dij(loop)*paramconfac,'|'
          endif
       end do
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       write(stdout,'(1x,a85)') '| Moriya-Vector   Fractional Coordinate           Cartesian Coordinate (Normalized) |'
       write(stdout,'(1x,a85)') '+-----------------------------------------------------------------------------------+'
       do loop=1,dij_num_nbors
          write(stdout,'(1x,a1,i6,4x,3F10.5,5x,a1,2x,3F10.5,5x,a1)') '|',loop,dij_vectors_frac(:,loop),'|',&
                dij_vectors_cart(:,loop),'|'
       end do
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       write(stdout,*) ' '
     end if

     if (have_singleion) then
       write(stdout,'(1x,a85)') '*---------------------------------- Single-ion-Axes --------------------------------*'
       if (paramconfac.eq.1.0_dp) then
          write(stdout,'(1x,a85)') '| Atoms   Fractional Coordinate     Cartesian Coordinate (Normalized)  \Delta(eV)   |'
       else
          write(stdout,'(1x,a85)') '| Atoms   Fractional Coordinate     Cartesian Coordinate (Normalized)  \Delta(Ryd)  |'
       endif
       write(stdout,'(1x,a85)') '+-----------------------------------------------------------------------------------+'
       do loop=1,num_atoms
          call utility_cart_to_frac (singleion_axes_cart(:,loop),vec_tmp1(:),recip_lattice)
          write(stdout,'(1x,a1,2x,a3,1x,3F8.4,5x,a1,2x,3F8.4,4x,a1,F15.8,1x,a1)') '|',singleion_atoms_symbol(loop)&
                                     ,vec_tmp1(:),'|',singleion_axes_cart(:,loop),'|',&
                                      singleion_parameters(loop)*paramconfac*k_B,'|'
       end do
       write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       write(stdout,*) ' '
     endif

     if (have_field) then
        write(stdout,'(1x,a85)') '*------------------------------- Magnetic Field ------------------------------------*'
        write(stdout,'(1x,a85)') '| Atoms   Fractional Coordinate     Cartesian Coordinate (Normalized)      H(T)     |'
        write(stdout,'(1x,a85)') '+-----------------------------------------------------------------------------------+'
        do loop=1,num_atoms
           call utility_cart_to_frac (field_axes_cart(:,loop),vec_tmp1(:),recip_lattice)
           write(stdout,'(1x,a1,2x,a3,1x,3F8.4,2x,a1,1x,3F8.4,7x,a1,1x,F15.8,1x,a1)') '|',field_atoms_symbol(loop)&
                                      ,vec_tmp1(:),'|',field_axes_cart(:,loop),'|',field_parameters(loop)*k_B/bohr_magn,'|'
        end do
        write(stdout,'(1x,a)') '*-----------------------------------------------------------------------------------*'
        write(stdout,*) ' '
     endif

     if (lorderparam) then
        write(stdout,'(1x,a85)') '*----------------------------- Order-Parameters-Axes -------------------------------*'
        write(stdout,'(1x,a85)') '| Atoms           Fractional Coordinate           Cartesian Coordinate (Normalized) |'
        write(stdout,'(1x,a85)') '+-----------------------------------------------------------------------------------+'
        do loop=1,num_atoms
           call utility_cart_to_frac (orderparam_axes_cart(:,loop),vec_tmp1(:),recip_lattice)
           write(stdout,'(1x,a1,2x,a3,5x,3F10.5,5x,a1,2x,3F10.5,5x,a1)') '|',orderparam_atoms_symbol(loop)&
                                      ,vec_tmp1(:),'|',orderparam_axes_cart(:,loop),'|'
        end do
        write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
        write(stdout,*) ' '
     endif

     if (lstaggered) then
        write(stdout,'(1x,a85)') '*----------------------------- Staggered Magnetization -----------------------------*'
        write(stdout,'(1x,a85)') '| Atom     Type          Fractional Coordinate           Staggered  Coefficient     |'
        write(stdout,'(1x,a85)') '+-----------------------------------------------------------------------------------+'
        i=0
        do nsp=1,num_species
           do nat=1,atoms_species_num(nsp)
              i = i+1
              if (staggered_coeff(i) > 0) then
                 write(stdout,'(1x,a1,1x,a3,6x,i3,7x,3F10.5,5x,a1,12x,a1,I1,13x,a1)') '|',atoms_symbol(nsp),nsp,&
                       atoms_pos_frac(:,nsp,nat),'|','+',staggered_coeff(i),'|'
              else
                 write(stdout,'(1x,a1,1x,a3,6x,i3,7x,3F10.5,5x,a1,12x,I2,13x,a1)') '|',atoms_symbol(nsp),nsp,&
                       atoms_pos_frac(:,nsp,nat),'|',staggered_coeff(i),'|'
              endif
           end do
        end do
        write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
       write(stdout,*) ' '
     endif
 
     if (lsfactor) then
        write(stdout,'(1x,a85)') '*-------------------------------- STRUCTURE-FACTOR ---------------------------------*'
        write(stdout,'(1x,a48,10x,I9,17x,a1)')    '|   Number of steps to take a SF measurement  : ',&
              sfactor_steps_measure,'|'
        write(stdout,'(1x,a85)') '|   2D q parameters (in reduced  coordinates):                                      |'
        call utility_cart_to_frac(sfactor_corner,vec_tmp1,real_lattice)
        write(stdout,'(1x,a14,2x,3F10.5,38x,a1)') '|     Corner: ', (vec_tmp1(i)/lenconfac,i=1,3),'|'
        call utility_cart_to_frac(sfactor_q1,vec_tmp1,real_lattice)
        write(stdout,'(1x,a14,2x,3F10.5,7x,a12,2x,i4,13x,a1)') &
             '|    Vector1: ', (vec_tmp1(i)/lenconfac,i=1,3),' Divisions:',sfactor_2dqmesh(1),'|'
        call utility_cart_to_frac(sfactor_q2,vec_tmp1,real_lattice)
        write(stdout,'(1x,a14,2x,3F10.5,7x,a12,2x,i4,13x,a1)') &
             '|    Vector2: ', (vec_tmp1(i)/lenconfac,i=1,3),' Divisions:',sfactor_2dqmesh(2),'|'
        if (lenconfac.eq.1.0_dp) then
           write(stdout,'(1x,a85)') '|   2D q parameters (in cartesian coordinates) (Ang^-1):                            |'
        else
           write(stdout,'(1x,a85)') '|   2D q parameters (in cartesian coordinates) (Bohr^-1):                           |'
        endif
        write(stdout,'(1x,a14,2x,3F10.5,38x,a1)') '|     Corner: ', (sfactor_corner(i),i=1,3),'|'
        write(stdout,'(1x,a14,2x,3F10.5,7x,a12,2x,i4,13x,a1)') &
             '|    Vector1: ', (sfactor_q1(i),i=1,3),' Divisions:',sfactor_2dqmesh(1),'|'
        write(stdout,'(1x,a14,2x,3F10.5,7x,a12,2x,i4,13x,a1)') &
             '|    Vector2: ', (sfactor_q2(i),i=1,3),' Divisions:',sfactor_2dqmesh(2),'|'

        write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
        write(stdout,*) ' '
     endif

   endif


101 format(24x,a3,2x,3F11.6)
102 format(1x,a1,3F11.6,1x,a1,3F11.6,1x,a,1F12.8,1x,a1)
103 format(1x,a1,3F11.6,1x,a1,3F11.6,1x,a,1i7,6x,a1)

  end subroutine param_write

  !==================================================================!
  subroutine param_write_seed
    !==================================================================!

    use mc_io, only : io_date
    use mc_comms, only : comms_reduce
    implicit none

    integer :: loop,i,j,k

     write(stdout,'(1x,a85)')      '*---------------------------- Random Number Generator ------------------------------*'
     write(stdout,'(1x,a85)')      '|  Random number generatorg algorithm         :          Mersenne Twister           |'
     if (.not.lseed) then
        write(stdout,'(1x,a85)')   '|  Initialized with seed derived from         :               Timer                 |'
     else
        write(stdout,'(1x,a85)')   '|  Initialized with seed derived from         :               *.mcin file           |'
     endif    
     write(stdout,'(1x,a85)')      '|  Seeds:                                                                           |'


     i = num_nodes/6
     j = mod(num_nodes,6)
     do loop=1,i
         write(stdout,'(1x,a1,5x,6I12,6x,a1)') '|',(seeds((loop-1)*6+k),k=1,6),'|'
     enddo
     if (j>0) then
        if ( j .eq. 1) then
           write(stdout,'(1x,a1,5x,1I12,66x,a1)') '|',(seeds(i*6+k),k=1,1),'|'
        elseif ( j .eq. 2) then
           write(stdout,'(1x,a1,5x,2I12,54x,a1)') '|',(seeds(i*6+k),k=1,2),'|'
        elseif ( j .eq. 3) then
           write(stdout,'(1x,a1,5x,3I12,42x,a1)') '|',(seeds(i*6+k),k=1,3),'|'
        elseif ( j .eq. 4) then
           write(stdout,'(1x,a1,5x,4I12,30x,a1)') '|',(seeds(i*6+k),k=1,4),'|'
        elseif ( j .eq. 5) then
           write(stdout,'(1x,a1,5x,5I12,18x,a1)') '|',(seeds(i*6+k),k=1,5),'|'
        endif
     endif

     if (spin_glass) then
        if (.not.lspin_glass_seed) then
           write(stdout,'(1x,a85)')   '|  Initialized spin glass seed                :               Timer                 |'
        else
           write(stdout,'(1x,a85)')   '|  Initialized spin glass seed                :  spin_glass_seed from *.mcin file   |'
        endif    
        write(stdout,'(1x,a48,10x,I12,14x,a1)') '|  Used seed for normal random generator      : ',spin_glass_seed,'|'
     endif    
     write(stdout,'(1x,a85)') '*-----------------------------------------------------------------------------------*'
     !write(stdout,*) ' '

  end subroutine param_write_seed

  !==================================================================!
  subroutine param_write_header
    !==================================================================!
    use mc_io, only : io_date
    implicit none


    character (len=9) :: cdate, ctime

    call io_date(cdate, ctime)

    write(stdout,*)
    write(stdout,*)  '              +-------------------------------------------------------+'
    write(stdout,*)  '              |                                                       |'
    write(stdout,*)  '              |                       ESpinS                          |'
    write(stdout,*)  '              |                                                       |'
    write(stdout,*)  '              +-------------------------------------------------------+'
    write(stdout,*)  '              |                                                       |'
    write(stdout,*)  '              |               Welcome to ESpinS Code                  |'
    write(stdout,*)  '              |        A Code for Thermodynamic properties of         |'
    write(stdout,*)  '              |        of Magnetic System based on Classical          |'
    write(stdout,*)  '              |               Monte Carlo simulation.                 |'
    write(stdout,*)  '              |                                                       |'
    write(stdout,*)  '              +-------------------------------------------------------+'
    write(stdout,*)  '              |      Execution started on ',cdate,' at ',ctime,'      |'
    write(stdout,*)  '              +-------------------------------------------------------+'

  end subroutine param_write_header

  !==================================================================!
  subroutine param_first_dealloc
    !==================================================================!
    !                                                                  !
    ! release memory from allocated parameters                         !
    !                                                                  !
    !===================================================================  
    use mc_io, only : io_error

    implicit none
    integer :: ierr

    if ( allocated ( atoms_label ) ) then
       deallocate (  atoms_label, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating atoms_label in param_first_dealloc')
    end if
    if ( allocated ( atoms_pos_cart ) ) then
       deallocate (  atoms_pos_cart, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating atoms_pos_cart in param_first_dealloc')
    end if

    if ( allocated( shells ) ) then
       deallocate( shells, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating shells in param_first_dealloc')
    end if

    if ( allocated( parameters_sigma ) ) then
       deallocate( parameters_sigma, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating parameters_sigma in param_first_dealloc')
    end if
    if ( allocated( parameters_jij ) ) then
       deallocate( parameters_jij, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating parameters_jij in param_first_dealloc')
    end if
    if ( allocated( jij_1atom_site ) ) then
       deallocate( jij_1atom_site, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating jij_1atom_site in param_first_dealloc')
    end if
    if ( allocated( jij_2atom_site ) ) then
       deallocate( jij_2atom_site, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating jij_2atom_site in param_first_dealloc')
    end if
    if ( allocated( jij ) ) then
       deallocate( jij, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating jij in param_first_dealloc')
    end if
    if ( allocated( jij_shell ) ) then
       deallocate( jij_shell, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating jij_shell in param_first_dealloc')
    end if
    if ( allocated( spin_glass_sigma ) ) then
       deallocate( spin_glass_sigma, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating spin_glass_sigma in param_first_dealloc')
    end if
    if ( allocated( mag_moments ) ) then
       deallocate( mag_moments, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating mag_moments in param_first_dealloc')
    end if

    if ( allocated( parameters_bij ) ) then
       deallocate( parameters_bij, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating parameters_bij in param_first_dealloc')
    end if
    if ( allocated( bij_1atom_site ) ) then
       deallocate( bij_1atom_site, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating bij_1atom_site in param_first_dealloc')
    end if
    if ( allocated( bij_2atom_site ) ) then
       deallocate( bij_2atom_site, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating bij_2atom_site in param_first_dealloc')
    end if
    if ( allocated( bij ) ) then
       deallocate( bij, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating bij in param_first_dealloc')
    end if

    if ( allocated( parameters_dij ) ) then
       deallocate( parameters_dij, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating parameters_dij in param_first_dealloc')
    end if

    if ( allocated( dij_1atom_site ) ) then
       deallocate( dij_1atom_site, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating dij_1atom_site in param_first_dealloc')
    end if
    if ( allocated( dij_2atom_site ) ) then
       deallocate( dij_2atom_site, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating dij_2atom_site in param_first_dealloc')
    end if
    if ( allocated( dij ) ) then
       deallocate( dij, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating dij in param_first_dealloc')
    end if
    if ( allocated( dij_vectors_frac ) ) then
       deallocate( dij_vectors_frac, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating dij_vectors_frac in param_first_dealloc')
    end if
    if ( allocated( dij_vectors_cart ) ) then
       deallocate( dij_vectors_cart, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating dij_vectors_cart in param_first_dealloc')
    end if

    if ( allocated( orderparam_atoms_symbol ) ) then
       deallocate( orderparam_atoms_symbol, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating order_parameter_atoms_symbol in param_first_dealloc')
    end if

    if ( allocated( singleion_atoms_symbol ) ) then
       deallocate( singleion_atoms_symbol, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating singleion_atoms_symbol in param_first_dealloc')
    end if

    if ( allocated( field_atoms_symbol ) ) then
       deallocate( field_atoms_symbol, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating field_atoms_symbol in param_first_dealloc')
    end if

    return

  end subroutine param_first_dealloc


  !==================================================================!
  subroutine param_second_dealloc
    !==================================================================!
    !                                                                  !
    ! release memory from allocated parameters                         !
    !                                                                  !
    !===================================================================  
    use mc_io, only : io_error

    implicit none
    integer :: ierr

    if ( allocated ( atoms_symbol ) ) then
       deallocate (  atoms_symbol, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating atoms_symbol in param_second_dealloc')
    end if
    if ( allocated ( atoms_pos_frac ) ) then
       deallocate (  atoms_pos_frac, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating atoms_pos_frac in param_second_dealloc')
    end if
    if ( allocated ( atoms_species_num ) ) then
       deallocate (atoms_species_num, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating atoms_species_num in param_second_dealloc')
    end if

    if ( allocated( orderparam_axes_cart ) ) then
       deallocate(orderparam_axes_cart, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating order_parameter_axes_cart in param_second_dealloc')
    end if

    if ( allocated( staggered_coeff ) ) then
       deallocate(staggered_coeff, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating staggered_coeff in param_second_dealloc')
    end if

    if ( allocated(binning_level) ) then
       deallocate(binning_level, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating binning_level in param_second_dealloc')
    end if

    if ( allocated(singleion_parameters) ) then
       deallocate(singleion_parameters, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating singleion_parameters in param_second_dealloc')
    end if
    
    if ( allocated(singleion_axes_cart) ) then
       deallocate( singleion_axes_cart, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating singleion_axes_cart in param_second_dealloc')
    end if

    if ( allocated(field_parameters) ) then
       deallocate( field_parameters, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating field_parameters in param_second_dealloc')
    end if

    if ( allocated(field_axes_cart) ) then
       deallocate( field_axes_cart, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating field_axes_cart in param_second_dealloc')
    end if

    if (allocated(seeds)) then
       deallocate(seeds, stat=ierr )
       if (ierr/=0) call io_error('Error in deallocating seeds in param_second_dealloc')
    endif

    if(allocated(tems))then
      deallocate(tems, stat=ierr )
      if (ierr/=0) call io_error('Error in deallocating tems in param_second_dealloc')
    endif

    return

  end subroutine param_second_dealloc

  !=======================================!
  subroutine param_in_file
    !=======================================!
    ! Load the *.mcin file into a character !
    ! array in_file, ignoring comments and  !
    ! blank lines and converting everything !
    ! to lowercase characters               !
    !=======================================!

    use mc_io,        only : io_file_unit,io_error,seedname,input1_file_flag,&
                             input2_file_flag,ham_file_flag
    use mc_utility,   only : utility_lowercase

    implicit none

    integer               :: in_unit,tot_num_lines,ierr,line_counter,loop,in1,in2
    character(len=maxlen) :: dummy

    if (input1_file_flag .or. ham_file_flag) then
       in_unit = io_file_unit( )
       open (in_unit, file=trim(seedname)//'.inp1.mcin',form='formatted',status='old',err=101)
    elseif (input2_file_flag) then
       in_unit = io_file_unit( )
       open (in_unit, file=trim(seedname)//'.inp2.mcin',form='formatted',status='old',err=102)
    else
       in_unit = io_file_unit( )
       open (in_unit, file=trim(seedname)//'.mcin',form='formatted',status='old',err=103)
    endif

    num_lines = 0; tot_num_lines = 0
    do
    if (input1_file_flag .or. ham_file_flag) then
       read(in_unit, '(a)', iostat = ierr, err= 201, end =210 ) dummy
    elseif (input2_file_flag) then
       read(in_unit, '(a)', iostat = ierr, err= 202, end =210 ) dummy
    else
       read(in_unit, '(a)', iostat = ierr, err= 203, end =210 ) dummy
    endif
       dummy = adjustl(dummy)
       tot_num_lines = tot_num_lines+1
       if (.not.dummy(1:1)=='!'  .and. .not.dummy(1:1)=='#') then
          if(len(trim(dummy)) > 0 ) num_lines = num_lines+1
       endif

    end do

101 call io_error('Error: Problem opening input file '//trim(seedname)//'.inp1.mcin')
102 call io_error('Error: Problem opening input file '//trim(seedname)//'.inp2.mcin')
103 call io_error('Error: Problem opening input file '//trim(seedname)//'.mcin')

201 call io_error('Error: Problem reading input file '//trim(seedname)//'.inp1.mcin')
202 call io_error('Error: Problem reading input file '//trim(seedname)//'.inp2.mcin')
203 call io_error('Error: Problem reading input file '//trim(seedname)//'.mcin')
210 continue
    rewind(in_unit)

    allocate(in_data(num_lines),stat=ierr)
    if (ierr/=0) call io_error('Error allocating in_data in param_in_file')

    line_counter=0
    do loop=1,tot_num_lines
       if (input1_file_flag .or. ham_file_flag) then
         read(in_unit, '(a)', iostat = ierr, err= 201 ) dummy
       elseif (input2_file_flag) then
         read(in_unit, '(a)', iostat = ierr, err= 202 ) dummy
       else
         read(in_unit, '(a)', iostat = ierr, err= 203 ) dummy
       endif
       dummy = utility_lowercase(dummy)
       dummy = adjustl(dummy)
       dummy = adjustl(dummy)
       if (dummy(1:1)=='!' .or.  dummy(1:1)=='#') cycle
       if (len(trim(dummy)) == 0 ) cycle
       line_counter = line_counter+1
       in1 = index(dummy,'!')
       in2 = index(dummy,'#')
       if (in1==0 .and. in2==0)  in_data(line_counter) = dummy
       if (in1==0 .and. in2>0)   in_data(line_counter) = dummy(:in2-1)
       if (in2==0 .and. in1>0)   in_data(line_counter) = dummy(:in1-1)
       if (in2>0 .and. in1>0)    in_data(line_counter) = dummy(:min(in1,in2)-1)
    end do

  end subroutine param_in_file

  !===========================================================================!
  subroutine param_get_keyword(keyword,found,c_value,l_value,i_value,r_value)
    !===========================================================================!
    !                                                                           !
    !             Finds the value of the required keyword.                      !
    !                                                                           !
    !===========================================================================!

    use mc_io,        only : io_error

    implicit none

    character(*),             intent(in)  :: keyword
    logical,                  intent(out) :: found
    character(*),optional,  intent(inout) :: c_value
    logical,optional,       intent(inout) :: l_value
    integer,optional,       intent(inout) :: i_value
    real(kind=dp),optional, intent(inout) :: r_value

    integer               :: kl,in,loop,itmp
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    do loop=1,num_lines
       in = index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1) cycle
       itmp = in+len(trim(keyword))
       if (in_data(loop)(itmp:itmp)/='=' &
            .and. in_data(loop)(itmp:itmp)/=':' &
            .and. in_data(loop)(itmp:itmp)/=' ') cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found = .true.
       dummy = in_data(loop)(kl+1:)
       in_data(loop)(1:maxlen) = ' '
       dummy = adjustl(dummy)
       if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy = dummy(2:)
          dummy = adjustl(dummy)
       end if
    end do

    if (found) then
       if (present(c_value)) c_value=dummy
       if (present(l_value)) then
          if (index(dummy,'t') > 0) then
             l_value = .true.
          elseif (index(dummy,'f') > 0) then
             l_value = .false.
          else
             call io_error('Error: Problem reading logical keyword '//trim(keyword))
          endif
       endif
       if (present(i_value)) read(dummy,*,err=220,end=220) i_value
       if (present(r_value)) read(dummy,*,err=220,end=220) r_value
    end if

    return

220 call io_error('Error: Problem reading keyword '//trim(keyword))


  end subroutine param_get_keyword


  !=========================================================================================!
  subroutine param_get_keyword_vector(keyword,found,length,c_value,l_value,i_value,r_value)
    !=========================================================================================!
    !                                                                                         !
    !                  Finds the values of the required keyword vector                        !
    !                                                                                         !
    !=========================================================================================!

    use mc_io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical,           intent(out) :: found
    integer,           intent(in)  :: length
    character(*),     optional, intent(inout) :: c_value(length)
    logical,          optional, intent(inout) :: l_value(length)
    integer,          optional, intent(inout) :: i_value(length)
    real(kind=dp),    optional, intent(inout) :: r_value(length)

    integer               :: kl,in,loop,i
    character(len=maxlen) :: dummy

    kl = len_trim(keyword)

    found = .false.

    do loop=1,num_lines
       in = index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1 ) cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found = .true.
       dummy = in_data(loop)(kl+1:)
       in_data(loop)(1:maxlen) = ' '
       dummy = adjustl(dummy)
       if (dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy = dummy(2:)
          dummy = adjustl(dummy)
       end if
    end do

    if (found) then
       if (present(c_value)) read(dummy,*,err=230,end=230) (c_value(i),i=1,length)
       if (present(l_value)) then
          ! I don't think we need this. Maybe read into a dummy charater
          ! array and convert each element to logical
       endif
       if (present(i_value)) read(dummy,*,err=230,end=230) (i_value(i),i=1,length)
       if (present(r_value)) read(dummy,*,err=230,end=230) (r_value(i),i=1,length)
    end if

    return

230 call io_error('Error: Problem reading keyword '//trim(keyword)//' in param_get_keyword_vector')


  end subroutine param_get_keyword_vector

  !========================================================!
  subroutine param_get_vector_length(keyword,found,length)
    !======================================================!
    !                                                      !
    !        Returns the length of a keyword vector        !
    !                                                      !
    !======================================================!

    use mc_io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical,           intent(out) :: found
    integer,           intent(out) :: length

    integer                        :: kl,in,loop,pos
    character(len=maxlen)          :: dummy

    kl = len_trim(keyword)

    found = .false.

    do loop=1,num_lines
       in = index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1) cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found = .true.
       dummy = in_data(loop)(kl+1:)
       dummy = adjustl(dummy)
       if (dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy = dummy(2:)
          dummy = adjustl(dummy)
       end if
    end do

    length = 0
    if (found) then
       if (len_trim(dummy)==0) call io_error('Error: keyword '//trim(keyword)//' is blank')
       length = 1
       dummy = adjustl(dummy)
       do
          pos = index(dummy,' ')
          dummy = dummy(pos+1:)
          dummy = adjustl(dummy)
          if(len_trim(dummy)>0) then
             length = length+1
          else
             exit
          endif

       end do

    end if

    return

  end subroutine param_get_vector_length


  !==============================================================================================!
  subroutine param_get_keyword_block(keyword,found,rows,columns,c_value,l_value,i_value,r_value)
    !==============================================================================================!
    !                                                                                              !
    !                           Finds the values of the required data block                        !
    !                           (The unit cell vectors and Dij vectors)                            !
    !                                                                                              !
    !==============================================================================================!

    use mc_constants, only : bohr
    use mc_io,        only : io_error
    use mc_utility,   only : utility_strip

    implicit none

    character(*),      intent(in)  :: keyword
    logical,           intent(out) :: found
    integer,           intent(in)  :: rows
    integer,           intent(in)  :: columns
    character(*),     optional, intent(inout) :: c_value(columns,rows)
    logical,          optional, intent(inout) :: l_value(columns,rows)
    integer,          optional, intent(inout) :: i_value(columns,rows)
    real(kind=dp),    optional, intent(inout) :: r_value(columns,rows)

    integer               :: in,ins,ine,loop,i,line_e,line_s,counter,blen
    logical               :: found_e,found_s,lconvert_length
    character(len=maxlen) :: dummy,end_st,start_st,ctemp

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)


    do loop=1,num_lines
       ins = index(in_data(loop),trim(keyword))
       if (ins==0) cycle
       in = index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s = loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s = .true.
    end do

    if (.not. found_s) then
       found = .false.
       return
    end if

    do loop=1,num_lines
       ine = index(in_data(loop),trim(keyword))
       if (ine==0) cycle
       in = index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e = loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e = .true.
    end do

    if (.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if (line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    ! number of lines of data in block
    blen = line_e-line_s-1

    if ((blen.ne.rows) .and. (blen.ne.rows+1)) &
         call io_error('Error: Wrong number of lines in block'//trim(keyword))          

    dummy = adjustl(in_data(line_s))
    ctemp = utility_strip(dummy)
    if (len(trim(ctemp)) .ne. (len(trim(keyword))+5)) &
       call io_error('Error: Unrecognised Begin block in '//trim(dummy)//' line')
    dummy = adjustl(in_data(line_e))
    ctemp = utility_strip(dummy)
    if(len(trim(ctemp)) .ne. (len(trim(keyword))+3)) & 
       call io_error('Error: Unrecognised End block in '//trim(dummy)//' line')
    
    found = .true.

    lconvert_length = .false.
    if (blen==rows+1) then
       dummy = in_data(line_s+1)
       dummy = adjustl(dummy)
       if (index(dummy,'ang').eq.1) then
          lconvert_length = .false.
          if (len(trim(dummy)) .ne. 3) call io_error('Error: Units in block '//trim(keyword)//' not recognised')
       elseif ( index(dummy,'bohr').eq.1) then
          lconvert_length = .true.
          if (len(trim(dummy)) .ne. 4) call io_error('Error: Units in block '//trim(keyword)//' not recognised')
       else
          call io_error('Error: Units in block '//trim(keyword)//' not recognised')
       endif
       in_data(line_s)(1:maxlen) = ' '
       line_s = line_s+1
    endif

!    r_value=1.0_dp
    counter = 0
    do loop=line_s+1,line_e-1
       dummy = in_data(loop)
       counter = counter+1
       if (present(c_value)) read(dummy,*,err=240,end=240) (c_value(i,counter),i=1,columns)
       if (present(l_value)) then
          ! I don't think we need this. Maybe read into a dummy charater
          ! array and convert each element to logical
       endif
       if (present(i_value)) read(dummy,*,err=240,end=240) (i_value(i,counter),i=1,columns)
       if (present(r_value)) read(dummy,*,err=240,end=240) (r_value(i,counter),i=1,columns)
    end do

    if (lconvert_length) then
       if (present(r_value)) then
          r_value = r_value*bohr
       endif
    endif

    in_data(line_s:line_e)(1:maxlen) = ' '


    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword))


  end subroutine param_get_keyword_block


  !=====================================================!
  subroutine param_get_block_length(keyword,found,rows,lunits)
    !=====================================================!
    !                                                     !
    !       Finds the length of the data block            !
    !                                                     !
    !=====================================================!

    use mc_io,        only : io_error
    use mc_utility,   only : utility_strip


    implicit none

    character(*),      intent(in)  :: keyword
    logical,           intent(out) :: found
    integer,           intent(out) :: rows
    logical, optional, intent(out) :: lunits

    integer               :: in,ins,ine,loop,line_e,line_s,pos
    logical               :: found_e,found_s
    character(len=maxlen) :: end_st,start_st,dummy,ctemp
    character(len=5)      :: c_punc=' ,;-:'
                         
    found_s = .false.
    found_e = .false.

    start_st= 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop=1,num_lines
       ins = index(in_data(loop),trim(keyword))
       if (ins==0) cycle
       in = index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s = loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s = .true.
    end do

    if (.not. found_s) then
       found = .false.
       return
    end if

    do loop=1,num_lines
       ine = index(in_data(loop),trim(keyword))
       if (ine==0) cycle
       in = index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e = loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e = .true.
    end do

    if (.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if (line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    rows = line_e-line_s-1

    dummy = adjustl(in_data(line_s))
    ctemp = utility_strip(dummy)
    if(len(trim(ctemp)) .ne. (len(trim(keyword))+5)) &
       call io_error('Error: Unrecognised Begin block in '//trim(dummy)//' line')
    dummy = adjustl(in_data(line_e))
    ctemp = utility_strip(dummy)
    if (len(trim(ctemp)) .ne. (len(trim(keyword))+3)) & 
       call io_error('Error: Unrecognised End block in '//trim(dummy)//' line')
    
    found = .true.

   if (present(lunits)) then
     dummy = in_data(line_s+1)
     dummy = adjustl(dummy)
     pos = scan(dummy,c_punc)
     if (pos > 1) ctemp=dummy(:pos-1)
     if ((index(ctemp,'ang').ne.0) .or. (index(ctemp,'bohr').ne.0) .or. & 
          (index(ctemp,'ev').ne.0) .or. (index(ctemp,'ryd').ne.0 )) go to 555
     lunits = .false.
   endif

    if(rows<=0) then !cope with empty blocks
       found = .false.
       in_data(line_s:line_e)(1:maxlen) = ' '
    end if

    return

555 lunits = .true.

    if (rows<=1) then !cope with empty blocks
       found = .false.
       in_data(line_s:line_e)(1:maxlen) = ' '
    end if

    return


  end subroutine param_get_block_length

    !====================================================================!
    subroutine param_get_range_vector(keyword,found,length,lcount,i_value)
    !====================================================================!
    !   Read a range vector eg. 1,2,3,4-10  or 1 3 400:100               !
    !   if(lcount) we return the number of states in length              !
    !====================================================================!
    use mc_io,        only : io_error

    implicit none

    character(*),      intent(in)    :: keyword
    logical          , intent(out)   :: found
    integer,           intent(inout) :: length
    logical,           intent(in)    :: lcount
    integer, optional, intent(out)   :: i_value(length)

    integer   :: kl, in,loop,num1,num2,i_punc
    integer   :: counter,i_digit,loop_r,range_size
    character(len=maxlen) :: dummy
    character(len=10), parameter :: c_digit="0123456789"
    character(len=2) , parameter :: c_range="-:"
    character(len=3) , parameter :: c_sep=" ,;"
    character(len=5) , parameter :: c_punc=" ,;-:"
    character(len=5)  :: c_num1,c_num2


    if(lcount .and. present(i_value) ) call io_error('param_get_range_vector: incorrect call')

    kl=len_trim(keyword)

    found=.false.

    do loop=1,num_lines
       in=index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1 ) cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found=.true.
       dummy=in_data(loop)(kl+1:)
       dummy=adjustl(dummy)
       if(.not. lcount) in_data(loop)(1:maxlen) = ' '
       if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy=dummy(2:)
          dummy=adjustl(dummy)
       end if
    end do

    if(.not. found) return

    counter=0
    if (len_trim(dummy)==0) call io_error('Error: keyword '//trim(keyword)//' is blank')
    dummy=adjustl(dummy)
    do
       i_punc=scan(dummy,c_punc)
       if(i_punc==0) call io_error('Error parsing keyword '//trim(keyword))
       c_num1=dummy(1:i_punc-1)
       read(c_num1,*,err=101,end=101) num1
       dummy=adjustl(dummy(i_punc:))
       !look for range
       if(scan(dummy,c_range)==1) then
          i_digit=scan(dummy,c_digit)
          dummy=adjustl(dummy(i_digit:))
          i_punc=scan(dummy,c_punc)
          c_num2=dummy(1:i_punc-1)
          read(c_num2,*,err=101,end=101) num2
          dummy=adjustl(dummy(i_punc:))
          range_size=abs(num2-num1)+1
          do loop_r=1,range_size
             counter=counter+1
             if(.not. lcount) i_value(counter)=min(num1,num2)+loop_r-1
          end do
       else
          counter=counter+1
          if(.not. lcount) i_value(counter)=num1
       end if

       if(scan(dummy,c_sep)==1) dummy=adjustl(dummy(2:))
       if(scan(dummy,c_range)==1) call io_error('Error parsing keyword '//trim(keyword)//' incorrect range')
       if(index(dummy,' ')==1) exit
    end do

    if(lcount) length=counter
    if(.not.lcount) then
       do loop=1,counter-1
          do loop_r=loop+1,counter
             if(i_value(loop)==i_value(loop_r)) &
                call io_error('Error parsing keyword '//trim(keyword)//' duplicate values')
          end do
        end do
    end if

    return

101 call io_error('Error parsing keyword '//trim(keyword))


   end  subroutine param_get_range_vector

  !===================================!
  subroutine param_get_atoms_moments(lunits)
    !===================================!
    !                                   !
    !   Fills the atom data block       !
    !                                   !
    !===================================!

    use mc_constants, only : bohr,zero,eps5
    use mc_utility,   only : utility_frac_to_cart,utility_cart_to_frac,&
                             utility_strip
    use mc_io,        only : io_error
    implicit none

    logical, intent(in)   :: lunits

    real(kind=dp)         :: atoms_pos_frac_tmp(3,num_atoms)
    real(kind=dp)         :: atoms_pos_cart_tmp(3,num_atoms)
    character(len=20)     :: keyword
    integer               :: in,ins,ine,loop,i,line_e,line_s,counter
    integer               :: counter1,blen
    integer               :: i_tmp,loop2,max_sites,ierr
    logical               :: found_e,found_s,found,frac
    character(len=maxlen) :: dummy,end_st,start_st,ctemp1
    character(len=maxlen) :: ctemp(num_atoms)
    character(len=maxlen) :: atoms_label_tmp(num_atoms)
    logical               :: lconvert_length
    !! for magnetic moments
    real(kind=dp)         :: mag_moments_tmp(num_atoms)

    keyword = "atoms_cart"
    frac = .false.
    call param_get_block_length("atoms_frac",found,i_tmp)
    if (found) then
       keyword = "atoms_frac"
       frac = .true.
    end if

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop=1,num_lines
       ins = index(in_data(loop),trim(keyword))
       if (ins==0) cycle
       in = index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s = loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s = .true.
    end do

    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0) cycle
       in = index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e = loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e = .true.
    end do

    if (.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    endif

    if (line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    dummy = adjustl(in_data(line_s))
    ctemp1 = utility_strip(dummy)
    if (len(trim(ctemp1)) .ne. (len(trim(keyword))+5)) &
       call io_error('Error: Unrecognised Begin block in '//trim(dummy)//' line')
    dummy = adjustl(in_data(line_e))
    ctemp1 = utility_strip(dummy)
    if (len(trim(ctemp1)) .ne. (len(trim(keyword))+3)) &
       call io_error('Error: Unrecognised End block in '//trim(dummy)//' line')

    blen = line_e-line_s-1

    lconvert_length = .false.
    if (blen==num_atoms+1) then  ! or    if (lunits) then
       dummy = in_data(line_s+1)
       dummy = adjustl(dummy)
       if (index(dummy,'ang').eq.1) then
          lconvert_length = .false.
          if (len(trim(dummy)) .ne. 3) call io_error('Error: Units in block '//trim(keyword)//' not recognised')
       elseif (index(dummy,'bohr').eq.1) then
          lconvert_length = .true.
          if (len(trim(dummy)) .ne. 4) call io_error('Error: Units in block '//trim(keyword)//' not recognised')
       else
          call io_error('Error: Units in block '//trim(keyword)//' not recognised')
       endif
       in_data(line_s)(1:maxlen) = ' '
       line_s = line_s+1
    endif

    counter = 0
    do loop=line_s+1,line_e-1
       dummy = in_data(loop)
       counter = counter+1
       if (frac) then
          read(dummy,*,err=340,end=340) atoms_label_tmp(counter),(atoms_pos_frac_tmp(i,counter),i=1,3),mag_moments_tmp(counter)
       else
          read(dummy,*,err=340,end=340) atoms_label_tmp(counter),(atoms_pos_cart_tmp(i,counter),i=1,3),mag_moments_tmp(counter)
       end if
    end do

    if (lconvert_length) atoms_pos_cart_tmp = atoms_pos_cart_tmp*bohr

    in_data(line_s:line_e)(1:maxlen) = ' '

    if (frac) then
      do loop=1,num_atoms
            !! Translate atoms to the home unit cell
          if (abs(atoms_pos_frac_tmp(1,loop)-nint(atoms_pos_frac_tmp(1,loop))) .lt. eps5)&
                  atoms_pos_frac_tmp(1,loop)=nint(atoms_pos_frac_tmp(1,loop))
          if (abs(atoms_pos_frac_tmp(2,loop)-nint(atoms_pos_frac_tmp(2,loop))) .lt. eps5)&
                  atoms_pos_frac_tmp(2,loop)=nint(atoms_pos_frac_tmp(2,loop))
          if (abs(atoms_pos_frac_tmp(3,loop)-nint(atoms_pos_frac_tmp(3,loop))) .lt. eps5)&
                  atoms_pos_frac_tmp(3,loop)=nint(atoms_pos_frac_tmp(3,loop))
          atoms_pos_frac_tmp(:,loop) = atoms_pos_frac_tmp(:,loop)-floor(atoms_pos_frac_tmp(:,loop)) 
          call utility_frac_to_cart (atoms_pos_frac_tmp(:,loop),atoms_pos_cart_tmp(:,loop),real_lattice)
      end do
    else
      do loop=1,num_atoms
          call utility_cart_to_frac (atoms_pos_cart_tmp(:,loop),atoms_pos_frac_tmp(:,loop),recip_lattice)
            !! Translate atoms to the home unit cell
          if (abs(atoms_pos_frac_tmp(1,loop)-nint(atoms_pos_frac_tmp(1,loop))) .lt. eps5)&
                  atoms_pos_frac_tmp(1,loop)=nint(atoms_pos_frac_tmp(1,loop))
          if (abs(atoms_pos_frac_tmp(2,loop)-nint(atoms_pos_frac_tmp(2,loop))) .lt. eps5)&
                  atoms_pos_frac_tmp(2,loop)=nint(atoms_pos_frac_tmp(2,loop))
          if (abs(atoms_pos_frac_tmp(3,loop)-nint(atoms_pos_frac_tmp(3,loop))) .lt. eps5)&
                  atoms_pos_frac_tmp(3,loop)=nint(atoms_pos_frac_tmp(3,loop))
          atoms_pos_frac_tmp(:,loop) = atoms_pos_frac_tmp(:,loop)-floor(atoms_pos_frac_tmp(:,loop)) 
          call utility_frac_to_cart (atoms_pos_frac_tmp(:,loop),atoms_pos_cart_tmp(:,loop),real_lattice)
      end do
    end if

    if (any(mag_moments_tmp .eq. zero)) &
        call io_error('Error: mag_moments can not be zero')

    ! Now we sort the data into the proper structures
    num_species = 1
    ctemp(1) = atoms_label_tmp(1)
    do loop=2,num_atoms
       do loop2=1,loop-1
          if (trim(atoms_label_tmp(loop))==trim(atoms_label_tmp(loop2))) exit
          if (loop2==loop-1) then 
             num_species = num_species+1
             ctemp(num_species) = atoms_label_tmp(loop)
          end if
       end do
    end do

    allocate(atoms_species_num(num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_species_num in param_get_atoms_moments')
    allocate(atoms_label(num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_label in param_get_atoms_moments')
    allocate(atoms_symbol(num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_symbol in param_get_atoms_moments')

    atoms_species_num(:) = 0

    do loop=1,num_species
       atoms_label(loop)=ctemp(loop)
       do loop2=1,num_atoms
          if (trim(atoms_label(loop))==trim(atoms_label_tmp(loop2))) then
             atoms_species_num(loop) = atoms_species_num(loop)+1
          end if
       end do
    end do

    max_sites = maxval(atoms_species_num)
    allocate(atoms_pos_frac(3,num_species,max_sites),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_pos_frac in param_get_atoms_moments')
    allocate(atoms_pos_cart(3,num_species,max_sites),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_pos_cart in param_get_atoms_moments')
    allocate(mag_moments(num_species,max_sites),stat=ierr)
       if (ierr/=0) call io_error('Error allocating mag_moments in param_get_atoms_moments')

    counter1 = 0
    convert_matrix_atoms = 0
    do loop=1,num_species
       counter = 0
       do loop2=1,num_atoms
          if (trim(atoms_label(loop))==trim(atoms_label_tmp(loop2))) then
             counter  = counter+1
             counter1 = counter1+1
             convert_matrix_atoms(loop2,counter1) = 1
             atoms_pos_frac(:,loop,counter) = atoms_pos_frac_tmp(:,loop2)
             atoms_pos_cart(:,loop,counter) = atoms_pos_cart_tmp(:,loop2)
             mag_moments(loop,counter)      = mag_moments_tmp(loop2)
          end if
       end do
    end do

    do loop=1,num_species    
       atoms_symbol(loop)(1:3) = atoms_label(loop)(1:3)
    end do

    return

340 call io_error('Error: Problem reading block keyword '//trim(keyword))

  end subroutine param_get_atoms_moments

  !=====================================================================================!
  subroutine param_get_parameters(keyword,found,rows,fatom_site,satom_site,parameters,shell,sigma)
    !=====================================================================================!
    !                                                                                     !
    !   Read the (jij,bij,dij) block and return the fractional                            !
    !   position of the atom site and its neighbors and the coupling parameters           !
    !                                                                                     !                                    
    !=====================================================================================!

    use mc_constants, only : bohr,hart
    use mc_utility,   only : utility_cart_to_frac,utility_frac_to_cart,&
                             utility_strip
    use mc_io,        only : io_error

    implicit none
    !!
    character(*)            , intent(in)     :: keyword
    integer                 , intent(in)     :: rows
    logical                 , intent(out)    :: found             
    real(kind=dp)           , intent(inout)  :: fatom_site(3,rows)
    real(kind=dp)           , intent(inout)  :: satom_site(3,rows)
    real(kind=dp)           , intent(inout)  :: parameters(rows)
    integer       , optional, intent(inout)  :: shell(rows)
    real(kind=dp) , optional, intent(inout)  :: sigma(rows)

    real(kind=dp)                    :: pos_frac1(3),pos_frac2(3)
    real(kind=dp)                    :: pos_cart1(3),pos_cart2(3)
    integer                          :: in,ins,ine,loop,line_e,line_s,counter
    integer                          :: line,pos1,pos2,blen,i,counter_sh,counter_sig
    logical                          :: found_e,found_s,found_length,found_energy
    character(len=maxlen)            :: dummy,end_st,start_st,ctemp
    !
    real(kind=dp)                    :: parameters_tmp                                
    real(kind=dp)                    :: shell_tmp,sigma_tmp                                
    logical                          :: lconvert_length,lconvert_energy
    character(len=5) , parameter     :: c_punc=" ,;-:"


    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop=1,num_lines
       ins = index(in_data(loop),trim(keyword))
       if (ins==0) cycle
       in = index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s = loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s = .true.
    end do

    if (.not. found_s) then
       found = .false.
       return
    end if

    do loop=1,num_lines
       ine = index(in_data(loop),trim(keyword))
       if (ine==0) cycle
       in = index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e = loop
       if (found_e) then
          call io_error('param_get_parameters: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e = .true.
    end do

    if (.not. found_e) then
       call io_error('param_get_parameters: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if (line_e<=line_s) then
       call io_error('param_get_parameters: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    dummy = adjustl(in_data(line_s))
    ctemp = utility_strip(dummy)
    if (len(trim(ctemp)) .ne. (len(trim(keyword))+5)) &
       call io_error('Error: Unrecognised Begin block in '//trim(dummy)//' line')
    dummy = adjustl(in_data(line_e))
    ctemp = utility_strip(dummy)
    if (len(trim(ctemp)) .ne. (len(trim(keyword))+3)) & 
       call io_error('Error: Unrecognised End block in '//trim(dummy)//' line')
    
    blen = line_e-line_s-1
    found = .true.

    lconvert_length = .false.
    lconvert_energy = .false.
    if (blen == rows+1) then
      dummy = adjustl(in_data(line_s+1))
      if (len_trim(dummy)==0) call io_error('Error: keyword '//trim(keyword)//' is blank')
      found_length = .false.
      found_energy = .false.
      do i=1,2
         pos1 = scan(dummy,c_punc)
         if(pos1==0 .or. pos1==1) call io_error('Error parsing keyword '//trim(keyword)) 
         ctemp = dummy(1:pos1-1)
         select case(ctemp)
         case ('ang')
              lconvert_length = .false.
              if (found_length) &
                 call io_error('Error: Found the length unit in block '//trim(keyword)//' more than once')
              found_length = .true.
         case ('bohr')
              lconvert_length = .true.
              if (found_length)  &
                 call io_error('Error: Found the length unit in block '//trim(keyword)//' more than once')
              found_length = .true.
         case ('ev')
              lconvert_energy = .false.
              if (found_energy) &
                 call io_error('Error: Found the energy unit in block '//trim(keyword)//' more than once')
              found_energy = .true.
         case ('ryd')
              lconvert_energy = .true.
              if (found_energy) &
                 call io_error('Error: Found the energy unit in block '//trim(keyword)//' more than once')
              found_energy = .true.
         case default
              call io_error('Error: Units in block '//trim(keyword)//' not recognised')
         end select
         dummy = adjustl(dummy(pos1:))
         if(scan(dummy,c_punc) == 1) dummy = adjustl(dummy(2:))
         if(index(dummy,' ')==1) exit
      enddo
       in_data(line_s)(1:maxlen) = ' '
       line_s = line_s + 1
    endif
 
    counter=0
    do line=line_s+1,line_e-1
       dummy = utility_strip(in_data(line))
       dummy = adjustl(dummy)
       pos1 = index(dummy,':')
       if (pos1==0) call io_error('param_read_parameters: malformed positions &
                    &definition: '//trim(dummy))
       ctemp = dummy(:pos1-1)
       ! Read the first atomic site
       if (index(ctemp,'c1=')>0) then
          ctemp = ctemp(4:)
          read(ctemp,*,err=105,end=105) (pos_cart1(i),i=1,3)
          if (lconvert_length) pos_cart1 = pos_cart1 * bohr
          call utility_cart_to_frac (pos_cart1(:),pos_frac1(:),recip_lattice)          
       elseif (index(ctemp,'f1=')>0) then
          ctemp = ctemp(4:)
          read(ctemp,*,err=105,end=105) (pos_frac1(i),i=1,3)
       else
          call io_error('param_get_parameters: Atom site not recognised '//trim(ctemp))
       end if

       !Now we know the sites for the first atom. Get the its neighbor
       
       dummy = dummy(pos1+1:)
       pos1 = index(dummy,':')
       if (pos1==0) call io_error('param_read_parameters: malformed postions &
                    &definition: '//trim(dummy))
       ctemp = dummy(:pos1-1)
       if (index(ctemp,'c2=')>0) then
          ctemp = ctemp(4:)
          read(ctemp,*,err=105,end=105) (pos_cart2(i),i=1,3)
          if (lconvert_length) pos_cart2 = pos_cart2 * bohr
          call utility_cart_to_frac (pos_cart2(:),pos_frac2(:),recip_lattice)
       elseif (index(ctemp,'f2=')>0) then
          ctemp = ctemp(4:)
          read(ctemp,*,err=105,end=105) (pos_frac2(i),i=1,3)
       else
          call io_error('param_get_parameters: Exchanges site not recognised '//trim(ctemp))
       end if

       dummy = dummy(pos1+1:)
       pos1=index(dummy,':')
       if (pos1 == 0) then
          ctemp = dummy
       else
          ctemp = dummy(:pos1 - 1)
       endif
       if (keyword .eq. 'jij_parameters') then
          pos2 = index(ctemp,'jij=')
       elseif (keyword .eq. 'bij_parameters') then
          pos2 = index(ctemp,'bij=')
       elseif (keyword .eq. 'dij_parameters') then
          pos2 = index(ctemp,'dij=')
       endif
       if(pos2==0) call io_error('param_read_parameters: malformed parameters' &
                                 //' definition: '//trim(dummy))
       read(ctemp(5:),*,err=105,end=105) parameters_tmp

       counter_sh=0
       counter_sig=0
       if (pos1 > 0) then
          dummy = dummy(pos1+1:)       
          do
            pos2 = index(dummy,':')
            if (pos2 == 0) then
               ctemp = dummy
            else
               ctemp = dummy(:pos2 - 1)
            endif
            if (index(ctemp, 'sh=') == 1) then
               read (ctemp(4:), *, err=105, end=105) shell_tmp
               counter_sh = counter_sh+1
            elseif (index(ctemp, 'sig=') == 1) then
               read (ctemp(5:), *, err=105, end=105) sigma_tmp
               counter_sig = counter_sig+1
            !else
            !   if (present(shell)) then
            !      call io_error('param_get_parameters: Problem reading shell '//trim(ctemp))
            !   elseif (present(sigma)) then
            !      call io_error('param_get_parameters: Problem reading sig '//trim(ctemp))
            !   else
            !      call io_error('param_get_parameters: Problem reading following part of ' &
            !                    //' line '//trim(ctemp)) 
            !   endif
            endif
            if (pos2 == 0) exit
            dummy = dummy(pos2+1:)
          enddo
       endif
       if (present(shell) .and. counter_sh .eq. 0) &
            call io_error('param_get_parameters: lspincorr is true but no shell is '&
                          //'found in '//trim(keyword)// ' block')
       if (present(sigma) .and. counter_sig .eq. 0) &
            call io_error('param_get_parameters: spin_glass is true but no sigma is '&
                          //'found in '//trim(keyword)// ' block')
       if (present(shell).and. counter_sh .ne. 1) &
          call io_error('param_get_parameters: Problem reading shell in '//trim(keyword)//' block')
       if (present(sigma).and. counter_sig .ne. 1) &
          call io_error('param_get_parameters: Problem reading sigma in '//trim(keyword)//' block')

!       if (present(shell)) then
!         ! Get the Jij or Bij or Dij parameters
!         pos1 = index(dummy,':')
!         if(pos1==0) call io_error('param_read_parameters: malformed parameters &
!                                    & definition: '//trim(dummy))
! 
!         if (keyword .eq. 'jij_parameters') pos2=index(dummy,'j=')
!         if (keyword .eq. 'bij_parameters') pos2=index(dummy,'b=')
!         if (keyword .eq. 'dij_parameters') pos2=index(dummy,'d=')
!         if(pos2==0) call io_error('param_read_parameters: malformed parameters &
!                                  & definition: '//trim(dummy))
!         ctemp = dummy(pos2+2:pos1-1)
!         read(ctemp,*,err=105,end=105) parameters_tmp      
! 
!        ! Get the shell number
!         dummy = dummy(pos1+1:)
!         pos1  = index(dummy,'sh=')
!         if (pos1==0) call io_error('param_read_parameters: malformed shells &
!                                    & definition: '//trim(dummy))
!         ctemp = dummy(pos1+3:)
!         read(ctemp,*,err=105,end=105) neigh                   
!       else
!         ! Get the Jij or Bij or Dij parameters
!
!         if (keyword .eq. 'jij_parameters') pos1=index(dummy,'j=')
!         if (keyword .eq. 'bij_parameters') pos1=index(dummy,'b=')
!         if (keyword .eq. 'dij_parameters') pos1=index(dummy,'d=')
!         if(pos1==0)  call io_error('param_read_parameters: malformed parameters &
!                                    & definition: '//trim(dummy))
!         ctemp = dummy(pos1+2:)
!         read(ctemp,*,err=105,end=105) parameters_tmp      
!       endif

       counter = counter+1
       fatom_site(:,counter) = pos_frac1
       satom_site(:,counter) = pos_frac2
       parameters(counter)   = parameters_tmp
       if (present(shell)) &
          shell(counter)     = shell_tmp       
       if (present(sigma)) &
          sigma(counter)     = sigma_tmp       
    end do !end loop over parameters block

    if (lconvert_energy) parameters = parameters*hart*0.5_dp
    if (lconvert_energy .and. present(sigma)) &
       sigma = sigma*hart*0.5_dp

    in_data(line_s:line_e)(1:maxlen) = ' '

    return

105 call io_error('param_get_parameters: Problem reading ' //trim(keyword) //trim(ctemp))         

  end subroutine param_get_parameters


  !=====================================================================================!
  subroutine param_get_inp2jbd(keyword,found,rows,param,species,shell,sigma)
    !=====================================================================================!
    !                                                                                     !
    !                                                                                     !                                    
    !=====================================================================================!

    use mc_constants, only : bohr,hart
    use mc_utility,   only : utility_strip                              
    use mc_io,        only : io_error

    implicit none
    !!
    character(*)           , intent(in)    :: keyword
    integer                , intent(in)    :: rows
    logical                , intent(out)   :: found             
    real(kind=dp)          , intent(inout) :: param(rows)
    real(kind=dp), optional, intent(inout) :: sigma(rows)
    integer                , intent(inout) :: species(2,rows)
    integer                , intent(inout) :: shell(rows)

    integer                          :: in,ins,ine,loop,line_e,line_s
    integer                          :: line,pos1,pos2,blen,i,counter
    logical                          :: found_e,found_s,found_length,found_energy
    character(len=maxlen)            :: dummy,end_st,start_st,ctemp
    logical                          :: lconvert_length,lconvert_energy
    character(len=5) , parameter     :: c_punc=" ,;-:"

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop=1,num_lines
       ins = index(in_data(loop),trim(keyword))
       if (ins==0) cycle
       in = index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s = loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s = .true.
    end do

    if (.not. found_s) then
       found = .false.
       return
    end if

    do loop=1,num_lines
       ine = index(in_data(loop),trim(keyword))
       if (ine==0) cycle
       in = index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e = loop
       if (found_e) then
          call io_error('param_get_parameters: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e = .true.
    end do

    if (.not. found_e) then
       call io_error('param_get_parameters: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if (line_e<=line_s) then
       call io_error('param_get_parameters: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    dummy = adjustl(in_data(line_s))
    ctemp = utility_strip(dummy)
    if (len(trim(ctemp)) .ne. (len(trim(keyword))+5)) &
       call io_error('Error: Unrecognised Begin block in '//trim(dummy)//' line')
    dummy = adjustl(in_data(line_e))
    ctemp = utility_strip(dummy)
    if (len(trim(ctemp)) .ne. (len(trim(keyword))+3)) & 
       call io_error('Error: Unrecognised End block in '//trim(dummy)//' line')
    
    blen = line_e-line_s-1
    found = .true.

    lconvert_length = .false.
    lconvert_energy = .false.
    if (blen == rows+1) then
      dummy = adjustl(in_data(line_s+1))
      if (len_trim(dummy)==0) call io_error('Error: keyword '//trim(keyword)//' is blank')
      found_length = .false.
      found_energy = .false.
      do i=1,2
         pos1 = scan(dummy,c_punc)
         if(pos1==0 .or. pos1==1) call io_error('Error parsing keyword '//trim(keyword)) 
         ctemp = dummy(1:pos1-1)
         select case(ctemp)
         case ('ang')
              lconvert_length = .false.
              if (found_length) &
                 call io_error('Error: Found the length unit in block '//trim(keyword)//' more than once')
              found_length = .true.
         case ('bohr')
              lconvert_length = .true.
              if (found_length)  &
                 call io_error('Error: Found the length unit in block '//trim(keyword)//' more than once')
              found_length = .true.
         case ('ev')
              lconvert_energy = .false.
              if (found_energy) &
                 call io_error('Error: Found the energy unit in block '//trim(keyword)//' more than once')
              found_energy = .true.
         case ('ryd')
              lconvert_energy = .true.
              if (found_energy) &
                 call io_error('Error: Found the energy unit in block '//trim(keyword)//' more than once')
              found_energy = .true.
         case default
              call io_error('Error: Units in block '//trim(keyword)//' not recognised')
         end select
         dummy = adjustl(dummy(pos1:))
         if(scan(dummy,c_punc) == 1) dummy = adjustl(dummy(2:))
         if(index(dummy,' ')==1) exit
      enddo
       in_data(line_s)(1:maxlen) = ' '
       line_s = line_s + 1
    endif

    counter = 0
    do line=line_s+1,line_e-1
       dummy = utility_strip(in_data(line))
       dummy = adjustl(dummy)
       counter = counter + 1
       pos1 = index(dummy,':')
       if (pos1==0) call io_error('param_read_param1: malformed types &
                    &definition: '//trim(dummy))
       ctemp = dummy(:pos1-1)
       ! Read the first atomic site
       if (index(ctemp,'t1=')>0) then
          ctemp = ctemp(4:)
          read(ctemp,*,err=105,end=105) species(1,counter)
       else
          call io_error('param_get_inp2jbd: first type not recognised '//trim(ctemp))
       end if

       !Now we know the type of first atom. Get the type of its neighbor
       
       dummy = dummy(pos1+1:)
       pos1 = index(dummy,':')
       if (pos1==0) call io_error('param_read_param1: malformed types &
                    &definition: '//trim(dummy))
       ctemp = dummy(:pos1-1)
       if (index(ctemp,'t2=')>0) then
          ctemp = ctemp(4:)
          read(ctemp,*,err=105,end=105) species(2,counter)
       else
          call io_error('param_get_inp2jbd: second type not recognised '//trim(ctemp))
       end if

       dummy = dummy(pos1+1:)
       pos1 = index(dummy,':')
       if (pos1==0) call io_error('param_get_inp2jbd: malformed types &
                    &definition: '//trim(dummy))
       ctemp = dummy(:pos1-1)
       if (index(ctemp,'sh=')>0) then
          ctemp = ctemp(4:)
          read(ctemp,*,err=105,end=105) shell(counter)
       else
          call io_error('param_get_inp2jbd: shell not recognised '//trim(ctemp))
       end if
       dummy = dummy(pos1+1:)
       if (keyword .eq. 'parameters_jij') then
           pos2 = index(dummy,'jij=')
       elseif (keyword .eq. 'parameters_bij') then
           pos2 = index(dummy,'bij=')
       elseif (keyword .eq. 'parameters_dij') then
           pos2 = index(dummy,'dij=')
       endif
       if (pos2==0) call io_error('param_read_param1: malformed parameters &
                    &definition: '//trim(dummy))
       pos1 = index(dummy,':')
       if (pos1>0) then
          ctemp = dummy(pos2+4:pos1-1)
       else
          ctemp = dummy(pos2+4:)
       endif
       read(ctemp,*,err=105,end=105) param(counter)

       if (pos1>0) then
          dummy = dummy(pos1+1:)
          if (present(sigma)) then
             pos2 = index(dummy,'sig=')
             if (pos2>0) then
                ctemp = dummy(pos2+4:)
                read(ctemp,*,err=105,end=105) sigma(counter)
             else
               call io_error('param_get_inp2jbd: sigma not recognised '//trim(ctemp))
             endif
          endif
       endif

    end do !end loop over parameters block
    if (lconvert_energy) param = param*hart*0.5_dp
    if (lconvert_energy .and. present(sigma)) sigma = sigma*hart*0.5_dp
 
    in_data(line_s:line_e)(1:maxlen) = ' '

    return

105 call io_error('param_get_inp2jbd: Problem reading ' //trim(keyword) //trim(ctemp))         

  end subroutine param_get_inp2jbd

  !===============================================================!
  subroutine param_get_axes(keyword,found,rows,c_value,r_value,param_value,lsort,lnorm)
    !================================================================!
    !                                                                !
    !   Return the axes in cartesian cordinate (Ang) and             !
    !   relative parameter without any converting the parameters     !
    !   in the case of lnorm=.true. norm of vector will be return    !
    !   if lsort = .true. change the ordeing according to the atoms  !
    !                                                                !
    !================================================================!

    use mc_constants, only : bohr,hart,eps4
    use mc_utility,   only : utility_frac_to_cart,utility_strip
    use mc_io,        only : io_error
    implicit none

    character(*),              intent(in) :: keyword
    logical,                  intent(out) :: found
    integer,                   intent(in) :: rows   
    logical,                   intent(in) :: lsort  
    logical,                   intent(in) :: lnorm  
    character(*) ,          intent(inout) :: c_value(rows)
    real(kind=dp),          intent(inout) :: r_value(3,rows)
    real(kind=dp),optional, intent(inout) :: param_value(rows)
    !
    real(kind=dp)          :: axes_frac_tmp(3,rows)
    real(kind=dp)          :: axes_cart_tmp(3,rows)
    real(kind=dp)          :: param_tmp(rows)
    character(len=maxlen)  :: atoms_label_tmp(rows)
    integer                :: in,ins,ine,loop,i,line_e,line_s,counter
    integer                :: loop2,ic,pos,blen
    character(len=maxlen)  :: dummy,end_st,start_st,ctemp
    logical                :: found_e,found_s,frac,found_length,found_energy
    logical                :: lconvert_length,lconvert_energy
    character(len=5) , parameter     :: c_punc=" ,;-:"


    frac = .false.
    dummy = trim(adjustl(keyword))
    pos = index(dummy,'frac')
    if(pos .ne. 0 ) frac = .true.

    found_s = .false.
    found_e = .false.

    start_st = 'begin '//trim(keyword)
    end_st = 'end '//trim(keyword)

    do loop=1,num_lines
       ins = index(in_data(loop),trim(keyword))
       if (ins==0) cycle
       in = index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s = loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s = .true.
    end do

    if(.not. found_s) then
       found = .false.
       return
    end if

    do loop=1,num_lines
       ine = index(in_data(loop),trim(keyword))
       if (ine==0) cycle
       in = index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e = loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e = .true.
    end do

    if (.not.found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if (line_e<=line_s) then
       call io_error('Param_get_axes: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    dummy = adjustl(in_data(line_s))
    ctemp = utility_strip(dummy)
    if (len(trim(ctemp)) .ne. (len(trim(keyword))+5)) &
       call io_error('Error: Unrecognised Begin block in '//trim(dummy)//' line')
    dummy = adjustl(in_data(line_e))
    ctemp = utility_strip(dummy)
    if (len(trim(ctemp)) .ne. (len(trim(keyword))+3)) & 
       call io_error('Error: Unrecognised End block in '//trim(dummy)//' line')
    
    blen = line_e-line_s-1
    found = .true.

    lconvert_length = .false.
    lconvert_energy = .false.
    if (blen == rows+1) then
      dummy = adjustl(in_data(line_s+1))
      if (len_trim(dummy)==0) call io_error('Error: keyword '//trim(keyword)//' is blank')
      found_length = .false.
      found_energy = .false.
      do i=1,2
         pos = scan(dummy,c_punc)
         if(pos==0 .or. pos==1) call io_error('Error parsing keyword '//trim(keyword)) 
         ctemp = dummy(1:pos-1)
         select case(ctemp)
         case ('ang')
              lconvert_length = .false.
              if (found_length) &
                 call io_error('Error: Found the length unit in block '//trim(keyword)//' more than once')
              found_length = .true.
         case ('bohr')
              lconvert_length = .true.
              if (found_length)  &
                 call io_error('Error: Found the length unit in block '//trim(keyword)//' more than once')
              found_length = .true.
         case ('ev')
              if (.not. present(param_value)) &
                   call io_error('Error: Units in block '//trim(keyword)//' not recognised')
              lconvert_energy = .false.
              if (found_energy) &
                 call io_error('Error: Found the energy unit in block '//trim(keyword)//' more than once')
              found_energy = .true.
         case ('ryd')
              if (.not. present(param_value)) &
                   call io_error('Error: Units in block '//trim(keyword)//' not recognised')
              lconvert_energy = .true.
              if (found_energy) &
                 call io_error('Error: Found the energy unit in block '//trim(keyword)//' more than once')
              found_energy = .true.
         case default
              call io_error('Error: Units in block '//trim(keyword)//' not recognised')
         end select
         dummy = adjustl(dummy(pos:))
         if(scan(dummy,c_punc) == 1) dummy = adjustl(dummy(2:))
         if(index(dummy,' ')==1) exit
      enddo
       in_data(line_s)(1:maxlen) = ' '
       line_s = line_s + 1
    endif

    axes_cart_tmp = 0
    counter = 0
    do loop=line_s+1,line_e-1
       dummy = in_data(loop)
       counter = counter+1
       if (present(param_value))then
         if (frac) then
            read(dummy,*,err=240,end=240) atoms_label_tmp(counter),&
                (axes_frac_tmp(i,counter),i=1,3),param_tmp(counter)
         else
            read(dummy,*,err=240,end=240) atoms_label_tmp(counter),&
                (axes_cart_tmp(i,counter),i=1,3),param_tmp(counter)
         end if
       else
         if (frac) then
            read(dummy,*,err=240,end=240) atoms_label_tmp(counter),(axes_frac_tmp(i,counter),i=1,3)
         else
            read(dummy,*,err=240,end=240) atoms_label_tmp(counter),(axes_cart_tmp(i,counter),i=1,3)
         end if
       endif
    end do

    if (lconvert_length)     axes_cart_tmp  = axes_cart_tmp*bohr
    if (lconvert_energy)     param_tmp      = param_tmp*hart*0.5_dp

    in_data(line_s:line_e)(1:maxlen) = ' '

    if (frac) then
      do loop=1,rows
          call utility_frac_to_cart (axes_frac_tmp(:,loop),axes_cart_tmp(:,loop),real_lattice)
      end do
    end if

    ! Now we sort the data according to the atoms 
    if (lsort) then
      counter=0
      do loop=1,num_species
         do loop2=1,num_atoms
            if( trim(atoms_symbol(loop))==trim(atoms_label_tmp(loop2) )) then
               counter = counter+1
               c_value(counter) = atoms_symbol(loop)
               r_value(:,counter) = axes_cart_tmp(:,loop2)
               if(present(param_value)) param_value(counter) = param_tmp(loop2)
            end if
         end do
      end do
      if (counter .ne. num_atoms) &
         call io_error('Error: Wrong atoms symbol at the block keyword '//trim(keyword))
      ! Atom labels (eg, si --> Si)
      do loop=1,num_atoms
         ic = ichar(c_value(loop)(1:1))
         if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
              c_value(loop)(1:1) = char(ic+ichar('Z')-ichar('z'))
      enddo
    else
      c_value = atoms_label_tmp
      if(present(param_value)) param_value = param_tmp
      do loop=1,rows
         r_value(:,loop) = axes_cart_tmp(:,loop)
      enddo
    endif

    ! Normalize Vectors
    if (lnorm) then
      do loop=1,rows
         if(all(abs(r_value(:,loop)).lt. eps4 )) &
            call io_error('Error: All component of an axe in block '//trim(keyword)//' are zero')
         r_value(:,loop) = r_value(:,loop)&
                            /sqrt(dot_product(r_value(:,loop),r_value(:,loop)))
      enddo
    endif

    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword))

  end subroutine param_get_axes

  !===========================================!
  subroutine param_memory_estimate
    !===========================================!
    !                                           !
    !! Estimate how much memory we will allocate
    !                                           !
    !===========================================!

    use mc_comms, only: on_root

    implicit none

    real(kind=dp), parameter :: size_log   = 1.0_dp
    real(kind=dp), parameter :: size_int   = 4.0_dp
    real(kind=dp), parameter :: size_int8  = 8.0_dp
    real(kind=dp), parameter :: size_real  = 8.0_dp
    real(kind=dp), parameter :: size_cmplx = 16.0_dp
    real(kind=dp) :: mem_param,mem_jij,mem_bij,mem_dij,mem_mc,mem_pt
    real(kind=dp) :: mem_param1,mem_jij1,mem_bij1,mem_dij1
    integer       :: ndim_data,ndim_sum,ndim_pm,num_freq

    mem_param  = 0.0_dp
    mem_param1 = 0.0_dp
    mem_jij    = 0.0_dp
    mem_jij1   = 0.0_dp
    mem_bij    = 0.0_dp
    mem_dij    = 0.0_dp
    mem_bij1   = 0.0_dp
    mem_dij1   = 0.0_dp
    mem_mc     = 0.0_dp
    mem_pt     = 0.0_dp

    ! First the data stored in the parameters module
    mem_param = mem_param + (num_nodes)*size_int                                     !seeds

    if (allocated(atoms_species_num)) then
      mem_param  = mem_param  + (num_species)*size_int                               !atoms_species_num
      mem_param1 = mem_param1 + (num_species)*size_real                              !atoms_label
      mem_param  = mem_param  + (num_species)*size_real                              !atoms_symbol
      mem_param  = mem_param  + (3*maxval(atoms_species_num)*num_species)*size_real  !atoms_pos_frac
      mem_param1 = mem_param1 + (3*maxval(atoms_species_num)*num_species)*size_real  !atoms_pos_cart
      mem_param1 = mem_param1 + maxval(atoms_species_num)*num_species*size_real      !mag_moments           
    endif

    mem_param1 = mem_param + 3*jij_num_nbors*size_real                              !jij_1atom_site
    mem_param1 = mem_param + 3*jij_num_nbors*size_real                              !jij_2atom_site
    mem_param1 = mem_param + jij_num_nbors*size_real                                !jij           
    mem_param1 = mem_param + jij_num_nbors*size_real                                !spin_glass_sigma           
    mem_param1 = mem_param + jij_num_nbors*size_int                                 !jij_shell

    if (have_biquad) then
       mem_param1 = mem_param1 + 3*bij_num_nbors*size_real                           !bij_1atom_site
       mem_param1 = mem_param1 + 3*bij_num_nbors*size_real                           !bij_2atom_site
       mem_param1 = mem_param1 + bij_num_nbors*size_real                             !bij           
    endif
 
 
    if (have_dm) then
       mem_param1 = mem_param1 + 3*dij_num_nbors*size_real                           !dij_1atom_site
       mem_param1 = mem_param1 + 3*dij_num_nbors*size_real                           !dij_2atom_site
       mem_param1 = mem_param1 + dij_num_nbors*size_real                             !dij           
       mem_param1 = mem_param1 + 3*dij_num_nbors*size_real                           !dij_vectors_frac
       mem_param1 = mem_param1 + 3*dij_num_nbors*size_real                           !dij_vectors_cart
    endif
 
    if (have_singleion) then
       mem_param  = mem_param  + num_atoms*size_real                                 !singleion_parameters
       mem_param  = mem_param  + 3*num_atoms*size_real                               !singleion_axes_cart
       mem_param1 = mem_param1 + num_atoms*size_real                                 !singleion_atoms_symbol
    endif
 
    if (have_field) then
       mem_param  = mem_param  + num_atoms*size_real                                 !field_parameters
       mem_param  = mem_param  + 3*num_atoms*size_real                               !field_axes_cart
       mem_param1 = mem_param1 + num_atoms*size_real                                 !field_atoms_symbol
    endif
 
    if (lorderparam) then
       mem_param  = mem_param  + 3*num_atoms*size_real                               !orderparam_axes_cart
       mem_param1 = mem_param1 + num_atoms*size_real                                 !orderparam_atoms_symbol
    endif
 
    if (lstaggered) &   
       mem_param = mem_param + num_atoms*size_real                                   !staggered_coeff        
                                                                                  
    if (lbinerror) &                                                              
       mem_param = mem_param + num_binning_level*size_real                           !binning_level          
                                                                                  
    mem_param = mem_param + tems_num*size_real                                       !tem                   
                                                                                
    mem_jij  = mem_jij + num_total_atoms*jij_num_nbors*size_real                                    !jij_matrix
    mem_jij  = mem_jij + num_total_atoms*size_real                                                  !mag_moments_matrix
    mem_jij  = mem_jij + num_total_atoms*jij_num_nbors*size_int                                     !jij_nbors_matrix
    mem_jij  = mem_jij + num_total_atoms*size_int                                                   !jij_num_nbors_matrix
    mem_jij1 = mem_jij1 + num_total_atoms*jij_num_nbors*size_real                                   !jij_tot_matrix      
    mem_jij1 = mem_jij1 + num_total_atoms*jij_num_nbors*size_int                                    !jij_tot_nbors_matrix      
    mem_jij1 = mem_jij1 + num_total_atoms*jij_num_nbors*size_int                                    !shell_matrix    
    mem_jij1 = mem_jij1 + 3*num_supercell*size_int                                                  !lmn  
    mem_jij1 = mem_jij1 + num_total_atoms*jij_num_nbors*size_real                                   !sigma_tot_matrix      
    mem_jij1 = mem_jij1 + num_total_atoms*jij_num_nbors*size_real                                   !sigma_matrix
    if (lspincorr) then
       mem_jij = mem_jij + jij_shell_max*num_total_atoms*maxval(jij_shell)*size_int                   !shell_nbors_matrix
       mem_jij = mem_jij + jij_shell_max*num_total_atoms*size_int                                     !shell_num_nbors_matrix
       mem_jij = mem_jij + jij_shell_max*num_species*num_species*size_int                             !shell_num_nbors_kind  
       mem_jij = mem_jij + num_total_atoms*size_int                                                   !kind_atoms  
    endif
 
    if (have_biquad) then 
       mem_bij  = mem_bij  + num_total_atoms*bij_num_nbors*size_real                                      !bij_matrix
       mem_bij  = mem_bij  + num_total_atoms*bij_num_nbors*size_int                                       !bij_nbors_matrix
       mem_bij  = mem_bij  + num_total_atoms*size_int                                                     !bij_num_nbors_matrix
       mem_bij1 = mem_bij1 + num_total_atoms*bij_num_nbors*size_real                                      !bij_tot_matrix      
       mem_bij1 = mem_bij1 + num_total_atoms*bij_num_nbors*size_int                                       !bij_tot_nbors_matrix      
       mem_bij1 = mem_bij1 + 3*num_supercell*size_int                                                     !lmn  
    endif
  
    if (have_dm) then 
       mem_dij  = mem_dij  + 3*num_total_atoms*dij_num_nbors*size_real                                         !dij_matrix
       mem_dij  = mem_dij  + num_total_atoms*dij_num_nbors*size_int                                            !dij_nbors_matrix
       mem_dij  = mem_dij  + num_total_atoms*size_int                                                      !dij_num_nbors_matrix
       mem_dij1 = mem_dij1 + 3.0_dp*num_total_atoms*dij_num_nbors*size_real                                  !dij_tot_matrix      
       mem_dij1 = mem_dij1 + num_total_atoms*dij_num_nbors*size_int                                          !dij_tot_matrix      
       mem_dij1 = mem_dij1 + 3*num_supercell*size_int                                                   !lmn  
    endif

    if (lorderparam.and.lstaggered) then
       ndim_data = 12; ndim_pm = 3; ndim_sum = 4
    elseif (lorderparam.or.lstaggered) then
       ndim_data = 9 ; ndim_pm = 2; ndim_sum = 3
    else
       ndim_data = 6 ; ndim_pm = 1; ndim_sum = 2
    endif
    num_freq = nint((2.0_dp*num_atoms*maxval(mag_moments))/0.01_dp)+1
  
    if (.not.pt) then
       mem_mc = mem_mc + 3*num_total_atoms*size_real                                                          !local_spin_matrix 
       if (index(initial_sconfig,'file')>0 ) &                                                                
          mem_mc = mem_mc + 3*num_total_atoms*tems_num*size_real                                              !local_spin_matrix_file 
!!         mem_mc = mem_mc + 3*num_total_atoms*int(tems_num/num_nodes)*size_real                               !local_spin_matrix_file 
       mem_mc = mem_mc + ndim_data*size_real                                                                  !local_data          
       if (lbinerror) mem_mc = mem_mc + ndim_sum*num_binning_level*size_real                                  !local_sum          
       if (lbinerror) mem_mc = mem_mc + ndim_sum*num_binning_level*size_real                                  !local_binning        
       if (lspincorr) mem_mc = mem_mc + jij_shell_max*num_species*num_species*size_real                       !local_spincorr
       if (lspincorr) mem_mc = mem_mc + jij_shell_max*num_species*num_species*size_real                       !local_spincorr_abs
       if (lsfactor) mem_mc = mem_mc + 3*sfactor_nqpts*size_real                                              !local_sfactor_matrix
       if (lsfactor) mem_mc = mem_mc + 3*num_total_atoms*size_real                                            !atoms_supercell_pos_cart
       if (energy_write) mem_mc = mem_mc + tems_num*size_int                                                  !energy_unit              
       if (energy_write) mem_mc = mem_mc + energy_num_print*size_real                                         !local_E              
       mem_mc = mem_mc + 3*num_total_atoms*num_nodes*size_real                                                !global_spin_matrix 
       mem_mc = mem_mc + num_nodes*size_real                                                                  !global_tem 
       mem_mc = mem_mc + 2*num_nodes*size_int8                                                                !global_rejected ?????????????? 
       mem_mc = mem_mc + 3*num_nodes*size_real                                                                !global_energy
       if (lspincorr) mem_mc = mem_mc + jij_shell_max*num_species*num_species*num_nodes*size_real             !global_spincorr
       if (lspincorr) mem_mc = mem_mc + jij_shell_max*num_species*num_species*num_nodes*size_real             !global_spincorr_abs
       mem_mc = mem_mc + ndim_data*num_nodes*size_real                                                        !global_data          
       if (lbinerror) mem_mc = mem_mc + ndim_sum*num_binning_level*num_nodes*size_real                        !global_binning        
       if (lsfactor) mem_mc = mem_mc + 3*sfactor_nqpts*num_nodes*size_real                                    !global_sfactor_matrix
       if (energy_write) mem_mc = mem_mc + energy_num_print*num_nodes*size_real                               !global_E              
       mem_mc = mem_mc + num_freq*size_real                                                                   !m_array      
       mem_mc = mem_mc + ndim_pm*num_freq*size_real                                                           !local_pm     
       mem_mc = mem_mc + ndim_pm*num_freq*num_nodes*size_real                                                 !global_pm     
    else
 
       mem_pt = mem_pt + 3*num_total_atoms*tems_num*size_real                                                 !local_spin_matrix 
       mem_pt = mem_pt + 3*tems_num*size_real                                                                 !tmp_spin_matrix 
       mem_pt = mem_pt + tems_num*size_int8                                                                   !local_rejected ?????????????? 
                                                                                                              
       mem_pt = mem_pt + 3*tems_num*size_real                                                              !local_energy          
       mem_pt = mem_pt + ndim_data*tems_num*size_real                                                      !local_data          
       if (lbinerror) mem_pt = mem_pt + ndim_sum*num_binning_level*tems_num*size_real                      !local_sum          
       if (lbinerror) mem_pt = mem_pt + ndim_sum*num_binning_level*tems_num*size_real                      !local_binning        
       if (lspincorr) mem_pt = mem_pt + jij_shell_max*num_species*num_species*tems_num*size_real           !local_spincorr
       if (lspincorr) mem_pt = mem_pt + jij_shell_max*num_species*num_species*tems_num*size_real           !local_spincorr_abs
       if (energy_write) mem_pt = mem_pt + tems_num*size_int                                               !energy_unit              
       if (energy_write) mem_pt = mem_pt + energy_num_print*tems_num*size_real                             !local_E              
       if (lsfactor) mem_pt = mem_pt + 3*num_total_atoms*size_real                                         !atoms_supercell_pos_cart
       if (lsfactor) mem_pt = mem_pt + 3*sfactor_nqpts*tems_num*size_real                                  !local_sfactor_matrix
                                                                                                              
       mem_pt = mem_pt + tems_num*size_int8                                                                !global_rejected ?????????????? 
       mem_pt = mem_pt + 3*num_total_atoms*tems_num*size_real                                              !global_spin_matrix 
       mem_pt = mem_pt + ndim_data*tems_num*size_real                                                      !global_data          
       mem_pt = mem_pt + 3*tems_num*size_real                                                              !global_energy
       mem_pt = mem_pt + tems_num*size_int                                                                 !map          
       mem_pt = mem_pt + 2*tems_num*size_int8                                                              !num_swaps_tot 
       if (lspincorr) mem_pt = mem_pt + jij_shell_max*num_species*num_species*num_nodes*tems_num*size_real !global_spincorr
       if (lspincorr) mem_pt = mem_pt + jij_shell_max*num_species*num_species*num_nodes*tems_num*size_real !global_spincorr_abs
       if (lbinerror) mem_pt = mem_pt + ndim_sum*num_binning_level*tems_num*size_real                      !global_binning        
       if (lsfactor) mem_pt = mem_pt + 3*sfactor_nqpts*tems_num*size_real                                  !global_sfactor_matrix
       if (energy_write) mem_pt = mem_pt + energy_num_print*tems_num*size_real                             !global_E              
       mem_pt = mem_pt + num_freq*size_real                                                                !m_array      
       mem_pt = mem_pt + ndim_pm*num_freq*int(tems_num/num_nodes)*size_real                                !local_pm     
       mem_pt = mem_pt + ndim_pm*num_freq*tems_num*size_real                                               !global_pm     
    endif
    if (on_root) then
      write (stdout, '(1x,a)') '*===================================================================================*'
      write (stdout, '(1x,a)') '|                                 MEMORY ESTIMATE                                   |'
      write (stdout, '(1x,a)') '|            Maximum RAM allocated during each phase of the calculation             |'
      write (stdout, '(1x,a)') '*===================================================================================*'
      write (stdout, '(1x,"|",24x,a17,f16.2,a,23x,"|")') 'Jij:',&
            (mem_param + mem_param1 + mem_jij + mem_jij1)/(1024**2), ' Mb'
      if (have_biquad) write (stdout, '(1x,"|",24x,a17,f16.2,a,23x,"|")') 'Bij:',&
                             (mem_param + mem_param1 + mem_jij + mem_bij + mem_bij1)/(1024**2), ' Mb'
      if (have_dm) write (stdout, '(1x,"|",24x,a17,f16.2,a,23x,"|")') 'Dij:',&
                  (mem_param + mem_param1 + mem_jij + mem_bij + mem_dij + mem_dij1)/(1024**2), ' Mb'
      if (.not. pt) then
        write (stdout, '(1x,"|",24x,a17,f16.2,a,23x,"|")') 'Monte-Carlo:', &
                        (mem_param + mem_jij + mem_bij + mem_dij + mem_mc)/(1024**2), ' Mb'
      else
        write (stdout, '(1x,"|",24x,a17,f16.2,a,23x,"|")') 'Parall-Tempering:',&
                        (mem_param + mem_jij + mem_bij + mem_dij + mem_pt)/(1024**2), ' Mb'
      endif
      write (stdout, '(1x,a)') '|                                                                                   |'
      write (stdout, '(1x,a)') '*-----------------------------------------------------------------------------------*'
      write (stdout, *) ' '
    endif

    return
  end subroutine param_memory_estimate


end module mc_parameters
