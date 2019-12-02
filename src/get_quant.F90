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

module mc_get_quant     

  use mc_constants,  only : dp,zero
  use mc_parameters, only : num_atoms,num_total_atoms,state
  use stdtypes
  use mtprng,        only : mtprng_init,mtprng_rand_real2,mtprng_state

  implicit none

  private 

  public :: metropolis 
  public :: energy_get
  public :: delta_energy_get
  public :: magnetization_get
  public :: staggered_get
  public :: orderparam_get
  public :: spincorr_get
  public :: Pm_get
  public :: mcarlo_read_sconfig
  public :: mcarlo_get_atoms_pos
  public :: mcarlo_structure_factor
  public :: mcarlo_structure_factor_write_mc
  public :: mcarlo_structure_factor_write_pt
  public :: open_E_unit
  public :: write_E 
  public :: close_E_unit

contains

    !==================================================================!
    subroutine metropolis(loop_at,beta,spin_matrix,energy_oc,rejected)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  

      use mc_constants,  only : twopi,eps8
      use mc_jij,        only : mag_moments_matrix
      use mc_parameters, only : mcarlo_mode,tilt_angles_max

      implicit none
      integer,       intent(in)    :: loop_at
      real(kind=dp), intent(in)    :: beta 
      real(kind=dp), intent(inout) :: spin_matrix(3,num_total_atoms) 
      real(kind=dp), intent(inout) :: energy_oc 
      integer(8),    intent(inout) :: rejected
      real(kind=dp)                :: random
      real(kind=dp)                :: costheta,sintheta, cosphi,sinphi
      real(kind=dp)                :: phi,theta,fac
      real(kind=dp)                :: spin_old(3),mag
      real(kind=dp)                :: delta_energy
      real(kind=dp)                :: delta_spin_matrix(3)
      real(kind=dp)                :: spin_old_norm(3),spinxy,spinx,spiny,spinz

      mag = mag_moments_matrix(loop_at)
      spin_old(:) = spin_matrix(:,loop_at)

      if (index(mcarlo_mode,'rand').eq. 1) then
        random = mtprng_rand_real2(state)
!        call random_number(random)
        costheta = 2.0_dp*random-1.0_dp
        sintheta = sqrt(1.0_dp-costheta**2)
     
        random = mtprng_rand_real2(state)
!        call random_number(random)
        phi = random*twopi

        spin_matrix(1,loop_at) = mag*sintheta*cos(phi)
        spin_matrix(2,loop_at) = mag*sintheta*sin(phi)
        spin_matrix(3,loop_at) = mag*costheta
        delta_spin_matrix(:) = spin_matrix(:,loop_at)-spin_old(:)
      elseif (index(mcarlo_mode,'const') > 0) then

        ! The idea of coonstrain mode is to increase acceptance ratio. 
        ! This can be done by restriction of angles theta and phi, but how we can do this without
        ! lossing the uniform distribution in a spherical!. 
        ! we can reach to this goole by first producting uniform distribution on sphere and 
        ! then transfer them to a restricted angle by multiply \theta to number less than 1, here is "tilt_angles_max". 
        ! (Please see Fig.3 of  Journal of Computational Science 5 (2014) 696–700).
        ! After generating the random phi and theta, then we rotate them to spin_old direction.
        ! In this way we generate new direction for spin which are near to spin_old. If tilt_angles_max incrases, then 
        ! the step change in each MC step increases. 

        random = mtprng_rand_real2(state)
        theta = tilt_angles_max*acos(2.0_dp*random-1.0_dp)

        random = mtprng_rand_real2(state)
        phi = random*twopi
       
        sintheta = sin(theta) 
        spinx = sintheta*cos(phi)
        spiny = sintheta*sin(phi)
        spinz = cos(theta)

        !phi0=atan2(spin_old(2),spin_old(1))
        !c=cos(phi0)
        !s=sin(phi0)
        spinxy = sqrt(spin_old(1)**2+spin_old(2)**2)
        if (spinxy > eps8) then
             cosphi = spin_old(1)/spinxy
             sinphi = spin_old(2)/spinxy
        else
            cosphi = 1.0_dp
            sinphi = 0.0_dp
        endif  

        !spin_old_norm(:)=1.d0/mag*spin_old(:) 
        ! we re-normilzed spin_old becuase after applying the rotation operators many times on spins 
        ! we should expect of round-off error in spin sizes. So to remove this effect we re-nomalize spins.
        spin_old_norm(:) = 1.0_dp/norm2(spin_old)*spin_old(:) 
        
            
        !The rotation matrix decompose to two rotation matrix: M=M(phi)_z * M(theta)_y
        ! M(phi)_z =  [[cos(phi), -sin(phi), 0] ,[sin(phi), cos(phi)], [0; 0,0,1]]
        ! M(theta)_y = [[cos(theta), 0, sin(theta)],[0, 1, 0], [-sin(theta), 0, cos(theta)]]
        ! M=[[cos(phi)*cos(theta), -sin(phi), cos(phi)*sin(theta)], [sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta)], [-sin(theta),0, cos(theta)]]
        ! we can summerize them as follows:
        !Mxx=cosphi*spin_old_norm(3)
        !Mxy=-sinphi
        !Mxz=spin_old_norm(1)

        !Myx=sinphi*spin_old_norm(3)
        !Myy=cosphi
        !Myz=spin_old_norm(2)

        !Mzx=-sqrt(spin_old_norm(1)**2+spin_old_norm(2)**2)
        !Mzy=0
        !Mzz=spin_old_norm(3)
         
        !rotated spin
        !sp_rx=Mxx*sp_x+Mxy*sp_y+Mxz*sp_z
        !sp_ry=Myx*sp_x+Myy*sp_y+Myz*sp_z
        !sp_rz=Mzx*sp_x+Mzy*sp_y+Mzz*sp_z


        spin_matrix(1,loop_at) = mag*(cosphi*spin_old_norm(3)*spinx-sinphi*spiny+spin_old_norm(1)*spinz)
        spin_matrix(2,loop_at) = mag*(sinphi*spin_old_norm(3)*spinx+cosphi*spiny+spin_old_norm(2)*spinz)
        spin_matrix(3,loop_at) = mag*(-sqrt(spin_old_norm(1)**2+spin_old_norm(2)**2)*spinx+spin_old_norm(3)*spinz)
        delta_spin_matrix(:) = spin_matrix(:,loop_at)-spin_old(:)


      endif

      call delta_energy_get(loop_at,spin_matrix,delta_spin_matrix,delta_energy)
      energy_oc = energy_oc+delta_energy

      if (delta_energy .gt. zero) then
        fac = exp(-beta*delta_energy)
        random = mtprng_rand_real2(state)
!        call random_number(random)
        if (fac .lt. random) then
           spin_matrix(:,loop_at) = spin_old(:)
           energy_oc = energy_oc-delta_energy
           rejected = rejected+1
        endif
      endif

      return

    end subroutine metropolis


   !==================================================================!
    subroutine energy_get(spin_matrix,energy_tot)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  

      use mc_jij,        only : jij_nbors_matrix,jij_num_nbors_matrix,&
                                jij_matrix
      use mc_bij,        only : bij_nbors_matrix,bij_num_nbors_matrix,bij_matrix
      use mc_dij,        only : dij_nbors_matrix,dij_num_nbors_matrix,dij_matrix
      use mc_utility,    only : utility_cross
      use mc_parameters, only : singleion_parameters,singleion_axes_cart,&
                                field_parameters,field_axes_cart,&
                                have_dm,have_biquad,have_singleion,have_field

      implicit none

      integer                    :: i,j,loop,axis
      real(kind=dp), intent(in)  :: spin_matrix(3,num_total_atoms)
      real(kind=dp), intent(out) :: energy_tot
      real(kind=dp)              :: sidotsj,sicrosssj(3)
      real(kind=dp)              :: sidotaxis,sidotfield

      energy_tot = zero
      do i=1,num_total_atoms
         do loop=1,jij_num_nbors_matrix(i)
            j=jij_nbors_matrix(i,loop)
            energy_tot=energy_tot-jij_matrix(i,loop)*dot_product(spin_matrix(:,i),spin_matrix(:,j))
         enddo
      enddo

      if (have_field) then
        do i=1,num_total_atoms
           axis = mod(i,num_atoms)
           if (axis .eq. 0) axis=num_atoms
           sidotfield = dot_product(spin_matrix(:,i),field_axes_cart(:,axis))
           energy_tot = energy_tot+2.0_dp*field_parameters(axis)*sidotfield
        enddo
      endif

      if (have_singleion) then
        do i=1,num_total_atoms
           axis = mod(i,num_atoms)
           if (axis .eq. 0) axis=num_atoms
           sidotaxis = dot_product(spin_matrix(:,i),singleion_axes_cart(:,axis))
           energy_tot = energy_tot+2.0_dp*singleion_parameters(axis)*sidotaxis*sidotaxis
        enddo
      endif

      if (have_biquad) then
         do i=1,num_total_atoms
            do loop=1,bij_num_nbors_matrix(i)
               j = bij_nbors_matrix(i,loop)
               sidotsj = dot_product(spin_matrix(:,i),spin_matrix(:,j))
               energy_tot = energy_tot+bij_matrix(i,loop)*sidotsj*sidotsj
            enddo
         enddo
      endif

      if (have_dm) then
         do i=1,num_total_atoms
            do loop=1,dij_num_nbors_matrix(i)
               j = dij_nbors_matrix(i,loop)
               call utility_cross(spin_matrix(:,i),spin_matrix(:,j),sicrosssj)
               energy_tot = energy_tot+dot_product(dij_matrix(:,i,loop),sicrosssj(:))
            enddo
         enddo
      endif


      energy_tot = energy_tot/2.0_dp

      return

    end subroutine energy_get

    !==================================================================!
    subroutine delta_energy_get(loop_at,spin_matrix,delta_spin_matrix,delta_energy)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  

      use mc_jij,        only : jij_nbors_matrix,jij_num_nbors_matrix,&
                                jij_matrix
      use mc_bij,        only : bij_nbors_matrix,bij_num_nbors_matrix,bij_matrix
      use mc_dij,        only : dij_nbors_matrix,dij_num_nbors_matrix,dij_matrix
      use mc_utility,    only : utility_cross
      use mc_parameters, only : singleion_parameters,singleion_axes_cart,&
                                field_parameters,field_axes_cart,&
                                have_dm,have_biquad,have_singleion,have_field


      implicit none

      integer,       intent(in)  :: loop_at
      real(kind=dp), intent(in)  :: spin_matrix(3,num_total_atoms)
      real(kind=dp), intent(in)  :: delta_spin_matrix(3)
      real(kind=dp), intent(out) :: delta_energy
      integer                    :: j,loop,axis
      real(kind=dp)              :: sidotdeltasj,si(3),sidotsj
      real(kind=dp)              :: sicrossdeltasj(3)
      real(kind=dp)              :: sidotaxis,deltasidotaxis

      delta_energy = zero

      do loop=1,jij_num_nbors_matrix(loop_at)
         j = jij_nbors_matrix(loop_at,loop)
         delta_energy = delta_energy-jij_matrix(loop_at,loop)*dot_product(delta_spin_matrix(:),spin_matrix(:,j))
      enddo

      if (have_field)then
        axis = mod(loop_at,num_atoms)
        if (axis.eq.0) axis = num_atoms
        delta_energy = delta_energy+field_parameters(axis)*&
                                  dot_product(delta_spin_matrix(:),field_axes_cart(:,axis))
      endif

      if (have_singleion) then
        axis = mod(loop_at,num_atoms)
        if (axis.eq.0) axis = num_atoms
        si = spin_matrix(:,loop_at)-delta_spin_matrix(:)
        sidotaxis = dot_product(si(:),singleion_axes_cart(:,axis))
        deltasidotaxis = dot_product(delta_spin_matrix(:),singleion_axes_cart(:,axis))
        delta_energy = delta_energy+singleion_parameters(axis)*(deltasidotaxis**2+&
                                    2.0_dp*sidotaxis*deltasidotaxis)
      endif

      if (have_biquad) then
        do loop=1,bij_num_nbors_matrix(loop_at)
           j = bij_nbors_matrix(loop_at,loop)
           si = spin_matrix(:,loop_at)-delta_spin_matrix(:)
           sidotsj = dot_product(si(:),spin_matrix(:,j))
           sidotdeltasj = dot_product(delta_spin_matrix(:),spin_matrix(:,j))
           delta_energy = delta_energy+bij_matrix(loop_at,loop)*(sidotdeltasj*sidotdeltasj+&
                                                            2.0_dp*sidotsj*sidotdeltasj)
        enddo
      endif

      if (have_dm) then
         do loop=1,dij_num_nbors_matrix(loop_at)
            j = dij_nbors_matrix(loop_at,loop)
            call utility_cross(delta_spin_matrix,spin_matrix(:,j),sicrossdeltasj)
            delta_energy = delta_energy+dot_product(dij_matrix(:,loop_at,loop),sicrossdeltasj(:))
         enddo
       endif

      return

    end subroutine delta_energy_get

    !==================================================================!
    subroutine magnetization_get(spin_matrix,m2_oc)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  

      implicit none

      real(kind=dp),intent(in)  :: spin_matrix(3,num_total_atoms)
      real(kind=dp),intent(out) :: m2_oc
      integer                   :: j
      real(kind=dp)             :: m_oc(3)

      m_oc = zero

      do j=1,num_total_atoms
         m_oc(:) = m_oc(:)+spin_matrix(:,j)/real(num_total_atoms,dp)
      enddo
      m2_oc = dot_product(m_oc,m_oc)

      return

    end subroutine magnetization_get

    !==================================================================!
    subroutine staggered_get(spin_matrix,staggered_m2)
      !==================================================================!
      !                                                                  !
      !            ST=|\sum{(S_i.op_axes_i)}/(#atoms)|                        ! 
      !                                                                  !
      !===================================================================  
      use mc_parameters, only : staggered_coeff

      implicit none

      real(kind=dp),intent(in)  :: spin_matrix(3,num_total_atoms)
      real(kind=dp),intent(out) :: staggered_m2
      integer                   :: j,staggered_natom
      real(kind=dp)             :: staggered_m(3)

      staggered_m = zero

      do j=1,num_total_atoms
         staggered_natom = mod(j,num_atoms)
         if (staggered_natom.eq. 0) staggered_natom = num_atoms
         staggered_m(:) = staggered_m(:)+spin_matrix(:,j)*staggered_coeff(staggered_natom)/real(num_total_atoms,dp)
      enddo
      staggered_m2 = dot_product(staggered_m,staggered_m)

      return

    end subroutine staggered_get


    !==================================================================!
    subroutine orderparam_get(spin_matrix,orderparam2)
      !==================================================================!
      !                                                                  !
      !            OP=|\sum_{i}{(S_i.op_axes_i)}/(#atoms)|               ! 
      !                                                                  !
      !===================================================================  
      use mc_parameters, only : orderparam_axes_cart

      implicit none

      real(kind=dp),intent(in)  :: spin_matrix(3,num_total_atoms)
      real(kind=dp),intent(out) :: orderparam2
      integer                   :: j,axis
      real(kind=dp)             :: orderparam 

      orderparam = zero

      do j=1,num_total_atoms
         axis = mod(j,num_atoms)
         if (axis.eq. 0) axis=num_atoms
         orderparam = orderparam+dot_product(spin_matrix(:,j),orderparam_axes_cart(:,axis))/real(num_total_atoms,dp)
      enddo
      orderparam2 = orderparam*orderparam

      return

    end subroutine orderparam_get

    !====================================================================!
    subroutine Pm_get(spin_matrix,num_freq,pdf_m,pdf_order,pdf_staggered)
      !====================================================================!
      !                                                                    !
      !         Probility distribution function magnetization              ! 
      !         & order parameter                                          !
      !                                                                    !
      !====================================================================!  
 
      use mc_jij,        only : mag_moments_matrix
      use mc_parameters, only : orderparam_axes_cart,staggered_coeff

      implicit none
 
      real(kind=dp),           intent(in)    :: spin_matrix(3,num_total_atoms)
      integer,                 intent(in)    :: num_freq
      real(kind=dp),           intent(inout) :: pdf_m(num_freq)
      real(kind=dp), optional, intent(inout) :: pdf_order(num_freq)
      real(kind=dp), optional, intent(inout) :: pdf_staggered(num_freq)
      integer                                :: j,axis,m_p,i
      real(kind=dp)                          :: orderparam,m2_oc,m_oc(3)
      real(kind=dp)                          :: m_max,d_m,m_min
      real(kind=dp)                          :: staggered_m(3),staggered_m2
 
      m_max = zero; d_m = 0.01_dp
      do i=1,num_atoms
         m_max = m_max+abs(mag_moments_matrix(i))
      enddo
      m_min = -m_max

      do j=1,num_total_atoms,num_atoms
         m_oc = zero
         orderparam = zero
         staggered_m = zero

         do i=0,num_atoms-1
            m_oc(:) = m_oc(:)+spin_matrix(:,j+i)
         enddo
         m2_oc = dot_product(m_oc,m_oc)
         m_p = max(nint(((sqrt(m2_oc)-m_min)/(m_max-m_min)) &
                 *real(num_freq-1,kind=dp)) + 1, 1)
         pdf_m(m_p) = pdf_m(m_p)+1.0_dp/d_m

         if (present(pdf_order)) then
           do i=0,num_atoms-1
              axis = mod(j+i,num_atoms)
              if(axis.eq.0) axis = num_atoms
              orderparam = orderparam+dot_product(spin_matrix(:,j+i),orderparam_axes_cart(:,axis))
           enddo
           m_p = max(nint(((orderparam-m_min)/(m_max-m_min)) &
                 *real(num_freq-1,kind=dp)) + 1, 1)
           pdf_order(m_p) = pdf_order(m_p)+1.0_dp/d_m
         endif

         if (present(pdf_staggered)) then
            do i=0,num_atoms-1
               axis = mod(j+i,num_atoms)
               if(axis.eq.0) axis = num_atoms
               staggered_m(:) = staggered_m(:)+spin_matrix(:,j+i)*staggered_coeff(axis)
            enddo
            staggered_m2 = dot_product(staggered_m,staggered_m)
            m_p = max(nint(((sqrt(staggered_m2)-m_min)/(m_max-m_min)) &
                  *real(num_freq-1,kind=dp)) + 1, 1)
            pdf_staggered(m_p) = pdf_staggered(m_p)+1.0_dp/d_m
         endif
      enddo

      return
 
    end subroutine Pm_get
 
   !==================================================================!
    subroutine spincorr_get(spin_matrix,spincorr,spincorr_abs,spincorr_tot_abs)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  

      use mc_parameters, only : jij_shell_max,num_species
      use mc_jij,        only : shell_nbors_matrix,shell_num_nbors_matrix,&
                                kind_atoms,shell_num_nbors_kind,mag_moments_matrix

      implicit none

      real(kind=dp), intent(in  )  :: spin_matrix(3,num_total_atoms)
      real(kind=dp), intent(inout) :: spincorr(jij_shell_max,num_species,num_species)
      real(kind=dp), intent(inout) :: spincorr_abs(jij_shell_max,num_species,num_species)
      real(kind=dp), intent(inout) :: spincorr_tot_abs(jij_shell_max,num_species,num_species)
      real(kind=dp)                :: spincorr_tmp(num_species,num_species)
      real(kind=dp)                :: spincorr_abs_tmp(num_species,num_species)
      integer                      :: i,j,loop,ishell,nsp_i,nsp_j

      do ishell=1,jij_shell_max
         spincorr_tmp = 0.0_dp
         spincorr_abs_tmp = 0.0_dp
         do i=1,num_total_atoms
            nsp_i = kind_atoms(i)
            do loop=1,shell_num_nbors_matrix(ishell,i)
               j = shell_nbors_matrix(ishell,i,loop)
               nsp_j = kind_atoms(j)
               spincorr_tmp(nsp_i,nsp_j) = spincorr_tmp(nsp_i,nsp_j)+&
                                     dot_product(spin_matrix(:,i),spin_matrix(:,j))/ &
                                     (shell_num_nbors_kind(ishell,nsp_i,nsp_j)*mag_moments_matrix(i)*mag_moments_matrix(j))                 
               spincorr_abs_tmp(nsp_i,nsp_j) = spincorr_abs_tmp(nsp_i,nsp_j)+&
                                   abs(dot_product(spin_matrix(:,i),spin_matrix(:,j)))/ &
                                   (shell_num_nbors_kind(ishell,nsp_i,nsp_j)*mag_moments_matrix(i)*mag_moments_matrix(j))                 
            enddo
         enddo
         spincorr(ishell,:,:)         = spincorr(ishell,:,:)         + spincorr_tmp
         spincorr_tot_abs(ishell,:,:) = spincorr_tot_abs(ishell,:,:) + abs(spincorr_tmp)
         spincorr_abs(ishell,:,:)     = spincorr_abs(ishell,:,:)     + spincorr_abs_tmp
      enddo

      return

    end subroutine spincorr_get
!
    !================================================================!
    subroutine mcarlo_read_sconfig(nxx,t,spin,spin_file,num_points)
      !================================================================!
      !                                                                !
      !                                                                !
      !                                                                !
      !================================================================!

      use mc_constants, only : eps4
      use mc_io,        only : stdout, io_file_unit, io_error, maxlen

      implicit none

      integer,           intent(in) :: nxx
      integer,           intent(in) :: num_points
      real(kind=dp),     intent(in) :: t(num_points)
      real(kind=dp),    intent(out) :: spin(3,nxx,num_points)
      character(len=50), intent(in) :: spin_file
      ! 
      integer               :: i,j,nw,file_unit
      integer               :: nt,k
      real(kind=dp)         :: tmp_t
      character(len=maxlen) :: dummy

      file_unit = io_file_unit()

      open(unit=file_unit, file=spin_file, form='formatted', &
           status='old', action='read',err=101)

      write(stdout,'(/a)',advance='no') ' Reading S matrix from '//trim(spin_file)//'  : '

      read(file_unit,'(a)',err=102,end=102) dummy
      write(stdout, '(a)') trim(dummy)
      read(file_unit,'(a)',err=102,end=102) dummy ! supercell size line

      read(file_unit,*,err=102,end=102) nw,nt ! number of total atoms & number of temperatures
      if (nw.ne.nxx) call io_error('Error: Wrong matrix size in mcarlo: read_sconfig')
      if (nt.ne.num_points) call io_error('Error: Wrong temperature number in mcarlo: read_sconfig')

      do k=1,nt
         read(file_unit,*,err=102,end=102) tmp_t
         if (abs(tmp_t-t(k)) .gt. eps4 ) call io_error('Error: Mismatch in temperatures: read_sconfig')
         do j=1,nxx
            read(file_unit,*,err=102,end=102) (spin(i,j,k),i=1,3)
         if(all(abs(spin(:,j,k)).lt. eps4 )) &
            call io_error('Error: All components of an spin in block are zero in '//trim(spin_file))
            spin(:,j,k) = spin(:,j,k)/sqrt(dot_product(spin(:,j,k),spin(:,j,k)))
         enddo
      enddo

      close(unit=file_unit)

      return

101    call io_error('Error: Problem opening input file '//spin_file)
102    call io_error('Error: Problem reading input file '//spin_file)
    end subroutine mcarlo_read_sconfig

    !==================================================================!
    subroutine mcarlo_get_atoms_pos(pos_supecell_cart)
      !==================================================================!
      ! This subroutine calculate the position of all atoms at supercell !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  
      use mc_io,         only : io_error
      use mc_utility,    only : utility_frac_to_cart
      use mc_parameters, only : atoms_pos_frac,num_supercell,&
                                supercell_size,real_lattice,num_species,&
                                atoms_species_num
      implicit none

      integer                    :: nnx,n,m,l
      integer                    :: nsp,nat,isupercell
      integer                    :: lmn(3,num_supercell)
      real(kind=dp)              :: pos_frac(3)
      real(kind=dp),intent(out)  :: pos_supecell_cart(3,num_total_atoms)


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

      nnx = 0
      do isupercell=1,num_supercell  ! Supercell
         do nsp=1,num_species
            do nat=1,atoms_species_num(nsp) ! these two loops determine the atoms at the home unit cell
               nnx = nnx+1
               pos_frac(:) = atoms_pos_frac(:,nsp,nat)+lmn(:,isupercell)
               call utility_frac_to_cart(pos_frac(:),pos_supecell_cart(:,nnx),real_lattice)
            enddo
         enddo
      enddo

      return
    end subroutine mcarlo_get_atoms_pos

    !==================================================================!
    subroutine mcarlo_structure_factor(spin_matrix,atoms_supercell_pos_cart,sfactor_matrix)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !,===================================================================  
      use mc_constants,  only : cmplx_0,cmplx_i
      use mc_parameters, only : sfactor_corner,sfactor_q1,sfactor_q2,&
                                sfactor_polar,sfactor_2dqmesh,sfactor_nqpts
      use mc_utility,    only : utility_cross

      implicit none

      real(kind=dp), intent(in)    :: spin_matrix(3,num_total_atoms)               
      real(kind=dp), intent(in)    :: atoms_supercell_pos_cart(3,num_total_atoms)               
      real(kind=dp), intent(inout) :: sfactor_matrix(3,sfactor_nqpts) 
      integer                      :: loop_xy,loop_x,loop_y,loop
      real(kind=dp)                :: qpt(3),q1,q2
      real(kind=dp)                :: multiple_q(3),qdotr
      real(kind=dp)                :: pcrossq(3),SF,NSF
      complex(kind=dp)             :: fac,Fq(3)

      ok:do loop_xy=0,sfactor_nqpts-1
         loop_x = loop_xy/(sfactor_2dqmesh(2) + 1)
         loop_y = loop_xy-loop_x*(sfactor_2dqmesh(2) + 1)
         ! q1 and q2 are the coefficients of the q-point in the basis
         ! (sfactor_q1,sfactor_q2)
         q1 = loop_x/real(sfactor_2dqmesh(1),dp)
         q2 = loop_y/real(sfactor_2dqmesh(2),dp)
         ! I need the qpt in catesian coordinate
         !qpt(:) = qvec(3,:)+q1*qvec(1,:)+q2*qvec(2,:)
         qpt(:) = sfactor_corner+q1*sfactor_q1+q2*sfactor_q2
         if (all(qpt(:) .eq. zero) ) then
           Fq = cmplx_0
           cycle ok
         endif
         !call utility_cross(pvec,qpt,pcrossq)
         call utility_cross(sfactor_polar,qpt,pcrossq)
         Fq = cmplx_0
         NSF = zero
         SF = zero
         do loop=1,num_total_atoms
            qdotr = dot_product(atoms_supercell_pos_cart(:,loop),qpt)
            fac = exp(cmplx_i*qdotr)
            ! Non polarize
            multiple_q = dot_product(spin_matrix(:,loop),qpt(:))/dot_product(qpt,qpt)
            Fq(:) = Fq(:)+(spin_matrix(:,loop)-multiple_q*qpt(:))*fac
            ! Spin Flip 
            !NSF = NSF+dot_product(pvec,spin_matrix(:,loop)-multiple_q*qpt(:))*fac/sqrt((dot_product(pvec,pvec)))
            NSF = NSF + dot_product(sfactor_polar,spin_matrix(:,loop)-multiple_q*qpt(:))*fac&
                        /sqrt((dot_product(sfactor_polar,sfactor_polar)))
!            SF=SF+dot_product(spin_matrix(:,loop)*fac,pcrossq(:))/ &
!                  (sqrt(dot_product(pvec,pvec))*sqrt(dot_product(qpt,qpt)))
            ! Non Spin Flip(NSF)
!            NSF=NSF+dot_product(spin_matrix(:,loop)*fac,pvec(:))/&
!                    sqrt(dot_product(pvec,pvec))
         enddo

         sfactor_matrix(1,loop_xy+1) = sfactor_matrix(1,loop_xy+1)+real(dot_product(Fq,Fq),dp)
!         sfactor_matrix(2,loop_xy+1)=sfactor_matrix(2,loop_xy+1)+SF*SF
         sfactor_matrix(2,loop_xy+1) = sfactor_matrix(2,loop_xy+1)+&
                                       real(dot_product(Fq,Fq),dp)-NSF*NSF
         sfactor_matrix(3,loop_xy+1) = sfactor_matrix(3,loop_xy+1)+NSF*NSF

      enddo ok

      return

    end subroutine mcarlo_structure_factor

    !==================================================================!
    subroutine mcarlo_structure_factor_write_mc(sfactor_matrix,start)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  
      use mc_io,         only : stdout,io_file_unit,seedname
      use mc_constants,  only : twopi                        
      use mc_utility,    only : utility_recip_lattice,utility_cross 
      use mc_parameters, only : sfactor_corner,sfactor_q1,sfactor_q2,&
                                sfactor_2dqmesh,sfactor_nqpts,&
                                tems,tems_num,lsfactor_polar
      use mc_comms,      only : num_nodes

      implicit none

      real(kind=dp), intent(in)  :: sfactor_matrix(3,sfactor_nqpts,0:num_nodes-1) 
      integer,       intent(in)  :: start
      integer                    :: loop_xy,loop_x,loop_y,i
      integer                    :: script_unit,n,n1,n2,n3
      integer,       allocatable :: sfactor_unit(:)
      real(kind=dp)              :: qvec(3,3)
      real(kind=dp)              :: avec_2d(3,3),areaq1q2,rdum
      real(kind=dp)              :: q1mod,q2mod,ymod,cosq1q2,cosyq2
      real(kind=dp)              :: yvec(3),zvec(3)
      real(kind=dp)              :: qpt_x,qpt_y
      real(kind=dp)              :: q1,q2
      character(len=40)          :: filename

      qvec(1,:) = sfactor_q1
      qvec(2,:) = sfactor_q2

      ! z_vec (orthogonal to the slice)
      call utility_cross(qvec(1,:),qvec(2,:),zvec)
      ! y_vec (orthogonal to q1=x_vec)
      call utility_cross(zvec,qvec(1,:),yvec)
      ! Moduli b1,b2,y_vec
      q1mod = sqrt(qvec(1,1)**2+qvec(1,2)**2+qvec(1,3)**2)
      q2mod = sqrt(qvec(2,1)**2+qvec(2,2)**2+qvec(2,3)**2)
      ymod = sqrt(yvec(1)**2+yvec(2)**2+yvec(3)**2)
      areaq1q2 = sqrt(zvec(1)**2+zvec(2)**2+zvec(3)**2)
      qvec(3,:) = zvec(:)/areaq1q2

      ! Now that we have qvec(3,:), we can compute the dual vectors
      ! avec_2d as in the 3D case
      call utility_recip_lattice(qvec,avec_2d,rdum)

      ! Cosine of the angle between y_vec and b2
      cosyq2 = yvec(1)*qvec(2,1)+yvec(2)*qvec(2,2)+yvec(3)*qvec(2,3)
      cosyq2 = cosyq2/(ymod*q2mod)
      ! Cosine of the angle between b1=x_vec and b2
      cosq1q2 = qvec(1,1)*qvec(2,1)+qvec(1,2)*qvec(2,2)+qvec(1,3)*qvec(2,3)
      cosq1q2 = cosq1q2/(q1mod*q2mod)

      allocate(sfactor_unit(0:num_nodes-1))
      i = 0
      do n=(start-1)*num_nodes+1,start*num_nodes
         n1 = n/100
         n2 = (n-n1*100)/10
         n3 = n-n1*100-n2*10
         sfactor_unit(i) = io_file_unit()
         filename = trim(seedname)//'_sf-T'&
                    //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat'
         open(sfactor_unit(i),file=filename,status='unknown',form='formatted',action='write')
          i = i+1
      enddo

      if(start*num_nodes .eq. tems_num) then
        do n=1,start*num_nodes
           n1 = n/100
           n2 = (n-n1*100)/10
           n3 = n-n1*100-n2*10
           filename = trim(seedname)//'_sf-T'&
                      //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat'
           write(stdout,'(3x,a,f8.3,3x,a,3x,a)') 'Output Structure factor at T = ',tems(n),'writted at',filename
        enddo
      endif

      do n=0,num_nodes-1
         do loop_xy=0,sfactor_nqpts-1
            loop_x = loop_xy/(sfactor_2dqmesh(2) + 1)
            loop_y = loop_xy-loop_x*(sfactor_2dqmesh(2) + 1)
            ! q1 and q2 are the coefficients of the q-point in the basis
            ! (sfactor_q1,sfactor_q2)
            q1 = loop_x/real(sfactor_2dqmesh(1),dp)
            q2 = loop_y/real(sfactor_2dqmesh(2), dp)
            ! Add to (q1,q2) the projection of sfactor_corner on the
            ! (sfactor_q1,sfactor_q2) plane, expressed as a linear
            ! combination of sfactor_q1 and sfactor_q2 
            q1 = q1+dot_product(sfactor_corner,avec_2d(1,:))/twopi
            q2 = q2+dot_product(sfactor_corner,avec_2d(2,:))/twopi

            qpt_x = q1*q1mod+q2*q2mod*cosq1q2
            qpt_y = q2*q2mod*cosyq2

            if (lsfactor_polar) then
              write(sfactor_unit(n),136) qpt_x,qpt_y,sfactor_matrix(1,loop_xy+1,n),&
                    sfactor_matrix(2,loop_xy+1,n),sfactor_matrix(3,loop_xy+1,n)
            else
              write(sfactor_unit(n),135) qpt_x,qpt_y,sfactor_matrix(1,loop_xy+1,n)
            endif
         enddo
         close(sfactor_unit(n))
      enddo

      filename = trim(seedname)//'_sf.py'
      open(script_unit,file=filename,status='unknown',form='formatted',action='write')
      write(script_unit,'(a)') 'import pylab as pl'
      write(script_unit,'(a)') 'import numpy as np'
      write(script_unit,'(a)') 'import matplotlib.mlab as ml'
      write(script_unit,'(a)') 'from collections import OrderedDict'
      write(script_unit,'(a)') 'from mpl_toolkits.axes_grid1 import'&
      //' make_axes_locatable '
      write(script_unit,'(a)') 'import matplotlib.cm as cm'
      write(script_unit,'(a)') ' '
      write(script_unit,'(a)') "points = np.loadtxt('"//trim(seedname)//"_sf-T001.dat')"
      write(script_unit,'(a)') 'points_x=points[:,0]'
      write(script_unit,'(a)') 'points_y=points[:,1]'
      write(script_unit,'(a)') 'val_array=points[:,2]'
      write(script_unit,'(a)') 'num_pt=len(points)'
      write(script_unit,'(a)') ' '

      write(script_unit,'(a)') '#points_x=points_x*??'
      write(script_unit,'(a)') '#points_y=points_y*??'
      write(script_unit,'(a)') ' '
      write(script_unit,'(a)')&
           'x_coord=list(OrderedDict.fromkeys(points_x))'
      write(script_unit,'(a)')&
           'y_coord=list(OrderedDict.fromkeys(points_y))'
      write(script_unit,'(a)') 'dimx=len(x_coord)'
      write(script_unit,'(a)') 'dimy=len(y_coord)'
      write(script_unit,'(a)') ' '
      write(script_unit,'(a)') "outfile ='"//trim(seedname)//"_sf-T001.pdf'"
      write(script_unit,'(a)') ' '
      write(script_unit,'(a)')&
           'vval=val_array.reshape(dimx,dimy).transpose()'
      write(script_unit,'(a)') "#pl.contourf(x_coord,y_coord,"&
       //"vval,1000,origin='lower')"
      write(script_unit,'(a)') 'pl.imshow(vval,origin="lower",'&
             //'extent=(min(x_coord),max(x_coord),min(y_coord),'&
         //'max(y_coord)),interpolation="kaiser")'
      write(script_unit,'(a)') '#pl.imshow(vval,origin="lower",'&
             //'extent=(min(x_coord),max(x_coord),min(y_coord),'&
         //'max(y_coord)),interpolation="kaiser",cmap=cm.hot)'

      write(script_unit,'(a)') ' '
      write(script_unit,'(a)') 'ax = pl.gca()'
      write(script_unit,'(a)') 'divider = make_axes_locatable(ax)'
      write(script_unit,'(a)') 'cax = divider.append_axes("right",'&
       //'"5%", pad="3%")'
      write(script_unit,'(a)') 'pl.colorbar(cax=cax)'
      write(script_unit,'(a)') 'ax.xaxis.set_visible(True)'
      write(script_unit,'(a)') 'ax.yaxis.set_visible(True)'
      write(script_unit,'(a)') 'ax.xaxis.set_major_locator'&
       //'(pl.MultipleLocator(1.0))'
      write(script_unit,'(a)') 'ax.yaxis.set_major_locator'&
       //'(pl.MultipleLocator(1.0))'
      write(script_unit,'(a)') "#ax.set_xlabel('(hh0)',fontsize = 25)"
      write(script_unit,'(a)') "#ax.set_ylabel('(00l)',fontsize = 25)"
      write(script_unit,'(a)') ' '
      write(script_unit,'(a)') 'pl.savefig(outfile)'
      write(script_unit,'(a)') 'pl.show()'
      close(script_unit)

      filename = trim(seedname)//'_sf.gnu'
      open(script_unit,file=filename,status='unknown',form='formatted',action='write')
      write(script_unit,'(a)') '#!/usr/bin/gnuplot'
      write(script_unit,'(a)') 'set term post eps enh color '
      write(script_unit,'(a)') "set output '"//trim(seedname)//"_sf-T001.eps'"
      write(script_unit,'(a)') 'set border lw 1.7'
      write(script_unit,'(a)') 'set size ratio -1'
      write(script_unit,'(a)') '#set xrange [:]'
      write(script_unit,'(a)') '#set yrange [:]'
      write(script_unit,'(a)') "#set xlabel 'h' offset 0,-1,0"
      write(script_unit,'(a)') "#set ylabel 'l'offset -3,0,0"
      write(script_unit,'(a)') '#set xtics font "Times-Roman,30"'
      write(script_unit,'(a)') '#set ytics font "Times-Roman, 30"'
      write(script_unit,'(a)') '#set cbtics font "Times-Roman, 30"'

      write(script_unit,'(a)') 'set xlabel font "Times-Roman-Bold-Italic, 45"'
      write(script_unit,'(a)') 'set ylabel font "Times-Roman-Bold-Italic, 45"'
      write(script_unit,'(a)') "#set ytics ('-1' -1,'0'  0,'1' 1,'2' 2,'3' 3,'4' 4)"
      write(script_unit,'(a)') "#set xtics ('0'  0,'1' 1,'2' 2,'3' 3,'4' 4)"
      write(script_unit,'(a)') "set colorbox"
      write(script_unit,'(a)') "# set cbtics 0.4"
      write(script_unit,'(a)') '#set title "T=0.11"  font "Times-Roman-Bold-Italic, 40"'

      write(script_unit,'(a)') 'set pm3d '
      write(script_unit,'(a)') 'set dgrid3d 100,100 ,3'
      write(script_unit,'(a)') 'set view map'
      write(script_unit,'(a)') 'set palette defined ( 0 "#000090",\'
      write(script_unit,'(a)') '     1 "#000fff",\'
      write(script_unit,'(a)') '     2 "#0090ff",\'
      write(script_unit,'(a)') '     3 "#0fffee",\'
      write(script_unit,'(a)') '     4 "#90ff70",\'
      write(script_unit,'(a)') '     5 "#ffee00",\'
      write(script_unit,'(a)') '     6 "#ff7000",\'
      write(script_unit,'(a)') '     7 "#ee0000",\'
      write(script_unit,'(a)') '     8 "#7f0000")'
      write(script_unit,'(a)') "splot  '"//trim(seedname)//"_sf-T001.dat' u 1:2:3 notitle palette "
      write(script_unit,'(a)') 'unset tics'
      write(script_unit,'(a)') 'unset key'
      close(script_unit)

135    format(3x,E15.8,3x,E15.8,3x,E15.8)
136    format(3x,E15.8,3x,E15.8,3x,E15.8,3x,E15.8,3x,E15.8)

      return

    end subroutine mcarlo_structure_factor_write_mc

    !==================================================================!
    subroutine mcarlo_structure_factor_write_pt(sfactor_matrix)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  
      use mc_io,         only : stdout,io_file_unit,seedname
      use mc_constants,  only : twopi                        
      use mc_utility,    only : utility_recip_lattice,utility_cross 
      use mc_parameters, only : sfactor_corner,sfactor_q1,sfactor_q2,&
                                sfactor_2dqmesh,sfactor_nqpts,&
                                tems,tems_num,lsfactor_polar

      implicit none

      real(kind=dp),intent(in)   :: sfactor_matrix(3,sfactor_nqpts,tems_num) 
      integer                    :: loop_xy,loop_x,loop_y
      integer                    :: script_unit,n,n1,n2,n3
      integer,       allocatable :: sfactor_unit(:)
      real(kind=dp)              :: qvec(3,3)
      real(kind=dp)              :: avec_2d(3,3),areaq1q2,rdum
      real(kind=dp)              :: q1mod,q2mod,ymod,cosq1q2,cosyq2
      real(kind=dp)              :: yvec(3),zvec(3)
      real(kind=dp)              :: qpt_x,qpt_y
      real(kind=dp)              :: q1,q2
      character(len=40)          :: filename

      qvec(1,:) = sfactor_q1
      qvec(2,:) = sfactor_q2

      ! z_vec (orthogonal to the slice)
      call utility_cross(qvec(1,:),qvec(2,:),zvec)
      ! y_vec (orthogonal to q1=x_vec)
      call utility_cross(zvec,qvec(1,:),yvec)
      ! Moduli b1,b2,y_vec
      q1mod = sqrt(qvec(1,1)**2+qvec(1,2)**2+qvec(1,3)**2)
      q2mod = sqrt(qvec(2,1)**2+qvec(2,2)**2+qvec(2,3)**2)
      ymod = sqrt(yvec(1)**2+yvec(2)**2+yvec(3)**2)
      areaq1q2 = sqrt(zvec(1)**2+zvec(2)**2+zvec(3)**2)
      qvec(3,:) = zvec(:)/areaq1q2

      ! Now that we have qvec(3,:), we can compute the dual vectors
      ! avec_2d as in the 3D case
      call utility_recip_lattice(qvec,avec_2d,rdum)

      ! Cosine of the angle between y_vec and b2
      cosyq2 = yvec(1)*qvec(2,1)+yvec(2)*qvec(2,2)+yvec(3)*qvec(2,3)
      cosyq2 = cosyq2/(ymod*q2mod)
      ! Cosine of the angle between b1=x_vec and b2
      cosq1q2 = qvec(1,1)*qvec(2,1)+qvec(1,2)*qvec(2,2)+qvec(1,3)*qvec(2,3)
      cosq1q2 = cosq1q2/(q1mod*q2mod)

      allocate(sfactor_unit(tems_num))
      do n=1,tems_num
         n1 = n/100
         n2 = (n-n1*100)/10
         n3 = n-n1*100-n2*10
         sfactor_unit(n) = io_file_unit()
         filename = trim(seedname)//'_sf-T'&
                    //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat'
         write(stdout,'(3x,a,f8.3,3x,a,3x,a)') 'Output Structure factor at T = ',tems(n),'writted at',filename
         open(sfactor_unit(n),file=filename,status='unknown',form='formatted',action='write')
      enddo


      do n=1,tems_num
         do loop_xy=0,sfactor_nqpts-1
            loop_x = loop_xy/(sfactor_2dqmesh(2) + 1)
            loop_y = loop_xy-loop_x*(sfactor_2dqmesh(2) + 1)
            ! q1 and q2 are the coefficients of the q-point in the basis
            ! (sfactor_q1,sfactor_q2)
            q1 = loop_x/real(sfactor_2dqmesh(1),dp)
            q2 = loop_y/real(sfactor_2dqmesh(2), dp)
            ! Add to (q1,q2) the projection of sfactor_corner on the
            ! (sfactor_q1,sfactor_q2) plane, expressed as a linear
            ! combination of sfactor_q1 and sfactor_q2 
            q1=q1+dot_product(sfactor_corner,avec_2d(1,:))/twopi
            q2=q2+dot_product(sfactor_corner,avec_2d(2,:))/twopi

            qpt_x = q1*q1mod+q2*q2mod*cosq1q2
            qpt_y = q2*q2mod*cosyq2
            if (lsfactor_polar) then
              write(sfactor_unit(n),136) qpt_x,qpt_y,sfactor_matrix(1,loop_xy+1,n),&
                    sfactor_matrix(2,loop_xy+1,n),sfactor_matrix(3,loop_xy+1,n)
            else
              write(sfactor_unit(n),135) qpt_x,qpt_y,sfactor_matrix(1,loop_xy+1,n)
            endif
         enddo
         close(sfactor_unit(n))
      enddo

      filename=trim(seedname)//'_sf.py'
      open(script_unit,file=filename,status='unknown',form='formatted',action='write')
      write(script_unit,'(a)') 'import pylab as pl'
      write(script_unit,'(a)') 'import numpy as np'
      write(script_unit,'(a)') 'import matplotlib.mlab as ml'
      write(script_unit,'(a)') 'from collections import OrderedDict'
      write(script_unit,'(a)') 'from mpl_toolkits.axes_grid1 import'&
      //' make_axes_locatable '
      write(script_unit,'(a)') 'import matplotlib.cm as cm'
      write(script_unit,'(a)') ' '
      write(script_unit,'(a)') "points = np.loadtxt('"//trim(seedname)//"_sf-T001.dat')"
      write(script_unit,'(a)') 'points_x=points[:,0]'
      write(script_unit,'(a)') 'points_y=points[:,1]'
      write(script_unit,'(a)') 'val_array=points[:,2]'
      write(script_unit,'(a)') 'num_pt=len(points)'
      write(script_unit,'(a)') ' '

      write(script_unit,'(a)') '#points_x=points_x*??'
      write(script_unit,'(a)') '#points_y=points_y*??'
      write(script_unit,'(a)') ' '
      write(script_unit,'(a)')&
           'x_coord=list(OrderedDict.fromkeys(points_x))'
      write(script_unit,'(a)')&
           'y_coord=list(OrderedDict.fromkeys(points_y))'
      write(script_unit,'(a)') 'dimx=len(x_coord)'
      write(script_unit,'(a)') 'dimy=len(y_coord)'
      write(script_unit,'(a)') ' '
      write(script_unit,'(a)') "outfile ='"//trim(seedname)//"_sf-T001.pdf'"
      write(script_unit,'(a)') ' '
      write(script_unit,'(a)')&
           'vval=val_array.reshape(dimx,dimy).transpose()'
      write(script_unit,'(a)') "#pl.contourf(x_coord,y_coord,"&
       //"vval,1000,origin='lower')"
      write(script_unit,'(a)') 'pl.imshow(vval,origin="lower",'&
             //'extent=(min(x_coord),max(x_coord),min(y_coord),'&
         //'max(y_coord)),interpolation="kaiser")'
      write(script_unit,'(a)') '#pl.imshow(vval,origin="lower",'&
             //'extent=(min(x_coord),max(x_coord),min(y_coord),'&
         //'max(y_coord)),interpolation="kaiser",cmap=cm.hot)'

      write(script_unit,'(a)') ' '
      write(script_unit,'(a)') 'ax = pl.gca()'
      write(script_unit,'(a)') 'divider = make_axes_locatable(ax)'
      write(script_unit,'(a)') 'cax = divider.append_axes("right",'&
       //'"5%", pad="3%")'
      write(script_unit,'(a)') 'pl.colorbar(cax=cax)'
      write(script_unit,'(a)') 'ax.xaxis.set_visible(True)'
      write(script_unit,'(a)') 'ax.yaxis.set_visible(True)'
      write(script_unit,'(a)') 'ax.xaxis.set_major_locator'&
       //'(pl.MultipleLocator(1.0))'
      write(script_unit,'(a)') 'ax.yaxis.set_major_locator'&
       //'(pl.MultipleLocator(1.0))'
      write(script_unit,'(a)') "#ax.set_xlabel('(hh0)',fontsize = 25)"
      write(script_unit,'(a)') "#ax.set_ylabel('(00l)',fontsize = 25)"
      write(script_unit,'(a)') ' '
      write(script_unit,'(a)') 'pl.savefig(outfile)'
      write(script_unit,'(a)') 'pl.show()'
      close(script_unit)

      filename = trim(seedname)//'_sf.gnu'
      open(script_unit,file=filename,status='unknown',form='formatted',action='write')
      write(script_unit,'(a)') '#!/usr/bin/gnuplot'
      write(script_unit,'(a)') 'set term post eps enh color '
      write(script_unit,'(a)') "set output '"//trim(seedname)//"_sf-T001.eps'"
      write(script_unit,'(a)') 'set border lw 1.7'
      write(script_unit,'(a)') 'set size ratio -1'
      write(script_unit,'(a)') '#set xrange [:]'
      write(script_unit,'(a)') '#set yrange [:]'
      write(script_unit,'(a)') "#set xlabel 'h' offset 0,-1,0"
      write(script_unit,'(a)') "#set ylabel 'l'offset -3,0,0"
      write(script_unit,'(a)') '#set xtics font "Times-Roman,30"'
      write(script_unit,'(a)') '#set ytics font "Times-Roman, 30"'
      write(script_unit,'(a)') '#set cbtics font "Times-Roman, 30"'

      write(script_unit,'(a)') 'set xlabel font "Times-Roman-Bold-Italic, 45"'
      write(script_unit,'(a)') 'set ylabel font "Times-Roman-Bold-Italic, 45"'
      write(script_unit,'(a)') "#set ytics ('-1' -1,'0'  0,'1' 1,'2' 2,'3' 3,'4' 4)"
      write(script_unit,'(a)') "#set xtics ('0'  0,'1' 1,'2' 2,'3' 3,'4' 4)"
      write(script_unit,'(a)') "set colorbox"
      write(script_unit,'(a)') "# set cbtics 0.4"
      write(script_unit,'(a)') '#set title "T=0.11"  font "Times-Roman-Bold-Italic, 40"'

      write(script_unit,'(a)') 'set pm3d '
      write(script_unit,'(a)') 'set dgrid3d 400,400 ,3'
      write(script_unit,'(a)') 'set view map'
      write(script_unit,'(a)') 'set palette defined ( 0 "#000090",\'
      write(script_unit,'(a)') '     1 "#000fff",\'
      write(script_unit,'(a)') '     2 "#0090ff",\'
      write(script_unit,'(a)') '     3 "#0fffee",\'
      write(script_unit,'(a)') '     4 "#90ff70",\'
      write(script_unit,'(a)') '     5 "#ffee00",\'
      write(script_unit,'(a)') '     6 "#ff7000",\'
      write(script_unit,'(a)') '     7 "#ee0000",\'
      write(script_unit,'(a)') '     8 "#7f0000")'
      write(script_unit,'(a)') "splot  '"//trim(seedname)//"_sf-T001.dat' u 1:2:3 notitle palette "
      write(script_unit,'(a)') 'unset tics'
      write(script_unit,'(a)') 'unset key'
      close(script_unit)

135    format(3x,E15.8,3x,E15.8,3x,E15.8)
136    format(3x,E15.8,3x,E15.8,3x,E15.8,3x,E15.8,3x,E15.8)

      return

    end subroutine mcarlo_structure_factor_write_pt

    !==================================================================!
    subroutine open_E_unit(nsteps,e_unit)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  
      use mc_io,         only : io_file_unit,seedname,io_date
      use mc_parameters, only : tems,tems_num

      implicit none
      integer, intent(in)    :: nsteps
      integer, intent(inout) :: e_unit(tems_num)
      integer                :: n,n1,n2,n3
      character(len=40)      :: filename
      character(len=9)       :: cdate, ctime


      do n=1,tems_num
         n1 = n/100
         n2 = (n-n1*100)/10
         n3 = n-n1*100-n2*10
         e_unit(n) = io_file_unit()
         filename = trim(seedname)//'_energy-T'&
                    //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat'
         open(e_unit(n),file=filename,status='unknown',form='formatted',action='write')
         call io_date(cdate,ctime)
         write(e_unit(n),*) '#written on '//cdate//' at '//ctime ! Date and time
         write(e_unit(n),'(3x,a,f9.4)') '# T = ',tems(n)
         write(e_unit(n),'(3x,i16)') nsteps
      enddo

      return

    end subroutine open_E_unit                

    !==================================================================!
    subroutine write_E(nmc,nt,nprint,energy,e_unit)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  
  
      implicit none

      integer,       intent(in) :: nmc
      integer,       intent(in) :: nt,nprint
      real(kind=dp), intent(in) :: energy(0:nprint-1,0:nt-1)               
      integer,    intent(inout) :: e_unit(nt)
      integer                   :: loop,loop_e,n

      do loop=0,nt-1
         n = loop+1
         do loop_e=1,nprint-1
            write(e_unit(n),'(I14,4x,1E16.8)')  nmc-nprint+loop_e,energy(loop_e,loop)/real(num_total_atoms,dp)
         enddo
         write(e_unit(n),'(I14,4x,1E16.8)')  nmc,energy(0,loop)/real(num_total_atoms,dp)
      enddo
 
      return

    end subroutine write_E                

    !==================================================================!
    subroutine close_E_unit(e_unit)
      !==================================================================!
      !                                                                  !
      !                                                                  ! 
      !                                                                  !
      !===================================================================  

      use mc_io,         only : stdout,io_file_unit,seedname
      use mc_parameters, only : tems,tems_num

      implicit none

      integer,    intent(inout) :: e_unit(tems_num)
      integer                   :: n,n1,n2,n3
      character(len=40)         :: filename

      do n=1,tems_num
         n1 = n/100
         n2 = (n-n1*100)/10
         n3 = n-n1*100-n2*10
         e_unit(n) = io_file_unit()
         filename = trim(seedname)//'_energy-T'&
                    //achar(48+n1)//achar(48+n2)//achar(48+n3)//'.dat'
         write(stdout,'(3x,a,f9.4,3x,a,3x,a)') 'Output Energy at T = ',tems(n),'writted at',filename
         close(e_unit(n))
      enddo
 
      return

    end subroutine close_E_unit                

    
    !==================================================================!
    function xrand(init)
      !==================================================================!
      use stdtypes
      use mtprng

      implicit none

      real*8 xrand
      integer init
      logical lfirst
      type(mtprng_state) :: state
      save state
      data lfirst/.true./
      save lfirst

      !  Initialize in first call:
      if (lfirst) then
         call mtprng_init(init,state)
         lfirst = .false.
      endif
      xrand =  mtprng_rand_real2(state)

      return

    end function xrand


end module mc_get_quant      
