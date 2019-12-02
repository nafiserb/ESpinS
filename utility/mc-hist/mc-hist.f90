!######################################################################
! This program is part of
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

program histogram

  implicit none

  integer, parameter         :: dp=selected_real_kind(15,300)
  integer, parameter         :: k=1500         ! for integration scheme
  real(dp), parameter        :: pi=3.141592653589793238462643383279_dp
  real(dp)                   :: d_omega,t1,t2
  real(dp)                   :: smearing_cutoff,degauss,arg
  integer                    :: mp_type
  real(dp)                   :: energy_min,energy_max,rdum,e
  integer                    :: num_energy,iene,n,num_freq,ifreq,min_e,max_e,loop_e
  logical                    :: smearing
  character(len=30)          :: filename,dummy
  real(kind=dp), allocatable :: energy(:),energy_array(:),hist(:)
  integer                    :: iostat

  call cpu_time(t1)
  call get_command_argument(1,filename)
 
  open(unit=1,file=adjustl(trim(filename)),action='read',iostat=iostat)
  if (iostat /= 0) then
     write(0,'(a)') "Error: Couldn't find file '" // trim(filename) // "' file."
     stop 
  end if
  open(unit=3,file=adjustl(trim(filename))//".hist",action='write')
  write(3,'(a)') '# Running mc-hist. Written by N Rezaei. (c) 2018.'
  write(3,'(a,a,a)') '# Reading input file ', adjustl(trim(filename)), ' ...'

  
  ! read input parameters
  read(1,*) dummy ! 
  read(1,*) dummy ! 
  read(1,*) num_energy  ! number of energies
  write(3,'(a,i12)') '# Number of Energies : ', num_energy
  write(3,'(a,i12)') '#      Energy               Histogram  '
  
  allocate(energy(num_energy))
  do iene=1,num_energy
     read(1,*) n,energy(iene)                                  
  enddo
  close (1)
  ! end read input parameters

  d_omega = 0.010_dp
  energy_min = minval(energy)
  energy_max = maxval(energy)
  num_freq = nint((energy_max-energy_min)/d_omega)+1
  if(num_freq==1) num_freq = 2
  d_omega = (energy_max-energy_min)/(num_freq-1)

  allocate(energy_array(num_freq))
  allocate(hist(num_freq))

  do ifreq=1,num_freq
     energy_array(ifreq) = energy_min+real(ifreq-1,dp)*d_omega
  enddo

  write(*,*) 'Enter gaussian broadening = '
  read (*,*) degauss

  hist = 0.0_dp
!  e_width=energy_array(2)-energy_array(1)
  mp_type = 0
  smearing_cutoff = 10.0_dp
     do iene=1,num_energy
       ! Faster optimization: I precalculate the indices
       if (degauss .le.  0.0_dp) then
          min_e = max(nint((energy(iene) - energy_array(1))/&
                  (energy_array(size(energy_array))-energy_array(1)) &
                  * real(size(energy_array)-1,kind=dp)) + 1, 1)
          max_e = min(nint((energy(iene) - energy_array(1))/&
                  (energy_array(size(energy_array))-energy_array(1)) &
                  * real(size(energy_array)-1,kind=dp)) + 1, size(energy_array))
          smearing = .false.
       else
          min_e = max(nint((energy(iene) - smearing_cutoff * degauss - energy_array(1))/&
                  (energy_array(size(energy_array))-energy_array(1)) &
                  * real(size(energy_array)-1,kind=dp)) + 1, 1)
          max_e = min(nint((energy(iene) + smearing_cutoff * degauss - energy_array(1))/&
                  (energy_array(size(energy_array))-energy_array(1)) &
                  * real(size(energy_array)-1,kind=dp)) + 1, size(energy_array))
          smearing = .true.
       end if

!min_e=1
!max_e=num_freq   

     do loop_e=min_e,max_e
!print*,loop_e
        if (smearing) then
           arg = (energy_array(loop_e)-energy(iene))/degauss
           rdum = methfessel(arg,mp_type)/degauss
        else
           rdum = 1._dp/(energy_array(2)-energy_array(1))
        end if

        hist(loop_e) = hist(loop_e)+rdum 
     enddo

   enddo

   do ifreq=1,num_freq
      write(3,'(2x,1E16.8,5x,1E16.8)'),energy_min+real(ifreq-1,dp)*d_omega,hist(ifreq) 
   enddo

   deallocate(energy)
   deallocate(energy_array)
   deallocate(hist)
   close(3)

    call cpu_time(t2)
print*,'time :',t2-t1,'(s)'

contains

   !============================================!
   !function methfessel
     !============================================!
     function methfessel(x,n)
      implicit none
       real(kind=dp),intent(in) :: x
       integer,intent(in)       ::   n
       integer                  :: i,ni
       real(kind=dp)            :: methfessel,arg
       real(kind=dp)            :: hd,hp,a
    
       arg = min (200.d0, x**2)
       methfessel = exp ( - arg) / sqrt (pi)
       if(n.eq.0) return
       hd = 0.0d0
       hp = exp ( - arg)
       ni = 0
       a = 1.0_dp / sqrt (pi)
       do i = 1, n
          hd = 2.0d0 * x * hp - 2.0d0 * DBLE (ni) * hd
          ni = ni + 1
          a = - a / (DBLE (i) * 4.0d0)
          hp = 2.0d0 * x * hd-2.0d0 * DBLE (ni) * hp
          ni = ni + 1
          methfessel = (methfessel+ a * hp)
       enddo

       return
     end function methfessel
   !===============================================

  
end program histogram
