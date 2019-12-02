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

module mc_constants

  implicit none

  private

  integer,          parameter, public :: dp      = selected_real_kind(15,300)
  real(kind=dp),    parameter, public :: pi      = 3.141592653589793238462643383279_dp
  real(kind=dp),    parameter, public :: twopi   = 2*pi
  real(kind=dp),    parameter, public :: zero    = 0.0_dp 
  complex(kind=dp), parameter, public :: cmplx_0 = (0.0_dp,0.0_dp)
  complex(kind=dp), parameter, public :: cmplx_1 = (1.0_dp,0.0_dp)
  complex(kind=dp), parameter, public :: cmplx_i = (0.0_dp,1.0_dp)

  real(kind=dp),    parameter, public :: bohr    = 0.5291772108_dp
                    
  real(kind=dp),    parameter, public :: eps2    = 1.0e-2_dp
  real(kind=dp),    parameter, public :: eps3    = 1.0e-3_dp
  real(kind=dp),    parameter, public :: eps4    = 1.0e-4_dp
  real(kind=dp),    parameter, public :: eps5    = 1.0e-5_dp
  real(kind=dp),    parameter, public :: eps6    = 1.0e-6_dp
  real(kind=dp),    parameter, public :: eps7    = 1.0e-7_dp
  real(kind=dp),    parameter, public :: eps8    = 1.0e-8_dp
  real(kind=dp),    parameter, public :: eps10   = 1.0e-10_dp
                    
  real(kind=dp),    parameter, public :: k_B       =  8.6173303e-5_dp     ! eV/K
  real(kind=dp),    parameter, public :: hart      = 27.21138602_dp       ! eV
  real(kind=dp),    parameter, public :: bohr_magn =  5.7883818012e-5_dp  ! eV/T

end module mc_constants
