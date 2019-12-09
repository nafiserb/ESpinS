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

! some parts of this module were taken from Wannier90 and modified according to our purposes.

module mc_utility

  use mc_constants, only : dp

  implicit none

  private

  public :: utility_recip_lattice
  public :: utility_cart_to_frac
  public :: utility_frac_to_cart
  public :: utility_lowercase
  public :: utility_strip
  public :: utility_cross          
  public :: utility_symmetric     
  public :: utility_asymmetric_pairs     
  public :: utility_sort     

contains

  !===================================================================
  subroutine utility_recip_lattice (real_lat,recip_lat,volume)  !
    !==================================================================!
    !                                                                  !
    !  Calculates the reciprical lattice vectors and the cell volume   !
    !                                                                  !
    !===================================================================

    use mc_constants,  only : dp,twopi,eps5
    use mc_io,         only : io_error

    implicit none
    real(kind=dp), intent(in)  :: real_lat (3, 3)
    real(kind=dp), intent(out) :: recip_lat (3, 3)  
    real(kind=dp), intent(out) :: volume

    recip_lat(1,1) = real_lat(2,2)*real_lat(3,3)-real_lat(3,2)*real_lat(2,3)
    recip_lat(1,2) = real_lat(2,3)*real_lat(3,1)-real_lat(3,3)*real_lat(2,1)
    recip_lat(1,3) = real_lat(2,1)*real_lat(3,2)-real_lat(3,1)*real_lat(2,2)
    recip_lat(2,1) = real_lat(3,2)*real_lat(1,3)-real_lat(1,2)*real_lat(3,3)
    recip_lat(2,2) = real_lat(3,3)*real_lat(1,1)-real_lat(1,3)*real_lat(3,1)
    recip_lat(2,3) = real_lat(3,1)*real_lat(1,2)-real_lat(1,1)*real_lat(3,2)
    recip_lat(3,1) = real_lat(1,2)*real_lat(2,3)-real_lat(2,2)*real_lat(1,3)
    recip_lat(3,2) = real_lat(1,3)*real_lat(2,1)-real_lat(2,3)*real_lat(1,1)
    recip_lat(3,3) = real_lat(1,1)*real_lat(2,2)-real_lat(2,1)*real_lat(1,2)

    volume = real_lat(1,1)*recip_lat(1,1) + &
             real_lat(1,2)*recip_lat(1,2) + &
             real_lat(1,3)*recip_lat(1,3)  


    if (abs(volume) < eps5 ) then
       call io_error(' Found almost zero Volume in utility_recip_lattice')
    end if

    recip_lat = twopi*recip_lat/volume
    volume = abs(volume)

    return

  end subroutine utility_recip_lattice


  !===================================================================
  subroutine utility_frac_to_cart(frac,cart,real_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================  
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3,3)
    real(kind=dp), intent(in)  :: frac(3)
    real(kind=dp), intent(out) :: cart(3)

    integer :: i

    do i=1,3
       cart(i) = real_lat(1,i)*frac(1) + real_lat(2,i)*frac(2) + real_lat(3,i)*frac(3) 
    end do

    return

  end subroutine utility_frac_to_cart


  !===================================================================
  subroutine utility_cart_to_frac(cart,frac,recip_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================  
    use mc_constants, only : twopi
    implicit none

    real(kind=dp), intent(in)  :: recip_lat(3,3)
    real(kind=dp), intent(out) :: frac(3)
    real(kind=dp), intent(in)  :: cart(3)

    integer :: i

    do i=1,3
       frac(i) = recip_lat(i,1)*cart(1) + recip_lat(i,2)*cart(2) + recip_lat(i,3)*cart(3) 
    end do

    frac = frac/twopi


    return

  end subroutine utility_cart_to_frac

  !=============================!
  function utility_strip(string)!
    !=============================!
    !                             !
    !    Strips string of all     !
    !        blank spaces         !
    !                             !
    !=============================!

    use mc_io, only : maxlen

    implicit none

    character(len=*), intent(in) :: string
    character(len=maxlen) :: utility_strip

    integer :: ispc,ipos,ilett,icount

    ! Initialise
    utility_strip=repeat(' ',maxlen)

    ispc = ichar(' ')
    icount = 0
    do ipos=1,len(string)
       ilett = ichar(string(ipos:ipos))
       if (ilett.ne.ispc) then
          icount = icount+1
          utility_strip(icount:icount) = string(ipos:ipos)
       endif
    enddo

    utility_strip = trim(utility_strip)

    return

  end function utility_strip


  !=================================!
  function utility_lowercase(string)!
    !=================================!
    !                                 !
    ! Takes a string and converts to  !
    !      lowercase characters       !
    !                                 !
    !=================================!

    use mc_io, only : maxlen

    implicit none

    character(len=*), intent(in) :: string
    character(len=maxlen) :: utility_lowercase

    integer :: iA,iZ,idiff,ipos,ilett

    iA = ichar('A')
    iZ = ichar('Z')
    idiff = iZ-ichar('z')

    utility_lowercase = string

    do ipos=1,len(string)
       ilett = ichar(string(ipos:ipos))
       if ((ilett.ge.iA).and.(ilett.le.iZ)) &
            utility_lowercase(ipos:ipos) = char(ilett-idiff)
    enddo

    utility_lowercase = trim(adjustl(utility_lowercase))

    return

  end function utility_lowercase

  !=============================================================!
  subroutine utility_cross(a,b,c)
    !=============================================================!
    !                                                             !
    ! Return cross  product of two vector a and b                 !
    !                                                             !
    !                       C = A X B                             !
    !                                                             !
    !=============================================================!


    implicit none

    real(kind=dp),  intent(in)  :: a(3)
    real(kind=dp),  intent(in)  :: b(3)
    real(kind=dp),  intent(out) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

    return

  end subroutine utility_cross

  !=============================================================!
  subroutine utility_symmetric(n,a,counter) 
    !=============================================================!
    !                                                             !
    !                                                             !
    !=============================================================!


    implicit none

    integer, intent(in)    :: n
    real(kind=dp),  intent(in)  :: a(n,n)
    integer,        intent(out) :: counter
    integer                     :: i,j            


    counter=0
    do i=1,n
       do j=i+1,n
         if (a(i,j) .ne. a(j,i)) then
            counter=counter+1
         endif
       enddo
    enddo 

    return

  end subroutine utility_symmetric

  !=============================================================!
  subroutine utility_asymmetric_pairs(n,counter,a,asym) 
    !=============================================================!
    !                                                             !
    !                                                             !
    !=============================================================!


    implicit none

    integer,        intent(in)  :: n,counter
    real(kind=dp),  intent(in)  :: a(n,n)
    integer,        intent(out) :: asym(n,counter)
    integer                     :: i,j,counter1            


    asym=0
    do i=1,n
       counter1=0
       do j=i+1,n
         if (a(i,j) .ne. a(j,i))   then
            counter1=counter1+1
            asym(i,counter1)=j
         endif
       enddo
    enddo 

    return

  end subroutine utility_asymmetric_pairs

  !==================================================================!
  subroutine utility_sort(n,r_array,i_array,conv)
    !==================================================================!
    !                                                                  !
    ! Sorting an posistive array and giving the convert matrix         !
    ! Doing the search in this order gives a dramatic speed up         !
    !                                                                  !
    !==================================================================!  
    implicit none

    integer,                 intent(in)    :: n
    real(kind=dp), optional, intent(inout) :: r_array(n)
    integer,       optional, intent(inout) :: i_array(n)
    integer,       optional, intent(out)   :: conv(n,n)
    integer                                :: loop,indx(1)
    real(kind=dp)                          :: a(n),a_cp(n)

    if(present(i_array)) a = i_array*1.0_dp
    if(present(r_array)) a = r_array
    if(present(conv)) conv = 0

    do loop=n,1,-1
       indx = internal_maxloc(a,n)
       a_cp(loop) = a(indx(1))
       if(present(conv)) conv(loop,indx(1)) = 1
       a(indx(1)) = -1.0_dp
    end do

    if(present(i_array)) i_array = int(a_cp)
    if(present(r_array)) r_array = a_cp

  end subroutine utility_sort

  !=========================================================================!
  function internal_maxloc(a,dim)
    !=========================================================================!
    !                                                                         !
    !  A predictable maxloc.                                                  !
    !                                                                         !
    !=========================================================================!

    use mc_constants, only : eps8
    implicit none


    integer,       intent(in) :: dim
    real(kind=dp), intent(in) :: a(dim)
    integer                   :: internal_maxloc
    integer                   :: guess(1),loop,counter
    integer                   :: list(dim)

    list=0
    counter=1

    guess=maxloc(a)
    list(1)=guess(1)
    ! look for any degenerate values
    do loop=1,dim
       if (loop==guess(1)) cycle
       if ( abs(a(loop)-a(guess(1))) < eps8 ) then
          counter=counter+1
          list(counter)=loop
       endif
    end do
    ! and always return the lowest index
    internal_maxloc=minval(list(1:counter))

  end function internal_maxloc

end module mc_utility
