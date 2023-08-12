module ModCompare

  implicit none

  include "constants.h" 
  integer, parameter :: i4 = kind(1)
  integer, parameter :: r4 = kind(1.0)
  private
  public :: compare

  interface compare
     module procedure &
          c0dr4, c1pr4, c2pr4, c3pr4, c4pr4, &
          c0di4, c1pi4, c2pi4, c3pi4, c4pi4, &
          c0dc,  c1dc,  c2dc,  c3dc,  c4dc,  &
          c0dl,  c1dl,  c2dl,  c3dl,  c4dl

  end interface

  integer, parameter :: sizepowers=62
  integer(i8),parameter :: powers(sizepowers)=(/ &
       2_i8,  4_i8,  8_i8,  16_i8,  32_i8,  &
       64_i8,  128_i8,  256_i8,  512_i8,  1024_i8, &
       2048_i8,  4096_i8,  8192_i8,  16384_i8,  32768_i8,  &
       65536_i8,  131072_i8,  262144_i8,  524288_i8,  1048576_i8, &
       2097152_i8,  4194304_i8,  8388608_i8,  16777216_i8,  33554432_i8,  &
       67108864_i8,  134217728_i8,  268435456_i8,  536870912_i8,  1073741824_i8, &
       2147483648_i8,  4294967296_i8,  8589934592_i8,  17179869184_i8,  34359738368_i8,  &
       68719476736_i8,  137438953472_i8,  274877906944_i8,  549755813888_i8,  1099511627776_i8, &
       2199023255552_i8,  4398046511104_i8,  8796093022208_i8,  17592186044416_i8,  35184372088832_i8,  &
       70368744177664_i8,  140737488355328_i8,  281474976710656_i8,  562949953421312_i8,  1125899906842624_i8, &
       2251799813685248_i8,  4503599627370496_i8,  9007199254740992_i8,  18014398509481984_i8,  36028797018963968_i8,  &
       72057594037927936_i8,  144115188075855872_i8,  288230376151711744_i8,  576460752303423488_i8,  1152921504606846976_i8, &
       2305843009213693952_i8,  4611686018427387904_i8/)

  integer, parameter :: maxbins=6
  integer, parameter :: lboundbins(maxbins) = (/ &
       0, 10, 20, 30, 40, 50 /)
  integer, parameter :: uboundbins(maxbins) = (/ &
       10, 20, 30, 40, 50, 60/)

contains


  ! compares two distinct r4 values



  subroutine TwoEntriesR4(a1, a2, &
       maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
    real(r4), intent(in) :: a1
    real(r4), intent(in) :: a2
    real(r4), intent(inout) :: maxdif
    real(r4), intent(inout) :: maxar1
    real(r4), intent(inout) :: maxar2
    integer,  intent(inout) :: maxbit
    integer,  intent(in)    :: maxbitsmantissa
    integer,  intent(inout) :: bins(maxbins)
    real(r4) :: difabs, spa, rintervals
    integer :: intervals, difbits, j


    ! how many floating point intervals at this entry

    difabs = abs(a1 - a2)
    spa = min(spacing(abs(a1)),spacing(abs(a2)))
    rintervals = difabs/spa
    intervals = ceiling(rintervals,i8)

    ! how many different bits among entries

    do difbits = 1, min(sizepowers, maxbitsmantissa)
       if (intervals < powers(difbits)) then
          exit
       end if
    end do

    ! put difference in corresponding bin

    do j = 1, maxbins
       if ( difbits> lboundbins(j) .and. &
            difbits<=uboundbins(j)        ) then
          bins(j) = bins(j) + 1
          exit
       end if
    end do

    ! get maximum difference

    if (difabs > maxdif) then
       maxbit = difbits
       maxdif = difabs
       maxar1 = a1
       maxar2 = a2
    end if
  end subroutine TwoEntriesR4



  ! rank independent single precision real output


  subroutine OutputR4(h, msg, verb, cntdif, sizein, zero1, zero2, &
       maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
    character(len=*), intent(in) :: h
    character(len=*), intent(in) :: msg
    logical,          intent(in) :: verb
    integer(i8),      intent(in) :: cntdif
    integer(i8),      intent(in) :: sizein
    logical,          intent(in) :: zero1
    logical,          intent(in) :: zero2
    real,             intent(in) :: maxdif
    real,             intent(in) :: maxar1
    real,             intent(in) :: maxar2
    integer,          intent(in) :: maxbit
    integer,          intent(in) :: maxbitsmantissa
    integer,          intent(in) :: bins(maxbins)
    integer :: i
    character(len=20) :: c0, c1

    if (cntdif == 0_i8) then

       if (verb) then

          ! no differences; verify if any array is null

          if (zero1 .and. zero2) then
             write(*,"(a,' both null')") msg
          else
             write(*,"(a,' matches')") msg
          end if
       end if

    else

       ! there are differences

       write(c0,"(i20)") cntdif
       write(c1,"(i20)") sizein
       write (*,"(a,1x,a,' differences in ',a,' entries; (',i3,'%)')") &
            h//" "//msg//":", trim(adjustl(c0)), trim(adjustl(c1)), (100*cntdif)/sizein

       ! case one array is null

       if (zero1) then
          write (*,"(10x,' first null; max abs second=',e11.3)") &
               maxar2
       else if (zero2) then
          write (*,"(10x,' second null; max abs first=',e11.3)") &
               maxar1
       else

          ! both arrays not null

          write (*,"(2x,' max rel dif: ',i3,' mantissa bits in ',i3)") &
               maxbit, maxbitsmantissa 
          write (*,"(2x,1p,' dif=',e10.3,', spacing=',e10.3,&
               &', entry1=',e10.3,', entry2=',e10.3)")&
               maxdif, spacing(maxdif), maxar1, maxar2
          do i = 1, maxbins
             write (*,"(2x,f6.2,'% differences from ',i2,' to ',i2,' bits')") &
                  100.0*real(bins(i))/real(cntdif), lboundbins(i), uboundbins(i)
             if (maxbit <= uboundbins(i)) then
                exit
             end if
          end do
       end if
    end if
  end subroutine OutputR4



  ! case single precision real scalars



  subroutine c0dr4 (a1, a2, msg, verb)
    real(kind=r4),    intent(in) :: a1
    real(kind=r4),    intent(in) :: a2
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: bins(maxbins)
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: maxbit
    integer :: maxbitsmantissa
    logical :: zero1
    logical :: zero2

    real(kind=r4) :: maxdif
    real(kind=r4) :: maxar1
    real(kind=r4) :: maxar2
    character(len=*), parameter :: h="**(c0dr4)**"



    ! local variables

    sizein = 1_i8
    cntdif = 0_i8       ! how many different entries
    maxbit = 0          ! how many different bits at maximum difference
    bins = 0            ! count differences in bins

    zero1  = a1==0.0_r4 ! first arg is null
    zero2  = a2==0.0_r4 ! second arg is null

    maxdif = 0.0_r4     ! maximum difference among entries of both arrays
    maxar1 = a1         ! value of a1 at entry with maximum difference
    maxar2 = a2         ! value of a2 at entry with maximum difference
    maxbitsmantissa=digits(maxdif)

    ! count differences 

    if (a1 /= a2) then
       cntdif = 1_i8
    else
       cntdif = 0_i8
    end if

    ! there are differences and both scalars not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then
       call TwoEntriesR4(a1, a2, &
            maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
    end if

    ! output

    call OutputR4(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
  end subroutine c0dr4



  ! case single precision real 1D arrays



  subroutine c1dr4 (a1, a2, msg, verb)
    real(kind=r4),    intent(in) :: a1(:)
    real(kind=r4),    intent(in) :: a2(:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: bins(maxbins)
    integer :: d1a1
    integer :: d1a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: i
    integer :: maxbit
    integer :: maxbitsmantissa
    logical :: zero1
    logical :: zero2

    real(kind=r4) :: maxdif
    real(kind=r4) :: maxar1
    real(kind=r4) :: maxar2
    character(len=*), parameter :: h="**(c1dr4)**"

    ! input arrays shape

    d1a1=size(a1); d1a2=size(a2)
    sizein=int(d1a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched sizes when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries
    maxbit = 0          ! how many different bits at maximum difference
    bins = 0            ! count differences in bins

    zero1  = .true.     ! first arg is null
    zero2  = .true.     ! second arg is null

    maxdif = 0.0_r4     ! maximum difference among entries of both arrays
    maxar1 = 0.0_r4     ! value of a1 at entry with maximum difference
    maxar2 = 0.0_r4     ! value of a2 at entry with maximum difference
    maxbitsmantissa=digits(maxdif)

    ! count differences 

    do i = 1, d1a1
       if (a1(i) /= a2(i)) then
          cntdif = cntdif + 1_i8
       end if
       zero1 = zero1 .and. a1(i)==0.0_r4
       zero2 = zero2 .and. a2(i)==0.0_r4
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do i = 1, d1a1
          if (a1(i) /= a2(i)) then
             call TwoEntriesR4(a1(i), a2(i), &
                  maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
          end if
       end do
    end if

    ! output

    call OutputR4(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
  end subroutine c1dr4



  ! case single precision real 2D arrays



  subroutine c2dr4 (a1, a2, msg, verb)
    real(kind=r4),    intent(in) :: a1(:,:)
    real(kind=r4),    intent(in) :: a2(:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: bins(maxbins)
    integer :: d1a1, d2a1
    integer :: d1a2, d2a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2
    integer :: maxbit
    integer :: maxbitsmantissa
    logical :: zero1
    logical :: zero2

    real(kind=r4) :: maxdif
    real(kind=r4) :: maxar1
    real(kind=r4) :: maxar2
    character(len=*), parameter :: h="**(c2dr4)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    sizein=int(d1a1,i8) * int(d2a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries
    maxbit = 0          ! how many different bits at maximum difference
    bins = 0            ! count differences in bins

    zero1  = .true.     ! first arg is null
    zero2  = .true.     ! second arg is null

    maxdif = 0.0_r4     ! maximum difference among entries of both arrays
    maxar1 = 0.0_r4     ! value of a1 at entry with maximum difference
    maxar2 = 0.0_r4     ! value of a2 at entry with maximum difference
    maxbitsmantissa=digits(maxdif)

    ! count differences 

    do ind2 = 1, d2a1
       do ind1 = 1, d1a1
          if (a1(ind1,ind2) /= a2(ind1,ind2)) then
             cntdif = cntdif + 1_i8
          end if
          zero1 = zero1 .and. a1(ind1,ind2)==0.0_r4
          zero2 = zero2 .and. a2(ind1,ind2)==0.0_r4
       end do
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do ind2 = 1, d2a1
          do ind1 = 1, d1a1
             if (a1(ind1,ind2) /= a2(ind1,ind2)) then
                call TwoEntriesR4(a1(ind1,ind2), a2(ind1,ind2), &
                     maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
             end if
          end do
       end do
    end if

    ! output

    call OutputR4(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
  end subroutine c2dr4



  ! case single precision real 3D arrays



  subroutine c3dr4 (a1, a2, msg, verb)
    real(kind=r4),    intent(in) :: a1(:,:,:)
    real(kind=r4),    intent(in) :: a2(:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: bins(maxbins)
    integer :: d1a1, d2a1, d3a1
    integer :: d1a2, d2a2, d3a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2, ind3
    integer :: maxbit
    integer :: maxbitsmantissa
    logical :: zero1
    logical :: zero2

    real(kind=r4) :: maxdif
    real(kind=r4) :: maxar1
    real(kind=r4) :: maxar2
    integer :: cntDifInd1
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(c3dr4)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    d3a1=size(a1,3); d3a2=size(a2,3)
    sizein=int(d1a1,i8) * int(d2a1,i8) * int(d3a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    else if (d3a1 /= d3a2) then
       call fatal_error(h//' unmatched third  dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries
    maxbit = 0          ! how many different bits at maximum difference
    bins = 0            ! count differences in bins

    zero1  = .true.     ! first arg is null
    zero2  = .true.     ! second arg is null

    maxdif = 0.0_r4     ! maximum difference among entries of both arrays
    maxar1 = 0.0_r4     ! value of a1 at entry with maximum difference
    maxar2 = 0.0_r4     ! value of a2 at entry with maximum difference
    maxbitsmantissa=digits(maxdif)

    ! count differences 

    do ind3 = 1, d3a1
       do ind2 = 1, d2a1
          do ind1 = 1, d1a1
             if (a1(ind1,ind2,ind3) /= a2(ind1,ind2,ind3)) then
                cntdif = cntdif + 1_i8
             end if
             zero1 = zero1 .and. a1(ind1,ind2,ind3)==0.0_r4
             zero2 = zero2 .and. a2(ind1,ind2,ind3)==0.0_r4
          end do
       end do
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       write(str(1),"(i8)") d1a1
       write(str(2),"(i8)") d2a1
       write(str(3),"(i8)") d3a1
       write(*,"(a)") h//" for field "//trim(adjustl(msg))//&
            " localy declared ("//&
            "1:"//trim(adjustl(str(1)))//","//&
            "1:"//trim(adjustl(str(2)))//","//&
            "1:"//trim(adjustl(str(3)))//"):"

       do ind3 = 1, d3a1
          do ind2 = 1, d2a1
             cntDifInd1=0
             do ind1 = 1, d1a1
                if (a1(ind1,ind2,ind3) /= a2(ind1,ind2,ind3)) then
                   cntDifInd1=cntDifInd1+1
                   call TwoEntriesR4(a1(ind1,ind2,ind3), a2(ind1,ind2,ind3), &
                        maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
                end if
             end do
!!$             if (cntDifInd1 /= 0) then
!!$                write(str(1),"(i8)") cntDifInd1
!!$                write(str(2),"(i8)") ind2
!!$                write(str(3),"(i8)") ind3
!!$                write(*,"(a)")h//" there are "//trim(adjustl(str(1)))//&
!!$                     " differences at (:,"//&
!!$                     trim(adjustl(str(2)))//","//&
!!$                     trim(adjustl(str(3)))//")"
!!$             end if
          end do
       end do
    end if

    ! output

    call OutputR4(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
  end subroutine c3dr4



  ! case single precision real 4D arrays



  subroutine c4dr4 (a1, a2, msg, verb)
    real(kind=r4),    intent(in) :: a1(:,:,:,:)
    real(kind=r4),    intent(in) :: a2(:,:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: bins(maxbins)
    integer :: d1a1, d2a1, d3a1, d4a1
    integer :: d1a2, d2a2, d3a2, d4a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2, ind3, ind4
    integer :: maxbit
    integer :: maxbitsmantissa
    logical :: zero1
    logical :: zero2

    real(kind=r4) :: maxdif
    real(kind=r4) :: maxar1
    real(kind=r4) :: maxar2
    character(len=*), parameter :: h="**(c4dr4)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    d3a1=size(a1,3); d3a2=size(a2,3)
    d4a1=size(a1,4); d4a2=size(a2,4)
    sizein=int(d1a1,i8) * int(d2a1,i8) * int(d3a1,i8) * int(d4a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    else if (d3a1 /= d3a2) then
       call fatal_error(h//' unmatched third  dimension when comparing '//msg)
    else if (d4a1 /= d4a2) then
       call fatal_error(h//' unmatched forth  dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries
    maxbit = 0          ! how many different bits at maximum difference
    bins = 0            ! count differences in bins

    zero1  = .true.     ! first arg is null
    zero2  = .true.     ! second arg is null

    maxdif = 0.0_r4     ! maximum difference among entries of both arrays
    maxar1 = 0.0_r4     ! value of a1 at entry with maximum difference
    maxar2 = 0.0_r4     ! value of a2 at entry with maximum difference
    maxbitsmantissa=digits(maxdif)

    ! count differences 

    do ind4 = 1, d4a1
       do ind3 = 1, d3a1
          do ind2 = 1, d2a1
             do ind1 = 1, d1a1
                if (a1(ind1,ind2,ind3,ind4) /= a2(ind1,ind2,ind3,ind4)) then
                   cntdif = cntdif + 1_i8
                end if
                zero1 = zero1 .and. a1(ind1,ind2,ind3,ind4)==0.0_r4
                zero2 = zero2 .and. a2(ind1,ind2,ind3,ind4)==0.0_r4
             end do
          end do
       end do
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do ind4 = 1, d4a1
          do ind3 = 1, d3a1
             do ind2 = 1, d2a1
                do ind1 = 1, d1a1
                   if (a1(ind1,ind2,ind3,ind4) /= a2(ind1,ind2,ind3,ind4)) then
                      call TwoEntriesR4(a1(ind1,ind2,ind3,ind4), a2(ind1,ind2,ind3,ind4), &
                           maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
                   end if
                end do
             end do
          end do
       end do
    end if

    ! output

    call OutputR4(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
  end subroutine c4dr4



  ! compares two distinct I4 values



  subroutine TwoEntriesI4(a1, a2, maxdif, maxar1, maxar2)
    integer, intent(in) :: a1
    integer, intent(in) :: a2
    integer, intent(inout) :: maxdif
    integer, intent(inout) :: maxar1
    integer, intent(inout) :: maxar2
    integer :: difabs

    difabs = abs(a1 - a2)

    ! get maximum difference

    if (difabs > maxdif) then
       maxdif = difabs
       maxar1 = a1
       maxar2 = a2
    end if
  end subroutine TwoEntriesI4



  ! rank independent single precision integer output


  subroutine OutputI4(h, msg, verb, cntdif, sizein, zero1, zero2, &
       maxdif, maxar1, maxar2)
    character(len=*), intent(in) :: h
    character(len=*), intent(in) :: msg
    logical,          intent(in) :: verb
    integer(i8),      intent(in) :: cntdif
    integer(i8),      intent(in) :: sizein
    logical,          intent(in) :: zero1
    logical,          intent(in) :: zero2
    integer,          intent(in) :: maxdif
    integer,          intent(in) :: maxar1
    integer,          intent(in) :: maxar2
    character(len=20) :: c0, c1

    if (cntdif == 0_i8) then

       if (verb) then

          ! no differences; verify if any array is null

          if (zero1 .and. zero2) then
             write(*,"(a,' both null')") msg
          else
             write(*,"(a,' matches')") msg
          end if
       end if

    else

       ! there are differences

       write(c0,"(i20)") cntdif
       write(c1,"(i20)") sizein
       write (*,"(a,1x,a,' differences in ',a,' entries; (',i3,'%)')") &
            h//" "//msg//":", trim(adjustl(c0)), trim(adjustl(c1)), (100*cntdif)/sizein

       ! case one array is null

       if (zero1) then
          write (*,"(10x,' first null; max abs second=',i10)") &
               maxar2
       else if (zero2) then
          write (*,"(10x,' second null; max abs first=',i10)") &
               maxar1
       else
          write (*,"(2x,1p,' dif=',i10, &
               &', entry1=',i10,', entry2=',i10)")&
               maxdif, maxar1, maxar2
       end if
    end if
  end subroutine OutputI4



  ! case single precision integer scalars



  subroutine c0di4 (a1, a2, msg, verb)
    integer,          intent(in) :: a1
    integer,          intent(in) :: a2
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer(i8) :: sizein
    integer(i8) :: cntdif
    logical :: zero1
    logical :: zero2

    integer :: maxdif
    integer :: maxar1
    integer :: maxar2
    character(len=*), parameter :: h="**(c0di4)**"

    ! input arrays shape

    sizein=1

    ! local variables

    zero1  = a1==0      ! first arg is null
    zero2  = a2==0      ! second arg is null

    maxdif = abs(a1-a2) ! maximum difference among entries of both arrays
    maxar1 = a1    ! value of a1 at entry with maximum difference
    maxar2 = a2    ! value of a2 at entry with maximum difference

    ! count differences 

    if (a1 /= a2) then
       cntdif = 1_i8
    else
       cntdif = 0_i8
    end if

    ! output

    call OutputI4(h, msg, verb, cntdif, 1_i8, zero1, zero2, &
         maxdif, maxar1, maxar2)
  end subroutine c0di4



  ! case single precision integer 1D arrays



  subroutine c1di4 (a1, a2, msg, verb)
    integer,          intent(in) :: a1(:)
    integer,          intent(in) :: a2(:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1
    integer :: d1a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1
    logical :: zero1
    logical :: zero2

    integer :: maxdif
    integer :: maxar1
    integer :: maxar2
    character(len=*), parameter :: h="**(c1di4)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    sizein=int(d1a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is null
    zero2  = .true.     ! second arg is null

    maxdif = 0     ! maximum difference among entries of both arrays
    maxar1 = 0     ! value of a1 at entry with maximum difference
    maxar2 = 0     ! value of a2 at entry with maximum difference

    ! count differences 

    do ind1 = 1, d1a1
       if (a1(ind1) /= a2(ind1)) then
          cntdif = cntdif + 1_i8
       end if
       zero1 = zero1 .and. a1(ind1)==0.0_r4
       zero2 = zero2 .and. a2(ind1)==0.0_r4
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do ind1 = 1, d1a1
          if (a1(ind1) /= a2(ind1)) then
             call TwoEntriesI4(a1(ind1), a2(ind1), &
                  maxdif, maxar1, maxar2)
          end if
       end do
    end if

    ! output

    call OutputI4(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2)
  end subroutine c1di4



  ! case single precision integer 2D arrays



  subroutine c2di4 (a1, a2, msg, verb)
    integer,          intent(in) :: a1(:,:)
    integer,          intent(in) :: a2(:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1, d2a1
    integer :: d1a2, d2a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2
    logical :: zero1
    logical :: zero2

    integer :: maxdif
    integer :: maxar1
    integer :: maxar2
    character(len=*), parameter :: h="**(c2di4)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    sizein=int(d1a1,i8) * int(d2a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is null
    zero2  = .true.     ! second arg is null

    maxdif = 0     ! maximum difference among entries of both arrays
    maxar1 = 0     ! value of a1 at entry with maximum difference
    maxar2 = 0     ! value of a2 at entry with maximum difference

    ! count differences 

    do ind2 = 1, d2a1
       do ind1 = 1, d1a1
          if (a1(ind1,ind2) /= a2(ind1,ind2)) then
             cntdif = cntdif + 1_i8
          end if
          zero1 = zero1 .and. a1(ind1,ind2)==0.0_r4
          zero2 = zero2 .and. a2(ind1,ind2)==0.0_r4
       end do
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do ind2 = 1, d2a1
          do ind1 = 1, d1a1
             if (a1(ind1,ind2) /= a2(ind1,ind2)) then
                call TwoEntriesI4(a1(ind1,ind2), a2(ind1,ind2), &
                     maxdif, maxar1, maxar2)
             end if
          end do
       end do
    end if

    ! output

    call OutputI4(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2)
  end subroutine c2di4



  ! case single precision integer 3D arrays



  subroutine c3di4 (a1, a2, msg, verb)
    integer,          intent(in) :: a1(:,:,:)
    integer,          intent(in) :: a2(:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1, d2a1, d3a1
    integer :: d1a2, d2a2, d3a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2, ind3
    logical :: zero1
    logical :: zero2

    integer :: maxdif
    integer :: maxar1
    integer :: maxar2
    character(len=*), parameter :: h="**(c3di4)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    d3a1=size(a1,3); d3a2=size(a2,3)
    sizein=int(d1a1,i8) * int(d2a1,i8) * int(d3a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    else if (d3a1 /= d3a2) then
       call fatal_error(h//' unmatched third  dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is null
    zero2  = .true.     ! second arg is null

    maxdif = 0     ! maximum difference among entries of both arrays
    maxar1 = 0     ! value of a1 at entry with maximum difference
    maxar2 = 0     ! value of a2 at entry with maximum difference

    ! count differences 

    do ind3 = 1, d3a1
       do ind2 = 1, d2a1
          do ind1 = 1, d1a1
             if (a1(ind1,ind2,ind3) /= a2(ind1,ind2,ind3)) then
                cntdif = cntdif + 1_i8
             end if
             zero1 = zero1 .and. a1(ind1,ind2,ind3)==0.0_r4
             zero2 = zero2 .and. a2(ind1,ind2,ind3)==0.0_r4
          end do
       end do
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do ind3 = 1, d3a1
          do ind2 = 1, d2a1
             do ind1 = 1, d1a1
                if (a1(ind1,ind2,ind3) /= a2(ind1,ind2,ind3)) then
                   call TwoEntriesI4(a1(ind1,ind2,ind3), a2(ind1,ind2,ind3), &
                        maxdif, maxar1, maxar2)
                end if
             end do
          end do
       end do
    end if

    ! output

    call OutputI4(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2)
  end subroutine c3di4



  ! case single precision integer 4D arrays



  subroutine c4di4 (a1, a2, msg, verb)
    integer,          intent(in) :: a1(:,:,:,:)
    integer,          intent(in) :: a2(:,:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1, d2a1, d3a1, d4a1
    integer :: d1a2, d2a2, d3a2, d4a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2, ind3, ind4
    logical :: zero1
    logical :: zero2

    integer :: maxdif
    integer :: maxar1
    integer :: maxar2
    character(len=*), parameter :: h="**(c4di4)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    d3a1=size(a1,3); d3a2=size(a2,3)
    d4a1=size(a1,4); d4a2=size(a2,4)
    sizein=int(d1a1,i8) * int(d2a1,i8) * int(d3a1,i8) * int(d4a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    else if (d3a1 /= d3a2) then
       call fatal_error(h//' unmatched third  dimension when comparing '//msg)
    else if (d4a1 /= d4a2) then
       call fatal_error(h//' unmatched forth  dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is null
    zero2  = .true.     ! second arg is null

    maxdif = 0     ! maximum difference among entries of both arrays
    maxar1 = 0     ! value of a1 at entry with maximum difference
    maxar2 = 0     ! value of a2 at entry with maximum difference

    ! count differences 

    do ind4 = 1, d4a1
       do ind3 = 1, d3a1
          do ind2 = 1, d2a1
             do ind1 = 1, d1a1
                if (a1(ind1,ind2,ind3,ind4) /= a2(ind1,ind2,ind3,ind4)) then
                   cntdif = cntdif + 1_i8
                end if
                zero1 = zero1 .and. a1(ind1,ind2,ind3,ind4)==0.0_r4
                zero2 = zero2 .and. a2(ind1,ind2,ind3,ind4)==0.0_r4
             end do
          end do
       end do
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do ind4 = 1, d4a1
          do ind3 = 1, d3a1
             do ind2 = 1, d2a1
                do ind1 = 1, d1a1
                   if (a1(ind1,ind2,ind3,ind4) /= a2(ind1,ind2,ind3,ind4)) then
                      call TwoEntriesI4(a1(ind1,ind2,ind3,ind4), a2(ind1,ind2,ind3,ind4), &
                           maxdif, maxar1, maxar2)
                   end if
                end do
             end do
          end do
       end do
    end if

    ! output

    call OutputI4(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2)
  end subroutine c4di4



  ! compares two distinct character values



  subroutine TwoEntriesC(a1, a2, maxdif, maxar1, maxar2)
    character(len=*), intent(in) :: a1
    character(len=*), intent(in) :: a2
    integer, intent(inout) :: maxdif
    character(len=*), intent(inout) :: maxar1
    character(len=*), intent(inout) :: maxar2
    integer :: difabs

    difabs = abs(len_trim(a1) - len_trim(a2))

    ! get maximum difference

    if (difabs >= maxdif) then
       maxdif = difabs
       maxar1 = a1
       maxar2 = a2
    end if
  end subroutine TwoEntriesC



  ! rank independent character output


  subroutine OutputC(h, msg, verb, cntdif, sizein, zero1, zero2, &
       maxdif, maxar1, maxar2)
    character(len=*), intent(in) :: h
    character(len=*), intent(in) :: msg
    logical,          intent(in) :: verb
    integer(i8),      intent(in) :: cntdif
    integer(i8),      intent(in) :: sizein
    logical,          intent(in) :: zero1
    logical,          intent(in) :: zero2
    integer,          intent(in) :: maxdif
    character(len=*), intent(in) :: maxar1
    character(len=*), intent(in) :: maxar2
    character(len=20) :: c0, c1

    if (cntdif == 0_i8) then

       if (verb) then

          ! no differences; verify if any array is empty

          if (zero1 .and. zero2) then
             write(*,"(a,' both empty')") msg
          else
             write(*,"(a,' matches')") msg
          end if
       end if

    else

       ! there are differences

       write(c0,"(i20)") cntdif
       write(c1,"(i20)") sizein
       write (*,"(a,1x,a,' differences in ',a,' entries; (',i3,'%)')") &
            h//" "//msg//":", trim(adjustl(c0)), trim(adjustl(c1)), (100*cntdif)/sizein

       ! case one array is null

       if (zero1) then
          write (*,"(10x,' first empty, second not empty - one entry is **',a,'**')") &
               trim(adjustl(maxar2))
       else if (zero2) then
          write (*,"(10x,' second empty; first not empty - one entry is **',a,'**')") &
               trim(adjustl(maxar1))
       else
          write (*,"(2x,1p,' max length dif=',i10,&
               &', entry1=**',a,'**, entry2=**',a,'**')")&
               maxdif, trim(adjustl(maxar1)), trim(adjustl(maxar2))
       end if
    end if
  end subroutine OutputC



  ! case scalar character



  subroutine c0dc (a1, a2, msg, verb)
    character(len=*), intent(in) :: a1
    character(len=*), intent(in) :: a2
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer(i8) :: sizein
    integer(i8) :: cntdif
    logical :: zero1
    logical :: zero2

    integer :: maxdif
    character(len=len(a1)) :: maxar1
    character(len=len(a2)) :: maxar2
    character(len=*), parameter :: h="**(c0dc)**"

    ! input arrays shape

    sizein=1_i8

    if ( trim(adjustl(a1)) /= &
         trim(adjustl(a2))     ) then
       cntdif = 1_i8
    else
       cntdif = 0_i8
    end if

    zero1  = len_trim(a1)==0 ! first arg is empty
    zero2  = len_trim(a2)==0 ! second arg is empty

    maxdif = 0     ! maximum difference among entries of both arrays
    maxar1 = " "   ! value of a1 at entry with maximum difference
    maxar2 = " "   ! value of a2 at entry with maximum difference

    ! count differences 


    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       if ( trim(adjustl(a1)) /= &
            trim(adjustl(a2))     ) then
          call TwoEntriesC(a1, a2, &
               maxdif, maxar1, maxar2)
       end if
    end if

    ! output

    call OutputC(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2)
  end subroutine c0dc



  ! case 1D character arrays



  subroutine c1dc (a1, a2, msg, verb)
    character(len=*), intent(in) :: a1(:)
    character(len=*), intent(in) :: a2(:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1
    integer :: d1a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1
    logical :: zero1
    logical :: zero2

    integer :: maxdif
    character(len=len(a1)) :: maxar1
    character(len=len(a2)) :: maxar2
    character(len=*), parameter :: h="**(c1dc)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    sizein=int(d1a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is empty
    zero2  = .true.     ! second arg is empty

    maxdif = 0     ! maximum difference among entries of both arrays
    maxar1 = " "   ! value of a1 at entry with maximum difference
    maxar2 = " "   ! value of a2 at entry with maximum difference

    ! count differences 

    do ind1 = 1, d1a1
       if ( trim(adjustl(a1(ind1))) /= &
            trim(adjustl(a2(ind1)))     ) then
          cntdif = cntdif + 1_i8
       end if
       zero1 = zero1 .and. len_trim(a1(ind1))==0
       zero2 = zero2 .and. len_trim(a2(ind1))==0
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do ind1 = 1, d1a1
          if ( trim(adjustl(a1(ind1))) /= &
               trim(adjustl(a2(ind1)))     ) then
             call TwoEntriesC(a1(ind1), a2(ind1), &
                  maxdif, maxar1, maxar2)
          end if
       end do
    end if

    ! output

    call OutputC(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2)
  end subroutine c1dc



  ! case 2D character arrays



  subroutine c2dc (a1, a2, msg, verb)
    character(len=*), intent(in) :: a1(:,:)
    character(len=*), intent(in) :: a2(:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1, d2a1
    integer :: d1a2, d2a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2
    logical :: zero1
    logical :: zero2

    integer :: maxdif
    character(len=len(a1)) :: maxar1
    character(len=len(a2)) :: maxar2
    character(len=*), parameter :: h="**(c2dc)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    sizein=int(d1a1,i8) * int(d2a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is empty
    zero2  = .true.     ! second arg is empty

    maxdif = 0     ! maximum difference among entries of both arrays
    maxar1 = " "   ! value of a1 at entry with maximum difference
    maxar2 = " "   ! value of a2 at entry with maximum difference

    ! count differences 

    do ind2 = 1, d2a1
       do ind1 = 1, d1a1
          if ( trim(adjustl(a1(ind1,ind2))) /= &
               trim(adjustl(a2(ind1,ind2)))     ) then
             cntdif = cntdif + 1_i8
          end if
          zero1 = zero1 .and. len_trim(a1(ind1,ind2))==0
          zero2 = zero2 .and. len_trim(a2(ind1,ind2))==0
       end do
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do ind2 = 1, d2a1
          do ind1 = 1, d1a1
             if ( trim(adjustl(a1(ind1,ind2))) /= &
                  trim(adjustl(a2(ind1,ind2)))     ) then
                call TwoEntriesC(a1(ind1,ind2), a2(ind1,ind2), &
                     maxdif, maxar1, maxar2)
             end if
          end do
       end do
    end if

    ! output

    call OutputC(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2)
  end subroutine c2dc



  ! case 3D character arrays



  subroutine c3dc (a1, a2, msg, verb)
    character(len=*), intent(in) :: a1(:,:,:)
    character(len=*), intent(in) :: a2(:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1, d2a1, d3a1
    integer :: d1a2, d2a2, d3a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2, ind3
    logical :: zero1
    logical :: zero2

    integer :: maxdif
    character(len=len(a1)) :: maxar1
    character(len=len(a2)) :: maxar2
    character(len=*), parameter :: h="**(c3dc)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    d3a1=size(a1,3); d3a2=size(a2,3)
    sizein=int(d1a1,i8) * int(d2a1,i8) * int(d3a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    else if (d3a1 /= d3a2) then
       call fatal_error(h//' unmatched third  dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is empty
    zero2  = .true.     ! second arg is empty

    maxdif = 0     ! maximum difference among entries of both arrays
    maxar1 = " "   ! value of a1 at entry with maximum difference
    maxar2 = " "   ! value of a2 at entry with maximum difference

    ! count differences 

    do ind3 = 1, d3a1
       do ind2 = 1, d2a1
          do ind1 = 1, d1a1
             if ( trim(adjustl(a1(ind1,ind2,ind3))) /= &
                  trim(adjustl(a2(ind1,ind2,ind3)))     ) then
                cntdif = cntdif + 1_i8
             end if
             zero1 = zero1 .and. len_trim(a1(ind1,ind2,ind3))==0
             zero2 = zero2 .and. len_trim(a2(ind1,ind2,ind3))==0
          end do
       end do
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do ind3 = 1, d3a1
          do ind2 = 1, d2a1
             do ind1 = 1, d1a1
                if ( trim(adjustl(a1(ind1,ind2,ind3))) /= &
                     trim(adjustl(a2(ind1,ind2,ind3)))     ) then
                   call TwoEntriesC(a1(ind1,ind2,ind3), a2(ind1,ind2,ind3), &
                        maxdif, maxar1, maxar2)
                end if
             end do
          end do
       end do
    end if

    ! output

    call OutputC(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2)
  end subroutine c3dc



  ! case 4D character arrays



  subroutine c4dc (a1, a2, msg, verb)
    character(len=*), intent(in) :: a1(:,:,:,:)
    character(len=*), intent(in) :: a2(:,:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1, d2a1, d3a1, d4a1
    integer :: d1a2, d2a2, d3a2, d4a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2, ind3, ind4
    logical :: zero1
    logical :: zero2

    integer :: maxdif
    character(len=len(a1)) :: maxar1
    character(len=len(a2)) :: maxar2
    character(len=*), parameter :: h="**(c4dc)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    d3a1=size(a1,3); d3a2=size(a2,3)
    d4a1=size(a1,4); d4a2=size(a2,4)
    sizein=int(d1a1,i8) * int(d2a1,i8) * int(d3a1,i8) * int(d4a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    else if (d3a1 /= d3a2) then
       call fatal_error(h//' unmatched third  dimension when comparing '//msg)
    else if (d4a1 /= d4a2) then
       call fatal_error(h//' unmatched forth  dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is empty
    zero2  = .true.     ! second arg is empty

    maxdif = 0     ! maximum difference among entries of both arrays
    maxar1 = " "   ! value of a1 at entry with maximum difference
    maxar2 = " "   ! value of a2 at entry with maximum difference

    ! count differences 

    do ind4 = 1, d4a1
       do ind3 = 1, d3a1
          do ind2 = 1, d2a1
             do ind1 = 1, d1a1
                if ( trim(adjustl(a1(ind1,ind2,ind3,ind4))) /= &
                     trim(adjustl(a2(ind1,ind2,ind3,ind4)))     ) then
                   cntdif = cntdif + 1_i8
                end if
                zero1 = zero1 .and. len_trim(a1(ind1,ind2,ind3,ind4))==0
                zero2 = zero2 .and. len_trim(a2(ind1,ind2,ind3,ind4))==0
             end do
          end do
       end do
    end do

    ! there are differences and both arrays not null

    if (cntdif /= 0_i8 .and. .not. (zero1 .and. zero2)) then

       ! for all distinct entries:

       do ind4 = 1, d4a1
          do ind3 = 1, d3a1
             do ind2 = 1, d2a1
                do ind1 = 1, d1a1
                   if ( trim(adjustl(a1(ind1,ind2,ind3,ind4))) /= &
                        trim(adjustl(a2(ind1,ind2,ind3,ind4)))     ) then
                      call TwoEntriesC(a1(ind1,ind2,ind3,ind4), a2(ind1,ind2,ind3,ind4), &
                           maxdif, maxar1, maxar2)
                   end if
                end do
             end do
          end do
       end do
    end if

    ! output

    call OutputC(h, msg, verb, cntdif, sizein, zero1, zero2, &
         maxdif, maxar1, maxar2)
  end subroutine c4dc



  ! rank independent character output


  subroutine OutputL(h, msg, verb, cntdif, sizein, zero1, zero2)
    character(len=*), intent(in) :: h
    character(len=*), intent(in) :: msg
    logical,          intent(in) :: verb
    integer(i8),      intent(in) :: cntdif
    integer(i8),      intent(in) :: sizein
    logical,          intent(in) :: zero1
    logical,          intent(in) :: zero2
    character(len=20) :: c0, c1

    if (cntdif == 0_i8) then

       if (verb) then

          ! no differences; verify if any array is empty

          if (zero1 .and. zero2) then
             write(*,"(a,' both false')") msg
          else
             write(*,"(a,' matches')") msg
          end if
       end if

    else

       ! there are differences

       write(c0,"(i20)") cntdif
       write(c1,"(i20)") sizein
       write (*,"(a,1x,a,' differences in ',a,' entries; (',i3,'%)')") &
            h//" "//msg//":", trim(adjustl(c0)), trim(adjustl(c1)), (100*cntdif)/sizein

       ! case one array is null

       if (zero1) then
          write (*,"(10x,' first false, second not false')")
       else if (zero2) then
          write (*,"(10x,' second false, first not false')") 
       end if
    end if
  end subroutine OutputL



  ! case logical scalars



  subroutine c0dl (a1, a2, msg, verb)
    logical,          intent(in) :: a1
    logical,          intent(in) :: a2
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer(i8) :: sizein
    integer(i8) :: cntdif
    logical :: zero1
    logical :: zero2
    character(len=*), parameter :: h="**(c0dl)**"

    ! input arrays shape

    sizein=1_i8

    ! local variables

    zero1  = .not. a1   ! first arg is false
    zero2  = .not. a2   ! second arg is false

    ! count differences 

    if (a1 .NEQV. a2) then
       cntdif = 1_i8
    else
       cntdif = 0_i8
    end if

    ! output

    call OutputL(h, msg, verb, cntdif, sizein, zero1, zero2)
  end subroutine c0dl



  ! case 1D logical arrays



  subroutine c1dl (a1, a2, msg, verb)
    logical,          intent(in) :: a1(:)
    logical,          intent(in) :: a2(:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1
    integer :: d1a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1
    logical :: zero1
    logical :: zero2
    character(len=*), parameter :: h="**(c1dl)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    sizein=int(d1a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is false
    zero2  = .true.     ! second arg is false

    ! count differences 

    do ind1 = 1, d1a1
       if ( a1(ind1) .NEQV. &
            a2(ind1)     ) then
          cntdif = cntdif + 1_i8
       end if
       zero1 = zero1 .and. .not. a1(ind1)
       zero2 = zero2 .and. .not. a2(ind1)
    end do

    ! output

    call OutputL(h, msg, verb, cntdif, sizein, zero1, zero2)
  end subroutine c1dl



  ! case 2D logical arrays



  subroutine c2dl (a1, a2, msg, verb)
    logical,          intent(in) :: a1(:,:)
    logical,          intent(in) :: a2(:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1, d2a1
    integer :: d1a2, d2a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2
    logical :: zero1
    logical :: zero2
    character(len=*), parameter :: h="**(c2dl)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    sizein=int(d1a1,i8) * int(d2a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is false
    zero2  = .true.     ! second arg is false

    ! count differences 

    do ind2 = 1, d2a1
       do ind1 = 1, d1a1
          if ( a1(ind1,ind2) .NEQV. &
               a2(ind1,ind2)     ) then
             cntdif = cntdif + 1_i8
          end if
          zero1 = zero1 .and. .not. a1(ind1,ind2)
          zero2 = zero2 .and. .not. a2(ind1,ind2)
       end do
    end do

    ! output

    call OutputL(h, msg, verb, cntdif, sizein, zero1, zero2)
  end subroutine c2dl



  ! case 3D logical arrays



  subroutine c3dl (a1, a2, msg, verb)
    logical,          intent(in) :: a1(:,:,:)
    logical,          intent(in) :: a2(:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1, d2a1, d3a1
    integer :: d1a2, d2a2, d3a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2, ind3
    logical :: zero1
    logical :: zero2
    character(len=*), parameter :: h="**(c3dl)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    d3a1=size(a1,3); d3a2=size(a2,3)
    sizein=int(d1a1,i8) * int(d2a1,i8) * int(d3a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    else if (d3a1 /= d3a2) then
       call fatal_error(h//' unmatched third  dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is false
    zero2  = .true.     ! second arg is false

    ! count differences 

    do ind3 = 1, d3a1
       do ind2 = 1, d2a1
          do ind1 = 1, d1a1
             if ( a1(ind1,ind2,ind3) .NEQV. &
                  a2(ind1,ind2,ind3)     ) then
                cntdif = cntdif + 1_i8
             end if
             zero1 = zero1 .and. .not. a1(ind1,ind2,ind3)
             zero2 = zero2 .and. .not. a2(ind1,ind2,ind3)
          end do
       end do
    end do

    ! output

    call OutputL(h, msg, verb, cntdif, sizein, zero1, zero2)
  end subroutine c3dl



  ! case 4D logical arrays



  subroutine c4dl (a1, a2, msg, verb)
    logical,          intent(in) :: a1(:,:,:,:)
    logical,          intent(in) :: a2(:,:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    integer :: d1a1, d2a1, d3a1, d4a1
    integer :: d1a2, d2a2, d3a2, d4a2
    integer(i8) :: sizein
    integer(i8) :: cntdif
    integer :: ind1, ind2, ind3, ind4
    logical :: zero1
    logical :: zero2
    character(len=*), parameter :: h="**(c4dl)**"

    ! input arrays shape

    d1a1=size(a1,1); d1a2=size(a2,1)
    d2a1=size(a1,2); d2a2=size(a2,2)
    d3a1=size(a1,3); d3a2=size(a2,3)
    d4a1=size(a1,4); d4a2=size(a2,4)
    sizein=int(d1a1,i8) * int(d2a1,i8) * int(d3a1,i8) * int(d4a1,i8)
    if (d1a1 /= d1a2) then
       call fatal_error(h//' unmatched first dimension when comparing '//msg)
    else if (d2a1 /= d2a2) then
       call fatal_error(h//' unmatched second dimension when comparing '//msg)
    else if (d3a1 /= d3a2) then
       call fatal_error(h//' unmatched third  dimension when comparing '//msg)
    else if (d4a1 /= d4a2) then
       call fatal_error(h//' unmatched forth  dimension when comparing '//msg)
    end if

    ! local variables

    cntdif = 0_i8       ! how many different entries

    zero1  = .true.     ! first arg is false
    zero2  = .true.     ! second arg is false

    ! count differences 

    do ind4 = 1, d4a1
       do ind3 = 1, d3a1
          do ind2 = 1, d2a1
             do ind1 = 1, d1a1
                if ( a1(ind1,ind2,ind3,ind4) .NEQV. &
                     a2(ind1,ind2,ind3,ind4)     ) then
                   cntdif = cntdif + 1_i8
                end if
                zero1 = zero1 .and. .not. a1(ind1,ind2,ind3,ind4)
                zero2 = zero2 .and. .not. a2(ind1,ind2,ind3,ind4)
             end do
          end do
       end do
    end do

    ! output

    call OutputL(h, msg, verb, cntdif, sizein, zero1, zero2)
  end subroutine c4dl



  ! case single precision real 1D pointer arrays



  subroutine c1pr4 (a1, a2, msg, verb)
    real(kind=r4),    pointer, intent(in) :: a1(:)
    real(kind=r4),    pointer, intent(in) :: a2(:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    logical :: assoc1
    logical :: assoc2
    character(len=*), parameter :: h="**(c1pr4)**"

    ! both pointers have the same association status

    assoc1=associated(a1)
    assoc2=associated(a2)

    if (assoc1 .and. assoc2) then
       call c1dr4 (a1, a2, msg, verb)
    else if ( &
         ((.not. assoc1) .and. assoc2) .or. &
         ((.not. assoc2) .and. assoc1)) then
       call fatal_error(h//" only one pointer is associated when comparing "//msg)
    else
       if (verb) then
          write(*,"(a)") h//" both pointers are dissasociated when comparing //msg"
       end if
    end if
  end subroutine c1pr4



  ! case single precision real 2D pointer arrays



  subroutine c2pr4 (a1, a2, msg, verb)
    real(kind=r4),    pointer, intent(in) :: a1(:,:)
    real(kind=r4),    pointer, intent(in) :: a2(:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    logical :: assoc1
    logical :: assoc2
    character(len=*), parameter :: h="**(c2pr4)**"

    ! both pointers have the same association status

    assoc1=associated(a1)
    assoc2=associated(a2)

    if (assoc1 .and. assoc2) then
       call c2dr4 (a1, a2, msg, verb)
    else if ( &
         ((.not. assoc1) .and. assoc2) .or. &
         ((.not. assoc2) .and. assoc1)) then
       call fatal_error(h//" only one pointer is associated when comparing "//msg)
    else
       if (verb) then
          write(*,"(a)") h//" both pointers are dissasociated when comparing //msg"
       end if
    end if
  end subroutine c2pr4



  ! case single precision real 3D pointer arrays



  subroutine c3pr4 (a1, a2, msg, verb)
    real(kind=r4),    pointer, intent(in) :: a1(:,:,:)
    real(kind=r4),    pointer, intent(in) :: a2(:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    logical :: assoc1
    logical :: assoc2
    character(len=*), parameter :: h="**(c3pr4)**"

    ! both pointers have the same association status

    assoc1=associated(a1)
    assoc2=associated(a2)

    if (assoc1 .and. assoc2) then
       call c3dr4 (a1, a2, msg, verb)
    else if ( &
         ((.not. assoc1) .and. assoc2) .or. &
         ((.not. assoc2) .and. assoc1)) then
       call fatal_error(h//" only one pointer is associated when comparing "//msg)
    else
       if (verb) then
          write(*,"(a)") h//" both pointers are dissasociated when comparing //msg"
       end if
    end if
  end subroutine c3pr4



  ! case single precision real 4D pointer arrays



  subroutine c4pr4 (a1, a2, msg, verb)
    real(kind=r4),    pointer, intent(in) :: a1(:,:,:,:)
    real(kind=r4),    pointer, intent(in) :: a2(:,:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    logical :: assoc1
    logical :: assoc2
    character(len=*), parameter :: h="**(c4pr4)**"

    ! both pointers have the same association status

    assoc1=associated(a1)
    assoc2=associated(a2)

    if (assoc1 .and. assoc2) then
       call c4dr4 (a1, a2, msg, verb)
    else if ( &
         ((.not. assoc1) .and. assoc2) .or. &
         ((.not. assoc2) .and. assoc1)) then
       call fatal_error(h//" only one pointer is associated when comparing "//msg)
    else
       if (verb) then
          write(*,"(a)") h//" both pointers are dissasociated when comparing //msg"
       end if
    end if
  end subroutine c4pr4



  ! case single precision integer 1D pointer arrays



  subroutine c1pi4 (a1, a2, msg, verb)
    integer(kind=i4),    pointer, intent(in) :: a1(:)
    integer(kind=i4),    pointer, intent(in) :: a2(:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    logical :: assoc1
    logical :: assoc2
    character(len=*), parameter :: h="**(c1pi4)**"

    ! both pointers have the same association status

    assoc1=associated(a1)
    assoc2=associated(a2)

    if (assoc1 .and. assoc2) then
       call c1di4 (a1, a2, msg, verb)
    else if ( &
         ((.not. assoc1) .and. assoc2) .or. &
         ((.not. assoc2) .and. assoc1)) then
       call fatal_error(h//" only one pointer is associated when comparing "//msg)
    else
       if (verb) then
          write(*,"(a)") h//" both pointers are dissasociated when comparing //msg"
       end if
    end if
  end subroutine c1pi4



  ! case single precision integer 2D pointer arrays



  subroutine c2pi4 (a1, a2, msg, verb)
    integer(kind=i4),    pointer, intent(in) :: a1(:,:)
    integer(kind=i4),    pointer, intent(in) :: a2(:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    logical :: assoc1
    logical :: assoc2
    character(len=*), parameter :: h="**(c2pi4)**"

    ! both pointers have the same association status

    assoc1=associated(a1)
    assoc2=associated(a2)

    if (assoc1 .and. assoc2) then
       call c2di4 (a1, a2, msg, verb)
    else if ( &
         ((.not. assoc1) .and. assoc2) .or. &
         ((.not. assoc2) .and. assoc1)) then
       call fatal_error(h//" only one pointer is associated when comparing "//msg)
    else
       if (verb) then
          write(*,"(a)") h//" both pointers are dissasociated when comparing //msg"
       end if
    end if
  end subroutine c2pi4



  ! case single precision integer 3D pointer arrays



  subroutine c3pi4 (a1, a2, msg, verb)
    integer(kind=i4),    pointer, intent(in) :: a1(:,:,:)
    integer(kind=i4),    pointer, intent(in) :: a2(:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    logical :: assoc1
    logical :: assoc2
    character(len=*), parameter :: h="**(c3pi4)**"

    ! both pointers have the same association status

    assoc1=associated(a1)
    assoc2=associated(a2)

    if (assoc1 .and. assoc2) then
       call c3di4 (a1, a2, msg, verb)
    else if ( &
         ((.not. assoc1) .and. assoc2) .or. &
         ((.not. assoc2) .and. assoc1)) then
       call fatal_error(h//" only one pointer is associated when comparing "//msg)
    else
       if (verb) then
          write(*,"(a)") h//" both pointers are dissasociated when comparing //msg"
       end if
    end if
  end subroutine c3pi4



  ! case single precision integer 4D pointer arrays



  subroutine c4pi4 (a1, a2, msg, verb)
    integer(kind=i4),    pointer, intent(in) :: a1(:,:,:,:)
    integer(kind=i4),    pointer, intent(in) :: a2(:,:,:,:)
    character(len=*), intent(in) :: msg   ! output header
    logical,          intent(in) :: verb  ! verboses output

    logical :: assoc1
    logical :: assoc2
    character(len=*), parameter :: h="**(c4pi4)**"

    ! both pointers have the same association status

    assoc1=associated(a1)
    assoc2=associated(a2)

    if (assoc1 .and. assoc2) then
       call c4di4 (a1, a2, msg, verb)
    else if ( &
         ((.not. assoc1) .and. assoc2) .or. &
         ((.not. assoc2) .and. assoc1)) then
       call fatal_error(h//" only one pointer is associated when comparing "//msg)
    else
       if (verb) then
          write(*,"(a)") h//" both pointers are dissasociated when comparing //msg"
       end if
    end if
  end subroutine c4pi4
end module ModCompare
