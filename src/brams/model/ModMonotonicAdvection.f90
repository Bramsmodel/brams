!**(JP)** Inicio Insercao
module ModUtils
  ! ModUtils: 
  !    Contains:
  !       (1) procedures to allocate, deallocate and dump multi-dimensional pointer arrays
  !       (2) procedures to find and read variables named __<name> at XXX-head.txt files
  !       (3) kind names parameters
  !       (4) error dealing procedure 
  ! Exports:
  !    Subroutine Alloc
  !       allocate memory to pointer array, generating fatal error 
  !       if already allocated. Generic interface for real 1,2,3 and 4D arrays,
  !       as well as 2D integer array.
  !    Subroutine Dealloc
  !       deallocate memory allocated to pointer array: no-op if already
  !       deallocated. Generic interface for real 1,2,3 and 4D arrays,
  !       as well as 2D integer array.
  !    Subroutine Dump
  !       dumps at stdout maximum, minimum and average value of pointer array, 
  !       for allocated pointer arrays. if array not allocated, dumps message.
  !       Generic interface for real 2,3 and 4D arrays, as well as 2D integer array.

  use iso_fortran_env, only: output_unit

  use ModParallelEnvironment, only: MsgDump
  implicit none

  private
  public :: i4, i8, r4, r8

  integer, parameter :: r4 = selected_real_kind(6)  ! kind for 32-bits real numbers
  integer, parameter :: i4 = selected_int_kind(9)   ! kind for 32-bits integer numbers
  integer, parameter :: r8 = selected_real_kind(15) ! kind for 64-bits real numbers
  integer, parameter :: i8 = selected_int_kind(14)  ! kind for 64-bits integer numbers

  interface Dump
     module procedure Dump1D_r, Dump2D_r, Dump3D_r, Dump4D_r, Dump1D_i, Dump2D_i
  end interface Dump

contains



  ! Dump1D_r: prints at stdout maximum, minimum and average values of 
  !           real pointer array ptr(:), if allocated. Prints message if
  !           ptr not allocated. Message contains caller procedure name 
  !           (caller) and ptr field name (fieldName)



  subroutine Dump1D_r(ptr, fieldName, caller)
    real, pointer                :: ptr(:)    ! intent(in)
    character(len=*), intent(in) :: fieldName
    character(len=*), intent(in) :: caller
    character(len=8) :: c0
    character(len=16) :: c3
    character(len=16) :: c4
    character(len=16) :: c5

    if (.not. associated(ptr)) then
       write(*,"(a)") trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//" not associated"
    else
       write(c0,"(i8)") size(ptr,1)
       write(c3,"(e16.7)") minval(ptr)
       write(c4,"(e16.7)") maxval(ptr)
       write(c5,"(e16.7)") sum(ptr)/size(ptr)
       call MsgDump (trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//"("//&
            trim(adjustl(c0))//") with "//&
            " minval="//trim(adjustl(c3))//&
            " maxval="//trim(adjustl(c4))//&
            " aveval="//trim(adjustl(c5)))
    end if
  end subroutine Dump1D_r



  ! Dump2D_r: prints at stdout maximum, minimum and average values of 
  !           real pointer array ptr(:,:), if allocated. Prints message if
  !           ptr not allocated. Message contains caller procedure name 
  !           (caller) and ptr field name (fieldName)



  subroutine Dump2D_r(ptr, fieldName, caller)
    real, pointer                :: ptr(:,:)  ! intent(in)
    character(len=*), intent(in) :: fieldName
    character(len=*), intent(in) :: caller
    character(len=8) :: c0
    character(len=8) :: c1
    character(len=16) :: c3
    character(len=16) :: c4
    character(len=16) :: c5

    if (.not. associated(ptr)) then
       write(*,"(a)") trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//" not associated"
    else
       write(c0,"(i8)") size(ptr,1)
       write(c1,"(i8)") size(ptr,2)
       write(c3,"(e16.7)") minval(ptr)
       write(c4,"(e16.7)") maxval(ptr)
       write(c5,"(e16.7)") sum(ptr)/size(ptr)
       call MsgDump(trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//"("//&
            trim(adjustl(c0))//","//&
            trim(adjustl(c1))//") with "//&
            " minval="//trim(adjustl(c3))//&
            " maxval="//trim(adjustl(c4))//&
            " aveval="//trim(adjustl(c5)))
    end if
  end subroutine Dump2D_r



  ! Dump3D_r: prints at stdout maximum, minimum and average values of 
  !           real pointer array ptr(:,:,:), if allocated. Prints message if
  !           ptr not allocated. Message contains caller procedure name 
  !           (caller) and ptr field name (fieldName)



  subroutine Dump3D_r(ptr, fieldName, caller)
    real, pointer                :: ptr(:,:,:) ! intent(in)
    character(len=*), intent(in) :: fieldName
    character(len=*), intent(in) :: caller
    character(len=8) :: c0
    character(len=8) :: c1
    character(len=8) :: c2
    character(len=16) :: c3
    character(len=16) :: c4
    character(len=16) :: c5

    if (.not. associated(ptr)) then
       write(*,"(a)") trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//" not associated"
    else
       write(c0,"(i8)") size(ptr,1)
       write(c1,"(i8)") size(ptr,2)
       write(c2,"(i8)") size(ptr,3)
       write(c3,"(e16.7)") minval(ptr)
       write(c4,"(e16.7)") maxval(ptr)
       write(c5,"(e16.7)") sum(ptr)/size(ptr)
       call MsgDump(trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//"("//&
            trim(adjustl(c0))//","//&
            trim(adjustl(c1))//","//&
            trim(adjustl(c2))//") with "//&
            " minval="//trim(adjustl(c3))//&
            " maxval="//trim(adjustl(c4))//&
            " aveval="//trim(adjustl(c5)))
    end if
  end subroutine Dump3D_r



  ! Dump4D_r: prints at stdout maximum, minimum and average values of 
  !           real pointer array ptr(:,:,:,:), if allocated. Prints message if
  !           ptr not allocated. Message contains caller procedure name 
  !           (caller) and ptr field name (fieldName)



  subroutine Dump4D_r(ptr, fieldName, caller)
    real, pointer                :: ptr(:,:,:,:) ! intent(in)
    character(len=*), intent(in) :: fieldName
    character(len=*), intent(in) :: caller
    character(len=8) :: c0
    character(len=8) :: c1
    character(len=8) :: c2
    character(len=8) :: c6
    character(len=16) :: c3
    character(len=16) :: c4
    character(len=16) :: c5

    if (.not. associated(ptr)) then
       write(*,"(a)") trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//" not associated"
    else
       write(c0,"(i8)") size(ptr,1)
       write(c1,"(i8)") size(ptr,2)
       write(c2,"(i8)") size(ptr,3)
       write(c6,"(i8)") size(ptr,4)
       write(c3,"(e16.7)") minval(ptr)
       write(c4,"(e16.7)") maxval(ptr)
       write(c5,"(e16.7)") sum(ptr)/size(ptr)
       call MsgDump(trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//"("//&
            trim(adjustl(c0))//","//&
            trim(adjustl(c1))//","//&
            trim(adjustl(c2))//","//&
            trim(adjustl(c6))//") with "//&
            " minval="//trim(adjustl(c3))//&
            " maxval="//trim(adjustl(c4))//&
            " aveval="//trim(adjustl(c5)))
    end if
  end subroutine Dump4D_r



  ! Dump1D_i: prints at stdout maximum, minimum and average values of 
  !           integer pointer array ptr(:), if allocated. Prints message if
  !           ptr not allocated. Message contains caller procedure name 
  !           (caller) and ptr field name (fieldName)



  subroutine Dump1D_i(ptr, fieldName, caller)
    integer, pointer             :: ptr(:)     ! intent(in)
    character(len=*), intent(in) :: fieldName
    character(len=*), intent(in) :: caller
    character(len=8) :: c0
    character(len=16) :: c3
    character(len=16) :: c4
    character(len=16) :: c5

    if (.not. associated(ptr)) then
       write(*,"(a)") trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//" not associated"
    else
       write(c0,"(i8)") size(ptr,1)
       write(c3,"(e16.7)") minval(ptr)
       write(c4,"(e16.7)") maxval(ptr)
       write(c5,"(e16.7)") sum(ptr)/size(ptr)
       call MsgDump(trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//"("//&
            trim(adjustl(c0))//") with "//&
            " minval="//trim(adjustl(c3))//&
            " maxval="//trim(adjustl(c4))//&
            " aveval="//trim(adjustl(c5)))
    end if
  end subroutine Dump1D_i



  ! Dump2D_i: prints at stdout maximum, minimum and average values of 
  !           integer pointer array ptr(:,:), if allocated. Prints message if
  !           ptr not allocated. Message contains caller procedure name 
  !           (caller) and ptr field name (fieldName)



  subroutine Dump2D_i(ptr, fieldName, caller)
    integer, pointer             :: ptr(:,:)   ! intent(in)
    character(len=*), intent(in) :: fieldName
    character(len=*), intent(in) :: caller
    character(len=8) :: c0
    character(len=8) :: c1
    character(len=16) :: c3
    character(len=16) :: c4
    character(len=16) :: c5

    if (.not. associated(ptr)) then
       write(*,"(a)") trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//" not associated"
    else
       write(c0,"(i8)") size(ptr,1)
       write(c1,"(i8)") size(ptr,2)
       write(c3,"(e16.7)") minval(ptr)
       write(c4,"(e16.7)") maxval(ptr)
       write(c5,"(e16.7)") sum(ptr)/size(ptr)
       call MsgDump(trim(adjustl(caller))//" "//&
            trim(adjustl(fieldName))//"("//&
            trim(adjustl(c0))//","//&
            trim(adjustl(c1))//") with "//&
            " minval="//trim(adjustl(c3))//&
            " maxval="//trim(adjustl(c4))//&
            " aveval="//trim(adjustl(c5)))
    end if
  end subroutine Dump2D_i
end module ModUtils
module ModCompare

  use ModUtils, only: i4
  use ModUtils, only: i8
  use ModUtils, only: r4
  use ModUtils, only: r8
  use ModParallelEnvironment, only: MsgDump
  implicit none

  private
  public :: compare

  interface compare
     module procedure &
          c0dr4, c1pr4, c2pr4, c3pr4, c4pr4, &
          c0di4, c1pi4, c2pi4, c3pi4, c4pi4, &
          c0dc,  c1dc,  c2dc,  c3dc,  c4dc,  &
          c0dl,  c1dl,  c2dl,  c3dl,  c4dl

  end interface compare

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
    character(len=12) :: c0, c1, c2, c3

    if (cntdif == 0_i8) then

       if (verb) then

          ! no differences; verify if any array is null

          if (zero1 .and. zero2) then
             call MsgDump(trim(msg)//' both null')
          else
             call MsgDump(trim(msg)//' matches')
          end if
       end if

    else

       ! there are differences

       write(c0,"(i8)") cntdif
       write(c1,"(i8)") sizein
       write(c2,"(f7.2)") real(100*cntdif)/real(sizein)
       call MsgDump(trim(msg)//": there are "//trim(adjustl(c0))//&
            " differences in "//trim(adjustl(c1))//" entries; ("//&
            trim(adjustl(c2))//" %)")

       ! case one array is null

       if (zero1) then
          write(c2,"(e12.3)") maxar2
          call MsgDump(trim(msg)//" first null; max abs second="//&
               trim(adjustl(c2)))
       else if (zero2) then
          write(c2,"(e12.3)") maxar1
          call MsgDump(trim(msg)//" second null; max abs second="//&
               trim(adjustl(c2)))
       else

          ! both arrays not null

          write(c0,"(i8)") maxbit
          write(c1,"(i8)") maxbitsmantissa 
          call MsgDump(trim(msg)//" differ; max rel dif: "//&
               trim(adjustl(c0))//" mantissa bits in "//&
               trim(adjustl(c1)))
          write(c0,"(e12.3)") maxdif
          write(c1,"(e12.3)") spacing(maxdif)
          write(c2,"(e12.3)") maxar1
          write(c3,"(e12.3)") maxar2
          call MsgDump(trim(msg)//" dif="//trim(adjustl(c0))//&
               "; spacing="//trim(adjustl(c1))//&
               "; entry1="//trim(adjustl(c2))//&
               "; entry2="//trim(adjustl(c3)))

          do i = 1, maxbins
             write(c0,"(f7.2)") 100.0*real(bins(i))/real(cntdif)
             write(c1,"(i8)") lboundbins(i)
             write(c2,"(i8)") uboundbins(i)
             call MsgDump(trim(msg)//" "//trim(adjustl(c0))//&
                  "% differences from "//trim(adjustl(c1))//&
                  " to "//trim(adjustl(c2))//" bits ")
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
    character(len=16) :: strReal(3)
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
       call MsgDump(h//" for field "//trim(adjustl(msg))//&
            " localy declared ("//&
            "1:"//trim(adjustl(str(1)))//","//&
            "1:"//trim(adjustl(str(2)))//","//&
            "1:"//trim(adjustl(str(3)))//"):")

       do ind3 = 1, d3a1
          do ind2 = 1, d2a1
             do ind1 = 1, d1a1
                if (a1(ind1,ind2,ind3) /= a2(ind1,ind2,ind3)) then
                   call TwoEntriesR4(a1(ind1,ind2,ind3), a2(ind1,ind2,ind3), &
                        maxdif, maxar1, maxar2, maxbit, maxbitsmantissa, bins)
                   if (verb) then
                      write(str(1),"(i8)") ind1
                      write(str(2),"(i8)") ind2
                      write(str(3),"(i8)") ind3
                      write(strReal(1),"(e16.7)") a1(ind1,ind2,ind3)
                      write(strReal(2),"(e16.7)") a2(ind1,ind2,ind3)
                      write(strReal(3),"(e16.7)") abs(a1(ind1,ind2,ind3)-a2(ind1,ind2,ind3))
                      call MsgDump(h//" diff at ("//&
                           trim(adjustl(str(1)))//","//&
                           trim(adjustl(str(2)))//","//&
                           trim(adjustl(str(3)))//")"//&
                           "; entry1="//trim(adjustl(strReal(1)))//&
                           "; entry2="//trim(adjustl(strReal(2)))//&
                           "; diff="//trim(adjustl(strReal(3))))
                   end if
                end if
             end do
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
    character(len=8) :: c0, c1, c2, c3

    if (cntdif == 0_i8) then

       if (verb) then

          ! no differences; verify if any array is null

          if (zero1 .and. zero2) then
             call MsgDump(trim(msg)//" both null")
          else
             call MsgDump(trim(msg)//" matches")
          end if
       end if

    else

       ! there are differences

       write(c0,"(i8)") cntdif
       write(c1,"(i8)") sizein
       write(c2,"(f7.2)") real(100*cntdif)/real(sizein)
       call MsgDump(trim(msg)//": there are "//trim(adjustl(c0))//&
            " differences in "//trim(adjustl(c1))//" entries; ("//&
            trim(adjustl(c2))//" %)")

       if (zero1) then
          write(c2,"(i8)") maxar2
          call MsgDump(trim(msg)//" first null; max abs second="//&
               trim(adjustl(c2)))
       else if (zero2) then
          write(c2,"(i8)") maxar1
          call MsgDump(trim(msg)//" second null; max abs second="//&
               trim(adjustl(c2)))
       else

          ! both arrays not null

          write(c0,"(i8)") maxdif
          write(c2,"(i8)") maxar1
          write(c3,"(i8)") maxar2
          call MsgDump(trim(msg)//" dif="//trim(adjustl(c0))//&
               "; entry1="//trim(adjustl(c2))//&
               "; entry2="//trim(adjustl(c3)))
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

!**(JP)** Fim Insercao

!----------------------------------------------------------------------!
! Optional advection scheme for CCATT-BRAMS/BRAMS models version 4.2+  !
! Based on Walcek, 2000 (JGR) and Walcek and Aleksic, 1998 (ATENV).    !
! The scheme is highly conservative, monotonic and keeps mass mixing   !
! ratio positive definite. 					       !
! Implemented by Saulo Freitas (saulo.freitas@cptec.inpe.br) @ Jun/2009!
! MPI/Paralelized by L. Flavio/J. Panneta                              !
!----------------------------------------------------------------------!

module ModMonotonicAdvection

  use ModCompare, only: Compare

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModMessageSet, only: &
       UpdateFieldAdressAtAdvMnt, &
       PostSendRecvMsgs, &
       WaitSendRecvMsgs

  use ModGridDims, only: &
       GridDims

  use ModDomainDecomp, only: &
       DomainDecomp

  use ModGrid, only: &
       Grid, &
       DumpGrid, &
       DeepCopyToScalarTabAtOneGrid, &
       DeepCopyFromScalarTabAtOneGrid       

  use mem_grid, only:        &
       dtlt,   & !intent(in)
       time,   &
       ngrids, & !intent(in)
       ngrid,  & !intent(in)
       dzt,    & !intent(in)
       dztn,   & !intent(in)
       grid_g, & !intent(in)
       grid_g, & !intent(in)
       naddsc, & !intent(in)
       hw4   , & !intent(in)
       if_adap,& !intent(in)
       dyncore_flag  !intent(in)

  use mem_basic, only: &
       basic_g  !intent(in)

  use micphys, only: &
       level !intent(in)

  use rconstants, only: &
       cp,p00,cv,rgas,cpi   !intent(in)

  use mem_aer1, only: &
       aerosol,    &       !intent(in)
       num_scalar_aer_1st !intent(in)

  use mem_chem1, only: &
       nspecies_transported !intent(in)

  use module_dry_dep, only: &
       sedim_type,          &
       dd_sedim,            &
       naer_transported

  use var_tables, only : scalar_tab & ! (var_p = IN, var_t = INOUT)
       ,num_scalar   ! (in)

  use advMessageMod, only: &
       SendMessageI, &
       RecvMessageI, &
       SendMessageJ, &
       RecvMessageJ, &
       newM2, &
       newM3, &
       newIa, &
       newIz, &
       newJa, &
       newJz, &
       nRecvI, &
       nRecvJ, &
       nSendI, &
       nSendJ, &
       totalrecvi, &
       totalsendi, &
       totalrecvj, &
       totalsendj

  use ParLib, only: &
       parf_send_noblock_real, &
       parf_get_noblock_real, &
       parf_wait_any_nostatus, &
       parf_wait_all_nostatus


  use ModNamelistFile, only: &
       NamelistFile

  use ccatt_start, only: &
       ccatt               ! (in)

  implicit none

  private

  type MonotonicAdvection
     real,pointer :: u3d(:,:,:)
     real,pointer :: v3d(:,:,:)
     real,pointer :: w3d(:,:,:)
     real,pointer :: vc3d_in(:,:,:)
     real,pointer :: vc3d_out(:,:,:)
     real,pointer :: vc3d_x(:,:,:)
     real,pointer :: vc3d_y(:,:,:)
     real,pointer :: dd0_3d(:,:,:)
     real,pointer :: dd0_3du(:,:,:)
     real,pointer :: dd0_3dv(:,:,:)
     real,pointer :: dd0_3dw(:,:,:)
     real,pointer :: den0_3d(:,:,:)
     real,pointer :: den1_3d(:,:,:)
     real,pointer :: den2_3d(:,:,:)
     real,pointer :: den3_3d(:,:,:)
     real,pointer :: l_dxtW(:,:,:)
     real,pointer :: l_dytW(:,:,:)
     real,pointer :: dxtW(:,:)
     real,pointer :: dytW(:,:)
     real,pointer :: dztW(:)
  end type MonotonicAdvection

  public :: MonotonicAdvection
  public :: CreateMonotonicAdvection
  public :: DestroyMonotonicAdvection
  public :: advmnt_driver  ! Subroutine
  public :: StoreNamelistFileAtRadvc_mnt ! Subroutine

  ! public names, set by StoreNamelistFileAtRadvc_mnt
  integer, public :: advmnt 
  integer, public :: GhostZoneLength 

  ! module private variables

  ! flow control flags
  integer, parameter :: ON=1,OFF=0
  integer, parameter :: use_true_density  = 1 ! 0= OFF, 1=ON

  ! for theoretical experiments
  integer, parameter :: theor_wind = 0        ! 0= OFF, 1=ON

  ! constants
  real, parameter :: c1 = cv/rgas
  real, parameter :: c2 = p00/rgas

  ! all fields with enlarged ghost zones
  type advmnt_vars
     real,pointer :: u3d(:,:,:)
     real,pointer :: v3d(:,:,:)
     real,pointer :: w3d(:,:,:)
     real,pointer :: vc3d_in(:,:,:)
     real,pointer :: vc3d_out(:,:,:)
     real,pointer :: vc3d_x(:,:,:)
     real,pointer :: vc3d_y(:,:,:)
     real,pointer :: dd0_3d(:,:,:)
     real,pointer :: dd0_3du(:,:,:)
     real,pointer :: dd0_3dv(:,:,:)
     real,pointer :: dd0_3dw(:,:,:)
     real,pointer :: den0_3d(:,:,:)
     real,pointer :: den1_3d(:,:,:)
     real,pointer :: den2_3d(:,:,:)
     real,pointer :: den3_3d(:,:,:)
     real,pointer :: l_dxtW(:,:,:)
     real,pointer :: l_dytW(:,:,:)
     real,pointer :: dxtW(:,:)
     real,pointer :: dytW(:,:)
     real,pointer :: dztW(:)
  end type advmnt_vars

  ! single variable containing all enlarged ghost zone fields
  ! for all grids
  type(advmnt_vars), allocatable :: advmnt_g(:)

  ! advmnt_g initialization flag
  ! advmnt_g should be initialized only once
!!$  integer :: mnt_adv_jnitialized=0 ! 0=not initialized; 1=initialized

  integer :: nSend_i
  integer :: nSend_j
  integer :: nRecv_i
  integer :: nRecv_j
  integer, parameter :: bigdump=1
  real, allocatable :: buffcomm(:,:,:)
  integer :: nRec_i
  integer :: nSnd_i
  integer :: nRec_j
  integer :: nSnd_j
  integer :: bufSendTotalLength_i
  integer :: bufSendTotalLength_j
  integer :: bufREcvTotalLength_i
  integer :: bufRecvTotalLength_j
  real, allocatable :: bufRecv(:)
  real, allocatable :: bufSend(:)


contains




  function CreateMonotonicAdvection(oneGrid) result(oneAdvMnt)
    type(Grid), pointer, intent(in) :: oneGrid
    type(MonotonicAdvection), pointer :: oneAdvMnt

    integer :: mzpAdvMnt
    integer :: mxpAdvMnt
    integer :: mypAdvMnt

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateMonotonicAdvection)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneGrid)) then
       call fatal_error(h//" oneGrid not associated")
    end if
    if (.not. associated(oneGrid%NodeDimsAdvMnt)) then
       call fatal_error(h//" oneGrid%NodeDimsAdvMnt not associated")
    end if

    mzpAdvMnt = oneGrid%NodeDimsAdvMnt%mzp
    mxpAdvMnt = oneGrid%NodeDimsAdvMnt%mxp
    mypAdvMnt = oneGrid%NodeDimsAdvMnt%myp

    allocate(oneAdvMnt, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%u3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%u3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%v3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%v3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%w3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%w3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dd0_3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dd0_3d  fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dd0_3du(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dd0_3du fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dd0_3dv(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dd0_3dv fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dd0_3dw(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dd0_3dw fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%den0_3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%den0_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%den1_3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%den1_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%den2_3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%den2_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%den3_3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%den3_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%l_dxtW(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%l_dxtW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%l_dytW(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%l_dytW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dxtW(mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dxtW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dytW(mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dytW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dztW(mzpAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dztW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%vc3d_in(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%vc3d_in  fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%vc3d_out(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%vc3d_out fails with stat="//&
            trim(adjustl(str(1))))
    end if

    if (dumpLocal) then
       write(str(2),"(i8)") mzpAdvMnt
       write(str(3),"(i8)") mxpAdvMnt
       write(str(4),"(i8)") mypAdvMnt
       call MsgDump(h//" allocated oneAdvMnt%u3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%v3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%w3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dd0_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dd0_3du("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dd0_3dv("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dd0_3dw("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%den0_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%den1_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%den2_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%den3_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%l_dxtW("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%l_dytW("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dxtW("//&
            trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dytW("//&
            trim(adjustl(str(3)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dztW("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%vc3d_in("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%vc3d_out("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
    end if
  end function CreateMonotonicAdvection





  subroutine DestroyMonotonicAdvection(oneAdvMnt)
    type(MonotonicAdvection), pointer, intent(inout) :: oneAdvMnt

    integer :: ierr
    character(len=8) :: str(10)
    logical :: dumpLocal=.false.
    character(len=*), parameter :: h="**(DestroyMonotonicAdvection)**"

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    if (associated(oneAdvMnt)) then
       deallocate(oneAdvMnt%u3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%u3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%v3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%v3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%w3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%w3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dd0_3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dd0_3d  fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dd0_3du, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dd0_3du fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dd0_3dv, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dd0_3dv fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dd0_3dw, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dd0_3dw fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%den0_3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%den0_3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%den1_3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%den1_3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%den2_3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%den2_3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%den3_3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%den3_3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%l_dxtW, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%l_dxtW fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%l_dytW, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%l_dytW fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dxtW, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dxtW fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dytW, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dytW fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dztW, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dztW fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%vc3d_in, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%vc3d_in  fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%vc3d_out, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%vc3d_out fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    nullify(oneAdvMnt)

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine DestroyMonotonicAdvection





  subroutine InitializeGridSpacings(&
       mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
       iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
       jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
       dxt, dyt, fmapt, rtgt, dztn, &
       dxtW, dytW, dztW)

    ! computes new grid spacing on x, y and z
    ! for Wal2cek monotonic advection

    integer, intent(in) :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer, intent(in) :: mxp
    ! x dimension of external fields 
    integer, intent(in) :: myp
    ! y dimension of external fields 
    integer, intent(in) :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer, intent(in) :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer, intent(in) :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    real, pointer, intent(in) :: dxt(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: dyt(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: fmapt(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: rtgt(:,:)
    ! external field, pointer and values are intent(in)
    real, intent(in) :: dztn(:)
    ! external field
    real, pointer, intent(in) :: dxtW(:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: dytW(:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: dztW(:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)

    ! local var

    integer :: iExtern
    integer :: jExtern
    integer :: i
    integer :: j
    integer :: k
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(InitializeGridSpacings)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") mzp
       write(str(2),"(i8)") i1ExternAtAdvMnt
       write(str(3),"(i8)") iMxpExternAtAdvMnt
       write(str(4),"(i8)") j1ExternAtAdvMnt
       write(str(5),"(i8)") jMypExternAtAdvMnt
       call MsgDump(h//" set values of"//&
            " dxtW, dytW both at ("//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")"//&
            " and dztW at (1:"//trim(adjustl(str(1)))//"), nullifying remaining positions")
       write(str(1),"(i8)") iOffset
       write(str(2),"(i8)") jOffset
       call MsgDump(h//" external fields index = monotonic fields index + ("//&
            trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//")")
    end if

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          dxtW(i,j) = 0.0
          dytW(i,j) = 0.0
       end do
    end do

    ! fill Monotonic Advection fields only where
    ! external fields are in range
    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       jExtern = j + jOffset
       ! set Monotonic Advection west ghost zone fields to zero 
       do i = 1, i1ExternAtAdvMnt - 1
          dxtW(i,j) = 0.0
          dytW(i,j) = 0.0
       end do
       ! fill where both Monotonic Advection and external fields
       ! are in range
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          rtgti = 1. / rtgt(iExtern,jExtern)

          !- at init/rams_grid.f90:
          !     dxt(i,j)=fmapt(i,j)/(xmn(i,ngrid)-xmn(i-1,ngrid))
          !     dyt(i,j)=fmapt(i,j)/(ymn(j,ngrid)-ymn(j-1,ngrid))

          dxtW(i,j) = 1./(dxt(iExtern,jExtern) * fmapt(iExtern,jExtern) * rtgti)
          dytW(i,j) = 1./(dyt(iExtern,jExtern) * fmapt(iExtern,jExtern) * rtgti)
       end do
       ! set Monotonic Advection east ghost zone fields to zero 
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          dxtW(i,j) = 0.0
          dytW(i,j) = 0.0
       end do
    end do

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          dxtW(i,j) = 0.0
          dytW(i,j) = 0.0
       end do
    end do

    ! z ranges of Monotonic Advection and external fields are identical
    do k = 1,mzp
       !- at init/gridset.f90:
       !  dztn(k,ifm) = 1. / (zmn(k,ifm) - zmn(k-1,ifm))
       ! Por que o Jacobiano nao depende de Z, o dztw depende somente
       ! de z.
       !dztW(k,i,j) = 1./ ( dzt(k) * rtgti * fmapt(i,j)**2 )
       dztW(k) = 1./ ( dztn(k) ) !
    end do
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine InitializeGridSpacings




  subroutine GetTrueDensities(&
       mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
       iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
       jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
       level,&
       rtp, rv, pp, pi0, theta, &
       dd0_3d, dd0_3du, dd0_3dv, dd0_3dw)

    integer, intent(in) :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer, intent(in) :: mxp
    ! x dimension of external fields 
    integer, intent(in) :: myp
    ! y dimension of external fields 
    integer, intent(in) :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer, intent(in) :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer, intent(in) :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    integer, intent(in) :: level
    real, pointer, intent(in) :: rtp(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: rv(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: pp(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: pi0(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: theta(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: dd0_3d(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: dd0_3du(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: dd0_3dv(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: dd0_3dw(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)

    ! local var
    integer :: iExtern
    integer :: jExtern
    integer :: i
    integer :: j
    integer :: k
    integer :: i1
    integer :: j1
    real :: c3

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(GetTrueDensities)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") lbound(dd0_3d,1)
       write(str(2),"(i8)") ubound(dd0_3d,1)
       write(str(3),"(i8)") lbound(dd0_3d,2)
       write(str(4),"(i8)") ubound(dd0_3d,2)
       write(str(5),"(i8)") lbound(dd0_3d,3)
       write(str(6),"(i8)") ubound(dd0_3d,3)
       call MsgDump(h//" dd0_3d dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") lbound(dd0_3du,1)
       write(str(2),"(i8)") ubound(dd0_3du,1)
       write(str(3),"(i8)") lbound(dd0_3du,2)
       write(str(4),"(i8)") ubound(dd0_3du,2)
       write(str(5),"(i8)") lbound(dd0_3du,3)
       write(str(6),"(i8)") ubound(dd0_3du,3)
       call MsgDump(h//" dd0_3du dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") lbound(dd0_3dv,1)
       write(str(2),"(i8)") ubound(dd0_3dv,1)
       write(str(3),"(i8)") lbound(dd0_3dv,2)
       write(str(4),"(i8)") ubound(dd0_3dv,2)
       write(str(5),"(i8)") lbound(dd0_3dv,3)
       write(str(6),"(i8)") ubound(dd0_3dv,3)
       call MsgDump(h//" dd0_3dv dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") lbound(dd0_3dw,1)
       write(str(2),"(i8)") ubound(dd0_3dw,1)
       write(str(3),"(i8)") lbound(dd0_3dw,2)
       write(str(4),"(i8)") ubound(dd0_3dw,2)
       write(str(5),"(i8)") lbound(dd0_3dw,3)
       write(str(6),"(i8)") ubound(dd0_3dw,3)
       call MsgDump(h//" dd0_3dw dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") mzp
       write(str(2),"(i8)") i1ExternAtAdvMnt
       write(str(3),"(i8)") iMxpExternAtAdvMnt
       write(str(4),"(i8)") j1ExternAtAdvMnt
       write(str(5),"(i8)") jMypExternAtAdvMnt
       call MsgDump(h//" set values of"//&
            " dd0_3d, dd0_3du, dd0_3dv, dd0_3dw, all at ("//&
            "1:"//trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")"//&
            ", nullifying remaining positions")
       write(str(1),"(i8)") iOffset
       write(str(2),"(i8)") jOffset
       call MsgDump(h//" external fields index = monotonic fields index + ("//&
            trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//")")
    end if

    ! dd0_3d computation

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
          end do
       end do
    end do

    !- true air density at points "T"

    c3 = c2 * (cpi**c1)
    if( level == 0 ) then
       do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
          jExtern = j + jOffset
          ! set Monotonic Advection west ghost zone fields to zero
          do i = 1, i1ExternAtAdvMnt-1
             do k = 1, mzp
                dd0_3d(k,i,j) = 0.0
             end do
          end do
          ! fill where both Monotonic Advection and external fields
          ! are in range
          do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
             iExtern = i + iOffset
             do k = 1, mzp
                dd0_3d(k,i,j) = (c3/theta(k,iExtern,jExtern))*&
                     (pi0(k,iExtern,jExtern)+pp(k,iExtern,jExtern))**c1
             end do
          end do
          ! set Monotonic Advection east ghost zone fields to zero
          do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
             do k = 1, mzp
                dd0_3d(k,i,j) = 0.0
             end do
          end do
       end do
    else
       do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
          jExtern = j + jOffset
          ! set Monotonic Advection west ghost zone fields to zero
          do i = 1, i1ExternAtAdvMnt-1
             do k = 1, mzp
                dd0_3d(k,i,j) = 0.0
             end do
          end do
          ! fill where both Monotonic Advection and external fields
          ! are in range
          do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
             iExtern = i + iOffset
             do k = 1, mzp
                dd0_3d(k,i,j) = (c3/theta(k,iExtern,jExtern))*&
                     (1. + rtp(k,iExtern,jExtern))/ &
                     (1. + 1.61*rv(k,iExtern,jExtern))*&
                     (pi0(k,iExtern,jExtern)+pp(k,iExtern,jExtern))**c1
             end do
          end do
          ! set Monotonic Advection east ghost zone fields to zero
          do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
             do k = 1, mzp
                dd0_3d(k,i,j) = 0.0
             end do
          end do
       end do
    end if

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
          end do
       end do
    end do


    ! use dd0_3d to compute dd0_3du and dd0_3dv

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
          end do
       end do
    end do

    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       j1 = min(j+1,jMypExternAtAdvMnt)
       ! set Monotonic Advection west ghost zone fields to zero
       do i = 1, i1ExternAtAdvMnt-1
          do k = 1, mzp
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
          end do
       end do
       !- true air density computation
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          i1 = min(i+1,iMxpExternAtAdvMnt)
          do k = 1,mzp
             dd0_3du(k,i,j) = .5 * (dd0_3d(k,i,j) + dd0_3d(k,i1,j))
             dd0_3dv(k,i,j) = .5 * (dd0_3d(k,i,j) + dd0_3d(k,i,j1))
          end do
       end do
       ! set Monotonic Advection east ghost zone fields to zero
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          do k = 1, mzp
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
          end do
       end do
    end do

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
          end do
       end do
    end do

    ! compute true air density for w

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do

    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       ! set Monotonic Advection west ghost zone fields to zero
       do i = 1, i1ExternAtAdvMnt-1
          do k = 1, mzp
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
       ! compute dd0_3dw 
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          do k = 1,mzp-1
             dd0_3dw(k,i,j) = 0.5*(dd0_3d(k,i,j) + dd0_3d(k+1,i,j))
          end do
          dd0_3dw(mzp,i,j)=dd0_3dw(mzp-1,i,j)
       end do
       ! set Monotonic Advection east ghost zone fields to zero
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          do k = 1, mzp
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine GetTrueDensities



  subroutine PrepareWinds(&
       ng, mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
       iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
       jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
       dtlt, uc, up, vc, vp, wc, wp, &
       fmapui, fmapvi, rtgt, rtgu, rtgv, f13t, f23t, &
       u3d, v3d, w3d, &
       aerosol, naer_transported, dd_sedim, dzt, ndt_z)

    integer, intent(in) :: ng ! grid number, should dissapear
    integer, intent(in) :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer, intent(in) :: mxp
    ! x dimension of external fields 
    integer, intent(in) :: myp
    ! y dimension of external fields 
    integer, intent(in) :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer, intent(in) :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer, intent(in) :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    real, intent(in) :: dtlt

    real, pointer, intent(in) :: uc(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: up(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: vc(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: vp(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: wc(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: wp(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: fmapui(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: fmapvi(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: rtgt(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: rtgu(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: rtgv(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: f13t(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: f23t(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: u3d(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: v3d(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: w3d(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    integer, intent(in) :: aerosol
    ! flag for aerosol computation
    integer, intent(in) :: naer_transported
    ! # transported aerosol
    type(sedim_type), intent(in) :: dd_sedim(:,:)
    real, intent(in) :: dzt(:)
    integer, intent(inout) :: ndt_z(:)
    ! aerosol sedimentation timestep control

    !- local var

    integer :: jm
    integer :: jp
    integer :: im
    integer :: ip
    integer :: ispc
    integer :: i
    integer :: j
    integer :: k
    integer :: iExtern
    integer :: jExtern
    real :: cx1
    real :: cx2
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(PrepareWinds)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       call MsgDump(h//" starts; computes u3d, v3d and w3d just at a section"//&
            " restricted to the original ghost zone of 1")
    end if

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             w3d(k,i,j) = 0.0
             u3d(k,i,j) = 0.0
             v3d(k,i,j) = 0.0
          end do
       end do
    end do

    ! u3d, u3d, and w3d are input as the velocity components (averaged
    ! between past and current time levels) times dtlt.
    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       jExtern = j + jOffset
       ! set Monotonic Advection west ghost zone fields to zero
       do i = 1, i1ExternAtAdvMnt-1
          do k = 1, mzp
             w3d(k,i,j) = 0.0
             u3d(k,i,j) = 0.0
             v3d(k,i,j) = 0.0
          end do
       end do
       ! u3d, u3d, and w3d are input as the velocity components (averaged
       ! between past and current time levels) times dtlt.
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          do k = 1, mzp
             w3d(k,i,j) = ( wc(k,iExtern,jExtern) + wp(k,iExtern,jExtern) )*0.5
             u3d(k,i,j) = ( uc(k,iExtern,jExtern) + up(k,iExtern,jExtern) )*0.5
             v3d(k,i,j) = ( vc(k,iExtern,jExtern) + vp(k,iExtern,jExtern) )*0.5
          end do
       end do
       ! set Monotonic Advection east ghost zone fields to zero
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          do k = 1, mzp
             w3d(k,i,j) = 0.0
             u3d(k,i,j) = 0.0
             v3d(k,i,j) = 0.0
          end do
       end do
    end do

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             w3d(k,i,j) = 0.0
             u3d(k,i,j) = 0.0
             v3d(k,i,j) = 0.0
          end do
       end do
    end do

    ! transform w3d from cartesian vertical velocity to sigma_z velocity

    ! Add contribution to w3d from horiz winds crossing sloping sigma surfaces,
    ! and include 1/rtgt factor in w3d
    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       jm = max(j1ExternAtAdvMnt,j-1)
       jp = min(jMypExternAtAdvMnt,j+1)
       jExtern = j + jOffset
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          im = max(i1ExternAtAdvMnt,i-1)
          ip = min(iMxpExternAtAdvMnt,i+1)
          rtgti = 1. / rtgt(iExtern,jExtern)
          do k = 1,mzp-1
             w3d(k,i,j) = &
                  ( &
                  (u3d(k,i,j)+u3d(k+1,i,j)+u3d(k,im,j)+u3d(k+1,im,j)) * f13t(iExtern,jExtern) + &
                  (v3d(k,i,j)+v3d(k+1,i,j)+v3d(k,i,jm)+v3d(k+1,i,jm)) * f23t(iExtern,jExtern)  &
                  ) * hw4(k) + w3d(k,i,j) * rtgti
          end do
       end do
    end do

    ! include map factors on u and v

    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       jExtern = j + jOffset
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          cx1 = fmapui(iExtern,jExtern) * rtgu(iExtern,jExtern)
          cx2 = fmapvi(iExtern,jExtern) * rtgv(iExtern,jExtern)
          do k = 1,mzp-1
             u3d(k,i,j) = u3d(k,i,j) * cx1
             v3d(k,i,j) = v3d(k,i,j) * cx2
          end do
       end do
    end do
    !-----------------------------------------
    !- control for aerosol sedimentation
    if(aerosol > 0 .and. naer_transported > 0) then
       ! very crude estimation of CFL violation and fix for the number of sub-timesteps
       ! for large particles
       do ispc=1,naer_transported
          ndt_z(ispc)=ceiling(maxval(abs(dd_sedim(ispc,ng)%v_sed_part))*dtlt*maxval(dzt(1:mzp)))
       end do
    end if
    !- end of aerosol sedimentation
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine PrepareWinds





  subroutine GetWalceksDensities(&
       mzp, dtlt, mxpAdvMnt, mypAdvMnt, &
       i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
       j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
       u3d, v3d, w3d, &
       dd0_3d, dd0_3du, dd0_3dv, dd0_3dw, &
       dxtW, dytW, dztW, &
       den0_3d, den1_3d, den2_3d, den3_3d)

    integer, intent(in) :: mzp
    real, intent(in) :: dtlt
    integer, intent(in) :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer, intent(in) :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    real, pointer, intent(in) :: u3d(:,:,:)
    real, pointer, intent(in) :: v3d(:,:,:)
    real, pointer, intent(in) :: w3d(:,:,:)
    real, pointer, intent(in) :: dd0_3d(:,:,:)
    real, pointer, intent(in) :: dd0_3du(:,:,:)
    real, pointer, intent(in) :: dd0_3dv(:,:,:)
    real, pointer, intent(in) :: dd0_3dw(:,:,:)
    real, pointer, intent(in) :: dztW(:)
    real, pointer, intent(in) :: dxtW(:,:)
    real, pointer, intent(in) :: dytW(:,:)
    real, pointer, intent(in) :: den0_3d(:,:,:)
    real, pointer, intent(in) :: den1_3d(:,:,:)
    real, pointer, intent(in) :: den2_3d(:,:,:)
    real, pointer, intent(in) :: den3_3d(:,:,:)


    ! local var
    integer i,j,k

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(GetWalceksDensities)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       call MsgDump(h//" starts; computes den0_3d, den1_3d, den2_3d and den3_3d"//&
            " just at a section restricted to the original ghost zone of 1")
    end if


    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             den0_3d(k,i,j) = 0.0
             den1_3d(k,i,j) = 0.0
             den2_3d(k,i,j) = 0.0
             den3_3d(k,i,j) = 0.0
          end do
       end do
    end do

    do  j = jMypExternAtAdvMnt, j1ExternAtAdvMnt+1, -1
       ! set Monotonic Advection west ghost zone fields to zero
       do i = 1, i1ExternAtAdvMnt
          do k = 1, mzp
             den0_3d(k,i,j) = 0.0
             den1_3d(k,i,j) = 0.0
             den2_3d(k,i,j) = 0.0
             den3_3d(k,i,j) = 0.0
          end do
       end do

       do  i = i1ExternAtAdvMnt+1, iMxpExternAtAdvMnt
          den0_3d(1,i,j) = 0.0
          den1_3d(1,i,j) = 0.0
          den2_3d(1,i,j) = 0.0
          den3_3d(1,i,j) = 0.0
          do k = 2, mzp
             den0_3d(k,i,j)=dd0_3d(k,i,j)
             den1_3d(k,i,j)=den0_3d(k,i,j)- dtlt/dxtW(i,j)*&
                  (dd0_3du(k,i,j)*u3d(k,i,j)-dd0_3du(k,i-1,j)*u3d(k,i-1,j))
             den2_3d(k,i,j)=den1_3d(k,i,j)- dtlt/dytW(i,j)*&
                  (dd0_3dv(k,i,j)*v3d(k,i,j)-dd0_3dv(k,i,j-1)*v3d(k,i,j-1))
             den3_3d(k,i,j)=den2_3d(k,i,j)- dtlt/dztW(k)  *&
                  (dd0_3dw(k,i,j)*w3d(k,i,j)-dd0_3dw(k-1,i,j)*w3d(k-1,i,j))
          end do
       end do
       ! set Monotonic Advection east ghost zone fields to zero
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          do k = 1, mzp
             den0_3d(k,i,j) = 0.0
             den1_3d(k,i,j) = 0.0
             den2_3d(k,i,j) = 0.0
             den3_3d(k,i,j) = 0.0
          end do
       end do
    end do

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             den0_3d(k,i,j) = 0.0
             den1_3d(k,i,j) = 0.0
             den2_3d(k,i,j) = 0.0
             den3_3d(k,i,j) = 0.0
          end do
       end do
    end do

    !srf- BC for den3_3d
    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       do k = 1, mzp
          den3_3d(k,i1ExternAtAdvMnt,j)=den3_3d(k,i1ExternAtAdvMnt+1,j)
       end do
    end do

    do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
       do k = 1, mzp
          den3_3d(k,i,j1ExternAtAdvMnt)=den3_3d(k,i,j1ExternAtAdvMnt+1)
       end do
    end do
  end subroutine GetWalceksDensities








  subroutine CheckBorders(m1, m2, m3, field, &
       nRecv, procRecv, tagRecv, iaRecv, izRecv, jaRecv, jzRecv, &
       bufRecvStart, bufRecvLength, bufRecvTotalLength, &
       nSend, procSend, tagSend, iaSend, izSend, jaSend, jzSend, &
       bufSendStart, bufSendLength, bufSendTotalLength,mynum,op,ie)
    integer, intent(in) :: m1,mynum,op,ie
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    real,    intent(inout) :: field(m1,m2,m3)
    integer, intent(in) :: nRecv
    integer, intent(in) :: procRecv(nRecv)
    integer, intent(in) :: tagRecv(nRecv)
    integer, intent(in) :: iaRecv(nRecv)
    integer, intent(in) :: izRecv(nRecv)
    integer, intent(in) :: jaRecv(nRecv)
    integer, intent(in) :: jzRecv(nRecv)
    integer, intent(in) :: bufRecvStart(nRecv)
    integer, intent(in) :: bufRecvLength(nRecv)
    integer, intent(in) :: bufRecvTotalLength
    integer, intent(in) :: nSend
    integer, intent(in) :: procSend(nSend)
    integer, intent(in) :: tagSend(nSend)
    integer, intent(in) :: iaSend(nSend)
    integer, intent(in) :: izSend(nSend)
    integer, intent(in) :: jaSend(nSend)
    integer, intent(in) :: jzSend(nSend)
    integer, intent(in) :: bufSendStart(nSend)
    integer, intent(in) :: bufSendLength(nSend)
    integer, intent(in) :: bufSendTotalLength

    integer :: fout,i
    character :: opc

    fout=80+mynum
    opc='Y'
    if(op==1) opc='X'

    if(ie==0) then
       write (fout,'("Borders updated, direction ",A)') opc
       return
    end if


    write (fout,'(" Updating borders, direction  ",A)') opc ; call flush(fout)
    write (fout,'("nRecv: ",I3.3," nSend: ",I3.3)') nRecv,nSend; call flush(fout)
    write (fout,'("TotRecv: ",I6.6," TotSend: ",I6.6)') bufRecvTotalLength,bufSendTotalLength; call flush(fout)
    write (fout,'(A)') '---------------------------------- Send ---------------------------'; call flush(fout)
    write (fout,'(7(A3,1X),2(A,1X))') 'nSn','prc','tag','ia ','iz ','ja ','jz ','Start','Length'; call flush(fout)
    do i=1,nRecv
       write (fout,'(7(I3.3,1X),2(I6.6,1X))') i,procRecv(i),tagRecv(i),iaRecv(i),izRecv(i),jaRecv(i),jzRecv(i), &
            bufRecvStart(i),bufRecvLength(i); call flush(fout)
    end do
    write (fout,'(A)') '------------------------------- Receive  --------------------------'; call flush(fout)
    write (fout,'(7(A3,1X),2(A,1X))') 'nRv','prc','tag','ia ','iz ','ja ','jz ','Start','Length'
    do i=1,nRecv
       write (fout,'(7(I3.3,1X),2(I6.6,1X))')	i,procSend(i),tagSend(i),iaSend(i),izSend(i),jaSend(i),jzSend(i), &
            bufSendStart(i),bufSendLength(i); call flush(fout)
    end do
  end subroutine CheckBorders





  subroutine AdvectTendency(mzp, mxp, &
       iOffset, jOffset, &
       ia, iz, ja, jz, dtl, &
       scalarp3D, AdvMntField, scalart1D)
    integer, intent(in) :: mzp ! external field dimension
    integer, intent(in) :: mxp ! external field dimension
    integer, intent(in) :: iOffset
    integer, intent(in) :: jOffset
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(in) :: dtl
    real, pointer, intent(in) :: scalarp3D(:,:,:)
    ! external field, regular ghost zone, values are intent(in)
    real, pointer, intent(in) :: AdvMntField(:,:,:)
    ! monotonic advection field, extended ghost zone, values are intent(in)
    real, pointer, intent(in) :: scalart1D(:)
    ! 3D external field with regular ghost zone mapped into a 1D field;
    ! values are intent(inout)

    integer :: i
    integer :: iAdvMnt
    integer :: j
    integer :: jAdvMnt
    integer :: k
    integer :: ind1D
    real :: dtli

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(AdvectTendency)**"
    character(len=8) :: str(10)
    if (dumpLocal) then
       write(str(1),"(i8)") mzp-1
       write(str(2),"(i8)") ia
       write(str(3),"(i8)") iz
       write(str(4),"(i8)") ja
       write(str(5),"(i8)") jz
       call MsgDump(h//" at (2:"//trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
    end if

    dtli = 1. / dtl

    ! indexing external field

    do j = ja,jz
       jAdvMnt=j-jOffset
       do i = ia,iz
          iAdvMnt=i-iOffset
          ind1D=(i-1+(j-1)*mxp)*mzp
          do k = 2,mzp-1
             scalart1D(k+ind1D) = scalart1D(k+ind1D) + &
                  (AdvMntField(k,iAdvMnt,jAdvMnt)-scalarp3D(k,i,j)) * dtli
          end do
       end do
    end do
  end subroutine AdvectTendency






  subroutine AdvectMnt(oneAdvMnt, oneGrid, ngrid,m1,m2,m3,ia,iz,ja,jz,dt,mynum,n,&
       current_aer_ispc,current_ndt_z,IsThisScalarAer)

    type(MonotonicAdvection), pointer, intent(in) :: oneAdvMnt
    type(Grid), pointer, intent(in) :: oneGrid
    integer , intent(in) :: m1,ngrid
    integer , intent(in) :: m2
    integer , intent(in) :: m3
    integer , intent(in) :: ia
    integer , intent(in) :: iz
    integer , intent(in) :: ja
    integer , intent(in) :: jz,n
    integer , intent(in) :: mynum
    real    , intent(in) :: dt
    integer , intent(in) :: current_ndt_z,current_aer_ispc
    logical , intent(in) :: IsThisScalarAer
    !- local var
    !REAL,DIMENSION(m1)               :: dxx
    !REAL,DIMENSION(m2,m3)            :: dxy
    real masscon,initialmass,vol
    integer nrec,itz
    integer ibegin,iend,jbegin,jend
    !- type of sedimentation scheme (0= Walcek, 1=upwind)
    integer , parameter :: iupwind = 0
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(AdvectMnt)**"
    character(len=8) :: str(10)

    iBegin= newIa(ngrid)-1
    iEnd  = newIz(ngrid)+1
    jBegin= newJa(ngrid)-1
    jEnd  = newJz(ngrid)+1

    !--- do X-advection
    if (dumpLocal) then
       call MsgDump(h//" starts; update borders of vc3d_in for x advection")
    end if
    call PostSendRecvMsgs(oneGrid%AdvMntScaSendX, oneGrid%AdvMntScaRecvX)
    call WaitSendRecvMsgs(oneGrid%AdvMntScaSendX, oneGrid%AdvMntScaRecvX)
    if (dumpLocal) then
       call MsgDump(h//" invoke Advec3DX to advect vc3d_in, storing result in vc3d_out")
    end if
    call Advec3DX(m1,newM2(ngrid),newM3(ngrid),2,newM2(ngrid)-1,2,newM3(ngrid)-1, &
         oneAdvMnt%vc3d_in,                             &
         oneAdvMnt%u3d,oneAdvMnt%den0_3d,         &
         oneAdvMnt%den1_3d,dt,oneAdvMnt%dxtW,     &
         oneAdvMnt%dd0_3du,                             &
         oneAdvMnt%vc3d_out, mynum )

    !--- do Y-advection

    if (dumpLocal) then
       call MsgDump(h//" update borders of vc3d_out for y advection")
    end if
    call PostSendRecvMsgs(oneGrid%AdvMntScaSendY, oneGrid%AdvMntScaRecvY)
    call WaitSendRecvMsgs(oneGrid%AdvMntScaSendY, oneGrid%AdvMntScaRecvY)
    if (dumpLocal) then
       call MsgDump(h//" invoke Advec3DY to advect vc3d_out, storing result in vc3d_in")
    end if
    call Advec3DY(m1,newM2(ngrid),newM3(ngrid),2,newM2(ngrid)-1,2,newM3(ngrid)-1, &
         oneAdvMnt%vc3d_out,                                        &
         oneAdvMnt%v3d,oneAdvMnt%den1_3d,                     &
         oneAdvMnt%den2_3d,dt,oneAdvMnt%dytW,                 &
         oneAdvMnt%dd0_3dv,                                         &
         oneAdvMnt%vc3d_in, mynum )

    !--- do k-advection
    if (dumpLocal) then
       call MsgDump(h//" invoke Advec3d_Z to advect vc3d_in, storing result in vc3d_out")
    end if
    call Advec3DZ(m1,newM2(ngrid),newM3(ngrid),ibegin,iend,jbegin,jend, &
         oneAdvMnt%vc3d_in,                                   &
         oneAdvMnt%w3d,oneAdvMnt%den2_3d,               &
         oneAdvMnt%den3_3d,dt,oneAdvMnt%dztW,           &
         oneAdvMnt%dd0_3dw,                                   &
         oneAdvMnt%vc3d_out, mynum )


    !- aerosol section to include sedimentation
    !- the sedimentation process is done using pure cartesian coordinates
    !- so, all sedimentation velocities are treat as cartesian vertical velocities
    !- which are positive downwards.
    if (dumpLocal) then
       write(str(1),"(i8)") aerosol
       write(str(2),"(l)") IsThisScalarAer
       call MsgDump(h//" aerosol="//trim(adjustl(str(1)))//&
            "; IsThisScalarAer="//trim(adjustl(str(2))))
    end if
    if(aerosol > 0 .and. IsThisScalarAer) then

       !-srf introducing a time-splitting for aerosol sedimentation

       if (dumpLocal) then
          write(str(1),"(i8)") iupwind
          call MsgDump(h//" iupwind="//trim(adjustl(str(1))))
       end if
       if(iupwind == 0 ) then
          ! - Walcek method
          ! this routine works _only_ for mass concentration or density (kg/m3)
          ! converting mixing ratio (kg/kg) to density (kg/m3)
          oneAdvMnt%vc3d_in(:,:,:)=oneAdvMnt%vc3d_out(:,:,:) * oneAdvMnt%den0_3d(:,:,:)

          !- do time splitting for aerosols with large fall velocities
          do itz=1,current_ndt_z
             call Advec3DZSedim(m1,m2,m3,ia,iz,ja,jz,                        &
                  oneAdvMnt%vc3d_in(:,iBegin:iEnd,jBegin:jEnd),	 &
                  dd_sedim(current_aer_ispc,ngrid)%v_sed_part,          & !fall velocity
                  dt/float(current_ndt_z),                              & !subtimestep
                  dzt(1:m1),grid_g(ngrid)%rtgt,	                 &
                  oneAdvMnt%vc3d_out(:,iBegin:iEnd,jBegin:jEnd),  &
                  mynum )

             ! copy output to input array for the next sup-timestep
             if(itz < current_ndt_z) oneAdvMnt%vc3d_in(:,:,:)=oneAdvMnt%vc3d_out(:,:,:)

          end do
          ! converting back mass concentration to mixing ratio
          oneAdvMnt%vc3d_out(:,:,:)=&
               oneAdvMnt%vc3d_out(:,:,:)/&
               oneAdvMnt%den0_3d(:,:,:)

       else if(iupwind == 1 ) then
          ! - upwind method
          !- do time splitting for aerosols with large fall velocities
          do itz=1,current_ndt_z
             call Advec3DZSedimUpw(m1,m2,m3,ia,iz,ja,jz,                          &
                  dd_sedim(current_aer_ispc,ngrid)%v_sed_part,          & !fall velocity
                  dt/float(current_ndt_z),                              & !subtimestep
                  dzt(1:m1),grid_g(ngrid)%rtgt,	                                 &
                  oneAdvMnt%vc3d_out(:,iBegin:iEnd,jBegin:jEnd),  &
                  mynum )

          end do
       end if
    end if
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine AdvectMnt


  subroutine Advec3DX(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,den0,&
       den1,dt,dxx,&
       dd0,&
       qn,mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: den0(m1,m2,m3)
    real   , intent(in) :: den1(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dxx(m2,m3)
    real   , intent(in) :: dd0(m1,m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    integer :: ii
    integer :: ji
    integer :: ii0
    integer :: ji0
    integer :: ie
    integer :: je
    integer :: ie0
    integer :: je0
    integer :: ipos
    integer :: iia
    integer :: iiz
    integer :: nvar
    integer :: nf
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3DX)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    ! copy input field to output field
    qn=q0
    imxmn=.false.

    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do j=ja,jz
       do k=2,m1-1
          if(u(k,1,j)>=zr0) flux(k,1,j)= q0(k,1,j)*u(k,1,j)*dt*dd0(k,1,j)
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
    do j=ja,jz
       do  i=2,m2-1 ! ia,iz-1 or 1,iz-1
          do k=2,m1-1
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k,i-1,j),q0(k,i+1,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k,i-1,j),q0(k,i+1,j))+eps)        !       extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k,i,j  )< zr0) ck1= q0(k,i+1,j)
             if(u(k,i-1,j)>=zr0) ck2= q0(k,i-1,j)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    do j=ja,jz
       do  i=2,m2-1 ! ia,iz-1 or 1,iz-1
          do k=2,m1-1
             if(u(k,i,j)<zr0) cycle
             if(u(k,i-1,j)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(i,j)               ! Courant number
                x1n= (1.-x1)*(q0(k,i+1,j)-q0(k,i-1,j))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k,i-1,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i+1,j)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i+1,j))  ), max(q0(k,i,j),q0(k,i+1,j)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j) = max(vcmin(k,i,j),min(vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k,i-1,j)/dxx(i,j))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(i,j)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k,i-1,j)
             end if                                                  !Eq-9a
          end do
       end do
    end do

    ! If periodic boundary conditions are assumed, it is necessary
    !   to recalculate the updated mixing ratio at cell 1 if there
    !   is inflow to that cell from the boundary between IDIM and 1
    !   Here these statements are commented out, but should be uncommented
    !   if this subroutine is needed for periodic boundary conditions,
    !   and then one of the calling arguements to the subroutine is IPERIOD
    !   which is set to "1" if you assume period boundary conditions
    !      IF(IPERIOD==1) THEN
    !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
    !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
    !      END IF
    !
    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    !
    do j=ja,jz
       do k=2,m1-1
          if(u(k,m2-1,j)<zr0) flux(k,m2-1,j)= &
               q0(k,m2,j)*u(k,m2-1,j)*dt*dd0(k,m2-1,j)
       end do
    end do

    do j=ja,jz
       do i=m2-1,2,-1 !iz,ia,-1
          do k=2,m1-1
             if(u(k,i-1,j)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j) + &
                     flux(k,i-1,j)/dxx(i,j))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k,i-1,j))/dxx(i,j)     ! Courant number
                x1n= (1.-x1)*(q0(k,i-1,j)-q0(k,i+1,j))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k,i+1,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i-1,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i-1,j)) ), max(q0(k,i,j),q0(k,i-1,j)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j)= max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j)-x1*cf1*dd0(k,i-1,j))/den1(k,i,j) ))
                flux(k,i-1,j)=dxx(i,j)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop y-z
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DX



  subroutine Advec3DY(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,den0,&
       den1,dt,dxx,&
       dd0,&
       qn,mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: den0(m1,m2,m3)
    real   , intent(in) :: den1(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dxx(m2,m3)
    real   , intent(in) :: dd0(m1,m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    integer :: ii
    integer :: ji
    integer :: ii0
    integer :: ji0
    integer :: ie
    integer :: je
    integer :: ie0
    integer :: je0
    integer :: ipos
    integer :: iia
    integer :: iiz
    integer :: nvar
    integer :: nf
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3DY)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    ! copy input field to output field
    qn= q0
    imxmn=.false.

    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do i=ia,iz
       do k=2,m1-1
          if(u(k,i,1)>=zr0) flux(k,i,1)= q0(k,i,1)*u(k,i,1)*dt*dd0(k,i,1)
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !	mixing ratio at t+dt. If these limits are ever violated,
    !	non-monotonic (oscillatory) behavior in solution results
    do i=ia,iz
       do  j=2,m3-1 ! ja,jz
          do k=2,m1-1
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k,i,j-1),q0(k,i,j+1))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k,i,j-1),q0(k,i,j+1))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k,i,j  )< zr0) ck1= q0(k,i,j+1)
             if(u(k,i,j-1)>=zr0) ck2= q0(k,i,j-1)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    do i=ia,iz
       do  j=2,m3-1 ! ja,jz
          do k=2,m1-1
             if(u(k,i,j)<zr0) cycle
             if(u(k,i,j-1)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(i,j)               ! Courant number
                x1n= (1.-x1)*(q0(k,i,j+1)-q0(k,i,j-1))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k,i,j-1)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i,j+1)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i,j+1))  ), max(q0(k,i,j),q0(k,i,j+1)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j)= max(  vcmin(k,i,j),   min(   vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k,i,j-1)/dxx(i,j))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(i,j)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k,i,j-1) !Eq-9a
             end if
          end do
       end do
    end do

    ! If periodic boundary conditions are assumed, it is necessary
    !   to recalculate the updated mixing ratio at cell 1 if there
    !   is inflow to that cell from the boundary between IDIM and 1
    !   Here these statements are commented out, but should be uncommented
    !   if this subroutine is needed for periodic boundary conditions,
    !   and then one of the calling arguements to the subroutine is IPERIOD
    !   which is set to "1" if you assume period boundary conditions
    !      IF(IPERIOD==1) THEN
    !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
    !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
    !      END IF
    !
    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    !
    do i=ia,iz
       do k=2,m1-1
          if(u(k,i,m3-1)<zr0) flux(k,i,m3-1)= &
               q0(k,i,m3)*u(k,i,m3-1)*dt*dd0(k,i,m3-1)
       end do
    end do

    do i=ia,iz
       do j=m3-1,2,-1 !jz,ja,-1
          do k=2,m1-1
             if(u(k,i,j-1)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j) + &
                     flux(k,i,j-1)/dxx(i,j))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k,i,j-1))/dxx(i,j)     ! Courant number
                x1n= (1.-x1)*(q0(k,i,j-1)-q0(k,i,j+1))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k,i,j+1)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i,j-1)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i,j-1)) ), max(q0(k,i,j),q0(k,i,j-1)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j)= max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j)-x1*cf1*dd0(k,i,j-1))/den1(k,i,j) ))
                flux(k,i,j-1)=dxx(i,j)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop x-z
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DY

  subroutine Advec3DZ(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,den0,&
       den1,dt,dxx,&
       dd0,&
       qn,mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: den0(m1,m2,m3)
    real   , intent(in) :: den1(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dxx(m1)
    real   , intent(in) :: dd0(m1,m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3DZ)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    ! copy input field to output field
    qn = q0
    imxmn=.false.


    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
    do j=ja,jz
       do i=ia,iz
          do k=2,m1-1 
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
             if(u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do


    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do j=ja,jz
       do i=ia,iz
          if(u(1,i,j)>=zr0) flux(1,i,j)= &
               q0(1,i,j)*u(1,i,j)*dt*dd0(1,i,j)
          do k=2,m1-1
             if(u(k,i  ,j)<zr0) cycle
             if(u(k-1,i,j)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(k)               ! Courant number
                x1n= (1.-x1)*(q0(k+1,i,j)-q0(k-1,i,j))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k-1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k+1,i,j)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k+1,i,j))  ), max(q0(k,i,j),q0(k+1,i,j)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j)= max(  vcmin(k,i,j),   min(   vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k-1,i,j)/dxx(k))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(k)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k-1,i,j)
             end if                                                  !Eq-9a
          end do
       end do
    end do

    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    do j=ja,jz
       do i=ia,iz
          if(u(m1-1,i,j)<zr0) flux(m1-1,i,j)=&
               q0(m1,i,j)*u(m1-1,i,j)*dt*dd0(m1-1,i,j)
          do k=m1-1,2,-1
             if(u(k-1,i,j)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(k) + flux(k-1,i,j)/dxx(k))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k-1,i,j))/dxx(k)     ! Courant number
                x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k+1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k-1,i,j)) ), max(q0(k,i,j),q0(k-1,i,j)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j) = max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(k)-x1*cf1*dd0(k-1,i,j))/den1(k,i,j) ))
                flux(k-1,i,j)=dxx(k)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop y-x
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DZ


  subroutine Advec3DZSedim(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,&
       dt,&
       dzt,rtgt,&
       qn,&
       mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dzt(m1)
    real   , intent(in) :: rtgt(m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3DZSedim)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    qn = q0
    imxmn=.false.

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
    do j=ja,jz
       do i=ia,iz
          do  k=2,m1-1 !
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(-u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
             if(-u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
             if(k==2) ck2= q0(k,i,j) !for sedim only
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    do j=ja,jz
       do i=ia,iz
          rtgti=1./rtgt(i,j)
          flux(m1-1,i,j)=q0(m1,i,j)*(-u(m1-1,i,j))*dt
          do k=m1-1,2,-1
             !srf       x1=  dt*ABS(u(k-1,i,j))/dxx(k)     ! Courant number
             x1=  dt*abs(u(k-1,i,j))*dzt(k)*rtgti     ! Courant number
             if(k==2) x1 = 0. ! no flux below sfc terrain,for sedim only
             x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
             cf= q0(k,i,j) + x1n                                       !Eq-4b
             if(imxmn(k+1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
             if(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
             cf1= min( max( cf, min(q0(k,i,j),q0(k-1,i,j)) ), max(q0(k,i,j),q0(k-1,i,j)) )
             if(k>2) then  !for sedim only
                qn(k,i,j) = max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                                !srf                 (q0(k,i,j)-flux(k,i,j)/dxx(k)      -x1*cf1) ))
                     (q0(k,i,j)-flux(k,i,j)*dzt(k)*rtgti-x1*cf1) ))
             else
                qn(k,i,j) = (q0(k,i,j)-flux(k,i,j)*dzt(k)*rtgti-x1*cf1)
             end if
             !srf	   flux(k-1,i,j)=dxx(k)             *(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
             flux(k-1,i,j)=(1./(dzt(k)*rtgti))*(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
          end do
       end do
    end do !- big loop y-x
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DZSedim

  subroutine Advec3DZSedimUpw(m1,m2,m3, ia,iz,ja,jz,u,dt,dzt,rtgt,qn,mynum)

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dzt(m1)
    real   , intent(in) :: rtgt(m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3DZSedimUpw)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    !- big loop y-x
    do j=ja,jz
       do i=ia,iz
          rtgti=1./rtgt(i,j)
          !srf dxx = dz = rtgti/dzt
          !srf qn(m1-1,i,j) = qn(m1-1,i,j) / (1.0 - dt*u(m1-1,i,j)/dxx(m1-1)      )
          qn(m1-1,i,j) = qn(m1-1,i,j) / (1.0 + dt*u(m1-1,i,j)*dzt(m1-1)*rtgti)
          do k=m1-2,2,-1 !
             !srf    qn(k,i,j)= 1.0/(1.0+dt*u(k,i,j)/dxx(k))&
             !srf               *( qn(k,i,j)+ dt*u(k,i,j) /dxx(k+1) * qn(k+1,i,j) )
             qn(k,i,j)= 1.0/(1.0 + dt*u(k,i,j)*dzt(k)*rtgti)&
                  *( qn(k,i,j) + dt*u(k+1,i,j)*dzt(k+1)*rtgti * qn(k+1,i,j) )
             !   tc(i,j,l,k) = 1.0/(1.0+dt_settl(k)*vd_cor/delz(i,j,l2))&
             !  	 *(tc(i,j,l,k) + dt_settl(k)*vd_cor /delz(i,j,l2-1) &
             !  	 * tc(i,j,l+1,k))
          end do
       end do
    end do !- big loop y-x
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DZSedimUpw


  subroutine InitializeDensities(mzp, mxp, myp, &
       mxpAdvMnt, mypAdvMnt, &
       iOffset, i1ExternAtAdvMnt, iMxpExternAtAdvMnt,  &
       jOffset, j1ExternAtAdvMnt, jMypExternAtAdvMnt,  &
       dn0, dn0u, dn0v, &
       dd0_3d ,dd0_3du,dd0_3dv,dd0_3dw )

    integer, intent(in) :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer, intent(in) :: mxp
    ! x dimension of external fields 
    integer, intent(in) :: myp
    ! y dimension of external fields 
    integer, intent(in) :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer, intent(in) :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer, intent(in) :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    real, pointer, intent(in) :: dn0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dn0u(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dn0v(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dd0_3d(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: dd0_3du(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: dd0_3dv(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: dd0_3dw(:,:,:)
    ! pointer intent(in), values intent(out)

    ! local var
    integer :: i
    integer :: j
    integer :: k
    integer :: iExtern
    integer :: jExtern

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(InitializeDensities)**"
    character(len=8) :: str(10)

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do


    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       jExtern = j + jOffset
       ! set Monotonic Advection west ghost zone fields to zero
       do i = 1, i1ExternAtAdvMnt-1
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
       ! fill where both Monotonic Advection and external fields
       ! are in range
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          do k = 1, mzp
             dd0_3d (k,i,j) = dn0 (k,iExtern,jExtern)
             dd0_3du(k,i,j) = dn0u(k,iExtern,jExtern)
             dd0_3dv(k,i,j) = dn0v(k,iExtern,jExtern)
          end do
          do k = 1,mzp-1
             dd0_3dw(k,i,j) = 0.5*&
                  (dn0(k,iExtern,jExtern) + dn0(k+1,iExtern,jExtern))
          end do
          dd0_3dw(mzp,i,j)=dd0_3dw(mzp-1,i,j)
       end do
       ! set Monotonic Advection east ghost zone fields to zero
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do


    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do
  end subroutine InitializeDensities
  



  subroutine PrepareTheorWinds(mzp, mxp, myp,&
       iOffset, i1ExternAtAdvMnt, iMxpExternAtAdvMnt,  &
       jOffset, j1ExternAtAdvMnt, jMypExternAtAdvMnt,  &
       dtlt, time,  &
       u3d, v3d, w3d, &
       dd0_3d, dd0_3du, dd0_3dv, dd0_3dw)

    integer, intent(in) :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer, intent(in) :: mxp
    ! x dimension of external fields 
    integer, intent(in) :: myp
    ! y dimension of external fields 
    integer, intent(in) :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    real, intent(in) :: dtlt
    real, intent(in) :: time
    real, pointer, intent(in) :: u3d(:,:,:)
    ! pointer intent(in), values intent(inout)
    real, pointer, intent(in) :: v3d(:,:,:)
    ! pointer intent(in), values intent(inout)
    real, pointer, intent(in) :: w3d(:,:,:)
    ! pointer intent(in), values intent(inout)
    real, pointer, intent(in) :: dd0_3d(:,:,:)
    ! pointer intent(in), values intent(inout)
    real, pointer, intent(in) :: dd0_3du(:,:,:)
    ! pointer intent(in), values intent(inout)
    real, pointer, intent(in) :: dd0_3dv(:,:,:)
    ! pointer intent(in), values intent(inout)
    real, pointer, intent(in) :: dd0_3dw(:,:,:)
    ! pointer intent(in), values intent(inout)

    !- local var
    integer :: i
    integer :: j
    integer :: k
    integer :: iExtern
    integer :: jExtern
    real :: dtlto2
    real :: ai0s  =  25.0
    real :: aj0s  =  50.0
    real :: umx   =   80.0
    real, parameter :: pii   =   3.141592653589793
    real    :: umax  =   0.0
    real    :: anrev,curnt,rx,xa,ilop,nrec,ya
    real    :: periodo  =   6.*3600.
!!$    real, parameter :: iwndty = 0  ! 0-rotating
!!$    real, parameter :: iwndty = 1  ! 1-divergent winds
    real, parameter :: iwndty = 2 

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(PrepareTheorWinds)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    dtlto2 =  10.!*dtlt


    if(iwndty==1) ai0s= 50.5
    ilop= ai0s-21.  ! needed for printouts
    ! Define wind fields (rotation or divergent) and initial mixing ratios
    ! Cone at (25,50) for rotating winds; Cone at (50,50) divergent winds
    do j = jMypExternAtAdvMnt, j1ExternAtAdvMnt, -1
       jExtern = j + jOffset
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          do k = 1, mzp
             dd0_3d (k,i,j)=1.
             dd0_3du(k,i,j)=1.
             dd0_3dv(k,i,j)=1.
             dd0_3dw(k,i,j)=1.

             u3d(k,i,j)=-2.*umx*(real(jExtern)-real(110)/2.-.5)/real(110)
             v3d(k,i,j)= 2.*umx*(real(iExtern)-real(100)/2.-.5)/real(100)
             w3d(k,i,j)= 0.
          end do
       end do
    end do
             
    if(iwndty==1) then
       do j = jMypExternAtAdvMnt, j1ExternAtAdvMnt, -1
          jExtern = j + jOffset
          do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
             iExtern = i + iOffset
             do k = 1, mzp
                xa=pii/25.
                u3d(k,i,j)=umx*&
                     sin(xa*real(iExtern))*&
                     sin(xa*(real(jExtern)))
                v3d(k,i,j)=umx*&
                     cos(xa*(real(iExtern)-.5))*&
                     cos(xa*(real(jExtern)+.5))
             end do
          end do
       end do

    else if (iwndty==2) then
       do j = jMypExternAtAdvMnt, j1ExternAtAdvMnt, -1
          jExtern = j + jOffset
          do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
             iExtern = i + iOffset
             do k = 1, mzp
                xa=pii/100. ! myp=mxp
                u3d(k,i,j)=umx*&
                     (sin(xa*real(iExtern)))**2*&
                     sin(2*xa*(real(jExtern)))*&
                     cos(pii*time/periodo)
                v3d(k,i,j)=-umx*&
                     (sin(xa*real(jExtern)))**2*&
                     sin(2*xa*(real(iExtern)))*&
                     cos(pii*time/periodo)
             end do
          end do
       end do

    else if (iwndty==3) then
       do j = jMypExternAtAdvMnt, j1ExternAtAdvMnt, -1
          jExtern = j + jOffset
          do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
             iExtern = i + iOffset
             do k = 1, mzp
                xa=pii/100. ! myp=mxp
                ya=50.
                u3d(k,i,j)=-umx*&
                     (sin(xa*real(iExtern)))**2*&
                     sin(2.*xa*(real(jExtern)-ya))*&
                     cos(pii*time/periodo)
                v3d(k,i,j)=0.5*&
                     umx*&
                     sin(2.*xa*real(iExtern))*&
                     cos(xa*(real(jExtern)-ya))*&
                     cos(pii*time/periodo)
             end do
          end do
       end do
    end if
  end subroutine PrepareTheorWinds


  
  subroutine advmnt_driver(oneGrid, varn, &
       m1 ,m2 ,m3 ,ia,iz,ja,jz,izu,jzv,&
       i0,j0,nodemyp,nodemxp,nodemzp,mynum)

    type(Grid), pointer, intent(in) :: oneGrid
    integer , intent(in) :: m1
    integer , intent(in) :: m2
    integer , intent(in) :: m3
    integer , intent(in) :: ia
    integer , intent(in) :: iz
    integer , intent(in) :: ja
    integer , intent(in) :: jz
    integer , intent(in) :: izu
    integer , intent(in) :: jzv
    integer , intent(in) :: i0
    integer , intent(in) :: j0
    integer, intent(in) :: nodemxp(:,:)
    integer, intent(in) :: nodemyp(:,:)
    integer, intent(in) :: nodemzp(:,:)
    integer , intent(in) :: mynum
    character(len=*),intent(in) :: varn

    !--- local vars
    integer :: n
    integer :: ng
    integer :: mxyzp
    integer :: i
    integer :: j
    integer :: k
    integer :: iExtern
    integer :: jExtern
    integer :: ierr
    integer :: procfile
    integer :: ibegin
    integer :: iend
    integer :: jbegin
    integer :: jend
    integer :: i_scl
    integer :: sori
    integer :: sorj
    integer :: sosi
    integer :: sosj
    integer :: current_aer_ispc
    integer :: current_ndt_z
    integer, target :: ndt_z(naer_transported)
    integer, target :: ndtZ(naer_transported)
    integer, pointer :: p1(:) => null()
    integer, pointer :: p2(:) => null()
    real, pointer :: scalarp
    real, pointer :: scalart
    logical  :: IsThisScalarAer =.false.

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(advmnt_driver)**"
    character(len=8) :: str(11)

    type(MonotonicAdvection), pointer :: oneAdvMnt

    integer :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer :: mxp
    ! x dimension of external fields 
    integer :: myp
    ! y dimension of external fields 
    integer :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection

    ! dimension of external fields (regular ghost zone width)

    mzp=oneGrid%NodeDims%mzp
    mxp=oneGrid%NodeDims%mxp
    myp=oneGrid%NodeDims%myp

    ! dimension of Monotonic Advection fields (wide ghost zone width)

    mxpAdvMnt=oneGrid%NodeDimsAdvMnt%mxp
    mypAdvMnt=oneGrid%NodeDimsAdvMnt%myp

    ! index external = index Monotonic Advection + offset

    iOffset = oneGrid%NodeDimsAdvMnt%i0 - oneGrid%NodeDims%i0 
    i1ExternAtAdvMnt = 1 - iOffset
    iMxpExternAtAdvMnt = mxp - iOffset
    jOffset = oneGrid%NodeDimsAdvMnt%j0 - oneGrid%NodeDims%j0 
    j1ExternAtAdvMnt = 1 - jOffset
    jMypExternAtAdvMnt = myp - jOffset

    ! current grid
    ! necessary while there are global variables outside oneGrid
    
    ng = OneGrid%Id
    
    if (dumpLocal) then
       write(str(1),"(i8)") mzp
       call MsgDump(h//"mzp="//trim(adjustl(str(1))))
       write(str(1),"(i8)") mxp
       call MsgDump(h//"mxp="//trim(adjustl(str(1))))
       write(str(1),"(i8)") myp
       call MsgDump(h//"myp="//trim(adjustl(str(1))))
       write(str(1),"(i8)") mxpAdvMnt
       call MsgDump(h//"mxpAdvMnt="//trim(adjustl(str(1))))
       write(str(1),"(i8)") mypAdvMnt
       call MsgDump(h//"mypAdvMnt="//trim(adjustl(str(1))))
       write(str(1),"(i8)") iOffset
       call MsgDump(h//"iOffset="//trim(adjustl(str(1))))
       write(str(1),"(i8)") i1ExternAtAdvMnt
       call MsgDump(h//"i1ExternAtAdvMnt="//trim(adjustl(str(1))))
       write(str(1),"(i8)") iMxpExternAtAdvMnt
       call MsgDump(h//"iMxpExternAtAdvMnt="//trim(adjustl(str(1))))
       write(str(1),"(i8)") jOffset
       call MsgDump(h//"jOffset="//trim(adjustl(str(1))))
       write(str(1),"(i8)") j1ExternAtAdvMnt
       call MsgDump(h//"j1ExternAtAdvMnt="//trim(adjustl(str(1))))
       write(str(1),"(i8)") jMypExternAtAdvMnt
       call MsgDump(h//"jMypExternAtAdvMnt="//trim(adjustl(str(1))))
    end if

!!$    if(mnt_adv_jnitialized == OFF) then
    if (mynum == 0) stop 'ADV MNT called with mynum = 0, try np = 2'
    if (dumpLocal) then
       call MsgDump(h//" enter initialization")
    end if
    
    if(use_true_density == OFF) then
       if (dumpLocal) then
          call MsgDump(h//" invokes InitializeDensities")
       end if
       call InitializeDensities(mzp, mxp, myp, &
            mxpAdvMnt, mypAdvMnt, &
            iOffset, i1ExternAtAdvMnt, iMxpExternAtAdvMnt,  &
            jOffset, j1ExternAtAdvMnt, jMypExternAtAdvMnt,  &
            basic_g(ng)%dn0, basic_g(ng)%dn0u, basic_g(ng)%dn0v, &
            oneAdvMnt%dd0_3d, oneAdvMnt%dd0_3du, &
            oneAdvMnt%dd0_3dv, oneAdvMnt%dd0_3dw)
    end if

!!$       mnt_adv_jnitialized= ON
!!$   end if

    oneAdvMnt => CreateMonotonicAdvection(oneGrid)

    call UpdateFieldAdressAtAdvMnt(&
         oneGrid%AdvMntUVSendX, oneGrid%AdvMntUVRecvX, &
         oneGrid%AdvMntUVSendY, oneGrid%AdvMntUVRecvY, &
         oneGrid%AdvMntDxDySendX, oneGrid%AdvMntDxDyRecvX, &
         oneGrid%AdvMntDxDySendY, oneGrid%AdvMntDxDyRecvY, &
         oneGrid%AdvMntDd0SendX, oneGrid%AdvMntDd0RecvX, &
         oneGrid%AdvMntDd0SendY, oneGrid%AdvMntDd0RecvY, &
         oneGrid%AdvMntDenSendX, oneGrid%AdvMntDenRecvX, &
         oneGrid%AdvMntDenSendY, oneGrid%AdvMntDenRecvY, &
         oneGrid%AdvMntScaSendX, oneGrid%AdvMntScaRecvX, &
         oneGrid%AdvMntScaSendY, oneGrid%AdvMntScaRecvY, &
         oneAdvMnt%u3d, oneAdvMnt%v3d, &
         oneAdvMnt%dxtW, oneAdvMnt%dytW, &
         oneAdvMnt%dd0_3d, oneAdvMnt%dd0_3du, &
         oneAdvMnt%dd0_3dv, oneAdvMnt%dd0_3dw, &
         oneAdvMnt%den0_3d, oneAdvMnt%den1_3d, &
         oneAdvMnt%den2_3d, oneAdvMnt%den3_3d, &
         oneAdvMnt%vc3d_in, oneAdvMnt%vc3d_out)

    if (dumpLocal) then
       call MsgDump(h//" grid after UpdateFieldAdressAtAdvMnt")
       call DumpGrid(oneGrid)
    end if

    call InitializeGridSpacings(&
         mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
         iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
         jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
         grid_g(ng)%dxt, &
         grid_g(ng)%dyt, &
         grid_g(ng)%fmapt, &
         grid_g(ng)%rtgt, &
         dztn(:,ng), &
         oneAdvMnt%dxtW, &
         oneAdvMnt%dytW, &
         oneAdvMnt%dztW)

!!$    call Compare(oneAdvMnt%dxtW, advmnt_g(ng)%dxtW, "dxtW", .true.)
!!$    call Compare(oneAdvMnt%dytW, advmnt_g(ng)%dytW, "dytW", .true.)
!!$    call Compare(oneAdvMnt%dztW, advmnt_g(ng)%dztW, "dztW", .true.)

    !- This scheme is not applied to advect  U, V, and W
    if (varn .eq. 'V' .or. varn .eq. 'ALL') then
       stop 'not using mnt to advect u,v,w'
    end if
    if (if_adap /= 0) then
       stop 'MNT advection not ready for shaved eta'
    end if

    !- get actual air densities, if using them instead of basic state fields
    if(use_true_density == ON) then
       call GetTrueDensities(&
            mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
            iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
            jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
            level,&
            basic_g(ng)%rtp, &
            basic_g(ng)%rv, &
            basic_g(ng)%pp, &
            basic_g(ng)%pi0, &
            basic_g(ng)%theta, &
            oneAdvMnt%dd0_3d, &
            oneAdvMnt%dd0_3du, &
            oneAdvMnt%dd0_3dv, &
            oneAdvMnt%dd0_3dw)

!!$       call Compare(oneAdvMnt%dd0_3d, advmnt_g(ng)%dd0_3d, "dd0_3d", .true.)
!!$       call Compare(oneAdvMnt%dd0_3du, advmnt_g(ng)%dd0_3du, "dd0_3du", .true.)
!!$       call Compare(oneAdvMnt%dd0_3dv, advmnt_g(ng)%dd0_3dv, "dd0_3dv", .true.)
!!$       call Compare(oneAdvMnt%dd0_3dw, advmnt_g(ng)%dd0_3dw, "dd0_3dw", .true.)
    end if

    !- prepare wind velocities including map factors

    ndtZ=0
    call PrepareWinds(&
         ng, mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
         iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
         jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
         dtlt, &
         basic_g(ng)%uc, basic_g(ng)%up, &
         basic_g(ng)%vc, basic_g(ng)%vp, &
         basic_g(ng)%wc, basic_g(ng)%wp, &
         grid_g(ng)%fmapui, grid_g(ng)%fmapvi, &
         grid_g(ng)%rtgt, grid_g(ng)%rtgu, grid_g(ng)%rtgv, &
         grid_g(ng)%f13t, grid_g(ng)%f23t, &
         oneAdvMnt%u3d, oneAdvMnt%v3d, oneAdvMnt%w3d, &
         aerosol, naer_transported, &
         dd_sedim, dzt, ndtZ)

    if(theor_wind == on) then
       if (dumpLocal) then
          call MsgDump (h//" invokes PrepareTheorWinds")
       end if

       call PrepareTheorWinds(mzp, mxp, myp,&
            iOffset, i1ExternAtAdvMnt, iMxpExternAtAdvMnt,  &
            jOffset, j1ExternAtAdvMnt, jMypExternAtAdvMnt,  &
            dtlt, time,  &
            oneAdvMnt%u3d, oneAdvMnt%v3d, oneAdvMnt%w3d, &
            oneAdvMnt%dd0_3d, oneAdvMnt%dd0_3du, &
            oneAdvMnt%dd0_3dv, oneAdvMnt%dd0_3dw)
    end if

    !- prepare Walcek's air densities

    call GetWalceksDensities(&
         mzp, dtlt, mxpAdvMnt, mypAdvMnt, &
         i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
         j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
         oneAdvMnt%u3d, oneAdvMnt%v3d, oneAdvMnt%w3d, &
         oneAdvMnt%dd0_3d, oneAdvMnt%dd0_3du, &
         oneAdvMnt%dd0_3dv, oneAdvMnt%dd0_3dw, &
         oneAdvMnt%dxtW, oneAdvMnt%dytW, oneAdvMnt%dztW, &
         oneAdvMnt%den0_3d, oneAdvMnt%den1_3d, &
         oneAdvMnt%den2_3d, oneAdvMnt%den3_3d)

!!$    call Compare(oneAdvMnt%den0_3d, advmnt_g(ng)%den0_3d, "den0_3d", .true.)
!!$    call Compare(oneAdvMnt%den1_3d, advmnt_g(ng)%den1_3d, "den1_3d", .true.)
!!$    call Compare(oneAdvMnt%den2_3d, advmnt_g(ng)%den2_3d, "den2_3d", .true.)
!!$    call Compare(oneAdvMnt%den3_3d, advmnt_g(ng)%den3_3d, "den3_3d", .true.)

    call PostSendRecvMsgs(oneGrid%AdvMntUVSendX, oneGrid%AdvMntUVRecvX)
    call PostSendRecvMsgs(oneGrid%AdvMntDxDySendX, oneGrid%AdvMntDxDyRecvX)
    call PostSendRecvMsgs(oneGrid%AdvMntDd0SendX, oneGrid%AdvMntDd0RecvX)
    call PostSendRecvMsgs(oneGrid%AdvMntDenSendX, oneGrid%AdvMntDenRecvX)
    call WaitSendRecvMsgs(oneGrid%AdvMntUVSendX, oneGrid%AdvMntUVRecvX)
    call WaitSendRecvMsgs(oneGrid%AdvMntDxDySendX, oneGrid%AdvMntDxDyRecvX)
    call WaitSendRecvMsgs(oneGrid%AdvMntDd0SendX, oneGrid%AdvMntDd0RecvX)
    call WaitSendRecvMsgs(oneGrid%AdvMntDenSendX, oneGrid%AdvMntDenRecvX)

    oneAdvMnt%l_dxtW=0.0
    oneAdvMnt%l_dytW=0.0

!!$    call Compare(oneAdvMnt%dd0_3d, advmnt_g(ng)%dd0_3d, "dd0_3d", .true.)
!!$    call Compare(oneAdvMnt%dd0_3du, advmnt_g(ng)%dd0_3du, "dd0_3du", .true.)
!!$    call Compare(oneAdvMnt%dd0_3dv, advmnt_g(ng)%dd0_3dv, "dd0_3dv", .true.)
!!$    call Compare(oneAdvMnt%dd0_3dw, advmnt_g(ng)%dd0_3dw, "dd0_3dw", .true.)
!!$    call Compare(oneAdvMnt%den0_3d, advmnt_g(ng)%den0_3d, "den0_3d", .true.)
!!$    call Compare(oneAdvMnt%den1_3d, advmnt_g(ng)%den1_3d, "den1_3d", .true.)
!!$    call Compare(oneAdvMnt%den2_3d, advmnt_g(ng)%den2_3d, "den2_3d", .true.)
!!$    call Compare(oneAdvMnt%den3_3d, advmnt_g(ng)%den3_3d, "den3_3d", .true.)
!!$    call Compare(oneAdvMnt%l_dxtW, advmnt_g(ng)%l_dxtW, "l_dxtW", .true.)
!!$    call Compare(oneAdvMnt%l_dytW, advmnt_g(ng)%l_dytW, "l_dytW", .true.)
!!$    call Compare(oneAdvMnt%dxtW, advmnt_g(ng)%dxtW, "dxtW", .true.)
!!$    call Compare(oneAdvMnt%dytW, advmnt_g(ng)%dytW, "dytW", .true.)

    if (dumpLocal) then
       call MsgDump (h//" invokes initial_fields_update exchanging borders on y")
    end if

    call PostSendRecvMsgs(oneGrid%AdvMntUVSendY, oneGrid%AdvMntUVRecvY)
    call PostSendRecvMsgs(oneGrid%AdvMntDxDySendY, oneGrid%AdvMntDxDyRecvY)
    call PostSendRecvMsgs(oneGrid%AdvMntDd0SendY, oneGrid%AdvMntDd0RecvY)
    call PostSendRecvMsgs(oneGrid%AdvMntDenSendY, oneGrid%AdvMntDenRecvY)
    call WaitSendRecvMsgs(oneGrid%AdvMntUVSendY, oneGrid%AdvMntUVRecvY)
    call WaitSendRecvMsgs(oneGrid%AdvMntDxDySendY, oneGrid%AdvMntDxDyRecvY)
    call WaitSendRecvMsgs(oneGrid%AdvMntDd0SendY, oneGrid%AdvMntDd0RecvY)
    call WaitSendRecvMsgs(oneGrid%AdvMntDenSendY, oneGrid%AdvMntDenRecvY)

!!$    call MsgDump(h//" depois de InitialFieldsUpdate on y:")
!!$    call Compare(oneAdvMnt%u3d, advmnt_g(ng)%u3d, "u3d", .true.)
!!$    call Compare(oneAdvMnt%v3d, advmnt_g(ng)%v3d, "v3d", .true.)
!!$    call Compare(oneAdvMnt%w3d, advmnt_g(ng)%w3d, "w3d", .true.)
!!$    call Compare(oneAdvMnt%dxtW, advmnt_g(ng)%dxtW, "dxtW", .true.)
!!$    call Compare(oneAdvMnt%dytW, advmnt_g(ng)%dytW, "dytW", .true.)
!!$    call Compare(oneAdvMnt%dd0_3d, advmnt_g(ng)%dd0_3d, "dd0_3d", .true.)
!!$    call Compare(oneAdvMnt%dd0_3du, advmnt_g(ng)%dd0_3du, "dd0_3du", .true.)
!!$    call Compare(oneAdvMnt%dd0_3dv, advmnt_g(ng)%dd0_3dv, "dd0_3dv", .true.)
!!$    call Compare(oneAdvMnt%dd0_3dw, advmnt_g(ng)%dd0_3dw, "dd0_3dw", .true.)
!!$    call Compare(oneAdvMnt%dd0_3d, advmnt_g(ng)%dd0_3d, "dd0_3d", .true.)
!!$    call Compare(oneAdvMnt%dd0_3du, advmnt_g(ng)%dd0_3du, "dd0_3du", .true.)
!!$    call Compare(oneAdvMnt%dd0_3dv, advmnt_g(ng)%dd0_3dv, "dd0_3dv", .true.)
!!$    call Compare(oneAdvMnt%dd0_3dw, advmnt_g(ng)%dd0_3dw, "dd0_3dw", .true.)
!!$    call Compare(oneAdvMnt%den0_3d, advmnt_g(ng)%den0_3d, "den0_3d", .true.)
!!$    call Compare(oneAdvMnt%den1_3d, advmnt_g(ng)%den1_3d, "den1_3d", .true.)
!!$    call Compare(oneAdvMnt%den2_3d, advmnt_g(ng)%den2_3d, "den2_3d", .true.)
!!$    call Compare(oneAdvMnt%den3_3d, advmnt_g(ng)%den3_3d, "den3_3d", .true.)

    !- ready to do advection, loop over all scalars
    if(advmnt == 1) then
       i_scl=1                                            !- all scalars
    else if(advmnt == 2) then
       i_scl=num_scalar(ng) - NSPECIES_TRANSPORTED +1  !- only chemical + aer species
    else if(advmnt == 3) then
       i_scl=2                                            !- all scalars, but not theta_il
    end if

    call DeepCopyToScalarTabAtOneGrid(oneGrid)
    
    if (dumpLocal) then
       write(str(1),"(i8)") num_scalar(ng)
       call MsgDump(h//" there are "//trim(adjustl(str(1)))//" scalar fields")
    end if

    !srf- do n=1,num_scalar(ng)     ! original
    do n=i_scl,num_scalar(ng)

       !- if RK or ABM3 scheme, THP/THC are not transported here

       if (dyncore_flag == 2) then
          if (oneGrid%ScalarTab(n)%name == 'THC' .or. &
               oneGrid%ScalarTab(n)%name == 'THP') cycle
       end if

       !srf - somente para gases e aerossois
       !     do n=num_scalar(ng) - NSPECIES_TRANSPORTED +1,num_scalar(ng)
       !      if (scalar_tab(n,ng)%name /= 'COP' .and. scalar_tab(n,ng)%name /= 'CH4P') cycle
       !          scalar_tab(n,ng)%name /= 'O3P'  ) cycle

       !- Aerosol sedimentation
       IsThisScalarAer  = .false.
       current_aer_ispc = 0
       current_ndt_z    = 1
       if(ccatt == 1 .and. aerosol > 0 .and. n >= num_scalar_aer_1st) then
          !srf-  We are going to include sedimentation of aerosols at
          !      vertical advection tendency. It is supposed that scalars
          !      with  N >= num_scalar_aer_1st are _all_ aerosols .
          !
          IsThisScalarAer=.true.
          current_aer_ispc = n - num_scalar_aer_1st + 1
          current_ndt_z    = ndt_z (current_aer_ispc)

       end if

       if (associated(oneGrid%ScalarTab(n)%var_p_3D)) then

          ! set oneAdvMnt%vc3d_in north border to zero
          do j = 1, j1ExternAtAdvMnt-1
             do i = 1, mxpAdvMnt
                do k = 1, mzp
                   oneAdvMnt%vc3d_in(k,i,j) = 0.0
                end do
             end do
          end do
          do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
             jExtern = j + jOffset
             ! set oneAdvMnt%vc3d_in west border to zero
             do i = 1, i1ExternAtAdvMnt-1
                do k = 1, mzp
                   oneAdvMnt%vc3d_in(k,i,j) = 0.0
                end do
             end do
             ! copy scalartab external field to the
             ! inner part of oneAdvMnt%vc3d_in
             do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
                iExtern = i + iOffset
                do k = 1, mzp
                   oneAdvMnt%vc3d_in(k,i,j) = oneGrid%ScalarTab(n)%var_p_3D(k,iExtern,jExtern)
                end do
             end do
             ! set oneAdvMnt%vc3d_in east border to zero
             do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
                do k = 1, mzp
                   oneAdvMnt%vc3d_in(k,i,j) = 0.0
                end do
             end do
          end do
          ! set oneAdvMnt%vc3d_in south border to zero
          do j = jMypExternAtAdvMnt+1, mypAdvMnt
             do i = 1, mxpAdvMnt
                do k = 1, mzp
                   oneAdvMnt%vc3d_in(k,i,j) = 0.0
                end do
             end do
          end do

          call AdvectMnt(oneAdvMnt, oneGrid, ng,m1,m2,m3,ia,iz,ja,jz,dtlt,mynum,n, &
               current_aer_ispc,current_ndt_z,IsThisScalarAer)


!!$          call Compare(advmnt_g(ng)%vc3d_in, oneAdvMnt%vc3d_in, &
!!$               "vc3d_in", .true.)

          call AdvectTendency(mzp, mxp, &
               iOffset, jOffset, &
               ia, iz, ja, jz, dtlt, &
               scalarp3D=oneGrid%ScalarTab(n)%var_p_3D, &
               AdvMntField=oneAdvMnt%vc3d_out, &
               scalart1D=oneGrid%ScalarTab(n)%var_t_1D)

!!$          call Compare(scalar_tab(n,ng)%var_t_1D, &
!!$               oneGrid%ScalarTab(n)%var_t_1D, &
!!$               "var_t_1D", .true.)
       else
          if (dumpLocal) then
             call MsgDump (h//" not associated "//&
                  trim(adjustl(oneGrid%ScalarTab(n)%name)))
          end if
       end if


    end do

    call DeepCopyFromScalarTabAtOneGrid(oneGrid)
    call DestroyMonotonicAdvection(oneAdvMnt)
    
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine advmnt_driver




















  subroutine Advec3d_X(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,den0,&
       den1,dt,dxx,&
       dd0,&
       qn,mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: den0(m1,m2,m3)
    real   , intent(in) :: den1(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dxx(m2,m3)
    real   , intent(in) :: dd0(m1,m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    integer :: ii
    integer :: ji
    integer :: ii0
    integer :: ji0
    integer :: ie
    integer :: je
    integer :: ie0
    integer :: je0
    integer :: ipos
    integer :: iia
    integer :: iiz
    integer :: nvar
    integer :: nf
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3d_X)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    ! copy input field to output field
    qn=q0
    imxmn=.false.

    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do j=ja,jz
       do k=2,m1-1
          if(u(k,1,j)>=zr0) flux(k,1,j)= q0(k,1,j)*u(k,1,j)*dt*dd0(k,1,j)
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
    do j=ja,jz
       do  i=2,m2-1 ! ia,iz-1 or 1,iz-1
          do k=2,m1-1
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k,i-1,j),q0(k,i+1,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k,i-1,j),q0(k,i+1,j))+eps)        !       extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k,i,j  )< zr0) ck1= q0(k,i+1,j)
             if(u(k,i-1,j)>=zr0) ck2= q0(k,i-1,j)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    do j=ja,jz
       do  i=2,m2-1 ! ia,iz-1 or 1,iz-1
          do k=2,m1-1
             if(u(k,i,j)<zr0) cycle
             if(u(k,i-1,j)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(i,j)               ! Courant number
                x1n= (1.-x1)*(q0(k,i+1,j)-q0(k,i-1,j))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k,i-1,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i+1,j)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i+1,j))  ), max(q0(k,i,j),q0(k,i+1,j)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j) = max(vcmin(k,i,j),min(vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k,i-1,j)/dxx(i,j))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(i,j)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k,i-1,j)
             end if                                                  !Eq-9a
          end do
       end do
    end do

    ! If periodic boundary conditions are assumed, it is necessary
    !   to recalculate the updated mixing ratio at cell 1 if there
    !   is inflow to that cell from the boundary between IDIM and 1
    !   Here these statements are commented out, but should be uncommented
    !   if this subroutine is needed for periodic boundary conditions,
    !   and then one of the calling arguements to the subroutine is IPERIOD
    !   which is set to "1" if you assume period boundary conditions
    !      IF(IPERIOD==1) THEN
    !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
    !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
    !      END IF
    !
    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    !
    do j=ja,jz
       do k=2,m1-1
          if(u(k,m2-1,j)<zr0) flux(k,m2-1,j)= &
               q0(k,m2,j)*u(k,m2-1,j)*dt*dd0(k,m2-1,j)
       end do
    end do

    do j=ja,jz
       do i=m2-1,2,-1 !iz,ia,-1
          do k=2,m1-1
             if(u(k,i-1,j)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j) + &
                     flux(k,i-1,j)/dxx(i,j))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k,i-1,j))/dxx(i,j)     ! Courant number
                x1n= (1.-x1)*(q0(k,i-1,j)-q0(k,i+1,j))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k,i+1,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i-1,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i-1,j)) ), max(q0(k,i,j),q0(k,i-1,j)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j)= max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j)-x1*cf1*dd0(k,i-1,j))/den1(k,i,j) ))
                flux(k,i-1,j)=dxx(i,j)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop y-z
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3d_X





  subroutine Advec3d_Y(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,den0,&
       den1,dt,dxx,&
       dd0,&
       qn,mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: den0(m1,m2,m3)
    real   , intent(in) :: den1(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dxx(m2,m3)
    real   , intent(in) :: dd0(m1,m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    integer :: ii
    integer :: ji
    integer :: ii0
    integer :: ji0
    integer :: ie
    integer :: je
    integer :: ie0
    integer :: je0
    integer :: ipos
    integer :: iia
    integer :: iiz
    integer :: nvar
    integer :: nf
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3d_Y)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    ! copy input field to output field
    qn= q0
    imxmn=.false.

    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do i=ia,iz
       do k=2,m1-1
          if(u(k,i,1)>=zr0) flux(k,i,1)= q0(k,i,1)*u(k,i,1)*dt*dd0(k,i,1)
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !	mixing ratio at t+dt. If these limits are ever violated,
    !	non-monotonic (oscillatory) behavior in solution results
    do i=ia,iz
       do  j=2,m3-1 ! ja,jz
          do k=2,m1-1
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k,i,j-1),q0(k,i,j+1))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k,i,j-1),q0(k,i,j+1))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k,i,j  )< zr0) ck1= q0(k,i,j+1)
             if(u(k,i,j-1)>=zr0) ck2= q0(k,i,j-1)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    do i=ia,iz
       do  j=2,m3-1 ! ja,jz
          do k=2,m1-1
             if(u(k,i,j)<zr0) cycle
             if(u(k,i,j-1)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(i,j)               ! Courant number
                x1n= (1.-x1)*(q0(k,i,j+1)-q0(k,i,j-1))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k,i,j-1)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i,j+1)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i,j+1))  ), max(q0(k,i,j),q0(k,i,j+1)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j)= max(  vcmin(k,i,j),   min(   vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k,i,j-1)/dxx(i,j))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(i,j)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k,i,j-1) !Eq-9a
             end if
          end do
       end do
    end do

    ! If periodic boundary conditions are assumed, it is necessary
    !   to recalculate the updated mixing ratio at cell 1 if there
    !   is inflow to that cell from the boundary between IDIM and 1
    !   Here these statements are commented out, but should be uncommented
    !   if this subroutine is needed for periodic boundary conditions,
    !   and then one of the calling arguements to the subroutine is IPERIOD
    !   which is set to "1" if you assume period boundary conditions
    !      IF(IPERIOD==1) THEN
    !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
    !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
    !      END IF
    !
    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    !
    do i=ia,iz
       do k=2,m1-1
          if(u(k,i,m3-1)<zr0) flux(k,i,m3-1)= &
               q0(k,i,m3)*u(k,i,m3-1)*dt*dd0(k,i,m3-1)
       end do
    end do

    do i=ia,iz
       do j=m3-1,2,-1 !jz,ja,-1
          do k=2,m1-1
             if(u(k,i,j-1)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j) + &
                     flux(k,i,j-1)/dxx(i,j))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k,i,j-1))/dxx(i,j)     ! Courant number
                x1n= (1.-x1)*(q0(k,i,j-1)-q0(k,i,j+1))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k,i,j+1)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i,j-1)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i,j-1)) ), max(q0(k,i,j),q0(k,i,j-1)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j)= max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j)-x1*cf1*dd0(k,i,j-1))/den1(k,i,j) ))
                flux(k,i,j-1)=dxx(i,j)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop x-z
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3d_Y





  subroutine Advec3d_Z(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,den0,&
       den1,dt,dxx,&
       dd0,&
       qn,mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: den0(m1,m2,m3)
    real   , intent(in) :: den1(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dxx(m1)
    real   , intent(in) :: dd0(m1,m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3d_Z)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    ! copy input field to output field
    qn = q0
    imxmn=.false.


    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
    do j=ja,jz
       do i=ia,iz
          do k=2,m1-1 
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
             if(u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do


    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do j=ja,jz
       do i=ia,iz
          if(u(1,i,j)>=zr0) flux(1,i,j)= &
               q0(1,i,j)*u(1,i,j)*dt*dd0(1,i,j)
          do k=2,m1-1
             if(u(k,i  ,j)<zr0) cycle
             if(u(k-1,i,j)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(k)               ! Courant number
                x1n= (1.-x1)*(q0(k+1,i,j)-q0(k-1,i,j))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k-1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k+1,i,j)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k+1,i,j))  ), max(q0(k,i,j),q0(k+1,i,j)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j)= max(  vcmin(k,i,j),   min(   vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k-1,i,j)/dxx(k))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(k)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k-1,i,j)
             end if                                                  !Eq-9a
          end do
       end do
    end do

    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    do j=ja,jz
       do i=ia,iz
          if(u(m1-1,i,j)<zr0) flux(m1-1,i,j)=&
               q0(m1,i,j)*u(m1-1,i,j)*dt*dd0(m1-1,i,j)
          do k=m1-1,2,-1
             if(u(k-1,i,j)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(k) + flux(k-1,i,j)/dxx(k))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k-1,i,j))/dxx(k)     ! Courant number
                x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k+1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k-1,i,j)) ), max(q0(k,i,j),q0(k-1,i,j)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j) = max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(k)-x1*cf1*dd0(k-1,i,j))/den1(k,i,j) ))
                flux(k-1,i,j)=dxx(k)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop y-x
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3d_Z





  subroutine Advec3d_Z_sedim(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,&
       dt,&
       dzt,rtgt,&
       qn,&
       mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dzt(m1)
    real   , intent(in) :: rtgt(m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3d_Z_sedim)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    qn = q0
    imxmn=.false.

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
    do j=ja,jz
       do i=ia,iz
          do  k=2,m1-1 !
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(-u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
             if(-u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
             if(k==2) ck2= q0(k,i,j) !for sedim only
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    do j=ja,jz
       do i=ia,iz
          rtgti=1./rtgt(i,j)
          flux(m1-1,i,j)=q0(m1,i,j)*(-u(m1-1,i,j))*dt
          do k=m1-1,2,-1
             !srf       x1=  dt*ABS(u(k-1,i,j))/dxx(k)     ! Courant number
             x1=  dt*abs(u(k-1,i,j))*dzt(k)*rtgti     ! Courant number
             if(k==2) x1 = 0. ! no flux below sfc terrain,for sedim only
             x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
             cf= q0(k,i,j) + x1n                                       !Eq-4b
             if(imxmn(k+1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
             if(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
             cf1= min( max( cf, min(q0(k,i,j),q0(k-1,i,j)) ), max(q0(k,i,j),q0(k-1,i,j)) )
             if(k>2) then  !for sedim only
                qn(k,i,j) = max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                                !srf                 (q0(k,i,j)-flux(k,i,j)/dxx(k)      -x1*cf1) ))
                     (q0(k,i,j)-flux(k,i,j)*dzt(k)*rtgti-x1*cf1) ))
             else
                qn(k,i,j) = (q0(k,i,j)-flux(k,i,j)*dzt(k)*rtgti-x1*cf1)
             end if
             !srf	   flux(k-1,i,j)=dxx(k)             *(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
             flux(k-1,i,j)=(1./(dzt(k)*rtgti))*(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
          end do
       end do
    end do !- big loop y-x
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3d_Z_sedim





  subroutine Advec3d_Z_sedim_upw(m1,m2,m3, ia,iz,ja,jz,u,dt,dzt,rtgt,qn,mynum)

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dzt(m1)
    real   , intent(in) :: rtgt(m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3d_Z_sedim_upw)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    !- big loop y-x
    do j=ja,jz
       do i=ia,iz
          rtgti=1./rtgt(i,j)
          !srf dxx = dz = rtgti/dzt
          !srf qn(m1-1,i,j) = qn(m1-1,i,j) / (1.0 - dt*u(m1-1,i,j)/dxx(m1-1)      )
          qn(m1-1,i,j) = qn(m1-1,i,j) / (1.0 + dt*u(m1-1,i,j)*dzt(m1-1)*rtgti)
          do k=m1-2,2,-1 !
             !srf    qn(k,i,j)= 1.0/(1.0+dt*u(k,i,j)/dxx(k))&
             !srf               *( qn(k,i,j)+ dt*u(k,i,j) /dxx(k+1) * qn(k+1,i,j) )
             qn(k,i,j)= 1.0/(1.0 + dt*u(k,i,j)*dzt(k)*rtgti)&
                  *( qn(k,i,j) + dt*u(k+1,i,j)*dzt(k+1)*rtgti * qn(k+1,i,j) )
             !   tc(i,j,l,k) = 1.0/(1.0+dt_settl(k)*vd_cor/delz(i,j,l2))&
             !  	 *(tc(i,j,l,k) + dt_settl(k)*vd_cor /delz(i,j,l2-1) &
             !  	 * tc(i,j,l+1,k))
          end do
       end do
    end do !- big loop y-x
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3d_Z_sedim_upw



  subroutine StoreNamelistFileAtRadvc_mnt(oneNamelistFile)

    ! import NameListFile values into module variables

    type(namelistFile), pointer :: oneNamelistFile

    advmnt = oneNamelistFile%advmnt
    GhostZoneLength=oneNamelistFile%GhostZoneLength
  end subroutine StoreNamelistFileAtRadvc_mnt
end module ModMonotonicAdvection
