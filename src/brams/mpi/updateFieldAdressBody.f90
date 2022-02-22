    integer :: iMsg
    type(MessageData), pointer :: oneMessageData => null()
    character(len=8) :: c0
    logical, parameter :: dumpLocal=.true.

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if
    if (.not. associated(oneMessageSet)) then
       if (dumpLocal) then
          call MsgDump(h//" null MessageSet")
       end if
       return
    else if (.not. associated(field)) then
       call fatal_error(h//" field not associated")
    end if
    if (dumpLocal) then
       write(c0,"(i8)") oneMessageSet%nMsgs
       call MsgDump(h//" will run through "//&
            trim(adjustl(c0))//" message data of message set "//&
            trim(adjustl(oneMessageSet%name)))
    end if

    do iMsg = 1, oneMessageSet%nMsgs
       oneMessageData => oneMessageSet%msgData(iMsg)
       if (.not. associated(oneMessageData)) then
          write(c0,"(i8)") iMsg
          call fatal_error(h//" msgData("//trim(adjustl(c0))//&
               ") not associated")
       end if

       call UpdateFieldAdress(oneMessageData, field, name)
    end do
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
