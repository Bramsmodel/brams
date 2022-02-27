    logical :: scalar
    real, pointer :: qz(:,:,:)
    real, pointer :: qx(:,:,:)
    real, pointer :: qy(:,:,:)
    real, pointer :: scr(:,:,:)
    real, pointer :: ufx_local(:,:,:)
    real, pointer :: vfx_local(:,:,:)
    real, pointer :: wfx_local(:,:,:)

    !external functions
    real, external :: flux_upwind,fq2, fq3, fq4, fq5, fq6, fq

    logical, parameter :: IsToDump=.false.
    !< Do a dump os communications?
    logical, parameter :: variable=.true.
    !<
    integer, parameter :: mzi=-2, myi=-2, mxi=-2

    integer :: mzpp3,mxpp3,mypp3
    integer :: mzppks,mxppis,myppjs

    logical, parameter :: dumpLocal=.true.
    character(len=8) :: str(10)

    if (dumpLocal) then
       call MsgDump(h//" to advect "//trim(adjustl(vname)))
    end if

    mzpp3=mzp+3; mxpp3=mxp+3; mypp3=myp+3
    mzppks=mzp+ks; mxppis=mxp+is; myppjs=myp+js

    allocate(qx(mzppks,mxppis,myppjs));qx=real_init
    allocate(qy(mzppks,mxppis,myppjs));qy=real_init
    allocate(qz(mzppks,mxppis,myppjs));qz=real_init
    allocate(scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3));scr=real_init
    allocate(ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3));ufx_local=real_init
    allocate(vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3));vfx_local=real_init
    allocate(wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3));wfx_local=real_init

!!$    call DumpGrid(OneGrid)
!!$    call parf_barrier(1)
!!$    call parf_exit_mpi()
!!$    if (mzppks > 0) stop
    
    call UpdateFieldAdress(OneGrid%LuFlaSendNorth, scr, "SCR")
    call UpdateFieldAdress(OneGrid%LuFlaSendNorth, ufx_local, "UFX")
    call UpdateFieldAdress(OneGrid%LuFlaSendNorth, vfx_local, "VFX")
    call UpdateFieldAdress(OneGrid%LuFlaSendNorth, wfx_local, "WFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvNorth, scr, "SCR")
    call UpdateFieldAdress(OneGrid%LuFlaRecvNorth, ufx_local, "UFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvNorth, vfx_local, "VFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvNorth, wfx_local, "WFX")

    call UpdateFieldAdress(OneGrid%LuFlaSendSouth, scr, "SCR")
    call UpdateFieldAdress(OneGrid%LuFlaSendSouth, ufx_local, "UFX")
    call UpdateFieldAdress(OneGrid%LuFlaSendSouth, vfx_local, "VFX")
    call UpdateFieldAdress(OneGrid%LuFlaSendSouth, wfx_local, "WFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvSouth, scr, "SCR")
    call UpdateFieldAdress(OneGrid%LuFlaRecvSouth, ufx_local, "UFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvSouth, vfx_local, "VFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvSouth, wfx_local, "WFX")

    call UpdateFieldAdress(OneGrid%LuFlaSendEast, scr, "SCR")
    call UpdateFieldAdress(OneGrid%LuFlaSendEast, ufx_local, "UFX")
    call UpdateFieldAdress(OneGrid%LuFlaSendEast, vfx_local, "VFX")
    call UpdateFieldAdress(OneGrid%LuFlaSendEast, wfx_local, "WFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvEast, scr, "SCR")
    call UpdateFieldAdress(OneGrid%LuFlaRecvEast, ufx_local, "UFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvEast, vfx_local, "VFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvEast, wfx_local, "WFX")

    call UpdateFieldAdress(OneGrid%LuFlaSendWest, scr, "SCR")
    call UpdateFieldAdress(OneGrid%LuFlaSendWest, ufx_local, "UFX")
    call UpdateFieldAdress(OneGrid%LuFlaSendWest, vfx_local, "VFX")
    call UpdateFieldAdress(OneGrid%LuFlaSendWest, wfx_local, "WFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvWest, scr, "SCR")
    call UpdateFieldAdress(OneGrid%LuFlaRecvWest, ufx_local, "UFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvWest, vfx_local, "VFX")
    call UpdateFieldAdress(OneGrid%LuFlaRecvWest, wfx_local, "WFX")

    !- flag to determine if a scalar is being advected
    scalar = .false.
    IF(is==0 .and. js==0 .and. ks==0) scalar = .true.

    if(IsToDump) &
         call dumpXYvar(scp,vname,'a',1,mxp,1,myp,1,mzp,600.0,660.0)

    !- copy input arrays to extended arrays
    call copyMyPart(scp,scr,ufx_local,vfx_local,wfx_local, &
         ufx,vfx,wfx,mxp,myp,mzp,is,js,ks, &
         mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
         vname)

    if(.not. variable) &
         call expandBorder(mxp,myp,mzp,is,js,ks, &
         mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
         scr,ufx_local,vfx_local,wfx_local)

    if(IsToDump) &
         call dumpXYvar(wfx_local,vname,'b',-2,mxp+3,-2,myp+3,-2,mzp+3,0.0,0.0)

    ! Set x & y boundary values in halo zones
    if (nmachs>1) then
       if (dumpLocal) then
          write(str(1),"(i8)") size(scr,1)
          write(str(2),"(i8)") size(scr,2)
          write(str(3),"(i8)") size(scr,3)
          call MsgDump(h//" exchange border of "//vname//" dimensioned ("//&
               trim(adjustl(str(1)))//","//&
               trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//")")
       end if
       call PostSendRecvMsgs(OneGrid%LuFlaSendNorth, OneGrid%LuFlaRecvNorth)
       call PostSendRecvMsgs(OneGrid%LuFlaSendSouth, OneGrid%LuFlaRecvSouth)
       call PostSendRecvMsgs(OneGrid%LuFlaSendEast, OneGrid%LuFlaRecvEast)
       call PostSendRecvMsgs(OneGrid%LuFlaSendWest, OneGrid%LuFlaRecvWest)
       call WaitSendRecvMsgs(OneGrid%LuFlaSendNorth, OneGrid%LuFlaRecvNorth)
       call WaitSendRecvMsgs(OneGrid%LuFlaSendSouth, OneGrid%LuFlaRecvSouth)
       call WaitSendRecvMsgs(OneGrid%LuFlaSendEast, OneGrid%LuFlaRecvEast)
       call WaitSendRecvMsgs(OneGrid%LuFlaSendWest, OneGrid%LuFlaRecvWest)
    end if

    if(IsToDump) &
         call dumpXYvar(wfx_local,vname,'c',-2,mxp+3,-2,myp+3,-2,mzp+3,0.0,0.0)

    select case (order_h)
    case (1)
       call compXYInterface_or1(mxp,myp,mzp,ks,is,js,&
            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
            mzppks,mxppis,myppjs, &
            scr,ufx_local,vfx_local,&
            border(mynum,:),qx,qy, &
            variable, vname)
    case (2)
       call compXYInterface_or2(mxp,myp,mzp,ks,is,js,&
            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
            mzppks,mxppis,myppjs, &
            scr,ufx_local,vfx_local,&
            border(mynum,:),qx,qy, &
            variable, vname)
    case (3)
       call compXYInterface_or3(mxp,myp,mzp,ks,is,js,&
            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
            mzppks,mxppis,myppjs, &
            scr,ufx_local,vfx_local,&
            border(mynum,:),qx,qy, &
            variable, vname)
    case (4)
       call compXYInterface_or4(mxp,myp,mzp,ks,is,js,&
            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
            mzppks,mxppis,myppjs, &
            scr,ufx_local,vfx_local,&
            border(mynum,:),qx,qy, &
            variable, vname)
    case(5,6)
       call compXYInterface_or56(mxp,myp,mzp,ks,is,js,&
            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
            mzppks,mxppis,myppjs, &
            scr,ufx_local,vfx_local,&
            border(mynum,:),qx,qy, &
            variable, vname, order_h)
    case default
       write (*,fmt='(A)') 'Advect Error : the order_h must be from 1 to 6'
       stop 'ERROR!'
    end select

    select case (order_v)
    case (1)
       call compZInterface_or1(mxp,myp,mzp,ks,is,js,&
            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
            mzppks,mxppis,myppjs, &
            scr,wfx_local,qz, &
            variable, vname)
    case (2)
       call compZInterface_or2(mxp,myp,mzp,ks,is,js,&
            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
            mzppks,mxppis,myppjs, &
            scr,wfx_local,qz, &
            variable, vname)
    case (3)
       call compZInterface_or3(mxp,myp,mzp,ks,is,js,&
            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
            mzppks,mxppis,myppjs, &
            scr,wfx_local,qz, &
            variable, vname)
    case (4)
       call compZInterface_or4(mxp,myp,mzp,ks,is,js,&
            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
            mzppks,mxppis,myppjs, &
            scr,wfx_local,qz, &
            variable, vname)
    case(5,6)
       call compZInterface_or56(mxp,myp,mzp,ks,is,js,&
            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
            mzppks,mxppis,myppjs, &
            scr,wfx_local,qz, &
            variable, vname,order_v)
    case default
       write (*,fmt='(A)') 'Advect Error : the order_v must be from 1 to 6'
       stop 'ERROR!'
    end select

    !print *, 'Is scalar? ',scalar,is,js,ks;call flush(6)
    !iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_notice,vname//' is scalar?: ',(/logical2Int(scalar),is,js,ks/),'I1')

    IF(pd_or_mnt_constraint > 0 .and. scalar .and. vname .ne. 'thc' .and. vname .ne. 'pc') THEN  

       !-- positivity/monotonicity constraints
       call PosMonConstraints(mxp,myp,mzp,is,js,ks,ia,iz,ja,jz,&
            pd_or_mnt_constraint, &
            dt,ufx,vfx,wfx,vt3dh,vt3dj,vt3dk,scp, &
            scr,ufx_local,vfx_local,wfx_local, &
            qx,qy,qz,mzi,mzpp3,mxi,mxpp3,myi, &
            mypp3,mzppks,mxppis,myppjs,mynum,vname,sct)

    ENDIF

    call CreateTendency(mxp,myp,mzp,is,js,ks,ia,iz,ja,jz, &
         mzppks,mxppis,myppjs, &
         dt,ufx,vfx,wfx, &
         vt3dh,vt3dj,vt3dk,scp, &
         qx,qy,qz,sct,vname,mynum)
    !-----check neg mass
    !IF(pd_or_mnt_constraint > 0 .and. scalar .and. vname .ne. 'thc' .and. vname .ne. 'pc') THEN  
    !IF( scalar .and. vname .ne. 'thc' .and. vname .ne. 'pc') THEN  
    !scp_new=scp+sct*dt
    !where(scp_new<0.) sct=-scp/dt
    !print*,"maxmin=",trim(vname),1.-abs(minval(scp_new)-maxval(scp_new))/(1.e-20+maxval(scp_new))
    !call flush(6)
    !ENDIF
    !-----check neg mass


    deallocate(qx ,qy,    qz,    scr,          ufx_local,    vfx_local ,   wfx_local)
