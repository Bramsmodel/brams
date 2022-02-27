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

    logical, parameter :: variable=.true.
    !<
    integer, parameter :: mzi=-2, myi=-2, mxi=-2

    integer :: mzpp3,mxpp3,mypp3
    integer :: mzppks,mxppis,myppjs

    logical, parameter :: dumpLocal=.false.
    character(len=8) :: str(10)

    call PostSendRecvMsgsVariableAdress(OneGrid%WideGhostZoneSend, OneGrid%WideGhostZoneRecv, &
         scp, ufx, vfx, wfx)

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

    !- flag to determine if a scalar is being advected
    scalar = .false.
    if(is==0 .and. js==0 .and. ks==0) scalar = .true.

    !- copy input arrays to extended arrays
    call copyMyPart(scp,scr,ufx_local,vfx_local,wfx_local, &
         ufx,vfx,wfx,mxp,myp,mzp,is,js,ks, &
         mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
         vname)

    if(.not. variable) &
         call expandBorder(mxp,myp,mzp,is,js,ks, &
         mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
         scr,ufx_local,vfx_local,wfx_local)

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
    end if

    call WaitSendRecvMsgs(OneGrid%WideGhostZoneSend, OneGrid%WideGhostZoneRecv, &
         1, mzp, scr, ufx_local, vfx_local, wfx_local)

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

    if(pd_or_mnt_constraint > 0 .and. scalar .and. vname .ne. 'thc' .and. vname .ne. 'pc') then  

       !-- positivity/monotonicity constraints
       call PosMonConstraints(mxp,myp,mzp,is,js,ks,ia,iz,ja,jz,&
            pd_or_mnt_constraint, &
            dt,ufx,vfx,wfx,vt3dh,vt3dj,vt3dk,scp, &
            scr,ufx_local,vfx_local,wfx_local, &
            qx,qy,qz,mzi,mzpp3,mxi,mxpp3,myi, &
            mypp3,mzppks,mxppis,myppjs,mynum,vname,sct)

    endif

    call CreateTendency(mxp,myp,mzp,is,js,ks,ia,iz,ja,jz, &
         mzppks,mxppis,myppjs, &
         dt,ufx,vfx,wfx, &
         vt3dh,vt3dj,vt3dk,scp, &
         qx,qy,qz,sct,vname,mynum)

    deallocate(qx ,qy,    qz,    scr,          ufx_local,    vfx_local ,   wfx_local)
