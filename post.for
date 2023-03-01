!#############################################################
! List of subroutines
!#############################################################

!       subroutine tecplot_T(ki)
!               prints temperature T of the fluid cells in tecplot format

!       subroutine tecplot_S(ki)
!               prints Entrophy S of the fluid cells in tecplot format

!       subroutine tecplot_p(ki)
!               prints pressure p of the fluid cells in tecplot format using a time-stamp

!        subroutine csvplot_p(ki,flag)
!               prints pressure p of the fluid cells in csv format

!       subroutine tec_turb(ki)
!               prints turbulent fields of the fluid cells in tecplot format

!	subroutine tecplot_phi(ki)
!               prints varaible phi of the fluid cells when the Level Set Method is applied  in tecplot format 
!               It uses a time-stamp

!	subroutine tecbin(ki)
!               prints relevant field values of the fluid cells in order to re-start the simulation

!       subroutine tec_instant(ki)

!        subroutine tec_instant_u(ki)
!               prints velocity component u of the fluid cells in tecplot format using a time-stamp

!        subroutine tec_instant_v(ki)
!               prints velocity component v of the fluid cells in tecplot format using a time-stamp

!        subroutine tec_instant_w(ki)
!               prints velocity component w of the fluid cells in tecplot format using a time-stamp

!        subroutine csv_instant_u(ki)
!               prints velocity component u of the fluid cells in csv format

!        subroutine csv_instant_v(ki)
!               prints velocity component v of the fluid cells in csv format

!        subroutine csv_instant_w(ki)
!               prints velocity component w of the fluid cells in csv format

!        subroutine csv_instant_ustar(ki,flag)
!               prints velocity component ustar of the fluid cells in csv format
!               ustar is the non-solenoindal velocity component before computing the pressure equation

!        subroutine csv_instant_vstar(ki,flag)
!               prints velocity component vstar of the fluid cells in csv format
!               vstar is the non-solenoindal velocity component before computing the pressure equation

!        subroutine csv_instant_wstar(ki,flag)
!               prints velocity component wstar of the fluid cells in csv format
!               wstar is the non-solenoindal velocity component before computing the pressure equation

!        subroutine save2vtk(ki,flag)
!               prints one of these 4 field values of the fluid cells in vtk format using sintax of an unstructured mesh
!               flag = 1 ---> it prints velocity component u
!               flag = 2 ---> it prints velocity component v
!               flag = 3 ---> it prints velocity component w
!               flag = 4 ---> it prints pressure p

!        subroutine save2vtk_rectilinearGrif(ki,flag)
!               prints one of these 4 field values of the fluid cells in vtk format using sintax of a rectilinear mesh
!               flag = 1 ---> it prints velocity component u
!               flag = 2 ---> it prints velocity component v
!               flag = 3 ---> it prints velocity component w
!               flag = 4 ---> it prints pressure p

!       subroutine tec_inst_plane(ki)


!##########################################################################
        subroutine tecplot_T(ki)
!##########################################################################
        use multidata
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecout_T_'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep
        toti=(ie+1)-(is-1)+1
        totj=(je+1)-(js-1)+1
        totk=(ke+1)-(ks-1)+1

        write (88,*) 'title = Temperature field'
        write (88,*)'variables=x,y,z,T,Tm,Ttm'
        write (88,*)'zone ', ' i=',toti,', ',
     &  ' j=',totj,', k= ',totk,' f=point'
        do k=ks-1,ke+1
           do j=js-1,je+1
              do i=is-1,ie+1
                 write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%zc(k)
     & ,dom(ib)%T(i,j,k),dom(ib)%Tm(i,j,k),dom(ib)%Ttm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do
88      format (15e25.8)

        end
!##########################################################################
        subroutine tecplot_S(ki)
!##########################################################################
        use multidata
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecout_S_'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep
        toti=(ie+1)-(is-1)+1
        totj=(je+1)-(js-1)+1
        totk=(ke+1)-(ks-1)+1

        write (88,*) 'title =  Scalar transport field'
        write (88,*)'variables=x,y,z,S,Sm,Stm,SUM,SVM,SWM'
        write (88,*)'zone  i=',toti,',  j=',totj,', k= ',totk,' f=point'

        do k=ks-1,ke+1
           do j=js-1,je+1
              do i=is-1,ie+1
                 write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%zc(k)
     & ,dom(ib)%S(i,j,k),dom(ib)%Sm(i,j,k),dom(ib)%Stm(i,j,k)
     & ,dom(ib)%SUtm(i,j,k),dom(ib)%SVtm(i,j,k),dom(ib)%SWtm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do
88      format (15e25.8)
        end
!##########################################################################
        subroutine tecplot_p(ki)
!##########################################################################
	use vars
        use multidata
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
	integer :: ntoti,ntotj,ntotk, itoti,itotj,itotk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecout_p'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

!        is=dom(ib)%isp; ie=dom(ib)%iep
!        js=dom(ib)%jsp; je=dom(ib)%jep
!        ks=dom(ib)%ksp; ke=dom(ib)%kep
!        toti=(ie+1)-(is-1)+1
!        totj=(je+1)-(js-1)+1
!        totk=(ke+1)-(ks-1)+1

         itoti=1 	; ntoti=dom(ib)%ttc_i
         itotj=1 	; ntotj=dom(ib)%ttc_j
         itotk=1 	; ntotk=dom(ib)%ttc_k

        write (88,*) 'title = Pressure field'
        write (88,*)'variables=x,y,z,P,PM,ppm,vis,visM,ksgs,eps,epsm'

!        write (88,*)'zone ', ' i=',toti,', ',
!     &  ' j=',totj,', k= ',totk,' f=point'

	! This way we add the time stamp
         write (88,*)'zone ','STRANDID=', 2, 'SOLUTIONTIME=', ctime,
     &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
     &    'zonetype=', 'ordered',', DATAPACKING=point'


!        do k=ks-1,ke+1
!           do j=js-1,je+1
!              do i=is-1,ie+1

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti 
                 write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),
     & dom(ib)%zc(k),dom(ib)%p(i,j,k),dom(ib)%pm(i,j,k),
     & dom(ib)%ppm(i,j,k),dom(ib)%vis(i,j,k),dom(ib)%vism(i,j,k),
     & dom(ib)%ksgs(i,j,k),dom(ib)%eps(i,j,k),dom(ib)%epsm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do
88      format (15e25.8)

        end
!##########################################################################
        subroutine csvplot_p(ki)
!##########################################################################
	use vars
        use multidata
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
	integer :: ntoti,ntotj,ntotk, itoti,itotj,itotk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp
        !ib = 8

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='csvout_p'//trim(adjustl(chb))//'.csv'
        !open (unit=88, file=gf)
87      format(1x, *(g0, ",")) 
        OPEN(unit = 88, access = "sequential", action = "write",
     &   status = "replace", file = gf, form = "formatted") 

!        is=dom(ib)%isp; ie=dom(ib)%iep
!        js=dom(ib)%jsp; je=dom(ib)%jep
!        ks=dom(ib)%ksp; ke=dom(ib)%kep
!        toti=(ie+1)-(is-1)+1
!        totj=(je+1)-(js-1)+1
!        totk=(ke+1)-(ks-1)+1

         itoti=1 	; ntoti=dom(ib)%ttc_i
         itotj=1 	; ntotj=dom(ib)%ttc_j
         itotk=1 	; ntotk=dom(ib)%ttc_k

         write(88,*) 'x coord, y coord, z coord, p, pm,'
   !     write (88,*) 'title = Pressure field'
   !     write (88,*)'variables=x,y,z,P,PM,ppm,vis,visM,ksgs,eps,epsm'
!  !      write (88,*)'zone ', ' i=',toti,', ',
!  !   &  ' j=',totj,', k= ',totk,' f=point'
	! This way we add the time stamp
    !     write (88,*)'zone ','STRANDID=', 2, 'SOLUTIONTIME=', ctime,
    ! &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
    ! &    'zonetype=', 'ordered',', DATAPACKING=point'


!        do k=ks-1,ke+1
!           do j=js-1,je+1
!              do i=is-1,ie+1

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti 
                 write (88,87) dom(ib)%xc(i),dom(ib)%yc(j),
     & dom(ib)%zc(k),dom(ib)%p(i,j,k),dom(ib)%pm(i,j,k)
   !              write (88,87) dom(ib)%xc(i),dom(ib)%yc(j),
   !  & dom(ib)%zc(k),dom(ib)%p(i,j,k),dom(ib)%pm(i,j,k),
   !  & dom(ib)%ppm(i,j,k),dom(ib)%vis(i,j,k),dom(ib)%vism(i,j,k)!,
   !  & dom(ib)%ksgs(i,j,k),dom(ib)%eps(i,j,k),dom(ib)%epsm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do

!88      format (15e25.8)

        end
!##########################################################################
        subroutine tec_turb(ki)
!##########################################################################
        use multidata
        use vars
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        double precision u_cn,v_cn,w_cn,um_cn,vm_cn,wm_cn
        double precision uum_cn,vvm_cn,wwm_cn,uvml,uwml,vwml
        integer :: itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecturb'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=pl+1; ie=dom(ib)%ttc_i-pl
        js=pl+1; je=dom(ib)%ttc_j-pl
        ks=pl+1; ke=dom(ib)%ttc_k-pl
        toti=ie-(is-1)+1
        totj=je-(js-1)+1
        totk=ke-(ks-1)+1

!Plot every X points in every direction: change itotX respectively
	  itoti=1 	; ntoti=toti/itoti+1	
	  itotj=1 	; ntotj=totj/itotj+1
	  itotk=1 	; ntotk=totk/itotk+1


	IF(L_LSM.EQ..FALSE.) THEN
        write (88,*) 'title = turb'
        write (88,*)'variables=x,y,z,U,V,W,UM,VM,WM'
     & ,'uuM,vvM,wwM,uvM,uwM,vwM'
        write (88,*)'zone ', ' i=',ntoti,', ',
     &  ' j=',ntotj,', k= ',ntotk,' f=point'
	ELSE !WHEN LSM IS ALSO ACTIVATED
        write (88,*) 'title =  turb'
        write (88,*)'variables=x,y,z,U,V,W,UM,VM,WM'
     & ,'uuM,vvM,wwM,uvM,uwM,vwM,phi,p'
        write (88,*)'zone ', ' i=',ntoti,', ',
     &  ' j=',ntotj,', k= ',ntotk,' f=point'
	ENDIF

        do k=ks-1,ke+1,itotk	
           do j=js-1,je+1,itotj
              do i=is-1,ie+1,itoti    

                 u_cn  =0.25*(dom(ib)%u(i,j,k)+
     &dom(ib)%u(i,j+1,k)+dom(ib)%u(i,j,k+1)+
     &dom(ib)%u(i,j+1,k+1))
                 um_cn  =0.25*(dom(ib)%um(i,j,k)+
     &dom(ib)%um(i,j+1,k)+dom(ib)%um(i,j,k+1)+
     &dom(ib)%um(i,j+1,k+1))
                 uum_cn  =0.25*(dom(ib)%uum(i,j,k)+
     &dom(ib)%uum(i,j+1,k)+dom(ib)%uum(i,j,k+1)+
     &dom(ib)%uum(i,j+1,k+1))

                 v_cn  =0.25*(dom(ib)%v(i,j,k)+
     &dom(ib)%v(i+1,j,k)+dom(ib)%v(i,j,k+1)+
     &dom(ib)%v(i+1,j,k+1))
                 vm_cn  =0.25*(dom(ib)%vm(i,j,k)+
     &dom(ib)%vm(i+1,j,k)+dom(ib)%vm(i,j,k+1)+
     &dom(ib)%vm(i+1,j,k+1))
                 vvm_cn  =0.25*(dom(ib)%vvm(i,j,k)+
     &dom(ib)%vvm(i+1,j,k)+dom(ib)%vvm(i,j,k+1)+
     &dom(ib)%vvm(i+1,j,k+1))

                 w_cn  =0.25*(dom(ib)%w(i,j,k)+
     &dom(ib)%w(i+1,j,k)+dom(ib)%w(i,j+1,k)+
     &dom(ib)%w(i+1,j+1,k)) 
                 wm_cn  =0.25*(dom(ib)%wm(i,j,k)+
     &dom(ib)%wm(i+1,j,k)+dom(ib)%wm(i,j+1,k)+
     &dom(ib)%wm(i+1,j+1,k)) 
                 wwm_cn  =0.25*(dom(ib)%wwm(i,j,k)+
     &dom(ib)%wwm(i+1,j,k)+dom(ib)%wwm(i,j+1,k)+
     &dom(ib)%wwm(i+1,j+1,k)) 

                 uvml  =0.125*(dom(ib)%uvm(i,j,k)+
     &dom(ib)%uvm(i+1,j,k)    +dom(ib)%uvm(i,j+1,k)+
     &dom(ib)%uvm(i+1,j+1,k)  +dom(ib)%uvm(i,j,k+1)+
     &dom(ib)%uvm(i+1,j,k+1)  +dom(ib)%uvm(i,j+1,k+1)+
     &dom(ib)%uvm(i+1,j+1,k+1))
                 uwml  =0.125*(dom(ib)%uwm(i,j,k)+
     &dom(ib)%uwm(i+1,j,k)    +dom(ib)%uwm(i,j+1,k)+
     &dom(ib)%uwm(i+1,j+1,k)  +dom(ib)%uwm(i,j,k+1)+
     &dom(ib)%uwm(i+1,j,k+1)  +dom(ib)%uwm(i,j+1,k+1)+
     &dom(ib)%uwm(i+1,j+1,k+1))
                 vwml  =0.125*(dom(ib)%vwm(i,j,k)+
     &dom(ib)%vwm(i+1,j,k)    +dom(ib)%vwm(i,j+1,k)+
     &dom(ib)%vwm(i+1,j+1,k)  +dom(ib)%vwm(i,j,k+1)+
     &dom(ib)%vwm(i+1,j,k+1)  +dom(ib)%vwm(i,j+1,k+1)+
     &dom(ib)%vwm(i+1,j+1,k+1))

	IF(L_LSM.EQ..FALSE.) THEN
          write (88,88) dom(ib)%x(i),dom(ib)%y(j),
     &  dom(ib)%z(k),u_cn,v_cn,w_cn,um_cn,vm_cn,wm_cn
     &  ,uum_cn,vvm_cn,wwm_cn,uvml,uwml,vwml
	else
          write (88,88) dom(ib)%x(i),dom(ib)%y(j),
     &  dom(ib)%z(k),u_cn,v_cn,w_cn,um_cn,vm_cn,wm_cn
     &  ,uum_cn,vvm_cn,wwm_cn,uvml,uwml,vwml
     &  ,dom(ib)%phi(i,j,k),dom(ib)%p(i,j,k)
	ENDIF
              end do
           end do
        end do

        close (88)

        end do
88      format (20e25.8)
        end
!##########################################################################
        subroutine tecplot_phi(ki)
!##########################################################################
        use vars
        use multidata
        implicit none

        integer :: sn,sn1,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb,chb1
        character*27 :: gf

        do ib=1,nbp

	  if (L_anim_phi) then
          write(chb,'(i4)') dom_id(ib)
          write(chb1,'(i6)') ki
          sn=len(trim(adjustl(chb)))
          sn1=len(trim(adjustl(chb1)))
          chb=repeat('0',(4-sn))//trim(adjustl(chb))
          chb1=repeat('0',(6-sn1))//trim(adjustl(chb1))
          gf='tecout_phi_'//trim(adjustl(chb))//'_'//
     & trim(adjustl(chb1))//'.plt'
	  else
          write(chb,'(i4)') dom_id(ib)
          sn=len(trim(adjustl(chb)))
          chb=repeat('0',(4-sn))//trim(adjustl(chb))
          gf='tecout_phi_'//trim(adjustl(chb))//'.plt'
	  endif      

	  open (unit=88, file=gf)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep
        toti=(ie+1)-(is-1)+1
        totj=(je+1)-(js-1)+1
        totk=(ke+1)-(ks-1)+1

        write (88,*)'title = Free-surface from LSM'
        write (88,*)'variables=x,y,z,phi,phim,dens,mu'
	  if (L_anim_phi) then
          write (88,*)'zone ','STRANDID=', 1, 'SOLUTIONTIME=', ctime,
     &    ' i=',toti,', ',' j=',totj,', k= ',totk,
     &    'zonetype=', 'ordered',', DATAPACKING=point'
	  else
          write(88,*)'zone  i=',toti,', j=',totj,', k= ',totk,' f=point'
	  endif

        do k=ks-1,ke+1
           do j=js-1,je+1
              do i=is-1,ie+1
                write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%zc(k)
     & ,dom(ib)%phi(i,j,k),dom(ib)%phim(i,j,k),dom(ib)%dens(i,j,k),
     & dom(ib)%mu(i,j,k)
              end do
           end do
        end do
        close (88)

        end do

88      format (8e20.8)

        end
!##########################################################################
        subroutine tecbin(ki)
!##########################################################################
        use multidata
        use vars
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: ki,ib,inind,jnind,knind
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecbin'//trim(adjustl(chb))//'.bin'
        open (unit=88, file=gf, form='unformatted')

        toti=dom(ib)%ttc_i
        totj=dom(ib)%ttc_j
        totk=dom(ib)%ttc_k

        write (88) toti,totj,totk
        write (88) pl
!====================================================================
        inind=0; jnind=0; knind=0
        if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5) inind=-1
        if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5) jnind=-1
        if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5) knind=-1
        write (88) inind,jnind,knind
!====================================================================

        do k=1,totk
           do j=1,totj
              do i=1,toti

        write (88) dom(ib)%x(i),dom(ib)%y(j),dom(ib)%z(k),
     & dom(ib)%p(i,j,k),dom(ib)%pm(i,j,k),dom(ib)%ppm(i,j,k),
     & dom(ib)%vis(i,j,k),dom(ib)%vism(i,j,k),
     & dom(ib)%u(i,j,k),dom(ib)%um(i,j,k),dom(ib)%uum(i,j,k),
     & dom(ib)%v(i,j,k),dom(ib)%vm(i,j,k),dom(ib)%vvm(i,j,k),
     & dom(ib)%w(i,j,k),dom(ib)%wm(i,j,k),dom(ib)%wwm(i,j,k),
     & dom(ib)%uvm(i,j,k),dom(ib)%uwm(i,j,k),dom(ib)%vwm(i,j,k),
     & dom(ib)%ksgs(i,j,k),dom(ib)%eps(i,j,k),dom(ib)%epsm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do
        end
!##########################################################################
        subroutine tec_instant(ki)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecinst'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
        open (unit=95, file=gf)        

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1



!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
	   itotk=1 	; ntotk=dom(ib)%ttc_k

        write(95,*)'variables=x,y,z,U,V,W'
!        write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'
          write (95,*)'zone ','STRANDID=', 1, 'SOLUTIONTIME=', ctime,
     &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
     &    'zonetype=', 'ordered',', DATAPACKING=point'

!        do k=ks-1,ke+1,itotk	
!           do j=js-1,je+1,itotj
!              do i=is-1,ie+1,itoti     

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti     
             write (95,88) dom(ib)%x(i),dom(ib)%y(j),dom(ib)%z(k),
     &  dom(ib)%u(i,j,k),dom(ib)%v(i,j,k),dom(ib)%w(i,j,k)
              end do
           end do
        end do

        close (95)
        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end


!##########################################################################
        subroutine tec_instant_u(ki)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecinst_u'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
        open (unit=95, file=gf)        

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1



!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
	   itotk=1 	; ntotk=dom(ib)%ttc_k

        write(95,*)'variables=x,y,z,U'
!        write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'
          write (95,*)'zone ','STRANDID=', 111, 'SOLUTIONTIME=', ctime,
     &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
     &    'zonetype=', 'ordered',', DATAPACKING=point'

!        do k=ks-1,ke+1,itotk	
!           do j=js-1,je+1,itotj
!              do i=is-1,ie+1,itoti     

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti     
             write (95,88) dom(ib)%x(i),dom(ib)%yc(j),dom(ib)%zc(k),
     &  dom(ib)%u(i,j,k)
              end do
           end do
        end do

        close (95)
        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end

!##########################################################################
        subroutine tec_instant_v(ki)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecinst_v'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
        open (unit=95, file=gf)        

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1



!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
	   itotk=1 	; ntotk=dom(ib)%ttc_k

        write(95,*)'variables=x,y,z,V'
!        write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'
          write (95,*)'zone ','STRANDID=', 222, 'SOLUTIONTIME=', ctime,
     &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
     &    'zonetype=', 'ordered',', DATAPACKING=point'

!        do k=ks-1,ke+1,itotk	
!           do j=js-1,je+1,itotj
!              do i=is-1,ie+1,itoti     

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti     
             write (95,88) dom(ib)%xc(i),dom(ib)%y(j),dom(ib)%zc(k),
     &  dom(ib)%v(i,j,k)
              end do
           end do
        end do

        close (95)
        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end

!##########################################################################
        subroutine tec_instant_w(ki)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecinst_w'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
        open (unit=95, file=gf)        

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1



!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
	   itotk=1 	; ntotk=dom(ib)%ttc_k

        write(95,*)'variables=x,y,z,W'
!        write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'
          write (95,*)'zone ','STRANDID=', 333, 'SOLUTIONTIME=', ctime,
     &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
     &    'zonetype=', 'ordered',', DATAPACKING=point'

!        do k=ks-1,ke+1,itotk	
!           do j=js-1,je+1,itotj
!              do i=is-1,ie+1,itoti     

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti     
             write (95,88) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%z(k),
     &  dom(ib)%w(i,j,k)
              end do
           end do
        end do

        close (95)
        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end
!##########################################################################
        subroutine csv_instant_u(ki)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp
        !ib = 8

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='csvinst_u'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.csv'
        !open (unit=95, file=gf)
87      format(1x, *(g0, ",")) 
        OPEN(unit = 95, access = "sequential", action = "write",
     &   status = "replace", file = gf, form = "formatted") 

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1



!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
	   itotk=1 	; ntotk=dom(ib)%ttc_k

         write(95,*) 'x coord, y coord, z coord, u'
  !      write(95,*)'variables=x,y,z,U'
! !       write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'
  !        write (95,*)'zone ','STRANDID=', 111, 'SOLUTIONTIME=', ctime,
  !   &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
  !   &    'zonetype=', 'ordered',', DATAPACKING=point'

!        do k=ks-1,ke+1,itotk	
!           do j=js-1,je+1,itotj
!              do i=is-1,ie+1,itoti     

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti     
   !          write (95,88) dom(ib)%x(i),dom(ib)%yc(j),dom(ib)%zc(k),
              write (95,87) dom(ib)%x(i) ,dom(ib)%yc(j), dom(ib)%zc(k),
     &  dom(ib)%u(i,j,k)
              end do
           end do
        end do

        close (95)

        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)
        !87      format (e11.4,a,e11.4,a,e11.4,a,e11.4)
        !87      format '(1x, F, 3(",", F))'
        
        end

!##########################################################################
        subroutine csv_instant_v(ki)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp
        !ib = 8

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='csvinst_v'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.csv'
        !open (unit=95, file=gf)
87      format(1x, *(g0, ",")) 
        OPEN(unit = 95, access = "sequential", action = "write",
     &   status = "replace", file = gf, form = "formatted") 

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1



!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
	   itotk=1 	; ntotk=dom(ib)%ttc_k

         write(95,*) 'x coord, y coord, z coord, v'
   !     write(95,*)'variables=x,y,z,V'
!  !      write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'
   !       write (95,*)'zone ','STRANDID=', 222, 'SOLUTIONTIME=', ctime,
   !  &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
   !  &    'zonetype=', 'ordered',', DATAPACKING=point'

!        do k=ks-1,ke+1,itotk	
!           do j=js-1,je+1,itotj
!              do i=is-1,ie+1,itoti     

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti     
   !          write (95,88) dom(ib)%xc(i),dom(ib)%y(j),dom(ib)%zc(k),
              write (95,87) dom(ib)%xc(i),dom(ib)%y(j),dom(ib)%zc(k),
     &  dom(ib)%v(i,j,k)
              end do
           end do
        end do

        close (95)

        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end

!##########################################################################
        subroutine csv_instant_w(ki)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp
        !ib = 8

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='csvinst_w'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.csv'
        !open (unit=95, file=gf)        
87      format(1x, *(g0, ",")) 
        OPEN(unit = 95, access = "sequential", action = "write",
     &   status = "replace", file = gf, form = "formatted") 

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1



!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
           itotk=1 	; ntotk=dom(ib)%ttc_k

         write(95,*) 'x coord, y coord, z coord, w'
  !      write(95,*)'variables=x,y,z,W'
! !       write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'
  !        write (95,*)'zone ','STRANDID=', 333, 'SOLUTIONTIME=', ctime,
  !   &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
  !   &    'zonetype=', 'ordered',', DATAPACKING=point'

!        do k=ks-1,ke+1,itotk	
!           do j=js-1,je+1,itotj
!              do i=is-1,ie+1,itoti     

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti     
  !           write (95,88) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%z(k),
  !   &  dom(ib)%w(i,j,k)
             write (95,87) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%z(k),
     &  dom(ib)%w(i,j,k)
              end do
           end do
        end do

        close (95)

        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end
!##########################################################################
        subroutine csv_instant_ustar(ki,flag)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer, intent(in):: flag
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp
        !ib = 8

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        if(flag.eq.1) then      ! 2 to indicate that is right before applying the ghost method
        gf='csvustar1'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.csv'
        elseif(flag.eq.2) then  ! 2 to indicate that is right after applying the ghost method
        gf='csvustar2'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.csv'
        endif
        !open (unit=95, file=gf)
87      format(1x, *(g0, ",")) 
        OPEN(unit = 95, access = "sequential", action = "write",
     &   status = "replace", file = gf, form = "formatted") 

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1



!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
	   itotk=1 	; ntotk=dom(ib)%ttc_k

         write(95,*) 'x coord, y coord, z coord, ustar'
  !      write(95,*)'variables=x,y,z,U'
! !       write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'
  !        write (95,*)'zone ','STRANDID=', 111, 'SOLUTIONTIME=', ctime,
  !   &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
  !   &    'zonetype=', 'ordered',', DATAPACKING=point'

!        do k=ks-1,ke+1,itotk	
!           do j=js-1,je+1,itotj
!              do i=is-1,ie+1,itoti     

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti     
              write (95,87) dom(ib)%x(i) ,dom(ib)%yc(j), dom(ib)%zc(k),
     &  dom(ib)%ustar(i,j,k)
              end do
           end do
        end do

        close (95)

        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)
        !87      format (e11.4,a,e11.4,a,e11.4,a,e11.4)
        !87      format '(1x, F, 3(",", F))'
        
        end

!##########################################################################
        subroutine csv_instant_vstar(ki,flag)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer, intent(in):: flag
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp
        !ib = 8

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        if(flag.eq.1) then      ! 2 to indicate that is right before applying the ghost method
        gf='csvvstar1'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.csv'
        elseif(flag.eq.2) then  ! 2 to indicate that is right after applying the ghost method
        gf='csvvstar2'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.csv'
        endif
        !open (unit=95, file=gf)
87      format(1x, *(g0, ",")) 
        OPEN(unit = 95, access = "sequential", action = "write",
     &   status = "replace", file = gf, form = "formatted") 

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1



!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
	   itotk=1 	; ntotk=dom(ib)%ttc_k

         write(95,*) 'x coord, y coord, z coord, vstar'
   !     write(95,*)'variables=x,y,z,V'
!  !      write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'
   !       write (95,*)'zone ','STRANDID=', 222, 'SOLUTIONTIME=', ctime,
   !  &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
   !  &    'zonetype=', 'ordered',', DATAPACKING=point'

!        do k=ks-1,ke+1,itotk	
!           do j=js-1,je+1,itotj
!              do i=is-1,ie+1,itoti     

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti     
              write (95,87) dom(ib)%xc(i),dom(ib)%y(j),dom(ib)%zc(k),
     &  dom(ib)%vstar(i,j,k)
              end do
           end do
        end do

        close (95)

        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end

!##########################################################################
        subroutine csv_instant_wstar(ki,flag)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer, intent(in):: flag
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp
        !ib = 8

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        if(flag.eq.1) then      ! 2 to indicate that is right before applying the ghost method
        gf='csvwstar1'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.csv'
        elseif(flag.eq.2) then  ! 2 to indicate that is right after applying the ghost method
        gf='csvwstar2'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.csv'
        endif
        !open (unit=95, file=gf)        
87      format(1x, *(g0, ",")) 
        OPEN(unit = 95, access = "sequential", action = "write",
     &   status = "replace", file = gf, form = "formatted") 

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1



!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
           itotk=1 	; ntotk=dom(ib)%ttc_k

         write(95,*) 'x coord, y coord, z coord, wstar'
  !      write(95,*)'variables=x,y,z,W'
! !       write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'
  !        write (95,*)'zone ','STRANDID=', 333, 'SOLUTIONTIME=', ctime,
  !   &    ' i=',ntoti,', ',' j=',ntotj,', k= ',ntotk,
  !   &    'zonetype=', 'ordered',', DATAPACKING=point'

!        do k=ks-1,ke+1,itotk	
!           do j=js-1,je+1,itotj
!              do i=is-1,ie+1,itoti     

         do k = 1,ntotk, itotk	
            do j = 1, ntotj, itotj
               do i = 1, ntoti,itoti     
             write (95,87) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%z(k),
     &  dom(ib)%wstar(i,j,k)
              end do
           end do
        end do

        close (95)

        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end
!##########################################################################
        subroutine save2vtk(ki,flag)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none

        integer, intent(in):: flag ! Indicates fields to write: u, v, w, p
        character*7 :: fileName ! file name
        character*70 :: heather ! Heather in the vtk file

        !integer, intent(in):: tf,pf,nne
        integer:: tf,pf,nne
        ! tf: number of rows in matrix t
        ! pf: number of rows in matrix p
        ! nne: number of points to read in each row of the element table
        !p: table of node coordinates
        !t: table element (rows) and their node inidexes (columns)
        !e: table of nodes per side

        ! change this to global variables
        DOUBLE PRECISION,allocatable,dimension(:,:):: p, t

        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        ! ki should have intent(in)
        character*8 :: chb,chb2
        character*35 :: gf ! file name
        character*35 :: message,message2,message3
        integer:: nds,n1,n2,n3,n4
        character*8 :: pf_str, tf_str, nds_str


!	if (nbp.gt.1) RETURN

        do ib=1,nbp
        !ib = 8

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        if (flag.eq.1) then
                fileName = 'vtk_u__'
                heather ='2D u velocity Distribution'
        elseif(flag.eq.2) then
                fileName = 'vtk_v__'
                heather ='2D v velocity Distribution'
        elseif(flag.eq.3) then
               fileName = 'vtk_w__'
               heather ='2D w velocity Distribution'
        elseif(flag.eq.4) then
                fileName = 'vtk_p'
                heather ='2D p pressure Distribution'
        elseif(flag.eq.5) then
                fileName = 'vtk_mu'
                heather ='2D mu Distribution'
        elseif(flag.eq.6) then
                fileName = 'vtk_nut'
                heather ='2D vis Distribution'
        elseif(flag.eq.7) then
                fileName = 'vtk_phi'
                heather ='2D phi Distribution'
        endif 


        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))

        gf=fileName//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.vtk'

        !open (unit=95, file=gf)        
87      format(1x, *(g0, " ")) 
        OPEN(unit = 95, access = "sequential", action = "write",
     &   status = "replace", file = gf, form = "formatted") 

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1

!Plot every X points in every direction: change itotX respectively
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j
           itotk=1 	; ntotk=dom(ib)%ttc_k

         !write(95,*) '# vtk DataFile Version 3.0'
         write(95,'(a)') '# vtk DataFile Version 3.0'
        ! check
         !write(95,*) trim(adjustl(heather)) !'2D Temperature Distribution\n'
         write(95,'(a)') adjustl(trim(adjustl(heather))) !'2D Temperature Distribution\n'

         write(95,'(a)') 'ASCII'

         write(95,'(a)') 'DATASET UNSTRUCTURED_GRID'

         ! lets get table p (nodes coordinates)
         !pf = 0.0d00
         !do k = 1,ntotk, itotk	
         !!do j = 1, ntotj, itotj
         !j = 8 ! lets just take the cross-section along the centerline
         !do i = 1, ntoti,itoti
         !       pf = pf + 1
         !end do
         !!end do
         !end do

         ! Alternative
         !pf = ntoti*ntotk
         pf = ntoti*ntotj*ntotk ! general form for subdomain ib
         allocate(p(pf,3))

         write(pf_str,'(i8)') pf
         message = 'POINTS '//trim(adjustl(pf_str))//' float'
         ! check
         !message = trim(adjustl(message))
         !write(95,*) message
         write(95,'(a)') trim(adjustl(message))


         !do k = 1,ntotk, itotk	
         !!do j = 1, ntotj, itotj
         !j = 8 ! lets just take the cross-section along the centerline
         !do i = 1, ntoti,itoti
         !      p( i + ntoti*(k-1),1 ) = dom(ib)%xc(i)
         !      p( i + ntoti*(k-1),2 ) = dom(ib)%yc(j)
         !      p( i + ntoti*(k-1),3 ) = dom(ib)%z(k)
         !      !p( i + ntoti*(j-1) + ntoti*ntotj(k-1),1 ) = dom(ib)%xc(i)
         !      !p( i + ntoti*(j-1) + ntoti*ntotj(k-1),2 ) = dom(ib)%yc(j)
         !      !p( i + ntoti*(j-1) + ntoti*ntotj(k-1),3 ) = dom(ib)%z(k)
         !end do
         !!end do
         !end do

         !do k = 1,pf	
         !        write (95,87) p(k,1),p(k,2),p(k,3)
         !&  dom(ib)%wstar(i,j,k)
         !   end do

         ! Alternative, this way we do not need matrix p
         do k = 1,ntotk, itotk	
         do j = 1, ntotj, itotj
         !j = 8 ! lets just take the cross-section along the centerline
         do i = 1, ntoti,itoti
         if (flag.eq.1) then
                ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                write (95,87) dom(ib)%x(i),dom(ib)%yc(j),dom(ib)%zc(k)
        elseif(flag.eq.2) then
                write (95,87) dom(ib)%xc(i),dom(ib)%y(j),dom(ib)%zc(k)
        elseif(flag.eq.3) then
                write (95,87) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%z(k)
        elseif(flag.eq.4) then
                write (95,87) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%zc(k)
        endif 

         end do
         end do
         end do



         ! lets get table t (elements and node inidexes)
         !tf = (ntoti-1)*(ntotk-1)
         tf = (ntoti-1)*(ntotj-1)*(ntotk-1) ! general form for subdomain ib
         allocate(t(tf,4))


        ! the orientation of the nodes inside the element is 1, 2, 4, 3
         !      2-------3
         !      |       |
         !      |       |
         !      0-------1
         !
         !do k = 1,ntotk-1, itotk	
         !!do j = 1, ntotj-1, itotj
         !j = 8 ! lets just take the cross-section along the centerline
         !do i = 1, ntoti-1,itoti 
         !      t( i+(ntoti-1)*(k-1),1 ) = i
         !      t( i+(ntoti-1)*(k-1),2 ) = i+1
         !      t( i+(ntoti-1)*(k-1),3 ) = i+1+ntoti*(k)
         !      t( i+(ntoti-1)*(k-1),4 ) = i+ntoti*(k)
         !      !t( i+(ntoti-1)*(j-1)+(ntoti-1)*(ntotj-1)*(k-1),1 ) = i
         !      !t( i+(ntoti-1)*(j-1)+(ntoti-1)*(ntotj-1)*(k-1),2 ) = i+1
         !      !t( i+(ntoti-1)*(j-1)+(ntoti-1)*(ntotj-1)*(k-1),3 ) = i+1+ntoti*(k)
         !      !t( i+(ntoti-1)*(j-1)+(ntoti-1)*(ntotj-1)*(k-1),4 ) = i+ntoti*(k)
         !end do
         !!end do
         !end do

         nne = 4 ! check this on vtk documentation
         write(tf_str,'(i8)') tf
         nds = (tf*(nne+1)) ! size of the matrix to read by the vtk reader, 
         !there is an extra column to indicate the number of points to read
         write(nds_str,'(i8)') nds
         message = 'CELLS '//trim(adjustl(tf_str))//
     &		' '//trim(adjustl(nds_str))

         write(95,'(a)') trim(adjustl(message))
         !do k = 1, tf
         !       write (95,87) nne,t(k,1),t(k,2),t(k,3),t(k,4)
         !enddo

        !Alternative
         do k = 0,ntotk-2, itotk	
         do j = 1, ntotj-1, itotj
         !j = 8 ! lets just take the cross-section along the centerline
         do i = 0, ntoti-2,itoti 
                ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                n1 = i+ntoti*(k)
                n2 = i+1+ntoti*(k)
                n3 = i+1+ntoti*(k+1) 
                n4 = i+ntoti*(k+1)
                write(95,87) nne, n1,n2,n3,n4
         end do
         end do
         end do

         message = 'CELL_TYPES '//trim(adjustl(tf_str))
         message = trim(adjustl(message))
         write(95,'(a)') message
         do k = 1, tf
                write(95,87) 9 ! for a quad element use index 9, , for an hexaedrum use index 11
         enddo

        !Alternative
         !do k = 1,ntotk-1, itotk	
         !!do j = 1, ntotj-1, itotj
         !j = 8 ! lets just take the cross-section along the centerline
         !do i = 1, ntoti-1,itoti 
         !       ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
         !       write(95,87) 8 ! for a quad element use index 8, for an hexaedrum use index 11
         !end do
         !!end do
         !end do


         message = 'CELL_DATA '//trim(adjustl(tf_str))
         !message = trim(adjustl(message))
         !write(95,*) message
         write(95,'(a)') trim(adjustl(message))

         message = 'NORMALS cell_normals float'
         !message = trim(adjustl(message))
         !write(95,*) message
         write(95,'(a)') trim(adjustl(message))         

        do k = 1, tf    ! lets write the normals of each surface element
                ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                write (95,87) 0,1,0
        enddo

        message = 'POINT_DATA '//trim(adjustl(pf_str))
         !message = trim(adjustl(message))
         !write(95,*) message
        write(95,'(a)') trim(adjustl(message))

        if (flag.eq.1) then
        message = 'SCALARS u_velocity float '//trim(adjustl('1')) ! '1' inidicates the number of components
        elseif(flag.eq.2) then
        message = 'SCALARS v_velocity float '//trim(adjustl('1'))
        elseif(flag.eq.3) then
        message = 'SCALARS w_velocity float '//trim(adjustl('1'))
        elseif(flag.eq.4) then
        message = 'SCALARS p_pressure float '//trim(adjustl('1'))
        elseif(flag.eq.5) then
        message = 'SCALARS mu_viscosity float '//trim(adjustl('1'))
        elseif(flag.eq.6) then
        message = 'SCALARS vis_viscosity float'//trim(adjustl('1'))
        elseif(flag.eq.7) then
        message = 'SCALARS phi_LSM float'//trim(adjustl('1'))
        endif 
         !message = trim(adjustl(message))
         !write(95,*) message
        write(95,'(a)') trim(adjustl(message))

        message = 'LOOKUP_TABLE default'
         !message = trim(adjustl(message))
         !write(95,*) message
        write(95,'(a)') trim(adjustl(message))

        ! lets write the field data
        !do k = 1, pf    
        !        write (95,87) dom(ib)%wstar(i,j,k) ! we can not just use wstar here since it does not match pf in the number of rows
        !enddo


         ! Alternative, this way we do not need matrix p
         do k = 1,ntotk, itotk	
         do j = 1, ntotj, itotj
         !j = 8 ! lets just take the cross-section along the centerline
         do i = 1, ntoti,itoti
         if (flag.eq.1) then
                ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                write (95,87) dom(ib)%u(i,j,k)
        elseif(flag.eq.2) then
                write (95,87) dom(ib)%v(i,j,k)
        elseif(flag.eq.3) then
                write (95,87) dom(ib)%w(i,j,k)
        elseif(flag.eq.4) then
                write (95,87) dom(ib)%p(i,j,k)
        elseif(flag.eq.5) then
                write (95,87) dom(ib)%mu(i,j,k)
        elseif(flag.eq.6) then
                write (95,87) dom(ib)%vis(i,j,k)
        elseif(flag.eq.7) then
                write (95,87) dom(ib)%phi(i,j,k)
        endif 

         end do
         end do
         end do


         deallocate(p,t)

        close (95)

        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end

!##########################################################################
        subroutine save2vtk_rectilinearGrid(ki,flag)!,flag)
!##########################################################################
                use multidata
                use vars
                use mpi
                implicit none
        
                integer, intent(in):: flag ! Indicates fields to write: u, v, w, p
                character*7 :: fileName ! file name
                character*70 :: heather ! Heather in the vtk file
        
                !integer, intent(in):: tf,pf,nne
                integer:: tf,pf,nne
                ! tf: number of rows in matrix t
                ! pf: number of rows in matrix p
                ! nne: number of points to read in each row of the element table
                !p: table of node coordinates
                !t: table element (rows) and their node inidexes (columns)
                !e: table of nodes per side
        
                ! change this to global variables
                DOUBLE PRECISION,allocatable,dimension(:,:):: p, t
        
                integer :: sn,i,j,k,toti,totj,totk
                integer :: is,ie,js,je,ks,ke
                integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
                ! ki should have intent(in)
                character*8 :: chb,chb2,nx_str,ny_str,nz_str
                character*35 :: gf ! file name
                character*35 :: message,message2,message3
                integer:: nds,n1,n2,n3,n4
                character*8 :: pf_str, tf_str, nds_str

                do ib = 1, nbp ! lets sweep only the subdomains of each procesor
        if (flag.eq.1) then
                fileName = 'vtk_u__'
                heather ='2D U Velocity Distribution'
        elseif(flag.eq.2) then
                fileName = 'vtk_v__'
                heather ='2D V Velocity Distribution'
        elseif(flag.eq.3) then
               fileName = 'vtk_w__'
               heather ='2D W Velocity Distribution'
        elseif(flag.eq.4) then
                fileName = 'vtk_p__'
                heather ='2D P Pressure Distribution'
        elseif(flag.eq.5) then
                fileName = 'vtk_nu'
                heather ='2D mu Distribution'
        elseif(flag.eq.6) then
                fileName = 'vtk_nut'
                heather ='2D vis Distribution'
        elseif(flag.eq.7) then
                fileName = 'vtk_phi'
                heather ='2D phi LSM Distribution'
        elseif(flag.eq.8) then
                fileName = 'vtk_rho'
                heather ='2D rho Density Distribution'
        elseif (flag.eq.-1) then
                fileName = 'vtk_uoo__'
                heather ='2D Uoo Velocity Distribution'
        elseif(flag.eq.-2) then
                fileName = 'vtk_voo__'
                heather ='2D Voo Velocity Distribution'
        elseif(flag.eq.-3) then
               fileName = 'vtk_woo__'
               heather ='2D Woo Velocity Distribution'
        elseif(flag.eq.-4) then
               fileName = 'vtk_p_o__'
               heather ='2D P_o Pressure Distribution'
        endif 

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))

        gf=trim(adjustl(chb))//fileName//
     &		'.'//trim(adjustl(chb2))//'.vtk'

        !open (unit=95, file=gf)
        ! CSV format
87      format(1x, *(g0, " ")) 
        OPEN(unit = 95, access = "sequential", action = "write",
     &   status = "replace", file = gf, form = "formatted") 

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1

!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj=dom(ib)%ttc_j   !ntotj = 1 
           itotk=1 	; ntotk=dom(ib)%ttc_k


         !write(95,*) '# vtk DataFile Version 3.0'
           write(95,'(a)') '# vtk DataFile Version 3.0'
           ! check
            !write(95,*) trim(adjustl(heather)) !'2D Temperature Distribution\n'
            write(95,'(a)') adjustl(trim(adjustl(heather))) !'2D Temperature Distribution\n'
   
            write(95,'(a)') 'ASCII'
   
            write(95,'(a)') 'DATASET RECTILINEAR_GRID'
   
   
            write(nx_str,'(i8)') ntoti
            write(ny_str,'(i8)') ntotj
            write(nz_str,'(i8)') ntotk
            message = 'DIMENSIONS '//adjustl(trim(adjustl(nx_str)))//' '
     &  //adjustl(trim(adjustl(ny_str)))//' '
     &  //adjustl(trim(adjustl(nz_str)))
            write(95,'(a)') trim(adjustl(message))
   
   
            message = 'X_COORDINATES '//trim(adjustl(nx_str))//' float'
            write(95,'(a)') trim(adjustl(message))
            do i = 1, ntoti,itoti
            if ((flag.eq.1).or.(flag.eq.-1)) then
                   ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                   write (95,87) dom(ib)%x(i)
            else
                   write (95,87) dom(ib)%xc(i)
            endif 
            enddo

            message = 'Y_COORDINATES '//trim(adjustl(ny_str))//' float'
            write(95,'(a)') trim(adjustl(message))
            do j = 1, ntotj
            !j = 8
            if ((flag.eq.2).or.(flag.eq.-2)) then
                   ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                   write (95,87) dom(ib)%y(j)
            else
                   write (95,87) dom(ib)%yc(j)
            endif
   
            enddo
   
            message = 'Z_COORDINATES '//trim(adjustl(nz_str))//' float'
            write(95,'(a)') trim(adjustl(message))
            do k = 1, ntotk
            if ((flag.eq.3).or.(flag.eq.-3)) then
                   ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                   write (95,87) dom(ib)%z(k)
            else
                   write (95,87) dom(ib)%zc(k)
            endif
            enddo
   
            ! lets get table t (elements and node inidexes)
            tf = (ntoti-1)*(ntotj-1)*(ntotk-1)
            !tf = (ntoti-1)*(ntotk-1)
            write(tf_str,'(i8)') tf
            message = 'CELL_DATA '//trim(adjustl(tf_str))
            !message = trim(adjustl(message))
            !write(95,*) message
            write(95,'(a)') trim(adjustl(message))


            pf = ntoti*ntotj*ntotk
            write(pf_str,'(i8)') pf
           message = 'POINT_DATA '//trim(adjustl(pf_str))
            !message = trim(adjustl(message))
            !write(95,*) message
           write(95,'(a)') trim(adjustl(message))
   
           if (flag.eq.1) then
           message = 'SCALARS u_velocity float '//trim(adjustl('1')) ! '1' inidicates the number of components
           elseif(flag.eq.2) then
           message = 'SCALARS v_velocity float '//trim(adjustl('1'))
           elseif(flag.eq.3) then
           message = 'SCALARS w_velocity float '//trim(adjustl('1'))
           elseif(flag.eq.4) then
           message = 'SCALARS p_pressure float '//trim(adjustl('1'))
           elseif(flag.eq.5) then
           message = 'SCALARS mu_viscosity float '//trim(adjustl('1'))
           elseif(flag.eq.6) then
           message = 'SCALARS vis_viscosity float'//trim(adjustl('1'))
           elseif(flag.eq.7) then
           message = 'SCALARS LSM float '//trim(adjustl('1'))
           elseif(flag.eq.8) then
           message = 'SCALARS dens_density float '//trim(adjustl('1'))
           elseif (flag.eq.-1) then
           message = 'SCALARS uoo_velocity float '//trim(adjustl('1')) ! '1' inidicates the number of components
           elseif(flag.eq.-2) then
           message = 'SCALARS voo_velocity float '//trim(adjustl('1'))
           elseif(flag.eq.-3) then
           message = 'SCALARS woo_velocity float '//trim(adjustl('1'))
           elseif(flag.eq.-4) then
           message = 'SCALARS p_o_pressure float '//trim(adjustl('1'))
           endif 
   
            !message = trim(adjustl(message))
            !write(95,*) message
           write(95,'(a)') trim(adjustl(message))
   
           message = 'LOOKUP_TABLE default'
            !message = trim(adjustl(message))
            !write(95,*) message
           write(95,'(a)') trim(adjustl(message))
   
           ! lets write the field data
           !do k = 1, pf    
           !        write (95,87) dom(ib)%wstar(i,j,k) ! we can not just use wstar here since it does not match pf in the number of rows
           !enddo

         ! Alternative, this way we do not need matrix p
         do k = 1,ntotk, itotk	
         do j = 1, ntotj, itotj
         !j = 8 ! lets just take the cross-section along the centerline
         do i = 1, ntoti,itoti
         if (flag.eq.1) then
                ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                write (95,87) dom(ib)%u(i,j,k)
        elseif(flag.eq.2) then
                write (95,87) dom(ib)%v(i,j,k)
        elseif(flag.eq.3) then
                write (95,87) dom(ib)%w(i,j,k)
        elseif(flag.eq.4) then
                write (95,87) dom(ib)%p(i,j,k)
        elseif(flag.eq.5) then
                write (95,87) dom(ib)%mu(i,j,k)
        elseif(flag.eq.6) then
                write (95,87) dom(ib)%vis(i,j,k)
        elseif(flag.eq.7) then
                write (95,87) dom(ib)%phi(i,j,k)!,dom(ib)%phim(i,j,k),
!     & dom(ib)%dens(i,j,k),dom(ib)%mu(i,j,k)
        elseif(flag.eq.8) then
            write (95,87) dom(ib)%dens(i,j,k)
         elseif (flag.eq.-1) then
                write (95,87) dom(ib)%uoo(i,j,k)
        elseif(flag.eq.-2) then
                write (95,87) dom(ib)%voo(i,j,k)
        elseif(flag.eq.-3) then
                write (95,87) dom(ib)%woo(i,j,k)
        elseif(flag.eq.-4) then
                write (95,87) dom(ib)%p(i,j,k)
        endif 

         end do
         end do
         end do


           close(95)
                  
                enddo
!##########################################################################
        end !subroutine save2vtk_rectilinearGrid
!##########################################################################

!##########################################################################
        subroutine save2vtk_rectilinearPlane(ki,flag,cellFlag)!,flag)
!##########################################################################
                use multidata
                use vars
                use mpi
                implicit none
        
                integer, intent(in):: flag ! Indicates fields to write: u, v, w, p
                integer, intent(in):: cellFlag ! indicate the y cell id for cryating a crossection plane
                character*7 :: fileName ! file name
                character*70 :: heather ! Heather in the vtk file
        
                !integer, intent(in):: tf,pf,nne
                integer:: tf,pf,nne
                ! tf: number of rows in matrix t
                ! pf: number of rows in matrix p
                ! nne: number of points to read in each row of the element table
                !p: table of node coordinates
                !t: table element (rows) and their node inidexes (columns)
                !e: table of nodes per side
        
                ! change this to global variables
                DOUBLE PRECISION,allocatable,dimension(:,:):: p, t
        
                integer :: sn,i,j,k,toti,totj,totk
                integer :: is,ie,js,je,ks,ke
                integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
                ! ki should have intent(in)
                character*8 :: chb,chb2,nx_str,ny_str,nz_str
                character*35 :: gf ! file name
                character*35 :: message,message2,message3
                integer:: nds,n1,n2,n3,n4
                character*8 :: pf_str, tf_str, nds_str

                do ib = 1, nbp ! lets sweep only the subdomains of each procesor
        if (flag.eq.1) then
                fileName = 'vtk_upl'
                heather ='2D U Velocity Distribution'
        elseif(flag.eq.2) then
                fileName = 'vtk_vpl'
                heather ='2D V Velocity Distribution'
        elseif(flag.eq.3) then
               fileName = 'vtk_wpl'
               heather ='2D W Velocity Distribution'
        elseif(flag.eq.4) then
                fileName = 'vtk_ppl'
                heather ='2D P Pressure Distribution'
        elseif(flag.eq.5) then
                fileName = 'vtk_npl'
                heather ='2D mu Distribution'
        elseif(flag.eq.6) then
                fileName = 'vtkntpl'
                heather ='2D vis Distribution'
        elseif(flag.eq.7) then
                fileName = 'vtkphpl'
                heather ='2D phi LSM Distribution'
        elseif(flag.eq.8) then
                fileName = 'vtkropl'
                heather ='2D rho Density Distribution'
        elseif (flag.eq.-1) then
                fileName = 'vtkoupl'
                heather ='2D Uoo Velocity Distribution'
        elseif(flag.eq.-2) then
                fileName = 'vtkovpl'
                heather ='2D Voo Velocity Distribution'
        elseif(flag.eq.-3) then
               fileName = 'vtkowpl'
               heather ='2D Woo Velocity Distribution'
        elseif(flag.eq.-4) then
               fileName = 'vtkoppl'
               heather ='2D P_o Pressure Distribution'
        endif 

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))

        gf=trim(adjustl(chb))//fileName//
     &		'.'//trim(adjustl(chb2))//'.vtk'

        !open (unit=95, file=gf)
        ! CSV format
87      format(1x, *(g0, " ")) 
        OPEN(unit = 95, access = "sequential", action = "write",
     &   status = "replace", file = gf, form = "formatted") 

!        is=pl+1; ie=dom(ib)%ttc_i-pl
!        js=pl+1; je=dom(ib)%ttc_j-pl
!        ks=pl+1; ke=dom(ib)%ttc_k-pl
!        toti=ie-(is-1)+1	!toti=nicell
!        totj=je-(js-1)+1
!        totk=ke-(ks-1)+1

!Plot every X points in every direction: change itotX respectively
!	   itoti=1 	 ntoti=toti/itoti+1
!	   itotj=1 	; ntotj=totj/itotj+1
!	   itotk=1 	; ntotk=totk/itotk+1
	
	   itoti=1 	; ntoti=dom(ib)%ttc_i
	   itotj=1 	; ntotj = 1 !ntotj=dom(ib)%ttc_j
           itotk=1 	; ntotk=dom(ib)%ttc_k


         !write(95,*) '# vtk DataFile Version 3.0'
           write(95,'(a)') '# vtk DataFile Version 3.0'
           ! check
            !write(95,*) trim(adjustl(heather)) !'2D Temperature Distribution\n'
            write(95,'(a)') adjustl(trim(adjustl(heather))) !'2D Temperature Distribution\n'
   
            write(95,'(a)') 'ASCII'
   
            write(95,'(a)') 'DATASET RECTILINEAR_GRID'
   
   
            write(nx_str,'(i8)') ntoti
            write(ny_str,'(i8)') ntotj
            write(nz_str,'(i8)') ntotk
            message = 'DIMENSIONS '//adjustl(trim(adjustl(nx_str)))//' '
     &  //adjustl(trim(adjustl(ny_str)))//' '
     &  //adjustl(trim(adjustl(nz_str)))
            write(95,'(a)') trim(adjustl(message))
   
   
            message = 'X_COORDINATES '//trim(adjustl(nx_str))//' float'
            write(95,'(a)') trim(adjustl(message))
            do i = 1, ntoti,itoti
            if ((flag.eq.1).or.(flag.eq.-1)) then
                   ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                   write (95,87) dom(ib)%x(i)
            else
                   write (95,87) dom(ib)%xc(i)
            endif 
            enddo

            message = 'Y_COORDINATES '//trim(adjustl(ny_str))//' float'
            write(95,'(a)') trim(adjustl(message))
            !do j = 1, ntotj
            j = cellFlag
            if ((flag.eq.2).or.(flag.eq.-2)) then
                   ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                   write (95,87) dom(ib)%y(j)
            else
                   write (95,87) dom(ib)%yc(j)
            endif
   
            !enddo
   
            message = 'Z_COORDINATES '//trim(adjustl(nz_str))//' float'
            write(95,'(a)') trim(adjustl(message))
            do k = 1, ntotk
            if ((flag.eq.3).or.(flag.eq.-3)) then
                   ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                   write (95,87) dom(ib)%z(k)
            else
                   write (95,87) dom(ib)%zc(k)
            endif
            enddo
   
            ! lets get table t (elements and node inidexes)
            !tf = (ntoti-1)*(ntotj-1)*(ntotk-1)
            tf = (ntoti-1)*(ntotk-1)
            write(tf_str,'(i8)') tf
            message = 'CELL_DATA '//trim(adjustl(tf_str))
            !message = trim(adjustl(message))
            !write(95,*) message
            write(95,'(a)') trim(adjustl(message))


            pf = ntoti*ntotj*ntotk
            write(pf_str,'(i8)') pf
           message = 'POINT_DATA '//trim(adjustl(pf_str))
            !message = trim(adjustl(message))
            !write(95,*) message
           write(95,'(a)') trim(adjustl(message))
   
           if (flag.eq.1) then
           message = 'SCALARS u_velocity float '//trim(adjustl('1')) ! '1' inidicates the number of components
           elseif(flag.eq.2) then
           message = 'SCALARS v_velocity float '//trim(adjustl('1'))
           elseif(flag.eq.3) then
           message = 'SCALARS w_velocity float '//trim(adjustl('1'))
           elseif(flag.eq.4) then
           message = 'SCALARS p_pressure float '//trim(adjustl('1'))
           elseif(flag.eq.5) then
           message = 'SCALARS mu_viscosity float '//trim(adjustl('1'))
           elseif(flag.eq.6) then
           message = 'SCALARS vis_viscosity float'//trim(adjustl('1'))
           elseif(flag.eq.7) then
           message = 'SCALARS LSM float '//trim(adjustl('1'))
           elseif(flag.eq.8) then
           message = 'SCALARS dens_density float '//trim(adjustl('1'))
           elseif (flag.eq.-1) then
           message = 'SCALARS uoo_velocity float '//trim(adjustl('1')) ! '1' inidicates the number of components
           elseif(flag.eq.-2) then
           message = 'SCALARS voo_velocity float '//trim(adjustl('1'))
           elseif(flag.eq.-3) then
           message = 'SCALARS woo_velocity float '//trim(adjustl('1'))
           elseif(flag.eq.-4) then
           message = 'SCALARS p_o_pressure float '//trim(adjustl('1'))
           endif 
   
            !message = trim(adjustl(message))
            !write(95,*) message
           write(95,'(a)') trim(adjustl(message))
   
           message = 'LOOKUP_TABLE default'
            !message = trim(adjustl(message))
            !write(95,*) message
           write(95,'(a)') trim(adjustl(message))
   
           ! lets write the field data
           !do k = 1, pf    
           !        write (95,87) dom(ib)%wstar(i,j,k) ! we can not just use wstar here since it does not match pf in the number of rows
           !enddo

         ! Alternative, this way we do not need matrix p
         do k = 1,ntotk, itotk	
         !do j = 1, ntotj, itotj
         j = cellFlag ! lets just take the cross-section along the centerline
         do i = 1, ntoti,itoti
         if (flag.eq.1) then
                ! check, there may by a trailing space that may cause isues for the vtk reader of paraview
                write (95,87) dom(ib)%u(i,j,k)
        elseif(flag.eq.2) then
                write (95,87) dom(ib)%v(i,j,k)
        elseif(flag.eq.3) then
                write (95,87) dom(ib)%w(i,j,k)
        elseif(flag.eq.4) then
                write (95,87) dom(ib)%p(i,j,k)
        elseif(flag.eq.5) then
                write (95,87) dom(ib)%mu(i,j,k)
        elseif(flag.eq.6) then
                write (95,87) dom(ib)%vis(i,j,k)
        elseif(flag.eq.7) then
                write (95,87) dom(ib)%phi(i,j,k)!,dom(ib)%phim(i,j,k),
!     & dom(ib)%dens(i,j,k),dom(ib)%mu(i,j,k)
        elseif(flag.eq.8) then
                write (95,87) dom(ib)%dens(i,j,k)
         elseif (flag.eq.-1) then
                write (95,87) dom(ib)%uoo(i,j,k)
        elseif(flag.eq.-2) then
                write (95,87) dom(ib)%voo(i,j,k)
        elseif(flag.eq.-3) then
                write (95,87) dom(ib)%woo(i,j,k)
        elseif(flag.eq.-4) then
                write (95,87) dom(ib)%p(i,j,k)
        endif 

         end do
         !end do
         end do


           close(95)
                  
                enddo
!##########################################################################
        end !subroutine save2vtk_rectilinearPlane
!##########################################################################

!##########################################################################
        subroutine tec_inst_plane(ki)
!##########################################################################
        use multidata
        use vars
	  use mpi
        implicit none
        integer :: sn,i,j,k,toti,totj,totk,is,ie,js,je,ks,ke,ki,ib
        integer :: itoti,itotj,itotk,ntoti,ntotj,ntotk
	  double precision :: xval,yval,zval
	  integer :: iplane
        character*8 :: chb,chb2
        character*35 :: gf

        do ib=1,nbp

        is=pl+1; ie=dom(ib)%ttc_i-pl
        js=pl+1; je=dom(ib)%ttc_j-pl
        ks=pl+1; ke=dom(ib)%ttc_k-pl
        toti=ie-(is-1)+1 ; itoti=1 	; ntoti=toti/itoti+1
        totj=je-(js-1)+1 ; itotj=1 	; ntotj=totj/itotj+1
        totk=ke-(ks-1)+1 ; itotk=1 	; ntotk=totk/itotk+1
	  
!!!! Set the plane orientation and location at xval, yval or zval
	  iplane=2 			!1=x-plane, 2=y-plane, 3=z-plane
	  xval=(xen-xst)/2 ;  yval=(yen-yst)/2  ;  zval=(zen-zst)/2 

	 if(iplane.eq.1) then
        do i=is-1,ie
	    if(dom(ib)%x(i).le.xval .and. dom(ib)%x(i+1).gt.xval) then

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecplane'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
           open (unit=95, file=gf)     
        write(95,*)'variables=y,z,u,v,w'
        write(95,89)'zone i= 1, j=',ntotj,', k=',ntotk,' f=point'		
           do k=ks-1,ke+1,itotk	!k=pl,nicell+pl
             do j=js-1,je+1,itotj
             write (95,88) dom(ib)%y(j),dom(ib)%z(k),
     &  dom(ib)%u(i,j,k),dom(ib)%v(i,j,k),dom(ib)%w(i,j,k)
             end do
           end do
           close (95)
	    endif ! within the interval
	  enddo
	 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	 if(iplane.eq.2) then
        do j=js-1,je
	    if(dom(ib)%y(j).le.yval .and. dom(ib)%y(j+1).gt.yval) then

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecplane'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
           open (unit=95, file=gf)     
        write(95,*)'variables=x,z,u,v,w'
        write(95,89)'zone i=',ntoti,', j=1 , k=',ntotk,' f=point'		
           do k=ks-1,ke+1,itotk	
             do i=is-1,ie+1,itoti
             write (95,88) dom(ib)%x(i),dom(ib)%z(k),
     &  dom(ib)%u(i,j,k),dom(ib)%v(i,j,k),dom(ib)%w(i,j,k)
!     &  ,dom(ib)%phi(i,j,k)
             end do
           end do
           close (95)
	    endif ! within the interval
	  enddo
	 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	 if(iplane.eq.3) then
        do k=ks-1,ke
	    if(dom(ib)%z(k).le.zval .and. dom(ib)%z(k+1).gt.zval) then

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecplane'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
           open (unit=95, file=gf)     
        write(95,*)'variables=x,y,u,v,w'
        write(95,89)'zone i=',ntoti,', j=',ntotj,', k=1 f=point'		
           do j=js-1,je+1,itotj
            do i=is-1,ie+1,itoti
             write (95,88) dom(ib)%x(i),dom(ib)%y(j),
     &  dom(ib)%u(i,j,k),dom(ib)%v(i,j,k),dom(ib)%w(i,j,k)
             end do
           end do
           close (95)
	    endif ! within the interval
	  enddo
	 endif

	enddo !ib

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end
