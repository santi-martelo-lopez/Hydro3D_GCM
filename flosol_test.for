!##########################################################################
        subroutine flosol
!##########################################################################
        use vars
        use mpi
        use multidata
	  use vars_pt
        implicit none
        real *8:: wtimedum,wtime_total,wtime_solver,wtime_ib,wtime_cd
        real*8:: wtime_gh
	  real *8:: wtime_lsm
        integer ib,i,j,k,kutta,kuttacond,jtime,ii
        double precision :: alfark(3)

! Set some constants ---------------------
! ..... 3-STEP RUNGE KUTTA
        if (conv_sch.eq.4) then
           alfark(1)=1./3.
           alfark(2)=0.5
           alfark(3)=1.0
           kuttacond=3
! ..... 2-STEP RUNGE KUTTA
        else if (conv_sch.eq.3) then
           alfark(1)=0.5
           alfark(2)=1.0
           alfark(3)=1.0
           kuttacond=2
        end if

        alfabc = 1
        alfapr = 1.0

        if (LRESTART.eq..false.) cnt_pt = 1 

        numfile1=1002; numfile2=1003; numfile3=1004
	  numfile4=1005; numfile5=1006
        if(myrank.eq.0) then
           if (pressureforce) then
            open (unit=numfile1, file='forcn.dat')
            write (numfile1,*)'variables=time,forcn,qsttp,flwsum'
           end if
           if (pressureforce_y) then
            open (unit=numfile4, file='forcn_y.dat')
            write (numfile4,*)'variables=time,forcn_y,qsttp_y,flwsum_y,'
           end if
           if (pressureforce_z) then
            open (unit=numfile5, file='forcn_z.dat')
            write (numfile5,*)'variables=time,forcn_z,qsttp_z,flwsum_z,'
           end if

           open (unit=numfile2, file='rms.dat')
           write (numfile2,*)'variables=iter,time,rmax,dt,Mdef'
           if (L_LSM) then
            open (unit=numfile3, file='normV.dat')
            write (numfile3,*)'variables=ntime,normV,steps,ctime,dt'
	     end if
        end if

        itime_start=ntime+1
	  jtime = 0
	  ireadinlet = 0
	  iaddinlet = 1

	if (myrank.eq.0) then
	   open(unit=203, file='worktime.dat')
	   write(203,*)'Variables=it,dt,CFL,C-D,PSolver,IBM,LSM,Total' !Add LPT if needed
	endif

!================== Start Time Loop ====================================
        do itime=itime_start,itime_end

	     if (myrank.eq.0) wtime_total = MPI_WTIME ( )
!Calculate time step
           call checkdt

           ctime=ctime+dt  ;  ntime = ntime + 1

!Reading inflow data, bc_west = 7
	   if (read_inflow.eq..true.) then  
           		if (ireadinlet.eq.ITMAX_PI.and.iaddinlet.eq.1) then
				iaddinlet=-1
          		elseif (ireadinlet.eq.1.and.iaddinlet.eq.-1) then
				iaddinlet=1
	     		endif
          		ireadinlet=ireadinlet+iaddinlet
!Reading inflow data from SEM files 
	   elseif ((bc_w.eq.8)) then 
           		if (ireadinlet.eq.ITMAX_SEM.and.iaddinlet.eq.1) then 
				iaddinlet=-1 
          		elseif (ireadinlet.eq.1.and.iaddinlet.eq.-1) then 
				iaddinlet=1 
	     		endif 
          		ireadinlet=ireadinlet+iaddinlet 
	   endif
!Internal mapping 
		if(bc_w.eq.77)	   call write_mapping
!-----------------------------------------------------------------------

           if (LENERGY) call boundT
           if (LSCALAR) call boundS
           !
           !
           ! ghost cell method
           !Store variables from previous time step:
           do ib=1,nbp
              do k=1,dom(ib)%ttc_k
              do i=1,dom(ib)%ttc_i
              do j=1,dom(ib)%ttc_j
                 dom(ib)%uoo(i,j,k)=dom(ib)%u(i,j,k)
                 dom(ib)%voo(i,j,k)=dom(ib)%v(i,j,k)
                 dom(ib)%woo(i,j,k)=dom(ib)%w(i,j,k)
           		if (LENERGY) dom(ib)%To(i,j,k)=dom(ib)%T(i,j,k)
	     	      if (LSCALAR) dom(ib)%So(i,j,k)=dom(ib)%S(i,j,k)
              end do
              end do
              end do
           end do
         !
         !
         !
           !Calculate energy equation (e.g. temperature)
           if (LENERGY) call energy
!Calculate sediment transport
           if (LSCALAR) call sediment_4thtest
!Calculate free-surface using level set method
           if (L_LSM)  then
			CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
		      if (myrank.eq.0) wtime_lsm = MPI_WTIME ( )
             call LSM_3D
             call heaviside 
			CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
			if (myrank.eq.0) wtime_lsm = MPI_WTIME( )-wtime_lsm
           end if
!Initial calculations in IBM to allocate particles to domains and cpus
           if (LIMB)  call IB_previous ! Here the location of the body and delta functions are updated
           !
           ! ghost cell method
           !
           ! If the solid can move, we need to serach for new cell at each time step
           if (LIMB)  then
            !if(rotating(M).eq..FALSE.) then
            !   if (itime.eq.itime_start) then
            !      ! If the solid is fixed, we only need to search for the solid cells at the first time step, 
            !      ! additionally we do not need to check for emerging cells
            !      call searchGhostCells ! Searches for ghost cells, if it finds any, the flag mucx = 1,2,3 or 4
            !   endif
            !elseif(rotating(M).eq..TRUE.) then
            !if (itime.eq.itime_start)then
               call searchGhostCells
               !After that, we nned to check for emerging cells, new solid cells, or update the info on the solid cells
               call checkEmergingCells ! it checks if mucx = mucxold ?
            !      ! if False and mucx = 0, we have found an emerging cell
            !      ! if False and mucxold = 0, we have found a new fluid cell inside the inmersed body
            !      ! mucxold = mucx
            !      ! mucx = 0
               call ghostCell(0) ! updating solid cells after the movement of the solid
               ! Write forces at play in the ghost cell method
               call ghostCell_forces
            if ((mod(itime,n_out).eq.0).and.(itime.ge.itime_start)) then
                  call save2vtk_rectilinearGrid(itime,-1)   ! -1 : flag for writting uoo velocities
             !     call save2vtk_rectilinearGrid(itime,-2)   ! -2 : flag for writting voo velocities
                  call save2vtk_rectilinearGrid(itime,-3)   ! -3 : flag for writting woo velocities
                  call save2vtk_rectilinearGrid(itime,-4)  ! -4 : flag for writting p_o pressure
               endif
            !
            !endif
           endif
           !
!Release and alloscation of particles in LPT
	     if (LPT) then
		if (mod(itime,tsnr).eq.0 .and. myrank.eq.0) call release_pt
		if (myrank.eq.0) call alloc_pt	
		call MPI_pt
	     endif
!Calculation of the sub-grid scale viscosity term
           if(SGS) then
              if(sgs_model.eq.1) then
                 call eddyv_smag
              else if(sgs_model.eq.2) then
                 call eddyv_wale
              else if(sgs_model.eq.3) then
                 call eddyv_1eqn
              else if(sgs_model.eq.4) then
                 call eddyv_keps
              end if
           end if
           !
! When streamwise periodic boundary conditions, the inlet velocities are
! corrected with the pressure gradient
           if (PERIODIC) call pressure_forcing
           !
! Resolve the convection and diffusion components
           CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
   ! Forces all taks in communiator to synchronize
   ! A tasks enters the barrier and waits for all other tasks to reach it
   ! Deadlock if not all required tasks reach it
	if (myrank.eq.0) wtime_cd = MPI_WTIME ( )
           if(conv_sch.eq.3 .or. conv_sch.eq.4) then
            do kutta=1,kuttacond
               if (kutta.eq.1) then
                  alfabc = 1
               else
                  alfabc = 0
               endif
               alfapr=alfark(kutta)
               select case (differencing)
                  case (1) 
                     call rungek_conv2nd
                  case (2) 
                     call rungek_conv4th
                  case (3)
                     call rungek_convWENO
               end select
               !
! Settig global velocities to zero after solving the Pressure Poisson Equation
               call ghostCell(1)
               !
               if (diff_sch.eq.3) then
                  select case (differencing)
                     case (1) 
                        call rungek_diff2nd
                     case (2) 
                        call rungek_diff4th
                     case (3)
                        call rungek_diff2nd
                        !              call ghostCell(1)
                  end select
                  !Write results 
                           !IF ((mod(itime,n_out).eq.0).and.(itime.ge.itime_start)) then           
                           ! csv format
                           !call csv_instant_ustar(itime,1) ! 1 to indicate that is right before applying the ghost method
                           !call csv_instant_vstar(itime,1)
                           !call csv_instant_wstar(itime,1)
                           !endif
                  ! Settig global velocities to zero after solving the Pressure Poisson Equation
                           call ghostCell(1)
                  !Write results
                           ! IF ((mod(itime,n_out).eq.0).and.(itime.ge.itime_start)) then           
                              ! csv format
                              !call csv_instant_ustar(itime,2) ! 2 to indicate that is right after applying the ghost method
                              !call csv_instant_vstar(itime,2)
                              !call csv_instant_wstar(itime,2)
                              ! vtk format
                           !call save2vtk_rectilinearGrif(itime,1) ! 1 : flag for writting u velocities
                           !call save2vtk_rectilinearGrif(itime,2) ! 2 : flag for writting v velocities
                           !call save2vtk_rectilinearGrif(itime,3) ! 3 : flag for writting w velocities
                           ! endif
               else
                  call diffusion
                  write(*,*) "diffusion"
               end if
               if (kutta.lt.kuttacond) call calvel
            end do
           else
              call convection
              write(*,*) "convection"
! Settig global velocities to zero after solving the Pressure Poisson Equation
!              call ghostCell(1)
              
              call diffusion
              write(*,*) "diffusion 2"
           end if

!Write results 
           !IF ((mod(itime,n_out).eq.0).and.(itime.ge.itime_start)) then           
           ! csv format
           !call csv_instant_ustar(itime,1) ! 1 to indicate that is right before applying the ghost method
           !call csv_instant_vstar(itime,1)
           !call csv_instant_wstar(itime,1)
           !endif
! Settig global velocities to zero after solving the Pressure Poisson Equation
           !call ghostCell(1)
!Write results
          ! IF ((mod(itime,n_out).eq.0).and.(itime.ge.itime_start)) then           
            ! csv format
            !call csv_instant_ustar(itime,2) ! 2 to indicate that is right after applying the ghost method
            !call csv_instant_vstar(itime,2)
            !call csv_instant_wstar(itime,2)
            ! vtk format
         !call save2vtk_rectilinearGrif(itime,1) ! 1 : flag for writting u velocities
         !call save2vtk_rectilinearGrif(itime,2) ! 2 : flag for writting v velocities
         !call save2vtk_rectilinearGrif(itime,3) ! 3 : flag for writting w velocities
           ! endif

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if (myrank.eq.0) wtime_cd=MPI_WTIME( )-wtime_cd

!Lagrangian particle tracking:
           	if (LPT) then								
        		call exchange(11)
        		call exchange(22)
        		call exchange(33)
			if (np_loc.gt.0) call particle_tracking			!Procs without particles do not enter
			call final_LPT
		endif
!Immersed boundary method
         IF (LIMB) then
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	     if (myrank.eq.0) wtime_ib = MPI_WTIME ( ) 
 			call IBM					
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	     if (myrank.eq.0) wtime_ib = MPI_WTIME ( ) - wtime_ib
	   ENDIF
!Correction of the outlet velocity in order to guarantee mass conservation
           call correctoutflux

!Resolve Poisson pressure correction equation
	     if (myrank.eq.0) wtime_solver = MPI_WTIME ( )
           if(solver.eq.1) then
              call pressure_1sweep	!SIP solver
           else if(solver.eq.2) then
              call newsolv_mg		!Multi-grid
           end if
	     if (myrank.eq.0) wtime_solver = MPI_WTIME ( ) - wtime_solver

! Settig global velocities to zero after solving the Pressure Poisson Equation
        !call ghostCell(2)

!!2nd iteration of the PPE in order to enforce the mass conservation equation
!!Resolve Poisson pressure correction equation
!	     if (myrank.eq.0) wtime_solver = MPI_WTIME ( )
!           if(solver.eq.1) then
!              call pressure_1sweep	!SIP solver
!           else if(solver.eq.2) then
!              call newsolv_mg		!Multi-grid
!           end if
!	     if (myrank.eq.0) wtime_solver = MPI_WTIME ( ) - wtime_solver
!
!! Settig global velocities to zero after solving the Pressure Poisson Equation
!        call ghostCell(2)
!
!!3rd iteration of the PPE in order to enforce the mass conservation equation
!!Resolve Poisson pressure correction equation
!	     if (myrank.eq.0) wtime_solver = MPI_WTIME ( )
!           if(solver.eq.1) then
!              call pressure_1sweep	!SIP solver
!           else if(solver.eq.2) then
!              call newsolv_mg		!Multi-grid
!           end if
!	     if (myrank.eq.0) wtime_solver = MPI_WTIME ( ) - wtime_solver
!
!! Settig global velocities to zero after solving the Pressure Poisson Equation
!        call ghostCell(2)
        
       ! ! obtaining pressure force around the body
       ! call pressureProve
       ! call shearStressProve

!Calculate averaged variables:
        if(time_averaging.and.
     & 		ctime.ge.t_start_averaging1)	call update_mean

!Write time series in selected points:
	  if (ctime.ge.t_start_averaging1) call timesig				

!Write outflow planes
	  if ((save_inflow.eq..true.) .and. itime.ge.itime_start)
     &     call write_inflow		       !set itime since when you'd like to write inflow planes

! Write forces at play in the ghost cell method
!      call ghostCell_forces

!Write results in tecplot format 
       IF ((mod(itime,n_out).eq.0).and.(itime.ge.itime_start)) then

	  	  IF (LTURB.EQ..TRUE.)     call tec_turb(itime)
!	  	  IF (LINST.EQ..TRUE.)     call tec_instant(itime)
	  	  IF (LINST.EQ..TRUE.) then     
			!call tec_instant_u(itime)
			!call tec_instant_v(itime)
         !call tec_instant_w(itime)
         
         ! csv format
			!call csv_instant_u(itime)
			!call csv_instant_v(itime)
         !call csv_instant_w(itime)

         !call save2vtk(itime,1) ! 1 : flag for writting u velocities
         !call save2vtk(itime,2) ! 2 : flag for writting v velocities
         !call save2vtk(itime,3) ! 3 : flag for writting w velocities
         !call save2vtk(itime,4) ! 4 : flag for writting p pressure values
         !call save2vtk_rectilinearGrid(itime,1) ! 1 : flag for writting u velocities
         !call save2vtk_rectilinearGrid(itime,2) ! 2 : flag for writting v velocities
         !call save2vtk_rectilinearGrid(itime,3) ! 3 : flag for writting w velocities

	 	  endif
		  IF (LPLAN.EQ..TRUE.)     call tec_inst_plane(itime)
        
        !IF (LTECP.EQ..TRUE.)     call tecplot_p(itime)
        !IF (LTECP.EQ..TRUE.)     call csvplot_p(itime)
        IF (LTECP.EQ..TRUE.) then
        !call save2vtk_rectilinearGrid(itime,4) ! 4 : flag for writting p pressure values
        endif


	  	  IF (LTECBIN.EQ..TRUE.)   call tecbin(itime)
		  !IF (L_LSM.EQ..TRUE.)     call tecplot_phi(itime)
        IF (L_LSM.EQ..TRUE.) then
   !      call save2vtk_rectilinearGrid(itime,7) ! 4 : flag for writting phi pressure values
         endif
        !IF (L_LSMbase.EQ..TRUE.) call tecplot_phi(itime)	!!
        IF (L_LSMbase.EQ..TRUE.) then
   !      call save2vtk_rectilinearGrid(itime,7) ! 4 : flag for writting phi pressure values
         endif
        IF (LENERGY.EQ..TRUE.)   call tecplot_T(itime)	!!
	 	  IF (LSCALAR.EQ..TRUE.)   call tecplot_S(itime)	!!
!File containing information needed for restarting
	  	  IF (LTECBIN.EQ..TRUE.)   then
         open (unit=101, file='final_ctime.dat')
         IF(myrank.eq.0) write(101,'(i8,3F15.6)')ntime,ctime,forcn,qstpn
         close(101)
		  ENDIF
!Lagrangian Particle Tracking final checks
	       if (LPT.eq..TRUE. .AND. myrank.eq.0) then          			
	    		open(30,file='final_particle.dat')
	   	 	write(30,*) np
	    		write(30,*) cnt_pt
	    		do k=1,np
	      	 write(30,*) xp_pt(k),yp_pt(k),zp_pt(k)
	      	 write(30,*) uop_pt(k),vop_pt(k),wop_pt(k)  
	 		 write(30,*) dp_pt(k)
	    		end do
	    		close(30)
		  endif
        ENDIF
        if (LPT) then									
		if ((mod(itime,tsteps_pt).eq.0).and.
     &		(itime.gt.itime_start)) then 			
 	           	if (myrank.eq.0) call TECPARTICLE(cnt_pt)
      	     	call TECPLOT(cnt_pt)
			cnt_pt = cnt_pt + 1
		endif
        end if

!Write on screen and output.dat:
           if (MOD(itime,20).eq.0 .and. myrank.eq.0)  then
              wtimedum = MPI_WTIME ( ) - wtime
              write (6,5000) itime,ctime,dt,dtavg
              write (6,5500) wtimedum
              write (6,*) ' '
              write (numfile,5000) itime,ctime,dt,dtavg
              write (numfile,*) ' '
              write (numfile,5500) wtimedum
              write (numfile,*) ' '
           end if

!OUTPUT THE TIME SPENT AT DIFFERENT SUBROUTINES  ---------------------
	 if (myrank.eq.0) then
	   wtime_total = MPI_WTIME ( ) - wtime_total
	   write(203,'(i6,1f12.7,6f12.5)')itime,dt,wtime_cfl
     &     ,wtime_cd,wtime_solver,wtime_ib,wtime_lsm,wtime_total
!STOP THE CODE IF:
	    if(wtime_cfl.gt.2.99 .AND. itime.gt.(itime_start+4)) then
		write(6,*) 'Code stopped due to LARGE CFL',wtime_cfl
		stop
	    endif
	 endif
!
 ! End Time Loop ---------------------------
!
       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        itime = itime - 1

	  close(203) !worktime.dat
!Output final results:
        if (mod(itime,n_out).ne.0) then
	  	  IF (LTURB.EQ..TRUE.)   call tec_turb(itime)
	  	  IF (LTECBIN.EQ..TRUE.) call tecbin(itime)
		  IF (L_LSM.EQ..TRUE.)   call tecplot_phi(itime)
		  IF (L_LSMbase.EQ..TRUE.)   call tecplot_phi(itime)	!!
	 	  IF (LENERGY.EQ..TRUE.)   call tecplot_T(itime)	!!
	 	  IF (LSCALAR.EQ..TRUE.)   call tecplot_S(itime)	!!
		  if (ctime.ge.t_start_averaging2)  call timesig
           open (unit=101, file='final_ctime.dat')
           if(myrank.eq.0) write (101,'(i8,3F15.6)') 
     &   ntime,ctime,forcn,qstpn
           close(101)
        end if

	if (myrank.eq.0) write (6,*) 'ctime=' , ctime
	if (myrank.eq.0) write (numfile,*) 'ctime=' , ctime

5000  format(/1x,10(1h=),' nrtstp=',i8,2x,'ctime=',e14.6,2x,
     & 'dt=',e14.6,'  dtavg=',e14.6)
5500  format(/1x,'Work took ',f18.8,2x,' seconds')
        end subroutine flosol
!##########################################################################

