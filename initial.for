!##########################################################################
        subroutine read_control
!##########################################################################
        use vars
        use multidata
        use mpi
        implicit none
        integer mgi,mgj,mgk,pow2,ib,i,j

        open (unit=12, file='control.cin')
!------DOMAIN SIZE AND DISCRETIZATION -------------------------------------    
	  read (12,*) 
	  read (12,*) 
        read (12,*) keyword,ubulk
        read (12,*) g_dx,g_dy,g_dz
        read (12,*) Re,Pr,Sc_t,beta
        read (12,*) conv_sch
        read (12,*) diff_sch
        read (12,*) differencing
        read (12,*) solver
        read (12,*) ngrid_input,mg_itrsch
        read (12,*) maxcy,irestr,iproln
        read (12,*) variTS,dt
        ! ramp function paramerts 
        read (12,*) tr_iter,dt_factor,maxcy_factor
        read (12,*) dt_rf,maxcy_rf
        read (12,*) uTrans
        read (12,*) lowerU,upperU
        ! tr_iter: number transition iterations to use default time step, dt; and maximum multigrid iterations, maxcy.
        ! dt_factor: maximum amplification factor for dt
        ! maxcy_factor: maximum amplification factor for maxcy
        ! dt_rf: type of rump function. 0: step, 1: linear
        ! maxcy_rf: type of rump function. 0: step, 1: linear
        !
        read (12,*) sweeps,safety_factor
        read (12,*) itime_end,LRESTART,reinitmean,n_out
        read (12,*) niter,eps,nswp(1),nswp(2),nswp(3),nswp(4)
	  read (12,*) 
        read (12,*) bc_w
        read (12,*) bc_e
        read (12,*) bc_s
        read (12,*) bc_n
        read (12,*) bc_b
        read (12,*) bc_t
	  read (12,*) L_n,fric
        read (12,*) save_inflow,ITMAX_PI
	  read (12,*) 
!===== Aristos 25/02/19 Waves ==========================================
        !theory = 2!Wave theory,(1=Linear, 2=2nd 0rder Stokes, 3=Solitary,11=cnoidal)
        !L_w = 2!L(wavelength)
        !omega = 0.5   ! omega: frequancy
        !H_w = 0.25   ! H_w: Ho(Wave Height)
        !tau = 0.0! tau(phase for solitary)
        !! length: initial water depth
        !u_prof = 1  ! u_prof: 1					current??? 
        !L_beach = 0 ! Wave Absorption Flag
        !beach   = 0 ! Wave Absorption(1=Artificial damping, 2=Relaxation method)
        !x_end_beach  = 0 ! end of domain
        read (12,*) theory !u_prof ! theory : Wave theory,(1=Linear, 2=2nd 0rder Stokes, 3=Solitary,11=cnoidal), u_prof : current??? 
        read (12,*) H_w,tau,L_w ! H_w: Ho(Wave Height), ! tau(phase for solitary), L(wavelength)
        read (12,*) L_beach,beach,x_end_beach ! Wave Absorption Flag, Wave Absorption(1=Artificial damping, 2=Relaxation method), end of domain
        read (12,*) 
!===== Aristos 25/02/19 Waves ==========================================
        read (12,*) UPROF		
        read (12,*) TI_SEM 
        read (12,*) ITMAX_SEM
	  read (12,*) 
        read (12,*) time_averaging,t_start_averaging1,
     & t_start_averaging2,noise
        read (12,*) SGS,sgs_model
        read (12,*) LMR
        read (12,*) pl_ex	
        read (12,*) LIMB,LROUGH,LSCALAR
        read (12,*) LPT,OMP_threads
	  read (12,*) 
        read (12,*) LENERGY,Th,Tc
        read (12,*) Tbc_w
        read (12,*) Tbc_e
        read (12,*) Tbc_s
        read (12,*) Tbc_n
        read (12,*) Tbc_b
        read (12,*) Tbc_t

! Starting to read parameters for the Level Set Method
        read (12,*)
        read (12,*) length,L_LSM,reinit	
	! length: inital depth of the free surface
	! L_LSM:
	! reinit:

        read (12,*) ntime_reinit,reldif_LSM,cfl_lsm
	  read (12,*) LENDS
	  read (12,*) L_LSMbase
	  read (12,*) L_anim_phi,L_anim_grd 
	  read (12,*) densl,densg,nul,nug
        read (12,*) grx,gry,grz
	  read (12,*) slope
	  read (12,*) 

! Starting to read parameters for creating costumed output files
        read (12,*) LTURB
        read (12,*) LINST	 
        read (12,*) LPLAN	 
        read (12,*) LTECP	                 	 
        read (12,*) LTECBIN	
	  read (12,*) 
	  read (12,*) n_unstpt
	  	allocate(id_unst(n_unstpt),i_unst(n_unstpt)
     &		  ,j_unst(n_unstpt),k_unst(n_unstpt))
		do i=1,n_unstpt
			read (12,*)id_unst(i),i_unst(i),j_unst(i),k_unst(i)
		enddo

	 if (bc_w.eq.7) read_inflow=.true.

        mul = nul * densl
        mug = nug * densg

	if (trim(L_n).eq.'n') fric=fric**0.33

!Pressure forcing for periodic conditions: PABLO 02/2018
            PERIODIC=.false. ; pressureforce=.false.
		pressureforce_y=.false. ; pressureforce_z=.false. 

        if(bc_w.eq.5 .or. bc_e.eq.5) then
           PERIODIC=.true.  ; pressureforce=.true.
        end if

        if(bc_s.eq.5 .or. bc_n.eq.5) then
           PERIODIC=.true.  ;  pressureforce_y=.true.
        end if

        if(bc_t.eq.5 .or. bc_b.eq.5)  then
           PERIODIC=.true.  ; pressureforce_z=.true.
        end if

	 dtinit=dt

        do ib=1,nbp
           dom(ib)%bc_west=bc_w
           dom(ib)%bc_east=bc_e
           dom(ib)%bc_south=bc_s
           dom(ib)%bc_north=bc_n
           dom(ib)%bc_bottom=bc_b
           dom(ib)%bc_top=bc_t
           dom(ib)%Tbc_west=Tbc_w
           dom(ib)%Tbc_east=Tbc_e
           dom(ib)%Tbc_south=Tbc_s
           dom(ib)%Tbc_north=Tbc_n
           dom(ib)%Tbc_bottom=Tbc_b
           dom(ib)%Tbc_top=Tbc_t
           dom(ib)%ngrid=ngrid_input
        end do

        if(differencing.eq.3 .and. pl_ex.ne.2) then
           pl_ex=2
           print*,'error: you select WENO but do not assign',
     &'  correct number of ghost planes,now it is corrected to 2'
        end if

        if(SGS .and. sgs_model.eq.3 .and. pl_ex.ne.2) then
           pl_ex=2
           print*,'error: you select 1-EQN model but do not assign',
     &'  correct number of ghost planes,now it is corrected to 2'
        end if

	  if (L_LSM .and. solver.eq.1) then
	   if (myrank.eq.0) then
          print*,'Error: SIP solver not presently compatible with LSM'
	   endif
	   stop
	  endif

	  if (L_LSM .and. differencing.ne.3) then
	   if (myrank.eq.0) then
          print*,'Error: WENO differencing must be used with LSM'
	   endif
	   stop
	  endif

	  if (L_LSMbase .and. L_LSM) then
	   if (myrank.eq.0) then
          print*,'Error: L_LSMbase and L_LSM cannot both be true!'
	   endif
	   stop
	  endif

	  if (L_anim_phi .and. L_LSM.eq..FALSE.) then
	   if (myrank.eq.0) then
          print*,'Error: L_anim_phi cannot be true if L_LSM is false!'
	   endif
	   stop
	  endif

	  if (L_LSM .and. bc_w.ne.1 .or. L_LSMbase .and. bc_w.ne.1) then
	   if (myrank.eq.0) then
          print*,'Warning: LSM might not work well with this inlet BC!'
	   endif
!	   stop
	  endif

        end 

!##########################################################################
        subroutine initial
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,ib,tti,ttj,ttk
        integer :: ii,jj,kk
	! id: subdomain id number
        integer :: glevel,gl,mgc_i,mgc_j,mgc_k,is,ie,js,je,ks,ke
	! glevel: an integer defined and used in multigrid used linked with multigrid step(ngrid)
	! mgc_i,mgc_j,mgc_k: integers/constants used in multigrid, calculated as a function of the total cells
	! is,ie,js,je,ks,ke,..., isp,iep,jsp,jep,ksp,kep : the first letter stands for the direction (i in x,j in y and k in z). e stands for the end of a sub-domain and s for the start 
	! and the final letter is used to know that you evaluate a variable at the centre of the cell, i.e pressure. In general you use this and the following one and also ie,je,ke, 
	!and is,js,ks to define the start and the end of a do loop in the code (basically you scan a sub-domain from the start to the end 


        double precision    :: ndx,ndy,ndz,nwxend,nwyend,nwzend

        do ib=1,nbp
           tti=dom(ib)%ttc_i; ttj=dom(ib)%ttc_j
           ttk=dom(ib)%ttc_k

           do i=1,26
              dom(ib)%tg(i)=i*10**5+dom_id(ib)
           end do

           allocate(dom(ib)%x(tti),dom(ib)%y(ttj),dom(ib)%z(ttk))
           allocate(dom(ib)%xc(tti),dom(ib)%yc(ttj),dom(ib)%zc(ttk))

           allocate(dom(ib)%u(tti,ttj,ttk))
           allocate(dom(ib)%v(tti,ttj,ttk))
           allocate(dom(ib)%w(tti,ttj,ttk))
           allocate(dom(ib)%p(tti,ttj,ttk))
           allocate(dom(ib)%pp(tti,ttj,ttk))
           allocate(dom(ib)%sup(tti,ttj,ttk))
           allocate(dom(ib)%ustar(tti,ttj,ttk))
           allocate(dom(ib)%vstar(tti,ttj,ttk))
           allocate(dom(ib)%wstar(tti,ttj,ttk))

           allocate(dom(ib)%uo(tti,ttj,ttk),dom(ib)%uoo(tti,ttj,ttk))
           allocate(dom(ib)%vo(tti,ttj,ttk),dom(ib)%voo(tti,ttj,ttk))  
           allocate(dom(ib)%wo(tti,ttj,ttk),dom(ib)%woo(tti,ttj,ttk))

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Aplied in the ghost cell method
           !Creating a copy of the current velocity field:
           ALLOCATE(dom(ib)%ucpy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%vcpy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%wcpy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%pcpy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%phicpy(tti,ttj,ttk))
           !!!ALLOCATE(dom(ib)%phimcpy(tti,ttj,ttk))
           !!!ALLOCATE(dom(ib)%rhocpy(tti,ttj,ttk))
           !!!ALLOCATE(dom(ib)%mucpy(tti,ttj,ttk))
           !!!
           ! variables for weno scheme
           ALLOCATE(dom(ib)%dudx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%dudy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%dudz(tti,ttj,ttk))
           !!!
           ALLOCATE(dom(ib)%dvdx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%dvdy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%dvdz(tti,ttj,ttk))
           !!!
           ALLOCATE(dom(ib)%dwdx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%dwdy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%dwdz(tti,ttj,ttk))
           !!!
           ! variables for marking cells
           ALLOCATE(dom(ib)%mucx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwcz(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mpc(tti,ttj,ttk))
           !!!
           ALLOCATE(dom(ib)%mucxold(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyold(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczold(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mpcold(tti,ttj,ttk))
           !!
           ALLOCATE(dom(ib)%mucxoldnn(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldnn(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldnn(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mpcoldnn(tti,ttj,ttk))
           !!
           ALLOCATE(dom(ib)%mucxoldib(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldib(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldib(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mpcoldib(tti,ttj,ttk))
           !!
           ALLOCATE(dom(ib)%mucxoldll(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldll(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldll(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldll(tti,ttj,ttk))
           !!                                      
           ALLOCATE(dom(ib)%mucxoldi(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldi(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldi(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldi(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxoldj(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldj(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldj(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldj(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxoldk(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldk(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldk(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldk(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxoldx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldx(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldx(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxoldy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldy(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldy(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxoldz(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldz(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldz(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldz(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxoldmi(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldmi(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldmi(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldmi(tti,ttj,ttk))
           !!                        m             
           ALLOCATE(dom(ib)%mucxoldmj(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldmj(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldmj(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldmj(tti,ttj,ttk))
           !!                        m             
           ALLOCATE(dom(ib)%mucxoldmk(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldmk(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldmk(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldmk(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxoldmx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldmx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldmx(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldmx(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxoldmy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldmy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldmy(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldmy(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxoldmz(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyoldmz(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczoldmz(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcoldmz(tti,ttj,ttk))
           !!                                     
           do kk=1,dom(ib)%ttc_k
            do ii=1,dom(ib)%ttc_i
            do jj=1,dom(ib)%ttc_j
               dom(ib)%mucxold(ii,jj,kk) = 0
               !dom(ib)%mvcyold(ii,jj,kk) = 0
               dom(ib)%mwczold(ii,jj,kk) = 0
               dom(ib)%mpcold(ii,jj,kk) = 0
            end do
            end do
            end do
           !!                                     
           ALLOCATE(dom(ib)%mucxnn(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcynn(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwcznn(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mpcnn(tti,ttj,ttk))
           !!
           ALLOCATE(dom(ib)%mucxib(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyib(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczib(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mpcib(tti,ttj,ttk))
           !!
           ALLOCATE(dom(ib)%mucxll(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyll(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczll(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcll(tti,ttj,ttk))
           !!                                      
           ALLOCATE(dom(ib)%mucxi(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyi(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczi(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpci(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxj(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyj(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczj(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcj(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxk(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyk(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczk(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpck(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczx(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcx(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczy(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcy(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxz(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcyz(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczz(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcz(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxmi(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcymi(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczmi(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcmi(tti,ttj,ttk))
           !!                        m             
           ALLOCATE(dom(ib)%mucxmj(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcymj(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczmj(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcmj(tti,ttj,ttk))
           !!                        m             
           ALLOCATE(dom(ib)%mucxmk(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcymk(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczmk(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcmk(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxmx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcymx(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczmx(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcmx(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxmy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcymy(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczmy(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcmy(tti,ttj,ttk))
           !!                                     
           ALLOCATE(dom(ib)%mucxmz(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mvcymz(tti,ttj,ttk))
           ALLOCATE(dom(ib)%mwczmz(tti,ttj,ttk))
           ALLOCATE( dom(ib)%mpcmz(tti,ttj,ttk))
           !!
           ! LU matrices for the bilinear least squared, q4
           ALLOCATE(dom(ib)%a3blsu(tti,ttj,ttk,6,4))
           ALLOCATE(dom(ib)%a3blsv(tti,ttj,ttk,6,4))
           ALLOCATE(dom(ib)%a3blsw(tti,ttj,ttk,6,4))
           ALLOCATE(dom(ib)%a3blsp(tti,ttj,ttk,6,4))
           ! permuation matrices for the bilinear least squared, q4
           ALLOCATE(dom(ib)%a2blsu(tti,ttj,ttk,4,4))
           ALLOCATE(dom(ib)%a2blsv(tti,ttj,ttk,4,4))
           ALLOCATE(dom(ib)%a2blsw(tti,ttj,ttk,4,4))
           ALLOCATE(dom(ib)%a2blsp(tti,ttj,ttk,4,4))
           ! pivot vectors for the bilinear least squared, q4
           ALLOCATE(dom(ib)%p2blsu(tti,ttj,ttk,4))
           ALLOCATE(dom(ib)%p2blsv(tti,ttj,ttk,4))
           ALLOCATE(dom(ib)%p2blsw(tti,ttj,ttk,4))
           ALLOCATE(dom(ib)%p2blsp(tti,ttj,ttk,4))
           !
           !
           ! LU matrices for the tseng  method, q4
           ALLOCATE(dom(ib)%a3tsengu(tti,ttj,ttk,6,6))
           ALLOCATE(dom(ib)%a3tsengv(tti,ttj,ttk,6,6))
           ALLOCATE(dom(ib)%a3tsengw(tti,ttj,ttk,6,6))
           ALLOCATE(dom(ib)%a3tsengp(tti,ttj,ttk,6,6))
           ! permuation matrices for the tseng  method, q4
           ALLOCATE(dom(ib)%a2tsengu(tti,ttj,ttk,4,6))
           ALLOCATE(dom(ib)%a2tsengv(tti,ttj,ttk,4,6))
           ALLOCATE(dom(ib)%a2tsengw(tti,ttj,ttk,4,6))
           ALLOCATE(dom(ib)%a2tsengp(tti,ttj,ttk,4,6))
           ! pivot vectors for the tseng  method, q4
           ALLOCATE(dom(ib)%p2tsengu(tti,ttj,ttk,6))
           ALLOCATE(dom(ib)%p2tsengv(tti,ttj,ttk,6))
           ALLOCATE(dom(ib)%p2tsengw(tti,ttj,ttk,6))
           ALLOCATE(dom(ib)%p2tsengp(tti,ttj,ttk,6))
           !
           !
           ! LU matrices for the quadratic least square method, q4
           ALLOCATE(dom(ib)%aqlsu(tti,ttj,ttk,6,6))
           ALLOCATE(dom(ib)%aqlsv(tti,ttj,ttk,6,6))
           ALLOCATE(dom(ib)%aqlsw(tti,ttj,ttk,6,6))
           ALLOCATE(dom(ib)%aqlsp(tti,ttj,ttk,6,6))
           !
           ! permuation matrices for the quadratic least square method, q4
           ALLOCATE(dom(ib)%a2qlsu(tti,ttj,ttk,4,6))
           ALLOCATE(dom(ib)%a2qlsv(tti,ttj,ttk,4,6))
           ALLOCATE(dom(ib)%a2qlsw(tti,ttj,ttk,4,6))
           ALLOCATE(dom(ib)%a2qlsp(tti,ttj,ttk,4,6))
           !
           ! pivot vectors for the quadratic least square method, q4
           ALLOCATE(dom(ib)%pqlsu(tti,ttj,ttk,6))
           ALLOCATE(dom(ib)%pqlsv(tti,ttj,ttk,6))
           ALLOCATE(dom(ib)%pqlsw(tti,ttj,ttk,6))
           ALLOCATE(dom(ib)%pqlsp(tti,ttj,ttk,6))
           !
           !
           ! determinatn fot the quadratic least square method, q16
           ALLOCATE(dom(ib)%da3u(tti,ttj,ttk))
           ALLOCATE(dom(ib)%da3v(tti,ttj,ttk))
           ALLOCATE(dom(ib)%da3w(tti,ttj,ttk))
           ALLOCATE(dom(ib)%da3p(tti,ttj,ttk))
           !
           ! Transpose matrices for the quadratic least square method, q16
           ALLOCATE(dom(ib)%aspu(tti,ttj,ttk,16,6))
           ALLOCATE(dom(ib)%aspv(tti,ttj,ttk,16,6))
           ALLOCATE(dom(ib)%aspw(tti,ttj,ttk,16,6))
           ALLOCATE(dom(ib)%aspp(tti,ttj,ttk,16,6))
           !
           ! Inverse matrices for the quadratic least square method, q16
           ALLOCATE(dom(ib)%ainvu(tti,ttj,ttk,6,6))
           ALLOCATE(dom(ib)%ainvv(tti,ttj,ttk,6,6))
           ALLOCATE(dom(ib)%ainvw(tti,ttj,ttk,6,6))
           ALLOCATE(dom(ib)%ainvp(tti,ttj,ttk,6,6))
           !
           !
           ! q4 HP FEM shape functions
           ALLOCATE(dom(ib)%np1(tti,ttj,ttk,12))
           ALLOCATE(dom(ib)%npu(tti,ttj,ttk,12))
           ALLOCATE(dom(ib)%npv(tti,ttj,ttk,12))
           ALLOCATE(dom(ib)%npw(tti,ttj,ttk,12))
           ALLOCATE(dom(ib)%npp(tti,ttj,ttk,12))
           !
           ! q4 bilinear shape functions
           ALLOCATE(dom(ib)%nblu(tti,ttj,ttk,4))
           ALLOCATE(dom(ib)%nblv(tti,ttj,ttk,4))
           ALLOCATE(dom(ib)%nblw(tti,ttj,ttk,4))
           ALLOCATE(dom(ib)%nblp(tti,ttj,ttk,4))
           !
           ! q9 Lagrange shape functions
           ALLOCATE(dom(ib)%nlgu(tti,ttj,ttk,9))
           ALLOCATE(dom(ib)%nlgv(tti,ttj,ttk,9))
           ALLOCATE(dom(ib)%nlgw(tti,ttj,ttk,9))
           ALLOCATE(dom(ib)%nlgp(tti,ttj,ttk,9))
           !!                      
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           allocate(dom(ib)%ap(tti,ttj,ttk),dom(ib)%su(tti,ttj,ttk))
           allocate(dom(ib)%ae(tti,ttj,ttk),dom(ib)%aw(tti,ttj,ttk))
           allocate(dom(ib)%an(tti,ttj,ttk),dom(ib)%as(tti,ttj,ttk))
           allocate(dom(ib)%at(tti,ttj,ttk),dom(ib)%ab(tti,ttj,ttk))

           allocate(dom(ib)%um(tti,ttj,ttk),dom(ib)%vm(tti,ttj,ttk))
           allocate(dom(ib)%wm(tti,ttj,ttk),dom(ib)%pm(tti,ttj,ttk))
           allocate(dom(ib)%uum(tti,ttj,ttk),dom(ib)%vvm(tti,ttj,ttk))
           allocate(dom(ib)%wwm(tti,ttj,ttk),dom(ib)%uvm(tti,ttj,ttk))
           allocate(dom(ib)%uwm(tti,ttj,ttk),dom(ib)%vwm(tti,ttj,ttk))
           allocate(dom(ib)%ppm(tti,ttj,ttk))

           allocate(dom(ib)%ntav1(tti,ttj,ttk))
           allocate(dom(ib)%ntav2(tti,ttj,ttk))
           allocate(dom(ib)%facp1(tti,ttj,ttk))
           allocate(dom(ib)%facp2(tti,ttj,ttk))
           allocate(dom(ib)%facm1(tti,ttj,ttk))
           allocate(dom(ib)%facm2(tti,ttj,ttk))

	     dom(ib)%ntav1=0
    	     dom(ib)%ntav2=0
    	     dom(ib)%facp1=0
    	     dom(ib)%facp2=0
    	     dom(ib)%facm1=1
    	     dom(ib)%facm2=1

           allocate(dom(ib)%vis(tti,ttj,ttk),dom(ib)%vism(tti,ttj,ttk))
           allocate(dom(ib)%epsm(tti,ttj,ttk))
           allocate(dom(ib)%ksgs(tti,ttj,ttk))
           allocate(dom(ib)%ksgso(tti,ttj,ttk))
           allocate(dom(ib)%eps(tti,ttj,ttk))
           allocate(dom(ib)%epso(tti,ttj,ttk))
           allocate (dom(ib)%stfcinf(6,pl,ngg))
           allocate(dom(ib)%T(tti,ttj,ttk),dom(ib)%To(tti,ttj,ttk))
           allocate(dom(ib)%Tm(tti,ttj,ttk),dom(ib)%Ttm(tti,ttj,ttk))
	     if (LSCALAR) then
            allocate(dom(ib)%S(tti,ttj,ttk),dom(ib)%Sm(tti,ttj,ttk))
           	allocate(dom(ib)%So(tti,ttj,ttk),dom(ib)%Stm(tti,ttj,ttk))
	      allocate(dom(ib)%SUtm(tti,ttj,ttk))
	      allocate(dom(ib)%SVtm(tti,ttj,ttk))
	      allocate(dom(ib)%SWtm(tti,ttj,ttk))
	   	allocate(dom(ib)%sfactor(tti,ttj,ttk))
	   	allocate(dom(ib)%ibfactor(tti,ttj,ttk))
	     endif

           if (L_LSM .or. L_LSMbase)  then
            allocate(dom(ib)%dens(tti,ttj,ttk),
     &      dom(ib)%mu(tti,ttj,ttk),dom(ib)%ijkp_lsm(0:dom(ib)%ngrid))
	      allocate(dom(ib)%phi(tti,ttj,ttk),dom(ib)%phim(tti,ttj,ttk))
	    endif
        
           if (L_LSM) then! .or. L_LSMbase) then
             dom(ib)%ijkp_lsm = 0
             dom(ib)%ijkp_lsm(1)=(dom(ib)%ttc_i-2*pl)*
     & (dom(ib)%ttc_j-2*pl)*(dom(ib)%ttc_k-2*pl) 
             do glevel=2,dom(ib)%ngrid
               mgc_i=(dom(ib)%iep-dom(ib)%isp+1)/2**(glevel-1)+2
               mgc_j=(dom(ib)%jep-dom(ib)%jsp+1)/2**(glevel-1)+2
               mgc_k=(dom(ib)%kep-dom(ib)%ksp+1)/2**(glevel-1)+2
               dom(ib)%ijkp_lsm(glevel)=dom(ib)%ijkp_lsm(glevel-1)+
     & (mgc_i-2)*(mgc_j-2)*(mgc_k-2)
             end do
             dom(ib)%tot=dom(ib)%ijkp_lsm(dom(ib)%ngrid)
             allocate (dom(ib)%dens_mg(dom(ib)%tot))
           end if

           if (differencing.eq.3) allocate(dom(ib)%d1(tti,ttj,ttk),
     & dom(ib)%dphi_dxplus(tti,ttj,ttk),
     & dom(ib)%dphi_dxminus(tti,ttj,ttk),
     & dom(ib)%dphi_dyplus(tti,ttj,ttk),
     & dom(ib)%dphi_dyminus(tti,ttj,ttk),
     & dom(ib)%dphi_dzplus(tti,ttj,ttk),
     & dom(ib)%dphi_dzminus(tti,ttj,ttk))

              allocate (dom(ib)%tauwe(
     &dom(ib)%jsp-1:dom(ib)%jep+1,dom(ib)%ksp-1:dom(ib)%kep+1))
              allocate (dom(ib)%tauww(
     &dom(ib)%jsp-1:dom(ib)%jep+1,dom(ib)%ksp-1:dom(ib)%kep+1))
              allocate (dom(ib)%tauwn(
     &dom(ib)%isp-1:dom(ib)%iep+1,dom(ib)%ksp-1:dom(ib)%kep+1))
              allocate (dom(ib)%tauws(
     &dom(ib)%isp-1:dom(ib)%iep+1,dom(ib)%ksp-1:dom(ib)%kep+1))
              allocate (dom(ib)%tauwt(
     &dom(ib)%isp-1:dom(ib)%iep+1,dom(ib)%jsp-1:dom(ib)%jep+1))
              allocate (dom(ib)%tauwb(
     &dom(ib)%isp-1:dom(ib)%iep+1,dom(ib)%jsp-1:dom(ib)%jep+1))
              allocate (dom(ib)%tauwe2(
     &dom(ib)%jsp-1:dom(ib)%jep+1,dom(ib)%ksp-1:dom(ib)%kep+1))
              allocate (dom(ib)%tauww2(
     &dom(ib)%jsp-1:dom(ib)%jep+1,dom(ib)%ksp-1:dom(ib)%kep+1))
              allocate (dom(ib)%tauwn2(
     &dom(ib)%isp-1:dom(ib)%iep+1,dom(ib)%ksp-1:dom(ib)%kep+1))
              allocate (dom(ib)%tauws2(
     &dom(ib)%isp-1:dom(ib)%iep+1,dom(ib)%ksp-1:dom(ib)%kep+1))
              allocate (dom(ib)%tauwt2(
     &dom(ib)%isp-1:dom(ib)%iep+1,dom(ib)%jsp-1:dom(ib)%jep+1))
              allocate (dom(ib)%tauwb2(
     &dom(ib)%isp-1:dom(ib)%iep+1,dom(ib)%jsp-1:dom(ib)%jep+1))

           if(solver.eq.2) then
              allocate (dom(ib)%faz(ngrd_gl))
           end if

           allocate (dom(ib)%sendb_m1(ngg)) 
           allocate (dom(ib)%sendb_p1(ngg))
           allocate (dom(ib)%recvb_m1(ngg))
           allocate (dom(ib)%recvb_p1(ngg))
           allocate (dom(ib)%sendb_m2(ngg)) 
           allocate (dom(ib)%sendb_p2(ngg))
           allocate (dom(ib)%recvb_m2(ngg))
           allocate (dom(ib)%recvb_p2(ngg))
           allocate (dom(ib)%sendb_m3(ngg)) 
           allocate (dom(ib)%sendb_p3(ngg))
           allocate (dom(ib)%recvb_m3(ngg))
           allocate (dom(ib)%recvb_p3(ngg))

           allocate (dom(ib)%sc1m(ngc),dom(ib)%sc1p(ngc))
           allocate (dom(ib)%rc1m(ngc),dom(ib)%rc1p(ngc))
           allocate (dom(ib)%sc2m(ngc),dom(ib)%sc2p(ngc))
           allocate (dom(ib)%rc2m(ngc),dom(ib)%rc2p(ngc))
           allocate (dom(ib)%sc3m(ngc),dom(ib)%sc3p(ngc))
           allocate (dom(ib)%rc3m(ngc),dom(ib)%rc3p(ngc))
           allocate (dom(ib)%sc4m(ngc),dom(ib)%sc4p(ngc))
           allocate (dom(ib)%rc4m(ngc),dom(ib)%rc4p(ngc))

           allocate (dom(ib)%se1m(nge),dom(ib)%se1p(nge))
           allocate (dom(ib)%re1m(nge),dom(ib)%re1p(nge))
           allocate (dom(ib)%se2m(nge),dom(ib)%se2p(nge))
           allocate (dom(ib)%re2m(nge),dom(ib)%re2p(nge))
           allocate (dom(ib)%se3m(nge),dom(ib)%se3p(nge))
           allocate (dom(ib)%re3m(nge),dom(ib)%re3p(nge))
           allocate (dom(ib)%se4m(nge),dom(ib)%se4p(nge))
           allocate (dom(ib)%re4m(nge),dom(ib)%re4p(nge))
           allocate (dom(ib)%se5m(nge),dom(ib)%se5p(nge))
           allocate (dom(ib)%re5m(nge),dom(ib)%re5p(nge))
           allocate (dom(ib)%se6m(nge),dom(ib)%se6p(nge))
           allocate (dom(ib)%re6m(nge),dom(ib)%re6p(nge))

           dom(ib)%x(1)=dom(ib)%xsl +(-pl+1)*dom(ib)%dx
           dom(ib)%y(1)=dom(ib)%ysl +(-pl+1)*dom(ib)%dy
           dom(ib)%z(1)=dom(ib)%zsl +(-pl+1)*dom(ib)%dz
           dom(ib)%xc(1)=dom(ib)%x(1)-0.5*dom(ib)%dx
           dom(ib)%yc(1)=dom(ib)%y(1)-0.5*dom(ib)%dy
           dom(ib)%zc(1)=dom(ib)%z(1)-0.5*dom(ib)%dz

           do i=2,dom(ib)%ttc_i
              dom(ib)%x(i)=dom(ib)%x(i-1)+dom(ib)%dx
           end do
           do i=2,dom(ib)%ttc_j
              dom(ib)%y(i)=dom(ib)%y(i-1)+dom(ib)%dy
           end do
           do i=2,dom(ib)%ttc_k
              dom(ib)%z(i)=dom(ib)%z(i-1)+dom(ib)%dz
           end do

           do i=1,dom(ib)%ttc_i
              dom(ib)%xc(i)=dom(ib)%x(i)-0.5*dom(ib)%dx
           end do
           do i=1,dom(ib)%ttc_j
              dom(ib)%yc(i)=dom(ib)%y(i)-0.5*dom(ib)%dy
           end do
           do i=1,dom(ib)%ttc_k
              dom(ib)%zc(i)=dom(ib)%z(i)-0.5*dom(ib)%dz
           end do

           if (dom(ib)%inext.ge.0) then
              if(abs(dom(ib)%x(dom(ib)%iep)-dom(ib)%xel).gt.1e-5) then
                 print*,'mycpu#:',myrank,' error-11'
                 stop
              end if
           end if

           if (dom(ib)%jnext.ge.0) then
              if(abs(dom(ib)%y(dom(ib)%jep)-dom(ib)%yel).gt.1e-5) then
                 print*,'mycpu#:',myrank,' error-12'
                 stop
              end if
           end if

           if (dom(ib)%knext.ge.0) then
              if(abs(dom(ib)%z(dom(ib)%kep)-dom(ib)%zel).gt.1e-5) then
                 print*,'mycpu#:',myrank,' error-13'
                 stop
              end if
           end if

           if(solver.eq.2 .and. ngrd_gl.ge.2) then

              is=dom(ib)%isp; ie=dom(ib)%iep
              js=dom(ib)%jsp; je=dom(ib)%jep
              ks=dom(ib)%ksp; ke=dom(ib)%kep

              do glevel=2,ngrd_gl
                 if(glevel.gt.dom(ib)%ngrid) then
                    gl=dom(ib)%ngrid
                 else
                    gl=glevel
                 end if
                 mgc_i=(ie-is+1)/2**(gl-1)+2
                 mgc_j=(je-js+1)/2**(gl-1)+2
                 mgc_k=(ke-ks+1)/2**(gl-1)+2

                 ndx=dom(ib)%dx*2**(gl-1)
                 ndy=dom(ib)%dy*2**(gl-1)
                 ndz=dom(ib)%dz*2**(gl-1)

                 nwxend=dom(ib)%x(dom(ib)%isp-1)+ndx*(mgc_i-2)
                 nwyend=dom(ib)%y(dom(ib)%jsp-1)+ndy*(mgc_j-2)
                 nwzend=dom(ib)%z(dom(ib)%ksp-1)+ndz*(mgc_k-2)

                 if((abs(dom(ib)%x(dom(ib)%iep)-nwxend).gt.1e-8)
     &  .or.(abs(dom(ib)%y(dom(ib)%jep)-nwyend).gt.1e-8)
     &  .or.(abs(dom(ib)%z(dom(ib)%kep)-nwzend).gt.1e-8)) then
                    print*,'==ERROR==> in multigrid: max ngrid value'
                    stop
                 end if
              end do 
           end if

        end do
        

        end subroutine initial
!##########################################################################
        subroutine iniflux
!##########################################################################
        use vars
        use multidata
        use mpi
        implicit none
        integer i,j,k,ib,ispr,iepr,jspr,jepr,kspr,kepr
        double precision buffer_flomas

        MPI_FLT = MPI_DOUBLE_PRECISION

        flomas=0.0

        do ib=1,nbp
           ispr=pl+1; iepr=dom(ib)%ttc_i-pl
           jspr=pl+1; jepr=dom(ib)%ttc_j-pl
           kspr=pl+1; kepr=dom(ib)%ttc_k-pl

           if(dom(ib)%iprev.lt.0) then
		if (dom(ib)%bc_west.lt.61 .and. dom(ib)%bc_west.ne.4) then
              do j=jspr,jepr 
                 do k=kspr,kepr
	             if (L_LSM) then
			   if (dom(ib)%phi(dom(ib)%isu-1,j,k) .ge. 0.0) then
                       flomas=flomas+dom(ib)%u(dom(ib)%isu-1,j,k)*
     &  dom(ib)%dy*dom(ib)%dz
			   end if
			 else
                     flomas=flomas+dom(ib)%u(dom(ib)%isu-1,j,k)*
     &  dom(ib)%dy*dom(ib)%dz
			 end if
                 end do
              end do
		endif
           end if

           if(dom(ib)%jprev.lt.0) then
		if (dom(ib)%bc_south.lt.61 .and. dom(ib)%bc_south.ne.4) then
              do k=kspr,kepr
                 do i=ispr,iepr
	             if (L_LSM) then
			   if (dom(ib)%phi(i,dom(ib)%jsv-1,k) .ge. 0.0) then
                       flomas=flomas+dom(ib)%v(i,dom(ib)%jsv-1,k)*
     &  dom(ib)%dx*dom(ib)%dz
			   end if
			 else
                     flomas=flomas+dom(ib)%v(i,dom(ib)%jsv-1,k)*
     &  dom(ib)%dx*dom(ib)%dz
			 end if
                 end do
              end do
		endif
           end if

           if(dom(ib)%kprev.lt.0) then
	     if (dom(ib)%bc_bottom.lt.61 .and. dom(ib)%bc_bottom.ne.4) then
              do j=jspr,jepr
                 do i=ispr,iepr
	             if (L_LSM) then
			   if (dom(ib)%phi(i,j,dom(ib)%ksw-1) .ge. 0.0) then
                       flomas=flomas+dom(ib)%w(i,j,dom(ib)%ksw-1)*
     &  dom(ib)%dx*dom(ib)%dy
			   end if
			 else
                     flomas=flomas+dom(ib)%w(i,j,dom(ib)%ksw-1)*
     &  dom(ib)%dx*dom(ib)%dy
			 end if
                 end do
              end do
		endif
           end if
        end do

        buffer_flomas = flomas
        call MPI_ALLREDUCE(buffer_flomas,flomas,1,MPI_FLT,MPI_SUM,
     & MPI_COMM_WORLD,ierr)
         
        return
        end 
!##########################################################################
        subroutine correctoutflux
!##########################################################################
        use vars
        use multidata
        use mpi
        implicit none
        integer i,j,k,ib,ispr,iepr,jspr,jepr,kspr,kepr
        double precision fmout,fct,buffer_fmout

        MPI_FLT = MPI_DOUBLE_PRECISION

        fmout=0.0

        do ib=1,nbp
           ispr=pl+1; iepr=dom(ib)%ttc_i-pl
           jspr=pl+1; jepr=dom(ib)%ttc_j-pl
           kspr=pl+1; kepr=dom(ib)%ttc_k-pl

           if(dom(ib)%inext.lt.0) then
		if (dom(ib)%bc_east.lt.61 .and. dom(ib)%bc_east.ne.4) then    
              do j=jspr,jepr 
                 do k=kspr,kepr
	             if (L_LSM) then
			   if (dom(ib)%phi(dom(ib)%ieu+1,j,k) .ge. 0.0) then  
                       fmout=fmout+dom(ib)%u(dom(ib)%ieu+1,j,k)*
     &  dom(ib)%dy*dom(ib)%dz
			   end if
			 else
                     fmout=fmout+dom(ib)%u(dom(ib)%ieu+1,j,k)*
     &  dom(ib)%dy*dom(ib)%dz
			 end if
                 end do
              end do    
		endif
           end if

           if(dom(ib)%jnext.lt.0) then
		if (dom(ib)%bc_north.lt.61 .and. dom(ib)%bc_north.ne.4) then    
              do k=kspr,kepr
                 do i=ispr,iepr
	             if (L_LSM) then
			   if (dom(ib)%phi(i,dom(ib)%jev+1,k) .ge. 0.0) then
                       fmout=fmout+dom(ib)%v(i,dom(ib)%jev+1,k)*
     &  dom(ib)%dx*dom(ib)%dz
			   end if
			 else
                     fmout=fmout+dom(ib)%v(i,dom(ib)%jev+1,k)*
     &  dom(ib)%dx*dom(ib)%dz
			 end if
                 end do
              end do   
		endif 
           end if

           if(dom(ib)%knext.lt.0) then
		if (dom(ib)%bc_top.lt.61 .and. dom(ib)%bc_top.ne.4) then    
              do j=jspr,jepr
                 do i=ispr,iepr
	             if (L_LSM) then
			   if (dom(ib)%phi(i,j,dom(ib)%kew+1) .ge. 0.0) then
                       fmout=fmout+dom(ib)%w(i,j,dom(ib)%kew+1)*
     &  dom(ib)%dx*dom(ib)%dy
			   end if
			 else
                     fmout=fmout+dom(ib)%w(i,j,dom(ib)%kew+1)*
     &  dom(ib)%dx*dom(ib)%dy
			 end if
                 end do
              end do   
		endif 
           end if
        end do

        buffer_fmout = fmout
        call MPI_ALLREDUCE(buffer_fmout,fmout,1,MPI_FLT,MPI_SUM,
     &MPI_COMM_WORLD,ierr)
!Ratio of inflow and outflow
        fct=flomas/(fmout+1.E-30)
!Mass deficit
        Mdef=flomas-fmout
!Update outflow regarding the mass deficit at the inlet
        do ib=1,nbp
           if(dom(ib)%inext.lt.0) then
		if (dom(ib)%bc_east.lt.61 .and. dom(ib)%bc_east.ne.4) then    
              do j=dom(ib)%jsu,dom(ib)%jeu 
                 do k=dom(ib)%ksu,dom(ib)%keu  
	             if (L_LSM) then							!PABLO 04/18
			   if (dom(ib)%phi(dom(ib)%ieu+1,j,k) .ge. 0.0) then  
           dom(ib)%u(dom(ib)%ieu+1,j,k)=dom(ib)%u(dom(ib)%ieu+1,j,k)*fct
			   end if
			 else
           dom(ib)%u(dom(ib)%ieu+1,j,k)=dom(ib)%u(dom(ib)%ieu+1,j,k)*fct
			 end if
 !          dom(ib)%u(dom(ib)%ieu+1,j,k)=dom(ib)%u(dom(ib)%ieu+1,j,k)*fct
                 end do
              end do

              do j=dom(ib)%jsv,dom(ib)%jev 
                 do k=dom(ib)%ksv,dom(ib)%kev  
	             if (L_LSM) then							!PABLO 04/18
			   if (dom(ib)%phi(dom(ib)%iev+1,j,k) .ge. 0.0) then  
           dom(ib)%v(dom(ib)%iev+1,j,k)=dom(ib)%v(dom(ib)%iev+1,j,k)*fct
			   end if
			 else
           dom(ib)%v(dom(ib)%iev+1,j,k)=dom(ib)%v(dom(ib)%iev+1,j,k)*fct
			 end if
!           dom(ib)%v(dom(ib)%iev+1,j,k)=dom(ib)%v(dom(ib)%iev+1,j,k)*fct
                 end do
              end do

              do j=dom(ib)%jsw,dom(ib)%jew 
                 do k=dom(ib)%ksw,dom(ib)%kew 
	             if (L_LSM) then							!PABLO 04/18
			   if (dom(ib)%phi(dom(ib)%iew+1,j,k) .ge. 0.0) then  
           dom(ib)%w(dom(ib)%iew+1,j,k)=dom(ib)%w(dom(ib)%iew+1,j,k)*fct
			   end if
			 else
           dom(ib)%w(dom(ib)%iew+1,j,k)=dom(ib)%w(dom(ib)%iew+1,j,k)*fct
			 end if 
!           dom(ib)%w(dom(ib)%iew+1,j,k)=dom(ib)%w(dom(ib)%iew+1,j,k)*fct
                 end do
              end do
		endif
           end if

           if(dom(ib)%jnext.lt.0) then
		if (dom(ib)%bc_north.lt.61 .and. dom(ib)%bc_north.ne.4) then    
              do i=dom(ib)%isu,dom(ib)%ieu 
                 do k=dom(ib)%ksu,dom(ib)%keu  
	             if (L_LSM) then							!PABLO 04/18
			   if (dom(ib)%phi(i,dom(ib)%jeu+1,k) .ge. 0.0) then  
           dom(ib)%u(i,dom(ib)%jeu+1,k)=dom(ib)%u(i,dom(ib)%jeu+1,k)*fct
			   end if
			 else
           dom(ib)%u(i,dom(ib)%jeu+1,k)=dom(ib)%u(i,dom(ib)%jeu+1,k)*fct
			 end if
!           dom(ib)%u(i,dom(ib)%jeu+1,k)=dom(ib)%u(i,dom(ib)%jeu+1,k)*fct
                 end do
              end do

              do i=dom(ib)%isv,dom(ib)%iev 
                 do k=dom(ib)%ksv,dom(ib)%kev  
	             if (L_LSM) then							!PABLO 04/18
			   if (dom(ib)%phi(i,dom(ib)%jev+1,k) .ge. 0.0) then  
           dom(ib)%v(i,dom(ib)%jev+1,k)=dom(ib)%v(i,dom(ib)%jev+1,k)*fct
			   end if
			 else
           dom(ib)%v(i,dom(ib)%jev+1,k)=dom(ib)%v(i,dom(ib)%jev+1,k)*fct
			 end if
!           dom(ib)%v(i,dom(ib)%jev+1,k)=dom(ib)%v(i,dom(ib)%jev+1,k)*fct
                 end do
              end do

              do i=dom(ib)%isw,dom(ib)%iew 
                 do k=dom(ib)%ksw,dom(ib)%kew
	             if (L_LSM) then							!PABLO 04/18
			   if (dom(ib)%phi(i,dom(ib)%jew+1,k) .ge. 0.0) then  
           dom(ib)%w(i,dom(ib)%jew+1,k)=dom(ib)%w(i,dom(ib)%jew+1,k)*fct
			   end if
			 else
           dom(ib)%w(i,dom(ib)%jew+1,k)=dom(ib)%w(i,dom(ib)%jew+1,k)*fct
			 end if  
!           dom(ib)%w(i,dom(ib)%jew+1,k)=dom(ib)%w(i,dom(ib)%jew+1,k)*fct
                 end do
              end do
		endif
           end if

           if(dom(ib)%knext.lt.0) then
		if (dom(ib)%bc_top.lt.61 .and. dom(ib)%bc_top.ne.4) then    
              do i=dom(ib)%isu,dom(ib)%ieu 
                 do j=dom(ib)%jsu,dom(ib)%jeu
	             if (L_LSM) then							!PABLO 04/18
			   if (dom(ib)%phi(i,j,dom(ib)%keu+1) .ge. 0.0) then  
           dom(ib)%u(i,j,dom(ib)%keu+1)=dom(ib)%u(i,j,dom(ib)%keu+1)*fct
			   end if
			 else
           dom(ib)%u(i,j,dom(ib)%keu+1)=dom(ib)%u(i,j,dom(ib)%keu+1)*fct
			 end if   
!           dom(ib)%u(i,j,dom(ib)%keu+1)=dom(ib)%u(i,j,dom(ib)%keu+1)*fct
                 end do
              end do

              do i=dom(ib)%isv,dom(ib)%iev 
                 do j=dom(ib)%jsv,dom(ib)%jev 
	             if (L_LSM) then							!PABLO 04/18
			   if (dom(ib)%phi(i,j,dom(ib)%kev+1) .ge. 0.0) then  
           dom(ib)%v(i,j,dom(ib)%kev+1)=dom(ib)%v(i,j,dom(ib)%kev+1)*fct
			   end if
			 else
           dom(ib)%v(i,j,dom(ib)%kev+1)=dom(ib)%v(i,j,dom(ib)%kev+1)*fct
			 end if  
!           dom(ib)%v(i,j,dom(ib)%kev+1)=dom(ib)%v(i,j,dom(ib)%kev+1)*fct
                 end do
              end do

              do i=dom(ib)%isw,dom(ib)%iew 
                 do j=dom(ib)%jsw,dom(ib)%jew 
	             if (L_LSM) then							!PABLO 04/18
			   if (dom(ib)%phi(i,j,dom(ib)%kew+1) .ge. 0.0) then  
           dom(ib)%w(i,j,dom(ib)%kew+1)=dom(ib)%w(i,j,dom(ib)%kew+1)*fct
			   end if
			 else
           dom(ib)%w(i,j,dom(ib)%kew+1)=dom(ib)%w(i,j,dom(ib)%kew+1)*fct
			 end if  
!           dom(ib)%w(i,j,dom(ib)%kew+1)=dom(ib)%w(i,j,dom(ib)%kew+1)*fct
                 end do
              end do
		endif
           end if
        end do

        return
        end 
!##########################################################################
        subroutine initflowfield
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib,tti,ttj,ttk,pll
        integer :: sn,sn2,jtime,idfile,inind,jnind,knind
        double precision dum,ubw,ube,ubs,ubn,ubt,ubb,vb,wb,lz,dummy
	  double precision ::ufric
        double precision, dimension(23) :: dm
        character*8   :: chb1,numpt,x,y,z
        character*25  :: gf,unst
        character*100 :: dummyline

        do ib=1,nbp

           tti=dom(ib)%ttc_i
           ttj=dom(ib)%ttc_j
           ttk=dom(ib)%ttc_k

           if (LRESTART) then

              qzero=ubulk 
              open (unit=700, file='final_ctime.dat')
               read (700,'(i8,3F18.12)') ntime,ctime,forcn,qstpn
               gh_startTime = ntime
              close (700)

	     if(pressureforce_y .eq. .true.) then				!PABLO 12/18
              open (unit=700, file='forcn_y.dat')
	       read (700,*) 
		   do i=1,ntime
                read (700,'(4F18.12)',END=200) dum,forcn_y,qstpn_y,dum
	 	   enddo
200	continue
              close (700)
		    endif

	     if(pressureforce_z .eq. .true.) then
              open (unit=700, file='forcn_z.dat')
		    do i=1,ntime
                 read (700,'(4F18.12)',END=201) dum,forcn_z,qstpn_z,dum
		    enddo
201	continue
              close (700)
		    endif

              write(chb1,'(i4)') dom_id(ib)
              sn=len(trim(adjustl(chb1)))
              chb1=repeat('0',(4-sn))//trim(adjustl(chb1))

!===============================================================
              gf='tecbin'//trim(adjustl(chb1))//'.bin'
              open (unit=700, file=gf, form='unformatted',status='old')

              read (700) tti,ttj,ttk
              read (700) pll
              read (700) inind,jnind,knind

              if(pll.ne.pl .and. myrank.eq.0) then
        print*,'&*&* different number of overlapping layers!!',pl,pll
        write(numfile,*) '&*&* different number of overlapping layers!!'
        stop
              end if

              do k=1,ttk
                 do j=1,ttj
                    do i=1,tti
                       read (700) dm(1),dm(2),dm(3),dm(4),dm(5),dm(6),
     & dm(7),dm(8),dm(9),dm(10),dm(11),dm(12),dm(13),dm(14),dm(15),
     & dm(16),dm(17),dm(18),dm(19),dm(20),dm(21),dm(22),dm(23)
                        dom(ib)%p  (i,j,k)=dm(4)
                        dom(ib)%pm (i,j,k)=dm(5)
                        dom(ib)%ppm(i,j,k)=dm(6)
                        dom(ib)%vis(i,j,k)=dm(7)
                        dom(ib)%vism(i,j,k)=dm(8)			!!
                		dom(ib)%u  (i,j,k)=dm(9)
                		dom(ib)%um (i,j,k)=dm(10)
                		dom(ib)%uum(i,j,k)=dm(11)
                		dom(ib)%v  (i,j,k)=dm(12)
                		dom(ib)%vm (i,j,k)=dm(13)
                		dom(ib)%vvm(i,j,k)=dm(14)
                		dom(ib)%w  (i,j,k)=dm(15)
                		dom(ib)%wm (i,j,k)=dm(16)
                		dom(ib)%wwm(i,j,k)=dm(17)
                		dom(ib)%uvm(i,j,k)=dm(18)
                		dom(ib)%uwm(i,j,k)=dm(19)
                		dom(ib)%vwm(i,j,k)=dm(20)
                        dom(ib)%ksgs(i,j,k)=dm(21)
              	      dom(ib)%eps(i,j,k) =dm(22)
              	      dom(ib)%epsm(i,j,k)=dm(23)
                    end do
                 end do
              end do
              close (700)

		IF (LENERGY.EQ..TRUE. )then
              gf='tecplot_T'//trim(adjustl(chb1))//'.plt'
              open (unit=700, file=gf,status='old')
              read (700) 
              read (700) 
              read (700) 

              do k=dom(ib)%ksp-1,dom(ib)%kep+1
                 do j=dom(ib)%jsp-1,dom(ib)%jep+1
                    do i=dom(ib)%isp-1,dom(ib)%iep+1
              read (700) dummy,dummy,dummy,
     &  dom(ib)%T(i,j,k),dom(ib)%Tm(i,j,k),dom(ib)%Ttm(i,j,k)
                   end do
                 end do
              end do
              close (700)
		ENDIF
		IF (LSCALAR.EQ..TRUE. )then
              gf='tecplot_S'//trim(adjustl(chb1))//'.plt'
              open (unit=700, file=gf,status='old')
              read (700) 
              read (700) 
              read (700) 
              do k=dom(ib)%ksp-1,dom(ib)%kep+1
                 do j=dom(ib)%jsp-1,dom(ib)%jep+1
                    do i=dom(ib)%isp-1,dom(ib)%iep+1
              read (700) dummy,dummy,dummy,
     &  dom(ib)%S(i,j,k),dom(ib)%Sm(i,j,k),dom(ib)%Stm(i,j,k),
     &  dom(ib)%SUtm(i,j,k),dom(ib)%SVtm(i,j,k),dom(ib)%SWtm(i,j,k)
                   end do
                 end do
              end do
              close (700)
		ENDIF

!===============================================================

              if (reinitmean) then
                 ctime=0.d0         ; ntime=0
                 dom(ib)%um   = 0.d0; dom(ib)%vm   = 0.d0
                 dom(ib)%wm   = 0.d0; dom(ib)%pm   = 0.d0
                 dom(ib)%uum  = 0.d0; dom(ib)%vvm  = 0.d0
                 dom(ib)%wwm  = 0.d0; dom(ib)%uvm  = 0.d0
                 dom(ib)%uwm  = 0.d0; dom(ib)%vwm  = 0.d0
                 dom(ib)%ppm  = 0.d0
                 dom(ib)%vism = 0.d0; dom(ib)%epsm = 0.d0
	     if (LENERGY) then
                 dom(ib)%Tm   = 0.d0; dom(ib)%Ttm  = 0.d0
		ENDIF
	     if (LSCALAR) then
 		     dom(ib)%Sm   = 0.d0; dom(ib)%Stm=0.d0	
                 dom(ib)%SUtm = 0.0 ; dom(ib)%SVtm  = 0.0
                 dom(ib)%SWtm = 0.0 !; dom(ib)%S   = 0.0	   
		ENDIF
 
		     if (L_LSM) dom(ib)%phim  = 0.d0
              end if
           else

              qzero=ubulk 								
              qstpn=qzero ; qstpn_y=0.d0 ;  qstpn_z=0.d0 			!Pablo 02/2018
	        forcn=2.0/(Re*qzero) ; forcn_y = 0.d0 ; forcn_z=0.d0
              ctime=0.D0
              ntime=0

	        if (L_LSMbase .or. L_LSM) then
                do k=1,ttk
                  do j=1,ttj
                    do i=1,tti
            if(abs(grz).gt.0.0001d00)then
			    if (dom(ib)%zc(k).le.length) then
                        dom(ib)%u(i,j,k)  =Ubulk  
                        dom(ib)%uo(i,j,k) =Ubulk 
			      dom(ib)%uoo(i,j,k)=Ubulk
			    else
				dom(ib)%u(i,j,k)=0.D0
                        dom(ib)%uo(i,j,k)=0.D0 
			      dom(ib)%uoo(i,j,k)=0.D0
			    end if
            elseif(abs(gry).gt.0.0001d00)then
			    if (dom(ib)%yc(j).le.length) then
                        dom(ib)%u(i,j,k)  =Ubulk  
                        dom(ib)%uo(i,j,k) =Ubulk 
			      dom(ib)%uoo(i,j,k)=Ubulk
			    else
				dom(ib)%u(i,j,k)=0.D0
                        dom(ib)%uo(i,j,k)=0.D0 
			      dom(ib)%uoo(i,j,k)=0.D0
			    end if
            endif ! gry, grz values
			  end do
		      end do
		    end do
		  else
		    dom(ib)%u=0.d0
                dom(ib)%uo=0.d0 
		    dom(ib)%uoo=0.d0	
		  end if

	   	  lz=zen-zst

              dom(ib)%v=0.d0		; dom(ib)%w=0.d0
              dom(ib)%p=0.d0		; dom(ib)%pm=0.d0
              dom(ib)%vo=0.d0		; dom(ib)%voo=0.d0
              dom(ib)%wo=0.d0		; dom(ib)%woo=0.d0
		  dom(ib)%wm   = 0.d0
              dom(ib)%um   = 0.d0	; dom(ib)%vm   = 0.d0
              dom(ib)%uum  = 0.d0	; dom(ib)%vvm  = 0.d0
              dom(ib)%wwm  = 0.d0	; dom(ib)%uvm  = 0.d0
              dom(ib)%uwm  = 0.d0	; dom(ib)%vwm  = 0.d0
              dom(ib)%ppm  = 0.d0   ; dom(ib)%vism = 0.d0
              dom(ib)%ksgs = 0.d0
              dom(ib)%eps  = 0.d0	; dom(ib)%epsm  =0.d0
	     if (LENERGY) then
              dom(ib)%T=0.d0		; dom(ib)%To=0.d0
              dom(ib)%Tm=0.d0		; dom(ib)%Ttm=0.d0
		ENDIF
	     if (LSCALAR) then
              dom(ib)%S=0.d0		; dom(ib)%So=0.d0
              dom(ib)%Sm   =0.d0	; dom(ib)%Stm=0.d0
              dom(ib)%SUtm =0.d0    ; dom(ib)%SVtm=0.d0 
              dom(ib)%SWtm=0.d0
		ENDIF
	    if (sgs_model.eq.4) then
              dom(ib)%ksgs = (3.d0/2.d0)*(ubulk*0.1)**2.0		  
              dom(ib)%eps  = 0.09**0.75*dom(ib)%ksgs**1.5/(0.07*lz) 	
	    endif

              if (trim(keyword).eq.'channel') then
                if (.not.L_LSM .and. .not.L_LSMbase) then
		       dom(ib)%u=ubulk
		    end if
                   ubw=ubulk; ube=ubulk; ubs=ubulk				
		       ubn=ubulk; ubt=ubulk; ubb=ubulk
                   vb=0.0; wb=0.0
              else if (trim(keyword).eq.'cavity') then
                 dom(ib)%u=0.0
                 ubw=0.0; ube=0.0; ubs=0.0; ubn=0.0; ubt=1.0; ubb=0.0
                 vb=0.0; wb=0.0
		  else if (trim(keyword).eq.'column') then
                 dom(ib)%u=0.0
                 ubw=0.0; ube=0.0; ubs=0.0; ubn=0.0; ubt=0.0; ubb=0.0
                 vb=0.0; wb=0.0
         else if (trim(keyword).eq.'wave') then
                 dom(ib)%u=0.0
                 ubw=0.0; ube=0.0; ubs=0.0 
                    ubn=0.0; ubt=0.0; ubb=0.0
                 vb=0.0; wb=0.0
		  else if (trim(keyword).eq.'waves') then
                 dom(ib)%u=0.0
                 ubw=0.0; ube=0.0; ubs=0.0; ubn=0.0; ubt=0.0; ubb=0.0
                 vb=0.0; wb=0.0
              else
                 write (6,*) ' wrong keyword '
              end if

!..............U=> West and East ...............
              if (dom(ib)%iprev.lt.0) then
                do k=1,ttk
                  do j=1,ttj
	              if (L_LSMbase .or. (L_LSM.and..not.lrestart)) then
                 if(abs(grz).gt.0.0001d00)then
		           if (dom(ib)%zc(k).gt.length) then   
                        dom(ib)%u(dom(ib)%isu-1,j,k) = 0.0 
		           endif
                 elseif(abs(gry).gt.0.0001d00)then
		           if (dom(ib)%yc(j).gt.length) then   
                        dom(ib)%u(dom(ib)%isu-1,j,k) = 0.0 
		           endif
                 endif ! gry, grz values
	              else
                      dom(ib)%u(dom(ib)%isu-1,j,k) = ubw
 	              end if
                  end do
                end do 
              end if
              if (dom(ib)%inext.lt.0) then
                do k=1,ttk
                  do j=1,ttj
	              if (L_LSMbase .or. (L_LSM.and..not.lrestart)) then
                 if(abs(grz).gt.0.0001d00)then
		           if (dom(ib)%zc(k).gt.length) then   
                        dom(ib)%u(dom(ib)%ieu+1,j,k) = 0.0 
		           endif
                 elseif(abs(gry).gt.0.0001d00)then
		           if (dom(ib)%yc(j).gt.length) then   
                        dom(ib)%u(dom(ib)%ieu+1,j,k) = 0.0 
		           endif
                 endif ! gry, grz values
	              else
                      dom(ib)%u(dom(ib)%ieu+1,j,k) = ube
	              end if 
                  end do
                end do 
              end if
!.............U=> South and North .................
              if (dom(ib)%jprev.lt.0) then
                do k=1,ttk
                  do i=1,tti
	              if (L_LSMbase .or. (L_LSM.and..not.lrestart)) then
                 if(abs(grz).gt.0.0001d00)then
		           if (dom(ib)%zc(k).gt.length) then   
                        dom(ib)%u(i,dom(ib)%jsu-1,k) = 0.0 
		           endif
                 elseif(abs(gry).gt.0.0001d00)then
		           if (dom(ib)%yc(dom(ib)%jsu-1).gt.length) then   
                        dom(ib)%u(i,dom(ib)%jsu-1,k) = 0.0 
		           endif
                 endif ! gry, grz values
	              else
                      dom(ib)%u(i,dom(ib)%jsu-1,k) = ubs
	              end if
                  end do
                end do 
              end if
              if (dom(ib)%jnext.lt.0) then
                do k=1,ttk
                  do i=1,tti
	              if (L_LSMbase .or. (L_LSM.and..not.lrestart)) then
                if(abs(grz).gt.0.0001d00)then
		          if (dom(ib)%zc(k).gt.length) then  
                        dom(ib)%u(i,dom(ib)%jeu+1,k) = 0.0 
		          endif
                elseif(abs(gry).gt.0.0001d00)then
		          if (dom(ib)%yc(dom(ib)%jeu+1).gt.length) then  
                        dom(ib)%u(i,dom(ib)%jeu+1,k) = 0.0 
		          endif
                endif ! gry, grz values
	              else
                       dom(ib)%u(i,dom(ib)%jeu+1,k) = ubn
                    end if
                  end do
                end do 
              end if
!.............U=> Bottom and Top .................
              if (dom(ib)%kprev.lt.0) then
                 do j=1,ttj
                    do i=1,tti
                       dom(ib)%u(i,j,dom(ib)%ksu-1) = ubb
                    end do
                 end do 
              end if
              if (dom(ib)%knext.lt.0) then
                do j=1,ttj
                  do i=1,tti
	              if (L_LSMbase.or. (L_LSM.and..not.lrestart)) then 
                      dom(ib)%u(i,j,dom(ib)%keu+1) = 0.0 
	              else
                      dom(ib)%u(i,j,dom(ib)%keu+1) = ubt
	              end if
                  end do
                end do 
              end if
!........... V=> West and East ....................
              if (dom(ib)%iprev.lt.0) then
                 do k=1,ttk
                    do j=1,ttj
                       dom(ib)%v(dom(ib)%isv-1,j,k)    = vb
                    end do
                 end do 
              end if
              if (dom(ib)%inext.lt.0) then
                 do k=1,ttk
                    do j=1,ttj
                       dom(ib)%v(dom(ib)%iev+1,j,k)  = vb
                    end do
                 end do 
              end if
!............V=> South and North ....................
              if (dom(ib)%jprev.lt.0) then
                 do k=1,ttk
                    do i=1,tti
                       dom(ib)%v(i,dom(ib)%jsv-1,k)    = vb
                    end do
                 end do 
              end if
              if (dom(ib)%jnext.lt.0) then
                 do k=1,ttk
                    do i=1,tti
                       dom(ib)%v(i,dom(ib)%jev+1,k)  = vb
                    end do
                 end do 
              end if
!.............V=> Bottom and Top .................
              if (dom(ib)%kprev.lt.0) then
                 do j=1,ttj
                    do i=1,tti
                       dom(ib)%v(i,j,dom(ib)%ksv-1)    = vb
                    end do
                 end do 
              end if
              if (dom(ib)%knext.lt.0) then
                 do j=1,ttj
                    do i=1,tti
                       dom(ib)%v(i,j,dom(ib)%kev+1)  = vb
                    end do
                 end do 
              end if
!........... W=> West and East ....................
              if (dom(ib)%iprev.lt.0) then
                 do k=1,ttk
                    do j=1,ttj
                       dom(ib)%w(dom(ib)%isw-1,j,k)    = wb
                    end do
                 end do 
              end if
              if (dom(ib)%inext.lt.0) then
                 do k=1,ttk
                    do j=1,ttj
                       dom(ib)%w(dom(ib)%iew+1,j,k)  = wb
                    end do
                 end do 
              end if
!............W=> South and North ....................
              if (dom(ib)%jprev.lt.0) then
                 do k=1,ttk
                    do i=1,tti
                       dom(ib)%w(i,dom(ib)%jsw-1,k)    = wb
                    end do
                 end do 
              end if
              if (dom(ib)%jnext.lt.0) then
                 do k=1,ttk
                    do i=1,tti
                       dom(ib)%w(i,dom(ib)%jew+1,k)  = wb
                    end do
                 end do 
              end if
!.............W=> Bottom and Top .................
              if (dom(ib)%kprev.lt.0) then
                 do j=1,ttj
                    do i=1,tti
                       dom(ib)%w(i,j,dom(ib)%ksw-1)    = wb
                    end do
                 end do 
              end if
              if (dom(ib)%knext.lt.0) then
                 do j=1,ttj
                    do i=1,tti
                       dom(ib)%w(i,j,dom(ib)%kew+1)  = wb
                    end do
                 end do 
              end if


!.######### U=> Power law inlet conditions .##########
!NEED TO ADAPT THESE TO LSM
          IF (dom(ib)%bc_west.eq.12 .or. 
     &  (dom(ib)%bc_west.eq.8 .and. UPROF.eq.12)) THEN	
           do i = 1,tti ; do j = 1,ttj ; 	 do k = 1,ttk 
	  if (dom(ib)%yc(j).lt.((yen-yst)/2)) then
             dom(ib)%u(i,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	      *(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1.d0/7.d0)
 	  else
             dom(ib)%u(i,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	 *(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
	  endif
            dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k)*(1.0d0+1.0d0/7.0d0)
     &	 *(DABS(dom(ib)%zc(k)/(zen-zst)))**(1.d0/7.d0)
	     enddo ; end do ;  end do
          END IF
          IF (dom(ib)%bc_west.eq.13 .or. 
     &  (dom(ib)%bc_west.eq.8 .and. UPROF.eq.13)) THEN		
           do i = 1,tti ; do j = 1,ttj ; 	 do k = 1,ttk 
	  if (dom(ib)%yc(j).lt.((yen-yst)/2)) then
             dom(ib)%u(i,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	      *(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1.d0/7.d0)
 	  else
             dom(ib)%u(i,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	 *(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
	  endif
	     enddo ; end do ;  end do
          END IF
          IF (dom(ib)%bc_west.eq.14 .or. 
     &  (dom(ib)%bc_west.eq.8 .and. UPROF.eq.14)) THEN		
           do i = 1,tti ; do j = 1,ttj ; 	 do k = 1,ttk 
            dom(ib)%u(i,j,k) = ubulk*(1.d0+1.d0/7.d0)
     &  *(DABS(dom(ib)%zc(k)/(zen-zst)))**(1.d0/7.d0)
	     enddo ; end do ;  end do
          END IF
          IF (dom(ib)%bc_west.eq.15.or. 
     &  (dom(ib)%bc_west.eq.8 .and. UPROF.eq.15)) THEN			
           do i = 1,tti ; do j = 1,ttj ; 	 do k = 1,ttk 
	      dom(ib)%u(i,j,k) = 0.0187d0*
     &  (1.d0/0.41d0*DLOG(DABS(dom(ib)%zc(k))*0.0187d0*10**6)+5.d0)
	     enddo ; end do ;  end do
          END IF

          IF (dom(ib)%bc_west.eq.17.or. 
     &  (dom(ib)%bc_west.eq.8 .and. UPROF.eq.17)) THEN
           do i = 1,tti ; do j = 1,ttj  ; do k = 1,ttk
 		ufric=0.103594d0*ubulk+0.00568d0 !Fit of a power law distribution
                 dom(ib)%u(i,j,k)=ufric*
     &  (1.d0/0.41d0*LOG(ABS(dom(ib)%zc(k))*ufric*10**6))/LOG(10.0d0)

	     enddo ; end do ;  end do
           end if

       end if    !!!!!!!!!!!non-restart

!.################################################
!.###########   Allocate time series    ##########
	 do i=1,n_unstpt
	  if(dom_id(ib).eq.id_unst(i)) then
		jtime=itime_end-ntime
	     if (ntime*dt.lt.t_start_averaging1) then
			jtime=itime_end-INT(t_start_averaging1/dt)	
	     endif
		write(numpt,'(i3)') i
		unst='unst_'//trim(adjustl(numpt))//'.dat'
		idfile=499+i
		open (unit=idfile, file=unst)
		write(x,'(F6.3)')dom(ib)%x(i_unst(i))
		write(y,'(F6.3)')dom(ib)%y(j_unst(i))
		write(z,'(F6.3)')dom(ib)%z(k_unst(i))
		write(idfile,*)'ZoneT = "Time series:',x,y,z,'"'
		write(idfile,*)'Variables = ctime,u,v,w'
	  endif
	 enddo
      end do!ib
!.################################################
!.###########  Synthetic Eddy Method   ##########
          IF (bc_w.eq.8 .and. myrank.eq.0) then
		print*,'Writing the SEM inlet'
!		call SEM  !Generate the files for the inlet turbulent field
		print*,'Finish the SEM inlet'
	    ENDIF

70      format (10e25.8)
71      format (3F15.6)

        end subroutine 
!##########################################################################
