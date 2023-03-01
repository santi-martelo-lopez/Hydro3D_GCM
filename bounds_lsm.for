!######################################################################
      subroutine bound_LSM(op)
!######################################################################
      use vars
      use LSM
      use multidata

      implicit none
      integer :: I,J,K,np,nq,nr,op,ib,ipl
      double precision, pointer, dimension(:,:,:) :: FI
      real, parameter :: PI = 4 * atan (1.0)	
      double precision :: omega,k_w,c,h_wave
      !
      !
        do ib=1,nbp

       select case (op)
        Case (14)  
          fi => dom(ib)%PHI
        Case (15)  
          fi => dom(ib)%PHI_REINIT
        Case (16)  
          fi => dom(ib)%phi_new
        Case (17)  
          fi => dom(ib)%phi_init
        Case (18)  
          fi => dom(ib)%dens
        Case (19)  
          fi => dom(ib)%mu
       end select 

          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k  
!
!..... Boundary Conditions for PHI, zero gradient..
!
!=== EAST ===
      IF (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5) THEN   
       DO k = dom(ib)%ksp-pl, dom(ib)%kep+pl
        DO j = dom(ib)%jsp-pl, dom(ib)%jep+pl
	  DO ipl = 1,pl
	    FI (dom(ib)%iep+ipl,j,k) = FI(dom(ib)%iep+ipl-1,j,k)
      !
      !
!Fix the level set funtion at the outlet of the domain:
          if (dom(ib)%bc_east.ne.4.or. L_beach) then
           if(op.le.17) then ! level set functions
          !
          if(abs(grz).gt.0.0001d00)then
                 if (dom(ib)%zc(k).lt.length)   then
	    FI(dom(ib)%iep+ipl,j,k) = 1.0*abs(dom(ib)%zc(k)-length)
                 else if (dom(ib)%zc(k).gt.length)   then
	    FI(dom(ib)%iep+ipl,j,k) = -1.0*abs(dom(ib)%zc(k)-length)
                 else if (dom(ib)%zc(k).eq.length)  then
	    FI(dom(ib)%iep+ipl,j,k) = 0.0
                 end if  
          elseif(abs(gry).gt.0.0001d00)then
            if (dom(ib)%yc(j).lt.length)   then
      FI(dom(ib)%iep+ipl,j,k) = 1.0*abs(dom(ib)%yc(j)-length)
                       else if (dom(ib)%yc(j).gt.length)   then
      FI(dom(ib)%iep+ipl,j,k) = -1.0*abs(dom(ib)%yc(j)-length)
                       else if (dom(ib)%yc(j).eq.length)  then
      FI(dom(ib)%iep+ipl,j,k) = 0.0
                       end if  
            endif
                 !
           endif
          endif
      !
      !
!Fix the level set funtion at the outlet of the domain:
	  if(op.le.17) then ! level set functions
         if (dom(ib)%bc_east.eq.2 .or.dom(ib)%bc_east.eq.21) then
          !
          !
          if(abs(grz).gt.0.0001d00)then
                 if (dom(ib)%zc(k).lt.length)   then
	    FI(dom(ib)%iep+ipl,j,k) = 1.0*abs(dom(ib)%zc(k)-length)
                 else if (dom(ib)%zc(k).gt.length)   then
	    FI(dom(ib)%iep+ipl,j,k) = -1.0*abs(dom(ib)%zc(k)-length)
                 else if (dom(ib)%zc(k).eq.length)  then
	    FI(dom(ib)%iep+ipl,j,k) = 0.0
                 end if  
          elseif(abs(gry).gt.0.0001d00)then
            if (dom(ib)%yc(j).lt.length)   then
      FI(dom(ib)%iep+ipl,j,k) = 1.0*abs(dom(ib)%yc(j)-length)
                       else if (dom(ib)%yc(j).gt.length)   then
      FI(dom(ib)%iep+ipl,j,k) = -1.0*abs(dom(ib)%yc(j)-length)
                       else if (dom(ib)%yc(j).eq.length)  then
      FI(dom(ib)%iep+ipl,j,k) = 0.0
                       end if  
            endif
                 !
                 !
  	  endif; endif

	  ENDDO
        END DO  
       END DO
      END IF
!
!=== WEST ===
      IF (dom(ib)%iprev.lt.0 .and. dom(ib)%bc_west.ne.5) THEN
       DO k = dom(ib)%ksp-pl, dom(ib)%kep+pl
        DO j = dom(ib)%jsp-pl, dom(ib)%jep+pl 
	  DO ipl = 1,pl
	    FI(dom(ib)%isp-ipl,j,k) = FI(dom(ib)%isp-ipl+1,j,k)

!=====Aristos CHristou Wave Theories ============================ 

	  if(op.le.17 .and. dom(ib)%bc_west.eq.31) then                        ! level set functions
            !
            !
                if(abs(grz).gt.0.0001d00)then
                !
	    if (theory.eq.1) then                                          !Linear waves
             k_w=2.0d0*PI/L_w
             omega = sqrt(9.81d0*k_w*tanh(k_w*length))                     !======== k_w: wavenumber, omega:freq
             h_wave = length+H_w/2.0d0*cos(omega*ctime-PI/2.0d0)


	    elseif (theory.eq.2) then                                      !2nd Order Stokes
             k_w=2.0d0*PI/L_w
             omega = sqrt(9.81d0*k_w*tanh(k_w*length))
             h_wave = length+H_w/2.0d0*cos(omega*ctime-PI/2.0d0)     
     & +(H_w**2.0*k_w/16.0d0)*(3.0d0/((tanh(k_w*length))**3.0)
     & -1.0d0/(tanh(k_w*length)))*cos(2.0d0*(omega*ctime-PI/2.0d0))

	    elseif (theory.eq.3) then                                      !Solitary Waves
             k_w=sqrt(3.0d0*H_w/(4.0d0*length**3.0))                
             c=sqrt(9.81*(H_w+length))			       
             h_wave = (H_w/(cosh(k_w*(-c*ctime)+tau))**2.0)+length	   
            endif
                !
                if (dom(ib)%zc(k).lt.h_wave)   then
            FI(dom(ib)%isp-ipl,j,k) = 1.0*abs(dom(ib)%zc(k)-h_wave)
		else if (dom(ib)%zc(k).gt.h_wave)   then
            FI(dom(ib)%isp-ipl,j,k) = -1.0*abs(dom(ib)%zc(k)-h_wave)
                else if (dom(ib)%zc(k).eq.h_wave)  then
            FI(dom(ib)%isp-ipl,j,k) = 0.0
                end if
            !
            elseif(abs(gry).gt.0.0001d00)then
            !
	    if (theory.eq.1) then                                          !Linear waves
             k_w=2.0d0*PI/L_w
             omega = sqrt(9.81d0*k_w*tanh(k_w*length))                     !======== k_w: wavenumber, omega:freq
             h_wave = length+H_w/2.0d0*cos(omega*ctime-PI/2.0d0)


	    elseif (theory.eq.2) then                                      !2nd Order Stokes
             k_w=2.0d0*PI/L_w
             omega = sqrt(9.81d0*k_w*tanh(k_w*length))
             h_wave = length+H_w/2.0d0*cos(omega*ctime-PI/2.0d0)     
     & +(H_w**2.0*k_w/16.0d0)*(3.0d0/((tanh(k_w*length))**3.0)
     & -1.0d0/(tanh(k_w*length)))*cos(2.0d0*(omega*ctime-PI/2.0d0))

	    elseif (theory.eq.3) then                                      !Solitary Waves
             k_w=sqrt(3.0d0*H_w/(4.0d0*length**3.0))                
             c=sqrt(9.81*(H_w+length))			       
             h_wave = (H_w/(cosh(k_w*(-c*ctime)+tau))**2.0)+length	   
            endif
            !
                    if (dom(ib)%yc(j).lt.h_wave)   then
                FI(dom(ib)%isp-ipl,j,k) = 1.0*abs(dom(ib)%yc(j)-h_wave)
            else if (dom(ib)%yc(j).gt.h_wave)   then
                FI(dom(ib)%isp-ipl,j,k) = -1.0*abs(dom(ib)%yc(j)-h_wave)
                    else if (dom(ib)%yc(j).eq.h_wave)  then
                FI(dom(ib)%isp-ipl,j,k) = 0.0
                    endif
            endif ! gry, grz values
          endif ! op.le.17 .and. dom(ib)%bc_west.eq.31
                !
	  ENDDO
        END DO
       END DO
      END IF
!
!=== BOTTOM ===
      IF (dom(ib)%kprev.lt.0 .and. dom(ib)%bc_bottom.ne.5) THEN 
       DO j = dom(ib)%jsp-pl, dom(ib)%jep+pl
        DO i = dom(ib)%isp-pl, dom(ib)%iep+pl
	  DO ipl = 1,pl
	    FI (i,j,dom(ib)%ksp-ipl) = FI(i,j,dom(ib)%ksp-ipl+1)
	  ENDDO
        END DO
       END DO
      END IF
!
!=== TOP ===   
      IF (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5) THEN 
       DO j = dom(ib)%jsp-pl, dom(ib)%jep+pl
        DO i = dom(ib)%isp-pl, dom(ib)%iep+pl  
	  DO ipl = 1,pl
	    FI (i,j,dom(ib)%kep+ipl) = FI(i,j,dom(ib)%kep+ipl-1)
	  ENDDO
        END DO
       END DO   
      END IF
!
!=== SOUTH ===
      IF (dom(ib)%jprev.lt.0 .and. dom(ib)%bc_south.ne.5) THEN
       DO k = dom(ib)%ksp-pl, dom(ib)%kep+pl   
        DO i = dom(ib)%isp-pl, dom(ib)%iep+pl    
	  DO ipl = 1,pl
	    FI (i,dom(ib)%jsp-ipl,k) = FI(i,dom(ib)%jsp-ipl+1,k)
	  ENDDO
        END DO
       END DO
      END IF
!
!=== NORTH ===  
      IF (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5) THEN 
       DO k = dom(ib)%ksp-pl, dom(ib)%kep+pl   
        DO i = dom(ib)%isp-pl, dom(ib)%iep+pl     
	  DO ipl = 1,pl
	    FI (i,dom(ib)%jep+ipl,k) = FI(i,dom(ib)%jep+ipl-1,k)
	  ENDDO
        END DO
       END DO   
      END IF

      end do

	call exchange(op) 

      RETURN
      end subroutine bound_LSM
!#####################################################################
