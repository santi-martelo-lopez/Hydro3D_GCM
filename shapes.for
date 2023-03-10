!######################################################################
      SUBROUTINE imb_hemisphere(numIB)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: thc,PI,thz(80),zccos,zcsin
      INTEGER      :: M,L,K,c,nxhem,nyhem,p,Lt,Lrow1,Lrow2
      INTEGER      :: strlen2,nzr(80),izr,maxnzr
      CHARACTER*8  :: char_block2
      CHARACTER*31 :: gridfile
      DOUBLE PRECISION,allocatable,dimension (:) :: ztemp_layer
      DOUBLE PRECISION,allocatable,dimension (:,:) :: Rtemp_layer      
      INTEGER,allocatable,dimension (:) :: ctot_layer,nodes_layer
      INTEGER,allocatable,dimension (:,:) :: nodes_percyl_layer

      PI = 4.D0*DATAN(1.D0)

         write(char_block2,'(I8)') myrank
         strlen2=LEN(TRIM(ADJUSTL(char_block2)))
         char_block2=REPEAT('0',(3-strlen2))//TRIM(ADJUSTL(char_block2))
         gridfile='geom_HemiSphere_'//TRIM(ADJUSTL(char_block2))//'.dat'
         open (unit=2, file=gridfile)

	maxnzr=0 ;   maxnode = 0 ; M=numIB
	
	nzr(M) = nint((2.d0*PI*R(M)/2.d0)/dxm) !Number of planes
	thz(M) = PI/(nzr(M))
	maxnzr=max(maxnzr,nzr(M))
      write (2,*) 'variables=x,y,z'

	allocate (ztemp_layer(maxnzr),ctot_layer(maxnzr))
      allocate (nodes_layer(10000),Rtemp_layer(maxnzr,100))
	allocate (nodes_percyl_layer(maxnzr,100))

         nodes_percyl_layer = 0 ;    nodes_layer = 0 ;nodes_layer=0
                  
	  M=numIB
	     nodes(M) = 0
	     
	 do izr = nzr(M)/2+1,nzr(M)	 	  !*****
	      c =1
         if (izr.eq.1) then
	   Rtemp_layer(izr,c) = 0.d0
         else
	   Rtemp_layer(izr,c) = R(M)*cos(thz(M)*(izr-1)-(PI/2.d0))
	   ztemp_layer(izr) = R(M)*sin(thz(M)*(izr-1)-(PI/2.d0))
     &		+dzm*0.49999d0
         end if
	      do while (Rtemp_layer(izr,c).ge.0.)		!gt!!!
                nodes_layer(izr) = nodes_layer(izr) + 
     &  NINT(2.0*PI*Rtemp_layer(izr,c)/dxm)
                nodes_percyl_layer(izr,c) = 
     &  NINT(2.0*PI*Rtemp_layer(izr,c)/dxm)
               if (Rtemp_layer(izr,c).eq.0.) then
                nodes_layer(izr) = nodes_layer(izr) + 1
                nodes_percyl_layer(izr,c) = 1
               end if
               c = c + 1
               Rtemp_layer(izr,c) = Rtemp_layer(izr,c-1)- dxm              
             if(c.ge.cmax(M)) goto 555  
	      end do

555	CONTINUE
              ctot_layer(izr) = 1 !c - 1
	      if(ctot_layer(izr) .gt. 100) then
	      print*, 'allocate problem in Rtemp_layer'
	      end if
          end do      	    

	M=numIB ; Cz(M)=0.d0
!!!!!!!!!!!!!!!!!!! Bed of hemispheres --> uncomment the following :
	Cx(M)=0.08d0	; nxhem=20
	Cy(M)=0.09d0	; nyhem=7

	 izr=nzr(M)/2+1   !*****
	    K = 1
	 do while(izr.le.nzr(M))
	   do c = 1,ctot_layer(izr)
	    thc = 2.d0*PI/nodes_percyl_layer(izr,c)
	     do L = 1,nodes_percyl_layer(izr,c)
	      zccos = cos(2.d0*PI*(L-1)/nodes_percyl_layer(izr,c)+thc)
	      zcsin = sin(2.d0*PI*(L-1)/nodes_percyl_layer(izr,c)+thc)
            nodex(M,K) = Cx(M) + Rtemp_layer(izr,c)*zccos
	      nodey(M,K) = Cy(M) + Rtemp_layer(izr,c)*zcsin
	      nodez(M,K) = dzm*0.499d0 + R(M)*sin(thz(M)*(izr-1)-(PI/2.d0))
	      write(2,89) nodex(M,K), nodey(M,K), nodez(M,K)
               K = K + 1     
	     end do
	    end do
            izr = izr + 1     
       end do


	 K=K-1 
 	      nodes(M) = K-1
  	      maxnode = max(maxnode,nodes(M))	  

	  L=0 ;	 Lt=K	!points of a hemisphere

! First row
	     do p=2,nxhem
		Do L=1,K
		 Lt=Lt+1
		 nodex(M,Lt)=nodex(M,L)+0.16d0*(p-1)		
		 nodey(M,Lt)=nodey(M,L)		
		 nodez(M,Lt)=nodez(M,L)	
	      write(2,89) nodex(M,Lt), nodey(M,Lt), nodez(M,Lt)	
		Enddo
	     enddo
	     Lrow1=Lt

! Second row
	     do p=1,nxhem-1
		Do L=1,K
		 Lt=Lt+1
		 nodex(M,Lt)=nodex(M,L)+0.16d0*(p-1)+0.08d0		
		 nodey(M,Lt)=nodey(M,L)+0.16d0+0.005d0		
		 nodez(M,Lt)=nodez(M,L)	
	      write(2,89) nodex(M,Lt), nodey(M,Lt), nodez(M,Lt)		
		Enddo
	     enddo
	     Lrow2=Lt

! 3rd row
	     do K=1,Lrow1
		 Lt=Lt+1
		 nodex(M,Lt)=nodex(M,K)
		 nodey(M,Lt)=nodey(M,K)+0.16d0*2.d0+0.005d0	
		 nodez(M,Lt)=nodez(M,K)
	      write(2,89) nodex(M,Lt), nodey(M,Lt), nodez(M,Lt)	
	     enddo
! 5th row
	     do K=1,Lrow1
		 Lt=Lt+1
		 nodex(M,Lt)=nodex(M,K)	
		 nodey(M,Lt)=nodey(M,K)+0.16d0*4.d0+0.005d0		
		 nodez(M,Lt)=nodez(M,K)
	      write(2,89) nodex(M,Lt), nodey(M,Lt), nodez(M,Lt)	
	     enddo
! 7th row
	     do K=1,Lrow1
		 Lt=Lt+1
		 nodex(M,Lt)=nodex(M,K)	
		 nodey(M,Lt)=nodey(M,K)+0.005d0+0.16d0*6.d0		
		 nodez(M,Lt)=nodez(M,K)
	      write(2,89) nodex(M,Lt), nodey(M,Lt), nodez(M,Lt)	
	     enddo

! 4th row
	     do K=Lrow1+1,Lrow2
		 Lt=Lt+1
		 nodex(M,Lt)=nodex(M,K)
		 nodey(M,Lt)=nodey(M,K)+0.16d0*2.d0	
		 nodez(M,Lt)=nodez(M,K)
	      write(2,89) nodex(M,Lt), nodey(M,Lt), nodez(M,Lt)	
	     enddo
! 6th row
	     do K=Lrow1+1,Lrow2
		 Lt=Lt+1
		 nodex(M,Lt)=nodex(M,K)	
		 nodey(M,Lt)=nodey(M,K)+0.16d0*4.d0		
		 nodez(M,Lt)=nodez(M,K)
	      write(2,89) nodex(M,Lt), nodey(M,Lt), nodez(M,Lt)	
	     enddo

 	      nodes(M) = Lt !K-1
  	      maxnode = max(maxnode,nodes(M))	  

	deallocate (ztemp_layer,ctot_layer,nodes_layer)
	deallocate (nodes_percyl_layer,Rtemp_layer)

      close (2)

	  write(6,*) 'Number of layers : ',maxnzr
        write(6,87) '** Hemisphere',numIB,', # of nodes : ',nodes(numIB)
            
   87 FORMAT (a,i2,a,i6)
   88 FORMAT (i15)
   89 FORMAT (3e20.5)

      RETURN
      end 
!######################################################################
      SUBROUTINE imb_dune(numIB)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: thc,PI,xtemp,ytemp
      INTEGER      :: P,M,L,I,K,c,nxlay,nylay,maxc,strlen2
      CHARACTER*8  :: char_block2
      CHARACTER*31 :: gridfile

         PI = 4.D0*DATAN(1.D0)       

         write(char_block2,'(I3)') numIB
         strlen2=LEN(TRIM(ADJUSTL(char_block2)))
         char_block2=REPEAT('0',(3-strlen2))//TRIM(ADJUSTL(char_block2))
         gridfile='geom_Dune_'//TRIM(ADJUSTL(char_block2))//'.dat'
         open (unit=2, file=gridfile)
         write (2,*) 'variables="x","y","z"'

	 maxnode = 0   ;  M=numIB   ; nodes(M) = 0 ; L = 1 
	 nylay=((yen-yst)/dym)  ;  nxlay=((xen-xst)/dxm) 
	
	write(6,*)'Dune, long:',nxlay,'trans:',nylay

	 Do P=1,nylay
        Do 	K=1,nxlay
	 xtemp=dxm*0.4999999d0+(K-1)*dxm ; ytemp=dym*0.4999999d0+(P-1)*dym

	 if (xtemp.le.0.08d0) then	!1 section
	      maxc=INT(0.16d0/dzm)+1
		do c=1,maxc
	       nodex(M,L)=xtemp
	 	 nodey(M,L)=ytemp
	 	 nodez(M,L)=dzm*0.4999999d0+(c-1)*dzm
         	 write(2,89) nodex(M,L),nodey(M,L),nodez(M,L)
		 L=L+1
		enddo
	 endif

	 if (xtemp.gt.0.08d0 .and. xtemp.le.0.40d0) then	!2 section
	      maxc=INT((0.20d0-(xtemp)*0.5d0)/dzm)+1
		do c=1,maxc
	       nodex(M,L)=xtemp
	 	 nodey(M,L)=ytemp
	 	 nodez(M,L)=(0.20d0-xtemp*0.5d0)-(c-1)*dzm
!	 	 nodez(M,L)=dzm*0.4999999d0+(c-1)*dzm
         	 write(2,89) nodex(M,L),nodey(M,L),nodez(M,L)
		 L=L+1
		enddo
	 endif

	 if (xtemp.gt.0.60d0 .and. xtemp.le.1.10d0) then !3 section
	      maxc=INT((xtemp*0.036d0-0.0216d0)/dzm)+1
		do c=1,maxc
	       nodex(M,L)=xtemp
	 	 nodey(M,L)=ytemp
	 	 nodez(M,L)=(xtemp*0.036d0-0.0216d0)-(c-1)*dzm
!	 	 nodez(M,L)=dzm*0.4999999d0+(c-1)*dzm
         	 write(2,89) nodex(M,L),nodey(M,L),nodez(M,L)
		 L=L+1
		enddo
	 endif

	 if (xtemp.gt.1.10d0 .and. xtemp.le.2.6d0) then	!4 section
	      maxc=INT((xtemp*0.088d0-0.0788d0)/dzm)+1
		do c=1,maxc
	       nodex(M,L)=xtemp
	 	 nodey(M,L)=ytemp
	 	 nodez(M,L)=(xtemp*0.088d0-0.0788d0)-(c-1)*dzm
!	 	 nodez(M,L)=dzm*0.4999999d0+(c-1)*dzm
         	 write(2,89) nodex(M,L),nodey(M,L),nodez(M,L)
		 L=L+1
		enddo
	 endif
	 if (xtemp.gt.2.60d0 .and. xtemp.le.3.12d0) then !5 section
	      maxc=INT((xtemp*0.01923d0+0.10d0)/dzm)+1
		do c=1,maxc
	       nodex(M,L)=xtemp
	 	 nodey(M,L)=ytemp
 	 	 nodez(M,L)=(xtemp*0.01923d0+0.10d0)-(c-1)*dzm
!	 	 nodez(M,L)=dzm*0.4999999d0+(c-1)*dzm
         	 write(2,89) nodex(M,L),nodey(M,L),nodez(M,L)
		 L=L+1
		enddo
	 endif
	 if (xtemp.gt.3.12d0) then	!6 section
	      maxc=INT((0.16d0-0.d0)/dzm)+1
		do c=1,maxc
	       nodex(M,L)=xtemp
	 	 nodey(M,L)=ytemp
	 	 nodez(M,L)=dzm*0.4999999d0+(c-1)*dzm
         	 write(2,89) nodex(M,L),nodey(M,L),nodez(M,L)
		 L=L+1
		enddo
	 endif
 	 Enddo
	Enddo
879	CONTINUE          
                     
           nodes(M) = L-1
           maxnode = max(maxnode,nodes(M))          
      close (2)
               
	  write(6,*) 'Number of layers : ',nylay
        write(6,87) '** Dune',numIB,', # of nodes  :  ',nodes(numIB)
            
   87 FORMAT (a,i2,a,i6)
   89 FORMAT (3e20.5)
      RETURN
      end
!######################################################################
      SUBROUTINE imb_cone(numIB)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: thc,PI,Rtemp
      INTEGER      ::M,L,I,K,c,nlay,maxc,strlen2,nodes_percyl
      CHARACTER*8  :: char_block2
      CHARACTER*31 :: gridfile

         PI = 4.D0*DATAN(1.D0)       

         write(char_block2,'(I3)') numIB
         strlen2=LEN(TRIM(ADJUSTL(char_block2)))
         char_block2=REPEAT('0',(3-strlen2))//TRIM(ADJUSTL(char_block2))
         gridfile='geom_Cone_'//TRIM(ADJUSTL(char_block2))//'.dat'
         open (unit=2, file=gridfile)
         write (2,*) 'variables="x","y","z"'

	 maxnode = 0   ;  M=numIB   
	 Rtemp = R(M); nodes(M) = 0 ; L = 1 
	 nlay=((zend(M)-zini(M))/dzm) 

        Do 	K=1,nlay
		Rtemp=(R(M)-R(M)/(zend(M)-zini(M))*K*dzm)/2.d0
		if(K*dzm.gt.zen) goto 879
         	maxc=INT(Rtemp/dxm)  

		Do c=1,maxc
            nodes_percyl = NINT(2.d0*PI*Rtemp/dxm)
            thc = 2.d0*PI/nodes_percyl

		 Do I=1,nodes_percyl
	 nodex(M,L)=Cx(M)+Rtemp*cos(2.d0*PI*(I-1)/nodes_percyl+thc+K/3)
       nodey(M,L)=Cy(M)+Rtemp*sin(2.d0*PI*(I-1)/nodes_percyl+thc+K/3)
       nodez(M,L)=zini(M)+dzm*(K-1)+0.49999d0*dzm
          write(2,89) nodex(M,L),nodey(M,L),nodez(M,L)
          L = L + 1     
		 Enddo

		Rtemp = Rtemp - dxm
		Enddo
	  Enddo

879	CONTINUE          
                     
           nodes(M) = L-1
           maxnode = max(maxnode,nodes(M))          
      close (2)
               
	  write(6,*) 'Number of layers : ',nlay
        write(6,87) '** CONE',numIB,', # of nodes  :  ',nodes(numIB)
            
   87 FORMAT (a,i2,a,i8)
   89 FORMAT (3e20.5)
      RETURN
      end
!#############################################################
      SUBROUTINE imb_square(numIB)
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB      
      DOUBLE PRECISION:: nodexmin,nodexmax,nodeymin,nodeymax
      INTEGER      :: M,L,nin,njn,I,J,k,strlen2,c,nlay,maxc,clay
      CHARACTER*8  :: char_block2
      CHARACTER*31 :: gridfile

         write(char_block2,'(I8)') myrank
         strlen2=LEN(TRIM(ADJUSTL(char_block2)))
         char_block2=REPEAT('0',(3-strlen2))//TRIM(ADJUSTL(char_block2))
         gridfile='geom_Squ_'//TRIM(ADJUSTL(char_block2))//'.dat'
         open (unit=2, file=gridfile)
         write (2,*) 'variables=x,y,z'
       
	  M=numIB   ;  maxnode = 0 ; L=0

         if (linfin(numIB).eq.1) then
           zini(M)=0.d0 ;   nlay=((zen-zst)/dzm)                             
         else if (linfin(numIB).eq.0) then
	   nlay=((zend(M)-zini(M))/dzm)               
         endif
         
         maxc=INT(R(M)/dxm)+1
         
         if(cmax(M).gt.maxc) cmax(M)=maxc
                  
	DO K=1,nlay
	 if(linfin(M).eq.0) then   !This adds the lids of the finite
	  if (K.le.cmax(M) .or. K.ge.(nlay-cmax(M)+1)) then	 
	   clay=maxc
	  else 
	   clay=cmax(M)
	  endif
	 else
	   clay=cmax(M)
	 endif
	 
	 

	if(axis(numIB).eq.2) then
	 Do c=1,clay
 	   nodexmin = Cx(M)+(c-1)*dxm  
 	   nodexmax = Cx(M)+2.0*(R(M)-(c-1)*dxm)
         nodeymin = Cz(M)+(c-1)*dzm 
         nodeymax = Cz(M)+2.0*(R(M)-(c-1)*dzm) 	

         nin=2*((R(M)-dxm*(c-1))/dxm)+1 !elements in x-direction
         njn=2*((R(M)-dzm*(c-1))/dzm)+1 !elements in y-direction		  
                  
          do I=1,nin				!Bottom horizontal edge
            L=L+1
            nodex(M,L)=nodexmin+dxm*(I-1)
            nodey(M,L)=zini(M) +dym*(K-1) 
            nodez(M,L)=nodeymin            
          enddo
          do J=2,njn				!Right vertical edge	
            L=L+1
            nodex(M,L)=nodexmax
            nodey(M,L)=zini(M) +dym*(K-1)  
            nodez(M,L)=nodeymin+dzm*(J-1)                         
          enddo
           do I=2,nin
            L=L+1
            nodex(M,L)=nodexmax-dxm*(I-1)
            nodey(M,L)=zini(M) +dym*(K-1)
            nodez(M,L)=nodeymax                        
           enddo
           do J=2,njn-1
            L=L+1
            nodex(M,L)=nodexmin
            nodey(M,L)=zini(M) +dym*(K-1) 
            nodez(M,L)=nodeymax-dzm*(J-1)                       
           enddo   
	 Enddo

	ELSE
	 Do c=1,clay	
 	   nodexmin = Cx(M)-R(M)+(c-1)*dxm  
 	   nodexmax = Cx(M)+R(M)-(c-1)*dxm 
         nodeymin = Cy(M)-R(M)+(c-1)*dym 
         nodeymax = Cy(M)+R(M)-(c-1)*dym 	
	
         nin=2*((R(M)-dxm*(c-1))/dxm)+1 
         njn=2*((R(M)-dym*(c-1))/dym)+1	 
                  
          do I=1,nin
            L=L+1
            nodex(M,L)=nodexmin+dxm*(I-1)
            nodey(M,L)=nodeymin
            nodez(M,L)=zini(M)+dzm*(K-1)            
          enddo
          do J=2,njn
            L=L+1
            nodex(M,L)=nodexmax
            nodey(M,L)=nodeymin+dym*(J-1)
            nodez(M,L)=zini(M)+dzm*(K-1)                         
          enddo
           do I=2,nin
            L=L+1
            nodex(M,L)=nodexmax-dxm*(I-1)
            nodey(M,L)=nodeymax
            nodez(M,L)=zini(M)+dzm*(K-1)                         
           enddo
           do J=2,njn-1
            L=L+1
            nodex(M,L)=nodexmin
            nodey(M,L)=nodeymax-dym*(J-1)
            nodez(M,L)=zini(M)+dzm*(K-1)                         
           enddo         
       Enddo !c

	ENDIF!axis

         ENDDO
    

         nodes(M)=L
         maxnode = max(maxnode,nodes(M))        
         
         Do L = 1,nodes(M)
            write (2,89) nodex(M,L),nodey(M,L),nodez(M,L)        
         End do  

      close (2)

   89 FORMAT (3e20.5)

      RETURN
      end
!######################################################################
      SUBROUTINE imb_cylinder(numIB)
!######################################################################
! This subroutine genereates a cylinder extruded along iether the OX, OY, or OZ axis.
! If the user want the cylinder completely filled, it is advides that cmax(M), the number of concentric layers of the sylinder defined "geo.cin" is 
! strictly bigger thant variable maxc = INT(R(M)/dxm) + 1   

      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: thc,PI,zccos,zcsin
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: Rtemp
	! Array that stores the radius of the exterior and inner layers

      INTEGER(kind = 8), DIMENSION(:), ALLOCATABLE::nodes_percyl
	! Array that stores the number markers in the exterior and inner layers

      INTEGER      ::M,L,I,K,c,ctot,nlay,maxc,strlen2
      INTEGER*8	   :: nodesPerLayer
      INTEGER(kind = 8)	   :: nodesPerStation
      CHARACTER*8  :: char_block2
      CHARACTER*31 :: gridfile

         PI = 4.D0*DATAN(1.D0)       

         write(char_block2,'(I8)') numIB
         strlen2=LEN(TRIM(ADJUSTL(char_block2)))
         char_block2=REPEAT('0',(3-strlen2))//TRIM(ADJUSTL(char_block2))
         gridfile='geom_Cyl_'//TRIM(ADJUSTL(char_block2))//'.dat'
         open (unit=2, file=gridfile)
         write (2,*) 'variables=x,y,z'

	! Initializing varaibles
	 maxnode = 0  ; M=numIB  
         nodes(M) = 0; c =1           
                    
	! Total number of layers possible to define
         maxc = INT(R(M)/dxm) + 1         		 
         if(cmax(M).gt.maxc) cmax(M)=maxc
		 
		 ! Case for infinie cylinder (this needs revision)
         if (linfin(M).eq.1) then
			zini(M)=zst 
			if(axis(M).eq.1) nlay=((xen-xst)/dxm)-1         
			if(axis(M).eq.2) nlay=((yen-yst)/dym)-1                             
			if(axis(M).eq.3) nlay=((zen-zst)/dzm)-1                                                
		endif
		
		! Case for user define finite cylinfer
		if (linfin(M).eq.0) then
			if(axis(M).eq.1) nlay=((zend(M)-zini(M))/dxm)+1           
			if(axis(M).eq.2) nlay=((zend(M)-zini(M))/dym)+1
			if(axis(M).eq.3) nlay=((zend(M)-zini(M))/dzm)+1
		endif		 		 
		 
	 write(*,*)  
	 write(*,*)                 
	 write(*,*) "Starting generation of cylinder geometry."

	 Allocate( MrkrsLy(M,nlay,maxc) )
	 MrkrsLy = 0
	 Allocate( Rtemp(cmax(M)), nodes_percyl(cmax(M)) )
	 Rtemp = R(M); nodes_percyl = 0
	 nodesPerStation = 0

            Do while (Rtemp(c).ge.0.d0 .and. c.le.maxc) 
               IF (c.eq.1) then
                nodes(M) = nodes(M) + NINT(2.d0*PI*Rtemp(c)/dxm)
                nodes_percyl(c) =     NINT(2.d0*PI*Rtemp(c)/dxm)
               else
                nodes(M) = nodes(M) + NINT(2.d0*PI*Rtemp(c)/(dxm))
                nodes_percyl(c)=NINT(2.0*PI*Rtemp(c)/(dxm))
               end if
			   
	      nodesPerLayer = NINT(2.d0*PI*Rtemp(c)/dxm)
		  do K = 1, nlay
			MrkrsLy(M,K,c) = nodesPerLayer
	      enddo
		!  write(*,*) M,K,c, MrkrsLy(M,nlay,c), nodesPerLayer
		  nodesPerStation = nodesPerStation + nodesPerLayer  
	      write(*,'(i11,A,i4)') nodesPerLayer," markers @ layer# ",c

              IF (Rtemp(c).eq.0.d0) nodes(M)        = nodes(M) + 1
              IF (Rtemp(c).eq.0.d0) nodes_percyl(c) = 1
              c = c + 1
              Rtemp(c) = Rtemp(c-1) - dxm                
           End do

		
		!!!!!!!!!!!!!!!!!!!!!
		! Lets initialize the global variable keeping track on the number of markers generated	 in any section
		allocate(NumStation(bodynum))
		NumStation(M) = nlay

		if(numIB .eq.1) then
			ibmSt_k = nlay
		else
			ibmSt_k = ibmSt_k + nlay
		endif
		ibmSt(M) = ibmSt_k
		write(*,*) "Number of Stations: ",ibmSt(M)
		write(*,*) "Total number of Stations: ",ibmSt_k
		 !!!!!!!!!!!!!!!!!!!!!
		 ! Lets initialize the global variable keemping track on the number of layers generated
		 allocate(NumLayers(bodynum,NumStation(M)))
		 do K = 1, nlay 
			NumLayers(M,k) = maxc
		 enddo
		 write(*,'(A,i11)') "   Number of layers: ", NumLayers(M,1)
		 !!!!!!!!!!!!!!!!!!!!!
		 ! Lets initialize the global variable keeping track on the number of markers generated	 in the outer most layer		   
		 allocate(MrkrsExLy(bodynum,NumStation(M)))
		 do K = 1, nlay
			MrkrsExLy(M,K) = nodes_percyl(1)
			if(M.eq.1) then
				ibmMkrsEL( k )  = nodes_percyl(1)
			else
				ibmMkrsEL( ibmSt(M-1) + k )  = nodes_percyl(1)
			endif
		 enddo
		 write(*,*) "Number of markers in exterior layer: ",ibmMkrsEL(ibmSt_k)
		 !!!!!!!!!!!!!!!!!!!!!
		 ! Lets initialize the global variable keeping track on the number of markers generated	 in any section
		 allocate(MrkrsSection (bodynum,NumStation(M)))		  
	   	 do K = 1, nlay
			MrkrsSection(M,K) = nodesPerStation
			if(M.eq.1) then
				ibmStMkrs( k )  = nodesPerStation
			else
				ibmStMkrs( ibmSt(M-1) + k )  = nodesPerStation
			endif
		 enddo
		  write(*,*) "Number of markers per station: ",ibmStMkrs(ibmSt_k)

           ctot = c - 1  
		   
           ! Calculating the number of stations
		   
	  L=1
	  Do K=1,nlay
		 if(linfin(M).eq.0) then   !This adds the lids of the finite
		 
		 ! ??? Check the purpose of this if conditionals
		  if (K.le.cmax(M) .or. K.ge.(nlay-cmax(M)+1)) then	 
		   ctot=maxc
		  else 
		   ctot=cmax(M)
		  endif
		 
		 else
		   ctot=cmax(M)
		endif
		
           Do c = 1,ctot
            thc = 2.d0*PI/nodes_percyl(c)
             Do I = 1,nodes_percyl(c)
	        zccos= cos(2.d0*PI*(I-1)/nodes_percyl(c)+thc)
	        zcsin= sin(2.d0*PI*(I-1)/nodes_percyl(c)+thc)
			if(axis(numIB).eq.1) then
			  nodez(M,L)= Cz(M) + Rtemp(c)*zccos
		        nodey(M,L)= Cy(M) + Rtemp(c)*zcsin
		        nodex(M,L)= zini(M)+dxm*(K-1)+dxm*0.499999d0
			elseif(axis(numIB).eq.2) then
			  nodex(M,L)= Cx(M) + Rtemp(c)*zccos
!		        nodey(M,L)= zini(M)+dym*(K-1)+dym*0.499999d0
		        nodey(M,L)= zini(M)+dym*(K-1)
		        nodez(M,L)= Cz(M) + Rtemp(c)*zcsin
			elseif(axis(numIB).eq.3) then
			  nodex(M,L)= Cx(M) + Rtemp(c)*zccos
		        nodey(M,L)= Cy(M) + Rtemp(c)*zcsin
		        nodez(M,L)= zini(M)+dzm*(K-1)+dzm*0.499999d0
			endif
              write(2,89) nodex(M,L),nodey(M,L),nodez(M,L)
              L = L + 1     
             End do  
            End do
			
			! Lets write the number of markers generated at each station ( crossection in the span-wise direction)
		    write(*,'(i11,A,i4)') (L-1)," markers @ station # ",K
			
          End Do
          
                     
           nodes(M) = L-1
           maxnode = max(maxnode,nodes(M))      
			write(*,*)                 
			write(*,*) "Markers generated for this particular cylinder: ",maxnode

      close (2)
      
           
        write(6,*) ' '
        write(6,87) '** CYLINDER',numIB,', # of nodes  :  ',nodes(numIB)

	Deallocate(Rtemp,nodes_percyl)
	 write(*,*)                 
	 write(*,*) "Generation of cylinder geometry completed."
	 write(*,*)  
	 write(*,*)  
            
   87 FORMAT (a,i2,a,i5)
   89 FORMAT (3e20.5)

      RETURN
      end

!######################################################################
      subroutine linespace(x, x_start, x_end, size)
!######################################################################
    ! =========================================================
    ! creates a 1-dimensional array with evenly spaced elements
    ! =========================================================
    ! This subroutine genereates a custom defined geometry
        use vars
        use multidata
        use imb
        use mpi
		implicit none
        integer:: nn, ii
		integer, intent(in) :: size ! number of nodes
        DOUBLE PRECISION, intent(in) :: x_start, x_end
        DOUBLE PRECISION, dimension(size), intent(out):: x
		DOUBLE PRECISION:: dh
		!
        dh = (x_end - x_start) / (size - 1)
		x = 0.0d00
		do nn = 1, size
			x(nn) = x_start + (nn-1)*dh
		enddo
        !x(1:x_length) = [(x_start + (ii-1)*dx,ii = 1, x_length)]
    
    
    	end

!######################################################################
      subroutine lines(x,z, x_start, x_end, 
     &				z_start, z_end, size)
!######################################################################
    ! =========================================================
    ! creates a 2 1-dimensional arrays with evenly spaced elements
    ! =========================================================
    ! This subroutine genereates a custom defined geometry
	 use vars
	 use multidata
	 use imb
	 use mpi
	 implicit none
	 integer:: nn, ii
	 integer, intent(in) :: size ! number of nodes
	 DOUBLE PRECISION, intent(in) :: x_start, x_end
	 DOUBLE PRECISION, intent(in) :: z_start, z_end
	 DOUBLE PRECISION, dimension(size), intent(out):: x,z
	 DOUBLE PRECISION:: dx, dz
	 dx = (x_end - x_start) / (size - 1)
	 dz = (z_end - z_start) / (size - 1)
	 x = 0.0d00
	 z = 0.0d00
	 do nn = 1, size
		 x(nn) = x_start + (nn-1)*dx
		 z(nn) = z_start + (z_end - z_start)*( (nn-1)*dx ) 
     &		/ (x_end - x_start) 
		 !z(nn) = z_start + (nn-1)*dz
	 enddo
	 !x(1:x_length) = [(x_start + (ii-1)*dx,ii = 1, x_length)]

	 end	subroutine	lines
!######################################################################
      SUBROUTINE imb_userDefinedGeo(numIB)
!######################################################################
! This subroutine genereates a custom defined geometry
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: thc,PI,zccos,zcsin
	  INTEGER:: nn, s_n, sk
      DOUBLE PRECISION:: x_l,z_l
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: x,z
      INTEGER(kind = 8), DIMENSION(:), ALLOCATABLE::nodes_percyl
	! Array that stores the number markers in the exterior and inner layers

      INTEGER      ::M,L,I,K,c,ctot,nlay,maxc,strlen2
      INTEGER*8	   :: nodesPerLayer
      INTEGER(kind = 8)	   :: nodesPerStation
      CHARACTER*8  :: char_block2
      CHARACTER*31 :: gridfile
         !
         PI = 4.D0*DATAN(1.D0)
         !
         write(char_block2,'(I8)') numIB
         strlen2=LEN(TRIM(ADJUSTL(char_block2)))
         char_block2=REPEAT('0',(3-strlen2))//TRIM(ADJUSTL(char_block2))
         gridfile='geom_udg_'//TRIM(ADJUSTL(char_block2))//'.dat'
         open (unit=2, file=gridfile)
         write (2,*) 'variables=x,y,z'
         !
	 !   Initializing varaibles
	 maxnode = 0  ; M=numIB
	 nodes(M) = 0; c =1
	 !   Total number of layers possible to define
	 maxc = 1
	 	! Case for user defined finite cylinfer
	 	if (linfin(M).eq.0) then
            ! lets compute the number of stations for this body
	 		if(axis(M).eq.1) nlay=((zend(M)-zini(M))/dxm)+1           
	 		if(axis(M).eq.2) nlay=((zend(M)-zini(M))/dym)+1
	 		if(axis(M).eq.3) nlay=((zend(M)-zini(M))/dzm)+1
        endif
     !
	 write(*,*)  
	 write(*,*)                 
	 write(*,*) "Starting generation of cylinder geometry."
	 Allocate( MrkrsLy(M,nlay,maxc) )
	 MrkrsLy = 0
	 !! lets specify the corner locations of each station
	 !ncrs = 5
	 Allocate( nodes_percyl(cmax(M)) )
	 !Allocate(cornersx(ncrs),cornersz(ncrs))
	 nodes_percyl = 0
	 nodesPerStation = 0
	 !cornersx(M,1) = 0.3
	 !cornersz(M,1) = 0.01
	 !cornersx(M,2) = 0.6
	 !cornersz(M,2) = 0.01
	 !cornersx(M,3) = 0.6
	 !cornersz(M,3) = 0.1
	 !cornersx(M,4) = 0.3
	 !cornersz(M,4) = 0.1
	 !cornersx(M,5) = 0.3
	 !cornersz(M,5) = 0.01
	 !
	 !debug
	 write(*,*) cornersx
	 write(*,*) cornersz
	 !
	 nodes(M) = 0
	 L = 0
	 do K=1,nlay ! sweeping stations
		nodesPerStation = 0
		!do c = 1,cmax(M) ! Sweeping layers
		nodes_percyl(c) = 0
		nodesPerLayer = 0
		do nn = 1, ncrs-1 ! sweeping corssection corners
			x_l = abs(cornersx(M,nn+1)-cornersx(M,nn))
			z_l = abs(cornersz(M,nn+1)-cornersz(M,nn))
			s_n = NINT(sqrt(x_l**2 + z_l**2 )/dxm)+1           ! warning, dxm
			!
			! checking if at least we will be creating 2 elements in x and z arrays
			if(s_n.ge.2)then
				!pass
			else
				s_n = s_n + 1
			endif
			!
			allocate(x(s_n),z(s_n))
			call linespace(x, cornersx(M,nn), cornersx(M,nn+1), s_n)
			call linespace(z, cornersz(M,nn), cornersz(M,nn+1), s_n)
			!			call lines(x,z, cornersx(M,nn), cornersx(M,nn+1), 
			!     &				cornersz(M,nn), cornersz(M,nn+1), s_n)
			!write(*,*) z
			do sk = 1,s_n-1   ! sweeping segement between current corners
				L = L + 1
				nodex(M,L)= Cx(M) + x(sk) !(sk-1)*x_l/(s_n-1) !
				nodey(M,L)= zini(M)+dym*(K-1)
				nodez(M,L)= Cz(M) + z(sk) !(sk-1)*z_l/(s_n-1) ! 
				nodes(M) = nodes(M) + 1
				nodes_percyl(c) = nodes_percyl(c) + 1
				nodesPerLayer = nodesPerLayer + 1
				write(2,89) nodex(M,L),nodey(M,L),nodez(M,L)
		  	enddo ! sweeping segement between current corners
		  deallocate(x,z)
		End do    ! sweeping corssection corners
		!
		nodesPerStation = nodesPerStation + nodesPerLayer  
		write(*,'(i11,A,i4)') nodesPerLayer," markers @ layer# ",c	  
		MrkrsLy(M,K,c) = nodesPerLayer
		!  write(*,*) M,K,c, MrkrsLy(M,nlay,c), nodesPerLayer
		!enddo ! sweeping layers
		! Lets write the number of markers generated at each station ( crossection in the span-wise direction)
		write(*,'(i11,A,i4)') (L)," markers @ station # ",K
	  enddo !sweeping stations
		!!!!!!!!!!!!!!!!!!!!!
		! Lets initialize the global variables keeping track on the number of markers generated	 in any section
		allocate(NumStation(bodynum))
		NumStation(M) = nlay
        !
		if(numIB .eq.1) then
			ibmSt_k = nlay
		else
			ibmSt_k = ibmSt_k + nlay    ! current number of stations for all the bodies so far
		endif
		ibmSt(M) = ibmSt_k  !array keeping track of the number of stations in each body
		write(*,*) "Number of Stations: ",ibmSt(M)
		write(*,*) "Total number of Stations: ",ibmSt_k
		 !!!!!!!!!!!!!!!!!!!!!
		 ! Lets initialize the global variable keemping track on the number of layers generated
		 allocate(NumLayers(bodynum,NumStation(M)))
		 do K = 1, nlay 
			NumLayers(M,k) = maxc
		 enddo
		 write(*,'(A,i11)') "   Number of layers: ", NumLayers(M,1)
		 !!!!!!!!!!!!!!!!!!!!!
		 ! Lets initialize the global variable keeping track on the number of markers generated	 in the outer most layer		   
		 allocate(MrkrsExLy(bodynum,NumStation(M)))
		 do K = 1, nlay
			MrkrsExLy(M,K) = nodes_percyl(1)
			if(M.eq.1) then
				ibmMkrsEL( k )  = nodes_percyl(1)
			else
				ibmMkrsEL( ibmSt(M-1) + k )  = nodes_percyl(1)
			endif
		 enddo
		 write(*,*) "Number of markers in exterior layer: ",ibmMkrsEL(ibmSt_k)
		 !!!!!!!!!!!!!!!!!!!!!
		 ! Lets initialize the global variable keeping track on the number of markers generated	 in any section
		 allocate(MrkrsSection (bodynum,NumStation(M)))		  
	   	 do K = 1, nlay
			MrkrsSection(M,K) = nodesPerStation
			if(M.eq.1) then
				ibmStMkrs( k )  = nodesPerStation
			else
				ibmStMkrs( ibmSt(M-1) + k )  = nodesPerStation
			endif
		 enddo
		 write(*,*) "Number of markers per station: ",ibmStMkrs(ibmSt_k)
		 maxnode = max(maxnode,nodes(M))
		 write(*,*)
	write(*,*) "Markers generated for this user defined geometry: ",maxnode
		 close (2)
		 write(6,*) ' '
		 write(6,87) '** CYLINDER',numIB,', # of nodes  :  ',nodes(numIB)
		 !
	Deallocate(nodes_percyl,cornersx,cornersz)
	 write(*,*)                 
	 write(*,*) "Generation of user defined geometry completed."
	 write(*,*)  
	 write(*,*)  
            
   87 FORMAT (a,i2,a,i5)
   89 FORMAT (3e20.5)

      RETURN
      end

!#############################################################
      SUBROUTINE imb_cube(numIB)
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: nodexmin,nodexmax,nodeymin,nodeymax
      DOUBLE PRECISION:: delta_L_x,delta_L_y,delta_L_z,nn
      DOUBLE PRECISION:: nodezmin,nodezmax
      INTEGER      :: M,L,nin,njn,nkn,I,J,k,strlen2
      CHARACTER*8  :: char_block2
      CHARACTER*31 :: gridfile

         write(char_block2,'(I8)') myrank
         strlen2=LEN(TRIM(ADJUSTL(char_block2)))
         char_block2=REPEAT('0',(3-strlen2))//TRIM(ADJUSTL(char_block2))
         gridfile='geom_Cube_'//TRIM(ADJUSTL(char_block2))//'.dat'
         open (unit=2, file=gridfile)
         write (2,*) 'variables="x","y","z"'
      
         maxnode = 0
! -------------------------------------------------
	 M=numIB
 	
!         nodexmin=Cx(M)-R(M); nodexmax=Cx(M)+R(M)
!         nodeymin=Cy(M)-R(M); nodeymax=Cy(M)+R(M)
!         nodezmin=Cz(M)-R(M); nodezmax=Cz(M)+R(M)

	 delta_L_x = R(M)
	 print*,delta_L_x
	 delta_L_y = 0.5d00*(zend(M) - zini(M)) 
	 print*,delta_L_y
	 delta_L_z = 5.0d00*R(M)
	 print*,delta_L_z

         nodexmin=Cx(M)-delta_L_x; nodexmax=Cx(M)+delta_L_x
	 print*,nodexmin
	 print*,nodexmax
         nodeymin=zini(M) ; nodeymax=zend(M)
	 print*,nodeymin
	 print*,nodeymax
         nodezmin=Cz(M)-delta_L_z; nodezmax=Cz(M)+delta_L_z
	 print*,nodezmin
	 print*,nodezmax

!         nin =2.d0*R(M)/dxm
	 nin = (nodexmax - nodexmin)/dxm
	 print*,nin
	 nn = (nodexmax - nodexmin)/dxm
	 print*,nn
         if(abs(nin-2.0*delta_L_x/dxm).ge.0.999999999)
     & then
            nin=nin+1
            print*,'absurd-x!!!!!'
         end if
         nin=nin+1
	 print*,'Number of nodes along the OX axis'
	 print*,nin


!         njn =2.d0*R(M)/dym
	 njn = NINT( (nodeymax - nodeymin)/dym)
	 print*,njn
	 nn = (nodeymax - nodeymin)/dym 
	 print*,nn
         if(abs(njn-2.0*delta_L_y/dym).ge.0.999999999)
     & then
            njn=njn+1
            print*,'absurd-y!!!!!'
         end if
         njn=njn+1
	 print*,'Number of nodes along the OY axis'
	 print*,njn


!         nkn =2.d0*R(M)/dzm
	 nkn =(nodezmax - nodezmin)/dzm
	 print*,nkn
	 nn = (nodezmax - nodezmin)/dzm
	 print*,nn
         if(abs(nkn-2.0*delta_L_z/dzm).ge.0.999999999)
     & then
            nkn=nkn+1
            print*,'absurd-z!!!!!'
         end if
         nkn=nkn+1
	 print*,'Number of nodes along the OZ axis'
	 print*,nkn

         nodes(M)= nkn * nin * njn
         maxnode = max(maxnode,nodes(M))

        Do k = 1,nkn
         Do I = 1, nin
          Do J = 1, njn 
             L = (k-1)*nin*njn + (I-1)*njn + J
              IF (I.eq.1) THEN 
               nodex(M,L) = nodexmin
              else IF (I.eq.nin) THEN 
               nodex(M,L) = nodexmax
              ELSE
               nodex(M,L) = nodex(M,L-njn) + dxm
              END IF
              IF (J.eq.1) THEN 
               nodey(M,L) = nodeymin
              else IF (J.eq.njn) THEN 
               nodey(M,L) = nodeymax
              ELSE
               nodey(M,L) = nodey(M,L-1) + dym
              END IF
              IF (K.eq.1) THEN 
               nodez(M,L) = nodezmin
              else IF (K.eq.nkn) THEN 
               nodez(M,L) = nodezmax
              ELSE
               nodez(M,L) = nodez(M,L-nin*njn) + dzm
              END IF
          End do         
         End do
        end do
         Do L = 1,nodes(M)
            write (2,89) nodex(M,L),nodey(M,L),nodez(M,L)        
         End do  

      close (2)

   88 FORMAT (i5)
   89 FORMAT (3e20.5)

      RETURN
      end
!######################################################################
      SUBROUTINE imb_sphere(numIB)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: thc,PI,thz(80),zccos,zcsin
      INTEGER      :: M,L,K,c
      INTEGER      :: strlen2,nzr(80),izr,maxnzr
      CHARACTER*8  :: char_block2
      CHARACTER*31 :: gridfile
      DOUBLE PRECISION,allocatable,dimension (:,:) :: ztemp_layer
      DOUBLE PRECISION,allocatable,dimension (:,:,:) :: Rtemp_layer      
      INTEGER,allocatable,dimension (:,:) :: ctot_layer,nodes_layer
      INTEGER,allocatable,dimension (:,:,:) :: nodes_percyl_layer


      PI = 4.D0*DATAN(1.D0)

         write(char_block2,'(I8)') myrank
         strlen2=LEN(TRIM(ADJUSTL(char_block2)))
         char_block2=REPEAT('0',(3-strlen2))//TRIM(ADJUSTL(char_block2))
         gridfile='geom_Sphere_'//TRIM(ADJUSTL(char_block2))//'.dat'
         open (unit=2, file=gridfile)

	maxnzr=0 ;   maxnode = 0
	M=numIB
	
	nzr(M) = nint((2.d0*PI*R(M)/2.d0)/dxm) !Number of planes
	thz(M) = PI/(nzr(M))
	maxnzr=max(maxnzr,nzr(M))
        write (2,*) 'variables="x","y","z"'

	allocate (ztemp_layer(1,maxnzr),ctot_layer(1,maxnzr))
        allocate (nodes_layer(1,10000),Rtemp_layer(1,maxnzr,100))
	allocate (nodes_percyl_layer(1,maxnzr,100))

         nodes_percyl_layer = 0 ;    nodes_layer = 0 ;nodes_layer=0
                  
	  M=numIB
	     nodes(M) = 0
	     
	 do izr = 1,nzr(M)	 
	      c =1
         if (izr.eq.1) then
	   Rtemp_layer(M,izr,c) = 0.d0
         else
	   Rtemp_layer(M,izr,c) = R(M)*cos(thz(M)*(izr-1)-(PI/2.d0))
	   ztemp_layer(M,izr) = R(M)*sin(thz(M)*(izr-1)-(PI/2.d0))+Cz(M)
         end if

	      do while (Rtemp_layer(M,izr,c).ge.0.)		!gt!!!
                nodes_layer(M,izr) = nodes_layer(M,izr) + 
     &  NINT(2.0*PI*Rtemp_layer(M,izr,c)/dxm)
                nodes_percyl_layer(M,izr,c) = 
     &  NINT(2.0*PI*Rtemp_layer(M,izr,c)/dxm)
               if (Rtemp_layer(M,izr,c).eq.0.) then
                nodes_layer(M,izr) = nodes_layer(M,izr) + 1
                nodes_percyl_layer(M,izr,c) = 1
               end if
               c = c + 1
               Rtemp_layer(M,izr,c) = Rtemp_layer(M,izr,c-1)- dxm
               
             if(c.ge.cmax(M)) goto 555 
               
	      end do

555	CONTINUE
              ctot_layer(M,izr) = c - 1
	      if(ctot_layer(M,izr) .gt. 100) then
	      print*, 'allocate problem in Rtemp_layer'
	      end if
          end do      	    

	M=numIB
	 izr=1
	    K = 1
	 do while(izr.le.nzr(M))
	   do c = 1,ctot_layer(M,izr)
	    thc = 2.d0*PI/nodes_percyl_layer(M,izr,c)
	     do L = 1,nodes_percyl_layer(M,izr,c)
	      zccos = cos(2.d0*PI*(L-1)/nodes_percyl_layer(M,izr,c)+thc)
	      zcsin = sin(2.d0*PI*(L-1)/nodes_percyl_layer(M,izr,c)+thc)
              nodex(M,K) = Cx(M) + Rtemp_layer(M,izr,c)*zccos
	      nodey(M,K) = Cy(M) + Rtemp_layer(M,izr,c)*zcsin
	      nodez(M,K) = Cz(M) + R(M)*sin(thz(M)*(izr-1)-(PI/2.d0))
	      write(2,89) nodex(M,K), nodey(M,K), nodez(M,K)
               K = K + 1     
	     end do
	    end do
            izr = izr + 1     
          end do
 	      nodes(M) = K-1
  	      maxnode = max(maxnode,nodes(M))	          

	deallocate (ztemp_layer,ctot_layer,nodes_layer)
	deallocate (nodes_percyl_layer,Rtemp_layer)

      close (2)

   88 FORMAT (i15)
   89 FORMAT (3e20.6)

      RETURN
      end 
!#############################################################
      SUBROUTINE imb_file(numIB)
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: PI,an
      INTEGER      :: L,nin,I,K,nlay,dummy,strlen
      CHARACTER*8  :: char_block2
      CHARACTER*31 :: gridfile
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION (:) :: xfile,yfile,zfile

       PI = 4.D0*DATAN(1.D0)

         write(char_block2,'(I3)') numIB
         strlen=LEN(TRIM(ADJUSTL(char_block2)))
         char_block2=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(char_block2))
         gridfile='geom_body_'//TRIM(ADJUSTL(char_block2))//'.dat'
         open (unit=2, file=gridfile)
!----      Load the file and proceed to interpolate     ----------
 	open(unit=1, name=filepoints(numIB))
	if (ibturbine(numIB).eq..true.)read(1,*)nin,nscatter(numIB)	!Mesh points	
	if (ibturbine(numIB).eq..false.)read(1,*)nin	!Mesh points	

	select case(axis(numIB))			!Brunho2015
	CASE (1) !-------------------------------------------
		allocate(yfile(nin),zfile(nin))
		    DO L=1,nin
			 read(1,*)dummy,yfile(L),zfile(L)
		    ENDDO
		close (1)
		maxnode=0  
		Cxor(numIB)=Cx(numIB) 
		Cyor(numIB)=Cy(numIB)
		Czor(numIB)=Cz(numIB)			           
           if (linfin(numIB).eq.1) then
             zini(numIB)=xst ;   nlay=((xen-xst)/dxm)!-1                             
           else if (linfin(numIB).eq.0) then
	     nlay=((zend(numIB)-zini(numIB))/dxm)!-1               
           endif							
		nodes(numIB)=nin*nlay*imbnumber(numIB)
		maxnode = max(maxnode,nodes(numIB))
  		K=1
  		do L=1,nlay  
   		 DO I=1,nin 
	 	  nodex(numIB,K)=zini(numIB)+dxm*(L)+1.d-7 
		  nodey(numIB,K)=yfile(I)+1.d-7
		  nodez(numIB,K)=zfile(I)+1.d-7 		 
     		  K=K+1	 
   		 ENDDO
   		Enddo
	CASE (2) !-------------------------------------------
		allocate(xfile(nin),zfile(nin))
		    DO L=1,nin
			 read(1,*)xfile(L),zfile(L)
		    ENDDO
		close (1)
		maxnode=0  
		Cxor(numIB)=Cx(numIB) 
		Cyor(numIB)=Cy(numIB)
		Czor(numIB)=Cz(numIB)
           if (linfin(numIB).eq.1) then
             zini(numIB)=yst ;   nlay=((yen-yst)/dym)!-1                             
           else if (linfin(numIB).eq.0) then
	     nlay=((zend(numIB)-zini(numIB))/dym)!-1               
           endif			
		nodes(numIB)=nin*nlay*imbnumber(numIB)
		 	maxnode = max(maxnode,nodes(numIB))
  		K=1
  		do L=1,nlay  
   		 DO I=1,nin 
	 	  nodex(numIB,K)=xfile(I)+1.d-7
		  nodey(numIB,K)=zini(numIB)+dym*(L)+1.d-7
		  nodez(numIB,K)=zfile(I)+1.d-7		 
     		  K=K+1	 
   		 ENDDO
   		Enddo
	CASE (3) !-------------------------------------------
		allocate(xfile(nin),yfile(nin))
		    DO L=1,nin
			 read(1,*)xfile(L),yfile(L)
		    ENDDO
		close (1)
		maxnode=0  
		Cxor(numIB)=Cx(numIB) 
		Cyor(numIB)=Cy(numIB)
		Czor(numIB)=Cz(numIB)
           if (linfin(numIB).eq.1) then
             zini(numIB)=zst ;   nlay=((zen-zst)/dzm)!-1  !it was-1                          
           else if (linfin(numIB).eq.0) then
	     nlay=((zend(numIB)-zini(numIB))/dzm)+1               
           endif
		nodes(numIB)=nin*nlay*imbnumber(numIB)
		maxnode = max(maxnode,nodes(numIB))
  		K=1
  		do L=1,nlay  
   		 DO I=1,nin 
	 	  nodex(numIB,K)=xfile(I) 
		  nodey(numIB,K)=yfile(I)
	 	  nodez(numIB,K)=zini(numIB)+dzm*(L)-0.49*dzm		 
     		  K=K+1	 
   		 ENDDO
   		Enddo
	 write(6,*)'Number of layers : ',nlay
	  CASE (-1) !-------------------------------------------
		allocate(xfile(nin),yfile(nin),zfile(nin))
		  DO L=1,nin
		   read(1,*)xfile(L),yfile(L),zfile(L)
		  ENDDO
		close (1)
		maxnode=0   ; nlay=1
		Cxor(numIB)=Cx(numIB) 
		Cyor(numIB)=Cy(numIB)
		Czor(numIB)=Cz(numIB)
		nodes(numIB)=nin*imbnumber(numIB)
		maxnode = max(maxnode,nodes(numIB))
   		DO I=1,nin 
	 	  nodex(numIB,I)=xfile(I)+1.d-12
		  nodey(numIB,I)=yfile(I)+1.d-12
		  nodez(numIB,I)=zfile(I)+1.d-12			  
   		ENDDO
	end select
 !--------Local and global coordinates  --------------------------------
	IF (ibturbine(numIB).eq..true.) then	!if body is a turbine 
	   an=pitch(numIB)*PI/180.d0	   !Angle of attack in radians
	   do i=1,nin*nlay		   !Rotate the body.
	    nodex(numIB,i)=nodex(numIB,i)-xaero(numIB)
	    nodey(numIB,i)=nodey(numIB,i)-yaero(numIB) 
	    nodez(numIB,i)=nodez(numIB,i)-zaero(numIB) !Gravity centre
!Rotated by the pitch angle
	 nodex(numIB,i)= nodex(numIB,i)*dcos(an)+nodey(numIB,i)*dsin(an)
	 nodey(numIB,i)=-nodex(numIB,i)*dsin(an)+nodey(numIB,i)*dcos(an)
!Set this as the local coordinates
	    nodexlocal(numIB,i)=nodex(numIB,i)    
	    nodeylocal(numIB,i)=nodey(numIB,i) 
	    nodezlocal(numIB,i)=nodez(numIB,i)
	   enddo
!From local to global coords applying the centre
	   do i=1,nin*nlay		
	    nodex(numIB,i)=nodex(numIB,i)+Cxor(numIB)
	    nodey(numIB,i)=nodey(numIB,i)+Cyor(numIB)
	    nodez(numIB,i)=nodez(numIB,i)+Czor(numIB)
	   enddo
	 ELSE	!All bodies but turbine
	  do i=1,nin*nlay		
	    nodexlocal(numIB,i)=nodex(numIB,i) !Set this as the local coordinates
	    nodeylocal(numIB,i)=nodey(numIB,i)
	    nodezlocal(numIB,i)=nodez(numIB,i)
	    nodex(numIB,i)=nodex(numIB,i)+Cxor(numIB)
	    nodey(numIB,i)=nodey(numIB,i)+Cyor(numIB)
	    nodez(numIB,i)=nodez(numIB,i)+Czor(numIB)
	   enddo   	   	    
	ENDIF

	if(imb_shape(numIB).eq.5) call imb_number(numIB)
	if(imb_shape(numIB).ne.5) print*,'subroutine not finished'
       close(2)

   88 FORMAT (i5)
   89 FORMAT (3e25.5)
      RETURN
      end        
!#############################################################
      SUBROUTINE imb_number(numIB)
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      INTEGER      :: I,K

	K=nodes(numIB)/imbnumber(numIB) !# of Lagrangians per unit
      write(6,*)' '
      write(6,'(a,i3,a)')'**********   THE BODY ',numIB,' :'
	write(6,*)'Total # ofnodes              :',nodes(numIB)
	write(6,*)'# of bodies                  :',imbnumber(numIB)
	IF(imb_shape(numIB).eq.5 .and. axis(numIB).eq.-1) 
     &		write(6,*)'No extrusion of the body !!!!'
!	IF(LSELFST(numIB))      write(6,*)'Turbine Self-Starting        :  YES'
!	IF(.not.LSELFST(numIB)) write(6,*)'Turbine Self-Starting        :  NO'
	write(2,*)'variables="x","y","z"'
!-----------------	1- body      --------------------------------
	if (imbnumber(numIB).eq.1) then
	  do i=1,K
	   nodex(numIB,i) = nodex(numIB,i) 
	   nodey(numIB,i) = nodey(numIB,i)
	   nodez(numIB,i) = nodez(numIB,i) 	
	   write (2,89) nodex(numIB,i),nodey(numIB,i),nodez(numIB,i)    
	  enddo
	endif
!-----------------	2- bodies      --------------------------------
      IF (imbnumber(numIB).eq.2) then
	 do i=1,K
	   nodex(numIB,i) = nodex(numIB,i) 
	   nodey(numIB,i) = nodey(numIB,i) + R(numIB)
	   nodez(numIB,i) = nodez(numIB,i) 		
	  write (2,89) nodex(numIB,i),nodey(numIB,i),nodez(numIB,i)    
	 enddo
	 do i=1,K
	  nodex(numIB,K+i)=-nodexlocal(numIB,i) 
	  nodey(numIB,K+i)=-nodeylocal(numIB,i)
	  nodex(numIB,K+i)=nodex(numIB,K+i) + Cxor(numIB) 
	  nodey(numIB,K+i)=nodey(numIB,K+i)+Cyor(numIB)-R(numIB)
	  nodez(numIB,K+i)=nodez(numIB,i) 	 
         write (2,89) nodex(numIB,K+i),nodey(numIB,K+i),nodez(numIB,K+i) 
	 enddo
      ENDIF
!-----------------	3- bodies      --------------------------------
	IF (imbnumber(numIB).eq.3) then
	do i=1,K
	  nodex(numIB,i) = nodex(numIB,i) 
	  nodey(numIB,i) = nodey(numIB,i) + R(numIB) 
	  nodez(numIB,i) = nodez(numIB,i) 		
	  write (2,89) nodex(numIB,i),nodey(numIB,i),nodez(numIB,i)  
	enddo
        do i=1,K
	  nodex(numIB,K+i)=nodexlocal(numIB,i)*(-0.5)-
     &                  nodeylocal(numIB,i)*(DSQRT(3.)/2)
	  nodey(numIB,K+i)=nodexlocal(numIB,i)*(DSQRT(3.)/2)+
     &                  nodeylocal(numIB,i)*(-0.5)
	  nodex(numIB,K+i)=nodex(numIB,K+i)+Cxor(numIB)-
     &			  R(numIB)*DSQRT(3.)/2 
	  nodey(numIB,K+i)=nodey(numIB,K+i)+Cyor(numIB)-R(numIB)*0.5 
	  nodez(numIB,K+i)=nodez(numIB,i) 	 
         write (2,89) nodex(numIB,K+i),nodey(numIB,K+i),nodez(numIB,K+i) 
        enddo
        do i=1,K
	 nodex(numIB,2*K+i)=nodexlocal(numIB,i)*(-0.5)-
     &                     nodeylocal(numIB,i)*(-DSQRT(3.)/2)
	 nodey(numIB,2*K+i)=nodexlocal(numIB,i)*(-DSQRT(3.)/2)+
     &                     nodeylocal(numIB,i)*(-0.5)
	 nodex(numIB,2*K+i)=nodex(numIB,2*K+i)+Cxor(numIB)+
     &				R(numIB)*DSQRT(3.)/2 
	 nodey(numIB,2*K+i)=nodey(numIB,2*K+i)+Cyor(numIB)-
     &				R(numIB)*0.5 
	 nodez(numIB,2*K+i)=nodez(numIB,i) 	 
         write (2,89) nodex(numIB,2*K+i),nodey(numIB,2*K+i)
     &   		,nodez(numIB,2*K+i) 
        enddo
       ENDIF
!-----------------	4- bodies      --------------------------------
	IF (imbnumber(numIB).eq.4) then
       do i=1,K
	nodex(numIB,i) = nodex(numIB,i)
	nodey(numIB,i) = nodey(numIB,i) + R(numIB)
	nodez(numIB,i) = nodez(numIB,i) 		
	write (2,89) nodex(numIB,i),nodey(numIB,i),nodez(numIB,i)   
       enddo
        do i=1,K
          nodex(numIB,K+i)=-nodeylocal(numIB,i)
          nodey(numIB,K+i)=nodexlocal(numIB,i)
          nodex(numIB,K+i)=nodex(numIB,K+i)+Cxor(numIB)-R(numIB) 
          nodey(numIB,K+i)=nodey(numIB,K+i) + Cyor(numIB) 
	  nodez(numIB,K+i)=nodez(numIB,i) 	 
	 write (2,89) nodex(numIB,K+i),nodey(numIB,K+i),nodez(numIB,K+i) 
       enddo
        do i=1,K
          nodex(numIB,2*K+i)=-nodexlocal(numIB,i)
          nodey(numIB,2*K+i)=-nodeylocal(numIB,i)
          nodex(numIB,2*K+i)=nodex(numIB,2*K+i) + Cxor(numIB) 
          nodey(numIB,2*K+i)=nodey(numIB,2*K+i) + Cyor(numIB)-R(numIB)
	  nodez(numIB,2*K+i)=nodez(numIB,i) 	            
         write (2,89) nodex(numIB,2*K+i),nodey(numIB,2*K+i)
     &   		,nodez(numIB,2*K+i) 
       enddo
        do i=1,K
          nodex(numIB,3*K+i)=nodeylocal(numIB,i)
          nodey(numIB,3*K+i)=-nodexlocal(numIB,i)
          nodex(numIB,3*K+i)=nodex(numIB,3*K+i) + Cxor(numIB)+R(numIB) 
          nodey(numIB,3*K+i)=nodey(numIB,3*K+i) + Cyor(numIB)  
         write (2,89) nodex(numIB,3*K+i),nodey(numIB,3*K+i)
     &   		,nodez(numIB,3*K+i) 
       enddo
	ENDIF

   88 FORMAT (i5)
   89 FORMAT (3e25.5)

	RETURN
	END SUBROUTINE 
!#############################################################
      SUBROUTINE imb_moved(numIB)
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      INTEGER, intent(in) :: numIB
      DOUBLE PRECISION:: PI
      INTEGER      :: L,I,K,strlen,strlen2,Geom_Time1,nplot    
      CHARACTER*8  :: char_block,num_ib
      CHARACTER*31 :: gridfile1     
       PI = 4.D0*DATAN(1.D0)

      IF (myrank.ne.master) RETURN

       Geom_Time1=301 

	nplot=1

	if(nplot.eq.2) then
         if (mod(itime,n_out).eq.0) then !Write in file when outputs are writen
           write(char_block,'(I6)') itime
           strlen=LEN(TRIM(ADJUSTL(char_block)))
           char_block=REPEAT('0',(6-strlen))//TRIM(ADJUSTL(char_block))
           write(num_ib,'(I2)') numIB
           strlen2=LEN(TRIM(ADJUSTL(num_ib)))
           num_ib=REPEAT('0',(2-strlen2))//TRIM(ADJUSTL(num_ib))
           gridfile1='Blade_'//TRIM(ADJUSTL(num_ib))//
     &      '_Time_'//TRIM(ADJUSTL(char_block))//'.plt'
           open (unit=Geom_Time1, file=gridfile1)
	   write(Geom_Time1,*)'title = points'
	   write(Geom_Time1,*)'variables=x,y,z'
	   write(Geom_Time1,*)
     &    'zone   ','i=  ',nodes(numIB),'DATAPACKING = POINT'
	 endif
 	endif

!- Distinguish between turbines (HA and VA) and non-turbines --------------
	IF (ibturbine(numIB) 	.eq..true.) then		! Turbine
!-------------------------------------------------------------------	
	IF (turax(numIB).eq.1) then	! Vertical Axis 
       do L=1,nodes(numIB)
      nodexlocal(numIB,L)=-R0(numIB,L)*dsin(rads(numIB)+alpha0(numIB,L))
      nodeylocal(numIB,L)= R0(numIB,L)*dcos(rads(numIB)+alpha0(numIB,L))                    
       enddo   

	  do L=1,nodes(numIB)
		nodex(numIB,L) = nodexlocal(numIB,L) + Cxor(numIB)
		nodey(numIB,L) = nodeylocal(numIB,L) + Cyor(numIB)
   	    if (mod(itime,n_out).eq.0) then
	if(nplot.eq.2) write(Geom_Time1,89) 
     &		nodex(numIB,L),nodey(numIB,L),nodez(numIB,L)   
    	    endif
	  enddo

	ENDIF !Vertical axis      
!-------------------------------------------------------------------		
	IF (turax(numIB).eq.2 .OR. turax(numIB).eq.3) then	! Horizontal Axis  ! AND Actuator line -Pablo 09/17
       do L=1,nodes(numIB)
       nodeylocal(numIB,L)=R0(numIB,L)*dsin(rads(numIB)+alpha0(numIB,L))
       nodezlocal(numIB,L)=R0(numIB,L)*dcos(rads(numIB)+alpha0(numIB,L))                    
       enddo   
	  do L=1,nodes(numIB)
		nodey(numIB,L) = nodeylocal(numIB,L) + Cyor(numIB)
		nodez(numIB,L) = nodezlocal(numIB,L) + Czor(numIB)
   	    if (mod(itime,n_out).eq.0) then
	if(nplot.eq.2) write(Geom_Time1,89) 
     &		nodex(numIB,L),nodey(numIB,L),nodez(numIB,L)   
    	    endif
	  enddo

	ENDIF !Horizontal axis          
	ENDIF !Turbine
!--------------------------------------------------
!_----- NON-TURBINE BODY which is moving   ------- 
       IF (ibturbine(numIB) .eq. .false.) then		! Non-Turbine
	do L=1,nodes(numIB)
	 nodex(numIB,L)=
     &   R0(numIB,L)*DSIN(rads(numIB)+alpha0(numIB,L))+Cxor(numIB)
	 nodey(numIB,L)=
     &   R0(numIB,L)*DCOS(rads(numIB)+alpha0(numIB,L))+Cyor(numIB)

!	 nodex(numIB,L)=nodexlocal(numIB,L)*dcos(rads(numIB))- 
!     &		 	nodeylocal(numIB,L)*dsin(rads(numIB))+Cxor(numIB)
!	 nodey(numIB,L)=nodexlocal(numIB,L)*dsin(rads(numIB))+ 
!     &		 	nodeylocal(numIB,L)*dcos(rads(numIB))+Cyor(numIB)
         if (mod(itime,n_out).eq.0) then
	if(nplot.eq.2)  write(Geom_Time1,89)
     &   nodex(numIB,L),nodey(numIB,L),nodez(numIB,L) 
	 endif
	enddo
	if(nplot.eq.2 .and. mod(itime,n_out).eq.0)  close(Geom_Time1)
       ENDIF	!Not a turbine
!---------------------------------------------------

   88 FORMAT (i5)
   89 FORMAT (3f25.5)
        RETURN
		END SUBROUTINE	



!#############################################################
      SUBROUTINE vertical_oscillating_cylinder(numIB)
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      INTEGER, intent(in) :: numIB
      DOUBLE PRECISION:: PI
      DOUBLE PRECISION:: amp_z, frec_z, disz
      INTEGER      :: L,I,K,strlen,strlen2,Geom_Time1,nplot    
      CHARACTER*8  :: char_block,num_ib
      CHARACTER*31 :: gridfile1     
	  !INTEGER:: gh_startTime = 0	  
	  PI = 4.D0*DATAN(1.D0)
   

      if(itime .gt. gh_startTime ) then        

      IF (myrank.eq.master) THEN

       Geom_Time1=301 

	nplot=1
		 
         if (mod(itime,n_out).eq.0) then !Write in file when outputs are writen
           write(char_block,'(I6)') itime
           strlen=LEN(TRIM(ADJUSTL(char_block)))
           char_block=REPEAT('0',(6-strlen))//TRIM(ADJUSTL(char_block))
           write(num_ib,'(I2)') numIB
           strlen2=LEN(TRIM(ADJUSTL(num_ib)))
           num_ib=REPEAT('0',(2-strlen2))//TRIM(ADJUSTL(num_ib))
           gridfile1='tecinscMotion'//TRIM(ADJUSTL(num_ib))//
     &      '_Time_'//TRIM(ADJUSTL(char_block))//'.plt'
           open (unit=Geom_Time1, file=gridfile1)
	   write(Geom_Time1,*)'title = points'
	   write(Geom_Time1,*)'variables=x,y,z'
           write(Geom_Time1,*) 'zone ','STRANDID=', 101,
     &    'SOLUTIONTIME=', ctime,
     &    ' i=',nodes(numIB),'zonetype=',
     & 	  'ordered',', DATAPACKING=point'

 
	 endif


	amp_z = 0.2d0
	frec_z = 2.0d0*PI*0.156
	disz = amp_z*DSIN( frec_z*(CTIME) )	! ! vertical displacement in oscillating cylinder case


	! Updating node velocities

	U_p = 0.0d0
	W_p = frec_z*amp_z*DCOS(frec_z*(CTIME) )
	V_p = 0.0d0

	if(itime .lt. (gh_startTime+2) ) then
        do L = 1, nodes(numIB)
	!allocate(nodezo(1,10000))
	!nodezo(numIB,L) = nodez(numIB,L)
	nodezlocal(numIB,L) = nodez(numIB,L)	
	enddo
	endif

	! Updating node location
	do L = 1, nodes(numIB)
	nodez(numIB,L) =  nodezlocal(numIB,L) + disz
	if (mod(itime,n_out).eq.0) then
		write(Geom_Time1,89) nodex(numIB,L),nodey(numIB,L),nodez(numIB,L)   
	endif						
	enddo

	if( mod(itime,n_out).eq.0)  close(Geom_Time1)



	if (myrank.eq.0) then
		 write(*,*) ctime, gh_startTime
		 write(*,*) U_p, W_p, disz
		 write(*,*) Cx(numIB),Cz(numIB) + disz
		 write(*,*) nodezlocal(numIB,1), nodez(numIB,1)
	endif

	ENDIF ! End of master processor

	else ! The immersed body is has not yet started to move

	if (myrank.eq.0) then
		U_p = 0.0d0
		W_p = 0.0d0
		V_p = 0.0d0
	 	write(*,*) 104, W_p, nodez(numIB,1)
	endif ! End of master processor

	endif

   88 FORMAT (i5)
   89 FORMAT (3f25.5)
		RETURN
!#############################################################
        END SUBROUTINE
!#############################################################

!#############################################################
      SUBROUTINE horizontal_oscillating_cylinder(numIB)
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      INTEGER, intent(in) :: numIB
      DOUBLE PRECISION:: PI
      DOUBLE PRECISION:: amp_z, frec_z, disz
      DOUBLE PRECISION:: amp_x, frec_x, disx
      INTEGER      :: L,I,K,strlen,strlen2,Geom_Time1,nplot    
      CHARACTER*8  :: char_block,num_ib
      CHARACTER*31 :: gridfile1     
	  !INTEGER:: gh_startTime = 0	  
	  PI = 4.D0*DATAN(1.D0)
   

      if(itime .gt. gh_startTime ) then        

      IF (myrank.eq.master) THEN

       Geom_Time1=301 

	nplot=1
		 
         if (mod(itime,n_out).eq.0) then !Write in file when outputs are writen
           write(char_block,'(I6)') itime
           strlen=LEN(TRIM(ADJUSTL(char_block)))
           char_block=REPEAT('0',(6-strlen))//TRIM(ADJUSTL(char_block))
           write(num_ib,'(I2)') numIB
           strlen2=LEN(TRIM(ADJUSTL(num_ib)))
           num_ib=REPEAT('0',(2-strlen2))//TRIM(ADJUSTL(num_ib))
           gridfile1='tecinscMotion'//TRIM(ADJUSTL(num_ib))//
     &      '_Time_'//TRIM(ADJUSTL(char_block))//'.plt'
           open (unit=Geom_Time1, file=gridfile1)
	   write(Geom_Time1,*)'title = points'
	   write(Geom_Time1,*)'variables=x,y,z'
           write(Geom_Time1,*) 'zone ','STRANDID=', 101,
     &    'SOLUTIONTIME=', ctime,
     &    ' i=',nodes(numIB),'zonetype=',
     & 	  'ordered',', DATAPACKING=point'

 
	 endif


	amp_x = 0.2d0
	frec_x = 2.0d0*PI*0.2 !0.156
	disx = amp_x*DSIN( frec_x*(CTIME) )	! ! vertical displacement in oscillating cylinder case


	! Updating node velocities

	U_p = frec_x*amp_x*DCOS(frec_x*(CTIME) ) !0.0d0
	W_p = 0.0d00 !frec_z*amp_z*DCOS(frec_z*(CTIME) )
	V_p = 0.0d0

	if(itime .lt. (gh_startTime+2) ) then
        do L = 1, nodes(numIB)
	!allocate(nodezo(1,10000))
	!nodezo(numIB,L) = nodez(numIB,L)
	nodexlocal(numIB,L) = nodex(numIB,L)	!nodezlocal(numIB,L) = nodez(numIB,L)	
	enddo
	endif

	! Updating node location
	do L = 1, nodes(numIB)
	!nodez(numIB,L) =  nodezlocal(numIB,L) + disz
	nodex(numIB,L) =  nodexlocal(numIB,L) + disx
	if (mod(itime,n_out).eq.0) then
		write(Geom_Time1,89) nodex(numIB,L),nodey(numIB,L),nodez(numIB,L)   
	endif						
	enddo

	if( mod(itime,n_out).eq.0)  close(Geom_Time1)



	if (myrank.eq.0) then
		 write(*,*) ctime, gh_startTime
		 write(*,*) U_p, W_p, disx !disz
		 write(*,*) Cx(numIB) + disx ,Cz(numIB)!Cx(numIB),Cz(numIB) + disz
		 write(*,*) nodezlocal(numIB,1), nodez(numIB,1)
	endif

	ENDIF ! End of master processor

	else ! The immersed body is has not yet started to move

	if (myrank.eq.0) then
		U_p = 0.0d0
		W_p = 0.0d0
		V_p = 0.0d0
	 	write(*,*) 104, W_p, nodez(numIB,1)
	endif ! End of master processor

	endif

   88 FORMAT (i5)
   89 FORMAT (3f25.5)
		RETURN
!#############################################################
        END SUBROUTINE
!#############################################################

!#############################################################
      SUBROUTINE imb_moved_shades(numIB)
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      INTEGER, intent(in) :: numIB
      DOUBLE PRECISION:: PI,xtop,ytop,ztop
      INTEGER      :: L,I,K,strlen,strlen2,GT1,numnodes,ntri,ii,P   
      CHARACTER*8  :: char_block,num_ib
      CHARACTER*31 :: gridfile1     
       PI = 4.D0*DATAN(1.D0)

      if (mod(itime,n_out).ne.0) RETURN

      IF (myrank.ne.master) goto 797

	IF (ibturbine(numIB).ne..true. .and. turax(numIB).ne.1) RETURN !Only works for VAT
	 	K=nodes(numIB)/imbnumber(numIB) 	! Vertical Axis Turbine
         GT1=301  
         write(char_block,'(I8)') itime
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(6-strlen))//TRIM(ADJUSTL(char_block))
           write(num_ib,'(I8)') numIB
           strlen2=LEN(TRIM(ADJUSTL(num_ib)))
           num_ib=REPEAT('0',(2-strlen2))//TRIM(ADJUSTL(num_ib))
           gridfile1='Blade_'//TRIM(ADJUSTL(num_ib))//
     &      '_Shade_'//TRIM(ADJUSTL(char_block))//'.plt'
         open (unit=GT1, file=gridfile1)
	   write(GT1,*)'TITLE = "Shaded blades"'
	   write(GT1,*)'VARIABLES="x","y","z"'
	   ntri=nscatter(numIB)*imbnumber(numIB)
	   write(GT1,*)'zone N=',(2*nscatter(numIB))*imbnumber(numIB),
     &    'E=',ntri,'DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL'
!Points writting
	Do P=1,imbnumber(numIB)
	  do L=1,nscatter(numIB)
	   ii=(P-1)*K+L
	   write(GT1,89)nodex(numIB,ii),nodey(numIB,ii),nodez(numIB,ii)   
	  enddo
	  do L=1,nscatter(numIB)
   	   ii=(P-1)*K+L
!		if (axis(numIB).eq.3) then
		xtop=nodex(numIB,ii)
		ytop=nodey(numIB,ii)
           	if (linfin(numIB).eq.1) ztop=zen
           	if (linfin(numIB).eq.0) ztop=zend(numIB)
!		endif
!		if (axis(numIB).eq.-1) then
!		xtop=nodex(numIB,ii)
!		ytop=nodey(numIB,ii)
!           	ztop=nodez(numIB,ii)
!		endif
	   write (GT1,89)xtop,ytop,ztop  
	  enddo
	Enddo !imbnumber
!Edges writting
	Do P=1,imbnumber(numIB)
	  do L=(2*(P-1))*nscatter(numIB)+1,(2*P-1)*nscatter(numIB)-1
	   write (GT1,88)L,L+1,nscatter(numIB)+L,nscatter(numIB)+L+1   
	  enddo
	   write (GT1,88)
     &     (2*P-1)*nscatter(numIB),(2*(P-1))*nscatter(numIB)+1,
     &	   (2*P)*nscatter(numIB),(2*P-1)*nscatter(numIB)+1  
	Enddo !imbnumber

	  close(GT1)

797	CONTINUE

   88 FORMAT (4i5)
   89 FORMAT (3e25.6)

        RETURN
        END SUBROUTINE	
