!######################################################################
      module imb
!######################################################################
	SAVE
!     DOUBLE PRECISION::xt(5),yt(5),xdt(5),ydt(5),xddt(5),yddt(5),yto
	!
	INTEGER:: interpolationScheme 	! selection of interpolation scheme
	INTEGER:: ncrs 					! number of corssection corners
	DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:):: cornersx,cornersz ! crossection corners coordinates
	!Forces of the ghost cell method
	DOUBLE PRECISION:: p1(3) !pressure p1(1) : pressure force in the OX axis direction
	DOUBLE PRECISION:: sh1(3) !shear stress sh1(1) : shear stress force in the OX axis direction
	DOUBLE PRECISION:: fv1(3) !volume force vf1(1) : volume force force in the OX axis direction due to the inmersed boundary forcing
	DOUBLE PRECISION:: tf1(3) !total force tf1(1) = pf(1) + sh1(1) : total force in the OX axis direction
	INTEGER:: forcefilegh ! file id for writing force values
	INTEGER:: deps ! Indicates the neighbourhood in which to look for cells when interpolating with delta functions
	INTEGER:: imax,jmax,kmax,imin,jmin,kmin

	! Variables for the case of Vertical oscillating cylinder
	DOUBLE PRECISION,allocatable,dimension(:) ::dis_z,dis_z_in, a_z, f_z
	DOUBLE PRECISION:: U_p, V_p, W_p ! Prescribed cartesian velocities
	INTEGER:: transverseFlag ! Flag for testing the bechmark case of a cylinder undergoing prescribe transverse oscillations

	! Variables for the case of Axial oscillating cylinder
	DOUBLE PRECISION,allocatable,dimension(:) ::dis_x,dis_x_in, a_x, fz_x	

	! copy of current velocity field
	INTEGER,allocatable,dimension(:,:,:)::ucpy, vcpy, wcpy

	! pi number definition
	DOUBLE PRECISION:: piNum = 4.D0*DATAN(1.D0)

	INTEGER :: imp_proc_master,imb_block_master,forcefilej
      INTEGER :: bodynum,maxnode,master,maxnodeIBS,mdfsteps,yangcase
	INTEGER :: ibm_out_forces,nIBslv,nl    
	
	! Serial Processing
	INTEGER,allocatable,dimension(:)::NumStation
	! NumSections: Number of sections created
	INTEGER,allocatable,dimension(:,:)::MrkrsExLy, NumLayers, MrkrsSection
	! MrkrsExLy: Number of markers in the most exterior layer (Ussually the first layer). Its dimensions are (Number of bodies, Number of crossections)
	! NumLayers: Number of layers created. Its dimensions are (Number of bodies, Number of crossections)
	! MrkrsSection: Number of markers per section. Its dimensions are (Number of bodies, Number of crossections)	
	INTEGER,allocatable,dimension(:,:,:)::MrkrsLy
	! MrkrsLy: Number of markers in each layer. Its dimensions are (Number of bodies, Number of crossections, Number of layers per section)

	! Parallel Processing
	INTEGER:: ibmSt_k	! counter of the total number of stations
	INTEGER,allocatable,dimension(:)::ibmSt !# Stations of each body, for each postion, it stores the 
												! number of stations defined in each body. Dimension: bodynum
	INTEGER,allocatable,dimension(:)::ibmStMkrs	!# IB points is each station of each body, for each postion jj, it stores the 
												! number of markers in station jj of body M. Dimension: imbSt_k
												! ibmStMkrs( (M-1)*ibmSt(M-1) + M*jj ) : number of markers at station jj for body M
	INTEGER,allocatable,dimension(:)::ibmMkrsEL	!# markers in the exterior layer in each station of each body
	! for each postion jj, it stores the number of markers in the exterior layer in station jj for body M
	! ibmMkrsEL( (M-1)*ibmSt(M-1) + M*jj ) : number of markers at exterior layer in station jj for body M
	! Dimension: imbSt_k




	LOGICAL,allocatable,dimension(:):: rotating,ibturbine
	
	DOUBLE PRECISION::lambda,sigma,dxm,dym,dzm,nxl
      DOUBLE PRECISION,allocatable::Cx(:),Cxor(:),Cy(:),Cyor(:)
      DOUBLE PRECISION,allocatable::Cz(:),Czor(:),l2norm(:)			
	DOUBLE PRECISION,allocatable,dimension(:)::xaero,yaero,zaero
	DOUBLE PRECISION,allocatable,dimension(:)::radsin,rads,pitch
	DOUBLE PRECISION,allocatable,dimension(:)::R,reddelta
      DOUBLE PRECISION,allocatable,dimension(:)::nodex_loc,nodey_loc
	DOUBLE PRECISION,allocatable,dimension(:)::nodez_loc
      DOUBLE PRECISION,allocatable,dimension(:)::FX1_loc,FX2_loc
      DOUBLE PRECISION,allocatable,dimension(:)::FX3_loc
      DOUBLE PRECISION,allocatable,dimension(:)::R0_loc,alpha0_loc
      DOUBLE PRECISION,allocatable,dimension(:)::U_Beta1_loc
	DOUBLE PRECISION,allocatable,dimension(:)::U_Beta2_loc
      DOUBLE PRECISION,allocatable,dimension(:)::U_Beta3_loc,zini,zend
      DOUBLE PRECISION,allocatable,dimension(:,:)::dh1_loc,dh2_loc
      DOUBLE PRECISION,allocatable,dimension(:,:)::dh3_loc
      DOUBLE PRECISION,allocatable,dimension(:,:)::nodexlocal
      DOUBLE PRECISION,allocatable,dimension(:,:)::nodeylocal
      DOUBLE PRECISION,allocatable,dimension(:,:)::nodezlocal
      DOUBLE PRECISION,allocatable,dimension(:)::FX1NF,FX2NF,FX3NF
      DOUBLE PRECISION,allocatable,dimension(:,:)::FX1,FX2,FX3
      DOUBLE PRECISION,allocatable,dimension(:,:)::FX1M,FX2M,FX3M
      DOUBLE PRECISION,allocatable,dimension(:,:)::nodex,nodey,nodez
      DOUBLE PRECISION,allocatable,dimension(:,:)::U_Beta1,U_Beta2
      DOUBLE PRECISION,allocatable,dimension(:,:)::U_Beta3
      DOUBLE PRECISION,allocatable,dimension(:,:)::alpha0,R0    
       
	INTEGER,allocatable,dimension(:,:) :: I_nr_V,J_nr_V,K_nr_V
	INTEGER,allocatable,dimension(:,:) :: I_nr_U,J_nr_U,K_nr_U
	INTEGER,allocatable,dimension(:,:) :: I_nr_W,J_nr_W,K_nr_W
	INTEGER,allocatable,dimension(:) :: imbnumber,imb_shape
	INTEGER,allocatable,dimension(:) :: kmaxU,kmaxV,kmaxW,nodes
	INTEGER,allocatable,dimension(:) :: imb_block,turax	
	INTEGER,allocatable,dimension(:) :: lag_bod_loc,cmax,linfin
      INTEGER,allocatable,dimension(:) :: imb_block_loc,axis
      INTEGER,allocatable,dimension(:) :: imbinblock_loc,rott_loc
	INTEGER,allocatable,dimension(:) :: nscatter	
               
      CHARACTER*32, allocatable, dimension (:) :: filepoints

	INTEGER,ALLOCATABLE,DIMENSION(:)::numIBslave,Lslave,Lslv,dsplc
	INTEGER,ALLOCATABLE,DIMENSION(:)::imb_proc_loc

	DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::xdh_U,ydh_U,zdh_U
	DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::xdh_V,ydh_V,zdh_V
	DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::xdh_W,ydh_W,zdh_W

	INTEGER,ALLOCATABLE,DIMENSION(:,:)::Irec_U,Jrec_U,Krec_U
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::Irec_V,Jrec_V,Krec_V
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::Irec_W,Jrec_W,Krec_W

	DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::USTAR_mas,VSTAR_mas
	DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::WSTAR_mas

!Self-starting 07_2017:
	LOGICAL,allocatable,dimension(:):: LSELFST
	DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::acc_ST,SUMtorque_ST	

!Actuator line 07_2017:
	DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::r_act,c_act,Pit_act
	INTEGER,ALLOCATABLE,DIMENSION(:) :: ppc_act
	DOUBLE PRECISION::Cl_act,Cd_act
	DOUBLE PRECISION,ALLOCATABLE :: F_X_MEAN(:,:),F_Q_MEAN(:,:)
	DOUBLE PRECISION,ALLOCATABLE :: F_Y_MEAN(:,:),F_Z_MEAN(:,:)

      end module imb
