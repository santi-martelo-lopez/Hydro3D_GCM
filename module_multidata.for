!##########################################################################
        module multidata
!##########################################################################
	  SAVE
        integer :: nbp,nbpmax,num_domains
        integer :: rdivmax,idom,jdom,kdom
        integer,allocatable,dimension(:) :: rdiv
        integer,allocatable,dimension(:,:) :: rdv
        integer,allocatable,dimension(:) :: dom_id,dom_ad,dom_indid
        integer,allocatable,dimension(:) :: i_unst,j_unst,k_unst
        integer,allocatable,dimension(:) :: id_unst
	  integer,allocatable,dimension(:) :: imbinblk	!Pablo
        double precision, allocatable,dimension(:,:)::xcor,ycor,zcor !Pablo

        type multidom
           integer :: inext,iprev,jnext,jprev,knext,kprev
           integer :: corprev1,corprev2,corprev3,corprev4
           integer :: cornext1,cornext2,cornext3,cornext4
           integer :: edgprev1,edgprev2,edgprev3
           integer :: edgprev4,edgprev5,edgprev6
           integer :: edgnext1,edgnext2,edgnext3
           integer :: edgnext4,edgnext5,edgnext6
           integer :: per_in,per_ip,per_jn,per_jp,per_kn,per_kp
           integer :: ttc_i,ttc_j,ttc_k,ttc_ijk,ngrid
		!sv	 ttc_i : number of cells of the subdomain along the OX axis  
		! check: localparameters.for:140:           dom(ib)%ttc_i=nicell+2*pl


           integer :: isp,iep,jsp,jep,ksp,kep
	   integer :: isu,ieu,jsu,jeu,ksu,keu
	   integer :: isv,iev,jsv,jev,ksv,kev
	   integer :: isw,iew,jsw,jew,ksw,kew
           integer :: nwork,nvars
           integer :: mximb
           integer :: rq_m1,rq_p1,rq_m2,rq_p2,rq_m3,rq_p3
           integer :: rq_c1m,rq_c2m,rq_c3m,rq_c4m
           integer :: rq_c1p,rq_c2p,rq_c3p,rq_c4p
           integer :: rq_e1m,rq_e2m,rq_e3m,rq_e4m,rq_e5m,rq_e6m
           integer :: rq_e1p,rq_e2p,rq_e3p,rq_e4p,rq_e5p,rq_e6p
           double precision    :: xsl,ysl,zsl,xel,yel,zel,dx,dy,dz
		! distances relevant to gridsize

	     double precision,pointer,dimension(:) :: tauw
           double precision, pointer, dimension(:,:,:) :: S,So,Sm,Stm
           double precision, pointer, dimension(:,:,:) :: SUtm,SVtm,SWtm
	     double precision, pointer, dimension(:,:,:) :: sfactor
           double precision, pointer, dimension(:) :: x,y,z,xc,yc,zc
		! xc, yc, zc : coordinates of the cell center
           double precision, pointer, dimension(:,:,:) :: u,v,w,p,pp

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Aplication to ghost cell method
           !!! Copy of field variables from the latest time step
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: ucpy
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: vcpy
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: wcpy
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: pcpy
           !!! gradient of velocities needed to calculate shear stress
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: dudx
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: dudy
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: dudz
           !!!
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: dvdx
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: dvdy
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: dvdz
           !!!
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: dwdx
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: dwdy
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:):: dwdz
           !!! Flags for keeping track of the moving cells
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mucx ! moving u velocity cells for current time step
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mvcy ! moving v velocity cells for current time step
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mwcz ! moving w velocity cells for current time step
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mpc ! moving w velocity cells for current time step
           !!! Flags for keeping track of the emerging cells
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mucxold ! moving u velocity cells for former time step
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mvcyold ! moving v velocity cells for former time step
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mwczold ! moving w velocity cells for former time step
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mpcold ! moving w velocity cells for former time step
           !!
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mucxoldnn,mucxnn ! marker id 
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mucxoldib,mucxib ! domain id to which the marker belongs
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mucxoldll,mucxll ! layer id to which the marker belongs
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mucxoldi,mucxi !  dom(ib)%x id location of the fluid cells where the marker belongs
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mucxoldj,mucxj
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mucxoldk,mucxk
           !!
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mvcyoldnn,mvcynn ! marker id 
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mvcyoldib,mvcyib
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mvcyoldll,mvcyll
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mvcyoldi,mvcyi
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mvcyoldj,mvcyj
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mvcyoldk,mvcyk
           !!
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mwczoldnn,mwcznn ! marker id 
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mwczoldib,mwczib
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mwczoldll,mwczll
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mwczoldi,mwczi
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mwczoldj,mwczj
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mwczoldk,mwczk
           !!
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mpcoldnn,mpcnn ! marker id 
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mpcoldib,mpcib 
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mpcoldll,mpcll 
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mpcoldi,mpci
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mpcoldj,mpcj
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:):: mpcoldk,mpck 
           !!
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mucxoldx,mucxx ! fluid cell location where the marker belongs
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mucxoldy,mucxy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mucxoldz,mucxz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mucxoldmx
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mucxmx ! mirror fluid cell location associated to the marker
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mucxoldmy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mucxmy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mucxoldmz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mucxmz
       !!
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mvcyoldx,mvcyx
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mvcyoldy,mvcyy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mvcyoldz,mvcyz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mvcyoldmx
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mvcymx
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mvcyoldmy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mvcymy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mvcyoldmz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mvcymz
       !!
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mwczoldx,mwczx
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mwczoldy,mwczy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mwczoldz,mwczz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mwczoldmx
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mwczmx
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mwczoldmy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mwczmy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mwczoldmz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mwczmz
       !!
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mpcoldx,mpcx
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mpcoldy,mpcy  
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mpcoldz,mpcz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mpcoldmx,mpcmx
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mpcoldmy,mpcmy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: mpcoldmz,mpcmz
           !!                                                     
          INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mucxoldmi,mucxmi !  dom(ib)%x id location of the mirror fluid cell associated to the marker 
          INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mucxoldmj,mucxmj
          INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mucxoldmk,mucxmk
           !!                                                     
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mvcyoldmi,mvcymi
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mvcyoldmj,mvcymj
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mvcyoldmk,mvcymk
           !!                                                     
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mwczoldmi,mwczmi
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mwczoldmj,mwczmj
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mwczoldmk,mwczmk
           !!                                                     
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mpcoldmi,mpcmi
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mpcoldmj,mpcmj  
           INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: mpcoldmk,mpcmk
           !!
           !
           ! LU matrices for the tseng method, q4
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a3tsengu
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a3tsengv
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a3tsengw
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a3tsengp
           !
           ! permuation matrices for the tseng method, q4
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2tsengu
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2tsengv
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2tsengw
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2tsengp
           !
           ! pivot vectors for the tseng method, q4
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: p2tsengu
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: p2tsengv
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: p2tsengw
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: p2tsengp
           !
           !
           ! LU matrices for the bilinear least squared, q4
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a3blsu
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a3blsv
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a3blsw
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a3blsp
           !
           ! permuation matrices for the bilinear least squared, q4
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2blsu
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2blsv
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2blsw
         DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2blsp
           !
           ! pivot vectors for the bilinear least squared, q4
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: p2blsu
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: p2blsv
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: p2blsw
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: p2blsp
           !
           !
           ! LU matrices for the quadratic leas square method, q4
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: aqlsu
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: aqlsv
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: aqlsw
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: aqlsp
           !
           ! permuation matrices for the quadratic leas square method, q4
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2qlsu
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2qlsv
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2qlsw
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: a2qlsp
           !
           ! pivot vectors for the quadratic leas square method, q4
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: pqlsu
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: pqlsv
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: pqlsw
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:):: pqlsp
           !
           ! Deternimant of the core marix a3 in the quadratic least square method
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: da3u
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: da3v
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: da3w
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: da3p
           ! Transpose matrices for the quadratic least square method
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: aspu
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: aspv
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: aspw
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: aspp
           ! Inverse matrices for the quadratic least square method
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: ainvu
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: ainvv
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: ainvw
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:,:):: ainvp
           !
           ! q4 HP FEM shape functions
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::np1
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::npu
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::npv
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::npw
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::npp
           !
           ! q4 bilinear shape functions
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::nblu
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::nblv
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::nblw
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::nblp
           !
           ! q9 Lagrange shape functions
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::nlgu
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::nlgv
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::nlgw
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) ::nlgp
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! ! Aplication to IBM
           double precision, pointer, dimension(:,:,:) :: ksgs,ksgso
           double precision, pointer, dimension(:,:,:) :: eps,epso
           double precision, pointer, dimension(:,:,:) :: T,To,Tm,Ttm
           double precision, pointer, dimension(:,:,:) :: su,ap,sup
           double precision, pointer, dimension(:,:,:) :: ae,aw,as,an
           double precision, pointer, dimension(:,:,:) :: at,ab
           double precision, pointer, dimension(:,:,:) :: um,vm,wm
           double precision, pointer, dimension(:,:,:) :: pm,ppm
           double precision, pointer, dimension(:,:,:) :: vis,vism,epsm
           double precision, pointer, dimension(:,:,:) :: uum,vvm,wwm
           double precision, pointer, dimension(:,:,:) :: uvm,uwm,vwm
           double precision, pointer, dimension(:,:,:) :: ustar,vstar
           double precision, pointer, dimension(:,:,:) :: wstar
           double precision, pointer, dimension(:,:,:) :: uo,uoo,vo
           double precision, pointer, dimension(:,:,:) :: voo,wo,woo
           double precision, pointer, dimension(:,:,:) :: stfcinf
           double precision, pointer, dimension(:,:,:) :: facp1,facp2
           double precision, pointer, dimension(:,:,:) :: facm1,facm2
           double precision, pointer, dimension (:) :: dh1,dh2,dh3 
           integer, pointer, dimension (:,:,:,:) :: ndimb,imbbodynum
           integer, pointer, dimension (:) :: faz,cntp
           integer, pointer, dimension (:,:,:) :: ntav1,ntav2
           integer, dimension (26) :: tg
	     integer, pointer, dimension (:,:,:) :: ibfactor 
           double precision, pointer, dimension(:)     :: cof
           double precision, pointer, dimension(:,:)   :: tauwe,tauww
           double precision, pointer, dimension(:,:)   :: tauws,tauwn
           double precision, pointer, dimension(:,:)   :: tauwt,tauwb
           double precision, pointer, dimension(:,:)   :: tauwe2,tauww2
           double precision, pointer, dimension(:,:)   :: tauws2,tauwn2
           double precision, pointer, dimension(:,:)   :: tauwt2,tauwb2
           double precision, pointer, dimension(:) :: sendb_m1,sendb_p1
           double precision, pointer, dimension(:) :: recvb_m1,recvb_p1
           double precision, pointer, dimension(:) :: sendb_m2,sendb_p2
           double precision, pointer, dimension(:) :: recvb_m2,recvb_p2
           double precision, pointer, dimension(:) :: sendb_m3,sendb_p3
           double precision, pointer, dimension(:) :: recvb_m3,recvb_p3
           double precision, pointer, dimension(:)::sc1m,sc1p,rc1m,rc1p
           double precision, pointer, dimension(:)::sc2m,sc2p,rc2m,rc2p
           double precision, pointer, dimension(:)::sc3m,sc3p,rc3m,rc3p
           double precision, pointer, dimension(:)::sc4m,sc4p,rc4m,rc4p
           double precision, pointer, dimension(:)::se1m,se1p,re1m,re1p
           double precision, pointer, dimension(:)::se2m,se2p,re2m,re2p
           double precision, pointer, dimension(:)::se3m,se3p,re3m,re3p
           double precision, pointer, dimension(:)::se4m,se4p,re4m,re4p
           double precision, pointer, dimension(:)::se5m,se5p,re5m,re5p
           double precision, pointer, dimension(:)::se6m,se6p,re6m,re6p
           double precision, pointer, dimension(:,:,:) ::d1,dphi_dxplus
           double precision, pointer, dimension(:,:,:) ::dphi_dyplus,
     &dphi_dzplus,dphi_dxminus,dphi_dyminus,dphi_dzminus
!============================== LSM VARIABLES =============================
           double precision, pointer, dimension(:,:,:) :: phi_init,
     &phi_new,phi_reinit,phi,dphi_dx,dphi_dy,dphi_dz,s_phi0,h_phi,
     &dens,mu,phim,abs_dphi_check
		!  Reinit should be turned on (T) and this is to apply a method called re-initilization when calculating the free-surface to make sure that distance function phi maintain 
		! its property. This method uses the variables presented in a following question of yours (s_phi0,phi_reinit) you can go through some papers regarding  level set method 
		! which eventually will make things more clear
          !
          ! Implementation of the ghost cell method for LSM
          double precision, pointer, dimension(:,:,:):: phicpy,phimcpy
          double precision, pointer, dimension(:,:,:):: rhocpy,mucpy

           double precision, pointer, dimension(:)     :: dens_mg
           real, pointer, dimension(:) :: sendb_m,sendb_p
           real, pointer, dimension(:) :: recvb_m,recvb_p
           integer, pointer, dimension (:) :: ijkp_lsm
	     integer :: tot
		! tot: again an integer linked with the total cells, different meaning in different subroutines 
           integer :: niul,njul,nkul,nivl,njvl,nkvl,niwl,njwl,nkwl
           integer :: nipl,njpl,nkpl,nigl,njgl,nkgl
           integer :: nipl2,njpl2,nkpl2
           double precision, pointer, dimension(:,:,:) :: resmax
           double precision, pointer, dimension(:,:,:) :: resfact
           double precision, pointer, dimension(:,:,:) :: resfact1
!==========================================================================
           integer bc_west,bc_east,bc_south,bc_north,bc_bottom,bc_top 
           integer Tbc_west,Tbc_east,Tbc_south,Tbc_north
           integer Tbc_bottom,Tbc_top 
	   logical :: coarse_ng,fine_ng
        end type multidom

        type (multidom), pointer, dimension(:) :: dom

        end module multidata
!##########################################################################
