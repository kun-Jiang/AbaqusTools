	!=======================================================
	!			 Define the global variables
	!=======================================================
      module global
 
      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.
 
      integer numElem,ElemOffset,err,NodeMax
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
	
      !real*8, allocatable :: globalSdv(:,:,:)
      real*8 globalDetFs(numElem,4)   ! to allocate the deformation
      real*8 globalCTau(numElem,4)    ! to allocate the concentration
      real*8 globalSdv(numElem,4,20)
      ! real*8 globalHydro(1000,numElem,4)   ! the hydrostatic pressure
      ! real*8 globalDetF(1000,numElem,4)    ! the total deformation
      integer test
 
	end module global

C***************************************************************************
	!!=======================================================
	!!   This solvTens contains the matrix multiplication
	!!=======================================================
      MODULE SolvTENS

      implicit none
      contains
      !----------------------------------------------------
       function Matmul422(A,B)
       integer i,j,k,L,M,N
       real*8 A(3,3,3,3),B(3,3),Matmul422(3,3)
       Matmul422 = 0.0
        do i=1,3
          do j=1,3
            do k=1,3
              do L=1,3
         Matmul422(i,j) = Matmul422(i,j) + A(i,j,k,L)*B(k,L)
              enddo
            enddo
          enddo  
        enddo
        end function Matmul422
      !----------------------------------------------------
       function Matmul220(A,B)
       integer i,j,k,L,M,N
       real*8 A(3,3),B(3,3),Matmul220
       Matmul220 = 0.0
        do i=1,3
          do j=1,3
          Matmul220 = Matmul220 + A(i,j)*B(i,j)
          enddo  
        enddo
       end function Matmul220
      !----------------------------------------------------
       function Matmul444(A,B)
       integer i,j,k,L,M,N
       real*8 A(3,3,3,3),B(3,3,3,3),Matmul444(3,3,3,3)        
        Matmul444 = 0.0
        do i=1,3
          do j=1,3
            do k=1,3
              do L=1,3
               do m=1,3
                do n=1,3
				Matmul444(i,j,k,L) = Matmul444(i,j,k,L) + A(i,j,m,n)*B(m,n,k,L)
                enddo
               enddo
              enddo
            enddo
          enddo  
        enddo
	 return
	 end function Matmul444
      !----------------------------------------------------
	 function Matmul244(A,B)
	 integer i,j,k,L,M,N
	 real*8 A(3,3),B(3,3,3,3),Matmul244(3,3,3,3)
        Matmul244 = 0.0
        do i=1,3
          do j=1,3
            do k=1,3
              do L=1,3
               do m=1,3
				Matmul244(i,j,k,L) = Matmul244(i,j,k,L) + A(i,m)*B(m,j,k,L)
               enddo
              enddo
            enddo
          enddo  
        enddo
	 return
       end function Matmul244
      !----------------------------------------------------
	 function Matmul242(A,B)
	 integer i,j,k,L,M,N
	 real*8 A(3,3),B(3,3,3,3),Matmul242(3,3)
        Matmul242 = 0.0
        do i=1,3
          do j=1,3
            do k=1,3
              do L=1,3
			   Matmul242(i,j) = Matmul242(i,j) + A(k,L)*B(k,L,i,j)
              enddo
            enddo
          enddo  
        enddo
	 return
	 end function Matmul242
      !----------------------------------------------------
	 function Matmul224A(A,B)
	 integer i,j,k,L,M,N        
	 real*8 A(3,3),B(3,3),Matmul224A(3,3,3,3)
        Matmul224A = 0.0
        do i=1,3
         do j=1,3
          do k=1,3
            do L=1,3      
			 Matmul224A(i,j,k,L)=Matmul224A(i,j,k,L)+A(i,k)*B(j,L)          
            enddo
           enddo
          enddo
	  enddo
	 return
	 end function Matmul224A 
      !----------------------------------------------------
	 function MatmulMV(A,B,N)
	 integer i,j,k,L,M,N        
	 real*8 A(N,N),B(N),MatmulMV(N)
        MatmulMV = 0.0
        do i=1,N
         do j=1,N
            MatmulMV(i)=MatmulMV(i)+A(i,j)*B(j)          
	   enddo
        enddo       
       return
       end function MatmulMV
      !----------------------------------------------------
       function Matmul224T(A,B)
       integer i,j,k,L,M,N        
       real*8 A(3,3),B(3,3),Matmul224T(3,3,3,3)
        Matmul224T = 0.0
        do i=1,3
         do j=1,3
           do k=1,3
             do L=1,3      
			  Matmul224T(i,j,k,L)= Matmul224T(i,j,k,L) +A(i,L)*B(j,k)          
             enddo
            enddo
          enddo
        enddo
       return
       end function Matmul224T  
      !----------------------------------------------------
       function Matmul224I(A,B)
       integer i,j,k,L
       real*8 A(3,3),B(3,3),Matmul224I(3,3,3,3)
        Matmul224I=0.0
        do i=1,3
         do j=1,3
           do k=1,3
             do L=1,3      
			  Matmul224I(i,j,k,L)= Matmul224I(i,j,k,L) +A(i,j)*B(k,L)       
             enddo
            enddo
          enddo
        enddo  
       return
       end function Matmul224I     
      end module SolvTENS
      
C****************************************************************************
      subroutine UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,U,DU,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)

      use global
      use SolvTENS
c
      IMPLICIT NONE
c
c     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
c
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
c
c     VARIABLES PASSED INTO UEL 
c
      REAL(8) :: PROPS,coords,U,DU,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
c
	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(*),ENERGY(8),
	1 PROPS(*),coords(MCRD,NNODE),U(NDOFEL),DU(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),
	2 TIME(2),PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
	3 PREDEF(2,NPREDF,NNODE),LFLAGS(*),DDLMAG(MDLOAD,*),JPROPS(*)

      integer lenJobName,lenOutDir
	character*256 jobName,outDir,fileName

      real*8 zero,one,two,Three
      parameter(zero=0.d0,one=1.d0,two=2.d0,Three=3.d0)
	
	! Open the debug/error message file
       call getJobName(jobName,lenJobName)
       call getOutDir(outDir,lenOutDir)
           fileName = outDir(1:lenOutDir) // '\Subroutine_' // jobName(1:lenJobName) // '.dat'
       open(unit=80,file=fileName,status='UNKNOWN') 
C	Call particular subroutine to calculate different layer
	 If (JTYPE .EQ. ONE) Then
C	Call UEL_Chem to calculate the chemical reaction and diffusion caused swell
		Call UEL_Chem(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,U,DU,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)
	 Else 
C	 Call UEL_Frac to calculate the fracture phase field
		Call UEL_Frac(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,Vel,Accn,JTYPE,
	2	 TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     3     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     4     NJPROP,PERIOD)
	 Endif

      return
	end subroutine UEL
	
      !////////////////////////////////////////////////////////////
      ! //        
      ! //              Chemical diffusion - reaction
      ! //        
      ! ////////////////////////////////////////////////////////////   
      subroutine UEL_Chem(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,U,DU,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)

      use global
      use SolvTENS

      IMPLICIT NONE

c     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS

      REAL(8) :: RHS,AMATRX,SVARS,ENERGY

C     VARIABLES PASSED INTO UEL 

      REAL(8) :: PROPS,coords,U,DU,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 thetaNew(nNode)
      real*8 thetaOld(nNode),dtheta(nNode),muNew(nNode)
      real*8 muOld(nNode),dMU(nNode),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,2)
      real*8 coordsC(mcrd,nNode),NCL
      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face,nInt,nIntS
      parameter(nDim=2,nInt=4,nIntS=1)

      integer ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,faceFlag,ndofn
	character*256 jobName,outDir,fileName

      real*8 Iden(3,3),theta0,phi0,Ru(2*nNode,1),Rc(nNode,1),body(3)
      real*8 Kuu(2*nNode,2*nNode),Kcc(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C,Vmol
      real*8 sh(nNode),detMapJ,phi_t,dsh(nNode,2),detMapJC,phiLmt,umeror
      real*8 dshC(nNode,2),mu_tau,mu_t,mutau,MUT(3),dMUdX(2,1),dMUdt
      real*8 SpTanMod(3,3,3,3),phi_tau,dPdt,DphiDmu,DphidotDmu,Mfluid
      real*8 Smat(3,1),Bmat(3,2*nNode),BodyForceRes(2*nNode,1),Qmat(4,4)
	real*8 sigHD,MuRef,CRate(3),JudgeReact,dCRdMu(3,3),Achem,ChemR,Mpfluid,
     +	DMPDMU,dMUPdx(3,1),DmDmu   
      real*8 xi(nInt,2),w(nInt),TrM0,detF,flux,wS(nIntS),ds,xLocal(nIntS),yLocal(nIntS)
      
      real*8 Nvec(1,nNode),ResFac,TanFac,detFs_t 

      real*8 D,C0,C_T,C_tau,MLi,dCdt,dCdMu,dCdotdMu,ResTemp,dt_M
      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)
	parameter(NCL=10.D0)
      !-----------------------------
      ! Do nothing if a ``dummy'' step
      !
      ! Check the procedure type, this should be a coupled
      !  temperature displacement or pore pressure displacement
      !  which are any of the following (64,65,72,73)
      !
      if((lflags(1).eq.64).or.(lflags(1).eq.65).or.
     +     (lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! all is good
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and chekc the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         !write(80,*) 'Abaqus does not have the right procedure'
         !write(80,*) 'go back and chekc the procedure type'
         !write(80,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem, I think this only matters
      !  when it comes to output in viewer
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         !write(80,*) 'Abaqus thinks you are doing'
         !write(80,*) 'a small displacement analysis'
         !write(80,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         !write(80,*) 'Abaqus thinks you are doing'
         !write(80,*) 'a linear purturbation step'
         call xit         
      endif
         
      if(dtime.eq.0.0) return

      ! Get element parameters
      !
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point
   
      ! Identity tensor
      !
      call onem(Iden,3)

      ! Obtain initial conditions
      !
      theta0   =  props(1)  
      Vmol     =  props(4)     
      D        =  props(6)
      C0       =  props(8)
      
      ! Initialize the residual and tangent matrices to zero
      Rc     = zero
      Kcc    = zero
      Energy = zero
      flux   = zero
   
      ! Body forces
      !
      body(1:3) = zero

      ! Obtain nodal displacements and chemical potentials
	
      k = 0
      do i=1,nNode
                 k =  k  + 1        
         MuNew(i)  =  U(k)
         dMu(i)    =  DU(k,1)
         MuOld(i)  =  MuNew(i)-dMu(i)  
      enddo
      
      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + U(j)
         enddo
      enddo      
      
       !
       !  chemical potential increment of polymer
       !
           do i=1,nNode
              if(dabs(dMu(i)).gt.1.d6)then
				WRITE(80,*) 'Increment half',427
               pnewdt = 0.5
               return
              endif
           enddo 
    
      !----------------------------------------------------------------
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         !write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif
          
      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
		   WRITE(80,*) 'increment half',451
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
		   WRITE(80,*) 'increment half',458
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         !write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
       endif
    

      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
		   WRITE(80,*) 'increment half',475
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
		   WRITE(80,*) 'increment half',482
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         !write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif         
         
         
      !----------------------------------------------------------------
      ! Begin the loop over body integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            !write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         !write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif
   
      ! Loop over integration points
      jj = 0   ! jj is used for tracking the state variables
      do intpt=1,nIntPt
         ! Obatin state variables from the previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then     
            C_T  =  C0          
         else   
            C_T  =  svars(1+jj)  ! concentration of the diffusant
         endif
         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            !write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
   
         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
		   WRITE(80,*) 'increment half',543
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
		   WRITE(80,*) 'increment half',550
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            !write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif

         ! Obtain the chemical potential and its derivative's at 
         !  this intPt at the begining and end of the incrment            
         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
		   WRITE(80,*) 'increment half',568
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
		   WRITE(80,*) 'increment half',575
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            !write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
          endif
            
         ! Obtain the chemical potential and its derivative's at 
         !  this intPt at the begining and end of the incrment            
             Mu_tau = zero
             Mu_t   = zero
             dMudt  = zero
             dMudX  = zero
   
         do k = 1,nNode
             Mu_tau = Mu_tau + MuNew(k)*sh(k)
             Mu_t   = Mu_t   + MuOld(k)*sh(k)

             if(Mu_tau.ge.0.0) then
                 Mu_tau = -1.e-6
             endif           
             
             do i=1,nDim
                 dMudX(i,1) = dMudX(i,1) + MuNew(k)*dshC(k,i)
             enddo
         enddo
            dMudt = (Mu_tau - Mu_t)/dtime   

         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point 
         !
         !===================================
         !  mass diffusion   
         !TrM0  = gloTrSig(jelem,intpt,1) * 1000.0   ! this is the hydro-static stress
         TrM0    = zero
         detF    = one
         sigHD   = zero
         detFs_t = zero
         !
         !call SolvFracInfTot(C_Tau,dCdotdMu,dCdt,props,nprops,Mu_Tau,MuRef,C_T,detF,sigHD,dtime,jelem,detFs_t)
         
         call SolvFracInfTot(C_Tau,dCdotdMu,dCdt,props,nprops,Mu_Tau,MuRef,C_T,detF,sigHD,dtime,jelem,detFs_t)

         
         if(stat.eq.0) then
		   WRITE(80,*) 'increment half',625
            pnewdt = 0.5
            return
         endif
  
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         svars(1+jj) = C_tau
         jj = jj + nlSdv ! setup for the next intPt

         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         
         if (kinc.ne.0) then        
            globalCTau(jelem-ElemOffset,intpt) = C_tau
            globalDetFs(jelem-ElemOffset,intpt) = one  + Vmol * C_tau
           ! globalSdv(jelem,intPt,1)   = C_tau   ! polymer volume fraction
           ! globalSdv(jelem,intpt,2)   = C0
           ! globalSdv_t(jelem,intPt,1) = C_t
         end if

	!=========================================================     
	!       chemical reactoin,    
	!  Here to avoid the inconvergence of the chemical reaction 
	!  where too small or negative concentrations
	!  When the concentrations are sufficiently small, 
	!  we neglect the chemical reaction
      !  the chemical potential for current step  
      ! to avoid computational convergency
              
	!call Umat_ChemReact(JudgeReact,props,nprops,C_T,MuT,CRate,dCRdMu,Achem,ChemR,Ncl,time,dtime,Jelem)
      
      ! (1) Compute/update the chemical potential residual vector
         
         do kk=1,nNode
             Nvec(1,kk)=sh(kk)
         enddo
         
      Mfluid  = D   ! here the diffusion coefficient
         
      ResFac  = zero
      ResFac  = -dCdt
      Rc = Rc + detmapJC*w(intpt)*
     1  (
     2    transpose(Nvec)*ResFac - Mfluid*matmul(dshC,dMudX)
      !3 + transpose(Nvec) * CRate(3)
     4   ) 
      

C============================================================
C           update the tangent matrix         
C============================================================ 

      !*********************
      !  (1) Kcc          Kcc
      !*********************
      !Compute/update the chemical potential tangent matrix
       TanFac = zero
       TanFac = DCdotDmu
       DmDmu  = zero   
            
       KCC  =  KCC  + detmapJC*w(intPt)*
     1  (
     2   TanFac * matmul(transpose(Nvec),Nvec) + Mfluid*matmul(dshC,transpose(dsh))
     4   + DmDmu*matmul(matmul(dsh,dMUdX),Nvec)
      !5   -dCRDmu(2,2)*matmul(transpose(Nvec),Nvec)     
     6   )
        
       enddo
           

      !
      ! End the loop over body integration points
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      ! Start loop over surface fluid flux terms here
      !
      !
      if(ndload.gt.0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i=1,ndload
            !
            ! based on my convention the face which the flux/traction
            !  acts on is the flux/traction ``label''
            !
            face = jdltyp(i,1) ! label
            flux = adlmag(i,1) ! flux magnitude
            
            if((face.ge.1).and.(face.le.4)) then
               !
               ! fluid flux applied
               !
               if(face.eq.1) then
                  faceFlag = 1
               elseif(face.eq.2) then
                  faceFlag = 2
               elseif(face.eq.3) then
                  faceFlag = 3
               elseif(face.eq.4) then
                  faceFlag = 4
               endif
               !
               if(nIntS.eq.1) then
                  call xintSurf2D1pt(faceFlag,xLocal,yLocal,wS)
               elseif(nIntS.eq.2) then
                  call xintSurf2D2pt(faceFlag,xLocal,yLocal,wS)
               elseif(nIntS.eq.3) then
                  call xintSurf2D3pt(faceFlag,xLocal,yLocal,wS)
               else
                  write(*,*) 'Invalid nIntS points, nIntS=',nIntS
                  !write(80,*) 'Invalid nIntS points, nIntS=',nIntS
                  call xit
               endif
               !
               ! loop over integ points
               !
               do ii=1,nIntS
                  !
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),faceFlag,coordsC,sh,ds)
                  !
                  ! Modify the residual, loop over nodes, recall
                  !  sh(n)=0 when n is not on this face
                  !
                  do n=1,nNode
                     ! Rc(n,1) = Rc(n,1) - wS(ii)*ds*sh(n)*flux
                  enddo
                  !
                  ! No change to the tangent matrix
                  !
               enddo ! loop over nIntS
               !
            else
               write(*,*) 'Unknown face=',face
               !write(80,*) 'Unknown face=',face
               call xit
            endif

         enddo ! loop over ndload
         endif ! ndload.gt.0 or not
      !
      ! End loop over flux and traction terms
      !----------------------------------------------------------------
      ! output the solved variables into the global variables
  
   
      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.  This
      !  is essentially giving Abaqus the residual and the tangent matrix.
      !
      ! Return Abaqus the right hand side vector
      !
      ! the sequence of the degree is U1,U2,U3,Phi,T,C 
         Rhs(:,1) = zero
         Amatrx   = zero
         nDofN    = nDofEl/nNode
         
      do i=1,nNode
         A11 = nDofN*(i-1)+1
         A12 = nDim*(i-1)+1
         !
         ! Chemical potential of polymer
         !
         rhs(A11,1)= RC(i,1)   
      enddo
      
      !
      ! Return Abaqus the tangent matrix
      !
      do i=1,nNode
         do j=1,nNode
            A11  = nDofN*(i-1)+1
            A12  = nDim*(i-1)+1
            B11  = nDofN*(j-1)+1
            B12  = nDim*(j-1)+1
         ! chemical mass - chemical mass
            amatrx(A11,B11) = KCC(i,j)  
         enddo
      enddo
      !
      return
      end subroutine UEL_Chem
C****************************************************************************
         
      !////////////////////////////////////////////////////////////
      ! //        
      ! //               phase fracture
      ! //        
      ! ////////////////////////////////////////////////////////////        
************************************************************************
      SUBROUTINE UEL_Frac(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     ==================================================================
      INCLUDE 'ABA_PARAM.INC'
C     ==================================================================
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,THREE=3.D0,
     1 TOLER=1.0D-8,FOUR=4.D0,RP25 = 0.25D0,HALF=0.5D0,SIX=6.D0,
     2 N_ELEM=1,NSTVTO=2,NSTVTT=14,NSTV=18)
C     ==================================================================
C     Initialization for all the element types
C     ==================================================================
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
     
       INTEGER I,J,L,K,K1,K2,K3,K4,IX,IY

       REAL*8 AINTW(NNODE),XII(NNODE,2),XI(2),dNdxi(NNODE,2),
     1 VJACOB(2,2),dNdx(NNODE,2),VJABOBINV(2,2),AN(4),BP(2,NDOFEL),
     2 DP(2),SDV(NSTV),BB(3,NDOFEL),CMAT(3,3),EPS(3),STRESS(3),
     3 VNI(2,NDOFEL),ULOC(2),PHASENOD(NNODE)
     
       REAL*8 DTM,THCK,HIST,CLPAR,GCPAR,EMOD,ENU,PARK,ENG
 
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,4)
      integer lenJobName,lenOutDir
	character*256 jobName,outDir,fileName	
	! Open the debug/error message file
       call getJobName(jobName,lenJobName)
       call getOutDir(outDir,lenOutDir)
           fileName = outDir(1:lenOutDir) // '\Subroutine_' // jobName(1:lenJobName) // '.dat'
       open(unit=80,file=fileName,status='UNKNOWN') 
C
C     ==================================================================
C     ******************************************************************
C     Constructing elemet TYPE 1 (phase)
C     ******************************************************************
C     ==================================================================
C      
       IF (JTYPE.EQ.Two) THEN
C     ==================================================================
C     Time an iteration variables
C     ==================================================================
       TIMEZ=USRVAR(JELEM-Two*N_ELEM,17,1)
       IF (TIMEZ.LT.TIME(2)) THEN
        USRVAR(JELEM-Two*N_ELEM,17,1)=TIME(2)
        USRVAR(JELEM-Two*N_ELEM,18,1)=ZERO
       ELSE
        USRVAR(JELEM-Two*N_ELEM,18,1)=USRVAR(JELEM-Two*N_ELEM,18,1)+ONE
       ENDIF
       STEPITER=USRVAR(JELEM-Two*N_ELEM,18,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
       CLPAR=PROPS(1)
       GCPAR =PROPS(2)
       THCK = PROPS(3)
C     ==================================================================
C     Initial preparations
C     ==================================================================
       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
        END DO
       END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
       XII(1,1) = -ONE/THREE**HALF
       XII(1,2) = -ONE/THREE**HALF
       XII(2,1) = ONE/THREE**HALF
       XII(2,2) = -ONE/THREE**HALF
       XII(3,1) = ONE/THREE**HALF
       XII(3,2) = ONE/THREE**HALF
       XII(4,1) = -ONE/THREE**HALF
       XII(4,2) = ONE/THREE**HALF
       DO I=1,NNODE
        AINTW(I) = ONE
       END DO
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
       DO INPT=1,NNODE
C     Initializing solution dependent variables (phase,history)
        DO I=1,NSTVTO
          SDV(I)=SVARS(NSTVTO*(INPT-1)+I)
        END DO
C
C     Local coordinates of the integration point
        XI(1) = XII(INPT,1)
        XI(2) = XII(INPT,2) 
C     Shape functions and local derivatives
        CALL SHAPEFUN(AN,dNdxi,XI)
C     Jacobian
        DO I = 1,2
         DO J = 1,2
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
          END DO
         END DO
        END DO
C       Determinant of Jacobian matrix
        DTM = ZERO
        DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
        IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT	
        END IF
C     Inverse of Jacobian
        VJABOBINV(1,1)=VJACOB(2,2)/DTM
        VJABOBINV(1,2)=-VJACOB(1,2)/DTM
        VJABOBINV(2,1)=-VJACOB(2,1)/DTM
        VJABOBINV(2,2)=VJACOB(1,1)/DTM
C        
C     Derivatives of shape functions respect to global ccordinates
        DO K = 1,NNODE
         DO I = 1,2
          dNdx(K,I) = ZERO
          DO J = 1,2
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
          END DO
         END DO
        END DO
C
C     Calculating B matrix (B=LN)
       DO INODE=1,NNODE
        BP(1,INODE)=dNdx(INODE,1)
        BP(2,INODE)=dNdx(INODE,2)
       END DO
C
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
        PHASE=ZERO
        DPHASE=ZERO
        DO I=1,4
         PHASE=PHASE+AN(I)*U(I)
        END DO
        DO I=1,4
         DPHASE=DPHASE+AN(I)*DU(I,1)
        END DO
C        
        IF (STEPITER.EQ.ZERO) THEN
          SDV(1)=PHASE-DPHASE
        ELSE
          SDV(1)=PHASE
        ENDIF
C
C     Gradient
        DO I=1,2
         DP(I)=ZERO
        END DO
        DO I=1,2
         DO J=1,4
          DP(I)=DP(I)+BP(I,J)*U(J)
         END DO
        END DO
C
C     ==================================================================
C     Calculating elastic ENERGY history
C     ==================================================================
	IF (STEPITER.EQ.ZERO) THEN
         ENGN=USRVAR(JELEM-Two*N_ELEM,13,INPT)
	ELSE
         ENGN=USRVAR(JELEM-Two*N_ELEM,16,INPT)
	ENDIF
C        
	HISTN=USRVAR(JELEM-Two*N_ELEM,16,INPT)
	IF (ENGN.GT.HISTN) THEN
         HIST=ENGN
	ELSE
         HIST=HISTN
	ENDIF
	SDV(2)=HIST
	
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
        DO I=1,NNODE
         DO K=1,NNODE
          DO J=1,2
           AMATRX(I,K)=AMATRX(I,K)+BP(J,I)*BP(J,K)*DTM*
     1      THCK*GCPAR*CLPAR*AINTW(INPT)
          END DO
          AMATRX(I,K)=AMATRX(I,K)+AN(I)*AN(K)*DTM*THCK*
     1     AINTW(INPT)*(GCPAR/CLPAR+TWO*HIST)
         END DO
        END DO
C        
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
        DO I=1,NDOFEL
         DO J=1,2
           RHS(I,1)=RHS(I,1)-BP(J,I)*DP(J)*GCPAR*CLPAR*
     1      AINTW(INPT)*DTM*THCK
         END DO
         RHS(I,1)=RHS(I,1)-AN(I)*AINTW(INPT)*DTM*THCK*
     1    ((GCPAR/CLPAR+TWO*HIST)*PHASE-TWO*HIST)
	END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C	NSTVTO=2,NSTVTT=14,NSTV=18
C     ==================================================================
        DO I=1,NSTVTO
         SVARS(NSTVTO*(INPT-1)+I)=SDV(I)
         USRVAR(JELEM-Two*N_ELEM,I+NSTVTT,INPT)=SVARS(NSTVTO*(INPT-1)+I)
        END DO
       END DO

      ENDIF
C      
      RETURN
	END
****************************************************************************

      !////////////////////////////////////////////////////////////
      ! //        
      ! //              Mechanical.
      ! //        Here we use UMAT to solve the mechanical 
      ! ////////////////////////////////////////////////////////////  
         
C========================================================
C           Umat  Mechanical 
C    the mechanical properites are solved by the umat
C========================================================
        subroutine UMAT(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      
      use global
      use SolvTENS
      
      include 'aba_param.inc'

      dimension stress(ntens),statev(nstatv),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      character*80 cmname,file1
      character*256 jobName,outDir,fileName

      integer i,j,k,l,iterError,lenJobName,lenOutDir,NELEMAN
	INTEGER NNPT

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),theta_t,theta_tau,T_tau(3,3)
      real*8 effStr,detF,dTRdF(3,3,3,3),SpTanMod(3,3,3,3),trB
      real*8 Finv(3,3),Bdis(3,3),trBdis,Bdis0(3,3),B_tau(3,3),B0(3,3)
      real*8 KS_TAU(3,3),dBetdF(3,3),dBetIdF(3,3),dInvLangdF(3,3),
     1 dBdis0dF(3,3,3,3),dKS0dF(3,3,3,3),dKS00dF(3,3,3,3),dKSdF(3,3,3,3)
      real*8 temp,dtemp,dXSRdF(3,3),dXSRIdF(3,3),alpha,theta_0
      
      real*8 SigStrs(3,3),KirStrs(3,3),KirStrV(6),Fe_tau(3,3),detFe

      real*8 Ndir0(3),Ndir_tau(3),detFg,SigHD
      
      real*8 LamdG_t,LamdG_tau,FG_t(3,3),FG_tau(3,3),SigStrV(6)
      real*8 Ftinv(3,3),DF_tau(3,3),DFe_tau(3,3),FgInv(3,3),SpUT(3,3)
      real*8 DFG_tau(3,3),DFgInv(3,3),DFeInv(3,3),Fe_t(3,3),Lamd
       
	REAL*8 detFs,Fs_inv(3,3),Fep_tau(3,3)
	
      real*8 MatTan(3,3,3,3),S_t,S_tau,Fp_t(3,3),Fp_tau(3,3),gBarP_t
      real*8 gBarP_tau,nuP_t,nuP_tau,gBarLmt,plasticwork
      real*8 StressTempJac(3,3),PlWrkDefJac(3,3)
      real*8 tao_tau,DESDT(6,6)
           
      real*8 Teq_tau(3,3),SpUUEQMod(3,3,3,3)
      real*8 Tne_tau(3,3),SpUUNeMod(3,3,3,3)
	
	REAL*8 RSlip(3,3,12),gBarP_t_Temp(12),gBarP_tau_Temp(12)
	REAL*8 PHASE,ENGE,STRAIN(3,3),T_tauD(3),T_tauV(3,3),STRAIN_D(3),STRAIN_V(3,3)
	REAL*8 STRS(3,3),STRS_D(3),STRS_V(3,3)
	REAL*8 Eyoung,poisson,Gshear,Kbulk,Lambda,Ramp_B
      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0,PARK = 1e-10)
	
      PARAMETER(N_ELEM=1,NSTV=18)

	 COMMON/KUSER/USRVAR(N_ELEM,NSTV,4)
       call onem(Iden,3)
	!Read material properties
       Eyoung   =  props(1)
       poisson  =  props(2)
C	Define two-order Identity matrix
	call onem(Iden,3)
	 
C	Obtain old and new deformation gradients
      F_t   = dfgrd0
      F_tau = dfgrd1	
C	call subroutine to bring the swelling deformation
	call mdet(F_tau,detF)
	!call Sub_Ch_Swell(kinc,noel,npt,DetFs)
	DetFs = One
C ***********************************************************************
C					solving the total stress
C ***********************************************************************
	call GetStrs(F_t,F_tau,props,nprops,SigStrs,KirStrs,
     1	SpTanMod,DetFs,theta_t,theta_tau)
	!T_tau = SigStrs           ! Cauchy stresses
	SpUT  = zero
	T_tau =  SigStrs
C	
	NELEMAN=NOEL
	NNPT = NPT
	IF (NNPT.EQ.3) THEN
		NNPT=4
	ELSEIF (NNPT.EQ.4) THEN
		NNPT=3
	ENDIF
      call SolvEnergy(F_tau,T_tau,ENGE)
	!ENGE = 100.d0
	IF(ENGE .GT. ZERO) THEN
		USRVAR(NELEMAN,13,NNPT) = ENGE
	ELSE
		USRVAR(NELEMAN,13,NNPT) = ZERO
	END IF
	PHASE = USRVAR(NELEMAN,15,NNPT)
	WRITE(80,*)'ENGE',ENGE
C	Stress and stiffness degration 
	!SpTanMod = SpTanMod*((ONE-PHASE)**TWO+PARK)
C 	Updating Phase field paramenter        
	DO I=1,18
		STATEV(I)=USRVAR(NELEMAN,I,NNPT)
	END DO
C	Return Abaqus/Standard the Cauchy stress
	if(ntens.eq.6) then
C
C	3D problem
C
		stress(1) = T_tau(1,1)
		stress(2) = T_tau(2,2)
		stress(3) = T_tau(3,3)
		stress(4) = T_tau(1,2)
		stress(5) = T_tau(1,3)
		stress(6) = T_tau(2,3)
         
	elseif(ntens.eq.4) then
C
C	2D problem
C   
		stress(1) = T_tau(1,1)
		stress(2) = T_tau(2,2)
		stress(3) = T_tau(3,3)
		stress(4) = T_tau(1,2)       
	endif

C ***********************************************************************
C		Return Abaqus/Standard the stress-deformation jacobian
C ***********************************************************************

      if(ntens.eq.6) then
         call jac3D(SpTanMod,ddsdde)
         call Jac31D(SpUT,ddsddt)
      elseif(ntens.eq.4) then
         call jac2D(SpTanMod,ddsdde)
         call Jac21D(SpUT,ddsddt) 
         !ddsddt = zero
	endif
	WRITE(80,*)'C22',DDSDDE(2,2)
	STATEV(19) = DDSDDE(2,2)
C
C	Return Abaqus the volumetric heat generated caused by plasticity
C
	rpl      =   0.0
	drpldt     =  0.d0

      return
      end subroutine umat

C************************************************************************
C						Material subroutine
C************************************************************************
      subroutine Umat_gent(props,nprops,F_tau,T_tau,
     +        SpUUMod,stat)

      implicit none
      !
      integer i,j,k,l,m,n,nprops,stat
      !
      real*8 props(nprops),F_tau(3,3),ER(3,1),T_tau(3,3),D_tau(3,1),
     +  SpUUMod(3,3,3,3),SpUPhiMod(3,3,3),SpPhiUMod(3,3,3),
     +  SpPhiPhiMod(3,3),Iden(3,3),Gshear,Kbulk,Imax,permit,detF,
     +  Finv(3,3),FT(3,3),FinvT(3,3),C_tau(3,3),Cinv(3,3),detC,trC,
     +  I1bar,fac,GShearGent,E(3,1),vec(3,1),I6,TR_tau(3,3),dGdF(3,3),
     +  dTRdF(3,3,3,3)
      !
      real*8 zero,one,two,three,fourth,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0)
 

      ! Identity matrix
      !
      call onem(Iden,3)
 

      ! Obtain relevant material properties
      !
      Gshear  = props(17) ! Shear modulus
      Kbulk   = props(18) ! Bulk modulus
      Imax    = props(19) ! Max I1bar
      permit  = props(20) ! Permittivity


      ! Compute the determinant of the deformation gradient
      !
      call mdet(F_tau,detF)

      ! Compute the determinant, the inverse, the transpose, 
      !  and the inverse transpose of the deformation gradient
      !
      call matInv3D(F_tau,Finv,detF,stat)
      FT     = transpose(F_tau)
      FinvT  = transpose(Finv)

      ! Compute the right Cauchy-Green tensor, its inverse, 
      !  its determinant, and its trace
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      call matInv3D(C_tau,Cinv,detC,stat)
      trC = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)
 
      ! Compute the trace of the distortional right Cauchy-Green tensor
      !
      I1bar = (detF**(-two/three))*trC
 
      ! Compute the ``I-3'' factor appearing in the Gent model
      !
      fac  = (I1bar - three)/Imax
      if(fac.gt.0.95d0) fac = 0.95d0
      fac = one/(one - fac)
 
      ! Compute the ``shear'' modulus. Note: this is not really the shear
      !  modulus, but it will help when computing the material tangent later
      !
      GShearGent = Gshear*fac
 
      ! Compute the derivative of the ``shear modulus'' with respect
      !  to the deformation gradient for use in the material tangent
      !
      dGdF   = two*(Gshear/Imax)*(detF**(-two/three))*
     +     fac*fac*(F_tau - third*trC*FinvT)
      
 

      ! Compute the 1st Piola stress
      !
      TR_tau = (detF**(-two/three))*GShearGent*(F_tau-third*trC*FinvT)
     +     + Kbulk*detF*(detF - one)*FinvT
 

      ! Compute the Cauchy stress
      !
      T_tau = (one/detF)*matmul(TR_tau,transpose(F_tau))


      ! Calculate the material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + (detF**(-two/three))*dGdF(k,l)*
     +                 (
     +                 F_tau(i,j) - third*trC*Finv(j,i)
     +                 )
     +                 + (detF**(-two/three))*GshearGent*
     +                 (
     +                 (-two/three)*F_tau(i,j)*Finv(l,k)
     +                 + (two/9.d0)*trC*Finv(j,i)*Finv(l,k)
     +                 + Iden(i,k)*Iden(j,l)
     +                 + third*trC*Finv(l,i)*Finv(j,k)
     +                 - (two/three)*Finv(j,i)*F_tau(k,l)
     +                 )
     +                 + detF*Kbulk*
     +                 (
     +                 (detF-one)*Finv(j,i)*Finv(l,k)
     +                 + detF*Finv(j,i)*Finv(l,k)
     +                 - (detF-one)*Finv(l,i)*Finv(j,k)
     +                 )
               enddo
            enddo
         enddo
      enddo
      
      
      ! Calculate the spatial tangent modulus
      !
      SpUUMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpUUMod(i,j,k,l) = SpUUMod(i,j,k,l) +
     +                      (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
C--------------------------------      
      return
	end subroutine Umat_gent


C****************************************************************************        
	SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
	1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
	2 JMAC,JMATYP,MATLAYO,LACCFLA)

    ! This subroutine is used to transfer SDV's from the UEL
    !  onto the dummy mesh for viewing.  Note that an offset of
    !  ElemOffset is used between the real mesh and the dummy mesh.
    !  If your model has more than ElemOffset UEL elements, then
    !  this will need to be modified.
   
	use global
   
	include 'ABA_PARAM.INC'

	CHARACTER*80 CMNAME,ORNAME
	CHARACTER*3 FLGRAY(15)
	DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
	DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
	PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,THREE=3.D0,
	1 TOLER=1.0D-8,FOUR=4.D0,RP25 = 0.25D0,HALF=0.5D0,SIX=6.D0,
	2 N_ELEM=1,NSTVTO=2,NSTVTT=14,NSTV=18)
C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.
	!if (kinc.eq.0) then
	!globalSdv(noel-ElemOffset*2,npt,1)=0.0
	!end if
C	Dummy mesh
	!IF (NPT .EQ. Three) NPT = Four
	!IF (NPT .EQ. Four)  NPT = Three
C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.
C	Concentration
	uvar(1) = globalCTau(noel-Four*N_ELEM,npt)
C	Swelling 
	uvar(2) = globalDetFs(noel-FOUR*N_ELEM,npt)
C	Phase variable
	!uvar(3) = globalSdv(noel-Three*N_ELEM,npt,15)  !dummy mesh
    
C      Do i=1,NUVARM
C       uvar(i) = globalSdv(noel-2*ElemOffset,npt,i)
C      enddo
    
	return
	end subroutine uvarm
    
C**********************************************************************
	SUBROUTINE MPROD(A,B,C)
 
C 	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(3,3),B(3,3),C(3,3)

	DO 2 I = 1, 3
	  DO 2 J = 1, 3
	    C(I,J) = 0.D0
	    DO 1 K = 1, 3
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE
C
	RETURN
	END

C**********************************************************************
	SUBROUTINE MPROD4(A,B,C)
 
C	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(4,4),B(4,4),C(4,4)

	DO 2 I = 1, 4
   	  DO 2 J = 1, 4
	    C(I,J) = 0.D0
	    DO 1 K = 1, 4
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE

	RETURN
	END
	
C**********************************************************************
	SUBROUTINE ZEROM(A)
C
C	THIS SUBROUTINE SETS ALL ENTRIES OF A 3 BY 3 MATRIX TO 0.D0.
C**********************************************************************

	REAL*8 A(3,3)

	DO 1 I=1,3
	  DO 1 J=1,3
	    A(I,J) = 0.D0
1	CONTINUE
C	
	RETURN
	END
	
c**********************************************************
      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(80,*) 'WARNING: SUBROUTINE matInv3:'
        write(80,*) 'WARNING: DET of MAT=',DET_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D
	
C///////////////////////////////////////////
C                  utiities
C///////////////////////////////////////////


      
C****************************************************************************
      subroutine SolvFracInfTot(C_Tau,dCdotdMu,dCdt,MatProp,nMatProp,MuTau,MuRef,C_T,
     1 detF,sigHD,dtime,noel,detFs_t)
      use global
      use SolvTENS
      implicit none
        
      integer I,J,K,L,M,N,nMatProp,IMax,KIter,Stat,noel,nMDiffs
      parameter(nMDiffs=10)
      real*8 Mdiffs(nMDiffs)
      
      real*8 C_PER,C_Tau,C_M,MatProp(nMatProp),MuTau,MuRef,C_t
      real*8 FerrS,JacbK,JacbInv,dPhi,detF,sigHD
      real*8 one,two,zero,three,half,ErrorC
      real*8 C_Temp,dPhiError,dtime,PerMu,detFs_t
      real*8 dCdotdMu,deltaMu,dCdt,dCdt_per,dCdt_m
      real*8 Mu_Per,Mu_M,Phi_Per,phi_M
      real*8 Phi0,DPDT_PER,DPDT_M
      real*8 theta0,theta_tau,Kbulk,chi,mu0,Vmol,Rgas
 
      parameter(IMax=50,zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5)
      parameter(PerMu = 1.e-8)
      
      ! ===================
      ! here we found using the following initial values can help to solving the phi
      ! obtain the material properties
       Kbulk = 1.e8
       chi = zero
       mu0 = MatProp(5)
       Vmol = MatProp(4)
       Rgas = MatProp(2)
       theta0 = MatProp(1)
       theta_Tau = theta0
      
      ! solve the phi by the provided chemical potential
       MDiffs(1) = muTau
       MDiffs(2) = mu0
       MDiffs(3) = Rgas
       MDiffs(4) = theta_tau
       MDiffs(5) = chi
       MDiffs(6) = Vmol
       MDiffs(7) = Kbulk
       Mdiffs(8) = detF
       Mdiffs(9) = detFs_t
      
       call solvePhi(C_t,C_tau,MDiffs,nMDiffs,stat)
 
      ! the phi rate
       dCdt = (C_Tau - C_t)/dtime
 
       ! Compute the perturbation on the chemical potential
       if(dabs(muTau).gt.one) then
            deltaMu = PerMu  * abs(MuTau)
       else
            deltaMU = 1.d-8
       endif
 
       ! solve the phi by the provided chemical potential
       MDiffs = zero
       MDiffs(1) = muTau + deltaMu
       MDiffs(2) = mu0
       MDiffs(3) = Rgas
       MDiffs(4) = theta_tau
       MDiffs(5) = chi
       MDiffs(6) = Vmol
       MDiffs(7) = Kbulk
       Mdiffs(8) = detF
       Mdiffs(9) = detFs_t  
       call solvePhi(C_t,C_Per,MDiffs,nMDiffs,stat)  
       dPdt_per = zero
       dPdt_per = (C_Per - C_t)/dtime
 
       ! solve the phi by the provided chemical potential
       MDiffs = zero
       MDiffs(1) = muTau - deltaMu
       MDiffs(2) = mu0
       MDiffs(3) = Rgas
       MDiffs(4) = theta_tau
       MDiffs(5) = chi
       MDiffs(6) = Vmol
       MDiffs(7) = Kbulk
       Mdiffs(8) = detF
       Mdiffs(9) = detFs_t  
       call solvePhi(C_t,C_M,MDiffs,nMDiffs,stat)
       dPdt_m = zero
       dPdt_m = (C_M - C_t)/dtime
    
       dCdotdMu = zero   
       dCdotdMu = (dPdt_per - dPdt_m)/(two * deltaMu)
       
      return
	end subroutine SolvFracInfTot
	
C****************************************************************************
!=========================
      !! subroutine SolvFracInfTot(phi_Tau,dPhidotdMu,dPhidt,MatProp,
      !!1 nMatProp,MuTau,MuRef,Phi_T,detF,sigHD,dtime,noel,detFs_t)
      !! use SolvTENS
      !! implicit none
      !!   
      !! integer I,J,K,L,M,N,nMatProp,IMax,KIter,Stat,noel
      !! real*8 phi_Tau,MatProp(nMatProp),MuTau,MuRef,phi_t
      !! real*8 FerrS,JacbK,JacbInv,dPhi,detF,sigHD
      !! real*8 one,two,zero,three,half,ErrorC
      !! real*8 Phi_Temp,dPhiError,dtime,PerMu,detFs_t
      !! real*8 dPhidotdMu,deltaMu,dPhidt,dPdt_per,dPdt_m
      !! real*8 Mu_Per,Mu_M,Phi_Per,phi_M
      !! real*8 Phi0
      !!
      !! parameter(IMax=50,zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5)
      !! parameter(PerMu = 1.e-8)
      !! 
      !! ! ===================
      !! ! here we found using the following initial values can help to solving the phi
      !!
      !! ! solve the phi by the provided chemical potential
      !!   call SubSolvFrac(Phi_Tau,Stat,MatProp,nMatProp,MuTau,MuRef,Phi_T,detF,sigHD,detFs_t)
      !!
      !! ! the phi rate
      !!   dPhidt = zero
      !!   dPhidt = (Phi_Tau - Phi_t)/dtime
      !!
      !! ! solve the dphidotdmu
      !!
      !!   deltaMu    = zero
      !!  if(abs(MuTau).gt.one) then        
      !!      deltaMu = PerMu * abs(MuTau)
      !!  else
      !!      deltaMu = PerMu
      !!  endif
      !!    
      !!   Mu_Per = zero
      !!   Mu_Per = MuTau + deltaMu
      !!   call SubSolvFrac(phi_Per,Stat,MatProp,nMatProp,Mu_Per,MuRef,Phi_T,detF,sigHD,detFs_t)
      !!     
      !!   Mu_M = zero
      !!   Mu_M = MuTau - deltaMu
      !!   call SubSolvFrac(phi_M,Stat,MatProp,nMatProp,Mu_M,MuRef,Phi_T,detF,sigHD,detFs_t)
      !! 
      !!   dPdt_per = zero
      !!   dPdt_m   = zero
      !!   dPdt_per = (Phi_Per - Phi_T)/dtime
      !!   dPdt_m = (phi_M - Phi_T)/dtime
      !!
      !!   dPhidotdMu = zero
      !!   dPhidotdMu = (dPdt_per - dPdt_m)/(two * deltaMu)
      !! 
      !! return
      !! end subroutine SolvFracInfTot      
      
      
      ! subroutine to solve the fractoin 
	subroutine SubSolvFrac(phi_Tau,Stat,MatProp,nMatProp,MuTau,MuRef,phi_t,detF,sigHD,detFs_t)
      use SolvTENS
      implicit none
      
      integer I,J,K,L,M,N,nMatProp,IMax,KIter,Stat
      real*8 phi_Tau,MatProp(nMatProp),MuTau,MuRef,phi_t,detF,sigHD 
      real*8 FerrS,JacbK,JacbInv,dPhi
      real*8 one,two,zero,three,half,ErrorC
      real*8 Phi_Temp,dPhiError,detFs_t

      parameter(IMax=50,zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5,ErrorC=1.d-4)
        ! loop to obtain the fraction
      
      Stat = 0
         
      Phi_Temp = phi_t
      
      Do 100 I=1,IMax
          
          call SolvFrac(FerrS,JacbK,Stat,MuTau,MatProp,nMatProp,MuRef,Phi_Temp,detF,sigHD,detFs_t)
              
          dPhi = FerrS/JacbK
          
          Phi_Temp = Phi_Temp - dPhi
          dPhiError = abs(dPhi)
          
           write(80,*) 'I = ', I,dPhiError
           write(80,*) 'FerrS,JacbK', FerrS,JacbK           
           write(80,*) 'MuTau,Phi_Temp',MuTau,Phi_Temp


          if (dPhiError.le.ErrorC) then  
              phi_tau = Phi_Temp
              ! write(*,*) 'a convergent result is obtained'
              return            
          endif
          
          if(I.ge.IMax) then               
              Stat = 1  
              phi_tau = Phi_Temp
              write(*,*)'too many tries on the solving of the phi' 
              write(*,*) phi_t,MuRef
              write(*,*) MuTau,detFs_t
              stop
              return
          endif    
100     continue    
        
        return
        end subroutine SubSolvFrac
        
****************************************************************************
        subroutine SolvFrac(FerrS,JacbK,Stat,Mu,MatProp,nMatProp,MuRef,Phi,detF,sigHD,detFs_t)
        use SolvTENS
        
        integer I,J,K,L,M,N,nMatProp,IMax,KIter,Stat,Ncl
        real*8 phi_Tau,MatProp(nMatProp),Mu,MuRef,phi,detF,sigHD
        real*8 Vd,Vf,VPcl,RT,KBulk,temp1,temp2,temp3
        real*8 FerrS,JacbK,Vmol,detFs_t
        
        real*8 zero,one,two,three,half,LangMult
        parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5,LangMult=1.e8)
        parameter(Ncl=10)
        
        Stat = 0
        
        ! material proprties
        
        Vmol = MatProp(4)
        Rgas = MatProp(2)
        Temp = MatProp(1)
        KBulk = 1.e8
        MuRef = zero
        RT = Rgas*Temp

        !
        sigHD = zero
        temp1 = sigHD/RT
        
        FerrS = zero
        FerrS = (MuRef - MU)/RT + dlog(Vmol*Phi) + one ! - dlog(detFs_t) 
      
        ! Jacobean matrix
        temp2 = KBulk/RT
        JacbK = zero
        JacbK = one/Phi

        return
	  end subroutine SolvFrac    
	  
****************************************************************************
      subroutine Umat_ChemReact(JudgeReact,props,nprops,C_T,MuTau,CRate,dCRdMuM,
     1 Achem,ChemR,Ncl,time,dtime,kelem)
      
      implicit none
       
      integer i,j,k,l,m,n,nprops,nNode,NMax,CTemp,CPTot,NTSeg
      integer Ncl,kelem
      real*8 props(nprops),CE_t,CP_t,CPcl_t,CP_temp
      real*8 CE_tau,CP_tau,CPcl_tau
      real*8 ReTot1,ReTot21,ReTot22,ReTot31,ReTot32,ReTot33
      real*8 Kbulk,VmolP,VmolE,detFs
      real*8 CEG_t,Kr,CEG_tau,dtime,time
      real*8 mu_tau,muP_tau,muPcl_tau,Rgas,Temp,RT
      real*8 PhiPcl,Achem,ChemR,kf,kb,VmolPcl,dCPclrRatedmu      
      real*8 CegRRate,CpRRate,CpclRRate,dCeRRatedmu,dCpRRatedmu
      real*8 mPcl,mPf,mPd,C_T(3),MuTau(3),CRate(3),dCRdMuM(3,3),RIndex
      real*8 JudgeReact

      real*8 zero,one,ten
      parameter(zero=0.0,one=1.0,ten=10.d0)
      !--------------------------
      ! material parameters
      Rgas  = props(4)
      Temp  = props(3) 
      VmolE = props(6)
      VmolP = props(7)
      VmolPcl = Ncl*VmolP
      RT = Rgas*Temp

      !-----------------------------
      ! chemical affinity and reaction rate
      ! MuTau = [MuPcl,MuPf,MuPd]
       Achem = MuTau(1) + (Ncl-1)*MuTau(3) - Ncl*MuTau(2)
      
      ! ---------------------
      ! to realize a complete dissolution, we take following assumption 
       if (JudgeReact.eq.1.0) then
           kf    = props(21)*1.e-2      ! the forward chemical reaction
           kb    = props(22)*1.e-2      ! the backward chemical reaction
       else
           kf    = props(21)      ! the forward chemical reaction
           kb    = props(22)      ! the backward chemical reaction
       endif
       
        ChemR = (kf-kb)*Achem       
      !-------------------------
      ! solve the reaction rate
        CRate    =  zero
        CRate(1) = -1*ChemR
        CRate(2) =  Ncl*ChemR
        CRate(3) = -(Ncl-1)*ChemR

      !---------------------------------
      ! derivative to the chemical potential
        mPcl  = -1.0
        mPf   =  1.0*Ncl
        mPd   = -1.0*(Ncl-1)      

        dCRdMuM = zero
        dCRdMuM(1,1) = -mPcl*mPcl*(kf-kb)
        dCRdMuM(1,2) = -mPcl*mPf*(kf-kb)
        dCRdMuM(1,3) = -mPcl*mPd*(kf-kb)
        
        dCRdMuM(2,1) = -mPf*mPcl*(kf-kb)
        dCRdMuM(2,2) = -mPf*mPf*(kf-kb)
        dCRdMuM(2,3) = -mPf*mPd*(kf-kb)        
        
        dCRdMuM(3,1) = -mPd*mPcl*(kf-kb)
        dCRdMuM(3,2) = -mPd*mPf*(kf-kb)
        dCRdMuM(3,3) = -mPd*mPd*(kf-kb)
      !------------
      return
	end subroutine Umat_ChemReact
      
		
C**********************************************************************
C				Calculate the swelling deformation
C**********************************************************************
      subroutine Sub_Ch_Swell(kinc,noel,npt,DetFs)
      use global
      implicit none
      
      integer kinc,I,J,K,noel,npt
      real*8 DetFs,Iden(3,3),one,zero
      real*8 Ce_tau,Cp_tau,CPcl_tau
      parameter(one=1.d0,zero=0.d0)
      
      call onem(Iden,3)
      
      if (kinc.le.1)then
          detFs = one
	else
          detFs = GlobalDetFs(noel,npt)
      endif
      
      return
	end subroutine Sub_Ch_Swell
	
C**********************************************************************
C				     solving the stress
C**********************************************************************
      subroutine GetStrs(F_t,F_tau,props,nprops,SigStrs,KirStrs,
     1 SpTanMod,DetFs,theta_t,theta_tau)
       use global
       use SolvTENS
       implicit none
       
       integer i,j,k,l,m,n,nprops
       real*8  SigStrs(3,3),KirStrs(3,3)
       real*8  F_t(3,3),F_tau(3,3),props(nprops)
       real*8  FG_tau(3,3),FgInv(3,3),FInv(3,3)
       real*8  Be_tau(3,3),trBe,Bedis(3,3),trBedis,Bedis0(3,3)
       real*8  B0(3,3),Ident2(3,3),Fe_tau(3,3),C(3,3)
       real*8 Gshear,Eyoung,poisson,detFe,SpTanMod(3,3,3,3),detFs,detF
       real*8 Kbulk,alpha,theta_t,theta_tau
       real*8 CV(3),SV(3)
       DOUBLE PRECISION,dimension(3,3):: LogC,EH,RE,EHD,CInv,MandStrs,LogCD,CR,SR,B,T_tau
       DOUBLE PRECISION,dimension(3,3,3,3) :: Ident4,Ident4S,Ident4D,Ident4A,Ident4T
       DOUBLE PRECISION,dimension(3,3,3,3) :: dCdF,PC,dCdFFT,LME,DSDEC,dTRdF
       real*8 one,two,three,half,zero,third
       parameter(one=1.0,half=0.5,two=2.d0,three=3.d0,zero=0.d0,third=1.d0/3.d0)

C	Read material properties
	Eyoung   =  props(1)
	poisson  =  props(2)
C	Calculate elastic modulus
	Gshear  =  Eyoung/(two*(one+poisson))
      Kbulk	=  Eyoung/(three*(one-two*poisson))

       ! define identity tensors
       call Ident4Def(Ident4,Ident4S,Ident4D,Ident4A,Ident4T,Ident2)
       
       ! compute some matrix
       C = matmul(transpose(F_tau),F_tau)
       B = matmul(F_tau,transpose(F_tau))
       call mdet(F_tau,detF) 
        
      T_tau = (Gshear*(B-Ident2) + Kbulk*dlog(detF)*Ident2)/detF
	SigStrs = T_tau
      WRITE(80,*)'F_tau',F_tau
	!WRITE(80,*)'C',C
	!WRITE(80,*)'B',B
	WRITE(80,*)'T_tau',T_tau
       ! Compute dTRdF, the so-called material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3          
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + Gshear*Ident2(i,k)*Ident2(j,l)
     +                 + Gshear*Finv(l,i)*Finv(j,k)
     +                 + Kbulk*Finv(j,i)*Finv(l,k)
     +                 - Kbulk*dlog(detF)*Finv(l,i)*Finv(j,k)
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate the so-called spatial tangent modulus, based
      !  on the push forward of the material tangent modulus
      !
      SpTanMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) + 
     +                       (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo


	return
	end subroutine GetStrs
C**********************************************************************
	
C**********************************************************************
      subroutine jac2D(SpTanMod,ddsdde)
      real*8 SpTanMod(3,3,3,3),ddsdde(4,4)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)

	end subroutine jac2D
C**********************************************************************
	
C**********************************************************************
      subroutine jac21D(SpUT,ddsddt)
      real*8 SpUT(3,3),ddsddt(4)
      
      ddsddt(1) = SpUT(1,1)
      ddsddt(2) = SpUT(2,2)
      ddsddt(3) = SpUT(3,3)
      ddsddt(4) = SpUT(1,2)
      
	end subroutine jac21D
C**********************************************************************
	
C**********************************************************************
      subroutine jac3D(SpTanMod,ddsdde)

      real*8 SpTanMod(3,3,3,3),ddsdde(6,6)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)
      ddsdde(1,5) = SpTanMod(1,1,1,3)
      ddsdde(1,6) = SpTanmod(1,1,2,3)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)
      ddsdde(2,5) = SpTanMod(2,2,1,3)
      ddsdde(2,6) = SpTanmod(2,2,2,3)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)
      ddsdde(3,5) = SpTanMod(3,3,1,3)
      ddsdde(3,6) = SpTanmod(3,3,2,3)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)
      ddsdde(4,5) = SpTanMod(1,2,1,3)
      ddsdde(4,6) = SpTanmod(1,2,2,3)

      ddsdde(5,1) = SpTanMod(1,3,1,1)
      ddsdde(5,2) = SpTanMod(1,3,2,2)
      ddsdde(5,3) = SpTanMod(1,3,3,3)
      ddsdde(5,4) = SpTanMod(1,3,1,2)
      ddsdde(5,5) = SpTanMod(1,3,1,3)
      ddsdde(5,6) = SpTanmod(1,3,2,3)

      ddsdde(6,1) = SpTanMod(2,3,1,1)
      ddsdde(6,2) = SpTanMod(2,3,2,2)
      ddsdde(6,3) = SpTanMod(2,3,3,3)
      ddsdde(6,4) = SpTanMod(2,3,1,2)
      ddsdde(6,5) = SpTanMod(2,3,1,3)
      ddsdde(6,6) = SpTanmod(2,3,2,3)

      return
	end subroutine jac3D    
C**********************************************************************
	
C**********************************************************************
      subroutine jac31D(SpUT,ddsddt)
      real*8 SpUT(3,3),ddsddt(6)
      ddsddt(1) = SpUT(1,1)
      ddsddt(2) = SpUT(2,2)
      ddsddt(3) = SpUT(3,3)
      ddsddt(4) = SpUT(1,2)
      ddsddt(5) = SpUT(1,3)
      ddsddt(6) = SpUT(2,3)

	end subroutine jac31D      
C**********************************************************************
	
C**********************************************************************
	SUBROUTINE M3INV(A,AINV)

C 	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C	[A] AND PLACES THE RESULT IN [AINV]. 
C 	IF DET(A) IS ZERO, THE CALCULATION
C 	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

	REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C	A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
C	DET		-- THE COMPUTED DETERMINANT OF [A].
C	ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C			   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C			   IS CALLED THE COFACTOR OF A(I,J).
C	AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C			   OBTAINED BY REPLACING EACH ELEMENT OF
C			   [A] BY ITS COFACTOR, AND THEN TAKING
C			   TRANSPOSE OF THE RESULTING MATRIX.
C	AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C			   [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------

	CALL MDET(A,DET)
	IF ( DET .EQ. 0.D0 ) THEN
	  write(*,10)
	  STOP
	ENDIF
	CALL MCOFAC(A,ACOFAC)
	CALL MTRANS(ACOFAC,AADJ)
	DO 1 I = 1,3
	DO 1 J = 1,3
	     AINV(I,J) = AADJ(I,J)/DET
1	CONTINUE
10	FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +         10X,'PROGRAM TERMINATED')

	RETURN
	END SUBROUTINE M3INV
C**********************************************************************
	SUBROUTINE SHAPEFUN(AN,dNdxi,xi)
	INCLUDE 'ABA_PARAM.INC'
	Real*8 AN(4),dNdxi(4,2)
	Real*8 XI(2)
	PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)
C
C     Values of shape functions as a function of local coord.
	AN(1) = ONE/FOUR*(ONE-XI(1))*(ONE-XI(2))
	AN(2) = ONE/FOUR*(ONE+XI(1))*(ONE-XI(2))
	AN(3) = ONE/FOUR*(ONE+XI(1))*(ONE+XI(2))
	AN(4) = ONE/FOUR*(ONE-XI(1))*(ONE+XI(2))
C
C     Derivatives of shape functions respect to local coordinates
       DO I=1,4
         DO J=1,2
             dNdxi(I,J) =  ZERO
         END DO
       END DO
       dNdxi(1,1) =  MONE/FOUR*(ONE-XI(2))
       dNdxi(1,2) =  MONE/FOUR*(ONE-XI(1))
       dNdxi(2,1) =  ONE/FOUR*(ONE-XI(2))
       dNdxi(2,2) =  MONE/FOUR*(ONE+XI(1))
       dNdxi(3,1) =  ONE/FOUR*(ONE+XI(2))
       dNdxi(3,2) =  ONE/FOUR*(ONE+XI(1))
       dNdxi(4,1) =  MONE/FOUR*(ONE+XI(2))
       dNdxi(4,2) =  ONE/FOUR*(ONE-XI(1))
       RETURN
	END

C****************************************************************************
C					DEFINE THE IDENTITY MATRIX
C****************************************************************************
	subroutine onem(A,N)
	! A: The matrix of N by N
	! N: The dimension of the matrix
	implicit none
	integer i,j,N

	real*8 A(N,N)
	do i=1,N
		do J=1,N
		if (i .eq. j) then
			A(i,j) = 1.0
		else
			A(i,j) = 0.0
		end if
		end do
	end do
	return
	end subroutine onem

****************************************************************************

	subroutine mdet(A,det)

	! This subroutine calculates the determinant
	! of the 3 by 3 matrix [A]

	implicit none

	real*8  A(3,3),det


	det = A(1,1)*A(2,2)*A(3,3) 
	1	  + A(1,2)*A(2,3)*A(3,1)
	2	  + A(1,3)*A(2,1)*A(3,2)
	3	  - A(3,1)*A(2,2)*A(1,3)
	4	  - A(3,2)*A(2,3)*A(1,1)
	5	  - A(3,3)*A(2,1)*A(1,2)


	return
	end subroutine mdet

C**********************************************************************
      subroutine SolvLogC(C,LogC,EH,CV)
      implicit none
        
      integer i,j,k,Noel
      real*8 LogC(3,3),C(3,3),CV(3,3),CE(3)
      real*8 eigE(3,3),U(3,3),EH(3,3),eigEH(3,3)
      real*8 zero,one,two,half
      parameter(zero=0.0,one=1.0,two=2.d0,half=0.5)

      call SPECTRAL(C,CE,CV)
      
      eigE = zero
      do i=1,3
          eigE(i,i)=dlog(CE(i))
      enddo

      eigEH = zero
      do i=1,3
         eigEH(i,i)=dlog(sqrt(CE(i)))
      enddo        
      LogC = matmul(CV,matmul(eigE,transpose(CV)))
      EH= matmul(CV,matmul(eigEH,transpose(CV)))

      return
	end subroutine SolvLogC  
C**********************************************************************
	
C**********************************************************************
      subroutine Ident4Def(Ident4,Ident4S,Ident4D,Ident4A,Ident4T,Iden)
      implicit none
      integer i,j,k,L,m,n
      real*8 Ident4S(3,3,3,3),Ident4D(3,3,3,3),Ident4A(3,3,3,3)
      real*8 Ident4T(3,3,3,3),Ident4(3,3,3,3),Iden(3,3)
      real*8 three,zero
      parameter(three=3.d0,zero=0.d0)
      !
      ! 2nd Identity matrix
      !
      call onem(Iden,3)    
      Ident4S = zero
      Ident4D = zero
      Ident4A = zero
      Ident4T = zero
      Ident4 = zero
      ! Ident4A
      do i=1,3
        do j=1,3
          do k=1,3
            do L=1,3
              Ident4A(i,j,k,L) = Iden(i,k)*Iden(j,L)
              enddo
            enddo
        enddo
      enddo
      ! Ident4T        
      do i=1,3
        do j=1,3
          do k=1,3
            do L=1,3
              Ident4T(i,j,k,L) = Iden(i,L)*Iden(j,k)
              enddo
            enddo
        enddo
      enddo
      ! Ident4 
      do i=1,3
        do j=1,3
          do k=1,3
            do L=1,3
              Ident4(i,j,k,L) = Iden(i,j)*Iden(k,L)
              enddo
            enddo
        enddo
      enddo   
      ! Ident4S,Ident4D      
      Ident4S = 0.5*(Ident4A + Ident4T)
      Ident4D = Ident4S - Ident4/three

      return
	end subroutine  Ident4Def    
	
C**********************************************************************
	SUBROUTINE MCOFAC(A,ACOFAC)
 
C 	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C 	AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************

	REAL*8  A(3,3), ACOFAC(3,3)

	ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
	ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
	ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
	ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
	ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
	ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
	ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
	ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
	ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

	RETURN
	END SUBROUTINE MCOFAC
	
C**********************************************************************
	SUBROUTINE MTRANS(A,ATRANS)
 
C	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

	REAL*8 A(3,3),ATRANS(3,3)

	DO 1 I=1,3
 	  DO 1 J=1,3
	    ATRANS(J,I) = A(I,J)
1	CONTINUE

	RETURN
	END SUBROUTINE MTRANS
	
C**********************************************************************
C	THE FOLLOWING SUBROUTINES CALCULATE THE SPECTRAL
C	DECOMPOSITION OF A SYMMETRIC THREE BY THREE MATRIX
C**********************************************************************
	SUBROUTINE SPECTRAL(A,D,V)
C
C	THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
C	A SYMMETRIC 3 BY 3 MATRIX [A]. 
C
C	THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
C	EIGENVALUES IN ASCENDING ORDER, AND
C	A MATRIX [V] WHOSE COLUMNS CONTAIN THE CORRESPONDING
C	EIGENVECTORS.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER(NP=3)
	DIMENSION D(NP),V(NP,NP)
	DIMENSION A(3,3),E(NP,NP)

	DO 2 I = 1,3
          DO 1 J= 1,3
            E(I,J) = A(I,J)
1	  CONTINUE
2	CONTINUE

	CALL JACOBI(E,3,NP,D,V,NROT)
	!CALL EIGSRT(D,V,3,NP)

	RETURN
      END
C**********************************************************************
	SUBROUTINE EIGSRT(D,V,N,NP)

C	GIVEN THE EIGENVALUES [D] AND EIGENVECTORS [V] AS OUTPUT FROM
C	JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO ASCENDING ORDER, 
C	AND REARRANGES THE COLUMNS OF [V] ACCORDINGLY.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", P. 348.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION D(NP),V(NP,NP)

	DO 13 I = 1,N-1
	  K = I
	  P = D(I)
	  DO 11 J = I+1,N
	    IF (D(J) .GE. P) THEN
	      K = J
	      P = D(J)
	    END IF
11	  CONTINUE
	  IF (K .NE. I) THEN
	    D(K) = D(I)
	    D(I) = P
	    DO 12 J = 1,N
	      P = V(J,I)
	      V(J,I) = V(J,K)
	      V(J,K) = P
12	    CONTINUE
  	  ENDIF
13	CONTINUE

	RETURN
      END
C**********************************************************************
	SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

C	COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C	MATRIX [A], WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL 
C	NP BY BP ARRAY. ON OUTPUT, ELEMENTS OF [A] ABOVE THE DIAGONAL 
C	ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
C	AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
C	VECTOR D RETURNS THE EIGENVALUES OF [A] IN ITS FIRST N ELEMENTS.
C	[V] IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
C	[A] WHOSE COLUMNS CONTAIN, ON OUTPUT, THE NORMALIZED
C	EIGENVECTORSOF [A]. NROT RETURNS THE NUMBER OF JACOBI ROTATIONS
C	WHICH WERE REQUIRED.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", PAGE 346.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NMAX =100)
	DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

C	INITIALIZE [V] TO THE IDENTITY MATRIX

	DO 12 IP = 1,N	
	  DO 11 IQ = 1,N
	    V(IP,IQ) = 0.D0
11        CONTINUE
          V(IP,IP) = 1.D0
12	CONTINUE

C	INITIALIZE [B] AND [D] TO THE DIAGONAL OF [A], AND Z TO ZERO.
C	THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
C	IN EQUATION (11.1.14)

	DO 13 IP = 1,N
	  B(IP) = A(IP,IP)
	  D(IP) = B(IP)
	  Z(IP) = 0.D0
13	CONTINUE
C
	NROT = 0
	DO 24 I = 1,50

C	SUM OFF-DIAGONAL ELEMENTS

          SM = 0.D0
          DO 15 IP = 1, N-1
            DO 14 IQ = IP + 1, N
	      SM = SM + DABS ( A(IP,IQ ))
14          CONTINUE
15        CONTINUE

C	IF SUM = 0., THEN RETURN. THIS IS THE NORMAL RETURN
C	WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE 
C	UNDERFLOW.

          IF ( SM .EQ. 0.D0) RETURN
C
C	  IF ( SM .LT. 1.0D-15) RETURN

C	IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
C	|A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE, 
C	SEE EQUATION (11.1.25). THEREAFTER TRESH = 0.

          IF ( I .LT. 4) THEN
            TRESH = 0.2D0*SM/N**2
          ELSE
            TRESH = 0.D0
          ENDIF
C
          DO 22 IP = 1, N-1
            DO 21 IQ = IP+1,N
              G = 100.D0*DABS(A(IP,IQ))

C	AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
C	OFF-DIAGONAL ELEMENT IS SMALL.

	      IF ((I .GT. 4) .AND. (DABS(D(IP))+G .EQ. DABS(D(IP)))
     +            .AND. ( DABS(D(IQ))+G .EQ. DABS(D(IQ)))) THEN
                A(IP,IQ) = 0.D0
              ELSE IF ( DABS(A(IP,IQ)) .GT. TRESH) THEN
                H = D(IQ) - D(IP)
                IF (DABS(H)+G .EQ. DABS(H)) THEN

C	T = 1./(2.*THETA), EQUATION(11.1.10)

	          T =A(IP,IQ)/H
	        ELSE
	          THETA = 0.5D0*H/A(IP,IQ)
	          T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
	          IF (THETA .LT. 0.D0) T = -T
	        ENDIF
	        C = 1.D0/DSQRT(1.D0 + T**2)
	        S = T*C
	        TAU = S/(1.D0 + C)
	        H = T*A(IP,IQ)
	        Z(IP) = Z(IP) - H
	        Z(IQ) = Z(IQ) + H
	        D(IP) = D(IP) - H
	        D(IQ) = D(IQ) + H
	        A(IP,IQ) = 0.D0

C	CASE OF ROTATIONS 1 <= J < P
				
	        DO 16 J = 1, IP-1
	          G = A(J,IP)
	          H = A(J,IQ)
	          A(J,IP) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
16	        CONTINUE

C	CASE OF ROTATIONS P < J < Q

	        DO 17 J = IP+1, IQ-1
	          G = A(IP,J)
	          H = A(J,IQ)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
17	        CONTINUE

C	CASE OF ROTATIONS Q < J <= N

	        DO 18 J = IQ+1, N
                  G = A(IP,J)
	          H = A(IQ,J)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(IQ,J) = H + S*(G - H*TAU)
18	        CONTINUE
	        DO 19 J = 1,N
	          G = V(J,IP)
	          H = V(J,IQ)
	          V(J,IP) = G - S*(H + G*TAU)
	          V(J,IQ) = H + S*(G - H*TAU)
19	        CONTINUE
	        NROT = NROT + 1
              ENDIF
21	    CONTINUE
22	  CONTINUE

C	UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z

	  DO 23 IP = 1, N
	    B(IP) = B(IP) + Z(IP)
	    D(IP) = B(IP)
	    Z(IP) = 0.D0
23	  CONTINUE
24	CONTINUE

C	IF THE ALGORITHM HAS REACHED THIS STAGE, THEN
C	THERE ARE TOO MANY SWEEPS, PRINT A DIAGNOSTIC
C	AND STOP.

	WRITE (*,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'

	RETURN
	END      
C**********************************************************************

	
C**********************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape2Da'
         call xit
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
	end subroutine mapShape2Da
C**********************************************************************

	
C**********************************************************************
      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape2D'
         call xit
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
	end subroutine mapShape2D
C**********************************************************************

C**********************************************************************
      subroutine xintSurf2D2pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(2),yLocal(2),w(2),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = -dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

	end subroutine xintSurf2D2pt
C**********************************************************************

C**********************************************************************
      subroutine xintSurf2D3pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(3),yLocal(3),w(3),zero,one,two,three,five,eight,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,five=5.d0,
     +     eight=8.d0,nine=9.d0)


      ! Gauss weights
      !
      w(1) = five/nine
      w(2) = eight/nine
      w(3) = five/nine
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = -one
         xLocal(2) = zero
         yLocal(2) = -one
         xLocal(2) = dsqrt(three/five)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(three/five)
         xLocal(2) = one
         yLocal(2) = zero
         xLocal(3) = one
         yLocal(3) = dsqrt(three/five)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = one
         xLocal(2) = zero
         yLocal(2) = one
         xLocal(3) = dsqrt(three/five)
         yLocal(3) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(three/five)
         xLocal(2) = -one
         yLocal(2) = zero
         xLocal(3) = -one
         yLocal(3) = -dsqrt(three/five)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

	end subroutine xintSurf2D3pt
C**********************************************************************

C**********************************************************************
      subroutine xintSurf2D1pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),w(1),zero,one,two
      parameter(zero=0.d0,one=1.d0,two=2.d0)


      ! Gauss weights
      !
      w(1) = two
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = zero
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

	end subroutine xintSurf2D1pt
C**********************************************************************

C**********************************************************************
      subroutine computeSurf(xLocal,yLocal,face,coords,sh,ds)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the length ds, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 4-node quadrilateral elements

      implicit none

      integer face

      real*8 xLocal,yLocal,ds,dshxi(4,2),sh(4),dXdXi,dXdEta,dYdXi
      real*8 dYdEta,one,coords(2,4),fourth,shape,normal(2,1)
      parameter(one=1.d0,fourth=1.d0/4.d0)

      sh(1) = fourth*(one - xLocal)*(one - yLocal)
      sh(2) = fourth*(one + xLocal)*(one - yLocal)
      sh(3) = fourth*(one + xLocal)*(one + yLocal)
      sh(4) = fourth*(one - xLocal)*(one + yLocal)
      
      dshxi(1,1) = -fourth*(one - yLocal)
      dshxi(1,2) = -fourth*(one - xLocal)
      dshxi(2,1) = fourth*(one - yLocal)
      dshxi(2,2) = -fourth*(one + xLocal)
      dshxi(3,1) = fourth*(one + yLocal)
      dshxi(3,2) = fourth*(one + xLocal)
      dshxi(4,1) = -fourth*(one + yLocal)
      dshxi(4,2) = fourth*(one - xLocal)

      dXdXi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)
     +     + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
      dXdEta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)
     +     + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
      dYdXi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)
     +     + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
      dYdEta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)
     +     + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)


      ! Jacobian of the mapping
      !
      if((face.eq.2).or.(face.eq.4)) then
         ds = dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
      elseif((face.eq.1).or.(face.eq.3)) then
         ds = dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
      else
         write(*,*) 'never should get here'
         call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.2).or.(face.eq.4)) then
         normal(1,1) = dYdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         normal(2,1) = -dXdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         if(face.eq.4) normal = -normal
      elseif((face.eq.1).or.(face.eq.3)) then
         normal(1,1) = dYdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         normal(2,1) = -dXdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         if(face.eq.3) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif

      return
	end subroutine computeSurf
C**********************************************************************

C**********************************************************************
      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element


      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
	end subroutine calcShape2DLinear
C**********************************************************************

C**********************************************************************
      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
	end subroutine xint2D4pt
C**********************************************************************

C**********************************************************************
      subroutine xint2D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2), w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0


      return
	end subroutine xint2D1pt
C**********************************************************************

C**********************************************************************
      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
	end subroutine matInv2D
		
C**********************************************************************
      subroutine solvePhi(con_t,con_tau,args,nargs,stat)
       
       implicit none
       
       integer nargs,NeoHookean,Langevin,material,stat
       parameter(NeoHookean=1,Langevin=2)

       real*8 args(nargs),f,df,muLi,muLi0,Rgas,theta,chi,VLimol,pnewdt
       real*8 Gshear,Kbulk,detF,con,RT,Jcon,Con_tau,Con_t,con_temp,ki

       real*8 zero,one,two,three,third,MuBound,error,errorC,con_tempt
       real*8 fifty,TrM0
       parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0
     +,errorC=1.d-6,fifty=50.d0)

      ! Obtain relevant quantities
      !
        muLi  =   args(1)
        muLi0 =   args(2)
        Rgas  =   args(3)
        theta =   args(4)
        chi   =   args(5)
        VLimol=   args(6)
        Kbulk =   args(7)
        detF  =   args(8)
        TrM0  =   args(9)

      ! Compute the useful quantity
      !
        RT = Rgas * theta

      ! Compute the residual
      !

      call phiFunc(con_t,f,df,args,nargs)
         
          if (df.eq.zero) then
              con_tau = con_t
          else
              con_temp = con_t - f/df
              error = abs(con_temp-con_t)
          endif
                
          Jcon =  zero   
          ki   =  zero
          
      ! itetrative to solve the convergent convergence when provided with the chemical pontential and 
      ! deformation
      do while (Jcon == zero) 
          ki  =  ki+1
          if (error.lt.errorC) then
              Jcon = one
              con_tau = con_temp
          else
              call phiFunc(con_temp,f,df,args,nargs)
             
              if (df.eq.zero) then
                con_tau = Con_temp
                Jcon=one
              else
                con_tempt= con_temp
                con_temp = Con_tempt - f/df

                error = abs(con_temp-con_tempt)
                Jcon=zero
              end if
          endif
         ! when the negative concentration is obtained, there must
         ! be something wrong here
         stat = 1 
         if (con_tau.le.zero) then
            stat = 0
            ! write(80,*)'here,it is because the solved concentration is less than zero, the Pnewdt = 0.5'
         endif 
      
         ! When too much iterative steps undergo, there must be something wrong here
         if(ki.eq.fifty) then
             write(80,*) 'Too much iterative steps for the 
     1                    convergent concentration'
             call xit
         endif
      
      end do
      
      return
      end subroutine solvePhi      
C----------------------------------------------------------------         
       subroutine phiFunc(con,f,df,args,nargs)
       
       implicit none
       
       integer nargs,NeoHookean,Langevin,material
       parameter(NeoHookean=1,Langevin=2)

       real*8 args(nargs),f,df,muLi,muLi0,Rgas,theta,chi,VLimol
       real*8 Gshear,Kbulk,detF,con,RT,TrM0

       real*8 zero,one,two,three,third,MuBound
       parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0)

      ! Obtain relevant quantities
      !
        muLi  = args(1)
        muLi0 = args(2)
        Rgas  = args(3)
        theta = args(4)
        chi   = args(5)
        VLimol= args(6)
        Kbulk = args(7)
        detF  = args(8)
        TrM0  = args(9)

      ! Compute the useful quantity
      !
        RT = Rgas*theta

      ! Compute the residual
      !
	  
      f = (muLi0 - muLi - TrM0 * VLimol)/RT
     + + dlog(one - (one/(one+VLimol*con))) + (one/(one+VLimol*con))
      !
       if(con.lt.0.00000001d0) then
          df = zero
       else
          df = one/con - (VLimol/(one+VLimol*con)) 
     + - (VLimol/(one+VLimol*con)**two)
       endif

       return
       end subroutine phiFunc 
****************************************************************************      
	SUBROUTINE RAMP(A,B)
C 	This is a Ramp fuction
	IMPLICIT NONE
	REAL*8 A,B,ZERO
	PARAMETER(ZERO=0.d0)
	IF(A .LT. ZERO) THEN
		B = ZERO
	ELSE
		B = A
      END IF
      END SUBROUTINE RAMP

****************************************************************************          
       subroutine SolvEnergy(F_tau,T_tau,energy)
      implicit none
      integer I,J,K,L,M,N
      double precision,dimension(3,3)::EVect,PKVect,Ident2,E,PK2,FInv,F_tau,T_tau
      double precision,dimension(3)::EVal,PKVal
      real*8 one,two,zero,three,half,detF,energy
      parameter(one=1.d0,zero=0.d0,half=0.5d0)
      
      call onem(Ident2,3)
      call mdet(F_tau,detF)
      call M3INV(F_tau,FInv)
      
      ! To sovle the Green strain
      
      E = half*(matmul(transpose(F_tau),F_tau)-Ident2)
      call SPECTRAL(E,EVal,EVect)
      
      ! To solve the 2nd PK
      
      PK2 = detF*matmul(FInv,matmul(T_tau,transpose(FInv)))

      call SPECTRAL(PK2,PKVal,PKVect)
      
      Energy = zero
      do I=1,3
          if(PKVal(I).gt.zero)then
              !if(PKVal(I).Lt.zero)then
              !    write(*,*)'it is a big problem here'
              !    write(*,*)'F_tau',F_tau
              !    write(*,*)'T_tau',T_Tau
              !    
              !    
                  !write(80,*)'EVal',EVal
              !    write(*,*)'EVect',EVect
                  !write(80,*)'PKVal',PKVal
              !    write(*,*)'PKVect',PKVect
              !    call exit
              !endif
              Energy = half*Eval(I)*PKVal(I) + Energy
          endif
      enddo
      
      end subroutine  SolvEnergy