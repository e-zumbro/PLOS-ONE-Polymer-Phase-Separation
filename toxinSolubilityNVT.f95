program toxinSolubilityNVT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Main simulation file
!
! Emiko Zumbro
! ezumbro@mit.edu
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE parameters
USE routinesMultTox
USE functions_Emi

IMPLICIT NONE

! Initiate all of the variables
INTEGER :: i,j,n,totSteps,iStart

! positions
REAL(DP), DIMENSION(:), ALLOCATABLE :: rx,ry,rz
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: posTox ! dim x maxTox
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dRTox ! beads x maxTox

! forces
REAL(DP), DIMENSION(:), ALLOCATABLE :: fljx,fljy,fljz,fljxTox,fljyTox,fljzTox ! beads
REAL(DP), DIMENSION(:), ALLOCATABLE :: fljxToxTox,fljyToxTox,fljzToxTox ! maxTox
REAL(DP), DIMENSION(:), ALLOCATABLE :: fljxToxTot,fljyToxTot,fljzToxTot ! beads
REAL(DP), DIMENSION(:), ALLOCATABLE :: fspringx,fspringy,fspringz,fspringxTox,fspringyTox,fspringzTox
REAL(DP), DIMENSION(:), ALLOCATABLE :: FbindMag,fbindx,fbindy,fbindz
REAL(DP), DIMENSION(:), ALLOCATABLE :: fwx,fwy,fwz ! polymer beads

! delta positions and cummulative potentials
REAL(DP), DIMENSION(:), ALLOCATABLE :: delXToxUnit,delYToxUnit,delZToxUnit
REAL(DP), DIMENSION(:), ALLOCATABLE :: dUdx,dUdy,dUdz,dUdxTox,dUdyTox,dUdzTox

! random multipliers
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: deltaWt,deltaWtTox
INTEGER :: nRun, nSeed
INTEGER*8, DIMENSION(:) :: readSeed(seedSize)
INTEGER*4, DIMENSION(:) :: seed(seedSize)

! Gaussian random number holders
REAL(DP) :: gausNum = 0.0_DP
INTEGER :: rowG,colG

! binding information
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: prevOmegaTot
INTEGER, DIMENSION(:,:), ALLOCATABLE :: prevOmega, newOmega
INTEGER, DIMENSION(:), ALLOCATABLE :: beadsBound, timeBoundUnbound, timeTypeBond
INTEGER, DIMENSION(:), ALLOCATABLE :: idxToxBnd, bdBndIdx
REAL(DP), DIMENSION(:), ALLOCATABLE :: delE_0
REAL(DP), DIMENSION(:), ALLOCATABLE :: delE_UB

! Variables for keeping toxin concentration constant
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dRToxTox, existTox

! Variables for assigning random indexes
INTEGER, DIMENSION(:), ALLOCATABLE :: old_idxs, rand_idxs
INTEGER :: nIdx = 0
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: delE0Tot

! Variables for neighbor lists
INTEGER, DIMENSION(:), ALLOCATABLE :: nNghbrsPP, nNghbrsToxTox, nNghbrsPolyTox
INTEGER, DIMENSION(:,:), ALLOCATABLE :: nghbrListPP, nghbrListToxTox, nghbrListPolyTox

! restart info
INTEGER :: restart

CHARACTER(80) :: restart_conf
CHARACTER(80) :: restart_rand
CHARACTER(80) :: configuration

restart_conf="dataOut/restart_conf.res"
restart_rand="dataOut/restart_rand.res"

OPEN (UNIT=20,FILE="dataOut/posToxX.dat",STATUS="REPLACE")
CLOSE (20)
OPEN (UNIT=20,FILE="dataOut/posToxY.dat",STATUS="REPLACE")
CLOSE (20)
OPEN (UNIT=20,FILE="dataOut/posToxZ.dat",STATUS="REPLACE")
CLOSE (20)
OPEN (UNIT=20,FILE="dataOut/confX.dat",STATUS="REPLACE")
CLOSE (20)
OPEN (UNIT=20,FILE="dataOut/confY.dat",STATUS="REPLACE")
CLOSE (20)
OPEN (UNIT=20,FILE="dataOut/confZ.dat",STATUS="REPLACE")
CLOSE (20)
OPEN (UNIT=20,FILE="dataOut/idxToxBnd.dat",STATUS="REPLACE")
CLOSE (20)
OPEN (UNIT=20,FILE="dataOut/numBnd.dat",STATUS="REPLACE")
CLOSE (20)

print *, "about to input.par"

OPEN(10,file="input.par",status="unknown")
REWIND (10)
READ (10,FMT=*)
READ (10,FMT=*) restart
READ (10,FMT=*)
READ (10,FMT=*) Nmer    
READ (10,FMT=*)
READ (10,FMT=*) beads
READ (10,FMT=*)
READ (10,FMT=*) delE_0_center
READ (10,FMT=*)
READ (10,FMT=*) bindingSites
READ (10,FMT=*)
READ (10,FMT=*) tot_t
READ (10,FMT=*)
READ (10,FMT=*) isPolymer
READ (10,FMT=*)
READ (10,FMT=*) randNewToxin
READ (10,FMT=*)
READ (10,FMT=*) stdDevPolyAff
READ (10,FMT=*)
READ (10,FMT=*) nRun
READ (10,FMT=*)
READ (10,FMT=*) nTox
READ (10,FMT=*)
READ (10,FMT=*) epsToxTox
READ (10,FMT=*)
READ (10,FMT=*) gammaw
READ (10,FMT=*)
READ (10,FMT=*) epsilon
READ (10,FMT=*)
READ (10,FMT=*) epsToxPoly
CLOSE (10)

totSteps = tot_t/deltaT
print *, totSteps
!This size was calculated to have the same concentration of beads as I had in all the data I ran with only 10 beads 1.1*beads; 
sizeBox = ((11.0_DP**3.0_DP)*(64.0_DP/10.0_DP))**(1.0/3.0_DP) !((11.0_DP**3.0_DP)*(beads/10.0_DP))**(1.0/3.0_DP) ! 200
print *, "sizeBox = ", sizeBox

! This is for NVT (we aren't going to change the number of tox, so don't need
!	to allocate extra space in advance)
maxTox = nTox

! Get out your random seed and seed your function
!print *, "about to rSeed.par"
OPEN(10,file="rSeed.par",status="unknown")
REWIND (10)
DO nSeed = 1,seedSize
READ (10,FMT=*) readSeed(nSeed)
END DO

CLOSE (10)

seed = INT(readSeed,4)

! Initiate your random seed
CALL RANDOM_SEED(PUT=seed)
!print *, "size of seed = ", SIZEOF(seed)

print *, "seeded"

ALLOCATE (rx(beads),ry(beads),rz(beads),posTox(dim,maxTox),dRTox(beads,maxTox))
ALLOCATE (fljx(beads),fljy(beads),fljz(beads),fljxTox(beads),fljyTox(beads),fljzTox(beads))
ALLOCATE (fljxToxTox(maxTox),fljyToxTox(maxTox),fljzToxTox(maxTox))
ALLOCATE (fljxToxTot(beads),fljyToxTot(beads),fljzToxTot(beads))

ALLOCATE (fspringx(beads),fspringy(beads),fspringz(beads))
ALLOCATE (fspringxTox(beads),fspringyTox(beads),fspringzTox(beads))
ALLOCATE (fwx(beads),fwy(beads),fwz(beads))

ALLOCATE (fbindx(beads),fbindy(beads),fbindz(beads),FbindMag(beads))
ALLOCATE (prevOmegaTot(beads,bindingSites,maxTox),prevOmega(beads,bindingSites), newOmega(beads,bindingSites))
ALLOCATE (timeBoundUnbound(maxTox),timeTypeBond(maxTox))
ALLOCATE (delXToxUnit(beads), delYToxUnit(beads), delZToxUnit(beads))
ALLOCATE (deltaWt(dim,beads), deltaWtTox(dim,maxTox),beadsBound(beads))
ALLOCATE (dUdx(beads),dUdy(beads),dUdz(beads),dUdxTox(maxTox),dUdyTox(maxTox),dUdzTox(maxTox))
ALLOCATE (delE_0(beads),delE_UB(beads))
ALLOCATE (dRToxTox(maxTox,maxTox))
! ALLOCATE (nToxAggTime(timeAvgInt),dRNewTox(maxTox),closeToxIdx(maxTox))
ALLOCATE (idxToxBnd(beads),bdBndIdx(bindingSites))
ALLOCATE (old_idxs(Nmer),rand_idxs(Nmer/2))
ALLOCATE (delE0Tot(50,beads))
ALLOCATE (nNghbrsPP(beads),nghbrListPP(beads,beads),nNghbrsToxTox(maxTox),nghbrListToxTox(maxTox,maxTox))
ALLOCATE (nNghbrsPolyTox(maxTox),nghbrListPolyTox(maxTox,beads))

delE_0 = delE_0_center

delE_UB = -delE_0 + 0.5_DP
posTox = 0.0

print *,"allocated"

!initialize the configuration
IF (restart == 0) THEN

	iStart = 0
	ry=0.0
	rz=0.0

	! Only works for polymers with 64 beads or less
	CALL setInitialPos(rx,ry,rz)
	!print *, rx

	CALL setInitialPosTox(posTox(1,1:nTox),posTox(2,1:nTox),posTox(3,1:nTox),rx,ry,rz)

	!print *, "posTox = ", posTox

ELSE
	! Start from a previous time interval - haven't implemented this

	print *,"Not implemented, start from time 0"

END IF

print *, "iStart = ",iStart

OPEN(10,file="dataOut/Output.dat",status="REPLACE")
! WRITE (10,*)
! WRITE (10,*) "This is the file that contains the parameters"
! WRITE (10,*) "used in the simulation."
! WRITE (10,*)
WRITE (10,*) "startingStep",iStart
WRITE (10,*) "restart",restart
WRITE (10,*) "beads",beads    
WRITE (10,*) "Nmer",Nmer
WRITE (10,*) "bindingSites",bindingSites
WRITE (10,*) "maxTox",maxTox
WRITE (10,*) "nTox0",nTox
WRITE (10,*) "tot_t",tot_t
WRITE (10,*) "delE_0_center",delE_0_center
WRITE (10,*) "deltaT",deltaT
WRITE (10,*) "tot_t",tot_t
WRITE (10,*) "checkBindingInterval",checkBindingInterval
WRITE (10,*) "D",D
WRITE (10,*) "DTox",DTox
WRITE (10,*) "dim",dim
WRITE (10,*) "sizeBox",sizeBox
WRITE (10,*) "d_wrt",d_wrt
WRITE (10,*) "nRun",nRun
WRITE (10,*) "epsilonLJPP",epsilon
WRITE (10,*) "epsLJToxTox",epsToxTox
WRITE (10,*) "epsLJToxPoly",epsToxPoly
WRITE (10,*) "gammaw",gammaw
CLOSE (10)

OPEN (UNIT=10,FILE="dataOut/iPos.dat", STATUS="REPLACE")
	DO j=1,beads
	   WRITE (UNIT=10,FMT=*) rx(j),ry(j),rz(j)
	END DO
CLOSE (10)

OPEN (UNIT=10,FILE="dataOut/iPosTox.dat", STATUS="REPLACE")
	DO j=1,maxTox
	   WRITE (UNIT=10,FMT=*) posTox(1,j),posTox(2,j),posTox(3,j)
	END DO
CLOSE (10)

OPEN (UNIT=10,FILE="dataOut/randSeed.dat", STATUS="REPLACE")
   WRITE (UNIT=10,FMT=*) seed
CLOSE (10)

OPEN (UNIT=10,FILE="dataOut/delE0.dat", STATUS="REPLACE")
	DO j=1,beads
	   WRITE (UNIT=10,FMT=*) delE_0(j)
	END DO
CLOSE (10)

OPEN (UNIT=10,FILE="dataOut/delEUB.dat", STATUS="REPLACE")
	DO j=1,beads
	   WRITE (UNIT=10,FMT=*) delE_UB(j)
	END DO
CLOSE (10)

! Initiate everything to start unbound and unaggregated
IF (restart==0) THEN
	prevOmegaTot = 0
	prevOmega = 0
	timeBoundUnbound = 0
	idxToxBnd = 0
	timeTypeBond = 0
ELSE
	prevOmega = 0
	timeBoundUnbound = 0
	timeTypeBond = 0
ENDIF

print *, "made it to function"

DO i=iStart+1,totSteps

	DO colG=1,dim
		DO rowG=1,beads
			CALL gasdev_s(gausNum)
			deltaWt(colG,rowG) = gausNum
		END DO
	END DO

	deltaWt = deltaWt * SQRT(deltaT)
	
	DO colG=1,dim
		DO rowG=1,nTox
			CALL gasdev_s(gausNum)
			deltaWtTox(colG,rowG) = gausNum
		END DO
	END DO
	!print *, "deltaWtTox = ", deltaWtTox

	deltaWtTox = deltaWtTox * SQRT(deltaT)

	! Calculate the spring forces due to connectivity
	IF (isPolymer) THEN
		CALL springForces(fspringx,fspringy,fspringz,rx,ry,rz)
		IF (gammaw > 0.0) THEN
			! Calculate the wormlike chain forces to control bending
			CALL worm_like(fwx,fwy,fwz,rx,ry,rz)
		ELSE
			fwx = 0.0
			fwy = 0.0
			fwz = 0.0
		END IF
	ELSE
		! there is no connectivity and all the spring forces are 0
		fspringx = 0.0
		fspringy = 0.0
		fspringz = 0.0

		fwx = 0.0
		fwy = 0.0
		fwz = 0.0
	END IF

	! print *, "Calculated spring forces"
	IF (beads > 1) THEN
		! fljxTox=0.0
	 !    fljyTox=0.0
	 !    fljzTox=0.0

		! Update the neighbor list
		IF ((i==(iStart+1)).OR.(MODULO(i,checkBindingInterval/10)==0)) THEN
			CALL updateNeighborListPP(rx,ry,rz,nNghbrsPP,nghbrListPP)
		END IF

		! Calculate the lennard-jones excluded volume forces between inhibitor beads
		CALL rljmodNLBC(fljx,fljy,fljz,rx,ry,rz,nNghbrsPP,nghbrListPP)

	END IF

	fbindx = 0.0
	fbindy = 0.0
	fbindz = 0.0

	fljxTox=0.0
    fljyTox=0.0
    fljzTox=0.0

    fljxToxTot=0.0
    fljyToxTot=0.0
    fljzToxTot=0.0

    dUdxTox=0.0
    dUdyTox=0.0
    dUdzTox=0.0

	idxToxBnd = 0


    ! Find LJ excluded volume between all the toxins and themselves if there are multiple toxins
	IF (nTox > 1) THEN
		

		! Update the neighbor list
		IF ((i==(iStart+1)).OR.(MODULO(i,checkBindingInterval/10)==0)) THEN
		 	CALL updateNeighborListToxTox(posTox(1,1:nTox),posTox(2,1:nTox),posTox(3,1:nTox),nNghbrsToxTox,nghbrListToxTox)
		END IF

		! Calculate the lennard-jones excluded volume forces between inhibitor beads
		! Neighbor List incompatible with changing toxin concentration dynamically
		CALL rljToxToxNLBC(fljxToxTox,fljyToxTox,fljzToxTox,posTox(1,1:nTox),posTox(2,1:nTox),posTox(3,1:nTox), &
			& nNghbrsToxTox,nghbrListToxTox)

	END IF
	
	IF (beads > 0) THEN

			IF ((i==(iStart+1)).OR.(MODULO(i,checkBindingInterval/10)==0)) THEN
			 	CALL updateNeighborListPolyTox(rx,ry,rz,posTox(1,1:nTox),posTox(2,1:nTox),posTox(3,1:nTox),&
			 		& nNghbrsPolyTox,nghbrListPolyTox)
			END IF

		DO j = 1,nTox
			CALL rljToxNLBC(fljxTox,fljyTox,fljzTox,dRTox(:,j),rx,ry,rz,posTox(:,j),delXToxUnit,delYToxUnit,&
				& delZToxUnit,nNghbrsPolyTox(j),nghbrListPolyTox(j,:))

			prevOmega = prevOmegaTot(:,:,j)
			
			! If it is a binding interval, check to see if a binding event happens
			IF ((MODULO(i,checkBindingInterval)==0).AND.(delE_0_center<0)) THEN
				CALL bound(prevOmegaTot, prevOmega, dRTox(:,j), FbindMag, newOmega, delE_0, delE_UB)
			ELSE
				newOmega = prevOmega
				beadsBound = SUM(newOmega,2) ! by this one toxin
				FbindMag = -k_bind*beadsBound*(dRTox(:,j)-l_bind)
			END IF

			! Calculate the bound forces with their direction
			DO n = 1,beads
	            !del_ToxUnit points from the toxin to the inhibitor 
	            ! I think this makes it correct for the forces on the bead, but the toxin will have to have a negative of the fbind_
	            ! In order to have fbindx be dimension n, add the previous fbindx as you iterate through the toxins, j
	            fbindx(n) = fbindx(n) + FbindMag(n)*delXToxUnit(n)
	            fbindy(n) = fbindy(n) + FbindMag(n)*delYToxUnit(n)
	            fbindz(n) = fbindz(n) + FbindMag(n)*delZToxUnit(n)
	        END DO

			! Update the Omega
			prevOmegaTot(:,:,j) = newOmega

			bdBndIdx = 0
			! Update the tracker for whether my toxins is overall bound or not
			IF (ANY(SUM(newOmega,2)==1)) THEN !newOmega(beads,bindingSites)
	            ! sum along binding sites to get a 1 for every bead that is bound and a 0 for every polymer bead that is unbound and then count the bead as bound if any of its binding sites is bound (ie, any of the bead slots = 1)
	            timeBoundUnbound(j) = 1 !timeBoundUnbound(ntox)
	            timeTypeBond(j) = SUM(SUM(newOmega,2),1)
	            bdBndIdx = FIND(SUM(newOmega,2)==1);
	            DO n = 1,bindingSites
	            	IF (bdBndIdx(n)>0) THEN
	            		idxToxBnd(bdBndIdx(n)) = j
	            	END IF
	            END DO

	        ELSE
	            timeBoundUnbound(j) = 0
	            timeTypeBond(j) = 0
	        END IF

	        ! Add together the forces on the toxins
	        dUdxTox(j) = SUM(FbindMag*delXToxUnit,1) + SUM(fljxTox,1) - fljxToxTox(j) !(maxTox)
	        dUdyTox(j) = SUM(FbindMag*delYToxUnit,1) + SUM(fljyTox,1) - fljyToxTox(j)
	        dUdzTox(j) = SUM(FbindMag*delZToxUnit,1) + SUM(fljzTox,1) - fljzToxTox(j)

	        ! Add together all of the lj forces from each toxin on the inhibitor beads
	        fljxToxTot = fljxToxTot+fljxTox
	        fljyToxTot = fljyToxTot+fljyTox
	        fljzToxTot = fljzToxTot+fljzTox
	    END DO
	ELSE
		DO j = 1,nTox
			dUdxTox(j) = - fljxToxTox(j) !(maxTox)
	        dUdyTox(j) = - fljyToxTox(j)
	        dUdzTox(j) = - fljzToxTox(j)
	    END DO
	END IF

	IF (beads > 0) THEN
		! Add together the forces on the polymer
		dUdx = -(fspringx + fwx + fljx + fbindx + fljxToxTot)
		dUdy = -(fspringy + fwy + fljy + fbindy + fljyToxTot)
		dUdz = -(fspringz + fwz + fljz + fbindz + fljzToxTot)

		! Calculate the next position of the polymer
		!print *, "rx = ", rx
		rx = rx + (Vf-dUdx/zeta)*deltaT + SQRT(2.0_DP*D)*deltaWt(1,:)
		ry = ry + (Vf-dUdy/zeta)*deltaT + SQRT(2.0_DP*D)*deltaWt(2,:)
		rz = rz + (Vf-dUdz/zeta)*deltaT + SQRT(2.0_DP*D)*deltaWt(3,:)
		!print *, "rx = ", rx
	END IF

	! Calculate the next position of the toxin
	!print *, "posTox = ", posTox
	posTox(1,:) = posTox(1,:) + (Vf-dUdxTox/zeta)*deltaT + SQRT(2.0_DP*DTox)*deltaWtTox(1,:)
	posTox(2,:) = posTox(2,:) + (Vf-dUdyTox/zeta)*deltaT + SQRT(2.0_DP*DTox)*deltaWtTox(2,:)
	posTox(3,:) = posTox(3,:) + (Vf-dUdzTox/zeta)*deltaT + SQRT(2.0_DP*DTox)*deltaWtTox(3,:)
	!print *, "posTox = ", posTox


	! Reestablish the boundary conditions for the polymer
    DO j = 1,beads

        IF (rx(j) > 0.5*sizeBox) THEN
            rx(j) = rx(j) - sizeBox
        ELSE IF (rx(j) < -0.5*sizeBox) THEN
            rx(j) = rx(j) + sizeBox
        END IF

        IF (ry(j) > 0.5*sizeBox) THEN
            ry(j) = ry(j) - sizeBox
        ELSE IF (ry(j) < -0.5*sizeBox) THEN
            ry(j) = ry(j) + sizeBox
        END IF

        IF (rz(j) > 0.5*sizeBox) THEN
            rz(j) = rz(j) - sizeBox
        ELSE IF (rz(j) < -0.5*sizeBox) THEN
            rz(j) = rz(j) + sizeBox
        END IF

    END DO
    !print *, "rx = ", rx

    ! Reestablish the boundary conditions for the toxin
    DO j = 1,nTox

        IF (posTox(1,j) > 0.5*sizeBox) THEN
            posTox(1,j) = posTox(1,j) - sizeBox
        ELSE IF (posTox(1,j) < -0.5*sizeBox) THEN
            posTox(1,j) = posTox(1,j) + sizeBox
        END IF

        IF (posTox(2,j) > 0.5*sizeBox) THEN
            posTox(2,j) = posTox(2,j) - sizeBox
        ELSE IF (posTox(2,j) < -0.5*sizeBox) THEN
            posTox(2,j) = posTox(2,j) + sizeBox
        END IF

        IF (posTox(3,j) > 0.5*sizeBox) THEN
            posTox(3,j) = posTox(3,j) - sizeBox
        ELSE IF (posTox(3,j) < -0.5*sizeBox) THEN
            posTox(3,j) = posTox(3,j) + sizeBox
        END IF

    END DO

    ! Save things to files at specified intervals
	IF (mod(i,checkBindingInterval) == 0) THEN

        OPEN (UNIT=10,FILE="dataOut/idxToxBnd.dat"&
             &,STATUS="unknown",POSITION="append")
        !DO j=1,maxTox
	    	WRITE (UNIT=10,FMT=*) idxToxBnd
		!END DO
        CLOSE (10)

        OPEN (UNIT=10,FILE="dataOut/numBnd.dat"&
             &,STATUS="unknown",POSITION="append")
        !DO j=1,maxTox
	    	WRITE (UNIT=10,FMT=*) SUM(timeBoundUnbound)
		!END DO
        CLOSE (10)


    END IF


	IF (mod(i,d_wrt) == 0) THEN

		! Write out the toxin position
		OPEN (UNIT=20,FILE="dataOut/posToxX.dat"&
             &,STATUS="unknown",POSITION="append")
        WRITE (UNIT=20,FMT=*) posTox(1,:)
        CLOSE (20)

        OPEN (UNIT=20,FILE="dataOut/posToxY.dat"&
             &,STATUS="unknown",POSITION="append")
        WRITE (UNIT=20,FMT=*) posTox(2,:)
        CLOSE (20)

        OPEN (UNIT=20,FILE="dataOut/posToxZ.dat"&
             &,STATUS="unknown",POSITION="append")
        WRITE (UNIT=20,FMT=*) posTox(3,:)
        CLOSE (20)

        ! Write out the polymer position
        OPEN (UNIT=20,FILE="dataOut/confX.dat"&
             &,STATUS="unknown",POSITION="append")
        WRITE (UNIT=20,FMT=*) rx
        CLOSE (20)

        OPEN (UNIT=20,FILE="dataOut/confY.dat"&
             &,STATUS="unknown",POSITION="append")
        WRITE (UNIT=20,FMT=*) ry
        CLOSE (20)

        OPEN (UNIT=20,FILE="dataOut/confZ.dat"&
             &,STATUS="unknown",POSITION="append")
        WRITE (UNIT=20,FMT=*) rz
        CLOSE (20)

    	!print *, "made it out of writing out files"

	END IF
     
    IF (mod(i,100000) == 0) THEN

        OPEN (UNIT=20,FILE="dataOut/counter.dat",STATUS="unknown")
        REWIND (20)
        WRITE (20,FMT='(I8)') i
        CLOSE (20)

        OPEN (UNIT=20,FILE=restart_conf,STATUS="unknown")
        REWIND (UNIT=20)
        DO j=1,beads
           WRITE (UNIT=20,FMT=*) rx(j),ry(j),rz(j)
        END DO
        CLOSE (20)

    END IF

	IF (mod(i,1000000) == 0) THEN
		print *, "Step ",i,"of ",totSteps
	END IF

	END DO


! Deallocate all of our variables, not sure why we do this
DEALLOCATE (rx,ry,rz,posTox)
DEALLOCATE (fljx,fljy,fljz,fljxTox,fljyTox,fljzTox)
DEALLOCATE (fljxToxTot,fljyToxTot,fljzToxTot,fljxToxTox,fljyToxTox,fljzToxTox)
DEALLOCATE (fspringx,fspringy,fspringz)
DEALLOCATE (fspringxTox,fspringyTox,fspringzTox)
DEALLOCATE (fbindx,fbindy,fbindz,FbindMag)
DEALLOCATE (prevOmegaTot,prevOmega,newOmega)
DEALLOCATE (timeBoundUnbound)
DEALLOCATE (dRTox,delXToxUnit,delYToxUnit,delZToxUnit)
DEALLOCATE (deltaWt,deltaWtTox,beadsBound)
DEALLOCATE (dUdx,dUdy,dUdz,dUdxTox,dUdyTox,dUdzTox)
DEALLOCATE (nNghbrsPP,nghbrListPP)
DEALLOCATE (fwx,fwy,fwz)

print *, "Done!"
  
END program toxinSolubilityNVT
