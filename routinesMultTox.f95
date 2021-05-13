MODULE routinesMultTox

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutines used in toxinSolubilityNVT.f95
!
! Emiko Zumbro
! ezumbro@mit.edu
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

CONTAINS

  SUBROUTINE worm_like(fwx,fwy,fwz,rx,ry,rz)

  !This subroutine uses a DOT product interaction between adjacent bond vectors
  !to implement the WLC model (from Alfredo, added 07/26/19)

  USE parameters

  IMPLICIT NONE

  INTEGER :: i,n

  REAL(DP) :: delta,constw
  REAL(DP) :: X_kj,Y_kj,Z_kj,X_ij,Y_ij,Z_ij
  REAL(DP) :: X_kj_Adj,Y_kj_Adj,Z_kj_Adj,X_ij_Adj,Y_ij_Adj,Z_ij_Adj
  REAL(DP), DIMENSION(:), INTENT(IN):: rx,ry,rz
  REAL(DP), DIMENSION(:), INTENT(OUT) :: fwx,fwy,fwz

  fwx=0.0
  fwy=0.0
  fwz=0.0

  constw=gammaw/(sigma**2)

  IF (Nmer == 1) THEN
    print *, fwx
    print *, "forgot to change isPolymer to 0"

  ELSE IF (Nmer == beads) THEN 

    DO i=2,SIZE(rx)-1


    X_kj = rx(i+1)-rx(i)
    Y_kj = ry(i+1)-ry(i)
    Z_kj = rz(i+1)-rz(i)

    X_ij = rx(i-1)-rx(i)
    Y_ij = ry(i-1)-ry(i)
    Z_ij = rz(i-1)-rz(i)

    ! Correct for periodic boundary conditions
    X_kj_Adj = X_kj - sizeBox*NINT(X_kj/sizeBox)
    Y_kj_Adj = Y_kj - sizeBox*NINT(Y_kj/sizeBox)
    Z_kj_Adj = Z_kj - sizeBox*NINT(Z_kj/sizeBox)

    X_ij_Adj = X_ij - sizeBox*NINT(X_ij/sizeBox)
    Y_ij_Adj = Y_ij - sizeBox*NINT(Y_ij/sizeBox)
    Z_ij_Adj = Z_ij - sizeBox*NINT(Z_ij/sizeBox)

    fwx(i) = fwx(i)+ X_kj_Adj + X_ij_Adj !rx(i+1)+rx(i-1)-2.0*rx(i) !x32-x12
    fwy(i) = fwy(i)+ Y_kj_Adj + Y_ij_Adj !ry(i+1)+ry(i-1)-2.0*ry(i) ! (derivative of r32.r12 wrt r2?)
    fwz(i) = fwz(i)+ Z_kj_Adj + Z_ij_Adj !rz(i+1)+rz(i-1)-2.0*rz(i)

    fwx(i-1) = fwx(i-1)-X_kj_Adj !(rx(i+1)-rx(i)) ! -x32
    fwy(i-1) = fwy(i-1)-Y_kj_Adj !(ry(i+1)-ry(i)) ! (derivative of r32.r12 wrt r1?)
    fwz(i-1) = fwz(i-1)-Z_kj_Adj !(rz(i+1)-rz(i))

    fwx(i+1) = -X_ij_Adj !-(rx(i-1)-rx(i)) ! -x12 (derivative of r32.r12 wrt r3?)
    fwy(i+1) = -Y_ij_Adj !-(ry(i-1)-ry(i))
    fwz(i+1) = -Z_ij_Adj !-(rz(i-1)-rz(i))

    END DO

  ELSE IF (Nmer < beads) THEN

    DO n = 1,beads/Nmer
      DO i=Nmer*(n-1)+2,Nmer*n-1

        X_kj = rx(i+1)-rx(i)
        Y_kj = ry(i+1)-ry(i)
        Z_kj = rz(i+1)-rz(i)

        X_ij = rx(i-1)-rx(i)
        Y_ij = ry(i-1)-ry(i)
        Z_ij = rz(i-1)-rz(i)

        ! Correct for periodic boundary conditions
        X_kj_Adj = X_kj - sizeBox*NINT(X_kj/sizeBox)
        Y_kj_Adj = Y_kj - sizeBox*NINT(Y_kj/sizeBox)
        Z_kj_Adj = Z_kj - sizeBox*NINT(Z_kj/sizeBox)

        X_ij_Adj = X_ij - sizeBox*NINT(X_ij/sizeBox)
        Y_ij_Adj = Y_ij - sizeBox*NINT(Y_ij/sizeBox)
        Z_ij_Adj = Z_ij - sizeBox*NINT(Z_ij/sizeBox)

        fwx(i) = fwx(i)+ X_kj_Adj + X_ij_Adj !rx(i+1)+rx(i-1)-2.0*rx(i) !x32-x12
        fwy(i) = fwy(i)+ Y_kj_Adj + Y_ij_Adj !ry(i+1)+ry(i-1)-2.0*ry(i) ! (derivative of r32.r12 wrt r2?)
        fwz(i) = fwz(i)+ Z_kj_Adj + Z_ij_Adj !rz(i+1)+rz(i-1)-2.0*rz(i)

        fwx(i-1) = fwx(i-1)-X_kj_Adj !(rx(i+1)-rx(i)) ! -x32
        fwy(i-1) = fwy(i-1)-Y_kj_Adj !(ry(i+1)-ry(i)) ! (derivative of r32.r12 wrt r1?)
        fwz(i-1) = fwz(i-1)-Z_kj_Adj !(rz(i+1)-rz(i))

        fwx(i+1) = -X_ij_Adj !-(rx(i-1)-rx(i)) ! -x12 (derivative of r32.r12 wrt r3?)
        fwy(i+1) = -Y_ij_Adj !-(ry(i-1)-ry(i))
        fwz(i+1) = -Z_ij_Adj !-(rz(i-1)-rz(i))
        
      END DO

    END DO

  END IF


  fwx=constw*fwx
  fwy=constw*fwy
  fwz=constw*fwz


  END SUBROUTINE worm_like

  SUBROUTINE updateNeighborListPP(rx,ry,rz,nNghbrs,nghbrList)

    USE parameters

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: rx,ry,rz
    INTEGER, DIMENSION(:), INTENT(OUT) :: nNghbrs(beads)
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: nghbrList(beads,beads)

    REAL(DP) :: deltaX, deltaY, deltaZ, deltaXAdj, deltaYAdj, deltaZAdj
    REAL(DP) :: delta, cutoff
    ! REAL(DP), DIMENSION(:,:) :: deltaTot
    INTEGER :: i,j

    nNghbrs = 0
    nghbrList = 0
    cutoff = (LJcutoff+1.0_DP)*sigma

    DO i=1,SIZE(rx)
       DO j=i+1,SIZE(rx)
          
          deltaX = rx(j)-rx(i)
          deltaY = ry(j)-ry(i)
          deltaZ = rz(j)-rz(i)

          ! Correct for periodic boundary conditions
          deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
          deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
          deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

          delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5

          ! if you are within the cutoff distance you are a neighbor
          IF (delta < cutoff) THEN
          ! add one to the number of neighbors beads i and j have
            nNghbrs(i) = nNghbrs(i)+1
            nNghbrs(j) = nNghbrs(j)+1
          ! then add bead i to j's neighbor list and vice versa
            nghbrList(i,nNghbrs(i)) = j
            nghbrList(j,nNghbrs(j)) = i
          END IF

      END DO
    END DO

  END SUBROUTINE updateNeighborListPP

  SUBROUTINE updateNeighborListToxTox(rx,ry,rz,nNghbrs,nghbrList)

    USE parameters

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: rx,ry,rz
    INTEGER, DIMENSION(:), INTENT(OUT) :: nNghbrs(maxTox)
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: nghbrList(maxTox,maxTox)

    REAL(DP) :: deltaX, deltaY, deltaZ, deltaXAdj, deltaYAdj, deltaZAdj
    REAL(DP) :: delta, cutoff
    INTEGER :: i,j

    nNghbrs = 0
    nghbrList = 0
    cutoff = (LJcutoff+1.0_DP)*sigmaToxTox

    DO i=1,SIZE(rx)
       DO j=i+1,SIZE(rx)
          
          deltaX = rx(j)-rx(i)
          deltaY = ry(j)-ry(i)
          deltaZ = rz(j)-rz(i)

          ! Correct for periodic boundary conditions
          deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
          deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
          deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

          delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5

          ! if you are within the cutoff distance you are a neighbor
          IF (delta < cutoff) THEN
          ! add one to the number of neighbors beads i and j have
            nNghbrs(i) = nNghbrs(i)+1
            nNghbrs(j) = nNghbrs(j)+1
          ! then add bead i to j's neighbor list and vice versa
            nghbrList(i,nNghbrs(i)) = j
            nghbrList(j,nNghbrs(j)) = i
          END IF

      END DO
    END DO

  END SUBROUTINE updateNeighborListToxTox

  SUBROUTINE updateNeighborListPolyTox(rx,ry,rz,tx,ty,tz,nNghbrs,nghbrList)

  USE parameters

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: rx,ry,rz,tx,ty,tz
    INTEGER, DIMENSION(:), INTENT(OUT) :: nNghbrs(maxTox)
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: nghbrList(maxTox,beads)

    REAL(DP) :: deltaX, deltaY, deltaZ, deltaXAdj, deltaYAdj, deltaZAdj
    REAL(DP) :: delta, cutoff
    INTEGER :: i,j

    nNghbrs = 0
    nghbrList = 0
    cutoff = (LJcutoff+1.0_DP)*sigmaToxPoly

    DO i=1,SIZE(tx)
       DO j=1,SIZE(rx)
          
          deltaX = rx(j)-tx(i)
          deltaY = ry(j)-ty(i)
          deltaZ = rz(j)-tz(i)

          ! Correct for periodic boundary conditions
          deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
          deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
          deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

          delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5

          ! if you are within the cutoff distance you are a neighbor
          IF (delta < cutoff) THEN
          ! add one to the number of neighbors toxin i has
            nNghbrs(i) = nNghbrs(i)+1
          ! then add bead j to i's neighbor list
            nghbrList(i,nNghbrs(i)) = j
          END IF

      END DO
    END DO

  END SUBROUTINE updateNeighborListPolyTox

  SUBROUTINE rljmodNLBC(fljx,fljy,fljz,rx,ry,rz,nNghbrs,nghbrList)
    ! Takes into account boundary conditions and neighbor list

    USE parameters

    IMPLICIT NONE

    ! Declare constants
    INTEGER :: i,j

    REAL(DP) :: delta,const
    REAL(DP) :: deltaX, deltaY, deltaZ, deltaXAdj, deltaYAdj, deltaZAdj
    REAL(DP), DIMENSION(:), INTENT(IN) :: rx,ry,rz
    INTEGER, DIMENSION(:), INTENT(IN) :: nNghbrs
    INTEGER, DIMENSION(:,:), INTENT(IN) :: nghbrList
    REAL(DP), DIMENSION(:), INTENT(OUT) :: fljx,fljy,fljz
    REAL(DP) :: fx,fy,fz

    fx=0.0
    fy=0.0
    fz=0.0

    fljx=0.0
    fljy=0.0
    fljz=0.0

    ! Loop through all the bead pairs and see what the LJ potential is between them
    DO i=1,SIZE(nNghbrs)
       DO j=1,nNghbrs(i)
          
          deltaX = rx(nghbrList(i,j))-rx(i)
          deltaY = ry(nghbrList(i,j))-ry(i)
          deltaZ = rz(nghbrList(i,j))-rz(i)

          ! Correct for periodic boundary conditions
          deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
          deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
          deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

          delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5

          IF (delta < LJcutoff*sigma) THEN

            const = ((sigma/delta)**12-(sigma/delta)**6)&
                 &*(epsilon/delta)/delta

            fx = - const*deltaXAdj 
            fy = - const*deltaYAdj
            fz = - const*deltaZAdj 
          
          ELSE

             fx=0.0
             fy=0.0
             fz=0.0

          END IF

          fljx(i) = fljx(i)+fx
          fljy(i) = fljy(i)+fy
          fljz(i) = fljz(i)+fz

       END DO
    END DO
       
  END SUBROUTINE rljmodNLBC

  SUBROUTINE rljToxToxNLBC(fljx,fljy,fljz,rx,ry,rz,nNghbrs,nghbrList)
    ! Takes into account boundary conditions and neighbor list
    ! Does not work with dynamically adjusting toxin concentration

    USE parameters

    IMPLICIT NONE

    ! Declare constants
    INTEGER :: i,j

    REAL(DP) :: delta,const
    REAL(DP) :: deltaX, deltaY, deltaZ, deltaXAdj, deltaYAdj, deltaZAdj
    REAL(DP), DIMENSION(:), INTENT(IN) :: rx,ry,rz
    INTEGER, DIMENSION(:), INTENT(IN) :: nNghbrs
    INTEGER, DIMENSION(:,:), INTENT(IN) :: nghbrList
    REAL(DP), DIMENSION(:), INTENT(OUT) :: fljx,fljy,fljz
    REAL(DP) :: fx,fy,fz

    fx=0.0
    fy=0.0
    fz=0.0

    fljx=0.0
    fljy=0.0
    fljz=0.0

    ! Loop through all the bead pairs and see what the LJ potential is between them
    DO i=1,SIZE(nNghbrs)

       DO j=1,nNghbrs(i)

          
          deltaX = rx(nghbrList(i,j))-rx(i)
          deltaY = ry(nghbrList(i,j))-ry(i)
          deltaZ = rz(nghbrList(i,j))-rz(i)

          ! Correct for periodic boundary conditions
          deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
          deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
          deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

          delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5

          IF (delta < LJcutoff*sigmaToxTox) THEN

            const = ((sigmaToxTox/delta)**12-(sigmaToxTox/delta)**6)&
                 &*(epsToxTox/delta)/delta

            fx = - const*deltaXAdj 
            fy = - const*deltaYAdj
            fz = - const*deltaZAdj
          
          ELSE

             fx=0.0
             fy=0.0
             fz=0.0
             ! print *, "outside Bounds made everything 0"

          END IF

          fljx(i) = fljx(i)+fx
          fljy(i) = fljy(i)+fy
          fljz(i) = fljz(i)+fz


       END DO
    END DO
       
  END SUBROUTINE rljToxToxNLBC


  SUBROUTINE rljToxNLBC(fljx,fljy,fljz,deltaTox,rx,ry,rz,posTox,delXToxUnit,delYToxUnit,delZToxUnit,nNghbrs,nghbrList)

    USE parameters

    IMPLICIT NONE

    ! Declare constants
    INTEGER :: i,j
    INTEGER, INTENT(IN) :: nNghbrs
    INTEGER, DIMENSION(:), INTENT(IN) :: nghbrList

    REAL(DP) :: delta,const
    REAL(DP) :: deltaX, deltaY, deltaZ, deltaXAdj, deltaYAdj, deltaZAdj
    REAL(DP), DIMENSION(:), INTENT(IN) :: rx,ry,rz,posTox
    REAL(DP), DIMENSION(:), INTENT(OUT) :: fljx,fljy,fljz,deltaTox(beads)
    REAL(DP), DIMENSION(:), INTENT(OUT) :: delXToxUnit(beads), delYToxUnit(beads), delZToxUnit(beads)
    REAL(DP) :: fx,fy,fz

    fx=0.0
    fy=0.0
    fz=0.0

    deltaTox=0.0

    delXToxUnit=0.0
    delYToxUnit=0.0
    delZToxUnit=0.0

    fljx=0.0
    fljy=0.0
    fljz=0.0

    ! Loop through all the bead pairs and see what the LJ potential is between them
    DO i=1,nNghbrs

      deltaX = rx(nghbrList(i))-posTox(1)
      deltaY = ry(nghbrList(i))-posTox(2)
      deltaZ = rz(nghbrList(i))-posTox(3)

      ! Correct for periodic boundary conditions
      deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
      deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
      deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

      delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5
      deltaTox(nghbrList(i)) = delta

      ! Make unit vectors for later use in bigger function
      delXToxUnit(nghbrList(i)) = deltaXAdj/delta
      delYToxUnit(nghbrList(i)) = deltaYAdj/delta
      delZToxUnit(nghbrList(i)) = deltaZAdj/delta

      IF (delta < LJcutoff*sigmaToxPoly) THEN

         const = ((sigmaToxPoly/delta)**12-(sigmaToxPoly/delta)**6)&
              &*(epsToxPoly/delta)/delta

         ! check to make sure the negative here is right!
         fx = - const*deltaXAdj 
         fy = - const*deltaYAdj
         fz = - const*deltaZAdj 
      
      ELSE

         fx=0.0
         fy=0.0
         fz=0.0

      END IF
      ! This part still needs work
      fljx(nghbrList(i)) = fljx(nghbrList(i))-fx
      fljy(nghbrList(i)) = fljy(nghbrList(i))-fy
      fljz(nghbrList(i)) = fljz(nghbrList(i))-fz

    END DO
  END SUBROUTINE rljToxNLBC


  SUBROUTINE springForces(fspringx,fspringy,fspringz,rx,ry,rz)

    USE parameters

    IMPLICIT NONE

    INTEGER :: i

    REAL(DP) :: delta,const
    REAL(DP) :: deltaX, deltaY, deltaZ, deltaXAdj, deltaYAdj, deltaZAdj
    REAL(DP), DIMENSION(:), INTENT(IN) :: rx,ry,rz
    REAL(DP), DIMENSION(:), INTENT(OUT) :: fspringx,fspringy,fspringz
    REAL(DP) :: fsxMid (beads), fsyMid (beads), fszMid (beads)
    REAL(DP) :: fsx (beads), fsy (beads), fsz (beads)
 
    fspringx = 0.0
    fspringy = 0.0
    fspringz = 0.0

    !print *, fspringx

    ! Compute spring forces from being connected to neighboring beads
    ! Compute End Forces (because only connected to one neighbor)
    IF (Nmer == 1) THEN
      print *, fspringx
      print *, "forgot to change isPolymer to 0"

    ELSE IF (Nmer == beads) THEN 

      ! Do first bead
      i=1
      
      deltaX = rx(i+1)-rx(i)
      deltaY = ry(i+1)-ry(i)
      deltaZ = rz(i+1)-rz(i)

      ! Correct for periodic boundary conditions
      deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
      deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
      deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

      delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5
      const = k_sp*(delta-l_sp)/delta

      fspringx(i) = const*(deltaXAdj)
      fspringy(i) = const*(deltaYAdj)
      fspringz(i) = const*(deltaZAdj)

      ! Do force on the last bead
      i=beads-1
      
      deltaX = rx(i+1)-rx(i)
      deltaY = ry(i+1)-ry(i)
      deltaZ = rz(i+1)-rz(i)

      ! Correct for periodic boundary conditions
      deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
      deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
      deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

      delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5
      const = k_sp*(delta-l_sp)/(delta)

      fspringx(i+1) = - const*(deltaXAdj)
      fspringy(i+1) = - const*(deltaYAdj)
      fspringz(i+1) = - const*(deltaZAdj)
        
    ELSE IF (Nmer < beads) THEN

      ! Calculate bonding of first beads
      DO i=1,beads-1,Nmer

        deltaX = rx(i+1)-rx(i)
        deltaY = ry(i+1)-ry(i)
        deltaZ = rz(i+1)-rz(i)

        ! Correct for periodic boundary conditions
        deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
        deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
        deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

        delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5
        const = k_sp*(delta-l_sp)/(delta)

        fspringx(i) = const*(deltaXAdj)
        fspringy(i) = const*(deltaYAdj)
        fspringz(i) = const*(deltaZAdj)

      END DO

      ! Calculate bonding of last beads
      DO i=(Nmer-1),(beads-1),Nmer

        deltaX = rx(i+1)-rx(i)
        deltaY = ry(i+1)-ry(i)
        deltaZ = rz(i+1)-rz(i)

        ! Correct for periodic boundary conditions
        deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
        deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
        deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

        delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5
        const = k_sp*(delta-l_sp)/(delta)

        fspringx(i+1) = - const*(deltaXAdj)
        fspringy(i+1) = - const*(deltaYAdj)
        fspringz(i+1) = - const*(deltaZAdj)

      END DO
        
    END IF
    
    
    ! Compute middle bead connectivity forces
    IF (Nmer > 2) THEN
      DO i=2,beads-1

        deltaX = rx(i+1)-rx(i)
        deltaY = ry(i+1)-ry(i)
        deltaZ = rz(i+1)-rz(i)

        ! Correct for periodic boundary conditions
        deltaXAdj = deltaX - sizeBox*NINT(deltaX/sizeBox)
        deltaYAdj = deltaY - sizeBox*NINT(deltaY/sizeBox)
        deltaZAdj = deltaZ - sizeBox*NINT(deltaZ/sizeBox)

        delta = (deltaXAdj**2+deltaYAdj**2+deltaZAdj**2)**0.5
        const = k_sp*(delta-l_sp)/(delta)

        ! If you are the end of an Nmer, make the force to the right 0
        ! Otherwise proceed normally
        IF (MOD(i,Nmer) /= 0) THEN
          fsxMid(i) = const*(deltaXAdj)
          fsyMid(i) = const*(deltaYAdj)
          fszMid(i) = const*(deltaZAdj)
        ELSE
          fsxMid(i) = 0.0
          fsyMid(i) = 0.0
          fszMid(i) = 0.0
        END IF

      END DO
      
      ! Add in the forces pulling back to the left
      ! This is 0 force pulling back to the left for left most end beads

        fsxMid(1) = fspringx(1)
        fsyMid(1) = fspringy(1)
        fszMid(1) = fspringz(1)

        fsxMid(beads) = fspringx(beads)
        fsyMid(beads) = fspringy(beads)
        fszMid(beads) = fspringz(beads)

        ! Reassign the middle forces to a temporary variable
        fsx = fsxMid
        fsy = fsyMid
        fsz = fszMid

      DO i=2,beads-1
        ! Without the temp variable, when I did -fsxMid(i-1) errors start to accummulate as you move right along the chain because you modify i-1 and then you use the modified value instead of the original value
       fsxMid(i) = fsxMid(i) - fsx(i-1)
       fsyMid(i) = fsyMid(i) - fsy(i-1)
       fszMid(i) = fszMid(i) - fsz(i-1)

      END DO

      DO i=1,beads,Nmer

        fsxMid(i) = fspringx(i)
        fsyMid(i) = fspringy(i)
        fszMid(i) = fspringz(i)

      END DO

      DO i=Nmer,beads,Nmer

        fsxMid(i) = fspringx(i)
        fsyMid(i) = fspringy(i)
        fszMid(i) = fspringz(i)

      END DO

      fspringx = fsxMid
      fspringy = fsyMid
      fspringz = fszMid

    END IF

  END SUBROUTINE springForces

! dR is size beads by 1 (cycle through this funciton for each toxin)
! prevOmegaTot is (beads,bindingSites,maxTox)
! prevOmega = (beads,bindingSites,existingToxins(k))
! newOmega = (beads,bindingSites)
  SUBROUTINE bound(prevOmegaTot, prevOmega, dRTox, FbindMag, newOmega,delE_0,delE_UB)

    USE parameters
    USE functions_Emi

    IMPLICIT NONE

    INTEGER :: i,j,k,cnt,beadToBind,beadBound,randIdx
    REAL(DP) :: deltaX, deltaY, deltaZ, deltaXAdj, deltaYAdj, deltaZAdj

    REAL(DP), DIMENSION(:), INTENT(IN) :: dRTox(beads),delE_0(beads),delE_UB(beads)
    INTEGER, DIMENSION(:,:), INTENT(IN) :: prevOmega(beads,bindingSites) ! how this toxin was previously bound
    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: prevOmegaTot(beads,bindingSites,maxTox) ! how all the other toxins are bound (to prevent binding to same bead)

    REAL(DP), DIMENSION(:), INTENT(OUT) :: FbindMag(beads)
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: newOmega(beads,bindingSites)

    REAL(DP), DIMENSION(:,:) :: oper(beads,bindingSites), alreadyBoundTemp(beads,bindingSites)
    INTEGER, DIMENSION(:) :: alreadyBound(beads), beadsBound(beads)
    INTEGER, DIMENSION(:), ALLOCATABLE :: alreadyBoundIdx, closeBds, possibleBinds
    REAL(DP), DIMENSION(:), ALLOCATABLE :: randPerm
    REAL(DP), DIMENSION(:) :: dR_modifier(beads), dRAvailableBeads(beads), delE_B(beads)
    LOGICAL :: mask(beads)

    delE_B = delE_0 + delE_UB
    newOmega = 0
    oper = 0.0
    dR_modifier = 0.0
    
    ! Generate uniform random numbers so that you can use them later to determine binding
    CALL RANDOM_NUMBER(oper)

    ! Keep track of the beads you've already bound so that you don't get more than one tox binding site bound to it
    alreadyBoundTemp = SUM(prevOmegaTot, 3)
    alreadyBound = SUM(alreadyBoundTemp, 2) ! This isn't affected by changing the storage of toxins because unbound and nonexistent toxins are both 0

    ! Check to make sure nothing has bound to the same bead twice
    mask = alreadyBound > 1
    IF (ANY(mask)) THEN
      print *, "Error: bound to the same bead twice"
    END IF
    
    ! Get the indexes of the beads that are bound so that you can keep track of them
    mask = alreadyBound==1
    alreadyBoundIdx = FIND(mask)

    ! Make the already bound beads too far to be picked up by the close beads search (so that they won't be available for new bonds to form)
    dR_modifier(alreadyBoundIdx) = l_bind*(1.0_DP+reach)*2.0_DP; ! 2 is just to guarantee bead is well outside reach of binding site
    dRAvailableBeads = dRTox + dR_modifier;

    ! Iterate through each binding site and assess whether you bind/unbind stay unbound to it
    DO j=1,bindingSites
      ! If this binding site is not already bound to something
      mask = prevOmega(:,j) == 0
      IF (ALL(mask)) THEN
        ! If your site is within binding range of any of the beads
        ! This will have a problem if there is no excluded volume between toxins and polymers
        mask = ((dRAvailableBeads<=(reach)).AND.(dRAvailableBeads>0))
        IF (ANY(mask)) THEN
          ! Find which beads are close
          closeBds = FIND(mask)
          !print *, "closeBds = ", closeBds
          !print *, "looking for new potential bond"
          ! Check if any of the close beads have a high enough energy to jump the binding energy barrier
          DO i = 1, SIZE(closeBds)
            IF (oper(i,j) < EXP(-delE_B(closeBds(i)))) THEN
              newOmega(closeBds(i),j) = 1;
              ! Make sure bead you bound to is out of reach
              dRAvailableBeads(closeBds(i)) = dRTox(closeBds(i)) + l_bind*(1.0_DP+reach)*2.0_DP;
              !print *, "found new bond"
            END IF
          END DO
          ! Check if you can bind to more than one bead
          IF (SUM(newOmega(:,j)) > 1) THEN
              ! Only pick one to bind to randomly
              !randomize the index of beads you could successfully bind to
              ALLOCATE (randPerm(SUM(newOmega(:,j))))
              CALL RANDOM_NUMBER(randPerm)

              randIdx = MAXLOC(randPerm,1)
              mask = newOmega(:,j) == 1
              possibleBinds = FIND(mask)

              beadToBind = possibleBinds(randIdx)
              newOmega(:,j) = 0
              newOmega(beadToBind,j) = 1
              ! Make sure the bead you just bound to is out of reach of other binding sites
              dRAvailableBeads(closeBds) = dRTox(closeBds)
              dRAvailableBeads(beadToBind) = dRTox(beadToBind) + l_bind*(1.0_DP+reach)*2.0_DP


              ! Deallocate randperm in case it is used again for another binding site
              DEALLOCATE (randPerm)
          END IF
        ! Otherwise, you're not within reach of the beads and everything should be unbound
        ELSE
          newOmega(:,j) = 0
        END IF
        
      ! If you're already bound to something
      ELSE
        beadBound = MAXLOC(prevOmega(:,j),1)
        !print *,"beadBound = ", beadBound
        ! If you don't have enough energy to unbind then stay bound to the same bead
        IF (oper(beadBound,j) > exp(-delE_UB(beadBound))) THEN
          newOmega(:,j) = prevOmega(:,j)
          !print *, "rebind to same bead"
        ELSE
          newOmega(:,j) = 0
          !print *, "unbound from bead"
        END IF
      END IF

    END DO

    ! Calculate the new force
    beadsBound = SUM(newOmega,2)
    FbindMag = -k_bind*beadsBound*(dRTox-l_bind)

    ! print *, "newOmega = ", newOmega
  END SUBROUTINE bound

! This isn't very elegant, creates inhibitor starting positions in s or ladder shape, only goes up to 92 beads
  SUBROUTINE setInitialPos(rx,ry,rz)
    USE parameters
    USE functions_Emi

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(OUT) :: rx(beads),ry(beads),rz(beads)
    REAL(DP) :: rightEdge,leftEdge
    INTEGER :: lengthRow,i
    REAL(DP), DIMENSION(:) :: iPosTmp(100)
    INTEGER :: maxBeads = 64

    rightEdge = (floor(sizeBox/2.0)-1.0)
    leftEdge = -(floor(sizeBox/2.0)-1.0)
    lengthRow = rightEdge-leftEdge+1

    ! Set up an array of x positions for your inhibitor to start from
    DO i=1,lengthRow
      iPosTmp(i) = leftEdge+l_sp*(i-1.0)
    END DO
    DO i=lengthRow+1,lengthRow+4
      iPosTmp(i) = rightEdge
    END DO
    DO i=1,lengthRow
      iPosTmp(i+lengthRow+4) = rightEdge-l_sp*(i-1.0)
    END DO
    DO i=2*lengthRow+4,2*lengthRow+8
      iPosTmp(i) = leftEdge
    END DO
    DO i=1,lengthRow
      iPosTmp(i+2*lengthRow+8) = leftEdge+l_sp*(i-1.0)
    END DO
    DO i=3*lengthRow+8,3*lengthRow+12
      iPosTmp(i) = rightEdge
    END DO
    DO i=1,lengthRow
      iPosTmp(i+3*lengthRow+12) = rightEdge-l_sp*(i-1.0)
    END DO
    DO i=4*lengthRow+12,4*lengthRow+16
      iPosTmp(i) = leftEdge
    END DO

    rx = iPosTmp(1:beads)
    ! print *, iPosTmp
    ! print *, rx

    ! Set y positions
    iPosTmp = 0.0
    DO i=1,lengthRow
      iPosTmp(i) = leftEdge
    END DO
    DO i=lengthRow+1,lengthRow+4
      iPosTmp(i) = leftEdge+(i-lengthRow)
    END DO
    DO i=1,lengthRow
      iPosTmp(i+lengthRow+4) = leftEdge+5
    END DO
    DO i=2*lengthRow+4,2*lengthRow+8
      iPosTmp(i) = leftEdge+(i+1-2*lengthRow)
    END DO
    DO i=1,lengthRow
      iPosTmp(i+2*lengthRow+8) = leftEdge+10
    END DO
    DO i=3*lengthRow+8,3*lengthRow+12
      iPosTmp(i) = leftEdge+(i+2-3*lengthRow)
    END DO
    DO i=1,lengthRow
      iPosTmp(i+3*lengthRow+12) = leftEdge+15
    END DO
    DO i=4*lengthRow+12,4*lengthRow+15
      iPosTmp(i) = leftEdge+(i+3-4*lengthRow)
    END DO

    ry = iPosTmp(1:beads)
    ! print *, iPosTmp
    ! print *, ry

    ! Set z positions
    rz = -4.0 !0.0

    IF (beads > maxBeads) THEN
     DO i = 1,(beads-maxBeads)
        rx(i+maxBeads) = rx(i)
        ry(i+maxBeads) = ry(i)
        rz(i+maxBeads) = rz(i)-2.0
     END DO
    END IF

  END SUBROUTINE setInitialPos

  SUBROUTINE setInitialPosTox(rx,ry,rz,px,py,pz)

    USE parameters
    USE functions_Emi

    IMPLICIT NONE

    INTEGER :: iTox, otherTox, ibds
    REAL(DP) :: dX,dY,dZ,dXAdj,dYAdj,dZAdj,dToxTmp
    REAL(DP), DIMENSION(:), INTENT(OUT) :: rx(nTox),ry(nTox),rz(nTox)
    REAL(DP), DIMENSION(:), INTENT(IN) :: px(beads),py(beads),pz(beads)


    DO iTox=1,nTox

  100   CALL RANDOM_NUMBER(rx(iTox))
        CALL RANDOM_NUMBER(ry(iTox))
        CALL RANDOM_NUMBER(rz(iTox))
        ! print *, "iTox = ",iTox
      rx(iTox) = rx(iTox)*sizeBox-sizeBox/2.0_DP
      ry(iTox) = ry(iTox)*sizeBox-sizeBox/2.0_DP
      rz(iTox) = rz(iTox)*sizeBox-sizeBox/2.0_DP

      ! Check to make sure you didn't just put your new toxin on top of other ones
      DO otherTox=1,iTox-1
            !print *, "iTox = ",iTox

            dX = rx(otherTox)-rx(iTox)
            dY = ry(otherTox)-ry(iTox)
            dZ = rz(otherTox)-rz(iTox)

            ! Correct for periodic boundary conditions
            dXAdj = dX - sizeBox*NINT(dX/sizeBox)
            dYAdj = dY - sizeBox*NINT(dY/sizeBox)
            dZAdj = dZ - sizeBox*NINT(dZ/sizeBox)

            dToxTmp = SQRT(dXAdj**2.0_DP+dYAdj**2.0_DP+dZAdj**2.0_DP)

            IF (dToxTmp < 1.5_DP*reach) THEN
              GOTO 100
            END IF
      END DO
      DO ibds=1,beads

            dX = px(ibds)-rx(iTox)
            dY = py(ibds)-ry(iTox)
            dZ = pz(ibds)-rz(iTox)

            ! Correct for periodic boundary conditions
            dXAdj = dX - sizeBox*NINT(dX/sizeBox)
            dYAdj = dY - sizeBox*NINT(dY/sizeBox)
            dZAdj = dZ - sizeBox*NINT(dZ/sizeBox)

            dToxTmp = (dXAdj**2+dYAdj**2+dZAdj**2)**0.5
            
            IF (dToxTmp < 1.5_DP*reach) THEN
              GOTO 100
            END IF

        END DO
    END DO

  END SUBROUTINE setInitialPosTox 
   
  SUBROUTINE init_random_seed(nRun,seed)

      INTEGER, INTENT(IN) :: nRun
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: seed
    
      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))
    
      CALL SYSTEM_CLOCK(COUNT=clock)
    
      seed = clock + 37 * (/ (i - 1, i = 1, n) /) + nRun
      print *, "seed = ", seed
      CALL RANDOM_SEED(PUT = seed)
    
      !DEALLOCATE(seed)

  END SUBROUTINE

  SUBROUTINE findAggToxins(dRToxTox,timeBoundUnbound,nToxAgg,closeToxIdx)

    USE parameters

    IMPLICIT NONE

    INTEGER :: i,j
    INTEGER, DIMENSION(:,:) :: closeMaskINT(nTox,nTox)
    INTEGER, DIMENSION(:), INTENT(IN) :: timeBoundUnbound(nTox)
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: dRToxTox
    INTEGER, DIMENSION(:), INTENT(OUT) :: closeToxIdx(nTox)
    INTEGER, INTENT(OUT) :: nToxAgg

    closeMaskINT = 0
    nToxAgg = 0
    closeToxIdx = 0

    DO i=1,nTox
      DO j = i+1,nTox
        IF (AND(dRToxTox(i,j) > 0.0, dRToxTox(i,j) < reach)) THEN
          closeMaskINT(i,j) = 1
          closeMaskINT(j,i) = 1
        END IF
      END DO
    END DO

    closeToxIdx = SUM(closeMaskINT,2)
    closeToxIdx = closeToxIdx + timeBoundUnbound

    DO i=1,nTox
      IF (closeToxIdx(i)>0) THEN
        nToxAgg = nToxAgg+1
        closeToxIdx(i)=1
      END IF
    END DO

    !print *, "nToxAgg = ", nToxAgg

  END SUBROUTINE findAggToxins

  SUBROUTINE adjNTox(nToxAgg,toxConc0,nToxNew)

    USE parameters

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: nToxAgg, toxConc0
    INTEGER, INTENT(OUT) :: nToxNew
    REAL(DP) :: toxConc, v0
    LOGICAL :: didntAdd = .true.

    v0 = sizeBox**3
    nToxNew = nTox
    ! Find the current free toxin concentration
    toxConc = (nTox-nToxAgg)/(v0-nToxAgg*(4.0/3.0)*pi*(l_sp/2.0)**3)

    ! Change the number of free toxins
      DO WHILE (toxConc < toxConc0)
        nToxNew = nToxNew+1
        toxConc = (nToxNew-nToxAgg)/(v0-nToxAgg*(4.0/3.0)*pi*(l_sp/2.0)**3)
        didntAdd = .false.
      END DO


      DO WHILE (AND(toxConc > toxConc0,didntAdd))
        nToxNew = nToxNew-1
        toxConc = (nToxNew-nToxAgg)/(v0-nToxAgg*(4.0/3.0)*pi*(l_sp/2.0)**3)
      END DO

      ! Make sure youre always slightly below the original concentration instead of oscillating between below and above because of adding and subtracting methods
      ! so if you added until higher than conc0, subtract one
      IF (.NOT.didntAdd) THEN
        nToxNew = nToxNew-1
      END IF
      
      !print *, "nToxNew", nToxNew
  END SUBROUTINE adjNTox

SUBROUTINE squeeze_to_string_long(prefix,inum1,middle,inum2,suffix,ss)
! ***************************************************************************
! Needs:
!    subroutine integer_to_string(,)
! ***************************************************************************
    IMPLICIT NONE
    
    INTEGER inum1,inum2
    CHARACTER(*) :: prefix, middle, suffix
    CHARACTER(60) :: strnum1, strnum2
    CHARACTER(*), INTENT(out) :: ss
    
    CALL integer_to_string(inum1,strnum1)
    CALL integer_to_string(inum2,strnum2)
    
    ss = prefix // TRIM(strnum1) // middle // TRIM(strnum2) // suffix
    
    RETURN
  END SUBROUTINE squeeze_to_string_long


  SUBROUTINE squeeze_to_string(prefix,inum,suffix,ss)
! ***************************************************************************
! Needs:
!    subroutine integer_to_string(,)
! ***************************************************************************
    IMPLICIT NONE
    
    INTEGER inum
    CHARACTER(*) :: prefix, suffix
    CHARACTER(60) :: strnum
    CHARACTER(*), INTENT(out) :: ss
    
    CALL integer_to_string(inum,strnum)
    
    ss = prefix // TRIM(strnum) // suffix
    
    RETURN
  END SUBROUTINE squeeze_to_string
  
  
  SUBROUTINE integer_to_string(jj,ss)
! ******************************************************************
! 11/02/00
! Converts an integer to a character variable of same value
! Character variable has left over space so it must be trimmed.
! ******************************************************************
    IMPLICIT NONE
    
    INTEGER :: ii, jj
    INTEGER :: maxdigits, idig, icount
    DOUBLE PRECISION :: mod1, div1, frame
    CHARACTER(*), INTENT(out) :: ss
    
    frame=DBLE(jj)
    maxdigits=INT(LOG10(frame))+1
    ss=' '
    mod1=10.0
    div1=1.0
    
    DO ii=1,maxdigits
       idig=INT(MOD(frame,mod1)/div1)
       ss = ACHAR(idig+48) // ss
       div1=mod1
       mod1=10.0*mod1
    END DO
    
    RETURN
  END SUBROUTINE integer_to_string


SUBROUTINE ARRAY_COPY(src,dest,n_copied,n_not_copied)
  
  USE parameters

  REAL(DP), DIMENSION (:), INTENT(IN) :: src
  REAL(DP), DIMENSION (:), INTENT(OUT) :: dest
  INTEGER, INTENT(OUT) :: n_copied,n_not_copied

  n_copied = MIN(SIZE(src),SIZE(dest))
  n_not_copied = SIZE(src)-n_copied
  dest(1:n_copied) = src(1:n_copied)


END SUBROUTINE ARRAY_COPY



! ******************************************************************
! From Numerical Recipes in Fortran 90 2nd Edition
! 01/18/18
! creates a vector of gaussian distributed numbers with zero mean
! and variance of 1, using RANDOM_NUMBER as the source of uniform
! random numbers
! ******************************************************************
! This is not working, not sure why - implemented one that does it one number at a time
SUBROUTINE gasdev_v(harvest)

  USE parameters

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(OUT) :: harvest
  REAL(DP), DIMENSION(SIZE(harvest)) :: rsq,v1,v2
  REAL(DP), ALLOCATABLE, DIMENSION(:), SAVE :: g
  INTEGER :: n,ng,nn,m
  INTEGER, SAVE :: last_allocated=0
  LOGICAL, SAVE :: gaus_stored=.false.
  LOGICAL, DIMENSION(SIZE(harvest)) :: mask

  n=SIZE(harvest)

  print *, "ENTERING FUNCTION...."
  print *, "last_allocated = ", last_allocated
  print *, "gaus_stored = ", gaus_stored

  IF (n /= last_allocated) THEN
    IF (last_allocated /= 0) THEN
      DEALLOCATE(g)
    END IF
    ALLOCATE(g(n))
    last_allocated = n
    gaus_stored=.false.
  END IF

  IF (gaus_stored) THEN
    harvest=g
    gaus_stored=.false.
  ELSE
    ng = 1

    DO

      IF (ng > n) THEN
        EXIT
      END IF
      print *, "ng = ", ng

      CALL RANDOM_NUMBER(v1(ng:n))
      CALL RANDOM_NUMBER(v2(ng:n))
      ! Make your uniform random numbers from -1 to 1
      v1(ng:n) = 2.0_DP*v1(ng:n)-1.0_DP
      v2(ng:n) = 2.0_DP*v2(ng:n)-1.0_DP
      rsq(ng:n) = v1(ng:n)**2+v2(ng:n)**2
      print *, "rsq = ", rsq
      print *, "v1 before Array copy = ", v1
      print *, "v2 before array copy = ", v2

      mask(ng:n) = (rsq(ng:n)>0.0 .and. rsq(ng:n)<1.0)
      print *, "mask = ", mask
      print *, "v1(ng:) = ", v1(ng:)
      CALL ARRAY_COPY(PACK(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
      print *, "result of v1 PACK = ", PACK(v1(ng:n),mask(ng:n))
      print *, "v1 after array copy = ", v1
      v2(ng:ng+nn-1) = PACK(rsq(ng:n),mask(ng:n))
      print *, "v2 after pack = ", v2
      print *, "n copied = ", nn
      print *, "not copied = ", m
      ng=ng+nn

    END DO
    ! Now make the Box-Muller transformation to get two normal deviates
    ! Return the amount needed and save the rest for next time
    rsq=sqrt(-2.0_DP*log(rsq)/rsq)
    harvest=v1*rsq
    g=v2*rsq
    gaus_stored=.true.

  END IF

  print *, "LEAVING FUNCTION...."
  print *, "last_allocated = ", last_allocated
  print *, "gaus_stored = ", gaus_stored
  print *, "LEFT FUNCTION."

END SUBROUTINE gasdev_v

SUBROUTINE gasdev_s(harvest)

  USE parameters

  IMPLICIT NONE

  REAL(DP), INTENT(OUT) :: harvest
  REAL(DP) :: rsq,v1,v2
  REAL(DP), SAVE :: g
  LOGICAL, SAVE :: gaus_stored=.false.

  !print *, "made it into gasdev_s"

  IF (gaus_stored) THEN
  !print *, "made it into if"
    harvest=g
    gaus_stored=.false.
  ELSE
  !print *, "made it into else"
    DO
      CALL RANDOM_NUMBER(v1)
      CALL RANDOM_NUMBER(v2)

      !print *, "called random_number"
      v1 = 2.0_DP*v1-1.0_DP
      v2 = 2.0_DP*v2-1.0_DP
      rsq = v1**2+v2**2
      IF (rsq > 0.0 .and. rsq < 1.0) THEN
        EXIT
      END IF
      !print *, "random number stuff is happening"
    END DO
    !print *, "about to do some box muller stuff"
    ! Now make the Box-Muller transformation to get two normal deviates
    ! Return the amount needed and save the rest for next time
    rsq = sqrt(-2.0_DP*log(rsq)/rsq)
    harvest = v1*rsq
    g = v2*rsq
    gaus_stored=.true.
    !print *, "made it through the box muller"
  END IF

  !print *, "made it through gasdev_s"

END SUBROUTINE gasdev_s

SUBROUTINE pick_random(old_idxs, rand_idxs, nIdx)
  ! pick a random index
      IMPLICIT NONE
 
    INTEGER, DIMENSION(:), INTENT(IN) :: old_idxs
    INTEGER, DIMENSION(:), INTENT(OUT) :: rand_idxs
    INTEGER, DIMENSION(:), ALLOCATABLE :: new_idxs
    INTEGER, INTENT(IN) :: nIdx
    REAL :: r ! random number from 0 to 1
    INTEGER :: r_int,i ! random index

    ALLOCATE (new_idxs(SIZE(old_idxs)))
    new_idxs = old_idxs
    DO i = 1,nIdx

      ! get a random number between 0 and 1 and turn it into an index
      call random_number(r)
      r_int = int(r*(size(new_idxs)-i+1)) + 1

      ! get the number associated with that index and return it (this is your actual random index you'll use in the your program)
      rand_idxs(i) = new_idxs(r_int)

      ! return a new matrix without the index you chose
      new_idxs(1:r_int-1) = new_idxs(1:r_int-1)
      new_idxs(r_int:SIZE(new_idxs)-1) = new_idxs((r_int+1):SIZE(new_idxs))

    END DO
    
  END SUBROUTINE pick_random

END MODULE routinesMultTox
