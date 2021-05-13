MODULE functions_Emi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! General homemade fortran functions - not specific to
! this simulation, but used in routinesMultTox and toxinSolubilityNVT
!
! Emiko Zumbro
! ezumbro@mit.edu
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE

CONTAINS

	REAL(DP) FUNCTION round(n)

		USE parameters

		IMPLICIT NONE

	  	REAL(DP), INTENT(IN) :: n
	  	REAL(DP) :: remainder
	  	REAL(DP) :: div = 1.0

	  	round = NINT(n)


	END FUNCTION round

	FUNCTION FIND(mask)

		USE parameters

		IMPLICIT NONE

		INTEGER, DIMENSION(:), ALLOCATABLE :: find
		LOGICAL, DIMENSION(:), INTENT(IN) :: mask
		INTEGER, DIMENSION(:), ALLOCATABLE :: indxs
		INTEGER :: i, cnt

		cnt = COUNT(mask)
	    ALLOCATE (indxs(cnt))
		indxs=0
	    cnt=1
		DO i=1,SIZE(mask,1)
		  IF (mask(i)) THEN
		    indxs(cnt)=i
		    cnt = cnt+1
		  END IF
		END DO

		find = indxs

	END FUNCTION FIND

END MODULE functions_Emi