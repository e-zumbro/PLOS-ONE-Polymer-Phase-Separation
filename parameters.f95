MODULE parameters
!
	IMPLICIT NONE
	!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! THIS MODULE HAS THE DECLARATION OF MOST OF THE VARIALBES USED
! IN THE SIMULATION. THE ONLY VARIABLES THAT ARE IN ANOTHER MODULE
! ARE THE RANDOM NUMBER GENERATOR VARIABLES, AND THOSE THAT ARE
! DECLARED AT RUNTIME IN input.par
!
! Emiko Zumbro
! ezumbro@mit.edu
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 
	!
	! This statement is simply to select the precision
	!	
	INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)

	REAL(DP), PARAMETER :: pi=3.141592653589793238462643383279502884197_DP

	REAL(DP), PARAMETER :: N_A = 6.0221409E23

	! Number of dimensions of space (ie: x,y,z)
	INTEGER, PARAMETER :: dim = 3

	! ! External field, ie: flow
	! REAL(DP), PARAMETER :: dUdx = 0.0

	! Equilibrium length/starting length of the springs
	REAL(DP), PARAMETER :: l_sp = 1.0 
	REAL(DP), PARAMETER :: kbT = 1.0 ! thermal energy
	REAL(DP), PARAMETER :: k_sp = 200.0 ! spring constant for potential force connecting polymer beads

	REAL(DP), PARAMETER :: k_bind = 10.0 !k_sp ! spring constant of bond between polymer and toxin

	REAL(DP), PARAMETER :: d_tox = l_sp
	
	REAL(DP), PARAMETER :: l_bind = l_sp/2.0_DP + d_tox/2.0_DP! equilibrium length of spring after binding
	REAL(DP), PARAMETER :: d_bead = l_sp ! diameter of a bead

	! How often should you check if a binding event happens
	INTEGER, PARAMETER :: checkBindingInterval = 100 !1000
      
	REAL(DP), PARAMETER :: D = 1.0 ! diffusion coefficient of Polymer
	REAL(DP), PARAMETER :: DTox = D ! diffusion coefficient of Toxin
	REAL(DP), PARAMETER :: Vf = 0.0 ! velocity of the fluid
	REAL(DP), PARAMETER :: deltaT = .0001 ! time step
	REAL(DP), PARAMETER :: zeta = 1.0 ! Not needed for now, so set to arbitrary number (for hydrodynamic forces/friction of fluid)

	INTEGER :: maxTox
	INTEGER, PARAMETER :: timeAvgInt = 100
	INTEGER :: nTox

	! Schmolochowsky limit for how far the toxin can see the inhibitor beads
	REAL(DP), PARAMETER :: reach = l_bind+SQRT(D*deltaT*100_DP)

	INTEGER :: Nmer

	INTEGER :: bindingSites

	REAL(DP) :: tot_t

	! The number of beads
	INTEGER :: beads

	! size of side of box for periodic boundary conditions (origin is at center of box)
	REAL(DP) :: sizeBox

	! Energy barriers
	REAL(DP) :: stdDevPolyAff

	REAL(DP) :: delE_0_center

	INTEGER :: randNumSeedMultiplier

	INTEGER, PARAMETER :: seedSize = 33

	! The effective step size

	REAL(DP) :: muo

	! The strength of the flow

	REAL(DP) :: wi

	! stretching constant

	REAL(DP), PARAMETER :: gamma=100.0


	! bk is the kuhn length in units of a
	! dmax is the maximum bond length

	REAL(DP), PARAMETER :: bk=2.0,dmax=2.2

	! sigma is the space constant in the lennard-jones force

	REAL(DP), PARAMETER :: sigma=l_sp

	REAL(DP), PARAMETER :: sigmaToxTox=d_tox

	REAL(DP), PARAMETER :: sigmaToxPoly=l_bind

	REAL(DP), PARAMETER :: LJcutoff=3.0_DP

	! epsilon is the strength of the lennard-jones potential

	REAL(DP) :: epsilon

	REAL(DP) :: epsToxTox

	REAL(DP) :: epsToxPoly ! this is now adjustable

	! constants to control the stiffness of the chain - wormlike chain model
	REAL(DP) :: gammaw 

	!
	! Options for changing 
	!

	! Decides whether the inhibitor is a polymer or just monomers floating around, true = is a polymer and experiences connectivity forces, 0 = free monomers in solution, don't experience connectivity forces
	LOGICAL :: isPolymer

	! This controls whether when a toxin unbinds, if it is replaced with a new toxin starting in a random place or if it just continues on
	LOGICAL :: randNewToxin


	! Cycles and writing blocks

	INTEGER :: block,d_wrt= 10 !10000!5000


END MODULE parameters
