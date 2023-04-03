!TAKES UNEDITED INPUT OF ALUMINUM AND GRAPHENE AND MAKES TWO VERSIONS OF
!A HYBRID NANOPARTICLE: ONE WITH A STEP IN THE GRAPHENE, AND ONE
!WITHOUT IT.


IMPLICIT NONE
INTEGER :: N_AL, N_C, I, J, TOTAL_ATOM
REAL*8 :: SEPARATION,OFFSET,DUMMY1,DUMMY2,DUMMY3

! COORDINATE ARRAYS PLAY MULTIPLE ROLES TO SAVE ROOM
REAL*8,ALLOCATABLE :: R1(:,:),R2(:,:)  ! STORES ORIGINAL COORDINATES
REAL*8,ALLOCATABLE :: RN1(:,:),RN2(:,:) ! STORES NEW COORDINATES


CHARACTER(2),ALLOCATABLE :: SPECIES1(:),SPECIES2(:)  ! ELEMENTS
REAL*8 :: EXT(2,3,2)  ! (SPECIES,DIRECTION,MIN=1 OR MAX=2)
REAL*8 :: DIM(2,3),CENTER(2,3),CENTER_NEW(2,3)















!----------------------------------------------------------------------
! WHERE FILES ARE READ
!----------------------------------------------------------------------

! STRUCTURE OF INPUT FILES SHOULD BE:

!  ------------------------>

! #NUMBER OF ATOMS
!
! X_COORD Y_COORD Z_COORD




OPEN(1,FILE = "Alrod.xyz",ACTION = 'READ')
OPEN(2,FILE = "graphene_sheet.xyz", ACTION = 'READ')

! READ IN NUMBER OF ATOMS
READ(1,*) N_AL
READ(1,*)
READ(2,*) N_C
READ(2,*)
TOTAL_ATOM = N_AL + N_C



! SET UP STORAGE ARRAYS
ALLOCATE(R1(3,N_AL),R2(3,N_C),SPECIES1(N_AL),SPECIES2(N_C))
ALLOCATE(RN1(3,N_AL),RN2(3,N_C))


R1 = 0.D0
R2 = 0.D0
RN1 = 0.D0
RN2 = 0.D0

! READ IN COORDINATES
DO I = 1,N_AL
 READ(1,*) SPECIES1(I),(R1(J,I),J=1,3)
ENDDO
DO I = 1,N_C
 READ(2,*) SPECIES2(I),(R2(J,I),J=1,3)
ENDDO


CLOSE(1);CLOSE(2)

!----------------------------------------------------------------------
! STRUCTURES ARE CENTERED
!----------------------------------------------------------------------



! EXTREMES ARE STORED HERE

	! FOR AL
EXT(1,1,1) = MINVAL(R1(1,:))
EXT(1,1,2) = MAXVAL(R1(1,:))
EXT(1,2,1) = MINVAL(R1(2,:))
EXT(1,2,2) = MAXVAL(R1(2,:))
EXT(1,3,1) = MINVAL(R1(3,:))
EXT(1,3,2) = MAXVAL(R1(3,:))


	! FOR C
EXT(2,1,1) = MINVAL(R2(1,:))
EXT(2,1,2) = MAXVAL(R2(1,:))
EXT(2,2,1) = MINVAL(R2(2,:))
EXT(2,2,2) = MAXVAL(R2(2,:))
EXT(2,3,1) = MINVAL(R2(3,:))
EXT(2,3,2) = MAXVAL(R2(3,:))


! LENGTH OF ALL SIDES OF EACH STRUCTURE
DIM(1,1) = ABS(EXT(1,1,1) + EXT(1,1,2)) ! AL - X
DIM(1,2) = ABS(EXT(1,2,1) + EXT(1,2,2)) ! AL - Y
DIM(1,3) = ABS(EXT(1,3,1) + EXT(1,3,2)) ! AL - Z
DIM(2,1) = ABS(EXT(2,1,1) + EXT(2,1,2)) ! C - X
DIM(2,2) = ABS(EXT(2,2,1) + EXT(2,2,2)) ! C - Y
DIM(2,3) = ABS(EXT(2,3,1) + EXT(2,3,2)) ! C - Z

! LOCATION OF CENTER IN BOTH STRUCTURES
DO I = 1,3
 CENTER(1,I) = 0.5D0*DIM(1,I) ! AL CENTER
 CENTER(2,I) = 0.5D0*DIM(2,I) ! C CENTER
ENDDO
  



! CENTER EXECUTION 

!====================================================================
!====================================================================
!====================================================================
! CENTERED STRUCTURES ARE WRITTEN OUT HERE WITH UNPHYSICAL CLIPPING
! ----------------------> FORT.3




WRITE(3,*) TOTAL_ATOM
WRITE(3,*)
DO I = 1,N_AL
 RN1(:,I) = R1(:,I) - CENTER(1,:)
 WRITE(3,*) SPECIES1(I),RN1(:,I)
ENDDO
DO I = 1,N_C
 RN2(:,I) = R2(:,I) !- CENTER(2,:)
 WRITE(3,*) SPECIES2(I),RN2(:,I)
ENDDO
WRITE(*,*) "CLIPPED STRUCTURE GENERATED"


!----------------------------------------------------------------------
! SHIFT AL BAR IN Z DIRECTION
!----------------------------------------------------------------------

! EXCHANGE VECTOR COORDINATES
R1= RN1; R2 = RN2
RN1 = 0.D0; RN2 = 0.D0


! IMPORTANT MESSAGE: CODY SUGGESTED I ALSO VARY THE DISTANCE BETWEEN
! THE BAR AND THE SHEET. I SHOULD TRY TO MAKE THIS DISTANCE EASILY
! ADJUSTABLE.

SEPARATION = 2.d0*6.4250683612  ! IN ATOMIC UNITS
! PRESENT VALUE IS BASED ON SPACING IN BILAYER GRAPHENE
! THIS VALUE IS CLOSE TO THE SPACINGS IN GRAPHITE
!====================================================================
!====================================================================
! CENTERED STRUCTURES ARE WRITTEN OUT HERE WITHOUT UNPHYSICAL CLIPPINg
! ----------------------> FORT.4




WRITE(4,*) TOTAL_ATOM
WRITE(4,*)
DO I = 1,N_AL
 WRITE(4,*) SPECIES1(I),R1(1,I),R1(2,I),R1(3,I)+0.5*SEPARATION
ENDDO
DO I = 1,N_C
 WRITE(4,*) SPECIES2(I),R2(1,I),R2(2,I),R2(3,I)-0.5*SEPARATION
ENDDO
WRITE(*,*) "UNCLIPPED STRUCTURE GENERATED"





!----------------------------------------------------------------------
! CLONE GRAPHENE LAYERS
!----------------------------------------------------------------------


! CHANGING OF THE GUARD
R1(3,:) = R1(3,:) + 0.5*SEPARATION
R2(3,:) = R2(3,:) - 0.5*SEPARATION



!====================================================================
!====================================================================
! LATEST INCARNATION OF STRUCTURE GENERATED WITH GRAPHENE CLONES
! ----------------------> FORT.7




WRITE(7,*) TOTAL_ATOM + N_C + N_C
WRITE(7,*)
DO I = 1,N_AL
 WRITE(7,*) SPECIES1(I),R1(1,I),R1(2,I),R1(3,I)
ENDDO
DO I = 1,N_C
 WRITE(7,*) SPECIES2(I),R2(1,I),R2(2,I),R2(3,I)
 WRITE(7,*) SPECIES2(I),R2(1,I),R2(2,I),R2(3,I)-.5*SEPARATION
 WRITE(7,*) SPECIES2(I),R2(1,I),R2(2,I),R2(3,I)-SEPARATION
ENDDO
WRITE(*,*) "CLONES GENERATED"

CLOSE(7) ! FOR LATER


!----------------------------------------------------------------------
! CUT FIRST LAYER
!----------------------------------------------------------------------



!!!	REQUIRES INPUT OF FILE LENGTH. CANNOT BE AUTOMATICALLY EXAMINED
! WITH JMOL


!====================================================================
!====================================================================
! sYSTEM WITH A GRAPHENE STEP
! ----------------------> FORT.8


! OFFSET APPLIED TO ENSURE PROPER LINE-UP BETWEEN AL AND STEP
OFFSET = 1


! FOR GETTING PROPER NUMBER OF ATOMS
J = 0
DO I = 1,N_C
 IF(R2(1,I).LT.MAXVAL(R1(1,:))+OFFSET) J = J + 1
ENDDO






WRITE(8,*) N_AL + N_C + N_C + J
WRITE(8,*)
DO I = 1,N_AL
 WRITE(8,*) SPECIES1(I),R1(1,I),R1(2,I),R1(3,I)
ENDDO

DO I = 1,N_C 
 IF(R2(1,I).LT.MAXVAL(R1(1,:))+OFFSET) WRITE(8,*) SPECIES2(I),R2(1,I),R2(2,I),R2(3,I)
 WRITE(8,*) SPECIES2(I),R2(1,I),R2(2,I),R2(3,I)-.5*SEPARATION
 WRITE(8,*) SPECIES2(I),R2(1,I),R2(2,I),R2(3,I)-SEPARATION
ENDDO
CLOSE(8) ! THIS WILL BE READ




!===================================================================
!===================================================================

!===================================================================
!===================================================================
!===================================================================
!===================================================================


! ROTATE AND PRINT BOTH HOMOGENOUS AND INHOMOGENOUS STRUCTURES
! FOUR TOTAL FILES TO BE GENERATED:
! INHOM.XYZ; INHOM.INP; HOM.XYZ; HOM.INP

OPEN(1000,FILE="INHOM.XYZ",ACTION='WRITE')
OPEN(1001,FILE="INHOM.INP",ACTION='WRITE')
OPEN(2000,FILE="HOM.XYZ",ACTION='WRITE')
OPEN(2001,FILE="HOM.INP",ACTION='WRITE')
OPEN(3000,FILE="fort.8",ACTION='READ')
OPEN(3001,FILE="fort.7",ACTION='READ')


DUMMY1 = MINVAL(R1(3,:))
DUMMY2 = MAXVAL(R2(3,:)) !GETS THE Z COORDINATES OF THE CENTERMOST
! ATOMS. TO BE USED FOR FINAL CENTERING FOR PROPER FLUX
DUMMY3 = ABS(0.5*(DUMMY1+DUMMY2))






DEALLOCATE(R1,SPECIES1)
TOTAL_ATOM = 0

READ(3000,*) TOTAL_ATOM
READ(3000,*)
ALLOCATE(SPECIES1(TOTAL_ATOM),R1(3,TOTAL_ATOM))


WRITE(1000,*) TOTAL_ATOM; WRITE(1001,*) TOTAL_ATOM
WRITE(1000,*)   ;WRITE(1001,*)
DO I = 1,TOTAL_ATOM   ! FOR INHOM.
 READ(3000,*) SPECIES1(I),(R1(J,I),J=1,3)


 WRITE(1000,*) SPECIES1(I),-1.D0*(R1(3,I)+DUMMY3),R1(2,I),R1(1,I)
 IF(SPECIES1(I)=="Al ") WRITE(1001,100) -1.D0*(R1(3,I)+DUMMY3),R1(2,I),R1(1,I),13,0
 IF(SPECIES1(I)=="C ")  WRITE(1001,100) -1.D0*(R1(3,I)+DUMMY3),R1(2,I),R1(1,I),6,0

ENDDO


100 format(3x,3f10.5,2i5)

DEALLOCATE(R1,SPECIES1)
TOTAL_ATOM = 0.D0

READ(3001,*) TOTAL_ATOM
READ(3001,*)
ALLOCATE(SPECIES1(TOTAL_ATOM),R1(3,TOTAL_ATOM))




WRITE(2000,*) TOTAL_ATOM; WRITE(2001,*) TOTAL_ATOM
WRITE(2000,*)   ;WRITE(2001,*)   
DO I = 1,TOTAL_ATOM   ! FOR HOM.
 READ(3001,*) SPECIES1(I),(R1(J,I),J=1,3)
 

 WRITE(2000,*) SPECIES1(I),-1.D0*(R1(3,I)+DUMMY3),R1(2,I),R1(1,I)
 IF(SPECIES1(I)=="Al ") WRITE(2001,*) -1.D0*(R1(3,I)+DUMMY3),R1(2,I),R1(1,I),"  13 0"
 IF(SPECIES1(I)=="C ")  WRITE(2001,*) -1.D0*(R1(3,I)+DUMMY3),R1(2,I),R1(1,I),"  6 0"

ENDDO


WRITE(*,*) "INP AND XYZ FILES GENERATE"

END
