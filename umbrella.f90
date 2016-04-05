MODULE ISING
	IMPLICIT NONE 
	REAL, PARAMETER               :: TEMP = 2.6
	INTEGER, PARAMETER            :: LSIZE = 20
	INTEGER, PARAMETER            :: NSCAN = 20, NSTEP = 5000
	INTEGER, PARAMETER            :: NSPIN = LSIZE*LSIZE 
	INTEGER, PARAMETER            :: MAG1  = -NSPIN, MAG2 = NSPIN
CONTAINS

SUBROUTINE INIT_TRANSITION_PROBABILITY(tprob, temperature) 
	IMPLICIT NONE 
	REAL, INTENT(OUT)             :: tprob(1:8) 
	REAL, INTENT(IN)              :: temperature 
	INTEGER                       :: dE

	! two possible value: exp(-4J) and exp(-8J)
	DO dE = 4, 8, 4
		tprob(dE) = EXP(-dE/temperature) 
	END DO 
END SUBROUTINE INIT_TRANSITION_PROBABILITY 

SUBROUTINE INIT_CONFIGURATION(spin, init_mag)
	IMPLICIT NONE 
	INTEGER, INTENT(OUT)         :: spin(LSIZE+2,LSIZE+2)
	INTEGER, INTENT(IN)          :: init_mag 
	INTEGER                      :: spin_up, spin_count
	INTEGER                    	 :: i, j 

	spin_count = 0 
	spin_up = (NSPIN + init_mag )/2
	DO i = 2, LSIZE + 1 
		DO j = 2, LSIZE + 1 
			IF ( spin_count < spin_up ) THEN 
				spin_count  = spin_count + 1	
				spin(i,j) = 1 
			ELSE 
				spin(i,j) = -1
			END IF 
		END DO 
	END DO 

	! Periodic boundary conditon 
	spin(1,:)       = spin(LSIZE+1,:) 
	spin(LSIZE+2,:) = spin(2,:)
	spin(:,1)       = spin(:,LSIZE+1)
	spin(:,LSIZE+2) = spin(:,2)
END SUBROUTINE INIT_CONFIGURATION

SUBROUTINE UMBRELLA(i, j, spin, tprob, total_mag, lmag, rmag)
	IMPLICIT NONE 
	INTEGER, INTENT(IN)         :: i, j
	INTEGER, INTENT(INOUT)      :: spin(LSIZE+2,LSIZE+2) 
	INTEGER, INTENT(IN)         :: lmag, rmag
	REAL, INTENT(IN)            :: tprob(1:8) 
	INTEGER, INTENT(INOUT)      :: total_mag 
	
	REAL                        :: random 
	INTEGER                     :: dE 
	LOGICAL                     :: move 
	
	move = .FALSE.

	dE = 2*spin(i,j)*(spin(i-1,j) + spin(i+1,j) + spin(i,j-1) + spin(i,j+1))
	
	IF ( dE <= 0 ) THEN 
		move = .TRUE. 
	ELSE 
		CAll RANDOM_NUMBER(random) 
		IF ( tprob(dE) > random ) THEN 
			move = .TRUE.
		END IF 
	END IF 

	IF ( move ) THEN 
		! flip the spin 
		spin(i,j) = -spin(i,j) 
		total_mag = total_mag + 2*spin(i,j)
		! if step outsite the window boundary, revert the step
		IF ( total_mag < lmag .OR. total_mag > rmag ) THEN 
			spin(i,j) = -spin(i,j) 
			total_mag = total_mag + 2*spin(i,j)
		END IF 
		! Periodic boundary condition
		IF ( i == 2 )       spin(LSIZE+2,j) = spin(i,j)
		IF ( i == LSIZE+1 ) spin(1,j)       = spin(i,j) 
		IF ( j == 2 )       spin(i,LSIZE+2) = spin(i,j)
		IF ( j == LSIZE+1 ) spin(i,1)       = spin(i,j)
	END IF 
END SUBROUTINE UMBRELLA

SUBROUTINE CONSTRUCT_HISTOGRAM(tprob, spin, mag, iscan, histogram) 
	IMPLICIT NONE 
	REAL, INTENT(IN)           :: tprob(1:8) 
	INTEGER, INTENT(INOUT)     :: spin(LSIZE+2, LSIZE+2), mag(NSCAN+1)
	INTEGER, INTENT(IN)        :: iscan 
	REAL, INTENT(OUT)          :: histogram(mag(iscan):mag(iscan+1))

	INTEGER                    :: i, j, istep, imag 
	INTEGER                    :: total_mag 

	total_mag = SUM(spin(2:LSIZE+1,2:LSIZE+1))

	DO istep = 1, NSTEP 
		DO i = 2, LSIZE + 1 
			DO j = 2, LSIZE + 1 
				CALL UMBRELLA(i, j, spin, tprob, total_mag, mag(iscan), mag(iscan+1))
				DO imag = mag(iscan), mag(iscan+1), 2 
					IF (total_mag == imag) histogram(imag) = histogram(imag) + 1 
				END DO 
			END DO 
		END DO 
	END DO 

	! normalization 
	histogram(:) = histogram(:)/(NSTEP*NSPIN) 
END SUBROUTINE CONSTRUCT_HISTOGRAM 

SUBROUTINE FREE_ENE_PROFILE(mag, iscan, histogram, free_ene, memory, output) 
	IMPLICIT NONE 
	INTEGER, INTENT(IN)        :: mag(NSCAN+1), iscan
	REAL, INTENT(IN)           :: histogram(mag(iscan):mag(iscan+1)) 
	REAL, INTENT(OUT)          :: free_ene(mag(iscan):mag(iscan+1)) 
	REAL, INTENT(INOUT)        :: memory
	INTEGER, INTENT(IN)        :: output 
	INTEGER                    :: imag

	! A(M) = -logP(M)
	free_ene(:) = -LOG(histogram(:))

	! make continuous transition 
	IF ( iscan /= 1 ) THEN 
		free_ene(:) = free_ene(:) + (memory - free_ene(mag(iscan)))
	END IF 

	DO imag = mag(iscan), mag(iscan+1), 2 
		WRITE (output, '(F15.8, 5X, E15.8)') REAL(imag)/NSPIN, free_ene(imag)
	END DO
	! store the last value to join with first point of next window
	memory = free_ene(mag(iscan+1))
END SUBROUTINE FREE_ENE_PROFILE
END MODULE ISING
