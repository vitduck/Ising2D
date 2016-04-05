PROGRAM MONTECARLO
	USE ISING
	IMPLICIT NONE 
	REAL                          :: tprob(1:8)
	REAL                          :: temperature(NTEMP) 
	INTEGER                       :: spin(LSIZE+2,LSIZE+2) 
	INTEGER                       :: itemp

	! scan NTEMP temperature
	temperature = (/ (TEMP - (itemp-1)*dTEMP, itemp = 1, NTEMP) /)

	! initialize spin configuration with total magnetization = 0	
	CALL INIT_CONFIGURATION(spin, 0)
	
	! output files
	OPEN (UNIT = 10, FILE = 'output.dat')
	OPEN (UNIT = 20, FILE = 'configuration.dat')
	OPEN (UNIT = 30, FILE = 'correlation.dat')

	! main loop 
	DO itemp = 1, NTEMP 
		! compute the transition probability
		CALL INIT_TRANSITION_PROBABILITY(tprob, temperature(itemp)) 
		! equilibrate system for ESTEP 
		CALL EQUILIBRATION(tprob, spin)
		! statistical analysis for NSTEP 
		CALL STATISTICAL_ANALYSIS(tprob, spin, temperature(itemp), 10, 20, 30)
	END DO 		
	
	! finalize 
	CLOSE (10)
	CLOSE (20)
	CLOSE (30)
END PROGRAM MONTECARLO

SUBROUTINE EQUILIBRATION(tprob, spin)
	USE ISING, ONLY: LSIZE, ESTEP, METROPOLIS
	IMPLICIT NONE
	REAL, INTENT(IN)            :: tprob(1:8)
	INTEGER, INTENT(INOUT)      :: spin(LSIZE+2,LSIZE+2) 
	INTEGER                     :: i, j, istep 

	DO istep = 1, ESTEP 
		DO i = 2, LSIZE + 1 
			DO j = 2, LSIZE + 1
				CALL METROPOLIS(tprob, spin,  i, j)
			END DO 
		END DO
	END DO
END SUBROUTINE EQUILIBRATION

SUBROUTINE STATISTICAL_ANALYSIS(tprob, spin, temperature, output, configuration, correlation)
	USE ISING, ONLY: LSIZE, NSTEP, NSPIN,& 
		             METROPOLIS, RESPONSE_FUNCTION,&
		             SPIN_CONFIGURATION, SPIN_CORRELATION
	IMPLICIT NONE 
	REAL, INTENT(IN)            :: tprob(1:8)
	INTEGER, INTENT(INOUT)      :: spin(LSIZE+2,LSIZE+2) 
	REAL, INTENT(IN)            :: temperature
	INTEGER, INTENT(IN)         :: output, configuration, correlation 
	INTEGER                     :: ene(NSTEP*NSPIN), mag(NSTEP*NSPIN) 
	REAL                        :: ene_per_spin(NSTEP*NSPIN), mag_per_spin(NSTEP*NSPIN)
	REAL                        :: ss1(NSTEP*NSPIN), ss2(NSTEP*NSPIN), ss3(NSTEP*NSPIN),&
		                           ss4(NSTEP*NSPIN), ss5(NSTEP*NSPIN), s(NSTEP*NSPIN) 
	INTEGER                     :: idata, istep, i, j, l, m 

	! use the center spin for correlation
	l = LSIZE/2 + 1 
	m = LSIZE/2 + 1
	idata = 0 
	DO istep = 1, NSTEP 
		DO i = 2, LSIZE + 1 
			DO j = 2, LSIZE + 1
				idata = idata + 1
				CALL METROPOLIS(tprob, spin, i, j, idata, ene(idata), mag(idata))
				s(idata)   = spin(l,m) 
				! 1st nearest neighbor
				ss1(idata) = 0.25*spin(l,m)*(spin(l-1,m) + spin(l+1,m) + spin(l,m-1) + spin(l,m+1))
				! 2nd nearest neighbor
				ss2(idata) = 0.25*spin(l,m)*(spin(l-1,m-1) + spin(l-1,m+1) + spin(l+1,m-1) + spin(l+1,m+1))
				! 3rd nearest neighbor
				ss3(idata) = 0.25*spin(l,m)*(spin(l-2,m) + spin(l+2,m) + spin(l,m-2) + spin(l,m+2))
				! 4th nearest neighbor
				ss4(idata) = 0.125*spin(l,m)*(spin(l+2,m+1) + spin(l+1,m+2) + spin(l-1,m+2) + spin(l-2,m+1) + & 
					spin(l+2,m-1) + spin(l+1,m-2) + spin(l-1,m-2) + spin(l-2,m-1))
				! 5th nearest neighbor
				ss5(idata) = 0.25*spin(l,m)*(spin(l-2,m-2) + spin(l-2,m+2) + spin(l+2,m-2) + spin(l+2,m+2))
			END DO 
		END DO 
	END DO
	ene_per_spin = REAL(ene)/NSPIN
	mag_per_spin = ABS(REAL(mag))/NSPIN
	CALL RESPONSE_FUNCTION(ene_per_spin, mag_per_spin, temperature, output)
	CALL SPIN_CONFIGURATION(spin, configuration)
	CALL SPIN_CORRELATION(s, ss1, ss2, ss3, ss4, ss5, temperature, correlation) 
END SUBROUTINE STATISTICAL_ANALYSIS	
