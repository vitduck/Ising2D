MODULE ISING
    IMPLICIT NONE 
    REAL, PARAMETER               :: TEMP = 3.5, dTEMP = 0.05
    INTEGER, PARAMETER            :: LSIZE = 20
    INTEGER, PARAMETER            :: NTEMP = 40, NSTEP = 10000, ESTEP = 5000
    INTEGER, PARAMETER            :: NSPIN = LSIZE*LSIZE 
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
    INTEGER                      :: i, j 

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

SUBROUTINE METROPOLIS(tprob, spin, i, j, idata, ene, mag)
    IMPLICIT NONE 
    REAL, INTENT(IN)                :: tprob(1:8) 
    INTEGER, INTENT(INOUT)          :: spin(LSIZE+2,LSIZE+2) 
    INTEGER, INTENT(IN)             :: i, j
    INTEGER, INTENT(IN), OPTIONAL   :: idata
    INTEGER, INTENT(OUT), OPTIONAL  :: ene, mag
    INTEGER, SAVE                   :: local_ene, local_mag
    INTEGER                         :: dE 
    REAL                            :: random 
    LOGICAL                         :: move 

    move = .FALSE. 
    dE = 2*spin(i,j)*(spin(i-1,j) + spin(i+1,j) + spin(i,j-1) + spin(i,j+1))

    ! Metropolis condition
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
        ! Periodic boundary condition
        IF ( i == 2 )       spin(LSIZE+2,j) = spin(i,j)
        IF ( i == LSIZE+1 ) spin(1,j)       = spin(i,j) 
        IF ( j == 2 )       spin(i,LSIZE+2) = spin(i,j)
        IF ( j == LSIZE+1 ) spin(i,1)       = spin(i,j)
    END IF 

    ! data collection if idata, ene, and mag are passed to subroutine
    IF ( PRESENT(idata) .AND. PRESENT(ene) .AND. PRESENT(mag) ) THEN 
        IF ( idata == 1 ) THEN 
            local_ene = TOTAL_ENE(spin)
            local_mag = TOTAL_MAG(spin)
        END IF 
        IF ( move ) THEN 
            ene = local_ene + dE 
            mag = local_mag + 2*spin(i,j)
        ELSE 
            ene = local_ene 
            mag = local_mag 
        END IF 
        ! save value of local_ene, and local_mag for next pass
        local_ene = ene 
        local_mag = mag 
    END IF  
END SUBROUTINE METROPOLIS

SUBROUTINE RESPONSE_FUNCTION(ene_per_spin, mag_per_spin, temperature, output) 
    IMPLICIT NONE 
    REAL, INTENT(IN)            :: ene_per_spin(NSTEP*NSPIN), mag_per_spin(NSTEP*NSPIN)
    REAL, INTENT(IN)            :: temperature
    INTEGER, INTENT(IN)         :: output
    REAL                        :: ave_ene, ave_ene2, ave_mag, ave_mag2, Cv, chi
    
    ! Heat capacity
    ave_ene      = SUM(ene_per_spin)/(NSTEP*NSPIN)
    ave_ene2     = SUM(ene_per_spin**2)/(NSTEP*NSPIN)
    Cv           = NSPIN*(ave_ene2 - ave_ene**2)/(temperature**2)

    ! Magnetic susceptability 
    ave_mag      = SUM(mag_per_spin)/(NSTEP*NSPIN)
    ave_mag2     = SUM(mag_per_spin**2)/(NSTEP*NSPIN)
    chi          = NSPIN*(ave_mag2 - ave_mag**2)/(temperature)

    WRITE (output, '(F7.3, 4(3X, F15.8))') temperature, ave_ene, ave_mag, Cv, chi
END SUBROUTINE RESPONSE_FUNCTION

SUBROUTINE SPIN_CORRELATION(s, ss1, ss2, ss3, ss4, ss5, temperature, output) 
    IMPLICIT NONE 
    REAL, INTENT(IN)            :: ss1(NSTEP*NSPIN), ss2(NSTEP*NSPIN), ss3(NSTEP*NSPIN),&
                                   ss4(NSTEP*NSPIN), ss5(NSTEP*NSPIN), s(NSTEP*NSPIN) 
    REAL, INTENT(IN)            :: temperature 
    INTEGER, INTENT(IN)         :: output 
    
    REAL                        :: G2_1, G2_2, G2_3, G2_4, G2_5
    REAL                        :: s2 

    s2 = (SUM(s)/(NSTEP*NSPIN))**2 
    
    G2_1 = SUM(ss1)/(NSTEP*NSPIN) - s2
    G2_2 = SUM(ss2)/(NSTEP*NSPIN) - s2
    G2_3 = SUM(ss3)/(NSTEP*NSPIN) - s2
    G2_4 = SUM(ss4)/(NSTEP*NSPIN) - s2
    G2_5 = SUM(ss5)/(NSTEP*NSPIN) - s2

    WRITE (output, '(F7.3, 5(3X, F15.8))') temperature, G2_1, G2_2, G2_3, G2_4, G2_5
END SUBROUTINE SPIN_CORRELATION 

SUBROUTINE SPIN_CONFIGURATION(spin, output) 
    IMPLICIT NONE 
    INTEGER, INTENT(INOUT)      :: spin(LSIZE+2,LSIZE+2)
    INTEGER, INTENT(IN)         :: output
    INTEGER                     :: i, j  

    DO i = 2, LSIZE + 1 
        DO j = 2, LSIZE + 1 
            IF ( spin(i,j) == 1 )  WRITE (output, '(I3,2(3X,I3))') i, j, 0 
            IF ( spin(i,j) == -1 ) WRITE (output, '(I3,2(3X,I3))') i, 0, j 
        END DO 
    END DO 
    ! two blank space to separate between data block
    WRITE (output,*)
    WRITE (output,*)
END SUBROUTINE SPIN_CONFIGURATION

INTEGER FUNCTION TOTAL_ENE(spin) 
    IMPLICIT NONE
    INTEGER, INTENT(INOUT)         :: spin(LSIZE+2,LSIZE+2) 
    INTEGER                        :: i, j  

    TOTAL_ENE = 0.0
    DO i = 2, LSIZE + 1 
        DO j = 2, LSIZE + 1 
            TOTAL_ENE = TOTAL_ENE - spin(i,j)*(spin(i-1,j) + spin(i+1,j) + spin(i,j-1) + spin(i,j+1))
            END DO 
    END DO 
    
    ! Double counting
    TOTAL_ENE = TOTAL_ENE/2
END FUNCTION TOTAL_ENE

INTEGER FUNCTION TOTAL_MAG(spin)
    IMPLICIT NONE 
    INTEGER, INTENT(INOUT)         :: spin(LSIZE+2,LSIZE+2) 

    TOTAL_MAG = SUM(spin(2:LSIZE+1,2:LSIZE+1))
END FUNCTION TOTAL_MAG
END MODULE ISING
