!gfortran dirichlet.f95 -fimplicit-none -O3 -o dirichlet && ./dirichlet 
      
      
      
      FUNCTION BOUNDARY(X,Y,H)
        !computes the value of (x*h,y*h). Valid only on the boundary.
        REAL BOUNDARY
        INTEGER :: X
        INTEGER :: Y
        REAL :: H
        BOUNDARY = (1/4.0)*(X*X*H*H+Y*Y*H*H)
        RETURN
      END FUNCTION
      
      FUNCTION ANALYTICAL(X,Y,H)
        REAL ANALYTICAL
        INTEGER :: X
        INTEGER :: Y
        REAL :: H
        ANALYTICAL = (1/4.0)*(X*X*H*H+Y*Y*H*H)
        RETURN
      END FUNCTION
      
      FUNCTION F(X,Y,H)
        !the right side of the poisson equation
        REAL F
        INTEGER X,Y
        REAL H
        F=1.0
        RETURN
      END FUNCTION
      
      
       FUNCTION NEWESTIMATE(X,Y,H,PHIS,LENGTH)
        !computes a new estimate for the value of (x*h,y*h). Not valid on the boundary
        INTEGER X,Y
        REAL H
        REAL NEWESTIMATE
        INTEGER LENGTH
        REAL F
        REAL :: PHIS(LENGTH,LENGTH)
        !WRITE (*,'(a,g12.4)') "F-PART IS ", H*H*F(X,Y,H)
        NEWESTIMATE = (1.0/4.0)*(PHIS(X+1,Y)+PHIS(X,Y+1)+PHIS(X-1,Y)+PHIS(X,Y-1)-H*H*F(X,Y,H))
        !WRITE (*,'(a,g12.4)') "NEWESTIMATE IS ", NEWESTIMATE
        RETURN
      END FUNCTION
      
      FUNCTION HIGHESTCHANGEFUN(OLD,NEW,PREVHIGHEST)
        !returns the gighest change in (phi-)value given a new old value, new assignment and the previous highest change
        REAL :: OLD
        REAL :: NEW
        REAL :: PREVHIGHEST
        REAL :: CHANGE
        REAL HIGHESTCHANGEFUN
        CHANGE = ABS(OLD-NEW)
        IF (CHANGE .GE. PREVHIGHEST) THEN
          HIGHESTCHANGEFUN = CHANGE
        ELSE
          HIGHESTCHANGEFUN = PREVHIGHEST
        END IF
        RETURN
      END FUNCTION
      
      FUNCTION LOWCHANGE(OLDLOWEST, LASTLOW)
        !returns the lesser of two variables. Intended to compare change magnitudes to determine the lowest one.
        REAL LOWCHANGE
        REAL OLDLOWEST
        REAL LASTLOW
        IF (OLDLOWEST .LT. LASTLOW) THEN
          LOWCHANGE = OLDLOWEST
        ELSE
          LOWCHANGE = LASTLOW
        END IF
      END FUNCTION
      
      PROGRAM SOLVER
        REAL, PARAMETER :: H = 0.003
        INTEGER, PARAMETER :: LENGTH = (1.0/H)+1
        INTEGER, PARAMETER :: SIZE = LENGTH*LENGTH
        REAL :: PHIS(LENGTH, LENGTH)
        INTEGER, PARAMETER :: out_unit=20
        REAL :: HIGHESTCHANGE
        REAL :: NEWVALUE
        REAL NEWESTIMATE
        REAL BOUNDARY
        INTEGER I,J
        REAL HIGHESTCHANGEFUN
        REAL LOWHIGHESTCHANGE
        REAL AVGCHANGE
        REAL LOWAVGCHANGE
        REAL LOWCHANGE
        REAL AVGERROR
        INTEGER NUMITERATIONS
        LOWAVGCHANGE = 20.0
        AVGCHANGE = 1.0
        LOWHIGHESTCHANGE = 20.0
        HIGHESTCHANGE = 10.0
        
        !random initialization
        DO I=1,LENGTH
          DO J=1,LENGTH
            PHIS(I,J) = RAND(0)*10
          END DO
        END DO
        
        !the computation loop. Computes the average change of phi and the highest change of phi. 
        !Remembers the lowest such values for any iteration. Terminates if both the average and highest change is higher than the historical lows.
        NUMITERATIONS = 0
        DO WHILE ((LOWHIGHESTCHANGE/HIGHESTCHANGE) > 1.0 .OR. (LOWAVGCHANGE/AVGCHANGE) > 1.0)
          NUMITERATIONS = NUMITERATIONS + 1
          LOWAVGCHANGE = LOWCHANGE(LOWAVGCHANGE, AVGCHANGE)
          LOWHIGHESTCHANGE = LOWCHANGE(LOWHIGHESTCHANGE, HIGHESTCHANGE)
          HIGHESTCHANGE = 0.0
          AVGCHANGE = 0.0
          DO I=1,LENGTH
            DO J=1,LENGTH
                IF (I==1 .OR. I==LENGTH .OR. J==1 .OR. J==LENGTH) THEN
                  NEWVALUE = BOUNDARY(I,J,H)
                ELSE
                  NEWVALUE = NEWESTIMATE(I,J,H,PHIS,LENGTH)
                END IF
                HIGHESTCHANGE = HIGHESTCHANGEFUN(PHIS(I,J),NEWVALUE,HIGHESTCHANGE)
                AVGCHANGE = AVGCHANGE + ABS(PHIS(I,J)-NEWVALUE)/SIZE
                !WRITE (*,*) "OLD VALUE IS ", PHIS(I,J)
                PHIS(I,J) = NEWVALUE
                !WRITE (*,*) "NEW VALUE IS ", PHIS(I,J)
            END DO
          END DO
        END DO
        
        !performance metrics number of iterations needed to compute result and the end average absolute error is printed 
        WRITE (*,*) "NUMITERATIONS IS ", NUMITERATIONS
        AVGERROR = 0.0
        DO I=1,LENGTH
          DO J=1,LENGTH
          IF (I==1 .OR. I==LENGTH .OR. J==1 .OR. J==LENGTH) THEN
            AVGERROR = AVGERROR + ABS(PHIS(I,J)-BOUNDARY(I,J,H))/SIZE
          ELSE
            AVGERROR = AVGERROR + ABS(PHIS(I,J)-NEWESTIMATE(I,J,H,PHIS,LENGTH))/SIZE
          END IF
          END DO
        END DO
        WRITE (*,*) "AVGERROR IS ", AVGERROR
        
        !write to file
        open (unit=out_unit,file="results.txt",action="write",status="replace")
        WRITE (OUT_UNIT,'(I3)') LENGTH
        WRITE (OUT_UNIT,'(I3)') LENGTH
        DO I=1,LENGTH
          DO J=1,LENGTH
              WRITE (OUT_UNIT,'(F10.5)') PHIS(I,J)
          END DO
        END DO
        close (out_unit)
      END
