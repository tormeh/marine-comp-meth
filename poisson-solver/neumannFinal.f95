!gfortran neumannFinal.f95 -fimplicit-none -O3 -o neumannFinal -fdefault-real-8 && time ./neumannFinal
      
      FUNCTION F(X,Y,H)
        !the right side of the poisson equation
        INTEGER X,Y
        REAL H
        REAL F
        REAL :: XC
        REAL :: YC
        XC = (X-1)*H
        YC = (Y-1)*H
        !F=12-12*XC-12*YC
        F=(6-12*XC)*(3*YC**2-2*YC**3) + (3*XC**2-2*XC**3)*(6-12*YC)
        RETURN
      END FUNCTION
           
      FUNCTION ANALYTICAL(X,Y,H)
        !computes the analytical solution
        REAL ANALYTICAL
        INTEGER :: X
        INTEGER :: Y
        REAL :: H
        REAL :: XC
        REAL :: YC
        XC = (X-1)*H
        YC = (Y-1)*H
        !ANALYTICAL = 3*XC**2 + 3*YC**2 - 2*XC**3 - 2*YC**3
        ANALYTICAL = (3*YC**2 - 2*YC**3)*(3*XC**2 - 2*XC**3)
        RETURN
      END FUNCTION
      
      FUNCTION NEWESTIMATE(X,Y,H,PHIS,LENGTH)
        !computes the value of ((x-1)*h,(y-h)*h). Valid everywhere.
        INTEGER X
        INTEGER  Y
        REAL H
        INTEGER LENGTH
        REAL NEWVALUE
        REAL :: PHIS(LENGTH,LENGTH)
        REAL NEWESTIMATE
        REAL SIMPLEESTIMATE
        REAL F
        REAL :: ORIGINDEMAND
        IF (X==1 .AND. Y==1) THEN
          PHIS(1,2)=0.0
          PHIS(2,1)=0.0
          PHIS(2,2)=0.0
          NEWESTIMATE=0.0
        ELSE IF (X==1 .AND. Y==LENGTH) THEN
          NEWESTIMATE = (PHIS(2,LENGTH-1) + PHIS(2,LENGTH) + PHIS(1,LENGTH-1))/3.0
        ELSE IF (X==LENGTH .AND. Y==1) THEN
          NEWESTIMATE = (PHIS(LENGTH-1,2) + PHIS(LENGTH,2) + PHIS(LENGTH-1,1))/3.0
        ELSE IF (X==LENGTH .AND. Y==LENGTH) THEN
          NEWESTIMATE = (PHIS(LENGTH-1,LENGTH-1) + PHIS(LENGTH-1,LENGTH) + PHIS(LENGTH,LENGTH-1))/3.0
        ELSE IF (X==1) THEN
          NEWVALUE = (1.0/4.0)*(2*PHIS(X+1,Y)+PHIS(X,Y+1)+PHIS(X,Y-1)-H*H*F(X,Y,H))
          NEWESTIMATE = NEWVALUE
        ELSE IF (X==LENGTH) THEN
          NEWVALUE = (1.0/4.0)*(PHIS(X,Y+1)+2*PHIS(X-1,Y)+PHIS(X,Y-1)-H*H*F(X,Y,H))
          NEWESTIMATE = NEWVALUE
        ELSE IF (Y==1) THEN
          NEWVALUE = (1.0/4.0)*(PHIS(X+1,Y)+2*PHIS(X,Y+1)+PHIS(X-1,Y)-H*H*F(X,Y,H))
          NEWESTIMATE = NEWVALUE
        ELSE IF (Y==LENGTH) THEN
          NEWVALUE = (1.0/4.0)*(PHIS(X+1,Y)+PHIS(X-1,Y)+2*PHIS(X,Y-1)-H*H*F(X,Y,H))
          NEWESTIMATE = NEWVALUE
        ELSE
          NEWESTIMATE = SIMPLEESTIMATE(X,Y,H,PHIS,LENGTH)
        END IF
        RETURN
      END FUNCTION
           
      
      FUNCTION SIMPLEESTIMATE(X,Y,H,PHIS,LENGTH)
        !computes a new estimate for the value of (x*h,y*h). Not valid on the boundary
        INTEGER X,Y
        REAL H
        REAL SIMPLEESTIMATE
        INTEGER LENGTH
        REAL F
        REAL :: PHIS(LENGTH,LENGTH)
        SIMPLEESTIMATE = (1.0/4.0)*(PHIS(X+1,Y)+PHIS(X,Y+1)+PHIS(X-1,Y)+PHIS(X,Y-1)-H*H*F(X,Y,H))
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
      
      FUNCTION ANALERROR(PHIS,LENGTH,H,SIZE)
        !computes the analytically determined error
        INTEGER SIZE
        INTEGER LENGTH
        REAL H
        REAL :: PHIS(LENGTH,LENGTH)
        REAL ANALYTICAL
        REAL AVGERROR
        INTEGER I,J
        REAL ANALERROR
        AVGERROR = 0.0
        DO J=1,LENGTH
          DO I=1,LENGTH
          AVGERROR = AVGERROR + ((PHIS(I,J)-ANALYTICAL(I,J,H))**2)/SIZE
          END DO
        END DO
        ANALERROR = AVGERROR
        RETURN
      END FUNCTION
      
      PROGRAM SOLVER
        REAL, PARAMETER :: H = 0.01
        INTEGER, PARAMETER :: LENGTH = (1.0/H)+1
        INTEGER, PARAMETER :: SIZE = LENGTH*LENGTH
        REAL :: PHIS(LENGTH, LENGTH)
        INTEGER, PARAMETER :: out_unit=20
        REAL :: HIGHESTCHANGE
        REAL :: NEWVALUE
        REAL SIMPLEESTIMATE
        REAL NEWESTIMATE
        REAL ANALYTICAL
        INTEGER I,J
        REAL HIGHESTCHANGEFUN
        REAL LOWHIGHESTCHANGE
        REAL AVGCHANGE
        REAL LOWAVGCHANGE
        REAL LOWCHANGE
        REAL AVGERROR
        REAL NUMAVGERROR
        INTEGER NUMITERATIONS
        REAL ANALERROR
        LOWAVGCHANGE = 20.0
        AVGCHANGE = 1.0
        LOWHIGHESTCHANGE = 20.0
        HIGHESTCHANGE = 10.0
        
        !a friendly reminder that what the real coordinate is (x-1)*h, not x*h
        WRITE (*,*) "(LENGTH-1)*H IS ", ((LENGTH-1)*H)
        
        DO I=1,LENGTH
          DO J=1,LENGTH
            PHIS(I,J) = RAND(0)*10
          END DO
        END DO
        
        NUMITERATIONS = 0
        !the actual computation is performed here
        DO WHILE  ((LOWHIGHESTCHANGE/HIGHESTCHANGE) > 1.00001 .OR. (LOWAVGCHANGE/AVGCHANGE) > 1.00001) !(NUMITERATIONS < 200000) !(ANALERROR(PHIS,LENGTH,H,SIZE)>0.000000006)! !((SIMPLEESTIMATE(2,2,H,PHIS,LENGTH)-0)**2>0.00001)!
          NUMITERATIONS = NUMITERATIONS + 1
          LOWAVGCHANGE = LOWCHANGE(LOWAVGCHANGE, AVGCHANGE)
          LOWHIGHESTCHANGE = LOWCHANGE(LOWHIGHESTCHANGE, HIGHESTCHANGE)
          HIGHESTCHANGE = 0.0
          AVGCHANGE = 0.0
          DO J=1,LENGTH
            DO I=1,LENGTH
                NEWVALUE = NEWESTIMATE(I,J,H,PHIS,LENGTH)
                HIGHESTCHANGE = HIGHESTCHANGEFUN(PHIS(I,J),NEWVALUE,HIGHESTCHANGE)
                AVGCHANGE = AVGCHANGE + ((PHIS(I,J)-NEWVALUE)**2)/SIZE
                PHIS(I,J) = NEWVALUE
            END DO
          END DO
        END DO
        
        
        
        WRITE (*,*) "NUMITERATIONS IS ", NUMITERATIONS
        !finds and prints the numerically determined error
        NUMAVGERROR = 0.0
        DO I=1,LENGTH
          DO J=1,LENGTH
          NUMAVGERROR = NUMAVGERROR + ((PHIS(I,J)-NEWESTIMATE(I,J,H,PHIS,LENGTH))**2)/SIZE
          END DO
        END DO
        WRITE (*,*) "AVERAGE SQUARE NUMERICALLY ESTIMATED ERROR IS ", NUMAVGERROR
        
        !performance metrics. Analytically determined error is printed
        AVGERROR = 0.0
        DO I=1,LENGTH
          DO J=1,LENGTH
          AVGERROR = AVGERROR + ((PHIS(I,J)-ANALYTICAL(I,J,H))**2)/SIZE
          END DO
        END DO
        WRITE (*,*) "AVERAGE SQUARE ANALYTICALLY DETERMINED ERROR IS ", AVGERROR
        
        !writes results to file
        open (unit=out_unit,file="results.txt",action="write",status="replace")
        WRITE (OUT_UNIT,'(I4)') LENGTH
        WRITE (OUT_UNIT,'(I4)') LENGTH
        DO I=1,LENGTH
          DO J=1,LENGTH
              WRITE (OUT_UNIT,'(F0.5)') PHIS(I,J)
          END DO
        END DO
        close (out_unit)
      END
