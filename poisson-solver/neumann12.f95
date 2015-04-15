!gfortran neumann12.f95 -fimplicit-none -O3 -o neumann12 && time ./neumann12
      
      FUNCTION F1(X,Y,H)
        INTEGER X,Y
        REAL H
        REAL F1
        REAL :: XC
        REAL :: YC
        XC = (X-1)*H
        YC = (Y-1)*H
        F1=12-12*XC-12*YC
        !F1=(6-12*XC)*(3*YC**2-2*YC**3) + (3*XC**2-2*XC**3)*(6-12*YC)
        RETURN
      END FUNCTION
           
      FUNCTION ANALYTICAL(X,Y,H)
        REAL ANALYTICAL
        INTEGER :: X
        INTEGER :: Y
        REAL :: H
        REAL :: XC
        REAL :: YC
        XC = (X-1)*H
        YC = (Y-1)*H
        ANALYTICAL = 3*XC**2 + 3*YC**2 - 2*XC**3 - 2*YC**3
        !ANALYTICAL = (3*YC**2 - 2*YC**3)*(3*XC**2 - 2*XC**3)
        RETURN
      END FUNCTION
      
      FUNCTION BOUNDARY(X,Y,H,PHIS,LENGTH)
        INTEGER X
        INTEGER  Y
        REAL H
        INTEGER LENGTH
        REAL NEWVALUE
        REAL :: PHIS(LENGTH,LENGTH)
        REAL BOUNDARY
        REAL NEWESTIMATE
        REAL F1
        IF (X==1 .AND. Y==1) THEN
          PHIS(1,2)=0.0
          PHIS(2,1)=0.0
          PHIS(2,2)=0.0
          BOUNDARY=0.0
        ELSE IF (X==1 .AND. Y==LENGTH) THEN
          BOUNDARY = (PHIS(2,LENGTH-1) + PHIS(2,LENGTH) + PHIS(1,LENGTH-1))/3.0
        ELSE IF (X==LENGTH .AND. Y==1) THEN
          BOUNDARY = (PHIS(LENGTH-1,2) + PHIS(LENGTH,2) + PHIS(LENGTH-1,1))/3.0
        ELSE IF (X==LENGTH .AND. Y==LENGTH) THEN
          BOUNDARY = (PHIS(LENGTH-1,LENGTH-1) + PHIS(LENGTH-1,LENGTH) + PHIS(LENGTH,LENGTH-1))/3.0
        ELSE IF (X==1) THEN
          NEWVALUE = (1.0/2.0)*(PHIS(X,Y+1)+PHIS(X,Y-1)-H*H*F1(X,Y,H))
          BOUNDARY = (PHIS(X+1,Y) + NEWVALUE)/2.0
        ELSE IF (X==LENGTH) THEN
          NEWVALUE = (1.0/2.0)*(PHIS(X,Y+1)+PHIS(X,Y-1)-H*H*F1(X,Y,H))
          BOUNDARY = (PHIS(X-1,Y) + NEWVALUE)/2.0
        ELSE IF (Y==1) THEN
          NEWVALUE = (1.0/2.0)*(PHIS(X+1,Y)+PHIS(X-1,Y)-H*H*F1(X,Y,H))
          BOUNDARY = (PHIS(X,Y+1) + NEWVALUE)/2.0
        ELSE IF (Y==LENGTH) THEN
          NEWVALUE = (1.0/2.0)*(PHIS(X+1,Y)+PHIS(X-1,Y)-H*H*F1(X,Y,H))
          BOUNDARY = (PHIS(X,Y-1) + NEWVALUE)/2.0
        ELSE
          BOUNDARY = NEWESTIMATE(X,Y,H,PHIS,LENGTH)
        END IF
        RETURN
      END FUNCTION
           
      
      FUNCTION NEWESTIMATE(X,Y,H,PHIS,LENGTH)
        INTEGER X,Y
        REAL H
        REAL NEWESTIMATE
        INTEGER LENGTH
        REAL F1
        REAL :: PHIS(LENGTH,LENGTH)
        !WRITE (*,'(a,g12.4)') "F-PART IS ", H*H*F1(X,Y,H)
        NEWESTIMATE = (1.0/4.0)*(PHIS(X+1,Y)+PHIS(X,Y+1)+PHIS(X-1,Y)+PHIS(X,Y-1)-H*H*F1(X,Y,H))
        !WRITE (*,'(a,g12.4)') "NEWESTIMATE IS ", NEWESTIMATE
        RETURN
      END FUNCTION
      
      FUNCTION HIGHESTCHANGEFUN(OLD,NEW,PREVHIGHEST)
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
        INTEGER SIZE
        INTEGER LENGTH
        REAL H
        REAL :: PHIS(LENGTH,LENGTH)
        REAL ANALYTICAL
        REAL AVGERROR
        INTEGER I,J
        REAL ANALERROR
        AVGERROR = 0.0
        DO I=1,LENGTH
          DO J=1,LENGTH
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
        REAL NEWESTIMATE
        REAL BOUNDARY
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
        
        WRITE (*,*) "(LENGTH-1)*H IS ", ((LENGTH-1)*H)
        
        DO I=1,LENGTH
          DO J=1,LENGTH
            PHIS(I,J) = RAND(0)*10
          END DO
        END DO
        
        NUMITERATIONS = 0
        DO WHILE (ANALERROR(PHIS,LENGTH,H,SIZE)>0.0000001) !((NEWESTIMATE(2,2,H,PHIS,LENGTH)-0)**2>0.00001) !((LOWHIGHESTCHANGE/HIGHESTCHANGE) > 1.0 .OR. (LOWAVGCHANGE/AVGCHANGE) > 1.0) !(NUMITERATIONS < 200000)
          NUMITERATIONS = NUMITERATIONS + 1
          LOWAVGCHANGE = LOWCHANGE(LOWAVGCHANGE, AVGCHANGE)
          LOWHIGHESTCHANGE = LOWCHANGE(LOWHIGHESTCHANGE, HIGHESTCHANGE)
          HIGHESTCHANGE = 0.0
          AVGCHANGE = 0.0
          DO I=1,LENGTH
            DO J=1,LENGTH
                !IF (I==1 .OR. I==LENGTH .OR. J==1 .OR. J==LENGTH) THEN
                NEWVALUE = BOUNDARY(I,J,H,PHIS,LENGTH)
                !ELSE IF (I==2 .OR. I==LENGTH-1 .OR. J==2 .OR. J==LENGTH-1) THEN
                  !NEWVALUE = BOUNDARYNEIGHBOR(I,J,H,PHIS,LENGTH)
                !ELSE
                  !NEWVALUE = NEWESTIMATE(I,J,H,PHIS,LENGTH)
                !END IF
                HIGHESTCHANGE = HIGHESTCHANGEFUN(PHIS(I,J),NEWVALUE,HIGHESTCHANGE)
                AVGCHANGE = AVGCHANGE + ((PHIS(I,J)-NEWVALUE)**2)/SIZE
                PHIS(I,J) = NEWVALUE
            END DO
          END DO
        END DO
        
        
        
        WRITE (*,*) "NUMITERATIONS IS ", NUMITERATIONS
        
        NUMAVGERROR = 0.0
        DO I=1,LENGTH
          DO J=1,LENGTH
          IF (I==1 .OR. I==LENGTH .OR. J==1 .OR. J==LENGTH) THEN
            NUMAVGERROR = NUMAVGERROR + ((PHIS(I,J)-BOUNDARY(I,J,H,PHIS,LENGTH))**2)/SIZE
          ELSE
            NUMAVGERROR = NUMAVGERROR + ((PHIS(I,J)-NEWESTIMATE(I,J,H,PHIS,LENGTH))**2)/SIZE
          END IF
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
        
        open (unit=out_unit,file="results.txt",action="write",status="replace")
        WRITE (OUT_UNIT,'(I4)') LENGTH
        WRITE (OUT_UNIT,'(I4)') LENGTH
        DO I=1,LENGTH
          DO J=1,LENGTH
              WRITE (OUT_UNIT,'(F10.5)') PHIS(I,J)
          END DO
        END DO
        close (out_unit)
      END
