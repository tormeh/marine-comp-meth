      
      FUNCTION F1(X,Y,H)
        INTEGER X,Y
        REAL H
        REAL F1
        F1=12-12*X*H-12*Y*H
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
          PHIS(0,1)=0.0
          PHIS(1,0)=0.0
          PHIS(1,1)=0.0
          BOUNDARY=0.0
        ELSE IF (X==1 .AND. Y==LENGTH) THEN
          PHIS(2,LENGTH) = PHIS(2,LENGTH-1)
          PHIS(1,LENGTH-1) = PHIS(2,LENGTH-1)
          BOUNDARY = PHIS(2,LENGTH-1)
        ELSE IF (X==LENGTH .AND. Y==1) THEN
          PHIS(LENGTH-1,1) = PHIS(LENGTH-1,2)
          PHIS(LENGTH,2) = PHIS(LENGTH-1,2)
          BOUNDARY = PHIS(LENGTH-1,2)
        ELSE IF (X==LENGTH .AND. Y==LENGTH) THEN
          PHIS(LENGTH-1,LENGTH) = PHIS(LENGTH-1,LENGTH-1)
          PHIS(LENGTH,LENGTH-1) = PHIS(LENGTH-1,LENGTH-1)
          BOUNDARY = PHIS(LENGTH-1,LENGTH-1)
        ELSE IF (X==1) THEN
          NEWVALUE = (1.0/2.0)*(PHIS(X,Y+1)+PHIS(X,Y-1)-H*H*F1(X,Y,H))
          PHIS(X+1,Y) = NEWVALUE
          BOUNDARY = NEWVALUE
        ELSE IF (X==LENGTH) THEN
          NEWVALUE = (1.0/2.0)*(PHIS(X,Y+1)+PHIS(X,Y-1)-H*H*F1(X,Y,H))
          PHIS(X-1,Y) = NEWVALUE
          BOUNDARY = NEWVALUE
        ELSE IF (Y==1) THEN
          NEWVALUE = (1.0/2.0)*(PHIS(X+1,Y)+PHIS(X-1,Y)-H*H*F1(X,Y,H))
          PHIS(X,Y+1) = NEWVALUE
          BOUNDARY = NEWVALUE
        ELSE IF (Y==LENGTH) THEN
          NEWVALUE = (1.0/2.0)*(PHIS(X+1,Y)+PHIS(X-1,Y)-H*H*F1(X,Y,H))
          PHIS(X,Y-1) = NEWVALUE
          BOUNDARY = NEWVALUE
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
        WRITE (*,'(a,g12.4)') "F-PART IS ", H*H*F1(X,Y,H)
        NEWESTIMATE = (1.0/4.0)*(PHIS(X+1,Y)+PHIS(X,Y+1)+PHIS(X-1,Y)+PHIS(X,Y-1)-H*H*F1(X,Y,H))
        WRITE (*,'(a,g12.4)') "NEWESTIMATE IS ", NEWESTIMATE
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
      
      PROGRAM SOLVER
        REAL, PARAMETER :: H = 0.1
        INTEGER, PARAMETER :: LENGTH = (1.0/H)+1
        REAL :: PHIS(LENGTH, LENGTH)
        INTEGER, PARAMETER :: out_unit=20
        REAL :: HIGHESTCHANGE
        REAL :: NEWVALUE
        REAL NEWESTIMATE
        REAL BOUNDARY
        INTEGER I,J
        REAL HIGHESTCHANGEFUN
        REAL PREVHIGHESTCHANGE
        PREVHIGHESTCHANGE = 20.0
        HIGHESTCHANGE = 10.0
        
        DO I=1,LENGTH
          DO J=1,LENGTH
            PHIS(I,J) = RAND(0)
          END DO
        END DO
        
        DO WHILE ((PREVHIGHESTCHANGE/HIGHESTCHANGE) > 1.0)
          PREVHIGHESTCHANGE = HIGHESTCHANGE
          HIGHESTCHANGE = 0.0
          DO I=1,LENGTH
            DO J=1,LENGTH
                IF (I==1 .OR. I==LENGTH .OR. J==1 .OR. J==LENGTH) THEN
                  NEWVALUE = BOUNDARY(I,J,H,PHIS,LENGTH)
                ELSE
                  NEWVALUE = NEWESTIMATE(I,J,H,PHIS,LENGTH)
                END IF
                HIGHESTCHANGE = HIGHESTCHANGEFUN(PHIS(I,J),NEWVALUE,HIGHESTCHANGE)
                WRITE (*,*) "OLD VALUE IS ", PHIS(I,J)
                PHIS(I,J) = NEWVALUE
                WRITE (*,*) "NEW VALUE IS ", PHIS(I,J)
            END DO
          END DO
        END DO
        
        open (unit=out_unit,file="results.txt",action="write",status="replace")
        DO I=1,LENGTH
          DO J=1,LENGTH
              WRITE (OUT_UNIT,'(F10.5)') PHIS(I,J)
          END DO
        END DO
        close (out_unit)
      END
