      
      
      
      
      FUNCTION BOUNDARY(X,Y,H)
        REAL BOUNDARY
        INTEGER :: X
        INTEGER :: Y
        REAL :: H
        BOUNDARY = (1/4.0)*(X*X*H*H+Y*Y*H*H)
        RETURN
      END FUNCTION
      
      
      
      FUNCTION F(X,Y)
        REAL F
        INTEGER X,Y
        F=1.0
        RETURN
      END FUNCTION
      
      
      FUNCTION NEWESTIMATE(X,Y,H,PHIS,LENGTH)
        INTEGER X,Y
        REAL H
        REAL NEWESTIMATE
        INTEGER LENGTH
        REAL F
        REAL :: PHIS(LENGTH,LENGTH)
        WRITE (*,'(a,g12.4)') "F-PART IS ", H*H*F(X,Y)
        NEWESTIMATE = (1.0/4.0)*(PHIS(X+1,Y)+PHIS(X,Y+1)+PHIS(X-1,Y)+PHIS(X,Y-1)-H*H*F(X,Y))
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
        HIGHESTCHANGE = 0.5
        
        DO I=1,LENGTH
          DO J=1,LENGTH
            PHIS(I,J) = 0.0
          END DO
        END DO
        
        DO WHILE (HIGHESTCHANGE > 0.00015)
          HIGHESTCHANGE = 0.0
          DO I=1,LENGTH
            DO J=1,LENGTH
                IF (I==1 .OR. I==LENGTH .OR. J==1 .OR. J==LENGTH) THEN
                  NEWVALUE = BOUNDARY(I,J,H)
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