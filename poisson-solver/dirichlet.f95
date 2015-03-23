      
      
      
      
      FUNCTION BOUNDARY(X,Y,H)
        INTEGER :: X
        INTEGER :: Y
        REAL :: H
        BOUNDARY = (1/4.0)*(X*X*H*H+Y*Y*H*H)
        RETURN
      END FUNCTION
      
      
      
      FUNCTION F(X,Y)
        INTEGER X,Y
        F=1.0
        RETURN
      END FUNCTION
      
      
      FUNCTION NEWESTIMATE(X,Y,H,PHIS,LENGTH)
        INTEGER X,Y
        REAL H
        REAL :: PHIS(LENGTH,LENGTH)
        NEWESTIMATE = (1/4.0)*(PHIS(X+1,Y)+PHIS(X,Y+1)+PHIS(X-1,Y)+PHIS(X,Y-1)-H*H*F(X,Y))
        RETURN
      END FUNCTION
      
      FUNCTION SOLVE(LENGTH)
        INTEGER :: LENGTH
      END FUNCTION
      
      PROGRAM SOLVER
        PARAMETER (H = 0.1)
        PARAMETER (LENGTH = (1.0/H)+1)
        PARAMETER (PHISSIZE = (LENGTH*LENGTH-4*LENGTH+2))
        REAL :: PHIS(LENGTH, LENGTH)
        integer, parameter :: out_unit=20
        REAL :: AVGCHANGE
        REAL :: NEWVALUE
        
        DO I=1,LENGTH
          DO J=1,LENGTH
            PHIS(I,J) = BOUNDARY(I,J,H)
          END DO
        END DO
        
        DO WHILE (AVGCHANGE > 0.00001)
          AVGCHANGE = 0.0
          DO I=1,LENGTH
            DO J=1,LENGTH
                IF (I==1 .OR. I==LENGTH .OR. J==1 .OR. j==LENGTH) THEN
                  PHIS(I,J) = BOUNDARY(I,J,H)
                ELSE
                  NEWVALUE = NEWESTIMATE(I,J,H,PHIS,LENGTH)
                  AVGCHANGE = AVGCHANGE + ABS(NEWVALUE-PHIS(I,J))/PHISSIZE
                  PHIS(I,J) = NEWVALUE
                END IF
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
