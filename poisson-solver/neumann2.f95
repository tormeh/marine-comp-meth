      
      
      
      
      FUNCTION BOUNDARY(X,Y,H,PHIS,LENGTH)
        INTEGER :: X
        INTEGER :: Y
        REAL :: H
        REAL :: PHIS(LENGTH,LENGTH)
        IF (X==0 .AND. Y==0) THEN
          BOUNDARY=0
        ELSE IF (X==0 .AND. Y==LENGTH) THEN
          BOUNDARY = NEWESTIMATE(X,Y,H,PHIS,LENGTH)
        ELSE IF (X==LENGTH .AND. Y==0) THEN
          BOUNDARY = NEWESTIMATE(X,Y,H,PHIS,LENGTH)
        ELSE IF (X==LENGTH .AND. Y==LENGTH) THEN
          BOUNDARY = NEWESTIMATE(X,Y,H,PHIS,LENGTH)
        ELSE IF (X==0) THEN
          BOUNDARY = PHIS(1,Y)
        ELSE IF (X==LENGTH) THEN
          BOUNDARY = PHIS(LENGTH-1,Y)
        ELSE IF (Y==0) THEN
          BOUNDARY = PHIS(X,1)
        ELSE IF (Y==LENGTH) THEN
          BOUNDARY = PHIS(X,LENGTH-1)
        END IF
        RETURN
      END FUNCTION
      
      
      FUNCTION F2(X,Y,H)
        INTEGER X,Y
        REAL H
        F2=(6-12*X*H)*(3*Y*Y*H*H-2*Y*Y*Y*H*H*H) + (3*X*X*H*H-2*X*X*X*H*H*H)*(6-12*Y*H)
        RETURN
      END FUNCTION      
      
      
      FUNCTION NEWESTIMATE(X,Y,H,PHIS,LENGTH)
        INTEGER X,Y
        REAL H
        REAL :: PHIS(LENGTH,LENGTH)
        NEWESTIMATE = (1/4.0)*(PHIS(X+1,Y)+PHIS(X,Y+1)+PHIS(X-1,Y)+PHIS(X,Y-1)-H*H*F2(X,Y,H))
        RETURN
      END FUNCTION
      
      
      PROGRAM SOLVER
        PARAMETER (H = 0.1)
        PARAMETER (LENGTH = (1.0/H)+1)
        PARAMETER (PHISSIZE = (LENGTH*LENGTH-4*LENGTH+2))
        REAL :: PHIS(LENGTH, LENGTH)
        INTEGER, PARAMETER :: out_unit=20
        REAL :: AVGCHANGE
        REAL :: NEWVALUE
        
        DO I=1,LENGTH
          DO J=1,LENGTH
            PHIS(I,J) = BOUNDARY(I,J,H,PHIS,LENGTH)
          END DO
        END DO
        
        DO WHILE (AVGCHANGE > 0.001)
          AVGCHANGE = 0.0
          DO I=1,LENGTH
            DO J=1,LENGTH
                IF (I==1 .OR. I==LENGTH .OR. J==1 .OR. j==LENGTH) THEN
                  PHIS(I,J) = BOUNDARY(I,J,H,PHIS,LENGTH)
                ELSE
                  NEWVALUE = NEWESTIMATE(I,J,H,PHIS,LENGTH)
                  AVGCHANGE = AVGCHANGE + ABS(NEWVALUE-PHIS(I,J))/PHISSIZE
                  PHIS(I,J) = NEWVALUE
                END IF
            END DO
          END DO
        END DO
        
        OPEN (unit=out_unit,file="results.txt",action="write",status="replace")
        DO I=1,LENGTH
          DO J=1,LENGTH
              WRITE (OUT_UNIT,*) PHIS(I,J)
          END DO
        END DO
        CLOSE (out_unit)
      END
