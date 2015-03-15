      
      
      
      
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
      
      
      FUNCTION NEWESTIMATE(X,Y,H)
        INTEGER X,Y
        REAL H
        NEWESTIMATE = (1/4.0)*(PHI(X+1,Y)+PHI(X,Y+H)+PHI(X-H,Y)+PHI(X,Y-H)-H*H*F(X,Y))
        RETURN
      END FUNCTION
      
      FUNCTION SOLVE(LENGTH)
        INTEGER :: LENGTH
      END FUNCTION
      
      PROGRAM SOLVER
        PARAMETER (H = 0.1)
        PARAMETER (LENGTH = (1.0/H)+1)
        REAL :: PHIS(LENGTH, LENGTH)
        
        DO I=1,LENGTH
          DO J=1,LENGTH
              IF (I==1 .OR. I==LENGTH .OR. J==1 .OR. j==LENGTH) THEN
                PHIS(I,J) = BOUNDARY(I,J,H)
              ELSE
                PHIS(I,J) = NEWESTIMATE(I,J,H)
              END IF
          END DO
        END DO
      END
