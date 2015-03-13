      PROGRAM SOLVER
      
      REAL FUNCTION F(X,Y)
      F=1.0
      RETURN
      
      REAL FUNCTION PHI(X,Y,H)
      REAL X,Y,H
      PHIESTIMATE = (1/4.0)*(PHI(X+H,Y)+PHI(X,Y+H)+PHI(X-H,Y)+PHI(X,Y-H)-H*H*F(X,Y))
      RETURN
      END
