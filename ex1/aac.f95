      PROGRAM MYPROG

      REAL X
      DO
         WRITE(*,*) 'Gi inn et tall (0 avslutter):' 
         READ(*,*)  A
         READ(*,*)  B
         READ(*,*)  N
         REAL, DIMENSION(N) :: ARR
         DO 10 I = 1, N
         ARR(I) = 
   10    CONTINUE
         WRITE(*,*) 'Tallet er ...........: ', X 
         WRITE(*,*) 'Tallet kvaårert er ..: ', X**2 
         WRITE(*,*) 'Roten av tallet er ..: ', SQRT(X) 
      END DO
      END
