      PROGRAM MYPROG

      REAL X
      DO
         WRITE(*,*) 'Gi inn et tall (0 avslutter):' 
         READ(*,*)  X 
         IF ( X .EQ. 0 ) EXIT
         WRITE(*,*) 'Tallet er ...........: ', X 
         WRITE(*,*) 'Tallet kvaårert er ..: ', X**2 
         WRITE(*,*) 'Roten av tallet er ..: ', SQRT(X) 
      END DO
      END
