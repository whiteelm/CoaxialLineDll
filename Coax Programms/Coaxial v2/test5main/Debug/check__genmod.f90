        !COMPILER-GENERATED INTERFACE MODULE: Wed May 05 20:37:30 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK__genmod
          INTERFACE 
            SUBROUTINE CHECK(EPS,N,W,BETAM)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: EPS
              COMPLEX(KIND=8) :: W(N)
              REAL(KIND=8) :: BETAM(N)
            END SUBROUTINE CHECK
          END INTERFACE 
        END MODULE CHECK__genmod
