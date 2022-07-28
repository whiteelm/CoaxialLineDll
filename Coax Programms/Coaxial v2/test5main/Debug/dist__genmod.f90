        !COMPILER-GENERATED INTERFACE MODULE: Wed May 05 20:37:30 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIST__genmod
          INTERFACE 
            FUNCTION DIST(ZZ,KS,N,Z)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ZZ
              INTEGER(KIND=4) :: KS
              COMPLEX(KIND=8) :: Z(N)
              REAL(KIND=8) :: DIST
            END FUNCTION DIST
          END INTERFACE 
        END MODULE DIST__genmod
