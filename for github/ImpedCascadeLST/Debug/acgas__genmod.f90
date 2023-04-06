        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 31 14:45:59 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ACGAS__genmod
          INTERFACE 
            SUBROUTINE ACGAS(AR,AI,N,BR,BI,L,JS)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: AR(N,N)
              REAL(KIND=8) :: AI(N,N)
              REAL(KIND=8) :: BR(N)
              REAL(KIND=8) :: BI(N)
              INTEGER(KIND=4) :: L
              INTEGER(KIND=4) :: JS(N)
            END SUBROUTINE ACGAS
          END INTERFACE 
        END MODULE ACGAS__genmod
