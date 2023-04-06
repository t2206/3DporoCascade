        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 31 14:45:59 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FINITE_RADIAL_MODE_EXPANSION__genmod
          INTERFACE 
            SUBROUTINE FINITE_RADIAL_MODE_EXPANSION(M,RH,RD,L,KKMN,B,INF&
     &,KM0N)
              INTEGER(KIND=4) :: L
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: RH
              REAL(KIND=8) :: RD
              REAL(KIND=8) :: KKMN(L)
              REAL(KIND=8) :: B(L,L)
              INTEGER(KIND=4) :: INF
              REAL(KIND=8) :: KM0N(L)
            END SUBROUTINE FINITE_RADIAL_MODE_EXPANSION
          END INTERFACE 
        END MODULE FINITE_RADIAL_MODE_EXPANSION__genmod
