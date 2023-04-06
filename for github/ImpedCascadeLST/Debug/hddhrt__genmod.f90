        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 31 14:45:59 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HDDHRT__genmod
          INTERFACE 
            SUBROUTINE HDDHRT(MM,RH,RD,A,B,H,EPS,X,N,M,HF)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: MM
              REAL(KIND=8) :: RH
              REAL(KIND=8) :: RD
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: H
              REAL(KIND=8) :: EPS
              REAL(KIND=8) :: X(N)
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: HF
              EXTERNAL HF
            END SUBROUTINE HDDHRT
          END INTERFACE 
        END MODULE HDDHRT__genmod
