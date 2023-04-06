!This program was created by SHEN, Zihan and was uploaded to Github as was used in the journal paper.
!Anyone who uses this program please refer to the author and the journal paper:
!Shen, Z., Wang, X., Sun, Y., Zhang, G., & Sun, X. (2022). Three-dimensional effects of cascade perforations on rotor–stator interaction noise. Journal of Fluid Mechanics, 952, A7. doi:10.1017/jfm.2022.871
    
    program ImpedCascadeLST
    
    USE omp_lib  !using parallel computation OpenMP
    use lapack95 
    use f95_precision    
    implicit none
    include 'mkl.fi'
    
    real*8,parameter:: pi=3.14159265358979323846d0
    integer,parameter:: NR=25,NZ=21,S_min=-160,S_max=160
    integer,parameter:: nr_s=5*NR,nz_s=9*NZ,NL=NR,INF=15001
    double complex,parameter::ixb=(0.0d0,1.0d0)
    
    integer ii,jj,kk,i,j,k,q,qq,Findex,m,n,ni,nj,nis,njs,KNUM,iil,jjn,Zindex

    ! NR Radial observation points number and expansion order of collocation method   ! NZ Axial points number and expansion order of collocation method 
    !  m=sB-qV q=(S_min,S_max)   
    ! NL  order of Namba's finite radial mode expansion, in this program NR>=NL
    !  nr_s=5*NR; nz_s=5*NZ;   !number of source points in the numerical integration！
    ! in this program nz_s should be (1,3,5,7...)*NZ in order to set observation points in the middle of two source points as in Whitehead     
    
    real*8 Rh,Rd,ChordLength,holeRadius,alphaH,f000,Finterval
    real*8 Mr,c0,rho0,beta_r,bxu,bxd,Mt,num_q
    double complex Pinlet
    integer m_in,n_in,num_B,num_V,qmin,qmax,InterP,Fnum,Zimflag,M_mode,Nsigma,flage,pindirect,num_S
    double complex Conductivity,Zimpedence1,Zimpedence2

    real*8 kmn(S_min:S_max,NL),kmn_pc(S_min:S_max,NR),km0n(NL)  !kmn: approximated eigenvalues using Namba's finite radial mode expansion !kmn_pc: 'precise' eigenvalues obtained from bisection method
    real*8 BB(S_min:S_max,NL,NL),KKmn_m(S_min:S_max,NL),BBB(S_min:S_max,NL,NL)           !Expansion coefficients and Kmn
    real*8 BB_tem(NL,NL),BB_inf(NL,NL),KKmn(NL),KKmn_inf(NL),BBB_inf(NL,NL),kmn_tem(NL)
    real*8 rr(NR),rr_s(nr_s),delta_ro(NR),delta_rs(nr_s)!radial Positions of the observation points and source points

    real*8 zxi_o(NZ),zxi_s(nz_s)     !Axial positions of the observation points and source points in z coordinates
    real*8  delta_zo(NZ),theta_o(NZ),delta_zs(nz_s),theta_s(nz_s),theta_cal(nz_s),zxi_cal(nz_s)
    double complex DeltaP(nz_s,nr_s),P_plus(nz_s,nr_s),P_minus(nz_s,nr_s),VV(NZ*NR),VVV(2*NZ*NR) 
    double complex KernelZplus(NZ*NR,NZ*NR),KernelZminus(NZ*NR,NZ*NR),Kernel(NZ*NR,NZ*NR),Kernel0(NZ*NR,NZ*NR)
    real*8 PSI0(NL,NR),PSI_inf(NL,NR),PSI0s(NL,nr_s),PSI_infs(NL,nr_s),dzmintest(NZ,nz_s),dzmin,dztest   

    !Kernel=zeros(NZ*NR,NZ*NR);      !matrix of the kernel function
    !VV=zeros(NZ*NR,1);            !RHS of the linear equations system
    real*8 f00,omega,k0,alpha3,mtest,test_1,test_2
    real*8 temp_o,zzo,zzs,zspancf(nz_s),zspancf_o(nz_s),exptemp1,exptemp2
    double complex tempinR,tempinS,temp_z,temp1,temp2,temp3,temp_K,temp_K2,V_sir(NZ,NR),V_sit(NZ,NR)!入射声波的径、周向声速度
    real*8 Ar1(NZ*NR,NZ*NR),Ai1(NZ*NR,NZ*NR),Br1(NZ*NR),Bi1(NZ*NR)
    INTEGER acgasflag,testflag,JS1(NZ*NR,NZ*NR)
    double complex knm(S_min:S_max,NR),kzin,Aijplus(NZ*NR),Aijminus(NZ*NR),temp_alpha
    REAL*8 phase_a,phase_b,dp_phase
    REAL*8 SPL_u,SPL_d,SoundEnergy_u,SoundEnergy_d,Eup(S_min:S_max,NL),Edown(S_min:S_max,NL)
    DOUBLE COMPLEX XXX_u,XXX_d,kappanm,LeanFactor1,Pin,Amn_u,Amn_d,alpha_1,alpha_2
    real*8 EnergyPin,Edtotal,Eutotal,Ct,Cr,Ct2,Cr2
    
    real*8 SSS,sum1,sum2
    double complex tempinS2
    
    !condition number
    double complex :: work(2*NZ*NR) ,awork(NZ*NR,NZ*NR)
    double complex :: kernel2(NZ*NR,NZ*NR)
    integer :: lda,infoa,ipiv(NZ*NR)
    real*8 :: rwork(2*NZ*NR),rcond,anormflag,anorm
    real*8 :: kernelabs(NZ*NR,NZ*NR)
    !call zgelss
    integer :: ldaf,nrhs,ldb,ldx

    DOUBLE COMPLEX :: a3(NZ*NR,NZ*NR)
    integer :: rank3
    integer,parameter :: lwork3=2*NR*NZ*NR*NZ
    DOUBLE COMPLEX :: work3(lwork3)
    real*8 :: rwork3(10*NR*NZ)
    real*8 :: s3(NR*NZ)
           
    !testing acgas norm
    real*8 :: anorm2,anorm3
    double complex :: alpha5
    DOUBLE COMPLEX :: VV1(NZ*NR)
    DOUBLE COMPLEX :: VV2(NZ*NR)
    DOUBLE COMPLEX :: VV3(NZ*NR)

    !normalization for comparison
    real*8 p0,rmany,phimany(10001),phimax,kmnmany,Pref
    integer phimaxpos(1)
      
    !comparision @ specific positions
    double complex:: deltaCp(nr_s)
    real*8 ::thetaCp
    
    !Position of hard-soft division line
    real*8 :: ZKzposhw
      
    REAL*8,EXTERNAL::PSI,GAMMA
            
    CALL omp_set_num_threads(20)!cores used in parallel computation
    !CALL omp_set_num_threads(max(NZ,NR)+1)

  
    OPEN(200, FILE='dP_amplitude.plt',status='replace')
    CLOSE(200)      
    OPEN(400, FILE='dP_phase.plt',status='replace')
    CLOSE(400)
    OPEN(1200, FILE='dP_real.plt',status='replace')
    CLOSE(1200)      
    OPEN(1400, FILE='dP_imag.plt',status='replace')
    CLOSE(1400)

    open(300,file='Amn_u.csv',status='replace')
    write(300,*) 'm,n,Amn_u'
    open(301,file='Amn_d.csv',status='replace')
    write(301,*) 'm,n,Amn_d'
    
    open(500,file='DeltaCp.csv',status='replace')
    write(500,*) 'r,Cp'
    
    !m_in=2;         !Circumferential mode of incident sound wave
    !n_in=2;         !Radial mode of incident sound wave
    num_B=16;         ! blade number of rotor
    num_V=24;           ! vane number of stator cascade
    num_q = 3.0d0           !q in theta(r)
    num_S = 1           !sB-qV, not the above q
    
    Rd = 1.0d0
    ChordLength=1.0d0*2.0d0*pi*Rd/dble(num_V)     ! bb Chord length
    !Rd=1.0d0/ChordLength;         ! Rd outer radius 
    Rh=0.5d0*Rd;          ! Rh inner radius of the annular duct
    !ChordLength = 1.0d0
    write(*,*) 'bb',ChordLength
    
    ZKzposhw = 1.5d0*ChordLength    !Position of hard-soft division line
    

    holeRadius=0.001d0;  !hole Radius of the perforations
    alphaH=0.02;!fractional open area of the perforations
    Conductivity = 1.0d0*(-0.00d0,-1.0d0)*2.0d0*holeRadius  !Rayleigh conductivity of the perforations
    !Conductivity=(0.0d0,0.0d0)
    Zimflag=1  !if=1,Perforation Conductivity case
    write(*,*) 'Zimflag=',Zimflag
    write(*,*) 'NR=',NR,'NZ=',NZ
    write(*,*) 'Smax=',S_max,'q=',num_q
    write(*,*)

    Mr=0.5d0;           ! Axial Mach number >0！
    c0=340.d0;          ! c0 speed of sound
    rho0=1.225d0;         ! rou0 density of air
    p0 = rho0*c0*c0/1.4d0         
    beta_r=dsqrt(1d0-Mr**2);
    write (*,*) beta_r
    
    Mt = 0.783d0            ! mt rotor tip mach number as in the benchmark problem of 3rd CAA Workshop
    omega = dble(num_S)*dble(num_B)*Mt*c0/Rd    
    write(*,*) 'omega,O',omega,omega/dble(num_B)

    bxu=-0.5d0*ChordLength
    bxd=1.5d0*ChordLength

    !can be used in sound wave incidence case, not used in vortical wave incidence case
    Pinlet = (0.0d0,0.0d0)
    pindirect = 1       !direction of incident sound wave，if>0,sound wave propagating downstream，if<0,going upstream (from downstream）
        
    write(*,*) 'solving B(m;n,l) and B(m_inf;n,l) and kmn(m,n),kmn(0,l)'
    CALL HQIUKMN(0,NL,Rh,Rd,km0n,1.0d-15)
    InterP = num_S*num_B;
    BB=0.0D0
    kmn_pc=0.0d0
    !启用并行
    CALL OMP_set_dynamic(.true.)            
    !$omp parallel private(M_mode,BB_tem,kmn_tem)
    !$OMP DO
      DO i=S_min,S_max,1
        M_mode=InterP-i*num_V
        write(*,*) 'M_mode=',M_mode
      IF(M_mode==0)THEN   
        DO j=1,NL
            BB(i,j,j)=1.0d0
        END DO
      ELSE
        CALL Finite_Radial_Mode_Expansion(M_mode,Rh,Rd,NL,KKmn,BB_tem,INF,km0n)
        DO j=1,NL
          DO k=1,NL
            BB(i,j,k)=BB_tem(j,k)
          END DO
          KKmn_m(i,j)=KKmn(j)
        END DO
      END IF
    
      IF(M_mode==0)THEN
        DO j=1,NL
          kmn(i,j)=km0n(j)
        END DO
      ELSE
        DO j=1,NL
          kmn(i,j)=DSQRT(KKmn_m(i,j)*(dble(M_mode)**2))
        END DO
      END IF
      if (abs(M_mode)<800) then
          CALL HQIUKMN(M_mode,NR,Rh,Rd,kmn_tem,1.0d-11)!the precision should not be too small to avoid huge computation effort，
          DO j=1,NL
            kmn_pc(i,j)=kmn_tem(j);
          END DO
      endif        
      END DO      
    !$OMP END DO
    !$omp end parallel
    CALL Finite_Radial_Mode_Expansion(INF,Rh,Rd,NL,KKmn_inf,BB_inf,INF-1,km0n)  !  
    do j=1,NL
        KKmn_inf(j) = dsqrt(KKmn_inf(j))
    enddo
    write(*,*) 'kmn_inf',KKmn_inf
    !---------------Coordinates for observation and source positions---------------------------------    
    delta_ro = (Rd - Rh)/dble(NR);
    do ii=1,NR
        rr(ii) = Rh + (Rd - Rh)/(2.0D0*dble(NR))*dble(2*ii-1)
    enddo
    delta_rs = (Rd - Rh)/dble(nr_s)
    do ii=1,nr_s
        rr_s(ii) = Rh + (Rd - Rh)/(2.0D0*dble(nr_s))*dble(2*ii-1)
    enddo

    do ii=1,NZ
        theta_o(ii) = pi*(dble(ii) - 0.5D0)/dble(NZ)
        zxi_o(ii) = ChordLength/2.0d0 - 0.5D0*cos(theta_o(ii))*ChordLength
        delta_zo(ii) = pi/dble(NZ)
    enddo
    do ii=1,nz_s
        theta_s(ii) = pi*(dble(ii) -1.0d0)/dble(nz_s)
        zxi_s(ii) =  ChordLength/2.0d0 - 0.5D0*cos(theta_s(ii))*ChordLength
        delta_zs(ii) = pi/dble(nz_s)
    enddo    
    do ii = 1,NZ
        do jj = 1,nz_s
            dzmintest(ii,jj) = dabs(zxi_o(ii)-zxi_s(jj))
        enddo
    enddo
    dzmin =MinVal(dzmintest)
    write(*,*) dzmin
    !----Solution of efficients BB，CB，DB，PSI(INF,l)-----------------------------
    BBB=0.0D0
      DO i=S_min,S_max,1       
        M_mode=InterP-i*num_V
        DO j=1,NL                            
          DO k=1,NL                    
            DO q=1,NL
              BBB(i,j,k)=BBB(i,j,k)+BB(i,q,k)*BB_inf(q,j)
            END DO
          END DO
        END DO
      END DO
    
     !对应展开法
    do ii=1,NR
        do jj=1,NL
            PSI0(jj,ii) = PSI(0,Rh,Rd,km0n(jj),rr(ii))/dsqrt(GAMMA(0,km0n(jj),Rh,Rd))!归一化的psi0
        enddo
    enddo    
      DO i=1,NR
        DO j=1,NL   !bessel function of infinite order in Namba's
          PSI_inf(j,i)=0.0D0
          DO k=1,NL
            PSI_inf(j,i)= PSI_inf(j,i)+BB_inf(k,j)*PSI0(k,i)
          END DO
        END DO
      END DO
      
    do ii=1,nr_s    !coordinates of recover DeltaPs
        do jj=1,NL
            PSI0s(jj,ii) = PSI(0,Rh,Rd,km0n(jj),rr_s(ii))/dsqrt(GAMMA(0,km0n(jj),Rh,Rd))!归一化的psi0
        enddo
    enddo  
    PSI_infs=0.0d0
    do ii=1,nr_s
        do j=1,NL
            do kk=1,NL
                PSI_infs(j,ii)=PSI_infs(j,ii)+BB_inf(kk,j)*PSI0s(kk,ii)
            enddo
        enddo
    enddo

    do Findex=0,0!Fnum-1
        !f00=f000+Findex*Finterval
        !write(*,*)
        !write(*,*) 'frequency is',f00,'Hz'
        !omega=2*pi*f00;
        k0=omega/c0;
        alpha3 = k0/Mr;
        temp_alpha = -ixb*alphaH*Conductivity/(pi*holeRadius**2*omega*rho0)

        write(*,*) 'Starting to solve Kernel'
        Kernel=(0.0d0,0.0d0)      
        KernelZplus=(0.0d0,0.0d0)
        KernelZminus=(0.0d0,0.0d0)
        DeltaP = (0.0d0,0.0d0)
        P_plus = (0.0d0,0.0d0)
        P_minus = (0.0d0,0.0d0)
        
        temp_o = -dble(num_V)/(2.0d0*pi*rho0*Mr*c0)*ChordLength/2.0D0    !z积分变换后的前置系数 
        knm=(0.0d0,0.0d0)
        do m=S_min,S_max      
            do n=1,NL
                if (k0**2>(beta_r**2*(kmn(m,n)**2))) then
                    knm(m,n) = dsqrt(k0**2 - beta_r**2*(kmn(m,n)**2));
                else
                    knm(m,n) = ixb*dsqrt(beta_r**2*(kmn(m,n)**2) - k0**2);
                endif
            enddo
        enddo
        do ii = -100,100   
            Nsigma = InterP+ii*num_V;
            if (Nsigma > 0) then
                exit
            endif
        enddo
        !write(*,*) 'NSIGMA',Nsigma
        testflag=0

        !CALL OMP_set_dynamic(.true.)            
        !$omp parallel 
        !$OMP DO private(M_mode,zspancf,tempinR,tempinS,KNUM,mtest,test_1,test_2,zzo,zzs,&
        &temp1,temp2,temp3,temp_z,exptemp1,exptemp2,temp_K,temp_K2,qmin,qmax,dztest)
        
        do ni=1,NZ
            do nj=1,NR
                zzo=zxi_o(ni)
                do nis=1,NZ
                    if (nis==1)then
                        do q=1,nz_s
                            zspancf(q)=2.0d0*dcos(theta_s(q)/2.0d0)**2
                        enddo
                    else
                        do q=1,nz_s   
                            zspancf(q)=dsin(dble(nis-1)*theta_s(q))*dsin(theta_s(q))
                        enddo
                    endif
                    do njs=1,NR    
                        tempinR = (0.0d0,0.0d0)
                        tempinS = (0.0d0,0.0d0)
                        mtest = dble(abs(InterP - S_max*num_V))
                        test_1 = (KKmn_inf(njs)/beta_r*dzmin*dble(NZ)/4.0d0)**2
                        test_2 = 1.0d0/mtest*dexp(-KKmn_inf(njs)/beta_r*mtest*dzmin*6.0d0)
                        if (testflag==0.and.test_1<test_2) then
                            write(*,*) 'Maybe it is better to adopt log|z-z| approximation of Namba',nis,ni
                            write(*,*) 'test1=',test_1,'test2=',test_2
                            test_2 = 0.0d0
                            testflag=1
                        endif

                        KNUM = nz_s    
                        do kk = 1,KNUM      
                            zzs=zxi_s(kk)
                            temp_z=cdexp(-ixb*Mr*k0/beta_r**2*(zzo-zzs))                            
                            !Judging the influence of z-z' to reduce the summation number of m！(In order to reduce the computation time)
                            mtest=dble(abs(InterP-S_max/8*num_V))
                            if (dexp(-KKmn_inf(njs)/beta_r*mtest*0.5d0*abs(zzs-zzo))/mtest<1.0d-14) then
                                qmin = S_min/16 -1
                                qmax = S_max/16 +1
                            elseif (dexp(-KKmn_inf(njs)/beta_r*mtest*abs(zzs-zzo))/mtest<1.0d-14) then
                                qmin = S_min/8 -1
                                qmax = S_max/8 +1     
                            elseif (dexp(-KKmn_inf(njs)/beta_r*2.0d0*mtest*abs(zzs-zzo))/2.0d0/mtest<1.0d-13)then
                                qmin = S_min/4 -1
                                qmax = S_max/4 +1
                            elseif (dexp(-KKmn_inf(njs)/beta_r*4.0d0*mtest*abs(zzs-zzo))/4.0d0/mtest<1.0d-11)then
                                qmin = S_min/2 -1
                                qmax = S_max/2 +1
                            else
                                qmin = S_min
                                qmax = S_max
                            endif
                            !qmin = S_min
                            !qmax = S_max
                            if (zzs<zzo) then
                                temp3=(0.0d0,0.0d0)
                                temp1=(0.0d0,0.0d0)
                                do m=qmin,qmax,1    !Summation of m 
                                    M_mode = InterP - m*num_V
                                    temp1=temp1 + 1.0d0*PSI_inf(njs,nj)/rr(nj)*&
                                        &(-dexp(-1.0d0*KKmn_inf(njs)*dble(abs(M_mode))/beta_r*(zzo-zzs))*(-0.5d0/KKmn_inf(njs)**2))                                                              
                                    do iil = 1,NL                         
                                        do jjn = 1,NL 
                                            temp1=temp1+1.0d0*BBB(m,iil,jjn)*BBB(m,njs,jjn)*dble(M_mode)**2*Mr*beta_r**2*PSI_inf(iil,nj)&
                                                &/rr(nj)*cdexp(+ixb*knm(m,jjn)/beta_r**2*(zzo-zzs))/(2.0d0*knm(m,jjn)*(Mr*knm(m,jjn)*1.0d0-k0)) 
                                        enddo
                                    enddo
                                enddo
                                do m=S_min,S_max,1    !对m求和 
                                    M_mode = InterP - m*num_V
                                    !for radially place straight vanes, phi-phi'=0，e**i0=1;
                                    temp3=temp3+1.0d0*PSI_inf(njs,nj)/rr(nj)*cdexp(ixb*alpha3*(zzo-zzs))*(-1.0d0/KKmn_inf(njs)**2)                                                                  
                                    do iil = 1,NL                              
                                        do jjn = 1,NL 
                                            temp3=temp3+1.0d0*BBB(m,iil,jjn)*BBB(m,njs,jjn)*cdexp(ixb*alpha3*(zzo-zzs))&
                                                &*dble(M_mode)**2*Mr**2/(k0**2+Mr**2*kmn(m,jjn)**2)*PSI_inf(iil,nj)/rr(nj)
                                        enddo
                                    enddo
                                enddo
          
                                tempinR = tempinR + (temp1*temp_z + temp3)*zspancf(kk)*delta_zs(kk)

                                exptemp1 = dexp(-KKmn_inf(njs)*dble(Nsigma)/beta_r*(zzo-zzs))
                                exptemp2 = dexp(KKmn_inf(njs)*dble(num_V)/beta_r*(zzo-zzs))
                                tempinS = tempinS - 0.5d0*zspancf(kk)*delta_zs(kk)*temp_z*&
                                    &(exptemp1/(exptemp2-1.0d0) + exptemp1 + 1.0d0/(exptemp1*(exptemp2-1.0d0)))
                            else
                                temp2=(0.0d0,0.0d0)
                                do m=qmin,qmax,1    !对m求和 
                                    M_mode = InterP - m*num_V
                                    temp2=temp2 + 1.0d0*PSI_inf(njs,nj)/rr(nj)*&
                                        &(-dexp(1.0d0*KKmn_inf(njs)*dble(abs(M_mode))/beta_r*(zzo-zzs))*(0.5d0/KKmn_inf(njs)**2))
                                    do iil = 1,NL                         
                                        do jjn = 1,NL  
                                            temp2=temp2+1.0d0*BBB(m,iil,jjn)*BBB(m,njs,jjn)*dble(M_mode)**2*Mr*beta_r**2*PSI_inf(iil,nj)&
                                                &/rr(nj)*cdexp(-ixb*knm(m,jjn)/beta_r**2*(zzo-zzs))/(2.0d0*knm(m,jjn)*(Mr*knm(m,jjn)*(-1.0d0)-k0)) 
                                        enddo
                                    enddo
                                enddo

                                tempinR = tempinR + (temp2)*temp_z*zspancf(kk)*delta_zs(kk)
                                exptemp1 = dexp(KKmn_inf(njs)*dble(Nsigma)/beta_r*(zzo-zzs))
                                exptemp2 = dexp(-KKmn_inf(njs)*dble(num_V)/beta_r*(zzo-zzs))
                                tempinS = tempinS + 0.5d0*zspancf(kk)*delta_zs(kk)*temp_z*&
                                    &(exptemp1/(exptemp2-1.0d0) + exptemp1 + 1.0d0/(exptemp1*(exptemp2-1.0d0)))
                            endif 
                            !trapezoidal rule ____difference is small from 0 order integration
                            if (kk == 1)then
                                tempinS=tempinS/2
                                tempinR=tempinR/2
                            endif
                        enddo
                        temp_K = temp_o*tempinR*KKmn_inf(njs)**2                     
                        temp_K2 = temp_o*tempinS*PSI_inf(njs,nj)/rr(nj)
                        !$omp atomic seq_cst,write,seq_cst
                        Kernel((nj-1)*NZ+ni,(njs-1)*NZ+nis) = temp_K+temp_K2 
                        !$omp end atomic
                    enddo
                enddo
            enddo
        enddo
        !$OMP END DO
        !$omp end parallel

        !the incident disturbance velocity as the RHS of the linear equations system
        do ii=1,NZ
            do jj=1,NR
                zzo=zxi_o(ii)
                VV((jj-1)*NZ+ii)=-1.0d0*0.1d0*Mr*c0*cdexp(ixb*dble(num_S)*dble(num_B)*&
                    &(omega/dble(num_S)/dble(num_B)*zzo/Mr/c0+0.0d0+2.0d0*pi*dble(num_q)/dble(num_B)*(rr(jj)-Rh)/(Rd-Rh) ))!
            enddo
        enddo 
        write(*,*) 'VV',VV,'omega',omega

        !------------------------------------------------------
        if (Zimflag == 1) then  
            
            do nis = 1,NZ
                if (nis==1)then
                    do q=1,NZ
                        zspancf_o(q)=dcos(theta_o(q)/2.0d0)/dsin(theta_o(q)/2.0d0)
                    enddo
                else
                    do q=1,NZ   
                        zspancf_o(q)=dsin(dble(nis-1)*theta_o(q))
                    enddo
                endif
                do njs = 1,NR
                    do ni=1,NZ
                        do nj=1,NR
                            zzo = zxi_o(ni)
                            if (zzo<ZKzposhw) then                           
                                KernelZplus((nj-1)*NZ+ni,(njs-1)*NZ+nis)=PSI_inf(njs,nj)*zspancf_o(ni)*temp_alpha
                            endif
                            !KernelZplus((nj-1)*NZ+ni,(njs-1)*NZ+nis)=PSI_inf(njs,nj)*zspancf_o(ni)*temp_alpha
                        enddo
                    enddo
                enddo
            enddo
            
            do i=1,NZ*NR
                do j=1,NZ*NR
                    !$omp atomic seq_cst,write,seq_cst
                    Kernel0(i,j)=Kernel(i,j)-KernelZplus(i,j)
                    !$omp end atomic
                enddo
            enddo
            
            
        DO i=1,NZ*NR
            DO j=1,NZ*NR
                Ar1(i,j)=DREAL(Kernel0(i,j))
                Ai1(i,j)=DIMAG(Kernel0(i,j))
            END DO
                Br1(i)=DREAL(VV(i))
                Bi1(i)=DIMAG(VV(i))
        END DO
        CALL ACGAS(Ar1,Ai1,NZ*NR,Br1,Bi1,acgasflag,JS1)
        IF(acgasflag==0)THEN
            WRITE(*,*) 'flag=',acgasflag,'We got an acgas problem!!'
        END IF
        WRITE(*,*) 'It is over to solve equation.'
        
        !testing___________
        VV3 = VV
        lda = NZ*NR
        anorm = 0
        do i=1,NR*NZ
            do j=1,NR*NZ
                kernel2(i,j) = Kernel0(i,j)
                kernelabs(i,j) = abs(kernel2(i,j))
            enddo
        enddo

        do j = 1,NZ*NR
            anormflag = sum(kernelabs(:,j))
            if (anormflag .GT. anorm) then
                anorm = anormflag
            end if
        end do
        
        write(*,*) 'norm of Matrix',anorm
        do i=1,NR*NZ
            do j=1,NR*NZ
                awork(i,j) = Kernel0(i,j)
            enddo
        enddo
        
        call zgetrf(NZ*NR,NZ*NR,awork,lda,ipiv,infoa)
        call zgecon('1',NZ*NR,awork,lda,anorm,rcond,work,rwork,infoa)
        write(*,*) 'infoa', infoa
        write(*,*) 'rcond',rcond
        !estimation the condition number of the kernel matrix
            !testing mkl zeglss SVD
        VV2 = VV
        infoa = 0
        rcond = 0.0d0
        !work3=(0.0d0,0.0d0)
        !ldaf = NZ*NR
        lda = NZ*NR
        ldb = NZ*NR
        !ldx = NZ*NR
        do i=1,NR*NZ
            do j=1,NR*NZ
                a3(i,j) = Kernel0(i,j)
            enddo
        enddo
        write(*,*) 'a3',a3(1,1)
        !write(*,*) 'boom'
        nrhs = 1        
        write(*,*) 'lda=',lda
        call zgelss(NZ*NR,NZ*NR,nrhs,a3,lda,VV,ldb,s3,rcond,rank3,work3,lwork3,rwork3,infoa)
        write(*,*) 'info3', infoa
        write(*,*) 'optimun lwork',work3(1)
      
        alpha5 = DCMPLX(1.0d0,0.0d0)      
        VV1 = VV
        call zgemv('n',NZ*NR,NZ*NR,alpha5,kernel2,lda,VV1,1,-alpha5,VV2,1)
        anorm2 = dznrm2(NZ*NR,VV2,1) 
        write(*,*) 'norm of svd VV2',anorm2
          
        !!testing acgas
        DO i=1,NZ*NR
            VV(i)=DCMPLX(Br1(i),Bi1(i))   !VV stores the solution of DeltaPs hereafter.
        END DO
      
        !testing acgas  
        alpha5 = DCMPLX(1.0d0,0.0d0)
        call zgemv('n',NZ*NR,NZ*NR,alpha5,kernel2,lda,VV,1,-alpha5,VV3,1) 
        anorm3 = dznrm2(NZ*NR,VV3,1) 
        write(*,*) 'norm of VV3',anorm3
      
        !use zeglss result
        if (anorm3 .GT. anorm2) then
            VV = VV1
            write(*,*) 'acgas is worse than zgelss'
        end if
        
        Aijplus=0.5d0*VV
        Aijminus=-0.5d0*VV
        
        do q=1,nz_s
            theta_cal(q)=theta_s(q) + pi/dble(nz_s)/2.0d0     
            zxi_cal(q)=ChordLength/2.0d0 - 0.5D0*dcos(theta_cal(q))*ChordLength
        enddo        
        do ni = 1,nz_s
            do nj = 1,nr_s
                zzs = zxi_cal(ni)
                do nis = 1,NZ    !i'，j'
                  if (nis ==1) then
                    do njs = 1,NR
                        P_plus(ni,nj)=P_plus(ni,nj)+Aijplus((njs-1)*NZ+nis)/tan(theta_cal(ni)/2.0d0)*PSI_infs(njs,nj);
                        P_minus(ni,nj)=P_minus(ni,nj)+Aijminus((njs-1)*NZ+nis)/tan(theta_cal(ni)/2.0d0)*PSI_infs(njs,nj);
                    enddo
                  else
                    do njs = 1,NR
                        P_plus(ni,nj)=P_plus(ni,nj)+Aijplus((njs-1)*NZ+nis)*dsin(dble(nis-1)*theta_cal(ni))*PSI_infs(njs,nj);
                        P_minus(ni,nj)=P_minus(ni,nj)+Aijminus((njs-1)*NZ+nis)*dsin(dble(nis-1)*theta_cal(ni))*PSI_infs(njs,nj);
                    enddo
                  endif
                enddo
            enddo
        enddo
        DeltaP = P_plus-P_minus
            
        else
            
        endif
        
        !=============Output of the amplitude and phase of the unsteady loading DeltaPs
        !WRITE (*,*) 'Output of the amplitude and phase of the unsteady loading DeltaPs'
        Pref = 0.5d0*rho0*Mr**2*c0**2
        deltaCp=0.0d0
        thetaCp=dacos(1-2.0d0*0.06d0)
        write(500,*) 'x/c=0.06'
        do nj=1,nr_s
            do nis = 1,NZ    !i'，j'
                if (nis == 1) then
                do njs = 1,NR
                    deltaCp(nj)=deltaCp(nj)+Aijplus((njs-1)*NZ+nis)/tan(thetaCp/2.0d0)*PSI_infs(njs,nj)
                enddo
                else
                do njs = 1,NR
                    deltaCp(nj)=deltaCp(nj)+Aijplus((njs-1)*NZ+nis)*dsin(dble(nis-1)*thetaCp)*PSI_infs(njs,nj)
                enddo
                endif
            enddo
            write(500,*) rr_s(nj),',',2*real(deltaCp(nj))/Pref,',',2*imag(deltaCp(nj))/Pref
        enddo
        
        deltaCp=0.0d0
        thetaCp=dacos(1-2.0d0*0.2d0)
        write(500,*) 'x/c=0.2'
        do nj=1,nr_s
            do nis = 1,NZ    !i'，j'
                if (nis == 1) then
                do njs = 1,NR
                    deltaCp(nj)=deltaCp(nj)+Aijplus((njs-1)*NZ+nis)/tan(thetaCp/2.0d0)*PSI_infs(njs,nj)
                enddo
                else
                do njs = 1,NR
                    deltaCp(nj)=deltaCp(nj)+Aijplus((njs-1)*NZ+nis)*dsin(dble(nis-1)*thetaCp)*PSI_infs(njs,nj)
                enddo
                endif
            enddo
            write(500,*) rr_s(nj),',',2*real(deltaCp(nj))/Pref,',',2*imag(deltaCp(nj))/Pref
        enddo
        
        deltaCp=0.0d0
        thetaCp=dacos(1-2.0d0*0.5d0)
        write(500,*) 'x/c=0.5'
        do nj=1,nr_s
            do nis = 1,NZ    !i'，j'
                if (nis == 1) then
                do njs = 1,NR
                    deltaCp(nj)=deltaCp(nj)+Aijplus((njs-1)*NZ+nis)/tan(thetaCp/2.0d0)*PSI_infs(njs,nj)
                enddo
                else
                do njs = 1,NR
                    deltaCp(nj)=deltaCp(nj)+Aijplus((njs-1)*NZ+nis)*dsin(dble(nis-1)*thetaCp)*PSI_infs(njs,nj)
                enddo
                endif
            enddo
            write(500,*) rr_s(nj),',',2*real(deltaCp(nj))/Pref,',',2*imag(deltaCp(nj))/Pref
        enddo
        
        deltaCp=0.0d0
        thetaCp=dacos(1-2.0d0*0.9d0)
        write(500,*) 'x/c=0.9'
        do nj=1,nr_s
            do nis = 1,NZ    !i'，j'
                if (nis == 1) then
                do njs = 1,NR
                    deltaCp(nj)=deltaCp(nj)+Aijplus((njs-1)*NZ+nis)/tan(thetaCp/2.0d0)*PSI_infs(njs,nj)
                enddo
                else
                do njs = 1,NR
                    deltaCp(nj)=deltaCp(nj)+Aijplus((njs-1)*NZ+nis)*dsin(dble(nis-1)*thetaCp)*PSI_infs(njs,nj)
                enddo
                endif
            enddo
            write(500,*) rr_s(nj),',',2*real(deltaCp(nj))/Pref,',',2*imag(deltaCp(nj))/Pref
        enddo
        
        
        OPEN(200,position='Append',FILE='dP_amplitude.plt')
        WRITE(200,*) 'title="unsteady surface pressure distribution_amplitude"'
        WRITE(200,*) 'variables=z,r,dP_amplitude'
        WRITE(200,*) 'zone I=,', nz_s, ',J=,',nr_s, ',f=POINT'      
        DO nj=1,nr_s
        DO ni=1,nz_s
            zzs=zxi_cal(ni)
            WRITE(200,*) zzs, rr_s(nj),ABS(DeltaP(ni,nj))/Pref
        END DO
        END DO 
        CLOSE(200) 
        
        OPEN(1200,position='Append',FILE='dP_real.plt')
        WRITE(1200,*) 'title="unsteady surface pressure distribution_real part"'
        WRITE(1200,*) 'variables=z,r,dP_real'
        WRITE(1200,*) 'zone I=,', nz_s, ',J=,',nr_s, ',f=POINT'      
        DO nj=1,nr_s
        DO ni=1,nz_s
            zzs=zxi_cal(ni)
            WRITE(1200,*) zzs, rr_s(nj),real(DeltaP(ni,nj))/Pref
        END DO
        END DO 
        CLOSE(1200)
      
        OPEN(400,position='Append',FILE='dP_phase.plt')
        WRITE(400,*) 'title="unsteady surface pressure distribution_phase"'
        WRITE(400,*) 'variables=z,r,dP_phase'
        WRITE(400,*) 'zone I=,', nz_s, ',J=,',nr_s, ',f=POINT'      
        DO nj=1,nr_s
        DO ni=1,nz_s
            zzs=zxi_cal(ni)          
            phase_a=DREAL(DeltaP(ni,nj))
            phase_b=DIMAG(DeltaP(ni,nj))
            if(phase_a>0.0D0 .and. phase_b>0.0D0)then
                dp_phase=ATAN(DABS(phase_b/phase_a))
            else if(phase_a<0.0D0 .and. phase_b>0.0D0)then
                dp_phase=PI-ATAN(DABS(phase_b/phase_a))
            else if(phase_a<0.0D0 .and. phase_b<0.0D0)then
                dp_phase=PI+ATAN(DABS(phase_b/phase_a))
            else if(phase_a>0.0D0 .and. phase_b<0.0D0)then
                dp_phase=2.0D0*PI-ATAN(DABS(phase_b/phase_a))
            else          
            end if
            WRITE(400,*) zzs, rr_s(nj),dp_phase       
        END DO
        END DO 
        CLOSE(400)
        
        OPEN(1400,position='Append',FILE='dP_imag.plt')
        WRITE(1400,*) 'title="unsteady surface pressure distribution_imaginary part"'
        WRITE(1400,*) 'variables=z,r,dP_imag'
        WRITE(1400,*) 'zone I=,', nz_s, ',J=,',nr_s, ',f=POINT'      
        DO nj=1,nr_s
        DO ni=1,nz_s
            zzs=zxi_cal(ni)          
            WRITE(1400,*) zzs, rr_s(nj),DIMAG(DeltaP(ni,nj))/Pref    
        END DO
        END DO 
        CLOSE(1400)
        
        
        SoundEnergy_u=0.0d0
        SoundEnergy_d=0.0d0
        do m=S_min,S_max,1
            M_mode=InterP-m*num_V            
            if (M_mode.ne.0 .and. abs(M_mode)<800)then
                do n=1,NL
                if (k0**2>(beta_r**2*(kmn_pc(m,n)**2)))then 
                    write(*,*) 'cut-on,M_mode,n=',M_mode,n
                    kappanm = dsqrt(k0**2-beta_r**2*(kmn_pc(m,n)**2))
                    alpha_1=(-Mr*k0-dsqrt(k0**2-beta_r**2*(kmn_pc(m,n)**2)))/beta_r**2
                    alpha_2=(-Mr*k0+dsqrt(k0**2-beta_r**2*(kmn_pc(m,n)**2)))/beta_r**2
                    XXX_u=(0.0d0,0.0d0)
                    XXX_d=(0.0d0,0.0d0)
                    do nis=1,nz_s
                        do njs=1,nr_s
                            zzs=zxi_cal(nis);
                            LeanFactor1 = -ixb*dble(M_mode)/rr_s(njs)*PSI(M_mode,Rh,Rd,kmn_pc(m,n),rr_s(njs))&
                                &/dsqrt(GAMMA(M_mode,kmn_pc(m,n),Rh,Rd))
                            XXX_u=XXX_u+LeanFactor1*cdexp(ixb*alpha_1*(-zzs))*DeltaP(nis,njs)&
                                &*ChordLength/2.0D0*delta_rs(njs)*delta_zs(nis)*dsin(theta_cal(nis))
                            XXX_d=XXX_d+LeanFactor1*cdexp(ixb*alpha_2*(-zzs))*DeltaP(nis,njs)&
                                &*ChordLength/2.0D0*delta_rs(njs)*delta_zs(nis)*dsin(theta_cal(nis))
                        enddo
                    enddo
                    SoundEnergy_u=SoundEnergy_u+dble(num_V)**2*omega/(2.0D0*pi*rho0*Mr**2*c0**2)&
                        &/(4.0D0*kappanm)*(beta_r**2*Mr/(-Mr*kappanm-k0))**2*(XXX_u*conjg(XXX_u))
                    SoundEnergy_d=SoundEnergy_d+dble(num_V)**2*omega/(2.0D0*pi*rho0*Mr**2*c0**2)&
                        &/(4.0D0*kappanm)*(beta_r**2*Mr/(Mr*kappanm-k0))**2*(XXX_d*conjg(XXX_d))
                endif
                enddo
            endif
        enddo
        !之前少了1/2！！！
        SoundEnergy_u = 0.5d0*SoundEnergy_u
        SoundEnergy_d = 0.5d0*SoundEnergy_d
        SPL_u=10.0D0*DLOG10(SoundEnergy_u*1.0D12/(0.5d0*(Rd**2-Rh**2)))
        SPL_d=10.0D0*DLOG10(SoundEnergy_d*1.0D12/(0.5d0*(Rd**2-Rh**2)))
        write(*,*) '============================================='
        write(*,*) 'SPL,Energy'
        WRITE(*,*)  SPL_u,SoundEnergy_u
        WRITE(*,*)  SPL_d,SoundEnergy_d
        write(*,*) '============================================='        

        open(2100,file='SPLresults.txt',status='replace')
        write(2100,*) 's=q=',num_S,num_q
        WRITE(2100,*)  SPL_u,SoundEnergy_u
        WRITE(2100,*)  SPL_d,SoundEnergy_d
        close(2100)
      
        write(*,*) '——————Amplitudes and Energy for each mode——————'        
        flage=0
        Eup = 0.0d0
        Edown = 0.0d0   
        do m=S_min,S_max,1
            M_mode=InterP-m*num_V
            if (abs(M_mode)<800) then
            do n=1,NL
                if(k0**2>beta_r**2*kmn_pc(m,n)**2)then
                    kappanm=dsqrt(k0**2.0d0-beta_r**2*(kmn_pc(m,n)**2))
                    flage=0
                    !write(*,*)  'M_mode, n = %d,%d',M_mode,n
                else
                    kappanm=ixb*dsqrt(beta_r**2*(kmn_pc(m,n)**2.0d0)-k0**2)
                    flage=1
                endif
                alpha_1 = (-Mr*k0-kappanm)/beta_r**2
                alpha_2 = (-Mr*k0+kappanm)/beta_r**2
                Amn_u=(0.0d0,0.0d0)
                Amn_d=(0.0d0,0.0d0)
                do nj=1,nr_s
                LeanFactor1=-ixb*dble(M_mode)/rr_s(nj)*PSI(M_mode,Rh,Rd,kmn_pc(m,n),rr_s(nj))/dsqrt(GAMMA(M_mode,kmn_pc(m,n),Rh,Rd))
                    do ni=1,nz_s
                        zzs=zxi_cal(ni)
                        Amn_u=Amn_u+LeanFactor1*1.0D0/kappanm*dble(num_V)/4.0D0/pi*cdexp(ixb*alpha_1*(bxu-zzs))&
                            &*DeltaP(ni,nj)*sin(theta_cal(ni))*ChordLength/2.0D0*delta_rs(nj)*delta_zs(ni)*(ixb)
                        Amn_d=Amn_d+LeanFactor1*1.0D0/kappanm*dble(num_V)/4.0D0/pi*cdexp(ixb*alpha_2*(bxd-zzs))&
                            &*DeltaP(ni,nj)*sin(theta_cal(ni))*ChordLength/2.0D0*delta_rs(nj)*delta_zs(ni)*(ixb)
                    enddo
                enddo
                if ((m==0 .or. m==1) .and. n<6 ) then
                    !normalization for namba comparison
                    kmnmany = kmn_pc(m,n)
                    do k=0,10000
                        rmany = dble(k)*(Rd-Rh)/10000.0d0 + Rh                        
                        phimany(k+1) = abs(PSI(M_mode,Rh,Rd,kmnmany,rmany)/dsqrt(GAMMA(M_mode,kmnmany,Rh,Rd)))!归一化的psi0
                    enddo
                    phimaxpos = maxloc(phimany)
                    !phimax = maxval(phimany)
                    rmany = dble(phimaxpos(1))*(Rd-Rh)/10000.0d0 + Rh
                    phimax = PSI(M_mode,Rh,Rd,kmnmany,rmany)/dsqrt(GAMMA(M_mode,kmnmany,Rh,Rd))
                    write(*,*) 'mn=',M_mode,n,'phimax=',phimax
                    write(*,*) 'Amn-ud',Amn_u/p0*phimax,Amn_d/p0*phimax
                    write(300,*) m,',',n,',',real(Amn_u/p0*phimax),',',imag(Amn_u/p0*phimax)
                    write(301,*) m,',',n,',',real(Amn_d/p0*phimax),',',imag(Amn_d/p0*phimax)
                endif
                                
                if (flage==0)then
                    if (M_mode==m_in.and.n==n_in)then
                        if (pindirect>0) then                            
                            Pin=Pinlet*cdexp(ixb*alpha_2*(bxd));
                            Amn_d = Pin+Amn_d;
                            EnergyPin=0.5d0*abs(Pin)**2*beta_r**4*omega*kappanm*2.0d0*pi/(rho0*c0**2)/(kappanm*Mr-k0)**2;
                        else
                            Pin=Pinlet*cdexp(ixb*alpha_1*(bxu));
                            Amn_u = Pin+Amn_u;
                            EnergyPin=0.5d0*abs(Pin)**2*beta_r**4*omega*kappanm*2.0d0*pi/(rho0*c0**2)/(kappanm*Mr+k0)**2;
                        endif
                        !之前少了1/2！！！
                        Eup(m,n)=0.5d0*abs(Amn_u)**2*beta_r**4*omega*kappanm*2.0d0*pi/(rho0*c0**2)/(kappanm*Mr+k0)**2;
                        Edown(m,n)=0.5d0*abs(Amn_d)**2*beta_r**4*omega*kappanm*2.0d0*pi/(rho0*c0**2)/(kappanm*Mr-k0)**2;                        
                    else
                        Eup(m,n)=0.5d0*abs(Amn_u)**2*beta_r**4*omega*kappanm*2.0d0*pi/(rho0*c0**2)/(kappanm*Mr+k0)**2;
                        Edown(m,n)=0.5d0*abs(Amn_d)**2*beta_r**4*omega*kappanm*2.0d0*pi/(rho0*c0**2)/(kappanm*Mr-k0)**2;
                    endif
                endif
            enddo
            endif
        enddo
        Eutotal=sum(Eup)
        Edtotal=sum(Edown)
        write(*,*) 'EuEd',Eutotal,Edtotal

    enddo
    close(300)
    close(301)


    print *, 'Hello World'
    read(*,*) 

    end program ImpedCascadeLST
    
    
    !-----------------------------------------------------------------------
!-----Solving the eigenvalue "KKmn==KKmn**2.0D0" and eigenvector "B(n,l,m)"
!-----------------------------------------------------------------------
      SUBROUTINE Finite_Radial_Mode_Expansion(m,Rh,Rd,L,KKmn,B,INF,km0n)
      IMPLICIT NONE
      INTEGER::I_r=5000       
      REAL*8::EPS=1.0d-12      !Precision for the solution of eigenvalues 
      INTEGER m,L,INF
      INTEGER i,j,ir
      INTEGER NUM(L)          
      REAL*8 Rh,Rd,km0n(L)     
      REAL*8 rr
      REAL*8 KKmn(L),B(L,L),VV(L,L)         !eigenvector(m,n:l)
      REAL*8,EXTERNAL::PSI,GAMMA

      B=0.0d0

      DO i=1,L
        DO j=1,L
          DO ir=1,2*I_r-1,2
            rr=Rh+ir*(Rd-Rh)/(2*I_r)
            B(i,j)=B(i,j)+ ((Rd-Rh)/I_r)/rr*PSI(0,Rh,Rd,km0n(i),rr)/DSQRT(GAMMA(0,km0n(i),Rh,Rd))&
                &*PSI(0,Rh,Rd,km0n(j),rr)/DSQRT(GAMMA(0,km0n(j),Rh,Rd))
          END DO
        END DO
      END DO
      
      DO i=1,L
        IF(m<INF)THEN       
          B(i,i)=B(i,i)+(km0n(i)/m)**2.0d0
        END IF
      END DO
      
      CALL CJCBJ(B,L,EPS,VV)

!     sorting the eigenvectors form smaller ones to larger ones
      NUM=1
      KKmn(1)=B(1,1)
      DO j=1,L
        DO i=1,L
          IF(B(num(j),num(j))>=B(i,i))THEN
            NUM(j)=i
            KKmn(j)=B(i,i)
          ELSE
          END IF
        END DO
        B(NUM(j),NUM(j))=1d10
      END DO

      DO i=1,L
        DO j=1,L
          B(i,j)=VV(i,num(j))
        END DO
      END DO
      
      END SUBROUTINE Finite_Radial_Mode_Expansion

!------------------------------------------------------------
!-----Jacobi method for the solution of the eigenvalues and eigenvectors of a matrix
!------------------------------------------------------------
      SUBROUTINE CJCBJ(A,N,EPS,V)
      DIMENSION A(N,N),V(N,N)
      DOUBLE PRECISION A,V,FF,FM,CN,SN,OMEGA,X,Y,EPS
      INTEGER P,Q
      DO 20 I=1,N
        V(I,I)=1.0
        DO 10 J=1,N
          IF (I.NE.J) V(I,J)=0.0
10      CONTINUE
20    CONTINUE
      FF=0.0
      DO 500 I=2,N
      DO 500 J=1,I-1
500     FF=FF+A(I,J)*A(I,J)
        FF=SQRT(2.0*FF)
205     FF=FF/(1.0*N)
25    DO 30 I=2,N
      DO 30 J=1,I-1
      IF (ABS(A(I,J)).GE.FF) THEN
        P=I
        Q=J
        GOTO 600
      END IF
30    CONTINUE
      IF (FF.GE.EPS) GOTO 205
      RETURN
600   X=-A(P,Q)
      Y=(A(Q,Q)-A(P,P))/2.0
      OMEGA=X/SQRT(X*X+Y*Y)
      IF (Y.LT.0.0) OMEGA=-OMEGA
      SN=1.0+SQRT(1.0-OMEGA*OMEGA)
      SN=OMEGA/SQRT(2.0*SN)
      CN=SQRT(1.0-SN*SN)
      FM=A(P,P)
      A(P,P)=FM*CN*CN+A(Q,Q)*SN*SN+A(P,Q)*OMEGA
      A(Q,Q)=FM*SN*SN+A(Q,Q)*CN*CN-A(P,Q)*OMEGA
      A(P,Q)=0.0
      A(Q,P)=0.0
      DO 60 J=1,N
        IF ((J.NE.P).AND.(J.NE.Q)) THEN
          FM=A(P,J)
          A(P,J)=FM*CN+A(Q,J)*SN
          A(Q,J)=-FM*SN+A(Q,J)*CN
        END IF
60    CONTINUE

      DO 70 I=1,N
        IF ((I.NE.P).AND.(I.NE.Q)) THEN
          FM=A(I,P)
          A(I,P)=FM*CN+A(I,Q)*SN
          A(I,Q)=-FM*SN+A(I,Q)*CN
        END IF
70    CONTINUE
      DO 80 I=1,N
        FM=V(I,P)
        V(I,P)=FM*CN+V(I,Q)*SN
        V(I,Q)=-FM*SN+V(I,Q)*CN
80    CONTINUE
      GOTO 25
      END
      
!-----------------------------------------------------------------------
!-----CALCULATING PSI(k_mn r)/sqrt(GAMMA): The radially normalized characteristic function 
!-----------------------------------------------------------------------
      FUNCTION GAMMA(M,Kmn,Rh,Rd)
      INTEGER M		
      REAL*8 Rh,Rd,Kmn
      REAL*8 GAMMA
      
      REAL*8,EXTERNAL::PSI
      
      IF(Rh==0.0d0)THEN
        IF(M==0)THEN
          GAMMA=0.5d0*(Rd**2.0d0)*PSI(M,Rh,Rd,Kmn,Rd)**2.0D0
        ELSE
          GAMMA=0.5D0*(Rd**2.0D0-(M**2.0D0)/(Kmn**2.0D0))*(PSI(M,Rh,Rd,Kmn,Rd)**2.0D0)
        END IF
      ELSE
        IF(M==0)THEN
          GAMMA=0.5d0*((Rd**2.0d0)*PSI(M,Rh,Rd,Kmn,Rd)**2.0D0-(Rh**2.0d0)*PSI(M,Rh,Rd,Kmn,Rh)**2.0D0)
          ELSE
          GAMMA=0.5D0*(Rd**2.0D0-(M**2.0D0)/(Kmn**2.0D0))*(PSI(M,Rh,Rd,Kmn,Rd)**2.0D0)-&
              &0.5D0*(Rh**2.0D0-(M**2.0D0)/(Kmn**2.0D0))*(PSI(M,Rh,Rd,Kmn,Rh)**2.0D0)
        END IF
      END IF
      
      RETURN
      END FUNCTION GAMMA    
      
!--------------------------------------------------------------
      FUNCTION PSI(M,Rh,Rd,Kmn,X)
      INTEGER M
      REAL*8 Rh,Rd,X,Kmn,A,B 
      REAL*8 PSI
      REAL*8,EXTERNAL::MBSL1,MBSL2
      
      IF(Rh==0.0D0)THEN
        PSI=MBSL1(M,Kmn*X)
      ELSE
        IF(M==0)THEN
          A=-MBSL2(1,Kmn*Rh)
          B=-MBSL1(1,Kmn*Rh)
          PSI=MBSL1(M,Kmn*X)-B/A*MBSL2(M,Kmn*X)
        ELSE
          A=0.5D0*(MBSL2(M-1,Kmn*Rd)-MBSL2(M+1,Kmn*Rd))
          B=0.5D0*(MBSL1(M-1,Kmn*Rd)-MBSL1(M+1,Kmn*Rd))
          PSI=MBSL1(M,Kmn*X)-B/A*MBSL2(M,Kmn*X)
        END IF
      END IF
      RETURN
      END FUNCTION PSI      

!-----------------------------------------------------------------------
!-----CALCULATING ROOT OF CHARACTERISTIC FUNCTION 
!-----------------------------------------------------------------------
      SUBROUTINE HQIUKMN(M,Nm,Rh,Rd,X,EPS)
      EXTERNAL hF
      INTEGER M,Nm
      DOUBLE PRECISION Rh,Rd,X(Nm),A,B,H,EPS,hF
      A=ABS(M)*0.9D0/Rd
      B=100000000.0d0  !B should be large since in our case Smax Smin is larger than Namba's LST
      H=1.0d0
      CALL hDDHRT(M,Rh,Rd,A,B,H,EPS,X,Nm,N,hF)
      END SUBROUTINE HQIUKMN
         
!-----------------------------------------------------------------------
!-----THE BISECTION METHOD OF FINDING ROOT OF EQUATION(STANDARD)
!-----(roots of hF(X(M))=0 )
!-----------------------------------------------------------------------
      SUBROUTINE hDDHRT(MM,Rh,Rd,A,B,H,EPS,X,N,M,hF)
      DOUBLE PRECISION A,B,H,EPS,hF,Z,Z0,Z1,Rh,Rd
      external hF
      DOUBLE PRECISION X(N)
      DOUBLE PRECISION Y,Y0,Y1
      M=0
      Z=A
      Y=hF(MM,Rh,Rd,Z)
10    IF ((Z.GT.B+H/2.0).OR.(M.EQ.N)) RETURN
      IF (ABS(Y).LT.EPS) THEN
        M=M+1
        X(M)=Z
        Z=Z+H/2.0d0
        Y=hF(mm,rh,rd,Z)
        GOTO 10
      END IF
      Z1=Z+H
      Y1=hF(MM,Rh,Rd,Z1)
      IF (ABS(Y1).LT.EPS) THEN
        M=M+1
        X(M)=Z1
        Z=Z1+H/2.0d0
        Y=hF(MM,Rh,Rd,Z)
        GOTO 10
      END IF
      IF (Y*Y1.GT.0.0d0) THEN
        Y=Y1
        Z=Z1
        GOTO 10
      END IF
20    IF (ABS(Z1-Z).LT.EPS) THEN
        M=M+1
        X(M)=(Z1+Z)/2.0d0
        Z=Z1+H/2.0d0
        Y=hF(MM,Rh,Rd,Z)
        GOTO 10
      END IF
      Z0=(Z1+Z)/2.0d0
      Y0=hF(MM,Rh,Rd,Z0)
      IF (ABS(Y0).LT.EPS) THEN
        M=M+1
        X(M)=Z0
        Z=Z0+H/2.0d0
        Y=hF(MM,Rh,Rd,Z)
        GOTO 10
      END IF
      IF (Y*Y0.LT.0.0d0) THEN
        Z1=Z0
        Y1=Y0
      ELSE
        Z=Z0
        Y=Y0
      END IF
      GOTO 20
      END SUBROUTINE hDDHRT
      
      
!-----------------------------------------------------------------------
!-----USER-DEFINED FUNCTION:hF
!-----in this program it refers to the zero determinant that calculates the characteristic function of hard-walled annular ducts
!-----------------------------------------------------------------------
      FUNCTION hF(MM,Rh,Rd,X)
      DOUBLE PRECISION X,hF,MBSL1,MBSL2,Rh,Rd
      external MBSL1,MBSL2
      IF(Rh==0.D0)THEN
        IF(MM==0)THEN
          hF=-MBSL1(MM+1,X*Rd)
          ELSE
          hF=0.5D0*(MBSL1(MM-1,X*Rd)-MBSL1(MM+1,X*Rd))
        END IF
      ELSE
        IF(MM==0)THEN
          hF=MBSL1(1,Rh*X)*MBSL2(1,Rd*X)-MBSL1(1,Rd*X)*MBSL2(1,Rh*X)
        ELSE
          hF=(MBSL1(MM-1,Rh*X)-MBSL1(MM+1,Rh*X))*(MBSL2(MM-1,Rd*X)-MBSL2(MM+1,Rd*X))&
     &       -(MBSL1(MM-1,Rd*X)-MBSL1(MM+1,Rd*X))*(MBSL2(MM-1,Rh*X)-MBSL2(MM+1,Rh*X))
        END IF
      END IF
      RETURN
      END FUNCTION hF

      
!-----------------------------------------------------------------------
!-----THE 1st BESSEL FUNCTION---of order NN: FUNCTION MBSL1(N,X)
!-----------------------------------------------------------------------     
      FUNCTION MBSL1(NN,X)
      DOUBLE PRECISION MBSL1,X
      DOUBLE PRECISION T,Y,Z,P,Q,S,B0,B1
      DOUBLE PRECISION A(6),B(6),C(6),D(6),E(5),F(5),G(5),H(5)
      DATA A/57568490574.0,-13362590354.0,651619640.7,-11214424.18,77392.33017,-184.9052456/
      DATA B/57568490411.0,1029532985.0,9494680.718,59272.64853,267.8532712,1.0/
      DATA C/72362614232.0,-7895059235.0,242396853.1,-2972611.439,15704.4826,-30.16036606/
      DATA D/144725228443.0,2300535178.0,18583304.74,99447.43394,376.9991397,1.0/
      DATA E/1.0,-0.1098628627D-02,0.2734510407D-04,-0.2073370639D-05,0.2093887211D-06/
      DATA F/-0.1562499995D-01,0.1430488765D-03,-0.6911147651D-05,0.7621095161D-06,-0.934935152D-07/
      DATA G/1.0,0.183105D-02,-0.3516396496D-04,0.2457520174D-05,-0.240337019D-06/
      DATA H/0.4687499995D-01,-0.2002690873D-03,0.8449199096D-05,-0.88228987D-06,0.105787412D-06/
      T=ABS(X)
      N=NN
      IF (N.LT.0) N=-N
      IF (N.NE.1) THEN
        IF (T.LT.8.0) THEN
          Y=T*T
          P=A(6)
          Q=B(6)
          DO 10 I=5,1,-1
            P=P*Y+A(I)
            Q=Q*Y+B(I)
10        CONTINUE
          P=P/Q
        ELSE
          Z=8.0/T
          Y=Z*Z
          P=E(5)
          Q=F(5)
          DO 20 I=4,1,-1
            P=P*Y+E(I)
            Q=Q*Y+F(I)
20        CONTINUE
          S=T-0.785398164
          P=P*COS(S)-Z*Q*SIN(S)
          P=P*SQRT(0.636619772/T)
        END IF
      END IF
      IF (N.EQ.0) THEN
        MBSL1=P
        RETURN
      END IF
      B0=P
      IF (T.LT.8.0) THEN
        Y=T*T
        P=C(6)
        Q=D(6)
        DO 30 I=5,1,-1
          P=P*Y+C(I)
          Q=Q*Y+D(I)
30      CONTINUE
        P=X*P/Q
      ELSE
        Z=8.0/T
        Y=Z*Z
        P=G(5)
        Q=H(5)
        DO 40 I=4,1,-1
          P=P*Y+G(I)
          Q=Q*Y+H(I)
40      CONTINUE
        S=T-2.356194491
        P=P*COS(S)-Z*Q*SIN(S)
        P=P*X*SQRT(0.636619772/T)/T
      END IF
      IF (N.EQ.1) THEN
        MBSL1=P
        RETURN
      END IF
      B1=P
      IF (X.EQ.0.0) THEN
        MBSL1=0.0
        RETURN
      END IF
      S=2.0/T
      IF (T.GT.1.0*N) THEN
        IF (X.LT.0.0) B1=-B1
        DO 50 I=1,N-1
          P=S*I*B1-B0
          B0=B1
          B1=P
50      CONTINUE
      ELSE
        M=SQRT(40.0*N)
        M=(N+M)/2
        M=M+M
        P=0.0
        Q=0.0
        B0=1.0
        B1=0.0
        DO 60 I=M-1,0,-1
          T=S*(I+1)*B0-B1
          B1=B0
          B0=T
          IF (ABS(B0).GT.1.0D+10) THEN
            B0=B0*1.0D-10
            B1=B1*1.0D-10
            P=P*1.0D-10
            Q=Q*1.0D-10
          END IF
          IF (MOD(I+2,2).EQ.0) Q=Q+B0
          IF (i+1.EQ.N) P=B1
60      CONTINUE
        Q=2.0*Q-B0
        P=P/Q
      END IF
      IF ((X.LT.0.0).AND.(MOD(N,2).EQ.1)) P=-P
      MBSL1=P
      RETURN
      END
      
      
!-----------------------------------------------------------------------
!-----THE 2nd BESSEL FUNCTION---of order NN: FUNCTION MBSL2(N,X)
!-----------------------------------------------------------------------
      FUNCTION MBSL2(NN,X)
      DOUBLE PRECISION MBSL2,X,MBSL1
      DOUBLE PRECISION Y,Z,P,Q,S,B0,B1
      DOUBLE PRECISION A(6),B(6),C(6),D(7),E(5),F(5),G(5),H(5)
      DATA A/-2.957821389D+09,7.062834065D+09,-5.123598036D+08,1.087988129D+07,-8.632792757D+04,2.284622733D+02/
      DATA B/4.0076544269D+10,7.452499648D+08,7.189466438D+06,4.74472647D+04,2.261030244D+02,1.0/
      DATA C/-4.900604943D+12,1.27527439D+12,-5.153438139D+10,7.349264551D+08,-4.237922726D+06,8.511937935D+03/
      DATA D/2.49958057D+13,4.244419664D+11,3.733650367D+09,2.245904002D+07,1.02042605D+05,3.549632885D+02,1.0/
      DATA E/1.0,-0.1098628627D-02,0.2734510407D-04,-0.2073370639D-05,0.2093887211D-06/
      DATA F/-0.1562499995D-01,0.1430488765D-03,-0.6911147651D-05,0.7621095161D-06,-0.934935152D-07/
      DATA G/1.0,0.183105D-02,-0.3516396496D-04,0.2457520174D-05,-0.240337019D-06/
      DATA H/0.4687499995D-01,-0.2002690873D-03,0.8449199096D-05,-0.88228987D-06,0.105787412D-06/
      N=NN
      IF (N.LT.0) N=-N
      IF (X.LT.0.0) X=-X
      IF (X+1.0.EQ.1.0) THEN
        MBSL2=-1.0D+35
        RETURN
      END IF
      IF (N.NE.1) THEN
        IF (X.LT.8.0) THEN
          Y=X*X
          P=A(6)
          Q=B(6)
          DO 10 I=5,1,-1
            P=P*Y+A(I)
            Q=Q*Y+B(I)
10        CONTINUE
          P=P/Q+0.636619772*MBSL1(0,X)*LOG(X)
        ELSE
          Z=8.0/X
          Y=Z*Z
          P=E(5)
          Q=F(5)
          DO 20 I=4,1,-1
            P=P*Y+E(I)
            Q=Q*Y+F(I)
20        CONTINUE
          S=X-0.785398164
          P=P*SIN(S)+Z*Q*COS(S)
          P=P*SQRT(0.636619772/X)
        END IF
      END IF
      IF (N.EQ.0) THEN
        MBSL2=P
        RETURN
      END IF
      B0=P
      IF (X.LT.8.0) THEN
        Y=X*X
        P=C(6)
        Q=D(7)
        DO 30 I=5,1,-1
          P=P*Y+C(I)
          Q=Q*Y+D(I+1)
30      CONTINUE
        Q=Q*Y+D(1)
        P=X*P/Q+0.636619772*(MBSL1(1,X)*LOG(X)-1.0/X)
      ELSE
        Z=8.0/X
        Y=Z*Z
        P=G(5)
        Q=H(5)
        DO 40 I=4,1,-1
          P=P*Y+G(I)
          Q=Q*Y+H(I)
40      CONTINUE
        S=X-2.356194491
        P=P*SIN(S)+Z*Q*COS(S)
        P=P*SQRT(0.636619772/X)
      END IF
      IF (N.EQ.1) THEN
        MBSL2=P
        RETURN
      END IF
      B1=P
      S=2.0/X
      DO 50 I=1,N-1
        P=S*I*B1-B0
        B0=B1
        B1=P
50    CONTINUE
      MBSL2=P
      RETURN
      END
      
      
!-----------------------------------------------------------------------
!-----solution of linear equations system
!-----------------------------------------------------------------------	
      SUBROUTINE ACGAS(AR,AI,N,BR,BI,L,JS)
      DIMENSION AR(N,N),AI(N,N),BR(N),BI(N),JS(N)
      DOUBLE PRECISION AR,AI,BR,BI,D,P,Q,S,W
      L=1
      DO 100 K=1,N-1
        D=0.0D0
        DO 10 I=K,N
        DO 10 J=K,N
          P=AR(I,J)*AR(I,J)+AI(I,J)*AI(I,J)
          IF (P.GT.D) THEN
            D=P
            JS(K)=J
            IS=I
          END IF
10      CONTINUE
        W=D
        IF (W+1.0.EQ.1.0) THEN
          WRITE(*,20)
          L=0
          RETURN
        END IF
20      FORMAT(1X,'  ERR**FAIL  ')
        DO 30 J=K,N
          P=AR(K,J)
          AR(K,J)=AR(IS,J)
          AR(IS,J)=P
          P=AI(K,J)
          AI(K,J)=AI(IS,J)
          AI(IS,J)=P
30      CONTINUE
        P=BR(K)
        BR(K)=BR(IS)
        BR(IS)=P
        P=BI(K)
        BI(K)=BI(IS)
        BI(IS)=P
        DO 50 I=1,N
          P=AR(I,K)
          AR(I,K)=AR(I,JS(K))
          AR(I,JS(K))=P
          P=AI(I,K)
          AI(I,K)=AI(I,JS(K))
          AI(I,JS(K))=P
50      CONTINUE
        DO 60 J=K+1,N
          P=AR(K,J)*AR(K,K)
          Q=-AI(K,J)*AI(K,K)
          S=(AR(K,K)-AI(K,K))*(AR(K,J)+AI(K,J))
          AR(K,J)=(P-Q)/D
          AI(K,J)=(S-P-Q)/D
60      CONTINUE
        P=BR(K)*AR(K,K)
        Q=-BI(K)*AI(K,K)
        S=(AR(K,K)-AI(K,K))*(BR(K)+BI(K))
        BR(K)=(P-Q)/D
        BI(K)=(S-P-Q)/D
        DO 90 I=K+1,N
          DO 80 J=K+1,N
            P=AR(I,K)*AR(K,J)
            Q=AI(I,K)*AI(K,J)
            S=(AR(I,K)+AI(I,K))*(AR(K,J)+AI(K,J))
            AR(I,J)=AR(I,J)-P+Q
            AI(I,J)=AI(I,J)-S+P+Q
80        CONTINUE
          P=AR(I,K)*BR(K)
          Q=AI(I,K)*BI(K)
          S=(AR(I,K)+AI(I,K))*(BR(K)+BI(K))
          BR(I)=BR(I)-P+Q
          BI(I)=BI(I)-S+P+Q
90      CONTINUE
100   CONTINUE
       D=AR(N,N)*AR(N,N)+AI(N,N)*AI(N,N)
       W=D
       IF (W+1.0.EQ.1.0) THEN
         L=0
         WRITE(*,20)
         RETURN
       END IF
       P=AR(N,N)*BR(N)
       Q=-AI(N,N)*BI(N)
       S=(AR(N,N)-AI(N,N))*(BR(N)+BI(N))
       BR(N)=(P-Q)/D
       BI(N)=(S-P-Q)/D
       DO 200 I=N-1,1,-1
         DO 150 J=I+1,N
           P=AR(I,J)*BR(J)
           Q=AI(I,J)*BI(J)
           S=(AR(I,J)+AI(I,J))*(BR(J)+BI(J))
           BR(I)=BR(I)-P+Q
           BI(I)=BI(I)-S+P+Q
150      CONTINUE
200    CONTINUE
      JS(N)=N
      DO 110 K=N,1,-1
        P=BR(K)
        BR(K)=BR(JS(K))
        BR(JS(K))=P
        P=BI(K)
        BI(K)=BI(JS(K))
        BI(JS(K))=P
110   CONTINUE
      RETURN
      END
  
    

    
