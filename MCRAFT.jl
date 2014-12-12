function MCRAFT(f,tfinal)

M0=5.0
I0 = 7.5e-3;
CTA0 = 1.0e-2;

TP = zeros(Int16, 70000000);
PTP1 = zeros(Int16, 40000000);
PTP2 = zeros(Int16, 40000000);
RTP = zeros(Int16, 3000000);
D = zeros(Int16, 1000000);
R = zeros(Int16, 5000);

Avogadro = 6.022141986e23;
Avogadrof = Avogadro * f;
t_f = tfinal;

        n_mon = ifloor(M0*Avogadrof);
        nmon_inic = n_mon;

        n_ini = ifloor(I0*Avogadrof);
        nini_inic = n_ini;

        n_cta = ifloor(CTA0*Avogadrof);
        ncta_inic = n_cta;
    
        print (n_cta+n_mon+n_ini\n)

        n_rad = 0.0;
        n_istar=0.0;
        n_rstar=0.0;
        n_raft1=0.0;
        n_raft2=0.0;
        n_pol = 0.0;

        n_nothing=0;

        t_r = 0.0;

        #CONTADORES DE REACCIONES
        #contadorR1 = 0;
        #contadorR2 = 0;
        #contadorR3 = 0;
        #contadorR4 = 0;
        #contadorR5 = 0;
        #contadorR6 = 0;
        #contadorR7 = 0;
        #contadorR8 = 0;
        #contadorR9 = 0;
        #contadorR10 = 0;
        #contadorR11 = 0;
        #contadorR12 = 0;
        #contadorR13 = 0;
        #contadorR14 = 0;

        #INICIALIZAMOS LAS CADENAS FINALES (MAXIMAS ALCANZADAS) EN 1#
        #nfinal_D = 1;
        #nfinal_R = 1;
        #nfinal_TP = 1;
        #nfinal_RTP = 1;
        #nfinal_PTP = 1;

        #Vnfinal = [nfinal_D, nfinal_R, nfinal_TP, nfinal_RTP, nfinal_PTP];

        last_D = 1;
        last_R = 1;
        last_TP = 1;
        last_RTP = 1;
        last_PTP = 1;

        #INICIALIZAMOS EL CONTEO DEL TIEMPO DE COMPUTO Y DE IMPRESION#
        t_print=0.1;

        #Las saqué afuera del LAZO DE MC, por T=cte en este caso particular
        #Ademas V=cte, kd y kf tienen un solo reactivo. Al resto las divido
        #por Avogrado*f

        eff=0.5;
        kd=0.036;
        kp=3600000/Avogadrof;
        ki=kp; #ASUMIMOS
        ktd=3.6e10/Avogadrof;
        ktc=3.6e10/Avogadrof;
        ka=3.6e9/Avogadrof;
        kf=36.0;
    
    #longitudes maximas de vectores D, R, TP, RTP
        #maxlong=[0,0,0,0,0,0]
    
        R1 = kd*n_ini;
        R2 = ki*n_mon*n_istar;
        R3 = kp*n_mon*n_rad;

        R4 = ka*n_cta*n_rad;
        R5 = ka*n_istar*n_raft1;
        R6 = 0.5*kf*n_rstar;
        R7 = 0.5*kf*n_rstar;

        R8 = 0.0;

        R9 = ka*n_rad*n_raft1;
        R10 = 0.0;
        R11 = 0.5*kf*n_raft2;
        R12 = 0.5*kf*n_raft2;

        R13 = 0.5*ktc*n_rad*(n_rad-1);
        R14 = 0.5*ktd*n_rad*(n_rad-1);


#INICIO LAZO MONTECARLO#
        while t_r < t_f

        #VELOCIDADES DE REACCION Y PROBABILIDADES# (en esta version actualizo las R directamente dentro del If)

            
            
            

            
            
            

            

            

            R_t = R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13+R14;

            r2 = rand();

            if r2 > (R1+R2+R3+R4)/(R_t)
                reacselec = 5;
                if r2 > (R1+R2+R3+R4+R5)/(R_t)
                    reacselec = 6;
                    if r2 > (R1+R2+R3+R4+R5+R6)/(R_t)
                        reacselec = 7;
                        if r2 > (R1+R2+R3+R4+R5+R6+R7)/(R_t)
                            reacselec = 8;
                            if r2 > (R1+R2+R3+R4+R5+R6+R7+R8)/(R_t)
                                reacselec = 9;
                                if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9)/(R_t)
                                    reacselec = 10;
                                    if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10)/(R_t)
                                        reacselec = 11;
                                        if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11)/(R_t)
                                            reacselec = 12;
                                            if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12)/(R_t)
                                                reacselec = 13;
                                                if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13)/(R_t)
                                                    reacselec = 14;
                                                end
                                            end
                                        end
                                    end
                                end
                            end

                        end
                    end
                end
            elseif r2 > (R1+R2)/(R_t)
                reacselec = 3;
                if r2 > (R1+R2+R3)/(R_t)
                    reacselec = 4;
                end
            elseif r2 < R1/(R_t) #CUIDADO ACA DI VUELTA EL SIGNO y la reaccion#
                reacselec = 1;
            else
                reacselec = 2;
            end


            #CALCULO DEL INTERVALO Dt#

            t_r = t_r + (-log(rand()))/(R_t);

            if reacselec == 1
                    #contadorR1 = contadorR1 + 1;
                    n_ini = n_ini - 1;
                    if rand() <= eff
                        n_istar = n_istar + 2;
                        R2 = ki*n_mon*n_istar;
                        R5 = ka*n_istar*n_raft1;
                    else
                        n_nothing = n_nothing + 2;
                    end
                    
                    R1 = kd*n_ini;                   
                    

            elseif reacselec == 2
                    #contadorR2 = contadorR2 + 1;

                    n_istar = n_istar - 1;
                    n_mon = n_mon - 1;
                    n_rad = n_rad + 1;

                    @inbounds R[last_R] = 1;
                    last_R = last_R + 1;
                    
                    R2 = ki*n_mon*n_istar;
                    R3 = kp*n_mon*n_rad;
                    R4 = ka*n_cta*n_rad;
                    R5 = ka*n_istar*n_raft1;
                    R9 = ka*n_rad*n_raft1;  
                    R13 = 0.5*ktc*n_rad*(n_rad-1);
                    R14 = 0.5*ktd*n_rad*(n_rad-1);
                    

            elseif reacselec == 3
                    #contadorR3 = contadorR3 + 1;

                    Rselec = ceil(n_rad*rand());

                    n_mon = n_mon - 1;

                    @inbounds R[Rselec] += 1;
            
                    R2 = ki*n_mon*n_istar;
                    R3 = kp*n_mon*n_rad;
                    

            elseif reacselec == 4
                    #contadorR4 = contadorR4 + 1;

                    Rselec = ceil(n_rad*rand());

                    n_rad = n_rad - 1;
                    n_cta = n_cta - 1;
                    n_rstar = n_rstar + 1;

                    @inbounds RTP[last_RTP] = R[Rselec];
                    last_RTP = last_RTP + 1;

                    @inbounds R[Rselec] = R[last_R-1];
                    @inbounds R[last_R-1] = 0;
                    last_R = last_R - 1;
            
                    R3 = kp*n_mon*n_rad;
                    R4 = ka*n_cta*n_rad;
                    R6 = 0.5*kf*n_rstar;
                    R7 = 0.5*kf*n_rstar;
                    R9 = ka*n_rad*n_raft1;
            R13 = 0.5*ktc*n_rad*(n_rad-1);
            R14 = 0.5*ktd*n_rad*(n_rad-1);

          elseif reacselec == 5

                  #contadorR5 = contadorR5 + 1;

                  TPselec = ceil(n_raft1*rand());

                  n_istar = n_istar - 1;
                  n_raft1 = n_raft1 - 1;
                  n_rstar = n_rstar + 1;

                  @inbounds RTP[last_RTP] = TP[TPselec];
                  last_RTP = last_RTP + 1;

                  @inbounds TP[TPselec] = TP[last_TP-1];
                  @inbounds TP[last_TP-1] = 0;
                  last_TP = last_TP - 1;
            
                    R2 = ki*n_mon*n_istar;
                    R5 = ka*n_istar*n_raft1;
                    R6 = 0.5*kf*n_rstar;
                    R7 = 0.5*kf*n_rstar;
                    R9 = ka*n_rad*n_raft1;
                    

          elseif reacselec == 6
                  #contadorR6 = contadorR6 + 1;

                  RTP_selec = ceil(n_rstar*rand());

                  n_rstar = n_rstar - 1;
                  n_rad = n_rad + 1;
                  n_cta = n_cta + 1;

                  @inbounds R[last_R] = RTP[RTP_selec];
                  last_R = last_R + 1;

                  @inbounds RTP[RTP_selec] = RTP[last_RTP-1];
                  @inbounds RTP[last_RTP-1] = 0;
                  last_RTP = last_RTP - 1;
            
                    R3 = kp*n_mon*n_rad;
                    R4 = ka*n_cta*n_rad;
                    R6 = 0.5*kf*n_rstar;
                    R7 = 0.5*kf*n_rstar;
                    R9 = ka*n_rad*n_raft1;
            R13 = 0.5*ktc*n_rad*(n_rad-1);
            R14 = 0.5*ktd*n_rad*(n_rad-1);

         elseif reacselec == 7 #preeq 4 - fragmentación 2
                #contadorR7 = contadorR7 + 1;

                RTP_selec = ceil(n_rstar*rand());

                n_rstar = n_rstar - 1;
                n_istar = n_istar + 1;
                n_raft1 = n_raft1 + 1;
            
                @inbounds TP[last_TP] = RTP[RTP_selec];
                last_TP = last_TP + 1;

                @inbounds RTP[RTP_selec] = RTP[last_RTP-1];
                @inbounds RTP[last_RTP-1] = 0;
                last_RTP = last_RTP - 1;
            
                R2 = ki*n_mon*n_istar;
                R5 = ka*n_istar*n_raft1;
                R6 = 0.5*kf*n_rstar;
                R7 = 0.5*kf*n_rstar;
                R9 = ka*n_rad*n_raft1;

        elseif reacselec == 8 #reiniciacion y propagacion del rstar liberado por cta
                #contadorR8 = contadorR8 + 1;

                n_istar = n_istar - 1;
                n_mon = n_mon - 1;
                n_rad = n_rad + 1;
            
                R2 = ki*n_mon*n_istar;
                R3 = kp*n_mon*n_rad;
                R4 = ka*n_cta*n_rad;
                R5 = ka*n_istar*n_raft1;
                R9 = ka*n_rad*n_raft1;
            R13 = 0.5*ktc*n_rad*(n_rad-1);
            R14 = 0.5*ktd*n_rad*(n_rad-1);

        elseif reacselec == 9 #core eq - adicion 1
            #contadorR9 = contadorR9 + 1;

                Rselec = ceil(n_rad*rand());
                TPselec = ceil(n_raft1*rand());

                n_rad = n_rad - 1;
                n_raft1 = n_raft1 - 1;
                n_raft2 = n_raft2 + 1;

                @inbounds PTP1[last_PTP] = R[Rselec];
                @inbounds PTP2[last_PTP] = TP[TPselec];
                last_PTP = last_PTP + 1;

                @inbounds R[Rselec] = R[last_R-1];
                @inbounds R[last_R-1] = 0;
                last_R = last_R - 1;

                @inbounds TP[TPselec] = TP[last_TP-1];
                @inbounds TP[last_TP-1] = 0;
                last_TP = last_TP - 1;
            
                R3 = kp*n_mon*n_rad;
                R4 = ka*n_cta*n_rad;
                R5 = ka*n_istar*n_raft1;
                R9 = ka*n_rad*n_raft1;
                R11 = 0.5*kf*n_raft2;
                R12 = 0.5*kf*n_raft2;
            R13 = 0.5*ktc*n_rad*(n_rad-1);
            R14 = 0.5*ktd*n_rad*(n_rad-1);

        elseif reacselec == 10
                  #contadorR10 = contadorR10 + 1;

                  n_rad = n_rad - 1;
                  n_raft1 = n_raft1 - 1;
                  n_raft2 = n_raft2 + 1;
            
                  R3 = kp*n_mon*n_rad;
                  R4 = ka*n_cta*n_rad;
                  R5 = ka*n_istar*n_raft1;
                  R9 = ka*n_rad*n_raft1;
                R11 = 0.5*kf*n_raft2;
                R12 = 0.5*kf*n_raft2;
            R13 = 0.5*ktc*n_rad*(n_rad-1);
            R14 = 0.5*ktd*n_rad*(n_rad-1);
    
        elseif reacselec == 11
                 #contadorR11 = contadorR11 + 1;

                 PTPselec = ceil(n_raft2*rand());

                 n_raft2 = n_raft2 - 1;
                 n_rad = n_rad + 1;
                 n_raft1 = n_raft1 + 1;

                 @inbounds R[last_R] = PTP1[PTPselec];
                 @inbounds TP[last_TP] = PTP2[PTPselec];

                 @inbounds PTP1[PTPselec]=PTP1[last_PTP-1];
                 @inbounds PTP2[PTPselec]=PTP2[last_PTP-1];
                 @inbounds PTP1[last_PTP-1] = 0;
                 @inbounds PTP2[last_PTP-1] = 0;
                 last_PTP = last_PTP - 1;

                 last_R = last_R + 1;
                 last_TP = last_TP + 1;
            
                    R3 = kp*n_mon*n_rad;
                    R4 = ka*n_cta*n_rad;
                    R5 = ka*n_istar*n_raft1;
                    R9 = ka*n_rad*n_raft1;
                R11 = 0.5*kf*n_raft2;
                R12 = 0.5*kf*n_raft2;
            R13 = 0.5*ktc*n_rad*(n_rad-1);
            R14 = 0.5*ktd*n_rad*(n_rad-1);

        elseif reacselec == 12
                 #contadorR12 = contadorR12 + 1;

                 PTPselec = ceil(n_raft2*rand());

                 n_raft2 = n_raft2 - 1;
                 n_rad = n_rad + 1;
                 n_raft1 = n_raft1 + 1;

                 @inbounds TP[last_TP] = PTP1[PTPselec];
                 @inbounds R[last_R] = PTP2[PTPselec];

                 @inbounds PTP1[PTPselec]=PTP1[last_PTP-1];
                 @inbounds PTP2[PTPselec]=PTP2[last_PTP-1];
                 @inbounds PTP1[last_PTP-1] = 0;
                 @inbounds PTP2[last_PTP-1] = 0;
                 last_PTP = last_PTP - 1;

                 last_R = last_R + 1;
                 last_TP = last_TP + 1;
            
                R3 = kp*n_mon*n_rad;
                R4 = ka*n_cta*n_rad;
                R5 = ka*n_istar*n_raft1;
                R9 = ka*n_rad*n_raft1;
            R11 = 0.5*kf*n_raft2;
            R12 = 0.5*kf*n_raft2;
            R13 = 0.5*ktc*n_rad*(n_rad-1);
            R14 = 0.5*ktd*n_rad*(n_rad-1);

        elseif reacselec == 13
                  #contadorR13 = contadorR13 + 1;

                  R1selec = ceil(rand()*(n_rad-1));
                  R2selec = ceil(rand()*(n_rad-1));

                  while R1selec == R2selec
                      if n_rad != 2
                          R2selec = ceil(rand()*(n_rad-1));
                      else
                          R2selec = 2;
                      end
                  end

                  n_rad = n_rad - 2;
                  n_pol = n_pol + 1;

                  @inbounds D[last_D] = R[R1selec] + R[R2selec];
                  last_D = last_D + 1;

                  @inbounds R[R1selec] = R[last_R-1];
                  @inbounds R[last_R-1] = 0;
                  @inbounds R[R2selec] = R[last_R-2];
                  @inbounds R[last_R-2] = 0;
                  last_R = last_R - 2;
            
                R3 = kp*n_mon*n_rad;
                R4 = ka*n_cta*n_rad;
                R9 = ka*n_rad*n_raft1;
            R13 = 0.5*ktc*n_rad*(n_rad-1);
            R14 = 0.5*ktd*n_rad*(n_rad-1);

        elseif reacselec == 14
                  #contadorR14 = contadorR14 + 1;

                  R1selec = ceil(rand()*(n_rad-1));
                  R2selec = ceil(rand()*(n_rad-1));

                  while R1selec == R2selec
                      if n_rad != 2
                          R2selec = ceil(rand()*(n_rad-1));
                      else
                          R2selec = 2;
                      end
                  end

                  n_rad = n_rad - 2;
                  n_pol = n_pol + 2;

                  @inbounds D[last_D] = R[R1selec];
                  last_D = last_D + 1;

                  @inbounds D[last_D] = R[R2selec];
                  last_D = last_D + 1;

                  @inbounds R[R1selec] = R[last_R-1];
                  @inbounds R[last_R-1] = 0;
                  @inbounds R[R2selec] = R[last_R-2];
                  @inbounds R[last_R-2] = 0;
                  last_R = last_R - 2;
            
                R3 = kp*n_mon*n_rad;
                R4 = ka*n_cta*n_rad;
                R9 = ka*n_rad*n_raft1;
            R13 = 0.5*ktc*n_rad*(n_rad-1);
            R14 = 0.5*ktd*n_rad*(n_rad-1);

          end
        
        
#         if n_pol > maxlong[1] #D
#             maxlong[1]=n_pol
#         end
#         if n_rad>maxlong[2] #R
#             maxlong[2]=n_rad
#         end
#         if n_istar>maxlong[3] #
#             maxlong[3]=n_istar
#         end
#         if n_rstar>maxlong[4] #RTP
#             maxlong[4]=n_rstar
#         end
#         if n_raft1>maxlong[5] #TP
#             maxlong[5]=n_raft1
#         end
#         if n_raft2>maxlong[6] #PTP
#             maxlong[6]=n_raft2
#         end


            end

    return D,PTP1,PTP2
           
end
