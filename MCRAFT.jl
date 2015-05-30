function initialize_arrays()
    TP = zeros(Int16, 8000000);
    PTP1 = zeros(Int16, 5000000);
    PTP2 = zeros(Int16, 5000000);
    RTP = zeros(Int16, 3000000);
    D = zeros(Int16, 1000);
    R = zeros(Int16, 200)
    return TP, PTP1, PTP2, RTP, D, R
end


function MCRAFTd(f,tfinal,TP,PTP1,PTP2,RTP,D,R)

    M0=5.0
    I0 = 7.5e-3;
    CTA0 = 1.0e-2;
    Avogadro = 6.022141986e23;
    Avogadrof = Avogadro * f;
    t_f = tfinal;
        
    n_mon = floor(Integer, M0*Avogadrof);
    nmon_inic = n_mon;
    n_ini = floor(Integer, I0*Avogadrof);
    nini_inic = n_ini;
    n_cta = floor(Integer, CTA0*Avogadrof);
    ncta_inic = n_cta;
    
    n_rad = 0.0;
    n_istar=0.0;
    n_rstar=0.0;
    n_raft1=0.0;
    n_raft2=0.0;
    n_pol = 0.0;
    n_nothing=0;

    t_r = 0.0;

    """
    #reaction counters... these are not needed anymore. We used them to identify the most common reaction
    contadorR1 = 0;
    contadorR2 = 0;
    contadorR3 = 0;
    contadorR4 = 0;
    contadorR5 = 0;
    contadorR6 = 0;
    contadorR7 = 0;
    contadorR8 = 0;
    contadorR9 = 0;
    contadorR10 = 0;
    contadorR11 = 0;
    contadorR12 = 0;
    contadorR13 = 0;
    contadorR14 = 0;
    """
    #These might not be needed anymore, I have to check if findlast() works and if its faster.
    last_D = 1;
    last_R = 1;
    last_TP = 1;
    last_RTP = 1;
    last_PTP = 1;

    #This was used to print results every 0.1 hrs of reaction time, not really needed anymore.
    t_print=0.1;

    """
    These kinetic constants can be outside the loop, *only* because for this particular reaction
    the temperature and pressure remain constant. Otherwise they go inside the loop.
    """

    eff=0.5;
    kd=0.036;
    kp=3600000/Avogadrof;
    ki=kp; #assumed
    ktd=3.6e10/Avogadrof;
    ktc=3.6e10/Avogadrof;
    ka=3.6e9/Avogadrof;
    kf=36.0;
    
    #Both, maxlong and iter were used for debugging only.
    #maxlong=[0,0,0,0,0,0]
    iter=0
    
    """
    These are also outside of the loop now, because we need to initialize them.
    We are now updating the value inside each reaction simulation (there was no need to update ALL of them
    in every iteration)
    There might be a better way to do this...
    """
    
    R1 = kd*n_ini;
    R2 = ki*n_mon*n_istar;
    R3 = kp*n_mon*n_rad;
    R4 = ka*n_cta*n_rad;
    R5 = ka*n_istar*n_raft1;
    R6 = 0.5*kf*n_rstar; # 0.5 goes because n_rstar fragments with .5 chance to R6 and .5 to R7
    R7 = 0.5*kf*n_rstar; # same as R6
    R8 = 0.0; #since we decided reactions 8 and 10 won't be simulated, they can be erased I think
    R9 = ka*n_rad*n_raft1;
    R10 = 0.0; # same as R8
    R11 = 0.5*kf*n_raft2; # same as R6
    R12 = 0.5*kf*n_raft2; # same as R6
    R13 = ktc*n_rad*(n_rad-1); #
    R14 = ktd*n_rad*(n_rad-1); #


    #MC Loop begins here
    while t_r < t_f
        """
        This is the reaction selection block.
        I'm sure there's a better way to do this... The way it is, it works but it's not clean.
        Also, for more complex reaction sets we might need to improve it anyway.
        """
        R_t = R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13+R14;

        r2 = rand();

        if r2 > (R1+R2+R3+R4)/(R_t)
            reacselec = 5;
            if r2 > (R1+R2+R3+R4+R5)/(R_t)
                reacselec = 6;
                if r2 > (R1+R2+R3+R4+R5+R6)/(R_t)
                    reacselec = 7;
                    #if r2 > (R1+R2+R3+R4+R5+R6+R7)/(R_t)
                        #reacselec = 8;
                        if r2 > (R1+R2+R3+R4+R5+R6+R7+R8)/(R_t)
                            reacselec = 9;
                            #if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9)/(R_t)
                                #reacselec = 10;
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
                            #end
                        end

                    #end
                end
            end
        elseif r2 > (R1+R2)/(R_t)
            reacselec = 3;
            if r2 > (R1+R2+R3)/(R_t)
                reacselec = 4;
            end
            elseif r2 < R1/(R_t) #Careful, I inverted the inequality here and the reactions too
            reacselec = 1;
        else
            reacselec = 2;
        end


        #time update.
        t_r = t_r + (-log(rand()))/(R_t);

        """
        Simulation of each reaction:
        Once the reaction was selected, we do 4 things in the worst case scenario:
        1. update the number of molecules
        2. randomly choose which molecule will react (this is now done with O(k) at expenses of big arrays)
        3. update the arrays accordingly
        4. update the reaction rates
        """
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
            
            Rselec = ceil(Integer, n_rad*rand());

            n_mon = n_mon - 1;

            @inbounds R[Rselec] += 1;

            R2 = ki*n_mon*n_istar;
            R3 = kp*n_mon*n_rad;


        elseif reacselec == 4
            #contadorR4 = contadorR4 + 1;

            Rselec = ceil(Integer, n_rad*rand());

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

            TPselec = ceil(Integer, n_raft1*rand());

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

            RTP_selec = ceil(Integer, n_rstar*rand());

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

            RTP_selec = ceil(Integer, n_rstar*rand());

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
            
            
        elseif reacselec == 9 #core eq - adicion 1
            #contadorR9 = contadorR9 + 1;

            Rselec = ceil(Integer, n_rad*rand());
            TPselec = ceil(Integer, n_raft1*rand());

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

        elseif reacselec == 11
            #contadorR11 = contadorR11 + 1;

            PTPselec = ceil(Integer, n_raft2*rand());

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

            PTPselec = ceil(Integer, n_raft2*rand());

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

            R1selec = ceil(Integer, rand()*(n_rad-1));
            R2selec = ceil(Integer, rand()*(n_rad-1));

            while R1selec == R2selec
              if n_rad != 2
                  R2selec = ceil(Integer, rand()*(n_rad-1));
              else
                  R2selec = 2;
              end
            end

            n_rad = n_rad - 2;
            n_pol = n_pol + 1;

            #@inbounds D[last_D] = R[R1selec] + R[R2selec];
            #last_D = last_D + 1;
            @inbounds D[R[R1selec]+R[R2selec]] += 1

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

            R1selec = ceil(Integer, rand()*(n_rad-1));
            R2selec = ceil(Integer, rand()*(n_rad-1));

            while R1selec == R2selec
              if n_rad != 2
                  R2selec = ceil(Integer, rand()*(n_rad-1));
              else
                  R2selec = 2;
              end
            end

            n_rad = n_rad - 2;
            n_pol = n_pol + 2;

            @inbounds D[R[R1selec]] += 1
            @inbounds D[R[R2selec]] += 1
            
            """
            We changed the representation for D.
            """
            #@inbounds D[last_D] = R[R1selec];
            #@inbounds D[last_D] = R[R2selec];
            #last_D = last_D + 2;
            

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

        end #if
        

"""
This is all debugging code
        
    if t_r > t_print;

        i+=1

        if t_r-t_print < 0.1
            t_print = t_print + 0.1;
        else
            t_print = t_r + 0.1;
        endint

        @inbounds X[i] = 100.0*(nmon_inic - n_mon)/(nmon_inic);
    end

     if n_pol > maxlong[1] #D
         maxlong[1]=n_pol
     end
     if n_rad>maxlong[2] #R
         maxlong[2]=n_rad
     end
     if n_istar>maxlong[3] #
         maxlong[3]=n_istar
     end
     if n_rstar>maxlong[4] #RTP
         maxlong[4]=n_rstar
     end
     if n_raft1>maxlong[5] #TP
         maxlong[5]=n_raft1
     end
     if n_raft2>maxlong[6] #PTP
         maxlong[6]=n_raft2
     end
"""        
        #iter += 1
        end #while

#counters = [contadorR1,contadorR2,contadorR3,contadorR4,contadorR5,contadorR6,contadorR7,contadorR8,contadorR9,contadorR10,contadorR11,contadorR12,contadorR13,contadorR14]       
conv = 100.0*(nmon_inic - n_mon)/(nmon_inic)
    
    #return everything for comparison purposes only. All postprocessing is done before plotting now.
    return D,PTP1,PTP2,TP,RTP,R,conv,iter

end #function
