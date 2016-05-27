function MCRAFT_SF(N,tfinal,M0,I0,CTA0,TP,PTP1,PTP2,RTP,D,R)

    """

include("C:/Users/epintos/Dropbox/Esteban/Tesis - privado/Julia/MCRAFT_SF_05.jl")

    Variable initialization starts here.
    When porting to Julia 0.4, pay extra attention to the integer conversion and integer functions.
    For example, instead of ifloor(x) --> change to floor(Integer, x)

    Defaults: (These are passed as arguments now)
    M0 = 5.0
    I0 = 7.5e-3
    CTA0 = 1e-2

    M0 = 5.0
    I0 = 5e-3
    CTA0 = 5e-3
    """

    """
    Arrays are initialzed with zeros and for the moment, using the highest value for the largest test
    set. Eventually, the size will depend on the value of the parameter f.
    Arrays are initialized outside the function now.
    """

    Avogadro = 6.022141986e23
    f=N/(M0*Avogadro+I0*Avogadro+CTA0*Avogadro)
    Avogadrof = Avogadro * f

    t_f = tfinal #this might be unnecesary

    n_mon = floor(Integer,M0*Avogadrof)
    n_ini = floor(Integer,I0*Avogadrof)
    n_cta = floor(Integer,CTA0*Avogadrof)
    nmon_inic = n_mon
    nini_inic = n_ini
    ncta_inic = n_cta

    n_total = n_mon+n_ini+n_cta

    n_rad = 0
    n_istar=0
    n_rstar=0
    n_raft1=0
    n_raft2=0
    n_pol = 0
    n_nothing=0

    t_r = 0.0

    #Reaction counters, these are for debugging purposes only.
    #Reacciones = zeros(Int,14)


    """
    These are used to track the last element added to the array, right now this is faster than
    doing find(x)[end]. Will have to test findlast(x) once implemented in 0.4.
    """
    last_D = 1
    last_R = 1
    last_TP = 1
    last_RTP = 1
    last_PTP = 1

    #This is to print every t_print hrs of reaction time - mainly for debugging purposes.
    t_print=0.001

    """
    Initially the MC rate constants were inside the MC loop, since for this particular reaction set
    the temperature and volume remains constant, the rate constants won't change so they can be
    outside of the loop.
    """
    eff=0.5
    kd=0.036
    kp=3600000.0/Avogadrof
    ki=3600000.0/Avogadrof #assumption
    ktd=3.6e10/Avogadrof
    ktc=3.6e10/Avogadrof
    ka=3.6e9/Avogadrof
    kf=36.0

    #This is just to keep track of the size of the arrays actually used - just for debugging
    #maxlong=[0,0,0,0,0,0]

    #Iteration counter - (for debugging)
    #iter=0

    """
    Rate constants are calculated outside the loop initially and update only when needed, inside
    the loop. Also, only the constants that were modified are calculated.
    """
    Re = zeros(Float64,14)
    Re[1] = kd*n_ini
    Re[2] = ki*n_mon*n_istar
    Re[3] = kp*n_mon*n_rad
    Re[4] = ka*n_cta*n_rad
    Re[5] = ka*n_istar*n_raft1
    Re[6] = 0.5*kf*n_rstar #the 0.5 is there because RAFT* reacts with 0.5 chance to Re[6] and 0.5 to Re[7]
    Re[7] = 0.5*kf*n_rstar #same than Re[6]
    Re[8] = 0.0; #Since we are sure this reaction doesnt go anymore, it could be removed
    Re[9] = ka*n_rad*n_raft1
    Re[10] = 0.0; #same than Re[8]
    Re[11] = 0.5*kf*n_raft2 #same than Re[6]
    Re[12] = 0.5*kf*n_raft2 #same than Re[6]
    Re[13] = ktc*n_rad*(n_rad-1)
    Re[14] = ktd*n_rad*(n_rad-1)

    """Debugging variables"""
    """
    species = ["tiempo","M","I","CTA","R","TP*","R*","TP","PTP","D","nothing"]
    Concentrations = zeros(ceil(Integer,t_f/t_print)+1,length(species))
    Concentrations[1,:] = [t_r,M0,I0,CTA0,n_rad,n_istar,n_rstar,n_raft1,n_raft2,n_pol,n_nothing]
    print_counter = 1
    """

    function linearwalk(R)
      partialcumsum = 0.0
      r2 = rand()*sum(R)
      i = 1
      while true
          partialcumsum += R[i]
          if r2 <= partialcumsum
              return i
          end
          i += 1
      end
    end

    #MC simluation starts here
    while t_r < t_f

        reacselec = linearwalk(Re)

        t_r = t_r + randexp()/R_t #this is even faster


        if reacselec == 1
            #Reacciones[1] += 1

            n_ini = n_ini - 1;

            if rand() <= eff
                n_istar = n_istar + 2;
                Re[2] = ki*n_mon*n_istar;
                Re[5] = ka*n_istar*n_raft1;
            else
                n_nothing = n_nothing + 2;
            end

            Re[1] = kd*n_ini;


        elseif reacselec == 2
            #Reacciones[2] += 1

            n_istar = n_istar - 1;
            n_mon = n_mon - 1;
            n_rad = n_rad + 1;

            @inbounds R[last_R] = 1;
            last_R = last_R + 1;

            Re[2] = ki*n_mon*n_istar;
            Re[3] = kp*n_mon*n_rad;
            Re[4] = ka*n_cta*n_rad;
            Re[5] = ka*n_istar*n_raft1;
            Re[9] = ka*n_rad*n_raft1;
            Re[13] = ktc*n_rad*(n_rad-1);
            Re[14] = ktd*n_rad*(n_rad-1);


        elseif reacselec == 3
            #Reacciones[3] += 1

            Rselec = unsafe_trunc(Int,n_rad*rand())+1
            """
            if Rselec == 0
                println("Rselec = 0 en reacselec == 3")
                Rselec = n_rad
            end
            """
            n_mon = n_mon - 1;

            @inbounds R[Rselec] += 1;

            Re[2] = ki*n_mon*n_istar;
            Re[3] = kp*n_mon*n_rad;


        elseif reacselec == 4
            #Reacciones[4] += 1

            Rselec = unsafe_trunc(Int,n_rad*rand())+1
            """
            if Rselec == 0
                Rselec = n_rad
            end
            """
            n_rad = n_rad - 1;
            n_cta = n_cta - 1;
            n_rstar = n_rstar + 1;

            @inbounds RTP[last_RTP] = R[Rselec];
            last_RTP = last_RTP + 1;

            @inbounds R[Rselec] = R[last_R-1];
            @inbounds R[last_R-1] = 0;
            last_R = last_R - 1;

            Re[3] = kp*n_mon*n_rad;
            Re[4] = ka*n_cta*n_rad;
            Re[6] = 0.5*kf*n_rstar;
            Re[7] = 0.5*kf*n_rstar;
            Re[9] = ka*n_rad*n_raft1;
            Re[13] = ktc*n_rad*(n_rad-1);
            Re[14] = ktd*n_rad*(n_rad-1);

        elseif reacselec == 5
            #Reacciones[5] += 1

            TPselec = unsafe_trunc(Int,n_raft1*rand())+1

            n_istar = n_istar - 1;
            n_raft1 = n_raft1 - 1;
            n_rstar = n_rstar + 1;

            @inbounds RTP[last_RTP] = TP[TPselec];
            last_RTP = last_RTP + 1;

            @inbounds TP[TPselec] = TP[last_TP-1];
            @inbounds TP[last_TP-1] = 0;
            last_TP = last_TP - 1;

            Re[2] = ki*n_mon*n_istar;
            Re[5] = ka*n_istar*n_raft1;
            Re[6] = 0.5*kf*n_rstar;
            Re[7] = 0.5*kf*n_rstar;
            Re[9] = ka*n_rad*n_raft1;

        elseif reacselec == 6
           #Reacciones[6] += 1

            RTP_selec = unsafe_trunc(Int,n_rstar*rand())+1

            n_rstar = n_rstar - 1;
            n_rad = n_rad + 1;
            n_cta = n_cta + 1;

            @inbounds R[last_R] = RTP[RTP_selec];
            last_R = last_R + 1;

            @inbounds RTP[RTP_selec] = RTP[last_RTP-1];
            @inbounds RTP[last_RTP-1] = 0;
            last_RTP = last_RTP - 1;

            Re[3] = kp*n_mon*n_rad;
            Re[4] = ka*n_cta*n_rad;
            Re[6] = 0.5*kf*n_rstar;
            Re[7] = 0.5*kf*n_rstar;
            Re[9] = ka*n_rad*n_raft1;
            Re[13] = ktc*n_rad*(n_rad-1);
            Re[14] = ktd*n_rad*(n_rad-1);

        elseif reacselec == 7
            #Reacciones[7] += 1

            RTP_selec = unsafe_trunc(Int,n_rstar*rand())+1

            n_rstar = n_rstar - 1;
            n_istar = n_istar + 1;
            n_raft1 = n_raft1 + 1;

            @inbounds TP[last_TP] = RTP[RTP_selec];
            last_TP = last_TP + 1;

            @inbounds RTP[RTP_selec] = RTP[last_RTP-1];
            @inbounds RTP[last_RTP-1] = 0;
            last_RTP = last_RTP - 1;

            Re[2] = ki*n_mon*n_istar;
            Re[5] = ka*n_istar*n_raft1;
            Re[6] = 0.5*kf*n_rstar;
            Re[7] = 0.5*kf*n_rstar;
            Re[9] = ka*n_rad*n_raft1;

        elseif reacselec == 8
            #Reacciones[8] += 1


        elseif reacselec == 9
            #Reacciones[9] += 1

            Rselec = unsafe_trunc(Int,rand()*n_rad)+1
            """
            if Rselec == 0
                Rselec = n_rad
            end
            """
            TPselec = unsafe_trunc(Int,n_raft1*rand())+1

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

            Re[3] = kp*n_mon*n_rad;
            Re[4] = ka*n_cta*n_rad;
            Re[5] = ka*n_istar*n_raft1;
            Re[9] = ka*n_rad*n_raft1;
            Re[11] = 0.5*kf*n_raft2;
            Re[12] = 0.5*kf*n_raft2;
            Re[13] = ktc*n_rad*(n_rad-1);
            Re[14] = ktd*n_rad*(n_rad-1);

        elseif reacselec == 10
            #Reacciones[10] += 1

        elseif reacselec == 11
            #Reacciones[11] += 1

            PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

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

            Re[3] = kp*n_mon*n_rad;
            Re[4] = ka*n_cta*n_rad;
            Re[5] = ka*n_istar*n_raft1;
            Re[9] = ka*n_rad*n_raft1;
            Re[11] = 0.5*kf*n_raft2;
            Re[12] = 0.5*kf*n_raft2;
            Re[13] = ktc*n_rad*(n_rad-1);
            Re[14] = ktd*n_rad*(n_rad-1);

        elseif reacselec == 12
            #Reacciones[12] += 1

            PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

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

            Re[3] = kp*n_mon*n_rad;
            Re[4] = ka*n_cta*n_rad;
            Re[5] = ka*n_istar*n_raft1;
            Re[9] = ka*n_rad*n_raft1;
            Re[11] = 0.5*kf*n_raft2;
            Re[12] = 0.5*kf*n_raft2;
            Re[13] = ktc*n_rad*(n_rad-1);
            Re[14] = ktd*n_rad*(n_rad-1);

        elseif reacselec == 13
            #Reacciones[13] += 1

            R1selec = unsafe_trunc(Int,rand()*(n_rad-1))+1
            R2selec = unsafe_trunc(Int,rand()*(n_rad-1))+1

            while R1selec == R2selec
              if n_rad != 2
                    R2selec = unsafe_trunc(Int,rand()*(n_rad-1))+1
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

            Re[3] = kp*n_mon*n_rad;
            Re[4] = ka*n_cta*n_rad;
            Re[9] = ka*n_rad*n_raft1;
            Re[13] = ktc*n_rad*(n_rad-1);
            Re[14] = ktd*n_rad*(n_rad-1);

        elseif reacselec == 14
            #Reacciones[14] += 1

            R1selec = unsafe_trunc(Int,rand()*(n_rad-1))+1
            R2selec = unsafe_trunc(Int,rand()*(n_rad-1))+1

            while R1selec == R2selec
              if n_rad != 2
                    R2selec = unsafe_trunc(Int,rand()*(n_rad-1))+1
              else
                  R2selec = 2;
              end
            end

            n_rad = n_rad - 2;
            n_pol = n_pol + 2;

            @inbounds D[R[R1selec]] += 1
            @inbounds D[R[R2selec]] += 1

            #@inbounds D[last_D] = R[R1selec];
            #@inbounds D[last_D] = R[R2selec];
            #last_D = last_D + 2;

            @inbounds R[R1selec] = R[last_R-1];
            @inbounds R[last_R-1] = 0;
            @inbounds R[R2selec] = R[last_R-2];
            @inbounds R[last_R-2] = 0;
            last_R = last_R - 2;

            Re[3] = kp*n_mon*n_rad;
            Re[4] = ka*n_cta*n_rad;
            Re[9] = ka*n_rad*n_raft1;
            Re[13] = ktc*n_rad*(n_rad-1);
            Re[14] = ktd*n_rad*(n_rad-1);

        end #if


"""
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

        #iter += 1
        conv = 100.0*(nmon_inic - n_mon)/(nmon_inic)
"""
"""
        if t_r > t_print
            println(t_print)
            t_print += 0.1
            print_counter += 1
            Concentrations[print_counter,:] = [t_r,n_mon,n_ini,n_cta,n_rad,n_istar,n_rstar,n_raft1,n_raft2,n_pol,n_nothing]./Avogadrof
            Concentrations[print_counter,1] = t_r
        end
"""
    end #while

#counters = [contadorRe[1],contadorR2,contadorR3,contadorR4,contadorR5,contadorR6,contadorR7,contadorR8,contadorR9,contadorRe[10],contadorRe[11],contadorRe[12],contadorRe[13],contadorRe[14]]

    conv = 100.0*(nmon_inic - n_mon)/(nmon_inic)
    return D,PTP1,PTP2,TP,RTP,R,conv

end #function
