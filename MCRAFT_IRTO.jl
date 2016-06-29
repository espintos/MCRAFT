module Tmp2

  include("MCRAFT_IRTO_prerun.jl")

  """
  Nested functions generate allocation.
  This is why the linearwalk() function its outside MCRAFT_SF()
  """

  function linearwalk(Re,suma,r2)
  partialcumsum = 0.0
  i = 1
  while true
      partialcumsum += Re[i]
      if r2 <= partialcumsum
          return i
      end
      i += 1
  end
end

  function linearwalk_sorted(Re, sorted_indices, suma, r)
    partialcumsum = 0.0
    i = 1
    while true
        @inbounds partialcumsum += Re[sorted_indices[i]]
        if r <= partialcumsum
            return sorted_indices[i]
        end
        i += 1
    end
  end

  function get_reac(Re,r2)
    if r2 > (Re[1]+Re[2]+Re[3])+(Re[4])
      reacselec = 5;
      if r2 > (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5])
        reacselec = 6;
        if r2 > (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5]+Re[6])
          reacselec = 7;
          #if r2 > (Re[1]+Re[2]+Re[3]+Re[4]+Re[5]+Re[6]+Re[7])/(R_t)
            #reacselec = 8;
            if r2 > (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5]+Re[6])+(Re[7]+Re[8])
              reacselec = 9;
              #if r2 > (Re[1]+Re[2]+Re[3]+Re[4]+Re[5]+Re[6]+Re[7]+Re[8]+Re[9])/(R_t)
                #reacselec = 10;
                if r2 > (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5]+Re[6])+(Re[7]+Re[8]+Re[9])+(Re[10])
                  reacselec = 11;
                  if r2 > (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5]+Re[6])+(Re[7]+Re[8]+Re[9])+(Re[10]+Re[11])
                    reacselec = 12;
                    if r2 > (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5]+Re[6])+(Re[7]+Re[8]+Re[9])+(Re[10]+Re[11]+Re[12])
                      reacselec = 13;
                      if r2 > (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5]+Re[6])+(Re[7]+Re[8]+Re[9])+(Re[10]+Re[11]+Re[12])+(Re[13])
                        reacselec = 14;
                        if r2 > (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5]+Re[6])+(Re[7]+Re[8]+Re[9])+(Re[10]+Re[11]+Re[12])+(Re[13]+Re[14])
                          reacselec = 15;
                          if r2 > (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5]+Re[6])+(Re[7]+Re[8]+Re[9])+(Re[10]+Re[11]+Re[12])+(Re[13]+Re[14]+Re[15])
                            reacselec = 16;
                            if r2 > (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5]+Re[6])+(Re[7]+Re[8]+Re[9])+(Re[10]+Re[11]+Re[12])+(Re[13]+Re[14]+Re[15])+(Re[16])
                                reacselec = 17;
                            end
                          end
                        end
                      end
                    end
                  end
                end
              #end
            end
          #end
        end
      end
    elseif r2 > (Re[1]+Re[2])
      reacselec = 3;
      if r2 > (Re[1]+Re[2]+Re[3])
        reacselec = 4;
      end
    elseif r2 < Re[1] #CUIDADO ACA DI VUELTA EL SIGNO DE LA DESIGUALDAD y la reaccion
      reacselec = 1;
    else
      reacselec = 2;
    end
    return reacselec
  end


  function MCRAFT_IRTO(N,tfinal,M0,I0,CTA0,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)
    """
    Variable initialization starts here.
    When porting to Julia 0.4, pay extra attention to the integer conversion and integer functions.
    For example, instead of ifloor(x) --> change to floor(Integer, x)

    Defaults:
    M0 = 5.0
    I0 = 7.5e-3
    CTA0 = 1e-2

    M0 = 5.0
    I0 = 5e-3
    CTA0 = 5e-3

    Arrays are initialzed with zeros and for the moment, using the highest value for the largest test
    set. Eventually, the size will depend on the value of the parameter f.
    Arrays are initialized outside the function now.
    """

    tiempos=times
    matrix = indices_matrix
    matrix_i = 2
    sorted_indices = sorted_indices = matrix[:,matrix_i]
    matrix_t = tiempos[2]

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
    n_rad1 = 0
    n_rad2 = 0
    n_istar=0
    n_rstar=0
    n_raft1=0
    n_raft2=0
    n_pol = 0
    n_nothing=0

    t_r = 0.0

    """
    These are used to track the last element added to the array, right now this is faster than
    doing find(x)[end]. Will have to test findlast(x) once implemented in 0.4.
    """
    last_D = 1
    last_R = 1
    last_TP = 1
    last_RTP = 1
    last_PTP = 1

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
    kf=3.6e7
    kct=3.6e10/Avogadrof

    """
    Rate constants are calculated outside the loop initially and update only when needed, inside
    the loop. Also, only the constants that were modified are calculated.
    """

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

    Re[15] = kct*n_istar*n_raft2
    Re[16] = kct*n_rad1*n_raft2
    Re[17] = kct*n_rad2*n_raft2

    """Debugging variables"""
    """
    species = ["tiempo","M","I","CTA","R","TP*","R*","TP","PTP","D","nothing"]
    Concentrations = zeros(unsafe_trunc(Int,t_f/t_print)+1,length(species))
    Concentrations[1,:] = [t_r,M0,I0,CTA0,n_rad,n_istar,n_rstar,n_raft1,n_raft2,n_pol,n_nothing]
    print_counter = 1
    """

    #MC simluation starts here
    while t_r < t_f
      #Reaction selection - this could be improved
      R_t = (Re[1]+Re[2]+Re[3])+(Re[4]+Re[5]+Re[6])+(Re[7]+Re[8]+Re[9])+(Re[10]+Re[11]+Re[12])+(Re[13]+Re[14]+Re[15])+(Re[16]+Re[17])

      r2 = rand()*R_t;

      reacselec = linearwalk_sorted(Re, sorted_indices, R_t, r2)

      #Reaction time update (This might change in 0.4 to randexp() )
      t_r = t_r + randexp()/R_t #this is even faster

      if reacselec == 1
          #Reacciones[1] += 1

          n_ini = n_ini - 1;

          if rand() <= eff
              n_istar = n_istar + 2;
              @inbounds Re[2] = ki*n_mon*n_istar;
              @inbounds Re[5] = ka*n_istar*n_raft1;
              @inbounds Re[15] = kct*n_istar*n_raft2
          else
              n_nothing = n_nothing + 2;
          end

          @inbounds Re[1] = kd*n_ini;


      elseif reacselec == 2
          #Reacciones[2] += 1

          n_istar = n_istar - 1;
          n_mon = n_mon - 1;
          n_rad = n_rad + 1;
          n_rad1 = n_rad1 + 1;

          R[last_R] = 1;
          last_R = last_R + 1;

          @inbounds Re[2] = ki*n_mon*n_istar;
          @inbounds Re[3] = kp*n_mon*n_rad;
          @inbounds Re[4] = ka*n_cta*n_rad;
          @inbounds Re[5] = ka*n_istar*n_raft1;
          @inbounds Re[9] = ka*n_rad*n_raft1;
          @inbounds Re[13] = ktc*n_rad*(n_rad-1)
          @inbounds Re[14] = ktd*n_rad*(n_rad-1)
          @inbounds Re[15] = kct*n_istar*n_raft2
          @inbounds Re[16] = kct*n_rad1*n_raft2

      elseif reacselec == 3
          #Reacciones[3] += 1

          Rselec = unsafe_trunc(Int,n_rad*rand())+1

          if R[Rselec] == 1
              n_rad1 = n_rad1 - 1
              n_rad2 = n_rad2 + 1
              @inbounds Re[16] = kct*n_rad1*n_raft2
              @inbounds Re[17] = kct*n_rad2*n_raft2
          elseif R[Rselec] == 2
              n_rad2 = n_rad2 - 1
              @inbounds Re[17] = kct*n_rad2*n_raft2
          end

          n_mon = n_mon - 1;

          R[Rselec] += 1;

          @inbounds Re[2] = ki*n_mon*n_istar;
          @inbounds Re[3] = kp*n_mon*n_rad;

      elseif reacselec == 4
          #Reacciones[4] += 1

          Rselec = unsafe_trunc(Int,n_rad*rand())+1

          if R[Rselec] == 1
              n_rad1 = n_rad1 - 1
              @inbounds Re[16] = kct*n_rad1*n_raft2
          elseif R[Rselec] == 2
              n_rad2 = n_rad2 - 1
              @inbounds Re[17] = kct*n_rad2*n_raft2
          end
          n_rad = n_rad - 1;
          n_cta = n_cta - 1;
          n_rstar = n_rstar + 1;

          RTP[last_RTP] = R[Rselec];
          last_RTP = last_RTP + 1;

          R[Rselec] = R[last_R-1];
          R[last_R-1] = 0;
          last_R = last_R - 1;

          @inbounds Re[3] = kp*n_mon*n_rad;
          @inbounds Re[4] = ka*n_cta*n_rad;
          @inbounds Re[6] = 0.5*kf*n_rstar;
          @inbounds Re[7] = 0.5*kf*n_rstar;
          @inbounds Re[9] = ka*n_rad*n_raft1;
          @inbounds Re[13] = ktc*n_rad*(n_rad-1);
          @inbounds Re[14] = ktd*n_rad*(n_rad-1);


      elseif reacselec == 5
          #Reacciones[5] += 1

          TPselec = unsafe_trunc(Int,n_raft1*rand())+1

          n_istar = n_istar - 1;
          n_raft1 = n_raft1 - 1;
          n_rstar = n_rstar + 1;

          RTP[last_RTP] = TP[TPselec];
          last_RTP = last_RTP + 1;

          TP[TPselec] = TP[last_TP-1];
          TP[last_TP-1] = 0;
          last_TP = last_TP - 1;

          @inbounds Re[2] = ki*n_mon*n_istar;
          @inbounds Re[5] = ka*n_istar*n_raft1;
          @inbounds Re[6] = 0.5*kf*n_rstar;
          @inbounds Re[7] = 0.5*kf*n_rstar;
          @inbounds Re[9] = ka*n_rad*n_raft1;
          @inbounds Re[15] = kct*n_istar*n_raft2


      elseif reacselec == 6
         #Reacciones[6] += 1

          RTP_selec = unsafe_trunc(Int,n_rstar*rand())+1

          if RTP[RTP_selec] == 1
              n_rad1 = n_rad1 + 1
              @inbounds Re[16] = kct*n_rad1*n_raft2
          elseif RTP[RTP_selec] == 2
              n_rad2 = n_rad2 + 1
              @inbounds Re[17] = kct*n_rad2*n_raft2
          end

          n_rstar = n_rstar - 1;
          n_rad = n_rad + 1;
          n_cta = n_cta + 1;

          R[last_R] = RTP[RTP_selec];
          last_R = last_R + 1;

          RTP[RTP_selec] = RTP[last_RTP-1];
          RTP[last_RTP-1] = 0;
          last_RTP = last_RTP - 1;

          @inbounds Re[3] = kp*n_mon*n_rad;
          @inbounds Re[4] = ka*n_cta*n_rad;
          @inbounds Re[6] = 0.5*kf*n_rstar;
          @inbounds Re[7] = 0.5*kf*n_rstar;
          @inbounds Re[9] = ka*n_rad*n_raft1;
          @inbounds Re[13] = ktc*n_rad*(n_rad-1);
          @inbounds Re[14] = ktd*n_rad*(n_rad-1);


      elseif reacselec == 7
          #Reacciones[7] += 1

          RTP_selec = unsafe_trunc(Int,n_rstar*rand())+1

          n_rstar = n_rstar - 1;
          n_istar = n_istar + 1;
          n_raft1 = n_raft1 + 1;

          TP[last_TP] = RTP[RTP_selec];
          last_TP = last_TP + 1;

          RTP[RTP_selec] = RTP[last_RTP-1];
          RTP[last_RTP-1] = 0;
          last_RTP = last_RTP - 1;

          @inbounds Re[2] = ki*n_mon*n_istar;
          @inbounds Re[5] = ka*n_istar*n_raft1;
          @inbounds Re[6] = 0.5*kf*n_rstar;
          @inbounds Re[7] = 0.5*kf*n_rstar;
          @inbounds Re[9] = ka*n_rad*n_raft1;
          @inbounds Re[15] = kct*n_istar*n_raft2


      elseif reacselec == 8
          #Reacciones[8] += 1


      elseif reacselec == 9
          #Reacciones[9] += 1

          Rselec = unsafe_trunc(Int,rand()*n_rad)+1
          TPselec = unsafe_trunc(Int,n_raft1*rand())+1

          if R[Rselec] == 1
              n_rad1 = n_rad1 - 1
          elseif R[Rselec] == 2
              n_rad2 = n_rad2 - 1
          end
          n_rad = n_rad - 1;
          n_raft1 = n_raft1 - 1;
          n_raft2 = n_raft2 + 1;

          PTP1[last_PTP] = R[Rselec];
          PTP2[last_PTP] = TP[TPselec];
          last_PTP = last_PTP + 1;

          R[Rselec] = R[last_R-1];
          R[last_R-1] = 0;
          last_R = last_R - 1;

          TP[TPselec] = TP[last_TP-1];
          TP[last_TP-1] = 0;
          last_TP = last_TP - 1;

          @inbounds Re[3] = kp*n_mon*n_rad;
          @inbounds Re[4] = ka*n_cta*n_rad;
          @inbounds Re[5] = ka*n_istar*n_raft1;
          @inbounds Re[9] = ka*n_rad*n_raft1;
          @inbounds Re[11] = 0.5*kf*n_raft2;
          @inbounds Re[12] = 0.5*kf*n_raft2;
          @inbounds Re[13] = ktc*n_rad*(n_rad-1);
          @inbounds Re[14] = ktd*n_rad*(n_rad-1);
          @inbounds Re[15] = kct*n_istar*n_raft2
          @inbounds Re[16] = kct*n_rad1*n_raft2
          @inbounds Re[17] = kct*n_rad2*n_raft2


      elseif reacselec == 10
          #Reacciones[10] += 1

      elseif reacselec == 11
          #Reacciones[11] += 1

          PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

          if PTP1[PTPselec] == 1
              n_rad1 = n_rad1 + 1
          elseif PTP1[PTPselec] == 2
              n_rad2 = n_rad2 + 1
          end
          n_raft2 = n_raft2 - 1;
          n_rad = n_rad + 1;
          n_raft1 = n_raft1 + 1;

          R[last_R] = PTP1[PTPselec];
          TP[last_TP] = PTP2[PTPselec];

          PTP1[PTPselec]=PTP1[last_PTP-1];
          PTP2[PTPselec]=PTP2[last_PTP-1];
          PTP1[last_PTP-1] = 0;
          PTP2[last_PTP-1] = 0;
          last_PTP = last_PTP - 1;

          last_R = last_R + 1;
          last_TP = last_TP + 1;

          @inbounds Re[3] = kp*n_mon*n_rad;
          @inbounds Re[4] = ka*n_cta*n_rad;
          @inbounds Re[5] = ka*n_istar*n_raft1;
          @inbounds Re[9] = ka*n_rad*n_raft1;
          @inbounds Re[11] = 0.5*kf*n_raft2;
          @inbounds Re[12] = 0.5*kf*n_raft2;
          @inbounds Re[13] = ktc*n_rad*(n_rad-1);
          @inbounds Re[14] = ktd*n_rad*(n_rad-1);
          @inbounds Re[15] = kct*n_istar*n_raft2
          @inbounds Re[16] = kct*n_rad1*n_raft2
          @inbounds Re[17] = kct*n_rad2*n_raft2


      elseif reacselec == 12
          #Reacciones[12] += 1

          PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

          if PTP2[PTPselec] == 1
              n_rad1 = n_rad1 + 1
          elseif PTP2[PTPselec] == 2
              n_rad2 = n_rad2 + 1
          end
          n_raft2 = n_raft2 - 1;
          n_rad = n_rad + 1;
          n_raft1 = n_raft1 + 1;

          TP[last_TP] = PTP1[PTPselec];
          R[last_R] = PTP2[PTPselec];

          PTP1[PTPselec]=PTP1[last_PTP-1];
          PTP2[PTPselec]=PTP2[last_PTP-1];
          PTP1[last_PTP-1] = 0;
          PTP2[last_PTP-1] = 0;
          last_PTP = last_PTP - 1;

          last_R = last_R + 1;
          last_TP = last_TP + 1;

          @inbounds Re[3] = kp*n_mon*n_rad;
          @inbounds Re[4] = ka*n_cta*n_rad;
          @inbounds Re[5] = ka*n_istar*n_raft1;
          @inbounds Re[9] = ka*n_rad*n_raft1;
          @inbounds Re[11] = 0.5*kf*n_raft2;
          @inbounds Re[12] = 0.5*kf*n_raft2;
          @inbounds Re[13] = ktc*n_rad*(n_rad-1);
          @inbounds Re[14] = ktd*n_rad*(n_rad-1);
          @inbounds Re[15] = kct*n_istar*n_raft2
          @inbounds Re[16] = kct*n_rad1*n_raft2
          @inbounds Re[17] = kct*n_rad2*n_raft2

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

          if R[R1selec] == 1
              n_rad1 = n_rad1 - 1
              @inbounds Re[16] = kct*n_rad1*n_raft2
          elseif R[R1selec] == 2
              n_rad2 = n_rad2 - 1
              @inbounds Re[17] = kct*n_rad2*n_raft2
          end
          if R[R2selec] == 1
              n_rad1 = n_rad1 - 1
              @inbounds Re[16] = kct*n_rad1*n_raft2
          elseif R[R2selec] == 2
              n_rad2 = n_rad2 - 1
              @inbounds Re[17] = kct*n_rad2*n_raft2
          end
          n_rad = n_rad - 2;
          n_pol = n_pol + 1;

          #D[last_D] = R[R1selec] + R[R2selec];
          #last_D = last_D + 1;
          D[R[R1selec]+R[R2selec]] += 1

          R[R1selec] = R[last_R-1];
          R[last_R-1] = 0;
          R[R2selec] = R[last_R-2];
          R[last_R-2] = 0;
          last_R = last_R - 2;

          @inbounds Re[3] = kp*n_mon*n_rad;
          @inbounds Re[4] = ka*n_cta*n_rad;
          @inbounds Re[9] = ka*n_rad*n_raft1;
          @inbounds Re[13] = ktc*n_rad*(n_rad-1);
          @inbounds Re[14] = ktd*n_rad*(n_rad-1);

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

          if R[R1selec] == 1
              n_rad1 = n_rad1 - 1
              @inbounds Re[16] = kct*n_rad1*n_raft2
          elseif R[R1selec] == 2
              n_rad2 = n_rad2 - 1
              @inbounds Re[17] = kct*n_rad2*n_raft2
          end
          if R[R2selec] == 1
              n_rad1 = n_rad1 - 1
              @inbounds Re[16] = kct*n_rad1*n_raft2
          elseif R[R2selec] == 2
              n_rad2 = n_rad2 - 1
              @inbounds Re[17] = kct*n_rad2*n_raft2
          end
          n_rad = n_rad - 2;
          n_pol = n_pol + 2;

          D[R[R1selec]] += 1
          D[R[R2selec]] += 1

          #D[last_D] = R[R1selec];
          #D[last_D] = R[R2selec];
          #last_D = last_D + 2;

          R[R1selec] = R[last_R-1];
          R[last_R-1] = 0;
          R[R2selec] = R[last_R-2];
          R[last_R-2] = 0;
          last_R = last_R - 2;

          @inbounds Re[3] = kp*n_mon*n_rad;
          @inbounds Re[4] = ka*n_cta*n_rad;
          @inbounds Re[9] = ka*n_rad*n_raft1;
          @inbounds Re[13] = ktc*n_rad*(n_rad-1);
          @inbounds Re[14] = ktd*n_rad*(n_rad-1);

      elseif reacselec == 15

          PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

          n_istar = n_istar - 1
          n_raft2 = n_raft2 -1
          n_pol = n_pol + 1

          #D[last_D] = R[R1selec] + R[R2selec];
          #last_D = last_D + 1;
          D[PTP1[PTPselec]+PTP2[PTPselec]] += 1

          PTP1[PTPselec]=PTP1[last_PTP-1]
          PTP2[PTPselec]=PTP2[last_PTP-1]
          PTP1[last_PTP-1] = 0
          PTP2[last_PTP-1] = 0
          last_PTP = last_PTP - 1

          @inbounds Re[2] = ki*n_mon*n_istar;
          @inbounds Re[5] = ka*n_istar*n_raft1;
          @inbounds Re[11] = 0.5*kf*n_raft2 #same than @inbounds Re[6]
          @inbounds Re[12] = 0.5*kf*n_raft2 #same than @inbounds Re[6]
          @inbounds Re[15] = kct*n_istar*n_raft2
          @inbounds Re[16] = kct*n_rad1*n_raft2
          @inbounds Re[17] = kct*n_rad2*n_raft2

      elseif reacselec == 16
          #Reacciones[16] += 1

          Rselec = findfirst(R,1)
          PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

          n_rad1 = n_rad1 - 1
          n_rad = n_rad - 1
          n_raft2 = n_raft2 -1
          n_pol = n_pol + 1

          if R[Rselec] != 1
              println("error en @inbounds Re[16]")
              break
          end

          #D[last_D] = R[R1selec] + R[R1selec];
          #last_D = last_D + 1;
          D[1+PTP1[PTPselec]+PTP2[PTPselec]] += 1

          PTP1[PTPselec]=PTP1[last_PTP-1]
          PTP2[PTPselec]=PTP2[last_PTP-1]
          PTP1[last_PTP-1] = 0
          PTP2[last_PTP-1] = 0
          last_PTP = last_PTP - 1

          R[Rselec] = R[last_R-1];
          R[last_R-1] = 0;
          last_R = last_R - 1;

          @inbounds Re[3] = kp*n_mon*n_rad
          @inbounds Re[4] = ka*n_cta*n_rad
          @inbounds Re[9] = ka*n_rad*n_raft1
          @inbounds Re[11] = 0.5*kf*n_raft2 #same than @inbounds Re[6]
          @inbounds Re[12] = 0.5*kf*n_raft2 #same than @inbounds Re[6]
          @inbounds Re[13] = ktc*n_rad*(n_rad-1)
          @inbounds Re[14] = ktd*n_rad*(n_rad-1)
          @inbounds Re[15] = kct*n_istar*n_raft2
          @inbounds Re[16] = kct*n_rad1*n_raft2
          @inbounds Re[17] = kct*n_rad2*n_raft2

      elseif reacselec == 17
          #Reacciones[17] += 1

          Rselec = findfirst(R,2)
          PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

          n_rad2 = n_rad2 - 1
          n_rad = n_rad - 1
          n_raft2 = n_raft2 -1
          n_pol = n_pol + 1

          if R[Rselec] != 2
              println("error en @inbounds Re[17]")
              break
          end

          #D[last_D] = R[R1selec] + R[R1selec];
          #last_D = last_D + 1;
          D[2+PTP1[PTPselec]+PTP2[PTPselec]] += 1

          PTP1[PTPselec]=PTP1[last_PTP-1]
          PTP2[PTPselec]=PTP2[last_PTP-1]
          PTP1[last_PTP-1] = 0
          PTP2[last_PTP-1] = 0
          last_PTP = last_PTP - 1

          R[Rselec] = R[last_R-1];
          R[last_R-1] = 0;
          last_R = last_R - 1;

          @inbounds Re[3] = kp*n_mon*n_rad
          @inbounds Re[4] = ka*n_cta*n_rad
          @inbounds Re[9] = ka*n_rad*n_raft1
          @inbounds Re[11] = 0.5*kf*n_raft2 #same than @inbounds Re[6]
          @inbounds Re[12] = 0.5*kf*n_raft2 #same than @inbounds Re[6]
          @inbounds Re[13] = ktc*n_rad*(n_rad-1)
          @inbounds Re[14] = ktd*n_rad*(n_rad-1)
          @inbounds Re[15] = kct*n_istar*n_raft2
          @inbounds Re[16] = kct*n_rad1*n_raft2
          @inbounds Re[17] = kct*n_rad2*n_raft2

      end #if


      if t_r > matrix_t
        matrix_i += 1
        sorted_indices = sorted_indices = matrix[:,matrix_i]
        matrix_t = tiempos[matrix_i]
      end

    end #while

    #counters = [contadorRe[1],contadorRe[2],contadorRe[3],contadorRe[4],contadorRe[5],contadorRe[6],contadorRe[7],contadorRe[8],contadorRe[9],contadorRe[10],contadorRe[11],contadorRe[12],contadorRe[13],contadorRe[14]]
    conv = 100.0*(nmon_inic - n_mon)/(nmon_inic)
    return Vector{Int32}[D,PTP1,PTP2,TP,RTP,R]

  end #function

end #module
