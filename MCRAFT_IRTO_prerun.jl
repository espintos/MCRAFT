  function MCRAFT_IRTO_prerun(N,tfinal,M0,I0,CTA0,TP,PTP1,PTP2,RTP,D,R,Re)

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

      aux=1

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
      #Re = zeros(Float64,14)
      R1 = kd*n_ini
      R2 = ki*n_mon*n_istar
      R3 = kp*n_mon*n_rad
      R4 = ka*n_cta*n_rad
      R5 = ka*n_istar*n_raft1
      R6 = 0.5*kf*n_rstar #the 0.5 is there because RAFT* reacts with 0.5 chance to R6] and 0.5 to R7]
      R7 = 0.5*kf*n_rstar #same than R6]
      R8 = 0.0; #Since we are sure this reaction doesnt go anymore, it could be removed
      R9 = ka*n_rad*n_raft1
      R10 = 0.0; #same than R8]
      R11 = 0.5*kf*n_raft2 #same than R6]
      R12 = 0.5*kf*n_raft2 #same than R6]
      R13 = ktc*n_rad*(n_rad-1)
      R14 = ktd*n_rad*(n_rad-1)

      R15 = kct*n_istar*n_raft2
      R16 = kct*n_rad1*n_raft2
      R17 = kct*n_rad2*n_raft2

      """Debugging variables"""
      time_conv_reacs=zeros(19)
      time_conv_indices=zeros(19)
      Reactions = zeros(17)
      interval_fixed = 0.01
      interval = interval_fixed
      iter = 0
      sorted_indices = sortperm(Reactions,rev=true)
      """
      #maxlong=[0,0,0,0,0,0]
      species = ["tiempo","M","I","CTA","R","TP*","R*","TP","PTP","D","nothing"]
      Concentrations = zeros(ceil(Integer,t_f/t_print)+1,length(species)) ## VERIFICAR ESTA LINEA
      Concentrations[1,:] = [t_r,M0,I0,CTA0,n_rad,n_istar,n_rstar,n_raft1,n_raft2,n_pol,n_nothing]
      print_counter = 1
      """

      #MC simluation starts here
      while t_r < t_f

          #iter += 1

          #Reaction selection - this could be improved
          R_t = (R1+R2+R3)+(R4+R5+R6)+(R7+R8+R9)+(R10+R11+R12)+(R13+R14+R15)+(R16+R17);

          r2 = rand()*R_t;

          if r2 > (R1+R2)+(R3+R4)
              reacselec = 5;
              if r2 > (R1+R2+R3)+(R4+R5)
                  reacselec = 6;
                  if r2 > (R1+R2+R3)+(R4+R5+R6)
                      reacselec = 7;
                      #if r2 > (R1+R2+R3+R4+R5+R6+R7)/(R_t)
                          #reacselec = 8;
                          if r2 > (R1+R2+R3+R4+R5+R6+R7+R8)
                              reacselec = 9;
                              #if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9)/(R_t)
                                  #reacselec = 10;
                                  if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10)
                                      reacselec = 11;
                                      if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11)
                                          reacselec = 12;
                                          if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12)
                                              reacselec = 13;
                                              if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13)
                                                  reacselec = 14;
                                                  if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13+R14)
                                                      reacselec = 15;
                                                      if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13+R14+R15)
                                                          reacselec = 16;
                                                              if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13+R14+R15+R16)
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
          elseif r2 > (R1+R2)
              reacselec = 3;
              if r2 > (R1+R2+R3)
                  reacselec = 4;
              end
          elseif r2 < R1 #CUIDADO ACA DI VUELTA EL SIGNO DE LA DESIGUALDAD y la reaccion
              reacselec = 1;
          else
              reacselec = 2;
          end


          #Reaction time update (This might change in 0.4 to randexp() )
          t_r = t_r + randexp()/R_t #this is even faster

          if reacselec == 1
              Reactions[1] += 1

              n_ini = n_ini - 1

              if rand() <= eff
                  n_istar = n_istar + 2
                  R2 = ki*n_mon*n_istar;
                  R5 = ka*n_istar*n_raft1;
                  R15 = kct*n_istar*n_raft2
              else
                  n_nothing = n_nothing + 2;
              end

              R1 = kd*n_ini;


          elseif reacselec == 2
              Reactions[2] += 1

              n_istar = n_istar - 1;
              n_mon = n_mon - 1;
              n_rad = n_rad + 1;
              n_rad1 = n_rad1 + 1;

              @inbounds R[last_R] = 1;
              last_R = last_R + 1;

              R2 = ki*n_mon*n_istar;
              R3 = kp*n_mon*n_rad;
              R4 = ka*n_cta*n_rad;
              R5 = ka*n_istar*n_raft1;
              R9 = ka*n_rad*n_raft1;
              R13 = ktc*n_rad*(n_rad-1)
              R14 = ktd*n_rad*(n_rad-1)
              R15 = kct*n_istar*n_raft2
              R16 = kct*n_rad1*n_raft2

          elseif reacselec == 3
              Reactions[3] += 1

              Rselec = unsafe_trunc(Int,n_rad*rand())+1

              if R[Rselec] == 1
                  n_rad1 = n_rad1 - 1
                  n_rad2 = n_rad2 + 1
                  R16 = kct*n_rad1*n_raft2
                  R17 = kct*n_rad2*n_raft2
              elseif R[Rselec] == 2
                  n_rad2 = n_rad2 - 1
                  R17 = kct*n_rad2*n_raft2
              end

              n_mon = n_mon - 1;

              @inbounds R[Rselec] += 1;

              R2 = ki*n_mon*n_istar;
              R3 = kp*n_mon*n_rad;

          elseif reacselec == 4
              Reactions[4] += 1

              Rselec = unsafe_trunc(Int,n_rad*rand())+1

              if R[Rselec] == 1
                  n_rad1 = n_rad1 - 1
                  R16 = kct*n_rad1*n_raft2
              elseif R[Rselec] == 2
                  n_rad2 = n_rad2 - 1
                  R17 = kct*n_rad2*n_raft2
              end
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
              R13 = ktc*n_rad*(n_rad-1);
              R14 = ktd*n_rad*(n_rad-1);


          elseif reacselec == 5
              Reactions[5] += 1

              TPselec = unsafe_trunc(Int,n_raft1*rand())+1

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
              R15 = kct*n_istar*n_raft2


          elseif reacselec == 6
             Reactions[6] += 1

              RTP_selec = unsafe_trunc(Int,n_rstar*rand())+1

              if RTP[RTP_selec] == 1
                  n_rad1 = n_rad1 + 1
                  R16 = kct*n_rad1*n_raft2
              elseif RTP[RTP_selec] == 2
                  n_rad2 = n_rad2 + 1
                  R17 = kct*n_rad2*n_raft2
              end

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
              R13 = ktc*n_rad*(n_rad-1);
              R14 = ktd*n_rad*(n_rad-1);


          elseif reacselec == 7
              Reactions[7] += 1

              RTP_selec = unsafe_trunc(Int,n_rstar*rand())+1

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
              R15 = kct*n_istar*n_raft2


          elseif reacselec == 8
              Reactions[8] += 1


          elseif reacselec == 9
              Reactions[9] += 1

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
              R13 = ktc*n_rad*(n_rad-1);
              R14 = ktd*n_rad*(n_rad-1);
              R15 = kct*n_istar*n_raft2
              R16 = kct*n_rad1*n_raft2
              R17 = kct*n_rad2*n_raft2


          elseif reacselec == 10
              Reactions[10] += 1

          elseif reacselec == 11
              Reactions[11] += 1

              PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

              if PTP1[PTPselec] == 1
                  n_rad1 = n_rad1 + 1
              elseif PTP1[PTPselec] == 2
                  n_rad2 = n_rad2 + 1
              end
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
              R13 = ktc*n_rad*(n_rad-1);
              R14 = ktd*n_rad*(n_rad-1);
              R15 = kct*n_istar*n_raft2
              R16 = kct*n_rad1*n_raft2
              R17 = kct*n_rad2*n_raft2


          elseif reacselec == 12
              Reactions[12] += 1

              PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

              if PTP2[PTPselec] == 1
                  n_rad1 = n_rad1 + 1
              elseif PTP2[PTPselec] == 2
                  n_rad2 = n_rad2 + 1
              end
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
              R13 = ktc*n_rad*(n_rad-1);
              R14 = ktd*n_rad*(n_rad-1);
              R15 = kct*n_istar*n_raft2
              R16 = kct*n_rad1*n_raft2
              R17 = kct*n_rad2*n_raft2

          elseif reacselec == 13
              Reactions[13] += 1

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
                  R16 = kct*n_rad1*n_raft2
              elseif R[R1selec] == 2
                  n_rad2 = n_rad2 - 1
                  R17 = kct*n_rad2*n_raft2
              end
              if R[R2selec] == 1
                  n_rad1 = n_rad1 - 1
                  R16 = kct*n_rad1*n_raft2
              elseif R[R2selec] == 2
                  n_rad2 = n_rad2 - 1
                  R17 = kct*n_rad2*n_raft2
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
              R13 = ktc*n_rad*(n_rad-1);
              R14 = ktd*n_rad*(n_rad-1);

          elseif reacselec == 14
              Reactions[14] += 1

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
                  R16 = kct*n_rad1*n_raft2
              elseif R[R1selec] == 2
                  n_rad2 = n_rad2 - 1
                  R17 = kct*n_rad2*n_raft2
              end
              if R[R2selec] == 1
                  n_rad1 = n_rad1 - 1
                  R16 = kct*n_rad1*n_raft2
              elseif R[R2selec] == 2
                  n_rad2 = n_rad2 - 1
                  R17 = kct*n_rad2*n_raft2
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

              R3 = kp*n_mon*n_rad;
              R4 = ka*n_cta*n_rad;
              R9 = ka*n_rad*n_raft1;
              R13 = ktc*n_rad*(n_rad-1);
              R14 = ktd*n_rad*(n_rad-1);

          elseif reacselec == 15
            Reactions[15] += 1

              PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

              n_istar = n_istar - 1
              n_raft2 = n_raft2 -1
              n_pol = n_pol + 1

              #@inbounds D[last_D] = R[R1selec] + R[R2selec];
              #last_D = last_D + 1;
              @inbounds D[PTP1[PTPselec]+PTP2[PTPselec]] += 1

              @inbounds PTP1[PTPselec]=PTP1[last_PTP-1]
              @inbounds PTP2[PTPselec]=PTP2[last_PTP-1]
              @inbounds PTP1[last_PTP-1] = 0
              @inbounds PTP2[last_PTP-1] = 0
              last_PTP = last_PTP - 1

              R2 = ki*n_mon*n_istar;
              R5 = ka*n_istar*n_raft1;
              R11 = 0.5*kf*n_raft2 #same than R6
              R12 = 0.5*kf*n_raft2 #same than R6
              R15 = kct*n_istar*n_raft2
              R16 = kct*n_rad1*n_raft2
              R17 = kct*n_rad2*n_raft2

          elseif reacselec == 16
              Reactions[16] += 1

              Rselec = findfirst(R,1)
              PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

              n_rad1 = n_rad1 - 1
              n_rad = n_rad - 1
              n_raft2 = n_raft2 -1
              n_pol = n_pol + 1

              if R[Rselec] != 1
                  println("error en R16")
                  break
              end

              #@inbounds D[last_D] = R[R1selec] + R[R2selec];
              #last_D = last_D + 1;
              @inbounds D[1+PTP1[PTPselec]+PTP2[PTPselec]] += 1

              @inbounds PTP1[PTPselec]=PTP1[last_PTP-1]
              @inbounds PTP2[PTPselec]=PTP2[last_PTP-1]
              @inbounds PTP1[last_PTP-1] = 0
              @inbounds PTP2[last_PTP-1] = 0
              last_PTP = last_PTP - 1

              @inbounds R[Rselec] = R[last_R-1];
              @inbounds R[last_R-1] = 0;
              last_R = last_R - 1;

              R3 = kp*n_mon*n_rad
              R4 = ka*n_cta*n_rad
              R9 = ka*n_rad*n_raft1
              R11 = 0.5*kf*n_raft2 #same than R6
              R12 = 0.5*kf*n_raft2 #same than R6
              R13 = ktc*n_rad*(n_rad-1)
              R14 = ktd*n_rad*(n_rad-1)
              R15 = kct*n_istar*n_raft2
              R16 = kct*n_rad1*n_raft2
              R17 = kct*n_rad2*n_raft2

          elseif reacselec == 17
              Reactions[17] += 1

              Rselec = findfirst(R,2)
              PTPselec = unsafe_trunc(Int,n_raft2*rand())+1

              n_rad2 = n_rad2 - 1
              n_rad = n_rad - 1
              n_raft2 = n_raft2 -1
              n_pol = n_pol + 1

              if R[Rselec] != 2
                  println("error en R17")
                  break
              end

              #@inbounds D[last_D] = R[R1selec] + R[R2selec];
              #last_D = last_D + 1;
              @inbounds D[2+PTP1[PTPselec]+PTP2[PTPselec]] += 1

              @inbounds PTP1[PTPselec]=PTP1[last_PTP-1]
              @inbounds PTP2[PTPselec]=PTP2[last_PTP-1]
              @inbounds PTP1[last_PTP-1] = 0
              @inbounds PTP2[last_PTP-1] = 0
              last_PTP = last_PTP - 1

              @inbounds R[Rselec] = R[last_R-1];
              @inbounds R[last_R-1] = 0;
              last_R = last_R - 1;

              R3 = kp*n_mon*n_rad
              R4 = ka*n_cta*n_rad
              R9 = ka*n_rad*n_raft1
              R11 = 0.5*kf*n_raft2 #same than R6
              R12 = 0.5*kf*n_raft2 #same than R6
              R13 = ktc*n_rad*(n_rad-1)
              R14 = ktd*n_rad*(n_rad-1)
              R15 = kct*n_istar*n_raft2
              R16 = kct*n_rad1*n_raft2
              R17 = kct*n_rad2*n_raft2

          end #if

          conv = 100.0*(nmon_inic - n_mon)/(nmon_inic)

          if t_r > interval
            new_sorted_indices = sortperm(Reactions, rev=true)
            if new_sorted_indices[1:4] != sorted_indices[1:4]
              time_conv_reacs=hcat(time_conv_reacs,vcat(t_r,conv,Reactions))
              time_conv_indices=hcat(time_conv_indices,vcat(t_r,conv,new_sorted_indices))
              interval += interval_fixed*0.5
              Reactions=zeros(17)
            else
              interval += interval_fixed*2
            end
            sorted_indices = new_sorted_indices
          end

          if t_r > aux
            aux +=1
            println(t_r)
          end


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

      #counters = [contadorR1,contadorR2,contadorR3,contadorR4,contadorR5,contadorR6,contadorR7,contadorR8,contadorR9,contadorR10,contadorR11,contadorR12,contadorR13,contadorR14]
      conv = 100.0*(nmon_inic - n_mon)/(nmon_inic)
      #return Vector{Int32}[D,PTP1,PTP2,TP,RTP,R]
      #return time_conv_reacs,time_conv_indices

      times = squeeze(time_conv_indices[1,:],1)
      push!(times,t_f+1.0)
      matrix_indices = int(time_conv_indices[3:end,:])
      matrix_indices = hcat(matrix_indices,matrix_indices[:,end]) #duplicates las column
      return times,matrix_indices

  end #function
