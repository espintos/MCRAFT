Program MCRAFT
  ! NOT SURE...
  Implicit none

  ! declarations
  double precision :: M0, I0, CTA0, N, f, avo, avof, tf, tfinal, tr, tprint
  double precision :: eff, kd, kp, ki, ktd, ktc, ka, kf
  double precision :: R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,Rtotal
  double precision :: ra1, ra2, ra3, ra4, conv
  integer :: TP(30), PTP1(20), PTP2(20), RTP(20), D(10), R(20)
  integer :: nmon, nini, ncta, nmoni, ninii, nctai, ntotal
  integer :: nrad, nistar, nrstar, nraft1, nraft2, npol, nnothing
  integer :: lastd, lastr, lasttp, lastrtp, lastptp
  integer :: reacselec
  parameter (avo = 6.022141986e23)

  !Initializations
  M0 = 5.0
  I0 = 7.5e-3
  CTA0 = 0.01
  N = 1e8
  tfinal = 35

  f = N/(M0*avo+I0*avo+CTA0*avo)
  avof = avo * f
  tf = tfinal

  nmon = floor(M0*avof)
  nini = floor(I0*avof)
  ncta = floor(CTA0*avof)
  nmoni = nmon
  ninii = nini
  nctai = ncta
  ntotal = nmon+nini+ncta

  nrad = 0
  nistar = 0
  nrstar = 0
  nraft1 = 0
  nraft2 = 0
  npol = 0
  nnothing = 0

  tr = 0

  lastd = 1
  lastr = 1
  lasttp = 1
  lastrtp = 1
  lastptp = 1

  tprint = 0.1

  data TP /30*0/
  data PTP1 /20*0/
  data PTP2 /20*0/
  data RTP /20*0/
  data D /10*0/
  data R /20*0/

  eff = 0.5
  kd = 0.036
  kp = 3600000/avof
  ki = kp
  ktd = 3.6e10/avof
  ktc = 3.6e10/avof
  ka = 3.6e9/avof
  kf = 36

  R1 = kd*nini
  R2 = ki*nmon*nistar
  R3 = kp*nmon*nrad
  R4 = ka*ncta*nrad
  R5 = ka*nistar*nraft1
  R6 = 0.5*kf*nrstar
  R7 = 0.5*kf*nrstar
  R8 = 0.0
  R9 = ka*nrad*nraft1
  R10 = 0.0
  R11 = 0.5*kf*nraft2
  R12 = 0.5*kf*nraft2
  R13 = ktc*nrad*(nrad-1)
  R14 = ktd*nrad*(nrad-1)


  call init_random_seed()

do while (tr<tf)

    Rtotal = R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13+R14

    call random_number(ra1)

    if (ra1 > (R1+R2+R3+R4)/(Rtotal)) then
       reacselec = 5
       if (ra1 > (R1+R2+R3+R4+R5)/(Rtotal)) then
          reacselec = 6
          if (ra1 > (R1+R2+R3+R4+R5+R6)/(Rtotal)) then
             reacselec = 7
!            #if r2 > (R1+R2+R3+R4+R5+R6+R7)/(R_t)
!               #reacselec = 8
                if (ra1 > (R1+R2+R3+R4+R5+R6+R7+R8)/(Rtotal)) then
                   reacselec = 9
!                  #if r2 > (R1+R2+R3+R4+R5+R6+R7+R8+R9)/(R_t)
!                     #reacselec = 10
                      if (ra1 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10)/(Rtotal)) then
                         reacselec = 11
                         if (ra1 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11)/(Rtotal)) then
                            reacselec = 12
                            if (ra1 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12)/(Rtotal)) then
                               reacselec = 13
                               if (ra1 > (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13)/(Rtotal)) then
                                  reacselec = 14
                               endif
                            endif
                         endif
                      endif
!                 #end
                endif
!           #end
          endif
       endif
    elseif (ra1 > (R1+R2)/(Rtotal)) then
       reacselec = 3
       if (ra1 > (R1+R2+R3)/(Rtotal)) then
          reacselec = 4
       endif
    elseif (ra1 < R1/(Rtotal)) then
       reacselec = 1
    else
       reacselec = 2
    endif

    call random_number(ra2)
    tr = tr + (-log(ra2))/(Rtotal)

    if (reacselec == 1) then
        nini = nini - 1
        call random_number(ra3)
        if (ra3<=eff) then
            nistar = nistar + 2
            R2 = ki*nmon*nistar
            R5 = ka*nistar*nraft1
        else
            nnothing = nnothing + 2
        endif
        R1 = kd*nini

    elseif (reacselec == 2) then
        nistar = nistar - 1
        nmon = nmon - 1
        nrad = nrad + 1

        R2 = ki*nmon*nistar
        R3 = kp*nmon*nrad
        R4 = ka*ncta*nrad
        R5 = ka*nistar*nraft1
        R9 = ka*nrad*nraft1
        R13 = ktc*nrad*(nrad-1)
        R14 = ktd*nrad*(nrad-1)

    elseif (reacselec == 3) then
        nmon = nmon - 1

        R2 = ki*nmon*nistar
        R3 = kp*nmon*nrad

    elseif (reacselec == 4) then
        nrad = nrad - 1
        ncta = ncta - 1
        nrstar = nrstar + 1

        R3 = kp*nmon*nrad
        R4 = ka*ncta*nrad
        R6 = 0.5*kf*nrstar
        R7 = 0.5*kf*nrstar
        R9 = ka*nrad*nraft1
        R13 = ktc*nrad*(nrad-1)
        R14 = ktd*nrad*(nrad-1)

    elseif (reacselec == 5) then
        nistar = nistar - 1
        nraft1 = nraft1 - 1
        nrstar = nrstar + 1

        R2 = ki*nmon*nistar
        R5 = ka*nistar*nraft1
        R6 = 0.5*kf*nrstar
        R7 = 0.5*kf*nrstar
        R9 = ka*nrad*nraft1

    elseif (reacselec == 6) then
        nrstar = nrstar - 1
        nrad = nrad + 1
        ncta = ncta + 1

        R3 = kp*nmon*nrad
        R4 = ka*ncta*nrad
        R6 = 0.5*kf*nrstar
        R7 = 0.5*kf*nrstar
        R9 = ka*nrad*nraft1
        R13 = ktc*nrad*(nrad-1)
        R14 = ktd*nrad*(nrad-1)

    elseif (reacselec == 7) then
        nrstar = nrstar - 1
        nistar = nistar + 1
        nraft1 = nraft1 + 1

        R2 = ki*nmon*nistar
        R5 = ka*nistar*nraft1
        R6 = 0.5*kf*nrstar
        R7 = 0.5*kf*nrstar
        R9 = ka*nrad*nraft1

    elseif (reacselec == 8) then
        !Nothing

    elseif (reacselec == 9) then
        nrad = nrad - 1
        nraft1 = nraft1 - 1
        nraft2 = nraft2 + 1

        R3 = kp*nmon*nrad
        R4 = ka*ncta*nrad
        R5 = ka*nistar*nraft1
        R9 = ka*nrad*nraft1
        R11 = 0.5*kf*nraft2
        R12 = 0.5*kf*nraft2
        R13 = ktc*nrad*(nrad-1)
        R14 = ktd*nrad*(nrad-1)

    elseif (reacselec == 10) then
        !Nothing

    elseif (reacselec == 11) then
        nraft2 = nraft2 - 1
        nrad = nrad + 1
        nraft1 = nraft1 + 1

        R3 = kp*nmon*nrad
        R4 = ka*ncta*nrad
        R5 = ka*nistar*nraft1
        R9 = ka*nrad*nraft1
        R11 = 0.5*kf*nraft2
        R12 = 0.5*kf*nraft2
        R13 = ktc*nrad*(nrad-1)
        R14 = ktd*nrad*(nrad-1)

    elseif (reacselec == 12) then
        nraft2 = nraft2 - 1
        nrad = nrad + 1
        nraft1 = nraft1 + 1

        R3 = kp*nmon*nrad
        R4 = ka*ncta*nrad
        R5 = ka*nistar*nraft1
        R9 = ka*nrad*nraft1
        R11 = 0.5*kf*nraft2
        R12 = 0.5*kf*nraft2
        R13 = ktc*nrad*(nrad-1)
        R14 = ktd*nrad*(nrad-1)

    elseif (reacselec == 13) then
        nrad = nrad - 2
        npol = npol + 1

        R3 = kp*nmon*nrad
        R4 = ka*ncta*nrad
        R9 = ka*nrad*nraft1
        R13 = ktc*nrad*(nrad-1)
        R14 = ktd*nrad*(nrad-1)

    elseif (reacselec == 14) then
        nrad = nrad - 2
        npol = npol + 2

        R3 = kp*nmon*nrad
        R4 = ka*ncta*nrad
        R9 = ka*nrad*nraft1
        R13 = ktc*nrad*(nrad-1)
        R14 = ktd*nrad*(nrad-1)
    endif

conv = 100.0*(nmoni - nmon)/(nmoni)

enddo




  print *, conv

  end program MCRAFT
