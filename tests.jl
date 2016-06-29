using DataFrames
using Gadfly

include("aux_functions.jl")
include("plotting_routines.jl")

function initialize_vectors()
  TP = zeros(Int32,Int(1e7))
  PTP1 = zeros(Int32,Int(3e7))
  PTP2 = zeros(Int32,Int(3e7))
  RTP = zeros(Int32,Int(1e7))
  D = zeros(Int32,100000)
  R = zeros(Int32,10000)
  return TP,PTP1,PTP2,RTP,D,R
end

N=[1e6,1e7,1e8,1e9]

TP, PTP1, PTP2, RTP, D, R = initialize_vectors()
@time MC1=MCRAFT_SF(1e6,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R);

TP, PTP1, PTP2, RTP, D, R = initialize_vectors()
@time MC2=MCRAFT_SF(1e7,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R);

TP, PTP1, PTP2, RTP, D, R = initialize_vectors()
@time MC3=MCRAFT_SF(1e8,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R)

TP, PTP1, PTP2, RTP, D, R = initialize_vectors()
@time MC4=MCRAFT_SF(1e9,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R)

nombres=["D","PTP","TP","RTP","R"]
DF1=equalize_length2(nombres,MC1)
DF2=equalize_length2(nombres,MC2)
DF3=equalize_length2(nombres,MC3)
DF4=equalize_length2(nombres,MC4)

multipleplot(N[2:end],DF2[1],DF3[1],DF4[1])

"""
Estos son convirtiendo a vector
1e6    .082439 seconds (6 allocations: 448 bytes)
1e7   0.797506 seconds (7 allocations: 496 bytes)
1e8   8.879239 seconds (7 allocations: 496 bytes)
1e9 178.099637 seconds (7 allocations: 496 bytes)


Estos usando la funcion get_reac
1e6   0.091892 seconds (6 allocations: 448 bytes)
1e7   0.866231 seconds (7 allocations: 496 bytes)
1e8   9.929670 seconds (7 allocations: 496 bytes)
1e9 182.578304 seconds (7 allocations: 496 bytes)
"""

@profile MCRAFT_SF(1e8,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R)

Profile.print()
