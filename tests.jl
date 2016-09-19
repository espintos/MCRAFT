using DataFrames
using Gadfly

include("aux_functions.jl")
include("plotting_routines.jl")

function initialize_vectors()
  TP = zeros(Int32,Int(3e7))
  PTP1 = zeros(Int32,Int(3e7))
  PTP2 = zeros(Int32,Int(3e7))
  RTP = zeros(Int32,Int(1e7))
  D = zeros(Int32,100000)
  R = zeros(Int32,10000)
  Re = zeros(Float64,17)
  return TP,PTP1,PTP2,RTP,D,R,Re
end

N=[1e6,1e7,1e8,1e9]

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
times,indices_matrix=MCRAFT_SF_prerun(1e8,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re)
indices_matrix
times

function MC_bench(N,samples)
  a=zeros(Float64,Int(samples))
  for i = 1:samples
    TP, PTP1, PTP2, RTP, D, R = initialize_vectors()
    a[i] = @elapsed MCRAFT_SF(N,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)
  end
  return a
end

results1e9 = MC_bench(1e9,100)
clipboard(results1e9)

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
@time MC1=MCRAFT_SF(1e6,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
@time MC2=MCRAFT_SF(1e7,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
@time MC3=MCRAFT_SF(1e8,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
@time MC4=MCRAFT_SF(1e9,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
@time MC5=MCRAFT_SF(1e10,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)

nombres=["D","PTP","TP","RTP","R"]
DF1=equalize_length2(nombres,MC1)
DF2=equalize_length2(nombres,MC2)
DF3=equalize_length2(nombres,MC3)
DF4=equalize_length2(nombres,MC4)

multipleplot(N[2:end],DF2[1],DF3[1],DF4[1])

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
Profile.clear()
@profile MC3=MCRAFT_SF(1e8,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)

Profile.print()


"""
1e6   0.094027 seconds (36 allocations: 3.250 KB)
1e7   0.986682 seconds (37 allocations: 3.297 KB)
1e8  12.263587 seconds (37 allocations: 3.297 KB)
1e9 169.239320 seconds (37 allocations: 3.297 KB)

1e6      0.083198 seconds (30 allocations: 2.688 KB)
1e7      0.866330 seconds (31 allocations: 2.734 KB)
1e8      9.294375 seconds (31 allocations: 2.734 KB)
1e9    146.686201 seconds (31 allocations: 2.734 KB)
1e10  2051.214733 seconds (25 allocations: 2.172 KB)
"""
