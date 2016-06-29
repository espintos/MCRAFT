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
times,indices_matrix=MCRAFT_IRTO_prerun(1e8,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re)
indices_matrix
times

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
@time MC1=MCRAFT_IRTO(1e6,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
@time MC2=MCRAFT_IRTO(1e7,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
@time MC3=MCRAFT_IRTO(1e8,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
@time MC4=MCRAFT_IRTO(1e9,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re,times,indices_matrix)

nombres=["D","PTP","TP","RTP","R"]
DF1=equalize_length2(nombres,MC1)
DF2=equalize_length2(nombres,MC2)
DF3=equalize_length2(nombres,MC3)
DF4=equalize_length2(nombres,MC4)

multipleplot(N[2:end],DF2[1],DF3[1],DF4[1])

@profile MC3=MCRAFT_SF(1e8,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re)

Profile.print()

a=7.2712999512e7
r=7.219919099983637e7
A=[35.855999999999995,0.0,3.5991667655999996e7,3.599985599999999e7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,360719.99999999994,360719.99999999994,0.0,0.0,0.0]

get_reac(A,r)
