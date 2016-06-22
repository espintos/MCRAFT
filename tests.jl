using DataFrames
using Gadfly

include("aux_functions.jl")
include("plotting_routines.jl")

using StatsBase
function initialize_vectors()
  TP = zeros(Int32,Int(1e6))
  PTP1 = zeros(Int32,Int(3e6))
  PTP2 = zeros(Int32,Int(3e6))
  RTP = zeros(Int32,Int(1e6))
  D = zeros(Int32,100000)
  R = zeros(Int32,10000)
  Re = zeros(Float64,14)
  return TP,PTP1,PTP2,RTP,D,R,Re
end

N=[1e6,1e7,1e8,1e9]

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
@time MC1=MCRAFT_SF(1e6,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re)

TP, PTP1, PTP2, RTP, D, R = initialize_vectors()
@time MC2=MCRAFT_SF(1e7,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re);

TP, PTP1, PTP2, RTP, D, R = initialize_vectors()
@time MC3=MCRAFT_SF(1e8,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re)

TP, PTP1, PTP2, RTP, D, R = initialize_vectors()
@time MC4=MCRAFT_SF(1e9,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re)

nombres=["D","PTP","TP","RTP","R"]
DF1=equalize_length2(nombres,MC1)
DF2=equalize_length2(nombres,MC2)
DF3=equalize_length2(nombres,MC3)
DF4=equalize_length2(nombres,MC4)

multipleplot(N[2:end],DF2[1],DF3[1],DF4[1])

Profile.clear()
@profile MC2=MCRAFT_SF(1e7,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re);

Profile.print()

@allocated MC2=MCRAFT_SF(1e7,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re);
Profile.clear_malloc_data()
