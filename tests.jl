using DataFrames
using Gadfly

include("aux_functions.jl")
include("plotting_routines.jl")

function initialize_vectors()
  TP = zeros(Int32,Int(1e6))
  PTP1 = zeros(Int32,Int(1e7))
  PTP2 = zeros(Int32,Int(1e7))
  RTP = zeros(Int32,Int(1e5))
  D = zeros(Int32,10000)
  R = zeros(Int32,1000)
  return TP,PTP1,PTP2,RTP,D,R
end

TP, PTP1, PTP2, RTP, D, R = initialize_vectors();

@time MC1=MCRAFT_SF(1e6,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R)
@time MC2=MCRAFT_SF(1e7,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R)
@time MC3=MCRAFT_SF(1e8,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R)

nombres=["D","PTP1","PTP2","TP","RTP","R"]
equalize_length2(nombres,MC3)
