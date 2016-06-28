function initialize_vectors()
  TP = zeros(Int32,Int(1e7))
  PTP1 = zeros(Int32,Int(3e7))
  PTP2 = zeros(Int32,Int(3e7))
  RTP = zeros(Int32,Int(1e7))
  D = zeros(Int32,100000)
  R = zeros(Int32,10000)
  Re = zeros(Float64,14)
  return TP,PTP1,PTP2,RTP,D,R,Re
end

N=[1e6,1e7,1e8,1e9]

TP, PTP1, PTP2, RTP, D, R, Re = initialize_vectors()
times,indices_matrix=MCRAFT_SF_prerun(1e7,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R,Re)

times
indices_matrix
