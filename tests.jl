function initialize_vectors()
  TP = zeros(Int32,10000000)
  PTP1 = zeros(Int32,100000000)
  PTP2 = zeros(Int32,100000000)
  RTP = zeros(Int32,1000000)
  D = zeros(Int32,10000)
  R = zeros(Int32,1000)
  return TP,PTP1,PTP2,RTP,D,R
end

TP, PTP1, PTP2, RTP, D, R = initialize_vectors();

@time MCRAFT_SF(1e6,34,5,5e-3,5e-3,TP,PTP1,PTP2,RTP,D,R);

Profile.print()
