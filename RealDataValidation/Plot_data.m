T = size(NMSE_BatchDSC,2);
plot(1:T,10*log10(NMSE_IncreDSC(1,:)),1:T,10*log10(NMSE_BatchDSC(1,:)));