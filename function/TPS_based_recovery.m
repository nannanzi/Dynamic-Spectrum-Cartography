% TPS interpolation to recover SLF
function Sest = TPS_based_recovery(Xgrid,sampleindex,SLFmat,Rest)
    I = size(Xgrid,1);
    J = size(Xgrid,2);
    dB=0;
    x_IND = real(Xgrid(sampleindex));
    y_IND = imag(Xgrid(sampleindex));
    xgrid = real(Xgrid(:));
    ygrid = imag(Xgrid(:));
    lambda_tps = 0e-4;
    for rr=1:Rest
        stps{rr} = TPS( x_IND,y_IND,SLFmat(:,rr),xgrid,ygrid,lambda_tps,dB );
        Sest{rr} = reshape(stps{rr},I,J);
    end
end