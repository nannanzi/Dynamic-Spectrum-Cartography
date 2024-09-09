clc;clear all; close all;
% v = 0.01; %moving speed
% load(['Param_R8_sigma8_K64_T600_v',num2str(v,'%3.3f'),'.mat']);
warning off;
load('Param_R8_sigma8_K64_T600_v0.01.mat');
p = 0.1; %sampling raio per time slot
I = size(X4DT,1);
J = size(X4DT,2);
K  = size(X4DT,3);
T  = size(X4DT,4);
gridLen = I-1;
snr = 100;%snr
gridResolution = 1;%
x_grid = [0:gridResolution:gridLen];
y_grid = [0:gridResolution:gridLen];
[Xmesh_grid, Ymesh_grid] = meshgrid(x_grid, y_grid);
Xgrid = Xmesh_grid + 1i*Ymesh_grid;

check_timeslot = 550;

check_frequencybin = 16;

Xcube = X4DT(:,:,:,check_timeslot);
Xcubemat = tens2mat(Xcube,[],3);
singvalue = svd(Xcubemat);
%% checking number of emitters via SVD
for tau = 1:length(singvalue)
    toptauvals = sum(singvalue(1:tau));
    occ_ratio = toptauvals/sum(singvalue);
    if occ_ratio > 0.99
        Rest = tau;
        break;
    end
end




NMSE_CBTD = zeros(1,T);
NMSE_NMF = zeros(1,T);
NMSE_TPS = zeros(1,T);

rho = p;
for tt = check_timeslot
%% sampling pattern: guarantee that each column/row has at least one sample
    SampleIndextt = randperm(I*J,round(I*J*rho));
    Xcubet = squeeze(X4DT(:,:,:,tt));


    Pn = Xcubet.^2*10^(-snr/10);

    if snr>=1e2
        Pn =0;
    end
    Xnoisy = Xcubet + sqrt(Pn).*randn(I,J,K);
    Xmatt = tens2mat(Xnoisy,[],3);


    Ytt = Xmatt(SampleIndextt,:);
    Wmattt = zeros(I,J);
    Wmattt(SampleIndextt) = 1;
    Wall(:,:,tt) = Wmattt;
%% baseline 1: LL1-TPS

    Lr = 6; %hyper parameter of multi-linear rank-(Lr,Lr,1) decomposition
    [Shat_ll1_tt,Chat_ll1_tt] = CBTD(I,J,Ytt,Lr,Rest,SampleIndextt);
     % use interpolation to recover SLFs 
    Srhat_ll1_tt = TPS_based_recovery(Xgrid,SampleIndextt,Shat_ll1_tt,Rest);
    Xhat_ll1_tt = zeros(I,J,K);
    for rr=1:Rest
        Xhat_ll1_tt = Xhat_ll1_tt + outprod(Srhat_ll1_tt{rr}, Chat_ll1_tt(:,rr));
    end
    NMSE_CBTD(tt) = frob(Xhat_ll1_tt - Xcubet).^2/frob(Xcubet).^2;
%% baseline 2: NMF-TPS
    SelectInd = SPA(Ytt,Rest);
    Shat_nmf_tt = Ytt(:,SelectInd);
    Chat_nmf_tt = (Shat_nmf_tt\Ytt)';
    Srhat_nmf_tt = TPS_based_recovery(Xgrid,SampleIndextt,Shat_nmf_tt,Rest);
    Xhat_nmf_tt = zeros(I,J,K);

    for rr=1:Rest
        Xhat_nmf_tt = Xhat_nmf_tt + outprod(Srhat_nmf_tt{rr}, Chat_nmf_tt(:,rr));
    end

    NMSE_NMF(tt) = frob(Xhat_nmf_tt - Xcubet).^2/frob(Xcubet).^2;



%% baseline 3: Interpolation
    Xhat_tps_tt = zeros(I,J,K);
    dB=0;
    x_IND = real(Xgrid(SampleIndextt));
    y_IND = imag(Xgrid(SampleIndextt));
    xgrid = real(Xgrid(:));
    ygrid = imag(Xgrid(:));
    lambda_tps = 0e-4;

    for kk = 1:K
        ykk = Ytt(:,kk);
        stps = TPS(x_IND,y_IND,ykk,xgrid,ygrid,lambda_tps,dB);
        Xhat_tps_tt(:,:,kk) = reshape(stps,[I,J]);
    end

    NMSE_TPS(tt) = frob(Xhat_tps_tt - Xcubet).^2/frob(Xcubet).^2;
    
 

end




