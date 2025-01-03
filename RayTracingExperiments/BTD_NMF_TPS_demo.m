clc;clear all; close all;
warning off;
addpath(genpath('function'));
load('RayTracingMaps.mat');
X4DT = Xtrue_tens_all;
svalue = 1e-16;
[I,J,K,T] = size(X4DT);
gridLen = I-1;
% p = 0.1:0.01:0.2; %sampling raio per time slot
rho = 0.1;
snr = 100;%snr
gridResolution = 1;%
x_grid = [0:gridResolution:gridLen];
y_grid = [0:gridResolution:gridLen];
[Xmesh_grid, Ymesh_grid] = meshgrid(x_grid, y_grid);
Xgrid = Xmesh_grid + 1i*Ymesh_grid;

check_timeslot = randperm(T,1);

vis_bin = 7;
vis_time = 190;

Xcube = X4DT(:,:,:,check_timeslot);
Xcubemat = tens2mat(Xcube,[],3);
singvalue = svd(Xcubemat);
%% checking number of emitters via SVD
for tau = 1:length(singvalue)
    toptauvals = sum(singvalue(1:tau));
    occ_ratio = toptauvals/sum(singvalue);
    if occ_ratio > 0.85
        Rest = tau;
        break;
    end
end

%initializing cache variables for DWCPD

%initializing cache variables for OptDSC
for rr = 1:Rest
    for ii = 1:I
        Rcache{rr}{ii} = 0;
        scache{rr}{ii} = 0;
    end
    for jj = 1:J
        Pcache{rr}{jj} = 0;
        qcache{rr}{jj} = 0;
    end
end


NMSE_CBTD = zeros(1,T);
NMSE_NMF = zeros(1,T);
NMSE_TPS = zeros(1,T);

for tt = 1:T
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


    if tt == vis_time
        valmin = min(min(Xcubet(:,:,vis_bin)));
        valmax = max(max(Xcubet(:,:,vis_bin)));

        figure(1)
       
        Xslice = squeeze(Xhat_ll1_tt(:,:,vis_bin));

        Xslice(Xslice<valmin) = valmin;
        Xslice(Xslice>valmax) = valmax;


        contourf(10*log10(Xslice),100,'linecolor','None');
        colormap jet;
        set(gca,'xtick',[],'xticklabel',[])
        set(gca,'ytick',[],'yticklabel',[])
        title(['BTD-TPS,NMSE=', num2str(10*log10(NMSE_CBTD(tt)),'%3.3f')])
        set(gca,'FontName','Times New Roman','FontSize',15,'LineWid',1);
        axes('position',[0.2,0.02,.6,.3])
        axis off
        my_handle = colorbar('east');
        my_handle.Title.String='dB';

        figure(2)
        
        Xslice = squeeze(Xhat_nmf_tt(:,:,vis_bin));
        Xslice(Xslice<valmin) = valmin;
        Xslice(Xslice>valmax) = valmax;


       
        contourf(10*log10(Xslice),100,'linecolor','None');
        colormap jet;
        set(gca,'xtick',[],'xticklabel',[])
        set(gca,'ytick',[],'yticklabel',[])
        title(['NMF-TPS,NMSE=', num2str(10*log10(NMSE_NMF(tt)),'%3.3f')])
        set(gca,'FontName','Times New Roman','FontSize',15,'LineWid',1);
        axes('position',[0.2,0.02,.6,.3])
        axis off
        my_handle = colorbar('east');
        my_handle.Title.String='dB';

    
    
        my_handle.Title.String='dB';

        figure(3)
        
        Xslice = squeeze(Xhat_tps_tt(:,:,vis_bin));
        Xslice(Xslice<valmin) = valmin;
        Xslice(Xslice>valmax) = valmax;


       
        contourf(10*log10(Xslice),100,'linecolor','None');
        colormap jet;
        set(gca,'xtick',[],'xticklabel',[])
        set(gca,'ytick',[],'yticklabel',[])
        title(['TPS,NMSE=', num2str(10*log10(NMSE_TPS(tt)),'%3.3f')])
        set(gca,'FontName','Times New Roman','FontSize',15,'LineWid',1);
        axes('position',[0.2,0.02,.6,.3])
        axis off
        my_handle = colorbar('east');
        my_handle.Title.String='dB';

    
    
        my_handle.Title.String='dB';


    end

end



