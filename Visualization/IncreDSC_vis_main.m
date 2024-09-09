clc;clear all; close all;
warning off;
addpath(genpath('function'));
addpath(genpath('spectrum_data'));
load('Param_R8_sigma8_K64_T600_v0.1.mat');
svalue = 1e-16;
I = size(X4DT,1);
J = size(X4DT,2);
K  = size(X4DT,3);
T  = size(X4DT,4);
gridLen = I-1;
gridResolution = 1;
snr = 100;%snr of additive noise
x_grid = [0:gridResolution:gridLen];
y_grid = [0:gridResolution:gridLen];
[Xmesh_grid, Ymesh_grid] = meshgrid(x_grid, y_grid);
Xgrid = Xmesh_grid + 1i*Ymesh_grid;
p = 0.1; %sampling raio per time slot
Fr = 5; %hyper-parameter: CP-rank
% Fr = 15;
lambda = 0.9; %hyper-parameter: forgetting factor
% maxiter = 10; %hyper-parameter: iteration number
mu = 1e-8; % hyper-parameter: regularization
check_timeslot = 550;
check_frequencybin = 16;

Xcube = X4DT(:,:,:,check_timeslot);
Xcubemat = tens2mat(Xcube,[],3);
singvalue = svd(Xcubemat);
%% checking number of emitters via SVD
for tau = 1:length(singvalue)
    toptauvals = sum(singvalue(1:tau));
    occ_ratio = toptauvals/sum(singvalue);
    if occ_ratio > 0.999
        Rest = tau;
        break;
    end
end


NMSE_IncreDSC = zeros(1,T);

%initializing cache variables for IncreDSC

for rr = 1:Rest
    for ii = 1:I
        Rcache{rr}{ii} = 0;
        scache{rr}{ii} = 0;
    end
    for jj = 1:J
        Pcache{rr}{jj} = 0;
        qcache{rr}{jj} = 0;
    end
        A{rr} = rand(I,Fr);
        B{rr} = rand(J,Fr);
        C{rr} = [];
end




    % Wall = sampling_pattern_exclude(I,J,p);
    % cyc_len = length(Wall);
rho = p;
SampleIndexall = sampling_pattern_retain(I,J,rho,0);
cyc_len = length(SampleIndexall);
for tt = 1:T
%% sampling pattern: guarantee that each column/row has at least one sample
    cyc_idx = mod(tt,cyc_len) + cyc_len*(mod(tt,cyc_len) == 0);
    Wmatt = zeros(I,J);
    SampleIndex = SampleIndexall{cyc_idx};
    Wmatt(SampleIndex) = 1;

    Wvect = Wmatt(:);
    SampleIndextt = find(Wvect);%re-check
    Xcubet = squeeze(X4DT(:,:,:,tt));

    Pn = Xcubet.^2*10^(-snr/10);

    if snr>=1e2
        Pn =0;
    end
    Xnoisy = Xcubet + sqrt(Pn).*randn(I,J,K);

    Xmatt = tens2mat(Xnoisy,[],3);
    Ytt = Xmatt(SampleIndextt,:);
%% Proposed IncreDSC
 % Initializing via NMF
    SelectInd = SPA(Ytt,Rest);
    Ginit_opt_tt = Ytt(:,SelectInd);
    %column-wise normalization 
    Ginit_opt_tt = Ginit_opt_tt / diag(max(Ginit_opt_tt));
    Cinit_opt_tt = (Ginit_opt_tt\Ytt)';
    [Ginit_opt_tt, Cinit_opt_tt] = NMF_HALS(Ytt,Rest,50,Ginit_opt_tt,Cinit_opt_tt);
    Cinit_opt{tt} = Cinit_opt_tt;
    if tt > 1 %unifying ambiguities
        Cinit_opt_historical = Cinit_opt{tt-1};
        [err,~,~,Cinit_opt_tt] = cpderr(Cinit_opt_historical,Cinit_opt_tt);
        Cinit_opt{tt} =  Cinit_opt_tt;
        Ginit_opt_tt = Ytt/(Cinit_opt{tt}');
    end

    Xhatt = zeros(I,J,K);
    for rr = 1:Rest
        Gmat{rr} = zeros(I,J);
        Gmat{rr}(SampleIndextt) = Ginit_opt_tt(:,rr);
        Rcacherr = Rcache{rr};
        scacherr = scache{rr};
        Pcacherr = Pcache{rr};
        qcacherr = qcache{rr};
        if tt == 1
            [A{rr},B{rr},Rcacherr,scacherr,Pcacherr,qcacherr,c{rr}] = UpdateAll_ALS(A{rr},B{rr},Rcacherr,scacherr,Pcacherr,qcacherr,Gmat{rr},Wmatt,lambda);
            c{rr} = c{rr}';
        else
            [A{rr},B{rr},Rcacherr,scacherr,Pcacherr,qcacherr,c{rr}] = UpdateAll_Grad(A{rr},B{rr},c{rr},Rcacherr,scacherr,Pcacherr,qcacherr,Gmat{rr},Wmatt,lambda);
        end
        Rcache{rr} = Rcacherr;
        scache{rr} = scacherr;
        Pcache{rr} = Pcacherr;
        qcache{rr} = qcacherr;

        Srt = A{rr} * diag(c{rr}) * B{rr}';
        Sr{tt} = Srt;
        Xhatt = Xhatt + outprod(Srt, Cinit_opt_tt(:,rr));
    end

    NMSE_IncreDSC(tt) = frob(Xhatt - Xcubet).^2/frob(Xcubet).^2;
    
    if tt == check_timeslot
        figure(1)
        valmin = min(min(min(Xcubet(:,:,check_frequencybin))));
        valmax = max(max(max(Xcubet(:,:,check_frequencybin))));
        Xhatt(Xhatt<valmin) = valmin;
        Xhatt(Xhatt>valmax) = valmax;
        Xslice = squeeze(Xhatt(:,:,check_frequencybin));
        contourf(10*log10(Xslice),100,'linecolor','None');
        colormap jet;
        set(gca,'xtick',[],'xticklabel',[])
        set(gca,'ytick',[],'yticklabel',[])
        title(['IncreDSC,NMSE=', num2str(NMSE_IncreDSC(tt),'%3.3f')])
        set(gca,'FontName','Times New Roman','FontSize',15,'LineWid',1);
        axes('position',[0.2,0.02,.6,.3])
        axis off
        my_handle = colorbar('east');
        my_handle.Title.String='dB';

        figure(2)
        Xvis = squeeze(Xcubet(:,:,check_frequencybin));
        contourf(10*log10(Xvis),100,'linecolor','None');
        colormap jet;
        set(gca,'xtick',[],'xticklabel',[])
        set(gca,'ytick',[],'yticklabel',[])
        title('Ground-truth')
        set(gca,'FontName','Times New Roman','FontSize',15,'LineWid',1);
        axes('position',[0.2,0.02,.6,.3])
        axis off
        my_handle = colorbar('east');
        my_handle.Title.String='dB';
    end


end
    


