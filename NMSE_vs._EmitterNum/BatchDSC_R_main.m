clc;clear all; close all;
warning off;
addpath(genpath('function'));
addpath(genpath('spectrum_data'));
svalue = 1e-16;
I = 51;
J = 51;
K = 64;
T = 600;
gridLen = I-1;
gridResolution = 1;
snr = 100;%snr of additive noise
x_grid = [0:gridResolution:gridLen];
y_grid = [0:gridResolution:gridLen];
[Xmesh_grid, Ymesh_grid] = meshgrid(x_grid, y_grid);
Xgrid = Xmesh_grid + 1i*Ymesh_grid;
p = 0.1; %sampling raio per time slot
rho = p;
Fr = ceil(J*rho); %hyper-parameter: CP-rank
InitActive = round(1.5/rho); %hyper_parameter: initialization batchsize
lambda = 0.9; %hyper-parameter: forgetting factor
% maxiter = 10; %hyper-parameter: iteration number
mu = 1e-8; % hyper-parameter: regularization

Batchsize = 20;
Batchactivate = Batchsize:1:T;


Rall = 2:2:6;
NMSE_BatchDSC = zeros(length(Rall),T);
%initializing cache variables for DWCPD
for rnum = 1:length(Rall)
%initializing cache variables for OptDSC
    R = Rall(rnum);
    load(['Param_R',num2str(R,'%d'),'_sigma8_K64_T600_v0.01.mat']);
    check_timeslot = randperm(T,1);
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
    Rallest{rnum} = Rest;

    for rr = 1:Rest
        A{rr} = randn(I,Fr);
        B{rr} = randn(J,Fr);
        C{rr} = randn(Batchsize,Fr);
    end


   

% Wall = sampling_pattern_exclude(I,J,p);
% cyc_len = length(Wall);
% SampleIndexall = sampling_pattern_retain(I,J,p,0);
% cyc_len = length(SampleIndexall);
    for tt = 1:T
    %% sampling pattern: guarantee that each column/row has at least one sample
    %     cyc_idx = mod(tt,cyc_len) + cyc_len*(mod(tt,cyc_len) == 0);
    %     Wmatt = Wall{cyc_idx};
    %     Wmatt = Wall{tt};
        SampleIndex = randperm(I*J,round(I*J*rho));
        Wmatt = zeros(I,J);
    %     SampleIndex = SampleIndexall{cyc_idx};
        Wmatt(SampleIndex) = 1;
        Wall(:,:,tt) = Wmatt;
        Wvect = Wmatt(:);
        SampleIndextt = find(Wvect);
        Xcubet = squeeze(X4DT(:,:,:,tt));

        Pn = Xcubet.^2*10^(-snr/10);

        if snr>=1e2
            Pn =0;
        end
        Xnoisy = Xcubet + sqrt(Pn).*randn(I,J,K);

        Xmatt = tens2mat(Xnoisy,[],3);
        Ytt = Xmatt(SampleIndextt,:);
    %% Proposed BatchDSC
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
            Cinit_opt{tt} = Cinit_opt_tt;
            Ginit_opt_tt = Ytt/(Cinit_opt_tt');
        end

        Xhatt = zeros(I,J,K);
        for rr = 1:Rest
            Gmatr = zeros(I,J);
            Gmatr(SampleIndextt) = Ginit_opt_tt(:,rr);
            Gtens{rr}(:,:,tt) = Gmatr;
            if ismember(tt,Batchactivate)
                [A{rr},B{rr},C{rr}] = LowCPrankTC(Gtens{rr},Wall,Fr,10,lambda,Batchsize,A{rr},B{rr});
                Sestr = A{rr}*diag(C{rr}(end,:))*(B{rr})';
                Xhatt =  Xhatt + outprod(Sestr,Cinit_opt_tt(:,rr));
            end
        end
        NMSE_BatchDSC(rnum,tt) = frob(Xhatt - Xcubet).^2/frob(Xcubet).^2;
    end
end
    


