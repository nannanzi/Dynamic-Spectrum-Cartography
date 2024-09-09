clc;clear all; close all;
warning off;
addpath(genpath('function'));
addpath(genpath('spectrum_data'));
load('Param_R8_sigma8_K64_T600_v0.01.mat');
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
p = 0.2; %sampling raio per time slot
Batchsize = 20;
Batchactivate = Batchsize:1:T;


lambda = 0.9; %hyper-parameter: forgetting factor
% maxiter = 10; %hyper-parameter: iteration number
mu = 1e-8; % hyper-parameter: regularization
check_timeslot = randperm(T,1);
check_frequencybin = 16;


Restall = 5:10;
%initializing cache variables for DWCPD
NMSE_BatchDSC = zeros(length(Restall),T);
for rnum = 1:length(Restall)
%initializing cache variables for OptDSC
    rho = p;
    Rest = Restall(rnum);
    Fr = ceil(J*rho); %hyper-parameter: CP-rank
    InitActive = round(1.5/rho); %hyper_parameter: initialization batchsize
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
    %     figure(1)
    %     valmin = min(min(min(Xcubet)));
    %     Xhatt(Xhatt<valmin) = valmin;
    %     Xslice = squeeze(Xhatt(:,:,check_frequencybin));
    %     contourf(10*log10(Xslice),100,'linecolor','None');
    %     colormap jet;
    %     set(gca,'xtick',[],'xticklabel',[])
    %     set(gca,'ytick',[],'yticklabel',[])
    %     title('Inst_{Recov}')
    %     set(gca,'FontName','Times New Roman','FontSize',15,'LineWid',1);
    %     axes('position',[0.2,0.02,.6,.3])
    %     axis off
    %     my_handle = colorbar('east');
    %     my_handle.Title.String='dB';
    %     
    %     figure(2)
    %     Xvis = squeeze(Xcubet(:,:,check_frequencybin));
    %     contourf(10*log10(Xvis),100,'linecolor','None');
    %     colormap jet;
    %     set(gca,'xtick',[],'xticklabel',[])
    %     set(gca,'ytick',[],'yticklabel',[])
    %     title('Ground-truth')
    %     set(gca,'FontName','Times New Roman','FontSize',15,'LineWid',1);
    %     axes('position',[0.2,0.02,.6,.3])
    %     axis off


    %     my_handle.Title.String='dB';
    %         

    end
end
    


