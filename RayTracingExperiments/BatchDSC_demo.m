clc;clear all; close all;
warning off;
addpath(genpath('function'));
load('RayTracingMaps.mat');
X4DT = Xtrue_tens_all;
svalue = 1e-16;
[I,J,K,T] = size(X4DT);
gridLen = I-1;
gridResolution = 1;
x_grid = [0:gridResolution:gridLen];
y_grid = [0:gridResolution:gridLen];
[Xmesh_grid, Ymesh_grid] = meshgrid(x_grid, y_grid);
Xgrid = Xmesh_grid + 1i*Ymesh_grid;
%p = 0.1:0.01:0.2; %sampling raio per time slot
rho = 0.1; %sampling raio per time slot
Fr = 5; %hyper-parameter: CP-rank


lambda = 1; %hyper-parameter: forgetting factor
% maxiter = 10; %hyper-parameter: iteration number
mu = 1e-8; % hyper-parameter: regularization
check_timeslot = randperm(T,1);
% check_frequencybin = 16;

vis_bin = 7;
vis_time = 190;

% Batchsize = 20;
% Batchactivate = Batchsize:1:T;

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


NMSE_BatchDSC = zeros(1,T);


%     Batchsize = round(2/rho);
%     Batchactivate = Batchsize:1:T;




   

% Wall = sampling_pattern_exclude(I,J,p);
% cyc_len = length(Wall);
 SampleIndexall = sampling_pattern_retain(I,J,rho,0.2);
 cyc_len = length(SampleIndexall);
 Batchsize = cyc_len;
 Batchactivate = Batchsize:1:T;

 for rr = 1:Rest
    A{rr} = randn(I,Fr);
    B{rr} = randn(J,Fr);
    C{rr} = randn(Batchsize,Fr);
 end


for tt = 1:T
%% sampling pattern: guarantee that each column/row has at least one sample
    cyc_idx = mod(tt,cyc_len) + cyc_len*(mod(tt,cyc_len) == 0);
    %     Wmatt = Wall{cyc_idx};
    %     Wmatt = Wall{tt};
    Wmatt = zeros(I,J);
    SampleIndex = SampleIndexall{cyc_idx};
    
    Wmatt(SampleIndex) = 1;
    Wall(:,:,tt) = Wmatt;
    
    Wvect = Wmatt(:);
    SampleIndextt = find(Wvect);
    Xcubet = squeeze(X4DT(:,:,:,tt));
    
    Xmatt = tens2mat(Xcubet,[],3);
    Ytt = Xmatt(SampleIndextt,:);
    %% Proposed BatchDSC
    % Initializing via NMF
    SelectInd = SPA(Ytt,Rest);
    Ginit_opt_tt = Ytt(:,SelectInd);
    %column-wise normalization 
    Ginit_opt_tt = Ginit_opt_tt / diag(max(Ginit_opt_tt));
    Cinit_opt_tt = (Ginit_opt_tt\Ytt)';
    %        [Ginit_opt_tt, Cinit_opt_tt] = NMF_HALS(Ytt,Rest,50,Ginit_opt_tt,Cinit_opt_tt);
    Cinit_opt{tt} = Cinit_opt_tt;
    if tt > 1 %unifying ambiguities
        Cinit_opt_historical = Cinit_opt{tt-1};
        [err,~,~,Cinit_opt_tt] = cpderr(Cinit_opt_historical,Cinit_opt_tt);
        Cinit_opt{tt} = Cinit_opt_tt;
    %             Cinit_opt{tt} = Cinit_opt_historical;   
        Ginit_opt_tt = Ytt/(Cinit_opt{tt}');
    end
    
    Xhatt = zeros(I,J,K);
    for rr = 1:Rest
        Gmatr = zeros(I,J);
        Gmatr(SampleIndextt) = Ginit_opt_tt(:,rr);
        Gtens{rr}(:,:,tt) = Gmatr;
        if ismember(tt,Batchactivate)
            [A{rr},B{rr},C{rr}] = BatchCPD(Gtens{rr},Wall,Fr,10,lambda,Batchsize,A{rr},B{rr},C{rr});
            Sestr = A{rr}*diag(C{rr}(end,:))*(B{rr})';
            Xhatt =  Xhatt + outprod(Sestr,Cinit_opt_tt(:,rr));
        end
    end
    
    NMSE_BatchDSC(tt) = frob(Xhatt - Xcubet).^2/frob(Xcubet).^2;

    if tt == vis_time

        
        figure(1)
        valmin = min(min(Xcubet(:,:,vis_bin)));
        valmax = max(max(Xcubet(:,:,vis_bin)));
 
        Xslice = squeeze(Xhatt(:,:,vis_bin));

        Xslice(Xslice<valmin) = valmin;
        Xslice(Xslice>valmax) = valmax;

        contourf(10*log10(Xslice),100,'linecolor','None');
        colormap jet;
        set(gca,'xtick',[],'xticklabel',[])
        set(gca,'ytick',[],'yticklabel',[])
        title(['BatchDSC,NMSE=', num2str(10*log10(NMSE_BatchDSC(tt)),'%3.3f')])
        set(gca,'FontName','Times New Roman','FontSize',15,'LineWid',1);
        axes('position',[0.2,0.02,.6,.3])
        axis off
        my_handle = colorbar('east');
        my_handle.Title.String='dB';

    end

    
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

    


