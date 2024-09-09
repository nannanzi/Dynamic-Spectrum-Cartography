clc;clear;close all;
warning off;
addpath(genpath('tensorlab_2016-03-28'))
addpath(genpath('functions'))
% addpath(genpath('G:\DZ_getdata'));
% addpath(genpath('GUI_mat'));
addpath(genpath('output_mat'));
Rcp = 15;
psd_plot = 0;
pre_switch = 1;
svalue = 1e-16;
tt_initial = [];
ICPD_flag = 0;
load('SensorInfo_all.mat');
SichuanTower = [104.094978,30.661968];
dist_threshold = 10; %10km 

dist = LA_Lo_dist(SensorInfo_all(:,2:3),SichuanTower);
Sensor_near_index = find(dist<=dist_threshold);
Sensor_near_ID = SensorInfo_all(Sensor_near_index,1);
fold_path = cd;
%save([cd,'\GUI_mat\Sensor_near_ID.mat'],"Sensor_near_ID");
dist_LaLo = dist_threshold/111;
figure(1)
longitude_min = 103.6;
longitude_max = 104.6;
latitude_min = 30.3;
latitude_max = 31;
gridResolution = 0.05;%
x_grid = [103.6:gridResolution:104.6];
y_grid = [30:gridResolution:31];
[Xmesh_grid, Ymesh_grid] = meshgrid(x_grid, y_grid);
Xgrid = Xmesh_grid + 1i*Ymesh_grid;
[I,J] = size(Xgrid);
plot(Xmesh_grid,Ymesh_grid,'b',Xmesh_grid',Ymesh_grid','b','HandleVisibility','off');
plot(real(Xgrid),imag(Xgrid),'b',real(Xgrid)',imag(Xgrid)','b','HandleVisibility','off');
xlabel('Longtitude'); %横轴label（对应含义为经度�?
ylabel('Latitude');   
title('成都市区');        %标题
axis([longitude_min longitude_max latitude_min latitude_max]); %限制坐标区显示范围，只看规定区域的图�?
hold on; %matlab图像设置，画新图不会覆盖原图
Sensor_location = SensorInfo_all(:,[2:3]);
scatter(Sensor_location(:,1),Sensor_location(:,2),50,...
                'MarkerFaceColor','g',...
                'MarkerEdgeColor','k',...
              'LineWidth',1.5,'DisplayName','传感器位�?');

hold on;
scatter(SichuanTower(1),SichuanTower(2),40,'LineWidth',1.5,'MarkerFaceColor','Y','MarkerEdgeColor','flat','DisplayName','成都市广播电视塔');
legend;
%
Sensor_ID_block = [];%[202102044,202102126];
loc_ob_X = [103.95,104.15,104.15,103.95,103.95];
loc_ob_Y = [30.7,30.7,30.6,30.6,30.7];
% loc_ob_X = [104.1,104.15,104.15,104.1,104.1];
% loc_ob_Y = [30.7,30.7,30.6,30.6,30.7];
loc_ob_la = [max(loc_ob_Y),min(loc_ob_Y)];
loc_ob_lo = [max(loc_ob_X),min(loc_ob_X)];
plot(loc_ob_X,loc_ob_Y,'R','LineWidth',1.25,'HandleVisibility','off'); % 观察范围
sensor_ob_ID_idx = find(SensorInfo_all(:,2)>=loc_ob_lo(2)&SensorInfo_all(:,2)<=loc_ob_lo(1)&SensorInfo_all(:,3)>=loc_ob_la(2)&SensorInfo_all(:,3)<=loc_ob_la(1));
sensor_ob_ID = SensorInfo_all(sensor_ob_ID_idx,1);
sensor_ob_ID = setdiff(sensor_ob_ID,Sensor_ID_block);

% read real data
load('Pall_1_88_108_2023_03_07_1700_1900.mat');
      %时间分辨�?
            %频段起始
                    %  时间    ，时间段
T = length(P_save);
T_test = T; %T

% 传感器�?�择

Mt_max_test = intersect(P_save{1}.AvailableSensorID,sensor_ob_ID);
for tt = 1:T_test
    Mt_ID_temp1 = intersect(P_save{tt}.AvailableSensorID,sensor_ob_ID);
    Mt_ID_temp1 = setdiff(Mt_ID_temp1,Sensor_ID_block);
    %Mt_tt_ID{tt} = Mt_ID_temp1;
    Mt(tt) = length(Mt_ID_temp1);
    Mt_ID_temp2 = intersect(P_save{min(tt+1,T_test)}.AvailableSensorID,sensor_ob_ID);
    Mt_test_ID{tt} = intersect(Mt_ID_temp1,Mt_ID_temp2);

    Mt_max_test = union(Mt_max_test,Mt_ID_temp2);
    M(tt) = length(Mt_test_ID{tt});
end
[M_max,M_sensoridx] = max(Mt);
Sensor_ID_time = intersect(P_save{M_sensoridx}.AvailableSensorID,sensor_ob_ID); %�?大工�?

M_init = length(Sensor_ID_time);

Rt = zeros(T,1);
% 选择的传感器绘图
[~,plot_idx] = intersect(SensorInfo_all(:,1),Sensor_ID_time);
Sensor_init = SensorInfo_all(plot_idx,[2:3]);
scatter(Sensor_init(:,1),Sensor_init(:,2),50,...
                'MarkerFaceColor','b',...
                'MarkerEdgeColor','r',...
              'LineWidth',1.5,'DisplayName','选择传感器位�?');
%% 频段
Freband_ob = [600:1:640];
Fre_true = [609,629]-[Freband_ob(1)]+1;
%Freband_true = [Fre_true(1)-3:Fre_true(1)+3;Fre_true(2)-3:Fre_true(2)+3];
%Fre_true = [585]-[Freband_ob(1)]+1;
Freband_true = [Fre_true(1)-3:Fre_true(1)+3;Fre_true(2)-3:Fre_true(2)+3];
%Freband_true = [Fre_true(1)-3:Fre_true(1)+3];
K = length(Freband_ob);

%% 频偏预处�? 
for tt = 1:T_test
    [Sensor_working_tt_ID,Sensor_data_tt_idx] = intersect(P_save{tt}.AvailableSensorID,Sensor_ID_time);
    [~,Sensor_working_tt_idx] = intersect(Sensor_ID_time,Sensor_working_tt_ID);
    Mt_tt_ID{tt} = Sensor_working_tt_ID;
    %radio map
    Xomega = P_save{tt}.data(Sensor_data_tt_idx,Freband_ob);
%     if psd_plot == 1
%         plot_psd_function(Xomega',0);
%     end
    [output,test1{tt},test2{tt}] = Xdata_Preprocess(Xomega,Freband_true);
    if mod(tt,10)==0 && psd_plot == 1
        plot_psd_function(output',[],0);
    end
    output = col_norm(output')';
    if pre_switch == 1
        Xomega_pre{tt} = output;
    else
        Xomega_pre{tt} = Xomega;
    end
%     if tt >1
%         [Sensor_test_ID, Sensor_test_idx] = intersect(P_save{tt}.AvailableSensorID,Mt_test_ID{tt-1});
%         Xtest = P_save{tt}.data(Sensor_test_idx,Freband_ob);
%         Xtest_temp = Xdata_Preprocess(Xtest,Freband_true);
%         Xomega_test{tt} = col_norm(Xtest_temp')'; 
%     end
end

%%
rho = 0.15;
dW = 10;
if dW < 10 %dW should at least be 10 to fit the order of ARMA
    dW = 10;
end
iter_inner = 1;
%%
Rrank = Rcp;
lambda = 0.09;

for tt = 1:T_test
    %% Samling
    [Sensor_working_tt_ID,Sensor_data_tt_idx] = intersect(P_save{tt}.AvailableSensorID,Sensor_ID_time);
    [~,Sensor_working_tt_idx] = intersect(Sensor_ID_time,Sensor_working_tt_ID);
    Sensor_sampling_idx{tt} = Sensor_working_tt_idx;
    Wtt = zeros(M_init,K); % Mask
    rng(9,'Twister');
    sampIndextt = randperm(M_init,round(M_init*rho));

    samp_test(tt,:) = sampIndextt;
    Wtt(sampIndextt,:) = 1;
    Wacc(:,:,tt) = Wtt; % Mask tensor
    % radio map
    Xtrue_tt_mat = zeros(M_init,K);
    Xtrue_tt_mat(Sensor_working_tt_idx,:) = Xomega_pre{tt};
    Xacc(:,:,tt) = Xtrue_tt_mat;
    Yacc = Wacc.* Xacc; % Mask tensor
    Ytt = Wtt.*Xtrue_tt_mat;
    %% NMF
    %Xomega = Xomega_pre{tt};
    Xomega = Ytt(sampIndextt,:);
    if isempty(Xomega)
        for rr = 1:max(Rt)
            cacheA{rr,tt} = [];
            cacheB{rr,tt} = [];
            %cacheC{r;
        end
        continue;
    end
    Rt(tt) = 2;%estRnum(Xomega,0.9);
    if tt >= 1       
        
        SelectInd = SPA(Xomega,Rt(tt));
        Sest{tt} = Xomega(:,SelectInd);
        Cest_temp = (Sest{tt}\Xomega)';
        Cest_temp(Cest_temp<svalue) = svalue;
        Cest{tt} = Cest_temp;




%         if psd_plot == 1
%             plot_psd_function(Cest{tt},1);
%         end
        if tt>1
            P_mat = cpderr_order(Cest{1},Cest{tt});
            Permutation_mat{tt} = P_mat;
            Cest{tt} = Cest{tt}*Permutation_mat{tt};
            Sest{tt} = Sest{tt}*Permutation_mat{tt};
        end
    end
    %% 置换不确定�??
    tt_initial = [tt_initial,tt];

    % 辐射源X_r排序
    Xnmf_cpd_ob_tt =[];
    if Rt(tt)~=0
        for rr = 1:Rt(tt)
            sctt = reshape(Sest{tt}(:,rr),[],1); % 构成X
            Xnmf_cpd_ob_tt{rr}(:,:) = outprod(sctt,Cest{tt}(:,rr));
        end
    else
        for rr = 1:Rt(tt)
            sctt = reshape(Sest{tt}(:,rr),[],1); % 构成X
            Xnmf_cpd_ob_tt{rr}(:,:) = Xomega_pre{tt};
        end
    end
    Xnmf_cpd_temp = zeros(M_init,K);
    for rr = 1:Rt(tt)
        for dd = 1:length(sampIndextt)
            Xnmf_cpd_temp(sampIndextt(dd),:) = Xnmf_cpd_ob_tt{rr}(dd,:);
        end
        Xnmf_cpd{tt,rr} = Xnmf_cpd_temp;
    end
    %Xnmf_cpd_T(:,:,tt,:) = Xnmf_cpd; 
    
    %% Factors & ARMA initializing via DWCPD.
    if length(tt_initial) == dW
        for rr = 1:Rt(tt)
            %Y = squeeze(Xnmf_cpd_T(:,:,1:tt,rr));
            Y_cell = Xnmf_cpd(tt_initial,rr);
            Y = mat2tens(cell2mat(Y_cell),[M_init,K,dW],[],2);
            Wacc_init = Wacc(:,:,tt_initial);
            Yacc = Wacc_init.*Y;
            AdW{rr} = randn(M_init,Rrank);
            BdW{rr} = randn(K,Rrank);
            [AdW{rr},BdW{rr},CdW{rr}] = DWCPD(Yacc,Wacc_init,Rrank,dW,AdW{rr},BdW{rr}); % Ycc、Wacc
            
            %% Factors cache for ICPD
            cacheA{rr,tt} = AdW{rr};
            cacheB{rr,tt} = BdW{rr};
            cacheC{rr} = CdW{rr};
            Yt_rr = squeeze(Xnmf_cpd{tt,rr});
            for ii = 1:M_init
                Sit{rr}(ii,:) = Yt_rr(ii,:)*diag(Wtt(ii,:))*BdW{rr}*diag(CdW{rr}(end,:));
                Rit{rr}{ii} = diag(CdW{rr}(end,:))*BdW{rr}'*diag(Wtt(ii,:))*BdW{rr}*diag(CdW{rr}(end,:));
            end
            
            for jj = 1:K
                Qjt{rr}(jj,:) = (diag(CdW{rr}(end,:))*AdW{rr}'*diag(Wtt(:,jj))*Yt_rr(:,jj) )';
                Pjt{rr}{jj} = diag(CdW{rr}(end,:))*AdW{rr}'*diag(Wtt(:,jj))*AdW{rr}*diag(CdW{rr}(end,:));
            end
        end
        ICPD_flag = 1;
        continue;
    end
    
    if ICPD_flag
        %% ICPD Updating
        for rr =1:Rt(tt)
            %Yt_rr = squeeze(Xnmf_cpd_T(:,:,tt,rr));
            tt_back = 1;
            if tt == 22
                test=1;
            end
            while(1)
                Yt_rr = Xnmf_cpd{tt,rr};
                Ait{rr} = cacheA{rr,tt-tt_back};
                Bit{rr} = cacheB{rr,tt-tt_back};
                c_incre{rr} = cacheC{rr}(tt_back,:);
                %c_incre{rr} = C_incre{rr}(end,:)';
                if ~isempty(Ait{rr})
                    break;
                end
                tt_back = tt_back+1;
            end
            for iter = 1:iter_inner
               [Ait{rr},Bit{rr},Rit{rr},Sit{rr},Pjt{rr},Qjt{rr},cals{rr}] = UpdateAll_Incre_Re(Ait{rr},Bit{rr},Rit{rr},Sit{rr},Pjt{rr},Qjt{rr},Yt_rr,Wtt,lambda);
            end
            cacheA{rr,tt} = Ait{rr};
            cacheB{rr,tt} = Bit{rr};
            cacheC{rr}(tt,:) = cals{rr};
        end
        % NMSE of completion
        X_comp = zeros(M_init,K);
        for rr = 1:Rt(tt)
            Anmse_comp = cacheA{rr,end};
            Bnmse_comp = cacheB{rr,end};
            c_nmse_comp = cacheC{rr}(end,:);
            S_comp = Anmse_comp*diag(c_nmse_comp)*Bnmse_comp';   
            X_comp = X_comp + S_comp;
        end
        nmsett_comp = frob(X_comp - Xtrue_tt_mat)^2/frob(Xtrue_tt_mat)^2;
        nmsett_comp
        NMSEcomp(tt) = nmsett_comp; 
    end
    
end
figure;
plot(NMSEcomp,'LineWidth',1.25,'DisplayName',['ICPD-ARMA completion']);
grid on;
set(gca,'GridLineStyle',':');
xlabel('t');
ylabel('NMSE');
legend('Interpreter','latex');
set(gca,'FontName','times new roman');
xlim([1,length(NMSEcomp)]);
set(gca,'YScale','log');
ylim([1e-3 5e0]);
mean(NMSEcomp(20:end))

%% data 
save(['..\Xtrue_rho',num2str(rho),'.mat'],"Xacc"); % (M*K*T) truth
save(['..\Y_rho',num2str(rho),'.mat'],"Yacc"); % (M*K*T)
save(['..\mask_rho',num2str(rho),'.mat'],"Wacc"); % (M*K*T)

%% save NMSEmat
save(['..\..\result\NMSE_IncreDSC_ALS_rho',num2str(rho),'.mat'],"NMSEcomp");

