function X4DT = RMGeneratorFun(R,v,sigma_s,directional)
    addpath(genpath('function'));
    loss_f = @(x,d,alpha) min(1,(x/d).^(-alpha)); %synthetic pathloss component
    alpha = 2;%path loss coefficient
    d0 = 2;

    K = 2^6; %number of frequency bins
    %% Generating spatial meshgrids
    T = 600;
    gridLen = 50;
    gridResolution = 1;%
    x_grid = [0:gridResolution:gridLen];
    y_grid = [0:gridResolution:gridLen];
    [Xmesh_grid, Ymesh_grid] = meshgrid(x_grid, y_grid);
    Xgrid = Xmesh_grid + 1i*Ymesh_grid;
    % X1 = transpose(X1);
    [I,J] = size(Xgrid); 


    %% Generating PSDs
    indK = [1:K]';
    Sx =@(f0,a) sinc((indK-f0)/a).^2.*( abs((indK-f0)/a)<=1);  % basis sinc wave function
    Ctrue = zeros(K,R);
    ind_psd = 2:2:K-2;

    for rr = 1:R
        am = 1.5 + 0.5*rand(3,1); % random amplitude
        centerf0 = ind_psd(randperm(length(ind_psd),3));
        cr = am(1)*Sx(centerf0(1),3*(1+rand)) + (rand>0.5)*am(2)*Sx(centerf0(2),3*(1+rand)) + (rand>0.5)*am(3)*Sx(centerf0(3),3*(1+rand));
    %     cr = ones(K,1); 
        Ctrue(:,rr) = cr;
    end




    %% Generating emitters' traces
    for rr = 1:R
        start_location_rr = rand*(I-1) + 1j*rand*(J-1);
        clock_wise_rr = (rand < .5);
        rect_length_rr = rand*(I - real(start_location_rr));
        rect_width_rr = rand*(J - imag(start_location_rr));
        for tt = 1:T
            [xt,yt] = RectCircle(v*tt,rect_length_rr,rect_width_rr,clock_wise_rr);
            xi{rr}(tt) = real(start_location_rr) + xt;
            yi{rr}(tt) = imag(start_location_rr) + yt;
        end
        if directional %generating radius angle
            Ang1{rr} = 0.5*rand*pi;
            Ang2{rr} = Ang1{rr} + 0.5*pi + 0.5*rand*pi;
            Ang2{rr} = max{Ang2{rr},pi};
        end       
        location_set{rr} = xi{rr} + 1i*yi{rr};
    end

    %% Generating dynamic SLFs
    % sigma_s = 8;
    shadow = Shadowing(Xgrid,sigma_s); %larger sigma_s for heavier shadowing
    shadow_linear = 10.^(shadow/10);
    for tt = 1:T
        Xt = zeros(I,J,K);
        for rr = 1:R
            location = location_set{rr}(tt);
            loss_mat = abs(Xgrid - location);
            Sc{rr} = loss_f(loss_mat,d0,alpha).*shadow_linear;
            Sc{rr} = Sc{rr}/norm(Sc{rr},'fro');
            if directional
                Angrr = Algemat(Xgrid,location,Ang1{rr},Ang2{rr});
                Sc{rr} = Sc{rr}.*Ang{rr};
            end
            Xt = Xt + outprod(Sc{rr}, Ctrue(:,rr));
        end
        X4DT(:,:,:,tt) = Xt;
    end
end


