% update all factors
function [A,B,R,s,P,q,C] =  UpdateAll_Incre(A,B,C,R,s,P,q,Yt,Wt,lambda,mu)
    if nargin<11
        mu = 1e-6;
%         max_iter = 10;
    end
    Frank = size(A,2); %CP-rank
    I = size(Wt,1);
    J = size(Wt,2);
    yvec = Yt(:);
    Wvec = Wt(:);
    SampleIndex = find(Wvec);
    yomega = yvec(SampleIndex);

    %% update C
    krAB = kr(B,A);
    PartkrAB = krAB(SampleIndex,:);
    ct = (PartkrAB'*PartkrAB + mu*eye(Frank))^(-1) * PartkrAB' * yomega;
    ct(ct<1e-16) = 1e-16;
    C = [C;ct'];

    %% update A
    BDC = B.* ct';
    for ii = 1:I
        s{ii} = lambda * s{ii} + Yt(ii,:) .* Wt(ii,:) * BDC;
        R{ii} = lambda * R{ii} + BDC' .* Wt(ii,:) * BDC;
        A(ii,:) = s{ii}*(R{ii} + mu*eye(Frank))^-1;
    end
    A(A<1e-16) = 1e-16;

        %% update B
    ADC = A.* ct';
    for jj = 1:J
        q{jj} = lambda * q{jj} + ADC' .* Wt(:,jj)' * Yt(:,jj);
        P{jj} = lambda * P{jj} + ADC' .* Wt(:,jj)' * ADC;
        B(jj,:) = q{jj}'*(P{jj} + mu*eye(Frank))^-1;
    end
    B(B<1e-16) = 1e-16;
    
end