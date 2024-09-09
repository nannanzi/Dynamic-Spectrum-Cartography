% update all factors
function [A,B,R,s,P,q,c] =  UpdateAll_Grad(A,B,c,R,s,P,q,Yt,Wt,lambda,mu,svalue,max_iter)
    if nargin<11
        mu = 1e-8;
        svalue = 1e-16;
        max_iter = 10;
    end
    Frank = size(A,2); %CP-rank
    [I,J] = size(Wt);
    yvec = Yt(:);
    Wvec = Wt(:);
    SampleIndex = find(Wvec);
    yomega = yvec(SampleIndex);

    for it = 1:max_iter
        %% update C
        
        krAB = kr(B,A);
        PartkrAB = krAB(SampleIndex,:);
        LC = 1/norm(PartkrAB'*PartkrAB);
        gC = (mu*c + (c*PartkrAB' - yomega') * PartkrAB);
        c = c - 1/LC*gC;
        c(c<svalue) = svalue;
        
        %% update A
        BDC = B.* c;
        for ii = 1:I
            sii = lambda * s{ii} + Yt(ii,:) .* Wt(ii,:) * BDC;
            Rii = lambda * R{ii} + BDC' .* Wt(ii,:) * BDC;
            A(ii,:) = sii*(Rii + mu*eye(Frank))^-1;
            if it == max_iter
                s{ii} = sii;
                R{ii} = Rii;
            end
        end
         A(A<svalue) = svalue;

        %% update B
        ADC = A.* c;
        for jj = 1:J
            qjj = lambda * q{jj} + ADC' .* Wt(:,jj)' * Yt(:,jj);
            Pjj = lambda * P{jj} + ADC' .* Wt(:,jj)' * ADC;
            B(jj,:) = qjj'*(Pjj + mu*eye(Frank))^-1;
            if it == max_iter
                q{jj} = qjj;
                P{jj} = Pjj;
            end
        end
         B(B<svalue) = svalue;
    end
    
end