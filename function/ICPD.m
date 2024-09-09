function [outA,outB,outc,outs,outR,outq,outP] = ICPD(A,B,s,R,q,P,Wt,Yt,Rrank,mu, inner_iter, lambda,svalue)
%input: CPD factors A,B,C; sampling mask Wt; observation Yt.
% auxiliary variables s,R,q,P; hyper parameters: CP-rank Rank;
%regularization paramter mu; inner_iter
    I = length(A(:,1));
    J = length(B(:,1));
    Wtvec = Wt(:);
    yt = Yt(:);
    %% update C
    for nn = 1:inner_iter
        DBA = Wtvec.*kr(B,A);
        ct = (DBA'*DBA + mu*eye(Rrank))^(-1) * DBA'*yt;
%     outc = ct1;
        ct(ct<svalue) = svalue;
        %% update A
        BDC = B.* ct';
        for ii = 1:I
            s{ii} = lambda*s{ii} + Yt(ii,:).*Wt(ii,:)*BDC;
            R{ii} = lambda*R{ii} + BDC'.*Wt(ii,:)*BDC;
            A(ii,:) = s{ii}*(R{ii} + mu*eye(Rrank))^-1;
        end
        A(A<svalue) = svalue;
        %% update B
        ADC = A.* ct';
        for jj = 1:J
            q{jj} = lambda*q{jj} + ADC'.*Wt(:,jj)'*Yt(:,jj);
            P{jj} = lambda*P{jj} + ADC'.*Wt(:,jj)'*ADC;
            B(jj,:) = q{jj}'*(P{jj} + mu*eye(Rrank))^-1;
        end
        B(B<svalue) = svalue;
    end
    outc = ct;
    outA = A;
    outB = B;
    outs = s;
    outR = R;
    outq = q;
    outP = P;
end