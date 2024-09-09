function [Si,Ci] = CBTD(I,J,Yob3,L,Rest,SampleIndex)
    M = length(Yob3(:,1));
    P = zeros(M,I*J);
    for mm = 1:M
        row_selection = SampleIndex(mm);
        P(mm,row_selection) = 1;
    end
    lambda = 1e-8;
    svalue = 1e-16;
%     svalue = 0;
    tic;
    iter = 30;
    [S_init,C] = NMF_HALS(Yob3,Rest,50);%Initialized with ALS
    S = zeros(I*J,Rest);
    S(SampleIndex,:) = S_init;
    %%Projection Gradient
    for ii = 1:iter    
        LS = 1/norm(S'*S);
        PS = P*S;
% % % %     if ii ==1
% % % %         PS = S_init;
% % % %     end
        SPPS = PS'*PS;
        C = C - rand*LS * (lambda*C + (C*PS' - Yob3') * PS);
        C(C<svalue) = svalue;
        LC = 1/norm(C'*C);
        S = S - rand*LC * (P'*(PS*C'- Yob3)*C); 
        for rr = 1:Rest
            sr = reshape(S(:,rr),[],1);
            SrMat = reshape(sr,I,J);
            [Us,Ss,Vs] = svds(SrMat,L);
            SrMat = Us*Ss*Vs';
            SrMat(SrMat<svalue) =svalue;
            sr = reshape(SrMat,[],1);
            S(:,rr) = sr;
        end
    end
    Ci = C;
    Si = Yob3/Ci';
    Si(Si<svalue) = svalue;
end