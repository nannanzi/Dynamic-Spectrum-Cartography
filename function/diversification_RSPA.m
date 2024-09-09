function [k,Kfin,err] = diversification_RSPA(R,options);

% Selection step for the robust successive projection algorithm (RSPA); see
% Successive Projection Algorithm Robust to Outliers, N. Gillis, June 2019.
%
% *** Description ***
% RSPA sequentially generates candidate data points by generating different
% strongly convex functions f_i's (1<=i<=d). For each function, there is an
% associated data point R(:,k(i)) that maximizes it.
% In the last step, RSPA returns the data point with smallest
% component-wise l_p residual.
%
% Remark. The implementation could be improved by not computing explicitely
% the residual, by using the formula
%            ||(I-uu^T)v||^2 = ||v||^2 - (u^T v)^2;
% See the SPA implementation FastSepNMF.m at
% https://sites.google.com/site/nicolasgillis/code
%
%
% ****** Input ******
% R          : a data matrix
% --- options ---
% d          : number of functions f_i's generated hence number of data
%               points added to the set of candidates.
%               - default: d=10
% p          : the component-wise l_p norm is used to compute the residual
%               error and choose among the candiate data points
%               - default: p=1
% beta       : > 1, prameter that determines the functions f_i's, see the
%               paper above for more details.
%               - default: beta=4
%
% ****** Output ******
% k        : an index that is a vertex of conv(X) and that minmizes the
%            component-wise l_p norm of the residual after orthogonal
%            projection, that is, that minimizes
%            || (I-X(:,k)X(:,k)^T/||X(:,k)||_2^2) X ||_p ()
% Kfin     : set of d candidate indices maximizing some strongly convex
%            functions f_i (1 <= i <= d).
% err      : d residuals (see above) corresponding to the d indices in Kfin.

% Options
if nargin <= 1
    options = [];
end
if ~isfield(options,'d')
    options.d = 10;
end
if ~isfield(options,'p')
    options.p = 1;
end
if ~isfield(options,'beta')
    options.beta = 4;
end
% Algorithm
Kfin = []; % Will keep in memory the candidate indices
Y = R;
for i = 1 : options.d
    normY = sum(Y.^2);
    [a,b] = max(normY);
    % Error computation
    Kfin(i) = b;
    if options.d == 1 % SPA
        err(i) = 0;
    else
        ui = R(:,Kfin(i)) / norm( R(:,Kfin(i)) );
        Ri = R - ui*(ui'*R);
        %  Criterion to evaluate the error of the residual
        err(i) = sum( sum( Ri.^2 ).^(options.p/2) );
        if i < options.d % At the last step, it is useless to update Y
            uyi = Y(:,Kfin(i)) / norm( Y(:,Kfin(i)) );
            Yi = Y - uyi*(uyi'*Y);
            normYi = sum(Yi.^2);
            [ap,bp] = max(normYi);
            if ap < 1e-12 % this may happen if X has rank one
                break;
            else
                % Diversification
                x = Y(:,Kfin(i));
                u = x/norm(x);
                y = Y(:,bp);
                alpha(i) = 1 - sqrt( 1 - (options.beta*norm(x)^2 - norm(y)^2)/(options.beta*(u'*x)^2 - (u'*y)^2) );
                if alpha(i) >= 1 || alpha(i) <= 0
                    % This should in theory not happen ~ safety procedure
                    error('alpha should be between 0 and 1');
                end
                Y = Y - alpha(i) * u*(u'*Y);
            end
        end
    end
end
% Keep the best candidate index w.r.t. err
[a,b] = min(err);
k = Kfin(b);
end