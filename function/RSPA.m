function K = RSPA(X,r,options)

% Robust successive projection algorithm for separable NMF; see 
% Successive Projection Algorithm Robust to Outliers, N. Gillis, July 2019.
% 
% *** Description ***
% At each step of the algorithm, the column of X maximizing a strongly
% convex function and minimizing the residual error is extracted (see
% diversification_RSPA.m) and X is updated by projecting its columns onto 
% the orthogonal complement of the extracted column. 
%
%
% ****** Input ******
% X = WH + N : a (normalized) noisy separable matrix, that is, W is full rank, 
%              H = [I,H']P where I is the identity matrix, H'>= 0 and its 
%              columns sum to at most one, P is a permutation matrix, and
%              N is sufficiently small. 
% r          : number of columns to be extracted.   
% --- options ---
% d          : number of functions f_i's generated hence number of data
%               points added to the set of candidates. 
%               - default: d=10 
% p          : the component-wise l_p norm is used to compute the residual 
%               error and choose among the candiate data points
%               - default: p=1
% beta       : > 1, prameter that determines the functions f_i's, see the
%               paper above for more details
%               - default: beta=4
% normalize  : normalize=1 will scale the columns of X so that they sum to 
%              one, so that the matrix H will satisfy the sum-to-one 
%              assumption. 
%              - default: normalize=0 (no scaling). 
%              For example, in hyperspectral imaging, this assumption is 
%              often already satisfied and normalization is not necessary. 
% 
% ****** Output ******
% K        : index set of the extracted columns. 

% Options 
if nargin <= 2
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
if ~isfield(options,'normalize')
    options.normalize = 0; 
end
[m,n] = size(X); 
if options.normalize == 1
    % Normalization of the columns of X so that they sum to one
    D = spdiags(((sum(X)+1e-16).^(-1))', 0, n, n); 
    X = X*D; 
end
R = X; 
normX = sum(X.^2); 
normR = normX; 
nX = max(normX); 
i = 1; 
% Perform r steps (unless the approximation error <= 10^-12)
while i <= r && max(normR) > 1e-12 
    % Select the column of R with the diversification procedure
    b = diversification_RSPA(R,options); 
    % Update the index set, and the residual 
    K(i) = b; 
    u = R(:,b)/norm(R(:,b)); 
    R = R - u*(u'*R); 
    normR = sum(R.^2); 
    i = i + 1; 
end

end 