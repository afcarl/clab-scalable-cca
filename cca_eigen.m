function [A, B, varargout] = cca_eigen(X, Y, rank)
%note: function assumes that X and Y are un-normalized (i.e., not scaled
%by number of examples like 1/(n-1))
nout = max(nargout,1)-2; %-2 because A and B are required outputs
n = size(X,1);

%method for computing cov matrices below preserves sparsity
%if using matlab's inbuilt cov, it does not
p1 = size(X,2);
p2 = size(Y,2); 
Cxx = X'*X;
Cyy = Y'*Y;
Cxy = X'*Y; 
Cyx = Cxy';
gamma = 10^(-8);
%regularizer: depends on sparse/dense
if issparse(X)
    Cxx = Cxx + gamma*speye(p1, p1);
    Cyy = Cyy + gamma*speye(p2, p2);
else
    Cxx = Cxx + gamma*eye(p1, p1);
    Cyy = Cyy + gamma*eye(p2, p2);
end

Rx = chol(Cxx); %Cxx = Rx'*Rx - upper triangular; feature space size so not too bad
invRx = inv(Rx); 
Z = invRx'*Cxy*(Cyy\Cyx)*invRx; %as per Lagrange multiplier result
Z = 0.5*(Z'+Z); %ensures symmetric matrix - may numerically not be the case before
[U, lambda_2] = eigs(Z, rank); %rank-reduced and eigenvalues are sorted in descending order

%note: why the sqrt(n-1) correction term?
%Matlab's canoncorr function outputs scaled versions of projection matrices
%(canonical coefficients) A and B such that the covariance matrices of the
%canonical variables u and v derived from these projection matrices (i.e.,
%u = XA and v = YB) are orthonormal (not orthogonal). The sample covariance
%matrices u^T*u and v^T*v are un-normalized, so they must be divided by
%1/(n-1). This means that the un-normalized covariance matrices should
%equate to not the identity matrix I, but rather (n-1)*I. So, we need to
%multiply our recovered U and V matrices from the SVD by a correction term
%to match the Matlab output. We can simply scale orthonormal bases
%recovered by SVD (U and V) by sqrt(n-1) as a result. 
U = U*sqrt(n-1); %need to make this adjustment because covariance matrices are unnormalized
lambda = diag(sqrt(real(lambda_2)))'; %original lambda obtained is lambda^2
A = Rx \ U; %actual A values are this; above is a computational trick for Cholesky
B = (Cyy\Cyx)*A; %as per Lagrangian
B = B./repmat(lambda,p2,1); 

if nout > 0
    varargout{1} = lambda;
end
if nout > 1
    u = X*A;   
    v = Y*B; 
    varargout{2} = u;
    varargout{3} = v;
end
   