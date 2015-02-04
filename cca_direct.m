function [A, B, varargout] = cca_direct(X, Y, rank)
%note: function assumes that X and Y are un-normalized (i.e., not scaled
%by number of examples like 1/(n-1))
nout = max(nargout,1)-2; %-2 because A and B are required outputs
n = size(X,1);

p1 = size(X,2);
p2 = size(Y,2);
Cxx = X'*X;
Cyy = Y'*Y;
Cxy = X'*Y;
gamma = 10^(-8);
if issparse(X)
    Cxx = Cxx + gamma*speye(p1, p1);
    Cyy = Cyy + gamma*speye(p2, p2);
else
    Cxx = Cxx + gamma*eye(p1, p1);
    Cyy = Cyy + gamma*eye(p1, p1);
end

%below code unsparsifies the matrix, but Cxx and Cyy are relatively smaller
%than X,Y so this may be preferred
Cxx_root = sqrtm(full(Cxx)); 
Cyy_root = sqrtm(full(Cyy)); 
Cxy_til = (Cxx_root \ Cxy) / Cyy_root; %not a symmetric matrix
[U, S, V] = svds(Cxy_til, rank); 
U_correct = U*sqrt(n-1);
V_correct = V*sqrt(n-1); 
A = Cxx_root \ U_correct; 
B = Cyy_root \ V_correct; 

if nout > 0
    varargout{1} = diag(S)';
end
if nout > 1
    u = X*A;   
    v = Y*B; 
    varargout{2} = u;
    varargout{3} = v;
end