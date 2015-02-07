function [A, B, varargout] = cca_ils(X, Y, rank)
%note: function assumes that X and Y are un-normalized (i.e., not scaled
%by number of examples like 1/(n-1))

nout = max(nargout,1)-2; %-2 because A and B are required outputs
n = size(X,1);

p1 = size(X,2);
p2 = size(Y,2); 
Cxx = X'*X;
Cyy = Y'*Y;
gamma = 1;
iters = 1000;

if issparse(X)
    Cxx = Cxx + gamma*speye(p1, p1);
    Cyy = Cyy + gamma*speye(p2, p2);
else
    Cxx = Cxx + gamma*eye(p1, p1);
    Cyy = Cyy + gamma*eye(p2, p2);
end

G = randn(p1, rank);
Hy = (Y / Cyy)*Y';
Hx = (X / Cxx)*X';
X_cur = X*G;
Y_cur = 0;
for it = 1:iters    
    Y_cur = Hy*X_cur;
    [Y_cur, r_y] = qr(Y_cur, 0); %QR decomp in every iteration for stability reasons (slower)
    X_cur = Hx*Y_cur;   
    [X_cur, r_x] = qr(X_cur, 0);
end

lambda = diag(X_cur'*Y_cur); %correlations should form diagonal of cross-covariance of projected data
A = X \ X_cur; %pseudo-inverse of X used since it is not square
B = Y \ Y_cur;

if nout > 0
    varargout{1} = lambda;
end

if nout > 1
    u = X_cur;
    v = Y_cur;
    varargout{2} = u;
    varargout{3} = v;
end
   