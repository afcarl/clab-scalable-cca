function [A, B, varargout] = cca_svd(X, Y, rank)
%note: function assumes that X and Y are un-normalized (i.e., not scaled
%by number of examples like 1/(n-1))
nout = max(nargout,1)-2; %-2 because A and B are required outputs
n = size(X,1);     

k_int_x = size(X,2); %can do reduced rank, but won't recover exact CCA
k_int_y = size(Y,2); 
[u_x, s_x, v_x] = svds(X, k_int_x); %if k_int_x << size(X,2), then fast
[u_y, s_y, v_y] = svds(Y, k_int_y);
[U, S, V] = svds(u_x'*u_y, rank);

U_correct = U*sqrt(n-1);
V_correct = V*sqrt(n-1);
A = v_x*(s_x \ U_correct); 
B = v_y*(s_y \ V_correct); 

if nout > 0
    varargout{1} = diag(S)';
end
if nout > 1
    u = X*A;   
    v = Y*B; 
    varargout{2} = u;
    varargout{3} = v;
end
   