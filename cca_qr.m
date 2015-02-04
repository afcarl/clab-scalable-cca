function [A, B, varargout] = cca_qr(X, Y, rank)
%note: function assumes that X and Y are un-normalized (i.e., not scaled
%by number of examples like 1/(n-1))
nout = max(nargout,1)-2; %-2 because A and B are required outputs
n = size(X,1); 

%QR not a good idea when X is tall - decomposition is expensive and
%destroys sparsity - Q is dense and very large!
r_x = qr(X,0); %QR on the full X matrix, but only gives R, hence speed
q_x = X / r_x; %using substitution to solve, hence speed
r_y = qr(Y,0);
q_y = Y / r_y;       
[U, S, V] = svds(q_x'*q_y, rank); 

U_correct = U*sqrt(n-1);
V_correct = V*sqrt(n-1); 
A = r_x \ U_correct;
B = r_y \ V_correct; 

if nout > 0
    varargout{1} = diag(S)';
end
if nout > 1
    u = q_x*U_correct;   
    v = q_y*V_correct; 
    varargout{2} = u;
    varargout{3} = v;
end
   