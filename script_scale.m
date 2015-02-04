%Testing script for speed/scalability
%generate the data
n = 10000;
p1 = 100; 
p2 = p1;
rank = 50;

%dense data
% X = randn(n, p1); 
% Y = randn(n, p2);
% mean_X = mean(X); %preprocessing: mean-center the data
% mean_Y = mean(Y);
% X = X-mean_X(ones(n,1),:);
% Y = Y-mean_Y(ones(n,1),:);

%sparse data
density = 0.01;
X = sprand(n, p1, density);
Y = sprand(n, p2, density); 
X(X ~= 0) = 1;
Y(Y ~= 0) = 1; 

tic;
[A_eig, B_eig, r_eig] = cca_eigen(X,Y,rank);
r_eig %print out
timeSpent = toc;
fprintf('Eigen decomposition CCA Time Taken: %.3f sec\n',timeSpent);

tic;
[A_qr, B_qr, r_qr] = cca_qr(X,Y,rank);
r_qr
timeSpent = toc;
fprintf('Inverse Square Root CCA Time Taken: %.3f sec\n',timeSpent);

tic;
[A_svd, B_svd, r_svd] = cca_svd(X,Y,rank);
r_svd
timeSpent = toc;
fprintf('Inverse Square Root CCA Time Taken: %.3f sec\n',timeSpent);

tic;
[A_dir, B_dir, r_dir] = cca_direct(X,Y,rank);
r_dir
timeSpent = toc;
fprintf('Inverse Square Root CCA Time Taken: %.3f sec\n',timeSpent);