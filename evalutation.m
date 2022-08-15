% evaluation 

l1_x= load('l1.mat');
l1_x_true = l1_x.x_true; 
l1_x_rec_random = l1_x.x_recovered_random;
l1_x_rec_fourier = l1_x.x_recovered_fourier;


%% ﻿plot the averaged normalized recovery error 
% ||x−x_hat||_2=||x||_2 of the different reconstruction algorithms for each s.
% Declare a recovery successful if ||x − x_hat||2=||x||2 ≤ 10−6

l1_x_rand_dif = l1_x_true - l1_x_rec_random;  % 64x32x100
l1_x_fourier_dif = l1_x_true - l1_x_rec_fourier; 

%vecnorm(l1_x_true, 2)% compute columnnorms of x_true for all 100 repetitions 1x32x100
%vecnorm(l1_x_rand_dif, 2)
%vecnorm(l1_x_fourier_dif, 2);

rec_error_random = squeeze(vecnorm(l1_x_rand_dif, 2)./vecnorm(l1_x_true, 2));
rec_error_fourier = squeeze(vecnorm(l1_x_fourier_dif, 2)./vecnorm(l1_x_true, 2));

mean_rec_error_random = mean(rec_error_random, 2);
mean_rec_error_fourier = mean(rec_error_fourier, 2);

%l1_x_rec_random_mean = mean(l1_x_rand_dif, 3);
%l1_x_rec_fourier_mean = mean(l1_x_fourier_dif, 3);

%%
s = 1:1:32;
figure
%plot(s,mean_rec_error_random, s, mean_rec_error_random + 2 )
plot(s,[mean_rec_error_random,mean_rec_error_random + 2] )
legend
yline(10^-6)
title('Averaged Normalized Reconstruction error as a function of sparsity s')
xlabel('sparsity s')
ylabel('reconstruction error')





%%
M = [[1 2 3]; [4 5 6]]
Mhat = [[2 5 7]; [4 5 7]]
MMhat = M - Mhat
N = [[-1 2 1]; [3 5 6]]
Nhat = [[0 2 0]; [3 4 5]]
NNhat = N - Nhat
Msparse = vecnorm(MMhat,2)./vecnorm(M,2)
Nsparse = vecnorm(NNhat,2)./vecnorm(N,2)
MMNNhatsparseConcat = cat(3, MMhat, NNhat)
NMmean = squeeze(mean(cat(3,Msparse,Nsparse),3))


%%
NM = squeeze(mean(cat(3,M,N),3))
NMhat = squeeze(mean(cat(3,MMhat,NNhat),3))
NMhatNorm = vecnorm(NMhat,2)
NMsparse = vecnorm(NMhat,2)./vecnorm(NM,2)
%%

%l1_rec_err_random = sum(abs(B).^2)



