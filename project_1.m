clear variables;
close all;
%% 
n = 2^6;
m = 2^5;
A_random = randn([m, n]);
A_random = normc(A_random);
%%
A_fourier = dftmtx(n);

A_fourier = A_fourier./abs(A_fourier);

index = randsample(1:length(A_fourier), m);
A_fourier = A_fourier(index, :);
%%
sparse_vector = zeros([n, m]);
for s = 1:m
    nnz_values = randn([s, 1]);
    index = randsample(1:length(sparse_vector), s, true);
    sparse_vector(index, s) = nnz_values;
end
b_fourier_arr = A_fourier * sparse_vector;
b_random_arr = A_random * sparse_vector;
%% l1-minimization, fourier
x_hat1_saver = zeros(n, m);
b_fourier_arr = b_fourier_arr';

for idx = 1:m
    cvx_begin
        variable x_hat1(n, 1)
        minimize(norm(x_hat1, 1))
        subject to
            A_fourier * x_hat1 == b_fourier_arr(:, idx)
    cvx_end
    x_hat1_saver(:, idx) = x_hat1;
end
%% l1-minimization, gaussian

x_hat2_saver = zeros(n, m);
b_random_arr = b_random_arr';

for idx = 1:m
    cvx_begin
        variable x_hat2(n, 1)
        minimize(norm(x_hat2, 1))
        subject to
            A_random * x_hat2 == b_random_arr(:, idx)
    cvx_end
    x_hat2_saver(:, idx) = x_hat2;
end
%% OMP
x_hat3_saver = zeros(n, m);
support_saver = zeros(m, m);
for idx = 1:m
    support=[];
    res=b_random_arr(:, idx);
    count = 0;
    while norm(res, 2) > 10^(-6) && count < m
        cor = abs(A_random' * res);
        [a, b] = max(cor);
        support = [support b];
        A_sup = zeros(m, n);
        A_sup(:, support) = A_random(:, support);

        cvx_begin
            variable x_hat3(n, 1)
            minimize(norm(b_random_arr(:, idx) - A_sup * x_hat3, 2))
            subject to
            x_hat3(support, :) ~= 0;
        cvx_end

        res = b_random_arr(:, idx) - A_sup * x_hat3;
        count = count + 1;
    end
    x_hat3_saver(:, idx) = x_hat3;
    support_saver(idx, :) = support;
end
%% Matching Pursuit
x_hat4_saver = zeros(n, m);

for idx = 1:m
    res = b_random_arr(:, idx);
    count = 0;
    x_hat_4 = zeros(n, 1);
    while norm(res, 2) > 10^(-6) && count < m
        cor = abs(A_random' * res);
        [a, b] = max(cor);
        A_sup = A_random(:, b);
        t = A_sup' * res;
        x_hat_4(b, :) = x_hat_4(b, :) + t;
        res = res - t * A_sup;
        count = count + 1;
    end
    x_hat4_saver(:, idx) = x_hat4;
end
%% Iterative Hard Thresholding
