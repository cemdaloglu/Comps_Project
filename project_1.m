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
    while norm(res, 2) > 10^(-6)
        cor = abs(A_random' * res);
        [a, b] = max(cor);
        support = [support b];

        cvx_begin
            variable x_hat3(size(support, 2), 1)
            minimize(norm(b_random_arr(:, idx) - A_random(:,support) * x_hat3, 2))
        cvx_end

        res = b_random_arr(:, idx) - A_random(:,support) * x_hat3;
    end
    x_hat3_saver(1:size(support, 2), idx) = x_hat3;
    support_saver(idx, 1:size(support, 2)) = support;
end
