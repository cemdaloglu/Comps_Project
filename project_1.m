clear variables;
close all;
cvx_quiet TRUE;
%% 
n = 2^6;
m = 2^5;
A_random = randn([m, n]);
A_random = normc(A_random);

A_fourier = dftmtx(n);

A_fourier = A_fourier./abs(A_fourier);

index = randsample(1:length(A_fourier), m);
A_fourier = A_fourier(index, :);

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

for idx = 1:m
    cvx_begin
        variable x_hat2(n, 1)
        minimize(norm(x_hat2, 1))
        subject to
            A_random * x_hat2 == b_random_arr(:, idx)
    cvx_end
    x_hat2_saver(:, idx) = x_hat2;
end

x_hat2_saver(abs(x_hat2_saver) < 10^(-4)) = 0;
%% Orthogonal matching pursuit (OMP)
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
            variable x_hat3(length(support), 1)
            minimize(norm(b_random_arr(:, idx) - A_random(:, support) * x_hat3, 2))
        cvx_end

        x_sup = zeros(n, 1);
        x_sup(support, :) = x_hat3;

        res = b_random_arr(:, idx) - A_random * x_sup;
        count = count + 1;
    end
    x_hat3_saver(:, idx) = x_sup;
    support_saver(idx, 1:length(support)) = support;
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
    x_hat4_saver(:, idx) = x_hat_4;
end
%% Iterative Hard Thresholding
x_hat5_saver = zeros(n, m);

for idx = 1:m
    res = b_random_arr(:, idx);
    count = 0;
    x_hat_5 = zeros(n, 1);
    while norm(res, 2) > 10^(-6) && count < m
        cor = abs(x_hat_5 + A_random' * res);
        [a, b] = maxk(cor, idx);
        
        cor_s = x_hat_5 + A_random' * res;
        x_hat_5(b, :) = cor_s(b, :);
        res = b_random_arr(:, idx) - A_random * x_hat_5;
        count = count + 1;
    end
    x_hat5_saver(:, idx) = x_hat_5;
end
%% Basic Thresholding
x_hat6_saver = zeros(n, m);

for idx = 1:m
    cor = abs(A_random' * b_random_arr(:, idx));
    [a, b] = maxk(cor, idx);

    cvx_begin
        variable x_hat6(idx, 1)
        minimize(norm(b_random_arr(:, idx) - A_random(:, b) * x_hat6, 2))
    cvx_end

    x_sup = zeros(n, 1);
    x_sup(b, :) = x_hat6;

    x_hat6_saver(:, idx) = x_sup;
end
%% Hard Threshholding Pursuit
x_hat7_saver = zeros(n, m);

for idx = 1:m
    res = b_random_arr(:, idx);
    count = 0;
    x_hat_7 = zeros(n, 1);
    while norm(res, 2) > 10^(-6) && count < m
        cor = abs(x_hat_7 + A_random' * res);
        [a, b] = maxk(cor, idx);

        cvx_begin
            variable x_hat7(idx, 1)
            minimize(norm(b_random_arr(:, idx) - A_random(:, b) * x_hat7, 2))
        cvx_end

        x_sup = zeros(n, 1);
        x_sup(b, :) = x_hat7;

        res = b_random_arr(:, idx) - A_random * x_sup;
        count = count + 1;
    end
    x_hat7_saver(:, idx) = x_sup;
end
%% Compressive sampling matching pursuit (CoSaMP)
x_hat8_saver = zeros(n, m);

for idx = 1:m
    res = b_random_arr(:, idx);
    count = 0;
    while norm(res, 2) > 10^(-6) && count < m
        cor = abs(A_random' * res);
        [a, b] = maxk(cor, 2*idx);

        cvx_begin
            variable x_hat8(2*idx, 1)
            minimize(norm(b_random_arr(:, idx) - A_random(:, b) * x_hat8, 2))
        cvx_end
        
        x_ = zeros(n, 1);
        x_(b, :) = x_hat8;
        [a2, b2] = maxk(abs(x_), idx);
        x_sup = zeros(n, 1);
        x_sup(b2, :) = x_(b2, :); 

        res = b_random_arr(:, idx) - A_random * x_sup;
        count = count + 1;
    end
    x_hat8_saver(:, idx) = x_sup;
end
%% Subspace Pursuit
x_hat9_saver = zeros(n, m);

for idx = 1:m
    res = b_random_arr(:, idx);
    count = 0;
    while norm(res, 2) > 10^(-6) && count < m
        cor = abs(A_random' * res);
        [a, b] = maxk(cor, idx);

        cvx_begin
            variable x_hat9(idx, 1)
            minimize(norm(b_random_arr(:, idx) - A_random(:, b) * x_hat9, 2))
        cvx_end
          
        x_ = zeros(n, 1);
        x_(b, :) = x_hat9;
        [a2, b2] = maxk(abs(x_), idx);


        cvx_begin
            variable x_hat9_2(idx, 1)
            minimize(norm(b_random_arr(:, idx) - A_random(:, b) * x_hat9_2, 2))
        cvx_end
          
        x_sup = zeros(n, 1);
        x_sup(b, :) = x_hat9_2;

        res = b_random_arr(:, idx) - A_random * x_sup;
        count = count + 1;
    end
    x_hat9_saver(:, idx) = x_sup;
end
