function [ x_reconstructed ] = MP_algorithm(A, b_arr)

    [m,n] = size(A);
    x_reconstructed = zeros(n, m);

    for idx = 1:m
        res = b_arr(:, idx);
        count = 0;
        x_hat_4 = zeros(n, 1);
        while norm(res, 2) > 10^(-6) && count < m
            cor = abs(A' * res);
            [a, b] = max(cor);
            A_sup = A(:, b);
            t = A_sup' * res;
            x_hat_4(b, :) = x_hat_4(b, :) + t;
            res = res - t * A_sup;
            count = count + 1;
        end
        x_reconstructed(:, idx) = x_hat_4;
    end
end