function [ x_reconstructed ] = HTP_algorithm(A, b_arr)
    [m,n] = size(A);
    x_reconstructed = zeros(n, m);
    
    for idx = 1:m
        res = b_arr(:, idx);
        count = 0;
        x_hat_7 = zeros(n, 1);
        while norm(res, 2) > 10^(-6) && count < m
            cor = abs(x_hat_7 + A' * res);
            [a, b] = maxk(cor, idx);
    
            cvx_begin
                variable x_rec(idx, 1)
                minimize(norm(b_arr(:, idx) - A(:, b) * x_rec, 2))
            cvx_end
    
            x_sup = zeros(n, 1);
            x_sup(b, :) = x_hat7;
    
            res = b_arr(:, idx) - A * x_sup;
            count = count + 1;
        end
        x_reconstructed(:, idx) = x_sup;
    end
end 