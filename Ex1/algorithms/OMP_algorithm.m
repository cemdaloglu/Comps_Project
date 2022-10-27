function [ x_reconstructed ] = OMP_algorithm(A, b_arr)

    [m,n] = size(A);
    x_reconstructed = zeros(n, m);
    support_saver = zeros(m, m);
    for idx = 1:m
        support=[];
        res=b_arr(:, idx);
        count = 0;
        while norm(res, 2) > 10^(-6) && count < m
            cor = abs(A' * res);
            [a, b] = max(cor);
            support = [support b];
            A_sup = zeros(m, n);
            A_sup(:, support) = A(:, support);

    
            cvx_begin
                variable x_rec(size(support, 2), 1)
                minimize(norm(b_arr(:, idx) - A(:, support) * x_rec, 2))
            cvx_end

            x_sup = zeros(n, 1);
            x_sup(support, :) = x_rec;
    
            res = b_arr(:, idx) - A * x_sup;
            count = count + 1;
        end
        x_reconstructed(:, idx) = x_sup;
        support_saver(idx, 1:length(support)) = support;
    end
end 


