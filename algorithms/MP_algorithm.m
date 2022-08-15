function [ x_reconstructed ] = MP_algorithm(A, b_arr)

    [m,n] = size(A);
    x_rec_saver = zeros(n, m);
    support_saver = zeros(m, m);
    eyemat = eye(n);
    for idx = 1:m
        support=[];
        res=b_arr(:, idx);
        while norm(res, 2) > 10^(-6)
            cor = abs(A' * res);
            [a, j] = max(cor);
            support = [support j];

            t = A(:,j)' * res * eyemat(:,j);

           

            x_rec =  x_rec + t * 
    
            cvx_begin
                variable x_rec(size(support, 2), 1)
                minimize(norm(b_arr(:, idx) - A(:,support) * x_rec, 2))
            cvx_end
    
            res = b_arr(:, idx) - A(:,support) * x_rec;
        end
        x_rec_saver(1:size(support, 2), idx) = x_rec;
        support_saver(idx, 1:size(support, 2)) = support;
    end
end 


