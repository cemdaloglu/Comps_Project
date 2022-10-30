clear variables;
close all;
cvx_quiet True;
%% 
S = load('data3_440.mat');
A_random = S.trainimages;
b = S.testimages;
b_random = b(:, 1);

[m,n] = size(A_random);
%% 

p = 0.1;
q = 0.5;
r = 4;
%gamma = 0.1;
gamma = 1 /norm(A_random)^2;

x0 = zeros(size(A_random, 2), 1);
x = x0;
y = x0;

% 4.1 variables
e = randn([size(A_random, 1), 1]);
e = normc(e);
w = [x0; e];

idtt = speye(size(A_random, 1));
B = [A_random, idtt];

%% Run that part to reconstruct all faces 
n_IDs = size(b,2);
fig = figure('Name', 'All reconstructed faces ');
tiledlayout(ceil(sqrt(n_IDs)),ceil(sqrt(n_IDs)))
for person_ID = 1:n_IDs
    nexttile
    sprintf('Reconstructing face of ID %d. ', person_ID)
    w_hat = fistamod(B, b(:, person_ID), w, gamma, p, q, r, n, m);
    x_hat = w_hat(1:n, 1);
    x_hat(x_hat < 0.001) = 0;
    recovered_face = A_random*x_hat + e; 
    face_plot = reshape(recovered_face, 200, 200);
    image(face_plot)
    set(gca,'XTick',[],'YTick',[])
    title(person_ID)
end
saveas(fig,'reconstruction.png')



%%
function x_sparse = fistamod(B, b, x_sprs, gamma, p, q, r, n, m)
    
    % gradient
    GradF = @(w) B' * (B * w - b);

    % Proxy
    lambda = 10;
    R = @(x) lambda .* norm(x, 1);
    FBO = @(x, y) prox(y - gamma * GradF(y), gamma, R(x), n, m);
    
    tol = 10^(-6);
    t = 1;
    x = x_sprs;
    y = x;
    
    its = 1;
    maxits = 100;
    while(its<maxits)
        
        x_old = x;
        
        x = FBO(x, y);
        
        t_old = t;
        t = (p + sqrt(q + r * t_old^2)) / 2;
        a = (t_old - 1) / t;
        
        y = x + a * (x - x_old);
        
        res = norm(x_old - x, 'fro');
        
        
        if res / n < tol
            break;
        end
        its = its + 1; 
    end
    x_sparse = x;
end
%%

function x = prox(input, gamma, R, n, m)
    l = n+m;
    cvx_begin
        variable x_hat1(l, 1)
        minimize(gamma * R + 1/2 * pow_pos(norm(x_hat1 - input, 2), 2))
    cvx_end
    x = x_hat1;
end