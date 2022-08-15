clear all
%% Project: exercise 1


recovery("l1")

function recovery(algorithm)


    % Determine where your m-file's folder is.
    folder = fileparts(which("project_ex1"));
    % Add that folder plus all subfolders to the path.
    addpath(genpath(folder));
    
    %%
    n = 2^6;
    m = n/2;
    
    A_random = randnc(m,n);
    vecnorm(A_random);
    
    A_fourier_tmp = dftmtx(n);
    A_fourier_rowindex = randi([1, n], 1, m);
    A_fourier_tmp = A_fourier_tmp(A_fourier_rowindex,:);
    A_fourier = A_fourier_tmp ./ ( sum( ( A_fourier_tmp .* conj(A_fourier_tmp) ).^2 ) ).^0.5; % normalize
    
    % Saving the sensors for reproducibility
    save([folder '/sensors'], "A_random", "A_fourier")

    %% Read in again
    A_loaded = load('sensors.mat');
    A_fourier = A_loaded.A_fourier;
    A_random = A_loaded.A_random;

    
    %%
    nrReps = 100;
    
    x_true = zeros(n, m, nrReps);
    x_recovered_random = zeros(n, m, nrReps);
    x_recovered_fourier = zeros(n, m, nrReps);
    %disp(size(x_recovered_random))
    
    %% ﻿Make 100 repetitions for each sparsity value s
    for nrep = 1:nrReps
     
        fprintf(1,'Iteration %3.0f/%3.0f -> ',nrep,nrReps);
        
        %% ﻿For each value of the sparsity 1 ≤ s ≤ m,
        % ﻿construct a random s-sparse vector x by ﻿drawing the nonzero locations
        % uniformly at random, and the nonzero values from a Gaussian
        x_sparse_arr = randnc(n,m);
        for s = 1:m
            r = randi([1, n], 1, s);
            x_sparse_arr(setdiff(1:end,r), s) = 0;
        end

        x_true(:,:,nrep) = x_sparse_arr;
    
        b_fourier_arr = A_fourier * x_sparse_arr;
        b_random_arr = A_random * x_sparse_arr;
        
    
        %% Run algorithm
        % for every problem instance defined by s and recover x^ from
        % measurements b := Ax. Use for all algorithms a
        % residual-based stopping criterion, e.g. kb−Ax^k2 ≤ 10−6
        cvx_quiet TRUE
    
        if strcmp(algorithm, "l1")
            x_rec_random = l1_algorithm(A_random, b_random_arr);
            x_rec_fourier = l1_algorithm(A_fourier, b_fourier_arr);
        elseif strcmp(algorithm, "OMP")
            x_rec_random  = OMP_algorithm(A_random, b_random_arr);
            x_rec_fourier = OMP_algorithm(A_fourier, b_fourier_arr);
        elseif strcmp(algorithm, "MP")
            x_rec_random = MP_algorithm(A_random, b_random_arr);
            x_rec_fourier = MP_algorithm(A_fourier, b_fourier_arr);
        elseif strcmp(algorithm, "IHT")
            x_rec_random = IHT_algorithm(A_random, b_random_arr, s); % TODO: loop over s
            x_rec_fourier = IHT_algorithm(A_fourier, b_fourier_arr, s);    
        elseif strcmp(algorithm, "CoSaMP")
            x_rec_random = CoSaMP_algorithm(A_random, b_random_arr, s);
            x_rec_fourier = CoSaMP_algorithm(A_fourier, b_fourier_arr, s);    
        elseif strcmp(algorithm, "BT")
            x_rec_random = BT_algorithm(A_random, b_random_arr, s);
            x_rec_fourier = BT_algorithm(A_fourier, b_fourier_arr, s); 
        elseif strcmp(algorithm, "HTP")
            x_rec_random = HTP_algorithm(A_random, b_random_arr);
            x_rec_fourier = HTP_algorithm(A_fourier, b_fourier_arr); 
        elseif strcmp(algorithm, "SP")
            x_rec_random = SP_algorithm(A_random, b_random_arr, s);
            x_rec_fourier = SP_algorithm(A_fourier, b_fourier_arr, s); 
        else
            warning('Chose one of "OMP", "MP", "IHT", "CoSaMP", "BT", "HTP", "SP"')
        end
    
        % save in matrix
        x_recovered_random(:,:,nrep) = x_rec_random;
        x_recovered_fourier(:,:,nrep) = x_rec_fourier;
        save(strcat(folder, '/results/' ,algorithm), "x_true", "x_recovered_random", "x_recovered_fourier")
   
    end
    
    save(strcat(folder, '/results/' ,algorithm),"x_true", "x_recovered_random", "x_recovered_fourier")
    disp("done")
end

%% ﻿plot the averaged normalized recovery error ﻿of the
% different reconstruction algorithms for each s