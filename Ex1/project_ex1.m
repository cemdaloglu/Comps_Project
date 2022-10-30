clear all
%% Project: exercise 1

%% TODO: choose which matrix/ matrices to use fo the algorithms.
% Only use one at once.
%useRandom = false;
%useFourier = true;

useRandom = true;
useFourier = false;

rep_start = 1;

n = 2^6; 
   
% Determine where your m-file's folder is.
folder = fileparts(which("project_ex1"));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% Call the algorithm 
recovery("IHT", n, useRandom, useFourier, folder, rep_start)
 
function recovery(algorithm, n, useRandom, useFourier, folder, rep_start)
    m = n/2;

    % ﻿Make 100 repetitions for each sparsity value s
    nrReps = 100;
    disp(strcat(folder, '/results/', algorithm, '_fourier.mat'))

    if useRandom
        [A, ~ , x_true] = generate_and_load_sensors_xtrue(n, folder);
        if isfile(strcat(folder, '/results/', algorithm, '_random.mat'))
            disp("Loading x_recovered random...")
            x_rec = load(strcat(algorithm, '_random.mat'));
            x_recovered = x_rec.x_recovered;
        else
            x_recovered = zeros(n, m, nrReps);
        end

    elseif useFourier
        [~, A , x_true] = generate_and_load_sensors_xtrue(n, folder);
        if isfile(strcat(folder, '/results/', algorithm, '_fourier.mat'))
            disp("Loading x_recovered fourier ...")
            x_rec = load(strcat(algorithm, '_fourier.mat'));
            x_recovered = x_rec.x_recovered;
        else 
            x_recovered = zeros(n, m, nrReps);
        end
    end

    b_arr = A * x_true;
  
    for nrep = rep_start:nrReps
     
        fprintf(1,'Iteration %3.0f/%3.0f -> ',nrep,nrReps);
  
        %% Run algorithm
        % for every problem instance defined by s and recover x^ from
        % measurements b := Ax. Use for all algorithms a
        % residual-based stopping criterion, e.g. kb−Ax^k2 ≤ 10−6
        cvx_quiet TRUE
    
        if strcmp(algorithm, "l1")
            x_rec = l1_algorithm(A, b_arr);
        elseif strcmp(algorithm, "OMP")  
            x_rec = OMP_algorithm(A, b_arr);
        elseif strcmp(algorithm, "MP")
            x_rec = MP_algorithm(A, b_arr);
        elseif strcmp(algorithm, "IHT")
            x_rec = IHT_algorithm(A, b_arr);   
        elseif strcmp(algorithm, "CoSaMP")
            x_rec = CoSaMP_algorithm(A, b_arr);  
        elseif strcmp(algorithm, "BT")
            x_rec = BT_algorithm(A, b_arr); 
        elseif strcmp(algorithm, "HTP")
            x_rec = HTP_algorithm(A, b_arr); 
        elseif strcmp(algorithm, "SP")
            x_rec = SP_algorithm(A, b_arr); 
        else
            warning('Chose one of "OMP", "MP", "IHT", "CoSaMP", "BT", "HTP", "SP"')
        end
    
        % save in matrix
        x_recovered(:,:,nrep) = x_rec;
        disp("DONE - saving iteration...")
        if useRandom
            save(strcat(folder, '/results/' ,algorithm,"_random"), "x_recovered")
        elseif useFourier
            save(strcat(folder, '/results/' ,algorithm,"_fourier"), "x_recovered")
        end
   
    end
    
    if useRandom
        save(strcat(folder, '/results/' ,algorithm,"_random"), "x_recovered")
    elseif useFourier
        save(strcat(folder, '/results/' ,algorithm,"_fourier"), "x_recovered")
    end
    disp("done")
end




%% Function to generate the sensor matrices and x_true or read them, 
% if they already exist 
% Input: n: number of columns
function [A_random, A_fourier, x_true] = generate_and_load_sensors_xtrue(n, folder)
    m = n/2;

    if isfile([folder '/x_true.mat'])
        %% Read in again
        disp("Loading x_true matrices...") 
        x_loaded = load('x_true.mat');
        x_true = x_loaded.x_true;
    else
    %% ﻿For each value of the sparsity 1 ≤ s ≤ m,
    % ﻿construct a random s-sparse vector x by ﻿drawing the nonzero locations
    % uniformly at random, and the nonzero values from a Gaussian
    x_true = randnc(n,m);
    for s = 1:m
        r = randi([1, n], 1, s);
        x_true(setdiff(1:end,r), s) = 0;
    end

    % Saving the x_true for reproducibility
    save([folder '/x_true'], "x_true")
    end

    if isfile([folder '/sensors.mat'])
        %% Read in again
        disp("Loading sensor matrices...") 
        A_loaded = load('sensors.mat');
        A_fourier = A_loaded.A_fourier;
        A_random = A_loaded.A_random;

    else
        disp("Creating and saving sensor matrices...") 
        A_rand = randnc(m,n);
        A_random = normc(A_rand);

        % Fourier matrix 
        A_four = dftmtx(n);
        A_fourier_norms = splitapply(@norm,A_four,1:size(A_four,2));
        A_fourier_tmp = A_four./ A_fourier_norms;
        
        index = randsample(1:length(A_fourier_tmp), m); % uniformly at random choose m rows 
        A_fourier = A_fourier_tmp(index, :);

        % Saving the sensors for reproducibility
        save([folder '/sensors'], "A_random", "A_fourier")
    end 
    
end