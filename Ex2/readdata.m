%% load data 
path = 'data';
load([path '/data3_' '440' '.mat']);

%%
A = trainimages;
B = testimages;
gnd_Train = trainids;
gnd_Test = testids;
imDims = [200, 200];

% Normalizing dictionary and observations
for i = 1 : size(A,2)
    A(:,i) = A(:,i)/norm(A(:,i)) ;
end
for i = 1 : size(B,2)
    B(:,i) = B(:,i)/norm(B(:,i)) ;
end
[m, n] = size(A);
if mod(n,2) ~= 0 % make sure there are even number of images
    A = A(:,1:end-1);
end

%%
f = fft2(B,1024,1024) % padding to get higher resolution in frequency
imshow(fftshift(f)); % plotting
