function B=sparsepca(X)
%load('Faces.mat');
% X=hyperConvert2D(im)';
%X=generatedata(im);
K=1;
delta=inf;
stop = -[250 250 250];
% stop=0;
maxiter = 3000;
convergenceCriterion = 1e-9;
verbose = false;
[B SD L D paths] = spca(X, [], K, delta, stop, maxiter, convergenceCriterion, verbose);
% T=SL';
end