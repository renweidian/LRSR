


clc;
clear;
addpath(genpath('LRSR'))
addpath(genpath('functions'))
addpath(genpath('superpixel'))
S=imread('.\data\original_rosis.tif');
S=S(50:177,50:177,11:end);
S=double(S);
S=S/max(S(:));




  [M N L]=size(S); 
  S1=hyperConvert2D(S);
R=load('.\data\ikonos_spec_resp.mat');
R=R.ikonos_sp;
[~, valid_ik_bands] = intersect(R(:,1), 430:860);
no_wa = length(valid_ik_bands);
xx  = linspace(1, no_wa, L);

x = 1:no_wa;
F = zeros(5, L);
for i = 1:5 % 1 - pan; 2 - blue; 3 - green; 4 - red; 5 - NIR
    F(i,:) = spline(x, R(valid_ik_bands,i+1), xx);
end
F = F(2:5,:);


  for band = 1:size(F,1)
        div = sum(F(band,:));
        for i = 1:size(F,2)
            F(band,i) = F(band,i)/div;
        end
 end



SNRm=40;
Y               =    F*S1;
kernel_type     =    {'Gaussian_blur'};
sz=[M N];
sf              =    4;
s0=2;
par             =    Parameters_setting( sf, kernel_type, sz,s0 );
 sigmah1 = sqrt(sum(Y (:).^2)/(10^(SNRm/10))/numel(Y ));
 rng(10,'twister')
Y = Y + sigmah1*randn(size(Y));
MSI=hyperConvert3D(Y,sz(1),sz(2));

   X               =    par.H(S1);
   SNRm=30;
 sigmah = sqrt(sum(X (:).^2)/(10^(SNRm/10))/numel(X ));
rng(10,'twister')
X = X + sigmah*randn(size(X));
HSI=hyperConvert3D(X,sz(1)/sf,sz(2)/sf);


 
 
 
 
 
 
%% LRSR
par.total_patches =200;
par.lambda1=1e-4;
par.lambda2=1e-3;
par.K=24;
par.P=F;
[M1]     =    NSSR_HSI_SR( S1, sf,par,sz,X,Y );
Z5=hyperConvert3D(M1, M, N);
[psnr5,rmse5, ergas5, sam5, uiqi5,ssim5,DD5,CC5] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z5)), 0, 1.0/sf);



