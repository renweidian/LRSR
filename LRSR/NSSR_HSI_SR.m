function    [HSI_res1 ]    =    NSSR_HSI_SR( Z_ori, sf, par,sz,X,Y )





    D            =    Nonnegative_DL( X, par ); 


Y1=normalize(Y');
B=sparsepca(Y1);
Y1=Y1*B;
Y1=reshape(Y1,sz(1),sz(2));

[labels,~]=suppixel(Y1,par.total_patches);
 


% 
 A         =   zeros( par.K, size(Y,2) );
% AS = sunsal(D, X,'POSITIVITY','yes','ADDONE','yes');
%  A = hyperConvert2d(imresize(hyperConvert3d(AS,sz(1)/sf),sf));
% 
for i=1
 A       =    estimate_A(D, par.P*D, X, Y, par, Z_ori, sf, sz,A,labels );
 a1 =double(im2uint8(D*A)); 
b1 =double(im2uint8(Z_ori)); 
 MSE2             =    mean( mean( (a1-b1).^2 ) );
rmse1(i)=  sqrt(MSE2)
%   D = update_D(D,par.H(A),A,X,Y,par.P);
a1 =double(im2uint8(D*A)); 
b1 =double(im2uint8(Z_ori)); 
MSE2             =    mean( mean( (a1-b1).^2 ) );
rmse(i)=  sqrt(MSE2)
end
HSI_res1=D*A;
% HSI_res2        =    alternating_back_projection( HSI_res1,X,Y,par.P,create_H(sz, sf), par );



