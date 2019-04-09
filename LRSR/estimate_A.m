function  [ A]    =   estimate_A( D, D0, X, Y, par, Z_ori, sf, sz,A,lables  )
FB=par.fft_B;  
FBC=par.fft_BT;
nl=sz(1);
V1=ConvC(A, FB, nl);
V2=A;
V3=A;
V4=A;
shift=floor(sf/2)-1;
mask = zeros(sz(1), sz(2));
mask(shift+1:sf:sz(1), shift+1:sf:sz(1)) = 1;
maskim = repmat(mask, [1, 1, size(D,2)]);
mask = im2mat(maskim);
D02=D0'*D0;
D2=D'*D;

mu=1e-3;
II_FB = 1./((abs(FB.^2)) + 3);
IBD_B= FBC./((abs(FB.^2)) + 3);
G1=zeros(size(V1));
G2=zeros(size(V2));
G3=zeros(size(V3));
G4=zeros(size(V4));
  X = hyperConvert3D(X, sz(1)/sf, sz(2)/sf);
 Yhim_up = upsamp_HS(X, sf, sz(1), sz(2), shift);
% X2=zeros(sz(1),sz(2),size(D,1));
% X2(shift+1:sf:end,shift+1:sf:end,:)=X;
X2 = im2mat(Yhim_up);
DTX=D'*X2;
D0Y=D0'*Y;
lamba=2;

for i=1:100
%% A
%  gg(i)=(length(find(abs(A)<1e-6)))/numel(A)  
   A = ConvC(V1+G1, IBD_B, nl) + ConvC(V2+G2, II_FB, nl) +  ConvC(V3+G3, II_FB, nl)+ConvC(V4+G4, II_FB, nl);

   
   
% A = ConvC(ConvC(V1+G1,FBC,nl) + (V2+G2) + (V3+G3),II_FB,nl);

 %% V1
NU1 =ConvC(A, FB, nl)-G1;
      
        V1_com =(D2+mu*eye(size(D,2)))\(DTX+mu*NU1); 
               
V1=NU1.*(1-mask)+  V1_com .*mask; 

%% V2
V2=(lamba*D02+mu*eye(size(D0,2)))\(lamba*D0Y+mu*(A-G2));

%% V3
 for ll=0:1:par.total_patches-1
      index=find(lables==ll);
      ggg=A-G3;
    temp=ggg(:,index);
    
  
 aa3  =   repmat(mean( temp, 2 ),1,size(temp,2));
            aa4    =   temp - aa3;    
 C=par.lambda2/(2*mu);
   aa5= prox_nuclear(aa4,C)+ aa3;
 V3(:,index)=aa5;
 end
  V4=soft(A-G4, par.lambda1/(2*mu));
%  V4=ClosedWL1(A-G4,1e-5/(2*mu),eps);  
%% G1 G2 G3
G1=G1+(V1-ConvC(A, FB, nl));
G2=G2+(V2-A);
G3=G3+(V3-A);
G4=G4+(V4-A);
% aa1(i)=norm(V2-A)
 a1 =double(im2uint8(D*A)); 
b1 =double(im2uint8(Z_ori)); 
 MSE             =    mean( mean( (a1-b1).^2 ) );
rmse=  sqrt(MSE);
        PSNR      =   10*log10(255^2/MSE);                
        disp( sprintf('Iter %d, RMSE = %3.3f, PSNR = %3.2f', i, rmse, PSNR) );  


end




                              

        

     
