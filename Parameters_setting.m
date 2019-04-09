function   par   =  Parameters_setting( sf, kernel_type, sz,s0 )
% par.h          =    sz(1);
% par.w          =    sz(2);

if strcmp(kernel_type, 'uniform_blur')
    psf        =    ones(sf)/(sf^2);
elseif strcmp(kernel_type, 'Gaussian_blur')
    psf        =    fspecial('gaussian',7,2);
end

   
 par.fft_B      =    psf2otf(psf,sz);
  par.B     =    ifft2(par.fft_B) ;
 
% psfSize=size(psf);
%   padSize = sz - size(psf);
%    psf     = padarray(psf, padSize, 'post');
% 
%    psf    = circshift(psf,-floor(psfSize/2));
%    par.B     =    psf;
%  par.fft_B      =    fft2(par.B);
   
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);


    