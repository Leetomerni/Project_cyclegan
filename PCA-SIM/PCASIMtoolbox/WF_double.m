function [out]=WF_double(Iraw)
N=size(Iraw,3);
NPixel=size(Iraw,1);

for I = 1:N
    IIrawFFT(:,:,I) = fftshift(fft2(Iraw(:,:,I)));
end

WF = zeros(NPixel,NPixel);                                                  % Widefield image sum                      
[x,y] = meshgrid(1:NPixel,1:NPixel);                                        % Image coordinate
Tdeconv = zeros(NPixel,NPixel);                                             % Deconvoluted imgae sum
WFdeconv = zeros(NPixel,NPixel,3);              	                        % Deconvoluted imgae
WFdeconvFFT = zeros(NPixel,NPixel,3);                                       % Deconvoluted imgae spectrum               

for i = 1:3
    for j =1:3
        Tdeconv = Tdeconv+Iraw(:,:,(i-1)*3+j);
        WF(:,:)=WF(:,:)+Iraw(:,:,(i-1)*3+j);
    end
    WFdeconv(:,:,i) = Tdeconv/3;
    WFdeconvFFT(:,:,i) = fftshift(fft2(WFdeconv(:,:,i)));
end
WF = WF/N;                                                                                             

%% Wide-field image reconstruction with double size
fftWF2 = zeros(2*NPixel,2*NPixel);                                          % Wide-field image spectrum with double size
fftWF = fftshift(fft2(WF));                                                 % Wide-field image raw
fftWF2(NPixel/2+1:NPixel/2+NPixel,NPixel/2+1:NPixel/2+NPixel)=fftWF;        % Wide-field image with double size
WF2 = real(ifft2(fftshift(fftWF2)));
out=WF2;
end