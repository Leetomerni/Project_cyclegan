function [out]=RLdeconv(Iraw,PSF,judge)
if judge==true
    IIraw=deconvlucy(Iraw,PSF,5);
else
    IIraw=Iraw;
end
out=IIraw;
end