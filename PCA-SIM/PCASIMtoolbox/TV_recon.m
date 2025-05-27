function [out]=TV_recon(Iraw,param,par,sub_optimization)
for i=1:9
IIrawFFT(:,:,i)=fftshift(fft2(Iraw(:,:,i)));
end
NPixel=size(Iraw,1);
Coarse_frequency=zeros(NPixel*2, NPixel*2);
Fine_frequency=zeros(NPixel*2, NPixel*2);
Shifted=zeros(NPixel*2, NPixel*2,3);    
siz=size(Iraw(:,:,1));
param.attStrength=0.9;
param.a=1;                
param.attFWHM=1.0;
for I=1:3                          
    param.fac=[1,1];
    if(sub_optimization)
        Separate=separateBands(IIrawFFT(:,:,(I-1)*3+...
        1:I*3),par(I).phaOff,param.nrBands,param.fac);          
        dist=0.15;
        overlap=dist;
        divideByOtf=true;   
        [b0,b1]=commonRegion(Separate(:,:,1),Separate(:,:,2),1,2,...
         param.OtfProvider,-par(I).px,-par(I).py,dist,overlap,divideByOtf);  
        b0=ifft2(fftshift(b0));       
        b1=ifft2(fftshift(b1));       
        scal=1/sum(sum(real(b0).^2+imag(b0).^2));    
        %b1s=fourierShift(b1,-par.px,-par.py);  
        siz=size(b1);
        [x,y]=meshgrid(0:siz(2)-1,0:siz(1)-1);
        x=x/siz(2);
        y=y/siz(1);
        b1s=b1.*exp(2*pi*1i*(-par(I).py*y+(-par(I).px)*x));
        b1s=b1s.*conj(b0);              
        corr=sum(sum(b1s))*scal;  
        par(I).phaOff=par(I).phaOff-phase(corr);
        resMag=abs(corr);        
        resMag(resMag<0.35) = 0.7;        
        param.fac(2:param.nrBands)=resMag;        
    end
    
    Separate=separateBands(IIrawFFT(:,:,(I-1)*3+...            % Spectrum separation
        1:I*3),par(I).phaOff,param.nrBands,param.fac);   
    
    Shifted=zeros(NPixel*2,NPixel*2,3);                        % Separated spectrum with double size      
    Shifted(:,:,1)=placeFreq(Separate(:,:,1));                              
    Shifted(:,:,2)=placeFreq(Separate(:,:,2));   
    Shifted(:,:,3)=placeFreq(Separate(:,:,3));   
    Shifted(:,:,2)=NfourierShift(Shifted(:,:,2),-(2-1)*par(I).px,...           % Spectrum shift
        -(2-1)*par(I).py); 
    Shifted(:,:,3)=NfourierShift(Shifted(:,:,3), (2-1)*par(I).px,...           % Spectrum shift
         (2-1)*par(I).py);    
    Coarse_frequency=Coarse_frequency+Shifted(:,:,1)+Shifted(:,:,2)...      % Spectrum merging
        +Shifted(:,:,3);
    
    Shifted(:,:,1)=otfToVector(Shifted(:,:,1),param.OtfProvider,1,0,0,1,0);    % Deconvolution      
    Shifted(:,:,2)=otfToVector(Shifted(:,:,2),param.OtfProvider,2,...          % Deconvolution
        -(2-1)*par(I).px,-(2-1)*par(I).py,1,0);    
    Shifted(:,:,3)=otfToVector(Shifted(:,:,3),param.OtfProvider,2,...          % Deconvolution
         (2-1)*par(I).px,(2-1)*par(I).py,1,0);      
    Fine_frequency=Fine_frequency+Shifted(:,:,1)+Shifted(:,:,2)+...         % Deconvolved spectrum merging
        Shifted(:,:,3);
end

w=0.1; 
K0=[par(1).K0,par(2).K0,par(3).K0];
K=max([ceil(K0)]);                             
cutoff=floor(1*K)/param.sampleLateral+1.0;     
otfHiFi=zeros(NPixel*2,NPixel*2);
otfHiFi=writeApoVector(otfHiFi,param.OtfProvider,cutoff);   
Mask=zeros(NPixel*2,NPixel*2);
Mask(otfHiFi~=0)=1;
param.Dir=par;
wFilter0=WienerFilterWiener_2D(param);                                               
Wk0=otfHiFi./(wFilter0.wDenom+w^2);
fftWiener=real(ifft2(fftshift((Fine_frequency.*Wk0.*Mask))));
Temp=fftWiener;
Temp(Temp<0)=0;
Wiener=255*Temp/max(max(Temp));

iter_Bregman = 1e2+50;                                                      % number of iteration
siranu=150;              
zbei=0;
lamda = 1;
gexiang=1;
y = Wiener;                                                         % observed data 				
tic
ymax=max(y(:));
y=y./ymax;              

%% initialization
[sx,sy,sz] = size(y);
sizex=[sx,sy,sz] ;
x = zeros(sizex);                                                           % start point
ztidu(:,:,1)=1;                    
ztidu(:,:,2)=-1;

%FFT of difference operator 
xfft=fftn([1 -1],sizex).*conj(fftn([1 -1],sizex));
yfft=fftn([1 -1]',sizex).*conj(fftn([1 -1]',sizex));
zfft=fftn(ztidu,sizex).*conj(fftn(ztidu,sizex));

%% iteration
for TH=1:1
    clear gap
    b7 = zeros(sizex);
    b8 = zeros(sizex);
    b9 = zeros(sizex);
    d7 = zeros(sizex);
    d8 = zeros(sizex);
    d9 = zeros(sizex);
    x = zeros(sizex);
    for ii = 1:iter_Bregman
%% renew x
        frac = (siranu/lamda)*(y);               
        frac = frac-back_diff(d7-b7,1,2);
        frac = frac-back_diff(d8-b8,1,1);
        frac = frac-(zbei)*back_diff(d9-b9,1,3);
        frac = fftn(frac);
        divide = (siranu/lamda);
        divide = divide + xfft;
        divide = divide + yfft;
        divide = divide + (zbei)*zfft;
        x = real(ifftn(frac./divide));

%% calculate the dirivative of x
        u_x = forward_diff(x,1,2);
        u_y = forward_diff(x,1,1);
        u_z = forward_diff(x,1,3);
%% renew d
    if gexiang==1
         signd7 = abs(u_x+b7)-1/lamda;
        signd7(signd7<0)=0;
        signd7=signd7.*sign(u_x+b7);
        d7=signd7;
        
        signd8 = abs(u_y+b8)-1/lamda;
        signd8(signd8<0)=0;
        signd8=signd8.*sign(u_y+b8);
        d8=signd8;
        
        signd9 = abs(u_z+b9)-1/lamda;
        signd9(signd9<0)=0;
        signd9=signd9.*sign(u_z+b9);
        d9=signd9;
    else
        s_T = sqrt((u_x+b7).^2+(u_y+b8).^2+(u_z+b9).^2*(zbei));
        d7 = lamda*s_T.*(u_x+b7)./(lamda*s_T+1+eps);
        d8 = lamda*s_T.*(u_y+b8)./(lamda*s_T+1+eps);
        d9 = lamda*s_T.*(u_z+b9)./(lamda*s_T+1+eps);
    end
    
%% renew frac
        b7 = b7+(u_x-d7);
        b8 = b8+(u_y-d8);
        b9 = b9+(u_z-d9);
        x(x<0) = 0;
    end
end

TV=x.*ymax;

out=TV;

end

function [ separate ] = separateBands( IrawFFT,phaOff,bands,fac )
phaPerBand=(bands*2)-1;      
phases=zeros(1,phaPerBand);
for p=1:phaPerBand
    phases(p)=(2*pi*(p-1))/phaPerBand+phaOff;  
end
separate=separateBands_final(IrawFFT,phases,bands,fac);
end

function [separate] = separateBands_final(IrawFFT,phases,bands,fac)
for I=2:bands  
    fac(I)=fac(I)*0.5;  
end
comp=zeros(1,bands*2-1);
comp(1)=0;
for I=2:bands
    comp((I-1)*2)=I-1;
    comp((I-1)*2+1)=-(I-1);
end
compfac=zeros(1,bands*2-1);
compfac(1)=fac(1);
for I=2:bands
    compfac((I-1)*2)=fac(I);
    compfac((I-1)*2+1)=fac(I);
end

W=exp(1i*phases'*comp);
for I=1:bands*2-1
    W(I,:)=W(I,:).*compfac;
end

length=size(phases,2);
siz=size(IrawFFT(:,:,1));
Sk=zeros(siz(1),siz(2),bands*2-1);

S=reshape(IrawFFT,[prod(siz),length])*pinv(W)';
Sk=reshape(S,[siz,bands*2-1]);
separate=Sk;
end

function [ newb0,newb1 ] = commonRegion( band0, band1, bn0, bn1, otf, kx, ky, dist,weightLimit, divideByOtf)

siz=size(band0);
w=siz(2);
h=siz(1);
cnt=[siz(1)/2+1,siz(2)/2+1];   

weight0=zeros(h,w);
weight1=zeros(h,w);
wt0=zeros(h,w);
wt1=zeros(h,w);

if bn0==1
    weight0=otfToVector(weight0,otf,bn0,0,0,0,1);
    weight1=otfToVector(weight1,otf,bn1,0,0,0,1);
else
    weight0=otfToVector(weight0,otf,bn0,0,0,0,1);
    weight1=otfToVector(weight1,otf,bn1,0,0,0,1);
end
wt0=otfToVector(wt0,otf,bn0, kx, ky,0,1);   
wt1=otfToVector(wt1,otf,bn1,-kx,-ky,0,1);     

x=1:w;
y=1:h;
[x,y]=meshgrid(x,y);
rad=sqrt((y-cnt(1)).^2+(x-cnt(2)).^2);
max=sqrt(kx*kx+ky*ky);
ratio=rad./max;                                              

mask1=(abs(weight0)<weightLimit) | (abs(wt0)<weightLimit);   
cutCount=length(find(mask1~=0));
band0(mask1)=0;                               
mask2=abs(weight1)<weightLimit | abs(wt1)<weightLimit;   
band1(mask2)=0;                                         

weight0(weight0==0)=1;
weight1(weight1==0)=1;
wt0(wt0==0)=1;
wt1(wt1==0)=1;

if divideByOtf==true
    band0=band0./weight0;              
    band1=band1./weight1;             
end


mask=ratio<dist | ratio>(1-dist);

band0(mask)=0;


idx = repmat({':'}, ndims(mask), 1);
n = size(mask, 1 );
if kx>0
    idx{2}=[round(kx)+1:n 1:round(kx)];
else
    idx{2}=[n-round(abs(kx))+1:n 1:n-round(abs(kx))];
end
if ky>0
    idx{1}=[round(ky)+1:n 1:round(ky)];
else
    idx{1}=[n-round(abs(ky))+1:n 1:n-round(abs(ky))];
end
mask0=mask(idx{:});

band1(mask0)=0;
newb0=band0;
newb1=band1;
end

function [ ret ] = otfToVector( vec,otf,band,kx,ky,useAtt,write )
    siz=size(vec);
    w=siz(2);
    h=siz(1);
    cnt=siz/2+1;
    kx=kx+cnt(2);
    ky=ky+cnt(1);

    x=1:w;
    y=1:h;  
    [x,y]=meshgrid(x,y);

    rad=hypot(y-ky,x-kx);
    cycl=rad.*otf.cyclesPerMicron;   

    mask=cycl>otf.cutoff;
    cycl(mask)=0;

    val=getOtfVal(otf,band,cycl,useAtt);
    if write==0                      
        vec=vec.*val;
    else
        vec=val;
    end
    vec(mask)=0;
    ret=vec;

end

function [ val ]= getOtfVal(otf,band,cycl,att)
    pos=cycl./otf.cyclesPerMicron;
    cpos=pos+1;

    lpos=floor(cpos);
    hpos=ceil(cpos);
    f=cpos-lpos;

    if att==1
        retl=otf.valsAtt(lpos).*(1-f);
        reth=otf.valsAtt(hpos).*f; 
        val=retl+reth;
    else
        retl=otf.vals(lpos).*(1-f);
        reth=otf.vals(hpos).*f;
        val=retl+reth;
    end

    mask=ceil(cpos)>otf.sampleLateral;
    val(mask)=0;

end

function [out] = placeFreq( in )
    siz=size(in);
    w=siz(2);
    h=siz(1);
    out=zeros(2*w,2*h);
    out(h/2+1:h+h/2,w/2+1:w+w/2)=in;
end

function [outv] = NfourierShift( inv,kx,ky )
    inv=(ifft2(fftshift(inv)));
    siz=size(inv);
    [x,y]=meshgrid(0:siz(2)-1,0:siz(1)-1);
    x=x/siz(2);
    y=y/siz(1);
    outv=inv.*exp(2*pi*1i*(ky*y+kx*x));
    outv=fftshift(fft2((outv)));
end

function [ ret ] = writeApoVector( vec,otf,cutoff)
    [h,w]=size(vec); 
    cnt=[h/2+1,w/2+1];
    for y=1:h
        for x=1:w
            rad=hypot(y-cnt(1),x-cnt(2));
            cycl=rad*otf.cyclesPerMicron;
            frac=cycl/(otf.cutoff*cutoff);
            if frac<0 || frac>1
                valIdealotf=0;
            else
                valIdealotf=(1/pi)*(2*acos(frac)-sin(2*acos(frac)));
            end
            vec(y,x)=valIdealotf;
        end
    end
    ret=vec;
end

function [wFilter]=WienerFilterWiener_2D( sp )
    wFilter.sp=sp;
    wFilter.wDenom=updateCache(sp);
end

function [ wDenom ] = updateCache( sp )
    Temp=sp.OtfProvider.otf;
    siz=size(Temp(:,:,1));
    w=siz(2);
    h=siz(1);
    wDenom=zeros(2*h,2*w);
    for d=1:3
        for b=1:2
            wd=wDenom;
            [wDenom]=addWienerDenominatorWiener_2D(wd,sp,d,b);
        end
    end
end

function [wDenom]=addWienerDenominatorWiener_2D( wd,sp,d,b)
    siz=size(wd);
    w=siz(2);
    h=siz(1);
    dir=sp.Dir(d);
    cyclMicron=sp.cyclesPerMicron;
    cnt=[h/2+1,w/2+1];

    x=1:w;
    y=1:h;
    [x,y]=meshgrid(x,y);

    rad1=hypot(x-cnt(2)-(b-1)*dir.px,y-cnt(1)-(b-1)*dir.py)*cyclMicron;
    rad2=hypot(x-cnt(2)+(b-1)*dir.px,y-cnt(1)+(b-1)*dir.py)*cyclMicron;

    otfVal1=abs(getOtfVval(sp.OtfProvider,rad1)).^2;   
    otfVal2=abs(getOtfVval(sp.OtfProvider,rad2)).^2;

    if b==1
        if sp.OtfProvider.attStrength==0
            otfVal1=otfVal1/2;
            otfVal2=otfVal2/2;
        else
            otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.00,1.0*sp.OtfProvider.attFWHM)/2.0;
            otfVal2=otfVal2.*valAttenuation(rad2,sp.OtfProvider.attStrength/1.00,1.0*sp.OtfProvider.attFWHM)/2.0;
        end
    else
            otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1*sp.OtfProvider.attFWHM)/1.0;
            otfVal2=otfVal2.*valAttenuation(rad2,sp.OtfProvider.attStrength/1.0,1*sp.OtfProvider.attFWHM)/1.0;
    end
    wd=wd+otfVal1+otfVal2;
    wDenom=wd;
end

function [ val ]= getOtfVval(ret,cycl)
    mask=cycl>ret.cutoff;
    cycl(mask)=0; 

    pos=cycl./ret.cyclesPerMicron;
    cpos=pos+1;
    lpos=floor(cpos);
    hpos=ceil(cpos);
    f=(cpos-lpos); 

    retl=ret.vals(lpos).*(1-f);
    reth=ret.vals(hpos).*f;
    val=retl+reth;
    val(mask)=0;
end

function [va]= valAttenuation(dist,str,fwhm)
    va=1-str*(exp(-power(dist,2)/(power(0.5*fwhm,2))));
end

function out = back_diff(data,step,dim)
    [m n r]=size(data);
    SIZE=[m,n,r];
    position = ones(1,3);
    temp1 = zeros(SIZE+1);
    temp2 = zeros(SIZE+1);
    temp1(position(1):SIZE(1),position(2):SIZE(2),position(3):SIZE(3))=data;
    temp2(position(1):SIZE(1),position(2):SIZE(2),position(3):SIZE(3))=data;
    SIZE(dim)=SIZE(dim)+1;    
    position(dim)=position(dim)+1;
    temp2(position(1):SIZE(1),position(2):SIZE(2),position(3):SIZE(3))=data;
    temp1 = (temp1-temp2)/step;
    SIZE(dim)=SIZE(dim)-1;    
    out = temp1(1:SIZE(1),1:SIZE(2),1:SIZE(3));
end

function  out = forward_diff(data,step,dim) 
    [m n r]=size(data);
    SIZE=[m,n,r];
    position = ones(1,3);
    temp1 = zeros(SIZE+1);
    temp2 = zeros(SIZE+1);
    SIZE(dim)=SIZE(dim)+1;    
    position(dim)=position(dim)+1;
    temp1(position(1):SIZE(1),position(2):SIZE(2),position(3):SIZE(3))=data;
    temp2(position(1):SIZE(1),position(2):SIZE(2),position(3):SIZE(3))=data;
    SIZE(dim)=SIZE(dim)-1;   
    temp2(1:SIZE(1),1:SIZE(2),1:SIZE(3))=data;
    temp1 = (temp1-temp2)/step;
    SIZE(dim)=SIZE(dim)+1;   
    out = temp1(position(1):SIZE(1),position(2):SIZE(2),position(3):SIZE(3));    
    out = -out;
end