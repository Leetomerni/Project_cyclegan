function [out1,out2]=Parameter_estimation_COR( IIraw,param)
for i=1:9
IIrawFFT(:,:,i)=fftshift(fft2(IIraw(:,:,i)));
end
NPixel=size(IIrawFFT(:,:,1),1);
[x,y] = meshgrid(1:NPixel,1:NPixel);                                        % Image coordinate
Center=[NPixel/2+1,NPixel/2+1];                                             % Image ceteral coordinate
R=sqrt((y-Center(1)).^2+(x-Center(2)).^2);                                  % frequency radius
Mask_2NA=double(R<=1.0*(param.cutoff/param.cyclesPerMicron+1));             % Mask 2NA
NotchFilter0=getotfAtt(NPixel,param.cyclesPerMicron,0.5*param.cutoff,0,0);  % Filter for zero frequency
NotchFilter=NotchFilter0.*Mask_2NA;                                         % Filter for zero frequency

for angle_num = 1:3
    Spectrum=separateBands(IIrawFFT(:,:,(angle_num-1)*3+...
        1:angle_num*3),param.phaOff,param.nrBands,param.fac);               % Spectrum separated
    Spectra{1,angle_num}=Spectrum;                                          % Spectra all    
temp = Spectrum(:,:,2).*NotchFilter;
[yPos,xPos] = find(temp==max(max(temp))); 
end

tic
for angle_num=1:3
    Spectrum = Spectra{1,angle_num};                                        % Get spectrum;
    temp = Spectrum(:,:,2).*NotchFilter;
    [yPos,xPos] = find(temp==max(max(temp)));                               % peak localization
    peak.xPos = xPos(1);
    peak.yPos = yPos(1);
    
    cntrl = zeros(10,30);                                                   % Subregion around peak
    overlap = 0.15;                                                         %
    step = 2.5;                                                             % matching step
    kx = (peak.xPos-Center(2));                                             % peak shift x (integral)
    ky = (peak.yPos-Center(1));                                             % peak shift y (integral)
    
    [peak,cntrl] = fitPeak(Spectrum(:,:,1)/(max(max(abs(Spectrum(:,:,1)...
        )))),Spectrum(:,:,2)/(max(max(abs(Spectrum(:,:,2))))),1,2,...
        param.OtfProvider,-kx,-ky,overlap,step,cntrl);                      % fitting subpixel

    p1=getPeak(Spectrum(:,:,1),Spectrum(:,:,2),...
        1,2,param.OtfProvider,peak.kx,peak.ky,overlap);     
    
    param.Dir(angle_num).px = -peak.kx;                                     % new peak kx
    param.Dir(angle_num).py = -peak.ky;                                     % new peak ky
    param.Dir(angle_num).phaOff = -phase(p1);                               % new phase offset   
    K0(angle_num) = sqrt((param.Dir(angle_num).px)^2+...
        (param.Dir(angle_num).py)^2);                                       % frequency vector
    param.Dir(angle_num).K0=K0(angle_num);
    
end
toc
  out1=param.Dir;
  out2=K0(angle_num);
end

function [ otfAtt ] = getotfAtt( imgSize,cyclesPerMicron,attfwhm,kx,ky )
    w=imgSize;
    h=imgSize;
    siz=[h w];
    cnt=siz/2+1;           
    kx=kx+cnt(2);      
    ky=ky+cnt(1);
    otfAtt=zeros(h,w);     
    y=1:h;
    x=1:w;
    [x,y]=meshgrid(x,y);
    rad=hypot(y-ky,x-kx);  
    cycl=rad.*cyclesPerMicron; 
    otfAtt=(1-exp(-power(cycl,4)/(2*power(attfwhm,4))));  
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

function [ peak ,cntrl] = fitPeak( band0, band1, bn0, bn1, otf, kx, ky, weightLimit, search, cntrl)



resPhase=0;
resMag=0;
for iter=1:3
    b0=band0;
    b1=band1;
    dist=0.15;
    divideByOtf=true;
    [b0,b1]=commonRegion(b0,b1,bn0,bn1,otf,kx,ky,dist,weightLimit,divideByOtf);  
    b0=ifft2(ifftshift(b0));                                                
    b1=ifft2(ifftshift(b1));                                                     
    corr=zeros(10,10);
    scal=1/sum(sum(real(b0).^2+imag(b0).^2));   
    cmax=0;                                     
    cmin=inf;                                   
    newKx=0;                                   
    newKy=0;                                   
    
    tkx=kx;                                    
    tky=ky;                                     
    ts=search;                                                              
    for yi=1:10
        for xi=1:10
            
            xpos=tkx+(((xi-1)-4.5)/4.5)*ts;        
            ypos=tky+(((yi-1)-4.5)/4.5)*ts;
            b1s=b1;  
            b1s=fourierShift(b1,xpos,ypos);  
            b1s=b1s.*conj(b0);               
            corr(xi,yi)=sum(sum(b1s))*scal;  
        end
    end
    
    for yi=1:10
        for xi=1:10
            if abs(corr(xi,yi))>cmax
                cmax=abs(corr(xi,yi));           
                newKx=tkx+(((xi-1)-4.5)/4.5)*ts; 
                newKy=tky+(((yi-1)-4.5)/4.5)*ts; 
                resPhase=phase(corr(xi,yi));     
                resMag=abs(corr(xi,yi));         
            end
            
            if abs(corr(xi,yi)<cmin)
                cmin=abs(corr(xi,yi));
            end
        end
    end
    
    if ~isempty(cntrl)
        for yi=1:10
            for xi=1:10
                cntrl(yi,xi+(iter-1)*10)=(abs(corr(xi,yi))-cmin)/(cmax-cmin);  
            end
        end
    end
    
    kx=newKx;
    ky=newKy;
    
    search=search/3;  
end

peak.kx=kx;
peak.ky=ky;
peak.resPhase=resPhase;
peak.resMag=resMag;
end

function [ comshift ] = fourierShift( vec, kx, ky  )
    siz=size(vec);
    [x,y]=meshgrid(0:siz(2)-1,0:siz(1)-1);
    x=x/siz(2);
    y=y/siz(1);
    comshift=vec.*exp(2*pi*1i*(ky*y+kx*x));
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

function [ ret] = getPeak( band0,band1,bn0,bn1,otf,kx,ky,weightLimit )
    b0=band0;
    b1=band1;
    dist=0.15;
    [b0,b1]=commonRegion(band0,band1,bn0,bn1,otf,kx,ky,dist,weightLimit,true); 

    b0=ifft2(ifftshift(b0));
    b1=ifft2(ifftshift(b1));

    b1=fourierShift(b1,kx,ky);
    b1=b1.*conj(b0);

    scal=1/sum(sum(real(b0).^2+imag(b0).^2));  
    ret=sum(sum(b1))*scal;  
end