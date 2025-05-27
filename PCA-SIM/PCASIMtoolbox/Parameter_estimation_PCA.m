function [out1,out2]=Parameter_estimation_PCA( IIraw,param,Filter_size,Mask_size)
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

PCA_num=1;         

tic
for angle_num=1:3
    
    param.fac=[1 1];
    param.phaOff=0;
     Spectrum=separateBands(IIrawFFT(:,:,(angle_num-1)*3+1:angle_num*3),...
         param.phaOff,param.nrBands,param.fac); 
     temp = Spectrum(:,:,2).*NotchFilter;
     [yPos,xPos] = find(temp==max(max(temp)));       

    peak.xPos = xPos;
    peak.yPos = yPos;                                                     
    
for pca_size=1:PCA_num
  if(pca_size==1)
      kx = (peak.xPos-Center(2)) ;                                         
      ky = (peak.yPos-Center(1));                                          
      old_kx=kx ;   
      old_ky=ky;   
  else
      kx = old_kx ;   
      ky = old_ky;           
  end
  
  Temp=NfourierShift(placeFreq(Spectrum(:,:,2)),-(2-1)*old_kx,-(2-1)*...
      old_ky);
  MASK = ones(size(Temp));
  MASK(NPixel-Mask_size+1:NPixel+Mask_size+1, NPixel-Mask_size+1:NPixel...
      +Mask_size+1) = 0;
  Temp(MASK==1) = 0;
  ROI = Temp(NPixel-Filter_size+1:NPixel+Filter_size+1, NPixel-...
      Filter_size+1:NPixel+Filter_size+1);
    
  SIZE = 2*Filter_size+1;
  Space_ROI=(ifft2(ifftshift(ROI)));
  Phase_ROI=exp(1i*angle(Space_ROI));

  [U,S,V] = svd(Phase_ROI);                                                  % Singular value decomposition
  Y = U*S*V';
  SS=zeros(SIZE,SIZE);
  SS(1,1)=S(1,1);
  Phase_ROI= U*SS*V';                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

  UNV1=unwrap(angle(V(:,1)));                                                % Phase unwrapping
  UNV_T=UNV1';
  UNV_T_fit=UNV_T(1,2:SIZE);

  FIT=polyfit(1:SIZE-1,UNV_T_fit,1);                                         % Linear fitting
  dx = - FIT(1)*SIZE/2/pi;
  NewUNV1=polyval(FIT,1:SIZE); 
  NEWV1=exp(1i.*((NewUNV1')));
  V(:,1)=NEWV1;

  UNU1=unwrap(angle(U(:,1)));
  UNN_T=UNU1';
  UNN_T_fit=UNN_T(1,2:SIZE);
  FIT=polyfit(1:SIZE-1,UNN_T_fit,1);
  dy = FIT(1)*SIZE/2/pi;
  NewUN1=polyval(FIT,1:SIZE);
  NEWU1=exp(1i.*(NewUN1'));
  U(:,1)=NEWU1;

  NEW_Phase_ROI= U*SS*V';

  old_kx= old_kx + dx ;   
  old_ky= old_ky + dy;  
  param.Dir(angle_num).px  = old_kx;
  param.Dir(angle_num).py  = old_ky;
  param.Dir(angle_num).phaOff = - angle(NEW_Phase_ROI(Filter_size+1,...
        Filter_size+1));
  K0(angle_num) = sqrt((param.Dir(angle_num).px)^2+(param.Dir...
        (angle_num).py)^2);  
  param.Dir(angle_num).K0=K0(angle_num);
end
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
