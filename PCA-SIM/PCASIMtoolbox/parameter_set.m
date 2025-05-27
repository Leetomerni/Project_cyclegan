function [out]=parameter_set(Iraw,Pixelsize,NA,lambda,mag)
    NPixel=size(Iraw,1);
    param.imgSize = NPixel;                                                % Image size
    param.micronsPerPixel = Pixelsize/mag;                                 % Pixel size 
    param.cyclesPerMicron = 1/(NPixel*param.micronsPerPixel);
    param.NA = NA;                                                         % Objective NA
    param.lambda = lambda*1000;                                            % Wavelength (nm)
    param.cutoff = 1000/(0.5*param.lambda/param.NA);                       % Cutoff frequency 2NA
    param.sampleLateral = ceil(param.cutoff/param.cyclesPerMicron)+1;      % Cutoff frequency radius
    param.nrBands = 2;                                                     % Bands
    param.phaOff=0;                                                        % Phase offset initial value
    param.fac=ones(1,param.nrBands);                                       % Modulation initial value 
    param.attStrength = 0;
    param.OtfProvider = SimOtfProvider(param,param.NA,param.lambda,1);     % Generate approximate OTF                  
    PSF = abs(otf2psf((param.OtfProvider.otf)));                           % Generate approximate PSF 
    param.OTF=param.OtfProvider.otf;
    param.psf=abs(otf2psf((param.OtfProvider.otf)));
    out=param;
end

function [ ret ] = SimOtfProvider( param,NA,lambda,a)
    ret.na=NA;
    ret.lambda=lambda;
    ret.cutoff=1000/(0.61*lambda/NA);
    %ret.cutoff=1000/(0.5*lambda/NA);
    ret.imgSize=param.imgSize;
    ret.cyclesPerMicron=param.cyclesPerMicron;
    ret.sampleLateral=ceil(ret.cutoff/ret.cyclesPerMicron)+1;
    ret.estimateAValue=a;
    ret.maxBand=2;
    ret.attStrength=param.attStrength;
    ret.attFWHM=1.0;
    ret.useAttenuation=1;
    ret=fromEstimate(ret);
    
    ret.otf=zeros(param.imgSize,param.imgSize);
    ret.otfatt=zeros(param.imgSize,param.imgSize);
    ret.onlyatt=zeros(param.imgSize,param.imgSize);

    ret.otf=otfToVector(ret.otf,ret,1,0,0,0,1);
    ret.onlyatt=getonlyatt(ret,0,0);
    ret.otfatt=ret.otf.*ret.onlyatt;
end

function [va]=valIdealOTF(dist)
    if dist<0 || dist>1
        va=0;
        return;
    end
    va=(1/pi)*(2*acos(dist)-sin(2*acos(dist)));
end

function [va]= valAttenuation(dist,str,fwhm)
    va=(1-str*(exp(-power(dist,2)/(power(0.5*fwhm,2)))).^1);
end

function [ret] = fromEstimate(ret)
    ret.isMultiband=0;
    ret.isEstimate=1;
    vals1=zeros(1,ret.sampleLateral);
    valsAtt=zeros(1,ret.sampleLateral);
    valsOnlyAtt=zeros(1,ret.sampleLateral);


    for I=1:ret.sampleLateral
        v=abs(I-1)/ret.sampleLateral;
        r1=valIdealOTF(v)*power(ret.estimateAValue,v);
        vals1(I)=r1;
    end

    for I=1:ret.sampleLateral
        dist=abs(I-1)*ret.cyclesPerMicron;
        valsOnlyAtt(I)=valAttenuation(dist,ret.attStrength,ret.attFWHM);
        valsAtt(I)=vals1(I)*valsOnlyAtt(I);
    end

    ret.vals=vals1;
    ret.valsAtt=valsAtt;
    ret.valsOnlyAtt=valsOnlyAtt;

end

function [ onlyatt ] = getonlyatt( ret,kx,ky )
    w=ret.imgSize;
    h=ret.imgSize;
    siz=[h w];
    cnt=siz/2+1;
    kx=kx+cnt(2);
    ky=ky+cnt(1);
    onlyatt=zeros(h,w);

    y=1:h;
    x=1:w;
    [x,y]=meshgrid(x,y);
    rad=hypot(y-ky,x-kx);
    cycl=rad.*ret.cyclesPerMicron;
    onlyatt=valAttenuation(cycl,ret.attStrength,ret.attFWHM);
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


