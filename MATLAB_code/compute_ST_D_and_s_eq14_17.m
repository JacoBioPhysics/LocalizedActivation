function [D,s,SNR, var, varnormal] = compute_ST_D_and_s_alt(track,pixelsize,varsigma)

    %this function computes the single track D (apparent diff. coeff.) and s
    %(localization error) with the method presented in 
    %Vestergaard et al., PRE (2014)

    if nargin <2
        pixelsize=0.106;
    end
    
    %compute the S(t) function
    framelength = 0.0650; % in s
    illumination = 0.035; %in s
    illumination_fraction=illumination/framelength;

    s=zeros(1,1000);
    dt=framelength/double(length(s));
    s(1:round(illumination_fraction*1000))=1;

    N=sum(s);
    s=s/(N*dt);  
    ss=sum(s);

    %compute R
    dt=0.065/double(length(s));
    R=0;
    for i =1:length(s)
        S_t=compute_S(i,s,dt);
        R = R + (S_t.*(1-S_t)*dt);
    end
    R= (1/framelength) *R;
    %R=1/12.;
    %compute D
    
    if ~isfield(track,'x')
       track.x=track.tracksCoordAmpCG(1:8:end);
       track.y=track.tracksCoordAmpCG(2:8:end);
    end
    
    L=length(track.x);
    
    x1=track.x(1:end-1)*pixelsize;
    x2=track.x(2:end)*pixelsize;
    
    y1=track.y(1:end-1)*pixelsize;
    y2=track.y(2:end)*pixelsize;
    
    dx_X=( (x2-x1) );
    dx_Y=(  (y1-y2) );
        
    x3=track.x(3:end)*pixelsize;
    y3=track.y(3:end)*pixelsize;
    
    dx2_X=( (x3-x2(1:end-1)) );
    dx2_Y=( (y3-y2(1:end-1)));
    
    D_X= ( nanmean(dx_X.^2)/(2*framelength)+nanmean(dx_X(1:end-1).*dx2_X)/framelength);
    D_Y= ( nanmean(dx_Y.^2)/(2*framelength)+nanmean(dx_Y(1:end-1).*dx2_Y)/framelength);
 
    
    s_x=R*nanmean(dx_X.^2) + (2*R-1)*nanmean(dx_X(1:end-1).*dx2_X);
    s_y=R*nanmean(dx_Y.^2) + (2*R-1)*nanmean(dx_Y(1:end-1).*dx2_Y);

    s=0.5*(s_x+s_y);
    
    sigma=0.023^2;
    sigma=s;
    
    D=0.5*(D_X+D_Y);

    %varsigma=9.4317e-07;
    
    SNR=sqrt(2*D*L*framelength/sigma);
    
    if ~(imag(SNR)==0)
        SNR=-1;
    end
    
    e=sigma/(D*framelength)-2*R;
    
    N=length(track.x);
    
    %var=D^2*(  (6 + 4*e +2*e^2)/(N) +  (4*(1+e)^2)/(N^2) ) ;
    
    varnormal=D^2* ( (6 + 4*e + 2*e^2)/N + (4*(1 + e)^2)/N^2 );
    
    var=sqrt(varnormal)/D;
    sqrtvar=0;
    
end%2-0.5