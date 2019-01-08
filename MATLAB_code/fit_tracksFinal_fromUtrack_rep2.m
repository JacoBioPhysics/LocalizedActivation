function [ tracks3D] = fit_tracksFinal_fromUtrack_rep2( tracks, images_path)
% 
%INPUT:
%
%tracks: is the output structure of utrack, named "tracksFinal" in the
%       workspace. It is a multidimensional structure with fields explained in the
%       "Readme" from the software developers. They can be filtered etc.
%
%OUTPUT:
%
%tracks3D: is a multidimensional structure containing the new 3D tracks. X
%and Y position are re-fitted using a maximum likelihood routine (MLEwG_elli_iter)
%and z position is inferred from the fitted widths of the asymmetrical
%fitted gaussians.

cd(images_path)

%***size of the data to be fed to the MLEwG_elli_iter function
datafitsize=5;
zdispl_thres=450;
%***filter out short and "double, merged, splitted" tracks from tracksFinal
%newtracks=filter_tracks(tracks,threshold, max);
tracks3D=tracks;
%***paths definitions and image and calibration data is here loaded

%images_path=uigetdir('D:');
%images_path='D:\ANALISYS\TRIALS\SAMPLE_MOVIE_FOR_UTRACK_902wstimulus_FC';
%BCK=im_read2photons('C:\Users\solari\Desktop\BCK_ANDOR_emGAIN100.tif',1,.17);

calibrationx=importdata('D:\ANALISYS\calibrationXlong.mat');
calibrationy=importdata('D:\ANALISYS\calibrationYlong.mat');
ascissacal=importdata('D:\ANALISYS\ascissacalLong.mat');

%***gets the name of the images in the target folder
%d=dir(strcat(images_path,'\*.tif'));
d=dir(strcat(images_path,'\*.Tiff'));


%***define some guess parameters for fitting
a=106;
PSFw=2.5;
widthguess=PSFw*a;
s2=[512,512];
%**************************************************************************
%-------------START THE ANALYSIS-------------------------------------------
%**************************************************************************
  
        clear images
        diffs=[];
        
        %***Cuts and stores the images according to datafisize. Images to be 
        %fitted are centered in the x-y coordinates obtained from utrack.
for j = 1:length(tracks)
        j
        amps=[];
        vars=[];
        bcks=[];
        XCoord=[];
        YCoord=[];                
        ZCoord=[];
        amp_bcks=[];        
        goodnesses=[];
        SXs=[];
        SYs=[];
        
        x=tracks(j).tracksCoordAmpCG(1:8:end);
        y=tracks(j).tracksCoordAmpCG(2:8:end);        
        first=tracks(j).seqOfEvents(1);
        last=tracks(j).seqOfEvents(2);

        try
            [data_frames,xrest,yrest,datafitsize]=cut_frames_rep(first,last,d,x,y,datafitsize);
        catch
            tracks3D(j)=[];
            'diocane'
            continue
        end
        
        
        for i=1:length(data_frames)
            
            if ~isnan(data_frames{i})
                data=data_frames{i};
                t=thresholdOtsu(data);
                ampp=sum(sum(data*(data>t)));
                p0=[(datafitsize+1)*a,(datafitsize+1)*a ,a*PSFw,a*PSFw,0,ampp];
                if i==1
                    pF_prev=p0;
                end
                sd=size(data);
                [xx,yy]=meshgrid(1:sd(1),1:sd(2));
                xx=xx*a;
                yy=yy*a;

                


                [pF, var, exitfl, fval]=MLEwG_elli_iter(medfilt2(data),p0,a,0,1,0);
                %[pF, var, exitfl, fval]=MLEwG_elli_iter(wiener2(data),pF,a,0,1,0); 
                %[pF, var, exitfl, fval]=MLEwG_elli_iter((data),pF,a,0,0,0); 


                pF(6)
               
                
                if  (fix(pF(2)/a)>datafitsize ... 
                            ||fix(pF(1)/a)>datafitsize ...
                                ||fix(pF(2)/a)<0 ... 
                                    ||fix(pF(1)/a)<0  ...
                                        ||(fix(abs(pF_prev(1)-abs(pF(1))))>zdispl_thres) ...
                                            ||(fix(abs(pF_prev(2)-abs(pF(2))))>zdispl_thres) ...
                                                || pF(6)<0 || pF(3)>400 || pF(4)>400 || pF(3)<0 || pF(4)<0 )
                    
                    %pF=[NaN,NaN,NaN,NaN,NaN,NaN]    
                    [pF, var, exitfl, fval]=MLEwG_elli_iter(wiener2(data),p0,a,0,1,0);
                    [pF, var, exitfl, fval]=MLEwG_elli_iter(data,pF,a,0,1,0);
                    
                    

                    if (pF(6)<10 || pF(6)>5000 || pF(3)>400 || pF(4)>400) || pF(3)<0 || pF(4)<0 ...
                                    fix(pF(2)/a)>datafitsize ... 
                                         ||fix(pF(1)/a)>datafitsize ...
                                             ||fix(pF(2)/a)<0 ... 
                                                ||fix(pF(1)/a)<0;
                                                
                        [pF, var, exitfl, fval]=MLEwG_elli_iter(medfilt2(data),pF_prev,a,0,1,0);
            
                    end
                    pF_prev=pF;

                end            
                
                [pF, var, exitfl, fval]=MLEwG_elli_iter((data),pF,a,0,1,0);
                pF_prev=pF;

                

    %compute z position
                if ~isnan(pF(3)) && ~isnan(pF(4))
                        
                        


                        [v,ind]=min(sqrt(  ((calibrationy)-pF(4)).^2 + ((calibrationx)-pF(3)).^2   ));
                        z=ascissacal(ind);                        

                        if i~=1 && (abs(z-z_prev)>zdispl_thres)
                            [pF, var, exitfl, fval]=MLEwG_elli_iter((data),pF_prev,a,0,1,0);
                            [v,ind]=min(sqrt(  ((calibrationy)-pF(4)).^2 + ((calibrationx)-pF(3)).^2   ));
                            z=ascissacal(ind); 
                        end

                        if i~=1 && (abs(z-z_prev)>zdispl_thres)
                            z=z_prev;
                        end
                        z_prev=z;

                else
                        z=NaN;
                end            
                
                %ampbck=pF(7);
                amps=vertcat(amps,pF(6));
                XCoord=vertcat(XCoord,x(i)*a);
                YCoord=vertcat(YCoord,y(i)*a);
                ZCoord=vertcat(ZCoord,z);
                vars=vertcat(vars,var);
                bcks=vertcat(bcks,pF(5));
                SXs=vertcat(SXs,pF(3));
                SYs=vertcat(SYs,pF(4));
                %amp_bcks=vertcat(amp_bcks,ampbck);
                
                diffs=horzcat(diffs,sqrt((pF(1)-(datafitsize)*a )^2 + (pF(2)-(datafitsize)*a)^2 ));
                try
                    goodnesses=horzcat(goodnesses,sum(sum(double(data)-double(expectedgaussian (xx,yy,pF,a)))));
                catch
                    goodnesses=horzcat(goodnesses,0);
                end
            else %if isnan(data_frams{i})
                pF=[NaN,NaN,NaN,NaN,NaN,NaN];
            end

      
        end % for i
        
        tracks3D(j).framename = [images_path '/'  d(j).name];
        tracks3D(j).AMPs=amps;
        tracks3D(j).LineCoord=YCoord;
        tracks3D(j).ColoumnCoord=XCoord;
        tracks3D(j).ZCoord=ZCoord;        
        tracks3D(j).VARs=vars;
        tracks3D(j).bcks=bcks;
        tracks3D(j).amp_bcks=amp_bcks;
        tracks3D(j).diffs=diffs;
        tracks3D(j).goodnesses=goodnesses;
        tracks3D(j).SXs=SXs;
        tracks3D(j).SYs=SYs;
        

        
    
end %for j
    
    
end