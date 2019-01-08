function [ tracks3D, indbadtracks ] = fit_tracksFinal_fromUtrack_rep( tracks, images_path, threshold , save_track_path )
%get 3D tracks from utrack 2.2.1" output structure 
% 
%INPUT:
%
%tracks: is the output structure of utrack, named "tracksFinal" in the
%workspace. It is a multidimensional structure with fields explained in the
%"Readme" from the software developers
%
%OUTPUT:
%
%tracks3D: is a multidimensional structure containing the new 3D tracks. X
%and Y position are re-fitted using a maximum likelihood routine (MLEwG_elli_iter)
%and z position is inferred from the fitted widths of the asymmetrical
%fitted gaussians.

%***size of the data to be fed to the MLEwG_elli_iter function
datafitsize=10;
maximum=100;
%***filter out short and "double, merged, splitted" tracks from tracksFinal
newtracks=filter_tracks(tracks,threshold, maximum);
%newtracks=tracks;
%***paths definitions and image and calibration data is here loaded

%images_path=uigetdir('D:');
%images_path='D:\ANALISYS\TRIALS\SAMPLE_MOVIE_FOR_UTRACK_902wstimulus_FC';
BCK=im_read2photons('C:\Users\solari\Desktop\BCK_ANDOR_emGAIN100.tif',1,.17);

calibrationx=importdata('D:\ANALISYS\calibrationXlong.mat');
calibrationy=importdata('D:\ANALISYS\calibrationYlong.mat');
ascissacal=importdata('D:\ANALISYS\ascissacalLong.mat');

%***gets the name of the images in the target folder
d=dir(strcat(images_path,'\*.tif'));


%***define some guess parameters for fitting
a=64;
widthguess=165;

%**************************************************************************
%-------------START THE ANALYSIS-------------------------------------------
%**************************************************************************
    for i=1:length(newtracks)
        clear images
        n=0;
        %***Cuts and stores the images according to datafisize. Images to be 
        %fitted are centered in the x-y coordinates obtained from utrack.
        
        framesN=[newtracks(i).seqOfEvents(1,1),newtracks(i).seqOfEvents(end,1)];%frames numbers of the track
        x=(newtracks(i).tracksCoordAmpCG(1,1:8:end));
        y=(newtracks(i).tracksCoordAmpCG(1,2:8:end));
        x=fix(x);
        y=fix(y);
        images=cut_frames_from_Utrack2(framesN(1),framesN(2),d,x,y,BCK,datafitsize);%function to cut and store the sub-images
        if sum(sum(images))==0
            tracks3D(i).x=nan;
            tracks3D(i).y=nan;
            tracks3D(i).z=nan;
            continue
        end
        
        
        
        
        %if the cutting function was not able to cut the imaghes, just
        %discard the track. remove Nan later.
        
        if images==0
            tracks3D(i).x=nan;
            tracks3D(i).x=nan;           
            tracks3D(i).y=nan;
            tracks3D(i).SXs(i)=nan;
            tracks3D(i).SYs(i)=nan;
            tracks3D(i).bck(i)=nan;
            tracks3D(i).AMP(i)=nan;
            
            continue
        end
        
        s=size(images);
        bck=mean(mean(BCK(100:130,100:130)));
        
        for j=1:s(3) %start the analysis od the track, cutting frames from images, conuting guesseew and applying ML fit

            j
            imagesc(localAvgStd2D(  images(:,:,j),3))
            std_im=localAvgStd2D(  images(:,:,j),3); 
            ampguess=sum(sum(images(:,:,j)));
            %***initialize the guess parameters for each ellipse to be fitted
            t=max(max(localAvgStd2D(images(:,:,j),3)))*0.6;
            ROI_mask=(localAvgStd2D(images(:,:,j),5)>t);
            %ROI_mask=bwmorph(ROI_mask,'open',Inf);
            ROI_mask=bwmorph(ROI_mask,'clean',Inf);    
            
            counter=0;
            
            while counter<5 && ~any(any(ROI_mask)) 
                counter=counter+1;
                t=max(max(localAvgStd2D(images(:,:,j),3)))*(0.65-0.05*counter);
                ROI_mask=(localAvgStd2D(images(:,:,j),5)>t);               
                ROI_mask=bwmorph(ROI_mask,'clean',Inf);
            end

            
            %take just the biggest feature in the mask
            [bw,n]=bwlabel(ROI_mask);
            
            if n>1
                for i=1:n
                    ind(i)=numel(find(bw==i));
                end
                ROI_mask=(bw==max(ind));
            end
                
            if any(any(ROI_mask))                
                
                [xx,yy]=find(ROI_mask==1);

                tosum=images(:,:,j);
                ampguess=sum(sum(tosum(ROI_mask==1)));

                Lmax=max(xx);
                Lmin=min(xx);
                Cmax=max(yy);
                Cmin=min(yy);

                ROI_sL=Lmax-Lmin;
                ROI_sC=Cmax-Cmin;

                p0=[fix(datafitsize/2)*a, fix(datafitsize/2)*a, ROI_sC*a/2, ROI_sL*a/2, bck, ampguess ];
                
            else
               
                p0=[fix(datafitsize/2)*a, fix(datafitsize/2)*a, widthguess, widthguess, bck, ampguess ];
            end
            %***call the fitting function
            
            
            try
                [pF, var, exitflag] = MLEwG_elli_iter_correct(images(:,:,j),p0,a,0,0,0);
            catch
                pF(1:6)=NaN;
            end
            
            %try to improve the result of the fit if they are bad
           
            if ( pF(1)>2*datafitsize*a || pF(2)>2*datafitsize*a || pF(6)<0 ) 
                pF(1:6)=NaN;
                'diocan'
                i
            end
                    
            %***saves the fitted x,y positions and x-y widths of the PSF
         
            tracks3D(i).x(j)=x(j) + pF(1)/64-datafitsize/2;            
            tracks3D(i).y(j)=y(j) + pF(2)/64-datafitsize/2;
            tracks3D(i).SXs(j)=pF(3);
            tracks3D(i).SYs(j)=pF(4);
            tracks3D(i).bck(j)=pF(5);
            tracks3D(i).AMP(j)=pF(6);

        end
        %***estimate z position
        z=zeros(1,s(3));
        if ~isnan(pF(3))
            for b=1:s(3)

                [v,ind]=min(sqrt(  ((calibrationy)-tracks3D(i).SYs(b)).^2 + ((calibrationx)-tracks3D(i).SXs(b)).^2   ));
                z(b)=ascissacal(ind);

            end
            tracks3D(i).z=z;
        end
    end
    mkdir(save_track_path);
    save(strcat(save_track_path,'tracks3D.mat'),'tracks3D');
    if ~exist('indbadtracks')
        indbadtracks=0;
        save(strcat(save_track_path,'indbadtracks.mat'),'indbadtracks');
    end
        save(strcat(save_track_path,'indbadtracks.mat'),'indbadtracks');
end

