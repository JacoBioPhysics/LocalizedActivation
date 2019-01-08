function [masked_tracks] = filter_tracks_on_mask_rep( tracks, mask, Ntracks )
%UNTITLED Summary of this function goes here
%   INPUT:
%       tracks: tracksFinal output from u-track
%       mask  : bright field mask from python function BrightField_Segmentation
%   OUTPUT:
%       masked_tracks: only those tracks in the input that falls at least 80% in the mask

    if nargin < 3 %the function is called to filter a single track set

        saveresult=0;


        %mask=((mask-1)~=0);  

        mask=((mask)~=0);  
        mask=bwmorph(mask,'dilate');  %enlarge a bit the features in the mask (bacteria)
        inds=[];
        
        for i=1:length(tracks) 
        
            x=round(tracks(i).tracksCoordAmpCG(1:8:end));
            y=round(tracks(i).tracksCoordAmpCG(2:8:end));
            x=x(~isnan(x));
            y=y(~isnan(y));
            l=length(x);

            pointsOnMask=sum(mask(sub2ind(size(mask),y,x)));

            if ~(pointsOnMask>l*0.5)
                inds=horzcat(inds,i);
            end
            
        end
                
            

        tracks(inds)=[];
        masked_tracks=tracks;
        







    else %if the function is called for filtering many tracks 
        
        saveresult=1;

        tracks={};
        mask={};
        %have to improve here...
        prompt1='choose tracks to be filtered';
        prompt2='choose the corresponding masks for filtering';
        [trackname,trackpath]=uigetfile('D:\ANALISYS\2014\CLUSTER_OF_RECEPTORS_HIGH_THROUGHPUT_TRACKING\UNFLITERED_TRACKS\',prompt1,'multiselect','on');
        [maskname,maskpath]=uigetfile('D:\ANALISYS\2014\CLUSTER_OF_RECEPTORS_HIGH_THROUGHPUT_TRACKING\FILTERED_TRACKS_ON_MASK\*',prompt2,'multiselect', 'on');
        
        for i=1:Ntracks %load all the tracks and masks
            
            tracks{i}=importdata([trackpath trackname{i}])
            tracksnames{i}=trackname{i};        
            mask{i}=importdata([maskpath maskname{i}]);
            mask{i}=mask{i}~=0;
        end

        for i=1:Ntracks
            inds=[];
            %tracks{i}=filter_tracks_rep(tracks{i},10,100)   
                %better to filter out nan and short tracks later with get_displ..function
            for j=1:length(tracks{i})
                pointsOnMask=[];                
                l=length(tracks{i}(j).tracksCoordAmpCG(1:8:end));
                x=round(tracks{i}(j).tracksCoordAmpCG(1:8:end));
                y=round(tracks{i}(j).tracksCoordAmpCG(2:8:end));
                x=x(~isnan(x));
                y=y(~isnan(y));
                for k=1:length(x)                    
                    pointsOnMask(k)=(mask{i}(y(k),x(k)));
                end
                
                if ~(sum(pointsOnMask)>l*0.8)
                    inds=horzcat(inds,j);
                end

            end %for j
            tracks{i}(inds)=[];
            masked_tracks{i}=tracks{i};


            

           % masked_tracks{i}.filename=trackname;
           % masked_tracks{i}.maskname=maskname

        end  %for i                  

    end




    if saveresult

        save_path=uigetdir('D:\','choose the path where to save the filtered tracks');
        for m=1:length(tracksnames)
            tmpmasked=masked_tracks{m};
            save(strcat(save_path,'\masked_',tracksnames{m}), 'tmpmasked')
        end

    end %if

end %function
