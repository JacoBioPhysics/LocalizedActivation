
function [mobility_struct]=get_mobility_infos_fromUtrack_rep(source_dir, maskdir, Ndelays, nonan, nonoverlap)
%***************************************************************************************************************************************************************
%get_mobility_infos_fromUtrack_rep: repository function to compute histogram of displacements, average displacement per track, track length histogram, MSDs,
% for a given frame delay from utrack output. The histograms are computed for subdatasets from the same day but at consecutive time intervals. tHIS IS MEANT TO SEE IF THERE
% ARE VARIATIONS DUE TO MEASUREMENT TIME PASSING IN THE MEASURED QUANTITIES.
%***************************************************************************************************************************************************************
%INPUT:
%          source_dir           : directory where all the outputs from utrack belonging to different subdataset from the same set but different times are stored
%          Ndelay                : Number of delays (in frames) for which you want to compute the histograms
%          nonan                : logical, 0 or 1. if 0 tracks containing NANs are not excluded from the histograms, if 1 only connected tracks are retained
%          nonoverlap           : logical,0 or 1. if 0 overlapping time intervals are considered to computethe mean square displacement. if 1 only non overlapping time segments for each delay are considered.
%          maskdir              : directory where the masks are stored.
%          
%OUTPUT:
%          mobility_struct : structure containing all the analysis result
%               .single_track_av_disp_tot    : cell array. Each entry correspond to a different frame delay. The complete set of average displacements per track from all subdataset for that delay are in each entry
%               .displ_tot                    : cell array. Each entry correspond to a different frame delay. The complete set of displacements for all subdataset for that delay are in each entry
%               .single_track_av_disp_perdataset : cell array of cell arrays. Each first entry correspond to a different subdataset. Each second entry correspond to the frame delay. The complete set of average displacements per track for that subdataset are in each entry
%               .displ_perdataset                : cell array of cell arrays. Each first entry correspond to a different subdataset. Each second entry correspond to the frame delay. The complete set of displacements for that subdataset are in each entry
%               .MSD_tot                     : cell array. Each entry correspond to a different frame delay. The average of squared displacements from all subdataset for that delay are in each entry     
%               .MSD_perdataset                  :
%               .source_dir                  : directory containing the files analyzed
%               .files                       : names of the files analyzed            
%***************************************************************************************************************************************************************
    


    if nargin<5
        source_dir=uigetdir
        Ndelays=10;
        nonan=1;
        nonoverlap=1;
    end

    single_track_av_disp_tot={};
    
    D=dir(strcat(source_dir,'\*.mat'));
    Dm=dir([maskdir '\*.tif']);    
            
    inds=find(source_dir=='\');
    dataset=source_dir(inds(end-1):end);

    single_track_av_disp_tot=cell(Ndelays,1);
    displ_tot=cell(Ndelays,1); 
    MSD_tot={};

    track_length=cell(length(D),1);
    infotrack=cell(length(D),1);

    MSD_perdataset=cell(length(D),Ndelays);
    displ_perdataset=cell(length(D),Ndelays);
    single_track_av_disp_perdataset=cell(length(D),Ndelays);

%***************************************************************************************************************************************************************
%START THE ANALYSIS OF THE TRACKS 
%***************************************************************************************************************************************************************       
    for i=1:length(D)
        
        tmptrack=importdata(strcat(source_dir,'\',D(i).name));    
        %if nonan=0 then the displacement will come also from tracks containing NANs, if nonan=1 then the tracks containing NANs will not be taken into account
        if nonan
            mytrack=filter_tracks_noNAN_rep(tmptrack,10,100);
        else
            mytrack=filter_tracks_rep(tmptrack,10,100);
        end

        if maskdir~=0

            mask=importdata([maskdir '\' Dm(i).name]);
                        
            mytrack=Spatial_analysis_from_Utrack_rep(mytrack, mask);

            infotrack{i}=mytrack;
            figure
            try
                plot_tracks(infotrack{i})
                title(strcat( dataset, num2str(i) ), 'fontsize',20 );
            catch

            end

        end



       
        for m=1:Ndelays
            
            mydisp=[];%displacements per delay
             
            %initialize arrays into the cell array initialized above for the _tot variables            
            
            

            for k=1:(length(mytrack))
                
                tmpdisp=[];%displacements per track

                x=mytrack(k).tracksCoordAmpCG(1:8:end)*64;
                y=mytrack(k).tracksCoordAmpCG(2:8:end)*64;

                if(m==1)
                    track_length{i}=vertcat(track_length{i},length(x));
                end

                if nonoverlap
                    for j=1:m:fix(length(x)/m)-1

                        if ~isnan(x(j)) && ~isnan(y(j)) ... %&&~isnan(track.z(j)) ...
                            && ~isnan(x(j+m)) && ~isnan(y(j+m)) % &&~isnan(track.z(j+1))                                    
                            tmpdisp(j)=( abs(x(j+m)-x(j))^2+abs(y(j+m)-y(j))^2 );
                            mydisp=vertcat(mydisp, (tmpdisp(j)));
                        end
                        
                    end 

                else

                    for j=1:length(x)-m-1

                        if ~isnan(x(j)) && ~isnan(y(j)) ... %&&~isnan(track.z(j)) ...
                            && ~isnan(x(j+m)) && ~isnan(y(j+m)) % &&~isnan(track.z(j+1))                                    
                            tmpdisp(j)=( abs(x(j+m)-x(j))^2+abs(y(j+m)-y(j))^2 );
                            mydisp=vertcat(mydisp, (tmpdisp(j)));
                        end
                        
                        
                    end
                end

                single_track_av_disp_tot{m}=vertcat(single_track_av_disp_tot{m},mean(tmpdisp(:)));                
                single_track_av_disp_perdataset{i,m}=vertcat(single_track_av_disp_perdataset{i,m},mean(tmpdisp(:)));

            end
            
            displ_perdataset{i,m}=(mydisp);
            MSD_perdataset{i,m}=mean(mydisp);
            displ_tot{m}=vertcat(displ_tot{m}, mydisp);            
            
        end
        
    end    

%***************************************************************************************************************************************************************
%FINAL PLOTTING MODULE
%***************************************************************************************************************************************************************    
  %{   
    MSD_tot_err={}; %compute the total MSD and its standard error of the mean
    for i=1:Ndelays
        tmp=[];
        for k=1:length(D)
            tmp=vertcat(tmp,MSD_perdataset{k,i});
        end
        MSD_tot_err{i}=std(tmp)/sqrt(length(tmp));
        MSD_tot{i}=mean(tmp);
    end    
 

    cc=jet(2*length(D));
    
    figure %plot single track average displacement in time
    for q=1:length(D)
        Histo_setted(sqrt(cat(1,single_track_av_disp_perdataset{q,1})),(0:15:250),'pmf','delay 1','y',cc(2*q,:),q);
    end
    title(strcat('average track step displacement',dataset))

    figure %plot single track average displacement in time
    for q=1:length(D)
        Histo_setted(sqrt(cat(1,displ_perdataset{q}{2})),(0:15:250),'pmf','delay 1','y',cc(2*q,:),q);
    end
    title(strcat('step displacement',dataset));
%}
%***************************************************************************************************************************************************************
%FINAL SAVING MODULE
%***************************************************************************************************************************************************************    
    

    mobility_struct=struct();

    mobility_struct.displ_tot=displ_tot;
    mobility_struct.displ_perdataset=displ_perdataset;
    mobility_struct.single_track_av_disp_tot=single_track_av_disp_tot;
    mobility_struct.single_track_av_disp_perdataset=single_track_av_disp_perdataset;
    mobility_struct.MSD_tot=MSD_tot;
    mobility_struct.MSD_tot_err=MSD_tot_err;
    mobility_struct.MSD_perdataset=MSD_perdataset;

    mobility_struct.infotracks=infotrack;

    mobility_struct.track_length=track_length;
    mobility_struct.files = source_dir;
    mobility_struct.files = D;


    save('C:\Users\solari\Desktop\provastruct.mat','mobility_struct')
      
end
