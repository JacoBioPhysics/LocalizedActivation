%	AUthor: Jacopo Solari
%	date: started in December 2014
%	
%	SYNOPSYS: this class allows to analyze dataset from localized activation. Once defined an instance, you can call methods to load paths of folders 
%			  to analyze, compute mobility structures out of those and perform




classdef localized_activation < handle

	properties

		datasetname='' 						%name of the complete dataset, could be a date like 2014-12-14
		mobility_structures={};				%loaded or computed mobility structures, see get_mobility_structures function or method
		target_folders={};					%target folders for utrack analisys containing the prepared converted .tif images
		tracksFinal={};						%tracksFinal variables ouput of utrack, you have one for each datapaths{i}
		subdataset={};						%localized_cativation classes, childrens of the first instance of the class
		childrens=0;						%index for initialization: the first time you define a localized_activation object, it will not ask for the datapaths but just for the complete datasetname
		prepareddatafolder=''				%path where to store the prepared converted .tif images
		frames_separators={'2','1500'};		%frames interval to be analyzed with utrack
		saveroot='D:\ANALISYS\2015\1.CLUSTER_OF_RECEPTORS_TRACKING\1.Tracks_and_masks_for_mobility_structures\' ;		%path root to where to save the results of the tracking
		save_track_path={}					%complete path for saving
		datapaths={};						%names of the folders in a subdataset could be datasetname/stimulus/solari20141214153020
		maskpaths={};						%paths of the segmantation masks corresponding to each datapaths{i} for each datapaths{i} you will have the corresponding maskpath{i}		
		dt='0.065';							%time interval between two frames;
		expdate='';							%date of the experiment
        tracksmaskpath='';
	end


	methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parent class constructor, it will ask just for the complete dataset name

		function obj=localized_activation_new_mingus(childrens)
			if nargin<1				
				obj.childrens=0;
			else
				obj.childrens=1;
			end
			%ask to give a name to the subdataset
			if obj.childrens
				obj.datasetname=char(inputdlg('insert the name of the subdataset'));
			else
				obj.datasetname=char(inputdlg('insert the name of the complete dataset'));
			end
			%ask to select the folder containing data in Frames000001 folders for subdatasets
			if obj.childrens
				dataroot=uigetdir('D:\','select the folder containing the datasets to analyze');
				dd=dir(dataroot);
				dd(1:2)=[];
				dd=dd(cat(1,dd.isdir));
			
				for i = 1:length(dd)
			    	folderlist{i}=[dataroot '/' dd(i).name '/Frames000001'];
				end

				obj.datapaths=folderlist;
			end
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%method for initialization of the childrens of the parent instance. it will ask for subdatasetname (e.g. nostimulus, stimulus ..)
%and for the datapaths of the main folder containing single expreiments run folders containing Frames000001 folders

		function initialize(obj)

			obj.expdate=inputdlg('insert the date when the experiment was performed')			
			n=inputdlg('how many subdataset you want?')

			for i=1:str2num(n{1})
				obj.subdataset{i}=localized_activation_new_mingus(1);				
				obj.subdataset{i}.prepareddatafolder=['D:\ANALISYS\TRIALS\' obj.datasetname,'/'];
				obj.subdataset{i}.expdate=obj.expdate;

			end



		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%method to load pre-computed mobility structures

		function load_mobility_struct(obj)

			%load, if you want, some mobility structures.
			mobstruct_folder=uigetdir('D:/ANALISYS/2015/1.CLUSTER_OF_RECEPTORS_TRACKING/1.Tracks_and_masks_for_mobility_structures','choose the folder containing all the mobility structures');	
			ms_dir=dir(mobstruct_folder);
			ms_dir(1:2)=[]
			%ms_dir=ms_dir(cat(1,ms_dir.isdir));
			[selection] = listdlg('PromptString','Select which mobility structures to load:','ListString',{ms_dir.name})
			cd(mobstruct_folder)
			for i=1:length(selection)
				obj.mobility_structures{i}=importdata(ms_dir(selection(i)).name)
				obj.mobility_structures{i}.name=ms_dir(selection(i)).name
				'loaded'
			end
		end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%method to load a single pre-computed mobility structures and append to the la_class

		function load_single_mobility_struct(obj,i)

			%load, if you want, some mobility structures.
			[ms_name,mobstruct_folder]=uigetfile('D:/ANALISYS/2015/1.CLUSTER_OF_RECEPTORS_TRACKING/1.Tracks_and_masks_for_mobility_structures/*.mat','choose the folder containing all the mobility structures');	

			obj.mobility_structures{i}=importdata([mobstruct_folder ms_name]);
			obj.mobility_structures{i}.name=ms_name;

		end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%method to perform all tracking on each subdataset

		function all_dotrack(obj)

			for i=1:length(obj.subdataset)

				obj.subdataset{i}.dotrack()

			end

		end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%method to launch the utrack routine on single subdatasets

		function [tracks]=dotrack(obj)

			%[selection] = listdlg('ListString',obj.datapaths)

			for i=1:length(obj.datapaths)

				destination= [ obj.prepareddatafolder, '/' , obj.datasetname ,'/'];
		    	obj.target_folders{i}=[destination] %makes the directory inside the sub_datasetname folder indicating the frame number of the images inside the folder
		    	
   				mkdir(obj.target_folders{i});  %creates the directory inside the directory dataset
   				mkdir([ obj.saveroot, obj.expdate{:},'/DATA\']);
    			obj.save_track_path{i}=[ obj.saveroot, obj.expdate{:},'/DATA\',obj.datasetname, '/TRACKS\' ];				
				

			end

			for h=1:length(obj.target_folders)			
        			        
        	    PrepareImagesForUtrack_rep ... %converts in .tif
          	 		 (obj.datapaths{h},obj.target_folders{h}, obj.frames_separators{1}, obj.frames_separators{2},'T');
		   		cd(char(obj.target_folders{h}))

		    	[movieInfo,exceptions,localMaxima,background,psfSigma]=obj.detection(h)

		    	obj.tracking(h,movieInfo,exceptions,localMaxima,background,psfSigma)

	    	    ind2=find(obj.datapaths{h}=='/')
			    timedatedata=obj.datapaths{h}(ind2(end-1)+7:ind2(end)-1)
			 
			    strcat(obj.save_track_path{h},'tracksFinal',timedatedata,'.mat')
			    mkdir(obj.save_track_path{h})
			    tracks=obj.tracksFinal{h}
			    save(strcat(obj.save_track_path{h},'tracksFinal',timedatedata(2:end),'.mat'),'tracks')

	    	end
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%method for computing mobility structures out of tracksFinal variables
%folders have to be formatted like this:
%
%	date\nostimulus 	\MASK
%						\TRACKS
%		\stimulus 		\MASK
%						\TRACKS


		function get_all_mobility_struct(obj, OnSpot)

			if nargin <2
				OnSpot = 0; %this flag decides wheter to filter all the tracks only for on-spot activation or to keep all the tracks
			end

	    	%source_dir='C:\Users\solari\Desktop\utrack_output_20141216\2014-10-21\NOSTIMULUS\tracks' %'D:\ANALISYS\TRIALS\tracksandmask_spattrial\track' 
	    	source_dirs=uigetdir('D:\ANALISYS\2015\1.CLUSTER_OF_RECEPTORS_TRACKING\1.Tracks_and_masks_for_mobility_structures','choose the folder containing all the tracks and masks for all datasets');	    	
	    	save_dir=uigetdir('D:\ANALISYS\2015\1.CLUSTER_OF_RECEPTORS_TRACKING\1.Tracks_and_masks_for_mobility_structures','choose the folder where to save all the mobility structures');	    	
	    	%maskdir='C:\Users\solari\Desktop\utrack_output_20141216\2014-10-21\NOSTIMULUS\masks' %'D:\ANALISYS\TRIALS\tracksandmask_spattrial\mask' 

	    	d=dir(source_dirs);
	    	d(1:2)=[];
	    	d(~[d.isdir])=[];

	    	if sum([d.isdir]) ~= length(obj.subdataset)
	    		warndlg('The number of dataset does not correspond to the number of folders!')
    			return
			end	

	    	if OnSpot

	    		im_disks = create_disks_Microlens_array();
			else
				im_disks=0;
			end

			for i=1:length(obj.subdataset)

			    tracksmaskpath=strcat(source_dirs,'/', d(i).name);
			    obj.tracksmaskpath=tracksmaskpath;
			    obj.mobility_structures{i}=obj.subdataset{i}.get_mobility_struct(tracksmaskpath, OnSpot, im_disks);
				obj.mobility_structures{i}.name=obj.subdataset{i}.datasetname;
				['finished analyzing' obj.subdataset{i}.datasetname]
				ms=obj.mobility_structures{i};
				save(strcat(save_dir,'/','mobility_struct',num2str(i),'.mat'),'ms')
			end


		end



		function [mobility_struct]=get_mobility_struct(obj,tracksmaskpath,OnSpot,im_disks)
		%***************************************************************************************************************************************************************
		%get_mobility_infos_fromUtrack_rep: repository function to compute histogram of displacements, average displacement per track, track length histogram, MSDs,
		% for a given frame delay from utrack output. The histograms are computed for subdatasets from the same day but at consecutive time intervals. tHIS IS MEANT TO SEE IF THERE
		% ARE VARIATIONS DUE TO MEASUREMENT TIME PASSING IN THE MEASURED QUANTITIES.
		%***************************************************************************************************************************************************************
		%INPUT:
		%          source_dir           : directory where all the outputs from utrack belonging to different subdataset from the same set but different times are stored
		%                                   for instance all the result of utrack from a "nostimulus" folder    
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




%initialize parameters
			pixelsize=106;
			%pixelsize=64;
	        Ndelays=10;
	        nonan=1;
	        nonoverlap=0;
	        timedelay=0.2;
	        %timedelay=0.3;
	        nPoints_forfit=3; %number of points in the singletrack MSD for fitting to extract D
	        nPoints_forfit_alpha=5; %number of points in the singletrack MSD for fitting to extract alpha
	        threshold_track_length=1;

	        act_dt=249;

%look for tracks and masks	    		    
		    d2=dir(tracksmaskpath)
            
		    d2(1:2)=[];
            d2(~[d2.isdir])=[];
		    D=dir(strcat(tracksmaskpath, '/', d2(2).name, '/*.mat')); %variable containing all the paths to the tracks
		    Dm=dir([tracksmaskpath, '/',d2(1).name, '/*.tif']);    
		    D.name
		    %inds=find(source_dir=='/');
		    %dataset=source_dir(inds(end-1):end);

%initialize arrays to be filled in the analysis		    
		    displ_perdataset=cell(length(D),Ndelays);
		    autocorr_perdataset=cell(length(D),Ndelays);

		    displ_tot=cell(Ndelays,1); 
		    POLEdisp_tot=cell(Ndelays,1); 
		    MIDDLEdisp_tot=cell(Ndelays,1); 

		    MSD_tot={};
		    autocorr_tot={};

		    track_length=cell(length(D),1);
		    infotrack=cell(length(D),1);

		    MSD_perdataset=cell(length(D),Ndelays);
		    MSD_perdataset_err=cell(length(D),Ndelays);

			n_independent=cell(Ndelays,1);
			n_independent_tracks=cell(Ndelays,1);
		    single_track_av_disp_tot=cell(Ndelays,1);
		    POLARsingle_track_av_disp_tot=cell(Ndelays,1);
		    MIDDLEsingle_track_av_disp_tot=cell(Ndelays,1);
		    single_track_av_disp_perdataset=cell(length(D),Ndelays);
		    single_track_MSD={}

	    	

		%***************************************************************************************************************************************************************
		%START THE ANALYSIS OF THE TRACKS 
		%***************************************************************************************************************************************************************       
			
			masks=cell(1,length(D)); 		
			masks_paths=cell(1,length(D));	%initialize the paths of the masks files and the masks files contatiner

	        n_independent_tracks=0;			
	        all_filtered_tracks=[];

		    for i=1:length(D)
		        i

		        if Dm(i).name~=0 & OnSpot

		            mask=importdata([strcat(tracksmaskpath, '/', d2(1).name) , '/' Dm(i).name]);		           
		            mask=mask~=0;

		            mask=create_LA_spot_mask(im_disks, mask);

		            mask=bwlabel(mask);
            		masks{i}=mask; 

        		end

		        tmptrack=importdata(strcat(tracksmaskpath, '/', d2(2).name, '/' ,D(i).name));    

		        %if nonan=0 then the displacement will come also from tracks containing NANs, if nonan=1 then the tracks containing NANs will not be taken into account
		        if OnSpot

		        	mytrack=filter_tracks_on_LA_spots(tmptrack, mask,threshold_track_length,1000, act_dt);
	        	else
		            mytrack=filter_tracks_noNAN_rep(tmptrack,threshold_track_length,1000);
	            end
		        %else
		         %   mytrack=filter_tracks_rep(tmptrack,threshold_track_length,100);
		        %end

		        if Dm(i).name~=0 & OnSpot

		            mytrack=Spatial_analysis_from_Utrack_rep(mytrack, mask, timedelay, pixelsize, nPoints_forfit, nPoints_forfit_alpha); %repository function that compute spatial resolved analysis of mobility of tracks based on the mask
	            	['Just finished analyzing' D(i).name]
	            	
		            infotrack{i}=mytrack;
		            all_filtered_tracks=[all_filtered_tracks,mytrack];

		            if i==1
		            	all_filtered_tracks=mytrack;
	            	else
	            		all_filtered_tracks=[all_filtered_tracks,mytrack];
            		end
        		end

	        	if Dm(i).name~=0

		            mask=importdata([strcat(tracksmaskpath, '/', d2(1).name) , '/' Dm(i).name]);		           
		            mask=mask~=0;
		            mask=bwlabel(mask);
            		masks{i}=mask;               		

		            mytrack=Spatial_analysis_from_Utrack_rep(mytrack, mask, timedelay, pixelsize, nPoints_forfit, nPoints_forfit_alpha); %repository function that compute spatial resolved analysis of mobility of tracks based on the mask
	            	['Just finished analyzing' D(i).name]
	            	
		            infotrack{i}=mytrack;
		            all_filtered_tracks=[all_filtered_tracks,mytrack];

		            if i==1
		            	all_filtered_tracks=mytrack;
	            	else
	            		all_filtered_tracks=[all_filtered_tracks,mytrack];
            		end



		        end

		        d=1 %time interval in which you compute the velocity. If d-->0 you have instantaneoud velocity
		        
		        for m=1:Ndelays

		        	n_independent{m}=0;
		        	

		            mydisp=[];%displacements per delay
		            POLEdisp=[];
		            MIDDLEdisp=[];
		            myautocorr=[];
		            %initialize arrays into the cell array initialized above for the _tot variables            
		            
		            %tmpEMSD=[];

		            for k=1:(length(mytrack))
		                
		                tmpdisp=[];%displacements per track
	           			tmpPOLEdisp=[];
		            	tmpMIDDLEdisp=[];


		                x=mytrack(k).tracksCoordAmpCG(1:8:end)*pixelsize;
		                y=mytrack(k).tracksCoordAmpCG(2:8:end)*pixelsize;

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

		                    for j=1:(length(x)-m-1)

	                        %compute the squared displacement

		                        if ~isnan(x(j)) && ~isnan(y(j)) ... %&&~isnan(track.z(j)) ...
		                            && ~isnan(x(j+m)) && ~isnan(y(j+m)) % &&~isnan(track.z(j+1))  

		                            tmpdisp(j)=( (x(j+m)-x(j))^2+(y(j+m)-y(j))^2 );

		                            if mytrack(k).position(1)=='p'
		                            	POLEdisp=vertcat(POLEdisp,tmpdisp(j));
		                            	tmpPOLEdisp=vertcat(tmpPOLEdisp,tmpdisp(j));	                            	
		                            elseif mytrack(k).position(1)=='m'
		                            	MIDDLEdisp=vertcat(MIDDLEdisp,tmpdisp(j));
		                            	tmpMIDDLEdisp=vertcat(tmpMIDDLEdisp,tmpdisp(j));
		                            end		                            	
		                            mydisp=vertcat(mydisp, (tmpdisp(j)));

		                        else
		                        	tempdisp(j)=NaN;
	                        	end

                			end

            			end %if/else

            			n_independent{m}=n_independent{m}+fix(length(x)/m);

						%compute the ensemble average MSD

						%tmpEMSD=vertcat(tmpEMSD,(x(1) - x(m+1))^2);	   						     			

        				%compute the velocity autocorrelation
                		for j=1:(length(x)-m-d-1)
                			m2=m-1;
	                        if ~isnan(x(j+d)) && ~isnan(y(j+d)) ... %&&~isnan(track.z(j)) ...
	                            && ~isnan(x(j+d+m2)) && ~isnan(y(j+d+m2)) % &&~isnan(track.z(j+1)) 

                    			vtemp=[x(j+m2+d)-x(j+m2),y(j+m2+d)-y(j+m2)]/(d*timedelay);                           
	                            vtemp_d=[x(j+d)-x(j),y(j+d)-y(j)]/(d*timedelay);
	                            tmpautocorr=dot(vtemp_d,vtemp);
	                           	myautocorr=vertcat(myautocorr,tmpautocorr);
                           	end

                       	end
	                        
		                        
		                    

		                single_track_av_disp_tot{m}=vertcat(single_track_av_disp_tot{m},nanmean(tmpdisp(:)));      
  		                POLARsingle_track_av_disp_tot{m}=vertcat(POLARsingle_track_av_disp_tot{m},nanmean(tmpPOLEdisp(:)));                
		                MIDDLEsingle_track_av_disp_tot{m}=vertcat(MIDDLEsingle_track_av_disp_tot{m},nanmean(tmpMIDDLEdisp(:)));                

		                single_track_av_disp_perdataset{i,m}=vertcat(single_track_av_disp_perdataset{i,m},nanmean(tmpdisp(:)));

		            end %for k
		            
		            displ_perdataset{i,m}=(mydisp);
		            MSD_perdataset{i,m}=vertcat(MSD_perdataset{i,m},nanmean(mydisp));
   			        MSD_perdataset_err{i,m}=vertcat(MSD_perdataset_err{i,m},std(mydisp)/sqrt(length(mydisp)));
		            displ_tot{m}=vertcat(displ_tot{m}, mydisp);  
		            POLEdisp_tot{m}=vertcat(POLEdisp_tot{m}, POLEdisp);  
		           	MIDDLEdisp_tot{m}=vertcat(MIDDLEdisp_tot{m}, MIDDLEdisp);            
		        	autocorr_perdataset{i,m}=nanmean(myautocorr);
		        	%EMSD{m} = mean(tmpEMSD);
		        	%EMSD_err{m}=std(tmpEMSD)/sqrt(length(tmpEMSD));

		        end %for m

				n_independent_tracks=n_independent_tracks+length(mytrack)		        
		    end    %for i		    

	        MSD_tot_err={}; %compute the total MSD and its standard error of the mean
	        MSD_tot_err_ind_steps={};
	        MSD_tot_err_ind_tracks={};

		    for i=1:Ndelays
		        tmp=[];
		        tmp2=[];
		        for k=1:length(D)
		            tmp=vertcat(tmp,MSD_perdataset{k,i});
		            tmp2=vertcat(tmp2,autocorr_perdataset{k,i});
		        end
		        MSD_tot_err{i}=std(displ_tot{i})/sqrt(length(displ_tot{i}));
		      	MSD_tot_err_ind_steps{i}=std(displ_tot{i})/sqrt((n_independent{i}));
		     	MSD_tot_err_ind_tracks{i}=std(displ_tot{i})/sqrt((n_independent_tracks));
		        MSD_tot{i}=mean(displ_tot{i});
		        autocorr_tot{i}=mean(tmp2);
		    end 

			%compute the velocity autocorrelation function	


			%***************************************************************************************************************************************************************
			%FINAL SAVING MODULE
			%***************************************************************************************************************************************************************    
	    

		    mobility_struct=struct();

		    mobility_struct.n_independent_tracks=n_independent_tracks;
		    mobility_struct.displ_tot=displ_tot;
		    mobility_struct.POLEdisp_tot=POLEdisp_tot;
		    mobility_struct.displ_perdataset=displ_perdataset;
		    mobility_struct.single_track_av_disp_tot=single_track_av_disp_tot;
		    mobility_struct.single_track_av_disp_perdataset=single_track_av_disp_perdataset;
		    mobility_struct.MIDDLEsingle_track_av_disp_tot=MIDDLEsingle_track_av_disp_tot;
			mobility_struct.POLARsingle_track_av_disp_tot=POLARsingle_track_av_disp_tot;
		    mobility_struct.VelocityAutocorr_perdataset=autocorr_perdataset;

		    mobility_struct.MSD_tot=MSD_tot;
		    mobility_struct.MSD_tot_err=MSD_tot_err;
		  	mobility_struct.MSD_tot_err_ind_steps=MSD_tot_err_ind_steps;
		  	mobility_struct.MSD_tot_err_ind_tracks=MSD_tot_err_ind_tracks;
		    mobility_struct.MSD_perdataset=MSD_perdataset;
    		mobility_struct.MSD_perdataset_err=MSD_perdataset_err;
		   % mobility_struct.EMSD=EMSD;
		    %mobility_struct.EMSD_err=EMSD_err;
		    mobility_struct.VelocityAutocorr=autocorr_tot;


		    mobility_struct.infotracks=infotrack;
		    mobility_struct.all_filtered_tracks=all_filtered_tracks;

		    mobility_struct.track_length=track_length;
		    mobility_struct.files = dir([tracksmaskpath '/' d2(2).name]);
		    mobility_struct.files(1:2)=[];

		    for i = 1:length(mobility_struct.files)
		    	mobility_struct.filespaths(i).name=[tracksmaskpath, filesep,  d2(2).name ,filesep,mobility_struct.files(i).name];
	    	end

		    mobility_struct.masks_paths = masks_paths;


		    delay=5;
		    pixelsize=106;
		    Ndelays=10;

			it=get_middlepolar_displ_rep (mobility_struct,delay,pixelsize,Ndelays);

		    mobility_struct.middlepoleInfo=it;



		    %['C:\Users\solari\Desktop\mobility_structures20141216\' obj.datasetname '/' obj.subdataset{1}.datasetname '.mat'],'mobility_struct'

		    %save(['C:\Users\solari\Desktop\mobility_structures20141216\' obj.datasetname '_' obj.subdataset{1}.datasetname '.mat'],'mobility_struct')
		      
		end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subfunction of the .dotrack method that does the localizadion detection part

		function [movieInfo,exceptions,localMaxima,background,psfSigma]=detection(obj, ind)

			movieParam.imageDir = char(obj.target_folders{ind}) %directory where images are

			d=dir('*.tif');

			firstim=d(1).name;
			lastim=d(end).name;
			movieParam.filenameBase = '/Image'; %image file name base
			movieParam.firstImageNum = str2num(firstim(regexp(firstim,'/d'))); %number of first image in movie
			movieParam.lastImageNum = str2num(lastim(regexp(lastim,'/d'))); %number of last image in movie
			movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

			%% detection parameters
			%*************************************************************************************************************************************
			%THIS PARAMETERS INFLUENCES THE QUALITY OF THE OUTPUT TRACKS, CHOOSE THEM CAREFULLY. A WORKING FINE SET OF THEM IS:
			%	detectionParam.psfSigma = 1.25; 
			%	detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); 
			%	detectionParam.visual = 0; 
			%	detectionParam.doMMF = 1; 
			%	detectionParam.bitDepth = 16; 
			%	detectionParam.alphaLocMax = 0.01; 
			%	detectionParam.numSigmaIter = 0; 
			%	detectionParam.integWindow = 2; 
			%*************************************************************************************************************************************

			detectionParam.psfSigma = 1.3; %point spread function sigma (in pixels)
			%detectionParam.psfSigma = 2.3; %point spread function sigma (in pixels) for 250x
			detectionParam.testAlpha = struct('alphaR',0.01,'alphaA',0.01,'alphaD',0.01,'alphaF',0); %THEHIGHER THE ALPHA, THE MORE FEATURE WILL BE RECOGNIZED alpha-values for detection statistical tests
			detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
			detectionParam.doMMF = 0; %1 if mixture-model fitting, 0 otherwise
			detectionParam.bitDepth = 16; %Camera bit depth
			detectionParam.alphaLocMax = 0.01; %alpha-value for initial detection of local maxima
			detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
			detectionParam.integWindow =1; %number of frames before and after a frame for time integration

			detectionParam.calcMethod = 'g';

			%absolute background info and parameters...
			background.imageDir = 'C:\Users\solari\Desktop\';

			background.filenameBase = 'BCK_ANDOR_emGAIN100.tif';
			background.alphaLocMaxAbs = 0.01;
			detectionParam.background = background;

			%% additional input

			%saveResults
			mkdir(strcat(movieParam.imageDir,'/detection_output'));
			saveResults.dir = strcat(movieParam.imageDir,'/detection_output'); %directory where to save input and output
			saveResults.filename = 'detectionTest1.mat'; %name of file where input and output are saved
			% saveResults = 0;

			%verbose state
			verbose = 1;

			%% run the detection function
			[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
			    detectSubResFeatures2D_rep(movieParam,detectionParam,saveResults,verbose);
		end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subfunction of the .dotrack method that does the linking-gapclosing part

		function tracking(obj,ind,movieInfo,exceptions,localMaxima,background,psfSigma)

			gapCloseParam.timeWindow = 1; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
			gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
			gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.

			%optional input:
			gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

			%% cost matrix for frame-to-frame linking

			%function name
			costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';

			%parameters

			parameters.linearMotion = 0; %use linear motion Kalman filter.

			parameters.minSearchRadius = 1; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
			parameters.maxSearchRadius = 4; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
			parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

			parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
			parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

			parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
			% parameters.kalmanInitParam.searchRadiusFirstIteration = 10; %Kalman filter initialization parameters.

			%optional input
			parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.

			costMatrices(1).parameters = parameters;
			clear parameters

			%% cost matrix for gap closingss

			%function name
			costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';

			%parameters

			%needed all the time
			parameters.linearMotion = 0; %use linear motion Kalman filter.

			parameters.minSearchRadius = 1; %minimum allowed search radius.
			parameters.maxSearchRadius = 4; %maximum allowed search radius.
			parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.

			parameters.brownScaling = [0.1 0.02]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
			% parameters.timeReachConfB = 3; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
			parameters.timeReachConfB = gapCloseParam.timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).

			parameters.ampRatioLimit = [0.7 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

			parameters.lenForClassify = 1; %minimum track segment length to classify it as linear or random.

			parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
			parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

			parameters.linStdMult = 1*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.

			parameters.linScaling = [0.1 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
			% parameters.timeReachConfL = 4; %similar to timeReachConfB, but for the linear part of the motion.
			parameters.timeReachConfL = gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.

			parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

			%optional; if not input, 1 will be used (i.e. no penalty)
			parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).

			%optional; to calculate MS search radius
			%if not input, MS search radius will be the same as gap closing search radius
			parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

			costMatrices(2).parameters = parameters;
			clear parameters

			%% Kalman filter function names

			kalmanFunctions.reserveMem  = 'kalmanResMemLM';
			kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
			kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
			kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

			%% additional input

			%saveResults
			mkdir(strcat(pwd,'/tracking_output'))
			saveResults.dir = strcat(pwd,'/tracking_output'); %directory where to save input and output
			% saveResults.filename = 'tracksTest1DetectionAll1.mat'; %name of file where input and output are saved
			% saveResults = 0; %don't save results

			%verbose state
			verbose = 1;

			%problem dimension
			probDim = 2;

			%% tracking function call

			% [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(1:300),...
			%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

			% for i = 1 : 12
			for i = 1
			    %movieInfoTmp((i-1)*1200+1:i*1200) = movieInfo((i-1)*1200+1:i*1200);
			    saveResults.filename = ['tracks1All_' sprintf('%02i',i) '.mat'];
			    movieInfoTmp=movieInfo;
			    [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse_rep(movieInfoTmp,...
			        costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);		    	

			    clear movieInfoTmp
			    obj.tracksFinal{ind}=tracksFinal;
			end

		end %tracking


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	3D tracking and PHOTONCOUNTING	 methods	%%%%%%%%%%%%%%%

	
		function [tracks3D,photoncount]=get_all_3Dtracks(obj)

			for i=1:length(obj.mobility_structures)
				for j=1:length(obj.mobility_structures{i}.infotracks)
			
					[a]=obj.get_3Dtracks(i,j)
			

					obj.mobility_structures{i}.tracks3D{j}=a;
					
				
				
				end
			end
		end

		function [tracks3D]=get_3Dtracks(obj,mob_str_index,infotr_index)

			t=obj.mobility_structures{mob_str_index}.infotracks{infotr_index};
			t=filter_tracks_noNAN_rep(t,10,1000);
			%t=filter_tracks_on_mask_rep(t,obj.mobility_structures{mob_str_index}.masks{infotr_index});
			t=filter_tracks_after_activation_rep(t,250,3);
			obj.subdataset{mob_str_index}.datapaths{infotr_index}			
			tracks3D=fit_tracksFinal_fromUtrack_rep2(t,obj.subdataset{mob_str_index}.datapaths{infotr_index});
		end		

		function [AMP_all,AMP_pole,AMP_middle]=get_all_AMPs(obj, SaveOn)

			if SaveOn
				mysavepath='D:\ANALISYS\2017\loc_pre_estimation_Mortensen_20160705';
                %uigetdir('D:/ANALISYS/2015/1.CLUSTER_OF_RECEPTORS_TRACKING/','choose the folder where to sae the AMP files')
			end 

			L=length(obj.mobility_structures);

			%basepath='/\storage01\data\AMOLF\users\solari\JACKSON_DATA\2016-07-24\902_prolonged_nostimulation_withSwitch_after_first_dataset\'
            
            basepath=uigetdir('/\storage01\data\AMOLF\groups\shimizu-group\shared-data\CHETAN','choose the folder containing the tracking data organizaed in subdatasets') %select the folder with all the tracking data separated in subdatasets
			datasetpaths=dir(basepath);
			datasetpaths(1:2)=[];
			datasetpaths(~[datasetpaths.isdir]')=[];

            loc_var_tot=[];
            
			AMP_all={};
            AMP_all_sum={};
            AMP_all_mean={};
			AMP_pole={};
			AMP_middle={};
            loc_prec_all={};
            
            D_all={};
			D_pole={};
			D_middle={};
           
            TL_all={};
			TL_pole={};
			TL_middle={};            
            
            D_all_MSD={};
			D_pole_MSD={};
			D_middle_MSD={};            
            
			for i=1:L %loop over different FOVs		

				subdatasetpaths=dir([basepath filesep datasetpaths(i).name]);
				subdatasetpaths(1:2)=[];
				subdatasetpaths(~[subdatasetpaths.isdir]')=[];

				AMP_all{i}(1)=0;
				AMP_pole{i}(1)=0;
				AMP_middle{i}(1)=0;
                AMP_all_sum{i}(1)=0;
                AMP_all_mean{i}(1)=0;
                loc_prec_all{i}(1)=0;
                
				D_all{i}(1)=0;
				D_pole{i}(1)=0;
				D_middle{i}(1)=0; 
                
                D_all_MSD{i}(1)=0;
				D_pole_MSD{i}(1)=0;
				D_middle_MSD{i}(1)=0; 
                
                TL_all{i}(1)=0;
                TL_pole{i}(1)=0;
                TL_middle{i}(1)=0;
                
				for j=1:length(obj.mobility_structures{i}.infotracks)
					
					framespath=dir([basepath filesep datasetpaths(i).name  filesep subdatasetpaths(j).name filesep]);
					framespath(1:2)=[];
					framespath(~[framespath.isdir]')=[];

					images_path=[basepath filesep datasetpaths(i).name filesep subdatasetpaths(j).name filesep framespath(1).name];
					
					cd(images_path)
					d=dir(strcat(images_path,'/*.Tiff'));
					
					count=0;
					reverseStr='';

					for k=obj.mobility_structures{i}.infotracks{j} %loop over infotracks

						%nice progress bar
						percentDone = 100 * count / length(obj.mobility_structures{i}.infotracks{j});
						msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
						fprintf([reverseStr, msg]);
						reverseStr = repmat(sprintf('/b'), 1, length(msg));
						%%%%% 

						count=count+1;

						[a,PSF_size,amp_sum,loc_var]=obj.get_AMP(k,d); %call the functiona that returns the amp value from the first frame of the track via precise fitting.
						%k.AMP_firstframe=a;
                        loc_var_tot=[loc_var_tot,loc_var];
						obj.mobility_structures{i}.infotracks{j}(count).AMP_firstframe=a; %saves the AMP value in each infotrack
						
                        TL_all{i}=[TL_all{i},length(k.x)];
                        
						AMP_all{i}=[AMP_all{i},a];
                        AMP_all_sum{i}=[AMP_all_sum{i},amp_sum];
                        
                        if length(amp_sum)>5
                            AMP_all_mean{i}=[AMP_all_mean{i},nanmean(amp_sum(1:5))];
                        else
                            AMP_all_mean{i}=[AMP_all_mean{i},nanmean(amp_sum)];
                        end
                        
                        D_all{i}=[D_all{i},k.D_ST];
                        D_all_MSD{i}=[D_all_MSD{i},k.D/4];
                        
                        loc_prec_all{i}=[loc_prec_all{i},PSF_size/sqrt(a)];
                        k.AMP=a;                          
                        
						if k.position(1)=='p'
							AMP_pole{i}=[AMP_pole{i},a];
                            D_pole{i}=[D_pole{i},k.D_ST];
                            D_pole_MSD{i}=[D_pole_MSD{i},k.D/4];
                            TL_pole{i}=[TL_pole{i},length(k.x)];
                            
						end
						if k.position(1)=='m'
							AMP_middle{i}=[AMP_middle{i},a];									
                            D_middle{i}=[D_middle{i},k.D_ST];
                            D_middle_MSD{i}=[D_middle_MSD{i},k.D/4];
                            TL_middle{i}=[TL_middle{i},length(k.x)];
						end
						

					end%for k		
					fprintf(['finished mobility structure ' num2str(i) ' subdataset ' num2str(j) '/n']	)
				end%for j
			
				AMP_all{i}(1)=[];
				AMP_pole{i}(1)=[];
				AMP_middle{i}(1)=[];
                AMP_all_sum{i}(1)=[];
                loc_prec_all{i}(1)=[];
                
                D_all{i}(1)=[];
				D_pole{i}(1)=[];
				D_middle{i}(1)=[];    

				obj.mobility_structures{i}.AMP_firstframe_all=AMP_all{i};
				obj.mobility_structures{i}.AMP_firstframe_pole=AMP_pole{i};
				obj.mobility_structures{i}.AMP_firstframe_middle=AMP_middle{i};			

				AMP_all_temp=AMP_all{i};
				AMP_pole_temp=AMP_pole{i};		
				AMP_middle_temp=AMP_middle{i};
                loc_prec_ALL=loc_prec_all{i};


				if SaveOn
					save([mysavepath filesep 'AMP_MobStruct_all_sum' num2str(i) '.mat' ], 'AMP_all_sum');
                    save([mysavepath filesep 'AMP_MobStruct_all' num2str(i) '.mat' ], 'AMP_all_temp');
                    save([mysavepath filesep 'AMP_MobStruct_all_mean' num2str(i) '.mat' ], 'AMP_all_mean');
					save([mysavepath filesep 'AMP_MobStruct_pole' num2str(i) '.mat' ], 'AMP_pole_temp');
					save([mysavepath filesep 'AMP_MobStruct_middle' num2str(i) '.mat' ], 'AMP_middle_temp');
                    save([mysavepath filesep 'locPrec_MobStruct_correct' num2str(i) '.mat' ], 'loc_prec_ALL');
                    
   					save([mysavepath filesep 'D_all' num2str(i) '.mat' ], 'D_all');
                    save([mysavepath filesep 'D_all_MSD' num2str(i) '.mat' ], 'D_all_MSD');
                    save([mysavepath filesep 'D_pole_MSD' num2str(i) '.mat' ], 'D_pole_MSD');
					save([mysavepath filesep 'D_pole' num2str(i) '.mat' ], 'D_pole');
					save([mysavepath filesep 'D_middle' num2str(i) '.mat' ], 'D_middle');
                    
                    save([mysavepath filesep 'TL_all' num2str(i) '.mat' ], 'TL_all');
                    save([mysavepath filesep 'TL_pole' num2str(i) '.mat' ], 'TL_pole');
					save([mysavepath filesep 'TL_middle' num2str(i) '.mat' ], 'TL_middle');
                    save([mysavepath filesep 'loc_var_tot_fromMLEwG' num2str(i) '.mat' ], 'loc_var_tot');
				end

			end%for i

		end %fun
	
		function [aa,PSF_size,amp_sum,loc_var]=get_AMP(obj,infotrack,d,PSFw)	

		
			PSFw=100;
            
        
            
            try
                a=106;
                datafitsize=3;

                ampforsum=[];
                loc_var=[];
                amps=[];
                vars=[];
                bcks=[];
                XCoord=[];
                YCoord=[];                

                amp_bcks=[];        
                goodnesses=[];
                SXs=[];
                SYs=[];

                x=infotrack.tracksCoordAmpCG(1:8:end);
                y=infotrack.tracksCoordAmpCG(2:8:end);        
                first=infotrack.seqOfEvents(1);
                last=infotrack.seqOfEvents(2);

                try
                    [data_frames,xrest,yrest,datafitsize]=cut_frames_rep(first,last,d,x,y,datafitsize);
                catch
                    pF=[];
                    aa=0;
                    PSF_size=0;
                    fprintf('/nthe track is too close to the edge of the FOV\n')
                                        
                    amp_sum=0;
                    aa=0;
                    PSF_size=0;
                    'cazzo something was wrong in the get_amp'
                    
                    return
                end

                for i =1:length(data_frames)
                    data=data_frames{i};
                    t=thresholdOtsu(data)*0.8;
                    ampp=sum(sum(data*(data>t)));
                    
                    
                    [y0,x0]=find(data==max(data(:)));
                    %p0=[(datafitsize+1)*a,(datafitsize+1)*a ,PSFw,0,ampp];
                    p0=[x0*a,y0*a ,PSFw,1,ampp];
                    if i==1
                        pF_prev=p0;
                    end
                    sd=size(data);
                    [xx,yy]=meshgrid(1:sd(1),1:sd(2));
                    xx=xx*a;
                    yy=yy*a;

                    [pF,loc_var_tmp]=MLEwG((data),p0,a,0,1,0);

                    if pF(5)>100 && pF(5)<20000
                        ampforsum=[ampforsum,pF(5)];   
                        loc_var=[loc_var,loc_var_tmp];
                    else
                        ampforsum=[ampforsum,0];
                    end

                    if i ==1
                        aa=pF(5);                        
                        PSF_size=pF(3);
                    end
                end

                amp_sum=sum(ampforsum);
                            
            catch  
                amp_sum=0;
                aa=0;
                PSF_size=0;
                loc_var=0;
                'cazzo something was wrong'
            end
            
                            
		end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	REARRANGEMENT OF THE LOC AACT CLASS	 methods	%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		function separate_pole_middle_tracks(obj)

		%first divide into polar and middle tracks 	
			for j=1:length(obj.mobility_structures)
				poleind=[];
				middleind=[];
				poleTlength=[];
				middleTlength=[];
				poletracklength=[];
				for k=1:length(obj.mobility_structures{j}.all_filtered_tracks)
					if obj.mobility_structures{j}.all_filtered_tracks(k).position(1)=='m'						
						middleind=[middleind,k];
						middleTlength=[middleTlength,length(obj.mobility_structures{j}.all_filtered_tracks(k).tracksCoordAmpCG(1:8:end))];
					elseif obj.mobility_structures{j}.all_filtered_tracks(k).position(1)=='s'

					else
						%poletracklength=[poletracklength,length(obj.mobility_structures{j}.all_filtered_tracks(k).tracksCoordAmpCG(1:8:end))];
						poleind=[poleind,k];
						poleTlength=[poleTlength,length(obj.mobility_structures{j}.all_filtered_tracks(k).tracksCoordAmpCG(1:8:end))];
					end
				end

				poleTracks=obj.mobility_structures{j}.all_filtered_tracks(poleind);
				middleTracks=obj.mobility_structures{j}.all_filtered_tracks(middleind);

				obj.mobility_structures{j}.middlepoleInfo.middleTracks=middleTracks;
				obj.mobility_structures{j}.middlepoleInfo.middleTlength=middleTlength;

				obj.mobility_structures{j}.middlepoleInfo.poleTracks=poleTracks;
				obj.mobility_structures{j}.middlepoleInfo.poleTlength=poleTlength;

				obj.mobility_structures{j}.poletracklength=poletracklength;

			end
		end

		function separate_in_equal_pole_middle_tracks(obj)

			obj.separate_pole_middle_tracks;

			Mlengths=[];
			Plengths=[];
			for i=1:length(obj.mobility_structures)
				Mlengths=[Mlengths,length(obj.mobility_structures{i}.middlepoleInfo.middleTracks)];
				Plengths=[Plengths,length(obj.mobility_structures{i}.middlepoleInfo.poleTracks)];
			end

			for i=1:length(obj.mobility_structures)
				[ms,minds]=sort(obj.mobility_structures{i}.middlepoleInfo.middleTlength);
				[ps,pinds]=sort(obj.mobility_structures{i}.middlepoleInfo.poleTlength);

				%filter out shorter tracks in order to have same number of polar-middle tracks in all mob. struct.
				minds=minds(end-min(Mlengths)+1:end);
				pinds=pinds(end-min(Plengths)+1:end);

				tmpmid=obj.mobility_structures{i}.middlepoleInfo.middleTracks(minds);
				tmppol=obj.mobility_structures{i}.middlepoleInfo.poleTracks(pinds);

				obj.mobility_structures{i}.all_filtered_tracks_sameN=[tmpmid,tmppol];

			end

		end						

		function compute_MSD_on_sameNtracks(obj)

			obj.separate_in_equal_pole_middle_tracks;
			Ndelays=5;
			pixelsize=106;
			for i=1:length(obj.mobility_structures)
				n_independent_sameN=zeros(length(obj.mobility_structures),1);
				displ_sameN=cell(1,Ndelays);
				for m=1:Ndelays	
					m
					for k=1:length(obj.mobility_structures{i}.all_filtered_tracks_sameN)
						
						x=obj.mobility_structures{i}.all_filtered_tracks_sameN(k).tracksCoordAmpCG(1:8:end)*pixelsize;
						y=obj.mobility_structures{i}.all_filtered_tracks_sameN(k).tracksCoordAmpCG(2:8:end)*pixelsize;	
			            for j=1:(length(x)-m-1)

			            %compute the MSD

			                if ~isnan(x(j)) && ~isnan(y(j)) ... %&&~isnan(track.z(j)) ...
			                    && ~isnan(x(j+m)) && ~isnan(y(j+m)) % &&~isnan(track.z(j+1))  

			                    tmpdisp(j)=( (x(j+m)-x(j))^2+(y(j+m)-y(j))^2 );
			                    displ_sameN{m}=horzcat(displ_sameN{m},tmpdisp);
			                    clear tmpdisp
			                else
			                	tempdisp_m(j)=NaN;
			            	end

						end				
						clear x y				
					end
					n_independent_sameN(i)=n_independent_sameN(i)+k;			
				end

		        MSD_tot_sameN={}; %compute the total MSD and its standard error of the mean
		        MSD_tot_err_ind_tracks_sameN={};

			    for h=1:Ndelays

			        MSD_tot_err_ind_tracks_sameN{h}=std(displ_sameN{h})/sqrt((n_independent_sameN(i)));
			        MSD_tot_sameN{h}=mean(displ_sameN{h});

			    end 
			    obj.mobility_structures{i}.MSD_tot_err_ind_tracks_sameN=MSD_tot_err_ind_tracks_sameN;
			    obj.mobility_structures{i}.MSD_tot_sameN=MSD_tot_sameN;

			end %for m on m.s.
		end %fun




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	OPERATIONS ON MULTIPLE LOCALIZED ACTIVATION TRACKS	 methods	%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		function combine_la_classes(obj, la_array)
			obj.mobility_structures{1}.MSD_tot_combined=[];
			obj.mobility_structures{1}.MSD_tot_err_ind_tracks_combined=[];
			for i=1:length(la_array)

				obj.mobility_structures{1}.MSD_tot_combined=[obj.mobility_structures{1}.MSD_tot{:}]+[la_array{i}.mobility_structures{1}.MSD_tot{:}];
				obj.mobility_structures{1}.MSD_tot_err_ind_tracks_combined=[obj.mobility_structures{1}.MSD_tot_err_ind_tracks{:}]+[la_array{i}.mobility_structures{i}.MSD_tot_err_ind_tracks{:}];

			end

			obj.mobility_structures{1}.MSD_tot=obj.mobility_structures{1}.MSD_tot_combined/(double(i+1));
			obj.mobility_structures{1}.MSD_tot_err_ind_tracks=obj.mobility_structures{1}.MSD_tot_err_ind_tracks_combined/double(i+1);

		end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	PLOTTING	 methods	%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		function [Ds_tot,Ds_bins]=plot_Ds(obj)
%-------------------------------------------------------------------------------------------------------------------------------
%compute the diffusion coefficients before and after stimulus

			for i=1:length(obj.mobility_structures)
				S{i}=obj.mobility_structures{i}.name
				labels{i}=obj.mobility_structures{i}.name
			end
			[selection] = listdlg('PromptString','Select which mobility structure you want to visualize','ListString',S)
			cc=lines(length(S)+30); %nice colormap
			Ds_tot=cell(1,length(S));
			Ds_tot_err=cell(1,length(S));			

			for i=1:length(S)
				i
				tempDs=[];
				tempDs_err=[];
				for j=1:length(obj.mobility_structures{i}.infotracks)
					tempDs=[obj.mobility_structures{i}.infotracks{j}.D_ST];
					tempDs=tempDs(tempDs > 0 & ~isnan(tempDs)  );
					min(tempDs)
					Ds_tot{i}=horzcat(Ds_tot{i},tempDs);
					tempDs_err=[obj.mobility_structures{i}.infotracks{j}.D_err];
					Ds_tot_err{i}=horzcat(Ds_tot_err{i},tempDs_err);
				end
				Ds_tot{i}(isnan(Ds_tot{i}))=[];
				Ds_tot{i}(Ds_tot{i}==0)=[];
				Ds_tot{i}(Ds_tot{i}<0)=[];
			end



			Ds_bins=(0:0.03:1.2);			

%-------------------------------------------------------------------------------------------------------------------------------
	%plot the histograms of the diffusion coefficients before and after stimulus


			figure
			box
			t=0;
			for i=selection

				hold on	


				[f,t]=Histo_setted_rep2(Ds_tot{i},Ds_bins,cc(i,:),'p','y',i,t);	
				figure

				[freq,bincenters]=hist(Ds_tot{i},Ds_bins);
			    norm=sum(freq);
    			step=diff(Ds_bins(end-1:end));
    			freqnorm=freq/double(norm)/step;

    			bar(Ds_bins,freqnorm)

				set(gca,'fontsize',20)
				ylabel('Frequency','fontsize',20)
				xlabel('D [um^2/s]','fontsize',20)	
				%xlim([Ds_bins(1) Ds_bins(end)])		
				title('single-track diffusion coefficient')
				set(gca,'linewidth',2)

				

			end		
		end %plot_Ds

		function []=plot_results(obj)

			save_results_path=uigetdir('D:/ANALISYS/2015/1.CLUSTER_OF_RECEPTORS_TRACKING/','choose the folder where to save plots of results');
%-------------------------------------------------------------------------------------------------------------------------------
%plot parameters and initializatoin
			pixelsize=106;
			%pixelsize=64;
			timedelay=0.2
			delay=4
			S={}
			bins=(0:15:500)
			Ndelays=10;
%get the  list of the names of each mobility structure contained in the class
%-------------------------------------------------------------------------------------------------------------------------------
			for i=1:length(obj.mobility_structures)
				S{i}=obj.mobility_structures{i}.name
				labels{i}=obj.mobility_structures{i}.name
			end
			[selection] = listdlg('PromptString','Select which mobility structure you want to visualize','ListString',S)
			%h=myGUI() gui to select which things to plot

%-------------------------------------------------------------------------------------------------------------------------------			
%plot global MSD of each mobility structure selceted

			figure
			hold on		

			set(gca,'fontsize',20)
			xlabel('time delay \Delta t [s]','fontsize',20)
			ylabel('MSD [um^2]','fontsize',20)
			cc=lines(length(S)+30); %nice colormap
			box
			set(gca,'LineWidth',2)

			pFs={};
			handles=[];
			fitDs=[];
			fitErr={};

			for i=selection

				time=(1:length(cell2mat(obj.mobility_structures{i}.MSD_tot)))*timedelay;		

			%create the plot handles
				handles=vertcat( handles,shadedErrorBar(time,1e-6*cell2mat(obj.mobility_structures{i}.MSD_tot), ...
					1e-6*cell2mat(obj.mobility_structures{i}.MSD_tot_err_ind_tracks),{'color',cc(i,:)},1) )
			%overlap datapoints
				hold on
				plot(time,1e-6*cell2mat(obj.mobility_structures{i}.MSD_tot),'o','MarkerFaceColor',cc(i,:),'MarkerEdgeColor','k')
			%fit the msd
				[pF,r]=fitLinearMSD(time,cell2mat(obj.mobility_structures{i}.MSD_tot),time,cell2mat(obj.mobility_structures{i}.MSD_tot));

				%plot(time,pF(1)*time+pF(2),'k')

				fitErr{i}=sqrt(diag(inv(r.R)*inv(r.R)')*(r.normr)^2/r.df);
			
				fitErr{i}=num2str(round(fitErr{i}(1)*10000)/10);

				appD=num2str(round(pF(1)/4.*10000)/10);

				pFs{i}=['/bf{', labels{i} ', D =  ' '(' appD setstr(177) fitErr{i} ')' ' ^{.}10^{-3} um^2/s}'];
				fitDs(i)=[pF(1)/4.];
			end %for i

			for k=1:length(pFs)
				if isempty(pFs{k})
					pFs{k}='';
				end
			end

			l=legend([handles.mainLine],pFs,'Location','best','fontsize', 15 );
			%set(l,'interpreter','tex')
			title('global MSD')

			set(gca, 'XScale','log');
			set(gca, 'YScale','log');
			xlim([0.05 1]);
			ylim([0.004 0.035]);

			saveas(gcf,[save_results_path '/global MSD.pdf'])
			saveas(gcf,[save_results_path '/global MSD.fig'])	

%-------------------------------------------------------------------------------------------------------------------------------
	%makes bar plot of the fitted diff constants

			figure
			h=[];
			
			for i=selection			
				h(i)=bar(i,fitDs(i));
				hold on
				set(h(i),'facecolor',cc(i,:))

				errD(i)=str2num(fitErr{i})*1e-3;

			end

			errorbar((1:length(fitErr)),fitDs,errD,'k','linestyle','none')

			set(gca,'fontsize',20)
			set(gca,'XTickLabel',labels,'XTick',1:numel(labels))
			set(gcf,'position',[560   457   800   491])
			ylim([0,0.015])
			ylabel('Apparent Diff Coeff [ \mum^2/s ]','fontsize',25)			
			
			saveas(gcf,[save_results_path '/Diff coeff barplot.pdf'])
			saveas(gcf,[save_results_path '/Diff coeff barplot.fig'])


%-------------------------------------------------------------------------------------------------------------------------------
	%PLOT THE DISTRIBUTION OF SINGLE-TRACK diffusion coefficients computed as in (Vestergaard 2014 PRE, sec III.B)
threshold_track_length=5
    try
		for i=selection
			obj.mobility_structures{i}.ALLsingle_track_Ds;		
		end
	catch
		[dP,dM,sP,sM,dP_perdataset,dM_perdataset,sP_perdataset,sM_perdataset,Ds_perdataset,ss_perdataset,D, D_perdataset] = computeDandS_mingus(obj,threshold_track_length,pixelsize); %compute the distriubtion of single-track D for polar and middle tracks
		
		for i=selection
			NtracksP(i)=length(dP{i});
			NtracksM(i)=length(dM{i});			
			obj.mobility_structures{i}.POLEsingle_track_D=dP{i};   %assignment
			obj.mobility_structures{i}.MIDDLEsingle_track_D=dM{i}; %assignment
			obj.mobility_structures{i}.ALLsingle_track_D=D{i}; %assignment
			obj.mobility_structures{i}.D_perdataset=D_perdataset;  %assignment

		end
	end

	NtracksPmin=min(NtracksP);
	NtracksMmin=min(NtracksM);	

	figure
	%plot the middle distribution and bar plots
	for i=selection
		dM{i}=obj.mobility_structures{i}.MIDDLEsingle_track_D(1:NtracksMmin);
		Histo_setted_rep2([dM{i}],(0:0.008:0.2),cc(i,:),'p','y',i,0);
	end
	%set(gca,'fontsize',20)
	title('middle apparent D distribution')
	ylabel('Apparent Diff Coeff [ \mum^2/s ]','fontsize',25)	

	figure
	for i=selection			
		h(i)=bar(i,mean([dM{i} dP{i}] ));
		hold on
		set(h(i),'facecolor',cc(i,:))

		errDst(i)=(std([dM{i} dP{i}])/sqrt(length(dM{i})));
		errorbar(i,mean([dM{i} dP{i}] ),errDst(i),'k','linestyle','none')
	end
	title('global apparent D from single-tracks')
	set(gca,'XTickLabel',labels,'XTick',1:numel(labels))
	ylabel('Apparent Diff Coeff [ \mum^2/s ]','fontsize',25)	
			saveas(gcf,[save_results_path '/SINGLETRACK Diff coeff barplot GLOBAL.pdf'])
			saveas(gcf,[save_results_path '/SINGLETRACK Diff coeff barplot GLOBAL.fig'])
	figure
	for i=selection			
		h(i)=bar(i,mean(dM{i}));
		hold on
		set(h(i),'facecolor',cc(i,:))

		errDst(i)=(std(dM{i}))/sqrt(length(dM{i}));
		errorbar(i,mean(dM{i}),errDst(i),'k','linestyle','none')
	end
	title('middle apparent D from single-tracks')
	set(gca,'XTickLabel',labels,'XTick',1:numel(labels))
	ylabel('Apparent Diff Coeff [ \mum^2/s ]','fontsize',25)	
			saveas(gcf,[save_results_path '/SINGLETRACK Diff coeff barplot MIDDLE.pdf'])
			saveas(gcf,[save_results_path '/SINGLETRACK Diff coeff barplot MIDDLE.fig'])

	%plot the differences in ST diff coeff between consecutive datasets
	figure
	diffs=[];
	for i = selection(1:end-1)
		diffs(i)=mean([dM{i} dP{i}] ) - mean([dM{i+1} dP{i+1}] )
		hh(i)=bar(i,diffs(i));		
		hold on
		set(hh(i),'facecolor',cc(i,:))
	end		
	title('differences of D from single-tracks')
	set(gca,'fontsize',20)
	ylabel('Apparent Diff Coeff [ \mum^2/s ]','fontsize',25)	
	%plot the polar distribution of D and bar plots

	figure
	for i=selection
		dP{i}=obj.mobility_structures{i}.POLEsingle_track_D(1:NtracksPmin);
		Histo_setted_rep2([dP{i}],(0:0.002:0.04),cc(i,:),'p','y',i,0);
	end
	%set(gca,'fontsize',20)
	title('polar apparent D distribution')


	figure
	for i=selection			
		h(i)=bar(i,mean(dP{i}));
		hold on
		set(h(i),'facecolor',cc(i,:))

		errDst(i)=(std(dP{i}))/sqrt(length(dP{i}));
		errorbar(i,mean(dP{i}),errDst(i),'k','linestyle','none')
	end

	title('polar apparent D from single-tracks')
	set(gca,'XTickLabel',labels,'XTick',1:numel(labels))
	saveas(gcf,[save_results_path '/SINGLETRACK Diff coeff barplot POLAR.pdf'])
	saveas(gcf,[save_results_path '/SINGLETRACK Diff coeff barplot POLAR.fig'])	
	%set(gca,'fontsize',20)





%-------------------------------------------------------------------------------------------------------------------------------
%plot the alpha (coefficient describing type of diffusion, sub-normal-super) distribution for each condition

binsa=(0:0.15:2);
figure
	for i=selection
		Nsubdatasets=length(obj.mobility_structures{i}.infotracks);
		a{i}(1)=NaN;
		for j=1:Nsubdatasets
			a{i}=[a{i},obj.mobility_structures{i}.infotracks{j}.a];

			%divide into polar and middle tracks

		end
		a{i}=a{i}(~isnan(a{i}));
		a{i}=a{i}(a{i}>0.001);
		Histo_setted_rep2(a{i},binsa,cc(i,:),'p','y',i,0)	
	end
	title('alpha')

	saveas(gcf,[save_results_path '/Aplha_singletrack_distribution.pdf'])
	saveas(gcf,[save_results_path '/Aplha_singletrack_distribution.fig'])

bins=(0:20:600);

%plot total displecements
			figure
			set(gca,'fontsize',20)
			xlabel('Total Displacement [nm]','fontsize',25)
			ylabel('PDF','fontsize',25)
			cc=lines(50) %nice colormap
			t=0;
			for i=selection
				hold on	
				[f,tt]=Histo_setted_rep2(sqrt((obj.mobility_structures{i}.displ_tot{delay})),bins,cc(i,:),'p','y',i,t);	
				xlim([0,600])		

			end
			l=legend(labels,'Location','northeast','fontsize', 15 );
			set(l,'interpreter','tex');
			title(['total displacement dt= ' num2str(delay*timedelay)]);
			saveas(gcf,[save_results_path '/global displacement.pdf'])
			saveas(gcf,[save_results_path '/global displacement.fig'])			


%-------------------------------------------------------------------------------------------------------------------------------
%plot the track length distributions

			tl=cell(length(obj.subdataset),1);
			figure
			tracklength_bins=(0:1:1000);
			t=0;
			hold on

			for i=selection
				
				tl{i}=vertcat(obj.mobility_structures{i}.track_length{:})						
				[f,tt]=Histo_setted_rep2(tl{i},tracklength_bins, cc(i,:),'p','y',i,t);				
				
			end	
			title('tracklength')
			l=legend(labels,'Location','northeast','fontsize', 15 );

			saveas(gcf,[save_results_path '/tracklength.pdf'])
			saveas(gcf,[save_results_path '/tracklength.fig'])			

%-------------------------------------------------------------------------------------------------------------------------------
%plot displacements for middle-polar tracks

			
			cc=lines(length(selection)*30);
			time=(1:Ndelays)*timedelay;	
			it=cell(1,selection(end));
			for i=selection
				it{i}=obj.mobility_structures{i}.middlepoleInfo
			end

			t1=0;
			t2=0;
			for i=selection			
				hold on
				figure	
				[f1,tt1]=Histo_setted_rep2(sqrt(it{i}.middledispl_tot),bins,cc(i,:),'p','y',i,t1);
				[f2,tt2]=Histo_setted_rep2(sqrt(it{i}.poledispl_tot),bins,cc(i+1,:),'p','y',i,t2);
				title(['middle-pole-displ ' obj.mobility_structures{i}.name])
			end
			
%-------------------------------------------------------------------------------------------------------------------------------
%plot MSDs for middle-polar tracks

			pFs={};
			handles=[];
			figure
            ii=0;
			for i=selection
                ii=ii+1;
				hold on
				handles=vertcat(handles,shadedErrorBar(time,[it{i}.poleMSD{:}]*1e-6,[it{i}.poleMSD_err{:}]*1e-6,{'color',cc(i,:)},1))	;
				plot(time,1e-6*cell2mat(it{i}.poleMSD),'o','MarkerFaceColor',cc(i,:),'MarkerEdgeColor','k')
				[pF,r]=fitLinearMSD(time,[it{i}.poleMSD{:}],time,[it{i}.poleMSD{:}],[0,0.04]);
				pFs{ii}=[labels{i} ' ' '  ' num2str(pF(1)/4.) ' um^2/s'];
								
			end %for i

				title(['pole MSD ' obj.mobility_structures{i}.name])
				l=legend([handles.mainLine],pFs,'Location','best','fontsize', 15 );
				set(l,'interpreter','none')	
			set(gca,'fontsize',20)
			xlabel('time delay \Delta t [s]','fontsize',20)
			ylabel('MSD [um^2]','fontsize',20)
			cc=lines(length(S)+30); %nice colormap
			box
			set(gca,'LineWidth',2)
			set(gca, 'XScale','log');
			set(gca, 'YScale','log');
			xlim([0.05 1]);
			ylim([0.004 0.035]);
			saveas(gcf,[save_results_path '/polar MSD.pdf'])
			saveas(gcf,[save_results_path '/polar MSD.fig'])			


			pFs={};
			handles=[];
			figure	
            ii=0;
			for i=selection
				
				hold on
				ii=ii+1;
				handles=vertcat(handles,shadedErrorBar(time,[it{i}.middleMSD{:}]*1e-6,[it{i}.middleMSD_err{:}]*1e-6,{'color',cc(i,:)},1))	
				plot(time,1e-6*cell2mat(it{i}.middleMSD),'o','MarkerFaceColor',cc(i,:),'MarkerEdgeColor','k')
				[pF,r]=fitLinearMSD(time,[it{i}.middleMSD{:}],time,[it{i}.middleMSD{:}],[0,0.01]);
				pFs{ii}=[labels{i} ' ' '  ' num2str(pF(1)/4.) ' um^2/s'];
								
			end %for i

				title(['middle MSD ' obj.mobility_structures{i}.name])
				l=legend([handles.mainLine],pFs,'Location','best','fontsize', 15 );
				set(l,'interpreter','none')	
			set(gca,'fontsize',20)
			xlabel('time delay \Delta t [s]','fontsize',20)
			ylabel('MSD [um^2]','fontsize',20)
			cc=lines(length(S)+30); %nice colormap
			box
			set(gca,'LineWidth',2)
			set(gca, 'XScale','log');
			set(gca, 'YScale','log');
			xlim([0.05 1]);
			ylim([0.004 0.035]);
			saveas(gcf,[save_results_path '/middle MSD.pdf'])
			saveas(gcf,[save_results_path '/middle MSD.fig'])

			handlesM=[];
			handlesP=[];
			pFsP={};
			pFsM={};





%{



%-------------------------------------------------------------------------------------------------------------------------------
	%plot the MSD for each dataset 

			figure('position', [50 50 750 900])

			for i=selection
				
				subplot('position', [ 0.1 (0.1 +(i-1)*(1/length(selection)))*0.9 (1/length(selection))*1.2 (1/length(selection))*0.6])
				hold on
				set(gca,'fontsize',10)
				xlabel('time delay \Delta t [s]','fontsize',10)
				ylabel('MSD [nm^2]','fontsize',10)	

				handles=[]	
				pFs={}
				for j=1:length([obj.mobility_structures{i}.MSD_perdataset{:,1}])
											
					time=(1:length(([obj.mobility_structures{i}.MSD_perdataset{j,:}])))*timedelay;	

				%create the plot handles
					handles=vertcat(handles,shadedErrorBar(time,[obj.mobility_structures{i}.MSD_perdataset{j,:}], ... 
						[obj.mobility_structures{i}.MSD_perdataset_err{j,:}],{'color',cc(j,:)},1))					
			
				%fit the msd					
					[pF,r]=fitLinearMSD(time,[obj.mobility_structures{i}.MSD_perdataset{j,:}],time,[obj.mobility_structures{i}.MSD_perdataset{j,:}],[0,0.04])
					pFs{j}=[labels{i} ' ' num2str(j) '  ' num2str(pF(1)/4.) ' um^2/s']
								
				end %for j

				title('ensemble-time averaged MSD per dataset')
				l=legend([handles.mainLine],pFs,'Location',[(1/2) (0.1 +(i-1)*(1/length(selection)))*0.9 (1/length(selection))*1.2 (1/length(selection))*0.6],'fontsize', 15 )
				%set(l,'interpreter','tex')	
				
			end %for i
%-------------------------------------------------------------------------------------------------------------------------------
	%plot the MSD for each dataset and each subdataset in the same figure (to see behaviour over time)
			
			figure		
			k=0;
			set(gca,'fontsize',20)				
			xlabel('time delay \Delta t [s]','fontsize',15)
			ylabel('MSD [nm^2]','fontsize',15)					
			hold on

			handles=[]	
			pFs={}	
			for i=selection				

				for j=1:length([obj.mobility_structures{i}.MSD_perdataset{:,1}])
					k=k+1;		
					time=(1:length(([obj.mobility_structures{i}.MSD_perdataset{j,:}])))*timedelay;		

				%create the plot handles				
					handles=vertcat(handles,errorbar(time,[obj.mobility_structures{i}.MSD_perdataset{j,:}], ... 
						[obj.mobility_structures{i}.MSD_perdataset_err{j,:}],'color',cc(k,:),'linewidth',2))					

				%fit the msd
						
					[pF,r]=fitLinearMSD(time,[obj.mobility_structures{i}.MSD_perdataset{j,:}],time,[obj.mobility_structures{i}.MSD_perdataset{j,:}],[0,0.04])
					pFs{k}=[num2str(pF(1)/4.) ' um^2/s']
					labels2{k}=['dataset' num2str(k) ,pFs{k}]
					
				end %for j

				l=legend(labels2,'Location','best','fontsize', 10 )
				set(l,'interpreter','tex')	
				ylim([0 4*1e4])												

			end %for i
				title('global time-averaged MSD for ALL dataset')


			saveas(gcf,[save_results_path '/all MSDs.pdf'])
			saveas(gcf,[save_results_path '/all MSDs.fig'])	


%-------------------------------------------------------------------------------------------------------------------------------
%plot global MSD of each FOV of each subdataset for each condition (e.g.FOV1-nostim-stim-remov, FOV2-nostim-stim-remov...)


			figure('position', [50 50 750 900])

			L=length([obj.mobility_structures{1}.MSD_perdataset{:,1}])
			for j=1:length([obj.mobility_structures{1}.MSD_perdataset{:,1}])

				subplot('position', [ 0.1 (0.1 +(j-1)*(1/L))*0.9 (1/L)*1.2 (1/L)*0.6])
				set(gca,'fontsize',10)
				xlabel('time delay \Delta t [s]','fontsize',10)
				ylabel('MSD [nm^2]','fontsize',10)				
				hold on	

				handles=[]	
				pFs={}
				for i=selection
											
					time=(1:10)*timedelay;	
					try
						handles=vertcat(handles,shadedErrorBar(time,[obj.mobility_structures{i}.MSD_perdataset{j,:}], ... 
								[obj.mobility_structures{i}.MSD_perdataset_err{j,:}],{'color',cc(i,:)},1))					
			

						
						%fit the msd

						[pF,r]=fitLinearMSD(time,[obj.mobility_structures{i}.MSD_perdataset{j,:}],time,[obj.mobility_structures{i}.MSD_perdataset{j,:}],[0,0.04])
						pFs{i}=[labels{i} ' '  '  ' num2str(pF(1)/4.) ' um^2/s']
					catch
					end			
				end
				for k=1:length(pFs)
					if isempty(pFs{k})
						pFs{k}='';
					end
				end
				title(['FOV' num2str(j)])
				l=legend([handles.mainLine],pFs,'Location',[(1/2) (0.1 +(j-1)*(1/L))*0.9 (1/L)*1.2 (1/L)*0.6],'fontsize', 15 )
				set(l,'interpreter','tex')	
			end			
	
			saveas(gcf,[save_results_path '/MSD per FOV.pdf'])
			saveas(gcf,[save_results_path '/MSD per FOV.fig'])	

%-------------------------------------------------------------------------------------------------------------------------------
%plot mean of the MSDs of the same subdataset for each FOV (e.g. "(nostimulus_1 + stimulus_1 + recover_1)/3" and "(nostimulus_2 + stimulus_2 + recover_2)/3" etc)
%this is useful to get the time-dependent behavior (e.g. the poly-L-lysine transient effect over time)
			figure
			handles=[]	
			pFs={}
			for j=1:length([obj.mobility_structures{1}.MSD_perdataset{:,1}])


				av_MSD_per_ith_ds=zeros(10,1)';
				av_MSD_per_ith_ds_err=zeros(10,1)';

				for i=selection	
                    try
                        av_MSD_per_ith_ds=(av_MSD_per_ith_ds+[obj.mobility_structures{i}.MSD_perdataset{j,:}])
                        av_MSD_per_ith_ds_err=(av_MSD_per_ith_ds_err+[obj.mobility_structures{i}.MSD_tot_err_ind_tracks{:}])
                    catch
                        'ops'
                    end
				end					

				av_MSD_per_ith_ds=av_MSD_per_ith_ds/double(i);
				av_MSD_per_ith_ds_err=av_MSD_per_ith_ds_err/double(i);

				time=(1:Ndelays)*timedelay;	
				
				handles=vertcat(handles,shadedErrorBar(time,av_MSD_per_ith_ds,av_MSD_per_ith_ds_err,{'color',cc(j,:)},1))					

				set(gca,'fontsize',20)
				xlabel('time delay \Delta t [s]','fontsize',15)
				ylabel('MSD [nm^2]','fontsize',15)

				%fit the msd
				hold on	
				[pF,r]=fitLinearMSD(time,av_MSD_per_ith_ds,time,av_MSD_per_ith_ds,[0,0.04])
				pFs{j}=['FOV' num2str(j) ' ' num2str(pF(1)/4.) 'um^2/s']
							

	
			end

			for k=1:length(pFs)
				if isempty(pFs{k})
					pFs{k}='';
				end
			end
			title('mean MSD per FOV for all conditions')
			l=legend([handles.mainLine],pFs,'Location','best','fontsize', 10 )
			set(l,'interpreter','tex')







%-------------------------------------------------------------------------------------------------------------------------------

%plot fitted numbers of photons of each mobility structure selceted 

			try 		
				figure
			    cc=lines(length(S)) %nice colormap
			    
			    t=0;
				for i=selection
					for j=1:length(obj.mobility_structures{i}.tracks3D)
						AMPs{i}=vertcat(obj.mobility_structures{i}.tracks3D{j}.AMPs);


					end
					hold on
					amps=[AMPs{i}];
					amps(amps>5000)=[];
					amps(amps<0)=[];
					[f,t]=Histo_setted_rep2(amps,(0:0.001:1),cc(i,:),'p','y',i,t);
				end			

			catch
				'this mobility structures does not have 3D tracks!'
			end	
			




%plot amps of each mobility structure selceted from utrack

	
			
			figure
			cc=lines(length(S)) %nice colormap
			AMPs={}
			for i=selection
				AMPs{i}=[0]
				for j=1:length(obj.mobility_structures{i}.infotracks)					
					for k=1:length(obj.mobility_structures{i}.infotracks{j})
						AMPs{i}=horzcat(AMPs{i},[obj.mobility_structures{i}.infotracks{j}(k).tracksCoordAmpCG(4:8:end)]);
					end
				end
				hold on
				amps=[AMPs{i}];
				amps(amps>5000)=[];
				amps(amps<0)=[];
				Histo_setted_rep2(amps,(0:0.0002:0.01),cc(i,:),'p','y',0,0,0);
			end

			saveas(gcf,[save_results_path '/Fitted Amplitude Histogram.pdf'])
			saveas(gcf,[save_results_path '/Fitted Amplitude Histogram.fig'])			
%}
%-------------------------------------------------------------------------------------------------------------------------------

%{

%-------------------------------------------------------------------------------------------------------------------------------
	%plot the cumulative distr func of the displacement
			
			figure
			hold on
			for i=selection
				cdf=[];			

				[displCDF,n]=histc(sqrt(obj.mobility_structures{i}.displ_tot{delay}),bins);
				displCDF_norm=sum(displCDF);

				for j=1:length(displCDF)
					cdf(j)=sum(displCDF(1:j))/double(displCDF_norm);
				end

				plot(cdf,'color',cc(i,:),'LineWidth',4)
			end

			set(gca,'xscale','log')
			l=legend(labels,'Location','northeast','fontsize', 10 );
			xlabel('Displacement [nm]','fontsize',15)
			ylabel('CDF','fontsize',15)
			set(gca,'fontsize',15)
			title('CDF of displacement')

%-------------------------------------------------------------------------------------------------------------------------------
	%plot the cumulative distr func of the singketrack displacement
			
			figure
			hold on
			for i=selection
				cdf=[];			

				[displCDF,n]=histc(sqrt(obj.mobility_structures{i}.single_track_av_disp_tot{delay}),bins);
				displCDF_norm=sum(displCDF);

				for j=1:length(displCDF)
					cdf(j)=sum(displCDF(1:j))/double(displCDF_norm);
				end

				plot(cdf,'color',cc(i,:),'LineWidth',4)
			end
			title('singletrack displ CDF')
			set(gca,'xscale','log')
			l=legend(labels,'Location','northeast','fontsize', 10 );
			xlabel('Displacement [nm]','fontsize',15)
			ylabel('CDF','fontsize',15)
			set(gca,'fontsize',15)
			title(['CDF singletrack av disp dt=' num2str(timedelay*delay)])		


%-------------------------------------------------------------------------------------------------------------------------------
%plot MSDs for middle-polar tracks on the same figure
			figure
			for i=selection
				
				hold on
				
				handlesM=vertcat(handlesM,shadedErrorBar(time,[it{i}.middleMSD{:}],[it{i}.middleMSD_err{:}],{'color',cc(i,:)},1))	

				[pF,r]=fitLinearMSD(time,[it{i}.middleMSD{:}],time,[it{i}.middleMSD{:}],[0,0.04])
				pFsM{i}=[labels{i} ' ' '  ' num2str(pF(1)/4.) ' um^2/s']
				
				handlesP=vertcat(handlesP,shadedErrorBar(time,[it{i}.poleMSD{:}],[it{i}.poleMSD_err{:}],{'color',cc(i,:)},1))	
				[pF,r]=fitLinearMSD(time,[it{i}.poleMSD{:}],time,[it{i}.poleMSD{:}],[0,0.04])
				pFsP{i}=[labels{i} ' ' '  ' num2str(pF(1)/4.) ' um^2/s']				
								
			end %for i

				title(['middle-polar MSD ' obj.mobility_structures{i}.name])
				l=legend([handlesM.mainLine handlesP.mainLine],[pFsM pFsP],'Location','best','fontsize', 10 )
				set(l,'interpreter','none')	

			saveas(gcf,[save_results_path '/middle and polar MSD.pdf'])
			saveas(gcf,[save_results_path '/middle and polar MSD.fig'])



%-------------------------------------------------------------------------------------------------------------------------------
	%compute the diffusion coefficients before and after stimulus

			Ds_tot=cell(1,length(S));
			Ds_tot_err=cell(1,length(S));			

			for i=1:length(S)
				tempDs=[];
				tempDs_err=[];
				for j=1:length(obj.mobility_structures{i}.infotracks)
					tempDs=[obj.mobility_structures{i}.infotracks{j}.D];
					Ds_tot{i}=horzcat(Ds_tot{i},tempDs);
					tempDs_err=[obj.mobility_structures{i}.infotracks{j}.D_err];
					Ds_tot_err{i}=horzcat(Ds_tot_err{i},tempDs_err);
				end
				Ds_tot{i}=Ds_tot{i}(Ds_tot_err{i}<0.5 & Ds_tot{i}>0); %filter out bad fits
			end



			Ds_bins=(0:0.02:1);			

%-------------------------------------------------------------------------------------------------------------------------------
	%plot the histograms of the diffusion coefficients before and after stimulus


			figure

			t=0;
			for i=selection

				hold on	


				[f,t]=Histo_setted_rep2(Ds_tot{i},Ds_bins,cc(i,:),'p','y',i,t)	
			

			end			

			saveas(gcf,[save_results_path '/Diffusion Constant Histogram.pdf'])
			saveas(gcf,[save_results_path '/Diffusion Constant Histogram.fig'])

%-------------------------------------------------------------------------------------------------------------------------------
	%plot the cumulative distr func of the diffusion coefficients before and after stimulus


			figure
			for i=1:length(S)
				
			

				[cdf,x]=ecdf(Ds_tot{i}/4.);
				plot(x,cdf,'color',cc(i,:),'linewidth',2.5)

				hold on
			end
			title('Apparent Diffusion Coefficient','fontsize',20)
			l=legend(labels,'Location','best','fontsize', 15 );
			xlabel('Apparent Diffusion Coeff [um^2/s]','fontsize',25)
			ylabel('CDF','fontsize',25)
			set(gca,'fontsize',25,'XScale','log')
			title('CDF of apparent D')

%}







%{

%-------------------------------------------------------------------------------------------------------------------------------
%plot only the displacement comeing from tracks thresholded with diff coeff

			cc=lines(length(S)) %nice colormap
			figure
			
			threshold=0.01
			displ=cell(1,length(obj.subdataset));
			

			for i=selection
				for j=1:length(obj.mobility_structures{i}.infotracks)
					for k=1:length(obj.mobility_structures{i}.infotracks{j})
						position=obj.mobility_structures{i}.infotracks{j}(k).position;
						if isnan(position)
							continue
						end
						if (obj.mobility_structures{i}.infotracks{j}(k).D > threshold) && (position(1)=='r' || position(1)=='l') 
							x=obj.mobility_structures{i}.infotracks{j}(k).x*pixelsize;
							y=obj.mobility_structures{i}.infotracks{j}(k).y*pixelsize;
							xs=circshift(x,[0,delay]);
							ys=circshift(y,[0,delay]);
							tempdisp=(x(1:end-1)-xs(2:end)).^2+(y(1:end-1)-ys(2:end)).^2;
							displ{i}=horzcat(displ{i},tempdisp);
							
						end
					end
				end

				Histo_setted_rep(sqrt(displ{i}), cc(i,:),bins,'y');
				title('polar displ')
				
			end
			
			displ=cell(1,length(obj.subdataset));
			l=legend(labels,'Location','northeast','fontsize', 10 );

			figure
			for i=selection
				for j=1:length(obj.mobility_structures{i}.infotracks)
					for k=1:length(obj.mobility_structures{i}.infotracks{j})
						position=obj.mobility_structures{i}.infotracks{j}(k).position;
						if isnan(position)
							continue
						end
						if (obj.mobility_structures{i}.infotracks{j}(k).D > threshold) && (position(1)=='m' ) 
							x=obj.mobility_structures{i}.infotracks{j}(k).x*pixelsize;
							y=obj.mobility_structures{i}.infotracks{j}(k).y*pixelsize;
							xs=circshift(x,[0,delay]);
							ys=circshift(y,[0,delay]);
							tempdisp=(x(1:end-1)-xs(2:end)).^2+(y(1:end-1)-ys(2:end)).^2;
							displ{i}=horzcat(displ{i},tempdisp);
							
						end
					end
				end

				Histo_setted_rep(sqrt(displ{i}), cc(i,:),bins,'y');
				title('middle displ')
				
			end

			displ=cell(1,length(obj.subdataset));			
			l=legend(labels,'Location','northeast','fontsize', 10 );
			figure
			for i=1:length(obj.subdataset)
				for j=1:length(obj.mobility_structures{i}.infotracks)
					for k=1:length(obj.mobility_structures{i}.infotracks{j})
						position=obj.mobility_structures{i}.infotracks{j}(k).position;
						if isnan(position)
							continue
						end
						if ((obj.mobility_structures{i}.infotracks{j}(k).D > threshold) && (position(1)=='s'  ))
							x=obj.mobility_structures{i}.infotracks{j}(k).x*pixelsize;
							y=obj.mobility_structures{i}.infotracks{j}(k).y*pixelsize;
							xs=circshift(x,[0,delay]);
							ys=circshift(y,[0,delay]);
							tempdisp=(x(1:end-1)-xs(2:end)).^2+(y(1:end-1)-ys(2:end)).^2;
							displ{i}=horzcat(displ{i},tempdisp);
							
						end
					end
				end

				Histo_setted_rep(sqrt(displ{i}), cc(i,:),bins,'y');
				title('spurius displ')
				
			end	
			displ=cell(1,length(obj.subdataset));			
			l=legend(labels,'Location','northeast','fontsize', 10 );	
			figure
			for i=1:length(obj.subdataset)
				for j=1:length(obj.mobility_structures{i}.infotracks)
					for k=1:length(obj.mobility_structures{i}.infotracks{j})
						position=obj.mobility_structures{i}.infotracks{j}(k).position;
						if isnan(position)
							continue
						end
						if (obj.mobility_structures{i}.infotracks{j}(k).D < threshold)
							x=obj.mobility_structures{i}.infotracks{j}(k).x*pixelsize;
							y=obj.mobility_structures{i}.infotracks{j}(k).y*pixelsize;
							xs=circshift(x,[0,delay]);
							ys=circshift(y,[0,delay]);
							tempdisp=(x(1:end-1)-xs(2:end)).^2+(y(1:end-1)-ys(2:end)).^2;
							displ{i}=horzcat(displ{i},tempdisp);
							
						end
					end
				end

				Histo_setted_rep(sqrt(displ{i}), cc(i,:),bins,'y');
				title(['total thresholded displ, thr= ' num2str(threshold) ] )
				
			end	
			l=legend(labels,'Location','northeast','fontsize', 10 );		




%-------------------------------------------------------------------------------------------------------------------------------
			%plot thresholded diffusion coefficients
			figure
			threshold=0
			cc=lines(length(S)) %nice colormap
3
			for i=1:length(obj.subdataset)
				for j=1:length(obj.mobility_structures{i}.infotracks)
					tempDs=[obj.mobility_structures{i}.infotracks{j}.D];
					Ds{i}=horzcat(Ds{i},tempDs);
				end
			end
			Ds_tot=Ds;
			title('is it useful 1')

			figure

			bins2=(threshold:0.0025:0.04);
			%bins=(0:.0002:threshold)
			for i=1:length(obj.subdataset)
				data=Ds{i};
				data=data/4.; % MSD=2dDt
				data=data(data>threshold);
				%data=data(data<threshold);
				data=data(data>0)
				data=data(~isinf(data));
				data=data(~isnan(data));				
				f=Histo_setted_rep(data, cc(i,:),bins2,'y');
			end
			title('diffusion coefficient for different stimuli')
			l=legend(labels,'Location','northeast','fontsize', 10 );
			text(3/4.*bins(end),max(f), ['threshold = ' num2str(threshold)])

%-------------------------------------------------------------------------------------------------------------------------------
%plot the diff constant coming from the pole for hte diffretnt stimuls conditions
			
			figure
			subplot(1,3,2)
			Ds=[cell(1,length(obj.subdataset));			]
			for i=1:length(obj.subdataset)
				for j=1:length(obj.mobility_structures{i}.infotracks)
					for k=1:length(obj.mobility_structures{i}.infotracks{j})
						position=obj.mobility_structures{i}.infotracks{j}(k).position;
						if isnan(position)
							continue
						end
						if (position(1)=='r' || position(1)=='l') 
							tempDs=[obj.mobility_structures{i}.infotracks{j}(k).D];
							Ds{i}=horzcat(Ds{i},tempDs);
							
						else

						end
					end
				end				

			end
			for i=1:length(obj.subdataset)
				data=Ds{i};
				data=data/4.; % MSD=2dDt
				data=data(data>0)
				data=data(~isinf(data));
				data=data(~isnan(data));				
				f=Histo_setted_rep(data, cc(i,:),bins2,'y');
			end
			title('D OF POLE for different stimuli')
			l=legend(labels,'Location','northeast','fontsize', 10 );
			text(3/4.*bins(end),max(f), ['threshold = ' num2str(threshold)])

%-------------------------------------------------------------------------------------------------------------------------------
%plot the diff constant coming from the MIDDLECELL for hte diffretnt stimuls conditions
	
			subplot(1,3,3)
			Ds=[cell(1,length(obj.subdataset))];
			for i=1:length(obj.subdataset)
				for j=1:length(obj.mobility_structures{i}.infotracks)
					for k=1:length(obj.mobility_structures{i}.infotracks{j})
						position=obj.mobility_structures{i}.infotracks{j}(k).position;
						if isnan(position)
							continue
						end
						if (position(1)=='m')
							tempDs=[obj.mobility_structures{i}.infotracks{j}(k).D];
							Ds{i}=horzcat(Ds{i},tempDs);
							
						else
							
						end
					end
				end				

			end
			figure
			for i=1:length(obj.subdataset)
				data=Ds{i};
				data=data/4.; % MSD=2dDt
				data=data(data>0)
				data=data(~isinf(data));
				data=data(~isnan(data));				
				f=Histo_setted_rep(data, cc(i,:),bins2,'y');
			end
			title('D OF MIDDLE for different stimuli')
			l=legend(labels,'Location','northeast','fontsize', 10 );
			text(3/4.*bins(end),max(f), ['threshold = ' num2str(threshold)])






%-------------------------------------------------------------------------------------------------------------------------------

	%plot the amplitude distribution coming from the tracks immediately after activation		
			figure
			amp_tot_per_ds=cell(length(selection),1);
			Dtot=[]
			amp_tot_first=[]

			for i=1:length(selection)
				for j=1:length(obj.mobility_structures{i}.tracks3D)
						a=filter_tracks_noNAN_rep(obj.mobility_structures{i}.tracks3D{j},10,200)					
						for k=1:length(a)%(obj.mobility_structures{i}.tracks3D{j})

								i
								j
								a(k)
								if ~isempty(a(k).AMPs) && ~isempty(a(k).D)
									amp_tot_per_ds{i}=vertcat(amp_tot_per_ds{i},a(k).AMPs(1:5));
									Dtot=vertcat(Dtot,obj.mobility_structures{i}.tracks3D{j}(k).D);
									amp_tot_first=vertcat(amp_tot_first,mean(a(k).AMPs(1:5)));
								end


								if ~isempty(obj.mobility_structures{i}.tracks3D{j}(k).AMPs) && ~isempty(obj.mobility_structures{i}.tracks3D{j}(k).D)
									amp_tot_per_ds{i}=vertcat(amp_tot_per_ds{i},obj.mobility_structures{i}.tracks3D{j}(k).AMPs(1:3));
									Dtot=vertcat(Dtot,obj.mobility_structures{i}.tracks3D{j}(k).D);
									amp_tot_first=vertcat(amp_tot_first,mean(obj.mobility_structures{i}.tracks3D{j}(k).AMPs(1:5)));
								end

									

							
							
						%amp_tot_per_ds{i}=horzcat(amp_tot_per_ds{i},obj.mobility_structures{i}.infotracks{j}(k).tracksCoordAmpCG(4:8:end));
						end
				


				
				end



				binsamp=(1:75:5500);
				amp_tot_first;
				amp_tot_per_ds{i}(amp_tot_per_ds{i}<0)=[];
				amp_tot_per_ds{i}(amp_tot_per_ds{i}>5500)=[];
				amp_tot_per_ds{i}(isnan(amp_tot_per_ds{i}))=[];
				Histo_setted_rep(amp_tot_per_ds{i}, cc(i,:),binsamp,0.008*0.03*i,'Photons');	
				title('porcoddio')
			end	

			


			length(amp_tot_first)
			length(Dtot)

			Dtot=Dtot/4;
			cond=amp_tot_first<0;
			Dtot(amp_tot_first>5000)=[];
			Dtot(cond)=[];
			amp_tot_first(amp_tot_first	>5000)=[];
			amp_tot_first(amp_tot_first	<0)=[];

			
			amp_tot_first(Dtot>0.5)=[];
			Dtot(Dtot>0.5)=[];
			
			

			
			polyfit(amp_tot_first,Dtot,1)



			
			figure
			scatter(amp_tot_first,Dtot,'.')	
%			title('Fitted Number of Photons')
			set(gca,'fontsize',15)
			l=legend(labels,'Location','northeast','fontsize', 10 );
			xlim([0 5000])	
			xlabel('#photons')
			ylabel('D [um2/s]')



			figure
			scatter_bins=(min(amp_tot_first):(max(amp_tot_first)-min(amp_tot_first))/30:max(amp_tot_first))
			[binned_photons,err,n]=binned_scatter_plot_rep(amp_tot_first,Dtot,scatter_bins)	
			step=scatter_bins(3)-scatter_bins(2);
			shadedErrorBar(scatter_bins(1:end-1)+step/2,binned_photons,err)
			title('photons - D correlation')
			set(gca,'fontsize',15)
			l=legend(labels,'Location','northeast','fontsize', 10 );	
			xlabel('#photons')
			ylabel('D [um2/s]')
			figure
			plot(n)


%}		





		end %plotting
		

%-------------------------------------------------------------------------------------------------------------------------------
%plot the percentage of polar and middle tracks
		function plot_percentage_middlepolar(obj)
			L=length(obj.mobility_structures);
			polar_tracks=cell(1,L);
			middle_tracks=cell(1,L);
			spurius_tracks=cell(1,L);
			hold on
			figure
			for i=1:L
				labels{i}=obj.mobility_structures{i}.name;
				polar_tracks{i}=0;
				middle_tracks{i}=0;
				spurius_tracks{i}=0;
				for j=1:length(obj.mobility_structures{i}.infotracks)
					Ntracks=length(obj.mobility_structures{i}.infotracks{j});
					for k=1:Ntracks

						if obj.mobility_structures{i}.infotracks{j}(k).position(1) == 'm' 
							middle_tracks{i}=middle_tracks{i}+1;						
						elseif obj.mobility_structures{i}.infotracks{j}(k).position(1) == 'p' 
							polar_tracks{i}=polar_tracks{i}+1;							
						elseif obj.mobility_structures{i}.infotracks{j}(k).position(1) == 's'
							spurius_tracks{i}=spurius_tracks{i}+1;

						end

					end	
				end
				N=sum([middle_tracks{i},	polar_tracks{i},spurius_tracks{i}]);
				subplot(1,L,i)
				h=bar([middle_tracks{i}/N,polar_tracks{i}/N,spurius_tracks{i}/N]);
				set(gca,'XTickLabel',{'middle', 'polar', 'spurius'},'fontsize',25)
				percM=mat2str(middle_tracks{i}/double(middle_tracks{i}+polar_tracks{i}+spurius_tracks{i})*100);
				percP=mat2str(polar_tracks{i}/double(middle_tracks{i}+polar_tracks{i}+spurius_tracks{i})*100);
				percS=mat2str(spurius_tracks{i}/double(middle_tracks{i}+polar_tracks{i}+spurius_tracks{i})*100);
				text(0.9,middle_tracks{i} + middle_tracks{i}*0.05 ,[ percM(1:2) '%'],'fontsize',20)
				text(1.9,polar_tracks{i} + polar_tracks{i}*0.05 ,[ percP(1:2) '%'],'fontsize',20)
				text(2.9,spurius_tracks{i} + spurius_tracks{i}*0.05 ,[ percS(1:2) '%'],'fontsize',20)

				obj.mobility_structures{i}.middlepoleInfo.percP=percP;
				obj.mobility_structures{i}.middlepoleInfo.percM=percM;
				obj.mobility_structures{i}.middlepoleInfo.percS=percS;

				title(labels{i},'fontsize',20)
				ylim([0 0.65])
				%saveas(gcf,[save_results_path '/polar-middle-percentage_' labels{i} '.pdf'])
				%saveas(gcf,[save_results_path '/polar-middle-percentage_' labels{i} '.fig'])
			end
		end %plot_percentage_middle-polar



%-------------------------------------------------------------------------------------------------------------------------------
	%PLOT the correlation between the tracklength and the fitted diffusion coefficient
		function plot_TrackLength_D_correlation(obj)

			[bins,y,err,n,ll]=correlate_amp_and_D(obj,0);

			for i=1:length(obj.mobility_structures)
				S{i}=obj.mobility_structures{i}.name
				labels{i}=obj.mobility_structures{i}.name
			end
			[selection] = listdlg('PromptString','Select which mobility structure you want to visualize','ListString',S)

			cc=lines(length(y));

		    figure
		    for i=selection
		        %bins_center3{i}=bins_center3{i}(1:end)-BIN/double(2)
		        errorbar(bins{i},y{i},err{i},'.','color',cc(i,:),'markersize',35,'linewidth',2)
		        hold on
		        %plot(AMPs,AMPs*pF2(1)+pF2(2),'r')
		        xlabel('Tracklength', 'fontsize',25)
		        ylabel('D [um^2/s]', 'fontsize',25)
		        xlim([0 100])
		        ylim([0 0.04])
                legend('aaaaaaaaa')
		    end	

		    set(gca,'fontsize',20)
%{
		    figure
		    for i=selection
		        %bins_center3{i}=bins_center3{i}(1:end)-BIN/double(2)
		        errorbar(bins2{i},y2{i},err2{i},'.','color',cc(i,:),'markersize',35,'linewidth',2)
		        hold on
		        %plot(AMPs,AMPs*pF2(1)+pF2(2),'r')        
		        ylabel('D [um^2/s]', 'fontsize',25)
		        xlabel('/alpha', 'fontsize',25)
		        xlim([0 2])
		        ylim([0 0.02])
		    end	

		    set(gca,'fontsize',20)

		    figure
		    for i=selection
		        %bins_center3{i}=bins_center3{i}(1:end)-BIN/double(2)
		        errorbar(bins3{i},y3{i},err3{i},'.','color',cc(i,:),'markersize',35,'linewidth',2)
		        hold on
		        %plot(AMPs,AMPs*pF2(1)+pF2(2),'r')        
		        xlabel('Tracklength [frames]', 'fontsize',25)
		        ylabel('/alpha', 'fontsize',25)
		        xlim([0 50])
		        ylim([0 2])
		    end		

		    set(gca,'fontsize',20)

		    
		    figure
		    for i=selection
		        %bins_center3{i}=bins_center3{i}(1:end)-BIN/double(2)
		        errorbar(bins3{i},y3{i},err3{i},'o','color',cc(i,:))
		        hold on
		        %plot(AMPs,AMPs*pF2(1)+pF2(2),'r')        
		        ylabel('D [um^2/s]', 'fontsize',15)
		        xlabel('/alpha', 'fontsize',15)
		        xlim([0 2])
		        ylim([0 0.02])
		    end			    
%}
			%saveas(gcf,[save_results_path '/Tracklength_D_correlation.pdf'])
			%saveas(gcf,[save_results_path '/Tracklength_D_correlation.fig'])
		
		%-------------------------------------------------------------------------------------------------------------------------------
			%PLOT the difference in tracklength distribution

			figure
		    


		    BIN=bins{1}(2)-bins{1}(1);
		    bins=bins{1}(1:end-1);
		    bins=[bins, bins(end)+BIN];
		    mystep=ones(1,length(n{1}))*BIN;
		    
		    
		    for i=selection
		        N=double(sum(n{i}.*mystep));
		        plot(bins,[n{i}/N],'color',cc(i,:));
		        hold on
		    end

			legend(labels,'fontsize',20)
			title('number of tracks in each track-length bin')
			set(gca,'fontsize',20)
			ii=0;
			for i=selection(2:end)
				ii=ii+1;
				subplot(2,round((length(selection)-1)/2.),ii)
				plot(bins,n{i}-n{i-1},'*','color',cc(i,:))
				title([labels{i} '-' labels{i-1} 'N=' num2str(sum(n{i}))])
			end	

			%saveas(gcf,[save_results_path '/Tracklength_difference.pdf'])
			%saveas(gcf,[save_results_path '/Tracklength_difference.fig'])

		end %plot_tracklengtyh_D_correlation

%-------------------------------------------------------------------------------------------------------------------------------
%plot global singletrack displacements 	
		function plot_single_track_av_displacemnt(obj)

			delay=6;
			timedelay=0.2;
			bins=(0:10:300);
			for i=1:length(obj.mobility_structures)
				S{i}=obj.mobility_structures{i}.name
				labels{i}=obj.mobility_structures{i}.name
			end
			[selection] = listdlg('PromptString','Select which mobility structure you want to visualize','ListString',S)


			set(gca,'fontsize',20)
			xlabel('Singletrack av. Displacement [nm]','fontsize',25)
			ylabel('PDF','fontsize',25)
			cc=lines(length(S)); %nice colormap
			t=0;

			figure
			for i=selection
				hold on
				[f,tt]=Histo_setted_rep2(sqrt((obj.mobility_structures{i}.single_track_av_disp_tot{delay})),bins,cc(i,:),'p','y',i,t);				

			end

			set(gca,'fontsize',20)
			xlabel('Singletrack av. Displacement [nm]','fontsize',25)
			ylabel('PDF','fontsize',25)			

			figure
			for i = selection
				[f,tt]=Histo_setted_rep2(sqrt((obj.mobility_structures{i}.MIDDLEsingle_track_av_disp_tot{delay})),bins,cc(i,:),'p','y',i,t);	
			end
			set(gca,'fontsize',20)
			xlabel('MIDDLE Singletrack av. Displacement [nm]','fontsize',25)
			ylabel('PDF','fontsize',25)

			figure
			for i = selection
				[f,tt]=Histo_setted_rep2(sqrt((obj.mobility_structures{i}.POLARsingle_track_av_disp_tot{delay})),bins,cc(i,:),'p','y',i,t);	
			end
			set(gca,'fontsize',20)
			xlabel('POLAR Singletrack av. Displacement [nm]','fontsize',25)
			ylabel('PDF','fontsize',25)
			l=legend(labels,'Location','northeast','fontsize', 15 );
			set(l,'interpreter','tex');
			title(['singletrack av disp dt=' num2str(timedelay*delay)])			

			%saveas(gcf,[save_results_path '/average singletrack displacement.pdf'])
			%saveas(gcf,[save_results_path '/average singletrack displacement.fig'])	

		end %plot_singe_trak_av_displ

%plot all the MSD and displacements and D for the same mobilit ystructure

		function [allAppD,ST_Ds,ST_Ds_pole,ST_Ds_middle,ST_Ds_spurius,expTime,absoluteTimes,ST_Ds_std,ST_Ds_pole_std, ST_Ds_middle_std, weigths ...
			 , weigthsP,weigthsM] = plot_MSD_from_single_MS(obj,MS_index,plot_on,flag)
            
            if nargin<4
                flag=1
            end
         
			L=length([obj.mobility_structures{MS_index}.MSD_perdataset{:,1}])
			Ntot=obj.mobility_structures{MS_index}.n_independent_tracks;

			if plot_on

			end
			cc=lines(L+30); %nice colormap
			timedelay=0.2;
			pFs={};
			handles=[];
			fitDs=[];
			fitErr={};	
			allAppD=[];		
			absoluteTimes=[];
			ST_Ds=[];

			offset=0;
			firstname=obj.mobility_structures{MS_index}.files(1).name
			if firstname(1)=='.' offset=2; end  

			for i = 1:L

				msd=[obj.mobility_structures{MS_index}.MSD_perdataset{i,:}];
				msd_err=[obj.mobility_structures{MS_index}.MSD_perdataset_err{i,:}];
				time=(1:length(cell2mat(obj.mobility_structures{MS_index}.MSD_tot)))*timedelay;		

				%get the time between datasets;
				nm=obj.mobility_structures{MS_index}.files(i+offset).name



				ind=find(nm=='_')
                try 
                    ind=ind(2);
                catch
                    ''                   
                end
                
				hrs=nm(ind+1:ind+2);
				mins=nm(ind+3:ind+4);

				if i==1
					expTime{i}=0;
					RefExpTime=str2num(hrs)*60+str2num(mins);

					absoluteTimes(1)=RefExpTime;
				else
					expTime{i}=(str2num(hrs)*60+str2num(mins) - RefExpTime);
					absoluteTimes(i)=str2num(hrs)*60+str2num(mins);
				end

			%create the plot handles
				if plot_on
					handles=vertcat( handles,shadedErrorBar(time,1e-6*msd, ...
						1e-6*msd_err,{'color',cc(i,:)},1) )
				end
				%overlap datapoints
				if plot_on
					hold on
					plot(time,1e-6*msd,'o','MarkerFaceColor',cc(i,:),'MarkerEdgeColor','k')
				end
			%fit the msd
				[pF,r]=fitLinearMSD(time,(msd),time,msd);

				%plot(time,pF(1)*time+pF(2),'k')

				fitErr{i}=sqrt(diag(inv(r.R)*inv(r.R)')*(r.normr)^2/r.df);
			
				fitErr{i}=num2str(round(fitErr{i}(1)*10000)/10);

				appD=num2str(round(pF(1)/4.*10000)/10);
				allAppD=horzcat(allAppD,str2num(appD));
				pFs{i}=['/bf{', num2str(i) ', D =  ' '(' appD setstr(177) fitErr{i} ')' ' ^{.}10^{-3} um^2/s}'];
				fitDs(i)=[pF(1)/4.];

				ST_Ds(i)=mean(obj.mobility_structures{MS_index}.D_perdataset{MS_index,i});
				ST_Ds_std(i)=std(obj.mobility_structures{MS_index}.D_perdataset{MS_index,i}/sqrt(length(obj.mobility_structures{MS_index}.D_perdataset{MS_index,i})));
				ST_Ds_pole(i)=mean(obj.mobility_structures{MS_index}.dP_perdataset{MS_index,i});
				ST_Ds_pole_std(i)=std(obj.mobility_structures{MS_index}.dP_perdataset{MS_index,i}/sqrt(length(obj.mobility_structures{MS_index}.dP_perdataset{MS_index,i})));
				ST_Ds_middle(i)=mean(obj.mobility_structures{MS_index}.dM_perdataset{MS_index,i});
				ST_Ds_middle_std(i)=std(obj.mobility_structures{MS_index}.dM_perdataset{MS_index,i}/sqrt(length(obj.mobility_structures{MS_index}.dM_perdataset{MS_index,i})));

				weigths(i)=length(obj.mobility_structures{MS_index}.D_perdataset{MS_index,i})/Ntot;
				weigthsP(i)=length(obj.mobility_structures{MS_index}.dP_perdataset{MS_index,i})/Ntot;
				weigthsM(i)=length(obj.mobility_structures{MS_index}.dM_perdataset{MS_index,i})/Ntot;
				%weigths(i)=length(obj.mobility_structures{MS_index}.D_perdataset{MS_index,i})/Ntot;


				ST_Ds_spurius(i)=mean(obj.mobility_structures{MS_index}.dS_perdataset{MS_index,i});



			end %for i

			for k=1:length(pFs)
				if isempty(pFs{k})
					pFs{k}='';
				end
			end

			if plot_on
				l=legend([handles.mainLine],pFs,'Location','best','fontsize', 15 );
				%set(l,'interpreter','tex')
				title('global MSD')

				set(gca, 'XScale','log');
				set(gca, 'YScale','log');

				figure
				plot(allAppD,'-o','color',cc(MS_index,:) )
				ylim([5,20])
			end

			expTime=[expTime{:}];

		end%plotsingleMSD

%-----------------------------------------------------------------------------------------------------------------------------

		function [T,DsTot,ST_Ds,ST_Dp,ST_Dm]=plot_all_MSD_from_ALL_MS(obj,plot_on,saveOn)
              threshold_track_length=5;
              pixelsize=0.106;
			if nargin==0
				L=length(obj.mobility_structures);
				plot_on=0;
			end
			L=length(obj.mobility_structures);
			myYlim=0.04;

			if saveOn
				save_results_path=uigetdir('D:/ANALISYS/2015/1.CLUSTER_OF_RECEPTORS_TRACKING/','choose the folder where to save plots of results');
			end

			try
				for i=1:L
					obj.mobility_structures{i}.Dp_perdatasett;		
				end
			catch
				[dP,dM,sP,sM,dP_perdataset,dM_perdataset,sP_perdataset,sM_perdataset, dS_perdataset,ss_perdataset,D, D_perdataset] = computeDandS_mingus(obj,threshold_track_length,pixelsize); %compute the distriubtion of single-track D for polar and middle tracks
				for i=1:L
					NtracksP(i)=length(dP{i});
					NtracksM(i)=length(dM{i});			
					obj.mobility_structures{i}.dM_perdataset=dM_perdataset;   %assignment
					obj.mobility_structures{i}.dP_perdataset=dP_perdataset; %assignment
					obj.mobility_structures{i}.dS_perdataset=dS_perdataset; %assignment					
					obj.mobility_structures{i}.D_perdataset=D_perdataset;  %assignment
				end
			end
			
			DsTot=[];			
			if plot_on

				cc=lines(L);
				figure
				hold on

			end
			Ltot=1;

%plots the apparent D overtime for all the tracked clusters
			for i = 1 :L
				[Ds,ST_Ds,ST_Dp,ST_Dm,ST_Dsp,times,absoluteTimes,ST_Ds_std,ST_Ds_pole_std, ST_Ds_middle_std,weigths,weigthsP,weigthsM]=obj.plot_MSD_from_single_MS(i,0); %call the method that extracts Ds and times
			
				if i==1
					times0=absoluteTimes(1);
				end

				L2=length(Ds);

				DsTot=[DsTot,Ds];

				t=absoluteTimes-times0;

				%plot((Ltot:Ltot+L2-1),Ds*1e-3,'-o','color',cc(i,:),'LineWidth',3)
				%plot((Ltot:Ltot+L2-1),mean(Ds)*1e-3*ones(length((Ltot:Ltot+L2-1))),'.-','color',cc(i,:),'LineWidth',2)

				if plot_on
					errorbar(absoluteTimes-times0,ST_Ds,ST_Ds_std,'o','color',cc(i,:),'LineWidth',3)
					plot(absoluteTimes-times0,sum(ST_Ds.*weigths)/sum(weigths)*ones(length((Ltot:Ltot+L2-1))),'-','color',cc(i,:),'LineWidth',2)

					%plot(linspace(t(1),t(end),length([obj.mobility_structures{i}.D_perdataset{i,:}])),[obj.mobility_structures{i}.D_perdataset{i,:}],'.','color',cc(i,:))					
				end

				T=absoluteTimes-times0;
				Ltot=Ltot+L2;
			end

			tEnd=absoluteTimes-times0
			tEnd=tEnd(end)
			tEnd=250;

			DsTot=[];

			if plot_on
				set(gca,'fontsize',20)
				xlabel('time [mins]','fontsize',20)
				ylabel('D [um^2/s]','fontsize',20)
				title('all clusters')
				ylim([0.005 myYlim])
				xlim([0,tEnd])
				box
				set(gca,'linewidth',2)
				if saveOn
					saveas(gcf,[save_results_path '/all_clusters_ST.pdf'])
					saveas(gcf,[save_results_path '/all_clusters.fig'])	
				end				
				cc=lines(L);
				figure
				hold on
			end

			Ltot=1;

%plots the apparent D overtime from MSD of each FOV
			for i = 1 :L

				[Ds,ST_Ds,ST_Dp,ST_Dm,ST_Dsp,times,absoluteTimes,ST_Ds_std,ST_Ds_pole_std, ST_Ds_middle_std,weigths,weigthsP,weigthsM]=obj.plot_MSD_from_single_MS(i,0); %call the method that extracts Ds and times

				if i==1
					times0=absoluteTimes(1);
				end

				L2=length(Ds);

				DsTot=[DsTot,Ds];

				%plot((Ltot:Ltot+L2-1),Ds*1e-3,'-o','color',cc(i,:),'LineWidth',3)
				%plot((Ltot:Ltot+L2-1),mean(Ds)*1e-3*ones(length((Ltot:Ltot+L2-1))),'.-','color',cc(i,:),'LineWidth',2)
				t=absoluteTimes-times0;
				if plot_on
					plot(absoluteTimes-times0,Ds*1e-3,'o','color',cc(i,:),'LineWidth',3)
					plot(absoluteTimes-times0,nanmean(Ds)*1e-3*ones(length((Ltot:Ltot+L2-1))),'-','color',cc(i,:),'LineWidth',2)
					%plot(linspace(t(1),t(end),length([obj.mobility_structures{i}.D_perdataset{i,:}])),[obj.mobility_structures{i}.D_perdataset{i,:}],'.')
				end
				
				Ltot=Ltot+L2;

			end

			if plot_on
				set(gca,'fontsize',20)
				xlabel('time [mins]','fontsize',20)
				ylabel('D [um^2/s]','fontsize',20)			
				title('all clusters from MSD')
				ylim([0.005 myYlim])
				xlim([0,tEnd])
				box
				set(gca,'linewidth',2)

				cc=lines(L);
				figure
				hold on
			end

			DsTot=[];
			Ltot=1;

%plots the apparent D overtime for polar clusters
			for i = 1 :L

				[Ds,ST_Ds,ST_Dp,ST_Dm,ST_Dsp,times,absoluteTimes,ST_Ds_std,ST_Ds_pole_std, ST_Ds_middle_std,weigths,weigthsP,weigthsM]=obj.plot_MSD_from_single_MS(i,0); %call the method that extracts Ds and times
			
				if i==1
					times0=absoluteTimes(1);
				end

				L2=length(Ds);

				DsTot=[DsTot,Ds];

				%plot((Ltot:Ltot+L2-1),Ds*1e-3,'-o','color',cc(i,:),'LineWidth',3)
				%plot((Ltot:Ltot+L2-1),mean(Ds)*1e-3*ones(length((Ltot:Ltot+L2-1))),'.-','color',cc(i,:),'LineWidth',2)

				if plot_on				
					errorbar(absoluteTimes-times0,ST_Dp,ST_Ds_pole_std, 'o','color',cc(i,:),'LineWidth',3)
					plot(absoluteTimes-times0,sum(ST_Dp.*weigthsP)/sum(weigthsP)*ones(length((Ltot:Ltot+L2-1))),'-','color',cc(i,:),'LineWidth',2)
				end

				Ltot=Ltot+L2;
			end

			DsTot=[];

			if plot_on
				set(gca,'fontsize',20)
				xlabel('time [mins]','fontsize',20)
				ylabel('D [um^2/s]','fontsize',20)	
				title('polar clusters')
				ylim([0.005 myYlim])
				xlim([0,tEnd])
				box
				set(gca,'linewidth',2)
				if save_results_path
					saveas(gcf,[save_results_path '/polar_clusters_ST.pdf'])
					saveas(gcf,[save_results_path '/polar_clusters.fig'])	
				end
				cc=lines(L);
				figure
				hold on
			end

			Ltot=1;

%plots the apparent D overtime for middel clusters
		for i = 1 :L

				[Ds,ST_Ds,ST_Dp,ST_Dm,ST_Dsp,times,absoluteTimes,ST_Ds_std,ST_Ds_pole_std, ST_Ds_middle_std,weigths,weigthsP,weigthsM]=obj.plot_MSD_from_single_MS(i,0); %call the method that extracts Ds and times
				
				if i==1
					times0=absoluteTimes(1);
				end

				L2=length(Ds);

				DsTot=[DsTot,Ds];

				%plot((Ltot:Ltot+L2-1),Ds*1e-3,'-o','color',cc(i,:),'LineWidth',3)
				%plot((Ltot:Ltot+L2-1),mean(Ds)*1e-3*ones(length((Ltot:Ltot+L2-1))),'.-','color',cc(i,:),'LineWidth',2)

				if plot_on
					hold on
					errorbar(absoluteTimes-times0,ST_Dm,ST_Ds_middle_std,'o','color',cc(i,:),'LineWidth',3)
					plot(absoluteTimes-times0,sum(ST_Dm.*weigthsM)/sum(weigthsM)*ones(length((Ltot:Ltot+L2-1))),'-','color',cc(i,:),'LineWidth',2)
					
				end

				Ltot=Ltot+L2;
			end

			if plot_on
				set(gca,'fontsize',20)
				xlabel('time [mins]','fontsize',20)
				ylabel('D [um^2/s]','fontsize',20)	
				title('middle clusters')
				ylim([0.005 myYlim])		
				xlim([0,tEnd])
				box
				set(gca,'linewidth',2)

				if save_results_path
					saveas(gcf,[save_results_path '/middle_clusters_ST.pdf'])
					saveas(gcf,[save_results_path '/middle_clusters.fig'])	
				end			
				figure
				hold on
			end
			Ltot=1;

%plots the apparent D overtime for spurius clusters
			for i = 1 :L

				[Ds,ST_Ds,ST_Dp,ST_Dm,ST_Dsp,times,absoluteTimes,ST_Ds_std,ST_Ds_pole_std, ST_Ds_middle_std,weigths,weigthsP,weigthsM]=obj.plot_MSD_from_single_MS(i,0); %call the method that extracts Ds and times
				

				if i==1
					times0=absoluteTimes(1);
				end

				L2=length(Ds);
				DsTot=[DsTot,Ds];

				%plot((Ltot:Ltot+L2-1),Ds*1e-3,'-o','color',cc(i,:),'LineWidth',3)
				%plot((Ltot:Ltot+L2-1),mean(Ds)*1e-3*ones(length((Ltot:Ltot+L2-1))),'.-','color',cc(i,:),'LineWidth',2)

				if plot_on
					plot(absoluteTimes-times0,ST_Dsp,'o','color',cc(i,:),'LineWidth',3)
					plot(absoluteTimes-times0,nanmean(ST_Dsp)*ones(length((Ltot:Ltot+L2-1))),'-','color',cc(i,:),'LineWidth',2)
				end

				Ltot=Ltot+L2;
			end
			if plot_on
				set(gca,'fontsize',20)
				xlabel('time [mins]','fontsize',20)
				ylabel('D [um^2/s]','fontsize',20)	
				title('spurius clusters')
				ylim([0.005 myYlim])
				xlim([0,tEnd])		
				box	
				set(gca,'linewidth',2)					
			end

		end%plot_all_MSD_from_ALL_MS			




		function plot_TrackLength_D_correlation_2D(obj,Dthreshold, plot_On)

			if nargin < 2
				Dthreshold=0.02;
				plot_On=1;
            end

            ref_ms=1
            
			Dthreshold=0.14%*10^6;
            D_threshold=0.14%*10^6;
			TL_threshold=80;
			D_step=15;
			TL_step=15;

'diocane'
			L=length(obj.mobility_structures);

			Ds={}; %all singletrack diffusion coefficients per mobility structure
			TL={}; %all corresponding tracklenghts

			P_TL_D={};
			TL_D_axes={};

			%fills the array with tracklengths and corresponding D
			for i=1:L
                i
				Ds{i}(1)=0;
				TL{i}(1)=0;
			
				for j = obj.mobility_structures{i}.infotracks;			

					for k=[j{:}]
                        [d,s,SNR,var,a]=compute_ST_D_and_s_alt(k, 0.106, 7.5e-07);
                        k.D_ST2=d;
						if k.D_ST2 > 0 && k.D_ST2 < Dthreshold && length(k.x)<TL_threshold && length(k.x)>3 &&~isnan(k.D_ST2)
                            
							Ds{i}=[Ds{i}, k.D_ST2];
							TL{i}=[TL{i}, length(k.x)];

						end %if

					end %for k

				end %for j
                save('Ds_1020_13bins.mat','Ds')
                save('TLs_1020_13bins.mat','TL')
                bins_linearTL=linspace(1,TL_threshold,TL_step);%0:TL_threshold/TL_step:TL_threshold;
                bins_linearD=linspace(0,D_threshold,D_step);%0:D_threshold/D_step:D_threshold;
                binslogTL=logspace(log10(6),2,TL_step)
                binslogD=logspace(-3,log10(0.25),D_step)


                %[P_TL_D{i},TL_D_axes{i}]=hist3([TL{i}; Ds{i}]','edges',{binslogTL,binslogD});
                [P_TL_D{i},TL_D_axes{i}]=hist3([TL{i}; Ds{i}]','edges',{bins_linearTL,bins_linearD});
                %[P_TL_D{i},TL_D_axes{i}]=hist3([TL{i}; Ds{i}]','edges',{bins_linearTL,binslogD});
				Norm=sum(sum(P_TL_D{i}));
				P_TL_D{i}=P_TL_D{i}'/Norm; %normalize and transpose

				if plot_On  && i>ref_ms

					figure('position',[100,100,350,250])			


					%imagesc(TL_D_axes{i}{1},TL_D_axes{i}{2}, (P_TL_D{i}))%-P_TL_D{1}))
					data=P_TL_D{i}-P_TL_D{ref_ms};
					%imagesc(TL_D_axes{i}{1},TL_D_axes{i}{2}, (data))
                    h = pcolor(TL_D_axes{i}{1},TL_D_axes{i}{2}, (data));
                    h.EdgeColor = 'none';
					
                    %hold on
                    %x=(1:100);
                    %y= 7^2*0.023^2./(x*0.065*2);
                    %plot(x,y,'r','linewidth',3)
                    
                    
					colorbar
                    xlim([1,80])
                    ylim([0.00,D_threshold])
                    
                    %set(gca,'Yscale','log')
                    %set(gca,'Xscale','log')
                    
					%caxis([-0.002 0.0025])
					caxis([-0.001 0.002])
                    %caxis([-0.001 0.0025])
					axis xy

					set(gca,'fontsize',15)					
					ylabel('D_{app} [\mu m^2/2]','fontsize',30)
					xlabel('Tracklength','fontsize',30)
                    
                    x=importdata('/Users/Copo1/x.mat');
                    y=importdata('/Users/Copo1/y.mat');
                    hold on
                    plot(x,y,'r','linewidth',3)
%{					
					figure

					plot([TL_D_axes{i}{1}],(P_TL_D{i}(end,:)-P_TL_D{1}(end,:)))
					xlabel('D')

					figure
					P_TL_D_T=P_TL_D{i}'
					P_TL_D_T2=P_TL_D{i-1}'
					plot(TL_D_axes{i}{2},((P_TL_D_T(end,:)-P_TL_D_T2(end,:))))

					xlabel('TL')

%}
                    pos=[   360   271   560   420];
                    set(gcf, 'Position',pos)
                    set(gca,'FontWeight','normal','linewidth',3,'fontsize',30)
                    xlim([3.8 80])
                    ylim([0.005,0.14])
				end %if plot_on
                
            %saveas(gcf,['stim-',num2str(i),'_15bins','.pdf'])
            %saveas(gcf,['stim-',num2str(i),'_15bins', '.fig'])
			end %for15i

		end %function


%-------------------------------------------------------------------------------------------------------------------------------
		function plot_global_MSD(obj,saveOn)

			if nargin<2
				saveOn=0;
			end

%plot parameters and initializatoin
			pixelsize=106;
			%pixelsize=64;
			timedelay=0.2;
			delay=5;
			S={};
			bins=(0:20:600);
			Ndelays=10;
%get the  list of the names of each mobility structure contained in the class
%-------------------------------------------------------------------------------------------------------------------------------
			for i=1:length(obj.mobility_structures)
				S{i}=obj.mobility_structures{i}.name;
			end
			[selection] = listdlg('PromptString','Select which mobility structure you want to visualize','ListString',S)
			%h=myGUI() gui to select which things to plot

%-------------------------------------------------------------------------------------------------------------------------------			
%plot global MSD of each mobility structure selceted

			figure('position',[500,100,350,250])
			set(gca,'fontsize',25)
			box
			set(gca,'linewidth',2)
			xlabel('time delay \Delta t [s]','fontsize',15)
			ylabel('MSD [um^2]','fontsize',15)

			ylim([0.004 0.06])
			xlim([0.05 1]);
			set(gca,'XMinortick','on')

			set(gca,'YMinortick','on')
			set(gca,'XScale','log')
			set(gca,'YScale','log')
			hold on		
			%set(gcf,'position',[269 310 845 638])

			cc=lines(length(S)); %nice colormap
			box
			pFs={};
			handles=[];
			fitDs=[];
			fitErr={};
			
			time=(1:length(cell2mat(obj.mobility_structures{i}.MSD_tot)))*timedelay;

			


			for i=selection
				labels{i}=obj.mobility_structures{i}.name;
			%create the plot handles
				handles=vertcat( handles,shadedErrorBar(time,1e-6*cell2mat(obj.mobility_structures{i}.MSD_tot), ...
					1e-6*cell2mat(obj.mobility_structures{i}.MSD_tot_err_ind_tracks),{'color',cc(i,:)},0) );
			%overlap datapoints
				hold on
				%plot(time,1e-6*cell2mat(obj.mobility_structures{i}.MSD_tot),'o','MarkerFaceColor',cc(i,:),'MarkerEdgeColor','k')
			%fit the msd
				[pF,r]=fitLinearMSD(time,cell2mat(obj.mobility_structures{i}.MSD_tot),time,cell2mat(obj.mobility_structures{i}.MSD_tot));

				%plot(time,pF(1)*time+pF(2),'k')

				fitErr{i}=sqrt(diag(inv(r.R)*inv(r.R)')*(r.normr)^2/r.df);
			
				fitErr{i}=num2str(round(fitErr{i}(1)*10000)/10);

				appD=num2str(round(pF(1)/4.*10000)/10);

				pFs{i}=['/bf{', labels{i} ', D =  ' '(' appD setstr(177) fitErr{i} ')' ' ^{.}10^{-3} um^2/s}'];
				fitDs(i)=[pF(1)/4.];
			end %for i

			pFs=pFs(selection);

			%[legh,objh,outh,outm]=legend([handles.mainLine],pFs,'fontsize', 20 );
			%set(objh,'linewidth',3);
			%set(l,'interpreter','tex')


			%ylim([-min(cell2mat(obj.mobility_structures{2}.MSD_tot)),0.7*max(cell2mat(obj.mobility_structures{2}.MSD_tot))]);
			
			if saveOn
				save_results_path=uigetdir('D:/ANALISYS/2015/1.CLUSTER_OF_RECEPTORS_TRACKING/','choose the folder where to save plots of results');
				saveas(gcf,[save_results_path '/global MSD.pdf'])
				saveas(gcf,[save_results_path '/global MSD.fig'])					
			end

		end%plot_global_MSD

		function [average_D,av_times,average_D_std]=plot_D_ST_overtime(obj,dxThresh_low,toffset,mycolor)
            mycolor2=mycolor;
			if nargin < 3
				dxThresh_low=0;
				mycolor='g'
			end

            %figure('position',[500,-50,550,5950])
            plot_on=1;
			%dxThresh_low=60;
            ST_Ds_temp={};   	
            cc=lines(40);

            myYlim=0.038;
            myYlim_low=0.002;
            tEnd=250;

            %save_results_path=uigetdir('D:/ANALISYS/2015/1.CLUSTER_OF_RECEPTORS_TRACKING/','choose the folder where to save plots of results');
            %save_results_path='D:/ANALISYS/2015/1.CLUSTER_OF_RECEPTORS_TRACKING/trash';
			try 
				obj.mobility_structures{1}.infotracks{1}(1).average_dx;
			catch
				['ERROR:BEFORE RUNNING THIS FUNCTION you should re-run the get-mobility-struct to get average dx for single tracks']
			end

			L=length(obj.mobility_structures);
            Ltot=1;


            average_D=[];
            average_D_std=[];
            av_times=[];
                fraction=[];
                mytimes=[];
                low_all=[];
                low_all_err=[];
                times_all=[];
			for i =1:L

				ST_av_dx=[];

				[Ds,ST_Ds2,ST_Dp,ST_Dm,ST_Dsp,times,absoluteTimes,ST_Ds_std,ST_Ds_pole_std, ST_Ds_middle_std,weigths,weigthsP,weigthsM]=obj.plot_MSD_from_single_MS(i,0); %call the method that extracts Ds and times

				ST_Ds_all=[];
                ST_Ds=[];	
                ST_Ds_std=[];
                ST_Ds_low=[];	
                ST_Ds_std_low=[];

                ST_Ds_p=[];	
                ST_Ds_std_p=[];
                ST_Ds_low_p=[];	
                ST_Ds_std_low_p=[];

            	ST_Ds_m=[];	
                ST_Ds_std_m=[];
                ST_Ds_low_m=[];	
                ST_Ds_std_low_m=[];

            	ST_Ds_s=[];	
                ST_Ds_std_s=[];
                ST_Ds_low_s=[];	
                ST_Ds_std_low_s=[];

				if i==1
					times0=absoluteTimes(1);
				end

				L2=length(Ds);
                    Ds_low_N=0;
                    Ds_N=0;
                    
                    Ds_low_N_p=0;
                    Ds_N_p=0;
                    
                    Ds_low_N_m=0;
                    Ds_N_m=0;

                    Ds_low_N_s=0;
                    Ds_N_s=0;                    


				for j=1:length(obj.mobility_structures{i}.infotracks)

            	nmiddle=0;
            	npole=0;
            	npole_low=0;
            	nmiddle_low=0;
            	nspur=0;
            	nspur_low=0;            	

                    j
                    ST_Ds_temp=[];
                    ST_Ds_temp_low=[];
                    
                    ST_Ds_temp_p=[];
                    ST_Ds_temp_low_p=[];

                    ST_Ds_temp_m=[];
                    ST_Ds_temp_low_m=[];                    

                    ST_Ds_temp_s=[];
                    ST_Ds_temp_low_s=[];                    


					for k=obj.mobility_structures{i}.infotracks{j}
                        k.average_dx=1000;
                        %ST_Ds_temp{j}(1)=0;
                        %snrth=100;
						if k.average_dx > dxThresh_low && k.D_ST > 0.001 && ~isnan(k.D_ST) && k.D_ST<1
                            %[d,s,SNR,a,a]=compute_ST_D_and_s(k,0.106);
                            %if SNR>snrth
                            ST_Ds_temp=[ST_Ds_temp,k.D_ST];
                            %end
							%ST_av_dx=[ST_av_dx,k.average_dx];
							
							if k.position(1)=='p'
								npole=npole+1;
								ST_Ds_temp_p=[ST_Ds_temp_p,k.D_ST];
							end
							if k.position(1)=='m'
								nmiddle=nmiddle+1;
								ST_Ds_temp_m=[ST_Ds_temp_m,k.D_ST];
							end

							if k.position(1)=='s'
								nspur=nspur+1;
								ST_Ds_temp_s=[ST_Ds_temp_s,k.D_ST];
							end

						elseif k.average_dx < dxThresh_low && k.D_ST > 0 && ~isnan(k.D_ST)

							ST_Ds_temp_low=[ST_Ds_temp_low,k.D_ST];

							if k.position(1)=='p'
								npole_low=npole_low+1;
								ST_Ds_temp_low_p=[ST_Ds_temp_low_p,k.D_ST];
							end
							
							if k.position(1)=='m'
								nmiddle_low=nmiddle_low+1;
								ST_Ds_temp_low_m=[ST_Ds_temp_low_m,k.D_ST];
							end		

							if k.position(1)=='s'
								nspur_low=nspur_low+1;
								ST_Ds_temp_low_s=[ST_Ds_temp_low_s,k.D_ST];
							end																	
						end
					end%for k

					%all D_ST for mobile clusters
					ST_Ds_all=[ST_Ds_all,ST_Ds_temp];
					%average D per FOV of mobile clusters
                    ST_Ds=[ST_Ds,nanmean(ST_Ds_temp)];
                    %average D per FOV of immobile clusters
                    ST_Ds_low=[ST_Ds_low,nanmean(ST_Ds_temp_low)];
                    %number of mobile clusters
               		Ds_N=Ds_N+length(ST_Ds_temp);
               		%number of immobile clusters
               		Ds_low_N=Ds_low_N+length(ST_Ds_temp_low);

                    ST_Ds_p=[ST_Ds_p,nanmean(ST_Ds_temp_p)];
                    ST_Ds_low_p=[ST_Ds_low_p,nanmean(ST_Ds_temp_low_p)];
               		Ds_N_p=Ds_N_p+length(ST_Ds_temp_p);
               		Ds_low_N_p=Ds_low_N_p+length(ST_Ds_temp_low_p);
                    
                    ST_Ds_m=[ST_Ds_m,nanmean(ST_Ds_temp_m)];
                    ST_Ds_low_m=[ST_Ds_low_m,nanmean(ST_Ds_temp_low_m)];
               		Ds_N_m=Ds_N_m+length(ST_Ds_temp_m);
               		Ds_low_N_m=Ds_low_N_m+length(ST_Ds_temp_low_m);

                    ST_Ds_s=[ST_Ds_s,nanmean(ST_Ds_temp_s)];
                    ST_Ds_low_s=[ST_Ds_low_s,nanmean(ST_Ds_temp_low_s)];
               		Ds_N_s=Ds_N_s+length(ST_Ds_temp_s);
               		Ds_low_N_s=Ds_low_N_s+length(ST_Ds_temp_low_s);               		
               		
					ST_Ds_std=[ST_Ds_std,nanstd(ST_Ds_temp)/sqrt( length([obj.mobility_structures{i}.infotracks{j}]) )];
					ST_Ds_std_low=[ST_Ds_std_low,nanstd(ST_Ds_temp_low)/sqrt( length([obj.mobility_structures{i}.infotracks{j}]) )];
	
					ST_Ds_std_p=[ST_Ds_std_p,nanstd(ST_Ds_temp_p)/sqrt( npole )];
					ST_Ds_std_low_p=[ST_Ds_std_low_p,nanstd(ST_Ds_temp_low_p)/sqrt( length([obj.mobility_structures{i}.infotracks{j}]) )];

					ST_Ds_std_m=[ST_Ds_std_m,nanstd(ST_Ds_temp_m)/sqrt( nmiddle )];
					ST_Ds_std_low_m=[ST_Ds_std_low_m,nanstd(ST_Ds_temp_low_m)/sqrt( length([obj.mobility_structures{i}.infotracks{j}]) )];

					ST_Ds_std_s=[ST_Ds_std_s,nanstd(ST_Ds_temp_s)/sqrt( nspur )];
					ST_Ds_std_low_s=[ST_Ds_std_low_s,nanstd(ST_Ds_temp_low_s)/sqrt( length([obj.mobility_structures{i}.infotracks{j}]) )];


				%plot((Ltot:Ltot+L2-1),Ds*1e-3,'-o','color',cc(i,:),'LineWidth',3)
				%plot((Ltot:Ltot+L2-1),mean(Ds)*1e-3*ones(length((Ltot:Ltot+L2-1))),'.-','color',cc(i,:),'LineWidth',2)
                end%for j
                t=absoluteTimes-times0 -toffset;
                %if any(t)<0, t=t-t(1),end
                	
%plot results for all clusters                
				if plot_on

					if i==1
						d0=0;%mean(ST_Ds);
                        mycolor=cc(1,:);
                    else
                        mycolor=mycolor;%'m';
                    end

					figure(10)

					t=absoluteTimes-times0;
                    hold on
                    if i==400
                    	errorbar(158,nanmean(ST_Ds),nanstd(ST_Ds_std),'-o','color',mycolor,'LineWidth',2)
                	else                	
                		average_D_std=[average_D_std,nanmean(ST_Ds_std)/sqrt(length(ST_Ds_std))];
						%errorbar(mean(absoluteTimes-times0)-toffset,mean(ST_av_dx),0,'o','color',mycolor,'LineWidth',2)
						%errorbar(mean(absoluteTimes-times0)-toffset,mean(ST_Ds)-d0,mean(ST_Ds_std)/sqrt(length(ST_Ds_std)),'-o','color',mycolor,'LineWidth',2)
                        errorbar(mean(absoluteTimes-times0)-toffset,nanmean(ST_Ds),nanmean(ST_Ds_all)/sqrt(length(ST_Ds_all)),'-s','color',mycolor,'LineWidth',2,'markersize',8)
						%plot(linspace(t(1),t(end),length(ST_av_dx)),ST_av_dx,'.','color',cc(i,:))
						average_D=[average_D,nanmean(ST_Ds)];
						av_times=[av_times,mean(absoluteTimes-times0)-toffset]

					end
					%errorbar(mean(t),sum(ST_Ds.*weigths)/sum(weigths),std(ST_Ds)/sqrt(length(ST_Ds)),'o','color',cc(i,:),'LineWidth',2)  
					%plot(absoluteTimes-times0,sum(ST_Ds.*weigths)/sum(weigths)*ones(length((Ltot:Ltot+L2-1))),'-','color',cc(i,:),'LineWidth',2)					 	
					
					set(gca,'fontsize',20)
					xlabel('time [mins]','fontsize',20)
					%ylabel('D [um^2/s]','fontsize',20)	
					ylabel('D [um^2/s]','fontsize',20)		
					title('all clusters')
					ylim([myYlim_low myYlim])
					xlim([0,tEnd])
					box
					set(gca,'linewidth',2)	
xlim([0 250])
ylim([0.01 0.06])


					figure(1)
					t=absoluteTimes-times0;
                    hold on


                    if i ==1
                    	m=nanmean(ST_Ds);
                	end

                    if i==400
                    	errorbar(158,nanmean(ST_Ds),nanstd(ST_Ds_std),'-o','color',mycolor,'LineWidth',2)
                	else                	
                		tm=mean(absoluteTimes-times0);
						

                		mytext_percentage=num2str(length(ST_Ds_all(ST_Ds_all>m))/double(length(ST_Ds_all)));

						plot(linspace(tm-5,tm+5,length(ST_Ds_all))-toffset,ST_Ds_all,'.','color',mycolor)
						%plot(linspace(tm-5,tm+5,length(ST_Ds_all(ST_Ds_all>m)))-toffset,ST_Ds_all(ST_Ds_all>m),'.','color',cc(i,:))
						%text(tm-10-toffset, 0.1+mean(ST_Ds),[mytext_percentage(3:4) '%'],'fontsize',20);
						%errorbar(mean(absoluteTimes-times0)-toffset,mean(ST_Ds),std(ST_Ds_std),'o','color',mycolor,'LineWidth',2)
						errorbar(mean(absoluteTimes-times0)-toffset,nanmean(ST_Ds),nanstd(ST_Ds)/sqrt(length(ST_Ds)),'o','color','k','LineWidth',2,'markerfacecolor',mycolor,'markersize',8)

						
						

					end

					%errorbar(mean(t),sum(ST_Ds.*weigths)/sum(weigths),std(ST_Ds)/sqrt(length(ST_Ds)),'o','color',cc(i,:),'LineWidth',2)  
					%plot(absoluteTimes-times0,sum(ST_Ds.*weigths)/sum(weigths)*ones(length((Ltot:Ltot+L2-1))),'-','color',cc(i,:),'LineWidth',2)					 	
					%errorbar(absoluteTimes-times0,ST_Ds_low,ST_Ds_std_low,'o','color','k','LineWidth',2)	
					set(gca,'fontsize',20)
					xlabel('time [mins]','fontsize',20)
					%ylabel('D [um^2/s]','fontsize',20)	
					ylabel('D [um^2/s]','fontsize',20)	
					title('all clusters')
					ylim([myYlim_low myYlim])
					xlim([0,tEnd])
					box
					set(gca,'linewidth',2)	
xlim([0 250])
ylim([0.02 0.19])
					%saveas(gcf,[save_results_path '/' 'all_clusters_hi-low' '.pdf'])
					%saveas(gcf,[save_results_path '/' 'all_clusters_hi-low' '.fig'])

					figure(3)	
					hold on			 					
					plot(mean(t)-toffset,Ds_low_N/(Ds_N + Ds_low_N) ,'k^', 'markersize',8,'linewidth',2)					
                    fraction=[fraction, Ds_low_N/(Ds_N + Ds_low_N)]
                    mytimes=[mytimes,mean(absoluteTimes-times0)-toffset]
					set(gca,'fontsize',25)
					xlabel('time [mins]','fontsize',30)
					%ylabel('D [um^2/s]','fontsize',20)	
					ylabel('immobile cluster %','fontsize',30)	
					%title('all clusters')
					ylim([myYlim_low myYlim])
					xlim([0,tEnd])
					box
					set(gca,'linewidth',2)	
	

xlim([0 210])
ylim([0.15 0.4])

figure(20)
                    %errorbar(nanmean(absoluteTimes-times0),nanmean(ST_Ds_low),nanmean(ST_Ds_std_low)/sqrt(length(ST_Ds_std_low)),'o-','color',mycolor,'LineWidth',2)
                    low_all=[low_all,nanmean(ST_Ds_low)];
                    low_all_err=[low_all_err,nanmean(ST_Ds_std_low)/sqrt(length(ST_Ds_std_low))];
                    hold on
                    %{
					saveas(gcf,[save_results_path '/' 'all_clusters_hi-low_fraction' '.pdf'])
					saveas(gcf,[save_results_path '/' 'all_clusters_hi-low_fraction' '.fig'])					
					%plot(absoluteTimes-times0,sum(ST_Ds_low.*weigths)/sum(weigths)*ones(length((Ltot:Ltot+L2-1))),'-','color','k','LineWidth',2)                   

                %}
            end

%plot results for polar clusters                
				if plot_on
					figure(3)

                    hold on
					errorbar(mean(absoluteTimes-times0)-toffset,nanmean(ST_Ds_p),nanmean(ST_Ds_std_p)/sqrt(length(ST_Ds_std_p)),'o','color',mycolor,'LineWidth',2)
					%plot(absoluteTimes-times0,sum(ST_Ds_p.*weigths)/sum(weigths)*ones(length((Ltot:Ltot+L2-1))),'-','color',cc(i,:),'LineWidth',2)  
					%errorbar(absoluteTimes-times0,ST_Ds_low_p,ST_Ds_std_low_p,'o','color','k','LineWidth',2)	
					set(gca,'fontsize',20)
					xlabel('time [mins]','fontsize',20)
					%ylabel('D [um^2/s]','fontsize',20)	
					ylabel('D [um^2/s]','fontsize',20)	
					title('polar clusters')
					ylim([myYlim_low myYlim])
					xlim([0,tEnd])
					box
					set(gca,'linewidth',2)	
xlim([0 250])
ylim([0.006 0.035])
                end
                    %{
					saveas(gcf,[save_results_path '/' 'polar_clusters_hi-low' '.pdf'])
					saveas(gcf,[save_results_path '/' 'polar_clusters_hi-low' '.fig'])	

					figure(4)		
					hold on		
					plot(mean(t),Ds_N_p/(Ds_N_p + Ds_low_N_p) ,'k.', 'markersize',30)
					set(gca,'fontsize',20)
					xlabel('time [mins]','fontsize',20)
					%ylabel('D [um^2/s]','fontsize',20)	
					ylabel('Nfast/N','fontsize',20)	
					title('polar clusters')
					ylim([myYlim_low myYlim])
					xlim([0,tEnd])
					box
					set(gca,'linewidth',2)	
xlim([0 250])
ylim([0.6 0.9])
					saveas(gcf,[save_results_path '/' 'polar_clusters_hi-low_fraction' '.pdf'])
					saveas(gcf,[save_results_path '/' 'polar_clusters_hi-low_fraction' '.fig'])					
					%plot(absoluteTimes-times0,sum(ST_Ds_low.*weigths)/sum(weigths)*ones(length((Ltot:Ltot+L2-1))),'-','color','k','LineWidth',2)                   
                end

%plot results for spurius clusters                
				if plot_on
					figure(5)

                    hold on
					errorbar(absoluteTimes-times0,ST_Ds_s,ST_Ds_std_s,'o','color',cc(i,:),'LineWidth',2)
					plot(absoluteTimes-times0,sum(ST_Ds_s.*weigths)/sum(weigths)*ones(length((Ltot:Ltot+L2-1))),'-','color',cc(i,:),'LineWidth',2)  
					errorbar(absoluteTimes-times0,ST_Ds_low_s,ST_Ds_std_low_s,'o','color','k','LineWidth',2)
					set(gca,'fontsize',20)
					xlabel('time [mins]','fontsize',20)
					%ylabel('D [um^2/s]','fontsize',20)	
					ylabel('D [um^2/s]','fontsize',20)	
					title('spurius clusters')
					ylim([myYlim_low myYlim])
					xlim([0,tEnd])
					box
					set(gca,'linewidth',2)	
xlim([0 250])
ylim([0.006 0.035])
					saveas(gcf,[save_results_path '/' 'spurius_clusters_hi-low' '.pdf'])
					saveas(gcf,[save_results_path '/' 'spurius_clusters_hi-low' '.fig'])	

					figure(6)	
					hold on				
					plot(mean(t),Ds_N_s/(Ds_N_s + Ds_low_N_s) ,'k.', 'markersize',30)
					set(gca,'fontsize',20)
					xlabel('time [mins]','fontsize',20)
					%ylabel('D [um^2/s]','fontsize',20)	
					ylabel('Nfast/N','fontsize',20)	
					title('spurius clusters')
					ylim([myYlim_low myYlim])
					xlim([0,tEnd])
					box
					set(gca,'linewidth',2)		
xlim([0 250])
ylim([0.6 0.9])
					saveas(gcf,[save_results_path '/' 'spurius_clusters_hi-low_fraction' '.pdf'])
					saveas(gcf,[save_results_path '/' 'spurius_clusters_hi-low_fraction' '.fig'])				
					%plot(absoluteTimes-times0,sum(ST_Ds_low.*weigths)/sum(weigths)*ones(length((Ltot:Ltot+L2-1))),'-','color','k','LineWidth',2)                   
                end
%}
%plot results for middle clusters                
				if plot_on
					figure(7)
                    hold on
					errorbar(nanmean(absoluteTimes-times0)-toffset,nanmean(ST_Ds_m),nanmean(ST_Ds_std_m)/sqrt(length(ST_Ds_std_m)),'o','color',mycolor,'LineWidth',2)
					%plot(absoluteTimes-times0,sum(ST_Ds_m.*weigths)/sum(weigths)*ones(length((Ltot:Ltot+L2-1))),'-','color',cc(i,:),'LineWidth',2)  
					%errorbar(absoluteTimes-times0,ST_Ds_low_m,ST_Ds_std_low_m,'o','color','k','LineWidth',2)
					set(gca,'fontsize',20)
					xlabel('time [mins]','fontsize',20)
					%ylabel('D [um^2/s]','fontsize',20)	
					ylabel('D [um^2/s]','fontsize',20)		
					title('middle clusters')
					ylim([myYlim_low myYlim])
					xlim([0,tEnd])
					box
					set(gca,'linewidth',2)
xlim([0 250])
ylim([0.006 0.035])
                end
%{
					saveas(gcf,[save_results_path '/' 'middle_clusters_hi-low' '.pdf'])
					saveas(gcf,[save_results_path '/' 'middle_clusters_hi-low' '.fig'])	

					figure(8)
					hold on
					plot(mean(t),Ds_N_m/(Ds_N_m + Ds_low_N_m) ,'k.', 'markersize',30)
					set(gca,'fontsize',20)
					xlabel('time [mins]','fontsize',20)
					%ylabel('D [um^2/s]','fontsize',20)	
					ylabel('Nfast/N','fontsize',20)		
					title('middle clusters')
					ylim([myYlim_low myYlim])
					xlim([0,tEnd])
					box
					set(gca,'linewidth',2)
xlim([0 250])
ylim([0.6 0.9])

					saveas(gcf,[save_results_path '/' 'middle_clusters_hi-low' '.pdf'])
					saveas(gcf,[save_results_path '/' 'middle_clusters_hi-low' '.fig'])	

					%plot(absoluteTimes-times0,sum(ST_Ds_low.*weigths)/sum(weigths)*ones(length((Ltot:Ltot+L2-1))),'-','color','k','LineWidth',2)                   
                end
             %}  
				T=absoluteTimes-times0;
				Ltot=Ltot+L2;
                
                times_all=[times_all,nanmean(absoluteTimes-times0)-toffset]
                
                            
			end%for i				

            
                    figure(5)
					hold on			 					
					plot(mytimes,fraction ,'ko-', 'markersize',8,'linewidth',2)	
                    figure(20)
                    errorbar(times_all,low_all,low_all_err,'o-','color',mycolor2,'LineWidth',2)
                    

            mytimes
            fraction
		end%function plot_D_ST_overtime






		function plot_tracklength_dX_correlation2D(obj,ref_MS)

			save_results_path=uigetdir('D:/ANALISYS/2015/1.CLUSTER_OF_RECEPTORS_TRACKING/','choose the folder where to save plots of results');

			cc=lines(50);
			if nargin < 2
				dXthreshold=250;
				dXthreshold_low=0;
				dx_step=15;
				ref_MS=1;	%mobility structure index (subdataset) to use as a reference for the differences
				TLthreshold=80;
				TL_step=15;
				plot_On=1;
				pixelsize=106;			
			end

			if nargin < 3
				dXthreshold=250;
				dXthreshold_low=0;
				dx_step=15;
				TLthreshold=80;
				TL_step=15;
				plot_On=1;
				pixelsize=106;			
			end			

			L=length(obj.mobility_structures);

			dX={}; %all singletrack step displacement
			dX_p={}; %singletrack step displacement for polar tracks
			dX_m={}; %singletrack step displacement for middle tracks
			TL={}; %all corresponding tracklenghts
			TL_p={}; %all corresponding polar tracklenghts
			TL_m={}; %all corresponding middle tracklenghts

			P_TL_dX={};
			TL_dX_axes={};
			P_TL_dX_p={};
			TL_dX_axes_p={};
			P_TL_dX_m={};
			TL_dX_axes_m={};

			%fills the array with tracklengths and corresponding dX
			myAxes=perfect_subplot(3,L-1);

			for i=1:L

				dX{i}(1)=0;
				dX_p{i}(1)=0;
				dX_m{i}(1)=0;
				TL{i}(1)=0;
				TL_p{i}(1)=0;
				TL_m{i}(1)=0;
			
				for j = obj.mobility_structures{i}.infotracks;			

					for k=[j{:}]

						dx=k.x-circshift(k.x,[0,-1]);
						dy=k.y-circshift(k.y,[0,-1]);

						dx=dx(1:end-1);
						dy=dy(1:end-1);

						dX_temp=(sqrt(dx.^2+dy.^2))*pixelsize;

						if mean(dX_temp) > 0 && mean(dX_temp) < dXthreshold && length(k.x)<100 && mean(dX_temp) > dXthreshold_low

							dX{i}=[dX{i}, mean(dX_temp)];
							TL{i}=[TL{i}, length(k.x)];
							%TL{i}=[TL{i}, length(k.x)*ones(1,length(k.x)-1)];
							if k.position(1)=='p'
								dX_p{i}=[dX_p{i}, mean(dX_temp)];
								TL_p{i}=[TL_p{i}, length(k.x)];	
							end							
							if k.position(1)=='m'
								dX_m{i}=[dX_m{i}, mean(dX_temp)];
								TL_m{i}=[TL_m{i}, length(k.x)];
							end					

						end %if

					end %for k

				end %for j

				[P_TL_dX{i},TL_dX_axes{i}]=hist3([TL{i}; dX{i}]','edges',{1:TLthreshold/TL_step:TLthreshold,1:dXthreshold/dx_step:dXthreshold});
				[P_TL_dX_p{i},TL_dX_axes_p{i}]=hist3([TL_p{i}; dX_p{i}]','edges',{1:TLthreshold/TL_step:TLthreshold,1:dXthreshold/dx_step:dXthreshold});
				[P_TL_dX_m{i},TL_dX_axes_m{i}]=hist3([TL_m{i}; dX_m{i}]','edges',{1:TLthreshold/TL_step:TLthreshold,1:dXthreshold/dx_step:dXthreshold});
				
				Norm=sum(sum(P_TL_dX{i}));
				P_TL_dX{i}=P_TL_dX{i}'/Norm; %normalize and transpose
				Norm_p=sum(sum(P_TL_dX_p{i}));
				P_TL_dX_p{i}=P_TL_dX_p{i}'/Norm_p;
				Norm_m=sum(sum(P_TL_dX_m{i}));
				P_TL_dX_m{i}=P_TL_dX_m{i}'/Norm_m;


	
			end%for	i

			kk=1;

			for ii=1:L	



				if plot_On  && ii ~= ref_MS


					kk=kk+1;
	
					%all
					imagesc(TL_dX_axes{ii}{1},TL_dX_axes{ii}{2}, (P_TL_dX{ii}-P_TL_dX{ref_MS}),'parent',myAxes(kk-1))
					set(myAxes(kk-1), 'YDir', 'normal');
					text(0.8*TLthreshold,0.85*dXthreshold,num2str(sum(sum(abs(P_TL_dX{ii}-P_TL_dX{ref_MS})))),'parent',myAxes(kk-1))
					if ii~=L
						set(myAxes(kk-1), 'XTicklabel',[]);
					end
					set(myAxes(kk-1), 'YTicklabel',[]);			
					ylim(myAxes(kk-1),[0,dXthreshold])
					if kk==2
						title(myAxes(kk-1),'all')
					end

					
					caxis(myAxes(kk-1),[-0.002 0.002])


					%polar
					imagesc(TL_dX_axes_p{ii}{1},TL_dX_axes_p{ii}{2}, (P_TL_dX_p{ii}-P_TL_dX_p{ref_MS}),'parent',myAxes(kk-1 + (L-1)))
					set(myAxes(kk-1 + (L-1)), 'YDir', 'normal');
					text(0.8*TLthreshold,0.85*dXthreshold,num2str(sum(sum(abs(P_TL_dX_p{ii}-P_TL_dX_p{ref_MS})))),'parent',myAxes(kk-1 + (L-1)))
					if ii~=L
					    set(myAxes(kk-1 + (L-1)), 'XTicklabel',[]);				    			    
			    	else
				    	xlabel(myAxes(kk-1 + (L-1)),'Tracklength [frames]','fontsize',10)	
			    	end
				    set(myAxes(kk-1 + (L-1)), 'YTicklabel',[]);
					ylim(myAxes(kk-1 + (L-1)),[0,dXthreshold])
					if kk==2
						title(myAxes(kk-1 + (L-1)),'polar')
					end					
					
					caxis(myAxes(kk-1 + (L-1)),[-0.002 0.002])

					%middle
					imagesc(TL_dX_axes_m{ii}{1},TL_dX_axes_m{ii}{2}, (P_TL_dX_m{ii}-P_TL_dX_m{ref_MS}),'parent', myAxes(kk-1 + (L-1)*2))
					set(myAxes(kk-1 + (L-1)*2), 'YDir', 'normal');
					text(0.8*TLthreshold,0.85*dXthreshold,num2str(sum(sum(abs(P_TL_dX_m{ii}-P_TL_dX_m{ref_MS})))),'parent',myAxes(kk-1 + (L-1)*2))
					if ii~=L
						set(myAxes(kk-1 + (L-1)*2), 'XTicklabel',[]);
					end
					ylim(myAxes(kk-1 + (L-1)*2),[0,dXthreshold])
					if kk==2
						title(myAxes(kk-1 + (L-1)*2),'middle')
					end					
					ylabel(myAxes(kk-1 + (L-1)*2),['P(dx,L)_',num2str(ii),'-','P(dx,L)_',num2str(ref_MS)],'fontsize',10)							
					caxis(myAxes(kk-1 + (L-1)*2),[-0.002 0.002])

					%plot marginal distributions of differences of P(dx, L) for all clusters 
					figure(4)
					ylim(gca,[-0.1,0.1])
					box
					plot(TL_dX_axes{ii}{2}, sum(P_TL_dX{ii}-P_TL_dX{ref_MS},2),'color',cc(ii,:),'linewidth',3)
					hold on
					plot(TL_dX_axes{ii}{2}, zeros(1,length(TL_dX_axes{ii}{2})),'-k','linewidth',2)
					set(gca,'fontsize',25,'linewidth',3)
					ylabel(gca,['P(dx)_',num2str(ii),'-','P(dx)_',num2str(ref_MS)],'fontsize',25)
					xlabel(gca,'dx [nm]','fontsize',25)
					ylim(gca,[-0.1,0.1])
					saveas(gcf,[save_results_path '/' ['P(dx)_',num2str(ii),'-','P(dx)_',num2str(ref_MS)] '.pdf'])
					saveas(gcf,[save_results_path '/' ['P(dx)_',num2str(ii),'-','P(dx)_',num2str(ref_MS)] '.fig'])	

				end	%if					
			end %for ii

%{

					%all
					figure('position',[500,100,350,250])	

					imagesc(TL_dX_axes{i}{1},TL_dX_axes{i}{2}, (P_TL_dX{i}-P_TL_dX{1}))
					set(gca, 'YDir', 'normal');
				
					ylim([0,dXthreshold])

					colormap(jet)
					colorbar

					caxis([-0.003 0.003])

					set(gca,'fontsize',15)					
					ylabel('dX [nm]','fontsize',15)
					xlabel('Tracklength','fontsize',15)

					title(['dataset ' num2str(i) '- dataset 1 '])	

					%polar
					figure('position',[500,100,350,250])	

					imagesc(TL_dX_axes_p{i}{1},TL_dX_axes_p{i}{2}, (P_TL_dX_p{i}-P_TL_dX_p{1}))
					set(gca, 'YDir', 'normal');
				
					ylim([0,dXthreshold])

					colormap(jet)
					colorbar

					caxis([-0.003 0.003])

					set(gca,'fontsize',15)					
					ylabel('dX [nm]','fontsize',15)
					xlabel('Tracklength','fontsize',15)

					title(['dataset ' num2str(i) '- dataset 1 polar'])	

					%middle
					figure('position',[500,100,350,250])	

					imagesc(TL_dX_axes_m{i}{1},TL_dX_axes_m{i}{2}, (P_TL_dX_m{i}-P_TL_dX_m{1}))
					set(gca, 'YDir', 'normal');
				
					ylim([0,dXthreshold])

					colormap(jet)
					colorbar

					caxis([-0.003 0.003])

					set(gca,'fontsize',15)					
					ylabel('dX [nm]','fontsize',15)
					xlabel('Tracklength','fontsize',15)

					title(['dataset ' num2str(i) '- dataset 1 middle'])	
	%}
				

%{
				if plot_On  

					figure('position',[500,100,350,250])	


					imagesc(TL_dX_axes{i}{1},TL_dX_axes{i}{2}, (P_TL_dX{i}))%-P_TL_dX{1}))
					axis xy
					
					ylim([0,dXthreshold])

					colormap(jet)
					colorbar

					caxis([-0.0 0.05])

					set(gca,'fontsize',15)					
					ylabel('dX [nm]','fontsize',15)
					xlabel('Tracklength','fontsize',15)		

					title(['dataset ' num2str(i)])	
				end	%if
%}

			

		end%plot_tracklength_dX_correlation2D


		function []=plot_AMP_dX_correlation2D(obj)

			nargin

			if nargin < 2
				dXthreshold=250;
				dx_step=20;
				ref_MS=1;	%mobility structure index (subdataset) to use as a reference for the differences
				AMPthreshold=5000;
				AMP_step=10;
				plot_On=1;
				pixelsize=106;			
			end

			if nargin < 3
				dXthreshold=250;
				dx_step=10;
				AMPthreshold=5000;
				AMP_step=20;
				plot_On=1;
				pixelsize=106;			
			end			

			L=length(obj.mobility_structures);

			dX={}; %all singletrack step displacement
			dX_p={}; %singletrack step displacement for polar tracks
			dX_m={}; %singletrack step displacement for middle tracks
			AMP={}; %all corresponding tracklenghts
			AMP_p={}; %all corresponding polar tracklenghts
			AMP_m={}; %all corresponding middle tracklenghts

			P_AMP_dX={};
			AMP_dX_axes={};
			P_AMP_dX_p={};
			AMP_dX_axes_p={};
			P_AMP_dX_m={};
			AMP_dX_axes_m={};

			%fills the array with tracklengths and corresponding dX
			myAxes=perfect_subplot(3,L-1);

			for i=1:L

				dX{i}(1)=0;
				dX_p{i}(1)=0;
				dX_m{i}(1)=0;
				AMP{i}(1)=0;
				AMP_p{i}(1)=0;
				AMP_m{i}(1)=0;
			
				for j = obj.mobility_structures{i}.infotracks;			

					for k=[j{:}]

						dx=k.x-circshift(k.x,[0,-1]);
						dy=k.y-circshift(k.y,[0,-1]);

						dx=dx(1:end-1);
						dy=dy(1:end-1);

						dX_temp=(sqrt(dx.^2+dy.^2))*pixelsize;

						if mean(dX_temp) > 0 && mean(dX_temp) < dXthreshold && length(k.x)<1000

							dX{i}=[dX{i}, mean(dX_temp)];
							AMP{i}=[AMP{i}, k.AMP_firstframe];
							%AMP{i}=[AMP{i}, length(k.x)*ones(1,length(k.x)-1)];
							if k.position(1)=='p'
								dX_p{i}=[dX_p{i}, mean(dX_temp)];
								AMP_p{i}=[AMP_p{i}, k.AMP_firstframe];	
							end							
							if k.position(1)=='m'
								dX_m{i}=[dX_m{i}, mean(dX_temp)];
								AMP_m{i}=[AMP_m{i}, k.AMP_firstframe];
							end					

						end %if

					end %for k

				end %for j

				[P_AMP_dX{i},AMP_dX_axes{i}]=hist3([AMP{i}; dX{i}]','edges',{1:AMPthreshold/AMP_step:AMPthreshold,1:dXthreshold/dx_step:dXthreshold});
				[P_AMP_dX_p{i},AMP_dX_axes_p{i}]=hist3([AMP_p{i}; dX_p{i}]','edges',{1:AMPthreshold/AMP_step:AMPthreshold,1:dXthreshold/dx_step:dXthreshold});
				[P_AMP_dX_m{i},AMP_dX_axes_m{i}]=hist3([AMP_m{i}; dX_m{i}]','edges',{1:AMPthreshold/AMP_step:AMPthreshold,1:dXthreshold/dx_step:dXthreshold});
				
				Norm=sum(sum(P_AMP_dX{i}));
				P_AMP_dX{i}=P_AMP_dX{i}'/Norm; %normalize and transpose
				Norm_p=sum(sum(P_AMP_dX_p{i}));
				P_AMP_dX_p{i}=P_AMP_dX_p{i}'/Norm_p;
				Norm_m=sum(sum(P_AMP_dX_m{i}));
				P_AMP_dX_m{i}=P_AMP_dX_m{i}'/Norm_m;
		
			end%for	i

			kk=1;

			for ii=1:L	



				if plot_On  && ii ~= ref_MS

					kk=kk+1;

					%all
					imagesc(AMP_dX_axes{ii}{1},AMP_dX_axes{ii}{2}, (P_AMP_dX{ii}-P_AMP_dX{ref_MS}),'parent',myAxes(kk-1))
					set(myAxes(kk-1), 'YDir', 'normal');
					text(0.8*AMPthreshold,0.85*dXthreshold,num2str(sum(sum(abs(P_AMP_dX{ii}-P_AMP_dX{ref_MS})))),'parent',myAxes(kk-1))
					if ii~=L
						set(myAxes(kk-1), 'XTicklabel',[]);
					end
					set(myAxes(kk-1), 'YTicklabel',[]);			
					ylim([0,dXthreshold])
					if kk==2
						title(myAxes(kk-1),'all')
					end

					colormap(jet)
					caxis(myAxes(kk-1),[-0.003 0.003])

					%polar
					imagesc(AMP_dX_axes_p{ii}{1},AMP_dX_axes_p{ii}{2}, (P_AMP_dX_p{ii}-P_AMP_dX_p{ref_MS}),'parent',myAxes(kk-1 + (L-1)))
					set(myAxes(kk-1 + (L-1)), 'YDir', 'normal');
					text(0.8*AMPthreshold,0.85*dXthreshold,num2str(sum(sum(abs(P_AMP_dX_p{ii}-P_AMP_dX_p{ref_MS})))),'parent',myAxes(kk-1 + (L-1)))
					if ii~=L
					    set(myAxes(kk-1 + (L-1)), 'XTicklabel',[]);				    			    
			    	else
				    	xlabel(myAxes(kk-1 + (L-1)),'Tracklength [frames]','fontsize',10)	
			    	end
				    set(myAxes(kk-1 + (L-1)), 'YTicklabel',[]);
					ylim([0,dXthreshold])
					if kk==2
						title(myAxes(kk-1 + (L-1)),'polar')
					end					
					colormap(jet)
					caxis(myAxes(kk-1 + (L-1)),[-0.002 0.002])

					%middle
					imagesc(AMP_dX_axes_m{ii}{1},AMP_dX_axes_m{ii}{2}, (P_AMP_dX_m{ii}-P_AMP_dX_m{ref_MS}),'parent', myAxes(kk-1 + (L-1)*2))
					set(myAxes(kk-1 + (L-1)*2), 'YDir', 'normal');
					text(0.8*AMPthreshold,0.85*dXthreshold,num2str(sum(sum(abs(P_AMP_dX_m{ii}-P_AMP_dX_m{ref_MS})))),'parent',myAxes(kk-1 + (L-1)*2))
					if ii~=L
						set(myAxes(kk-1 + (L-1)*2), 'XTicklabel',[]);
					end
					ylim([0,dXthreshold])
					if kk==2
						title(myAxes(kk-1 + (L-1)*2),'middle')
					end					
					ylabel(myAxes(kk-1 + (L-1)*2),['P(dx)_',num2str(ii),'-','P(dx)_',num2str(ref_MS)],'fontsize',10)
					colormap(jet)					
					caxis(myAxes(kk-1 + (L-1)*2),[-0.003 0.003])

				end	%if					
			end %for ii		

		end%plot_amp_correlation_2D

		function [av_D_FOVs]=hist_average_D_from_FOVs(obj)

			av_D_FOVs=[];

			L= length(obj.mobility_structures);

			for i = 1:L

				N_FOVs=length(obj.mobility_structures{i}.infotracks)

				[Ds,ST_Ds,ST_Dp,ST_Dm,ST_Dsp,times,absoluteTimes,ST_Ds_std,ST_Ds_pole_std, ST_Ds_middle_std,weigths,weigthsP,weigthsM]=obj.plot_MSD_from_single_MS(i,0);

				av_D_FOVs=[av_D_FOVs, ST_Ds];

			end%for i

		end %hist_average_D_from_FOVs

		function plot_D_cdf(obj, Nbins)
threshold_track_length=5;
pixelsize=0.106
			if nargin < 2
				Nbins=100;
			end

			L=length(obj.mobility_structures);

			try Dall=obj.mobility_structures{1}.ALLsingle_track_D;

			catch
				[a,a,a,a,a,a,a,a,a,a,D, a]=computeDandS(obj,threshold_track_length,pixelsize);
				for k=1:L
					obj.mobility_structures{k}.ALLsingle_track_D=D{k};
				end
			end

			Dall={};
			Dm={};
			Dp={};

			%P=[0    0   435   378];

			h1=figure(1);
			ax1=gca(h1);
			%set(h1,'position',P)
			box

			h2=figure(2)
			ax2=gca(h2);
			%set(h2,'position',P)
			box

			h3=figure(3)
			ax3=gca(h3);
			%set(h3,'position',P)
			box			

			for i = 1:L

				Dall=obj.mobility_structures{i}.ALLsingle_track_D;
				%Dp=obj.mobility_structures{i}.POLEsingle_track_D;
				%Dm=obj.mobility_structures{i}.MIDDLEsingle_track_D;

				[Nall,binsAll]=hist(Dall,Nbins);				
				%[Np,binsp]=hist(Dp,Nbins);
				%[Nm,binsm]=hist(Dm,Nbins);

				Nall=Nall/sum(Nall)/(binsAll(3)-binsAll(2));
				%Np=Np/sum(Np)/(binsp(3)-binsp(2));
				%Nm=Nm/sum(Nm)/(binsm(3)-binsm(2));

				for j=1:Nbins

					cdfAll(j)=sum(Nall(1:j)*(binsAll(2)-binsAll(1)));
					%cdfp(j)=sum(Np(1:j)*(binsp(2)-binsp(1)));
					%cdfm(j)=sum(Nm(1:j)*(binsm(2)-binsm(1)));

				end %for j

				 %plotcdf of all D
				figure(1)
				hold on
				plot(ax1,binsAll,1-cdfAll,'linewidth',3)
				set(ax1,'FontWeight','normal','linewidth',2,'fontsize',25)
				xlim([0 0.15])
				ylim([0 0.9])
				ylabel('1 - cdf','fontsize',20)
				xlabel('D [um^2/s]','fontsize',20)
				title('all','fontsize',20)			
				
%{

				 %plotcdf of polar D
				figure(2)
				hold on
				plot(ax2,binsp,1-cdfp,'linewidth',3)
				set(ax2,'FontWeight','normal','linewidth',2,'fontsize',25)
				xlim([0 0.15])
				ylim([0 0.9])
				ylabel('1 - cdf','fontsize',20)
				xlabel('D [um^2/s]','fontsize',20)		
				title('polar','fontsize',20)			


				 %plotcdf of middle D
				figure(3)
				hold on
				plot(ax3,binsm,1-cdfm,'linewidth',3)				
				set(ax3,'FontWeight','normal','linewidth',2,'fontsize',25)
				xlim([0 0.15])
				ylim([0 0.9])
				ylabel('1 - cdf','fontsize',20)
				xlabel('D [um^2/s]','fontsize',20)	
				title('lateral','fontsize',20)			

				length(Dall)
				length(Dp)
				length(Dm)
%}
			end %for i

		end %plot_D_cdf

	end %methods

end %class

%{
		function plot_tracklength_dX_correlation2D(obj)

			if nargin < 2
				dXthreshold=250;
				plot_On=1;
				pixelsize=106;
			end

			L=length(obj.mobility_structures);

			dX={}; %all singletrack average displacement
			TL={}; %all corresponding tracklenghts

			P_TL_dX={};
			TL_dX_axes={};
			namess={}
			%fills the array with tracklengths and corresponding dX
			for i=1:L

				dX{i}(1)=0;
				TL{i}(1)=0;
			
				mydir=uigetdir(obj.saveroot);
				filess=dir(mydir);
				filess(1:2)=[];

				for c =1:length(filess)
					namess{c}=[mydir,filesep, filess(i).name]	;				
				end 



				for j = 1:length(obj.mobility_structures{i}.infotracks);	


					alltracks=importdata(namess{j});

					for m=1:length(alltracks)

						k=alltracks(m);

						k.x=k.tracksCoordAmpCG(1:8:end);
						k.y=k.tracksCoordAmpCG(2:8:end);

						dx=k.x-circshift(k.x,[0,-1]);
						dy=k.y-circshift(k.y,[0,-1]);

						dx=dx(1:end-1);
						dy=dy(1:end-1);

						dX_temp=mean(sqrt(dx.^2+dy.^2))*pixelsize;

						if dX_temp > 0 && dX_temp < dXthreshold && length(k.x)<100

							dX{i}=[dX{i}, dX_temp];
							TL{i}=[TL{i}, length(k.x)];

						end %if

					end %for k

				end %for j

				[P_TL_dX{i},TL_dX_axes{i}]=hist3([TL{i}; dX{i}]',[25,25])
			

				if plot_On  %&& i >1

					figure('position',[500,100,350,250])	

					hist(TL{i},100)

					%imagesc(TL_dX_axes{i}{2},TL_dX_axes{i}{1}, (P_TL_dX_temp{i})')

					figure('position',[500,100,350,250])			

					Norm1=sum(sum(P_TL_dX{i}));
					%Norm2=sum(sum(P_TL_dX{i-1}));

					P_TL_dX{i}=P_TL_dX{i}/Norm1;
					%P_TL_dX{i-1}=P_TL_dX{i-1}/Norm2;

					imagesc(TL_dX_axes{i}{1},TL_dX_axes{i}{2}, (P_TL_dX{i}))%-P_TL_dX{1}))
					axis xy
					
					ylim([0,dXthreshold])

					colormap(jet)
					colorbar

					%caxis([-0.01 0.01])

					axis xy

					set(gca,'fontsize',15)					
					ylabel('dX [nm]','fontsize',15)
					xlabel('Tracklength','fontsize',15)		
				end	%if

			end%for

		end%plot_tracklength_dX_correlation2D

%}