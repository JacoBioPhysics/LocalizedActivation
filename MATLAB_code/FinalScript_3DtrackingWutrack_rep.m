%this script allow the user to select target folder to analyze with the
%utrack software and then fit the result of the utrack analysis with
%MLEw_elli_iter for 3D trakcing and max likelihood fitting. 


%*************************************************************************
%--------------------DATA PREPARATION FOR FIT-----------------------------
%*************************************************************************

Nfolder=1;

%Initialize paths and names for:

%   frames_separators:  first and last frames to analyize in the sub_dataset
%   datasetname:        name of the complete dataset to analyze
%   saveroot:           root path where to save all the different datasets tracking results
%   prepareddatafolder: path where all the all_target_folders lie

%   folderlist:         getting and copying images in .tif format for u-track (it does not like .Tiff ..) from data folders
%   all_target_folders: '.tif' copied images directories paths for analysis with u-track
%   save_track_path:    paths where to save results of the tracking
%   sub_datasetname:    name of the sub_dataset like 'solari_20141002_190012'

frames_separators={'2','2000'};
datasetname='fakedata_symm';
saveroot=strcat('C:\Users\solari\Desktop\utrack_output\',datasetname,'\');
prepareddatafolder=['D:\ANALISYS\TRIALS\' datasetname];


folderlist={};
all_target_folders={};
save_track_path={};
sub_datasetname=[];

mkdir(prepareddatafolder);


%threshold for track length, only those tracks longer than threshold are kept
threshold=10;

%user selection of the data folders (have to be Frames00000x folders)


dataroot=uigetdir('D:','select the folder containing the datasets to analyze');
dd=dir(dataroot);
dd(1:2)=[];
dd=dd(cat(1,dd.isdir));

for i = 1:Nfolder
    folderlist{i}=[dataroot '\' dd(i).name '\Frames000001' ]
end
        



%prepare the folders and the '.tif' copied images 
for j=1:Nfolder
    
    ind=find(folderlist{j}=='\');
    sub_datasetname=folderlist{j}(ind(end-1):ind(end));
    destination=strcat(prepareddatafolder,sub_datasetname);
    mkdir(destination)   %make the directory to store converted data readable by utrack in sub_datasetname folders  

    all_target_folders{j}=strcat(destination,frames_separators(1),'-',frames_separators(2)); %makes the directory inside the sub_datasetname folder indicating the frame number of the images inside the folder
    mkdir(all_target_folders{j}{:});  %creates the directory inside the directory dataset
    save_track_path{j}=strcat(saveroot(1:end-1),sub_datasetname,frames_separators(1),'-',frames_separators(2),'\');
    folderlist{j}
    all_target_folders{j}

    ifdir=dir([folderlist{j} '\*.Tiff'])
    if ~isempty(ifdir)
        TiffORtif='T';
    else
        TiffORtif='t'
    end
        
    PrepareImagesForUtrack_rep ... %call the actual function that does the job
            (folderlist{j},all_target_folders{j}, frames_separators{1}, frames_separators{2},TiffORtif);
    
        
    
end

%*************************************************************************
%--------------------UTRACK TRACKING OF THE PREPARED DATA-----------------
%*************************************************************************

%***reshape the call array to the same level

all_target_folders=cat(1,all_target_folders{:});
all_target_folders=cat(1,all_target_folders(:));
Ntarget=length(all_target_folders);
save_track_path=cat(1,save_track_path{:});
save_track_path=cat(1,save_track_path(:));

for h =1:Ntarget
    
    h
        
    cd(all_target_folders{h})

    myscriptDetectGeneral_rep
    myscriptTrackGeneral_rep
    
    tracksFinal(1).datafolder=folderlist{h};
    tracksFinal(1).prepareddatafolder=all_target_folders{h};
    tracksFinal(1).save_track_path=save_track_path{h};

    ind2=find(all_target_folders{h}=='\');
    timedatedata=all_target_folders{h}(ind2(end-1)+7:ind2(end)-1);
 
    strcat(save_track_path,'tracksFinal',timedatedata,'.mat')
    mkdir(save_track_path{h})
    save(strcat(save_track_path{h},'tracksFinal',timedatedata(2:end),'.mat'),'tracksFinal')

    %{
    try
        fit_tracksFinal_fromUtrack(tracksFinal, all_target_folders{h}, threshold, save_track_path{h})
    catch
        mkdir('porcoddio ha crashato fit elli')
    end
    %}
end

    
    
    
    