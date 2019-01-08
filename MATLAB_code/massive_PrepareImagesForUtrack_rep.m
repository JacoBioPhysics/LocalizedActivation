function [] = massive_PrepareImagesForUtrack_rep(N1,N2,TiffORtif)

    %destination=uigetdir('H:\Jacopo','select the directory where to duplicate your dataset')
    massive_source=uigetdir('H:\Jacopo','select the directory containing datasets to duplicate your dataset');
    
    cd(massive_source);
    D=dir();
    
    D(1:2)=[];
    D(~[D.isdir])=[];
    
    mkdir('tif_images');
    cd('tif_images');
    
    for i = 1:length(D)        
        mkdir(D(i).name)
        cd(D(i).name)
        PrepareImagesForUtrack_rep([massive_source,'\', D(i).name '\Frames000001\' ],'.',N1,N2,'T');
        cd('..')
    end

end