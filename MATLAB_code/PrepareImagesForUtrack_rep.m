%select a folder from where to get images and a folder where to save them
%in order to make them recognizable by u-track2.1.1

function [] = PrepareImagesForUtrack_rep(source,destination,N1,N2,TiffORtif)


%source=uigetdir('D:','select the folder containing the images to be converted');
%destination=uigetdir('D:','Select the folder where to save the images');

    'Preparing Images..'
    if TiffORtif=='T'
        source
        d=dir(strcat(source,'\*.Tiff'))
    else
        d=dir(strcat(source,'\*.tif'))
    end
    

    N1=str2num(N1);
    N2=str2num(N2);

    if N2>length(d)
        N2=length(d)
    end

    for i=N1:N2
        
        
        img=imread(strcat(source,'\', d(i).name));        
        imwrite(uint16(img),char(strcat(destination,'\Image',sprintf('%4.4d',i),'.tif')));
        clear img
        
    end
    ['Images Prepared in' destination '!']

end