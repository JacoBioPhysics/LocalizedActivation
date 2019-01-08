function [angle,box,mask2, position, bin_edges,maskbox, msd2, msd2_err,pa,p, D_err,velocity_angles,D_ST,average_dx]=singletrack_spat_analysis_rep(singletrack, mask,i,timedelay,pixelsize, NpointsFit,NpointsFitAlpha)% tr_name,savepath)
    
% singletrack_spat_analysis: assign a cellular position to a single track 

%   INPUTS:     

%       singletrack: single track from a track set output from Utrack tracking routine
%       mask       : bright field mask from python function BrightField_Segmentation

% This function takes a track and the labeled mask and find the feature (bacterium) on which the track fall onto.
% then it rotates the cut mask corresponding to the identified bacterium in order to make it horizontal.
% after that it histograms the x position over the cell and if the frequency of localization in a particular cellular 
% compartment is > fraction, it assigns that position to the track.
%
%it also computes the MSD for single tracks and fit it to a straight line considering only the first n poinbts


    mask=bwlabel(mask);

    
    %identify the bacterium on which the tracks falls on by finding the
    %intersecition between the mask and the points of the track
    
    x=singletrack.tracksCoordAmpCG(1:8:end);
    y=singletrack.tracksCoordAmpCG(2:8:end);
    
    singletrack.x=x;
    singletrack.y=y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       CELLULAR POSITION ASSIGNMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fraction=.85; %fraction of point to be found in the middle or polar refgion for position assignment

    cols=fix(y);
    lines=fix(x);
    s=size(mask);
    inds=lines*s(1)+cols;
    inds=inds(~isnan(inds));
       
    value=mask(inds);
    value=nonzeros(value);
    value=fix(mean(value));

    if isnan(value)
        mask=bwmorph(mask,'dilate');
        value=mask(inds);
        value=nonzeros(value);
        value=fix(mean(value));

    end

    if isnan(value)
        angle=0;
        box=0;
        mask2=0;
        position='cazzo';
        bin_edges=0;
        maskbox=0;
        msd2=0;
        msd2_err=0;
        p=[NaN,NaN];
        pa=[NaN,NaN];
        D_err=0;
        velocity_angles=[];
        'porciodio'
        D_ST=0;
        average_dx=0;
        return
    end

    mask=mask==value;
    %mask=bwmorph(mask,'dilate'); 

    mask=int16(mask);

    mask(inds)=value+1;
 

    %isolate the bacterium from the mask
    rp=regionprops(mask);    
    box=rp.BoundingBox;
    box=fix(box);
    maskbox=(mask(box(2)+1:box(2)+box(4),box(1)+1:box(1)+box(3)));

    bwmaskbox=maskbox~=0;
    rp2=regionprops(bwmaskbox);



    %rotate the bacterium to an horizontal position
    %
    %find the couple of points that maximize the distance from the centroid
    boundary=bwboundaries(bwmaskbox);
    boundary=boundary{1};
    c=rp2.Centroid;
    distance=sqrt((boundary(:,1) - c(2)).^2 + (boundary(:,2) - c(1)).^2);
    
    
    myind=(find(distance==max(distance)));
    myind=myind(1);
    
    distance2=sqrt((boundary(:,1) - boundary(myind,1)).^2 + (boundary(:,2) - boundary(myind,2)).^2);
    myind2=find(distance2==max(distance2));
    myind2=myind2(1);




    ym=boundary(myind,2);
    yM=boundary(myind2,2);
    xm=boundary(myind,1);
    xM=boundary(myind2,1);

    angle=atan((xM-xm)/(yM-ym))/2./pi*360;

    try
        mask2=imrotate(maskbox , angle);
    catch
        'diocan'
        angle=0;
        box=0;
        mask2=0;
        position='cazzo';
        bin_edges=0;
        maskbox=0;
        msd2=0;
        msd2_err=0;
        p=[NaN,NaN];
        pa=[NaN,NaN]
        D_err=0;
        velocity_angles=[];
        D_ST=0;
        average_dx=0;        
        return
    end        

    rp3=regionprops(mask2);
    box=rp3.BoundingBox;
    box=floor(box);


    mask2=(mask2(box(2)+1:box(2)+box(4),box(1)+1:box(1)+box(3)));
   

    %define cellular regions (left pole, middel, right pole) 
    s=size(mask2);
    bin_edges=[0,s(2)/6,5*s(2)/6,s(2)];
    
    %isolate the points of the track and histogram using the cellular
    %region as bins
    [lrot,crot]=find(mask2==value+1);
    hist_freq=histc(crot,bin_edges);
    freq_tot=sum(hist_freq(1:3));


    %assign definite position if the track spend 80% of its time in a polar
    %or lateral region. Otherwise it is classified as "spurius".
    if hist_freq(1)/freq_tot>fraction
       position='pole';    %polar track
    end
    if hist_freq(3)/freq_tot>fraction
        position='pole';   %polar track
    end
    if hist_freq(2)/freq_tot>fraction 
        position='middle';
    end
    
    if ~exist('position') 
        position='s';   %spurius track
    end

    %compute the angles of the velocity over a certain threshold with the cell axis vector
    cell_axis=[cos(angle),sin(angle)];
    velocity_angles=[];
    for j=1:(length(x)-2)

        if ~isnan(x(j+1)) && ~isnan(y(j+1)) ... %&&~isnan(track.z(j)) ...
            && ~isnan(x(j)) && ~isnan(y(j)) % &&~isnan(track.z(j+1)) 

            vtemp=[x(j+1)-x(j),y(j+1)-y(j)]*pixelsize;                    
            if norm(vtemp)>150
                tmp_angle=dot(cell_axis,vtemp);
                tmp_angle=asind(tmp_angle/double(norm(vtemp)));            
                velocity_angles=vertcat(velocity_angles,tmp_angle);
            end
        else
            velocity_angles=vertcat(velocity_angles,0);
        end

    end

%compute single-track apparent diffusion coefficeint from displacements as in Vestergaars2014

    D_ST=compute_ST_D_and_s(singletrack);

%compute average displacement

    dx=x-circshift(x,[0,-1]);
    dy=y-circshift(y,[0,-1]);

    dx=dx(1:end-1);
    dy=dy(1:end-1);

    dX_temp=(sqrt(dx.^2+dy.^2))*pixelsize;

    if mean(dX_temp) > 0 &&  length(x)<100 %mean(dX_temp) < dXthreshold &&

        average_dx=mean(dX_temp);            

    end %if

 %{
    %redifine new tracks with additional information from fitting MSDs. r2
    %is kind of a measure of the goodness of fit
    track_info=struct('track',result_struct,'position',position,'hist',hist_freq, 'DiffConst', pF(1),'alpha',pF(2),'name',tr_name);%,'baseline',pF(3));
    newtrack.info=track_info;    
    newtrack.msds.r2=r2;
    newtrack.msds.resids=resids;
    newtrack.msds.cov=cov;
    newtrack.info.nPointsForFit=factor;
    newtrack.info.method=method;
    
    newtrack.mask=mask;
    
    %save the track to the savepath specified
    strcat(savepath,'INFOTRACK',tr_name)
    save(strcat(savepath,'INFOTRACK',tr_name),'newtrack');
    'ANDATOOOOOO'
    %save(strcat('..\track_info_140306_no_baseline4\','INFO',tr_name),'newtrack');
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MSD FITTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    msdfun=msdsFunctions_rep;

    [msd2,msd2_err]=msdfun.Overlapping(x,y,pixelsize);

    time=(1:length(msd2))*timedelay;

    if length(x)>10

        [p, r]=fitLinearMSD(x,y,time,msd2,[0,0.09],NpointsFit,0);
        D_err{i}=sqrt(diag(inv(r.R)*inv(r.R'))*(r.normr)^2/r.df);

        [pa]=fit_powerlaw_MSD(x,y,time,msd2,[0.01,0.5],NpointsFitAlpha,0);

    else
        p=[NaN, NaN];
        pa=[NaN,NaN];
        D_err=NaN;
    end

    

end