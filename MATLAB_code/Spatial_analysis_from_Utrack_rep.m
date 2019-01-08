function [ tracks ] = Spatial_analysis_from_Utrack_rep( tracks, mask, timedelay,pixelsize,NpointsFit,NpointsFitAlpha )

%Spatial_analysis_from_Utrack_rep: this function assigns a cellular position (pole, middle or spurius) to each track in a track set

% 	INPUTS: 

%		tracks: track set output of Utrack tracking routine
%		mask  : segmented bright field image from BRIGHT_FIELD_SEGMENTATION Wtershed function (Python)

%	OUTPUTS:

%		tracks: tracks with additional information on position 		
	tracks
	
	tracks=filter_tracks_on_mask_rep(tracks,mask)

%{
	for m=1:max(max(mask))
		m
		[a,b]=find(mask==m)
		if any(a) == 511 || any(b)==511 || any(a)==1 || any(b)==1
			mask(mask==m)=0;
		end
	end
	imagesc(mask)
%}
	for i=1:length(tracks)

	
		if mod(i, 100)==0
			i	
		end

		[myangle, box, mask2, position, bins,maskbox,msd2, msd_err,pa,p,D_err,velocity_angles,D_ST,average_dx] = singletrack_spat_analysis_rep(tracks(i),mask,i,timedelay,pixelsize,NpointsFit,NpointsFitAlpha);
		
		if isnan(mask2)
			tracks(i).position=NaN;
			continue
		end

		tracks(i).position=position;
		tracks(i).bins=bins;
		tracks(i).x=tracks(i).tracksCoordAmpCG(1:8:end);
		tracks(i).y=tracks(i).tracksCoordAmpCG(2:8:end);
		tracks(i).FrameNumbers=[ tracks(i).seqOfEvents(1,1) , tracks(i).seqOfEvents(2,1) ];
		tracks(i).MSD2=msd2;
		tracks(i).MSD2_err=msd_err;
		tracks(i).Da=pa(1);
		tracks(i).a=pa(2);		
		tracks(i).D=p(1);
		tracks(i).D_err=D_err;
		tracks(i).velocity_angles=velocity_angles;
		tracks(i).D_ST=D_ST;
		tracks(i).average_dx=average_dx;


	end



end

