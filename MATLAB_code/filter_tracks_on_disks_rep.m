%This function takes a set of tracks and a masked image of the microlens array (where each spots is a disk of diameter ~1um) as inpput. Then checks if the starting point of each track is within a disk of the masked ML image. If so the track is retained otherwise discarde.
%
%INPUTS
%	tracksin : input tracks
%	im_disk  : disk image coming from the output of the 'create_disks_Microlens_array' function


function [tracksout] = filter_tracks_on_disks_rep(tracksin, im_disks)

	inds=[];

	for i = 1:length(tracksin)

		tracksin(1)

		x = floor(tracksin(i).tracksCoordAmpCG(1:8:end));
		y = floor(tracksin(i).tracksCoordAmpCG(2:8:end));
        if x==0 | x<0
            x==1;
        end
        if y==0 | y<0
            y==1;
        end
        
        x
        y
        
		if ~im_disks(y(1),x(1))
			inds=horzcat(inds,i);
		end

	end

	tracksin(inds)=[];
	tracksout=tracksin;

end