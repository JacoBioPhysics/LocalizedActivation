

%this function takes a utrack output track structure and returns only those tracks that are produced after the activation events.

%	INPUTS:

%			tracks: 	utrack output (struct)
%			ActEvent:   frames interval at which an activation event occurs (int)
%			tolerance:  allowed number of frames after activation to keep the track (int)
%	OUTPUT:
%
%			filtered_tracks (struct)

function [filtered_tracks]=filter_tracks_after_activation_rep(tracks, ActEvent, tolerance)

		if nargin<3
			ActEvent=249;
			tolerance=20;
        end
    
		events=[tracks(:).seqOfEvents];
		begins=events(1,1:4:end);
		filtered_tracks=tracks(mod(begins,ActEvent)<tolerance);

end



