        
function [displacements,single_track_av_disp]=get_singlestep_displ_fromUtrack_outputnoNAN_rep(source_dir, Ndelays)

    displacements={};
    single_track_av_disp={};
    
    D=dir(strcat(source_dir,'\*.mat'));
    
    
    mytrack=[];    
    for i=1:length(D)
        
        tmptrack=importdata(strcat(source_dir,'\',D(i).name))
        
        if isfield(tmptrack,'tracksFinal')
            mytrack=[mytrack, tmptrack.tracksFinal];        
        else
            mytrack=[mytrack; tmptrack];
        end
        
    end
    
    track=filter_tracks(mytrack,5,100);
    figure
    for m=1:Ndelays
        disp=[];
        tmpdisp=[];
        for k=1:(length(track))
            
            x=track(k).tracksCoordAmpCG(1:8:end)*64;
            y=track(k).tracksCoordAmpCG(2:8:end)*64;
            if any(isnan(x))
                
                continue
            end
            for j=1:length(x)-m-1
                if ~isnan(x(j)) && ~isnan(y(j)) ... %&&~isnan(track.z(j)) ...
                    && ~isnan(x(j+m)) && ~isnan(y(j+m)) % &&~isnan(track.z(j+1))
                
                    %disp=sqrt( abs(x(j+1)-x(j))^2+abs(y(j+1)-y(j))^2 + abs(track.z(j+1)-track.z(j))^2 )
                    tmpdisp(j)=( abs(x(j+m)-x(j))^2+abs(y(j+m)-y(j))^2 );
                    disp=vertcat(disp, tmpdisp(j));
                end
            end
            single_track_av_disp{m}{k}=mean(tmpdisp(:));
            
        end
        displacements{m}=disp;
        cc=jet(m);
        Histo_setted(sqrt(cat(1,disp)),(0:25:550),'pmf',num2str(m),'y',cc(end,:),0)
    end
    displacements{end+1}=source_dir;
    single_track_av_disp{end+1}=source_dir;
end
