
function [Dp,Dm, sp, sm, Dp_perdataset, Dm_perdataset, sp_perdataset, sm_perdataset, Ds_perdataset,ss_perdataset,D, D_perdataset ] = computeDandS( loc_act_class, track_length_thresh )
%This function takes a localized activation class and for each infotracks
%in each mobility structure it computes the total distribution of apparent diff
%coefficients and localization errors.

    if nargin < 2
        track_length_thresh=5;
    end  
    
    SNRs=[];

    Ndataset=length(loc_act_class.mobility_structures);

    
    D={};%all tracks diff coeff

    Dp={};%polar tracks diff coeff
    Dm={};%middle tracks diff coeff
    Ds={};%spurius tracks diff coeff

    sp={};%polar tracks error on d
    sm={};%middle track error on D
    ss={};%spurius tracks error on D

    for i=1:Ndataset
        
        Dp{i}(1)=0;
        Dm{i}(1)=0;
        sp{i}(1)=0;
        sm{i}(1)=0;
        ss{i}(1)=0;
        D{i}(1)=0;
        Ds{i}(1)=0;
        
        Nsubdataset=length(loc_act_class.mobility_structures{i}.infotracks);

        for j=1:Nsubdataset

            Dp_perdataset{i,j}(1)=0;
            Dm_perdataset{i,j}(1)=0;
            Ds_perdataset{i,j}(1)=0;           
            sp_perdataset{i,j}(1)=0;
            ss_perdataset{i,j}(1)=0;                               
            sm_perdataset{i,j}(1)=0;            
            
            D_perdataset{i,j}(1)=0;
            
            mytracks=filter_tracks_noNAN_rep(loc_act_class.mobility_structures{i}.infotracks{j},track_length_thresh,1000);
            Ntracks=length(mytracks);
            
            for k=1:Ntracks
                
                track=loc_act_class.mobility_structures{i}.infotracks{j}(k);

                %computes the ST diffusion coeff and error
                [d,s,SNR]=compute_ST_D_and_s(track);
                SNRs=[SNRs,SNR];    

                %assigns d and the error to the proper variable (pole, middle  or spurius)
                if track.position(1)=='p' && d>0 
                    Dp{i}=[Dp{i},d];
                    sp{i}=[sp{i},s];
                    Dp_perdataset{i,j}=[Dp_perdataset{i,j},d];
                    sp_perdataset{i,j}=[sp_perdataset{i,j},s]; 
                    loc_act_class.mobility_structures{i}.infotracks{j}(k).D_ST=d;
                elseif track.position(1)=='m' && d>0
                    Dm{i}=[Dm{i},d];
                    sm{i}=[sm{i},s];
                    Dm_perdataset{i,j}=[Dm_perdataset{i,j},d];
                    sm_perdataset{i,j}=[sm_perdataset{i,j},s];
                    loc_act_class.mobility_structures{i}.infotracks{j}(k).D_ST=d;
                elseif track.position(1)=='s' && d>0
                    Ds{i}=[Ds{i},d];
                    ss{i}=[ss{i},s];
                    Ds_perdataset{i,j}=[Ds_perdataset{i,j},d];
                    ss_perdataset{i,j}=[ss_perdataset{i,j},s];                   
                    loc_act_class.mobility_structures{i}.infotracks{j}(k).D_ST=d;
                else
                    loc_act_class.mobility_structures{i}.infotracks{j}(k).D_ST=NaN;
                end

                if d>0
                    D{i}=[D{i},d];
                    D_perdataset{i,j}=[D_perdataset{i,j},d];
                end
                                
            end%for k
            

        end%for Nsubdataset
        

 
    end %for i Ndataset

    

    figure
    hist(SNRs,50);


end%function

function [Dp,Dm, sp, sm, Dp_perdataset, Dm_perdataset, sp_perdataset, sm_perdataset, Ds_perdataset,ss_perdataset,D_perdataset ] = computeDandS( loc_act_class, track_length_thresh )
%This function takes a localized activation class and for each infotracks
%in each mobility structure it computes the total distribution of apparent diff
%coefficients and localization errors.

    if nargin < 2
        track_length_thresh=5;
    end  
    
    Ndataset=length(loc_act_class.mobility_structures);
    Dp={};
    Dm={};
    D={};
    sp={};
    sm={};
    Ds={};
    ss={};

    for i=1:Ndataset
        
        Dp{i}(1)=0;
        Dm{i}(1)=0;
        sp{i}(1)=0;
        sm{i}(1)=0;
        ss{i}(1)=0;
        D{i}(1)=0;
        Ds{i}(1)=0;
        
        Nsubdataset=length(loc_act_class.mobility_structures{i}.infotracks);
        
        for j=1:Nsubdataset

            Dp_perdataset{i,j}(1)=0;
            Dm_perdataset{i,j}(1)=0;
            Ds_perdataset{i,j}(1)=0;           
            sp_perdataset{i,j}(1)=0;
            ss_perdataset{i,j}(1)=0;                               
            sm_perdataset{i,j}(1)=0;            
            
            D_perdataset{i,j}(1)=0;
            
            mytracks=filter_tracks_noNAN_rep(loc_act_class.mobility_structures{i}.infotracks{j},track_length_thresh,1000);
            Ntracks=length(mytracks);
            
            for k=1:Ntracks
                
                track=loc_act_class.mobility_structures{i}.infotracks{j}(k);
                [d,s]=compute_ST_D_and_s(track);
                
                if track.position(1)=='p' && d>0 
                    Dp{i}=[Dp{i},d];
                    sp{i}=[sp{i},s];
                    Dp_perdataset{i,j}=[Dp_perdataset{i,j},d];
                    sp_perdataset{i,j}=[sp_perdataset{i,j},s];                    
                elseif track.position(1)=='m' && d>0
                    Dm{i}=[Dm{i},d];
                    sm{i}=[sm{i},s];
                    Dm_perdataset{i,j}=[Dm_perdataset{i,j},d];
                    sm_perdataset{i,j}=[sm_perdataset{i,j},s];
                elseif track.position(1)=='s' && d>0
                    Ds{i}=[Ds{i},d];
                    ss{i}=[ss{i},s];
                    Ds_perdataset{i,j}=[Ds_perdataset{i,j},d];
                    ss_perdataset{i,j}=[ss_perdataset{i,j},s];                    
                end

                if d>0
                    D{i}=[D{i},d];
                    D_perdataset{i,j}=[D_perdataset{i,j},d];
                end
                                
            end
            

        end
        

 
    end %for i Ndataset

    
end%function


