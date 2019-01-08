
%% movie information
%
% Copyright (C) 2014 LCCB 
%
% This file is part of u-track.ss
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
movieParam.imageDir = all_target_folders{h}; %directory where images are

d=dir('*.tif');

firstim=d(1).name;
lastim=d(end).name;
movieParam.filenameBase = '\Image'; %image file name base
movieParam.firstImageNum = str2num(firstim(regexp(firstim,'\d'))); %number of first image in movie
movieParam.lastImageNum = str2num(lastim(regexp(lastim,'\d'))); %number of last image in movie
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

detectionParam.psfSigma = 2.3; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.1,'alphaA',0.1,'alphaD',0.1,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 0; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.1; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

detectionParam.calcMethod = 'g';

%absolute background info and parameters...
background.imageDir = 'C:\Users\solari\Desktop\';

background.filenameBase = 'BCK_ANDOR_emGAIN100.tif';
background.alphaLocMaxAbs = 0.1;
detectionParam.background = background;

%% additional input

%saveResults
mkdir(strcat(movieParam.imageDir,'\detection_output'));
saveResults.dir = strcat(movieParam.imageDir,'\detection_output'); %directory where to save input and output
saveResults.filename = 'detectionTest1.mat'; %name of file where input and output are saved
% saveResults = 0;

%verbose state
verbose = 1;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_rep(movieParam,detectionParam,saveResults,verbose);
