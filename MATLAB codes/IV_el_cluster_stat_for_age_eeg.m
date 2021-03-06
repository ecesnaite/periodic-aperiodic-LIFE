% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in 'Alterations in rhythmic and non-rhythmic resting-state
% EEG activity and their link to cognition in older age' paper.
% The code estimates correlations between age and PSD parameters. Then runs
% cluster statistics to with permutation tests.
% The code was addapted from the version created together with Dr Keyvan
% Mahjoory
% Last updated 22.06.2021

clear all, close all

saveDir = '';
resDir = ''; % data directory

eeglab

load(''); % laod 2D channel locations, channel labels, and channel neighbour information

cfg.chanLocs2D = new_2Dlocs_LIFE(:,2:3); % required ONLY for plotting

%%-------------------------------------------------------------------------
%% Initial settings of analysis

resTypeAll= {'power'};%,{'slope','power','iaf', 'theta'}

%regress other parameters
regOutPower = 0;
regOutSlope = 1;
regOutTheta = 1;
regOutIAF = 1;
regOutSex = 1;
regOutEdu = 1;


cfg.nRnd = 1000;   % number of randomization times
cfg.sigThresh = 0.05; % significance level

makeTopoPlot =1;
rmlH = 1; %Right handed subjects

%% =============  Loading Data  ===========================
load([resDir,''])%

if strcmp(resTypeAll,'power')
    resMatrix = alpha_pow; 
elseif strcmp(resTypeAll,'iaf')
    resMatrix =  freq; 
elseif strcmp(resTypeAll,'slope') 
    resMatrix = slope;
elseif strcmp(resTypeAll,'theta')  
    resMatrix = theta_pow;   
end

%% =============  Loading Factors, age, sex and education info  ==================

questData0 = str2double(IDs_matched(:,2)); % participant age information

% if there is data to regress out

datD = [];
regOutLabel = ''
if regOutPower
    datD = [datD alpha_pow];
    regOutLabel = 'Power';
    fprintf('Power was regressed out. \n')
end

if regOutSlope
    datD = [datD slope];
    regOutLabel = [regOutLabel, ', ','Slope'];
    fprintf('Slope was regressed out. \n')
end
if regOutIAF
    datD = [datD freq];
    regOutLabel = [regOutLabel, ', ','IAF'];
    fprintf('IAF was regressed out. \n')
end

if regOutTheta
    datD = [datD theta_pow];
    regOutLabel = [regOutLabel, ', ','Theta'];
    fprintf('Theta was regressed out. \n')
end

datC = []
regAgeLabel = ''

if regOutEdu
    datC = [datC edu'];
    regAgeLabel = [regAgeLabel, ', ', 'Edu'];
    fprintf('Education was regressed out. \n')
end
if regOutSex
    datC = [datC sex_cog_code];
    regAgeLabel = [regAgeLabel, ', ','Sex'];
    fprintf('Sex was regressed out. \n')
end

datA = questData0
datB.vals = resMatrix;
datB.chanLabels = chanLabels;
datB.chanNeighbours = {neighbours.label};
datB.chanNeighbours(2,:) = {neighbours.neighblabel};

%% =============  Correlation & cluster analysis  =========

[clNoPerm, CL, indxPerm] = el_k1_cluster_calc_sensor_LIFE(datA, datB, normalize(datC), normalize(datD), cfg);
plotTitle = ['corr(',char(resTypeAll), ' - age | ', char(regAgeLabel),', ',char(regOutLabel),') '];

if ~isempty(CL) % if cluster is found
    CLCell = struct2cell(CL');
    pvalClust = length(find( abs(cell2mat(CLCell(2,:))) > abs(clNoPerm.tmax) ))/cfg.nRnd;
    plotTitle = [plotTitle,'  P= ' num2str(pvalClust)];
end

