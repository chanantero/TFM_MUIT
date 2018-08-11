%% Experiment15

% Created on 18th July 2018.
% First use of SVTtoolbox 
% Documentation: http://matlab.sfstoolbox.org/en/2.4.2/

%% Preamble
pathSetUp;

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% Examples of use of SFStoolbox with native secondary source arrays

%% Creation of an array with the GTAC's loudspeaker array geometry and simulation

% "conf" will be the most important variable. It is a structure that stores
% all relevant information about how simulations are carried out. It must
% be passed as an input argument to all simulating functions. SFS_config is
% a function that returns a default version of this structure. The user can
% change the values of the fields of conf. However, the user should not
% edit SFS_config.
conf = SFS_config; 

% Create WFS array positions, orientations and weights
s = WFSToolSimple.generateScenario(96, 'originReference', 'roomCorner');
conf.secondary_sources.geometry = 'custom';
conf.secondary_sources.x0 = [s.loudspeakersPosition, s.loudspeakersOrientation, ones(96, 1)];
x0 = secondary_source_positions(conf);
ax = axes(figure);
figsize(540,404,'px');
conf.plot.realloudspeakers = true;
draw_loudspeakers(x0,conf);
ax.XLim = [s.roomPosition(1), s.roomPosition(1) + s.roomPosition(3)];
ax.YLim = [s.roomPosition(2), s.roomPosition(2) + s.roomPosition(4)];

