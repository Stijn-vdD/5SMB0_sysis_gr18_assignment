%% SETUP FOR PLOTS

% Declare this top header once.

%Sets the units of your root object (screen) to centimeters
set(0,'units','centimeters');
%Obtains this inch information
fig_setup = struct;

fig_setup.CM_SS = get(0,'screensize');
fig_setup.width = fig_setup.CM_SS(3);
fig_setup.height= fig_setup.CM_SS(4);
fig_setup.fig_wd = 9.8; %width of figure in cm
fig_setup.fig_wd_wide = 15;
fig_setup.fig_hgt= 9.8; %height of figure in cm
fig_setup.fntsize = 9; %fontsize in px
 
fig_setup.export_figures = 1;
% 'vector' for pdf
% 'image'  for png
fig_setup.img_format = 'vector';

color(1) = "#0047AB"; %cobalt blue
color(2) = "#FF5733"; %tomato
color(3) = "#13BF1D"; %dark pastel green
color(4) = "#A41495"; %purple
color(5) = "#E9E16D"; %dark yellow
color(6) = "#838383"; %dark grey
color(7) = "#87A5CF"; %blue-light

if strcmp(fig_setup.img_format,'vector')
    fig_setup.img_ext = ".pdf";
else
    fig_setup.img_ext = ".png";
end

if ~isfolder('output-figures')
    mkdir('output-figures');
end