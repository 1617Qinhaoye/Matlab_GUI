% gui_diff_circle.m

% GUI for the simulation of the diffraction pattern for a circular aperture
% INPUTS:
%       wavelegnth in namometers (nm) 
%       aperture radius in millimeters (mm)
%       distance between the aperture & observation planes (mm)
%       max radial distance (mm)


%clear all
clc
close all

f1 = figure('Color',[1 1 1],'Name','Diffraction', ...
     'NumberTitle','off');
% set(gcf,'units','normalized','position',[0.05 0.05 0.85 0.85]); 
set(gcf,'units','normalized') 
set(gcf,'position',[0.01 0.05 0.85 0.85]); 

Text_heading = uicontrol(gcf,'Style','text','units','normalized','Position',[0.05 0.92 0.4 0.05], ...
          'String','CIRCULAR','FontSize',20,'HorizontalAlignment','righ', ...
          'BackgroundColor',[1 1 1],'FontWeight','bold', ...
          'ForegroundColor',[0 0 0]);

% Input Initial Data
boxA = 632.8; boxB = 0.10; boxC = 1000; boxD = 15;

fs = 10;

Text_a = uicontrol(gcf,'Style','text','units','normalized','Position',[0.01 0.82 0.4 0.1], ...
          'String','    Wavelength (400-750)  ','FontSize',fs,'HorizontalAlignment','left', ...
          'BackgroundColor',[1 1 1]);
      
Edit_boxA = uicontrol(gcf,'Style','edit','units','normalized','Position',[0.03 0.82 0.05 0.05], ...
        'String',boxA,'FontSize',fs,'BackgroundColor',[1 1 1], ...
        'Callback','boxA = str2num(get(Edit_boxA,''String''));');
      
Text_b = uicontrol(gcf,'Style','text','units','normalized','Position',[0.18 0.82 0.4 0.1], ...
          'String','      Radius (0.1-10)  ','FontSize',fs,'HorizontalAlignment','left', ...
          'BackgroundColor',[1 1 1]);
      
Edit_boxB = uicontrol(gcf,'Style','edit','units','normalized','Position',[0.20 0.82 0.05 0.05], ...
        'String',boxB,'FontSize',fs,'BackgroundColor',[1 1 1], ...
        'Callback','boxB = str2num(get(Edit_boxB,''String''));');
    
Text_c = uicontrol(gcf,'Style','text','units','normalized','Position',[0.36 0.82 0.4 0.1], ...
          'String','   Distance (0.1-1000)  ','FontSize',fs,'HorizontalAlignment','left', ...
          'BackgroundColor',[1 1 1]);
      
Edit_boxC = uicontrol(gcf,'Style','edit','units','normalized','Position',[0.38 0.82 0.05 0.05], ...
        'String',boxC,'FontSize',fs,'BackgroundColor',[1 1 1], ...
        'Callback','boxC = str2num(get(Edit_boxC,''String''));');    
    
Text_d = uicontrol(gcf,'Style','text','units','normalized','Position',[0.62 0.82 0.4 0.1], ...
          'String','        Scale(0.1-100)  ','FontSize',fs,'HorizontalAlignment','left', ...
          'BackgroundColor',[1 1 1]);
      
Edit_boxD = uicontrol(gcf,'Style','edit','units','normalized','Position',[0.64 0.82 0.05 0.05], ...
        'String',boxD,'FontSize',fs,'BackgroundColor',[1 1 1], ...
        'Callback','boxD = str2num(get(Edit_boxD,''String''));');    
    
   
pushbutton_run = uicontrol(gcf,'Style','pushbutton','units','normalized','Position',...
     [0.7 0.93 0.05 0.05], 'FontSize',12,'FontWeight','bold', 'String','Start', ...
     'CallBack', 'gui_diff_circle_cal');
 
 
 pushbutton_close = uicontrol(gcf,'Style','pushbutton','units','normalized','Position',...
     [0.8 0.93 0.08 0.05], 'FontSize',12,'FontWeight','bold', 'String','Back to Main', ...
     'CallBack', 'close1');
 
 
 
      