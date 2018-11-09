function varargout = all(varargin)
% ALL MATLAB code for all.fig
%      ALL, by itself, creates a new ALL or raises the existing
%      singleton*.
%
%      H = ALL returns the handle to a new ALL or the handle to
%      the existing singleton*.
%
%      ALL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALL.M with the given input arguments.
%
%      ALL('Property','Value',...) creates a new ALL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before all_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to all_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help all

% Last Modified by GUIDE v2.5 03-Nov-2018 23:48:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @all_OpeningFcn, ...
                   'gui_OutputFcn',  @all_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before all is made visible.
function all_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to all (see VARARGIN)

% Choose default command line output for all
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes all wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = all_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

gui_diff_circle;



% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)


diffraction_1slit;


setup1;
figure(1);
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)

diffraction_2slit;
setup1;
figure(1);

% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)

% op_rs_xy_02.m

% Rectangular aperture
% Numerical integration of the Rayleigh-Sommerfield diffraction integral of
% the first kind
% Uses cartesain coordinates / S.I. units
% Calculate of energy enclosed in circles
% Calls   fn_distancePQ.m
% SYMBOLS:  irradiance = intensity = energy density u [W.m^-2]
%           energy aperture --> observation screen  U [W or J/s]


% Uses functions
%     simpson2d.m  distancePQ.m

% 12 oct 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

% ************************************************************************
%clear all
close all
clc
all
tic


% INPUT PARAMETERS -------------------------------------------------------
num = 30;             % number for observation space
nP = num*4+1;         % observation points for P  format  integer * 4 + 1
nQ = 59;              % aperture points for Q  must be ODD

wL = 632.8e-9;        % wavelength [m]
 
ax = 1e-4;             % width - aperture space  [m]
ay = 1e-4;  

uQmax = 1e-3;                % incident energy density  [W.m^-2]

xPmax = 2e-2;         % half width - observation space   [m] 
yPmax = 2e-2;
zP = 1;                % distance: aperture plan to observation plane [m]
%zP = 3e-2;

% default values
% num = 40;             
% nP = num*4+1;         % observation points for P  format  integer * 4 + 1
% nQ = 99;              % aperture points for Q  must be ODD
% wL = 632.8e-9;        
% xQmax = 1e-4;        
% yQmax = 2e-4;  
% uQmax = 1e-3;               
% xPmax = 2e-2;          
% yPmax = 2e-2;
% zP = 1;               


% SETUP -----------------------------------------------------------------
cL = 2.99792458e8;      % speed of light
eps0 = 8.854187e-12;    % permittivity of free space
nR = 1;                 % refractive index

k = 2*pi/wL;            % propagation constant
ik = 1i*k;              % j k

% Aperture space
xQmin = -ax/2;  xQmax = ax/2;
yQmin =  -ay/2; yQmax =  ay/2;

unit = ones(nQ,nQ);    % unit matrix
rPQ = zeros(nQ,nQ); rPQ3 = zeros(nQ,nQ);
MP1 = zeros(nQ,nQ); MP2 = zeros(nQ,nQ); kk = zeros(nQ,nQ); 
MP = zeros(nQ,nQ);

EQmax = sqrt(2*uQmax/(cL*nR*eps0));    % Electric field within aperture
EQ = EQmax .* ones(nQ,nQ);  

zQ = 0;
xQ1 = linspace(xQmin,xQmax,nQ);    
yQ1 = linspace(yQmin,yQmax,nQ);
[xQ, yQ] = meshgrid(xQ1,yQ1);
a = max([ax ay]);
d_RL = a^2/wL;                        % Rayleigh distance                  

uQ = (cL*nR*eps0/2) .* EQ .*EQ;       % energy density within aperture
UQ = simpson2d(uQ,xQmin,xQmax,yQmin,yQmax);    % energy transmitted from aperture per s [J/s]
UQtheory = uQmax * (xQmax-xQmin) * (yQmax - yQmin);

% Observation space
xPmin = -xPmax;
yPmin =  -yPmax;

EP = zeros(nP,nP);         % electric field

xP1 = linspace(xPmin,xPmax,nP);     
yP1 = linspace(yPmin,yPmax,nP);
[xP, yP] = meshgrid(xP1,yP1);
indexXY = num*2+1;               % defines X and Y axes  y = 0 and x = 0;

%          optical coordinates
vPx = (pi*(xQmax - xQmin)/wL) .* xP ./ sqrt(xP.^2 + zP^2);  
vPy = (pi*(yQmax - yQmin)/wL) .* yP ./ sqrt(yP.^2 + zP^2);


%          first min xP and yP values
thetax1 = asin(wL/ax); thetay1 = asin(wL/ay);
x1 = zP * tan(thetax1); y1 = zP * tan(thetay1);
 


% COMPUATION OF DIFFRACTION INTEGRAL ------------------------------------

% Electric Field EP - observation space   P
for c1 = 1 : nP
    for c2 = 1 :nP
       rPQ = fn_distancePQ(xP(c1,c2),yP(c1,c2),zP,xQ,yQ,zQ);
       rPQ3 = rPQ .* rPQ .* rPQ;
       kk = ik .* rPQ;
       MP1 = exp(kk);
       MP1 = MP1 ./ rPQ3;
       MP2 = zP .* (ik .* rPQ - unit);
             
       MP = EQ .* MP1 .* MP2;
       EP(c1,c2) = simpson2d(MP,xQmin,xQmax,yQmin,yQmax) ./(2*pi) ;
   end 
end

% Irradiance (energy density) u /  energy U  -   observation space   P 
uP = (cL*nR*eps0/2) .* abs(EP).^2;
uPmax = max(max(uP));
uPdB = 10 .* log10(uP./uPmax);

UP = simpson2d(uP,xPmin,xPmax,yPmin,yPmax);


% percentage of energy enclosed within a circle of radius of 1st min in x direction
Rc = 3.831*2;
uPenclosed = uP;
uPenclosed(vPx > Rc) = 0;
UPenclosed = simpson2d(uPenclosed,xPmin,xPmax,yPmin,yPmax);
UPenclosed_perc = 100*UPenclosed/UP;

%theta_c = asin(Rc/(k*ax));     % ax effective diameter of circular aperture
%xPc = zP * tan(theta_c);

% Draw circle
nc = 500;
theta = linspace(0,2*pi,nc);
xc = x1 .* cos(theta); yc = x1 .* sin(theta);


% GRAPHICS --------------------------------------------------------------
figure(1)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.6 0.6]);
fs = 12;
%tx = 'x_P (blue)  or y_P (red)  [m]';
%ty = 'irradiance  (a.u.)';

% energy density   uP vs Xp & yP
subplot(2,2,1);
tx = 'x_P (blue)  or y_P (red)  [m]';
ty = 'energy density  u  [W.m^{-2}]';
x = xP(indexXY,:); y = uP(indexXY,:);
plot(x,y ,'linewidth',2);
hold on
x = yP(:,indexXY); y = uP(:,indexXY);
plot(x,y,'linewidth',2);
xlabel(tx);   ylabel(ty);

% energy density   uPdB  vs xP & yP
subplot(2,2,3);   
tx = 'x_P (blue)  or y_P (red)  [m]';
ty = 'energy density  u  [dB]';
x = xP(indexXY,:);
y = uPdB(indexXY,:);
plot(x,y ,'linewidth',2);
hold on
x = yP(:,indexXY);
y = uPdB(:,indexXY);
plot(x,y,'linewidth',2)
xlabel(tx);   ylabel(ty);

% energy density   uP vs vPx & vPy
subplot(2,2,2);
tx = 'optical coordinate  v_P / \pi';
ty = 'energy density  u  [a.u.]';
x = vPx(indexXY,indexXY:end)./pi;
y = uP(indexXY,indexXY:end)./uPmax;
plot(x,y ,'linewidth',2);
xlabel(tx);   ylabel(ty);
grid on

% energy density   uPdB  vs xP & yP
subplot(2,2,4);   
tx = 'optical coordinate   v_P / \pi';
ty = 'energy density  u  [dB]';
x = vPx(indexXY,indexXY:end)./pi;
y = uPdB(indexXY,indexXY:end)./uPmax;
plot(x,y ,'linewidth',2);
xlabel(tx);   ylabel(ty);
grid on

figure(2)      % photograph plot  /  [3D] plot
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.3 0.3]);
pcolor(xP,yP,10.*(uP).^0.3);
shading interp
axis equal
colormap(gray)
axis off
set(gcf,'color','k')

hold on
x = xc; y = yc;
plot(x,y,'y','lineWidth',2);

figure(3)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.3 0.3 0.3 0.3]);
surf(xP,yP,10.*(uP).^0.3);
shading interp
axis equal
colormap(jet)
axis off
%set(gcf,'color','k')
xlabel('x');
ylabel('y')


figure(4)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.3 0.3 0.15 0.2]);
%tx = 'x_P (blue)  or y_P (red)  [m]';
tx = 'x_P  [m]';
%ty = 'energy density  u  [W.m^{-2}]';
ty = 'u  [W.m^{-2}]';
tm1 ='    z_P  =   ';
tm2 = num2str(zP,3);
tm3 = '  m';
tm = [tm1 tm2 tm3];
x = xP(indexXY,:); y = uP(indexXY,:);
plot(x,y ,'b','linewidth',2);
%hold on
%x = yP(:,indexXY); y = uP(:,indexXY);
%plot(x,y,'r','linewidth',2);
xlabel(tx);   ylabel(ty);
title(tm);
grid on
%set(gca,'Ylim',[0 1e-6]);
%set(gca,'Xlim',[-4e-4 4e-4]);


% COMMAND WINDOW OUTPUTS ------------------------------------------------
disp('Parameter summary  [SI units]');
fprintf('wavelength [m]  =  %3.3g \n',wL);
fprintf('nQ  =  %3.3d \n',nQ);
fprintf('nP  =  %3.3d \n',nP);
disp('  ')
disp('Aperature Space');
fprintf('X width [m] =  %3.3e \n',ax);
fprintf('Y width [m]  =  %3.3e \n',ay);
fprintf('energy density [W/m2] uQmax  =  %3.3e \n',uQmax);
fprintf('energy from aperture [J/s]   UQ(theory) = %3.3e \n',UQtheory);
fprintf('energy from aperture [J/s]   UQ(calculated) = %3.3e \n',UQ);
disp('  ')
disp('Observation Space');
fprintf('X width [m] =  %3.3e \n',2*xPmax);
fprintf('Y width [m]  =  %3.3e \n',2*yPmax);
fprintf('distance aperture to observation plane [m]   zP = %3.3e \n',zP);
fprintf('Rayleigh distance  [m]   d_RL = %3.3e \n',d_RL);
fprintf('Fraunhofer: position of 1st min in X direction [m]  =  %3.3e \n',x1);
fprintf('Fraunhofer: position of 1st min in Y direction [m]  =  %3.3e \n',y1);
disp('  ');
fprintf('energy to aperture [J/s]   UP = %3.3e \n',UP);
fprintf('max energy density  [W./m2]   uPmax = %3.3e \n',uPmax);
disp('   ');
fprintf('percentage energy enclosed within circle of radius (vPx = 3.14159)  = %3.1f \n',UPenclosed_perc);
disp('   ');


toc




% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
clc
clear
close all;
M=256;
N=256;
sizeA=100;
sizeB=50;
binPic1=zeros(M,N);
binPic2=zeros(M,N);
binPic1((M-sizeA)/2+1:(M-sizeA)/2+1+sizeA,(M-sizeA)/2+1:(M-sizeA)/2+1+sizeA)=1;
binPic2((N-sizeB)/2+1:(N-sizeB)/2+1+sizeB,(N-sizeB)/2+1:(N-sizeB)/2+1+sizeB)=1;
binPic1=uint8(binPic1);
binPic2=uint8(binPic2);
figure,
subplot(221),imshow(binPic1,[]),title('Slit 1');
subplot(222),imshow(binPic2,[]),title('Slit 2');
binPic1fft=mySpectrum(binPic1);             
binPic2fft=mySpectrum(binPic2);
subplot(223),imshow(binPic1fft),title('Diffraction pattern of 1');
subplot(224),imshow(binPic2fft),title('Diffraction pattern of 2');

[z,Fx,Fy]=mesh3D(binPic1fft);       
figure,mesh(Fx,Fy,z);title('3D plot of 1');
set(gcf,'unit','centimeters','position',[1,2,20,15])
z1=mesh3D(binPic2fft);              
figure,mesh(Fx,Fy,z1);title('3D plot of 2');
set(gcf,'unit','centimeters','position',[30,2,20,15])


% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
A = imread('aperture','jpg');
s=size(A);
A=A(:,:,1);
A=I;
		[M,N] = size(I);
		lambda = 600e-9;
		f = 16.5e-3;
		W = 4e-3;
		I = (fftshift(fft2(A))).^2;
		x = linspace(-lambda*f*N/(2*W),  lambda*f*N/(2*W), N);
		y = linspace(-lambda*f*M/(2*W),  lambda*f*M/(2*W), M);
		pcolor(x, y, abs(I));
        
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
%% edge diffraction


%% plot edge function
% a = \lambda z

% normalize to 1
a = 1;
N=401;

x = linspace(-2,10,N+1);

[cc,ss] = fresnel(sqrt(2/a)*x);

z = sqrt(1/(2*j))*complex(0.5+cc,0.5+ss);


f = z.*conj(z);
figure;
subplot(4,1,1)
plot(x,f,'k','LineWidth',1.5);
%axis([-xrange xrange 0 1.5]);
xlabel('x');
%ylabel('E(x)');

%% edge image

x = linspace(-2,8,512);
edge = edge_fcn(x);
edge = edge/max(max(edge));
edge = ones(size(x))'*edge;

subplot(4,1,2)
imshow(edge);
%imwrite(edge,'edge.jpg');


%% edge wavefront


a = 1;
N=400;

x = linspace(-1,5,N+1);

[cc,ss] = fresnel(sqrt(2/a)*x);

z = sqrt(1/(2*j))*complex(0.5+cc,0.5+ss);


g = z.*conj(z);

f =  - unwrap(angle(z))/(2*pi) - 1;

subplot(4,1,3);
plot(x,g,'k','LineWidth',1.5);
ylabel('irradiance');

subplot(4,1,4);
plot(x,f,'k','LineWidth',1.5);
%axis([-xrange xrange 0 1.5]);
xlabel('x');
ylabel('optical phase (waves)');
%ylabel('E(x)');
set(gcf,'unit','centimeters','position',[1,2,40,40])

% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
op_rs_rxy_cross;
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
