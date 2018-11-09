% fd_02.m

% [1D] FDTD simulation of an electromagnetic wave produced by a point
% source.
% The Finite Difference Time Domain method uses a centre-difference
% representation of the continuous partial differential equations to create
% an iterative numerical model of the progagation of electromagnetic waves.
% [ ... ]  used to state default values and units
% Units for time:   tStep   convert to m:  1 m = (1/zStep) zSteps
% Units for Z grid: zStep   convert to s:  1 s = (1/tStep) tSteps

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 180308

tic
close all
clear all
clc

% =======================================================================
% INPUTS   [defaults  units]
% =======================================================================
% Number of time steps   [1200 tSteps]
Nt = 500;
% Number of grid points for Z space [400 zSteps]
Nz = 150;
% Wavelength EM wave to determine cell size  [3e-2 m]
lambda = 80;

% SOURCE 
% Type: 1 (pulse)   2  (sinusoidal)   3 (modulated sinusoidal)
flagS = 1;
% Location [10 zSteps    zS > 2]
%   zS = round(Nz/2);
zS = 2;
% Amplitude   [1]   
A = 1;
% Width of Gaussian pulse   [pulse 12 modulated  25    tSteps]   
width = 12;
% Time centre for Gasussian pulse   [40 tSteps]
centre = 100;

% MEDIUM 1 
% Relative permittivity   [1]
eR1 = 1; 
% Conductivity   [0  S/m (siemens / metre]     
S1  = 0;  

% MEDIUM 2 
% Relative permittivity   [1]
eR2 = 4; 
% Conductivity   [0  S/m]     
S2  = 0.0;
% Start index Medium 2   [round(Nz/2)]   
M2 = round(Nz/2); 
% set Z grid points for dielectric medium
indexR = 50:50+20+10  ;

% Probe positions to monitor Esx as a function of time  [zSteps]
%cP = round([1 2 3 4 5].*Nz./6);
%cP = [60 80 100 120 150];
cP = [ 36  66 132];
% Scaling Y axis for Esx [1.1] and Hy [1.1] 
limE = 2.1;
limH = 2.1;
 
% BOUNDARY CONITIONS -
% flagBC = 1   Absorbing at both boundaries
% flagBC = 2   Perfect electric conductor PEC at both boundaries
% flagBC = 3   Perfect magnetic conductor PMC at both boundaries
flagBC = 1;  
   
% SETUP FOR SAVING ANIMATION 
% flag: 0 animated gif NOT saved or 1 saved /  file name
f_gif = 1;
ag_name = 'agFD01A.gif';
% Delay in seconds before displaying the next image  
delay = 0.15;  
% Frame counter start
nF = 1; 

% ======================================================================
% CALCULATIONS
% ======================================================================
% CONSTANTS 
   eps0 = 8.85e-12;     % permittivity of free space
   mu0 = 4*pi*1e-7;     % permeability of free space
   c0 = 3e8;            % speed of light
   D = 4;               % Stability factor

 % SETUP 
  dz = lambda / 40;        % cell size or grid spacing [m]
  dt = dz / (sqrt(D)*c0);  % time step (time increment)  [s]
  t = 1:Nt;                % time grid
  z = 1:Nz;                % Z space grid;
  f = c0/lambda;           % frequency of EM wave   [Hz]
  T = 1/f;                 % period of EM wave   [s]
  omega = 2*pi*f;          % angular frequency of EM wave   [1/s]
  
  % SOURCE    
  switch flagS
      case 1     % pulse 
          if zS == 1; zS = 2; end
        source = A.*exp(-0.5.*((centre - t)/width).^2);
      case 2     % sinusoidal
        source = A.*sin(2*pi*f*dt*t);
      case 3     % modulated
        source = A.*sin(2*pi*f*dt*t).*exp(-0.5.*((centre - t)/width).^2);
  end  % switch

 % ELECTRICAL PROPERTIES OF MEDIA
   eR = ones(1,Nz).* eR1;  % Relative permittivity
   
   eR(indexR) = eR2;        
   %eR(M2:end) = eR2;
   
  % eR(88:end) = 1;
   S = ones(1,Nz) .* S1;   % Conductivity 
   S(M2:Nz) = S2;    

% Medium 1
    f1 = f;              % EM wave frequency   [Hz]
    c1 = c0/sqrt(eR1);   % speed of light   [m/s]      
    T1 = 1/f1;           % period   [s]
    w1 = 2*pi*f1;        % angular frequency   [1/s]
    wL1 = c1/f1;         % wavelength   [m]
    k1 = 2*pi/wL1;       % angular wave number   [1/m]

% Medium 2 dielectric
    f2 = f1;             % EM wave frequency   [Hz]
    c2 = c0/sqrt(eR2);    % speed of light   [m/s]
    T2 = 1/f2;            % period   [s]
    w2 = 2*pi*f2;         % angular frequency   [1/s]
    wL2 = c2/f2;          % wavelength   [m]
    k2 = 2*pi/wL2;        % angular wave number   [1/m] 
   
% K constants 
  K1 = 1 - (dt .* S) ./ (2 .* eR .* eps0);
  K2 = 1 + (dt .* S) ./ (2 .* eR .* eps0);
  K3 = K1 ./K2;
  K4 = dt * c0 / dz;
  K5 = K4 ./ (eR .* K2);

% Initialize fields 
  E = zeros(Nt,Nz);
  H = zeros(Nt,Nz);
  E(1,zS) = source(1);

% FDTD CALCULATIONS ------------------------------------------------------
for ct = 2 : Nt
      
for cz = 2 : Nz-1
      E(ct,cz) = K3(cz)*E(ct-1,cz) - K5(cz) * (H(ct-1,cz) - H(ct-1,cz-1));
     % if ct < 100
     E(ct,zS) = source(ct);
     %E(ct,Nz-2) = -source(ct);
     % end
end   % cz 

% E absorbing B.C.
   if flagBC == 1
     if ct > 2
         E(ct,1) = E(ct-2,2); E(ct,Nz) = E(ct-2,Nz-1);
     end
   end
% E PEC B.C.
   if flagBC == 2
     E(ct,1) = 0; E(ct,Nz) = 0;
   end 

   for cz = 1 : Nz-1
     H(ct,cz) = H(ct-1,cz) - K4*(E(ct,cz+1) - E(ct,cz));
   end   % cz

% H absorbing B.C.
if flagBC == 1
   if ct > 2
       H(ct,1) = H(ct-2,2); H(ct,Nz) = H(ct-2,Nz-1);
  end
end

% H PMC B.C.
   if flagBC == 3
     H(ct,1) = 0; H(ct,Nz) = 0;
   end 

end   % ct

% END FDTD CALCULATIONS =================================================

% % Probe measurements 
%        flagP = zeros(5,1); th = 0.2; v12 = 0; v45 = 0;
%        ctA1(1) = 0; ctA2(1) = 0; ctA3(1) = 0; ctA4(1) = 0; ctA5(1)=0;
%     if flagS < 3
%        if max(abs(E(:,cP(1)))) > th
%         [ctA1, ct1] = findpeaks(abs(E(:,cP(1))),'MinPeakProminence',th);
%         flagP(1) = 1;
%        end
%        
%        if max(abs(E(:,cP(2)))) > th
%          [ctA2, ct2] = findpeaks(abs(E(:,cP(2))),'MinPeakProminence',th);
%          flagP(2) = 1;
%        end
%        
%        if max(abs(E(:,cP(3)))) > th
%         [ctA3, ct3] = findpeaks(abs(E(:,cP(3))),'MinPeakProminence',th);
%         flagP(3) = 1;
%        end
%        
%        if max(abs(E(:,cP(4)))) > th
%         [ctA4, ct4] = findpeaks(abs(E(:,cP(4))),'MinPeakProminence',th);
%          flagP(4) = 1;
%        end
%        
%        if max(abs(E(:,cP(5)))) > th
%          [ctA5, ct5] = findpeaks(abs(E(:,cP(5))),'MinPeakProminence',th);
%          flagP(5) = 1;
%        end
%        
%        if flagP(1)*flagP(2)  > 0 
%           deltaZ = (cP(2)-cP(1))*dz;
%           deltaT = (ct2(1) - ct1(1))*dt;
%           v12 = deltaZ / deltaT;
%        end
%        
%       if flagP(4)*flagP(5) > 0 
%          deltaZ = (cP(5)-cP(4))*dz;
%          deltaT = (ct5(1) - ct4(1))*dt;
%          v45 = deltaZ / deltaT;
%       end
%     end

% Sinusoidal source
   lambda1 = c0/(sqrt(eR1)*f);              % wavelenth medium 1
   lambda2 = c0/(sqrt(eR2)*f);  % wavelength medium 2
   
% Lossy Medium  
if S2 > 0
% Decay constant  
    alpha = (w2/c0) * sqrt(eR2/2) * sqrt( sqrt(1+(S2/(w2*eps0*eR2))^2) - 1 );
% angular wave number [1/m] and wavelength [m] in Medium 2 
    k = (w2/c0) * sqrt(eR2/2) * sqrt( sqrt(1+(S2/(w2*eps0*eR2))^2) + 1 );
    lambdaN = 2*pi/k;
    
% Find peaks in transmitted wave a: peak values / b: indices for peaks
% Input range for indices to findpeaks in Medium 2
    index2 = 200:260; 
    [a, b] = findpeaks(E(end,index2), z(index2),'MinPeakProminence',0.05);
% Decay constant from exponetial decay of EM wave in Medium 2  [m]  
    alphaN = -log(a(1)/a(end))/(b(1)-b(end))/dz;
end


% =======================================================================
% GRAPHICS
% =======================================================================

figure(1)   % Animation of E(t,z)and H(t,z) -----------------------------
   pos = [0.02 0.07 0.30 0.62];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   step = 10 ;
   for ct = 10 : step : Nt
   subplot(2,1,1) 
     xP = z; yP = E(ct,:);
     plot(xP,yP,'b','linewidth',2)
     hold on
     xxP = cP; yyP = zeros(length(xxP),1);
     hPlot = plot(xxP,yyP,'mo');
     set(hPlot,'markerfacecolor','m')
     xxP = [indexR(1) indexR(1)]; yyP = [-limE limE];
     plot(xxP,yyP,'m','linewidth',1);
     xxP = [indexR(end) indexR(end)]; yyP = [-limE limE];
     plot(xxP,yyP,'m','linewidth',1);
     grid on
     hold off
     set(gca,'yLim',[-limE limE])
  %  set(gca,'xLim',[0 5e-6])
     xlabel('z','fontsize',14)
     ylabel(' E_{sx} ','fontsize',14)
     tm1 ='width = ';
     tm2 = num2str(indexR(end)-indexR(1),'%3.0f  \n');
     tm3 = '   time step = ';
     tm4 = num2str(ct,'%4.0f  \n');
     tm = [tm1 tm2 tm3 tm4];
     title(tm)
     set(gca,'fontsize',12)
     
   subplot(2,1,2) 
     xP = z; yP = H(ct,:);
     plot(xP,yP,'r','linewidth',2)
     grid on
     xlabel('z  ','fontsize',14)
     ylabel(' H_{y} ','fontsize',14)
     set(gca,'yLim',[-limH limH])
     set(gca,'fontsize',12)
     % set(gca,'xLim',[0 5e-6])
     pause(0.1)
     
     if f_gif > 0 
       frame = getframe(1);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
     %  On the first loop, create the file. In subsequent loops, append.
       if nF == 1
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
       end
       nF = nF+1;
    end
     
   end

   
   
%  Probe measurements ----------------------------------------------------
figure(2)
  pos = [0.34 0.07 0.30 0.62];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w')
  LW = 2;
  
  subplot(3,1,1)
  xP = t; yP = E(:,cP(1));
  plot(xP,yP,'b','linewidth',LW)
  grid on
  ylabel('#1   E_{sx} ','fontsize',14)
  tm1 ='Width = ';
  tm2 = num2str(indexR(end)-indexR(1),'%3.0f  \n');
  tm3 = '     Probe Measurements';
  tm = [tm1 tm2 tm3];
  title(tm,'fontweight','normal') 
  
  subplot(3,1,2)
  yP = E(:,cP(2));
  plot(xP,yP,'b','linewidth',LW)
  grid on
  ylabel('#2   E_{sx} ','fontsize',14)
    
  subplot(3,1,3)
  yP = E(:,cP(3));
  plot(xP,yP,'b','linewidth',LW)
%   yP = E(:,cP(4));
%   plot(xP,yP,'k','linewidth',LW)
%   yP = E(:,cP(5));
%   plot(xP,yP,'c','linewidth',LW)
  grid on
%  legend('#1','#2','#3','#4','#5','Orientation','horizontal', ...
%     'location','south');
  xlabel('t ','fontsize',14)
  ylabel('#3   E_{sx} ','fontsize',14)
  
%  set(gca,'yLim',[-limE limE])
  set(gca,'fontsize',12)

  
% figure(3)
%    pos = [0.34 0.47 0.30 0.42];
%    set(gcf,'Units','normalized');
%    set(gcf,'Position',pos);
%    set(gcf,'color','w');
%    set(gca,'xLim',[0 100]);
%    set(gca,'yLim',[-50 100]);
%    
%    yd = 99;
%    d = 12;
%    
%    text(0,yd,'time step  dt =','fontsize',12,'color','k')
%    tm = num2str(dt,'%3.2e s  \n');
%    text(25,yd,tm,'fontsize',12,'color','k')
%    
%    text(50,yd,'cell size dz =','fontsize',12,'color','k')
%    tm = num2str(dz,'%3.2e m  \n');
%    text(72,99,tm,'fontsize',12,'color','k');
%    
%    yd = yd - d;
%    text(0,yd,'PROBE positions [grid points]','fontsize',12,'color','k')
%   
%    yd = yd-d;
%    tm1 = 'z_1 = ';
%    tm2 = num2str((cP(1)),'%4.0f  \n');
%    tm3 = '   z_2 = ';
%    tm4 = num2str((cP(2)),'%4.0f  \n');
%    tm5 = '   z_3 = ';
%    tm6 = num2str((cP(3)),'%4.0f  \n');
%    tm7 = '   z_4 = ';
%    tm8 = num2str((cP(4)),'%4.0f  \n');
%    tm9 = '   z_5 = ';
%    tm10 = num2str((cP(5)),'%4.0f  \n');
%    tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7 tm8 tm9 tm10];
%    text(0,yd,tm,'fontsize',12,'color','k')
%    
%    yd = yd - d;
%    tm1 = 'Velocity between Probes #1 and #2';
%    tm2 = '   v_{12} = ';
%    tm3 = num2str(v12,'%4.2e  m/s \n');
%    tm = [tm1 tm2 tm3];
%    text(0,yd,tm,'fontsize',12,'color','k')
%    
%    yd = yd - d;
%    tm1 = 'Velocity between Probes #4 and #5';
%    tm2 = '   v_{45} = ';
%    tm3 = num2str(v45,'%4.2e  m/s \n');
%    tm = [tm1 tm2 tm3];
%    text(0,yd,tm,'fontsize',12,'color','k')
%   
%   if flagS == 1
%    yd = yd - d;
%    tm = 'Amplitudes';
%    text(0,yd,tm,'fontsize',12,'color','k')
%   
%    yd = yd - d;
%    tm1 = '   A_1 = ';
%    tm2 = num2str(ctA1(1),'%2.3f   \n');
%    tm3 = '   A_2 = ';
%    tm4 = num2str(ctA2(1),'%2.3f   \n');
%    tm5 = '   A_3 = ';
%    tm6 = num2str(ctA3(1),'%2.3f   \n');
%    tm7 = '   A_4 = ';
%    tm8 = num2str(ctA4(1),'%2.3f   \n');
%    tm9 = '   A_5 = ';
%    tm10 = num2str(ctA5(1),'%2.3f   \n');
%    
%    tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7 tm8 tm9 tm10];
%    text(0,yd,tm,'fontsize',12,'color','k')
%   end
%  
%   if flagS == 2
%    yd = yd - d;
%    tm = 'Wavelengths';
%    text(0,yd,tm,'fontsize',12,'color','k')
%    
%    yd = yd - d;
%    tm1 = '   \lambda_1 = ';
%    tm2 = num2str(lambda1,'%2.3e   m   \n');
%    tm3 = '        \lambda_2 = ';
%    tm4 = num2str(lambda2,'%2.3e   m   \n');
%    tm = [tm1 tm2 tm3 tm4];
%    text(0,yd,tm,'fontsize',12,'color','k')
%   end
%   
%   
%  if S2 > 0
%    yd = yd - 2*d;
%    tm1 = 'Lossy Medium 2:';
%    tm2 = '   conductivity \sigma_2 = ';
%    tm3 = num2str(S2,'%2.3e   \n');
%    tm = [tm1 tm2 tm3];
%    text(0,yd,tm,'fontsize',12,'color','k')
%    
%    yd = yd - d;
%    tm = 'Decay constant  Theory (\alpha) and Numerical (\alpha_N)';
%    text(0,yd,tm,'fontsize',12,'color','k') 
%    
%    yd = yd - d;
%    tm1 = '\alpha = ';
%    tm2 = num2str(alpha,'%2.3e  m   \n');
%    tm3 = '     \alpha_N = ';
%    tm4 = num2str(alphaN,'%2.3e  m   \n');
%    tm = [tm1 tm2 tm3 tm4];
%    text(0,yd,tm,'fontsize',12,'color','k')
%    
%    yd = yd - d;
%    tm1 = 'wavelength (numerical prediction)   ';
%    tm2 = '\lambda_N = ';
%    tm3 = num2str(lambdaN,'%2.3e  m   \n');
%    tm = [tm1 tm2 tm3];
%    text(0,yd,tm,'fontsize',12,'color','k')
%  end
% 
%    axis off
   

figure(4)    % source  Esource(t)  --------------------------------------
   pos = [0.65 0.07 0.30 0.30];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');

   xP = t; yP = source;
   plot(xP,yP,'b','linewidth',2)
   grid on
   xlabel('t  ','fontsize',14)
   ylabel(' E_{source} ','fontsize',14)
   set(gca,'yLim',[-1.1 1.1])
   set(gca,'fontsize',12)

  
 toc


