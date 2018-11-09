%% direct convolution
% use convn instead of Fourier Transform

%% impulse response (propagation)

smoothed = true;
Nf = 50;

Nv = [1600 3200 4800];
for k=1:3
    N = Nv(k);
    D = 15;
    x = linspace(-1,1,N+1)*D;
    dx = x(2)-x(1);

    [cc2,ss2] = fresnel(sqrt(2*Nf)*(x+dx/2));
    [cc1,ss1] = fresnel(sqrt(2*Nf)*(x-dx/2));


    % this version is smoothed
    if smoothed
        f = (1/sqrt(2*Nf))*complex(cc2-cc1,ss2-ss1)/dx;
    else
        % qchirp is exact
        f = qchirp(x,Nf);
    end

    xrange = D;
    subplot(3,1,k);
    %plot(x,z,'k','LineWidth',2);
    plot(x,imag(f),'b',x,real(f),'k','LineWidth',2);
    %axis([-xrange xrange -1.5 1.5]);
    xlabel('x');
    ylabel('f(x)');
end

%% prop_response

Nf = 50;
smoothed = true;


N=3200;
D=4;

x=linspace(-1,1,N+1)*D;
dx = x(2)-x(1);

[cc2,ss2] = fresnel(sqrt(2*Nf)*(x+dx/2));
[cc1,ss1] = fresnel(sqrt(2*Nf)*(x-dx/2));


if smoothed
    f = (1/sqrt(2*Nf))*complex(cc2-cc1,ss2-ss1)/dx;
else
    f = qchirp(x,Nf);
end


xrange = D;
subplot(2,1,1);
%plot(x,z,'k','LineWidth',2);
plot(x,imag(f),'b',x,real(f),'k','LineWidth',2);
%axis([-xrange xrange -1.5 1.5]);
xlabel('x');
ylabel('f(x)');

g = rect(x/2);
z = convn(g,f,'same')*(sqrt(Nf)*dx);
y = z.*conj(z);

if (Nf<1)
    y = y/4;
end

subplot(2,1,2);
xrange = 2;
idx = find(abs(x)<xrange);
plot(x(idx),y(idx),'k'); %,'LineWidth',2);


%% compare to exact

subplot(2,1,1);
slit(Nf);
ylabel('irradiance');

%% Nf = 100

Nf = 100;
prop_response(Nf);

%%

subplot(2,1,1);
slit(Nf);
ylabel('irradiance');

%% Nf = 1000;

Nf = 1000;
prop_response(Nf);

%%

direct_compare(Nf);

%% Nf = 10;

Nf = 10;
prop_response(Nf);

%%

direct_compare(Nf);

%% Nf = 0.9

Nf = 0.9

% body of pcompare
% near field calculation
[fy fx] = prop_response(Nf);
% exact calculation
[y x] = slit(Nf);
idx = find(abs(fx)<2.0);
subplot(1,1,1);
plot(fx(idx),fy(idx),'b',x,y,'k','LineWidth',2);
if (Nf<1)
    xlabel('x/b');
    ylabel('normalized irradiance');
else
    xlabel('x/w');
    ylabel('irradiance');
end