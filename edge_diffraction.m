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

plot(x,f,'k','LineWidth',1.5);
%axis([-xrange xrange 0 1.5]);
xlabel('x');
%ylabel('E(x)');

%% edge image

x = linspace(-2,8,512);
edge = edge_fcn(x);
edge = edge/max(max(edge));
edge = ones(size(x))'*edge;
imshow(edge);
figure;
%imwrite(edge,'edge.jpg');


%% edge wavefront


a = 1;
N=400;

x = linspace(-1,5,N+1);

[cc,ss] = fresnel(sqrt(2/a)*x);

z = sqrt(1/(2*j))*complex(0.5+cc,0.5+ss);


g = z.*conj(z);

f =  - unwrap(angle(z))/(2*pi) - 1;

subplot(2,1,1);
plot(x,g,'k','LineWidth',1.5);
ylabel('irradiance');

subplot(2,1,2);
plot(x,f,'k','LineWidth',1.5);
%axis([-xrange xrange 0 1.5]);
xlabel('x');
ylabel('optical phase (waves)');
%ylabel('E(x)');


