%% exact slit diffraction

%% far field (Nf < 1);

Nf = [1e-8 0.05 0.25];
for k=1:3
    subplot(3,1,k);
    slit(Nf(k));
    if (k==2)
        ylabel('relative irradiance');
    end
end
xlabel('x/b');

%% near field (Nf > 1);

Nf = [1 8 16];
for k=1:3
    subplot(3,1,k);
    slit(Nf(k));
    if (k==2)
        ylabel('relative irradiance');
    end
end
xlabel('x/w');

%% very near field (Nf >> 1);

Nf = [10 100 1000];
for k=1:3
    subplot(3,1,k);
    slit(Nf(k));
    if (k==2)
        ylabel('relative irradiance');
    end
end
xlabel('x/w');