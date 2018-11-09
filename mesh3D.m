function [z,Fx,Fy]=mesh3D(f)
    [M,N]=size(f);
    dr = 0.5; dc = 0.5;
    Fx = ((0:M-1)-floor(M/2))/(M*dc);   %对Fx平移,并归一化
    Fy = ((0:N-1)-floor(N/2))/(N*dr);   %对Fy平移,并归一化
    [Fx,Fy] = meshgrid(Fx, Fy);         %画网格
    z=double(f(:,:));
end