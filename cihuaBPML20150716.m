%   无等离子体，纯电波BPML，反射较小
%   有等离子体时，能看出明显变化
clc;clear;close all;
c=3e8;
s0=8.85*10^(-12);
e=1.6*10^(-19);
me=9.11*10^(-31);
u0=4*pi*1e-7;
f0=1e9;
fv=2e10*0;
wp=19e9*0;
wc=7.04e10*0;
nx=50;      % 网格数
ny=50;
nz=50;
dx=c/f0/20;
dy=dx;
dz=dx;
dt=dx/c/2;
Jx(nx,ny,nz)=0;
Jy(nx,ny,nz)=0;
Jz(nx,ny,nz)=0;
Ex(nx,ny+1,nz+1)=0;
Ey(nx+1,ny,nz+1)=0;
Ez(nx+1,ny+1,nz)=0;
Hx(nx,ny,nz)=0;
Hy(nx,ny,nz)=0;
Hz(nx,ny,nz)=0;
%   BPML中参数初始化
sigmam=0;
sigmax(nx,ny,nz)=0;
sigmay(nx,ny,nz)=0;
sigmaz(nx,ny,nz)=0;
toux(nx,ny,nz)=0;
touy(nx,ny,nz)=0;
touz(nx,ny,nz)=0;
Ahx(nx,ny,nz)=0;
Ahy(nx,ny,nz)=0;
Ahz(nx,ny,nz)=0;
Bhx(nx,ny,nz)=0;
Bhy(nx,ny,nz)=0;
Bhz(nx,ny,nz)=0;
Aex(nx,ny,nz)=0;
Aey(nx,ny,nz)=0;
Aez(nx,ny,nz)=0;
Bex(nx,ny,nz)=0;
Bey(nx,ny,nz)=0;
Bez(nx,ny,nz)=0;
Exy(nx,ny,nz)=0;
Exz(nx,ny,nz)=0;
Eyx(nx,ny,nz)=0;
Eyz(nx,ny,nz)=0;
Ezx(nx,ny,nz)=0;
Ezy(nx,ny,nz)=0;
Hxy(nx,ny,nz)=0;
Hxz(nx,ny,nz)=0;
Hyx(nx,ny,nz)=0;
Hyz(nx,ny,nz)=0;
Hzx(nx,ny,nz)=0;
Hzy(nx,ny,nz)=0;

%    BPML参数设定
npml=8;
sigmam=(4+1)/(150*pi*dx);
for i=1:npml    
    sigmax(i,:,:)=sigmam*(npml+1-i)^4/npml^4; 
    sigmax(nx+1-i,:,:)=sigmam*(npml+1-i)^4/npml^4;
    sigmay(:,i,:)=sigmam*(npml+1-i)^4/npml^4; 
    sigmay(:,nx+1-i,:)=sigmam*(npml+1-i)^4/npml^4;
    sigmaz(:,:,i)=sigmam*(npml+1-i)^4/npml^4; 
    sigmaz(:,:,nx+1-i)=sigmam*(npml+1-i)^4/npml^4;
end
toux=u0/s0.*sigmax;
touy=u0/s0.*sigmay;
touz=u0/s0.*sigmaz;
Ahx=(2*u0-toux.*dt)./(2*u0+toux.*dt);
Ahy=(2*u0-touy.*dt)./(2*u0+touy.*dt);
Ahz=(2*u0-touz.*dt)./(2*u0+touz.*dt);
Bhx=2*dt./(2*u0+toux.*dt);
Bhy=2*dt./(2*u0+touy.*dt);
Bhz=2*dt./(2*u0+touz.*dt);
Aex=(2*s0-sigmax.*dt)./(2*s0+sigmax.*dt);
Aey=(2*s0-sigmay.*dt)./(2*s0+sigmay.*dt);
Aez=(2*s0-sigmaz.*dt)./(2*s0+sigmaz.*dt);
Bex=2*dt./(2*s0+sigmax.*dt);
Bey=2*dt./(2*s0+sigmay.*dt);
Bez=2*dt./(2*s0+sigmaz.*dt);
%  方便画图
I=1:1:nx;
J=1:1:ny;
K=1:1:nz;
Ek(nx,nz)=0;
%   等离子体参数初始化
npx1=round(0.5/6*nx);
npx2=round(5.5/6*nx);
npy1=round(0.5/6*ny);
npy2=round(5.5/6*ny);
npz1=round(2.5/5*nz);
npz2=round(3/5*nz);
%A(nx,ny,nz)=1;
%B(nx,ny,nz)=0;
%C(nx,ny,nz)=0;
%D(nx,ny,nz)=1;
%E(nx,ny,nz)=0;

A=(4-wc^2*dt^2)*exp(-fv*dt)/(4+wc^2*dt^2*exp(-fv*dt));
B=4*s0*wp^2*dt*exp(-fv*dt/2)/(4+wc^2*dt^2*exp(-fv*dt));
C=2*wc*dt*exp(-fv*dt/2)/(4+wc^2*dt^2*exp(-fv*dt));
D=1+exp(-fv*dt);
E=s0*wp^2*dt*exp(-fv*dt/2);
N=80;
for n=1:1:N 
    for i=npml+1:1:nx-npml
        for j=npml+1:1:ny-npml
            for k=npml+1:1:nz-npml
                Ex(i,j,k)=Ex(i,j,k)+dt/s0/dx*(Hz(i,j,k)-Hz(i,j-1,k)-Hy(i,j,k)+Hy(i,j,k-1))-dt/s0*Jx(i,j,k);    
                Ey(i,j,k)=Ey(i,j,k)+dt/s0/dx*(Hx(i,j,k)-Hx(i,j,k-1)-Hz(i,j,k)+Hz(i-1,j,k))-dt/s0*Jy(i,j,k);
                Ez(i,j,k)=Ez(i,j,k)+dt/s0/dx*(Hy(i,j,k)-Hy(i-1,j,k)-Hx(i,j,k)+Hx(i,j-1,k))-dt/s0*Jz(i,j,k);
            end
        end
    end
    for i=npx1:1:npx2;%npml+1:1:nx-npml
        for j=npy1:1:npy2;%npml+1:1:ny-npml
            for k=npz1:1:npz2;%npml+1:1:nz-npml
                Jx(i,j,k)=A*Jx(i,j,k)+B*Ex(i,j,k)-C*(D*Jy(i,j,k)+E*Ey(i,j,k));
                Jy(i,j,k)=A*Jy(i,j,k)+B*Ey(i,j,k)+C*(D*Jx(i,j,k)+E*Ex(i,j,k));
                Jz(i,j,k)=(D-1)*Jz(i,j,k)+E*Ez(i,j,k);
            end
        end
    end
    for i=npml+1:1:nx-npml
        for j=npml+1:1:ny-npml
            for k=npml+1:1:nz-npml
                Hx(i,j,k)=Hx(i,j,k)-dt/u0/dx*(Ez(i,j+1,k)-Ez(i,j,k)-Ey(i,j,k+1)+Ey(i,j,k));
                Hy(i,j,k)=Hy(i,j,k)+dt/u0/dx*(Ez(i+1,j,k)-Ez(i,j,k)-Ex(i,j,k+1)+Ex(i,j,k));
                Hz(i,j,k)=Hz(i,j,k)-dt/u0/dx*(Ey(i+1,j,k)-Ey(i,j,k)-Ex(i,j+1,k)+Ex(i,j,k));
            end
        end
    end
    %BPML边界条件 
    for i=2:nx-1     
        for j=2:ny-1   
            for k=2:nz-1   
                if (i<=npml)||(i>=nx-npml+1)||(j<=npml)||(j>=nx-npml+1)||(k<=npml)||(k>=nx-npml+1)                
                    Hxy(i,j,k)=Ahy(i,j,k)*Hxy(i,j,k)-Bhy(i,j,k)*(Ez(i,j+1,k)-Ez(i,j,k))/dy;
                    Hxz(i,j,k)=Ahz(i,j,k)*Hxz(i,j,k)-Bhz(i,j,k)*(-Ey(i,j,k+1)+Ey(i,j,k))/dz;
                    Hx(i,j,k)=Hxy(i,j,k)+Hxz(i,j,k);
                    Hyz(i,j,k)=Ahz(i,j,k)*Hyz(i,j,k)-Bhz(i,j,k)*(Ex(i,j,k+1)-Ex(i,j,k))/dz;
                    Hyx(i,j,k)=Ahx(i,j,k)*Hyx(i,j,k)-Bhx(i,j,k)*(-Ez(i+1,j,k)+Ez(i,j,k))/dx;
                    Hy(i,j,k)=Hyz(i,j,k)+Hyx(i,j,k);
                    Hzx(i,j,k)=Ahx(i,j,k)*Hzx(i,j,k)-Bhx(i,j,k)*(Ey(i+1,j,k)-Ey(i,j,k))/dx;
                    Hzy(i,j,k)=Ahy(i,j,k)*Hzy(i,j,k)-Bhy(i,j,k)*(-Ex(i,j+1,k)+Ex(i,j,k))/dy;
                    Hz(i,j,k)=Hzx(i,j,k)+Hzy(i,j,k);  
                    Exy(i,j,k)=Aey(i,j,k)*Exy(i,j,k)+Bey(i,j,k)*(Hz(i,j,k)-Hz(i,j-1,k))/dy;
                    Exz(i,j,k)=Aez(i,j,k)*Exz(i,j,k)+Bez(i,j,k)*(-Hy(i,j,k)+Hy(i,j,k-1))/dz;
                    Ex(i,j,k)=Exy(i,j,k)+Exz(i,j,k);
                    Eyz(i,j,k)=Aez(i,j,k)*Eyz(i,j,k)+Bez(i,j,k)*(Hx(i,j,k)-Hx(i,j,k-1))/dz;
                    Eyx(i,j,k)=Aex(i,j,k)*Eyx(i,j,k)+Bex(i,j,k)*(-Hz(i,j,k)+Hz(i-1,j,k))/dx;
                    Ey(i,j,k)=Eyz(i,j,k)+ Eyx(i,j,k);
                    Ezx(i,j,k)=Aex(i,j,k)*Ezx(i,j,k)+Bex(i,j,k)*(Hy(i,j,k)-Hy(i-1,j,k))/dx;
                    Ezy(i,j,k)=Aey(i,j,k)*Ezy(i,j,k)+Bey(i,j,k)*(-Hx(i,j,k)+Hx(i,j-1,k))/dy; 
                    Ez(i,j,k)=Ezx(i,j,k)+Ezy(i,j,k);
                end
            end 
        end
    end
    
    pause(0.0001);
    Ey(npml:nx-npml,npml:ny-npml,10)=sin(2*pi*f0*n*dt)+Ey(npml:nx-npml,npml:ny-npml,10);
    %figure(1) 
    for i=1:nx
        for j=1:nz
           Ek(i,j)=Ey(i,25,j);  
        end
    end
    figure(1)   
    contourf(I,J,Ek); %surf(Ez(:,:,25));%pcolor(Ez(:,:,25));
    colorbar;
    %caxis([-1 1]);
    figure(2)
    plot(Ek(25,:));
    axis([0 50 -1.5 1.5])
end







