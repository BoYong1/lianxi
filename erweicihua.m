%   无等离子体，纯电波BPML，反射较小
%   有等离子体时，能看出明显变化
clc;clear;close all;
c=3e8;
s0=8.85*10^(-12);
e=1.6*10^(-19);
me=9.11*10^(-31);
u0=4*pi*1e-7;
f0=1e9;
%   网格信息
dx=c/f0/30;
dz=dx;
dt=dx/c/2;
x=1.5;                 %  x=1m
z=1.5;                 %  z=1m
nx=round(x/dx);      % 网格数
nz=round(z/dz);
Jx(nx,nz)=0;
Jy(nx,nz)=0;
Jz(nx,nz)=0;
Ex(nx,nz+1)=0;
Ey(nx+1,nz+1)=0;
Ez(nx+1,nz)=0;
Hx(nx,nz)=0;
Hy(nx,nz)=0;
Hz(nx,nz)=0;
%   BPML中参数初始化
sigmam=0;
sigmax(nx,nz)=0;
sigmay(nx,nz)=0;
sigmaz(nx,nz)=0;
toux(nx,nz)=0;
touy(nx,nz)=0;
touz(nx,nz)=0;
Ahx(nx,nz)=0;
Ahy(nx,nz)=0;
Ahz(nx,nz)=0;
Bhx(nx,nz)=0;
Bhy(nx,nz)=0;
Bhz(nx,nz)=0;
Aex(nx,nz)=0;
Aey(nx,nz)=0;
Aez(nx,nz)=0;
Bex(nx,nz)=0;
Bey(nx,nz)=0;
Bez(nx,nz)=0;
Exy(nx,nz)=0;
Exz(nx,nz)=0;
Eyx(nx,nz)=0;
Eyz(nx,nz)=0;
Ezx(nx,nz)=0;
Ezy(nx,nz)=0;
Hxy(nx,nz)=0;
Hxz(nx,nz)=0;
Hyx(nx,nz)=0;
Hyz(nx,nz)=0;
Hzx(nx,nz)=0;
Hzy(nx,nz)=0;

%    BPML参数设定
npml=4;
sigmam=(4+1)/(150*pi*dx);
for i=1:npml    
    sigmax(i,:)=sigmam*(npml+1-i)^4/npml^4; 
    sigmax(nx+1-i,:)=sigmam*(npml+1-i)^4/npml^4;
    sigmaz(:,i)=sigmam*(npml+1-i)^4/npml^4; 
    sigmaz(:,nx+1-i)=sigmam*(npml+1-i)^4/npml^4;
end
toux=u0/s0.*sigmax;
touz=u0/s0.*sigmaz;
Ahx=(2*u0-toux.*dt)./(2*u0+toux.*dt);
Ahy=1;%(2*u0-touy.*dt)./(2*u0+touy.*dt);
Ahz=(2*u0-touz.*dt)./(2*u0+touz.*dt);
Bhx=2*dt./(2*u0+toux.*dt);
Bhy=2*dt./(2*u0);
Bhz=2*dt./(2*u0+touz.*dt);
Aex=(2*s0-sigmax.*dt)./(2*s0+sigmax.*dt);
Aey=1;%(2*s0-sigmay.*dt)./(2*s0+sigmay.*dt);
Aez=(2*s0-sigmaz.*dt)./(2*s0+sigmaz.*dt);
Bex=2*dt./(2*s0+sigmax.*dt);
Bey=2*dt./(2*s0);
Bez=2*dt./(2*s0+sigmaz.*dt);
%  方便画图
I=1:1:nx;
K=1:1:nz;
Ek(nx,nz)=0;
%   等离子体参数初始化
B0=0.8;
fv=2e10;
wp=19e9;
wc=e*B0/me;
npx1=round(5/15*nx);
npx2=round(10/15*nx);
npz1=round(5/15*nz);
npz2=round(10/15*nz);
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
N=250;
for n=1:1:N 
    for i=npml+1:1:nx-npml
            for k=npml+1:1:nz-npml
                Ex(i,k)=Ex(i,k)+dt/s0/dx*(-Hy(i,k)+Hy(i,k-1))-dt/s0*Jx(i,k);    
                Ey(i,k)=Ey(i,k)+dt/s0/dx*(Hx(i,k)-Hx(i,k-1)-Hz(i,k)+Hz(i-1,k))-dt/s0*Jy(i,k);
                Ez(i,k)=Ez(i,k)+dt/s0/dx*(Hy(i,k)-Hy(i-1,k))-dt/s0*Jz(i,k);
            end
    end
    for i=npx1:1:npx2;%npml+1:1:nx-npml
            for k=npz1:1:npz2;%npml+1:1:nz-npml
                Jx(i,k)=A*Jx(i,k)+B*Ex(i,k)-C*(D*Jy(i,k)+E*Ey(i,k));
                Jy(i,k)=A*Jy(i,k)+B*Ey(i,k)+C*(D*Jx(i,k)+E*Ex(i,k));
                Jz(i,k)=(D-1)*Jz(i,k)+E*Ez(i,k);
            end
    end
    for i=npml+1:1:nx-npml
            for k=npml+1:1:nz-npml
                Hx(i,k)=Hx(i,k)-dt/u0/dx*(-Ey(i,k+1)+Ey(i,k));
                Hy(i,k)=Hy(i,k)+dt/u0/dx*(Ez(i+1,k)-Ez(i,k)-Ex(i,k+1)+Ex(i,k));
                Hz(i,k)=Hz(i,k)-dt/u0/dx*(Ey(i+1,k)-Ey(i,k));
            end
    end
    %BPML边界条件 
    for i=2:nx-1     
            for k=2:nz-1   
                if (i<=npml)||(i>=nx-npml+1)||(k<=npml)||(k>=nx-npml+1)                
                    Hxy(i,k)=Ahy*Hxy(i,k);
                    Hxz(i,k)=Ahz(i,k)*Hxz(i,k)-Bhz(i,k)*(-Ey(i,k+1)+Ey(i,k))/dz;
                    Hx(i,k)=Hxy(i,k)+Hxz(i,k);
                    Hyz(i,k)=Ahz(i,k)*Hyz(i,k)-Bhz(i,k)*(Ex(i,k+1)-Ex(i,k))/dz;
                    Hyx(i,k)=Ahx(i,k)*Hyx(i,k)-Bhx(i,k)*(-Ez(i+1,k)+Ez(i,k))/dx;
                    Hy(i,k)=Hyz(i,k)+Hyx(i,k);
                    Hzx(i,k)=Ahx(i,k)*Hzx(i,k)-Bhx(i,k)*(Ey(i+1,k)-Ey(i,k))/dx;
                    Hzy(i,k)=Ahy*Hzy(i,k);
                    Hz(i,k)=Hzx(i,k)+Hzy(i,k);  
                    Exy(i,k)=Aey*Exy(i,k);
                    Exz(i,k)=Aez(i,k)*Exz(i,k)+Bez(i,k)*(-Hy(i,k)+Hy(i,k-1))/dz;
                    Ex(i,k)=Exy(i,k)+Exz(i,k);
                    Eyz(i,k)=Aez(i,k)*Eyz(i,k)+Bez(i,k)*(Hx(i,k)-Hx(i,k-1))/dz;
                    Eyx(i,k)=Aex(i,k)*Eyx(i,k)+Bex(i,k)*(-Hz(i,k)+Hz(i-1,k))/dx;
                    Ey(i,k)=Eyz(i,k)+ Eyx(i,k);
                    Ezx(i,k)=Aex(i,k)*Ezx(i,k)+Bex(i,k)*(Hy(i,k)-Hy(i-1,k))/dx;
                    Ezy(i,k)=Aey*Ezy(i,k); 
                    Ez(i,k)=Ezx(i,k)+Ezy(i,k);
                end
            end 
    end
    
    pause(0.0001);
    Ey(:,10)=sin(2*pi*f0*n*dt)+Ey(:,10);%round(nx/2)
    %figure(1) 
    for i=1:nx
        for j=1:nz
           Ek(i,j)=Ey(i,j);  
        end
    end
    figure(1)   
    contourf(I,K,Ek); %surf(Ez(:,:,25));%pcolor(Ez(:,:,25));
    colorbar;
    %caxis([-1 1]);
    %figure(2)
    %plot(Ek(round(nx/2),:));
    %axis([0 nx -1.5 1.5])
end







