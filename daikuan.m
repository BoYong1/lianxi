
clc
clear
e=1.6e-19;
me=9.11e-31;
mi=1800*me;
B=0:0.00001:1;
w=0:1e4:5.269e9;
wp=2.84e9;
wle=e.*B/me/(2*pi);
wli=e.*B/mi/(2*pi);

%NR=1-wp^2./((w-wle).*(w+wli));
%NR=1-(wp./w).^2./(1-wle./w);

figure(2)
plot(B,wle)
