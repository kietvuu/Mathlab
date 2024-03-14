%%
clear
close all
Rk = 0.7;
Rt =3.6;
E = 10;
x=0:0.001:15;
f= (E*x)./(x+Rk)-(E*x)./(x+Rt);
[fmax,indmax]=max(f); %largest value of vector f is its element nu
% mber indmin 
xmax=x(indmax); 

plot(x,f,'linewidth',1.5)
hold on
plot(xmax,fmax,'r.','markersize',20);
hold off
grid
xlabel('x')
ylabel('f(x)')
%%
clear
close all
kv = 50;
ke = 100;
r = 10000;

%minimum at x=x0
x0=sqrt((2*r*ke)/kv);

%minimum value
L0= kv * x0./2 + ke * r./x0;

%graph
x=0:x0/100:800*x0;
L=kv * x./2 + ke * r./x;

plot(x,L,'linewidth',1.5)
hold
plot(x0,L0,'r.','markersize',20)
hold off
grid
ylim([0,2*L0])
xlabel('x')
ylabel('L(x)')

%%
clear
close all
A=1;
% h = (A - pi*r.^2)./(2*pi*r)
% V = (r.*A-r.^3*pi)./2

r0=(A/(3*pi))^(1/2);
V0=(r0.*A-r0.^3*pi)./2;

r=0:r0/100:1.7*r0;
V = (r.*A-r.^3*pi)./2;

plot(r,V,'linewidth',1.5)
hold on
plot(r0,V0,'r.','markersize',20)
hold off
grid

%%
clear
close all
a = 5;
b = 9;
%MA = MB = b^2/4 + h.^2
%(MC + h).^2 = a^2- b^2/4
%s = -h + 2*sqrt(b^2/4+h.^2) + sqrt(a^2-b^2/4)
if 1> b/(2*a) > (sqrt(3)/2)
    h0 = sqrt(a^2-b^2/4);
    s0= -h0 + 2*sqrt(b^2/4+h0.^2) + sqrt(a^2-b^2/4);
    h=0:h0/100:h0;
    s =-h + 2*sqrt(b^2/4 + h.^2) + sqrt(a^2-b^2/4);
    plot(h,s,'linewidth',1.5)
    hold on
    plot(h0,s0,'r.','markersize',20)
    hold off
    grid
    title('case 2')
else
    h0 = b/sqrt(12);
    s0= -h0 + 2*sqrt(b^2/4+h0.^2) + sqrt(a^2-b^2/4);
    h=0:h0/100:4;
    s =-h + 2*sqrt(b^2/4 + h.^2) + sqrt(a^2-b^2/4);
    plot(h,s,'linewidth',1.5)
    hold on
    plot(h0,s0,'r.','markersize',20)
    hold off
    grid
    title('case 1')
end

%%
clear
close all
L = 10;
h = 4;
v = 1.05;
%PQ=1*t1, QB=v*t2
%t1+t1= PQ+QB/v = sqrt(h^2+x.^2) + (L-x)./v (x = AQ)
x=0:0.001:10;
f=sqrt(h^2+x.^2) + (L-x)./v;
[fmin,indmin]=min(f);%smallest value of vector f is its element number indmin 
xmin=x(indmin); 
plot(x,f,'linewidth',1.5)
hold on
plot(xmin,fmin,'r.','markersize',20);
hold off
grid
xlabel('x')
ylabel('f(x)')
%%
clear
close all
m = 2;
R = 1;
k = 3;
M = 10;
b = 5;

if b^2-2*k*M>0
    w=0:0.01:10;
    C=(m*R*w.^2)./sqrt((k-M*w.^2).^2+(b*w).^2);
    [Cmax,indmax]=max(C); %largest value of vector f is its element nu
% mber indmin 
    wmax=w(indmax); 
    plot(w,C,'linewidth',1.5)
    hold on 
    plot(wmax,Cmax,'r.','markersize',20);
    hold off
    grid
    title('case 1')
else
    w=0:0.001:3;
    C=(m*R*w.^2)./sqrt((k-M*w.^2).^2+(b*w).^2);
    [Cmax,indmax]=max(C); %largest value of vector f is its element nu
% mber indmin 
    wmax=w(indmax); 
    plot(w,C,'linewidth',1.5)
    hold on 
    plot(wmax,Cmax,'r.','markersize',20);
    hold off
    grid
    title('case 2')
end
    
