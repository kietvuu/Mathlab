% kiet vu
%contain:
%derivative 1 2 3 4
%spline_and_bezier_curves 1 2 3 4
%position_velocity_acceleration 1 2 4
%extremum_problems 123456
%integral 1  4
%taylor_polynomial 1

%% derivative
%ex1
clear
close all
%inline-function
f=@(x) 0.3*x.^3-0.5*x.^2-3*x+2; % f(x)
df=@(x) 0.9*x.^2-x-3; %f'(x)

%values at x0
x0= 4.5;
f0=f(x0); 
df0=df(x0);
x2 = 1.6;
f2=f(x2); 
df2=df(x2);
%graphs, when x=xmin...xmax
xmin=0;
xmax=5;
x=xmin:(xmax-xmin)/100:xmax;

figure(1)
plot(x,f(x),'b','linewidth',1.5)
hold on
L=1; %width of the tangent = 2L
x1 = x0 - f0/df0;
plot([x1,x0+L],[0,f0+df0*L],'k','linewidth',1)
plot(x0,f0,'r.','markersize',20)
plot(x1,0,'g.','markersize',20)

hold off
grid
ylabel('function f(x)')
title([' x_0 = ',num2str(x0),', f(x_0) = ',num2str(f0)])
axis([0,5,-4,12])

figure(2)
plot(x,f(x),'b','linewidth',1.5)
hold on
L=1; %width of the tangent = 2L
x1 = x0 - f0/df0;
plot([x1,x0+L],[0,f0+df0*L],'r','linewidth',1)
plot(x0,f0,'r.','markersize',20)
plot(x1,0,'g.','markersize',20)

hold off
grid
ylabel('function f(x)')
title([' x_0 = ',num2str(x0),', f(x_0) = ',num2str(f0)])
axis([0,5,-4,12])

figure(2)
plot(x,f(x),'b','linewidth',1.5)
hold on
L=1; %width of the tangent = 2L
x3 = x2 - f2/df2;
plot([x3,x2+L],[0,f2+df2*L],'r','linewidth',1)
plot(x2,f2,'r.','markersize',20)
plot(x3,0,'g.','markersize',20)

hold off
grid
ylabel('function f(x)')
title([' x_0 = ',num2str(x0),', f(x_0) = ',num2str(f0)])
axis([0,5,-4,12])


%%
%ex2
clear
close all
L = 5;
H = 3;
B = -20;

x=0:L/100:5;


%syms L H B a b
%S = solve(a*(L.^2)+b*L == H, 2*a*L+b == tand(B),a,b)
%a=S.a
%b=S.b

a = -(H - L*tand(B))/L^2; 
b = (2*H - L*tand(B))/L;
alpha = atand(b);
%a = -9,81/[2(v*cosd(alpha)).^2
%syms a b L H beta g v0 alpha
%S2= solve(a==-g/(2*(v0*cosd(alpha))^2),v0)
v = (2^(1/2)*(-9.81/a)^(1/2))/(2*cosd(alpha));

f=@(x)a.*(x.^2)+b.*x; % f(x)
df=@(x) 2*a.*x+b; %f'(x)
S=1;
plot(x,f(x),'b','linewidth',2)
hold on
plot([L,S+L],[H,f(L)+df(L)*S],'r','linewidth',1.5)
hold on
plot(L,H, 'b.', 'markersize',20)
hold off
grid

%%
%ex3 

clear
close all
x1=1;
y1=2;
k1=-0.25;
m1=1;
x2=5;
y2=4;
k2=0.5;
m2=-2;
%coefficient matrix
A=[x1^5,x1^4,x1^3,x1^2,x1,1
   5*x1^4,4*x1^3,3*x1^2,2*x1,1,0
   20*x1^3,12*x1^2,6*x1,2,0,0

  x2^5,x2^4,x2^3,x2^2,x2,1
   5*x2^4,4*x2^3,3*x2^2,2*x2,1,0
   20*x2^3,12*x2^2,6*x2,2,0,0];

%right hand side
B=[y1
   k1
   m1
   y2
   k2
   m2];

%solution i.e the coefficients a,b,c,d
sol=A\B;

a=sol(1);
b=sol(2);
c=sol(3);
d=sol(4);
e=sol(5);
f=sol(6);

f=@(x) a*x.^5 + b*x.^4 + c*x.^3 + d*x.^2 + e*x + f; % f(x)
df=@(x) 5*a*x.^4 + 4*b*x.^3 + 3*c*x.^2 + 2*d*x + e; %f'(x)
d2f=@(x) 20*a*x.^3 + 12*b*x.^2 + 6*c*x + 2*d; %f''(x)

%curvature
k=@(x) d2f(x)./(sqrt(1+df(x).^2)).^3;

%values
x1=1;
y1=f(x1);
df1=df(x1);
d2f1=d2f(x1);
k0=k(x1);

R1=abs(1/k0);

%polar angle of the tangent i.e vector [1,f'(x0)]
th1=atan2d(df1,1);
%center of the circle of curvature
if k0>=0
    kx0=x1+R1*cosd(th1+90);
    ky0=y1+R1*sind(th1+90);
else
    kx0=x1+R1*cosd(th1-90);
    ky0=y1+R1*sind(th1-90);
end
%points on the circle of curvature
t=0:360;
kx1=kx0+R1*cosd(t);
ky1=ky0+R1*sind(t);

%values
x2=5;
y2=f(x2);
df2=df(x2);
d2f2=d2f(x2);
k3=k(x2);

R2=abs(1/k3);

%polar angle of the tangent i.e vector [1,f'(x0)]
th2=atan2d(df2,1);
%center of the circle of curvature
if k3>=0
    kx3=x2+R2*cosd(th2+90);
    ky3=y2+R2*sind(th2+90);
else
    kx3=x2+R2*cosd(th2-90);
    ky3=y2+R2*sind(th2-90);
end
%points on the circle of curvature
t=0:360;
kx2=kx3+R2*cosd(t);
ky2=ky3+R2*sind(t);

%graphs, when x=xmin...xmax

dx=(x2-x1)/100;
x=x1:dx:x2;

figure(1)
plot(x1,y1,'r.','markersize',15)
hold on
plot(x2,y2,'g.','markersize',15)
plot(x,f(x),'b','linewidth',1)
plot(kx1,ky1,'r','linewidth',1)
L=1; %horizontal width =2L
plot([x1,x1+L],[y1,y1+df1*L],'r','linewidth',1)
plot([x1,kx0],[y1,ky0],'r','linewidth',1)
plot(kx0,ky0,'r.','markersize',10)

plot(kx2,ky2,'g','linewidth',1)
plot([x2,x2+L],[y2,y2+df2*L],'g','linewidth',1)
plot([x2,kx3],[y2,ky3],'g','linewidth',1)
plot(kx3,ky3,'g.','markersize',10)
hold off
grid
axis equal
xlabel('x')
ylabel('y')

title(['k1 = ',num2str(k1),', m1 = ',num2str(m1),...
       ', k2 = ',num2str(k2),', m2 = ',num2str(m2)])


%%
%ex4

clear
close all
A = 1;
a = 1.5;
Delta = 1.57;
B = 1.5;
b = 1;
T = 12.57;
t0 = 4.5;
t=0:T/100:T;

x=@(t) A*sin(a*t+Delta);
dx=@(t) A*a*cos(a*t+Delta);
d2x=@(t) -A*a^2*sin(a*t+Delta);
y=@(t) B*sin(b*t);
dy=@(t) B*b*cos(b*t);
d2y=@(t) -B*b^2*sin(b*t);
%length of tangent vector
c=@(t) sqrt(dx(t).^2+dy(t).^2);
%curvature
k=@(t) (dx(t).*d2y(t) - d2x(t).*dy(t))./(sqrt(dx(t).^2+dy(t).^2)).^3;

%values

x1=x(t0);
y1=y(t0);
df1=dy(t0);
df2=dx(t0);
d2f1=d2y(t0);
k0=k(t0);

R1=abs(1/k0);

%polar angle of the tangent i.e vector [1,f'(x0)]
th1=atan2d(df1,df2);
%center of the circle of curvature
if k0>=0
    kx0=x1+R1*cosd(th1+90);
    ky0=y1+R1*sind(th1+90);
else
    kx0=x1+R1*cosd(th1-90);
    ky0=y1+R1*sind(th1-90);
end
%points on the circle of curvature
m=0:360;
kx1=kx0+R1*cosd(m);
ky1=ky0+R1*sind(m);
L = 1;
figure(1)
plot(x(t),y(t),'b','linewidth',2)
hold
plot(x1,y1,'r.','markersize',15)
plot(kx1,ky1,'g','linewidth',1)
plot([x1,kx0],[y1,ky0],'r','linewidth',1)
plot([x1,x1+df2],[y1,y1+df1],'r','linewidth',2)
plot(kx0,ky0,'g.','markersize',10)
hold
axis equal
xlabel('x')
ylabel('y')
%
figure(2)
subplot(2,1,1)
plot(t,k(t),'b','linewidth',1.5)
grid
title('curvature \kappa(x)')
subplot(2,1,2)
plot(t,c(t),'b','linewidth',1.5)
grid
title('length of the tangent vector')
axis equal
xlabel('x')
ylabel('y')

%% spline_and_bezier_curves
%ex1
clear 
close all 
x=[1,4,3,1,0];
y=[1,2,4,5,3];
n=length(x);
t= 1:1:n;
%tapaus=case in finnish
tapaus=1; %1=natural, 2=clamped, 3=periodic

%tapaus=2, clamped spline
k1=1; %Y'(x1)
kn=2; %Y'(xn)
m1=-2;
mn=-1;

%if tapaus = 3, then requirement y(1)=y(n)
  

%coefficients from a function spline_curves_coefficients.m
%which has to be on the current folder 
%see spline_curves_coefficients.pdf
abcdx= spline_curves_coefficients(t,x,tapaus,k1,kn);
abcdy= spline_curves_coefficients(t,y,tapaus,m1,mn);
%Y_k=a_k*x^3+b_k*x^2+c_k*x+d_k
%coefficients a_k,b_k,c_k,d_k
%as rows in matrix abcd 
ax=abcdx(:,1);
bx=abcdx(:,2);
cx=abcdx(:,3);
dx=abcdx(:,4);

ay=abcdy(:,1);
by=abcdy(:,2);
cy=abcdy(:,3);
dy=abcdy(:,4);
% graph
figure(1)
plot(x,y,'r.','markersize',18)
hold on 
for k=1:n-1
  dt=(t(k+1)-t(k))/100;  
  tt=t(k):dt:t(k+1);
  xx=ax(k)*(tt-t(k)).^3+bx(k)*(tt-t(k)).^2+cx(k)*(tt-t(k))+dx(k);
  yy=ay(k)*(tt-t(k)).^3+by(k)*(tt-t(k)).^2+cy(k)*(tt-t(k))+dy(k);
  plot(xx,yy,'linewidth',1.5)
end 
hold off
grid
axis equal
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
if tapaus==1
  title('natural spline','fontsize',10)
elseif tapaus==2
  title('clamped spline','fontsize',10)
else
  title('periodic spline','fontsize',10)
end

%% 
%ex2
clear 
close all

%syms x2 y2 x1 y1 r s a b ca sa cb sb
%S3 = solve(x1+ r*ca == x4+ s*cb,y1+ r*sa == y4+ s*sb,r, s)
%S3.r
%S3.s
a= 80;
b = -50;
r = -(2* cosd(b) - 5*sind(b) - cosd(b)*1 + sind(b))/(cosd(a)*sind(b) - cosd(b)*sind(a));
x2 = 1 + r*cosd(a);
y2 = 1 + r*sind(a);
x3 = x2;
y3 = y2;

%control points

P1=[1,1];
P2=[x2,y2];
P3=[x3,y3];
P4=[5,2];

x1=P1(1);
y1=P1(2);
x2=P2(1);
y2=P2(2);
x3=P3(1);
y3=P3(2);
x4=P4(1);
y4=P4(2);




x=@(t) (1-t).^3*x1+3*(1-t).^2.*t*x2+3*(1-t).*t.^2*x3+t.^3*x4;
y=@(t) (1-t).^3*y1+3*(1-t).^2.*t*y2+3*(1-t).*t.^2*y3+t.^3*y4;

t=0:0.01:1;

figure(1)
plot([x1,x2,x3,x4],[y1,y2,y3,y4],'k.-','linewidth',1.5,'markersize',20)
hold 
plot(x(t),y(t),'b','linewidth',1.5)
hold off
grid
axis equal
xlabel('x','fontsize',14)
ylabel('y  ','rotation',0,'fontsize',14)

dx=-3*(1-t).^2*x1+...
   3*(1-4*t+3*t.^2)*x2+...
   3*(2*t-3*t.^2)*x3+...
   3*t.^2*x4;
d2x=6*(1-t)*x1+...
   3*(-4+6*t)*x2+...
   3*(2-6*t)*x3+...
   6*t*x4;
dy=-3*(1-t).^2*y1+...
   3*(1-4*t+3*t.^2)*y2+...
   3*(2*t-3*t.^2)*y3+...
   3*t.^2*y4;
d2y=6*(1-t)*y1+...
   3*(-4+6*t)*y2+...
   3*(2-6*t)*y3+...
   6*t*y4;

k=(dx.*d2y-d2x.*dy)./(dx.^2+dy.^2).^(3/2);
o=atan2d(dy, dx);
figure(2)
subplot(2,1,1)
plot(t,k,'linewidth',1.5)
grid
title('kaarevuus K(t)')
xlabel('t')
ylabel('curvature \kappa(t)')
subplot(2,1,2)
plot(t,o,'linewidth',1.5)
grid
title('suuntakulma 0(t)')
xlabel('t')
ylabel('y')
%% 
%ex3
clear 
close all

a= 2.54;
b = 5.2;
%control points
r = 4;
x0 = 6;
y0 = 5;
delta = -30;
beta = delta + 90;

x4 = x0 + r* cosd(delta);
y4 = y0 + r* sind(delta);
c = y4/sind(beta);
b = (3*c^2)/(2*r*sind(beta));
x3 = x4 - y4/tand(beta);
a = x3 - b; 
x2 = x3 - b;
x1 = x2 -a;

P1=[x1,0];
P2=[x2,0];
P3=[x3,0];
P4=[x4,y4];

x1=P1(1);
y1=P1(2);
x2=P2(1);
y2=P2(2);
x3=P3(1);
y3=P3(2);
x4=P4(1);
y4=P4(2);

x5 = x4;
y5 = y4;
delta = 0:360;
x5 = x0 + r*cosd(delta);
y5 = y0 + r*sind(delta);


x=@(t) (1-t).^3*x1+3*(1-t).^2.*t*x2+3*(1-t).*t.^2*x3+t.^3*x4;
y=@(t) (1-t).^3*y1+3*(1-t).^2.*t*y2+3*(1-t).*t.^2*y3+t.^3*y4;

t=0:0.01:1;

figure(1)
plot([x1,x2,x3,x4],[y1,y2,y3,y4],'k.-','linewidth',1.5,'markersize',20)
hold 
plot(x0,y0,'r.','markersize',15)
plot(x4,y4,'k.','markersize',20)
plot(x(t),y(t),'b','linewidth',1.5)
plot([x0,x4],[y0,y4],'r','linewidth',1.5)
plot(x5,y5,'r','linewidth',1)
hold off
grid
axis equal
xlabel('x','fontsize',14)
ylabel('y  ','rotation',0,'fontsize',14)

dx=-3*(1-t).^2*x1+...
   3*(1-4*t+3*t.^2)*x2+...
   3*(2*t-3*t.^2)*x3+...
   3*t.^2*x4;
d2x=6*(1-t)*x1+...
   3*(-4+6*t)*x2+...
   3*(2-6*t)*x3+...
   6*t*x4;
dy=-3*(1-t).^2*y1+...
   3*(1-4*t+3*t.^2)*y2+...
   3*(2*t-3*t.^2)*y3+...
   3*t.^2*y4;
d2y=6*(1-t)*y1+...
   3*(-4+6*t)*y2+...
   3*(2-6*t)*y3+...
   6*t*y4;

k=(dx.*d2y-d2x.*dy)./(dx.^2+dy.^2).^(3/2);

figure(2)

plot(t,k,'linewidth',1.5)
grid
title('kaarevuus K(t)')
xlabel('t')
ylabel('curvature \kappa(t)')

xlabel('t')
ylabel('y')

%%
%ex4
clear 
close all

B1=[0,0];
B2=[1,2];
B3=[2,1];
B4=[3,3];
B5=[4,0];
B6=[2,-1];


x1=B1(1);
y1=B1(2);
x2=B2(1);
y2=B2(2);
x3=B3(1);
y3=B3(2);
x4=B4(1);
y4=B4(2);
x5=B5(1);
y5=B5(2);
x6=B6(1);
y6=B6(2);


S6= 1/6 * (B6 + 4*B1 + B2);
P1=2/3*B1+1/3*B2;
P2=1/3*B1+2/3*B2;
S1= 1/6 * (B1 + 4*B2 + B3);

S1;
P3=2/3*B2+1/3*B3;
P4=1/3*B2+2/3*B3;
S2= 1/6 * (B2 + 4*B3 + B4);

S2;
P5=2/3*B3+1/3*B4;
P6=1/3*B3+2/3*B4;
S3= 1/6 * (B3 + 4*B4 + B5);

S3;
P7=2/3*B4+1/3*B5;
P8=1/3*B4+2/3*B5;
S4= 1/6 * (B4 + 4*B5 + B6);

S4;
P9=2/3*B5+1/3*B6;
P10=1/3*B5+2/3*B6;
S5= 1/6 * (B5 + 4*B6 + B1);

S5;
P11=2/3*B6+1/3*B1;
P12=1/3*B6+2/3*B1;
S6;

X1=S6(1);
Y1=S6(2);
X2=P1(1);
Y2=P1(2);
X3=P2(1);
Y3=P2(2);
X4=S1(1);
Y4=S1(2);

Z1=S1(1);
T1=S1(2);
Z2=P3(1);
T2=P3(2);
Z3=P4(1);
T3=P4(2);
Z4=S2(1);
T4=S2(2);

M1=S2(1);
N1=S2(2);
M2=P5(1);
N2=P5(2);
M3=P6(1);
N3=P6(2);
M4=S3(1);
N4=S3(2);

O1=S3(1);
P1=S3(2);
O2=P7(1);
P2=P7(2);
O3=P8(1);
P3=P8(2);
O4=S4(1);
P4=S4(2);

L1=S4(1);
K1=S4(2);
L2=P9(1);
K2=P9(2);
L3=P10(1);
K3=P10(2);
L4=S5(1);
K4=S5(2);

G1=S5(1);
H1=S5(2);
G2=P11(1);
H2=P11(2);
G3=P12(1);
H3=P12(2);
G4=S6(1);
H4=S6(2);


t=0:0.01:1;
x=@(t) (1-t).^3*x1+3*(1-t).^2.*t*x2+3*(1-t).*t.^2*x3+ 3*(1-t).*t.^2*x4+ 3*(1-t).*t.^2*x5+ t.^3*x6;
y=@(t) (1-t).^3*y1+3*(1-t).^2.*t*y2+3*(1-t).*t.^2*y3+ 3*(1-t).*t.^2*y4+ 3*(1-t).*t.^2*y5+t.^3*y6;

X=(1-t).^3*X1+3*(1-t).^2.*t*X2+3*(1-t).*t.^2*X3+t.^3*X4;
Y=(1-t).^3*Y1+3*(1-t).^2.*t*Y2+3*(1-t).*t.^2*Y3+t.^3*Y4;

Z=(1-t).^3*Z1+3*(1-t).^2.*t*Z2+3*(1-t).*t.^2*Z3+t.^3*Z4;
T=(1-t).^3*T1+3*(1-t).^2.*t*T2+3*(1-t).*t.^2*T3+t.^3*T4;

M=(1-t).^3*M1+3*(1-t).^2.*t*M2+3*(1-t).*t.^2*M3+t.^3*M4;
N=(1-t).^3*N1+3*(1-t).^2.*t*N2+3*(1-t).*t.^2*N3+t.^3*N4;

O=(1-t).^3*O1+3*(1-t).^2.*t*O2+3*(1-t).*t.^2*O3+t.^3*O4;
P=(1-t).^3*P1+3*(1-t).^2.*t*P2+3*(1-t).*t.^2*P3+t.^3*P4;

L=(1-t).^3*L1+3*(1-t).^2.*t*L2+3*(1-t).*t.^2*L3+t.^3*L4;
K=(1-t).^3*K1+3*(1-t).^2.*t*K2+3*(1-t).*t.^2*K3+t.^3*K4;

G=(1-t).^3*G1+3*(1-t).^2.*t*G2+3*(1-t).*t.^2*G3+t.^3*G4;
H=(1-t).^3*H1+3*(1-t).^2.*t*H2+3*(1-t).*t.^2*H3+t.^3*H4;
figure(2)

plot([x1,x2],[y1,y2],'k','linewidth',1.5)
hold on
plot([x2,x3],[y2,y3],'k','linewidth',1.5)
hold on
plot([x3,x4],[y3,y4],'k','linewidth',1.5)
hold on
plot([x4,x5],[y4,y5],'k','linewidth',1.5)
hold on
plot([x5,x6],[y5,y6],'k','linewidth',1.5)
hold on
plot([x6,x1],[y6,y1],'k','linewidth',1.5)
hold on

plot(X,Y,'-','linewidth',1.5)
hold on
%plot([X1,X2,X3,X4],[Y1,Y2,Y3,Y4],'b.-','markersize',20)

plot(Z,T,'-','linewidth',1.5)
hold on
%plot([Z1,Z2,Z3,Z4],[T1,T2,T3,T4],'b.-','markersize',20)

plot(M,N,'-','linewidth',1.5)
hold on
%plot([M1,M2,M3,M4],[N1,N2,N3,N4],'b.-','markersize',20)

plot(O,P,'-','linewidth',1.5)
hold on
%plot([O1,O2,O3,O4],[P1,P2,P3,P4],'b.-','markersize',20)

plot(L,K,'-','linewidth',1.5)
hold on
%plot([L1,L2,L3,L4],[K1,K2,K3,K4],'b.-','markersize',20)

plot(G,H,'-','linewidth',1.5)
hold on
%plot([G1,G2,G3,G4],[H1,H2,H3,H4],'b.-','markersize',20)
grid
axis equal
xlabel('x','fontsize',14)
ylabel('y  ','rotation',0,'fontsize',14)

%% position_velocity_acceleration
% ex1
clear
close all
%position
s=@(t) sqrt((5+2*t).^2-4^2);
%velocity
v=@(t) (10+4*t)./sqrt(4*t.^2+20*t+9);
%acceleration
a=@(t) -64./(4*t.^2+20*t+9).^(3/2);
%values at time t0
t0=3;
s(t0)
v(t0)
a(t0)

%graphs when t=0...T
%AP=9->3+2*T=sqrt(9^2-4^2)=sqrt(65)->
T=2;
t=0:T/100:T;

subplot(311)
plot(t,s(t),'linewidth',1.5)
grid
xlim([0,T])
title('position')
subplot(312)
plot(t,v(t),'linewidth',1.5)
grid
xlim([0,T])
title('velocity')
subplot(313)
plot(t,a(t),'linewidth',1.5)
grid
xlim([0,T])
title('acceleration')
xlabel('time t')

%%
%ex2
clear
close all
h = 5;
x = pi/6;
w = 2*pi;
alpha = -5*pi;
%postion
s2 =  h*tan(x);
% velocity
v2 = h*w*1/(cos(x))*1/(cos(x));
% acceleration
a2 = h*1/(cos(x))*1/(cos(x))*(2*w^2*tan(x)+alpha);

%%
%ex4
clear
r=2;
V=3;
T=4.1888;
t0=0.3*T;
th= V*t0/r;
x0=r*th;
y0=r;

 

t=0:T/100:T;
a=V*t./r;
d=0:1:360;

 

%center of point of the the circle
xc=x0+r*cosd(d);
yc=y0+r*sind(d);

 

x2 = r*a;
xk=x0-r*sin(th);
yk=y0-r*cos(th);
x=@(t) V*t-r*sin(V*t./r);
y=@(t) r-r*cos(V*t./r);

 

vx=@(t) V - V*cos(V*t./r);
vy=@(t) V*sin(V*t./r);
v=@(t) sqrt(vx(t).^2+vy(t).^2);

 

ax=@(t) V^2*1/r*sin(V*t./r);
ay=@(t) V.^2*1/r*cos(V*t./r);
a= @(t) sqrt(ax(t).^2+ay(t).^2);
aT= @(t) (vx(t).*ax(t)+vy(t).*ay(t))./v(t);
aN= @(t) (vx(t).*ay(t)-vy(t).*ax(t))./v(t);

 


figure (1)
plot(x(t),y(t),'b','linewidth',1.5)
hold
plot(xc,yc,'r','linewidth',0.5)
plot(x0,y0,'r.','markersize',10)
plot(xk,yk,'b.','markersize',10)
plot([0,2*pi*r],[0,0],'k')
p2=plot([x(t0),x(t0)+ax(t0)],[y(t0),y(t0)+ay(t0)],'g','linewidth',2);
p1=plot([x(t0),x(t0)+vx(t0)],[y(t0),y(t0)+vy(t0)],'r','linewidth',2);
p4=plot([x(t0),x(t0)+aN(t0)*(-vy(t0))/v(t0)],[y(t0),y(t0)+aN(t0)*vx(t0)/v(t0)],'m','linewidth',2);
p3=plot([x(t0),x(t0)+aT(t0)*vx(t0)/v(t0)],[y(t0),y(t0)+aT(t0)*vy(t0)/v(t0)],'c','linewidth',2);
hold off
grid
axis equal
legend([p1,p2,p3,p4],'v','a','a_T','a_N')

 

figure(2)
subplot(2,1,1)
plot(t,v(t),'linewidth',1.5)
grid
title('|v|')
xlim([0,T])
subplot(2,1,2)
plot(t,a(t),'linewidth',1.5)
grid
title('|a|')
xlabel('time t')
xlim([0,T])
ylim([2,7])

 

figure(3)
subplot(2,1,1)
plot(t,aT(t),'linewidth',1.5)
grid
title('tangential acceleration a_T')
xlim([0,T])
subplot(2,1,2)
plot(t,aN(t),'linewidth',1.5)
grid
title('normal acceleration a_N')
xlabel('time t')
xlim([0,T])
%% extremum_problems
%ex1
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
%ex2
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
%ex3
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
%ex4
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
%ex5
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
%ex6
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
%% integral
%ex1
clear
close all
a = 0.7;
h = 4;
f=@(x) h;
g=@(x) a*x.^2;
%upper limit x0
x0= h^(1/2)/a^(1/2);
%lower limit at x1
x1 = 0;
%area
A= (a*x1^3-3*h*x1)/3+(2*h^(3/2))/(3*sqrt(a));
% Centerpoint 
xp = 1/A*((a^2*x1^4-2*a*h*x1^2+h^2)/(4*a));
yp = 1/A*(((a^2*x1^5-5*h^2*x1)/5+ (4*h^(5/2))/(5*sqrt(a)))/2);
x = 0:1/100:x0;
plot([x1,x0], [h,h], 'b','linewidth',1.5)
hold on
plot (x, g(x), 'r', 'linewidth', 1.5)
plot ([0,0], [0,h], 'k', 'linewidth',1.5)
plot (xp, yp, 'k.', 'markersize', 20)
hold off
grid
axis equal
%% 
%ex4
clear
close all
%ex4
h = 3;
L = 4;
V = -(sqrt(L)*(3*L^2-20*L))/10;

%% taylor_polynomial

clear
close all
%degree n
n=5;
%graphs for x=a...b
a=-1;
b=1;
x=a:(b-a)/1000:b;
f=@(x) log(1+x);
x0=0;

Tn=0;

%factorial(k)=k!=k*(k-1)*...*2*1
for k=1:n
    Tn=Tn+(-1)^(k-1)/(k)*x.^(k); 
end

plot(x,f(x),'b','linewidth',1.5)
hold
plot(x,Tn,'r','linewidth',1.2)
plot(x0,f(x0),'k.','markersize',20)
hold off
grid
xlim([-1,1])
ylim([-2,2])
xlabel('x')
legend({'log(1+x)','T_n(x)'},'fontsize',12)
title(['degree n = ',num2str(n)])
