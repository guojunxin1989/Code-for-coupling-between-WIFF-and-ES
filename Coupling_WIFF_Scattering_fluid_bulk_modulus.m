clear;
close;
clc;

Km=26;%bulk modulus of the background medium,GPa
Gm=31;%shear modulus of the background medium
por=0.1;%porosity of the background medium
Ks=37;%bulk modulus of the solid 
Gs=44;%shear modulus of the solid
L=Km+4/3*Gm;
g=Gm/L;
a=1-Km/Ks;
den=2.49;%density of the fluid saturated background medium,g/cc
vis=0.001*10^(-9);%viscosity of fluid, Gpa*s
perm=10^(-15);%permeability of background, m^2
b=vis/perm;

d=1;%diameter of the fracture,m
r=d/2;%fracture radius
h=0.01;%fracture thickness,m
cden=0.05;
n0=cden/r^3;
 
f=-3:0.1:5;
f=10.^f;
Df=length(f);
Kfa=[3,2,1,0.00001];
N=length(Kfa);
v=zeros(N,Df);
inQ=zeros(N,Df);

for dN=1:N
    
    Kf=Kfa(dN);
    M=((a-por)/Ks+por/Kf)^(-1);
    H=L+a^2*M;

for FN=1:Df

FN
ft=f(FN);
w=2*pi*ft;
k0=10^(-3)*sqrt(den*w^2/L);%1/m
k1=10^(-3)*sqrt(den*w^2/H);%1/m
k2=sqrt(1i*w*b*H/(L*M));%1/m
k3=10^(-3)*sqrt(den*w^2/Gm);%1/m

abk=(2*pi)/d*2;
DM=20000;
y=abk/DM:abk/DM:abk;

%y=0.01:0.01:200;

YN=length(y);

q1=zeros(1,YN);
q2=zeros(1,YN);
q3=zeros(1,YN);

for lk=1:YN
    
yy=y(lk);
if yy<=k1
    q1(lk)=-1i*sqrt(k1^2-yy^2);
else
    q1(lk)=sqrt(yy^2-k1^2);
end
    q2(lk)=-1i*sqrt(k2^2-yy^2);
if yy<=k3
    q3(lk)=-1i*sqrt(k3^2-yy^2);
else
    q3(lk)=sqrt(yy^2-k3^2);
end

end

%T11=(k0^2-k1^2)*((M*L*h)./(q1*H)-2*Kf/k2^2*(1+(h*y.^2)./(2*q1))).*(2*y.^2-k3^2)+2*k3^2*Kf*a*M/H;
T11=(k0^2-k1^2)*((M*L*h)./(q1*H)-2*Kf/k2^2).*(2*y.^2-k3^2)+2*k3^2*Kf*a*M/H;
T12=(4*a*g*y.^2).*(y.^2-q2.*q3)-k2^2*(2*y.^2-k3^2);
%T21=q2.*(k2^2*(2*Kf/k2^2*(1+(h*y.^2)./(2*q2))-(M*L*h)./(q2*H)).*(2*y.^2-k3^2)+2*k3^2*Kf*a*M/H);
T21=q2.*(k2^2*(2*Kf/k2^2-(M*L*h)./(q2*H)).*(2*y.^2-k3^2)+2*k3^2*Kf*a*M/H);
T22=2*a*g*(1-g)*k3^2*y;
T1=(T11.*T12)./(T21.*T22);        
T33=a*g*((2*y.^2-k3^2).^2-4*((y.^2).*q1).*q3)+(2*y.^2-k3^2)*(k0^2-k1^2);
T44=2*a*g*(1-g)*k3^2*q1.*y;
T2=T33./T44;
T=(1+(a*M*k3^2)./(H*(2*y.^2-k3^2))).*(T1-T2)-1;

N=50;
ds=r/N;
Ma=zeros(N,N);
V=zeros(N,1);
for i=1:N
    for j=1:N
       si=(i-1)*ds;
       sj=(j-1)*ds;
       
       TS=(T.*sin(si*y)).*sin(sj*y);
       
       My=trapz(y,TS)*2/pi;
       if(i==j)
           Ma(i,j)=My*ds+1;
       else
           Ma(i,j)=My*ds;
       end
       
    end
    
    p0=H-a*M;
    V(i)=-p0*si;
end

P=Ma\V;
sum=0;
for i=1:N
    si=(i-1)*ds;
    sum=sum+si*P(i)*ds;
end

sum=sum*2/pi;
g=Gm/L;
cof=-2*a*Gm*(1-g)/(L*(1-a*M/H));
A1=sum/cof;
f0=a*k1^2/L*A1;
ke=k1*(1+(2*pi*n0)/k1^2*f0);
v(dN,FN)=w/real(ke);
inQ(dN,FN)=2*imag(ke)/real(ke);

end

end


figure(1)
semilogx(f,v(1,:),'c','linewidth',1.5);
hold on;
semilogx(f,v(2,:),'r','linewidth',1.5);
hold on;
semilogx(f,v(3,:),'g','linewidth',1.5);
hold on;
semilogx(f,v(4,:),'b','linewidth',1.5);

xlabel('Frequency (Hz)','FontSize',12);
ylabel('{\itP}-wave velocity (m/s)','FontSize',12);
legend('{\itK_{f}}=3 GPa','{\itK_{f}}=2 GPa','{\itK_{f}}=1 GPa','{\itK_{f}}=0 GPa');
set(gca, 'FontSize', 12);
xlim([0.001 100000]);
set(gca,'XTick',[0.001 0.01 0.1 1 10 100 1000 10000 100000]);

figure(2)
loglog(f,inQ(1,:),'c','linewidth',1.5);
hold on;
loglog(f,inQ(2,:),'r','linewidth',1.5);
hold on;
loglog(f,inQ(3,:),'g','linewidth',1.5);
hold on;
loglog(f,inQ(4,:),'b','linewidth',1.5);

xlabel('Frequency (Hz)','FontSize',12);
ylabel('Q^{-1}','FontSize',12);
legend('{\itK_{f}}=3 GPa','{\itK_{f}}=2 GPa','{\itK_{f}}=1 GPa','{\itK_{f}}=0 GPa');
set(gca, 'FontSize', 12);
xlim([0.001 100000]);
ylim([0.001 1]);
set(gca,'XTick',[0.001 0.01 0.1 1 10 100 1000 10000 100000]);







 
 

