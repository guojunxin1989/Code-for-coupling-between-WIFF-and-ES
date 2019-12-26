clear;
close;
clc;

Km=26;%bulk modulus of the background medium,GPa
Gm=31;%shear modulus of the background medium
por=0.1;%porosity of the background medium
Ks=37;%bulk modulus of the solid 
Gs=44;%shear modulus of the solid
Kf=2.25;%bulk modulus of the fluid
L=Km+4/3*Gm;
g=Gm/L;
a=1-Km/Ks;
M=((a-por)/Ks+por/Kf)^(-1);
H=L+a^2*M;
den=2.49;%density of the fluid saturated background medium,g/cc
vis=0.001*10^(-9);%viscosity of fluid, Gpa*s
perm=10^(-15);%permeability of background, m^2
b=vis/perm;

d=1;%diameter of the fracture,m
r=d/2;%radius of the fracture
h=0.01;%fracture thickness,m
cden=0.05;
n0=cden/r^3;
 
f=-3:0.1:5;
f=10.^f;
Df=length(f);
v=zeros(Df,1);
inQ=zeros(Df,1);

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
v(FN)=w/real(ke);
inQ(FN)=2*imag(ke)/real(ke);

end

Kmf=0.00000001;
Gmf=0.00000001;
a1=1;
a2=1;
a3=0.01;
pf=4/3*pi*a3*cden;
opt=2;
Sf=General_Eshelby_model(Km,Gm,Kmf,Gmf,a1,a2,a3,pf,opt);
ZN=Sf(3,3);
ZT=Sf(4,4);
La=L-2*Gm;
C0=[L,La,La,0,0,0;
    La,L,La,0,0,0;
    La,La,L,0,0,0;
    0,0,0,Gm,0,0;
    0,0,0,0,Gm,0;
    0,0,0,0,0,Gm;
    ];
S0=inv(C0);
S1=S0+Sf;
C1=inv(S1);
Cs1=anisotropy_Gassmann(C1,por,Ks,Kf);
Vp1=sqrt(Cs1(3,3)/den)*1000;

Kms=H-4/3*Gm;
Sf=General_Eshelby_model(Kms,Gm,Kmf,Gmf,a1,a2,a3,pf,opt);
Ha=H-2*Gm;
C0=[H,Ha,Ha,0,0,0;
    Ha,H,Ha,0,0,0;
    Ha,Ha,H,0,0,0;
    0,0,0,Gm,0,0;
    0,0,0,0,Gm,0;
    0,0,0,0,0,Gm;
    ];
S0=inv(C0);
S1=S0+Sf;
C1=inv(S1);
Cs2=anisotropy_Gassmann(C1,pf,Kms,Kf);
Vp2=sqrt(Cs2(3,3)/den)*1000;

T=(2*(H-a*M)^2*(2-4*a*g+3*a^2*g^2)*r^2*cden*vis)/(15*Gm*g*(1-g)^2*H*L*perm);
Sc=(pi*cden)/r;
C2=Cs2(3,3);
C1=Cs1(3,3);
Lc=pf/ZN;
Gc=pf/ZT;
Kc=Lc-4/3*Gc;
ac=1-Kc/Ks;
porf=1;
Mc=((ac-porf)/Ks+porf/Kf)^(-1);
Hc=Lc+ac^2*Mc;
G=(2*Sc*C2*(a*M/H-ac*Mc/Hc)^2)/sqrt((M*L*vis)/(H*perm));
% G=2*Sc*(H-a*M)^2*C2/H*sqrt(perm/(vis*H*M*L));
F1=((C2-C1)/(C1*G))^2;
F2=(C2-C1)^3/(2*C2*C1^2*T*G^2);
F3=(C2-C1)/C1;
s=1/C2*(1+F3./(1-F2+F2*sqrt(1-1i*(2*pi*f*F1)/F2^2)));
c=1./s;
Vpg=(real(1./sqrt(c/den))).^(-1)*1000;
inQg=abs(imag(c)./real(c));

Vp3=sqrt(H/den)*1000;

N=length(v);
Vp1a=zeros(N,1)+Vp1;
Vp2a=zeros(N,1)+Vp2;
Vp3a=zeros(N,1)+Vp3;

figure(1)
semilogx(f,v,'r','linewidth',1.5);
hold on;
semilogx(f,Vpg,'r--','linewidth',1.5);
hold on;
semilogx(f,Vp1a,'g','linewidth',1.5);
hold on;
semilogx(f,Vp2a,'b','linewidth',1.5);
hold on;
semilogx(f,Vp3a,'k','linewidth',1.5);

xlabel('Frequency (Hz)','FontSize',12);
ylabel('{\itP}-wave velocity (m/s)','FontSize',12);
legend('results of this paper','results of Guo et al.(2018a)');
set(gca, 'FontSize', 12);
xlim([0.001 100000]);
set(gca,'XTick',[0.001 0.01 0.1 1 10 100 1000 10000 100000]);

figure(2)
loglog(f,inQ,'r','linewidth',1.5);
hold on;
loglog(f,inQg,'r--','linewidth',1.5);

xlabel('Frequency (Hz)','FontSize',12);
ylabel('Q^{-1}','FontSize',12);
legend('results of this paper','results of Guo et al.(2018a)');
set(gca, 'FontSize', 12);
xlim([0.001 100000]);
set(gca,'XTick',[0.001 0.01 0.1 1 10 100 1000 10000 100000]);











 
 

