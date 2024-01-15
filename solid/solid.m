clear all
clc
%%input the single event[Events]
%%Time interval of outputs in second(1/6h,10mins;1/4,15mins)
%%output
%%[./result/stress33.92N133.72E30 km250strike.txt]:
%%normal stress,shear stress,first stress invariant
%%[./result/strain33.92N133.72E30 km250strike.txt]
%%strain tensor includes 6 independent elements,total 9
%%[./result/stresstensor33.92N133.72E30 km250strike.txt]
%%stress tensor includes 6 independent elements,total 9
tic %%it may need 20s to calculate the 4465 data
% Events=[ lat, lon, depth, strike, dip, rake, year, month, day, hour，min,sec];
% Events=[39.03,140.88,8,209,39,101,2008, 5,31, 15,0,0];%%case1
Events=[33.92,133.72,30,250,10,100,2008, 5,31, 15,0,0];%%case2
Ra=6371000.000;
%%read miu and lambda at r.
earthmodel=load("./model/earthmodel.dat");
miu=earthmodel(:,3);
lambda=earthmodel(:,4);
%%read love number or y1~y4 at depth and calcuate the derivatives of y1 and y3 with respect to r.
TLNs=ReadLoveNumber("./model/Tide_Love.dat",earthmodel);
%%read parameters from tide generating potential table.
[Dood,fre,co_K]=ReadPotential('./model/potential.dat');

lat=Events(1);       lon=Events(2);       depth=Events(3);
strike=Events(4);  dip=Events(5);       rake=Events(6);
year= Events(7);   month=Events(8); day=Events(9);     hour=Events(10);
min=Events(11);sec=Events(12);
min=min-15;
julday=juliandate(year,month,day,hour,min,sec);%%calculate Juliandate from 12:00 January 1, 4713 BC
hour2=hour+(min+sec/60)/60;
phase_t= phase(lon,julday, hour2,Dood, co_K);
[G0,G1,G2,Gr]=Cal_Geo(lat, depth);
N=length(phase_t);
Rx=Ra-depth*1000;
ID0=(Ra-Rx)/100.0;
ID=floor(ID0+1);%index of array of earthmodel parameters. r in meter
for h=1:4465           %the number of outputs
    stress2=zeros(9,1);
    strain2=zeros(9,1);
    phase_t=phase_t+fre*1/4;%phi=phi0+w*t(1/6h,10mins;1/4,15mins)Time interval of outputs in second
    for n=1:N
        degree=Dood(n,1);
        order=Dood(n,2);
        strain=Cal_strain( degree, order, lat, Rx, phase_t,TLNs, G0,G1,G2,Gr,co_K(n,1),ID,n);
        strain2=strain2+strain;
        stress=Cal_stress(strain,lambda,miu,ID);
        stress2=stress2+stress;
    end
    CC(h,:)=stress2;%%stress tensor includes 6 independent elements,total 9
    DD(h,:)=strain2;%%strain tensor includes 6 independent elements,total 9
    %%rotates to the fault plane
    Series(:)=Rotation( strike, dip, rake, stress2);
    S(h,:)= Series;
end
S2(:,1)=S(:,1);
S2(:,2)=S(:,2);
S2(:,3)=(S(:,1)+S(:,5)+S(:,9))/3;
newDir='./results';%%the results save in this file
%mkdir(newDir);
cd(newDir);
filename=['stress' num2str(lat) 'N' num2str(lon) 'E' num2str(depth) ' km' num2str(strike) 'strike' '.txt'];
fid=fopen(filename,'w');
fprintf(fid,'%16.6E  %16.6E  %16.6E\n',S2');%%normal stress,shear stress,first stress invariant
fclose(fid);

filename=['strain' num2str(lat) 'N' num2str(lon) 'E' num2str(depth) ' km' num2str(strike) 'strike' '.txt'];
fid=fopen(filename,'w');
fprintf(fid,'%16.6E  %16.6E  %16.6E  %16.6E  %16.6E  %16.6E  %16.6E  %16.6E  %16.6E\n',DD');

filename=['stresstensor' num2str(lat) 'N' num2str(lon) 'E' num2str(depth) ' km' num2str(strike) 'strike' '.txt'];
fid=fopen(filename,'w');
fprintf(fid,'%16.6E  %16.6E  %16.6E  %16.6E  %16.6E  %16.6E  %16.6E  %16.6E  %16.6E\n',CC');
toc


function  love=ReadLoveNumber(filename1,model)
% read love number or y1~y4 at depth and calcuate the derivatives of y1 and y3 with respect to r.
% convert to y1,y2,y3,y4,y1',y3'
Ra=6371000;
tlns=load(filename1);
Rr=model(:,1).*1000;
miu=model(:,3);
lambda=model(:,4);
N=length(Rr);
Rn=zeros(N,3);
Rn(:,1)=Rr.*Rr/Ra/Ra;
Rn(:,2)=Rn(:,1).*Rr/Ra;
Rn(:,3)=Rn(:,2).*Rr/Ra;
lam2miu=lambda+2*miu;
LMR=lambda./lam2miu./Rr;
love2=zeros(8001,3,4);
love=zeros(8001,3,6);
for j=1:3 %%degree 2,3,4
    for k=1:4 %%y1,y2,y3,y4
        love2(:,:,k)=reshape(tlns(:,k),3,8001)';
        love(:,j,k)=love2(:,j,k).*Rn(:,j);
    end
    love(:,j,5)=-2*LMR(:).*love(:,j,1)+love(:,j,2)./lam2miu(:)+(j+1)*(j+2).*LMR(:).*love(:,j,3);  % dy1/dr
    love(:,j,6)=-love(:,j,1)./Rr(:)+love(:,j,3)./Rr(:)+love(:,j,4)./miu(:);    % dy3/dr
end
end

function   [Dood,fre,co_K]=ReadPotential(filename2)
% read parameters from tide generating potential table.
fid = fopen(filename2, 'r');
podata = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines',69, 'CollectOutput', 1);
data = podata{1};
fclose(fid);
potential=data(1:end-1,:);
N=length(potential);
Dood=potential(:,2:8);
fre=potential(:,14);
COS=potential(:,15);
SIN=potential(:,16);
trans=zeros(5,5);
% transfer coefficient from HW's tidal potential expansion to Doodson's
% reference: Xi Q W. Chinese Journal of Geophysics,2007, 50(2):111-114.
% see table 2 but with modification on 20 30 31 41 42 with minus sign
trans(3,1)=-1.17778;
trans(3,2)=1.35998;
trans(3,3)=1.35998;
trans(4,1)=-2.22578;
trans(4,2)=-1.18040;
trans(4,3)=1.33547;
trans(4,4)=1.25909;
trans(5,1)=0.87786;
trans(5,2)=-1.05150;
trans(5,3)=-1.22140;
trans(5,4)=1.29233;
trans(5,5)=1.18708;
co_K=zeros(N,2);
for n=1:N
    if (abs(COS(n))<0.000001)&&(abs(SIN(n))<0.000001)
        error('error may occur in reading potential table!');
    end
    if abs(COS(n))<0.000001
        co_K(n,1)=SIN(n)/1.0e+10;
        co_K(n,2)=-90.0;
    else
        co_K(n,1)=COS(n)/1.0e+010;
        co_K(n,2)=0.0;
    end
    D1=Dood(n,1)+1;
    D2=Dood(n,2)+1;
    co_K(n,1)=co_K(n,1)/trans(D1,D2 );
end
end

function  phase_t= phase(lon,julday, hour,Dood, co_K)
%calculate astronomic argument according to Doodson codes
T=zeros(6,1);
DETMUT=51.0;
TimeSys=0;%time zone
TD = julday - TimeSys/24.0 + DETMUT/86400.0-2415020;
% 2415020 is the juliandate of 12:00 December 31, 1899
FS = TD/36525.0;
FS2 = FS*FS;
FS3 = FS2*FS;
T(2) = 270.434162 + 481267.883142*FS - 0.0011333*FS2 +0.000001889*FS3;
T(3) = 279.696681 + 36000.768925*FS + 0.0003027*FS2;
T(4) = 334.329568 + 4069.034034*FS - 0.0103249*FS2 - 0.0000125*FS3;
T(5) = 100.843202 + 1934.142008*FS - 0.002078 *FS2 - 0.000002 *FS3;
T(6) = 281.220868 + 1.719175*FS + 0.0004527*FS2 + 0.0000033*FS3;
T(1) = T(3) - T(2) + (hour-TimeSys)*15.0 + lon;
phase_t=zeros(2934,1);
for i=1:2934
    R = 0.0;
    for j=1:6
        R =R +T(j)*Dood(i,j+1);
    end
    phase_t(i) = R - 360.0*floor(R/360.0);
    phase_t(i)=phase_t(i)+co_K(i,2);         % convert sin to cos in coefficients part
end
end

function   [G0,G1,G2,Gr]=Cal_Geo(lat, depth)
Ra=6371000.000;
G0=zeros(5,5);
G1=G0;
G2=G1;
Gr=G2;
Rr=Ra-depth*1000;
dr=atan(1.0)/45.0;
p=dr*lat;
r2=(Rr/Ra)*(Rr/Ra);
r3=r2*(Rr/Ra);
r4=r2*r2;
s1=sin(p);
s2=s1*s1;
s3=s2*s1;
s4=s2*s2;
c1=cos(p);
c2=c1*c1;
c3=c2*c1;
c4=c2*c2;
Doodson=2.6206;
r2=r2*Doodson;
r3=r3*Doodson;
r4=r4*Doodson;
% geodetic coefficients
G0(3,1)=1.5*(1.0/3-s2 )*r2;
G0(3,2)=sin(2*p)*r2;
G0(3,3)=c2*r2;

G0(4,1)=sqrt(5.0)/2.0*s1*(3.0-5*s2)*r3;
G0(4,2)=3*sqrt(15.0)/16.0*c1*(1-5*s2)*r3;
G0(4,3)=3*sqrt(3.0)/2.0*s1*c2*r3;
G0(4,4)=c3*r3;

G0(5,1)=1.0/8*(3-30*s2+35*s4)*r4;
G0(5,2)=224.0/((3+sqrt(393.0))*sqrt(390.0+2*sqrt(393.0)))*sin(2*p)*(3-7.0*s2)*r4;
G0(5,3)=7.0/9.0*c2*(1-7*s2)*r4;
G0(5,4)=16.0/sqrt(27.0)*s1*c3*r4;
G0(5,5)=c4*r4;
% first order derivatives of geodetic coefficient with respect to r (multiplied by n/r)
Gr(3,1)=1.5*(1.0/3-s2 )*r2 *2/Rr;
Gr(3,2)=sin(2*p)*r2 *2/Rr;
Gr(3,3)=c2*r2 *2/Rr;

Gr(4,1)=sqrt(5.0)/2.0*s1*(3.0-5*s2)*r3 *3/Rr;
Gr(4,2)=3*sqrt(15.0)/16.0*c1*(1-5*s2)*r3 *3/Rr;
Gr(4,3)=3*sqrt(3.0)/2.0*s1*c2*r3 *3/Rr;
Gr(4,4)=c3*r3 *3/Rr;

Gr(5,1)=1.0/8*(3-30*s2+35*s4)*r4 *4/Rr;
Gr(5,2)=224.0/((3+sqrt(393.0))*sqrt(390.0+2*sqrt(393.0)))*sin(2*p)*(3-7.0*s2)*r4 *4/Rr;
Gr(5,3)=7.0/9.0*c2*(1-7*s2)*r4 *4/Rr;
Gr(5,4)=16.0/sqrt(27.0)*s1*c3*r4 *4/Rr;
Gr(5,5)=c4*r4 *4/Rr;
% first order derivatives of geodetic coefficient with respect to sita (represented by latitude)
G1(3,1)=1.5*( -2*s1*c1 ) *r2;
G1(3,2)=( 2*cos(2*p) ) *r2;
G1(3,3)=-( 2*c1*s1 )    *r2;

G1(4,1)=sqrt(5.0)/2.0*( 3.0*c1-15.0*s2*c1 )*r3;
G1(4,2)=3*sqrt(15.0)/16.0*( -s1+5*s3-10*s1*c2 )*r3;
G1(4,3)=3*sqrt(3.0)/2.0*( c3-2*c1*s2 ) *r3;
G1(4,4)=(-3.0*c2*s1 )*r3;

G1(5,1)=1.0/8*( -60.0*s1*c1+140*s3*c1 ) *r4;
G1(5,2)=224.0/((3+sqrt(393.0))*sqrt(390.0+2*sqrt(393.0)))*( 2*cos(2*p)*(3-7*s2)-sin(2*p)*7*sin(2*p) )  *r4;%%-
G1(5,3)=7.0/9.0*( -sin(2*p)*(1-7*s2)-c2*7*sin(2*p) ) *r4;
G1(5,4)=16.0/sqrt(27.0)*( c4-3*c2*s2 )*r4;
G1(5,5)=( -4*c3*s1 )*r4;

% second order derivatives of geodetic coefficient with respect to sita (represented by latitude)
G2(3,1)=1.5*( -2*cos(2*p) ) *r2;
G2(3,2)=( -4*sin(2*p) )   *r2;
G2(3,3)=( -2*cos(2*p) )   *r2;

G2(4,1)=sqrt(5.0)/2.0*( -3*s1+15*s3-30*s1*c2 ) *r3;
G2(4,2)=3*sqrt(15.0)/16.0*( -c1+35*s2*c1-10*c3 ) *r3;
G2(4,3)=3*sqrt(3.0)/2.0*( -7*s1*c2+2*s3 )*r3;
G2(4,4)=(-3.0*c3+6*s2*c1 )*r3;

G2(5,1)=1.0/8*( -60.0*cos(2*p)-140*s4+420*s2*c2 ) *r4;
G2(5,2)=224.0/((3+sqrt(393.0))*sqrt(390.0+2*sqrt(393.0)))*( -4*sin(2*p)*(3-7*s2)...
    -14*cos(2*p)*sin(2*p) -14*sin(4*p) )  *r4;
G2(5,3)=7.0/9.0*( -cos(2*p)*(1-7*s2)+14*sin(2*p)*sin(2*p)-14*c2*cos(2*p) ) *r4;
G2(5,4)=16.0/sqrt(27.0)*( -4*c3*s1-1.5*sin(4*p) ) *r4;
G2(5,5)=( -4*c4+12*c2*s2 ) *r4;

end

function   strain=Cal_strain( degree,  order, lat, r,phase_t,Love, G0,G1,G2,~,K, ID,index)
%includes 6 independent elements,rr,tt,ll,rt,rl,tl. r=radial,t=theta,l=lambda
dr=atan(1.0)/45;
p=90-lat;%%geocenter
c1=cos(dr*p);
s1=sin(dr*p);
s2=s1*s1;
l=degree-1;
n1=degree+1;
n2=order+1;
strain=zeros(1,9);
%rr
strain(1)=Love(ID,l,5)*G0(n1,n2)*K;
strain(1)=strain(1)*cos(dr*phase_t(index));          %cos()
%tt
strain(5)=Love(ID,l,1)/r*G0(n1,n2)*K+Love(ID,l,3)/r*G2(n1,n2)*K;
strain(5)=strain(5)*cos(dr*phase_t(index));          %cos()
%ll
strain(9)=Love(ID,l,1)/r*G0(n1,n2)*K+Love(ID,l,3)*(1.0/r/s2*G0(n1,n2)*(-order*order)...
    -1.0/r*c1/s1*G1(n1,n2))*K;
strain(9)=strain(9)*cos(dr*phase_t(index));          %cos()
%rt
strain(2)=Love(ID,l,1)/r*G1(n1,n2)*K-Love(ID,l,3)/r*G1(n1,n2)*K+Love(ID,l,6)*G1(n1,n2)*K;
strain(2)=strain(2)/2.0;
strain(2)=strain(2)*cos(dr*phase_t(index));          %cos()
strain(4)=strain(2);
%rl
strain(3)=Love(ID,l,1)/r/s1*G0(n1,n2)*(-order)*K+Love(ID,l,3)*(-1.0/r/s1*G0(n1,n2)*(-order) )*K...
    +Love(ID,l,6)/s1*G0(n1,n2)*(-order)*K;
strain(3)=strain(3)/2.0;
strain(3)=strain(3)*sin(dr*phase_t(index));          %sin()
strain(7)=strain(3);
%tl
strain(6)=Love(ID,l,3)*2*(-1/r/s1*G1(n1,n2)*(-order)-c1/s2/r*G0(n1,n2)*(-order))*K;
strain(6)=strain(6)/2.0;
strain(6)=strain(6)*sin(dr*phase_t(index));         %sin()
strain(8)=strain(6);
strain=strain';
end

function   stress=Cal_stress(strain,lambda,miu,ID)
%according to strain-stress constitution
%ID represents the depth
%stress includes 6 independent elements,rr,tt,ll,rt,rl,tl. r=radial,t=theta,l=lambda
stress=zeros(9,1);
dil=strain(1)+strain(5)+strain(9);
stress(1)=lambda(ID)*dil+2*miu(ID)*strain(1);
stress(5)=lambda(ID)*dil+2*miu(ID)*strain(5);
stress(9)=lambda(ID)*dil+2*miu(ID)*strain(9);
stress(2)=2*miu(ID)*strain(2);
stress(4)=stress(2);
stress(3)=2*miu(ID)*strain(3);
stress(7)=stress(3);
stress(6)=2*miu(ID)*strain(6);
stress(8)=stress(6);
end

function   Series=Rotation(strike, dip, rake,stress2)
dr=atan(1.0)/45.0;
s=dr*strike;
d=dr*dip;
r=dr*rake;
Cs=cos(-s);
Ss=sin(-s);
Cd=cos(d);
Sd=sin(d);
Cr=cos(r);
Sr=sin(r);
A=zeros(9,1);          %rotates along with r with angle of rake
B=zeros(9,1);            %rotates along with sita with angle of dip
C=zeros(9,1);            %rotates along with r with angle of -strike
A(1)=1.0;
A(5)=Cr;
A(6)=-Sr;
A(8)=Sr;
A(9)=Cr;
B(1)=Cd;
B(3)=Sd;
B(5)=1.0;
B(7)=-Sd;
B(9)=Cd;
C(1)=1.0;
C(5)=Cs;
C(6)=-Ss;
C(8)=Ss;
C(9)=Cs;
A=reshape(A,3,3);
B=reshape(B,3,3);
C=reshape(C,3,3);
stress2=reshape(stress2,3,3);
Series=(A*B*C)*stress2*(A*B*C)'; % Series=ABC*stress2*(ABC)'.
Series=reshape(Series,9,1);
Series=Series';
end




