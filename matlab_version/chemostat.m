function chemostat

clear
clc
warning('off')

%%% Initial conditions  

xa=0.003;    % Bh consuming threhalose (0.006773 for monoculture)
xb=0;        % Bh consuming glucose
xc=0;        % Bh inactive
xd=0;        % Bh dead
xe=0.003;    % Bt consuming glucose (0.0677 for monoculture)
xf=0;        % Bt consuming mannose
xg=0;        % Bt inactive
xh=0;        % Bt dead
xi=0.003;    % Ri consuming glucose (0.004 for monoculture)
xj=0;        % Ri consuming lactate + acetate
xk=0;        % Ri consuming inactive
xl=0;        % Ri consuming dead

s1=0.6845;    % trehalose  % 0.65541 ?
s2=8.1368;    % pyruvate 
s3=7.4126;    % glucose
s4=1.0000;    % glutamate
s5=0.2723;    % lactate
s6=1.9912;    % acetate
s7=1.0000;    % mannose
s8=0.4450;    % succinate
s9=0.6998;    % formate
s10=0.0037;   % butyrate

S0=[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10];
X0=[xa,xb,xc,xd,xe,xf,xg,xh,xi,xj,xk,xl];
v0=[X0,S0];

%%% Dilution rate and pH

D=0.0666;
pH=5.5;    % pH=0 -> uncontrolled

%%% Time

trans=0;
tend=300;
tstep=0.01;

ttt=[trans tend tstep];


%%% Simulation

[t,v]=run(v0,ttt,D,pH,S0);


%%% Figure

figure(1)
set(figure(1),'position',[600 300 1000 800])
clf

subplot(3,2,1)

plot(t,v(:,1:4),'LineWidth',2); % blautia
L=legend('xa','xb','xc','xd');

%xlim([0 tend])
%ylim([0 1])

H = gca; H.XAxis.FontSize = 14; H.YAxis.FontSize = 14;
set(L,'fontsize',12)

%xlabel('Time (h)','FontSize',20)
ylabel('Blautia','FontSize',20)

subplot(3,2,3)

plot(t,v(:,5:8),'LineWidth',2); % bacteroides
L=legend('xe','xf','xg','xh');

%xlim([0 tend])
%ylim([0 1])

H = gca; H.XAxis.FontSize = 14; H.YAxis.FontSize = 14;
set(L,'fontsize',12)

%xlabel('Time (h)','FontSize',20)
ylabel('Bacteriodes','FontSize',20)


subplot(3,2,5)

plot(t,v(:,9:12),'LineWidth',2); % roseburia
L=legend('xi','xj','xk','xl');

%xlim([0 tend])
%ylim([0 1])

H = gca; H.XAxis.FontSize = 14; H.YAxis.FontSize = 14;
set(L,'fontsize',12)

xlabel('Time (h)','FontSize',20)
ylabel('Roseburia','FontSize',20)


subplot(3,2,2)

plot(t,v(:,13:17),'LineWidth',2)
L=legend('s1=threhalose','s2=pyruvate','s3=glucose','s4=glutamate','s5=lactate');

%xlim([0 tend])
%ylim([0 1])

H = gca; H.XAxis.FontSize = 14; H.YAxis.FontSize = 14;
set(L,'fontsize',12)

%xlabel('Time (h)','FontSize',20)
ylabel('Metabolites','FontSize',20)


subplot(3,2,4)

plot(t,v(:,18:22),'LineWidth',2)
%L=legend('s1','s2','s3','s4','s5','s6','s7','s8','s9','s10');
L=legend('s6=acetate','s7=mannose','s8=succinate','s9=formate','s10=butyrate');

%xlim([0 tend])
%ylim([0 1])

H = gca; H.XAxis.FontSize = 14; H.YAxis.FontSize = 14;

set(L,'fontsize',12)

%xlabel('Time (h)','FontSize',20)
ylabel('Metabolites','FontSize',20)


subplot(3,2,6)

if pH==0
ph=fph(v(:,18),v(:,20),v(:,22));
else
ph=0*t+pH;
end

plot(t,ph,'k','LineWidth',2)

%xlim([0 tend])
%ylim([0 1])

H = gca; H.XAxis.FontSize = 14; H.YAxis.FontSize = 14;
set(L,'fontsize',12)

xlabel('Time (h)','FontSize',20)
ylabel('pH','FontSize',20)


figure(2)
clf


bh=v(:,1)+v(:,2)+v(:,3);  % total living bh
bt=v(:,5)+v(:,6)+v(:,7);  % total living bt
ri=v(:,9)+v(:,10)+v(:,11);  % total living ri

plot(t,ri,t,bh,t,bt,'linewidth',2)

%hold on
%fill([220 240 240 220 220],[0 0 10 10 0],[0.8 0.8 0.8],'LineStyle','none','facealpha',0.5,'HandleVisibility','off')
%hold on
%plot([100 100],[0 10],'k')

H = gca; H.XAxis.FontSize = 14; H.YAxis.FontSize = 14;

L=legend('Ri','Bh','Bt');
set(L,'fontsize',16)

xlabel('Time (h)','FontSize',20)
ylabel('Bacteria','FontSize',20)

title(sprintf('pH=%g, D=%g',pH,D),'fontsize',16)



% ====================================================================
% Run
% ====================================================================

function [t,x]=run(v0,ttt,D,pH,S0)

trans=ttt(1);
tend=ttt(2);
tstep=ttt(3);

ttrans = [0:tstep:trans];
tspan = [0:tstep:tend];
option = odeset('RelTol',1e-6);

if trans > 0 
    [t x] = ode45(@myequa,ttrans,v0,option,D,pH,S0);
    v0=x(end,:);
end

[t x] = ode45(@myequa,tspan,v0,option,D,pH,S0);


% ====================================================================
% Phi
% ====================================================================

function f = phi(ph,a,b)

f=(b^a/gamma(a))*ph^(a-1)*exp(-b*ph);


% ====================================================================
% Ph
% ====================================================================

function ph = fph(s6,s8,s10)

ph=-0.0491*s6-0.0972*s8-0.0436*s10+6.6369;

ph=max(3,ph);
ph=min(10,ph);


% ====================================================================
% Equations
% ====================================================================

function dv = myequa(t,vv,D,ph,S0)

%%% variables

xa=vv(1);
xb=vv(2);
xc=vv(3);
xd=vv(4);
xe=vv(5);
xf=vv(6);
xg=vv(7);
xh=vv(8);
xi=vv(9);
xj=vv(10);
xk=vv(11);
xl=vv(12);

s1=vv(13);
s2=vv(14);
s3=vv(15);
s4=vv(16);
s5=vv(17);
s6=vv(18);
s7=vv(19);
s8=vv(20);
s9=vv(21);
s10=vv(22);

%%% parameters

if ph==0;  
ph=fph(s6,s8,s10);     
end

% ph=ph-0.1*(t>100);   % pH perturbation

% D=D*(t>12);              % start dilution after a certain time
% D=D-D*(t>220)*(t<240);   % temporarily stop dilution

mumaxxa=0.19279348;     % xa_mumax	
mumaxxb=0.9781041;      % xb_mumax	
pxa=6.8180188;          % xa_pHopt
pxb=6.539596;           % xb_pHopt
axa=46.2739563;         % xa_pHalpha	
axb=62.518901;          % xb_pHalpha	
Kxas1=0.287608;         % xa_k_s1	
Kxas2=2.391651;         % xa_k_s2	 
Kxbs3=0.1039681;        % xb_k_s3	
Kxbs4=1.2700485;        % xb_k_s4	
Kxbs2=0.5;              % xb_k_s2	
gxas1=1.2417776;        % xa_g_s1	
gxas2=1.0;              % xa_g_s2	
gxas6s1=1.60E-10;       % xa_g_s6_s1 
gxas5s1=3;              % xa_g_s5_s1
gxas6s2=1.86995;        % xa_g_s6_s2
gxbs3=3.76388;          % xb_g_s3	  
gxbs4=0.5;              % xb_g_s4   
gxbs2=10.0;             % xb_g_s2  
gxbs6s3s4=3;            % xb_g_s6_s3_s4
gxbs6s2=5;              % xb_g_s6_s2	
r1=0.0250945;           % z1_r	
r2=1.5;                 % z2_r  
r3=0.00001;             % z3_r  
r4=0.29603745;          % z4_r	
r5=0.03555697;          % z5_r	
l1s1=4.57E-08;          % z1_l_s1	 
l2s1=0.21;              % z2_l_s1  
l4s3s4=0.1E-3;          % z4_l_s3_s4
l3s3s4=0.000129335;     % z4_l_s3_s4	

h1s1=1.0;               % z1_h_s1	  
h2s1=50;                % z2_h_s1  
h4s3s4=2.9055457;       % z4_h_s3_s4	
mumaxxe=0.9215459;      % xe_mumax	
mumaxxf=1.1965073;      % xf_mumax	
pxe=7.507812;           % xe_pHopt	
pxf=6.984095;           % xf_pHopt
axe=59.82343;           % xe_pHalpha
axf=71.8164;            % xf_pHalpha	
Kxes3=0.35302;          % xe_k_s3	 
Kxes2=10.0;             % xe_k_s2	 
Kxfs7=0.10598211;       % xf_k_s7	   
gxes2=2.9010309;        % xe_g_s2	   
gxes3=1.0269754;        % xe_g_s3	   
gxes5s2=1.71819;        % xe_g_s5_s2
gxes6s2=0.260295;       % xe_g_s6_s2
gxes6s3=0.9966;         % xe_g_s6_s3
gxes8s3=0.62527;        % xe_g_s8_s3
gxes9s2=0.61589;        % xe_g_s9_s2
gxfs7=0.4040234;        % xf_g_s7	   
gxfs6s7=0.853242;       % xf_g_s6_s7
gxfs8s7=2.025524;       % xf_g_s8_s7
r6=1.4959078;           % z6_r 
r7=0.9575532;           % z7_r 
r8=0.0765445;           % z8_r 
r9=0.0044568;           % z9_r 
r10=0.0001;             % z10_r
l6s3=0.009919;          % z6_l_s3	 
l6s7=0.509946;          % z6_l_s7	
l7s3=0.0065072;         % z7_l_s3	
l7ph=5.5;               % z7_l_pH 

l8s7=7.162E-07;         % z8_l_s7	 
l8ph=5.5;               % z8_l_pH 
l10s3=0.5;              % z10_l_s3  
h6s3=1.0006644;         % z6_h_s3	
h6s7=1.00221;           % z6_h_s7	 
h7s3=1.459916;          % z7_h_s3	 
h7ph=10;                % z7_h_pH  
h8s7=30;                % z8_h_s7	  
h8ph=10;                % z8_h_pH  
h10s3=10;               % z10_h_s3
mumaxxi=0.705967;       % xi_mumax
mumaxxj=0.0153452;      % xj_mumax
pxi=7.768276;           % xi_pHopt
pxj=7.6243042;          % xj_pHopt
axi=39.525597;          % xi_pHalpha
axj=20;                 % xj_pHalpha
Kxis2=0.547523;         % xi_k_s2
Kxis3=6.8653727;        % xi_k_s3	
Kxjs5=9.999458;         % xj_k_s5
Kxjs6=0.2756123;        % xj_k_s6
gxis2=2.78375132;       % xi_g_s2
gxis3=1.7594803;        % xi_g_s3
gxjs5=1.9604511;        % xj_g_s5
gxjs6=0.853387;         % xj_g_s6
gxis5s3=0.377444;       % xi_g_s5_s3 
gxis6s2=0.305128;       % xi_g_s6_s2
gxis6s3=0.727664;       % xi_g_s6_s3
gxis10s2=1.994832;      % xi_g_s10_s2
gxis10s3=6.865E-12;     % xi_g_s10_s3
gxjs10s5=2.988196;      % xj_g_s10_s5
gxjs10s6=1.994E-09;     % xj_g_s10_s6
r11=0.0219468;          % z11_r

r12=2.0;                % z12_r
r13=0.01;               % z13_r
r14=0.00001;            % z14_r
r15=0.1;                % z15_r
l11s5=2.331191;         % z11_l_s5
l12s3s2=0.003491;       % z12_l_s3_s2
l15s3s2=8;              % z15_l_s3_s2
h11s5s6=30;             % z11_h_s5
h12s3s2=1.11238;        % z12_h_s3_s2
h15s3s2=50;             % z15_h_s3_s2



%%% equations

Exa=phi(ph,axa,(axa-1)/pxa)/phi(pxa,axa,(axa-1)/pxa);
Exb=phi(ph,axb,(axb-1)/pxb)/phi(pxb,axb,(axb-1)/pxb);
Exe=phi(ph,axe,(axe-1)/pxe)/phi(pxe,axe,(axe-1)/pxe);
Exf=phi(ph,axf,(axf-1)/pxf)/phi(pxf,axf,(axf-1)/pxf);
Exi=phi(ph,axi,(axi-1)/pxi)/phi(pxi,axi,(axi-1)/pxi);
Exj=phi(ph,axj,(axj-1)/pxj)/phi(pxj,axj,(axj-1)/pxj);

Z1=(l1s1^h1s1/(s1^h1s1+l1s1^h1s1))*r1;
Z2=(s1^h2s1/(s1^h2s1+l2s1^h2s1))*r2;
Z3=r3;
Z4=(l4s3s4^h4s3s4/((min(s3,s4))^h4s3s4+l4s3s4^h4s3s4))*r4;
Z5=r5;
Z6=(l6s3^h6s3/(s3^h6s3+l6s3^h6s3))*(s7^h6s7/(s7^h6s7+l6s7^h6s7))*r6;
Z7=(l7s3^h7s3/(s3^h7s3+l7s3^h7s3))*(l7ph^h7ph/(l7ph^h7ph+ph^h7ph))*r7;
Z8=(l8s7^h8s7/(s7^h8s7+l8s7^h8s7))*(l8ph^h8ph/(l8ph^h8ph+ph^h8ph))*r8;
Z9=r9;
Z10=(s3^h10s3/(s3^h10s3+l10s3^h10s3))*r10;
Z11=((s5+s6)^h11s5s6/((s5+s6)^h11s5s6+l11s5^h11s5s6))*r11;
Z12=(l12s3s2^h12s3s2/((s3+s2)^h12s3s2+l12s3s2^h12s3s2))*r12;
Z13=r13;
Z14=r14;
Z15=((s3+s2)^h15s3s2/((s3+s2)^h15s3s2+l15s3s2^h15s3s2))*r15;

dxa=xa*(Exa*mumaxxa*(s1/(s1+Kxas1)+s2/(s2+Kxas2))-(Z1+Z3))+xb*Z2-D*xa;
dxb=xb*(Exb*mumaxxb*(s3/(s3+Kxbs3)*s4/(s4+Kxbs4)+s2/(s2+Kxbs2))-(Z2+Z4))+xa*Z1-D*xb;
dxc=xa*Z3+xb*Z4-xc*Z5-D*xc;
dxd=0; % xc*Z5-D*xd;
dxe=xe*(Exe*mumaxxe*(s2/(s2+Kxes2)+s3/(s3+Kxes3))-(Z6+Z7))+xf*Z10-D*xe;
dxf=xf*(Exf*mumaxxf*(s7/(s7+Kxfs7))-(Z8+Z10))+xe*Z6-D*xf;
dxg=xe*Z7+xf*Z8-xg*Z9-D*xg;
dxh=0; % xg*Z9-D*xh;
dxi=xi*(Exi*mumaxxi*(s3/(s3+Kxis3)+s2/(s2+Kxis2))-(Z11+Z12))+xj*Z15-D*xi;
dxj=xj*(Exj*mumaxxj*(s5/(s5+Kxjs5)+s6/(s6+Kxjs6))-(Z15+Z13))+xi*Z11-D*xj;
dxk=xj*Z13-xk*Z14-D*xk;
dxl=0; % xk*Z14-D*xl;

ds1=-gxas1*(s1/(s1+Kxas1))*Exa*mumaxxa*xa...
    -D*s1+D*S0(1);

ds2=-gxas2*(s2/(s2+Kxas2))*Exa*mumaxxa*xa ...
    -gxbs2*(s2/(s2+Kxbs2))*Exb*mumaxxb*xb ...
    -gxes2*(s2/(s2+Kxes2))*Exe*mumaxxe*xe ...
    -gxis2*(s2/(s2+Kxis2))*Exi*mumaxxi*xi ...
    -D*s2+D*S0(2);

ds3=-gxbs3*(s3/(s3+Kxbs3))*(s4/(s4+Kxbs4))*Exb*mumaxxb*xb ...
    -gxes3*(s3/(s3+Kxes3))*Exe*mumaxxe*xe ...
    -gxis3*(s3/(s3+Kxis3))*Exi*mumaxxi*xi ...
    -D*s3+D*S0(3);

ds4=-gxbs4*(s3/(s3+Kxbs3))*(s4/(s4+Kxbs4))*Exb*mumaxxb*xb ...
    -D*s4+D*S0(4);

ds5=-gxjs5*(s5/(s5+Kxjs5))*Exj*mumaxxj*xj ...
    +gxas5s1*(s1/(s1+Kxas1))*Exa*mumaxxa*xa...
    +gxes5s2*(s2/(s2+Kxes2))*Exe*mumaxxe*xe ...
    +gxis5s3*(s3/(s3+Kxis3))*Exi*mumaxxi*xi ... 
    -D*s5+D*S0(5);

ds6=-gxjs6*(s6/(s6+Kxjs6))*Exj*mumaxxj*xj...
    +((gxas6s1*s1/(s1+Kxas1))+(gxas6s2*s2/(s2+Kxas2)))*Exa*mumaxxa*xa...
    +(gxbs6s3s4*(s3/(s3+Kxbs3))*(s4/(s4+Kxbs4))+gxbs6s2*(s2/(s2+Kxbs2)))*Exb*mumaxxb*xb...  
    +((gxes6s2*s2/(s2+Kxes2))+(gxes6s3*s3/(s3+Kxes3)))*Exe*mumaxxe*xe...
    +(gxfs6s7*s7/(s7+Kxfs7))*Exf*mumaxxf*xf...
    +((gxis6s2*(s2/(s2+Kxis2)))+(gxis6s3*(s3/(s3+Kxis3))))*Exi*mumaxxi*xi ...
    -D*s6+D*S0(6);

ds7=-gxfs7*(s7/(s7+Kxfs7))*Exf*mumaxxf*xf ...
    -D*s7+D*S0(7);

ds8=gxes8s3*(s3/(s3+Kxes3))*Exe*mumaxxe*xe...
    +(gxfs8s7*s7/(s7+Kxfs7))*Exf*mumaxxf*xf ...
    -D*s8+D*S0(8);

ds9=gxes9s2*(s2/(s2+Kxes2))*Exe*mumaxxe*xe ...
    -D*s9+D*S0(9);

ds10=((gxis10s2*(s2/(s2+Kxis2)))+(gxis10s3*(s3/(s3+Kxis3))))*Exi*mumaxxi*xi...
     +((gxjs10s5*(s5/(s5+Kxjs5)))+(gxjs10s6*(s6/(s6+Kxjs6))))*Exj*mumaxxj*xj ...
     -D*s10+D*S0(10);

dv=10*[dxa;dxb;dxc;dxd;dxe;dxf;dxg;dxh;dxi;dxj;dxk;dxl;ds1;ds2;ds3;ds4;ds5;ds6;ds7;ds8;ds9;ds10];



