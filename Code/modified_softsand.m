function [Vp,Vs,Ip,PR,RHOB,Phi]=modified_softsand(Kf,RHOf,Phi,Quartz,Clay,Feldspar,Limestone,Dolomite,Pressure,PhiC)
% [Vp,Vs,Ip,PR,RHOB,Phi]=modified_softsand(Kf,RHOf,Phi,Quartz,Clay,Feldspar,Limestone,Dolomite,Pressure,PhiC)
%
% Modified softsand model for unconsolidated sediments.
%
%               INPUTS
%Kf             Fluid bulk modulus
%RHOf			Fluid density
%Phi			Porosity (fraction)
%Quartz, etc.	Volume mineral content in solid phase (fraction)
%Pressure		Effective pressure (MPa)
%PhiC			Critical porosity ~0.4
%               OUTPUTS
%Vp,Vs          Effective P and S-wave velocities
%Ip, PR         P-impedance and Poisson's ratio
%RHOB, Phi		Effective density and porosity
%
% Example
% [Vp,Vs,Ip,PR,RHOB,Phi]=modified_softsand(2.25,1,[0.01:0.01:0.4]',1,0,0,0,0,20,0.4);
% 
% see also SoftSand, waltonv

% At the critical porosity effective moduli are computed by 
% mixing rough and smooth grains of Walton (1987), combined with
% empirical coordination number vs. pressure relations inverted from 
% Zimmer's measurements. 
% The effective moduli between critical porosity and mineral endmember
% are interpolated using the lower Hashin-Shtrikman bound (same principle as
% in Dvorkin's softsand model
%
% written by Tanima Dutta, 2009

%Step-1A
%Balancing mineralogy
Dolomite=1-(Quartz+Clay+Feldspar+Limestone);
%Solid-phase elastic moduli
KsV=Quartz*36.6+Clay*21+Limestone*76.8+Dolomite*94.9+Feldspar.*75.6;
KsR=1./(Quartz/36.6+Clay/21+Limestone/76.8+Dolomite/94.9+Feldspar./75.6);
Ks=0.5*(KsV+KsR);
GsV=Quartz*45+Clay*7+Limestone*32+Dolomite*45+Feldspar.*25.6;
GsR=1./(Quartz/45+Clay/7+Limestone/32+Dolomite/45+Feldspar./25.6);
Gs=0.5*(GsV+GsR);
Ms=Ks+(4/3)*Gs;
NUs=0.5*(Ms./Gs-2)./(Ms./Gs-1);
%Solid-phase density
RHOs=Quartz*2.65+Clay*2.58+Limestone*2.71+Dolomite*2.87+Feldspar.*2.63;

%Step-1B: get C from C-pressure relation (This relation was obtained by inverting
%Zimmer's data using the combined Walton model for 60% rough and 40% smooth grains)
    pr=Pressure;
    
    a= 777.1;b= 0.001545;c=-770.7;
    Cp_pr= a*pr.^b+c;
         
    a= -4.457;b=-0.2724 ;c=9.401;
    Cs_pr= a*pr.^b+c;
    
   % (a,b,c are constants,input pressure in MPa, 
   % Cp_pr= coordination # for P-wave,Cs_pr= coordination # for S-wave 
   
%Step-2: run walton model using Cp to get Vp_dry
   
    %mineral moduli 
    Kmin= Ks*10^9;Gmin=Gs*10^9; Rhomin= RHOs;
    
    % get the critical porosity and compute density 
    porY=PhiC;
    den=(1-porY).*Rhomin * 10^3;

    % run walton  model
    [Vp_w_S,Vs_w_S]= waltonv ('s',Kmin,Gmin,Rhomin*10^3,pr*10^6,porY,Cp_pr);
    [K_WS,G_WS]=v2ku(Vp_w_S,Vs_w_S,den);
    lambda_WS= K_WS - (2/3)*G_WS;

    % use fractional value of alpha in walton model
    % alpha= fraction of rough contacts
    
    alpha=0.6;

    G_WF= (1+ 3*alpha* ((1-.06)./(2-.06))).* G_WS;
    lambda_WF = (1- 2*alpha* ((1-.06)./(2-.06))).* lambda_WS;
    
    K_WF = lambda_WF + (2/3)*G_WF;
    [Vp_w_F,Vs_w_F]= ku2v (K_WF,G_WF,den);
    Vp_dry=Vp_w_F;
    
 %Step-3: run walton model again using Cs to get Vs_dry
     
 
    [Vp_w_S,Vs_w_S]= waltonv ('s',Kmin,Gmin,Rhomin*10^3,pr*10^6,porY,Cs_pr);
    [K_WS,G_WS]=v2ku(Vp_w_S,Vs_w_S,den);
    lambda_WS= K_WS - (2/3)*G_WS;
    
    G_WF= (1+ 3*alpha* ((1-.06)./(2-.06))).* G_WS;
    lambda_WF = (1- 2*alpha* ((1-.06)./(2-.06))).* lambda_WS;
    
    
    K_WF = lambda_WF + (2/3)*G_WF;
    [Vp_w_F,Vs_w_F]= ku2v (K_WF,G_WF,den);
    Vs_dry=Vs_w_F;
    
    [Khat,Ghat]=v2ku (Vp_dry,Vs_dry,den);
    Khat=Khat/10^9;Ghat=Ghat/10^9;
    
 %Effective bulk and shear moduli at porosity Phi<=PhiC
    KDry1=1./((Phi./PhiC)./(Khat+4.*Ghat./3)+((PhiC-Phi)./PhiC)./(Ks+4.*Ghat./3))-4.*Ghat./3;
    ZZ1=(Ghat./6).*(9.*Khat+8.*Ghat)./(Khat+2.*Ghat);
    GDry1=1./((Phi./PhiC)./(Ghat+ZZ1)+((PhiC-Phi)./PhiC)./(Gs+ZZ1))-ZZ1; 
    MDry1 = KDry1+(4./3).*GDry1;
    NuDry1=0.5.*(MDry1./GDry1-2)./(MDry1./GDry1-1);
%Effective bulk and shear moduli at porosity Phi>PhiC
    KDry2 = 1./(((1-Phi)./(1-PhiC))./(Khat+4.*Ghat./3)+((Phi-PhiC)./(1-PhiC))./(4.*Ghat./3))-4.*Ghat./3;
    ZZ2 = (Ghat./6).*(9.*Khat+8.*Ghat)./(Khat+2.*Ghat);
    GDry2 = 1./(((1-Phi)./(1-PhiC))./(Ghat+ZZ2)+((Phi-PhiC)./(1-PhiC))./(ZZ2))-ZZ2;
    MDry2 = KDry2+(4./3).*GDry2;
    NuDry2=0.5.*(MDry2./GDry2-2)./(MDry2./GDry2-1);
    MDry=MDry1.*(Phi<=PhiC)+MDry2.*(Phi>PhiC);
    GDry=GDry1.*(Phi<=PhiC)+GDry2.*(Phi>PhiC);
    NuDry=NuDry1.*(Phi<=PhiC)+NuDry2.*(Phi>PhiC);
    KDry=MDry-(4/3)*GDry;
% ==============================================
% ================== Saturated Rock ============
    KSat=Ks.*(Phi.*KDry-(1+Phi).*Kf.*KDry./Ks+Kf)./((1-Phi).*Kf+Phi.*Ks-Kf.*KDry./Ks);
    MSat=KSat+(4/3)*GDry;
    RHOB=(1-Phi).*RHOs+Phi.*RHOf;
    Vp=sqrt(MSat./RHOB);
    Vs=sqrt(GDry./RHOB);
    Ip=Vp.*RHOB;
    PR=0.5*(MSat./GDry-2)./(MSat./GDry-1);