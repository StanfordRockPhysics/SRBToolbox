function [vp_sat,vs_sat,density_sat]= kataharaflsub(vp1,vs1,vp2,density1,density2,ROFL1,KFL1,ROFL2,KFL2,K0,PHI)
%[vp_sat,vs_sat,density_sat]= kataharaflsub(vp1,vs1,vp2,density1,density2,ROFL1,KFL1,ROFL2,KFL2,K0,PHI)
%Fluid substitution in laminated shaly-sand using Katahara's formulation (2004) 
% Example
% [vp_sat,vs_sat,density_sat]= kataharaflsub(2.3,1,2.8,2.1,2.3,1.0,2.25,0.7,.8,36,.35);
%
% Input
%
% vp1,vs1:      Vp and Vs of sand end-member with initial fluid(km/sec)
% vp2:          Vp of shale end-member(km/sec)
% density1:     Density of sand end-member(g/cc)
% density2:     Density of shale end-member(g/cc)
% ROFL1, KFL1:  Density (g/cc) and bulk modulus (GPa) of initial fluid
% ROFL2, KFL2:  Density (g/cc) and bulk modulus (GPa) of new fluid
% K0, PHI:      Mineral bulk modulus for sandstone(GPa), and sand porosity in fraction
%
% Output
%
% vp_sat,vs_sat,density_sat:    Vp, Vs and density of saturated shaly sand
% see also gassmnv

% written by Tanima Dutta, 2007

% compliance from velocity and density
    sand_point_den=density1;sh_point_den= density2;
    sand_point_compl= 1./(sand_point_den.*vp1.^2);
    sh_point_compl= 1./(sh_point_den.*vp2.^2);

%apply gassmn at sand end-member
    
    [vp_sd_sat,vs_sd_sat,rhob_sd_sat,k_sd_sat]=gassmnv (vp1,vs1,density1,ROFL1,KFL1,ROFL2,KFL2,K0,PHI)   

%join fluid substituted sand end-member with shale end-member
    Pcompliance_sat= linspace (sh_point_compl, (1/(rhob_sd_sat*(vp_sd_sat)^2)));
    Scompliance_sat= linspace (sh_point_compl, (1/(rhob_sd_sat*(vs_sd_sat)^2)));
    density_sat= linspace (sh_point_den, rhob_sd_sat);
    vp_sat=sqrt ( 1./ ((Pcompliance_sat.*density_sat)));
    vs_sat=sqrt ( 1./ ((Scompliance_sat.*density_sat)));

    compliance= linspace (sh_point_compl,sand_point_compl);
    density= linspace (sh_point_den, sand_point_den);
    
% graphical plot of fluid substitution in compliance-density domain
    figure;  
    plot (density,compliance,'k');
    hold on;
    plot(sh_point_den,sh_point_compl,'o','markerfacecolor',[0 .7 0],'markeredgecolor',[0 .7 0]);
    text(sh_point_den+.02,sh_point_compl,'shale point','fontsize',14,'color',[0 .7 0]);
    plot(sand_point_den,sand_point_compl,'o','markerfacecolor',[0 .7 0],'markeredgecolor',[0 .7 0]);
    text(sand_point_den+.02,sand_point_compl,'sand point','fontsize',14,'color',[0 .7 0]);
    
    plot (density_sat,Pcompliance_sat,'b');
    plot(rhob_sd_sat,(1/(rhob_sd_sat*(vp_sd_sat)^2)),'o','markerfacecolor',[0 .7 0],'markeredgecolor',[0 .7 0]);
    text(rhob_sd_sat+.02,(1/(rhob_sd_sat*(vp_sd_sat)^2)),'Fluid-sub sand point','fontsize',14,'color',[0 .7 0]);
    xlabel ('Density, g/cc'); ylabel ('Compliance, 1/GPa');




