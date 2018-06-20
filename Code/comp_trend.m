function [phi_sd vp_sd vs_sd phi_sh vp_sh vs_sh z_m]=comp_trend(z_m, plot_ind)
% [phi_sd vp_sd vs_sd phi_sh vp_sh vs_sh z_m]=comp_trend(z_m, plot_ind)
% 
% provides depth trends of porosity and velocity for shale and clean brine
% sands at shallow section (< 2000 m).
% Example
% [phi_sd vp_sd vs_sd phi_sh vp_sh vs_sh]=comp_trend([0:1:1800],1);
% Input 
% z_m:      depth (vector) in "meters"
% plot_ind: '1' for plotting, '0' otherwise
% Outputs
% phi_sd, vp_sd, vs_sd: porosity, Vp, Vs depth trends (vectors) in sand
% phi_sh, vp_sh, vs_sh: porosity, Vp, Vs depth trends in shale
%                       output velocities in m/s

% Written by Tanima Dutta, 2009


% Reference: Dutta, T., Mukerji, T., Mavko, G., Lane, T., 2009, Compaction
% trends for shale and clean sandstone in shallow sediments, 
% Gulf of Mexico:  The Leading Edge, 28, 260-266.
%
% rock property =a*exp(b*z)+c*exp(d*z); %General depth trend equation

z=z_m.*3.2808399; %Depth in ft

% coefficients
a_phi_sd= 0.01132; b_phi_sd= -0.004874; c_phi_sd= 0.3923; d_phi_sd=-1.678e-005;
a_phi_sh=0.2875 ; b_phi_sh=-0.00774; c_phi_sh=0.4384 ; d_phi_sh=-0.0001761;

a_vp_sd=6350 ; b_vp_sd= 2.532e-005; c_vp_sd= -500.3; d_vp_sd=-0.003647;
a_vs_sd=1951 ; b_vs_sd= 7.962e-005; c_vs_sd=-1788 ; d_vs_sd=-0.001942;

a_vp_sh=6917 ; b_vp_sh=4.633e-005 ; c_vp_sh=-1652 ; d_vp_sh=-0.0003646;
a_vs_sh=3540 ; b_vs_sh= 1.2e-005; c_vs_sh= -3536; d_vs_sh=-0.000421;

% porosity depth trends

phi_sd=a_phi_sd*exp(b_phi_sd*z)+c_phi_sd*exp(d_phi_sd*z);
phi_sh=a_phi_sh*exp(b_phi_sh*z)+c_phi_sh*exp(d_phi_sh*z);

% velocity depth trends

vp_sd=a_vp_sd*exp(b_vp_sd*z)+c_vp_sd*exp(d_vp_sd*z); vp_sd=0.3048*vp_sd; %Convert Vp, Vs from ft/s to m/s
vs_sd=a_vs_sd*exp(b_vs_sd*z)+c_vs_sd*exp(d_vs_sd*z); vs_sd=0.3048*vs_sd;

vp_sh=a_vp_sh*exp(b_vp_sh*z)+c_vp_sh*exp(d_vp_sh*z); vp_sh=0.3048*vp_sh;
vs_sh=a_vs_sh*exp(b_vs_sh*z)+c_vs_sh*exp(d_vs_sh*z); vs_sh=0.3048*vs_sh;

% plotting
if plot_ind==1
    figure, subplot(1,3,1), plot(phi_sd, z_m, phi_sh, z_m), axis ij, xlabel('Porosity'), ylabel('Depth (m)'), legend('Sand', 'Shale')
    subplot(1,3,2), plot(vp_sd, z_m, 'b', vs_sd, z_m, 'b:', vp_sh, z_m, 'r', vs_sh, z_m, 'r:'), axis ij, xlabel('Velocity (m/s)'), ylabel('Depth (m)'), legend('Sand Vp', 'Sand Vs', 'Shale Vp', 'Shale Vs');
    subplot(1,3,3),  plot(vs_sd./vp_sd, z_m, 'b', vs_sh./vp_sh, z_m, 'r'), axis ij, xlabel('Vs/Vp ratio'), ylabel('Depth (m)'), legend('Sand Vp/Vs', 'Shale Vp/Vs');
end
    