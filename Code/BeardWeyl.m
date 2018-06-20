function [sorting,phi_coarse,phi_med,phi_fine,phi_very_fine]= BeardWeyl()
% [sorting,phi_coarse,phi_med,phi_fine,phi_very_fine]= BeardWeyl()
% sorting versus porosity relations in sandstones for four
% different grain sizes using Beard and Weyl (1973)
% Output
% sorting:      sorting coefficient
% phi_coarse:   mean porosity for coarse sands
% phi_coarse:   mean porosity for medium sands
% phi_coarse:   mean porosity for fine sands
% phi_coarse:   mean porosity for very-fine sands
%

%written by Tanima Dutta, 2006

% sorting coefficient  for six different sorting classes
% A higher value indicates poor sorting
class={'Extremely well sorted', 'Very Well sorted','Well Sorted', ...
    'Moderately Sorted','Poorly sorted','Very poorly sorted'};
sorting= [0.072 0.207 0.389 0.787 1.267 2.128];

%Grain size-1: coarse sands
phi_up_coarse= [0.431 0.408 0.380 0.324 0.271 0.286];
phi_low_coarse= [0.428 0.415 0.384 0.333 0.298 0.252];
phi_coarse= 0.5* (phi_up_coarse+phi_low_coarse);

L_coarse= phi_coarse-phi_low_coarse;
U_coarse= phi_up_coarse- phi_coarse;

%Grain size-2: medium sands
phi_up_med= [0.417 0.402 0.381 0.342 0.315 0.258];
phi_low_med= [0.428 0.415 0.384 0.333 0.298 0.252];
phi_med= 0.5* (phi_up_med+phi_low_med);

L_med= phi_med- phi_low_med;
U_med= phi_up_med- phi_med;

%Grain size-3: fine sands
phi_up_fine= [0.413 0.398 0.391 0.339 0.304 0.285];
phi_low_fine=[0.435 0.408 0.397 0.343 0.310 0.290];
phi_fine= 0.5* (phi_up_fine+phi_low_fine);

L_fine= phi_fine- phi_low_fine;
U_fine= phi_up_fine- phi_fine;

%Grain size-4: very fine sands
phi_up_very_fine= [0.423 0.412 0.402 0.356 0.305 0.301];
phi_low_very_fine=[0.430 0.418 0.398 0.331 0.342 0.326];
phi_very_fine= 0.5* (phi_up_very_fine+phi_low_very_fine);

L_very_fine= phi_very_fine- phi_low_very_fine;
U_very_fine= phi_up_very_fine- phi_very_fine;

% Plot sorting versus porosity for different grain sizes
errorbar(sorting,phi_coarse,L_coarse,U_coarse);
hold,
errorbar(sorting,phi_med,L_med,U_med,'g');
errorbar(sorting,phi_fine,L_fine,U_fine,'r');
errorbar(sorting,phi_very_fine,L_very_fine,U_very_fine,'k');

legend ('coarse', 'medium','fine','very fine');
text (sorting,phi_very_fine,class);

xlabel ('Decreasing sorting');
ylabel ('Porosity');
