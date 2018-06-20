function TransMod()

%Transducer impulse response modeling code for transducer
%with P and S crystals, backing, single matching layer and fron load
%All inputs via user interface

% Written by Kevin Wolf, June 2009

%Get stcking order from user
ButtonName = questdlg('Transducer Stacking Order:', ...
    'Select transducer stack order',...
    'Backing  |  P-Crystal  |  S-Crystal  |  Matching Layer  |  Load',...
    'Backing  |  S-Crystal  |  P-Crystal  |  Matching Layer  |  Load',...
    'Backing  |  P-Crystal  |  S-Crystal  |  Matching Layer  |  Load');

switch ButtonName,
    case 'Backing  |  P-Crystal  |  S-Crystal  |  Matching Layer  |  Load',
        order = 1;
    case 'Backing  |  S-Crystal  |  P-Crystal  |  Matching Layer  |  Load',
        order = 2;
end

%Get crytal properties from user
prmt={'Crystal Area (m^2)', ...
    'Crystal Characteristic P Impedance (Rayl)', ...
    'Crystal P-Wave Resonant Frequency (Hz)', ...
    'Crystal P-Wave Velocity (m/s)', ...
    'Crystal P-Wave Quality Factor', ...
    'Crystal P-Wave Dielectic Constant', ...
    'Crystal P-Wave Coupling Coefficient', ...
    'Crystal Characteristic S Impedance (Rayl)', ...
    'Crystal S-Wave Resonant Frequency (Hz)', ...
    'Crystal S-Wave Velocity (m/s)', ...
    'Crystal S-Wave Quality Factor', ...
    'Crystal S-Wave Dielectic Constant', ...
    'Crystal S-Wave Coupling Coefficient'};
    
defans={2.32e-4,19e6,0.25e6,3050,15,270,0.3,12e6,0.25e6,1940,15,240,0.28};
for k=1:length(defans),defans{k}=num2str(defans{k});end;

xtalprop=inputdlg(prmt,'Crystal Properties',1,defans);

for k=1:length(xtalprop),
    xtalnum(k)=str2num(xtalprop{k});
end;
area=xtalnum(1);
Z_c_char_P=xtalnum(2);
fo_P=xtalnum(3);
V_c_P=xtalnum(4);
Q_c_P=xtalnum(5);
RelEps_P=xtalnum(6);
k_t=xtalnum(7); 
Z_c_char_S=xtalnum(8);
fo_S=xtalnum(9);
V_c_S=xtalnum(10);
Q_c_S=xtalnum(11);
RelEps_S=xtalnum(12);
k_15=xtalnum(13);

%Get backing layer properties from user
prmt={'Backing Layer Characteristic P Impedance (Rayl)', ...
    'Backing Layer Characteristic S Impedance (Rayl)'};
    
defans={19e6,12e6};
for k=1:length(defans),defans{k}=num2str(defans{k});end;

blprop=inputdlg(prmt,'Backing Layer Properties',1,defans);

for k=1:length(blprop),
    blnum(k)=str2num(blprop{k});
end;
Z_bl_char_P=blnum(1);
Z_bl_char_S=blnum(2);

%Get matching layer properties from user
prmt={'Matching Layer Thickness (m)' ...
    'Matching Layer Characteristic P Impedance (Rayl)', ...
    'Matching Layer P-Wave Velocity (m/s)', ...
    'Matching Layer P-Wave Quality Factor', ...
    'Matching Layer Characteristic S Impedance (Rayl)', ...
    'Matching Layer S-Wave Velocity (m/s)', ...
    'Matching Layer S-Wave Quality Factor'};

defans={0.00254,6.15e6,3820,10,4.34e6,2700,10};
for k=1:length(defans),defans{k}=num2str(defans{k});end;

mlprop=inputdlg(prmt,'Matching Layer Properties',1,defans);

for k=1:length(mlprop),
    mlnum(k)=str2num(mlprop{k});
end;
t_ml=mlnum(1);
Z_ml_char_P=mlnum(2);
V_ml_P=mlnum(3);
Q_ml_P=mlnum(4);
Z_ml_char_S=mlnum(5);
V_ml_S=mlnum(6);
Q_ml_S=mlnum(7);

%Get front load properties from user
prmt={'Front Load Characteristic P Impedance (Rayl)', ...
    'Front Load Characteristic S Impedance (Rayl)'};
    
defans={5.94e6,2.2e6};
for k=1:length(defans),defans{k}=num2str(defans{k});end;

flprop=inputdlg(prmt,'Front Load Properties',1,defans);

for k=1:length(flprop),
    flnum(k)=str2num(flprop{k});
end;
Z_fl_char_P=flnum(1);
Z_fl_char_S=flnum(2);

%Get pulse generator properties from user
prmt={'Pulse Generator Volatage (V)', ...
    'Pulse Generator Resistance (ohms)'};
    
defans={500,50};
for k=1:length(defans),defans{k}=num2str(defans{k});end;

pgprop=inputdlg(prmt,'Pulser Properties',1,defans);

for k=1:length(pgprop),
    pgnum(k)=str2num(pgprop{k});
end;
V_s=pgnum(1);
res=pgnum(2);



%General
j = sqrt(-1);
wo_P = 2*pi*fo_P;           %P resonant frequency of crystal (rad)
f_P = (1:1000:2*fo_P)';    %P frequency range (Hz)
w_P = 2*pi*f_P;             %P frequency range(rad)
wo_S = 2*pi*fo_S;           %S resonant frequency of crystal (rad)
f_S = (1:1000:2*fo_S)';    %S frequency range (Hz)
w_S = 2*pi*f_S;             %S frequency range(rad)

%Crystal Properties
Z_c_P = Z_c_char_P*area;
Z_c_S = Z_c_char_S*area;
d_P = V_c_P/(2*fo_P);     %crystal thickness for frequency given above (m)
d_S = V_c_S/(2*fo_S);     %crystal thickness for frequency given above (m)
E0 = 8.854e-12;         %free space permittivity
Eps_P = RelEps_P*E0;        %crystal permittivity (F/m)
Eps_S = RelEps_S*E0;        %crystal permittivity (F/m)
Co_P = Eps_P*area/d_P;           %clamped capacitance
Co_S = Eps_S*area/d_S;           %clamped capacitance

Cprime_P = -Co_P./(k_t^2*sinc(w_P/wo_P));
phi1_P = k_t/sqrt(2*fo_P*Co_P*Z_c_P);
phi2_P = sinc(f_P/(2*fo_P));
phi_P = phi1_P.*phi2_P;
Cprime_S = -Co_S./(k_15^2*sinc(w_S/wo_S));
phi1_S = k_15/sqrt(2*fo_S*Co_S*Z_c_S);
phi2_S = sinc(f_S/(2*fo_S));
phi_S = phi1_S.*phi2_S;

%Backing Load
Z_bl_P = Z_bl_char_P*area;
Z_bl_S = Z_bl_char_S*area;

%Matching Layer
Z_ml_P = Z_ml_char_P*area;
Z_ml_S = Z_ml_char_S*area;

%Front Load
Z_fl_P = Z_fl_char_P*area;
Z_fl_S = Z_fl_char_S*area;

%Calculate propagation constants
Gam_c_P = (j*2*pi*f_P/V_c_P)*(1-j/(2*Q_c_P));
Gam_c_S = (j*2*pi*f_S/V_c_S)*(1-j/(2*Q_c_S));
arg1_P = Gam_c_P*d_P/2;
arg1_S = Gam_c_S*d_S/2;

Gam_ml_P = (j*2*pi*f_P/V_ml_P)*(1-j/(2*Q_ml_P));
Gam_ml_S = (j*2*pi*f_S/V_ml_S)*(1-j/(2*Q_ml_S));
arg2_P = Gam_ml_P.*t_ml;
arg2_S = Gam_ml_S.*t_ml;

%Calculate ABCD matrices
if order==1
    [A_P,B_P,C_P,D_P]=TransConfig1(f_P,Gam_c_P,Gam_ml_P,Cprime_P,...
        Co_P,phi_P,d_P,d_S,t_ml,Z_c_P,Z_bl_P,Z_ml_P);
    [A_S,B_S,C_S,D_S]=TransConfig2(f_S,Gam_c_S,Gam_ml_S,Cprime_S,...
        Co_S,phi_S,d_S,d_P,t_ml,Z_c_S,Z_bl_S,Z_ml_S);
elseif order==2
    [A_P,B_P,C_P,D_P]=TransConfig2(f_P,Gam_c_P,Gam_ml_P,Cprime_P,...
        Co_P,phi_P,d_P,d_S,t_ml,Z_c_P,Z_bl_P,Z_ml_P);
    [A_S,B_S,C_S,D_S]=TransConfig1(f_S,Gam_c_S,Gam_ml_S,Cprime_S,...
        Co_S,phi_S,d_S,d_P,t_ml,Z_c_S,Z_bl_S,Z_ml_S);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Transducer Properties for P wave %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate transducer electrical impedance
Z_e_P = (A_P*Z_fl_P+B_P)./(C_P*Z_fl_P+D_P);
%Calculate velocity at front transducer
denom_P = res*(C_P*Z_fl_P+D_P)+A_P*Z_fl_P+B_P;
v_P = -V_s./denom_P;
%Displacement at front transducer
u_P = v_P./(j*w_P);
%Force applied to front load
denom2_P = (A_P*Z_fl_P+B_P).*(res+Z_e_P);
F_P = V_s.*Z_fl_P.*Z_e_P./denom2_P;
%Power delivered to load
P_P = 0.5*Z_fl_P*v_P.^2;
P_max = V_s^2/(8*res);
%Insertion Loss
IL_P = 10*log10(P_max./P_P);
IL2_P = -10*log10((4*res*Z_fl_P)/((denom_P).^2));
t1_P = 4*res*Z_fl_P*Z_e_P.^2;
t2_P = (res+Z_e_P).^2;
V_P = V_s*(Z_e_P./(res+Z_e_P));
t3_P = (v_P./V_P).^2;
IL3_P = -10*log10(t1_P./t2_P.*t3_P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate P Impulse Response %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pad = 10000;
%Calculate current and phase from current
I_P = V_P./Z_e_P;
psiI_P = atan(imag(I_P)./(real(I_P)));

psiI_P_pad = [psiI_P' zeros(1,pad)];
P_P_pad = [abs(P_P') zeros(1,pad)];

P_psiI_P =P_P_pad.*exp(j*psiI_P_pad);
P_allI_P = [P_psiI_P conj(fliplr(P_psiI_P(2:end)))];
IR_P = ifft(P_allI_P);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Transducer Properties for S wave %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate transducer electrical impedance
Z_e_S = (A_S*Z_fl_S+B_S)./(C_S*Z_fl_S+D_S);
%Calculate velocity at front transducer
denom_S = res*(C_S*Z_fl_S+D_S)+A_S*Z_fl_S+B_S;
v_S = -V_s./denom_S;
%Displacement at front transducer
u_S = v_S./(j*w_S);
%Force applied to front load
denom2_S = (A_S*Z_fl_S+B_S).*(res+Z_e_S);
F_S = V_s.*Z_fl_S.*Z_e_S./denom2_S;
%Power delivered to load
P_S = 0.5*Z_fl_S*v_S.^2;
P_max = V_s^2/(8*res);
%Insertion Loss
IL_S = 10*log10(P_max./P_S);
IL2_S = -10*log10((4*res*Z_fl_S)/((denom_S).^2));
t1_S = 4*res*Z_fl_S*Z_e_S.^2;
t2_S = (res+Z_e_S).^2;
V_S = V_s*(Z_e_S./(res+Z_e_S));
t3_S = (v_S./V_S).^2;
IL3_S = -10*log10(t1_S./t2_S.*t3_S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate S Impulse Response %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate current and phase from current
I_S = V_S./Z_e_S;
psiI_S = atan(imag(I_S)./(real(I_S)));

psiI_S_pad = [psiI_S' zeros(1,pad)];
P_S_pad = [abs(P_S') zeros(1,pad)];

P_psiI_S =P_S_pad.*exp(j*psiI_S_pad);
P_allI_S = [P_psiI_S conj(fliplr(P_psiI_S(2:end)))];
IR_S = ifft(P_allI_S);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT IMPLUSE RESPONSES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt_P = (1/1000)/length(P_allI_P);
temp_P = ones(length(P_allI_P),1);temp_P(1)=0;
t_P = cumsum(temp_P)*dt_P*1e6;
dt_S = (1/1000)/length(P_allI_S);
temp_S = ones(length(P_allI_S),1);temp_S(1)=0;
t_S = cumsum(temp_S)*dt_S*1e6;
figure;
subplot(2,1,1),plot(t_P,fftshift(real(IR_P)));
xlabel('Time (\mus) (relative time important, not absolute)');
ylabel('Volts (V)');
subplot(2,1,2),plot(t_S,fftshift(real(IR_S)));
xlabel('Time (\mus) (relative time important, not absolute)');
ylabel('Volts (V)');
