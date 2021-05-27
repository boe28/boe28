%% Condenser - Heat pump to Warm side system - Propane (R290)
%Finding the dynamic system that describe the energy transfer.
clc; clear all; close all; format shortEng;

display('Welcome to an interactive presentation of our system!')
display(' ')
display('This program is supplemental to our report: Control and system engineering to')
display('optimise energy transfer from natural refrigerant heat pump systems.')
display('Please see chapter 7 for background and additional information.')
display(' ')
display(' ')
display(' ')
display('The system presented here is based upon a mathematical representation taken from:')
display('State-space model for dynamic behavior of vapor compression liquid chiller')
display('https://www.sciencedirect.com/science/article/abs/pii/S0140700713001187#fig2')

x = input('Press "Enter" to continue!');
clc;
%% variable
%
display('There are a lot of variables in this model! Everything from the diameter of the pipes, area of the')
display('condenser walls, temperature of the substances in different states, for all the inputs and outputs.')
display('There is also enthalpy, mass, specific heat at a specific pressure, void fractions, and mass flowrate.')
display('Some are easy to obtain from heat pump datasheets, however some variables like')
display('inner physical dimensions are not public, as these often are company secrets.')
display('We have proposed the closest estimates to our ability.')
display('The variables can be modified according to your needs. New variables will not “ruin” the system,')
display('they will only change the outcome. However, if one of them has the same value (or very close value)')
display('as an another variable, there is a risk that it might cancel out one or more of the equations. ')
display('Some mathematical knowledge is required.')
display(' ')
x = input('Press "Enter" to continue!');
clc;
%diameter of pipe [m]
d_evaporator = 0.00137;     % Y Table 2 
d = 0.059;                  % Y Table 2
%density [kg m^3]
p2 = 2;                     % N (2 = compressor outlet/condenser inlet)
p_refrigerantVaporC = 1.5;    % N 
p_liquid = 550;               % N hva slags "liquid"?
p_refrigerantLiquid = 400; % Y 
p_heatExchangerWallC = 2700;   % N
p_coolantMedium = 1000;        % N hva er forskjellen fra liquid?
%enthalpy [1/(j kg)]
h1 = 450;                     % Y ? https://www.researchgate.net/figure/Ideal-vapour-compression-refrigeration-cycle-of-propane-on-the-p-h-diagram_fig3_335869667
h2 = 650;                     % Y ?
h_refrigerantVaporC = 95;   % N
h_latentC = 70;              % N
%convective heat transfer coefficient [W/(m^2 C)]
alfa_refrigerant = 650;       
alfa_coolantMedium = 1000;
%area of heat exchanger [m^2]
Ac = 8;                     % N
%length of flow
Lc = 4;                     % N
L_OverheatedRegionC = 0.85;    % N
L_2phaseRegionC = 0.4;        % N
%temperature [C]
t_kondensasjon = 87.4;         % N
t_heatExchangerWall_OverheatedRegion = 91.0;   % N
t_heatExchangerWall_twoPhaseRegionC = 88.2;    % N
t_heatExchangerWall_OverCoolingRegionC = 1; % N
t_compressorOutlet = 92.5;                  % N
t_condenserOutlet = 83.0;                   % N
t_coolantMedium_InletC = 45.0;              % N
t_coolantMedium_twoPhaseRegionC = 50.2;     % N
t_coolantMedium_overheatedRegionC = 55;   % N
t_coolantMedium_OutletC = 60.0;             % N
%d/dt
dp_dt = 1; %p_refrigerantVaporC derivert med tanke på t_k
% void fraction - 0 til 1 "% andel veske til gass"
voidFraction = 0.6;         % N
% specific heat at constant pressure [J / kg C]
c_p_Liquid = 1900;             
c_p_condenserOutlet = 1200;    
c_p_heatExchangerWallC = 900; 
c_p_coolantMedium = 3800; 
% volume [m^3]
V_refrigerantC = 0.8;         % N
V_heatExchangerWallC = 0.08;   % N
V_coolantMediumC = 0.8;       % N
% mass flow rate [kg s^−1]
G_refrigerant_condenserOutlet = 0.2;      
G_refrigerant_CompressorOutlet = 0.2;      
G_coolantMediumC = 3.5;  

%% Inputs and outputs
%
display('We have a MIMO system, with 10 inputs and 3 outputs')
display(' ')
display('This means that we have multiple Transfer functions,')
display('and aspects of interest in variation of all these inputs.')
display(' ')
display('We will now print all the results of the 72 equations that are in the ABC-Matrix.')
display(' ')
x = input('Press "Enter" to continue!');
clc;

deltaX = [L_OverheatedRegionC,t_kondensasjon,L_2phaseRegionC,t_condenserOutlet,t_heatExchangerWall_OverheatedRegion,t_heatExchangerWall_twoPhaseRegionC,t_heatExchangerWall_OverCoolingRegionC,t_coolantMedium_twoPhaseRegionC,t_coolantMedium_overheatedRegionC,t_coolantMedium_OutletC];  
deltaU = [G_refrigerant_CompressorOutlet,G_refrigerant_condenserOutlet,t_condenserOutlet,t_coolantMedium_InletC,G_coolantMediumC];
deltaY = [t_kondensasjon,t_condenserOutlet,t_coolantMedium_OutletC];

%% The elements to be Partial derivated
%
T1 = (1.0/8.0)*pi*(d_evaporator^2.0)*p2*(h2-h_refrigerantVaporC)
X11 = alfa_refrigerant*(Ac/Lc)*(t_heatExchangerWall_OverheatedRegion - ((t_kondensasjon + t_compressorOutlet) / 2.0))
X12 = -alfa_refrigerant*L_OverheatedRegionC*(Ac/(2.0*Lc))
X13 = alfa_refrigerant*L_OverheatedRegionC*(Ac/Lc)
X14 = h2-h_refrigerantVaporC
X15 = -alfa_refrigerant*L_OverheatedRegionC*(Ac/(2.0*Lc))
T2 = (1.0/8.0)*pi*(d_evaporator^2.0)*L_OverheatedRegionC*(dp_dt)
X21 = -alfa_refrigerant*(Ac/Lc)*(t_heatExchangerWall_OverheatedRegion - ((t_kondensasjon + t_compressorOutlet) / 2.0))*((p2-p_refrigerantVaporC)/p2*(h2-h_refrigerantVaporC))
X22 = alfa_refrigerant*L_OverheatedRegionC*(Ac/(2.0*Lc))*(p2-p_refrigerantVaporC)/p2*(h2-h_refrigerantVaporC) - alfa_refrigerant*(Ac/Lc)*(L_2phaseRegionC/h_latentC)
X23 = - alfa_refrigerant*(Ac/Lc)*(1/h_latentC)*(t_kondensasjon - t_heatExchangerWall_twoPhaseRegionC)
X24 = -alfa_refrigerant*L_OverheatedRegionC*(Ac/Lc)*(p2-p_refrigerantVaporC)/p2*(h2-h_refrigerantVaporC)
X25 = alfa_refrigerant*(Ac/Lc)*(L_2phaseRegionC/h_latentC)
X26 = p_refrigerantVaporC/p2
X27 = alfa_refrigerant*L_OverheatedRegionC*(Ac/Lc)*(p2-p_refrigerantVaporC)/p2*(h2-h_refrigerantVaporC)
T3 = (pi/4)*(d_evaporator^(2.0))*p_liquid*h_latentC*(1-voidFraction)
X31 = alfa_refrigerant*(Ac/Lc)*(L_2phaseRegionC)
X32 = alfa_refrigerant*(Ac/Lc)*(t_kondensasjon - t_heatExchangerWall_twoPhaseRegionC)
X33 = -alfa_refrigerant*(Ac/Lc)*(L_2phaseRegionC)
X34 = -h_latentC;
T4 = p_refrigerantLiquid*c_p_Liquid*V_refrigerantC*((Lc-L_OverheatedRegionC-L_2phaseRegionC)/(2.0*Lc))
X41 = ((G_refrigerant_condenserOutlet*c_p_condenserOutlet)/(Lc-L_OverheatedRegionC-L_2phaseRegionC))*(t_kondensasjon - t_condenserOutlet)-((T4/T2)*X21)
X42 = (G_refrigerant_condenserOutlet*c_p_condenserOutlet)-alfa_refrigerant*(Ac/(2.0*Lc))*(Lc-L_OverheatedRegionC-L_2phaseRegionC)-((T4/T2)*X22)
X43 = ((G_refrigerant_condenserOutlet*c_p_condenserOutlet)/(Lc-L_OverheatedRegionC-L_2phaseRegionC))*(t_kondensasjon - t_condenserOutlet)-((T4/T2)*X23)
X44 = -(G_refrigerant_condenserOutlet*c_p_condenserOutlet+alfa_refrigerant*(Ac/Lc)*(Lc-L_OverheatedRegionC-L_2phaseRegionC))
X45 = -((T4/T2)*X24)
X46 = -((T4/T2)*X25)
X47 = alfa_refrigerant*(Ac/Lc)*(Lc-L_OverheatedRegionC-L_2phaseRegionC)
X48 = -((T4/T2)*X26)
X49 = (t_kondensasjon - t_condenserOutlet)*c_p_condenserOutlet
X40 = -((T4/T2)*X27)
T5 = V_heatExchangerWallC*p_heatExchangerWallC*c_p_heatExchangerWallC*(L_OverheatedRegionC/Lc)
X51 = alfa_refrigerant*(Ac*L_OverheatedRegionC/(2.0*Lc))
X52 = (alfa_refrigerant*(Ac*L_OverheatedRegionC/Lc)-alfa_coolantMedium*(Ac*L_OverheatedRegionC/Lc))
X53 = alfa_coolantMedium*(Ac*L_OverheatedRegionC/(2.0*Lc))
X54 = alfa_coolantMedium*(Ac*L_OverheatedRegionC/(2.0*Lc))
X55 = alfa_refrigerant*(Ac*L_OverheatedRegionC/(2.0*Lc))
T6 = V_heatExchangerWallC*p_heatExchangerWallC*c_p_heatExchangerWallC*L_2phaseRegionC/Lc
X61 = alfa_refrigerant*L_2phaseRegionC*(Ac/Lc)
X62 = -(alfa_refrigerant*(Ac*L_2phaseRegionC/Lc)-alfa_coolantMedium*(Ac*L_2phaseRegionC/Lc))
X63 = alfa_refrigerant*L_2phaseRegionC*(Ac/(2.0*Lc))
X64 = alfa_refrigerant*L_2phaseRegionC*(Ac/(2.0*Lc))
T7 = V_heatExchangerWallC*p_heatExchangerWallC*c_p_heatExchangerWallC*(Lc-L_OverheatedRegionC-L_2phaseRegionC)/Lc
X71 = alfa_refrigerant*(Ac/(2.0*Lc))*(Lc-L_OverheatedRegionC-L_2phaseRegionC)
X72 = alfa_refrigerant*(Ac/(2.0*Lc))*(Lc-L_OverheatedRegionC-L_2phaseRegionC)
X73 = (alfa_refrigerant*(Ac/(Lc))*(Lc-L_OverheatedRegionC-L_2phaseRegionC)+alfa_coolantMedium*(Ac/(Lc))*(Lc-L_OverheatedRegionC-L_2phaseRegionC))
X74 = alfa_coolantMedium*(Ac/(2.0*Lc))*(Lc-L_OverheatedRegionC-L_2phaseRegionC)
X75 = alfa_coolantMedium*(Ac/(2.0*Lc))*(Lc-L_OverheatedRegionC-L_2phaseRegionC)
T8 = V_coolantMediumC*p_coolantMedium*c_p_coolantMedium*(Lc-L_OverheatedRegionC-L_2phaseRegionC)
X81 = G_coolantMediumC*c_p_coolantMedium*(t_coolantMedium_InletC-t_coolantMedium_twoPhaseRegionC)/(Lc-L_OverheatedRegionC-L_2phaseRegionC)
X82 = G_coolantMediumC*c_p_coolantMedium*(t_coolantMedium_InletC-t_coolantMedium_twoPhaseRegionC)/(Lc-L_OverheatedRegionC-L_2phaseRegionC)
X83 = alfa_coolantMedium*(Ac/Lc)*(Lc-L_OverheatedRegionC-L_2phaseRegionC)
X84 = -(G_coolantMediumC*c_p_coolantMedium+alfa_coolantMedium*(Ac/(2.0*Lc))*(Lc-L_OverheatedRegionC-L_2phaseRegionC))
X85 = G_coolantMediumC*c_p_coolantMedium-alfa_coolantMedium*(Ac/(2.0*Lc))*(Lc-L_OverheatedRegionC-L_2phaseRegionC)
X86 = c_p_coolantMedium*(t_coolantMedium_InletC-t_coolantMedium_twoPhaseRegionC)
T9 = V_coolantMediumC*p_coolantMedium*c_p_coolantMedium*L_2phaseRegionC/(2.0*Lc)
X91 = -(T9/-T8)*X81
X92 = (G_coolantMediumC*c_p_coolantMedium/L_2phaseRegionC)*(t_coolantMedium_twoPhaseRegionC-t_coolantMedium_overheatedRegionC) - (T9/T8)*X82
X93 = alfa_coolantMedium*(Ac/Lc)*L_2phaseRegionC
X94 = -(T9/T8)*X83
X95 = (G_coolantMediumC*c_p_coolantMedium-alfa_coolantMedium*(Ac/2.0*Lc)*L_2phaseRegionC) - (T9/T8)*X84
X96 = -(G_coolantMediumC*c_p_coolantMedium+alfa_coolantMedium*(Ac/2.0*Lc)*L_2phaseRegionC)
X97 = -(T9/T8)*X85
X98 = c_p_coolantMedium*(t_coolantMedium_twoPhaseRegionC-t_coolantMedium_overheatedRegionC)-(T9/T8)*X86
T0 = V_coolantMediumC*p_coolantMedium*c_p_coolantMedium*L_OverheatedRegionC/(2.0*Lc)
X01 = c_p_coolantMedium*(G_coolantMediumC/L_OverheatedRegionC)*(t_coolantMedium_overheatedRegionC-t_coolantMedium_OutletC)-(T0/T9)*X91
X02 = -(T0/T9)*X92
X03 = alfa_coolantMedium*(Ac/Lc)*L_OverheatedRegionC
X04 = -(T0/T9)*X93
X05 = -(T0/T9)*X94
X06 = -(T0/T9)*X95
X07 = c_p_coolantMedium*G_coolantMediumC - alfa_coolantMedium*(Ac/(2.0*Lc))*L_OverheatedRegionC -(T0/T9)*X96
X08 = -c_p_coolantMedium*G_coolantMediumC - alfa_coolantMedium*(Ac/(2.0*Lc))*L_OverheatedRegionC
X09 = -(T0/T9)*X97
X00 = c_p_coolantMedium*(t_coolantMedium_overheatedRegionC-t_coolantMedium_OutletC)-(T0/T9)*X98

x = input('Press "Enter" to continue!');
clc;
%% A-matrix
%
display('The A matrix, often known as the transition matrix.')
display('It explains the dynamic proparties of the system!')

A = [...
   X11/T1,X12/T1,0,0,X13/T1,0,0,0,0,0; ...
   X21/T2,X22/T2,X23/T2,0,X24/T2,X25/T2,0,0,0,0; ...
   0,X31/T3,X31/T3,0,0,X33/T3,0,0,0,0; ...
   X41/T4,X42/T4,X43/T4,X44/T4,X45/T4,X46/T4,X47/T4,0,0,0; ...
   0,X51/T5,0,0,X52/T5,0,0,0,X53/T5,X54/T5; ...
   0,X61/T6,0,0,0,X62/T6,0,X63/T6,X64/T6,0; ...
   0,X71/T7,0,X72/T7,0,0,X73/T7,X74/T7,0,0; ...
   X81/T8,0,X82/T8,0,0,0,X83/T8,X84/T8,0,0; ...
   X91/T9,0,X92/T9,0,0,X93/T9,X94/T9,X95/T9,X96/T9,0; ...
   X01/T0,0,X02/T0,0,X03/T0,X04/T0,X05/T0,X06/T0,X07/T0,X08/T0]


x = input('Press "Enter" to continue!');
clc;
%% B-Matrix
%
display('The B matrix, often known as the routing matrix/ input matrix.')
display('It explains the "routing" from input values to the state variables!')


B = [...
    X14/T1,0,X15/T1,0,0; ...
    X26/T2,0,X27/T2,0,0; ...
    0,X34/T3,0,0,0; ...
    X48/T4,X49/T4,X40/T4,0,0; ...
    0,0,X55/T5,0,0; ...
    0,0,0,0,0; ...
    0,0,0,X75/T7,0; ...
    0,0,0,X85/T8,X86/T8; ...
    0,0,0,X97/T9,X98/T9; ...
    0,0,0,X09/T0,X00/T0]


x = input('Press "Enter" to continue!');
clc;
%% C-Matrix
%
display('C matrix - "output matrix."')
display('The output linear-combination matrix.')
C = [...
    0,1,0,0,0,0,0,0,0,0; ...
    0,0,0,1,0,0,0,0,0,0; ...
    0,0,0,0,0,0,0,0,1,0]

%% D-Matrix
%
display('D matrix.')
display('Direct feedforward coefficient.')
D = [...
    0,0,0,0,0; ...
    0,0,0,0,0; ...
    0,0,0,0,0]
%D_1 = eye(3,5)


x = input('Press "Enter" to continue!');
clc;
%% system eq.
%
display('State system is based upon A, B, C and D matrix. Since this a large MIMO system')
display('with 5 control points and 3 outputs, we will convert it to 15 transfer functions.')
display('For an easier overview, input and output variables will be displayed on all the matrices.')
display('Afterwards we will display all the 15 transfer functions.')
display(' ')
x = input('Press "Enter" to continue!');
clc;

delta_X = {'L_OverheatedRegionC' 't_kondensasjon' 'L_2phaseRegionC' 't_condenserOutlet' 't_heatExchangerWall_OverheatedRegion' 't_heatExchangerWall_twoPhaseRegionC' 't_heatExchangerWall_OverCoolingRegionC' 't_coolantMedium_twoPhaseRegionC' 't_coolantMedium_overheatedRegionC' 't_coolantMedium_OutletC'};  
delta_U = {'G_refrigerant_CompressorOutlet' 'G_refrigerant_condenserOutlet' 't_condenserOutlet' 't_coolantMedium_InletC' 'G_coolantMediumC'};
delta_Y = {'t_kondensasjon' 't_condenserOutlet' 't_coolantMedium_OutletC'};

Sys_ss = ss(A,B,C,D,'statename',delta_X,'inputname',delta_U,'outputname',delta_Y)

x = input('Press "Enter" to continue!');
clc;

Sys_tf = tf(Sys_ss)

display(' ')

A_a = deltaX;
[deltaX] = gradient(A);
B_b = deltaU;
[deltaU] = gradient(B);
C_c = deltaY;
[deltaY] = gradient(C);

%deltaXdot = deltaX*A_a + deltaU*B_b
%deltaYdot = deltaX*C_c

x = input('Press "Enter" to continue!');
clc;

%% Transfer function's
%
display('We have 10 inputs, 5 control points and 3 outputs for our system.')
display('This means that we have 15 transfer functions for this system.')
display('For an example, here we have the transfer function for')
display('the output temperature from regulation of the mass flow')
display('rate by the the compressor:')
display(' ')
x = input('Press "Enter" to continue!');
clc;

TF_massFlowRateCompressorControl = Sys_tf(3,1)

%[num,den] = ss2tf(A,B,C,D,3);
%abc_Output1 = tf(num(1,:),den)
%abc_Output2 = tf(num(2,:),den)
%abc_Output3 = tf(num(3,:),den)

prompt = 'Do you want like to se the graphs of this TF? Y/N [Y]: ';
x = input(prompt,'s');
if strcmp(x,'y') | strcmp(x,'Y')
    figure(1)
    bode(TF_massFlowRateCompressorControl)
    grid on

    pause
    x = input('Press "Enter" to continue!');
    clc; close all;

    figure(1)
    step(TF_massFlowRateCompressorControl)
    grid on

    pause
    x = input('Press "Enter" to continue!');
    clc; close all;

    pole(TF_massFlowRateCompressorControl)
    figure(1)
    pzplot(TF_massFlowRateCompressorControl)
    grid on
end

x = input('Press "Enter" to continue!');
clc; close all;



%% Stability and reachability
%

display('It is important that all eigenvalues have negative real part.')
display('If at least one eigenvalue has a positive part,')
display('the equilibrium will be unstable.')
System_eig = eig(A);
System_det = det(A); %If stable, than everyting is reachable




display('The eigenvalue or characteristic polynomial is:')
System_eig

display(' ')
x = input('Press "Enter" to continue!');
clc;

display('For all negative determinants the system is stable,')
display('on the other hand, all the positives are unstable')
display('and need further adjustments to change pole placements.')

display('The determinants of the systems inputs are:')
System_det

%     display('The A matrix has a zero determinant and is "called singular or non-invertible"')
%     display('There exists a vector v =/= 0 such that Av = 0')
%     display('This means that NOT all components/ inputs are reachable')
%     display('We also know that this system is not stable!')

display(' ')

x = input('Press "Enter" to continue!');
clc; 

Stability_from_trace = trace(A)
display('By calculating the sum of the diagonal line of the A matrix')
display('we can calculate the trace of the system,')
display('which again will describe stability of the system - ')
display('whether it converges or diverges. It preferably should be negative.')

x = input('Press "Enter" to continue!');
clc; 

%% Observability
%
display('Observability depends on the A and C matrix, and not B!')
display('The A matrix has the dynamic components of the system.')
display('The C matrix is linear-combination of the output to the state variables!')
display('The A is a 10X10 matrix and C is a 3X10 matrix')
display('The matrix we will get is a 30X10 matrix of')
display('values that will determine observability!')
display(' ')



prompt = 'Do you want to see the observability matrix? Y/N [Y]: ';
x = input(prompt,'s');
if strcmp(x,'y') | strcmp(x,'Y')
    observability =  obsv(A,C)
    rank = rank(observability)
end

x = input('Press "Enter" to continue!');
clc;

%% Controlabilety
%
display('Controllability depends on the A and B matrix, and not C!')
display('The A matrix has the dynamic components of the system.')
display('The B matrix is routing input values to the state variables!')
display('The A is a 10X10 matrix and B is a 10X5 matrix')
display('We will therfore get 10X50 values that will determine controllability!')
display(' ')
prompt = 'Do you want to se the controlability matrix? Y/N [Y]: ';
x = input(prompt,'s');
if strcmp(x,'y') | strcmp(x,'Y')
    controlabilety =  ctrb(A,B)
end

x = input('Press "Enter" to continue!');
clc;

%%

%% Thank you for
%
display('Thank you for attention.')
display(' ')
x = input('The program is done. Press "Enter" to finish!');
clc; clear all; close all;


