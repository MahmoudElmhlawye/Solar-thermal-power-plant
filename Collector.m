% -------------------------------------------------------------------------------------------
% Parabolic trough Simple Design (PT) algorithm by Mahmoud Elmhlawye                 %
% Programmed by Mahmoud Elmhlawye at The University of                        %
% Science and Technology at Zewail City                                       %              
% Programming dates: Dec 2016-Jan 2017                                        %         %  
%                                               
% ------------------------------------------------------------------------------------------
% -- Citation Details:
% 1)Duffie, J. and Beckman, W. (2013). Solar engineering of thermal
% processes 
%. 1st ed. Hoboken: Wiley, pp.327-334.
% 2)Geankoplis, C. and Geankoplis, C. (2003). Transport
% processes and separation process principles. 1st ed. 
%Upper Saddle River, NJ: Prentice Hall Professional Technical Reference.
% ----------------------------------------------------------------%
% This Algorithm determines the plant area needed to provide      %
% Specific Power output; The collector used is Parabolic trough.  %
% Iterative method is used to determine the loss cofficient and   %
% Outlet Temperature.                                             %
% Averge values for the hear trasnfer fluid  should be used where %
% The Upper limit is the exit temperature and lower limit the     % 
% Inlet temperature.                                              %
% Efficiency of the plant power generation part was used as       %
% 35% for many rankine cycles.                                    % 
% The plant Area is assumed double the collector area             %
% --------------------------------------------------------------- %

% =============================================================== %
% Notes:                                                          %
% % Special thanks to Eng.Islam Ibrahim who taught the author 
% how to calculate the loss coefficient

%------------------------------------------------------------------

% The tube is made of stainless steel
% All unit SI except mentioned else
% Irradiance can be altered to fit your problem
%-----------------------------------------------------
%-----------------------------------------------------
%System Data 
function f=Collector(Q,Q_f,T_in,Texit,Cp,k3,density1)

%Q   Total thermal energy supplied by the collector system
%T_in  degree C;Inlet  T of thermal oil
%Texit degree c;outlet T of thermal oil
%Q_f Kg/s flow rate of heating fluid
%Cp averge heat capacity between Tin:Texit of thermanol (SI units)
%k3 averge Thermal Conductivity between Tin:Texit of thermanol (SI units)
%density1 averge density between Tin:Texit of thermanol (SI units)
L=100; %L   Length of one unit
W=5; % width of the parabolic concentrator
S=1000; % radiation per unit area of aperture
epsilon_r=0.31; %emittance of the reciever
D_ab=0.06; %diameter of absorber
D_o=0.09;  %outer diameter of absorber
t_o=0.004; %Thickness of outer layer

Tr=50+273; % reciever temperature **assumed**
%Cp=2.195*1000; %j/Kg.K
k1=16;% Conductivity of stainless steel W/m.C
k2=0.02046; %Conductivity of air
%k3= 0.106;   %Conductivity of thermanol
kc=1.4; %cover conductivity;
t_t=0.005; %Thickness of the tube
Tsky=5+273;
Tair=10+273;
vwind=4; %wind speed m/s
segma=5.67*10^-8; %Stefan–Boltzmann constant
%density1=885.1;
density2=1.232; %density of air
Viscosity1= 0.86*10^-3;
Viscosity2=1.794*10^-5;

%---------------------------------------------------
%---------------------------------------------------
%Clcution
%Convective heat transfer cofficients 
%wind heat transfer coefficient
Re1=density2*vwind*D_o/Viscosity2;
Nu1=0.3*Re1^0.6;
hw=Nu1*k2/D_o;
%---------------------------------------------------
%In pipe
v2=4*Q_f/((pi*D_ab^2)*density1);
 Re2=density1*v2*D_ab/Viscosity1;
 Pr=Viscosity1*Cp/k3;
 Nu2=0.023*(Re2^0.8)*(Pr^(0.3));
 h2=Nu2*k3/D_ab;
%loss coefficient calculation
Tr_new=400;
while abs(Tr-Tr_new)>0.1
    if Tr>Tr_new
        Tr=Tr-0.1*(Tr-Tr_new);
    else Tr=Tr+0.1*(Tr_new-Tr);
    end
    Tco=280; %initial guess surface Temperature
    Qloss1=100;
    Qloss2=200;
    while abs(Qloss1-Qloss2)>1
        Qloss1=pi*D_o*(hw*(Tco-Tair)+0.88*segma*(Tco^4-Tsky^4));
        Tci=Tco+Qloss1*log(D_o/(D_o-2*t_o))/(2*pi*kc*1);
        Qloss2=pi*D_ab*1*segma*(Tr^4-Tci^4)/(1/epsilon_r+(1-0.88)/(0.88)*(D_ab/(D_o-2*t_o)));
        Tco=Tco+0.01;
    end
    U_L=Qloss2/(pi*D_ab*1*(Tr-Tair));
    %----------------------------------------------------
    %Calculation of useful energy gain
    Ar=pi*D_ab*L; %outer area of reciever
    Aa=L*(W-D_o); %outer area of Collector
    C=Aa/Ar;
    F_dash=(1/U_L)/((1/U_L)+(D_ab/(h2*(D_ab-2*t_t)))+((D_ab/2*k1)*log(D_ab/(D_ab-2*t_t))));
    factor=Q_f*Cp/(Ar*U_L*F_dash);
    F_doubleDash=factor*(1-exp(-1/factor));
    F_R=F_doubleDash*F_dash;
    Qu=F_R*Aa*(S-(Ar/Aa)*U_L*(Tr-Tair));
    deltaT=Qu/(Q_f*Cp);
    Texit=T_in+deltaT;
    T_drop=Qu*((1/(pi*(D_ab-2*t_t)*L*h2))+(log(D_ab/(D_ab-2*t_t))/(2*pi*k1*L)));
    Tr_new=(T_in+Texit)/2 +T_drop;
end
efficiency=Q_f*Cp*(Texit-T_in)*100/(S*W*L);
power=Q_f*Cp*(Texit-T_in); % for each single module

N.modules=Q/(power); 
TotalArea=N.modules*Aa*2;

f=[N.modules efficiency TotalArea];


% You may usethose for presenting data 

% fprintf('        ***************String properties***************  \n');
% fprintf('UnitLength(m)  CollectorArea(m2)  ConcentrationRatio   UnitsNumber   StringPower(KW) \n');
% fprintf('%1f  %17f  %17.6f %18.6f %13.6f    \n',L/10,Aa/10,C,10,ceil(power/1000));
% 
% fprintf('        ***************Plant properties***************  \n');
% fprintf('T_inlet(C)  T_exit(C)   PlantPower(MW) ThermalEfficiency   TotalEfficiency  StringsNumber   Total Area(m2)  \n');
% fprintf('%1f  %10f  %11.6f %18.6f %15f %16.6f %17.6f   \n',T_in-273,round(Texit-273),StationPower/10^6,round(efficiency),round(efficiency*0.35),ceil(N.modules),ceil(TotalArea));
% 
% fprintf('        ***************Flow properties***************  \n');
% fprintf('Flow rate(Kg/s)  Speed(m/s)   Total thermal oil used in plant(kg)  \n');
% fprintf('%1f  %15f  %20.6f    \n',Q_f,v2,Q_f*L*N.modules/v2);
% 






