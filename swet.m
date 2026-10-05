%% Regmi Pramishan % MAE-3303 (Supersonic Airfoil Project)

clc;clear;close all;

%% INPUT VARIABLE

for aoa =6      % Angle of attack in degrees
aoa_radian = deg2rad(aoa);
for xu1 =0.9    %X-axis position for the maximum thickness of first top panel
for xl1 = 0.6   %X-axis position for the maximum thickness of first lower panel 
for tu = 0.03   % Maximum Y-axis position for the top panels

%% Constraints

M_inf = 2.5;    % Incoming Mach number of the airflow
c = 1;           % Chord length
for t_max = 0.1 % Maximum thickness of an airfoil 5% of the chord
tl = t_max - tu;%Maximum Y-axis position for the lower panels
gamma = 1.4;    %Specific heat ratio of an air
theta_max=thetamax(M_inf,gamma);
q_inf = 0.5*gamma*M_inf^2;

%% Geometrical Calculations

    %Upper surface
theta_u1 = atan(tu/xu1);   %Angle of first top panel
xu2 = c - xu1;  %X-axis position for the maximum thickness of second top panel
theta_u2 = atan(tu/xu2);%Angle of second top panel
   
    % Lower Surface
theta_l1 = atan(tl/xl1);        %Angle of first lower panel
xl2 = c-xl1;    %X-axis position for the maximum thickness of second lower panel
theta_l2 = atan(tl/xl2);    %Angle of second lower panel

%% SHOCK WAVE EXPANSION THEORY
    
    %Calculation on upper surface
if aoa_radian>theta_u1

    disp('expansionfan')

    theta_u = rad2deg(aoa_radian - theta_u1);    
    Pr_inf = (1+ ((gamma-1)/2)*(M_inf^2))^(gamma/(gamma-1)); 
    [M1_u] = pm(M_inf,theta_u,gamma); 
    Pr_u1 = (1+ ((gamma-1)/2)*(M1_u^2))^(gamma/(gamma-1)); 
    theta_2 = rad2deg(theta_u1) + rad2deg(theta_u2);
    [M2_u] = pm(M1_u,theta_2,gamma); 
    Pr_u2 = (1+ ((gamma-1)/2)*(M2_u^2))^(gamma/(gamma-1));
 
elseif aoa_radian<theta_u1
    
    disp('obliqueshockwave')

theta_u = rad2deg(theta_u1);
[beta] = oswbeta(M_inf,theta_u,gamma);
beta = deg2rad(beta);
Mi_n = M_inf*sin(beta);
M1_n = sqrt((1+(0.5*(gamma-1)*(Mi_n)^2))/((gamma*(Mi_n)^2)-(0.5*(gamma-1))));
theta_u = deg2rad(theta_u);
M1_bu = M1_n/sin(beta-theta_u);
Pr_iu =1 + (((2*gamma)/(gamma+1))*((Mi_n^2)-1)) ; 
Pr_u1 = (1+ ((gamma-1)/2)*(M1_bu^2))^(gamma/(gamma-1));

theta_2 = rad2deg(theta_u1) + rad2deg(theta_u2);
[M2_u] = pm(M1_bu,theta_2,gamma);
Pr_u2 = (1+ ((gamma-1)/2)*(M2_u^2))^(gamma/(gamma-1));

else

disp('No shock or expansion')

 Pr_inf = (1+ ((gamma-1)/2)*(M_inf^2))^(gamma/(gamma-1)); 
 theta_2 = rad2deg(theta_u1) + rad2deg(theta_u2);
 [M2_u] = pm(M_inf,theta_2,gamma); 
 Pr_u2 = (1+ ((gamma-1)/2)*(M2_u^2))^(gamma/(gamma-1));
 
end

% Calculation on the lower surface edge
theta_l = rad2deg(aoa_radian + theta_l1);
[beta] = oswbeta(M_inf,theta_l,gamma);
beta = deg2rad(beta);
Mi_n = M_inf*sin(beta);
M1_n = sqrt((1+(.5*(gamma-1)*(Mi_n)^2))/((gamma*(Mi_n)^2)-(.5*(gamma-1))));
theta_l = deg2rad(theta_l);
M1_l = M1_n/sin(beta-theta_l);

Pr_il =1 + (((2*gamma)/(gamma+1))*((Mi_n^2)-1)) ; 
Pr_l1 = (1+ ((gamma-1)/2)*(M1_l^2))^(gamma/(gamma-1));
 
theta_2l = rad2deg(theta_l1)+ rad2deg(theta_l2);
[M2_l] = pm(M1_l,theta_2l,gamma); 
Pr_l2 = (1+ ((gamma-1)/2)*(M2_l^2))^(gamma/(gamma-1));

%% Pressure Coefficients 

   % Upper Surface
   if aoa_radian>theta_u1
       CP_u1 = ((Pr_inf/Pr_u1) - 1)/ q_inf; 
       CP_u2 = ((Pr_inf/Pr_u2) - 1)/ q_inf;

   elseif aoa_radian<theta_u1
       CP_u1 = (Pr_iu-1)/q_inf;
       CP_u2 = ((Pr_iu*Pr_u1/Pr_u2) - 1)/ q_inf;

   else
       CP_u1 = 0;
       CP_u2 = ((Pr_inf*Pr_u2) - 1)/ q_inf ;
   end

    %Lower Surface
 CP_l1 = (Pr_il - 1)/ q_inf ;
 CP_l2 = ((Pr_il*Pr_l1/Pr_l2) - 1)/ q_inf ;
 
 %% Normal and Axial Coefficient
 
 C_n  = CP_l1 * (xl1) + CP_l2*xl2 - CP_u1 * xu1 - CP_u2*xu2;
 C_a  = CP_u1 *tan((theta_u1))* xu1+ CP_u2*tan((-theta_u2))*xu2 - CP_l1*tan((-theta_l1))* xl1 - CP_l2*(tan(theta_l2))*xl2;
 
 %% Lift and Drag Coefficient
 
 C_l = C_n*cos(aoa_radian) - C_a * sin(aoa_radian);
 C_d = C_n*sin(aoa_radian) + C_a * cos(aoa_radian);
 
 %% Pitching Moment about leading edge

 Cm_LE =  CP_u1*(0.5*xu1^2)+CP_u2*(0.5*(1-xu1^2))...
        - CP_l1*(0.5*xl1^2)-CP_l2*(0.5*(1-xl1^2))...
        + CP_u1*tan(theta_u1)*(0.5*tu*xu1)+CP_u2*tan(theta_u2)*(0.5*(tu)*(1-xu1))...
        - CP_l1*tan(-theta_l1)*(0.5*tl*xl1)- CP_l2*tan(theta_u2)*(0.5*(tl)*(1-xl1));
    
   %% LIFT TO DRAG RATIO

   Max_Lift = C_l/C_d;
   T=[aoa,xu1,xl1,tu,tl,Cm_LE,C_d,C_l,Max_Lift]
  
   %% Plot of an airfoil

   au=[0 xu1 1];
   al=[0 xl1 1];
   bu=[0 tu 0];
   bl=[0 -tl 0];
   plot (au,bu,al,bl)
end
end
end
end
end


