function [x_dot] = functional_electrical_simulation(x, x_ext, u)

% Inputs
%  x: State variables
%    x(1): Tibialis anterior muscle activation (0 to 1)
%    x(2): Angle of the foot to the horizontal axis (ground)
%    x(3): Angular velocity of the foot
%  -----------------------------------------------------------------------
%  x_ext: External state variables
%    x_ext(1): Acceleration of the ankle joint in the vertical direction
%    x_ext(2): Acceleration of the ankle joint in the horizontal direction
%    x_ext(3): Angle from the shank to the vertical axis
%    x_ext(4): Angular velocity of the shank
%  -----------------------------------------------------------------------
%  u: Tibialis anterior muscle excitation (0 to 1)
%  =======================================================================
% Output
%  x_dot: State equations
%    x_dot(1): R.O.C of the tibialis anterior muscle activation
%    x_dot(2): x[3] (Angular velocity of the foot)
%    x_dot(3): Angular acceleration of the foot
%  -----------------------------------------------------------------------


% constants
T_act = 0.01; %seconds
T_deact = 0.04; %seconds
J = 0.0197; %kg*m^2
d = 3.7; %cm
B = 0.82;
c_f = 11.45; %cm
m_f = 1.0275; %kg
a_v = 1.33;
f_v1 = 0.18;
f_v2 = 0.023;
v_max = -0.9; %m/s
f_max = 600; %N
W = 0.56;
l_T = 22.3; %cm
l_MT0 = 32.1; %cm
a = [2.10, -0.08, 7.97, 0.19, -1.79];
g = -9.81; %m/s^2 


% Determined using literature and trig since it wasn't included in the
% paper
l_CE_opt = 11.0; %cm

% Other funcs
l_MT = l_MT0 + d*(x_ext(3) - x(2));
l_CE = l_MT - l_T;
v_CE = d*(x_ext(4) - x(3));
f_fv = (1 - (v_CE/v_max)) / (1 + v_CE/(v_max*f_v1));
if v_CE < 0
    f_fv = (1 + (a_v*v_CE/f_v2)) / (1 + v_CE/f_v2);
end 
f_fl = exp(-((l_CE - l_CE_opt)/(W*l_CE_opt)));

T_grav = -m_f * c_f * cos(x(2)) * g;
T_acc = m_f * c_f * (x_ext(1) * sin(x(2)) - x_ext(2) * cos(x(2)));
T_ela = exp(a(1) + a(2)*x(2)) - exp(a(3) + a(4)*x(2)) + a(5);
f_m = x(1) * f_max * f_fl * f_fv;

% Outputs
x_dot(1) = (u-x(1)) * ((u/T_act) - ((1-u)/T_deact));
x_dot(2) = x(3);
x_dot(3) = (1/J) * ((f_m*d) + T_grav + T_acc + T_ela + B*(x_ext(4) - x(3)));

end