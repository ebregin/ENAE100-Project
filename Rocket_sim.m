% ENAE100
% Simulation for rocket 


global thrust_time T0 mf m0 A Cd;

mass_rocket = .693; %Kg
m0=mass_rocket+.202; %Kg
mf=mass_rocket+.09; %Kg
T0=165; %N
thrust_time=1; %burn time
A=10.36358/1550; %m^2
Cd=.5;
x0=[0 0];
tspan=[0 40];

[tout,xout] = ode45(@xdotsimplerocket,tspan,x0);



function xdot = xdotsimplerocket(t,x)
global thrust_time T0 mf m0 A Cd;

if t<thrust_time
    T=T0;
    m=(mf-m0)/thrust_time*t+m0;
    m_dot=(mf-m0)/thrust_time;
else
    T=0;
    m=mf;
    m_dot=0;
end
    
xdot(1,1) =  x(2);

xdot(2,1) = (T-sign(x(2))*.5*rho(x(1))*x(2)^2*Cd*A-9.8*m-2*m_dot*x(2))/m;

end

function density=rho(h)
p0=101325; %Pa
T0=288.15; %K
L=.0065; %Temp lapse rate
R=8.31446; %Ideal gas constant
M=.0289652; %molar mass of air

density = p0*M/R/T0*(1-L*h/T0)^((9.8*M/R/L)-1);
end


