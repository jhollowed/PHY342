% This program solves ordinary differential equations by the 
%   Runge-Kutta-Fehlberg adaptive step algorithm. 
%
function MyRKF(t0,tf,S,hmax) 

[t,y]=rkf(t0,tf,S,hmax);
plot(t,y, '.r', 'markersize', 15)
figure
plot(t(1: end-1), diff(t), 'linewidth', 2);
xlabel('t', 'fontsize', 20)
ylabel('v,y', 'fontsize', 20')

%figure
%plot(y(1, :), y(2, :), 'linewidth', 2, 'color', 'r')
%xlabel('x', 'fontsize', 20);
%ylabel('v', 'fontsize', 20);
end

function [tt,SS]=rkf(t0,tf,S0,hmax)
% Specify relative error tolerance.
epsilon = 1.e-5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify coefficients used in the R-K-F algorithm: 

a1 = 0.25; a2 = 0.375; a3 = 12.0/13.0; a4 = 1.0; a5 = 0.5; 

%   The coefficients Bmn are used to determine the 'y' at 
%   which the derivative is evaluated. 
b10 =   0.25;         b20 =     3.0/  32.0; b21 =     9.0/  32.0; 
b30 =1932.0/2197.0;   b31 = -7200.0/2197.0; b32 =  7296.0/2197.0; 
b40 = 439.0/ 216.0;   b41 =    -8.0;        b42 =  3680.0/ 513.0;  
                      b43 = -845.0/4104.0; 
b50 =  -8.0  /27.0;   b51 =     2.0;        b52 = -3544.0/2565.0; 
                      b53 =1859.0/4104.0;   b54 =   -11.0/  40.0; 

%   The Cn are used to evaluate the solution YHAT. 
c0  =    16.0/  135.0; c2 =  6656.0/12825.0; c3 = 28561.0/56430.0;
                       c4 =    -0.18;        c5 = 2.0/55.0; 

%   The Dn are used to evaluate the error. 
d0  =     1.0/  360.0; d2 = -128.0/4275.0;  d3  = -2197.0/75240.0;
          d4 =    0.02;         d5 = 2.0/55.0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
          
h = hmax;
hmin= 0.0001;  %redefine as desired
index =1;
tt(1) = t0;
SS(:,1) = S0;

%%%%%%%%%%%%%%%%% Runge-Kutta-Fehlberg takes place now %%%%%%%%%%%%%%%
while (t0 < tf) 

    f0 = derivs(t0,S0); 
    t  = t0 + a1*h; 
    S  = S0 + b10*h*f0; 
    f1 = derivs(t,S); 
    t  = t0 + a2*h; 
    S  = S0 + b20*h*f0 + b21*h*f1; 
    f2 = derivs(t,S); 
    t  = t0 + a3*h; 
    S  = S0 + h*(b30*f0+b31*f1+b32*f2); 
    f3 = derivs(t,S); 
    t  = t0 + a4*h; 
    S  = S0 + h*(b40*f0+b41*f1+b42*f2+b43*f3); 
    f4 = derivs(t,S); 
    t  = t0 + a5*h; 
    S  = S0 + h*(b50*f0+b51*f1+b52*f2+b53*f3+b54*f4); 
    f5 = derivs(t,S); 

    Shat = S0 + h*(c0*f0 + c2*f2 + c3*f3 + c4*f4 + c5*f5);
    
%%%%%%%%%%%%%%%%%%%% Checking to see if tolerance met %%%%%%%%%%%%%%%%   

    if length(S0) == 1   % use this if ODE is 1D
      RelErr = abs (h*(d0*f0 + d2*f2 + d3*f3 + d4*f4 + d5*f5)/Shat) 
      hnew   = h * (epsilon/RelErr)^0.2;
    else                 % otherwise use this for error
      RelErr = sqrt(sum((h*(d0*f0 + d2*f2 + d3*f3 + d4*f4 + d5*f5)./Shat).^2)/length(Shat));
      hnew   = h * (epsilon/RelErr)^0.2;
    end

%%%%%% Making sure h stays within prescribe bounds--see page 225 %%%%%%%%
    if hnew>hmax
        hnew = hmax;
    end
    if hnew< hmin
       error(strcat('\n Possible Problem: "hnew" is too small \n',...
	't0=  %11.8f\n hnew = %8.4gf < hmin = \%8.4g'),t0,hnew,hmin)
    end 

%%%%%%%%%%% step size appropriate, storing values in tt and SS %%%%%%%%%%

    if( hnew >= h )   
        t0 = t0 + h;     
        S0 = Shat;        

        index = index +1;
        tt(index) = t0;
        SS(:,index) = Shat;
    end
%%%%%%%%%%% step size wasn't good, resetting %%%%%%%%%%%%%%%%%%%%%%%%%%

        h = 0.9 * hnew;
       
end 
end


function [der] = derivs(t,S) 
% This function evaluates the derivative of the function 
% der = exp(t)*sin(S);
der = [100*9.8 - 0.1*(S/100)^2];
end 
