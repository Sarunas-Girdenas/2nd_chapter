% This file needs to be ran after dynare finished its computations
% First thing what we do is calculate steady state & then specify parameters
% needed for second order welfare approximation.
% This code calculates welfare. It corresponds to the model in the second
% chapter.
% before running this file dont clear the workspace. It needs variable 'M_'
% to obtain the order of variables in dynare code.
% also some other files should be pre-created (from looping dynare), they are loaded
% in the loop by name

%Parameters

beta      = 0.99; % discount factor
eta       = 0.5;  % relative bargaining power (workers)
alpha     = 0.5; %0.5;  % matching parameter
xi        = 0.76; % matching parameter (constant)
lambda    = 0.05; %0.05; % probability of match seperation
k         = 0.598; %0.711; % cost of posting vacancy
a         = 0.5; % Workers unemployment benefit
z         = 1.0;   % average productivity
phist     = 0.868; %enforcement parameter (shock)
epsilon   = 6; % 6 price elasticity paremeter
omega     = 0.25; %0.25; %0.75; %price adjustment probability
Xst       = epsilon/(epsilon-1); %price ratio/markup
sigma     = 2; %1; %risk aversion, comes from household utility function
xcoef     = (1/omega)*(1-omega)*(1-omega*beta); %Phillips curve coefficient


%Steady state variables

syms qi Si
eq1 = a-(z/Xst)+Si*(phist*(1-eta)*(1-beta)+1+eta*beta*(qi^((-alpha)/(1-alpha)))*(xi^(1/(1-alpha)))-beta*(1-lambda));
eq2 = k/qi-beta*(phist+1)*(1-eta)*Si;
sol = solve(eq1,eq2);
qi  = sol.qi;
Si  = sol.Si;
S1  = double(Si);
q1  = double(qi);
Sst = S1(1,1)
qst = q1(1,1)

%Since we solved for S q, we can compute all the other variables

brFst = phist*(1-eta)*Sst*beta
pEst = (qst^((-alpha)/(1-alpha)))*(xi^(1/(1-alpha)))
Nst = pEst/(lambda+pEst*(1-lambda))
Kst = (Nst*z)/(1-omega*beta)
Zst = Kst/Xst
thetast = (pEst/xi)^(1/alpha)
Cst = a- k*thetast+Nst*(z-a+k*thetast*(1-lambda))
intrst = 1/beta
Yst = Nst*z
inflst = 1;
ust = 1-(1-lambda)*Nst

%Lagrange Multipliers
%some parameters for computations of multipliers

f1 = (lambda+eta*pEst)/((1+phist)*(1-eta));
f2 = (brFst/(beta*Cst))*sigma*(1-beta)*(1-f1);
f3 = (beta*(1-lambda)*k*thetast-a+z)/(1-beta*(1-lambda)*(1-pEst));
f4 = k*lambda*Nst*thetast/pEst-a*lambda*Nst*f3;
f5 = alpha*pEst*eta*beta*Sst+(1-alpha)*f1*k/qst;

%Multipliers

Gs1  = Cst^(-sigma)*f4/(f4*f2+f5) %
Gs12 =-f3*(Cst^(-sigma)+f2*Gs1) %
Gs3  = f1*Gs1 %
Gs8  = Cst^(-sigma)+brFst/beta/Cst*(Gs1-Gs3)*sigma*(1-beta) %
Gs9  = Gs3*k/qst/pEst/alpha*(1-alpha) %
Gs10 = -Gs9-Gs1*eta*beta*Sst %
Gs7  = (Gs3-Gs1)*brFst %
Gs11 = -Gs8 %
Gs5  = Gs1/Xst/Nst %
Gs2  = Gs3*beta %
Gs6  = -Gs1/Nst %
Gs4  = -Gs1*(1-omega)/Nst/Xst %

%Some parameters we use for welfare computations

%Coeficcient for 4th equation

S41 = omega/(1-omega)*Kst*((epsilon*omega)/(1-omega)+epsilon-2);
S42 = 2*omega/(1-omega)*Kst;

%Coeficients for 5th equation
S51 = (-omega/intrst)*(epsilon-1)*(epsilon-2)*Kst;
S52 = 2*(omega/intrst)*(epsilon-1)*Kst;

%Coeficients for 6th equation
S61 = (-omega/intrst)*epsilon*(epsilon-1)*Zst;

%Coeficients for 7th equation
S71 = sigma*(-sigma-1);

%Coeficient for 9th equation
S91 = (-alpha/(1-alpha))*(1/(1-alpha))*(qst^(((alpha-2)/(1-alpha))+2))*(xi^(1/(1-alpha)));

%Coeficient for 10th equation
S101 = (-1)*xi*alpha*(alpha-1)*(thetast^alpha);
%Coeficients for 12 equation

S121 = (-xi)*alpha*(alpha-1)*(1-(1-lambda)*Nst)*(thetast^alpha);
S122 = 2*xi*alpha*(thetast^alpha)*Nst*(1-lambda);


for i = 1:41^2; %Number of steps (length of interval)
    name=sprintf('results_%d.mat',i); %load some data from dynare

    load(name);

%Now calculate welfare

%Extract matrices A & B
format short;
%The following piece of code extracts matrix A (variables) in the same
%order as variables are declared
% Let's get the variables list and call it C
C = oo_.dr.order_var;
%Now add it to our ghx matrix
A_1 = [C oo_.dr.ghx];
%Now sort the matrix in ascending order
A_2 = sortrows(A_1);
%Now delete the first column (which is sorting column)
A_2(:,1)=[];
%Call the matrix which is now is exactly the same is in dynare output
A = transpose(A_2);
n_pr = rows(A);
%Here we take predetermined variables
k1     = find(oo_.dr.kstate(:,2) <= M_.maximum_lag+1);
klag   = oo_.dr.kstate(k1, [1,2]);
state  = [ oo_.dr.order_var(klag(:,1)), klag(:,2)-M_.maximum_lag-2 ];
state1 = state(:,1);
A_sm   = A(:,state1);

%Now we do the same for shock matrix (ghu)

% Let's get the variables list and call it C
C = oo_.dr.order_var;
%Now add it to our ghu matrix
B_1 = [C oo_.dr.ghu];
%Now sort the matrix in ascending order
B_2 = sortrows(B_1);
%Now delete the first column (which is sorting column)
B_2(:,1) = [];
%Call the matrix which is now the same as in dynare output
B = transpose(B_2);
B;
%Taking predetermined variables
B_sm = B(:,state1);
%computing covariance matrix for predetermined variables
C1 = (beta^(1/2))*A_sm';
omega1 = [0.01 0;0 0.01]; % Shock correlation matrix, we call it omega1 because omega already exists
D_1 = (transpose(B_sm))*omega1*(B_sm)*(beta/(1-beta));
S_1_0 = zeros(n_pr,n_pr);
S_1 = D_1;
H_1 = C1;

format long;
%Sum the matrix

%t=0;

while max(max(abs(S_1-S_1_0)))>0.0000000000000000000000000001;
     S_1;
     S_1_0 = S_1;
     S_1 = S_1_0+H_1*S_1_0* H_1';
     H_1 = H_1^2;
     
     S_big=beta*A'*S_1*A+(beta/(1-beta))*B'*omega1*B;   
     
     %Calculate welfare for each step
     
     welfare1=   (-sigma*(Cst^(-1-sigma))*S_big(13,13))+...
                 Gs1*((-2/Xst)*z*S_big(6,6)+(2/Xst)*z*S_big(1,6)+2*brFst*intrst*S_big(11,11)+2*eta*beta*Sst*pEst*S_big(22,5)+2*intrst*brFst*(S_big(21,20)-S_big(20,11)-S_big(21,11)))+...    
                 Gs2*((2*brFst*intrst)*(S_big(23,23)+S_big(4,10)-S_big(4,23)-S_big(10,23))-2*phist*Sst*(1-eta)*S_big(9,22))+...
                 Gs3*2*k/qst*S_big(8,8)+...
                 Gs4*(S41*(S_big(11,11))+S42*S_big(16,11))+...
                 Gs5*(S51*S_big(23,23)-2*(omega/intrst)*Kst*S_big(10,10)+S52*S_big(10,23)+2*(omega/intrst)*Kst*S_big(10,26)-S52*S_big(23,26))+...          
                 Gs6*(S61*S_big(23,23)-2*(omega/intrst)*Zst*S_big(10,10)+2*(Yst/Xst)*S_big(12,6)-2*(Yst/Xst)*(S_big(6,6))+2*(omega/intrst)*epsilon*Zst*S_big(10,23)-2*(omega/intrst)*epsilon*Zst*S_big(25,23)+2*Zst*(omega/intrst)*S_big(25,10))+...
                 Gs7*(S71*(S_big(13,13))-sigma*(sigma-1)*S_big(24,24)+2*(sigma^2)*S_big(13,24)+2*sigma*S_big(13,23)-2*sigma*S_big(24,23))+...
                 Gs8*(2*k*(1-lambda)*Nst*thetast*S_big(7,19))+...
                 Gs9*(S91*S_big(8,8))+...
                 Gs10*(S101*S_big(7,7))+...
                 Gs11*(-2*Nst*z*S_big(2,1))+...
                 Gs12*(S121*S_big(7,7)+S122*S_big(7,19));
               %t=t+1     
                  x(i)=transpose(welfare1);
               
                                  
     end;


continue;

welfare2(i)=welfare1;
%i=i+1
%t=t+1

end;
     
% Plot 3D figure using gridfit

% The same length as in dynare file (parameters which we are looping) 
%alpha_u
%x=1.1:0.098:5.1;
%alpha_infl
%y=1.1:0.025:2.1;
%welfare which we want to plot
%z=welfare2;

% Plot the actual figure

xI = floor(min(x)):0.1:ceil(max(x));
yI = floor(min(y)):0.1:ceil(max(y));
%Using gridfit to interpolate surafe (data from dynare)
g = gridfit(x,y,z,xI,yI); % function gridfit could be obtained from the matlab central exchange
%Restrictin surface from being above zero (welfare)
g(g > 0) = 0;
%Plotting surface
V=surf(xI,yI,g);
%Setting axis
xlabel('Response to Unemployment');
ylabel('Response to Inflation');
zlabel('Unconditional Welfare');
%Rotate title of axis to make it look nicer
set(get(gca,'xlabel'),'rotation',15); %where angle is in degrees
set(get(gca,'ylabel'),'rotation',-25); %where angle is in degrees
set(get(gca,'zlabel'),'rotation',90); %where angle is in degrees
%Choosing nice looking color of the surface
colormap(bone);
