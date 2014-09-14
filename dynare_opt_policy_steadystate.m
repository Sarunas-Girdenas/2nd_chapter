% Steady State file for optimal policy

function [ys,check]=dynare_opt_policy(ys,exe)

global M_ lgy_
global beta xcoef eta alpha xi lambda k a z epsilon omega Xst sigma Sst brFst inflst intrst pEst qst Kst Zst Yst Cst Nst thetast phist ust
global Gs1 Gs2 Gs3 Gs4 Gs5 Gs6  Gs7 Gs8 Gs9 Gs10 Gs11 Gs12 Cmin Smin Nmin brFmin intrmin Splus inflplus Cplus Zplus Kplus

if isfield(M_,'param_nbr') == 1
NumberOfParameters = M_.param_nbr;
for i = 1:NumberOfParameters
  paramname = deblank(M_.param_names(i,:));
  eval([ paramname ' = M_.params(' int2str(i) ');']);
end
check = 0;
end

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


%Here we solve for q & S steady state

syms qi Si

eq1 = a-(z/Xst)+Si*(phist*(1-eta)*(1-beta)+1+eta*beta*(qi^((-alpha)/(1-alpha)))*(xi^(1/(1-alpha)))-beta*(1-lambda));
eq2 = k/qi-beta*(phist+1)*(1-eta)*Si;

sol = solve(eq1,eq2);

qi = sol.qi;
Si = sol.Si;

S1 = double(Si);
q1 = double(qi);

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

inflst = 1

ust = 1-(1-lambda)*Nst

% Computing lagrange multipliesrs

f1=(lambda +eta*pEst)/(1+phist)/(1-eta);
f2=brFst/beta/Cst*sigma*(1-beta)*(1-f1);
f3=(beta*(1-lambda)*k*thetast-a+z)/(1-beta*(1-lambda)*(1-pEst));
f4=k*lambda*Nst*thetast/pEst-a*lambda*Nst*f3;
f5=alpha*pEst*eta*beta*Sst+(1-alpha)*f1*k/qst;
Gs1=Cst^(-sigma)*f4/(f4*f2+f5);
Gs12=-f3*(Cst^(-sigma)+f2*Gs1);
Gs3=f1*Gs1;
Gs8=Cst^(-sigma)+brFst/beta/Cst*(Gs1-Gs3)*sigma*(1-beta);
Gs9=Gs3*k/qst/pEst/alpha*(1-alpha);
Gs10=-Gs9-Gs1*eta*beta*Sst;
Gs7=(Gs3-Gs1)*brFst;
Gs11=-Gs8;
Gs5=Gs1/Xst/Nst;
Gs2=Gs3*beta;
Gs6=-Gs1/Nst;
Gs4=-Gs1*(1-omega)/Nst/Xst;

S        = 0;
brF      = 0;
infl     = 0;
intr     = 0;
pE       = 0;
q        = 0;
X        = 0;
Y        = 0;
C        = 0;
N        = 0;
theta    = 0;
phi      = 0;
u        = 0;
zz       = 0;
Z        = 0;
K        = 0;
G1       = 0;
G2       = 0;
G3       = 0;
G4       = 0;
G5       = 0;
G6       = 0;
G7       = 0;
G8       = 0;
G9       = 0;
G10      = 0;
G11      = 0;
G12      = 0;
Cmin     = 0;
Smin     = 0;
Nmin     = 0;
brFmin   = 0;
intrmin  = 0;
Splus    = 0;
inflplus = 0;
Cplus    = 0;
Zplus    = 0;
Kplus    = 0;

for iter = 1:length(M_.params)
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

if isfield(M_,'param_nbr') == 1

if isfield(M_,'orig_endo_nbr') == 1
NumberOfEndogenousVariables = M_.orig_endo_nbr;
else
NumberOfEndogenousVariables = M_.endo_nbr;
end
ys = zeros(NumberOfEndogenousVariables,1);
for i = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(i,:));
  eval(['ys(' int2str(i) ') = ' varname ';']);
end
else
ys=zeros(length(lgy_),1);
for i = 1:length(lgy_)
    ys(i) = eval(lgy_(i,:));
end
check = 0;
end
end