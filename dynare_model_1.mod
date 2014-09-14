// this is the main dynare file for the model

var zz N S brF pE X theta q phi intr infl Y C u Z K Cmin Smin Nmin brFmin intrmin Splus inflplus Cplus Zplus Kplus;

varexo shock1 shock2;

parameters alpha_u alpha_i alpha_infl alpha_b beta eta alpha xi lambda k a z phist epsilon omega Xst sigma qst Sst Nst Cst Zst pEst Kst brFst intrst Yst inflst thetast xcoef ust;

// some parameters

alpha_infl=1.59;
alpha_i=0;
alpha_u=0.54;

model(linear);

Cmin=C(-1);
Smin=S(-1);
Nmin=N(-1);
brFmin=brF(-1);
intrmin=intr(-1);
Splus=S(+1);
inflplus=infl(+1);
Cplus=C(+1);
Zplus=Z(+1);
Kplus=K(+1);

// policy, taylor rule

intr=(1-alpha_i)*(alpha_infl*infl+alpha_u*u)+alpha_i*intr(-1); 

// Unemployment

u*ust=-(1-lambda)*N(-1)*Nst;

// First equation

Sst*S(+1)*((1-lambda)*beta-eta*beta*pEst)=Sst*S+brFst*(intrst*(brF(-1)+intr(-1)-infl)-brF)+eta*beta*Sst*pEst*pE-z/Xst*(zz-X);

// Second equation

S(+1)+phi=brF+intr-infl(+1);

// Third equation

brFst*brF=(-k/qst)*q-beta*(1-eta)*Sst*S(+1);

// Three equations instead of Phillips curve block

//(omega/(1-omega))*infl=Z-K;

K=(1-omega/intrst)*Y-(omega/intrst)*intr+omega*(epsilon-1)*(infl(+1)/intrst)+(omega/intrst)*K(+1);

Z=(1-omega/intrst)*(Y-X)+(omega/intrst)*(Z(+1)+epsilon*infl(+1)-intr);

//Phillips Curve

infl-beta*infl(+1)=xcoef*X;

// Fifth equation

intr=sigma*(C(+1)-C)+infl(+1);

// Sixth equation

Cst*C=Y*Yst-a*Nst*N-k*thetast*theta*(1-Nst*(1-lambda))+k*thetast*Nst*N(-1)*(1-lambda);

// Seventh equation

pE=-alpha*q/(1-alpha);

// Eight equation

pE=alpha*theta;

// Ninth equation

Y=N+zz;

// Tenth equation

Nst*N=Nst*N(-1)*(1-lambda)*(xi*(thetast^alpha)-1)-xi*alpha*(thetast^alpha)*theta*((1-lambda)*Nst-1);

// Eleventh equation (shock to firm's repayment probability)

phi=0.75*phi(-1)-shock1;//0.75

zz=0.75*zz(-1)+shock2;//0.75

end;

shocks;
var shock1;
stderr 1;
var shock2;
stderr 1;

end;

resid(1);

steady;

//options_.noprint=1; // supress output
check;

stoch_simul(irf=0,order = 1);

conditional_variance_decomposition=1;

// saving results
save results_1 oo_;


