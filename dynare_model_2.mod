// This is the model file which was used for grid search algorithm. Therefore it contains dynare loop.
// Purpose of the loop: simulate model for all the values in the vectors A_infl, B_i, C_i.
// To run the file we have to create arrays which we want to load and place them in the same
// directory


var zz N S brF pE X theta q phi intr infl Y C u Z K Cmin Smin Nmin brFmin intrmin Splus inflplus Cplus Zplus Kplus;

varexo shock1 shock2;

parameters alpha_infl alpha_i alpha_u beta eta alpha xi lambda k a z phist epsilon omega Xst sigma qst Sst Nst Cst Zst pEst Kst brFst intrst Yst inflst thetast xcoef ust;

// dynare loop

load A_infl
load B_i
load C_i
alpha_infls= [A_infl;];
alpha_is= [B_i;];
alpha_us= [C_i;];
@#for i in 1:3375
alpha_infl = alpha_infls(@{i});
alpha_i = alpha_is(@{i});
alpha_u = alpha_us(@{i});
stoch_simul(order=1);
save results_@{i} oo_;
@#endfor

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

intr=(1-alpha_i)*(alpha_infl*infl+alpha_u*Y)+alpha_i*intr(-1); 

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

phi=0.65*phi(-1)-shock1;//0.75

zz=0.65*zz(-1)+shock2;//0.75

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

save results_1 oo_;


