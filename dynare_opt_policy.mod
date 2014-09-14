// this is the first dynare file for the 2nd chapter


var zz N S brF pE X theta q phi intr infl Y C u Z K Cmin Smin Nmin brFmin intrmin Splus inflplus Cplus Zplus Kplus G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12;
varexo shock1 shock2;

parameters beta eta alpha xi ust lambda k a z phist epsilon omega Xst sigma qst Sst Nst Cst Zst pEst Kst brFst intrst Yst inflst thetast xcoef Gs1 Gs2 Gs3 Gs4 Gs5 Gs6 Gs7 Gs8 Gs9 Gs10 Gs11 Gs12;

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
//intr=1.5*infl;

// Unemployment

u*ust=-(1-lambda)*N(-1)*Nst;

// First equation

Sst*S(+1)*((1-lambda)*beta-eta*beta*pEst)=Sst*S+brFst*(intrst*(brF(-1)+intr(-1)-infl)-brF)+eta*beta*Sst*pEst*pE-z/Xst*(zz-X);

// Second equation

S(+1)+phi=brF+intr-infl(+1);

// Third equation

brFst*brF=(-k/qst)*q-beta*(1-eta)*Sst*S(+1);

// Fourth equation (Phillips curve)

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

phi=0.01*phi(-1)-shock1;

//Productivity shock

zz=0.01*zz(-1)+shock2;

//K
K=Y*(1-omega*beta)+beta * omega*(K(+1)-intr+infl(+1)*(epsilon-1));
//Z
Z=(Y-X)*(1-omega*beta)+beta * omega*(Z(+1)-intr+infl(+1)*(epsilon));
// dynamics of lagrange multipliers 13 equations from FOC
//1,i
-Gs1* brFst*(G1(+1)*brF-infl(+1)*intr)=brFst* Gs3*(G2+brF-infl(+1)+intr)+beta*Gs6*omega*Zst*(G6-intr+epsilon*infl(+1)+Z(+1))+Gs7*(G7+intr)+beta*Gs5*omega*Kst*(G5-intr+(epsilon-1)*infl(+1)+K(+1));
//2 Z
(1- omega)*G4=G6-omega*(G6(-1)-intr(-1)+epsilon*infl);
//3 N
Gs8*a*G8=Gs8*beta*(1-lambda)*k*(theta(+1)+G8(+1))-z*Gs11*(G11+zz)+G12*Gs12-beta *Gs12*(1-lambda)*(1-pEst)*G12(+1)+Gs12*beta*(1-lambda)*pEst*pE(+1);
//4 S
-Gs1*G1=-Gs1*((1-lambda)-eta*pEst)*G1(-1)+eta *pEst*Gs1*pE(-1)-Gs3*(1-eta)*(phist*(G2(-1)+phi)+G3(-1));
//5 C
-Cst^(1-sigma)*(1-sigma)*C= Gs7*sigma*(G7+sigma*(C(+1)-C)+infl(+1))-Gs7*sigma/beta*(G7(-1)+sigma*(C-C(-1))+infl)-Cst*Gs8*(G8+C);
// 6 X
G1+zz=G6+Y;
//7 b
Gs1*(G1(+1)+intr-infl(+1))-G1*Gs1+Gs3*(G2-G3-infl(+1)+intr)=0;
//8 pE
G1 *eta *beta*Sst*(G1+S(+1))+Gs9*G9+Gs10*G10=0;
//
G3-q=G9+pE;
//10 Y
Gs5*G5+Gs6/Xst*(G6-X)=Gs8*G8+Gs11*G11;
//11 K
Gs4*(G4+omega/(1-omega)*infl)+Gs5*G5=Gs5*omega*(G5(-1)-intr(-1)+(epsilon-1)*infl);
//12 infl
Gs1*brFst*intrst*(G1+brF(-1)+intr(-1)-infl)=-Gs3*brFst*intrst*(G2(-1)+brF(-1)+intr(-1)-infl) -(epsilon-1)*Gs5*Kst*omega*(G5(-1)-intr(-1)+(epsilon-1)*infl+ K)-epsilon *Zst*omega*Gs6*(G6(-1)-intr(-1)+epsilon*infl*Z)-1/beta*Gs7*(G7(-1)*sigma*(C-C(-1))+infl)+Gs4*omega/(1-omega)*Kst*(G4+(epsilon-1)*infl+omega/(1-omega)*epsilon*infl+K);
//13 theta
-Gs8*k*thetast*((1-(1-lambda)*Nst)*(G8*theta)-Nst*(1-lambda)*N(-1))-alpha*Gs10*pEst*(G10+pE)- Gs12*alpha*pEst*((1-(1-lambda)*Nst)*(G12+pE)-Nst*(1-lambda)*N(-1));
end;

shocks;
var shock1;
stderr 1;
var shock2;
stderr 1;

end;

resid(1);

steady;

check;

stoch_simul(irf=0,order = 1);



