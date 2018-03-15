*DICE123 with put file.
* Lines which begin with an asterisk (*) are remarks and
* will not be read by the computer using the GAMS program.

* This is an optimal growth model to calculate the optimal control
* rate and timing for the abatement of CO2 and other greenhouse gases.

* This is the revised version of the model as of December 1992.
* It is version DICE 1.2.3.
* The "1" indicates that it is a one-region model; the "2" that
* it is the second major version; and the "3" indicates that it
* uses the third-round estimate of the data.

* Documentation for this version is contained in W. D. Nordhaus,
* "Explaining the 'DICE'," Cowles Foundation Discussion Paper,
* January 1991.

* The calibration is to a 60-period run for the transversality

SETS     T         Time periods                                 /1*40/
         TFIRST(T) First period
         TLAST(T)  Last period

SCALARS  BET       Elasticity of marginal utility                /0/
         R         Rate of social time preference per year       /.03/
         GL0       Growth rate of population per decade          /.223/
         DLAB      Decline rate of pop growth per decade         /.195/
         DELTAM    Removal rate on carbon per decade             /.0833/
         GA0       Initial growth rate for technology per decade /.15/
         DELA      Decline rate of technol. change per decade    /.11  /
         SIG0      CO2-equivalent emissions-GNP ratio            /.519/
         GSIGMA    Growth of sigma per decade                    /-.1168/
         DK        Depreciation rate on capital per year         /.10/
         GAMA      Capital elasticity in production function     /.25/
         M0        CO2-equiv concentrations 1965 billions t C    /677/
         TL0       Lower stratum temperature (C) 1965            /.10/
         T0        Atmospheric temperature (C) 1965              /.2/
         ATRET     Marginal atmospheric retention rate           /.64/
         Q0        1965 world gross output trill 89 US dollars   /8.519/
         LL0       1965 world population millions                /3369/
         K0        1965 value capital trill 1989 US dollars      /16.03/
         C1        Climate-equation coefficient for upper level  /.226/
         LAM       Climate feedback factor                       /1.41/
         C3        Transfer coeffic. upper to lower stratum      /.440/
         C4        Transfer coeffic for lower level              /.02/
         A0        Initial level of total factor productivity    /.00963/
         A1        Damage coeff for co2 doubling (fraction GWP)  /.0133/
         B1        Intercept control cost function               /.0686/
         B2        Exponent of control cost function             /2.887/
         PHIK      Transversality coeff capital ($ per unit)     /140 /
         PHIM      Transversality coeff carbon ($ per ton)       /-9/
         PHITE     Transversalit coeff temper (bill $ per deg C) /-7000 /

PARAMETERS   L(T)          Level of population and labor
             AL(T)         Level of total factor productivity
             SIGMA(T)      CO2-equivalent-emissions output ratio
             RR(T)         Discount factor
             GA(T)         Growth rate of productivity from 0 to T
             FORCOTH(T)    Exogenous forcing for other greenhouse gases
             GL(T)         Growth rate of labor 0 to T
             GSIG(T)       Cumulative improvement of energy efficiency
             DUM(T)        Dummy variable 0 except 1 in last period;

TFIRST(T) = YES$(ORD(T) EQ 1);
TLAST(T)  = YES$(ORD(T) EQ CARD(T));
DISPLAY TFIRST, TLAST;

GL(T) = (GL0/DLAB)*(1-exp(-DLAB*(ord(t)-1)));
L(T)=LL0*exp(GL(t));
GA(T)  = (GA0/DELA)*(1-exp(-DELA*(ord(t)-1)));

AL(T) =  a0*exp(GA(t));
GSIG(T) = (GSIGMA/DELA)*(1-exp(-DELA*(ord(t)-1)));
SIGMA(T)=SIG0*exp(GSIG(t));

DUM(T)=1$(ord(T)  eq card(T));
RR(T) = (1+R)**(10*(1-ord(t)));

FORCOTH(T) = 1.42;
FORCOTH(T)$(ord(t) lt 15) = .2604+.125*ord(T)-.0034*ord(t)**2;


VARIABLES              MIU(T)         Emission control rate GHGs
                       FORC(T)        Radiative forcing, W per m2
                       TE(T)          Temperature, atmosphere C
                       TL(T)          Temperature, lower ocean C
                       M(T)           CO2-equivalent concentration bill t
                       E(T)           CO2-equivalent emissions bill t
                       C(T)           Consumption trill US dollars
                       K(T)           Capital stock trill US dollars
                       CPC(T)         Per capita consumption thousands US dol
                       PCY(t)         Per capita income thousands US dol
                       I(T)           Investment trill US dollars
                       S(T)           Savings rate as fraction of GWP
                       RI(T)          Real interest rate per annum
                       TRANS(T)       Transversality variable last period
                       Y(T)           Output

                       UTILITY;

POSITIVE VARIABLES MIU, E, TE, M, Y, C, K, I;

EQUATIONS     UTIL      Objective function
              YY(T)     Output equation
              CC(T)     Consumption equation
              KK(T)     Capital balance equation
              KK0(T)    Initial condition for K
              KC(T)     Terminal condition for K
              CPCE(t)   Per capita consumption definition
              PCYE(T)   Per capita income definition
              EE(T)     Emissions process
              SEQ(T)    Savings rate equation
              RIEQ(T)   Interest rate equation
              FORCE(T)  Radiative forcing equation
              MM(T)     CO2 distribution equation
              MM0(T)    Initial condition for M
              TTE(T)    Temperature-climate equation for atmosphere
              TTE0(T)   Initial condition for atmospheric temperature
              TLE(T)    Temperature-climate equation for lower oceans
              TRANSE(t) Transversality condition
              TLE0(T)   Initial condition for lower ocean;

* Equations of the model

KK(T)..        K(T+1) =L= (1-DK)**10 *K(T)+10*I(T);
KK0(TFIRST)..  K(TFIRST) =E= K0;
KC(TLAST)..    R*K(TLAST) =L= I(TLAST);

EE(T)..        E(T)=G=10*SIGMA(T)*(1-MIU(T))*AL(T)*L(T)**(1-GAMA)*K(T)**GAMA;
FORCE(T)..     FORC(T) =E=  4.1*(log(M(T)/590)/log(2))+FORCOTH(T);
MM0(TFIRST)..  M(TFIRST) =E= M0;
MM(T+1)..      M(T+1) =E= 590+ATRET*E(T)+(1 - DELTAM)*(M(T)-590);

TTE0(TFIRST).. TE(TFIRST) =E= T0;
TTE(T+1)..     TE(T+1) =E= TE(t)+C1*(FORC(t)-LAM*TE(t)-C3*(TE(t)-TL(t)));
TLE0(TFIRST).. TL(TFIRST) =E= TL0;
TLE(T+1)..     TL(T+1) =E= TL(T)+C4*(TE(T)-TL(T));

YY(T)..        Y(T) =E= AL(T)*L(T)**(1-GAMA)*K(T)**GAMA*(1-B1*(MIU(T)**B2))
                  /(1+(A1/9)*SQR(TE(T)));
SEQ(T)..       S(T) =e= I(T)/(.001+Y(T));
RIEQ(T)..      RI(T) =E= GAMA*Y(T)/K(T)- (1-(1-DK)**10)/10  ;

CC(T)..        C(T) =E= Y(T)-I(T);
CPCE(T)..      CPC(T) =e= C(T)*1000/L(T);
PCYE(T)..      PCY(T) =e= Y(T)*1000/L(T);

TRANSE(TLAST).. TRANS(TLAST)=E=RR(TLAST)
  *(PHIK*K(TLAST)+PHIM*M(TLAST)+PHITE*TE(TLAST));

UTIL..         UTILITY =E=
               SUM(T, 10 *RR(T)*L(T)*LOG(C(T)/L(T))/.55 +TRANS(T)*DUM(T));

*  Upper and Lower Bounds: General conditions imposed for stability

MIU.up(T) = 0.99;
MIU.lo(T) = 0.01;
K.lo(T) = 1;
TE.up(t) = 20;
M.lo(T) = 600;
C.lo(T) = 2;

* Upper and lower bounds for historical constraints

MIU.fx('1')=0.;
MIU.fx('2')=0.;
MIU.fx('3')=0.;

* Solution options

option iterlim = 99900;
option reslim = 99999;
option solprint = off;
option limrow = 0;
option limcol = 0;
model CO2 /all/;
solve CO2 maximizing UTILITY using nlp ;

* Display of results

display Y.l, C.l, S.l, K.l, MIU.l, E.l, M.l, TE.l, FORC.l, RI.l;
display CC.m, EE.m, KK.m, MM.m, TTE.m, CPC.l, TL.l, PCY.l, I.l;
display SIGMA, RR, L, AL, DUM, FORCOTH;


PARAMETERS  Tax(T);
tax(t) = -1*ee.m(t)*10000/cc.m(t);
FILE dice123a;
PUT dice123a;
dice123a.LJ = 1;
PUT "SIMULATION", @20, SYSTEM.TITLE //;
PUT "DATE",       @20, SYSTEM.DATE//;
PUT // @10, "RESULTS OF DICE MODEL";
PUT ///;
PUT // @10, "ECONOMIC RESULTS";
PUT // @5, "Y";
LOOP(T, PUT / T.TL:<2;
     PUT     @5,     Y.L(T):10:4);
PUT // @5 "C":
LOOP(T, PUT / T.TL:<2;
     PUT     @5,     C.L(T):10:4);
PUT // @10, "CLIMATIC RESULTS";
PUT // @5, "C tax";
LOOP(T, PUT / T.TL:<2;
     PUT     @5,     tax(T):8:5);
PUT // @10, "CLIMATIC RESULTS";

PUT // @10, "CLIMATIC RESULTS";
PUT // @5,  "Emiss";
LOOP(T, PUT / T.TL:<2;
     PUT     @5,      E.L(T):10:4);


PUT // @5, "M Conc";
LOOP(T, PUT / T.TL:<2;
     PUT     @5,    M.L(T):10:4);
PUT // @10, "CLIMATIC RESULTS";
PUT // @5, "TE";
LOOP(T, PUT / T.TL:<2;
     PUT     @5,   TE.L(T):8:5);

PUT // @10, "MARGINAL VALUES";

PUT // @10, "parameters";
PUT // @5, "AL";
LOOP(T, PUT / T.TL:<2;
     PUT     @5,    AL(T):10:4);




PUT // @5 "L";
LOOP(T, PUT / T.TL:<2;
     PUT     @5,  L(T):10:1);



PUT // @5, "Y", @15, "C", @25, "CPC", @35, "K", @45, "I";
LOOP(T, PUT / T.TL:<2;
     PUT     @5,     Y.L(T):10:4,
     PUT     @15,    C.L(T):10:4,
     PUT     @25,    CPC.L(T):10:5,
     PUT     @35,    K.L(T):10:4,
     PUT     @45,    I.L(T):10:4);

PUT // @10, "CLIMATIC RESULTS";
PUT // @5, "tax", @15, "miu", @25, "Emiss", @37, "M Conc", @49, "TE";
LOOP(T, PUT / T.TL:<2;
     PUT     @5,     tax(T):8:5,
     PUT     @15,    miu.l(T):10:4,
     PUT     @25,    E.L(T):10:4,
     PUT     @37,    M.L(T):10:4,
     PUT     @48,    TE.L(T):8:5);
PUT // @10, "MARGINAL VALUES";
PUT // @5, "TTE.M", @15, "M.M", @25, "K.M";
LOOP(T, PUT / T.TL:<2;

     PUT     @5,     TTE.M(T):10:4,
     PUT     @15,    MM.M(T):8:4,
     PUT     @25,    KK.M(T):8:4);
PUT // @10, "parameters";
PUT // @5, "AL", @15, "R", @25, "L", @35, "FORC";
LOOP(T, PUT / T.TL:<2;
     PUT     @5,     AL(T):10:4,
     PUT     @15,    RI.L(T):10:4,
     PUT     @25,    L(T):10:1,
     PUT     @37,    FORC.L(T):10:4);

execute_unload "DICE123.gdx"
