$ontext
DICE delta version 8
July 17, 2008.
This version is used for the DICE book, A Question of Balance (YUP, 2008).
We have included only the base, Hotelling, and optimal runs.
Exclude statements are removed so that it can run as a self-contained program.
Created September 5, 2008.

Note that this can be loaded into a data reading program,
$offtext

SETS  T                 Time periods                     /1*60/ ;

SCALARS

** Preferences
 B_ELASMU   Elasticity of marginal utility of consumption     /  2.0    /
 B_PRSTP    Initial rate of social time preference per year   / .015    /

** Population and technology
 POP0     2005 world population millions                  /6514     /
 GPOP0    Growth rate of population per decade            /.35      /
 POPASYM  Asymptotic population                           / 8600    /
 A0       Initial level of total factor productivity      /.02722   /
 GA0      Initial growth rate for technology per decade   /.092      /
 DELA     Decline rate of technol change per decade       /.001     /
 DK       Depreciation rate on capital per year           /.100     /
 GAMA     Capital elasticity in production function       /.300     /
 Q0       2005 world gross output trill 2005 US dollars   /61.1     /
 K0       2005 value capital trill 2005 US dollars        /137.     /

** Emissions
 SIG0     CO2-equivalent emissions-GNP ratio 2005         /.13418    /
 GSIGMA   Initial growth of sigma per decade              /-.0730    /
 DSIG     Decline rate of decarbonization per decade      /.003   /
 DSIG2    Quadratic term in decarbonization               / .000   /
 ELAND0   Carbon emissions from land 2005(GtC per decade) / 11.000  /

** Carbon cycle
 MAT2000  Concentration in atmosphere 2005 (GtC)          /808.9   /
 MU2000   Concentration in upper strata 2005 (GtC)        /1255     /
 ML2000   Concentration in lower strata 2005 (GtC)        /18365    /
 b11      Carbon cycle transition matrix                  /0.810712 /
 b12      Carbon cycle transition matrix                  /0.189288 /
 b21      Carbon cycle transition matrix                  /0.097213 /
 b22      Carbon cycle transition matrix                  /0.852787 /
 b23      Carbon cycle transition matrix                  /0.05     /
 b32      Carbon cycle transition matrix                  /0.003119 /
 b33      Carbon cycle transition matrix                  /0.996881 /

** Climate model
 T2XCO2   Equilibrium temp impact of CO2 doubling oC      / 3    /
 FEX0     Estimate of 2000 forcings of non-CO2 GHG        / -.06   /
 FEX1     Estimate of 2100 forcings of non-CO2 GHG        / 0.30   /
 TOCEAN0  2000 lower strat. temp change (C) from 1900     /.0068   /
 TATM0    2000 atmospheric temp change (C)from 1900       /.7307   /
 C1       Climate-equation coefficient for upper level    /.220    /
 C3       Transfer coeffic upper to lower stratum         /.300    /
 C4       Transfer coeffic for lower level                /.050    /
 FCO22X   Estimated forcings of equilibrium co2 doubling  /3.8     /

** Climate damage parameters calibrated for quadratic at 2.5 C for 2105
 A1       Damage intercept                                / 0.00000    /
 A2       Damage quadratic term                           /  0.0028388 /
 A3       Damage exponent                                 / 2.00       /

** Abatement cost
 EXPCOST2   Exponent of control cost function               /2.8   /
 PBACK      Cost of backstop 2005 000$ per tC 2005          /1.17  /
 BACKRAT    Ratio initial to final backstop cost            / 2    /
 GBACK      Initial cost decline backstop pc per decade     /.05   /
 LIMMIU     Upper limit on control rate                     / 1    /

** Participation
 PARTFRACT1  Fraction of emissions under control regime 2005 /1      /
 PARTFRACT2  Fraction of emissions under control regime 2015 /1      /
 PARTFRACT21 Fraction of emissions under control regime 2205 /1      /
 DPARTFRACT  Decline rate of participation                   /0      /

** Availability of fossil fuels
 FOSSLIM  Maximum cumulative extraction fossil fuels         / 6000  /

** Scaling and inessential parameters
  scale1 Scaling coefficient in the objective function       /194    /
  scale2 Scaling coefficient in the objective function       /381800 / ;

* Definitions for outputs of no economic interest
SETS
      TFIRST(T)
      TLAST(T)
      TEARLY(T)
      TLATE(T);

PARAMETERS
  L(T)          Level of population and labor
  AL(T)         Level of total factor productivity
  SIGMA(T)      CO2-equivalent-emissions output ratio
  R(T)          Instantaeous rate of social time preference
  RR(T)         Average utility social discount rate
  GA(T)         Growth rate of productivity from 0 to T
  FORCOTH(T)    Exogenous forcing for other greenhouse gases
  GL(T)         Growth rate of labor 0 to T
  GCOST1        Growth of cost factor
  GSIG(T)       Cumulative improvement of energy efficiency
  ETREE(T)      Emissions from deforestation
  COST1(t)      Adjusted cost for backstop
  PARTFRACT(T)  Fraction of emissions in control regime
  AA1           Variable A1
  AA2           Variable A2
  AA3           Variable A3
  ELASMU        Variable elasticity of marginal utility of consumption
  PRSTP         Variable nitial rate of social time preference per year
  LAM           Climate model parameter
  Gfacpop(T)    Growth factor population ;

PARAMETERS
  L(T)          Level of population and labor
  AL(T)         Level of total factor productivity
  SIGMA(T)      CO2-equivalent-emissions output ratio
  RR(T)         Average utility social discount factor
  GA(T)         Growth rate of productivity from 0 to T
  FORCOTH(T)    Exogenous forcing for other greenhouse gases
  GL(T)         Growth rate of labor 0 to T
  GCOST1        Growth of cost factor
  GSIG(T)       Cumulative improvement of energy efficiency
  ETREE(T)      Emissions from deforestation
  COST1(t)      Adjusted cost for backstop
  PARTFRACT(T)  Fraction of emissions in control regime
  AA1           Variable A1
  AA2           Variable A2
  AA3           Variable A3
  ELASMU        Variable elasticity of marginal utility of consumption
  PRSTP         Variable nitial rate of social time preference per year
  LAM           Climate model parameter
  Gfacpop(T)    Growth factor population ;

* Unimportant definitions to reset runs
TFIRST(T) = YES$(ORD(T) EQ 1);
TLAST(T)  = YES$(ORD(T) EQ CARD(T));
TEARLY(T) = YES$(ORD(T) LE 20);
TLATE(T)  = YES$(ORD(T) GE 21);
AA1 = A1;
AA2 = A2;
AA3 = A3;
ELASMU = B_ELASMU;
PRSTP  = B_PRSTP;

b11 = 1 - b12;
b21 = 587.473*B12/1143.894;
b22 = 1 - b21 - b23;
b32 = 1143.894*b23/18340;
b33 = 1 - b32 ;


* Important parameters for the model
LAM     = FCO22X/ T2XCO2;
Gfacpop(T) =   (exp(gpop0*(ORD(T)-1))-1)/exp(gpop0*(ORD(T)-1));
L(T)=POP0* (1- Gfacpop(T))+Gfacpop(T)*popasym;
ga(T)=ga0*EXP(-dela*10*(ORD(T)-1));
al("1") = a0;
LOOP(T, al(T+1)=al(T)/((1-ga(T))););
gsig(T)=gsigma*EXP(-dsig*10*(ORD(T)-1)-dsig2*10*((ord(t)-1)**2));sigma("1")=sig0;LOOP(T,sigma(T+1)=(sigma(T)/((1-gsig(T+1)))););
cost1(T) = (PBACK*SIGMA(T)/EXPCOST2)* ( (BACKRAT-1+ EXP (-gback* (ORD(T)-1) ) )/BACKRAT);
ETREE(T) = ELAND0*(1-0.1)**(ord(T)-1);
RR(t)=1/((1+prstp)**(10*(ord(T)-1)));
FORCOTH(T)= FEX0+ .1*(FEX1-FEX0)*(ORD(T)-1)$(ORD(T) LT 12)+ 0.36$(ORD(T) GE 12);
partfract(t) = partfract21;
PARTFRACT(T)$(ord(T)<25) = Partfract21 + (PARTFRACT2-Partfract21)*exp(-DPARTFRACT*(ORD(T)-2));
partfract("1")= PARTFRACT1;


VARIABLES
 MIU(T)          Emission control rate GHGs
 FORC(T)         Radiative forcing in watts per m2
 TATM(T)         Temperature of atmosphere in degrees C
 TOCEAN(T)       Temperatureof lower oceans degrees C
 MAT(T)          Carbon concentration in atmosphere GtC
 MATAV(T)        Average concentrations
 MU(T)           Carbon concentration in shallow oceans Gtc
 ML(T)           Carbon concentration in lower oceans GtC
 E(T)            CO2-equivalent emissions GtC
 C(T)            Consumption trillions US dollars
 K(T)            Capital stock trillions US dollars
 CPC(T)          Per capita consumption thousands US dollars
 PCY(t)          Per capita income thousands US dollars
 I(T)            Investment trillions US dollars
 S(T)            Gross savings rate as fraction of gross world product
 RI(T)           Real interest rate per annum
 Y(T)            Gross world product net of abatement and damages
 YGROSS(T)       Gross world product GROSS of abatement and damages
 YNET(T)         Output net of damages equation
 DAMAGES(T)      Damages
 ABATECOST(T)    Cost of emissions reductions
 CCA(T)          Cumulative industrial carbon emissions GTC
 PERIODU(t)      One period utility function
 UTILITY;

POSITIVE VARIABLES MIU, TATM, TOCE, E, MAT, MATAV, MU, ML, Y, YGROSS, C, K, I, CCA ;

EQUATIONS

 CCTFIRST(T)      First period cumulative carbon
 CCACCA(T)        Cumulative carbon emissions
 UTIL             Objective function
 YY(T)            Output net equation
 YNETEQ(T)        Output net of damages equation
 YGROSSEQ(T)      Output gross equation
 DAMEQ(T)         Damage equation
 ABATEEQ(T)       Cost of emissions reductions equation
 CC(T)            Consumption equation
 KK(T)            Capital balance equation
 KK0(T)           Initial condition for capital
 KC(T)            Terminal condition for capital
 CPCE(t)          Per capita consumption definition
 PCYE(T)          Per capita income definition
 EE(T)            Emissions equation
 SEQ(T)           Savings rate equation
 RIEQ(T)          Interest rate equation
 FORCE(T)         Radiative forcing equation
 MMAT0(T)         Starting atmospheric concentration
 MMAT(T)          Atmospheric concentration equation
 MMATAVEQ(t)      Average concentrations equation
 MMU0(T)          Initial shallow ocean concentration
 MMU(T)           Shallow ocean concentration
 MML0(T)          Initial lower ocean concentration
 MML(T)           Lower ocean concentration
 TATMEQ(T)        Temperature-climate equation for atmosphere
 TATM0EQ(T)       Initial condition for atmospheric temperature
 TOCEANEQ(T)      Temperature-climate equation for lower oceans
 TOCEAN0EQ(T)     Initial condition for lower ocean temperature
 PERIODUEQ(t)     Instantaneous utility function equation  ;

** Equations of the model

CCTFIRST(TFIRST).. CCA(TFIRST)=E=0;
CCACCA(T+1)..      CCA(T+1)=E=CCA(T)+ E(T);
KK(T)..            K(T+1) =L= (1-DK)**10 *K(T)+10*I(T);
KK0(TFIRST)..      K(TFIRST) =E= K0;
KC(TLAST)..        .02*K(TLAST) =L= I(TLAST);
EE(T)..            E(T)=E=10*SIGMA(T)*(1-MIU(T))*AL(T)*L(T)**(1-GAMA)*K(T)**GAMA + ETREE(T);
FORCE(T)..         FORC(T) =E=  FCO22X*((log((Matav(T)+.000001)/596.4)/log(2)))+FORCOTH(T);
MMAT0(TFIRST)..    MAT(TFIRST) =E= MAT2000;
MMU0(TFIRST)..     MU(TFIRST)  =E= MU2000;
MML0(TFIRST)..     ML(TFIRST)  =E= ML2000;
MMAT(T+1)..        MAT(T+1)    =E= MAT(T)*b11+MU(T)*b21 + E(T);
MMATAVEQ(t)..      MATAV(T)    =e= (MAT(T)+MAT(T+1))/2 ;
MML(T+1)..         ML(T+1)     =E= ML(T)*b33+b23*MU(T);
MMU(T+1)..         MU(T+1)     =E= MAT(T)*b12+MU(T)*b22+ML(T)*b32;
TATM0EQ(TFIRST)..  TATM(TFIRST) =E= TATM0;
TATMEQ(T+1)..      TATM(T+1) =E= TATM(t)+C1*(FORC(t+1)-LAM*TATM(t)-C3*(TATM(t)-TOCEAN(t)));
TOCEAN0EQ(TFIRST)..  TOCEAN(TFIRST) =E= TOCEAN0;
TOCEANEQ(T+1)..    TOCEAN(T+1) =E= TOCEAN(T)+C4*(TATM(T)-TOCEAN(T));
YGROSSEQ(T)..   YGROSS(T) =e= AL(T)*L(T)**(1-GAMA)*K(T)**GAMA;
DAMEQ(T)..      DAMAGES(t) =E= YGROSS(T)- YGROSS(T)/(1+aa1*TATM(T)+ aa2*TATM(T)**aa3);
YNETEQ(T)..     YNET(T) =E=  YGROSS(T)/(1+aa1*TATM(T)+ aa2*TATM(T)**aa3);
ABATEEQ(T)..    ABATECOST(T) =E= (PARTFRACT(T)**(1-expcost2))*YGROSS(T)*(cost1(t)*(MIU(T)**EXPcost2));
YY(T)..         Y(T) =E= YGROSS(T)*((1-(PARTFRACT(T)**(1-expcost2))*cost1(t)*(MIU(T)**EXPcost2)))/(1+aa1*TATM(T)+ aa2*TATM(T)**aa3);
SEQ(T)..        S(T)    =E= I(T)/(.001+Y(T));
RIEQ(T)..       RI(T)   =E= GAMA*Y(T)/K(T)- (1-(1-DK)**10)/10  ;
CC(T)..         C(T)    =E= Y(T)-I(T);
CPCE(T)..       CPC(T)  =E= C(T)*1000/L(T);
PCYE(T)..       PCY(T)  =E= Y(T)*1000/L(T);
PERIODUEQ(T)..  PERIODU(T)  =E=   ((C(T)/L(T))**(1-ELASMU)-1)/(1-ELASMU);
UTIL..          UTILITY =E= SUM(T, 10 *RR(T)*L(T)*(PERIODU(T))/scale1)+ scale2 ;

**  Upper and Lower Bounds: General conditions for stability

K.lo(T)         = 100;
MAT.lo(T)       = 10;
MU.lo(t)        = 100;
ML.lo(t)        = 1000;
C.lo(T)         = 20;
TOCEAN.up(T)    = 20;
TOCEAN.lo(T)    = -1;
TATM.up(t)      = 20;
miu.up(t)       = LIMMIU;
partfract("1")= 0.25372;

* First period predetermined by Kyoto Protocol
miu.fx("1")     = 0.005;

** Fix savings assumption for standardization if needed
*s.fx(t)=.22;

** Cumulative limits on carbon use at 6000 GtC
CCA.up(T) = FOSSLIM;

** Solution options
option iterlim = 99900;
option reslim = 99999;
option solprint = on;
option limrow = 0;
option limcol = 0;
model CO2 /all/;

* Optimal run
* Solution for optimal run

solve CO2 maximizing UTILITY using nlp ;
solve CO2 maximizing UTILITY using nlp ;
solve CO2 maximizing UTILITY using nlp ;
solve CO2 maximizing UTILITY using nlp ;
solve CO2 maximizing UTILITY using nlp ;
solve CO2 maximizing UTILITY using nlp ;

* Definition of opt results

Parameters
Year(t)         Date
opt_y(t)
opt_cpc(t)
opt_s(t)
opt_indem(t)
opt_sigma(t)
opt_tatm(t)
opt_mat(t)
opt_tax(t)
opt_ri(t)
opt_rr(t)
opt_al(t)
opt_forcoth(t)
opt_l(t)
opt_etree(t)
opt_yy(t)
opt_cc(t)
opt_miu(t)
opt_wem(t)
opt_ri(t)
opt_dam(t)
opt_abate(t)
opt_mcemis(t)
opt_utility
opt_scc(t) ;

Year(t)         = 2005 +10*(ord(t)-1);
opt_y(t)=y.l(t);
opt_cpc(t)=cpc.l(t);
opt_s(t)=s.l(t)     ;
opt_indem(t)= e.l(t)-etree(t);;
opt_sigma(t)=sigma(t) ;
opt_tatm(t)=tatm.l(t)  ;
opt_mat(t)=mat.l(t)     ;
opt_tax(t)=-1*ee.m(t)*1000/(kk.m(t)+.00000000001)       ;
opt_ri(t)=ri.l(t);
opt_rr(t)=rr(t)   ;
opt_al(t)=al(t)    ;
opt_forcoth(t)=forcoth(t);
opt_l(t)=l(t);
opt_etree(t)=etree(t);
opt_yy(t)=yy.m(t)     ;
opt_cc(t)=cc.m(t)      ;
opt_miu(t)=miu.l(t)     ;
opt_wem(t)= e.l(t);
opt_ri(t)=ri.l(t)         ;
opt_dam(t)= damages.l(t);
opt_abate(t) = abatecost.l(t);
opt_mcemis(t)= expcost2*cost1(t)*miu.l(t)**(expcost2-1)/sigma(t)*1000;
opt_utility=utility.l  ;
opt_scc(t)=-ee.m(t)/cc.m(t)*(1000)    ;

option decimals=6;
display opt_scc,y.l,opt_miu,ee.m, cc.m, yy.m;

execute_unload "DICE2007.gdx"
