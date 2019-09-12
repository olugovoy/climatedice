
# This is a draft of Julia version of Willam Nordhaus' DICE2016 model, obtained from
# the model author's website <https://sites.google.com/site/williamdnordhaus/dice-rice>,
# translated from original GAMS code by Oleg Lugovoy (https://github.com/olugovoy/).
#
# For any inconsistencies of this version with original GAMS version, please
# blame the the interpreter here: <https://github.com/olugovoy/climatedice/issues>
#
# Current version of DICE2016inJulia reads parameters from RData file ("DICE2016_params.rda"),
# solves the Model with IPOPT solver, produces basic plots of DICE variables,
# and saves results to CSV files.
# Standalone Julia version (with declaration of all the parameters in Julia),
# as well as several common DICE scenarios is in production and coming soon.
#
# Run "install.jl" file to add all required packages.

using RData
#help(load)
ds = load("DICE2016_params.rda", convert=true)
if true # PARAMETERS
    typeof(ds)
    ds["a0"]
    typeof(ds["t"])
    #println(ds["L"])
    pars = keys(ds)
    length(pars)
    typeof(pars)

    # DICE2016 Model
    t = ds["t"]
    #t = collect(1:1:100)
    l = ds["L"]
    al = ds["al"]
    gama = ds["gama"]
    k0 = ds["k0"]
    dk = ds["dk"]
    tstep = ds["tstep"]
    rr_orig = ds["rr"]
    #rr = rr_orig
    #elasmu = ds["elasmu"]
    #elasmu = 2 # as in DICE2007

    scale1 = ds["scale1"]
    scale2 = ds["scale2"]
    etree = ds["etree"]
    sigma = ds["sigma"]
    cumetree = ds["cumetree"]
    fco22x = ds["fco22x"]
    forcoth = ds["forcoth"]
    a1 = ds["a1"]
    a2 = ds["a2"]
    a3 = ds["a3"]
    cost1 = ds["cost1"]
    expcost2 = ds["expcost2"]
    pbacktime = ds["pbacktime"]
    b11 = ds["b11"]
    b12 = ds["b12"]
    b21 = ds["b21"]
    b22 = ds["b22"]
    b23 = ds["b23"]
    b32 = ds["b32"]
    b33 = ds["b33"]
    c1 = ds["c1"]
    c3 = ds["c3"]
    c4 = ds["c4"]
    t2xco2 = ds["t2xco2"]
    tatm0 = ds["tatm0"]
    tocean0 = ds["tocean0"]
    mu0 = ds["mu0"]
    mat0 = ds["mat0"]
    ml0 = ds["ml0"]
    fosslim = ds["fosslim"]
    optlrsav = ds["optlrsav"]
    miu0 = ds["miu0"]
    prstp = ds["prstp"]

    miu_up = vcat(repeat(1:1, inner = 30), fill(1.2, 70, 1))
end

using JuMP, Ipopt

dice = Model()

#if true
# Declare parameters which can be changes inside the Model
@NLparameter(dice, rr[i in t] == rr_orig[i])
@NLparameter(dice, elasmu == ds["elasmu"]) # Check TATM when rr = 0
# Declare variables
@variable(dice, 0.01 <= YGROSS[i in t] <= Inf, start = 100)
@variable(dice, 0.01 <= K[t] <= Inf)
@variable(dice, 0.01 <= C[t] <= Inf)
@variable(dice, 0.01 <= I[t] <= Inf)
@variable(dice, 0.01 <= S[t] <= .9) #, start = optlrsav
@variable(dice, 0.01 <= CPC[t] <= Inf)
#         RI(t)           Real interest rate (per annum)
@variable(dice, 0.0 <= RI[t] <= 1)#, start = 0)
#@variable(dice, RI[t], start = .05)
@variable(dice, 0 <= RI2[t] <= 1)
#         DAMAGES(t)      Damages (trillions 2005 USD per year)
@variable(dice, 0.0 <= DAMAGES[t] <= Inf)
#         DAMFRAC(t)      Damages as fraction of gross output
@variable(dice, 0.0 <= DAMFRAC[t] <= 1)
#         ABATECOST(t)    Cost of emissions reductions  (trillions 2005 USD per year)
@variable(dice, 0.0 <= ABATECOST[t] <= Inf)
#         MCABATE(t)      Marginal cost of abatement (2005$ per ton CO2)
@variable(dice, 0.01 <= MCABATE[t] <= Inf)
#         YNET(t)         Output net of damages equation (trillions 2005 USD per year)
@variable(dice, 0.01 <= YNET[t] <= Inf)
#         Y(t)            Gross world product net of abatement and damages (trillions 2005
@variable(dice, 0.01 <= Y[t] <= Inf)
# Emissions & climate
#         MIU(t)          Emission control rate GHGs
@variable(dice, 0.0 <= MIU[t] <= 1.2, start = 0.)
#@variable(dice, fill(0, 100, 1) <= MIU[t] <= MIU_UP[t], start = 0)
#         FORC(t)         Increase in radiative forcing (watts per m2 from 1900)
@variable(dice, 0.0 <= FORC[t] <= Inf, start = 1.)
#         TATM(t)         Increase temperature of atmosphere (degrees C from 1900)
@variable(dice, 0.0 <= TATM[t] <= 12., start = 1.)
#         TOCEAN(t)       Increase temperatureof lower oceans (degrees C from 1900)
@variable(dice, 0.0 <= TOCEAN[t] <= 20., start = 1.)
#         MAT(t)          Carbon concentration increase in atmosphere (GtC from 1750)
@variable(dice, 10.0 <= MAT[t] <= Inf, start = 10.)
#         MU(t)           Carbon concentration increase in shallow oceans (GtC from 1750)
@variable(dice, 100.0 <= MU[t] <= Inf, start = 100.)
#         ML(t)           Carbon concentration increase in lower oceans (GtC from 1750)
@variable(dice, 1000.0 <= ML[t] <= Inf, start = 1000.)
#         E(t)            Total CO2 emissions (GtCO2 per year)
@variable(dice, E[t])
#@variable(dice, 0.0 <= E[t] <= Inf, start = 1)
#         EIND(t)         Industrial emissions (GtCO2 per year)
@variable(dice, EIND[t], start = 0)
#@variable(dice, 0.0 <= EIND[t] <= Inf, start = 1)
#         CCA(t)          Cumulative industrial carbon emissions (GTC)
@variable(dice, 0.0 <= CCA[t] <= fosslim, start = 400.)
#         CCATOT(t)       Total carbon emissions (GtC)
@variable(dice, 1. <= CCATOT[t] <= Inf, start = 400.)
#         CPRICE(t)       Carbon price (2005$ per ton of CO2)
@variable(dice, 0.0 <= CPRICE[t] <= Inf, start = 1.)
#         PERIODU(t)      One period utility function
@variable(dice, PERIODU[t])
#         CEMUTOTPER(t)   Period utility
@variable(dice, CEMUTOTPER[t])
#         UTILITY         Welfare function;
@variable(dice, UTILITY)

#ygrosseq(t)..        YGROSS(t)      =E= (al(t)*(L(t)/1000)**(1-GAMA))*(K(t)**GAMA)
@NLconstraint(dice, YGROSSEQ[i in t], YGROSS[i] == al[i]*((l[i]/1000)^(1-gama))*K[i]^gama)
#typeof(YGROSSEQ)
#typeof(K)
#  kk(t+1)..            K(t+1)         =L= (1-dk)**tstep * K(t) + tstep * I(t);
@NLconstraint(dice, KK[i in t[1:end-1]], K[i+1] == (1-dk)^tstep * K[i] + tstep * I[i])
##JuMP.fix(K[1], k0; force=true)
#@NLconstraint(dice, KK0[1], K[1] == k0)
# cc(t)..  C(t)           =E= Y(t) - I(t)
@NLconstraint(dice, CC[i in t], C[i] == Y[i] - I[i])
# seq(t)..             I(t)           =E= S(t) * Y(t);
@NLconstraint(dice, SEQ[i in t], I[i] == S[i] * Y[i])
@NLconstraint(dice, SEQ1[i in t[(end-10):end]], S[i] == optlrsav)
#JuMP.fix(S[i=t[(end-10):end]], optlrsav) # Doesn't work this way
#  cpce(t)..            CPC(t)         =E= 1000 * C(t) / L(t);
@NLconstraint(dice, CPCE[i in t], CPC[i] == 1000 * C[i] / l[i])
#  periodueq(t)..       PERIODU(t)     =E= ((C(T)*1000/L(T))**(1-elasmu)-1)/(1-elasmu)-1;
@NLconstraint(dice, PERIODUEQ[i in t], PERIODU[i] == ((C[i]*1000/l[i])^(1-elasmu)-1)/(1-elasmu)-1)
#cemutotpereq(t)..    CEMUTOTPER(t)  =E= PERIODU(t) * L(t) * rr(t)
@NLconstraint(dice, CEMUTOTPEREQ[i in t], CEMUTOTPER[i] == PERIODU[i] * l[i] * rr[i])
# eeq(t)..             E(t)           =E= EIND(t) + etree(t);
@NLconstraint(dice, EEQ[i in t], E[i] == EIND[i] + etree[i])
# eindeq(t)..          EIND(t)        =E= sigma(t) * YGROSS(t) * (1-(MIU(t)));
@NLconstraint(dice, EINDEQ[i in t], EIND[i] == sigma[i] * YGROSS[i] * (1-(MIU[i])))
#ccacca(t+1)..        CCA(t+1)       =E= CCA(t)+ EIND(t)*5/3.666;
@NLconstraint(dice, CCACCA[i in t[1:end-1]], CCA[i+1] == CCA[i] + EIND[i] * 5/3.666)
##JuMP.fix(CCA[1], 400.; force= true)
#@NLconstraint(dice, CCACCA0[1], CCA[1] == 400)
# ccatoteq(t)..        CCATOT(t)      =E= CCA(t)+cumetree(t);
@NLconstraint(dice, CCATOTEQ[i in t], CCATOT[i] == CCA[i] + cumetree[i])
# force(t)..           FORC(t)        =E= fco22x * ((log((MAT(t)/588.000))/log(2))) + forcoth(t);
@NLconstraint(dice, FORCE[i in t], FORC[i] == fco22x * log(MAT[i]/588.000)/log(2) + forcoth[i])
# damfraceq(t) ..      DAMFRAC(t)     =E= (a1*TATM(t))+(a2*TATM(t)**a3) ;
@NLconstraint(dice, DAMFRACEQ[i in t], DAMFRAC[i] == (a1 * TATM[i]) + (a2 * TATM[i]^a3))
# dameq(t)..           DAMAGES(t)     =E= YGROSS(t) * DAMFRAC(t);
@NLconstraint(dice, DAMEQ[i in t], DAMAGES[i] == YGROSS[i] * DAMFRAC[i])
# abateeq(t)..         ABATECOST(t)   =E= YGROSS(t) * cost1(t) * (MIU(t)**expcost2);
@NLconstraint(dice, ABATEEQ[i in t], ABATECOST[i] == YGROSS[i] * cost1[i] * (MIU[i]^expcost2))
# mcabateeq(t)..       MCABATE(t)     =E= pbacktime(t) * MIU(t)**(expcost2-1);
@NLconstraint(dice, MCABATEEQ[i in t], MCABATE[i] == pbacktime[i] * MIU[i]^(expcost2-1))
# carbpriceeq(t)..     CPRICE(t)      =E= pbacktime(t) * (MIU(t))**(expcost2-1);
@NLconstraint(dice, CARBPRICEEQ[i in t], CPRICE[i] == pbacktime[i] * (MIU[i]^expcost2))
# mmat(t+1)..          MAT(t+1)       =E= MAT(t)*b11 + MU(t)*b21 + (E(t)*(5/3.666));
@NLconstraint(dice, MMAT[i in t[1:end-1]], MAT[i+1] == MAT[i]*b11 + MU[i]*b21 + E[i] * 5/3.666)
##JuMP.fix(MAT[1], mat0;force=true)
#@NLconstraint(dice, MMAT0[1], MAT[1] == mat0)
# mml(t+1)..           ML(t+1)        =E= ML(t)*b33  + MU(t)*b23;
@NLconstraint(dice, MML[i in t[1:end-1]], ML[i+1] == ML[i]*b33 + MU[i]*b23)
##JuMP.fix(ML[1], ml0;force=true)
#@NLconstraint(dice, MML0[1], ML[1] == ml0)
# mmu(t+1)..           MU(t+1)        =E= MAT(t)*b12 + MU(t)*b22 + ML(t)*b32;
@NLconstraint(dice, MMU[i in t[1:end-1]], MU[i+1] == MAT[i]*b12 + MU[i]*b22 + ML[i]*b32)
##JuMP.fix(MU[1], mu0;force=true)
#@NLconstraint(dice, MMU0[1], MU[1] == mu0)
# tatmeq(t+1)..        TATM(t+1)      =E= TATM(t) + c1 * ((FORC(t+1)-(fco22x/t2xco2)*TATM(t))-(c3*(TATM(t)-TOCEAN(t))));
@NLconstraint(dice, TATMEQ[i in t[1:end-1]], TATM[i+1] == TATM[i] +
                c1 * ((FORC[i+1]-(fco22x/t2xco2)*TATM[i])-(c3*(TATM[i]-TOCEAN[i]))))
##JuMP.fix(TATM[1], tatm0;force=true)
#@NLconstraint(dice, TATM0[1], TATM[1] == tatm0)
# toceaneq(t+1)..      TOCEAN(t+1)    =E= TOCEAN(t) + c4*(TATM(t)-TOCEAN(t));
@NLconstraint(dice, TOCEANEQ[i in t[1:end-1]], TOCEAN[i+1] == TOCEAN[i] + c4*(TATM[i]-TOCEAN[i]))
##JuMP.fix(TOCEAN[1], tocean0;force=true)
#@NLconstraint(dice, TOCEANEQ0[1], TOCEAN[1] == tocean0)
# yneteq(t)..          YNET(t)        =E= YGROSS(t)*(1-damfrac(t));
@NLconstraint(dice, YNETEQ[i in t], YNET[i] == YGROSS[i] * (1-DAMFRAC[i]))
# yy(t)..              Y(t)           =E= YNET(t) - ABATECOST(t);
@NLconstraint(dice, YY[i in t], Y[i] == YNET[i] - ABATECOST[i])
#rieq(t)..          RI(t)          =E=  (1+ (tstep*GAMA*Y(t)/K(t)-1+(1-dk)**tstep))**(1/tstep)-1 ;
@NLconstraint(dice, RIEQ2[i in t], RI2[i] == (1 + (tstep*gama*Y[i]/K[i]-1 + (1-dk)^tstep))^(1/tstep)-1)
#setvalue(RI2[100], .01)
#rieq(t+1)..          RI(t)          =E= (1+prstp) * (CPC(t+1)/CPC(t))**(elasmu/tstep) - 1;
@NLconstraint(dice, RIEQ[i in t[1:end-1]], RI[i] == (1+prstp) * (CPC[i+1]/CPC[i])^(elasmu/tstep) - 1)
##JuMP.fix(RI[t[end]], .01;force=true)
#@NLconstraint(dice, RIEQ[i in t[1:end-1]], RI[i] == 0.5)
# MIU constrains
#@NLconstraint(dice, MIU0[1], MIU[1] == miu0)
##JuMP.fix(MIU[1], miu0;force=true)
#@NLconstraint(dice, MIU1[i in t[1:30]], MIU[i] <= 1)
#@NLconstraint(dice, MIU2[i in t], MIU[i] <= 1.2)
@NLconstraint(dice, MIU1[i in t], MIU[i] <= miu_up[i])


### OBJECTIVE ###
#util..               UTILITY        =E= tstep * scale1 * sum(t,  CEMUTOTPER(t)) + scale2 ;
@NLobjective(dice, Max, tstep * scale1 * sum(CEMUTOTPER[i] for i in t) + scale2)
#print(dice)
#writeLP(dice, "dice_eq_test.txt", genericnames=false)
#end
#println("done")
# Try to solve DICE roughly
optimize!(dice, with_optimizer(Ipopt.Optimizer,
                               tol = 0.1,
                               print_level = 5,
                               print_frequency_iter = 250,
                               start_with_resto = "yes",
                               expect_infeasible_problem = "yes",
                               max_iter = 5000
                               ))


println("UTILITY = ", getobjectivevalue(dice), ", $dice_status")

# Try to resolve if solution is not Optimal
if dice_status != :Optimal
    setsolver(dice, IpoptSolver(tol = 1e-1,  print_frequency_iter = 250))
    for i = 1:5
        if dice_status != :Optimal
            dice_status = solve(dice)
            println("UTILITY = ", getobjectivevalue(dice), ", $dice_status")
        end
    end
end

# Solve with higher precision
setsolver(dice, IpoptSolver(tol = 1e-10,  print_frequency_iter = 250))
for i = 1:5
    dice_status = solve(dice)
    println("UTILITY = ", getobjectivevalue(dice), ", Status: $dice_status")
    if dice_status == :Optimal break end#i = 1e3 end
end
# println("UTILITY = ", getobjectivevalue(dice))
# dice_status

# Export results in CSV fieles
using DataFrames
using CSV

function var2csv(variable, t, folder = "Output")
    vv = eval(Symbol(variable))
    #println("vv")
    if (t == nothing) |  (t == "") | (typeof(getvalue(vv)) == Float64)
        df = DataFrame(
            l = getvalue(vv),
            m = getdual(vv),
            lo = getlowerbound(vv),
            up = getupperbound(vv))
    else
        if typeof(t) == String t = eval(Symbol(idx)) end
        df = DataFrame(
            t = t,
            l = getvalue(vv[:]),
            m = getdual(vv[:]),
            lo = getlowerbound(vv[:]),
            up = getupperbound(vv[:]))
    end
    # if t != nothing
    #     if typeof(t) == String t = eval(Symbol(idx)) end
    #     t = DataFrame(t = t)
    #     df = hcat(t, df)
    # end
    mkpath(folder)
    CSV.write(string(folder, "/", variable, ".csv"), df)
    println(string("saving.. ", folder, "/", variable, ".csv"))
end
var2csv("UTILITY", "")
#var2csv("ABATECOST", t)

VARS =["YGROSS", "K", "C", "I", "S", "CPC", "RI", "RI2",
    "DAMAGES", "DAMFRAC", "ABATECOST", "MCABATE", "YNET", "Y", "MIU",
    "FORC", "TATM", "TOCEAN", "MAT", "MU", "ML",
    "E", "EIND", "CCA", "CCATOT",
    "CPRICE", "PERIODU", "CEMUTOTPER"]

for i in VARS
    #println(i)
    var2csv(i, t)
end

using Gadfly
plot(x = t, y = getvalue(TATM[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(MIU[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(RI[:]), Geom.point, Geom.line)
plot(layer(x = t, y = getvalue(RI[:]), Geom.line),
    layer(x = t, y = getvalue(RI2[:]), Geom.point))
#
plot(layer(x = t[1:50], y = getvalue(RI[:])[1:50], Geom.line),
    layer(x = t[1:50], y = getvalue(RI2[:])[1:50], Geom.point))
#
plot(x = t, y = getvalue(CCA[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(CCATOT[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(EIND[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(E[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(MAT[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(ML[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(MU[:]), Geom.point, Geom.line)

plot(x = t, y = getvalue(S[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(K[:]), Geom.point, Geom.line)

plot(layer(x = t, y = getvalue(YGROSS[:]), Geom.line),
     layer(x = t, y = getvalue(YNET[:]), Geom.line),
     layer(x = t, y = getvalue(Y[:]), Geom.line))

plot(x = t, y = getvalue(DAMAGES[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(DAMFRAC[:]), Geom.point, Geom.line)
plot(x = t, y = getvalue(ABATECOST[:]), Geom.point, Geom.line)

plot(x = t, y = getvalue(rr[:]))
