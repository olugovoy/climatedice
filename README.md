<body style="background-color:#F4F6F7;">

# Climate DICE

[Dynamic Integrated Climate-Economy model (DICE)](https://en.wikipedia.org/wiki/DICE_model) model by [William Nordhaus](http://www.econ.yale.edu/~nordhaus/homepage/homepage.htm) is one of the mainstream Integrated Assessment Models (IAMs) for evaluation of *optimal* climate policy under certain assumptions. It is broadly cited and used in climate change economics discussions. Despite its simplicity, the model is a very helpful starting point to understand relations between climate change, economy, and policy with very long run effects.


<a title="By WolfgangRieger [Public domain], via Wikimedia Commons" href="https://commons.wikimedia.org/wiki/File%3APompeii_-_Osteria_della_Via_di_Mercurio_-_Dice_Players.jpg"><center><img width="512" alt="Pompeii - Osteria della Via di Mercurio - Dice Players" src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/16/Pompeii_-_Osteria_della_Via_di_Mercurio_-_Dice_Players.jpg/512px-Pompeii_-_Osteria_della_Via_di_Mercurio_-_Dice_Players.jpg"/></div> <div class="caption"></center></div>
<figcaption><center>Fig1. Dice Players in Pompeii (wall-painting).</figcaption>
</a><left>

## Motivation
Original DICE model is written in [**GAMS**](https://www.gams.com/), an efficient and powerful language and software for formulating and solving a broad set of optimization problems. DICE model is quite simple and can be translated and solved in various alternative programming languages. This project provides such alternative formulations of DICE model in [**R**](https://www.r-project.org/), [**Julia**](https://julialang.org/), and [**Python**](https://www.python.org/). All of them are free, open source software, which also provide a vast of tools for data analysis and graphical representation.

## Disclosure
The aim of this project is not to oppose nor endorse DICE model, but rather educational, to lower barriers of access to IAMs (like DICE) by translating them into alternative, preferably free and open source software; to draw attention of broader research community to economy-climate models (and hopefully, make our chances for fortune higher than the dice-players on the figure above had). 

## Original versions of DICE model in GAMS
Several versions of DICE, starting from **DICE123** (1991) to **DICE2016**, are openly available from William Nordhaus, websites (see http://www.econ.yale.edu/~nordhaus/homepage/homepage.htm and https://sites.google.com/site/williamdnordhaus/dice-rice). All of them are formulated in GAMS (there are some simplified Excel versions which do not offer optimization). The original GAMS source codes of the models are also stored in the "GAMS/" folder of this project for reference. The most recent version (*DICE2016*) has been translated to **R** and **Julia**. Earlier versions of DICE in R, as well as Python interpretation, are forthcoming. 

## DICE in R
R (see [The R Project for Statistical Computing](https://www.r-project.org/)) provides broad set of features for computing and graphics. This project aims to provide two R versions of the model. The first one is a "direct translation", standalone version with a step-by-step conversion of original GAMS code to R. This script could be helpful to those who are familiar with GAMS code of the model and/or want to understand the translation strategy, solver options, and the model itself more deeply. All variables and parameters are stored in R global environment, i.e. fully accessible for review and analysis, though this feature doesn't allow simultaneous work with several versions of DICE model family, thus models should be run sequentially to generate results for following analysis and comparison.

The second version of R code (*forthcoming*) is an R-package-style with all mechanics of optimization hidden under the hood of one function __solveDICE()__(*forthcoming*). It is more suitable for comparative analysis of scenarios and results from different DICE versions. All parameters to the model are supplied as parameters to the function and can be manipulated by users.


### The solution strategy
The original DICE model is formulated as a system of non-linear equations. The 2016 version has 26 blocks of equations and 28 blocks of variables, "blocks" mean grouping according to sets (dimensions), i.e. time periods ( _t=100_ ) in DICE. Removing predetermined values of variables for particular time periods, the model has _2,493 SINGLE EQUATIONS_ and _2,700 SINGLE VARIABLES_. However, the problem can be reduced to only two _control_ variables which define the system solution, i.e. _2*100_ less fixed. Though the excessiveness of the algebraic formulation, while probably trading some efficiency of solution for readablity and transparency of the problem.

Unlike GAMS (and other algebraic programming languages, f.i. [AMPL](https://ampl.com/)), __R__ (to my best knowledge) doesn't offer many options for modeling macro languages. (There are some existing initiatives, which are mostly suitable for very small problems, because of low time performance.) Due to the known time issues of R (while providing a number of other benefits), the reduction of the problem dimensionality is highly desirable, hopefully without loosing clarify of the problem formulation. Therefore the strategy of the following R code is as follows:  
    - choose two control variables (emissions control rate (__MIU__) and savings rate (__S__) in DICE are already refered as control variables by W.Nordhaus);  
    - represent every GAMS-equation of DICE model as an R-function of other variables and parameters;  
    - define an objective (global welfare (__UTILITY__)) as a function of control variables (__MIU__ and __S__);  
    - maximize the objective function using solvers for NLP problems, subject to the control variables (and constraints for particular scenarios).  

This formulation of the problem doesn't require constraints on non-control variables 
for __baseline__ scenario, and can be solved with several solvers, available in base R and other packages. Additional constrains, f.i. on temperature (__TATM__) require solvers which allow non-linear constraints on variables. So far the most flexible and robust solver (for this particular problem, from several tested) is offered in [NlcOptim](https://CRAN.R-project.org/package=NlcOptim) package by Xianyan Chen and Xiangrong Yin. The solution time with arbitrary starting values varies from 1 to 10 minutes on my machines, depending on the scenario (vs. less than a minute in GAMS). The results are identical to the original GAMS version except for some divergence to the end of the model horizon (2500s), which likely can be cured by scaling the problem, ignored by now due low importance of the problem.

## DICE in Julia
[Julia](https://julialang.org/) is a modern high-level programming language for numerical computing with the main focus on the code-compactness and performance. With [JuMP](https://github.com/JuliaOpt/JuMP.jl) package and [IPOPT](https://projects.coin-or.org/Ipopt) open-source solver (embedded in [ipopt](https://github.com/JuliaOpt/Ipopt.jl) package for *Julia*), the translation from GAMS is pretty straightforward. Reformulation of the problem is not needed. The version of DICE2016 in Julia produces identical results with GAMS version for baseline scenario (others are not tested).


## DICE in Python
(forthcoming)

## Comparative runs of different software versions of DICE model
(forthcoming)
