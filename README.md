<body style="background-color:#F4F6F7;">

# Climate DICE

[Dynamic Integrated Climate-Economy (DICE)](https://en.wikipedia.org/wiki/DICE_model) model by [William Nordhaus](http://www.econ.yale.edu/~nordhaus/homepage/homepage.htm) is one of the mainstream Integrated Assessment Models (IAMs) for evaluation of *optimal* climate policy under certain assumptions. It is broadly cited and used in climate change economics discussions. Despite its simplicity, the model is a very helpful starting point to understand relations between climate change, economy, and policy with very long run effects.


<a title="By WolfgangRieger [Public domain], via Wikimedia Commons" href="https://commons.wikimedia.org/wiki/File%3APompeii_-_Osteria_della_Via_di_Mercurio_-_Dice_Players.jpg"><center><img width="512" alt="Pompeii - Osteria della Via di Mercurio - Dice Players" src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/16/Pompeii_-_Osteria_della_Via_di_Mercurio_-_Dice_Players.jpg/512px-Pompeii_-_Osteria_della_Via_di_Mercurio_-_Dice_Players.jpg"/></div> <div class="caption"></center></div>
<figcaption><center>Fig.1 Dice Players in Pompeii (wall-painting).</figcaption>
</a><left>

## Motivation
Original DICE model is written in [**GAMS**](https://www.gams.com/), an efficient and powerful language and software for formulating and solving a broad set of optimization problems. DICE model is quite simple and can be translated and solved in various alternative programming languages. This project provides such alternative formulations of DICE model in [**R**](https://www.r-project.org/), [**Julia**](https://julialang.org/), and [**Python**](https://www.python.org/). All of them are free, open source software, which also provide a vast of tools for data analysis and graphical representation.

## Disclosure
The aim of this project is not to oppose nor endorse DICE model, but rather educational, to lower barriers of access to IAMs (like DICE) by translating them into alternative, preferably free and open source software; to draw attention of broader research community to economy-climate models (and hopefully, make our chances for fortune higher than the dice-players on the Fig.1 had). 

## Original versions of DICE model in GAMS
Several versions of DICE, starting from **DICE123** (1991) to **DICE2016**, are openly available from William Nordhaus' websites (see http://www.econ.yale.edu/~nordhaus/homepage/homepage.htm and https://sites.google.com/site/williamdnordhaus/dice-rice). All of them are formulated in GAMS (there are some simplified Excel versions which do not offer optimization). The original GAMS source codes of the models are also stored in the "GAMS/" folder of this project for reference. The most recent version (*DICE2016*) has been translated to **R** and **Julia**. Earlier versions of DICE in R, as well as Python interpretation, are forthcoming. 

## DICE in R
R (see [The R Project for Statistical Computing](https://www.r-project.org/)) provides broad set of features for computing and graphics. This project aims to provide two R versions of the model. The first one is a "direct translation", standalone version with a step-by-step conversion of original GAMS code to R. This script could be helpful to those who are familiar with GAMS code of the model and/or want to understand the translation strategy, solver options, and the model itself more deeply. All variables and parameters are stored in R global environment, i.e. fully accessible for review and analysis, though this feature doesn't allow simultaneous work with several versions of DICE model family, thus models should be run sequentially to generate results for following analysis and comparison.

The second version of R code (*forthcoming*) is an R-package-style with all mechanics of optimization hidden under the hood of one function __solveDICE__(*forthcoming*). It is more suitable for comparative analysis of scenarios and results from different DICE versions. All parameters to the model are supplied as parameters to the function and can be manipulated by users.

## DICE in Julia
[Julia](https://julialang.org/) is a modern high-level programming language for numerical computing with the main focus on the code-compactness and performance. With [JuMP](https://github.com/JuliaOpt/JuMP.jl) package and [IPOPT](https://projects.coin-or.org/Ipopt) open-source solver (embedded in [ipopt](https://github.com/JuliaOpt/Ipopt.jl) package for *Julia*), the translation from GAMS is pretty straightforward. Reformulation of the problem is not needed. This version of DICE2016 in Julia is tested with [JuliaPro v06.2.1 Personal](https://juliacomputing.com/products/juliapro.html) assembly, and produces identical results with GAMS version for baseline scenario (others are not tested).

## DICE in Python
(forthcoming)

## Comparative runs of different software versions of DICE model
(forthcoming)
