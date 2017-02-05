DenHaag is a package of functions for simulating populations using various options for managing inbreeding. The package provides a function for establising base population `new_pop()`, and the simplest option is then to add generations with random pairing of males and females using the function `add_gen()`. Both functions output a data frame containing data on pedigree (`"id"`, `"sire"`, `"dam"` and inbreeding coefficient `"f"`), a phenotype (`"ptype"`) with true breeding value (`"tbv"`), an estimate of breeding value ("ebv"). The populations are assumed to have two sexes so the data frame also contains data on sex (`"sex"`), and also a data column called `"noff"` which is set by the user to determine the number of offspring each parent has in the next generation. The EBV is calculated within `add_gen()` using a call to `blup()`, which in turn makes a call to function `a_inv()`.

The base is set up by calling `new_pop(nn,hh)` where `nn` is the number of individuals and `hh` is the heritability of the trait. The functions in the package assume an equal number of males and females in the population, and so `new_pop()` checks that the number in the base is an even positive integer and flags an error if not. It also flags an error if the heritability is outside the open interval (0,1).

``` r
library("DenHaag", lib.loc="~/R/win-library/3.3")
# new_pop() accepts numbers ...
my.df <- new_pop(7,0.25)
```

    ## [1] "Census size is expected to be even not odd! 7!"

``` r
my.df <- new_pop(8,1.25)
```

    ## [1] "Heritability out of bounds! 1.25!"

``` r
# ... or assigned variables
my.nn <- 8; my.hh <- 0.25
my.df <- new_pop(my.nn,my.hh)
```

    ## [1] "Creating 8 individuals for base generation"
    ##  id sex noff      ptype         ebv sire dam        tbv f
    ##   1   1    0  1.1115616  0.27789041    0   0  1.4538470 0
    ##   2   1    0 -0.6444710 -0.16111776    0   0 -0.4767513 0
    ##   3   1    0  0.5282433  0.13206082    0   0 -0.1035248 0
    ##   4   1    0 -0.4727104 -0.11817761    0   0  0.4150990 0
    ##   5   2    0 -0.8769355 -0.21923388    0   0  0.5052995 0
    ##   6   2    0  2.0394407  0.50986017    0   0  0.7150695 0
    ##   7   2    0  0.2724047  0.06810118    0   0  0.0370429 0
    ##   8   2    0  1.7381912  0.43454781    0   0  0.5985120 0

Note `"noff"` is set to 0. The user selects which individuals become parents and the number of offspring for each parent by assigning `"noff"` values.

``` r
my.df$noff[2]=4; my.df$noff[4]=4; my.df$noff[5]=4; my.df$noff[7]=4
```

The next generation is then formed by calling `add_gen(my.df,my.hh)`. The output from `add_gen()` is an extended data frame. The function `add_gen()` detects errors in setting `"noff"` such as total offspring not being a even positive integer, or numbers of offspring from males not being equal to the number of offsring from females.

``` r
my.df$noff[2]=4; my.df$noff[4]=3; my.df$noff[5]=4; my.df$noff[7]=3
my.df <- add_gen(my.df,my.hh)
```

    ## [1] "Offspring numbers are expected to be even not odd! 7!"

``` r
my.df$noff[2]=4; my.df$noff[4]=4; my.df$noff[5]=4; my.df$noff[7]=3
my.df <- add_gen(my.df,my.hh)
```

    ## [1] "Numbers of male and female parents unequal in matings!"

``` r
# the heritability can be entered as a number
my.df$noff[2]=4; my.df$noff[4]=4; my.df$noff[5]=4; my.df$noff[7]=4
my.df <- add_gen(my.df,0.25)
```

    ## [1] "Creating 8 offspring"
    ##  id sex noff      ptype          ebv sire dam         tbv f
    ##   1   1    0  1.1115616  0.277890408    0   0  1.45384703 0
    ##   2   1    0 -0.6444710 -0.034343695    0   0 -0.47675133 0
    ##   3   1    0  0.5282433  0.132060817    0   0 -0.10352478 0
    ##   4   1    0 -0.4727104 -0.167030941    0   0  0.41509903 0
    ##   5   2    0 -0.8769355 -0.068519477    0   0  0.50529954 0
    ##   6   2    0  2.0394407  0.509860174    0   0  0.71506947 0
    ##   7   2    0  0.2724047 -0.004692485    0   0  0.03704290 0
    ##   8   2    0  1.7381912  0.434547807    0   0  0.59851201 0
    ##   9   1    0 -0.9563188 -0.210212732    4   7 -0.03663645 0
    ##  10   1    0 -0.5640533 -0.181529215    4   5  0.54092939 0
    ##  11   1    0  0.9407434  0.090307702    2   5  0.61727108 0
    ##  12   1    0  0.7876834  0.011576024    4   5  0.56567800 0
    ##  13   2    0 -0.1624631 -0.124159193    4   5  0.69237522 0
    ##  14   2    0  0.8863580  0.109892779    2   7  0.13421549 0
    ##  15   2    0 -1.5201154 -0.233889136    2   7  0.34568413 0
    ##  16   2    0  0.7662527  0.092734885    2   7  0.03717500 0

Note that on output `my.df` has `"noff"` set to 0. Successive generations are then produced by repeated assigment of `"noff"` and calls to `add_gen()`.

Recommended contributions can be obtained for achieving a target group coancestry by setting a group coancestry, say `my.gc`, followed by a call to `oc_sel(my.df,my.tp,my.tc,my.gc)`. The parameter `my.tp` is a number which gives the number of eligible candidates for selection, which are assumed to be the **most recent** individuals. The parameter `my.tc` is the total number of offspring required in the next generation, which is only used to scale the optimum contributions to a projected number of offspring - although these projections are not integer! The function `oc_sel()` calls `a_mat()` to produce the numerator relationships among the candidates. The function prints a table listing selected parents and contributions but the function only returns `"TRUE"` or `"FALSE"` indicating the success of the algorithm. The user can then set `"noff"` guided by the recommendations.

``` r
my.gc <- 0.2
oc_sel(my.df,8,8,my.gc)
```

    ## [1] "Recommendations for 8 offspring from the most recent 8 parents"
    ##  id sex         ebv           c       noff
    ##  11   1  0.09030770 0.221437169 3.54299471
    ##  12   1  0.01157602 0.278562831 4.45700529
    ##  13   2 -0.12415919 0.001382103 0.02211365
    ##  14   2  0.10989278 0.267062447 4.27299916
    ##  16   2  0.09273488 0.231555450 3.70488719

    ## [1] TRUE

``` r
# NOTE due to authors inexperience of 'SWEAVE", the following may not agree with
# the recommendations!
my.df$noff[9]=1; my.df$noff[11]=4; my.df$noff[12]=3
my.df$noff[14]=3; my.df$noff[15]=2; my.df$noff[16]=3
my.df <- add_gen(my.df,0.25)
```

    ## [1] "Creating 8 offspring"
    ##  id sex noff       ptype          ebv sire dam          tbv     f
    ##   1   1    0  1.11156163  0.277890408    0   0  1.453847029 0.000
    ##   2   1    0 -0.64447104  0.063665337    0   0 -0.476751331 0.000
    ##   3   1    0  0.52824327  0.132060817    0   0 -0.103524785 0.000
    ##   4   1    0 -0.47271044 -0.170069089    0   0  0.415099031 0.000
    ##   5   2    0 -0.87693550 -0.019282829    0   0  0.505299538 0.000
    ##   6   2    0  2.03944069  0.509860174    0   0  0.715069469 0.000
    ##   7   2    0  0.27240472  0.041041750    0   0  0.037042896 0.000
    ##   8   2    0  1.73819123  0.434547807    0   0  0.598512012 0.000
    ##   9   1    0 -0.95631884 -0.187532041    4   7 -0.036636447 0.000
    ##  10   1    0 -0.56405325 -0.161729858    4   5  0.540929392 0.000
    ##  11   1    0  0.94074343  0.234962920    2   5  0.617271077 0.000
    ##  12   1    0  0.78768342  0.035891547    4   5  0.565677997 0.000
    ##  13   2    0 -0.16246309 -0.104359835    4   5  0.692375218 0.000
    ##  14   2    0  0.88635799  0.159595980    2   7  0.134215495 0.000
    ##  15   2    0 -1.52011541 -0.085890714    2   7  0.345684130 0.000
    ##  16   2    0  0.76625273  0.170294494    2   7  0.037174999 0.000
    ##  17   1    0  0.05761065 -0.003742504    9  14  0.779866416 0.125
    ##  18   1    0  0.89620687  0.106601338   12  15  0.722394743 0.000
    ##  19   1    0  0.20228744  0.197994877   11  14  0.088293309 0.125
    ##  20   1    0  0.56443435  0.144521567   11  15  0.729130167 0.125
    ##  21   2    0  0.19934688  0.116843571   12  16  0.833518634 0.000
    ##  22   2    0  1.31067520  0.360921063   11  16  0.603970399 0.125
    ##  23   2    0 -0.84060325 -0.031720732   12  16  0.201733473 0.000
    ##  24   2    0 -0.07369115  0.158569364   11  14  0.009346568 0.125

An alternative to setting a target group coancestry is to set a cost of inbreeding and use penalised contributions. Analagous to `oc_sel()`, recommended contributions for penalised contributions can be obtained by using the function `cost_f(my.df,my.tp,my.tc,my.cost)` where the first three arguments are as with `oc_sel()` and the final argument is set by the user.

``` r
my.cost <- 1000
cost_f(my.df,8,8,my.cost)
```

    ## [1] "Recommendations for 8 offspring from the most recent 8 parents"
    ## [1] "Recommendations have group coancestry 0.213942324423893"
    ##  id sex          ebv          c     noff
    ##  17   1 -0.003742504 0.16981576 2.717052
    ##  18   1  0.106601338 0.16029243 2.564679
    ##  19   1  0.197994877 0.07376138 1.180182
    ##  20   1  0.144521567 0.09613043 1.538087
    ##  21   2  0.116843571 0.15385620 2.461699
    ##  22   2  0.360921063 0.08990035 1.438406
    ##  23   2 -0.031720732 0.15370763 2.459322
    ##  24   2  0.158569364 0.10253581 1.640573

    ## [1] TRUE

``` r
# NOTE due to authors inexperience of 'SWEAVE", the following may not agree with
# the recommendations!
my.df$noff[17]=2; my.df$noff[18]=2; my.df$noff[19]=2; my.df$noff[20]=2
my.df$noff[21]=4; my.df$noff[22]=1; my.df$noff[23]=1; my.df$noff[24]=2
my.df <- add_gen(my.df,0.25)
```

    ## [1] "Creating 8 offspring"
    ##  id sex noff       ptype          ebv sire dam          tbv      f
    ##   1   1    0  1.11156163  0.277890408    0   0  1.453847029 0.0000
    ##   2   1    0 -0.64447104  0.057267161    0   0 -0.476751331 0.0000
    ##   3   1    0  0.52824327  0.132060817    0   0 -0.103524785 0.0000
    ##   4   1    0 -0.47271044 -0.122813103    0   0  0.415099031 0.0000
    ##   5   2    0 -0.87693550 -0.022494369    0   0  0.505299538 0.0000
    ##   6   2    0  2.03944069  0.509860174    0   0  0.715069469 0.0000
    ##   7   2    0  0.27240472  0.085111101    0   0  0.037042896 0.0000
    ##   8   2    0  1.73819123  0.434547807    0   0  0.598512012 0.0000
    ##   9   1    0 -0.95631884 -0.113950043    4   7 -0.036636447 0.0000
    ##  10   1    0 -0.56405325 -0.142853667    4   5  0.540929392 0.0000
    ##  11   1    0  0.94074343  0.190787354    2   5  0.617271077 0.0000
    ##  12   1    0  0.78768342  0.099294487    4   5  0.565677997 0.0000
    ##  13   2    0 -0.16246309 -0.085483644    4   5  0.692375218 0.0000
    ##  14   2    0  0.88635799  0.156106723    2   7  0.134215495 0.0000
    ##  15   2    0 -1.52011541 -0.061092742    2   7  0.345684130 0.0000
    ##  16   2    0  0.76625273  0.236332346    2   7  0.037174999 0.0000
    ##  17   1    0  0.05761065  0.111669857    9  14  0.779866416 0.1250
    ##  18   1    0  0.89620687  0.175162167   12  15  0.722394743 0.0000
    ##  19   1    0  0.20228744  0.124068561   11  14  0.088293309 0.1250
    ##  20   1    0  0.56443435  0.130563157   11  15  0.729130167 0.1250
    ##  21   2    0  0.19934688  0.292357774   12  16  0.833518634 0.0000
    ##  22   2    0  1.31067520  0.408833979   11  16  0.603970399 0.1250
    ##  23   2    0 -0.84060325  0.001641232   12  16  0.201733473 0.0000
    ##  24   2    0 -0.07369115  0.058652094   11  14  0.009346568 0.1250
    ##  25   1    0  0.88799838  0.348391271   20  22  0.891376060 0.2500
    ##  26   1    0  1.94193172  0.437137856   17  21  0.609946093 0.1250
    ##  27   1    0 -0.30129316  0.008284098   17  23  0.250287505 0.1250
    ##  28   1    0  0.15207180  0.200626496   19  21  0.204757260 0.1250
    ##  29   2    0  1.21734653  0.374272337   18  21  0.807455395 0.1875
    ##  30   2    0 -0.50397995  0.114779328   20  21  0.907517735 0.1250
    ##  31   2    0 -0.71120974 -0.010784954   19  24 -0.095289601 0.3125
    ##  32   2    0 -0.35996282  0.052465245   18  24 -0.477802959 0.1250

The function `add_ma_gen()` is a modification of add-gen to produce the next generation using maximum avoidance (minimum coancestry) mating. The function call to `add_ma_gen()` is identical to `add_gen()`. Function `add_ma_gen()` calls `a_mat()` to obtain the numerator relationships for the parents. Prior to printing the updated population data frame it prints out the group coancestry, the average inbreeding coefficient achieved and the estimate of alpha - which is expected to be negative!

``` r
my.df$noff[25]=2; my.df$noff[26]=2; my.df$noff[27]=2; my.df$noff[28]=2
my.df$noff[29]=2; my.df$noff[30]=2; my.df$noff[31]=2; my.df$noff[32]=2
my.df <- add_ma_gen(my.df,0.25)
```

    ## [1] "Creating 8 offspring with maximum avoidance"
    ## [1] "Group coancestry 0.25048828125 and average offspring F 0.1640625 with alpha -0.115309446254072"
    ##  id sex noff       ptype          ebv sire dam          tbv       f
    ##   1   1    0  1.11156163  0.277890408    0   0  1.453847029 0.00000
    ##   2   1    0 -0.64447104  0.043141635    0   0 -0.476751331 0.00000
    ##   3   1    0  0.52824327  0.132060817    0   0 -0.103524785 0.00000
    ##   4   1    0 -0.47271044 -0.129088743    0   0  0.415099031 0.00000
    ##   5   2    0 -0.87693550 -0.026334029    0   0  0.505299538 0.00000
    ##   6   2    0  2.03944069  0.509860174    0   0  0.715069469 0.00000
    ##   7   2    0  0.27240472  0.068549595    0   0  0.037042896 0.00000
    ##   8   2    0  1.73819123  0.434547807    0   0  0.598512012 0.00000
    ##   9   1    0 -0.95631884 -0.134596734    4   7 -0.036636447 0.00000
    ##  10   1    0 -0.56405325 -0.147188796    4   5  0.540929392 0.00000
    ##  11   1    0  0.94074343  0.175824616    2   5  0.617271077 0.00000
    ##  12   1    0  0.78768342  0.093652392    4   5  0.565677997 0.00000
    ##  13   2    0 -0.16246309 -0.089818773    4   5  0.692375218 0.00000
    ##  14   2    0  0.88635799  0.096030002    2   7  0.134215495 0.00000
    ##  15   2    0 -1.52011541 -0.054600128    2   7  0.345684130 0.00000
    ##  16   2    0  0.76625273  0.231032018    2   7  0.037174999 0.00000
    ##  17   1    0  0.05761065  0.045969684    9  14  0.779866416 0.12500
    ##  18   1    0  0.89620687  0.190536073   12  15  0.722394743 0.00000
    ##  19   1    0  0.20228744  0.045032698   11  14  0.088293309 0.12500
    ##  20   1    0  0.56443435  0.157215912   11  15  0.729130167 0.12500
    ##  21   2    0  0.19934688  0.270309900   12  16  0.833518634 0.00000
    ##  22   2    0  1.31067520  0.435020280   11  16  0.603970399 0.12500
    ##  23   2    0 -0.84060325 -0.005251553   12  16  0.201733473 0.00000
    ##  24   2    0 -0.07369115 -0.021505021   11  14  0.009346568 0.12500
    ##  25   1    0  0.88799838  0.446004679   20  22  0.891376060 0.25000
    ##  26   1    0  1.94193172  0.330042351   17  21  0.609946093 0.12500
    ##  27   1    0 -0.30129316 -0.032831826   17  23  0.250287505 0.12500
    ##  28   1    0  0.15207180  0.106468148   19  21  0.204757260 0.12500
    ##  29   2    0  1.21734653  0.455097781   18  21  0.807455395 0.18750
    ##  30   2    0 -0.50397995  0.107046253   20  21  0.907517735 0.12500
    ##  31   2    0 -0.71120974 -0.145378089   19  24 -0.095289601 0.31250
    ##  32   2    0 -0.35996282 -0.025995577   18  24 -0.477802959 0.12500
    ##  33   1    0 -1.32442290 -0.133447610   28  32 -0.372036380 0.18750
    ##  34   1    0  0.89431519  0.146206411   27  30  0.478694586 0.15625
    ##  35   1    0  0.01095697  0.399906730   25  29  0.756313986 0.15625
    ##  36   1    0  0.54940756  0.105039902   28  32  0.015727326 0.18750
    ##  37   2    0 -0.06454976  0.074258180   26  31  0.123151700 0.15625
    ##  38   2    0  2.21429761  0.653747817   25  29  1.091035242 0.15625
    ##  39   2    0 -0.98500995 -0.092980426   27  30  0.748245575 0.15625
    ##  40   2    0 -0.90747056 -0.022852511   26  31  0.638566927 0.15625

At any time the population data frame `my.df` can be saved and reloaded.

``` r
save(my.df,file="my_df.Rda")
# ... and later ...
load("my_df.Rda")
```
