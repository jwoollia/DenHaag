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
    ##  id sex noff       ptype          ebv sire dam        tbv f
    ##   1   1    0 -0.01896647 -0.004741616    0   0 -0.1097632 0
    ##   2   1    0 -2.97697931 -0.744244828    0   0 -1.1328339 0
    ##   3   1    0  0.76592142  0.191480355    0   0 -0.5010496 0
    ##   4   1    0  0.88931125  0.222327813    0   0  0.5368046 0
    ##   5   2    0  0.83839101  0.209597752    0   0  0.4477924 0
    ##   6   2    0  0.28563374  0.071408434    0   0  0.2438594 0
    ##   7   2    0 -0.61102920 -0.152757299    0   0  0.4680659 0
    ##   8   2    0  0.11916895  0.029792236    0   0  0.3289373 0

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
    ##  id sex noff       ptype          ebv sire dam         tbv f
    ##   1   1    0 -0.01896647 -0.004741616    0   0 -0.10976324 0
    ##   2   1    0 -2.97697931 -0.687085878    0   0 -1.13283391 0
    ##   3   1    0  0.76592142  0.191480355    0   0 -0.50104955 0
    ##   4   1    0  0.88931125  0.223780772    0   0  0.53680463 0
    ##   5   2    0  0.83839101  0.167214438    0   0  0.44779241 0
    ##   6   2    0  0.28563374  0.071408434    0   0  0.24385942 0
    ##   7   2    0 -0.61102920 -0.051762077    0   0  0.46806595 0
    ##   8   2    0  0.11916895  0.029792236    0   0  0.32893730 0
    ##   9   1    0 -0.56513725  0.086835482    4   5  0.10263844 0
    ##  10   1    0  0.46668267  0.234238328    4   5 -0.04455996 0
    ##  11   1    0  0.10360119  0.088522467    4   7  0.04063140 0
    ##  12   1    0  1.62974224 -0.083828803    2   7 -0.26051227 0
    ##  13   2    0 -0.65148242 -0.315870963    2   5 -0.59192732 0
    ##  14   2    0 -0.10756679 -0.332015808    2   7 -0.48379569 0
    ##  15   2    0  0.68091650  0.264843162    4   5  0.12334870 0
    ##  16   2    0 -1.70541715 -0.560280145    2   7 -0.18987525 0

Note that on output `my.df` has `"noff"` set to 0. Successive generations are then produced by repeated assigment of `"noff"` and calls to `add_gen()`.

Recommended contributions can be obtained for achieving a target group coancestry by setting a group coancestry, say `my.gc`, followed by a call to `oc_sel(my.df,my.tp,my.tc,my.gc)`. The parameter `my.tp` is a number which gives the number of eligible candidates for selection, which are assumed to be the **most recent** individuals. The parameter `my.tc` is the total number of offspring required in the next generation, which is only used to scale the optimum contributions to a projected number of offspring - although these projections are not integer! The function `oc_sel()` calls `a_mat()` to produce the numerator relationships among the candidates. The function prints a table listing selected parents and contributions but the function only returns `"TRUE"` or `"FALSE"` indicating the success of the algorithm. The user can then set `"noff"` guided by the recommendations.

``` r
my.gc <- 0.2
oc_sel(my.df,8,8,my.gc)
```

    ## [1] "Recommendations for 8 offspring from the most recent 8 parents"
    ##  id sex         ebv          c      noff
    ##   9   1  0.08683548 0.02382861 0.3812578
    ##  10   1  0.23423833 0.17656011 2.8249617
    ##  11   1  0.08852247 0.12867334 2.0587734
    ##  12   1 -0.08382880 0.17093794 2.7350071
    ##  13   2 -0.31587096 0.01092368 0.1747789
    ##  14   2 -0.33201581 0.09729196 1.5566713
    ##  15   2  0.26484316 0.39178436 6.2685497

    ## [1] TRUE

``` r
# NOTE due to authors inexperience of 'SWEAVE", the following may not agree with
# the recommendations!
my.df$noff[9]=1; my.df$noff[11]=4; my.df$noff[12]=3
my.df$noff[14]=3; my.df$noff[15]=2; my.df$noff[16]=3
my.df <- add_gen(my.df,0.25)
```

    ## [1] "Creating 8 offspring"
    ##  id sex noff        ptype          ebv sire dam         tbv     f
    ##   1   1    0 -0.018966465 -0.004741616    0   0 -0.10976324 0.000
    ##   2   1    0 -2.976979312 -0.586946494    0   0 -1.13283391 0.000
    ##   3   1    0  0.765921419  0.191480355    0   0 -0.50104955 0.000
    ##   4   1    0  0.889311253  0.219119258    0   0  0.53680463 0.000
    ##   5   2    0  0.838391007  0.109374468    0   0  0.44779241 0.000
    ##   6   2    0  0.285633737  0.071408434    0   0  0.24385942 0.000
    ##   7   2    0 -0.611029196  0.101555764    0   0  0.46806595 0.000
    ##   8   2    0  0.119168945  0.029792236    0   0  0.32893730 0.000
    ##   9   1    0 -0.565137252 -0.020875003    4   5  0.10263844 0.000
    ##  10   1    0  0.466682667  0.207451978    4   5 -0.04455996 0.000
    ##  11   1    0  0.103601188  0.230733853    4   7  0.04063140 0.000
    ##  12   1    0  1.629742243  0.118149631    2   7 -0.26051227 0.000
    ##  13   2    0 -0.651482418 -0.297742642    2   5 -0.59192732 0.000
    ##  14   2    0 -0.107566792 -0.200941067    2   7 -0.48379569 0.000
    ##  15   2    0  0.680916504  0.231489199    4   5  0.12334870 0.000
    ##  16   2    0 -1.705417149 -0.376606918    2   7 -0.18987525 0.000
    ##  17   1    0  1.502822889  0.103921576   12  16  0.63555289 0.250
    ##  18   1    0 -0.002848589 -0.062923969   11  16  0.05455064 0.125
    ##  19   1    0 -0.736008389  0.044701157   12  15 -0.02531604 0.000
    ##  20   1    0  0.978376206  0.289613242   12  15  0.33028066 0.000
    ##  21   2    0  0.566047591  0.093632278   11  14  0.26555832 0.125
    ##  22   2    0 -0.549280538 -0.140985676   11  16 -0.36199080 0.125
    ##  23   2    0 -1.432668915 -0.299731018    9  14 -0.84439193 0.000
    ##  24   2    0  1.152192933  0.177367327   11  14 -0.13374069 0.125

An alternative to setting a target group coancestry is to set a cost of inbreeding and use penalised contributions. Analagous to `oc_sel()`, recommended contributions for penalised contributions can be obtained by using the function `cost_f(my.df,my.tp,my.tc,my.cost)` where the first three arguments are as with `oc_sel()` and the final argument is set by the user.

``` r
my.cost <- 1000
cost_f(my.df,8,8,my.cost)
```

    ## [1] "Recommendations for 8 offspring from the most recent 8 parents"
    ## [1] "Recommendations have group coancestry 0.21719045198587"
    ##  id sex         ebv          c     noff
    ##  17   1  0.10392158 0.06538747 1.046200
    ##  18   1 -0.06292397 0.09186800 1.469888
    ##  19   1  0.04470116 0.17124981 2.739997
    ##  20   1  0.28961324 0.17149472 2.743916
    ##  21   2  0.09363228 0.08175398 1.308064
    ##  22   2 -0.14098568 0.11147033 1.783525
    ##  23   2 -0.29973102 0.22493799 3.599008
    ##  24   2  0.17736733 0.08183771 1.309403

    ## [1] TRUE

``` r
# NOTE due to authors inexperience of 'SWEAVE", the following may not agree with
# the recommendations!
my.df$noff[17]=2; my.df$noff[18]=2; my.df$noff[19]=2; my.df$noff[20]=2
my.df$noff[21]=4; my.df$noff[22]=1; my.df$noff[23]=1; my.df$noff[24]=2
my.df <- add_gen(my.df,0.25)
```

    ## [1] "Creating 8 offspring"
    ##  id sex noff        ptype           ebv sire dam         tbv      f
    ##   1   1    0 -0.018966465 -0.0047416163    0   0 -0.10976324 0.0000
    ##   2   1    0 -2.976979312 -0.6389678364    0   0 -1.13283391 0.0000
    ##   3   1    0  0.765921419  0.1914803547    0   0 -0.50104955 0.0000
    ##   4   1    0  0.889311253  0.1524379635    0   0  0.53680463 0.0000
    ##   5   2    0  0.838391007  0.1133684626    0   0  0.44779241 0.0000
    ##   6   2    0  0.285633737  0.0714084342    0   0  0.24385942 0.0000
    ##   7   2    0 -0.611029196 -0.0211408678    0   0  0.46806595 0.0000
    ##   8   2    0  0.119168945  0.0297922363    0   0  0.32893730 0.0000
    ##   9   1    0 -0.565137252 -0.0570541077    4   5  0.10263844 0.0000
    ##  10   1    0  0.466682667  0.1805859922    4   5 -0.04455996 0.0000
    ##  11   1    0  0.103601188  0.0452416961    4   7  0.04063140 0.0000
    ##  12   1    0  1.629742243  0.0438636765    2   7 -0.26051227 0.0000
    ##  13   2    0 -0.651482418 -0.3183257914    2   5 -0.59192732 0.0000
    ##  14   2    0 -0.107566792 -0.3056576388    2   7 -0.48379569 0.0000
    ##  15   2    0  0.680916504  0.2023981402    4   5  0.12334870 0.0000
    ##  16   2    0 -1.705417149 -0.5324736672    2   7 -0.18987525 0.0000
    ##  17   1    0  1.502822889 -0.0045788601   12  16  0.63555289 0.2500
    ##  18   1    0 -0.002848589 -0.3770366923   11  16  0.05455064 0.1250
    ##  19   1    0 -0.736008389  0.0266894421   12  15 -0.02531604 0.0000
    ##  20   1    0  0.978376206  0.1990561077   12  15  0.33028066 0.0000
    ##  21   2    0  0.566047591 -0.0680820469   11  14  0.26555832 0.1250
    ##  22   2    0 -0.549280538 -0.3637788839   11  16 -0.36199080 0.1250
    ##  23   2    0 -1.432668915 -0.3919094666    9  14 -0.84439193 0.0000
    ##  24   2    0  1.152192933  0.0009828418   11  14 -0.13374069 0.1250
    ##  25   1    0 -1.597539683 -0.5265882110   18  22  0.09620270 0.3125
    ##  26   1    0 -2.380771512 -0.2650884484   20  21  0.13733013 0.1250
    ##  27   1    0  1.763029987  0.2948847272   20  21  0.95530189 0.1250
    ##  28   1    0 -0.231416766 -0.0491720407   19  21 -0.51606115 0.1250
    ##  29   2    0 -1.652944152 -0.3744709359   18  24 -0.13912052 0.2500
    ##  30   2    0 -0.708282335 -0.2631581125   17  23  0.01040783 0.1250
    ##  31   2    0  0.650169995  0.0998272032   19  24  0.71717914 0.1250
    ##  32   2    0  0.312281461  0.0052471143   17  21 -0.38014876 0.1875

The function `add_ma_gen()` is a modification of add-gen to produce the next generation using maximum avoidance (minimum coancestry) mating. The function call to `add_ma_gen()` is identical to `add_gen()`. Function `add_ma_gen()` calls `a_mat()` to obtain the numerator relationships for the parents. Prior to printing the updated population data frame it prints out the group coancestry, the average inbreeding coefficient achieved and the estimate of alpha - which is expected to be negative!

``` r
my.df$noff[25]=2; my.df$noff[26]=2; my.df$noff[27]=2; my.df$noff[28]=2
my.df$noff[29]=2; my.df$noff[30]=2; my.df$noff[31]=2; my.df$noff[32]=2
my.df <- add_ma_gen(my.df,0.25)
```

    ## [1] "Creating 8 offspring with maximum avoidance"
    ## [1] "Group coancestry 0.271484375 and average offspring F 0.20703125 with alpha -0.0884718498659518"
    ##  id sex noff        ptype          ebv sire dam         tbv        f
    ##   1   1    0 -0.018966465 -0.004741616    0   0 -0.10976324 0.000000
    ##   2   1    0 -2.976979312 -0.660886825    0   0 -1.13283391 0.000000
    ##   3   1    0  0.765921419  0.191480355    0   0 -0.50104955 0.000000
    ##   4   1    0  0.889311253  0.127376123    0   0  0.53680463 0.000000
    ##   5   2    0  0.838391007  0.095131460    0   0  0.44779241 0.000000
    ##   6   2    0  0.285633737  0.071408434    0   0  0.24385942 0.000000
    ##   7   2    0 -0.611029196 -0.049884694    0   0  0.46806595 0.000000
    ##   8   2    0  0.119168945  0.029792236    0   0  0.32893730 0.000000
    ##   9   1    0 -0.565137252 -0.072747582    4   5  0.10263844 0.000000
    ##  10   1    0  0.466682667  0.162029345    4   5 -0.04455996 0.000000
    ##  11   1    0  0.103601188  0.012107364    4   7  0.04063140 0.000000
    ##  12   1    0  1.629742243 -0.015381017    2   7 -0.26051227 0.000000
    ##  13   2    0 -0.651482418 -0.335535502    2   5 -0.59192732 0.000000
    ##  14   2    0 -0.107566792 -0.375465344    2   7 -0.48379569 0.000000
    ##  15   2    0  0.680916504  0.144515707    4   5  0.12334870 0.000000
    ##  16   2    0 -1.705417149 -0.511509093    2   7 -0.18987525 0.000000
    ##  17   1    0  1.502822889 -0.019533557   12  16  0.63555289 0.250000
    ##  18   1    0 -0.002848589 -0.348420934   11  16  0.05455064 0.125000
    ##  19   1    0 -0.736008389 -0.077294845   12  15 -0.02531604 0.000000
    ##  20   1    0  0.978376206  0.094153102   12  15  0.33028066 0.000000
    ##  21   2    0  0.566047591 -0.204427954   11  14  0.26555832 0.125000
    ##  22   2    0 -0.549280538 -0.309169607   11  16 -0.36199080 0.125000
    ##  23   2    0 -1.432668915 -0.427979320    9  14 -0.84439193 0.000000
    ##  24   2    0  1.152192933 -0.084515857   11  14 -0.13374069 0.125000
    ##  25   1    0 -1.597539683 -0.362833215   18  22  0.09620270 0.312500
    ##  26   1    0 -2.380771512 -0.414952150   20  21  0.13733013 0.125000
    ##  27   1    0  1.763029987  0.083830872   20  21  0.95530189 0.125000
    ##  28   1    0 -0.231416766 -0.253336663   19  21 -0.51606115 0.125000
    ##  29   2    0 -1.652944152 -0.455982508   18  24 -0.13912052 0.250000
    ##  30   2    0 -0.708282335 -0.287499472   17  23  0.01040783 0.125000
    ##  31   2    0  0.650169995 -0.028573710   19  24  0.71717914 0.125000
    ##  32   2    0  0.312281461 -0.068739483   17  21 -0.38014876 0.187500
    ##  33   1    0  0.792419325 -0.103763496   25  32 -0.39945323 0.250000
    ##  34   1    0 -2.028072537 -0.234006107   27  31  0.56009859 0.203125
    ##  35   1    0 -1.609640272 -0.504336181   28  29 -0.39581145 0.203125
    ##  36   1    0  0.057428384 -0.157034032   27  29 -0.22153760 0.203125
    ##  37   2    0  0.044768031 -0.300826595   26  30  0.50713802 0.171875
    ##  38   2    0  1.054786182 -0.051637835   25  31  0.25765328 0.187500
    ##  39   2    0 -1.404809245 -0.385224870   26  32 -0.57723360 0.265625
    ##  40   2    0 -0.703794095 -0.325575017   28  30 -0.44170743 0.171875
