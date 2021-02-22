# Meet the Data {#data}

## WASH Benefits Example Dataset {#wash}

The data come from a study of the effect of water quality, sanitation, hand
washing, and nutritional interventions on child development in rural Bangladesh
(WASH Benefits Bangladesh): a cluster randomized controlled trial
[@luby2018effect]. The study enrolled pregnant women in their first or second
trimester from the rural villages of Gazipur, Kishoreganj, Mymensingh, and
Tangail districts of central Bangladesh, with an average of eight women per
cluster. Groups of eight geographically adjacent clusters were block randomized,
using a random number generator, into six intervention groups (all of which
received weekly visits from a community health promoter for the first 6 months
and every 2 weeks for the next 18 months) and a double-sized control group (no
intervention or health promoter visit). The six intervention groups were:

1. chlorinated drinking water;
2. improved sanitation;
3. hand-washing with soap;
4. combined water, sanitation, and hand washing;
5. improved nutrition through counseling and provision of lipid-based nutrient
   supplements; and
6. combined water, sanitation, handwashing, and nutrition.

In the handbook, we concentrate on child growth (size for age) as the outcome of
interest. For reference, this trial was registered with ClinicalTrials.gov as
NCT01590095.


```r
library(readr)
# read in data via readr::read_csv
dat <- read_csv(
  paste0("https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
         "wash-benefits/washb_data.csv")
)
head(dat)
#> # A tibble: 6 x 28
#>     whz tr     fracode month  aged sex    momage momedu momheight hfiacat  Nlt18
#>   <dbl> <chr>  <chr>   <dbl> <dbl> <chr>   <dbl> <chr>      <dbl> <chr>    <dbl>
#> 1  0    Contr… N05265      9   268 male       30 Prima…      146. Food Se…     3
#> 2 -1.16 Contr… N05265      9   286 male       25 Prima…      149. Moderat…     2
#> 3 -1.05 Contr… N08002      9   264 male       25 Prima…      152. Food Se…     1
#> 4 -1.26 Contr… N08002      9   252 female     28 Prima…      140. Food Se…     3
#> 5 -0.59 Contr… N06531      9   336 female     19 Secon…      151. Food Se…     2
#> # … with 1 more row, and 17 more variables: Ncomp <dbl>, watmin <dbl>,
#> #   elec <dbl>, floor <dbl>, walls <dbl>, roof <dbl>, asset_wardrobe <dbl>,
#> #   asset_table <dbl>, asset_chair <dbl>, asset_khat <dbl>, asset_chouki <dbl>,
#> #   asset_tv <dbl>, asset_refrig <dbl>, asset_bike <dbl>, asset_moto <dbl>,
#> #   asset_sewmach <dbl>, asset_mobile <dbl>
```
For the purposes of this handbook, we start by treating the data as independent and
identically distributed (i.i.d.) random draws from a very large target
population. We could, with available options, account for the clustering of the
data (within sampled geographic units), but, for simplification, we avoid these
details in the handbook, although modifications of our methodology for biased
samples, repeated measures, etc., are available.

We have 28 variables measured, of which 1 variable is set to be the outcome of
interest. This outcome, $Y$, is the weight-for-height Z-score (`whz` in `dat`);
the treatment of interest, $A$, is the randomized treatment group (`tr` in
`dat`); and the adjustment set, $W$, consists simply of *everything else*. This
results in our observed data structure being $n$ i.i.d. copies of $O_i = (W_i,
A_i, Y_i)$, for $i = 1, \ldots, n$.

Using the [`skimr` package](https://CRAN.R-project.org/package=skimr), we can
quickly summarize the variables measured in the WASH Benefits data set:


```r
library(skimr)
# optionally disable sparkline graphs for PDF output
skim_no_sparks <- skim_with(numeric = sfl(hist = NULL),
                            ts = sfl(line_graph = NULL))
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(dat))
} else {
  skim(dat)
}
```


Table: (\#tab:skim_washb_data)Data summary

|                         |     |
|:------------------------|:----|
|Name                     |dat  |
|Number of rows           |4695 |
|Number of columns        |28   |
|_______________________  |     |
|Column type frequency:   |     |
|character                |5    |
|numeric                  |23   |
|________________________ |     |
|Group variables          |None |


**Variable type: character**

|skim_variable | n_missing| complete_rate| min| max| empty| n_unique| whitespace|
|:-------------|---------:|-------------:|---:|---:|-----:|--------:|----------:|
|tr            |         0|             1|   3|  15|     0|        7|          0|
|fracode       |         0|             1|   2|   6|     0|       20|          0|
|sex           |         0|             1|   4|   6|     0|        2|          0|
|momedu        |         0|             1|  12|  15|     0|        3|          0|
|hfiacat       |         0|             1|  11|  24|     0|        4|          0|


**Variable type: numeric**

|skim_variable  | n_missing| complete_rate|   mean|    sd|     p0|    p25|   p50|    p75|   p100|hist  |
|:--------------|---------:|-------------:|------:|-----:|------:|------:|-----:|------:|------:|:-----|
|whz            |         0|          1.00|  -0.59|  1.03|  -4.67|  -1.28|  -0.6|   0.08|   4.97|▁▆▇▁▁ |
|month          |         0|          1.00|   6.45|  3.33|   1.00|   4.00|   6.0|   9.00|  12.00|▇▇▅▇▇ |
|aged           |         0|          1.00| 266.32| 52.17|  42.00| 230.00| 266.0| 303.00| 460.00|▁▂▇▅▁ |
|momage         |        18|          1.00|  23.91|  5.24|  14.00|  20.00|  23.0|  27.00|  60.00|▇▇▁▁▁ |
|momheight      |        31|          0.99| 150.50|  5.23| 120.65| 147.05| 150.6| 154.06| 168.00|▁▁▆▇▁ |
|Nlt18          |         0|          1.00|   1.60|  1.25|   0.00|   1.00|   1.0|   2.00|  10.00|▇▂▁▁▁ |
|Ncomp          |         0|          1.00|  11.04|  6.35|   2.00|   6.00|  10.0|  14.00|  52.00|▇▃▁▁▁ |
|watmin         |         0|          1.00|   0.95|  9.48|   0.00|   0.00|   0.0|   1.00| 600.00|▇▁▁▁▁ |
|elec           |         0|          1.00|   0.60|  0.49|   0.00|   0.00|   1.0|   1.00|   1.00|▆▁▁▁▇ |
|floor          |         0|          1.00|   0.11|  0.31|   0.00|   0.00|   0.0|   0.00|   1.00|▇▁▁▁▁ |
|walls          |         0|          1.00|   0.72|  0.45|   0.00|   0.00|   1.0|   1.00|   1.00|▃▁▁▁▇ |
|roof           |         0|          1.00|   0.99|  0.12|   0.00|   1.00|   1.0|   1.00|   1.00|▁▁▁▁▇ |
|asset_wardrobe |         0|          1.00|   0.17|  0.37|   0.00|   0.00|   0.0|   0.00|   1.00|▇▁▁▁▂ |
|asset_table    |         0|          1.00|   0.73|  0.44|   0.00|   0.00|   1.0|   1.00|   1.00|▃▁▁▁▇ |
|asset_chair    |         0|          1.00|   0.73|  0.44|   0.00|   0.00|   1.0|   1.00|   1.00|▃▁▁▁▇ |
|asset_khat     |         0|          1.00|   0.61|  0.49|   0.00|   0.00|   1.0|   1.00|   1.00|▅▁▁▁▇ |
|asset_chouki   |         0|          1.00|   0.78|  0.41|   0.00|   1.00|   1.0|   1.00|   1.00|▂▁▁▁▇ |
|asset_tv       |         0|          1.00|   0.30|  0.46|   0.00|   0.00|   0.0|   1.00|   1.00|▇▁▁▁▃ |
|asset_refrig   |         0|          1.00|   0.08|  0.27|   0.00|   0.00|   0.0|   0.00|   1.00|▇▁▁▁▁ |
|asset_bike     |         0|          1.00|   0.32|  0.47|   0.00|   0.00|   0.0|   1.00|   1.00|▇▁▁▁▃ |
|asset_moto     |         0|          1.00|   0.07|  0.25|   0.00|   0.00|   0.0|   0.00|   1.00|▇▁▁▁▁ |
|asset_sewmach  |         0|          1.00|   0.06|  0.25|   0.00|   0.00|   0.0|   0.00|   1.00|▇▁▁▁▁ |
|asset_mobile   |         0|          1.00|   0.86|  0.35|   0.00|   1.00|   1.0|   1.00|   1.00|▁▁▁▁▇ |
A convenient summary of the relevant variables is given just above, complete
with a small visualization describing the marginal characteristics of each
covariate. Note that the *asset* variables reflect socio-economic status of the
study participants. Notice also the uniform distribution of the treatment groups
(with twice as many controls); this is, of course, by design.

## International Stroke Trial Example Dataset {#ist}

The International Stroke Trial database contains individual patient data from
the International Stroke Trial (IST), a multi-national randomized trial
conducted between 1991 and 1996 (pilot phase between 1991 and 1993) that aimed
to assess whether early administration of aspirin, heparin, both aspirin and
heparin, or neither influenced the clinical course of acute ischaemic stroke
[@sandercock1997international]. The IST dataset includes data on 19,435 patients
with acute stroke, with 99\% complete follow-up. De-identified data are
available for download at https://datashare.is.ed.ac.uk/handle/10283/128.  This
study is described in more detail in @sandercock2011international. The example
data for this handbook considers a sample of 5,000 patients and the binary
outcome of recurrent ischemic stroke within 14 days after randomization.  Also
in this example data, we ensure that we have subjects with a missing outcome.


```r
# read in data
ist <- read_csv(
  paste0("https://raw.githubusercontent.com/tlverse/tlverse-handbook/master/",
         "data/ist_sample.csv")
)
head(ist)
#> # A tibble: 6 x 26
#>   RDELAY RCONSC SEX     AGE RSLEEP RATRIAL RCT   RVISINF RHEP24 RASP3  RSBP
#>    <dbl> <chr>  <chr> <dbl> <chr>  <chr>   <chr> <chr>   <chr>  <chr> <dbl>
#> 1     46 F      F        85 N      N       N     N       Y      N       150
#> 2     33 F      M        71 Y      Y       Y     Y       N      Y       180
#> 3      6 D      M        88 N      Y       N     N       N      N       140
#> 4      8 F      F        68 Y      N       Y     Y       N      N       118
#> 5     13 F      M        60 N      N       Y     N       N      N       140
#> # … with 1 more row, and 15 more variables: RDEF1 <chr>, RDEF2 <chr>,
#> #   RDEF3 <chr>, RDEF4 <chr>, RDEF5 <chr>, RDEF6 <chr>, RDEF7 <chr>,
#> #   RDEF8 <chr>, STYPE <chr>, RXHEP <chr>, REGION <chr>,
#> #   MISSING_RATRIAL_RASP3 <dbl>, MISSING_RHEP24 <dbl>, RXASP <dbl>,
#> #   DRSISC <dbl>
```

We have 26 variables measured, and the outcome of interest, $Y$, indicates
recurrent ischemic stroke within 14 days after randomization (`DRSISC` in
`ist`); the treatment of interest, $A$, is the randomized aspirin vs. no aspirin
treatment allocation (`RXASP` in `ist`); and the adjustment set, $W$, consists
of all other variables measured at baseline. In this data, the outcome is
occasionally missing, but there is no need to create a variable indicating this
missingness (such as $\Delta$) for analyses in the `tlverse`, since it is
automatically detected when `NA` are present in the outcome. This observed data
structure can be denoted as $n$ i.i.d. copies of $O_i = (W_i, A_i, \Delta_i,
\Delta Y_i)$, for $i = 1, \ldots, n$, where $\Delta$ denotes the binary
indicator that the outcome is observed.

Like before, we can summarize the variables measured in the IST sample data set
with `skimr`:


```r
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(ist))
} else {
  skim(ist)
}
```


Table: (\#tab:skim_ist_data)Data summary

|                         |     |
|:------------------------|:----|
|Name                     |ist  |
|Number of rows           |5000 |
|Number of columns        |26   |
|_______________________  |     |
|Column type frequency:   |     |
|character                |19   |
|numeric                  |7    |
|________________________ |     |
|Group variables          |None |


**Variable type: character**

|skim_variable | n_missing| complete_rate| min| max| empty| n_unique| whitespace|
|:-------------|---------:|-------------:|---:|---:|-----:|--------:|----------:|
|RCONSC        |         0|             1|   1|   1|     0|        3|          0|
|SEX           |         0|             1|   1|   1|     0|        2|          0|
|RSLEEP        |         0|             1|   1|   1|     0|        2|          0|
|RATRIAL       |         0|             1|   1|   1|     0|        3|          0|
|RCT           |         0|             1|   1|   1|     0|        2|          0|
|RVISINF       |         0|             1|   1|   1|     0|        2|          0|
|RHEP24        |         0|             1|   1|   1|     0|        3|          0|
|RASP3         |         0|             1|   1|   1|     0|        3|          0|
|RDEF1         |         0|             1|   1|   1|     0|        3|          0|
|RDEF2         |         0|             1|   1|   1|     0|        3|          0|
|RDEF3         |         0|             1|   1|   1|     0|        3|          0|
|RDEF4         |         0|             1|   1|   1|     0|        3|          0|
|RDEF5         |         0|             1|   1|   1|     0|        3|          0|
|RDEF6         |         0|             1|   1|   1|     0|        3|          0|
|RDEF7         |         0|             1|   1|   1|     0|        3|          0|
|RDEF8         |         0|             1|   1|   1|     0|        3|          0|
|STYPE         |         0|             1|   3|   4|     0|        5|          0|
|RXHEP         |         0|             1|   1|   1|     0|        4|          0|
|REGION        |         0|             1|  10|  26|     0|        7|          0|


**Variable type: numeric**

|skim_variable         | n_missing| complete_rate|   mean|    sd| p0| p25| p50| p75| p100|hist  |
|:---------------------|---------:|-------------:|------:|-----:|--:|---:|---:|---:|----:|:-----|
|RDELAY                |         0|             1|  20.14| 12.43|  1|   9|  19|  29|   48|▇▆▆▃▂ |
|AGE                   |         0|             1|  71.93| 11.65| 16|  65|  74|  81|   99|▁▁▃▇▂ |
|RSBP                  |         0|             1| 160.62| 27.84| 71| 140| 160| 180|  290|▁▇▇▁▁ |
|MISSING_RATRIAL_RASP3 |         0|             1|   0.05|  0.22|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|MISSING_RHEP24        |         0|             1|   0.02|  0.13|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|RXASP                 |         0|             1|   0.50|  0.50|  0|   0|   0|   1|    1|▇▁▁▁▇ |
|DRSISC                |        10|             1|   0.02|  0.15|  0|   0|   0|   0|    1|▇▁▁▁▁ |

## NHANES I Epidemiologic Follow-up Study (NHEFS) {#NHEFS}

This data is from the National Health and Nutrition Examination Survey (NHANES)
Data I Epidemiologic Follow-up Study. More coming soon.


```r
# read in data
nhefs_data <- read_csv(
  paste0("https://raw.githubusercontent.com/tlverse/tlverse-handbook/master/",
         "data/NHEFS.csv")
)
head(nhefs_data)
#> # A tibble: 6 x 64
#>    seqn  qsmk death yrdth modth dadth   sbp   dbp   sex   age  race income
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
#> 1   233     0     0    NA    NA    NA   175    96     0    42     1     19
#> 2   235     0     0    NA    NA    NA   123    80     0    36     0     18
#> 3   244     0     0    NA    NA    NA   115    75     1    56     1     15
#> 4   245     0     1    85     2    14   148    78     0    68     1     15
#> 5   252     0     0    NA    NA    NA   118    77     0    40     0     18
#> # … with 1 more row, and 52 more variables: marital <dbl>, school <dbl>,
#> #   education <dbl>, ht <dbl>, wt71 <dbl>, wt82 <dbl>, wt82_71 <dbl>,
#> #   birthplace <dbl>, smokeintensity <dbl>, smkintensity82_71 <dbl>,
#> #   smokeyrs <dbl>, asthma <dbl>, bronch <dbl>, tb <dbl>, hf <dbl>, hbp <dbl>,
#> #   pepticulcer <dbl>, colitis <dbl>, hepatitis <dbl>, chroniccough <dbl>,
#> #   hayfever <dbl>, diabetes <dbl>, polio <dbl>, tumor <dbl>,
#> #   nervousbreak <dbl>, alcoholpy <dbl>, alcoholfreq <dbl>, alcoholtype <dbl>,
#> #   alcoholhowmuch <dbl>, pica <dbl>, headache <dbl>, otherpain <dbl>,
#> #   weakheart <dbl>, allergies <dbl>, nerves <dbl>, lackpep <dbl>,
#> #   hbpmed <dbl>, boweltrouble <dbl>, wtloss <dbl>, infection <dbl>,
#> #   active <dbl>, exercise <dbl>, birthcontrol <dbl>, pregnancies <dbl>,
#> #   cholesterol <dbl>, hightax82 <dbl>, price71 <dbl>, price82 <dbl>,
#> #   tax71 <dbl>, tax82 <dbl>, price71_82 <dbl>, tax71_82 <dbl>
```

A snapshot of the data set is shown below:


```r
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(nhefs_data))
} else {
  skim(nhefs_data)
}
```


Table: (\#tab:skim_nhefs_data)Data summary

|                         |           |
|:------------------------|:----------|
|Name                     |nhefs_data |
|Number of rows           |1629       |
|Number of columns        |64         |
|_______________________  |           |
|Column type frequency:   |           |
|numeric                  |64         |
|________________________ |           |
|Group variables          |None       |


**Variable type: numeric**

|skim_variable     | n_missing| complete_rate|     mean|      sd|     p0|      p25|      p50|      p75|     p100|hist  |
|:-----------------|---------:|-------------:|--------:|-------:|------:|--------:|--------:|--------:|--------:|:-----|
|seqn              |         0|          1.00| 16552.36| 7498.92| 233.00| 10607.00| 20333.00| 22719.00| 25061.00|▂▂▂▂▇ |
|qsmk              |         0|          1.00|     0.26|    0.44|   0.00|     0.00|     0.00|     1.00|     1.00|▇▁▁▁▃ |
|death             |         0|          1.00|     0.20|    0.40|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▂ |
|yrdth             |      1311|          0.20|    87.57|    2.66|  83.00|    85.00|    88.00|    90.00|    92.00|▅▇▆▆▆ |
|modth             |      1307|          0.20|     6.26|    3.62|   1.00|     3.00|     6.00|    10.00|    12.00|▇▅▅▃▇ |
|dadth             |      1307|          0.20|    15.87|    8.91|   1.00|     8.00|    15.50|    24.00|    31.00|▇▇▇▆▇ |
|sbp               |        77|          0.95|   128.71|   19.05|  87.00|   116.00|   126.00|   140.00|   229.00|▃▇▂▁▁ |
|dbp               |        81|          0.95|    77.74|   10.63|  47.00|    70.00|    77.00|    85.00|   130.00|▁▇▅▁▁ |
|sex               |         0|          1.00|     0.51|    0.50|   0.00|     0.00|     1.00|     1.00|     1.00|▇▁▁▁▇ |
|age               |         0|          1.00|    43.92|   12.17|  25.00|    33.00|    44.00|    53.00|    74.00|▇▆▇▅▂ |
|race              |         0|          1.00|     0.13|    0.34|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|income            |        62|          0.96|    17.95|    2.66|  11.00|    17.00|    19.00|    20.00|    22.00|▂▂▂▇▅ |
|marital           |         0|          1.00|     2.50|    1.08|   2.00|     2.00|     2.00|     2.00|     8.00|▇▁▁▁▁ |
|school            |         0|          1.00|    11.14|    3.09|   0.00|    10.00|    12.00|    12.00|    17.00|▁▁▅▇▂ |
|education         |         0|          1.00|     2.70|    1.19|   1.00|     2.00|     3.00|     3.00|     5.00|▃▅▇▂▂ |
|ht                |         0|          1.00|   168.74|    9.05| 142.88|   161.78|   168.28|   175.38|   198.09|▁▇▇▅▁ |
|wt71              |         0|          1.00|    71.05|   15.73|  36.17|    59.65|    69.40|    79.95|   169.19|▅▇▂▁▁ |
|wt82              |        63|          0.96|    73.47|   16.16|  35.38|    61.69|    72.12|    83.46|   136.53|▂▇▅▁▁ |
|wt82_71           |        63|          0.96|     2.64|    7.88| -41.28|    -1.48|     2.60|     6.69|    48.54|▁▁▇▁▁ |
|birthplace        |        92|          0.94|    31.60|   14.50|   1.00|    22.00|    34.00|    42.00|    56.00|▃▅▆▇▅ |
|smokeintensity    |         0|          1.00|    20.55|   11.80|   1.00|    10.00|    20.00|    30.00|    80.00|▆▇▂▁▁ |
|smkintensity82_71 |         0|          1.00|    -4.74|   13.74| -80.00|   -10.00|    -1.00|     1.00|    50.00|▁▁▇▇▁ |
|smokeyrs          |         0|          1.00|    24.87|   12.20|   1.00|    15.00|    24.00|    33.00|    64.00|▅▇▆▃▁ |
|asthma            |         0|          1.00|     0.05|    0.21|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|bronch            |         0|          1.00|     0.09|    0.28|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|tb                |         0|          1.00|     0.01|    0.12|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|hf                |         0|          1.00|     0.00|    0.07|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|hbp               |         0|          1.00|     1.05|    0.96|   0.00|     0.00|     1.00|     2.00|     2.00|▇▁▁▁▇ |
|pepticulcer       |         0|          1.00|     0.10|    0.31|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|colitis           |         0|          1.00|     0.03|    0.18|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|hepatitis         |         0|          1.00|     0.02|    0.13|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|chroniccough      |         0|          1.00|     0.05|    0.23|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|hayfever          |         0|          1.00|     0.09|    0.29|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|diabetes          |         0|          1.00|     0.98|    1.00|   0.00|     0.00|     0.00|     2.00|     2.00|▇▁▁▁▇ |
|polio             |         0|          1.00|     0.01|    0.12|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|tumor             |         0|          1.00|     0.02|    0.15|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|nervousbreak      |         0|          1.00|     0.03|    0.17|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|alcoholpy         |         0|          1.00|     0.88|    0.34|   0.00|     1.00|     1.00|     1.00|     2.00|▁▁▇▁▁ |
|alcoholfreq       |         0|          1.00|     1.92|    1.31|   0.00|     1.00|     2.00|     3.00|     5.00|▇▇▅▃▁ |
|alcoholtype       |         0|          1.00|     2.48|    1.21|   1.00|     1.00|     3.00|     4.00|     4.00|▇▂▁▇▆ |
|alcoholhowmuch    |       417|          0.74|     3.29|    2.98|   1.00|     2.00|     2.00|     4.00|    48.00|▇▁▁▁▁ |
|pica              |         0|          1.00|     0.98|    1.00|   0.00|     0.00|     0.00|     2.00|     2.00|▇▁▁▁▇ |
|headache          |         0|          1.00|     0.63|    0.48|   0.00|     0.00|     1.00|     1.00|     1.00|▅▁▁▁▇ |
|otherpain         |         0|          1.00|     0.25|    0.43|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▂ |
|weakheart         |         0|          1.00|     0.02|    0.15|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|allergies         |         0|          1.00|     0.06|    0.24|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|nerves            |         0|          1.00|     0.14|    0.35|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▂ |
|lackpep           |         0|          1.00|     0.05|    0.22|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|hbpmed            |         0|          1.00|     1.01|    0.98|   0.00|     0.00|     1.00|     2.00|     2.00|▇▁▁▁▇ |
|boweltrouble      |         0|          1.00|     1.03|    0.97|   0.00|     0.00|     1.00|     2.00|     2.00|▇▁▁▁▇ |
|wtloss            |         0|          1.00|     0.03|    0.16|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▁ |
|infection         |         0|          1.00|     0.15|    0.36|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▂ |
|active            |         0|          1.00|     0.65|    0.65|   0.00|     0.00|     1.00|     1.00|     2.00|▇▁▇▁▂ |
|exercise          |         0|          1.00|     1.20|    0.74|   0.00|     1.00|     1.00|     2.00|     2.00|▃▁▇▁▇ |
|birthcontrol      |         0|          1.00|     1.08|    0.95|   0.00|     0.00|     1.00|     2.00|     2.00|▆▁▂▁▇ |
|pregnancies       |       903|          0.45|     3.69|    2.21|   1.00|     2.00|     3.00|     5.00|    15.00|▇▅▂▁▁ |
|cholesterol       |        16|          0.99|   219.97|   45.44|  78.00|   189.00|   216.00|   245.00|   416.00|▁▇▇▂▁ |
|hightax82         |        92|          0.94|     0.17|    0.37|   0.00|     0.00|     0.00|     0.00|     1.00|▇▁▁▁▂ |
|price71           |        92|          0.94|     2.14|    0.23|   1.51|     2.04|     2.17|     2.24|     2.69|▂▂▇▅▁ |
|price82           |        92|          0.94|     1.81|    0.13|   1.45|     1.74|     1.81|     1.87|     2.10|▁▂▇▅▂ |
|tax71             |        92|          0.94|     1.06|    0.22|   0.52|     0.94|     1.05|     1.15|     1.52|▁▂▇▂▂ |
|tax82             |        92|          0.94|     0.51|    0.11|   0.22|     0.44|     0.51|     0.57|     0.75|▁▃▇▇▂ |
|price71_82        |        92|          0.94|     0.33|    0.16|  -0.20|     0.20|     0.34|     0.44|     0.61|▁▂▆▇▅ |
|tax71_82          |        92|          0.94|     0.55|    0.15|   0.04|     0.46|     0.54|     0.62|     0.88|▁▂▇▆▃ |
