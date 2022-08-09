# Example Datasets

## WASH Benefits Bangladesh Study {#wash}

The example data come from a study of the effect of water quality, sanitation,
hand washing, and nutritional interventions on child development in rural
Bangladesh (WASH Benefits Bangladesh), a cluster randomized controlled trial
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
3. handwashing with soap;
4. combined water, sanitation, and handwashing;
5. improved nutrition through counseling and provision of lipid-based nutrient
   supplements; and
6. combined water, sanitation, handwashing, and nutrition.

In the handbook, we concentrate on child growth (size for age) as the outcome of
interest. For reference, this trial was registered with ClinicalTrials.gov under
registration number NCT01590095.


```r
library(readr)
# read in data via readr::read_csv
dat <- read_csv(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  )
)
```
For instructional purposes, we start by treating the data as independent and
identically distributed (i.i.d.) random draws from a large target population. We
could account for the clustering of the data (within sampled geographic units),
but, we avoid these details in this handbook for the sake of clarity of
illustration.  Modifications of TL methodology for biased samples, repeated
measures, and related complications, are readily available.

We have 28 variables measured, of which a single variable is set to
be the outcome of interest. This outcome, $Y$, is the weight-for-height Z-score
(`whz` in `dat`); the treatment of interest, $A$, is the randomized treatment
group (`tr` in `dat`); and the adjustment set (potential baseline confounders),
$W$, consists simply of *everything else*. This results in our observed data
structure being $n$ i.i.d.  copies of $O_i = (W_i, A_i, Y_i)$, for $i = 1,
\ldots, n$.

Using the [`skimr` package](https://CRAN.R-project.org/package=skimr), we can
quickly summarize the variables measured in the WASH Benefits data set:


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
A convenient summary of the relevant variables appears above, complete with a
sparkline visualizations describing the marginal characteristics of each
covariate. Note that the *asset* variables reflect socioeconomic status of the
study participants. Notice also the uniform distribution of the treatment groups
(with twice as many controls) -- this is, of course, by design.
