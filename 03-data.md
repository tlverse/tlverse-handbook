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


\begin{tabular}{l|l|r|r|r|r|r|r|r|r|r|r|r|r|r|r}
\hline
skim\_type & skim\_variable & n\_missing & complete\_rate & character.min & character.max & character.empty & character.n\_unique & character.whitespace & numeric.mean & numeric.sd & numeric.p0 & numeric.p25 & numeric.p50 & numeric.p75 & numeric.p100\\
\hline
character & tr & 0 & 1.0000 & 3 & 15 & 0 & 7 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & fracode & 0 & 1.0000 & 2 & 6 & 0 & 20 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & sex & 0 & 1.0000 & 4 & 6 & 0 & 2 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & momedu & 0 & 1.0000 & 12 & 15 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & hfiacat & 0 & 1.0000 & 11 & 24 & 0 & 4 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
numeric & whz & 0 & 1.0000 & NA & NA & NA & NA & NA & -0.5861 & 1.0321 & -4.67 & -1.28 & -0.6 & 0.08 & 4.97\\
\hline
numeric & month & 0 & 1.0000 & NA & NA & NA & NA & NA & 6.4547 & 3.3321 & 1.00 & 4.00 & 6.0 & 9.00 & 12.00\\
\hline
numeric & aged & 0 & 1.0000 & NA & NA & NA & NA & NA & 266.3150 & 52.1746 & 42.00 & 230.00 & 266.0 & 303.00 & 460.00\\
\hline
numeric & momage & 18 & 0.9962 & NA & NA & NA & NA & NA & 23.9059 & 5.2405 & 14.00 & 20.00 & 23.0 & 27.00 & 60.00\\
\hline
numeric & momheight & 31 & 0.9934 & NA & NA & NA & NA & NA & 150.5041 & 5.2267 & 120.65 & 147.05 & 150.6 & 154.06 & 168.00\\
\hline
numeric & Nlt18 & 0 & 1.0000 & NA & NA & NA & NA & NA & 1.6047 & 1.2473 & 0.00 & 1.00 & 1.0 & 2.00 & 10.00\\
\hline
numeric & Ncomp & 0 & 1.0000 & NA & NA & NA & NA & NA & 11.0432 & 6.3504 & 2.00 & 6.00 & 10.0 & 14.00 & 52.00\\
\hline
numeric & watmin & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.9487 & 9.4812 & 0.00 & 0.00 & 0.0 & 1.00 & 600.00\\
\hline
numeric & elec & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.5951 & 0.4909 & 0.00 & 0.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & floor & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.1067 & 0.3088 & 0.00 & 0.00 & 0.0 & 0.00 & 1.00\\
\hline
numeric & walls & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.7150 & 0.4515 & 0.00 & 0.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & roof & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.9853 & 0.1203 & 0.00 & 1.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & asset\_wardrobe & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.1672 & 0.3732 & 0.00 & 0.00 & 0.0 & 0.00 & 1.00\\
\hline
numeric & asset\_table & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.7344 & 0.4417 & 0.00 & 0.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & asset\_chair & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.7344 & 0.4417 & 0.00 & 0.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & asset\_khat & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.6132 & 0.4871 & 0.00 & 0.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & asset\_chouki & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.7813 & 0.4134 & 0.00 & 1.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & asset\_tv & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.3039 & 0.4600 & 0.00 & 0.00 & 0.0 & 1.00 & 1.00\\
\hline
numeric & asset\_refrig & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.0794 & 0.2705 & 0.00 & 0.00 & 0.0 & 0.00 & 1.00\\
\hline
numeric & asset\_bike & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.3191 & 0.4662 & 0.00 & 0.00 & 0.0 & 1.00 & 1.00\\
\hline
numeric & asset\_moto & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.0660 & 0.2484 & 0.00 & 0.00 & 0.0 & 0.00 & 1.00\\
\hline
numeric & asset\_sewmach & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.0647 & 0.2461 & 0.00 & 0.00 & 0.0 & 0.00 & 1.00\\
\hline
numeric & asset\_mobile & 0 & 1.0000 & NA & NA & NA & NA & NA & 0.8586 & 0.3485 & 0.00 & 1.00 & 1.0 & 1.00 & 1.00\\
\hline
\end{tabular}
A convenient summary of the relevant variables appears above, complete with a
sparkline visualizations describing the marginal characteristics of each
covariate. Note that the *asset* variables reflect socioeconomic status of the
study participants. Notice also the uniform distribution of the treatment groups
(with twice as many controls) -- this is, of course, by design.
