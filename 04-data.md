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
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  )
)
```
For the purposes of this handbook, we start by treating the data as independent
and identically distributed (i.i.d.) random draws from a very large target
population. We could, with available options, account for the clustering of the
data (within sampled geographic units), but, for simplification, we avoid these
details in the handbook, although modifications of our methodology for biased
samples, repeated measures, and related complications, are available.

We have 28 variables measured, of which a single variable is set to
be the outcome of interest. This outcome, $Y$, is the weight-for-height Z-score
(`whz` in `dat`); the treatment of interest, $A$, is the randomized treatment
group (`tr` in `dat`); and the adjustment set, $W$, consists simply of
*everything else*. This results in our observed data structure being $n$ i.i.d.
copies of $O_i = (W_i, A_i, Y_i)$, for $i = 1, \ldots, n$.

Using the [`skimr` package](https://CRAN.R-project.org/package=skimr), we can
quickly summarize the variables measured in the WASH Benefits data set:


\begin{tabular}{l|l|r|r|r|r|r|r|r|r|r|r|r|r|r|r}
\hline
skim\_type & skim\_variable & n\_missing & complete\_rate & character.min & character.max & character.empty & character.n\_unique & character.whitespace & numeric.mean & numeric.sd & numeric.p0 & numeric.p25 & numeric.p50 & numeric.p75 & numeric.p100\\
\hline
character & tr & 0 & 1.00000 & 3 & 15 & 0 & 7 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & fracode & 0 & 1.00000 & 2 & 6 & 0 & 20 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & sex & 0 & 1.00000 & 4 & 6 & 0 & 2 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & momedu & 0 & 1.00000 & 12 & 15 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & hfiacat & 0 & 1.00000 & 11 & 24 & 0 & 4 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
numeric & whz & 0 & 1.00000 & NA & NA & NA & NA & NA & -0.58608 & 1.03212 & -4.67 & -1.28 & -0.6 & 0.08 & 4.97\\
\hline
numeric & month & 0 & 1.00000 & NA & NA & NA & NA & NA & 6.45474 & 3.33214 & 1.00 & 4.00 & 6.0 & 9.00 & 12.00\\
\hline
numeric & aged & 0 & 1.00000 & NA & NA & NA & NA & NA & 266.31502 & 52.17465 & 42.00 & 230.00 & 266.0 & 303.00 & 460.00\\
\hline
numeric & momage & 18 & 0.99617 & NA & NA & NA & NA & NA & 23.90592 & 5.24055 & 14.00 & 20.00 & 23.0 & 27.00 & 60.00\\
\hline
numeric & momheight & 31 & 0.99340 & NA & NA & NA & NA & NA & 150.50407 & 5.22667 & 120.65 & 147.05 & 150.6 & 154.06 & 168.00\\
\hline
numeric & Nlt18 & 0 & 1.00000 & NA & NA & NA & NA & NA & 1.60469 & 1.24726 & 0.00 & 1.00 & 1.0 & 2.00 & 10.00\\
\hline
numeric & Ncomp & 0 & 1.00000 & NA & NA & NA & NA & NA & 11.04324 & 6.35044 & 2.00 & 6.00 & 10.0 & 14.00 & 52.00\\
\hline
numeric & watmin & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.94867 & 9.48125 & 0.00 & 0.00 & 0.0 & 1.00 & 600.00\\
\hline
numeric & elec & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.59510 & 0.49092 & 0.00 & 0.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & floor & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.10671 & 0.30878 & 0.00 & 0.00 & 0.0 & 0.00 & 1.00\\
\hline
numeric & walls & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.71502 & 0.45145 & 0.00 & 0.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & roof & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.98530 & 0.12035 & 0.00 & 1.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & asset\_wardrobe & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.16720 & 0.37319 & 0.00 & 0.00 & 0.0 & 0.00 & 1.00\\
\hline
numeric & asset\_table & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.73440 & 0.44170 & 0.00 & 0.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & asset\_chair & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.73440 & 0.44170 & 0.00 & 0.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & asset\_khat & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.61321 & 0.48707 & 0.00 & 0.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & asset\_chouki & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.78126 & 0.41344 & 0.00 & 1.00 & 1.0 & 1.00 & 1.00\\
\hline
numeric & asset\_tv & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.30394 & 0.46001 & 0.00 & 0.00 & 0.0 & 1.00 & 1.00\\
\hline
numeric & asset\_refrig & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.07945 & 0.27046 & 0.00 & 0.00 & 0.0 & 0.00 & 1.00\\
\hline
numeric & asset\_bike & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.31906 & 0.46616 & 0.00 & 0.00 & 0.0 & 1.00 & 1.00\\
\hline
numeric & asset\_moto & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.06603 & 0.24836 & 0.00 & 0.00 & 0.0 & 0.00 & 1.00\\
\hline
numeric & asset\_sewmach & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.06475 & 0.24611 & 0.00 & 0.00 & 0.0 & 0.00 & 1.00\\
\hline
numeric & asset\_mobile & 0 & 1.00000 & NA & NA & NA & NA & NA & 0.85857 & 0.34850 & 0.00 & 1.00 & 1.0 & 1.00 & 1.00\\
\hline
\end{tabular}
A convenient summary of the relevant variables is given just above, complete
with a small visualization describing the marginal characteristics of each
covariate. Note that the *asset* variables reflect socio-economic status of the
study participants. Notice also the uniform distribution of the treatment groups
(with twice as many controls); this is, of course, by design.
