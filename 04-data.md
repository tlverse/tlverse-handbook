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


```r
library(skimr)
# optionally disable sparkline graphs for PDF output
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(dat), format = "latex")
} else {
  skim(dat)
}
```


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
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-handbook/master/",
    "data/ist_sample.csv"
  )
)
```

We have 26 variables measured, and the outcome of interest, $Y$,
indicates recurrent ischemic stroke within 14 days after randomization (`DRSISC`
in `ist`); the treatment of interest, $A$, is the randomized aspirin vs. no
aspirin treatment allocation (`RXASP` in `ist`); and the adjustment set, $W$,
consists of all other variables measured at baseline. In this data, the outcome
is occasionally missing, but there is no need to create a variable indicating
this missingness (such as $\Delta$) for analyses in the `tlverse`, since it is
automatically detected when `NA` are present in the outcome. This observed data
structure can be denoted as $n$ i.i.d. copies of $O_i = (W_i, A_i, \Delta_i,
\Delta Y_i)$, for $i = 1, \ldots, n$, where $\Delta$ denotes the binary
indicator that the outcome is observed.

Like before, we can summarize the variables measured in the IST sample data set
with `skimr`:


```r
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(ist), format = "latex")
} else {
  skim(ist)
}
```


\begin{tabular}{l|l|r|r|r|r|r|r|r|r|r|r|r|r|r|r}
\hline
skim\_type & skim\_variable & n\_missing & complete\_rate & character.min & character.max & character.empty & character.n\_unique & character.whitespace & numeric.mean & numeric.sd & numeric.p0 & numeric.p25 & numeric.p50 & numeric.p75 & numeric.p100\\
\hline
character & RCONSC & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & SEX & 0 & 1.000 & 1 & 1 & 0 & 2 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RSLEEP & 0 & 1.000 & 1 & 1 & 0 & 2 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RATRIAL & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RCT & 0 & 1.000 & 1 & 1 & 0 & 2 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RVISINF & 0 & 1.000 & 1 & 1 & 0 & 2 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RHEP24 & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RASP3 & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RDEF1 & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RDEF2 & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RDEF3 & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RDEF4 & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RDEF5 & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RDEF6 & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RDEF7 & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RDEF8 & 0 & 1.000 & 1 & 1 & 0 & 3 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & STYPE & 0 & 1.000 & 3 & 4 & 0 & 5 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & RXHEP & 0 & 1.000 & 1 & 1 & 0 & 4 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
character & REGION & 0 & 1.000 & 10 & 26 & 0 & 7 & 0 & NA & NA & NA & NA & NA & NA & NA\\
\hline
numeric & RDELAY & 0 & 1.000 & NA & NA & NA & NA & NA & 20.14400 & 12.43485 & 1 & 9 & 19 & 29 & 48\\
\hline
numeric & AGE & 0 & 1.000 & NA & NA & NA & NA & NA & 71.93460 & 11.65016 & 16 & 65 & 74 & 81 & 99\\
\hline
numeric & RSBP & 0 & 1.000 & NA & NA & NA & NA & NA & 160.61560 & 27.84196 & 71 & 140 & 160 & 180 & 290\\
\hline
numeric & MISSING\_RATRIAL\_RASP3 & 0 & 1.000 & NA & NA & NA & NA & NA & 0.05000 & 0.21797 & 0 & 0 & 0 & 0 & 1\\
\hline
numeric & MISSING\_RHEP24 & 0 & 1.000 & NA & NA & NA & NA & NA & 0.01840 & 0.13441 & 0 & 0 & 0 & 0 & 1\\
\hline
numeric & RXASP & 0 & 1.000 & NA & NA & NA & NA & NA & 0.49780 & 0.50005 & 0 & 0 & 0 & 1 & 1\\
\hline
numeric & DRSISC & 10 & 0.998 & NA & NA & NA & NA & NA & 0.02365 & 0.15196 & 0 & 0 & 0 & 0 & 1\\
\hline
\end{tabular}

## NHANES I Epidemiologic Follow-up Study (NHEFS) {#NHEFS}

This data is from the National Health and Nutrition Examination Survey (NHANES)
Data I Epidemiologic Follow-up Study. More coming soon.


```r
# read in data
nhefs_data <- read_csv(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-handbook/master/",
    "data/NHEFS.csv"
  )
)
```

A snapshot of the data set is shown below:


```r
if (knitr::is_latex_output()) {
  knitr::kable(skim_no_sparks(nhefs_data), format = "latex")
} else {
  skim(nhefs_data)
}
```


\begin{tabular}{l|l|r|r|r|r|r|r|r|r|r}
\hline
skim\_type & skim\_variable & n\_missing & complete\_rate & numeric.mean & numeric.sd & numeric.p0 & numeric.p25 & numeric.p50 & numeric.p75 & numeric.p100\\
\hline
numeric & seqn & 0 & 1.00000 & 16552.36464 & 7498.91820 & 233.00000 & 10607.00000 & 20333.00000 & 2.2719e+04 & 2.5061e+04\\
\hline
numeric & qsmk & 0 & 1.00000 & 0.26274 & 0.44026 & 0.00000 & 0.00000 & 0.00000 & 1.0000e+00 & 1.0000e+00\\
\hline
numeric & death & 0 & 1.00000 & 0.19521 & 0.39649 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & yrdth & 1311 & 0.19521 & 87.56918 & 2.65941 & 83.00000 & 85.00000 & 88.00000 & 9.0000e+01 & 9.2000e+01\\
\hline
numeric & modth & 1307 & 0.19767 & 6.25776 & 3.61530 & 1.00000 & 3.00000 & 6.00000 & 1.0000e+01 & 1.2000e+01\\
\hline
numeric & dadth & 1307 & 0.19767 & 15.87267 & 8.90549 & 1.00000 & 8.00000 & 15.50000 & 2.4000e+01 & 3.1000e+01\\
\hline
numeric & sbp & 77 & 0.95273 & 128.70941 & 19.05156 & 87.00000 & 116.00000 & 126.00000 & 1.4000e+02 & 2.2900e+02\\
\hline
numeric & dbp & 81 & 0.95028 & 77.74483 & 10.63486 & 47.00000 & 70.00000 & 77.00000 & 8.5000e+01 & 1.3000e+02\\
\hline
numeric & sex & 0 & 1.00000 & 0.50952 & 0.50006 & 0.00000 & 0.00000 & 1.00000 & 1.0000e+00 & 1.0000e+00\\
\hline
numeric & age & 0 & 1.00000 & 43.91529 & 12.17043 & 25.00000 & 33.00000 & 44.00000 & 5.3000e+01 & 7.4000e+01\\
\hline
numeric & race & 0 & 1.00000 & 0.13198 & 0.33858 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & income & 62 & 0.96194 & 17.94767 & 2.66328 & 11.00000 & 17.00000 & 19.00000 & 2.0000e+01 & 2.2000e+01\\
\hline
numeric & marital & 0 & 1.00000 & 2.50338 & 1.08237 & 2.00000 & 2.00000 & 2.00000 & 2.0000e+00 & 8.0000e+00\\
\hline
numeric & school & 0 & 1.00000 & 11.13505 & 3.08960 & 0.00000 & 10.00000 & 12.00000 & 1.2000e+01 & 1.7000e+01\\
\hline
numeric & education & 0 & 1.00000 & 2.70350 & 1.19010 & 1.00000 & 2.00000 & 3.00000 & 3.0000e+00 & 5.0000e+00\\
\hline
numeric & ht & 0 & 1.00000 & 168.74096 & 9.05313 & 142.87500 & 161.78125 & 168.28125 & 1.7538e+02 & 1.9809e+02\\
\hline
numeric & wt71 & 0 & 1.00000 & 71.05213 & 15.72959 & 36.17000 & 59.65000 & 69.40000 & 7.9950e+01 & 1.6919e+02\\
\hline
numeric & wt82 & 63 & 0.96133 & 73.46922 & 16.15805 & 35.38020 & 61.68856 & 72.12119 & 8.3461e+01 & 1.3653e+02\\
\hline
numeric & wt82\_71 & 63 & 0.96133 & 2.63830 & 7.87991 & -41.28047 & -1.47840 & 2.60381 & 6.6896e+00 & 4.8538e+01\\
\hline
numeric & birthplace & 92 & 0.94352 & 31.59532 & 14.50050 & 1.00000 & 22.00000 & 34.00000 & 4.2000e+01 & 5.6000e+01\\
\hline
numeric & smokeintensity & 0 & 1.00000 & 20.55126 & 11.80375 & 1.00000 & 10.00000 & 20.00000 & 3.0000e+01 & 8.0000e+01\\
\hline
numeric & smkintensity82\_71 & 0 & 1.00000 & -4.73788 & 13.74136 & -80.00000 & -10.00000 & -1.00000 & 1.0000e+00 & 5.0000e+01\\
\hline
numeric & smokeyrs & 0 & 1.00000 & 24.87109 & 12.19807 & 1.00000 & 15.00000 & 24.00000 & 3.3000e+01 & 6.4000e+01\\
\hline
numeric & asthma & 0 & 1.00000 & 0.04850 & 0.21488 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & bronch & 0 & 1.00000 & 0.08533 & 0.27946 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & tb & 0 & 1.00000 & 0.01412 & 0.11802 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & hf & 0 & 1.00000 & 0.00491 & 0.06993 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & hbp & 0 & 1.00000 & 1.05095 & 0.95821 & 0.00000 & 0.00000 & 1.00000 & 2.0000e+00 & 2.0000e+00\\
\hline
numeric & pepticulcer & 0 & 1.00000 & 0.10374 & 0.30502 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & colitis & 0 & 1.00000 & 0.03376 & 0.18067 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & hepatitis & 0 & 1.00000 & 0.01719 & 0.13001 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & chroniccough & 0 & 1.00000 & 0.05402 & 0.22613 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & hayfever & 0 & 1.00000 & 0.08963 & 0.28573 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & diabetes & 0 & 1.00000 & 0.97974 & 0.99579 & 0.00000 & 0.00000 & 0.00000 & 2.0000e+00 & 2.0000e+00\\
\hline
numeric & polio & 0 & 1.00000 & 0.01412 & 0.11802 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & tumor & 0 & 1.00000 & 0.02333 & 0.15099 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & nervousbreak & 0 & 1.00000 & 0.02885 & 0.16744 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & alcoholpy & 0 & 1.00000 & 0.87600 & 0.33887 & 0.00000 & 1.00000 & 1.00000 & 1.0000e+00 & 2.0000e+00\\
\hline
numeric & alcoholfreq & 0 & 1.00000 & 1.92020 & 1.30714 & 0.00000 & 1.00000 & 2.00000 & 3.0000e+00 & 5.0000e+00\\
\hline
numeric & alcoholtype & 0 & 1.00000 & 2.47575 & 1.20816 & 1.00000 & 1.00000 & 3.00000 & 4.0000e+00 & 4.0000e+00\\
\hline
numeric & alcoholhowmuch & 417 & 0.74401 & 3.28713 & 2.98470 & 1.00000 & 2.00000 & 2.00000 & 4.0000e+00 & 4.8000e+01\\
\hline
numeric & pica & 0 & 1.00000 & 0.97545 & 0.99785 & 0.00000 & 0.00000 & 0.00000 & 2.0000e+00 & 2.0000e+00\\
\hline
numeric & headache & 0 & 1.00000 & 0.62983 & 0.48300 & 0.00000 & 0.00000 & 1.00000 & 1.0000e+00 & 1.0000e+00\\
\hline
numeric & otherpain & 0 & 1.00000 & 0.24616 & 0.43091 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & weakheart & 0 & 1.00000 & 0.02210 & 0.14705 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & allergies & 0 & 1.00000 & 0.06200 & 0.24123 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & nerves & 0 & 1.00000 & 0.14426 & 0.35146 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & lackpep & 0 & 1.00000 & 0.05095 & 0.21997 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & hbpmed & 0 & 1.00000 & 1.00552 & 0.98295 & 0.00000 & 0.00000 & 1.00000 & 2.0000e+00 & 2.0000e+00\\
\hline
numeric & boweltrouble & 0 & 1.00000 & 1.03499 & 0.96722 & 0.00000 & 0.00000 & 1.00000 & 2.0000e+00 & 2.0000e+00\\
\hline
numeric & wtloss & 0 & 1.00000 & 0.02578 & 0.15854 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & infection & 0 & 1.00000 & 0.14794 & 0.35515 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & active & 0 & 1.00000 & 0.65193 & 0.65274 & 0.00000 & 0.00000 & 1.00000 & 1.0000e+00 & 2.0000e+00\\
\hline
numeric & exercise & 0 & 1.00000 & 1.19521 & 0.73935 & 0.00000 & 1.00000 & 1.00000 & 2.0000e+00 & 2.0000e+00\\
\hline
numeric & birthcontrol & 0 & 1.00000 & 1.08471 & 0.94775 & 0.00000 & 0.00000 & 1.00000 & 2.0000e+00 & 2.0000e+00\\
\hline
numeric & pregnancies & 903 & 0.44567 & 3.69146 & 2.20560 & 1.00000 & 2.00000 & 3.00000 & 5.0000e+00 & 1.5000e+01\\
\hline
numeric & cholesterol & 16 & 0.99018 & 219.97396 & 45.44420 & 78.00000 & 189.00000 & 216.00000 & 2.4500e+02 & 4.1600e+02\\
\hline
numeric & hightax82 & 92 & 0.94352 & 0.16591 & 0.37212 & 0.00000 & 0.00000 & 0.00000 & 0.0000e+00 & 1.0000e+00\\
\hline
numeric & price71 & 92 & 0.94352 & 2.13875 & 0.22902 & 1.50659 & 2.03662 & 2.16797 & 2.2417e+00 & 2.6929e+00\\
\hline
numeric & price82 & 92 & 0.94352 & 1.80610 & 0.13064 & 1.45190 & 1.73999 & 1.81494 & 1.8677e+00 & 2.1030e+00\\
\hline
numeric & tax71 & 92 & 0.94352 & 1.05858 & 0.21623 & 0.52490 & 0.94495 & 1.04980 & 1.1548e+00 & 1.5225e+00\\
\hline
numeric & tax82 & 92 & 0.94352 & 0.50598 & 0.11189 & 0.21997 & 0.43994 & 0.50598 & 5.7190e-01 & 7.4792e-01\\
\hline
numeric & price71\_82 & 92 & 0.94352 & 0.33274 & 0.15504 & -0.20270 & 0.20099 & 0.33600 & 4.4379e-01 & 6.1206e-01\\
\hline
numeric & tax71\_82 & 92 & 0.94352 & 0.55261 & 0.15032 & 0.03600 & 0.46100 & 0.54395 & 6.2195e-01 & 8.8440e-01\\
\hline
\end{tabular}
