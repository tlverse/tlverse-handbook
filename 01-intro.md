# Introduction {-}

<!--
\begin{shortbox}
\Boxhead{Test}
test shortbox
\end{shortbox}

\begin{VT1}
\VH{Test}
test VT1
\end{VT1}
-->

> "One enemy of robust science is our humanity -- our appetite for
> being right, and our tendency to find patterns in noise, to see supporting
> evidence for what we already believe is true, and to ignore the facts that do
> not fit."
>
> --- @naturenews_2015

Scientific research is at a unique point in its history. The need to improve
rigor and reproducibility is greater than ever. Corroboration moves science
forward, yet there is growing alarm that results cannot be reproduced or
validated, suggesting that many discoveries may be false or less robust than
initially claimed [@baker2016there]. A failure to meet this need will result in
further declines in the rate of scientific progress, tarnishing of the
reputation of the scientific enterprise as a whole, and erosion of the public's
trust in scientific findings [@munafo2017manifesto; @naturenews2_2015].

> "The key question we want to answer when seeing the results of any scientific
> study is whether we can trust the data analysis."
>
> --- @peng2015reproducibility

Unfortunately, in its current state, the culture of statistical data analysis
enables, rather than precludes, the manner in which human bias may affect the
results of (ideally objective) data analytic efforts. A significant degree of
human bias enters into statistical data analysis efforts in the form of improper
model selection. All procedures for estimation and hypothesis testing are
derived based on a choice of a statistical model; thus, obtaining valid
estimates and statistical inference relies critically on the chosen model
containing an accurate representation of the process that generated the data.

Consider, for example, a hypothetical study in which a treatment is assigned to
a group of patients -- was the treatment assigned randomly, or were
characteristics of the individuals (i.e., "baseline covariates") taken into
account in making the treatment assignment decision? What's more, in light of
patient characteristics and heterogeneity in clinician decision-making, are
patients assigned to the treatment arm uniformly receiving the same treatment?
Any such knowledge can -- indeed, _must_ -- be incorporated into the choice of
statistical model.  Alternatively, the data could arise from an observational
study (or "quasi-experiment"), in which there is no (or very limited) control
over the treatment assignment mechanism. In such cases, available knowledge
about the data-generating process (DGP) is even more limited. In these
situations, the statistical model should contain *all* possible distributions of
the data. In practice, however, models are not selected based on scientific
knowledge available about the DGP; instead, models are often selected based on
(1) the philosophical leanings of the analyst, (2) the relative convenience of
implementation of statistical methods admissible within the choice of model, and
(3) the results of significance testing (i.e., p-values).

This practice of "cargo-cult statistics --- the ritualistic miming of statistics
rather than conscientious practice," [@stark2018cargo] is characterized by
arbitrary modeling choices, even when these choices often result in different
answers to the same research question. As opposed to its original purpose of
safeguarding the scientific process -- by providing formal techniques for
evaluating the veracity of a claim using properly designed experiments and data
collection procedures -- Statistics is increasingly often used "to aid and abet
weak science, a role it can perform well when used mechanically or
ritual[istically]" [@stark2018cargo].  The current trend in deriving scientific
discoveries by way of abusing statistical methods helps to explain the modern
epidemic of false findings from which scientific research is suffering
[@vdl2014entering].

> "We suggest that the weak statistical understanding is probably due to
> inadequate "statistics lite" education. This approach does not build up
> appropriate mathematical fundamentals and does not provide scientifically
> rigorous introduction into statistics. Hence, students' knowledge may remain
> imprecise, patchy, and prone to serious misunderstandings. What this approach
> achieves, however, is providing students with false confidence of being able
> to use inferential tools whereas they usually only interpret the p-value
> provided by black box statistical software. While this educational problem
> remains unaddressed, poor statistical practices will prevail regardless of
> what procedures and measures may be favored and/or banned by editorials."
>
> --- @szucs2017null

Our team at the University of California, Berkeley is uniquely positioned to
provide such an education. Spearheaded by Professor Mark van der Laan, and
now spreading rapidly through his students and colleagues who have greatly
enriched the field, the aptly named "Targeted Learning" paradigm emphasizes a
focus upon (i.e., "targeting of") the scientific question motivating a given
study or dataset. The philosophy of Targeted Learning runs counter to the
current cultural problem of "convenience statistics," which opens the door to
biased estimation, misleading data analytic results, and erroneous discoveries.
Targeted Learning (TL) embraces the fundamentals that formalized the field of
Statistics, notably including the dual notions that a statistical model must
represent real knowledge about the data-generating experiment and that a target
parameter (a particular feature of the data-generating probability distribution)
represents what we seek to learn from the data [@vdl2014entering]. In this way,
TL defines a ground truth and establishes a principled standard for inference,
thereby curtailing opportunities for our all-too-human biases (e.g., hindsight
bias, confirmation bias, and outcome bias) to infiltrate our efforts at
objective data analysis.

> "The key for effective classical [statistical] inference is to have
> well-defined questions and an analysis plan that tests those questions."
>
> --- @nosek2018preregistration

Inspired loosely by R.A. Fisher's influential classic _Statistical Methods for
Research Workers_ [@fisher1946statistical], this handbook aims to provide
practical training for students, researchers, industry professionals, and
academicians in the sciences (broadly considered -- whether biological,
physical, economic, or social), medicine and public health, statistics, and
numerous other allied disciplines, equipping them with the necessary knowledge
and skills to utilize the methodological developments of TL. The Targeted
Learning paradigm encompasses a principled set of techniques, united by a single
philosophy, for developing answers to queries with confidence, utilizing
advances in causal inference, state-of-the-art non/semi-parametric statistical
theory, and machine learning --- so that each and every data analysis is
realistic, reflecting appropriately what is known (and unknown) about the
process that generated the data,
<!--via pre-specified and data-adaptive machines -- nh: idk what this means-->
while remaining fully compatible with the guiding principles of computational
reproducibility.

Just as the conscientious use of modern statistical methodology is necessary to
ensure that scientific practice thrives --- robust, well-tested software plays a
critical role in allowing practitioners to access the published results of a
given scientific investigation. We concur with the view put forth by
@buckheit1995wavelab that "an article...in a scientific publication is not the
scholarship itself, it is merely advertising of the scholarship. The actual
scholarship is the complete software development environment and the complete
set of instructions which generated the figures," making the availability and
adoption of robust statistical software key to enhancing the transparency that
is an inherent (and assumed) aspect of the scientific process.

For a statistical methodology to be readily accessible in practice, it is
crucial that it is accompanied by user-friendly software
[@pullenayegum2016knowledge; @stromberg2004write]. The `tlverse` software
ecosystem, composed of a set of packages for the `R` language and environment for
statistical computing [@R], was developed to fulfill this need for the TL
methodological framework. Not only does this suite of software tools
facilitate computationally reproducible and efficient analyses, it is also a
tool for TL education. Rather than focusing on implementing a specific estimator
or a small set of related estimators, the design paradigm of the `tlverse`
ecosystem focuses on exposing the statistical framework of Targeted Learning
itself: all software packages in the `tlverse` ecosystem directly model the key
objects defined in the mathematical and theoretical framework of Targeted
Learning. What's more, the `tlverse` software packages share a core set of
design principles centered on extensibility, allowing for them to be used in
conjunction with each other and used cohesively as building blocks for
formulating sophisticated statistical analyses. For an introduction to the TL
framework, we recommend @coyle2021targeted's [recent review
paper](https://arxiv.org/abs/2006.07333).

In this handbook, the reader will embark on a journey through the `tlverse`
ecosystem. Guided by `R` programming exercises, case studies, and
intuition-building explanations, readers will learn to use this toolbox for
applying the TL statistical methodology, which will translate to real-world
causal analyses. Some preliminaries are required prior to this learning endeavor
-- for this, we provide a list of [recommended learning resources](#learn).
