# Robust Statistics and Reproducible Science {#robust}

<!--

test shortbox



test VT1

-->

> "One enemy of robust science is our humanity -- our appetite for
> being right, and our tendency to find patterns in noise, to see supporting
> evidence for what we already believe is true, and to ignore the facts that do
> not fit."
>
> --- @naturenews_2015

Scientific research is at a unique point in its history. The need to improve
rigor and reproducibility in our field is greater than ever; corroboration moves
science forward, yet there is growing alarm that results cannot be reproduced or
validated, suggesting the possibility that many discoveries may be false
[@baker2016there]. Consequences of not meeting this need will result in further
decline in the rate of scientific progress, the reputation of the sciences, and
the public's trust in scientific findings [@munafo2017manifesto;
@naturenews2_2015].

> "The key question we want to answer when seeing the results of any scientific
> study is whether we can trust the data analysis."
>
> --- @peng2015reproducibility

Unfortunately, in its current state, the culture of statistical data analysis
enables, rather than precludes, the manner in which human bias may affect the
results of (ideally objective) data analytic efforts. A significant degree of
human bias enters statistical analysis efforts in the form improper model
selection. All procedures for estimation and hypothesis testing are derived
based on a choice of statistical model; thus, obtaining valid estimates and
statistical inference relies critically on the chosen statistical model
containing an accurate representation of the process that generated the data.
Consider, for example, a hypothetical study in which a treatment was assigned to
a group of patients: Was the treatment assigned randomly or were characteristics
of the individuals (i.e., baseline covariates) used in making the treatment
decision? Such knowledge can should be incorporated in the statistical model.
Alternatively, the data could be from an observational study, in which there is
no control over the treatment assignment mechanism. In such cases, available
knowledge about the data-generating process (DGP) is more limited still.  If
this is the case, then the statistical model should contain *all* possible
distributions of the data. In practice, however, models are not selected based
on scientific knowledge available about the DGP; instead, models are often
selected based on (1) the philosophical leanings of the analyst, (2) the
relative convenience of implementation of statistical methods admissible within
the choice of model, and (3) the results of significance testing (i.e.,
p-values) applied within the choice of model.

This practice of "cargo-cult statistics --- the ritualistic miming of statistics
rather than conscientious practice," [@stark2018cargo] is characterized by
arbitrary modeling choices, even though these choices often result in different
answers to the same research question.  That is, "increasingly often,
[statistics] is used instead to aid and abet weak science, a role it can perform
well when used mechanically or ritually," as opposed to its original purpose of
safeguarding against weak science by providing formal techniques for evaluating
the veracity of a claim using properly collected data [@stark2018cargo]. This
presents a fundamental drive behind the epidemic of false findings from which
scientific research is suffering [@vdl2014entering].

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
spreading rapidly by many of his students and colleagues who have greatly
enriched the field, the aptly named "Targeted Learning" methodology emphasizes a
focus of (i.e., "targeting of") the scientific question at hand, running counter
to the current culture problem of "convenience statistics," which opens the door
to biased estimation, misleading analytic results, and erroneous discoveries.
Targeted Learning embraces the fundamentals that formalized the field of
statistics, notably including the notions that a statistical model must
represent real knowledge about the experiment that generated the data and that a
target parameter represents what we are seeking to learn from the data as a
feature of the distribution that generated it [@vdl2014entering]. In this way,
Targeted Learning defines a truth and establishes a principled standard for
estimation, thereby curtailing our all-too-human biases (e.g., hindsight bias,
confirmation bias, and outcome bias) from infiltrating our objective analytic
efforts.

> "The key for effective classical [statistical] inference is to have
> well-defined questions and an analysis plan that tests those questions."
>
> --- @nosek2018preregistration

This handbook aims to provide practical training to students, researchers,
industry professionals, and academicians in the sciences (whether biological,
physical, economic, or social), public health, statistics, and numerous other
fields, to equip them with the necessary knowledge and skills to utilize the the
methodological developments of Targeted Learning --- a technique that provides
tailored pre-specified machines for answering queries --- taking advantage of
estimators that are efficient, minimally biased, and that provide formal
statistical inference --- so that each and every data analysis incorporates
state-of-the-art statistical methodology, all while ensuring compatibility with
the guiding principles of computational reproducibility.

Just as the conscientious use of modern statistical methodology is necessary to
ensure that scientific practice thrives, robust, well-tested software plays a
critical role in allowing practitioners to direct access the published results
of a given scientific investigation.  In fact, "an article...in a scientific
publication is not the scholarship itself, it is merely advertising of the
scholarship. The actual scholarship is the complete software development
environment and the complete set of instructions which generated the figures,"
thus making the availability and adoption of robust statistical software key to
enhancing the transparency that is an inherent (and assumed) aspect of the
scientific process [@buckheit1995wavelab].

For a statistical methodology to be readily accessible in practice, it is
crucial that it is accompanied by user-friendly software
[@pullenayegum2016knowledge; @stromberg2004write]. The `tlverse` software
ecosystem, composed of a set of package for the `R` language and environment for
statistical computing [@R], was developed to fulfill this need for the Targeted
Learning methodological framework. Not only does this suite of software tools
facilitate computationally reproducible and efficient analyses, it is also a
tool for Targeted Learning education, since its workflow mirrors the central
aspects of the statistical methodology. In particular, the programming paradigm
central to the `tlverse` ecosystem does not focus on implementing a specific
estimator or a small set of related estimators. Instead, the focus is on
exposing the statistical framework of Targeted Learning itself --- all software
packages in the `tlverse` ecosystem directly model the key objects defined in
the mathematical and theoretical framework of Targeted Learning. What's more,
the `tlverse` software packages share a core set of design principles centered
on extensibility, allowing for them all to be used in conjunction with each
other and even used cohesively as building blocks for formulating sophisticated
statistical analyses. For an introduction to the Targeted Learning framework, we
recommend a [recent review paper](https://arxiv.org/abs/2006.07333) from
@coyle2021targeted.

In this handbook, the reader will embark on a journey through the `tlverse`
ecosystem. Guided by `R` programming exercises, case studies, and
intuition-building explanations, readers will learn to use a toolbox for
applying the Targeted Learning statistical methodology, which will translate to
real-world causal inference analyses. Some preliminaries are required prior to
this learning endeavor -- we have made available a list of [recommended learning
resources](#learn).
