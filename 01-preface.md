# Robust Statistics and Reproducible Science {#robust}

> "One enemy of robust science is our humanity — our appetite for
> being right, and our tendency to find patterns in noise, to see supporting
> evidence for what we already believe is true, and to ignore the facts that do
> not fit."
>
> --- @naturenews_2015

Scientific research is at a unique point in history. The need to improve rigor
and reproducibility in our field is greater than ever; corroboration moves
science forward, yet there is a growing alarm about results that cannot be
reproduced and that report false discoveries [@baker2016there]. Consequences of
not meeting this need will result in further decline in the rate of scientific
progression, the reputation of the sciences, and the public’s trust in its
findings [@munafo2017manifesto; @naturenews2_2015].

> "The key question we want to answer when seeing the results of any scientific
> study is whether we can trust the data analysis."
>
> --- @peng2015reproducibility

Unfortunately, at its current state the culture of data analysis and statistics
actually enables human bias through improper model selection. All hypothesis
tests and estimators are derived from statistical models, so to obtain valid
estimates and inference it is critical that the statistical model contains the
process that generated the data. Perhaps treatment was randomized or only
depended on a small number of baseline covariates; this knowledge should and
can be incorporated in the model. Alternatively, maybe the data is
observational, and there is no knowledge about the data-generating process (DGP).
If this is the case, then the statistical model should contain *all* data
distributions. In practice; however, models are not selected based on knowledge
of the DGP, instead models are often selected based on (1) the p-values they
yield, (2) their convenience of implementation, and/or (3) an analysts loyalty
to a particular model. This practice of "cargo-cult statistics --- the
ritualistic miming of statistics rather than conscientious practice,"
[@stark2018cargo] is characterized by arbitrary modeling choices, even though
these choices often result in different answers to the same research question.
That is, "increasingly often, [statistics] is used instead to aid and
abet weak science, a role it can perform well when used mechanically or
ritually," as opposed to its original purpose of safeguarding against weak
science [@stark2018cargo]. This presents a fundamental drive behind the epidemic
of false findings that scientific research is suffering from [@vdl2014entering].

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


Our team at The University of California, Berkeley, is uniquely positioned to
provide such an education. Spearheaded by Professor Mark van der Laan, and
spreading rapidly by many of his students and colleagues who have greatly
enriched the field, the aptly named "Targeted Learning" methodology targets the
scientific question at hand and is counter to the current culture of
"convenience statistics" which opens the door to biased estimation, misleading
results, and false discoveries. Targeted Learning restores the fundamentals that
formalized the field of statistics, such as the that facts that a statistical
model represents real knowledge about the experiment that generated the data,
and a target parameter represents what we are seeking to learn from the data as
a feature of the distribution that generated it [@vdl2014entering]. In this way,
Targeted Learning defines a truth and establishes a principled standard for
estimation, thereby inhibiting these all-too-human biases (e.g., hindsight bias,
confirmation bias, and outcome bias) from infiltrating analysis.

> "The key for effective classical [statistical] inference is to have
> well-defined questions and an analysis plan that tests those questions."
>
> --- @nosek2018preregistration

The objective for this handbook is to provide training to students, researchers,
industry professionals, faculty in science, public health, statistics, and other
fields to empower them with the necessary knowledge and skills to utilize the
sound methodology of Targeted Learning --- a technique that provides tailored
pre-specified machines for answering queries, so that each data analysis is
completely reproducible, and estimators are efficient, minimally biased, and
provide formal statistical inference.

Just as the conscientious use of modern statistical methodology is necessary to
ensure that scientific practice thrives, it remains critical to acknowledge the
role that robust software plays in allowing practitioners direct access to
published results. We recall that "an article...in a scientific publication is
not the scholarship itself, it is merely advertising of the scholarship. The
actual scholarship is the complete software development environment and the
complete set of instructions which generated the figures," thus making the
availability and adoption of robust statistical software key to enhancing the
transparency that is an inherent aspect of science [@buckheit1995wavelab].

For a statistical methodology to be readily accessible in practice, it is
crucial that it is accompanied by robust user-friendly software
[@pullenayegum2016knowledge; @stromberg2004write]. The `tlverse` software
ecosystem was developed to fulfill this need for the Targeted Learning
methodology. Not only does this software facilitate computationally reproducible
and efficient analyses, it is also a tool for Targeted Learning education since
its workflow mirrors that of the methodology. In particular, the `tlverse`
paradigm does not focus on implementing a specific estimator or a small set of
related estimators. Instead, the focus is on exposing the statistical framework
of Targeted Learning itself --- all `R` packages in the `tlverse` ecosystem
directly model the key objects defined in the mathematical and theoretical
framework of Targeted Learning. What's more, the `tlverse` `R` packages share a
core set of design principles centered on extensibility, allowing for them to be
used in conjunction with each other and built upon one other in a cohesive
fashion. For an introduction to Targeted Learning, we recommend the [recent
review paper](https://arxiv.org/abs/2006.07333) from @coyle2021targeted.

In this handbook, the reader will embark on a journey through the `tlverse`
ecosystem. Guided by `R` programming exercises, case studies, and
intuitive explanation readers will build a toolbox for applying the Targeted
Learning statistical methodology, which will translate to real-world causal
inference analyses. Some preliminaries are required prior to this learning
endeavor -- we have made available a list of [recommended learning
resources](#learn).
