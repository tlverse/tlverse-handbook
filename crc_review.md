<!--
> I usually begin by reviewing a proposal. Most authors find these comments
> useful and revise their proposal, and we usually decide whether to offer
> a contract at this stage. If we agree on a contract, I review chapters as
> written in sets of about 100 pages. You want to learn early if you have the
> essential topics, are writing at the intended level, and have enough
> motivation and examples so you can adjust your style if necessary.
-->

# CHAPMAN & HALL/CRC STATISTICS BOOK PROPOSAL

## TITLE AND AUTHOR/EDITOR(S)

1. Provisional title of your book.

   * The Hitchhiker's Guide to the tlverse: A Targeted Learning Practitioner's
      Handbook

2. Names, titles, affiliations, and email addresses for all authors/editors

      * Mark van der Laan, UC Berkeley, laan@berkeley.edu
      * Alan Hubbard, UC Berkeley, hubbard@berkeley.edu
      * Jeremy Coyle, UC Berkeley, jeremy.coyle@gmail.com
      * Nima Hejazi, UC Berkeley, nhejazi@berkeley.edu
      * Ivana Malenica, UC Berkeley, imalenica@berkeley.edu
      * Rachael Phillips, UC Berkeley, rachaelvphillips@berkeley.edu


## CONTENTS

3. Please include a complete table of contents including chapter and section
    headings.

     - The Roadmap for Targeted Learning

           - Learning Objectives
           - Introduction
           - The Roadmap
           - Summary of the Roadmap
           - Causal Target Parameters

     - Welcome to the [`tlverse`](https://tlverse.org)

           - Learning Objectives
           - What is the `tlverse`?
           - `tlverse` Components
           - Installation

     - Example Datasets and Case Studies

           - WASH Benefits Example Dataset
           - International Stroke Trial Example Dataset
           - Veterans' Administration Lung Cancer Trial Dataset

    - Cross-validation

          - Learning Objectives?
          - Why Split Your Sample?
          - Cross-validating arbitrary models
          - Exercises

    - Super (Machine) Learning

          - Learning Objectives
          - Motivation
          - Introduction
          - `sl3` "Microwave Dinner" Implementation
          - Cross-validated Super Learner
          - Variable Importance Measures with `sl3`
          - Exercises
          - Concluding Remarks
          - Appendix

    - The TMLE Framework

          - Learning Objectives
          - Introduction
          - Easy-Bake Example: `tmle3` for ATE
          - `tmle3` Components
          - Fitting `tmle3` with Multiple Parameters
          - Exercises
          - Summary

    - Optimal Individualized Treatments Regimes

          - Learning Objectives
          - Introduction to Optimal Individualized Interventions
          - Data Structure and Notation
          - Defining the Causal Effect of an Optimal Individualized Intervention
          - Interpreting the Causal Effect of Optimal Individualized Interventions
          - Evaluating the Causal Effect of an OIT with Binary Treatment
          - Evaluating the Causal Effect of an Optimal ITR with Categorical
               Treatment
          - Extensions to Causal Effect of an OIT
          - Variable Importance Analysis with OIT
          - Real-World Data and `tmle3mopttx`
          - Exercises

    - Stochastic Treatment Regimes

          - Learning Objectives
          - Introduction to Stochastic Interventions
          - Data Structure and Notation
          - Defining the Causal Effect of a Stochastic Intervention
          - Interpreting the Causal Effect of a Stochastic Intervention
          - Evaluating the Causal Effect of a Stochastic Intervention
          - Extensions: Variable Importance Analysis with Stochastic Interventions
          - Exercises

    - Appendix: A Primer on the `R6` Class System

          - Classes, Fields, and Methods
          - Object Oriented Programming: Python and R


## SUBJECT/AUDIENCE

4. Please describe in detail the subject of your book. Why will this book be
    important, who will find it useful, and what is new? What background will
    you assume and can you name any books that will suffice?

  *The Hitchhiker's Guide to the `tlverse`: A Targeted Learning Practitioner's
  Handbook* is an open-source and fully-reproducible handbook for applying the
  Targeted Learning framework and statistical methodology in practice using the
  `tlverse` software ecosystem (https://github.com/tlverse). The materials have
  been developed, tested, and refined across a series of short courses focused
  on the applications of Targeted Learning, all delivered at major academic
  conferences.

  The contents of this handbook are meant to serve as a reference guide for
  applying the Targeted Learning methodology in research settings, throughout
  the biomedical, health, and social sciences. The book will empower all manner
  of researchers to use this state-of-the-art methodology in rigorously
  developing statistical answers to scientific questions grounded in causal
  inference. Each section introduces a set of distinct causal questions, each
  motivated by a set of case studies, alongside statistical methodology and
  software for assessing the causal claims of interest.

  Necessary background for engaging with the handbook includes intermediate
  knowledge of statistical concepts and methodology (e.g., null hypotheses,
  hypothesis testing, confidence intervals), a basic working knowledge of
  statistical estimation and machine learning (e.g., what is a linear model?,
  what are prediction algorithms?), and some background in R programming. The
  handbook is intended to be self-contained -- as such, much background in
  causal inference and Targeted Learning are included. Introductory reference
  material will be provided in an appendix.

5. Will your book be primarily a reference/monograph or textbook? If a
    textbook, for which courses will it be the primary text and at what level
    is the course taught? Will you include exercises sets and supply a
    solutions manual?

  The book is intended as a reference, aiming to support researchers in using
  the Targeted Learning framework to develop rigorous answers to questions
  rooted in causal inference and arising in a vast array of application areas.
  Principally, the handbook aims to provide a gentle but thorough introduction
  to the Targeted Learning methodology, demonstrating its use in case studies
  rather than delving into theoretical detail at the level of existing reference
  texts (see below) and relevant primary literature. Secondarily, the handbook
  also serves as long-form documentation of the rich set of tools available in
  the `tlverse` software ecosystem. While not intended as a textbook, this work
  has been developed and refined through teaching a series of workshops; to aid
  in this, a limited set of exercises were developed (solutions to exercises
  will be made available online and possibly in an appendix).

6. What related books are available and what advantages will your book have?

    * 2011, Targeted Learning: Causal Inference for Observational and
        Experimental Data (Mark van der Laan and Sherri Rose)
    * 2018, Targeted Learning in Data Science: Causal Inference for Complex
        Longitudinal Studies (Mark van der Laan and Sherri Rose)

  Th work is focus __not__ on providing in-depth technical descriptions or in
  communicating the most recent developments in state-of-the-art statistical
  methodology, in stark contrast to the two existing books on Targeted Learning.
  Instead, the goal is to convey both the underlying unifying concepts and key
  details of the broad set of state-of-the-art Targeted Learning techniques in
  a manner that is both clear and complete, without burdening the reader with
  extraneous technical details. We aim for the book to serve as a comprehensive
  reference and handbook for researchers -- methodologists and domain
  specialists alike -- who wish to employ the central statistical tools of the
  Targeted Learning framework efficiently. The material covered is a mix of key
  ideas in causal inference, machine learning, and non/semi-parametric
  statistical theory, all presented in the context of developing solutions to
  complex and carefully considered real-world data analysis questions.
  Throughout the text, all aspects of the methodology are demonstrated using the
  accompanying free and open source `tlverse` software ecosystem
  (https://github.com/tlverse).


## PRODUCTION

7. Approximately how many printed pages will your book contain? Are color
    figures essential to your book? If so, about how many would have to be in
    color? Color printing is still very expensive and color figures will
    increase the price so black and white should be used unless color is
    essential.

    * Approximate page count: 200-250 pages.
    * Color figures are *not* essential to our book.

8. When would you hope to be able to submit the final draft of the book to us?
    Will you use LaTeX, bookdown, or Word? We will supply a style file for
    LaTeX authors and request an unformatted file from Word authors.

    * A final draft of the book will be available by 01 January 2021.
    * Our book has been prepared using the bookdown R package.


## REVIEWS

9. Please give the names and e-mail addresses of four people who would be
     qualified to give an opinion on your proposed book.

    * Ashley Naimi, University of Pittsburgh
    * Linda Valeri, Columbia University
    * Constantin Frangakis, Johns Hopkins University
    * Michael Hudgens, University of North Carolina
    * Peter Gilbert, Fred Hutchinson Cancer Research Center

<!--
  More ideas for the type of person we'd like to have as a reviewer:
  1. writer of TMLE software (e.g. Eric Polley, Susan Gruber)
  2. writer of well-documented / easy-to-use R software
  3. applied statistician or data scientist who has worked with TMLE/SL
-->

## MARKETING

10. If your book is aimed at a professional/research market, are there key
      societies outside of statistics to which it should be marketed? (All
      major statistical societies will be covered.)

   Society for Epidemiologic Research

11. Please list up to six key features of your proposed book that we can use
      in bulleted form.

    * Learn to leverage the powerful Targeted Learning framework and methodology
         without a mathematical deep-dive
    * Introduces both classical and cutting-edge topics in causal inference and
        non/semi-parametric statistics
    * Explores predictive and prescriptive data analytics with the `tlverse`
        software ecosystem
    * See Targeted Learning in action by engaging with real-world case studies
        and hands-on examples

12. Please list up to six key words or phrases that people interested in this
      topic may use to search Amazon or the web. Do not repeat words in the
      title as these will already be found.

    * Causal inference, machine learning, statistical computing,
        biostatistics, data science, predictive/prescriptive analytics

13. Please select the three most important markets for your book. Other
      categories are available including education, psychology, and economics
      so please mention other important disciplines.

   Selected market sections in __bold__ below.

> STA01A-Statistics-Statistics for Life Sciences

> STA02A-Statistics-Introductory Statistics & General References

> STA04A-Statistics-Probability Theory & Applications

> STA06A-Statistics-SPC/Reliability/Quality Control

> **STA07A-Statistics-Statistical Theory & Methods**

> **STA07J-Statistics-Statistical Learning & Data Mining**

> STA08A-Statistics-Computational Statistics

> STA09A-Statistics-Statistics for Business, Finance & Economics

> STA10A-Statistics-Statistics for Engineering and Physical Science

> **STA11A-Statistics-Biostatistics and Epidemiology**

> STA12A-Statistics-Statistics for the Social and Behavioral Sciences

> STA13A-Statistics-Environmental Statistics

> STA14A-Statistics-Statistical Genetics & Bioinformatics
