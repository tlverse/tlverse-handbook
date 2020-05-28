> I usually begin by reviewing a proposal. Most authors find these comments
> useful and revise their proposal, and we usually decide whether to offer
> a contract at this stage. If we agree on a contract, I review chapters as
> written in sets of about 100 pages. You want to learn early if you have the
> essential topics, are writing at the intended level, and have enough
> motivation and examples so you can adjust your style if necessary.

# CHAPMAN & HALL/CRC STATISTICS BOOK PROPOSAL

## TITLE AND AUTHOR/EDITOR(S)

1. Provisional title of your book.

   * The Hitchhiker's Guide to the tlverse: A Targeted Learning Practitioner's
      Handbook

2. Names, titles, affiliations, and email addresses for all authors/editors

      * Mark van der Laan, UC Berkeley, laan@berkeley.edu
      * Alan Hubbard, UC Berkeley, hubbard@berkeley.edu
      * Jeremy Coyle, jeremy.coyle@gmail.com
      * Nima Hejazi, UC Berkeley, nhejazi@berkeley.edu
      * Ivana Malenica, UC Berkeley, imalenica@berkeley.edu
      * Rachael Phillips, UC Berkeley, rachaelvphillips@berkeley.edu


## CONTENTS

3. Please include a complete table of contents including chapter and section
    headings.

    - [Why we need a statistical
      revolution](https://senseaboutscienceusa.org/super-learning-and-the-revolution-in-knowledge/)
    - The Roadmap and introductory case study: the WASH Benefits data
    - Introduction to the [`tlverse` software ecosystem](https://tlverse.org)
    - Cross-validation with the [`origami`](https://github.com/origami) package
    - Ensemble machine learning with the [`sl3`](https://github.com/tlverse/sl3)
       package
    - Targeted learning for causal inference with the
      [`tmle3`](https://github.com/tlverse/tmle3) package
    - Optimal treatments regimes and the
      [`tmle3mopttx`](https://github.com/tlverse/tmle3mopttx) package
    - Stochastic treatment regimes and the
      [`tmle3shift`](https://github.com/tlverse/tmle3shift) package


## SUBJECT/AUDIENCE

4. Please describe in detail the subject of your book. Why will this book be
    important, who will find it useful, and what is new? What background will
    you assume and can you name any books that will suffice?

  *The Hitchhiker's Guide to the `tlverse`: A Targeted Learning Practitioner's
  Handbook* is an open-source and fully-reproducible handbook for applying the
  Targeted Learning framework and statistical methodology in practice using the
  `tlverse` software ecosystem (https://github.com/tlverse).

  The contents of this handbook are meant to serve as a reference guide for
  applied research as well as materials that can be taught in a series of short
  courses focused on the applications of Targeted Learning. Each section
  introduces a set of distinct causal questions, motivated by a case study,
  alongside statistical methodology and software for assessing the causal claim
  of interest.

  [TO FILL IN]

5. Will your book be primarily a reference/monograph or textbook? If a
    textbook, for which courses will it be the primary text and at what level
    is the course taught? Will you include exercises sets and supply a
    solutions manual?

  The book is intended as a reference, aiming to help applied researchers get up
  to speed in using the Targeted Learning framework, and the tlverse software
  ecosystem, to bring statistical rigor to thee analysis of their data. As the
  book has been used as a reference for teaching workshops, it includes
  a limited set of exercises as well as solutions to those exercises. Since the
  exercises are entirely programming/coding-based, a solutions manual will be
  made freely available online.

6. What related books are available and what advantages will your book have?

    * 2011, Targeted Learning: Causal Inference for Observational and
        Experimental Data (Mark van der Laan and Sherri Rose)
    * 2018, Targeted Learning in Data Science: Causal Inference for Complex
        Longitudinal Studies (Mark van der Laan and Sherri Rose)

  The focus of this work is __not__ on providing in-depth technical descriptions
  or in communicating recent developments in state-of-the-art statistical
  methodology, in stark contrast to the two prior books on Targeted Learning.
  Instead, the goal is to convey key details of state-of-the-art techniques in
  an manner that is both clear and complete, without burdening the reader with
  extraneous technical information/details. We aim for the book to serve as
  a comprehensive reference and handbook for researchers -- methodologists and
  domain specialists alike -- who wish to employ the central tools the central
  statistical tools of the Targeted Learning framework efficiently. The material
  covered is a mix of key ideas in causal inference, machine learning, and
  non/semi-parametric statistical theory, all presented in the context of
  developing solutions to complex and carefully considered real-world data
  analysis questions using the accompanying free and open source `tlverse`
  software ecosystem (https://github.com/tlverse).


## PRODUCTION

7. Approximately how many printed pages will your book contain? Are color
    figures essential to your book? If so, about how many would have to be in
    color? Color printing is still very expensive and color figures will
    increase the price so black and white should be used unless color is
    essential.

    * Approximate page count: [TO FILL IN]
    * Color figures are *not* essential to our book.

8. When would you hope to be able to submit the final draft of the book to us?
    Will you use Latex, bookdown, or Word? We will supply a style file for
    LaTeX authors and request an unformatted file from Word authors.

    * A final draft of the book should be available by 01 September 2020.
    * Our book has been prepared using bookdown.


## REVIEWS

9. Please give the names and e-mail addresses of four people who would be
     qualified to give an opinion on your proposed book.

  [TO FILL IN]


## MARKETING

10. If your book is aimed at a professional/research market, are there key
      societies outside of statistics to which it should be marketed? (All
      major statistical societies will be covered.)

  [TO FILL IN]

11. Please list up to six key features of your proposed book that we can use
      in bulleted form.

  [TO FILL IN]

12. Please list up to six key words or phrases that people interested in this
      topic may use to search Amazon or the web. Do not repeat words in the
      title as these will already be found.

    * Causal inference, machine learning, applied statistics, statistical
        computing, statistical software (reproducible research?), biostatistics

13. Please select the three most important markets for your book. Other
      categories are available including education, psychology, and economics
      so please mention other important disciplines.

*Statistics*
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

*Biomedical Science*
> BIO01D-Biomedical Science-Bioinformatics
> BIO01I-Biomedical Science-Epidemiology
> BIO01O-Biomedical Science-Neuroscience
> BIO01V-Biomedical Science-Medical Devices

*Business*
> BUS08A-Business & Management-Business Information Systems

*Computer Science*
> CMS05-Computer Science & Engineering-Computational Biology
> CMS07-Computer Science & Engineering-Computer Graphics
> CMS08-Computer Science & Engineering-Visualization
> CMS11-Computer Science & Engineering-Data Mining and Knowledge Discovery
> CMS13-Computer Science & Engineering-Discrete Systems
> CMS15-Computer Science & Engineering-Machine Learning and Pattern Recognition
> CMS18-Computer Science & Engineering-High Performance Computing

*Environmental Science*
> ENV11A-Environmental Science-Environmental & Ecological Risk Assessment
> ENV12A-Environmental Science-Environmental  Modeling & Systems Analysis
> ENV13A-Environmental Science-Ecology

*Geoscience*
> GEO03-Geoscience-Geology

*Life Science*
> LIF05C-Life Science-Genetics

*Mathematics*
> MTH14A-Mathematics-Operations Research
> MTH15A-Mathematics-Probability Theory & Applications

*Medicine*
> MED02I-Medicine-Clinical Medicine-Microbiology & Infectious Diseases
> MED02J5-Medicine-Clinical Medicine-Medical Genetics
> PHR-Pharmaceutical Science & Regulation
> PHR03G-Pharmaceutical Science & Regulation-Clinical Trials
> PHR03L-Pharmaceutical Science & Regulation-Drug Development
