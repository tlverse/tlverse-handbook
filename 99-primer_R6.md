# `R6` Class System Primer {#r6}

A central goal of the Targeted Learning statistical paradigm is to estimate
scientifically relevant parameters in realistic (usually nonparametric) models.

The `tlverse` is designed using basic OOP principles and the `R6` OOP framework.
While we've tried to make it easy to use the `tlverse` packages without worrying
much about OOP, it is helpful to have some intuition about how the `tlverse` is
structured. Here, we briefly outline some key concepts from OOP. Readers
familiar with OOP basics are invited to skip this section.

## Classes, Fields, and Methods

The key concept of OOP is that of an object, a collection of data and functions
that corresponds to some conceptual unit. Objects have two main types of
elements:

1. _fields_, which can be thought of as nouns, are information about an object,
   and
2. _methods_, which can be thought of as verbs, are actions an object can
   perform.

Objects are members of classes, which define what those specific fields and
methods are. Classes can inherit elements from other classes (sometimes called
base classes) -- accordingly, classes that are similar, but not exactly the
same, can share some parts of their definitions.

Many different implementations of OOP exist, with variations in how these
concepts are implemented and used. R has several different implementations,
including `S3`, `S4`, reference classes, and `R6`. The `tlverse` uses the `R6`
implementation. In `R6`, methods and fields of a class object are accessed using
the `$` operator. For a more thorough introduction to `R`'s various OOP systems,
see http://adv-r.had.co.nz/OO-essentials.html, from Hadley Wickham's _Advanced
R_ [@wickham2014advanced].

## Object Oriented Programming: `Python` and `R`

OO concepts (classes with inherentence) were baked into Python from the first
published version (version 0.9 in 1991). In contrast, `R` gets its OO "approach"
from its predecessor, `S`, first released in 1976. For the first 15 years, `S`
had no support for classes, then, suddenly, `S` got two OO frameworks bolted on
in rapid succession: informal classes with `S3` in 1991, and formal classes with
`S4` in 1998. This process continues, with new OO frameworks being periodically
released, to try to improve the lackluster OO support in `R`, with reference
classes (`R5`, 2010) and `R6` (2014). Of these, `R6` behaves most like Python
classes (and also most like OOP focused languages like C++ and Java), including
having method definitions be part of class definitions, and allowing objects to
be modified by reference.
