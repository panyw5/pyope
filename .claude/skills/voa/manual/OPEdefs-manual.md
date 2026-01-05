# A Mathematica™ Package for Computing Operator Product Expansions (OPErefs 3.1)

K. Thielemans*

Theoretical Physics Group Imperial College London SW7 2BZ (UK)

April 1995

# Abstract

A general purpose MathematicaTM package for computing Operator Product Expansions of composite operators in meromorphic conformal field theory is described. Given the OPEs for a set of "basic" fields, OPEs of arbitrarily complicated composites can be computed automatically. Normal ordered products are always reduced to a standard form. As an explicit example, the conformal anomaly for superstrings is computed.

The most important extensions with respect to the first version of the package are the ability the check the Jacobi identities, and to compute Poisson brackets ("classical OPEs").

# 1 Introduction

Operator Product Expansions (OPEs) are extensively used in conformal field theories, for example in string theory and statistical physics. They are used to evaluate expectation values of several fields with arguments in neighbouring points. One then writes

$$
A (z) B (w) = \frac {C (w)}{(z - w) ^ {2}} + \frac {D (w)}{z - w} + \dots \tag {1}
$$

where the dots represent the regular terms in the Laurent expansion. The above expression is valid when evaluating expectation values for  $z \to w$ .

In handling OPEs, the problem of normal ordering requires special attention. In the so called point splitting regularization scheme, the normal ordered product of two operators is defined as the zero'th order term in their OPE. The OPE of an operator with a composite follows from a variation on Wick's theorem [1]. However, calculations quickly become long and error prone when composites of composites are involved.

The definition of normal ordering is noncommutative and nonassociative. As a consequence, it is recommended to use a standard order for the operators and introduce a standard way of normal ordering composites of several operators. To this end, there exists a set of prescriptions and again these formulas are conceptually quite simple, but have to be applied recursively in complicated cases.

OPErefs is written in MathematicaTM $^{1}$ , an interactive environment for performing symbolic computations. The advantages of writing the package in Mathematica are numerous. Its programming language has very powerful pattern matching capabilities. The result of a computation can be transformed with the help of built-in functions (expanding, factoring, collecting terms in specified variables ...) and one may even obtain output in TEX-form. Mathematica is running on a lot of machines, from PC's to supercomputers and all versions are completely compatible. And last but not least, it is considerably less time-consuming to program this problem in Mathematica than, for instance, in  $C$ .

In using the package, all "basic" operators of the theory have to be declared to be bosonic or fermionic (parafermions are not supported), and the OPEs of the basic operators must be given. Fields with any (also negative) conformal dimension may be used, the only restriction is that fractional

powers in the Laurent expansion are not supported. With this input, the package is able to compute OPEs of arbitrarily complicated composites, given enough memory and time. Also, normal ordered products are automatically reduced to standard form (putting operators in the order they are declared and normal ordering them from right to left).

Version 2.0 of OPErefs is already presented in [3]. Version 3.0 is described in [4]. In this paper, version 3.1 is briefly described. A more complete overview of the OPE formalism and of the internals of OPErefs is given in [5] to which I refer for further details.

The paper is organised as follows. In section 2, the necessary rules for computing OPEs are given and the algorithm is discussed briefly. The next section explains how to use the package. Finally, some runtimes are given.

# Notation

Input for and output from Mathematica is written in typeset font. Input lines are preceded by "In[n] :=", and corresponding output statements by "Out[n] =", as in Mathematica.

We use the following notation for OPEs

$$
A (z) B (w) = \sum_ {n <   = h (A, B)} \frac {[ A B ] _ {n} (w)}{(z - w) ^ {n}} , \tag {2}
$$

where  $h(A, B)$  is some finite number, and is usually given by the sum of the conformal dimensions of  $A$  and  $B$ .  $[AB]_0$  will be called the normal ordered product of  $A$  and  $B$ .

For Poisson (or in fact Dirac) brackets we use

$$
\{A (x), B \left(x _ {0}\right) \} _ {\mathrm {P B}} = \sum_ {n > 0} \frac {(- 1) ^ {n - 1}}{(n - 1) !} \left\{A B \right\} _ {n} \left(x _ {0}\right) \partial^ {n - 1} \delta^ {(2)} \left(x - x _ {0}\right), \tag {3}
$$

where  $\delta^{(2)}(x - x_0)$  is the Dirac delta-function on the plane. The derivative is with respect to the  $x$ -coordinate. We choose the normalisation factors such that:

$$
\{A B \} _ {n} (y) = \int d x ^ {2} (x - y) ^ {n - 1} \{A (x), B (y) \} _ {\mathrm {P B}}. \tag {4}
$$

# 2 Formulas

In [3], all formulas (extracted from [1]) needed to compute OPEs are given. The more complicated rules can be derived in the following way, see [5] for

more details.

Using the definition (2) for the OPEs, and Cauchy's residue formula for contour integrals, we can isolate the contribution of a certain part of the OPE by taking appropriate contour integrals:

$$
[ A [ B C ] _ {p} ] _ {q} (u) = \oint_ {C _ {u}} \frac {d z}{2 \pi i} (z - u) ^ {q - 1}
$$

$$
\oint_ {C _ {u}} \frac {d w}{2 \pi i} (w - u) ^ {p - 1} A (z) B (w) C (u), \tag {5}
$$

where  $C_u$  denotes a contour which encircles  $u$  once anti-clockwise. We can now use a contour deformation argument relating the contour integral in eq. (5) to a contour integral where the integration over  $w$  is performed last. This integral has two terms: one where the  $z$  contour is around  $u$ , and one where it is around  $w$ . We find:

$$
[ A [ B C ] _ {p} ] _ {q} = (- 1) ^ {| A | | B |} [ B [ A C ] _ {q} ] _ {p} + \sum_ {l > 0} \binom {q - 1} {l - 1} \left[ [ A B ] _ {l} C \right] _ {p + q - l}. \tag {6}
$$

This equation has to hold (inside correlators) for consistency of the OPE-formalism. It is valid for any integers  $p, q$ , i.e. also negative numbers. However, in practice we use the Jacobi identities eq. (6) for positive  $p, q$  as equations for the singular part of the OPEs, while for  $p$  or  $q$  zero, they define how one should calculate with composites.

OPErefs 3.1 can calculate OPEs and Poisson brackets. We concentrate on OPEs here and briefly comment which changes are needed for the other case. However, note that the Jacobi identities eq. (6) for positive  $p, q$  remain exactly the same.

# 2.0.1 Computing OPEs

The package should compute OPEs of arbitrarily complicated composites when a set of generators and their OPEs is given. For easy reference, we list the rules we need, aside from linearity. First, there are the rules involving derivatives:

$$
[ \partial A B ] _ {q} = - (q - 1) [ A B ] _ {q - 1} \tag {7}
$$

$$
[ A \partial B ] _ {q} = (q - 1) [ A B ] _ {q - 1} + \partial [ A B ] _ {q}. \tag {8}
$$

Next, we need the OPE of  $B$  with  $A$ , given the OPE of  $A$  with  $B$ :

$$
[ B A ] _ {q} = (- 1) ^ {| A | | B |} \sum_ {l \geq q} \frac {(- 1) ^ {l}}{(l - q) !} \partial^ {(l - q)} [ A B ] _ {l}. \tag {9}
$$

OPEs with composites can be calculated using eq. (6):

$$
\begin{array}{l} [ A [ B C ] _ {0} ] _ {q} = (- 1) ^ {| A | | B |} [ B [ A C ] _ {q} ] _ {0} + [ [ A B ] _ {q} C ] _ {0} \\ + \sum_ {l = 1} ^ {q - 1} \binom {q - 1} {l} \left[ \left[ A B \right] _ {q - l} C \right] _ {l} \tag {10} \\ \end{array}
$$

where  $q \geq 1$ . Furthermore, we will use:

$$
[ A A ] _ {0} = - \sum_ {l > 0} \frac {(- 1) ^ {l}}{2 l !} \partial^ {l} [ A A ] _ {l}, \quad \text {f o r} A \text {f e r m i o n i c}. \tag {11}
$$

These rules are the only ones needed to compute every OPE in the OPA. Indeed, when computing the OPE of  $A$  with  $B$ , we apply the following procedure:

- if  $A$  and  $B$  are generators whose OPE we know, return it as the result.  
- apply linearity if necessary.  
- if  $A$  is an operator with derivatives, use eq. (7).  
- if  $B$  is an operator with derivatives, use eq. (8).  
- if  $A$  or  $B$  contains a composite, apply eq. (11) if necessary.  
- if  $B$  is a composite, use eq. (10).  
- if  $A$  is a composite, use eq. (9).  
- if the OPE  $B(z) A(w)$  is known, compute the OPE  $A(z) B(w)$  using eq. (9).

This list should be used recursively until none of the rules applies, which means that the OPE has been calculated. The order in which we check the rules is in this case not important, but we will check them in a "top-down" order. Note that to compute an OPE of a composite with a generator, first eq. (9) is used, and eq. (10) in the next step.

The algorithm used to compute Poisson brackets is essentially the same. Some small changes in the rule eq. (10) are needed, namely all double contractions have to be dropped. Also, eq. (11) is changed to  $[AA]_0 = 0$  when  $A$  is fermionic.

# 2.0.2 Simplifying composites

We now discuss how we can reduce normal ordered products to a standard form. We define an order on the generators and their derivatives, e.g. lexicographic ordering. Given a composite, we apply the rules given below until all composites are normal ordered from right to left and the operators are ordered, i.e.  $[A[B[C\ldots]_0]_0]_0$ , and  $A \leq B \leq C$ . The relevant formulas are:

$$
\partial [ A B ] _ {0} = [ A \partial B ] _ {0} + [ \partial A B ] _ {0} \tag {12}
$$

$$
[ A B ] _ {- q} = \frac {1}{q !} \left[ \left(\partial^ {q} A\right) B \right] _ {0}, \quad q \geq 1 \tag {13}
$$

$$
[ B A ] _ {0} = (- 1) ^ {| A | | B |} [ A B ] _ {0} + (- 1) ^ {| A | | B |} \sum_ {l \geq 1} \frac {(- 1) ^ {l}}{l !} \partial^ {l} [ A B ] _ {l} \tag {14}
$$

$$
[ A [ B C ] _ {0} ] _ {0} = (- 1) ^ {| A | | B |} [ B [ A C ] _ {0} ] _ {0} + \tag {15}
$$

$$
[ \left([ A B ] _ {0} - (- 1) ^ {| A | | B |} [ B A ] _ {0}\right) C ] _ {0} \tag {16}
$$

For the case of Poisson brackets, the rules (14) and (16) drastically simplify, only the first terms remain.

# 2.0.3 Improvements

The rules given in this section up to now are sufficient to compute any OPE, and to reorder any composite into a standard form. However, some shortcuts exist  $^2$ , and where introduced in OPErefs 3.x.

$[A[BC]_0]_q$  can be computed using eq. (10), but an alternative is:

$$
\begin{array}{l} [ A [ B C ] _ {0} ] _ {q} = (- 1) ^ {| A | | B |} \left([ B [ A C ] _ {q} ] _ {0} + \sum_ {l \geq 0} \frac {(- 1) ^ {l + q}}{l !} [ \partial^ {l} [ B A ] _ {l + q} C ] _ {0} \right. \\ \left. + \sum_ {l = 1} ^ {q - 1} (- 1) ^ {l} \left[ \left[ B A \right] _ {l} C \right] _ {q - l}\right). \tag {17} \\ \end{array}
$$

This rule is more convenient when we know the OPE  $B(z)A(w)$  while  $A(z)B(w)$  has to be computed using eq. (9).

Similarly, to compute an OPE where the first operator is a composite, the algorithm as presented above uses eqs. (9) and (10). However, the next rule implements this in one step:

$$
\begin{array}{l} [ [ A B ] _ {0} C ] _ {q} = \sum_ {l \ge 0} \frac {1}{l !} [ \partial^ {l} A [ B C ] _ {l + q} ] _ {0} + (- 1) ^ {| A | | B |} \sum_ {l \ge 0} \frac {1}{l !} [ \partial^ {l} B [ A C ] _ {l + q} ] _ {0} \\ + (- 1) ^ {| A | | B |} \sum_ {l = 1} ^ {q - 1} [ B [ A C ] _ {q - l} ], \tag {18} \\ \end{array}
$$

where  $q\geq 1$  and:

$$
\begin{array}{l} [ [ A B ] _ {0} C ] _ {0} = [ A [ B C ] _ {0} ] _ {0} + \\ \sum_ {l > 0} \frac {1}{l !} [ \partial^ {l} A [ B C ] _ {l} ] _ {0} + (- 1) ^ {| A | | B |} \sum_ {l > 0} \frac {1}{l !} [ \partial^ {l} B [ A C ] _ {l} ] _ {0}. \tag {19} \\ \end{array}
$$

In OPErefs 3.1 these two last rules are applied when  $C$  is not a composite itself.

# 3 User's Guide

This section is intended as a user's guide to the package OPErefs 3.1. Explicit examples are given for most operations. Note that OPErefs 3.1 requires Mathematica 1.2 or later.

We introduce some special notations. Input for and output from Mathematica is written in typeset font. Input lines are preceded by "In[n] :=", and corresponding output statements by "Out[n] =", as in Mathematica.

As OPErefs is implemented as a Mathematica package, it has to be loaded before any of its global symbols is used. Loading the package a second time will clear all previous definitions of operators and OPEs, as well as all stored intermediate results. Assuming that the package is located in the Mathematica-path, e.g. in your current directory, issue:

$$
I n [ 1 ] := \quad <   <   O P E d e f s. m
$$

After loading OPErefs into Mathematica, help for all the global symbols is provided using the standard help-mechanism, e.g. ?OPE.

Now, you need to declare the operators that will be used. If you want to define bosonic operators  $\mathbf{T}$  and  $\mathbf{J}[\mathbf{i}]$  (any index could be used), and fermionic operators  $\mathbf{psi}[\mathbf{i}]$ , the corresponding statements are:

$$
\left. \operatorname {I n} / 2 \right] := \quad \text {B o s o n i c} [ T, J [ i _ {-} ] ]
$$

$$
I n [ 3 ] := \quad \text {F e r m i o n i c} [ \mathrm {p s i} [ \mathrm {i} _ {-} ] ]
$$

The order of the declarations fixes also the ordering of operators used by the program:

$$
T <   J [ 1 ] ^ {\prime} <   J [ 1 ] <   J [ 2 ] <   J [ i ] <   p s i [ 1 ] <   \dots \tag {20}
$$

By default, derivatives of an operator are considered "smaller" than the operator itself. This can be reversed using the global options N00ordering (see below).

Finally, the non-regular OPEs between the basic operators have to be given. An OPE can be specified in two different ways.

The first way is by listing the operators that occur at the poles, the first operator in the list is the one at the highest non-zero pole, the last operator has to be the one at the first order pole, e.g.:

$$
I n [ 4 ] := \quad O P E [ T, T ] = M a k e O P E [ \{c / 2 O n e, 0, 2 T, T ^ {\prime} \} ];
$$

Note the operator One which specifies the unit-operator.

The second way is by giving the OPE as a Laurent series expansion, adding the symbol  $\mathbf{0}$ rd which specifies the (implicit) arguments of the operators for which the OPE is defined<sup>3</sup>. The arguments for the operators can be any Mathematica expression.

Warning: it is important that the operators occurring as arguments of OPE in a definition should be given in standard order (20), otherwise wrong results will be generated.

The following statements define a  $SU(2)_k$ -Kac-Moody algebra:

$$
\begin{array}{l} I n [ 5 ] := \quad O P E [ J [ i _ {-} ], J [ i _ {-} ] ] := \\ \text {M a k e O P E} [ - k / 2 (z - w) ^ {\wedge} - 2 + \operatorname {O r d} [ z, w, 0 ] ] \\ \end{array}
$$

$$
\begin{array}{l} I n [ 6 ] := \quad O P E [ J [ 1 ], J [ 2 ] ] = \\ \text {M a k e O P E} [ J [ 3 ] [ w ] (z - w) ^ {\sim} - 1 + 0 r d [ z, w, 0 ] ]; \\ \end{array}
$$

$$
\begin{array}{l} I n [ 7 ] := \quad O P E [ J [ 2 ], J [ 3 ] ] = \\ \text {M a k e O P E} [ J [ 1 ] [ w ] (z - w) ^ {\wedge} - 1 + \operatorname {O r d} [ z, w, 0 ] ]; \\ \end{array}
$$

$$
\begin{array}{l} I n [ 8 ] := \quad O P E [ J [ 1 ], J [ 3 ] ] = \\ \text {M a k e O P E} [ - J [ 2 ] [ w ] (z - w) ^ {\wedge} - 1 + \operatorname {O r d} [ z, w, 0 ] ]; \\ \end{array}
$$

In fact, with the above definitions, one has to use always the explicit indices 1,2,3 for the currents  $J$ . If we would compute an OPE with current J[i]

where the index  $i$  is not 1,2 or 3, wrong results will be given. One can circumvent this peculiarity by reformulating the definitions.

A normal ordered product  $[AB]_0$  is entered in the form  $\mathsf{NO}[\mathsf{A},\mathsf{B}]$ . Multiple composites can be entered using only one NO head, e.g.  $\mathsf{NO}[\mathsf{A},\mathsf{B},\mathsf{C}]$ . This input is effectively translated into  $\mathsf{NO}[\mathsf{A},\mathsf{NO}[\mathsf{B},\mathsf{C}]]$ . All output is normal ordered with the same convention, i.e. from right to left (input can be in any order). Also, the operators in composites will always be ordered according to the standard order (20).

As an example, we can define the Sugawara energy-momentum tensor for  $\widehat{SU(2)}_k$ . The Mathematica output of an OPE is a list of the operators at the poles.

$$
\begin{array}{c} I n [ 9 ] := \hskip 1 4. 2 2 6 3 7 8 p t T s = - 1 / (k + 2) (N O [ J [ 1 ], J [ 1 ] ] + \\ \hskip 1 4. 2 2 6 3 7 8 p t N O [ J [ 2 ], J [ 2 ] ] + N O [ J [ 3 ], J [ 3 ] ]) ; \end{array}
$$

$$
I n [ 1 0 ] := \text {O P E S i m p l i f y} [ O P E [ T s, J [ 1 ] ] ]
$$

$$
O u t [ 1 0 ] = \ll 2 | | J [ 1 ] | | 1 | | J [ 1 ] ^ {\prime} > >
$$

Warning: when computing OPEs with composites, or when reordering composites, OPEdefines remembers by default some intermediate results. Thus, it is dangerous to change the definition of the basic OPEs after some calculations have been performed. For example, consider a constant  $a$  in an OPE. If calculations are performed after assigning a value to  $a$ , the intermediate results are stored with this value. Changing  $a$  afterwards will give wrong results.

The other globally defined functions available from the package are:

- OPEOperator[operator_, parity_] provides a more general way to declare an operator than Bosonic and Fermionic. The second argument is the parity of the operator such that  $(-1)^{\text{parity}}$  is  $+1$  for a boson, and  $-1$  for a fermion. It can be a symbolic constant. This is mainly useful for declaring a bc-system of unspecified parity, or a Kac-Moody algebra based on a super-Lie algebra. In such cases, the operator can contain a named pattern:

$$
I n [ 1 1 ] := \quad O P E O p e r a t o r [ J [ i _ {-} ], p a r i t y [ i ] ]
$$

If one wants to declare more operators, one can group each operator and its parity in a list:

$$
I n [ 1 2 ] := O P E O p e r a t o r \left[ \left\{b [ i _ {-} ], p a r i t y [ i ] \right\}, \left\{c [ i _ {-} ], p a r i t y [ i ] \right\} \right]
$$

See also SetOPEOptions[ParityMethod, _].

- OPEPole[n_][ope_] gets a single pole term of an OPE:

$$
I n [ 1 3 ] := \quad O P E P o l e [ 2 ] [ O u t [ 1 0 ] ]
$$

Out[13] = J[1]

OPEPole[n_][A_,B_] can also be used to compute only one pole term of an OPE:

In[14] := Factor[OPEPole[4][Ts, Ts]]

Out[14] = (3 k 0ne)/(2 (2 + k))

OPEPole can also give terms in the regular part of the OPE:

$\operatorname{In}[15] := \quad \text{OPEPole}[-1][T, T]$

Out[15] = NO[T', T]

- MaxPole[ope_] gives the order of the highest pole in the OPE.  
- OPEParity[A] returns an even (odd) integer of  $A$  is bosonic (fermionic).  
- OPESimplify[ope_, function_] "collects" all terms in ope with the same operator and applies function on the coefficients. When no second argument is given, the coefficients are Expanded.

In[16] := OPESimplify[OPE[J[1], NO[J[2], J[1]]]]

Out[16] = << 2 || (1 - k/2) J[2] || 1 ||

NO[J[1]，J[3]]+J[2]’>>

OPESimplify[pole_, function_] does the same simplifications on sums of operators.

- OPEMap[function_, ope_] maps function to all poles of ope.  
- GetCoefficients[expr_] returns a list of all coefficients of operators in expr which can be (a list of) OPEs or poles.  
- OPEJacobi[op1_, op2_, op3_] computes the Jacobi-identities (6) for the singular part of the OPEs of the three arguments. Due to the nature of eq. (6), the computing time will be smallest (in most cases) when op1 ≤ op2 ≤ op3 in the order (20). Note that it is sufficient to check the Jacobi-identities with the generators of the OPA.

The result of OPEJacobi is a double list of operators. It is generated by

Table[OPEPole[n][A,OPEPole[m][B,C]] +

corrections,  $\{\mathfrak{m},\max \} ,\{\mathfrak{n},\max \}\}$

All elements of the list should be zero up to null operators for the OPA to be associative.

- Delta[i_,j_] is the Kronecker delta symbol  $\delta_{ij}$ .

- ClearOPESavedValues[] clears all stored intermediate results, but not the definition of the operators and their OPEs. To clear everything, reload the package.  
- OPEToSeries[ope_] converts an OPE to a Laurent series expansion in  $\mathbf{z}$  and  $\mathbf{w}$ . The arguments can be set to  $\mathbf{x}$  and  $\mathbf{y}$  with:

$$
I n [ 1 7 ] := \text {S e t O P E O p t i o n s} [ \text {S e r i e s A r g u m e n t s}, \{x, y \} ]
$$

- TeXForm[ope_] gives TExoutput for an OPE. The same arguments are used as in OPETOseries.  
- OPESave[filename_] (with filename a string between double quotes) saves the intermediate results that OPErefs remembers to file (see the option OPESaving below).  
- SetOPEOptions is a function to set the global options of the package. The current options are:

- SetOPEOptions[SeriesArguments, {arg1_, arg2_}]: sets arguments to be used by TeXForm and OPETOSeries. One can use any Mathematica expression for arg1 and arg2.  
- SetOPEOptions[NOOrdering, n_]: if n is negative, order higher derivatives to the left (default), if n is positive, order them to the right.  
- SetOPEOptions[ParityMethod, 0|1]: makes it possible to use operators of an unspecified parity. When the second argument is 0 (default), all operators have to be declared to be bosonic or fermionic. When the argument is 1, OPEOperator can be used with a symbolic parity. Note that in this case, powers of  $-1$  are used to compute signs, which is slightly slower than the boolean function which is used by the first method.

This option is not normally needed as the use of OPEOperator with a non-integer second argument sets this option automatically.

- SetOPEOptions[OPESaving, boolean_]: if boolean evaluates to True (default), OPErefs stores the intermediate results when computing OPEs of composites and when reordering composites. This option is useful if Mathematica runs short of memory in a

large calculation, or when computing with dummy indices<sup>4</sup>.

- SetOPEOptions[OPEMethod, method_]: with the parameter method equal to QuantumOPEs enables normal OPE computations (default setting), while ClassicalOPEs enables Poisson bracket computations. Using this option implicitly calls ClearOPESavedValues[].

# 4 Example : The conformal anomaly in superstring theory

We consider only one free boson field  $X$  and one free fermion field  $\psi$  because additional free fields will have exactly the same OPEs and commute with each other. We denote  $\partial X$  with  $J$  and  $\psi$  with  $p s i$  (we normalise them such that they have a  $+1$  in their OPEs. The ghosts are a fermionic  $b, c$  system (operators  $b, c$ ) and a bosonic  $\beta, \gamma$  system (operators  $B, G$ ).  $b$  has conformal dimension 2 and  $\beta$  has  $3/2$ . It is now a trivial task to compute the conformal anomaly:

```latex
$\begin{array}{rl} & {\mathrm{In}[1]\coloneqq \mathrm{<  }\mathrm{<  }0\mathrm{E}\mathrm{d}\mathrm{e}\mathrm{s}. \mathrm{m}}\\ & {\mathrm{In}[2]\coloneqq \mathrm{Bosonic}[\mathrm{J},\mathrm{B},\mathrm{G}];\mathrm{Fermionic}[\mathrm{b},\mathrm{c},\mathrm{psi}];}\\ & {\mathrm{OPE}[\mathrm{J},\mathrm{J}] = \mathrm{MakeOPE}[\{\mathrm{One},\mathrm{0}\} ];}\\ & {\mathrm{OPE}[\mathrm{psi},\mathrm{psi}] = \mathrm{MakeOPE}[\{\mathrm{One}\} ];}\\ & {\mathrm{OPE}[\mathrm{b},\mathrm{c}] = \mathrm{MakeOPE}[\{\mathrm{One}\} ];}\\ & {\mathrm{OPE}[\mathrm{B},\mathrm{G}] = \mathrm{MakeOPE}[\{\mathrm{One}\} ];}\\ & {\mathrm{Tb} = 1 / 2\mathrm{NO}[\mathrm{J},\mathrm{J}];\mathrm{Tf} = -1 / 2\mathrm{NO}[\mathrm{psi},\mathrm{psi}^{\prime}];}\\ & {\mathrm{Tbc} = -2\mathrm{NO}[\mathrm{b},\mathrm{c}^{\prime}] - \mathrm{NO}[\mathrm{b}^{\prime},\mathrm{c}];}\\ & {\mathrm{TBG} = 3 / 2\mathrm{NO}[\mathrm{B},\mathrm{G}^{\prime}] + 1 / 2\mathrm{NO}[\mathrm{B}^{\prime},\mathrm{G}];}\\ & {\mathrm{In}[3]\coloneqq \mathrm{OPESimplify}[\mathrm{OPE}[\mathrm{Tb},\mathrm{Tb}]]}\\ & {\mathrm{Out}[3] = \ll 4||\mathrm{One}/2||3||0||2||\mathrm{NO}[\mathrm{J},\mathrm{J}]||1||}\\ & {\qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \text{} NO[\mathrm{J}^{\prime},\mathrm{J}] >>}\\ & {\mathrm{In}[4]\coloneqq \mathrm{OPESimplify}[\mathrm{OPE}[\mathrm{Tf},\mathrm{Tf}]]}\\ & {\mathrm{Out}[4] = \ll 4||\mathrm{One}/4||3||0||2||\mathrm{NO}[\mathrm{psi}^{\prime},\mathrm{psi}]||1||}\\ & {\qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad NO[\mathrm{psi}^{\prime}],\mathrm{psi}] / 2>>}\\ & {\mathrm{In}[5]\coloneqq \mathrm{OPESimplify}[\mathrm{OPE}[\mathrm{Tbc},\mathrm{Tbc}]]}\\ & {\mathrm{Out}[5] = \ll 4||-13\mathrm{One}||3||0||2||}\\ & {\qquad \qquad \qquad \qquad -4\mathrm{NO}[\mathrm{b},\mathrm{c}^{\prime}] - 2\mathrm{NO}[\mathrm{b}^{\prime},\mathrm{c}]||1||}\\ & {\qquad \qquad \qquad -2\mathrm{NO}[\mathrm{b},\mathrm{c}^{\prime}) - 3\mathrm{NO}[\mathrm{b}^{\prime},\mathrm{c}^{\prime}] - \mathrm{NO}[\mathrm{b}^{\prime}),\mathrm{c}] >>} \end{array}$
```

$$
I n [ 6 ] := \quad O P E S i m p l i f y [ O P E [ T B G, T B G ] - M a k e O P E [ \{2 T B G, T B G ^ {\prime} \} ]
$$

$$
O u t [ 6 ] = \ll 4 | | 1 1 0 n e / 2 | | 3 | | 0 | | 2 | | 0 | | 1 | | 0 > >
$$

We see that each bosonic (fermionic) field will contribute a central charge 1  $(1/2)$  to the total central charge of the theory. The  $b$ ,  $c$  system contributes  $-26$ , and the  $\beta$ ,  $\gamma$  system 11. This gives the well-known relation for the critical dimensions of the bosonic string  $D_{b} - 26 = 0$  and the superstring  $3/2D_{s} - 26 + 11 = 0$ . Moreover, we can easily verify that the energy-momentum tensors obey the Virasoro algebra.

The reader without experience in CFT is invited at this point to take out some time and compute the OPE for  $T_{BG}$ , for instance, by hand. Although this computation is rather trivial with OPErefs, the same calculation was attempted in [6] using the mode-algebra. There it proved not to be possible to compute the Virasoro algebra automatically due to difficulties with the infinite sums in the normal ordered products.

# 5 Performance

In [3], a free field realization for  $\widehat{B_2}$  level  $k$  using (bosonic)  $\beta$ ,  $\gamma$  systems was constructed. In Table 1, we tabulate CPU times for computing an OPE of two of the currents, and the Sugawara tensor for this realization. The first time given in the table is the time for evaluating the statement after loading the package and defining the realization. The time between brackets is measured when the statement is repeated. Note that version 2.0 of Mathematica is roughly 1.4 times slower than version 1.2!

Table 1: CPU time for the computation of the OPE of the currents corresponding to the positive simple root of  $\widehat{B_2}$  (statement 9) and the computation of the Sugawara tensor (statement 11) (see Ref. [3]) for Mathematica running on a PC 386 (25 Mhz).

<table><tr><td>Mathematica-version
OPErefs-version</td><td>1.2
2.0</td><td>1.2
3.1</td><td>2.0
3.1</td></tr><tr><td>In[9]
In[11]</td><td>23.5 (4.5) s
43.2 (11.6) s</td><td>14.9 (2.8) s
31.3 (9.4) s</td><td>19.3 (3.8) s
40.7 (12.1) s</td></tr></table>

# 6 How to get it, and the future

If you are interested in OPErefs, you can get it by Email from the author. Please put a reference to [3] in your paper when you use it. Questions, remarks and improvements are welcome. Already in testing-phase is OPEconf, a package which enables you to work with (quasi-)primaries and conformal blocks.

# 7 Acknowledgements

Many improvements implemented in OPErefs 3.0 are suggested by K. Hornfeck, who also helped me testing this version.

# References

[1] F. Bais, P. Bouwknegt, M. Surridge, K. Schoutens, Nucl. Phys. B304 (1988) 348.  
A. Sevrin, W. Troost, A. Van Proeyen, P. Spindel, Nucl. Phys. B311 465 (1988).  
[2] Mathematica, A system for Doing Mathematics by Computer, S. Wolfram, Addison-Wesley Publishing Company, Inc.  
[3] K. Thielemans, Int. J. Mod. Phys. C Vol. 2, No. 3, 787 (1991).  
[4] K. Thielemans, New computing techniques in Physics Research II, proceedings of the Second International Workshop on Software Engineering, Artificial Intelligence and Expert Systems in High Energy and Nuclear Physics, ed. D. Perret-Gallix, World Scientific (1992).  
[5] K. Thielemans, An Algorithmic Approach to Operator Product Expansions, W-algebras and W-strings, PhD thesis KU Leuven, june 1994.  
[6] W.M. Seiler, SUPERCALC, a REDUCE Package for commutator calculations, Karlsruhe preprint KA-THEP-20/90.