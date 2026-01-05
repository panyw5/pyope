# An Algorithmic Approach to Operator Product Expansions, W-Algebras and W-Strings

Promotor: Dr. W. Troost

Proefschrift ingediend tot het behalen van de graad van Doctor in de Wetenschappen door Kris Thielemans

Leuven, Juni 1994.

# Notations

# Abstract

String theory is currently the most promising theory to explain the spectrum of the elementary particles and their interactions. One of its most important features is its large symmetry group, which contains the conformal transformations in two dimensions as a subgroup. At quantum level, the symmetry group of a theory gives rise to differential equations between correlation functions of observables. We show that these Ward-identities are equivalent to Operator Product Expansions (OPEs), which encode the short-distance singularities of correlation functions with symmetry generators. The OPEs allow us to determine algebraically many properties of the theory under study. We analyse the calculational rules for OPEs, give an algorithm to compute OPEs, and discuss an implementation in Mathematica.

There exist different string theories, based on extensions of the conformal algebra to so-called  $W$ -algebras. These algebras are generically nonlinear. We study their OPEs, with as main results an efficient algorithm to compute the  $\beta$ -coefficients in the OPEs, the first explicit construction of the  $WB_{2}$ -algebra, and criteria for the factorisation of free fields in a  $W$ -algebra.

An important technique to construct realisations of  $W$ -algebras is Drinfeld's-Sokolov reduction. The method consists of imposing certain constraints on the elements of an affine Lie algebra. We quantise this reduction via gauged WZNW-models. This enables us in a theory with a gauged  $W$ -symmetry, to compute exactly the correlation functions of the effective theory.

Finally, we investigate the (critical)  $W$ -string theories based on an extension of the conformal algebra with one extra symmetry generator of dimension  $N$ . We clarify how the spectrum of this theory forms a minimal model of the  $W_N$ -algebra.

$$
\partial , \bar {\partial}: \frac {d}{d z}, \frac {d}{d \bar {z}}.
$$

$\delta_{\varepsilon}$  : infinitesimal transformation with  $\varepsilon (x)$  an infinitesimal parameter.

:  $T_{1}T_{2}$  : : normal ordening.

OPA : Operator Product Algebra, subsection 2.3.3.

OPE : Operator Product Expansions, section 2.3.

$[T_{1}(z)T_{2}(z_{0})]:$  the complete OPE.

$T_{1}(z)T_{2}(z_{0})$  : singular part of an OPE.

$[T_1T_2]_n$  : coefficient of  $(z - z_0)^{-n}$  in the OPE  $[T_{1}(z)T_{2}(z_{0})]$ , eq. (2.3.3).

$\ll \frac{c}{2} |0|2T|\partial T\gg$  : list of the operators at the singular poles in an OPE. Highest order pole is given first. First order pole occurs last in the list.

$\{T_{1}T_{2}\}_{n}(x_{0})$  : coefficient of  $\frac{(-1)^{n-1}}{(n-1)!}$ $\delta^{(2)}(x-x_0)$  in the Poisson bracket  $\{T_{1}(x)T_{2}(x_{0})\}$  subsection 2.3.5.

$\widehat{A}_m$  : mode of the operator  $A$ , section 2.4.

$L_{\{n\}}$  : is  $L_{n_p}\dots L_{n_1}$ $L_{\{-n\}}$  reversed order, chapter 4.

$\bar{g}$  : Lie algebra, see appendix B for additional conventions.

$\hat{g}:\mathrm{Ka}\check{\mathrm{c}}$  -Moody algebra.

$e_0, e_+$ ,  $e_-$ : generators of  $sl(2)$ .

$\mathcal{K}_{\pm}$  : kernel of the adjoint of  $e_{\pm}$ .

$\Pi$  : projection operator in  $\bar{g}$  , appendix B.

$\Pi_{hw}$  : projection on  $\kappa_{+}$ .

$\underline{c}$ ,  $\overline{c}$ : index limitation to generators of strictly negative, resp. positive grading, section 6.4.

$(a)_n$  : Pochhammer symbol,  $(a)_n = \frac{\Gamma(a + n)}{\Gamma(a)} = \prod_{j=0}^{n-1}(a + j)$ ;  $a \in \mathbb{R}$  and  $n \in \mathbb{N}$ , appendix 2.A.

# Contents

# 1 Introduction and outline 2

# 2 Conformal Field Theory and OPEs 7

2.1 Conformal transformations 7

2.2 Correlation functions and symmetry 8

2.2.1 Global conformal invariance 8

2.2.2 Ward identities 9

2.3Operator Product Expansions 11

2.3.1 OPEs and Ward identities 11

2.3.2 Consistency conditions for OPEs 12

2.3.3Operator Product Algebras 15

2.3.4 OPA-terminology 16

2.3.5 Poisson brackets 17

2.4 Mode algebra 18

2.5 Generating functionals 19

2.6 A few examples 21

2.6.1 The free massless scalar 21

2.6.2 The free Majorana fermion 23

2.6.3 Other first order systems 24

2.6.4 Wess-Zumino-Novikov-Witten models 24

2.A Appendix: Combinatorics 26

# 3 OPEs in Mathematica 28

3.1 Intention and history 28  
3.2 Design considerations 28  
3.3 Implementation 29  
3.3.1 Algorithm 29  
3.3.2 The internals of OPErefs 32  
3.3.3Operator handling 33  
3.3.4 Performance 37  
3.4 User's Guide 37  
3.5 Example : The conformal anomaly in superstring theory 39  
3.6 Future developments 40

# 3.7 Other packages 40

# 4  $\mathcal{W}$  -algebras 43

4.1 Introduction 43  
4.2 Highest weight representations and minimal models 44  
4.3 Consequences of the global conformal group 45

4.3.1 Consequences for OPEs 46  
4.3.2 Finding the quasiprimaries in an OPE 46

4.4 Consequences of the full conformal group 48

4.4.1 Restrictions of conformal covariance on the OPEs 48  
4.4.2 Virasoro descendants in a quasiprimary basis 49  
4.4.3 Virasoro descendants of the unit operator 51  
4.4.4 Finding the primaries in an OPE 52

4.5 An overview of  $\mathcal{W}$ -algebras 52

4.5.1 Direct construction 52  
4.5.2 Subalgebras of known  $\mathcal{W}$  -algebras 53  
4.5.3 Constructing a  $\mathcal{W}$ -algebra via a realization 53  
4.5.4 Superconformal algebras 54  
4.5.5 Attempts towards a classification 54

4.6 An example:  $\mathcal{W}_cB_2$  54

4.6.1 The  $\mathcal{W}_cB_2$ -algebra 54  
4.6.2 Coulomb Gas Realisation 55  
4.6.3 Highest weight representations 57

4.7 Discussion 57  
4.A Appendix 57

# 5 Factoring out Free Fields 62

5.1 Algorithms for factorisation 62

5.1.1 Free fermions 62  
5.1.2 Symplectic bosons 63  
5.1.3  $U(1)$  currents 64

5.2 Generating functionals 64  
5.3 Examples 65

5.3.1  $N = 3$  superconformal algebras 65  
5.3.2  $N = 4$  superconformal algebras 68

5.4 Discussion 70

# 6  $\mathcal{W}$  -algebras and gauged WZNW models 72

6.1 Classical Drinfeld-Sokolov reduction 72  
6.2 Quantum reduction 74

6.2.1 Gauged WZNW Model 74  
6.2.2 The induced action 75

6.3 An example:  $osp(N|2)$  76  
6.4 Cohomology 78

6.4.1 Computing the quantum cohomology 78  
6.4.2 Classical cohomology 80  
6.4.3 The quantum  $\mathcal{W}$  -algebra 81

6.5 An example:  $osp(N|2)$ , continued 82  
6.6 Quantum corrections to the extended action 82  
6.7 Discussion 83

# 7 Renormalisation factors in  $\mathcal{W}$  -Gravity 85

7.1 Introduction 85  
7.2 All-order results 86  
7.3 An example:  $so(N)$ -supergravity 89  
7.4 Semiclassical evaluation for the linear superconformal algebras 91

7.4.1  $N = 1$  93  
7.4.2  $N = 2$  94  
7.4.3  $N = 3,4$  95

7.5 Discussion 95

# 8 Critical  $\mathcal{W}$  -strings 98

8.1 The bosonic string 98  
8.2  $\mathcal{W}$  -strings 99  
8.3  $\mathcal{W}_{2,s}$  -strings 101  
8.4 Minimal models and  $\mathcal{W}_{2,s}$  -strings 103

8.4.1 The  $\mathcal{W}_{2,4}$  -string 104  
8.4.2 The  $\mathcal{W}_{2,5}$  string 105  
8.4.3 The  $\mathcal{W}_{2,6}$  string 106

8.5 Hierarchies of string embeddings 107  
8.6 Conclusion and discussion 110

# A Green's function for the Laplacian in 2d 112

# B Superconventions 114

C A Mathematica primer

Bibliography

# Chapter 1

# Introduction and outline

The first chapter of a recent PhD. thesis in contemporary high-energy physics necessarily stresses the importance of symmetry [42, 51, 45]. The reason for this is that symmetry is the most powerful organising principle available, and a theoretical physicist wants to assume as little as possible. This has the peculiar consequence that he or she ends up making the far-reaching assumption that "nature" has the largest symmetry we are able to find.

A striking example is provided by string theory. The universe seems to contain a large number of "elementary" particles. It is an appealing idea to think of these particles as different states of one single object. This would enable us to treat them in a symmetric way. The simplest objects in every-day experience which have such different eigenstates are (violin) strings. One then has to find which action governs a string-object moving through space-time. The simplest (in a certain sense) action, was found by Polyakov [158]. It is a generalisation of the action for a free relativistic particle. For a bosonic string in  $D$  dimensions the action is given by:

$$
S [ X ^ {\mu}, g ^ {i j} ] = - \frac {1}{4 \pi T} \int d x ^ {2} \sqrt {g (x)} g ^ {i j} (x) \partial_ {i} X ^ {\mu} (x) \partial_ {j} X _ {\mu} (x), \tag {1.1}
$$

where the fields  $X^{\mu}$  describe the position of the string, and  $x$  are coordinates which parametrise the two-dimensional surface ("world-sheet") which is swept out by the string as it moves in  $D$ -dimensional (flat) space-time.  $g^{ij}$  is the metric on the world-sheet, with inverse determinant  $g$ .  $T$  is a parameter that is related to the string tension.

The action (1.1) has (classically) a very large symmetry group, corresponding to reparametrisations of the world-sheet, and rescalings of the metric  $g^{ij}$ . These invariances are quite natural from the point of view of string theory. When viewing the theory defined by eq. (1.1) as a field theory in two dimensions, a first surprise awaits us. The field theory has an infinite dimensional symmetry group, which was quite uncommon in those days. A second surprise arises when we quantise the bosonic string theory. Requiring that the symmetry survives quantisation fixes the

number of space-time dimensions to 26. Somehow, this makes one hope that a more realistic string theory would "explain" why we are living in a four-dimensional world. The third surprise is that, while we started with a free theory, interactions seem also to be fixed by the action (1.1), with only one parameter  $T$ . This is in contrast with the grand-unified theories where interactions have to be put in by hand, requiring the introduction of a number of parameters that have to be fixed by comparing with experiments. A last surprise which we wish to mention, is that the spectrum of the physical states contains a particle with the correct properties for a graviton. String theory thus seems to incorporate quantum gravity. This is particularly fortunate because no other theory has been found yet which provides a consistent quantisation of gravity (in four dimensions).

These four features - symmetry, fixing the number of dimensions, "automatic" interactions and quantum gravity - were so attractive that many physicists decided to put the book of Popper [166] back on the shelf for a while. Indeed, although string theory certainly looks like a "good" theory, it still does not produce any results which are falsifiable, i.e. which can be contradicted by an experiment.

To make any chance of being a realistic theory, a number of flaws of the original bosonic string theory (like apparently giving the wrong space-time dimension) had to be resolved. Several roads can be followed, but we will concentrate here on the one which is most related to this thesis: enlarging the symmetry of the theory. In fact, the Polyakov action has many more symmetries than we alluded to. However, they are only global symmetries, i.e. generated by transformations with constant parameters. To make some of these symmetries local, one has to introduce extra gauge fields, which can be viewed as generalisations of the metric  $g^{ij}$  in eq. (1.1). In some cases, extra fields comparable to  $X^{\mu}(x)$  are added to the theory. For instance, in superstrings one uses fermionic fields describing coordinates in a Grassmann manifold. The resulting string theories are called "W-strings", and the (infinite dimensional) algebra formed by the infinitesimal transformations of the enlarged group is called an extended conformal algebra or, loosely speaking, a "W-algebra".

Unfortunately, during the process of enhancing the original bosonic string, one of its attractive features has been lost, namely its uniqueness. This is due to a number of reasons, but we will only mention two. Because an infinite number of  $\mathcal{W}$ -algebras exist, an infinite number of  $\mathcal{W}$ -string theories can be found (although certainly not all of them are candidates for a realistic theory). A second reason is that by adding an extra action for the metric to eq. (1.1), one can make a consistent quantum theory for other dimensions of space-time than 26, the "noncritical" strings. As it is well-known since Einstein that the metric is related to gravity, the study of consistent quantum actions for the metric provides a quantisation of gravity in two dimensions. Two-dimensional  $\mathcal{W}$ -gravity is interesting in its own respect because one hopes to gain some insight in how to construct consistent quantum gravity in four dimensions.

Although uniqueness has been lost, the other attractive features of string theory still survive. In particular, the symmetry group of the Polyakov action has even been enlarged. The study of this new kind of symmetry has influenced, and has been influenced by, many other branches of physics and mathematics. This has happened quite often in the history of string theory, and is sometimes regarded as an important motivation for studying strings. To be able to discuss the relation of  $\mathcal{W}$ -algebras to other fields in physics, we have to be somewhat more precise.

The Polyakov action (1.1) is invariant under general coordinate transformations  $x^i \to f^i(x)$  and local Weyl rescalings of the metric  $g_{ij}(x) \to \Lambda(x)g_{ij}(x)$ . A generalisation of the latter is to allow other fields  $\Phi(x)$  to rescale as  $\Phi(x) \to \Lambda(x)^h\Phi(x)$ .  $h$  is called the scaling dimension of the field  $\Phi(x)$ . The combination of the general coordinate and Weyl transformations can be used to gauge away components of the metric. In two dimensions one has exactly enough parameters to put the metric equal, at least locally, to the flat metric  $\eta_{ij}$  (the conformal gauge). However, this does not yet completely fix the gauge. Obviously, conformal transformations (coordinate transformations which scale the metric) combined with the appropriate Weyl rescaling form a residual symmetry group. Therefore, field theories which have general coordinate and Weyl invariance are called conformal field theories.

The situation in two dimensions is rather special. In light-cone coordinates,  $x_{\pm} = x_0 \pm x_1$ , every transformation  $x_{\pm} \rightarrow f^{\pm}(x_{\pm})$  is conformal. We see that the group formed by the conformal transformations is infinite dimensional in two dimensions. This makes clear why the symmetry group of string theory is so exceptionally large. The conformal transformations are generated by the energy-momentum tensor  $T^{ij}$  of the theory, which has scaling dimension  $h = 2$ . In fact, the algebra splits in two copies of the Virasoro algebra, related to the  $x_{+}$  and  $x_{-}$  transformations, and generated by  $T^{++}$  and  $T^{- - }$ . Similarly, an extended conformal algebra is formed by two copies of what is called a  $\mathcal{W}$ -algebra.

The symmetries of a theory have direct consequences for its correlation functions. The Ward identities are relations between  $N$ -point correlation functions where one of the fields is a symmetry generator, and  $(N - 1)$ -point functions. Usually, this does not yet fix the  $N$ -point function, and correlation functions have to be calculated tediously.

In two-dimensional conformal field theory, the consistency conditions imposed by the symmetries are so strong that the Ward identities determine all correlation functions with a symmetry generator in terms of those without any symmetry generators. This means that once the Ward identities have been found (which requires some regularisation procedure), all correlation functions with symmetry generators can be recursively computed.

In renormalisable field theories, "Operator Product Expansions" (OPEs) are introduced to calculate the short-distance behaviour of correlation functions [208]. This formalism has been extended by Belavin, Polyakov and Zamolodchikov [13] to two-dimensional conformal field theory. The Ward identities fix OPEs. Moreover, they impose a set of consistency conditions on the OPEs such that computing with OPEs amounts to applying a set of algebraic rules. They even almost determine the form of the OPEs. As an example, we give the OPE of one of the components of the energy-momentum tensor, which for any conformal field theory is:

$$
T ^ {+ +} (x) T ^ {+ +} (y) = \frac {c / 2}{\left(x ^ {-} - y ^ {-}\right) ^ {4}} + \frac {2 T ^ {+ +} (y)}{\left(x ^ {-} - y ^ {-}\right) ^ {2}} + \frac {\partial_ {-} T ^ {+ +} (y)}{\left(x ^ {-} - y ^ {-}\right)} + O \left(x ^ {-} - y ^ {-}\right) ^ {0}. \tag {1.2}
$$

Here, the "central charge"  $c$  is a number which can be determined by computing the two-point function  $< T^{++}(x)T^{++}(y)>$ , and depends on the theory we are considering. The important point here is that once the central charge is known, all correlation functions of  $T^{++}$  can be algebraically computed. From the OPE eq. (1.2), the Virasoro algebra can be derived and vice versa. Similarly, if the symmetry algebra of the theory forms an extended conformal algebra, the generators form an Operator Product Algebra. This contains exactly the same information as the  $\mathcal{W}$ -algebra, and indeed is often called a  $\mathcal{W}$ -algebra.

There exists by now a wealth of examples of  $\mathcal{W}$ -algebras. Among the best known are the affine Lie algebras and the linear superconformal algebras. When a  $\mathcal{W}$ -algebra contains a generator with scaling dimension larger than two, the  $\mathcal{W}$ -algebra is (in most cases) nonlinear. Some examples of such  $\mathcal{W}$ -algebras with only one extra generator are  $\mathcal{W}_3$  [211], the spin 4 algebra [29, 108] and the spin 6 algebra [74]. The Bershadsky-Knizhnik algebras [20, 133] have  $N$  supersymmetry generators and an affine  $so(N)$  subalgebra. Many other examples exist and no classification of  $\mathcal{W}$ -algebras seems as yet within reach.

The fields of a conformal field theory form a representation of its symmetry algebra. Belavin, Polyakov and Zamolodchikov [13] showed that under certain assumptions all fields of the theory are descendants of a set of primary fields. For a certain subclass of conformal field theories, the "minimal" models, the Ward identities fix all correlation functions. They also showed that the simplest minimal model corresponds to the Ising model at criticality. This connection with statistical mechanics of two-dimensional systems is due to the fact that a system becomes invariant under scaling transformations at the critical point of a phase transition. By hypothesising local conformal invariance, various authors (see [4, 118]) found the critical exponents

of many two-dimensional models. Some examples of statistical models are the Ising  $(m = 3)$ , tricritical Ising  $(m = 4)$ , 3-state Potts  $(m = 5)$ , tricritical 3-state Potts  $(m = 6)$ , and the Restricted Solid-on-Solid (any  $m$ ) models, where we denoted the number of the corresponding unitary Virasoro minimal model in brackets.

Another important connection was found, not in statistical mechanics, but in the study of integrable models in mathematics and physics. These models have two different Hamiltonian structures, whose Poisson brackets form examples of classical  $\mathcal{W}$ -algebras. For example, the Korteweg-de Vries (KdV) equation gives rise to a Virasoro Poisson bracket, while the Boussinesq equation has a Hamiltonian structure which corresponds to the classical  $\mathcal{W}_3$  algebra. Moreover, the relation between the two Hamiltonian structures gives rise to a powerful method of constructing classical  $\mathcal{W}$ -algebras. Drinfeld and Sokolov [60] found a hierarchy of equations of the KdV type based on the Lie algebras  $sl(N)$ . The Hamiltonian structures of these equations provide explicit realisations of the classical  $\mathcal{W}_N$  algebras, whose generators have scaling dimensions  $2, 3, \ldots, N$ . An extension of this method is still the most powerful way at our disposal to find (realisations of)  $\mathcal{W}$ -algebras.

This thesis is organised as follows. Chapter 2 gives an introduction to conformal field theory. We discuss how the Ward identities are derived. For this purpose, the Operator Product Expansion formalism is introduced. We determine the complete set of consistency conditions on the OPEs. We then show how an infinite dimensional Lie algebra, corresponding to the symmetry algebra of the conformal theory, can be found using OPEs. We define the generating functional of the correlation functions of symmetry generators. The Ward identities can be used to find functional equations for the generating functional (or induced action). We conclude the chapter with some important examples of conformal field theories: free-field theories and WZNW-models.

In chapter 3, the consistency conditions on OPEs are converted to a set of algorithms to compute with OPEs, suitable for implementation in a symbolic manipulation program. We then describe the Mathematica package OPErefs we developed. This package completely automates the computation of OPEs (and thus of correlation functions), given the set of OPEs of the generators of the  $\mathcal{W}$ -algebra.

The next chapter discusses  $\mathcal{W}$ -algebras using the Operator Product Expansion formalism. We first give some basic notions on highest weight representations of  $\mathcal{W}$ -algebras, of which minimal models are particular examples. The fields of a conformal field theory assemble themselves in highest weight representations. We then analyse the structure of  $\mathcal{W}$ -algebras using the consistency conditions found in chapter 2. The global conformal transformations fix the form of OPEs of quasiprimary fields, which are special examples of highest weight fields with respect to the global conformal algebra. A similar analysis is done for the local conformal transformations and primary fields. We then discuss the different methods which are used to construct  $\mathcal{W}$ -algebras and comment on the classification of the  $\mathcal{W}$ -algebras. Finally, as an example of the ideas in this chapter, the  $\mathcal{W}_cB_2$  algebra is studied in detail. The

complexity of the calculations shows the usefulness of OPErefs.

Goddard and Schwimmer [103] proved that free fermions can always be factored out of a  $\mathcal{W}$ -algebra. Chapter 5 extends this result to arbitrary free fields. This is an important result as it shows that a classification of  $\mathcal{W}$ -algebras need not be concerned with free fields. We provide explicit algorithms to perform this factorisation. We then show how the generating functionals of the  $\mathcal{W}$ -algebra obtained via factorisation are related to the generating functional of the original  $\mathcal{W}$ -algebra. The  $N = 3$  and  $N = 4$  linear superconformal algebras are treated as examples.

The Drinfeld-Sokolov method constructs a realisation for a classical  $\mathcal{W}$ -algebra via imposing constraints on the currents of a Kac-Moody algebra. In particular, for any embedding of  $sl(2)$  in a semi-simple (super)Lie algebra a different  $\mathcal{W}$ -algebra results. In chapter 6, Drinfeld-Sokolov reduction is extended to the quantum case. The reduction is implemented in an action formalism using a gauged WZNW-model, which enables us to find a path integral formulation for the induced action of the  $\mathcal{W}$ -algebra. The gauge fixing is performed using the Batalin-Vilkovisky [11] method. In a special gauge, the BV procedure reduces to a BRST approach [138]. Operator Product Expansions are used to perform a BRST quantisation. The cohomology of the BRST operator is determined in both the classical and quantum case. The results are then used to show that we indeed constructed a realisation of a quantum  $\mathcal{W}$ -algebra. The  $N$ -extended  $so(N)$  superconformal algebras [20, 133] are used as an illustration of the general method.

The results of chapter 6 are then used in chapter 7 to study  $\mathcal{W}$ -gravity theories in the light-cone gauge. The gauged WZNW-model is used as a particular matter sector for the coupling to  $\mathcal{W}$ -gravity. Using the path integral formulation of the previous chapter, the effective action can be computed. It is shown that the effective action can be obtained from its classical limit by simply inserting some renormalisation factors. Explicit expressions for these renormalisation factors are given. They contain the central charge of the  $\mathcal{W}$ -algebra and some parameters related to the (super)Lie algebra and the particular  $sl(2)$ -embedding for which a realisation of the  $\mathcal{W}$ -algebra can be found. The example of the previous chapter is used to construct the effective action of  $so(N)$  supergravity. We then check the results using the correspondence between the linear superconformal algebras and the  $so(N)$ $\mathcal{W}$ -algebras for  $N \leq 4$ , established in chapter 5. This is done using a semiclassical evaluation of the effective action for the linear superconformal algebras.

The last chapter contains a discussion of critical  $\mathcal{W}$ -strings. After a short review of  $\mathcal{W}$ -string theory, we concentrate to the case where the classical  $\mathcal{W}$ -algebra is formed by the energy-momentum tensor and a dimension  $s$  generator. These  $\mathcal{W}_{2,s}$  strings provide examples which can be analysed in considerably more detail than string theories based on more complicated  $\mathcal{W}$ -algebras. In particular, operator product expansions (and OPErefs) are used to provide some insight in the appearance of  $\mathcal{W}$ -minimal models in the spectrum of  $\mathcal{W}$ -strings. Finally, some comments are made on the recent developments initiated by Berkovits and Vafa [14]. They

showed how the bosonic string can be viewed as an  $N = 1$  superstring with a particular choice of vacuum, and a similar embedding of  $N = 1$  into  $N = 2$  superstrings. The hope arises that a hierarchy of string embeddings exists, restoring the uniqueness of string theory in some sense.

# Chapter 2

# Conformal Field Theory and Operator Product Expansions

This chapter gives an introduction to conformal field theory with special attention to Operator Product Expansions (OPEs). It is of course not complete as conformal field theory is a very wide subject, and many excellent reviews exist, e.g. [97, 122, 190]. Because it forms an introduction to the subject, some points are probably trivial for someone who feels at home in conformal field theory. However, some topics are presented from a new standpoint, a few new results (on the associativity of OPEs) are given, and notations are fixed for the rest of the work.

We start by introducing the conformal transformations. Then, we study the consequences of a symmetry of a conformal field theory on its correlation functions. In the quantum case, this information is contained in the Ward identities. These identities are then used in the third section to develop the OPE formalism. We study the consistency conditions for OPEs in detail and introduce the notion of an Operator Product Algebra (OPA). We then draw attention to the close analogy between OPEs and Poisson brackets. Section 2.4 defines the mode algebra of the symmetry generators. In the next section, we define the generating functionals of the theory and show how they are determined by the Ward identities. We conclude with some important examples of conformal field theories, free field theories and WZNW-models.

# 2.1 Conformal transformations

A conformal field theory is a field theory which is invariant under general coordinate transformations  $x^{i} \to f^{i}(x)$  and under the additional symmetry of Weyl invariance. The latter transformations correspond to local scale transformations of the metric  $g_{ij}(x) \to \Lambda(x)g_{ij}(x)$  and fields  $\Phi(x) \to \Lambda(x)^{h}\Phi(x)$ , where  $h$  is the scaling dimension of the field  $\Phi$ . The combination of these symmetries can be used to gauge away components of the metric. In two dimensions one has exactly enough parameters to put the metric equal, at least locally, to the Minkowski metric  $\eta_{ij}$  (the conformal

gauge). Obviously, conformal transformations (coordinate transformations which scale the metric) combined with the appropriate Weyl rescaling form a residual symmetry group. Therefore, we will study this conformal group first.

Conformal transformations are coordinate transformations which change the metric with a local scale factor. In a space-time of signature  $(p,q)$  they form a group isomorphic with  $SO(p + 1,q + 1)$ . However, in the complex plane it is well-known that all (anti-)analytic transformations are conformal. This extends to the Minkowski plane where in light-cone coordinates, the conformal transformations are given by:

$$
x ^ {\pm} \longrightarrow x ^ {\prime \pm} \left(x ^ {\pm}\right). \tag {2.1.1}
$$

In a space-time of signature  $(-1,1)$ , it is customary to perform a Wick rotation. We will always assume this has been done, and treat only the Euclidean case.

For a space of Euclidean signature, it is advantageous to use a complex basis  $(\tau + i\sigma, \tau - i\sigma)$ . In string theory, the space in which these coordinates live is a cylinder, as  $\sigma$  is used as a periodic coordinate. This also applies to two-dimensional statistical systems with periodic boundary conditions in one dimension. One then maps this cylinder to the full complex plane with coordinates

$$
(z, \bar {z}) = (\exp (\tau + i \sigma), \exp (\tau - i \sigma), \tag {2.1.2}
$$

where we will take a flat metric proportional to  $\delta_{ij}$  in the real coordinates, or in the complex coordinates:

$$
d s ^ {2} = 2 \sqrt {g} d z d \bar {z}. \tag {2.1.3}
$$

We will use the notation  $x$  for a coordinate of a point in the complex plane. Note that the points with fixed time  $\tau$  lie on a circle in the  $(z, \bar{z})$  plane.

It is often convenient to consider  $z$  and  $\bar{z}$  as independent coordinates (i.e. not necessarily complex conjugate). We can then restrict to the Euclidean plane by imposing  $\bar{z} = z^{*}$ . The plane with a Minkowski metric corresponds to  $z, \bar{z} \in \mathbb{R}$ .

As we are working in Euclidean space, the complex plane can be compactified to the Riemann sphere. Conformal field theory can also be defined on arbitrary Riemann surfaces, but we will restrict ourselves in this work to the sphere.

In the complex coordinates, a conformal transformation is given by:

$$
z \rightarrow z ^ {\prime} = f (z), \quad \bar {z} \rightarrow \bar {z} ^ {\prime} = \bar {f} (\bar {z}), \tag {2.1.4}
$$

where  $f$  is an analytic function, and  $\bar{f}$  is antianalytic.

Definition 2.1.1 A primary field transforms under the conformal transformation (2.1.4) as:

$$
\Phi (z, \bar {z}) \rightarrow \Phi^ {\prime} (f (z), \bar {f} (\bar {z})) = (\partial f (z)) ^ {- h} (\bar {\partial} \bar {f} (\bar {z})) ^ {- \bar {h}} \Phi (z, \bar {z}), \tag {2.1.5}
$$

where  $\partial \equiv \frac{d}{dz}$  and  $\bar{\partial} \equiv \frac{d}{d\bar{z}}$ . The numbers  $h$  and  $\bar{h}$  are called the conformal dimensions of the field  $\Phi$ .

For infinitesimal transformations of the coordinates  $f(z) = z - \varepsilon (z)$ , we see that the primary fields transform as:

$$
\delta_ {\varepsilon} \Phi (z, \bar {z}) = \varepsilon (z) \partial \Phi (z, \bar {z}) + h \partial \varepsilon (z) \Phi (z, \bar {z}). \tag {2.1.6}
$$

By choosing for  $\varepsilon(z)$  any power of  $z$  we see that the conformal transformations form an infinite algebra generated by:

$$
l _ {m} = - z ^ {m + 1} \partial \quad \text {a n d} \quad \bar {l} _ {m} = - \bar {z} ^ {m + 1} \bar {\partial}, \quad m \in \mathbb {Z}, \tag {2.1.7}
$$

which consists of two commuting copies of the Virasoro algebra, but without central extension (see further):

$$
\left[ l _ {m}, l _ {n} \right] = (m - n) l _ {m + n}, \tag {2.1.8}
$$

and analogous commutators for the  $\bar{l}_m$ .

Clearly,  $l_0$  corresponds to scaling transformations in  $z$ . The combination  $l_0 + \bar{l}_0$  generates scaling transformations in the complex plane  $x \to \lambda x$ , while  $i(l_0 - \bar{l}_0)$  generates rotations. This means that a field  $\Phi(x)$  with conformal dimensions  $h$  and  $\bar{h}$  has scaling dimension  $h + \bar{h}$  and spin  $|h - \bar{h}|$ .  $l_{-1}$  generates translations and  $l_1$  "special" conformal transformations. The subalgebra formed by  $\{l_{-1}, l_0, l_1\}$  corresponds to the globally defined and invertible conformal transformations:

$$
z \rightarrow \frac {a z + b}{c z + d}, \tag {2.1.9}
$$

where  $a, b, c, d \in \mathbb{C}$  and  $ac - bd = 1$ . These transformations form a group isomorphic to  $SL(2,\mathbb{C}) \approx SO(3,1)$ .<sup>1</sup>

Definition 2.1.2 A quasiprimary field transform as eqs. (2.1.5,2.1.6) for the global transformations. A scaling field transforms this way for translations and scaling transformations.

# 2.2 Correlation functions and symmetry

The physics of a quantum field theory is contained in the correlation functions. In a path integral formalism, a correlation function of fields  $\Phi_{i}$  corresponding to observables can be symbolically written as:

$$
<   \mathcal {R} \left\{\Phi_ {1} \left(x _ {1}\right) \Phi_ {2} \left(x _ {2}\right) \dots \right\} > = \frac {1}{\mathcal {N}} \int \left[ d \varphi_ {k} \right] \exp \left(- S \left[ \varphi_ {j} \right]\right) \Phi_ {1} \left(x _ {1}\right) \Phi_ {2} \left(x _ {2}\right) \dots , \tag {2.2.1}
$$

where  $S$  is a suitable action, a functional of the fields  $\varphi_{j}$  in the theory, and  $[d\varphi_k]$  denotes an appropriate measure.  $\mathcal{R}$  denotes time-ordering in a radial quantisation scheme [92], i.e.  $|x_i| > |x_{i + 1}|$ . We will drop this symbol for ease of notation. In fact, we will assume that for a correlation function  $< \mathcal{R}\{\Phi_1(x_1)\Phi_1(x_2)\ldots \} >$ , the expression for  $|x_{1}| < |x_{2}|$  can be obtained by analytic continuation from the expression for  $|x_{2}| < |x_{1}|$ . This property is commonly called "crossing symmetry". It can be checked for free fields, and we restrict ourselves to theories where it is true.

The normalisation constant in eq. (2.2.1) is given by:

$$
\mathcal {N} = \int \left[ d \varphi_ {k} \right] \exp \left(- S \left[ \varphi_ {j} \right]\right). \tag {2.2.2}
$$

We suppose that  $\mathcal{N}$  is different from zero (vacuum normalised to one). In general, the pathintegral (2.2.1) has to be computed perturbatively.

Symmetries of the theory put certain restrictions on the form of the correlation functions. We will investigate this in the next subsections.

# 2.2.1 Global conformal invariance

In this subsection, we consider the consequences of the invariance of the correlation functions of quasiprimary fields under the global conformal transformations (2.1.9). The results can be extended to conformal field theories in an arbitrary number of dimensions, but we treat only the two-dimensional case.

In conformal field theory one generally requires translation, rotation and scaling invariance of the correlation functions. This restricts all one-point functions of fields with zero conformal dimensions ( $h = \bar{h} = 0$ ) to be constant, and all others to be zero. Two-point functions have the following form:

$$
<   \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) > = \frac {\operatorname {c s t} _ {1 2}}{\left(z _ {1} - z _ {2}\right) ^ {h _ {1} + h _ {2}} \left(\bar {z} _ {1} - \bar {z} _ {2}\right) ^ {\bar {h} _ {1} + \bar {h} _ {2}}} , \tag {2.2.3}
$$

with  $\mathrm{cst}_{12}$  a constant. If one requires global conformal invariance (i.e. invariance under the special transformations as well), only fields with the same conformal dimensions can have non-zero two-point functions, see also eq. (2.2.9). Under this

# 2.2. Correlation functions and symmetry

assumption, three-point functions are restricted to:

$$
\begin{array}{l} <   \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) \Phi_ {3} \left(z _ {3}, \bar {z} _ {3}\right) > = \\ \operatorname{cst}_{123}\prod_{\substack{\{ijk\} \\ i <   j}}\frac{1}{(z_{i} - z_{j})^{h_{ijk}}(\bar{z}_{i} - \bar{z}_{j})^{\bar{h}_{ijk}}}, \tag{2.2.4} \\ \end{array}
$$

where  $h_{ijk} = h_i + h_j - h_k$  and the product goes over all permutations of  $\{1,2,3\}$ . It is also possible to find the restrictions from global conformal invariance on four (or more)-point functions. It turns out that the correlation function is a function of the harmonic ratios of the coordinates involved, apart from  $(z_i - z_j)$  factors like in (2.2.4) to have the correct scaling law. Crossing symmetry has to be checked for a four-point function, while it is automatic for two- and three-point functions.

Note that imposing invariance under the local conformal transformations would make all correlation functions zero.

# 2.2.2 Ward identities

In this subsection, we will discuss the consequences of a global symmetry of the action using Ward identities.

We start by considering an infinitesimal transformation  $\delta_{\varepsilon}$  of the fields, where  $\varepsilon(x)$  is an infinitesimal parameter. Throughout this work, we will only consider local transformations, i.e.  $\delta_{\varepsilon} \Phi(x)$  depends on a finite number of fields and their derivatives. If the transformation leaves the measure of the path integral (2.2.1) invariant, we find the following identity:

$$
\begin{array}{l} \int \left[ d \varphi_ {k} \right] \exp (- S [ \varphi_ {j} ]) \delta_ {\varepsilon} (\Phi_ {1} (x _ {1}) \Phi_ {2} (x _ {2}) \dots) = \\ \int \left[ d \varphi_ {k} \right] \exp \left(- S \left[ \varphi_ {j} \right]\right) \left(\delta_ {\varepsilon} S [ \varphi ]\right) \Phi_ {1} \left(x _ {1}\right) \Phi_ {2} \left(x _ {2}\right) \dots . \tag {2.2.5} \\ \end{array}
$$

This follows by considering a change of variables in the pathintegral. Using eq. (2.2.1), eq. (2.2.5) implies an identity between correlation functions. This Ward identity is useful in many cases, but at this moment we are interested in transformations  $\delta_{\varepsilon}$  which leave the action invariant if  $\delta_{\varepsilon}$  is constant, i.e. global symmetries.

As an example, we will derive the Ward identity for the energy-momentum tensor  $T^{ij}$ . Consider an infinitesimal transformation of the fields of the form:

$$
\delta_ {\varepsilon} \Phi (x) = \varepsilon^ {i} (x) \partial_ {t} \Phi (x) + \dots , \tag {2.2.6}
$$

where the dots denote corrections according to the tensorial nature of the field  $\Phi$ , see (2.1.6). The action is supposed to be invariant under these transformations if

the  $\varepsilon^i (x)$  are constants. Using Noether's law, one has a classically conserved current associated to this symmetry of the action  $S$ :

$$
\delta_ {\varepsilon} S = - \frac {1}{\pi} \int d ^ {2} x \sqrt {g} T _ {j} ^ {i} (x) \partial_ {i} \varepsilon^ {j} (x), \tag {2.2.7}
$$

where  $g$  is the absolute value of the determinant of  $g_{ij}$ . This equation defines the energy-momentum tensor  $T^{ij}$  up to functions whose divergence vanishes. One can show that the alternative definition

$$
T _ {i j} = \frac {- 2 \pi}{\sqrt {g}} \frac {\delta S}{\delta g ^ {i j}} \tag {2.2.8}
$$

satisfies (2.2.7). It is obvious that the tensor (2.2.8) is symmetric. If the action is invariant under Weyl scaling, we immediately see that it is traceless. Hence, in a conformal field theory in two dimensions, the energy-momentum tensor has only two independent components. In the complex basis, the components that remain are  $T^{zz}$  and  $T^{\bar{z}\bar{z}}$ . Using the metric eq. (2.1.3), we see from eq. (2.2.7) that Weyl invariance implies that the action is not only invariant under a global transformation, but also under a conformal transformation,  $\bar{\partial}\varepsilon^{z} = 0$ ,  $\varepsilon^{\bar{z}} = 0$ . This transformation is sometimes called semi-local.

Let us now consider the expectation value of some fields  $\Phi_k$ . Combining eqs. (2.2.5) and (2.2.7), we find:

$$
\begin{array}{l} <   \left(\delta_ {\varepsilon} \Phi_ {1} (x _ {1})\right) \Phi_ {2} (x _ {2}) \dots > + <   \Phi_ {1} (x _ {1}) \left(\delta_ {\varepsilon} \Phi_ {2} (x _ {2})\right) \dots > + \dots \\ = - \frac {\sqrt {g}}{\pi} \int d ^ {2} x (\partial_ {i} \varepsilon^ {j} (x)) <   T _ {j} ^ {i} (x) \Phi_ {1} \left(x _ {1}\right) \Phi_ {2} \left(x _ {2}\right) \dots >. \tag {2.2.9} \\ \end{array}
$$

The equation (2.2.9) is the Ward identity for the energy-momentum tensor. It shows that  $T^{ij}$  is the generator of general coordinate transformations. Note that  $T^{ij}$  has to be symmetric and traceless for the correlation functions to be rotation and scaling invariant.

We derived the Ward identity (2.2.9) in a formal way. In a given theory, it has to be checked using a regularised calculation. We will not do this in this work, and assume that eq. (2.2.9) holds, possibly with some quantum corrections. This enables us to make general statements for every theory where eq. (2.2.9) is valid. Similar Ward identities can be derived for any global symmetry of the action which leaves the measure invariant.

An important corollary of the Ward identity for a symmetry generator, is that it is conserved "inside" correlation functions. We will again treat  $T^{ij}$  as an example. We take  $\varepsilon^i(x)$  functions which go to zero when  $x$  goes to infinity, and which have no singularities. In this case, we can use partial integration. Now, the lhs of (2.2.9)

depends only on the value of  $\varepsilon$  (and its derivatives) in the  $x_{i}$ . This means that, after a partial integration, the coefficient of  $\varepsilon$  in the integrandum in the rhs of (2.2.9) has to be zero except at these points. We find:

$$
\frac {\partial}{\partial x ^ {i}} <   T ^ {i j} (x) \Phi_ {1} \left(x _ {1}\right) \Phi_ {2} \left(x _ {2}\right) \dots > = 0 \tag {2.2.10}
$$

where  $x\neq x_{k}$

From now on, we consider the case of a conformal field theory in two dimensions. We already showed that  $T^{ij}$  is symmetric and traceless, in the sense that all correlation functions vanish if one of the fields is the antisymmetric part of  $T^{ij}$  or its trace. In the complex basis the two remaining components are  $T^{\bar{z}\bar{z}}$  and  $T^{zz}$ . The conservation of  $T^{ij}$  eq. (2.2.10) gives:

$$
\frac {d}{d \bar {z}} <   T ^ {\bar {z} \bar {z}} (z, \bar {z}) \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) \dots > = 0, \quad z \neq z _ {i}, \tag {2.2.11}
$$

and an analogous equation for  $T^{zz}$ . This shows that  $< T^{\bar{z}\bar{z}}(z,\bar{z})\Phi_1\Phi_2\dots >$  is a holomorphic function of  $z$  with possible singularities in  $z_{i}$ . Using the assumption that the transformation of the fields  $\delta_{\varepsilon}\Phi_{i}(x_{i})$  is local, we see that the correlation function is a meromorphic function, with possible poles in  $z = z_{i}$ . We define:

$$
T \equiv g T ^ {\bar {z} \bar {z}} = T _ {z z} \quad \text {a n d} \quad \bar {T} \equiv g T ^ {z z} = T _ {\bar {z} \bar {z}} \tag {2.2.12}
$$

and we will write  $T$  as a function of  $z$  only.

For  $\varepsilon^i (x)$  as specified above eq. (2.2.10), only the singularities  $z = z_{i}$  contribute to the Ward identity (2.2.9). We now take  $\varepsilon^z = \varepsilon (z)$ ,  $\varepsilon^{\bar{z}} = 0$  in a neighbourhood of these points. Using eq. (A.5), the Ward identity (2.2.9) becomes for a conformal transformation:

$$
\begin{array}{l} \sum_ {j = 1} ^ {N} <   \Phi_ {1} (z _ {1}, \bar {z} _ {1}) \dots (\delta_ {\varepsilon} \Phi_ {j} (z _ {j}, \bar {z} _ {j})) \dots > \\ = \oint_ {C} \frac {d z}{2 \pi i} \varepsilon (z) <   T (z) \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) \dots >, \tag {2.2.13} \\ \end{array}
$$

where the contour  $C$  encircles all  $z_{j}$  exactly once anticlockwise (and an analogous formula for antianalytic variations).

If all  $\Phi_j$  are primary fields, we can use (2.1.6) to compute the  $lhs$  of the Ward identity (2.2.13):

$$
\begin{array}{l} \oint_ {C} \frac {d z}{2 \pi i} \varepsilon (z) <   T (z) \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) \dots > \\ = \sum_ {j = 1} ^ {N} \left(\varepsilon \left(z _ {j}\right) \partial_ {j} + h _ {j} \partial_ {j} \varepsilon \left(z _ {j}\right)\right) <   \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) \dots >. \tag {2.2.14} \\ \end{array}
$$

This fixes the singular part in  $(z = z_{i})$  of the meromorphic  $(N + 1)$ -point function  $< T\Phi_1\dots >$ . We will assume that all correlation functions have no singularity at infinity. Hence, the regular terms in  $(z - z_{i})$  vanish, except for a possible constant.

To know the transformation law of  $T(z)$  itself, we observe that since it is classically a rank two conformal tensor, its dimensions are  $h = 2$ ,  $\bar{h} = 0$ . However, due to quantum anomalies, a Schwinger term can arise in its transformation law:

$$
\delta_ {\varepsilon} T (z) = \varepsilon (z) \partial T (z) + 2 \partial \varepsilon (z) T (z) + \frac {c}{1 2} \partial^ {3} \varepsilon (z). \tag {2.2.15}
$$

In the path integral formalism, the Schwinger term is non-zero if the measure is not invariant under the symmetry transformation generated by  $T$ , see eq. (2.2.9). From dimensional arguments, it follows that  $c$  should have dimension zero. In fact,  $c$  is in general a complex number. When  $c$  is different from zero,  $T$  is not a primary field, but only quasiprimary, see definition 2.1.2.

The transformation law (2.2.15) involves only  $T$  and no other components of the energy-momentum tensor. In the classical case, this is a direct consequence of the conformal symmetry of the theory. We will assume that the same property holds in the quantum case.

Recall that eq. (2.2.14) determines the correlation function  $< T\Phi_1 \cdots >$  up to a constant. When we require invariance of the correlation function under scaling transformations  $(\varepsilon(z) \sim z)$ , this constant has to be zero. Under these assumptions, the Ward identity completely fixes all correlation functions with one of the fields equal to  $T$  and all other fields being primary:

$$
<   T (z) \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) \dots > =
$$

$$
\sum_ {j = 1} ^ {N} \left(\frac {h _ {j}}{(z - z _ {j}) ^ {2}} + \frac {1}{(z - z _ {j})} \partial_ {j}\right) <   \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) \dots >. \tag {2.2.16}
$$

We can now look at correlation functions of  $T$  only. Using eq. (2.2.15), we find in a completely analogous way to eq. (2.2.16):

$$
<   T (z) T \left(z _ {1}\right) T \left(z _ {2}\right) \dots > =
$$

$$
\sum_ {j = 1} ^ {N} <   T \left(z _ {1}\right) \dots \left(\frac {c / 2}{\left(z - z _ {j}\right) ^ {4}} + \frac {2 T \left(z _ {j}\right)}{\left(z - z _ {j}\right) ^ {2}} + \frac {\partial_ {j} T \left(z _ {j}\right)}{\left(z - z _ {j}\right)}\right) \dots T \left(z _ {N}\right) >. \tag {2.2.17}
$$

Let us consider a few examples. Because of scaling invariance, one-point functions of fields with non-zero dimension vanish. Together with the fact that  $< \mathbb{1} > = 1$ , we see that eq. (2.2.15) implies:

$$
<   \delta_ {\varepsilon} T > = \frac {c}{1 2} \partial^ {3} \varepsilon (z). \tag {2.2.18}
$$

# 2.3. Operator Product Expansions

The Ward identity (2.2.13) fixes then the two-point function to:

$$
<   T (z) T (w) > = \frac {c / 2}{(z - w) ^ {4}}, \tag {2.2.19}
$$

in accordance with eq. (2.2.3). Similarly, the constant in the three-point function (2.2.4) of three  $T$ 's is fixed to  $c$ .

The analysis of the Ward identity for the energy-momentum tensor can be repeated for every generator of a global symmetry of the theory. In two dimensions, all tensors can be decomposed in symmetric tensors. If they are in addition traceless, only two components remain  $W^{zzz\dots}$  and  $W^{\bar{z}\bar{z}\bar{z}\dots}$ . We again find that  $\bar{\partial} W^{\bar{z}\bar{z}\bar{z}\dots} = 0$ , and the correlation functions of  $W^{\bar{z}\bar{z}\bar{z}\dots}$  are meromorphic functions. For generators with a spinorial character, or generators with non-integer conformal dimensions, we can only conclude that the correlation functions are holomorphic, i.e. fractional powers of  $z - z_{i}$  could occur. We restrict ourselves in this work to the meromorphic sector of the conformal field theory. The symmetry algebra of the theory is a direct product $^3$ $\mathcal{A} \otimes \bar{\mathcal{A}}$  where  $\mathcal{A}$  is the algebra formed by the "chiral" generators (with holomorphic correlation functions), and  $\bar{\mathcal{A}}$  is its antichiral counterpart. From now on, we will concentrate on the chiral generators.

To conclude this section, we wish to stress that the Ward identities are regularisation dependent. In this work, we will not address the question of deriving the Ward identities for a given theory. However, once they have been determined, the Ward identities enable us to compute the singular part of any correlation function containing a symmetry generator in terms of (differential polynomials of) correlation functions of the fields of the theory. This determines those correlation functions up to a constant, because we required that a correlation function has no singularity at infinity. Furthermore, if the correlation functions are invariant under scale transformations, this extra constant can only appear when the sum of the conformal dimensions of all fields in the correlator is zero.

# 2.3 Operator Product Expansions

In this section, we discuss Operator Product Expansions (OPEs). In a first step, we show how they encode the information contained in the Ward identities. As such, OPEs provide an algebraic way of computing correlation functions. In the next subsection, we derive the consistency conditions on the OPEs. Subsection 2.3.3 introduces the concept of an Operator Product Algebra (OPA).

# 2.3.1 OPEs and Ward identities

To every field in a conformal field theory, we assign a unique element of a vectorspace  $\mathcal{V}$  of "operators". For the moment, we leave the precise correspondence open, but our goal is to define a bilinear operation in the vectorspace which enables us to compute correlation functions. We simply write the same symbol for the field and the corresponding operator. Similarly, if a field has conformal dimension  $h$ , we say that the corresponding operator has conformal dimension  $h$ .

Let us look at an example to clarify what we have in mind. Consider the Ward identities of the chiral component of the energy-momentum tensor  $T(z)$ . Eq. (2.2.16) suggests that we assign the OPE:

$$
T (z) \Phi \left(z _ {0}, \bar {z} _ {0}\right) = \frac {h \Phi \left(z _ {0} , \bar {z} _ {0}\right)}{\left(z - z _ {0}\right) ^ {2}} + \frac {\partial \Phi \left(z _ {0} , \bar {z} _ {0}\right)}{z - z _ {0}} + O \left(z - z _ {0}\right) ^ {0} \tag {2.3.1}
$$

for a primary field  $\Phi$  with conformal dimension  $h$ . We do not specify the regular part yet. Similarly, eq. (2.2.17) imposes the following OPE for  $T$  with itself:

$$
T (z) T \left(z _ {0}\right) = \frac {c / 2}{\left(z - z _ {0}\right) ^ {4}} + \frac {2 T \left(z _ {0}\right)}{\left(z - z _ {0}\right) ^ {2}} + \frac {\partial T \left(z _ {0}\right)}{z - z _ {0}} + O \left(z - z _ {0}\right) ^ {0}. \tag {2.3.2}
$$

We see that an OPE is a bilinear map from  $\mathcal{V} \otimes \mathcal{V}$  to the space of formal Laurent expansions in  $\mathcal{V}$ . We will use the following notation for OPEs:

$$
T _ {1} (z) T _ {2} \left(z _ {0}\right) = \sum_ {n <   = h \left(T _ {1}, T _ {2}\right)} \frac {\left[ T _ {1} T _ {2} \right] _ {n} \left(z _ {0}\right)}{\left(z - z _ {0}\right) ^ {n}}, \tag {2.3.3}
$$

where  $h(T_1, T_2)$  is some finite number. Because of the correspondence between OPEs and Ward identities, we see that all terms in the sum (2.3.3) have the same conformal dimension:  $\dim([T_1T_2]_n) = \dim(T_1) + \dim(T_2) - n$ . If no negative dimension fields are present in the field theory, we immediately infer that  $h(T_1, T_2)$  is less than or equal to the sum of the conformal dimensions of  $T_1$  and  $T_2$ . As noted before, we restrict ourselves in this work to the case where the correlation functions of the chiral symmetry generators  $T_i(z_i)$  are meromorphic functions in  $z_i$ . In other words, the sum in eq. (2.3.3) runs over the integer numbers. We denote the singular part of an OPE by  $T_1(z)T_2(z_0)$ . It is determined via:

$$
\oint_ {C _ {z _ {0}}} \frac {d z}{2 \pi i} \varepsilon (z) T _ {1} (z) T _ {2} (z _ {0}) = \delta_ {\varepsilon} ^ {T _ {1}} T _ {2} (z _ {0}), \tag {2.3.4}
$$

where we denoted the transformation generated by  $T_{1}$  with parameter  $\varepsilon$  as  $\delta_{\varepsilon}^{T_1}$ . Note that we can assign an OPE only when  $T_{1}$  is a symmetry generator.

The prescription we use for the moment to calculate correlation functions with OPEs, is simply an application of conformal Ward identities like eq. (2.2.13). As

an example, to compute a correlation function with a symmetry generator  $T_{1}$ , we substitute the contractions of  $T_{1}$  with the other fields in the correlator:

$$
\begin{array}{l} <   T _ {1} (z) \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) \dots > = \\ <   T _ {1} (z) \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) \dots > + \\ <   T _ {1} (z) \Phi_ {1} \left(z _ {1}, \bar {z} _ {1}\right) \Phi_ {2} \left(z _ {2}, \bar {z} _ {2}\right) \dots > + \dots \tag {2.3.5} \\ \end{array}
$$

In this way, we can always reduce any correlation function, where symmetry generators are present to a differential polynomial of correlation functions of fields only. Of course, if only generators were present from the start, we will end up with a function of all arguments  $z, w, \ldots$  containing one-point functions and  $< 1 > = 1$ .

We now consider the map from the space of fields to the vectorspace  $\mathcal{V}$  in more detail. The normal ordered product of two fields  $T_{1}$  and  $T_{2}$  is defined by considering a correlation function  $< T_{1}(z)T_{2}(w)X>$ , where  $X$  denotes an arbitrary sequence of operators, and taking the limit of  $z$  going to  $w$  after subtracting all singular terms:

$$
<   : T _ {1} T _ {2}: (w) X > =
$$

$$
\lim  _ {z \rightarrow w} <   \left(T _ {1} (z) T _ {2} (w) - \sum_ {n > 0} \frac {\left[ T _ {1} T _ {2} \right] _ {n} (w)}{(z - w) ^ {n}}\right) X >. \tag {2.3.6}
$$

This definition corresponds to a point-splitting regularisation prescription. In the case of non-meromorphic correlation functions, extra fractional powers of  $(z - w)$  are included in this definition. We will not treat this case here. In the notation of eq. (2.3.3), we see that the map from the space of the fields to  $\mathcal{V}$  is in this case given by:

$$
: T _ {1} T _ {2} \colon \rightarrow [ T _ {1} T _ {2} ] _ {0}. \tag {2.3.7}
$$

In fact, we can use the same procedure to define all operators in the regular part of the OPE<sup>5</sup>, i.e.

$$
\begin{array}{l} <   \sum_ {n \leq 0} \frac {[ T _ {1} T _ {2} ] _ {n} (w)}{(z - w) ^ {n}} X > = \\ <   \left(T _ {1} (z) T _ {2} (w) - \sum_ {n > 0} \frac {\left[ T _ {1} T _ {2} \right] _ {n} (w)}{(z - w) ^ {n}}\right) X >. \tag {2.3.8} \\ \end{array}
$$

The definition (2.3.8) implies that correlation functions can be computed by substituting for two operators their complete OPE. We write:

$$
<   T _ {1} \left(z _ {1}\right) T _ {2} \left(z _ {2}\right) X > = <   \left[ T _ {1} \left(z _ {1}\right) T _ {2} \left(z _ {2}\right) \right] X >. \tag {2.3.9}
$$

In this way, an OPE is now an identity between two bilocal operators.

Because one-point functions are constants and we assumed that correlators have no singularity at infinity, the definition (2.3.8) implies:

$$
<   \left[ T _ {1} (z) T _ {2} (w) \right] _ {n} > = 0, \quad \forall n <   0. \tag {2.3.10}
$$

If all fields in the correlation function eq. (2.3.9) are symmetry generators, we end up with one-point functions  $< T_{i}(z) >$ . These can of course not be computed in the OPE formalism. As infinite sums are involved, one should pay attention to convergence, which will in general only be assured when certain inequalities are obeyed involving the  $z_{i}$ . However, the result should be a meromorphic function with poles in  $(z_{i} - z_{j})$ . So analytic continuation can be used.

The definition (2.3.8) gives the same results as the straightforward application of the Ward identities eq. (2.3.5). This will be proved in the next subsection after determining the consistency requirements on the OPEs.

# 2.3.2 Consistency conditions for OPEs

In this subsection, we determine the consistency requirements on the OPEs by considering (contour integrals of) correlators. We first list the properties of correlation functions used in subsection 2.2.2:

# Assumption 2.3.1 The correlators have the following properties:

- translation and scaling invariance  
- no singularity at infinity  
- the correlation functions involving a chiral symmetry generator  $W(z)$  are meromorphic functions in  $z$  
crossing symmetry

We now suppose that there is a map from the space of fields of the conformal field theory to a vectorspace  $\mathcal{V}$ . A  $\mathbb{Z}_2$  grading in  $\mathcal{V}$  should exist, corresponding to bosonic and fermionic operators. We denote it with  $|A|$ . There is an even linear map  $\partial$  from  $\mathcal{V}$  to  $\mathcal{V}$ , and a sequence of bilinear operations for every  $n\in \mathbb{Z}$  which we denote by  $[AB]_n$ , see eq. (2.3.3). For  $n > 0$  the bilinear operations are defined in eq. (2.3.4). We require that correlation functions can be computed by substituting the complete OPE as done for eq. (2.3.9). From this definition, we determine the properties of

# 2.3. Operator Product Expansions

these maps. In the remainder of this subsection  $X$  and  $Y$  denote arbitrary sequences of operators.

We first determine  $[\partial A B]_n$  by requiring that  $\partial$  corresponds to the derivative on the fields. Because we have:

$$
<   X \partial A (z) Y > = \frac {d}{d z} <   X A (z) Y >, \tag {2.3.11}
$$

we find that the OPE  $\partial A(z)B(w)$  is given by taking the derivative of the OPE  $A(z)B(w)$  with respect to  $z$ :

$$
\partial A (z) B (w) = \sum_ {n <   = h (A, B)} - n \frac {[ A B ] _ {n} (w)}{(z - w) ^ {n + 1}}, \tag {2.3.12}
$$

and hence:

$$
[ \partial A B ] _ {n + 1} = - n [ A B ] _ {n}. \tag {2.3.13}
$$

This equation enables us to write all the terms in the regular part in terms of normal ordered operators:

$$
[ A B ] _ {- n} = \frac {1}{n !}: (\partial^ {n} A) B:, \quad n \in \mathbb {N}. \tag {2.3.14}
$$

By applying derivatives with respect to  $w$  we find:

$$
[ A \partial B ] _ {n + 1} = n [ A B ] _ {n} + \partial [ A B ] _ {n + 1}. \tag {2.3.15}
$$

We now investigate the consequences of the crossing symmetry of the correlators.

First, consider a correlator  $< A(z)B(w)X >$ . As  $X$  is completely arbitrary, we see that the OPE  $A(z)B(w)$  must be equal to the OPE  $(-1)^{|A||B|}B(w)A(z)$ . Looking at eq. (2.3.3) for both OPEs, one sees that a Taylor expansion of  $[AB]_n(z)$  around  $w$  has to be performed to compare the OPEs. The result is:

$$
[ B A ] _ {q} = (- 1) ^ {| A | | B |} \sum_ {l \geq q} \frac {(- 1) ^ {l}}{(l - q) !} \partial^ {(l - q)} [ A B ] _ {l} \tag {2.3.16}
$$

for arbitrary  $q$ .

# Intermezzo 2.3.1

Eq. (2.3.16) leads to a consistency equation when  $A = B$ :

$$
[ A A ] _ {q} = - \sum_ {l > 0} \frac {(- 1) ^ {l}}{2 l !} \partial^ {l} [ A A ] _ {q + l} \quad \text {i f} | A | + q \text {o d d}, \tag {2.3.17}
$$

where we shifted the summation index  $l$  over  $q$  with respect to eq. (2.3.16). If  $A$  is bosonic (fermionic), this relation determines an odd (even) pole of the OPE  $A(z)A(w)$

in terms of derivatives of the higher poles, and thus in terms of the higher even (odd) poles. We can write:

$$
[ A A ] _ {q} = \sum_ {k \geq 0} a _ {k} \partial^ {2 k + 1} [ A A ] _ {q + 2 k + 1},
$$

where the constants  $a_{k}$  satisfy the following recursion relation:

$$
a _ {k} = \frac {1}{2} \left(\frac {1}{(2 k + 1) !} - \sum_ {0 \leq l <   k} \frac {a _ {l}}{(2 k - 2 l) !}\right),
$$

which can be solved by using the generating function:

$$
f (x) = \sum_ {k \geq 0} a _ {k} x ^ {2 k} = \frac {\tanh  (x / 2)}{x},
$$

for which the Taylor expansion is well known. Hence:

$$
a _ {k} = 2 B _ {2 k + 2} \frac {2 ^ {2 k + 2} - 1}{(2 k + 2) !},
$$

where  $B_{n}$  are the Bernoulli numbers. The first values in this series are:

$$
a _ {0} = \frac {1}{2}, a _ {1} = - \frac {1}{2 4}, a _ {2} = \frac {1}{2 4 0}, a _ {3} = - \frac {1 7}{4 0 3 2 0}.
$$

Applying the rule (2.3.17) for  $q = 0$  and  $A$  fermionic, expresses the normal ordered product of a fermionic field with itself in terms of the operators occurring in the singular part of its OPE with itself. An example is found in the  $N = 1$  superconformal algebra. This algebra contains a supersymmetry generator  $G$  which has the following OPE:

$$
G (z) G (w) = \frac {2 c / 3}{(z - w) ^ {3}} + \frac {2 T (w)}{(z - w)} + O (z - w) ^ {0}.
$$

Applying the above formulas, we find:

$$
[ G G ] _ {0} = \partial T.
$$

It is clear that eq. (2.3.16) shows that (dropping sign factors):

$$
\begin{array}{l} <   \left[ [ A (z) B (w) ] C (u) \right] X > _ {\text {O P E}} \\ = <   \left[ B (w) A (z) \right] C (u) ] X > _ {\text {O P E}} \\ = <   [ C (u) [ A (z) B (w) ] ] X > _ {\text {O P E}}. \tag {2.3.18} \\ \end{array}
$$

This does not yet prove that the correlator is crossing symmetric. Indeed, for this we also need (dropping arguments as well):

$$
<   \left[ [ A B ] C \right] X > _ {\mathrm {O P E}} = <   \left[ [ C A ] B \right] X > _ {\mathrm {O P E}}. \tag {2.3.19}
$$

We see that crossing symmetry implies "associativity" of the OPEs: the order in which the OPEs are substituted should be irrelevant. This puts very stringent conditions on the OPEs. These conditions are most easily derived by using contour integrals of correlators. Indeed, we can isolate the contribution of a certain part of the OPE by taking appropriate contour integrals:

$$
\begin{array}{l} <   [ A [ B C ] _ {p} ] _ {q} (u) X > = \oint_ {C _ {u}} \frac {d z}{2 \pi i} (z - u) ^ {q - 1} \\ \oint_ {C _ {u}} \frac {d w}{2 \pi i} (w - u) ^ {p - 1} <   A (z) B (w) C (u) X >, \tag {2.3.20} \\ \end{array}
$$

where  $C_u$  denotes a contour which encircles  $u$  once anti-clockwise, not including any other points involved in the correlator. We can now use a contour deformation argument relating the contour integral in eq. (2.3.20) to a contour integral where the integration over  $w$  is performed last, see fig. 2.1. This integral has two terms: one where the  $z$  contour is around  $u$  (corresponding to the correlator  $(-1)^{|A| |B|} < [B[AC]]X >$ , and one where it is around  $w (< [[AB]C]X>)$ . Using the definition

![](images/a762b8d360dd00f65fee47c1f2e8623b46fbdc8e564b9e38c2d1f75986888795.jpg)  
Figure 2.1: Contour deformations

(2.3.3) for the OPEs, and Cauchy's residue formula for contour integrals, we arrive at:

$$
[ A [ B C ] _ {p} ] _ {q} = (- 1) ^ {| A | | B |} [ B [ A C ] _ {q} ] _ {p} + \sum_ {l > 0} \binom {q - 1} {l - 1} \left[ [ A B ] _ {l} C \right] _ {p + q - l}. \tag {2.3.21}
$$

A second equation follows in the same manner by interchanging the role of  $z$  and  $w$  in fig. 2.1 and bringing the second contour integral of the  $\text{rhs}$  of eq. (2.3.20) to the left:

$$
\begin{array}{l} [ A [ B C ] _ {p} ] _ {q} \\ = (- 1) ^ {| A | | B |} \left(\left[ B [ A C ] _ {q} \right] _ {p} - \sum_ {l > 0} \binom {p - 1} {l - 1} \left[ \left[ B A \right] _ {l} C \right] _ {p + q - l}\right). \tag {2.3.22} \\ \end{array}
$$

Both eqs. (2.3.21),(2.3.22) have to be satisfied for any  $p,q\in \mathbb{Z}$ .<sup>6</sup>

Similarly, by starting with:

$$
\begin{array}{l} <   \left[ \left[ A B \right] _ {p} C \right] _ {q} (u) X > = \oint_ {C _ {u}} \frac {d w}{2 \pi i} (w - u) ^ {q - 1} \\ \oint_ {C _ {w}} \frac {d z}{2 \pi i} (z - w) ^ {p - 1} <   A (z) B (w) C (u) X >, \tag {2.3.23} \\ \end{array}
$$

and using the same contour deformation of fig. 2.1, but now with the first term of the  $\text{rhs}$  brought to the left, we get:

$$
\begin{array}{l} [ [ A B ] _ {p} C ] _ {q} = \sum_ {l \geq q} (- 1) ^ {q - l} \binom {p - 1} {l - q} [ A [ B C ] _ {l} ] _ {p + q - l} - \tag {2.3.24} \\ (- 1) ^ {| A | | B |} \sum_ {l > 0} (- 1) ^ {p - l} \left( \begin{array}{c} p - 1 \\ l - 1 \end{array} \right) [ B [ A C ] _ {l} ] _ {p + q - l}. \\ \end{array}
$$

which has again to be satisfied for all  $p, q \in \mathbb{Z}$ .

Eq. (2.3.16) and eq. (2.3.21) for  $q = 0$  first appeared in [5] where they were derived using the mode algebra associated to OPEs (see section 2.4). In [186], contour deformation arguments were used to find eq. (2.3.21) and an equation related to eq. (2.3.22) for restricted ranges of  $p, q$ . Reference [186] also contains eq. (2.3.24) for  $q > 1$ .

At this point, we found the conditions on the OPEs of the operators in  $\mathcal{V}$  such that we find "correlation functions" which satisfy the assumptions 2.3.1. We still need to show that this procedure indeed gives the correlation functions of the conformal field theory we started with, i.e. that we find correlation functions which satisfy the Ward identities of the theory. Before proceeding, we prove the following lemma:

Lemma 2.3.1 The associativity condition eq. (2.3.21) for strictly positive  $q$  can be rewritten as:

$$
A (z) [ B (w) C (u) ] = (- 1) ^ {| A | | B |} \left[ B (w) A (z) C (u) \right] + \left[ A (z) B (w) C (u) \right]
$$

$$
f o r \left| w - u \right| <   \left| z - u \right|.
$$

# Proof :

The proof of this lemma is in some sense the reverse of the derivation of eq. (2.3.21). We multiply eq. (2.3.21) with  $(z - u)^q (w - u)^p$  and sum for  $p$  and  $q$  over the appropriate ranges. After subtracting the result from the proposition of this lemma, we need to

# 2.3. Operator Product Expansions

prove:

$$
\begin{array}{l} \left[ \begin{array}{c} A (z) B (w) C (u) \end{array} \right] = \\ \sum_ {q > 0} \sum_ {p} (z - u) ^ {- q} (w - u) ^ {- p} \sum_ {l > 0} \left( \begin{array}{c} q - 1 \\ l - 1 \end{array} \right) [   [ A B ] _ {l} C ] _ {p + q - l} (u). \\ \end{array}
$$

Due to the binomial coefficient, and because  $q$  is strictly positive, the sum over  $l$  is from one to  $q$ . Call  $r = p + q - l$  and  $s = q - l$ . We find:

$$
\begin{array}{l} r h s = \sum_ {r} \sum_ {l > 0} \left[ \left[ A B \right] _ {l} C \right] _ {r} (u) \\ \left(\sum_ {s \geq 0} \binom {l + s - 1} {l - 1} (z - u) ^ {- l - s} (w - u) ^ {- r + s}\right) \\ = \sum_ {r} \sum_ {l > 0} (z - u) ^ {- l} (w - u) ^ {- r} \left(1 - \frac {w - u}{z - u}\right) ^ {- l} \left[ [ A B ] _ {l} C ] _ {r} (u), \right. \\ \end{array}
$$

which is exactly the lhs. In the last step we used the Taylor expansion of  $(1 - x)^{-l}$  which only converges for  $|x| < 1$ .

The Ward identities have the form (2.2.13):

$$
\begin{array}{l} \oint_ {C} \frac {d z}{2 \pi i} \varepsilon (z) <   \left[ T _ {0} (z) T _ {1} \left(z _ {1}\right) \right] \dots T _ {n} \left(z _ {n}\right) > = \\ \oint_ {C _ {1}} \frac {d z}{2 \pi i} \varepsilon (z) <   T _ {0} (z) T _ {1} (z _ {1}) \dots T _ {n} (z _ {n}) > + \dots \\ + \oint_ {C _ {n}} \frac {d z}{2 \pi i} \varepsilon (z) <   T _ {0} (z) T _ {1} \left(z _ {1}\right) \dots T _ {n} \left(z _ {n}\right) >, \tag {2.3.25} \\ \end{array}
$$

where the contractions correspond by definition to the transformation generated by  $T_{0}$  of the field  $T_{i}$ , eq. (2.3.4). The contour  $C$  surrounds  $z_{1}\dots z_{n}$  once in the anticlockwise direction, while the contours  $C_i$  encircle only  $z_{i}$ .  $\varepsilon (z)$  is analytic in a region containing  $C$ .

Theorem 2.3.1 The Ward identities eq. (2.3.25) are satisfied for correlation functions which we compute by substituting two operators with their complete OPE.

# Proof :

The  $(n + 1)$ -point function  $< T_0(z)T_1(z_1)\dots T_n(z_n)>$  has only poles in  $z$  when  $z = z_i$ . We will use this to deform the contour  $C$ .

Suppose that the theorem holds for  $n$ -point functions. We rearrange the operators  $T_{i}$

such that  $|z_{1} - z_{2}|$  is smaller than all other  $|z_{i} - z_{j}|$ . We substitute  $T_{1}(z_{1})T_{2}(z_{2})$  by their OPE. Then we have to compute a  $n$ -point function.

We now split the contour  $C$  in a contour  $C_2^\prime$  and the contours  $C_i$ ,  $i > 2$ , used in eq. (2.3.25).  $C_2^\prime$  encircles  $z_{1}$  and  $z_{2}$ , but no other  $z_{i}$ . Because of the rearrangement, we can take  $C_2^\prime$  such that for any point  $z$  on  $C_2^\prime$ ,  $|z - z_2| > |z_1 - z_2|$ . By applying the recursion assumption we find:

$$
\begin{array}{l} \oint_ {C} \frac {d z}{2 \pi i} \varepsilon (z) <   T _ {0} (z) T _ {1} (z _ {1}) \dots T _ {n} (z _ {n}) > = \\ \oint_ {C _ {2} ^ {\prime}} \frac {d z}{2 \pi i} \varepsilon (z) <   T _ {0} (z) [ T _ {1} (z _ {1}) T _ {2} (z _ {2}) ] \dots T _ {n} (z _ {n}) > \\ + \oint_ {C _ {3}} \frac {d z}{2 \pi i} \varepsilon (z) <   T _ {0} (z) [ T _ {1} (z _ {1}) T _ {2} (z _ {2}) ] T _ {3} (z _ {3}) \dots T _ {n} (z _ {n}) > + \dots \\ \end{array}
$$

We can now use lemma 2.3.1 for the first term of the  $\text{rhs}$ . We indeed find the  $\text{rhs}$  of eq. (2.3.25).

We comment on when a correlation function can be computed by taking all contractions, which was the prescription we used to define OPEs in the previous subsection, see eq. (2.3.5). Computing a correlation function in this way, means that we drop the integrals in the Ward identity eq. (2.3.25). This can only be done if all one-point functions  $< [T_1 T_2]_0 >$  vanish (except if  $T_i = \mathbf{1}$ ). Because of scaling invariance, this is true when all operators (except  $\mathbf{1}$ ) have strictly positive dimension. Bowcock argued in [34] that correlation functions can be computed by substituting the complete OPE, or by using contractions. His argument is based on the claim that eq. (2.3.25) is true, but no proof is given.

To conclude this subsection, we wish to mention that Wilson Operator Product Expansions were already used outside the scope of two-dimensional conformal field theory [208, 191]. However, because the consistency requirements on the OPEs are especially strong in two-dimensional conformal field theory (due to the fact that the conformal algebra has infinite dimension), it is there where the full power of the formalism comes to fruition.

# 2.3.3 Operator Product Algebras

We can assemble the consistency conditions on OPEs in a definition, see also [28, 96].

Definition 2.3.1 An Operator Product Algebra (OPA) is a  $\mathbb{Z}_2$  graded vectorspace  $\mathcal{V}$  with elements  $1, A, B, C \cdots$ , an even linear map  $\partial$ , and a bilinear binary operation:

$$
[.. ] _ {l}: \mathcal {V} \otimes \mathcal {V} \rightarrow \mathcal {V}, \quad l \in \mathbb {Z},
$$

which is zero for  $l$  sufficiently large. The following properties hold:

- unity:

$$
[ \mathbf {1} A ] _ {l} = \delta_ {l} A
$$

- commutation : (eq. (2.3.16))

$$
[ B A ] _ {n} = (- 1) ^ {| A | | B |} \sum_ {l \geq n} \frac {(- 1) ^ {l}}{(l - n) !} \partial^ {(l - n)} [ A B ] _ {l} \quad \forall n \in \mathbb {Z}
$$

- associativity : (eq. (2.3.21))

$$
\begin{array}{l} [ A [ B C ] _ {p} ] _ {q} = (- 1) ^ {| A | | B |} [ B [ A C ] _ {q} ] _ {p} \\ + \sum_ {l > 0} \binom {q - 1} {l - 1} \left[ [ A B ] _ {l} C ] _ {p + q - l} \right. \qquad \forall p, q \in \mathbb {Z} \\ \end{array}
$$

The properties in this definition are sufficient to recover all consistency conditions of the previous subsection. Using  $\partial A = [A\mathbf{1}]_1$  and the associativity condition, the equations for OPEs of derivatives (2.3.13,2.3.15) follow. Similarly, the other associativity conditions (2.3.22) and (2.3.24) can be derived from the properties of an OPA.

An alternative definition would be to impose eqs. (2.3.13, 2.3.15), requiring the associativity condition eq. (2.3.21) only for  $p, q \geq 0$ . One then usually considers eq. (2.3.21) for  $p$  or  $q$  equal to zero as being the definition of how we can calculate OPEs of normal ordered products. The set of consistency conditions eq. (2.3.21) for  $p, q > 0$  has then to be checked for all generators of the OPA.

In some cases, the definition 2.3.1 is too strong. OPEs are intended to compute correlation functions. It is possible that there are some fields in the theory which have vanishing correlators with all other fields. These null fields should be taken into account in the definition of an OPA. Indeed, they could occur in every consistency equation for OPEs without affecting the results for the correlation functions. We need an algebraic definition of a null field. From the way we compute correlators using OPEs, we see that if  $N$  is a null field,  $[NA]_n$  should be again a null field, for any  $A \in \mathcal{V}$  and  $n \in \mathbb{Z}$ . If this would not be true, we could write down a nonvanishing correlator with  $N$ . We see that the null fields form an ideal in the OPA.

Definition 2.3.2 Consider an OPA as defined in def. 2.3.1. Suppose there is an ideal  $\mathcal{N}$  in the OPA, whose elements we call null operators (or null fields). We extend the definition of an OPA to algebras where the defining properties are only satisfied up to elements of  $\mathcal{N}$ .

This extension of the definition of an OPA was not considered in [28, 96].

It is in general difficult to check if we can consider an operator  $N$  to be null. We do not want that  $\mathbf{1} \in \mathcal{N}$ , because then all operators are null. Hence, a necessary

condition for  $N$  to be null is that we can find no operator  $A$  in  $\mathcal{V}$  such that  $[NA]_n \sim \mathbb{1}$  (for some  $n$ ). Usually, this is regarded as a sufficient condition, because one generally works with fields of strictly positive conformal dimension. Due to the scaling invariance of the correlation functions, all one-point functions are then zero, except  $< \mathbb{1}>$ . So, in this case, any field which does not produce the identity operator in some OPE has zero correlation functions.

# 2.3.4 OPA-terminology

In this subsection, we introduce some terminology which is continually used in the rest of this work.

Definition 2.3.3 We call  $A$  a composite operator if it is equal to  $[BC]_0$ , for some  $B, C \neq \mathbb{1}$ , except when  $B$  is equal to  $C$  and fermionic.

The condition when  $B = C$  in this definition comes from the considerations in intermezzo 2.3.1.

Definition 2.3.4 A set of operators is said to generate the OPA if all elements of  $\mathcal{V}$  can be constructed by using addition, scalar multiplication, derivation and taking composites.

Definition 2.3.5 A Virasoro operator  $T$  has the following non-zero poles in its OPE (2.3.2):

$$
[ T T ] _ {4} = \frac {c}{2} \mathbb {1}, \qquad [ T T ] _ {2} = 2 T, \qquad [ T T ] _ {1} = \partial T.
$$

Definition 2.3.6 An operator  $A$  is a scaling operator with respect to  $T$ , with (conformal) dimension  $h$  if:

$$
[ T A ] _ {2} = h A, \qquad [ T A ] _ {1} = \partial A.
$$

A quasiprimary operator is a scaling operator with in addition  $[TA]_3 = 0$ .

A primary operator is a scaling operator with in addition  $[TA]_n = 0$  for  $n \geq 3$ .

A conformal OPA is an OPA with a Virasoro  $T$ . All other operators of the OPA (except 1) are required to be scaling operators with respect to  $T$ .

# Intermezzo 2.3.2

As an application of the above definitions, we wish to show that when  $A$  and  $B$ , with dimension  $a$  and  $b$ , are scaling operators with respect to  $T$  then  $[AB]_n$  is a scaling

# 2.3. Operator Product Expansions

operator with dimension  $a + b - n$ . Let us compute the first order pole of this operator with  $T$  using eq. (2.3.21):

$$
\begin{array}{l} [ T [ A B ] _ {n} ] _ {1} = [ A [ T B ] _ {1} ] _ {n} + \left[ \left(T A \right] _ {1} B \right] _ {n} \\ = [ A \partial B ] _ {n} + [ \partial A B ] _ {n} \\ = \partial [ A B ] _ {n}, \\ \end{array}
$$

where we used the sum of eq. (2.3.13) and eq. (2.3.15) in the last step. For the second order pole we have:

$$
\begin{array}{l} [ T [ A B ] _ {n} ] _ {2} = [ A [ T B ] _ {2} ] _ {n} + \left[ [ T A ] _ {2} B \right] _ {n} + \left[ [ T A ] _ {1} B \right] _ {n + 1} \\ = (a + b) [ A B ] _ {n} + [ \partial A B ] _ {n + 1} \\ = (a + b - n) [ A B ] _ {n}. \\ \end{array}
$$

Definition 2.3.7 A map from the OPA to the half-integer numbers is called a dimension if it has the properties:

$$
\begin{array}{l} \dim (\mathbb {1}) = 0 \\ \dim (\partial A) = \dim (A) + 1 \\ \dim ([ A B ] _ {l}) = \dim (A) + \dim (B) - l. \\ \end{array}
$$

If such a map exists, we call the OPA graded.

Definition 2.3.8 A  $\mathcal{W}$ -algebra is a conformal OPA where one can find a set of generators which are quasiprimary.

Different definitions of a  $\mathcal{W}$ -algebra exist in the literature. Sometimes one requires that the generators are primary (except  $T$  itself). In this work, we will mainly consider  $\mathcal{W}$ -algebras of this subclass, and for which the number of generators is finite. The importance of  $\mathcal{W}$ -algebras lies in the fact that the chiral symmetry generators of a conformal field theory form a  $\mathcal{W}$ -algebra<sup>7</sup>. We will treat  $\mathcal{W}$ -algebras in more detail in chapter 4.

Finally, we introduce a notation for OPEs which lists only the operators in the singular terms, starting with the highest order pole. As an example, we will write a Virasoro OPE (2.3.2) as:

$$
T \times T = \ll \frac {c}{2} | 0 | 2 T | \partial T \gg . \tag {2.3.26}
$$

# 2.3.5 Poisson brackets

To conclude this section, we want to show the similarity between Poisson bracket calculations and OPEs.

In a light-cone quantisation scheme (choosing  $\bar{z}$  as "time") [209], the symmetry generators in classical conformal field theory obey Poisson brackets of the form:

$$
\left\{A (z), B \left(z _ {0}\right) \right\} _ {\mathrm {P B}} = \sum_ {n > 0} \frac {(- 1) ^ {n - 1}}{(n - 1) !} \left\{A B \right\} _ {n} \left(z _ {0}\right) \partial^ {n - 1} \delta \left(z - z _ {0}\right), \tag {2.3.27}
$$

where  $\{AB\}_{n}$  are also symmetry generators of the theory. The derivative is with respect to the  $z$ -coordinate. We choose the normalisation factors such that:

$$
\left\{A B \right\} _ {n} (z) = \int d z (z - w) ^ {n - 1} \left\{A (z), B (w) \right\} _ {\mathrm {P B}}. \tag {2.3.28}
$$

For convenience, we drop the subscript PB in the rest of this subsection. The Poisson brackets satisfy:

$$
\{A (z), B \left(z _ {0}\right) \} = (- 1) ^ {| A | | B |} \left\{B \left(z _ {0}\right), A (z) \right\}
$$

$$
\{\partial A (z), B (z _ {0}) \} = \frac {d}{d z} \{A (z), B (z _ {0}) \}
$$

$$
\begin{array}{l} \{\{A (z _ {1}), B (z _ {2}) \}, C (z _ {3}) \} = \{A (z _ {1}), \{B (z _ {2}), C (z _ {3}) \} \} - \\ (- 1) ^ {| A | | B |} \left\{B \left(z _ {2}\right), \left\{A \left(z _ {1}\right), C \left(z _ {3}\right) \right\} \right\} \tag {2.3.29} \\ \end{array}
$$

These relations imply identities for the  $\{AB\}_n$ . We do not list the consequences of the first two, as they are exactly the same as eqs. (2.3.16) and (3.3.1). The Jacobi identities give:

$$
\begin{array}{l} \{A \left\{B C \right\} _ {p} \} _ {q} (z _ {3}) - (- 1) ^ {| A | | B |} \left\{B \left\{A C \right\} _ {q} \right\} _ {p} (z _ {3}) \\ = \int d z _ {2} (z _ {2} - z _ {3}) ^ {p - 1} \int d z _ {1} \sum_ {l = 1} ^ {q} \binom {q - 1} {l - 1} \\ \left. \left(z _ {1} - z _ {2}\right) ^ {l - 1} \left(z _ {2} - z _ {3}\right) ^ {q - l} \{\{A \left(z _ {1}\right), B \left(z _ {2}\right) \}, C \left(z _ {3}\right) \} \right. \\ = \sum_ {l = 1} ^ {q} \binom {q - 1} {l - 1} \left\{\left\{A B \right\} _ {l} C \right\} _ {p + q - l} \left(z _ {3}\right) \}. \tag {2.3.30} \\ \end{array}
$$

These equations are of exactly the same form as the associativity conditions for OPEs (see (2.3.21) with  $q, p > 0$ ).

An important difference with OPEs is that no "regular" part is defined for Poisson brackets. In particular, normal ordering is not necessary. A Poisson bracket where a product of fields is involved, is simply:

$$
\{A (z), B (w) C (w) \} = \{A (z), B (w) \} C (w) + (- 1) ^ {| A | | B |} B (w) \{A (z), C (w) \}. \tag {2.3.31}
$$

When using the notation:

$$
\left\{A B \right\} _ {- n} (z) = \frac {1}{n !} \left(\partial^ {n} A (z)\right) B (z), \quad n \geq 0, \tag {2.3.32}
$$

we see that eq. (2.3.31) corresponds to eq. (2.3.21) for  $p = 0$  and with the double contractions  $l < q$  dropped. This is also true for a classical version of eq. (2.3.24). We can conclude that when using the correspondence:

$$
\left\{A B \right\} _ {n} \leftrightarrow [ A B ] _ {n}, \tag {2.3.33}
$$

computing with Poisson brackets follows almost the same rules as used for OPEs: one should drop double contractions and use a graded-commutative and associative normal ordering. In particular, as the Jacobi-identities are the same, any linear PB-algebra corresponds to an operator product algebra and vice-versa. For nonlinear algebras, this is no longer true because of normal ordering.

Due to this correspondence, we will often write "classical OPEs" for Poisson brackets.

# 2.4 Mode algebra

In this section, we show that there is an infinite dimensional algebra with a graded-symmetric bracket associated to every OPA. For every operator  $A$  we define the  $m$ -th mode of  $A$  by specifying how it acts on an operator:

$$
\widehat {A} _ {m} B \equiv [ A B ] _ {m + a}, \tag {2.4.1}
$$

where  $a$  is the conformal dimension of  $A$ . The shift in the index is made such that  $\widehat{A}_mB$  has dimension  $b - m$ , independent of  $A^8$ . Hence for operators with (half-)integer dimension,  $m$  is (half-)integer. An immediate consequence of this definition follows by considering (2.3.13):

$$
\widehat {(\partial A)} _ {m} = - (m + a) \widehat {A} _ {m}. \tag {2.4.2}
$$

We can now compute the graded commutator of two modes:

$$
[ \widehat {A} _ {m}, \widehat {B} _ {n} ] C = [ A [ B C ] _ {n + b} ] _ {a + m} - (- 1) ^ {| A | | B |} [ B [ A C ] _ {m + a} ] _ {n + b}. \tag {2.4.3}
$$

Using the associativity condition (2.3.21) we see that:

$$
\left[ \widehat {A} _ {m}, \widehat {B} _ {n} \right] = \sum_ {l > 0} \binom {m + a - 1} {l - 1} \left(\left[ \widehat {A B} \right] _ {l}\right) _ {m + n}. \tag {2.4.4}
$$

Note that this commutator is determined by the singular part of the OPE.

Theorem 2.4.1 Eq. (2.4.4) defines a graded commutator.

# Proof :

We can use the relation between the OPEs  $AB$  and  $BA$  (2.3.16) and eq. (2.4.2) to show that:

$$
\begin{array}{l} - (- 1) ^ {| A | | B |} [ \widehat {B} _ {n}, \widehat {A} _ {m} ] \\ = \sum_ {l > 0} \binom {n + b - 1} {l - 1} \sum_ {p > l} \frac {(- 1) ^ {p}}{(p - l) !} (\widehat {\partial^ {p - l} [ A B ] _ {p}}) _ {m + n} \\ = \sum_ {l > 0} \sum_ {p > l} \binom {n + b - 1} {l - 1} \binom {m + n + a + b - l - 1} {p - l} (\widehat {[ A B ] _ {p}}) _ {m + n}. \\ \end{array}
$$

Using eq. (2.A.4), we see that this is equal to (2.4.4).

A first example of a mode algebra is given by the modes of a Virasoro operator  $T$  with OPE (2.3.2). For historic reasons, we denote these modes with  $L$ . We find using eq. (2.4.4):

$$
\left[ \widehat {L} _ {m}, \widehat {L} _ {n} \right] = (m - n) \widehat {L} _ {m + n} + \frac {c}{2} \binom {m + 1} {3} \delta_ {m + n}, \tag {2.4.5}
$$

where we used that the modes of the unit operator (which is implicit in the fourth order pole of (2.3.2)) are given by:

$$
\widehat {\mathbf {1}} _ {m} = \delta_ {m}. \tag {2.4.6}
$$

The infinite dimensional Lie algebra with commutator (2.4.5) is called the Virasoro algebra. We see that this algebra is a central extension of the algebra of the classical generators of conformal transformations (2.1.8). The modes corresponding to the global conformal transformations  $\widehat{L}_1, \widehat{L}_0, \widehat{L}_{-1}$  form a finite dimensional subalgebra where the central extension drops out.

Similarly, the OPE of  $T$  with a primary field  $\Phi$  with dimension  $h$  (2.3.1) gives the following commutator:

$$
\left[ \hat {L} _ {m}, \hat {\Phi} _ {n} \right] = ((h - 1) m - n) \hat {\Phi} _ {m + n}. \tag {2.4.7}
$$

We have the following important theorem.

# 2.5. Generating functionals

Theorem 2.4.2 For a given OPA, the commutator of the corresponding mode algebra satisfies graded Jacobi identities modulo modes of null fields:

$$
[ \widehat {A} _ {k}, [ \widehat {B} _ {l}, \widehat {C} _ {m} ] ] = (- 1) ^ {| A | | B |} [ \widehat {B} _ {l}, [ \widehat {A} _ {k}, \widehat {C} _ {m} ] ] + [ [ \widehat {A} _ {k}, \widehat {B} _ {l} ], \widehat {C} _ {m} ] ]. \tag {2.4.8}
$$

Proof :

The  $lhs$  of eq. (2.4.8) is by definition (2.4.4) equal to:

$$
[ \widehat {A} _ {k}, [ \widehat {B} _ {l}, \widehat {C} _ {m} ] ] =
$$

$$
\sum_ {q > 0} \sum_ {p > 0} \binom {k + a - 1} {q - 1} \binom {l + b - 1} {p - 1} \left(\left[ A [ \widehat {B C} ] _ {p} \right] _ {q}\right) _ {k + l + m}. \tag {2.4.9}
$$

We can now use the associativity condition (2.3.21). Calling the summation index in (2.3.21)  $r$  and renaming  $q = s + r - p$ , it follows that the Jacobi identity (2.4.8) will be satisfied if:

$$
\sum_ {p > 0} \binom {k + a - 1} {r + s - p - 1} \binom {l + b - 1} {p - 1} \binom {r + s - p - 1} {r - 1} =
$$

$$
\binom {k + a - 1} {r - 1} \binom {k + l + a + b - r - 1} {s - 1}. \tag {2.4.10}
$$

After cancelling out factors, one sees that this equation follows from eq. (2.A.5).

Moreover, from the above proof it is clear that the reverse is also true:

Theorem 2.4.3 If the mode algebra satisfies the Jacobi-identities (2.4.8) up to null fields, the associativity conditions (2.3.21) are satisfied.

Modes of normal ordered operators are given by eq. (2.3.24):

$$
\left(\left[ \widehat {A B} \right] _ {0}\right) _ {m} = \sum_ {l}: \widehat {A} _ {l} \widehat {B} _ {m - l}:, \tag {2.4.11}
$$

where

$$
: \widehat {A} _ {l} \widehat {B} _ {m}: \equiv \left\{ \begin{array}{l l} \widehat {A} _ {l} \widehat {B} _ {m} & \text {i f} l \leq - a \\ (- 1) ^ {| A | | B |} \widehat {B} _ {m} \widehat {A} _ {l} & \text {i f} l > - a \end{array} \right. \tag {2.4.12}
$$

Consider an OPA where the generators have OPEs whose singular part contains composite operators. In this case, the mode algebra is only an infinite dimensional (super-)Lie algebra when those composites are viewed as new elements of the algebra, such that the commutators close linearly. Otherwise, the commutators close only in the enveloping algebra of the modes of the generators.

As the mode algebra contains the same information as the OPEs, one can always choose which one uses in a certain computation. For linear algebras (where the OPEs

close on a finite number of noncomposite operators) modes are very convenient. However, for nonlinear algebras the infinite sums in the modes of a composite are more difficult to handle.

The definition eq. (2.4.1) provides a realisation of the mode algebra. In a canonical quantisation scheme, another representation is found in terms of the creation and annihilation operators. For different periodicity conditions in the coordinate  $\sigma$  (relating  $z$  and  $\exp(2i\pi)z$ ) of the symmetry generators, formally the same algebra arises. However, the range of the indices differs. As an example, the  $N = 1$  superconformal algebra consists of a Virasoro operator  $T$  and a fermionic dimension  $3/2$  primary operator  $G$  with OPE:

$$
G \times G = \ll \frac {2 c}{3} \mathbb {1} | 0 | 2 T \gg . \tag {2.4.13}
$$

This gives for the anticommutator of the modes:

$$
\left[ \widehat {G} _ {m}, \widehat {G} _ {n} \right] = \frac {2 c}{3} \binom {m + 1 / 2} {2} \delta_ {m + n} + 2 \widehat {L} _ {m + n}. \tag {2.4.14}
$$

In the representation eq. (2.4.1) defined via the OPEs,  $m$  and  $n$  in this commutator are half-integer numbers. However, the algebra is also well-defined if  $m$  and  $n$  are integer. This corresponds to different boundary conditions on  $G(z)$ . The relation between the different modings of the linear superconformal algebras is studied in [178, 50]. We will use the notation  $\widehat{A}_m$  in the representation eq. (2.4.1) of the mode algebra, and drop the hats otherwise.

# 2.5 Generating functionals

In this section, we will define the generating functionals for the correlation functions of a conformal field theory and show that the Ward identities give a set of functional equations for these functionals.

Consider a conformal field theory with fields  $\phi_{i}$  and action  $S[\phi_i]$ . We denote the generators of the chiral symmetries of the theory with  $T_{k}$ . The partition function  $Z$  is defined by:

$$
\begin{array}{l} Z [ \mu ] = \frac {1}{\mathcal {N}} \int \left[ d \varphi_ {i} \right] \exp \left(- S \left[ \varphi_ {j} \right] - \frac {1}{\pi} \int d ^ {2} x \mu^ {k} (x) T _ {k} (x)\right) (2.5.1) \\ = <   \exp \left(- \frac {1}{\pi} \int d ^ {2} x \mu^ {k} (x) T _ {k} (x)\right) > _ {\mathrm {O P E}}. (2.5.2) \\ \end{array}
$$

Here the normalisation constant  $\mathcal{N}$  was defined in eq. (2.2.2) and  $\mu^k$  (the "sources") are non-fluctuating fields.  $Z$  is the generating functional for the correlation functions

of the generators  $T_{k}$ . Indeed, by functional derivation with respect to the sources we can determine every correlation function:

$$
<   T _ {1} (x _ {1}) T _ {2} (x _ {2}) T _ {3} (x _ {3}) > = (- \pi) ^ {3} \frac {\delta}{\delta \mu^ {1} (x _ {1})} \frac {\delta}{\delta \mu^ {2} (x _ {2})} \frac {\delta}{\delta \mu^ {3} (x _ {3})} Z [ \mu ] \Big | _ {\mu = 0}. \qquad (2. 5. 3)
$$

If a generator  $T_{k}$  is fermionic,  $\mu^{k}$  is Grasmann-odd, i.e. anticommuting. We will always use left-functional derivatives in this work.

It is often useful to define the induced action  $\Gamma$  as the generating functional of all "connected" diagrams:

$$
Z [ \mu ] = \exp (- \Gamma [ \mu ]). \tag {2.5.4}
$$

We can view the sources  $\mu^k$  as gauge fields. By assigning appropriate transformation rules for  $\mu^k$ , we can try to make the global symmetries local at the classical level. As an example, in a conformal field theory, the action is invariant under the conformal transformations generated by  $T$ , when the parameter  $\varepsilon^z$  is analytic. If we want to make the theory invariant for any  $\varepsilon^z$ , we have to couple the generator  $T$  to a gauge field. We find that:

$$
\delta_ {\varepsilon} (S + \frac {1}{\pi} \int \mu T) = \frac {1}{\pi} \int (- \bar {\partial} \varepsilon T + \delta \mu T + \mu \delta T), \tag {2.5.5}
$$

where we used the definition of  $T$  as a Noether current (eq. (2.2.7)). Applying the transformation rule for  $T$  (eq. (2.2.15)) in the classical case ( $c = 0$ ), we see that if the source  $\mu$  transforms as:

$$
\delta_ {\varepsilon} \mu = \bar {\partial} \varepsilon + \varepsilon \partial \mu - \partial \varepsilon \mu , \tag {2.5.6}
$$

the action  $S + \frac{1}{\pi} \int \mu T$  is classically invariant. When gauging not only the conformal symmetry, we would expect that higher order terms in the sources have to be added to eq. (2.5.1) to obtain invariance. However, Hull [115] proved that minimal coupling (addition of only linear terms) is sufficient to gauge chiral symmetries.

It is possible that the resulting local symmetry does not survive at the quantum level. For example, the Schwinger term in the transformation law of  $T$  (eq. (2.2.15)) breaks gauge invariance. In general, central terms in the transformation laws of the symmetry generators give rise to "universal" anomalies, which cannot be canceled by changing the transformation law of the gauge fields<sup>9</sup>. In this case, the induced action  $\Gamma$  is a (in general nonlocal) functional of the gauge fields  $\mu^k$  and is used in  $\mathcal{W}$ -gravity theories (see chapter 7) and non-critical  $\mathcal{W}$ -strings.

The Ward identities can be used to derive functional equations for the generating functionals  $Z$  and  $\Gamma$ . As an example, consider the induced action where only the energy-momentum tensor is coupled to a source (i.e. the generating functional for

correlation functions with only  $T$ ). Under the variation of the source (eq. (2.5.6)), the partition function transforms as:

$$
\begin{array}{l} \delta Z [ \mu ] = <   \exp \left(- \frac {1}{\pi} \int (\mu + \delta \mu) T\right) > - Z [ \mu ] \\ = - \frac {1}{\pi} \int (\bar {\partial} \varepsilon + \varepsilon \partial \mu - (\partial \varepsilon) \mu) <   T \exp \left(- \frac {1}{\pi} \int \mu T\right) > \tag {2.5.7} \\ \end{array}
$$

To compute the first integral with  $\bar{\partial}\varepsilon T$ , we note that in the complex basis, when  $\varepsilon^{\bar{z}} = 0$  and all  $\Phi_i = T$ , eq. (2.2.9) becomes:

$$
\begin{array}{l} - \frac {1}{\pi} \int d ^ {2} x \bar {\partial} \varepsilon (x) <   T (x) T (x _ {1}) \dots T (x _ {N}) > \\ = \sum_ {j = 1} ^ {N} <   T \left(x _ {1}\right) \dots \delta_ {\varepsilon} T \left(x _ {j}\right) \dots T \left(x _ {N}\right) >. \tag {2.5.8} \\ \end{array}
$$

We can now multiply this equation with  $\mu(x_1) \cdots \mu(x_N)$  and integrate over all  $x_j$ . Using crossing symmetry in the rhs of eq. (2.5.8), we get (in an obvious notation):

$$
\begin{array}{l} \int d ^ {2} x \bar {\partial} \varepsilon (x) <   T (x) \left(- \frac {1}{\pi} \int \mu T\right) ^ {N} > \\ = N \int d ^ {2} x \mu (x) <   \delta_ {\varepsilon} T (x) \left(- \frac {1}{\pi} \int \mu T\right) ^ {N - 1} >. \tag {2.5.9} \\ \end{array}
$$

This gives the following result:

$$
\begin{array}{l} \int d ^ {2} x \bar {\partial} \varepsilon (x) <   T (x) \exp \left(- \frac {1}{\pi} \int \mu T\right) > = \\ \int d ^ {2} x \mu (x) <   \delta_ {\varepsilon} T (x) \exp \left(- \frac {1}{\pi} \int \mu T\right) >. \tag {2.5.10} \\ \end{array}
$$

Using eq. (2.2.15), the variation of the partition function (2.5.7) becomes:

$$
\delta Z [ \mu ] = \frac {c}{1 2 \pi} \left(\int \varepsilon \partial^ {3} \mu\right) Z [ \mu ]. \tag {2.5.11}
$$

On the other hand, we can rewrite the  $r h s$  of eq. (2.5.7) using:

$$
<   T (x) \exp \left(- \frac {1}{\pi} \int \mu T\right) > = - \pi \frac {\delta Z}{\delta \mu (x)} [ \mu ] = \pi Z [ \mu ] \frac {\delta \Gamma}{\delta \mu (x)} [ \mu ]. \tag {2.5.12}
$$

# 2.6. A few examples

Combining (2.5.7) and (2.5.12) gives us a functional equation for  $Z$ , and thus for  $\Gamma$ :

$$
\frac {c}{1 2 \pi} \partial^ {3} \mu = \bar {\partial} t - \mu \partial t - 2 (\partial \mu) t \tag {2.5.13}
$$

where:

$$
t = \frac {\delta \Gamma}{\delta \mu}. \tag {2.5.14}
$$

We see that the Ward identities fix the generating functionals. This is of course no surprise, as we already knew that the Ward identities determine the correlation functions. The underlying theory  $S[\varphi_j]$  determines the central charge  $c$ .

By rescaling  $t$  in eq. (2.5.13), we can make the functional equation  $c$ -independent. We see that  $\Gamma = c\Gamma^{(0)}$ , where  $\Gamma^{(0)}$  is  $c$ -independent. In fact, the solution of the Ward identity (2.5.13) was given by Polyakov [159]:

$$
\Gamma [ \mu ] = \frac {c}{2 4 \pi} \int \mu \partial^ {2} (1 - \bar {\partial} ^ {- 1} \mu \partial) ^ {- 1} \frac {\partial}{\bar {\partial}} \mu , \tag {2.5.15}
$$

where the inverse derivative is defined in eq. (A.7). For a general set of symmetry generators, it will not be possible to solve the functional equation in closed form.

We now show how OPEs can be used to derive the functional equation on the partition function  $Z$  (eq. (2.5.1)). We start by computing:

$$
- \pi \bar {\partial} \frac {\delta Z}{\delta \mu^ {i} (z , \bar {z})} = \bar {\partial} <   T _ {i} (z) \exp \left(- \frac {1}{\pi} \int d ^ {2} x \mu^ {k} (x) T _ {k} (x)\right) >. \tag {2.5.16}
$$

Let us look at order  $N$  in the sources:

$$
\begin{array}{l} \bar {\partial} <   T _ {i} (z) \left(- \int \mu^ {k} T _ {k}\right) ^ {N} > \\ = - N \int d ^ {2} x _ {0} \mu^ {j} (x _ {0}) \bar {\partial} <   \sum_ {n > 0} \frac {[ T _ {i} T _ {j} ] _ {n} (z _ {0})}{(z - z _ {0}) ^ {n}} \left(- \int \mu^ {k} T _ {k}\right) ^ {N - 1} > \\ = \pi N \sum_ {n > 0} \frac {(- 1) ^ {n}}{(n - 1) !} \\ \partial^ {n - 1} \left(\mu^ {j} (z, \bar {z}) <   [ T _ {i} T _ {j} ] _ {n} (z) \left(- \int \mu^ {k} T _ {k}\right) ^ {N - 1} >\right) \\ = \pi N \sum_ {n > 0} \frac {\partial^ {n - 1} \mu^ {j} (z , \bar {z})}{(n - 1) !} <   \left[ T _ {j} T _ {i} \right] _ {n} (z) \left(- \int \mu^ {k} T _ {k}\right) ^ {N - 1} >, \tag {2.5.17} \\ \end{array}
$$

where we used eq. (A.3) in the second step, and the last step [42] relies on eq. (2.3.16). Hence,

$$
- \pi \bar {\partial} \frac {\delta Z}{\delta \mu^ {i} (z , \bar {z})} = \sum_ {n > 0} \frac {\left(\partial^ {n - 1} \mu^ {j} (z , \bar {z})\right)}{(n - 1) !} <   \left[ T _ {j} T _ {i} \right] _ {n} (z) \exp \left(- \frac {1}{\pi} \int \mu^ {k} T _ {k}\right) >, \tag {2.5.18}
$$

Note that only the singular parts of the OPEs contribute. One can easily check that eq. (2.5.18) reproduces the functional equation (2.5.13).

Consider now the case of a linear algebra, i.e. the singular parts of the OPEs contain only the generators  $T_{j}$  or their derivatives. The result (2.5.18) can then be expressed in terms of functional derivatives with respect to  $\mu^{j}$  of  $Z$  as in eq. (2.5.12). Central extension terms in eq. (2.5.18) are simply proportional to  $Z$ . This means an overall factor of  $Z$  can be divided out. If all central extension terms are proportional to  $c$ , we again infer that  $\Gamma = c\Gamma^{(0)}$ .

When the OPEs close nonlinearly, i.e. contain also normal ordered expressions of the generators, it is still possible to write down a functional equation for  $Z$  or  $\Gamma$  by using eq. (2.3.6). An explicit calculation of such a case is given in section 5.3. The resulting functional equations have not been solved up to now. We will treat this problem at several points in this work, but especially in chapter 7.

# 2.6 A few examples

# 2.6.1 The free massless scalar

The action for a massless scalar propagating in a space with metric  $g_{ij}$  is given by:

$$
S _ {s} [ X, g ^ {i j} ] = - \frac {1}{4 \pi \lambda} \int d x ^ {2} \sqrt {g} g ^ {i j} \partial_ {i} X (x) \partial_ {j} X (x), \tag {2.6.1}
$$

with  $g$  the absolute value of the determinant of  $g_{ij}$  and  $\lambda$  a normalisation constant. For  $D$  scalars, this gives an action for the bosonic string in  $D$  dimensions with flat target space, see chapter 8.

Before discussing the Ward identity of the energy-momentum tensor for this action, we first observe that the action (2.6.1) has a symmetry:

$$
\delta X (x) = \varepsilon (x). \tag {2.6.2}
$$

We will follow the reasoning of subsection 2.2.2 to find the Ward identity corresponding to this symmetry. The Noether current for this symmetry is  $\frac{1}{2\lambda}\partial_iX(x)$  (see eq. (2.2.7)). It satisfies the Ward identity (see eq. (2.2.9)):

$$
\begin{array}{l} - \frac {1}{2 \pi \lambda} \int d ^ {2} x \sqrt {g} g ^ {i j} \partial_ {i} \varepsilon (x) <   \partial_ {j} X (x) X (x _ {1}) X (x _ {2}) \dots > = \\ <   \varepsilon \left(x _ {1}\right) X \left(x _ {2}\right) \dots > + <   X \left(x _ {1}\right) \varepsilon \left(x _ {2}\right) \dots > + \dots \tag {2.6.3} \\ \end{array}
$$

This symmetry has a chiral component  $\partial X$  for which the equation corresponding to (2.2.16) is:

$$
\begin{array}{l} <   \partial X (z, \bar {z}) X (z _ {1}, \bar {z} _ {1}) X (z _ {2}, \bar {z} _ {2}) \dots > = \\ \frac {\lambda}{z - z _ {1}} <   X \left(z _ {2}, \bar {z} _ {2}\right) \dots > + \frac {\lambda}{z - z _ {2}} <   X \left(z _ {1}, \bar {z} _ {1}\right) \dots > + \dots \tag {2.6.4} \\ \end{array}
$$

corresponding to the OPE<sup>10</sup>:

$$
\partial X (z) X \left(z _ {0}, \bar {z} _ {0}\right) = \frac {\lambda}{z - z _ {0}} + O \left(z - z _ {0}\right) ^ {0}. \tag {2.6.5}
$$

This agrees with the propagator for  $X$  which is given by the Green's function for the Laplacian (see appendix A):

$$
<   X (z, \bar {z}) X \left(z _ {0}, \bar {z} _ {0}\right) > = \lambda \log \left(z - z _ {0}\right) \left(\bar {z} - \bar {z} _ {0}\right). \tag {2.6.6}
$$

Note that  $X$  itself is a field whose two-point function is not a holomorphic function of  $z$  or  $\bar{z}$ . This in fact makes it impossible to use the OPE formalism with the field  $X$ . This need not surprise us, as  $X$  not a symmetry generator. We will only work with the chiral symmetry generator  $\partial X$ . For easy reference we give its OPE, which follows from eq. (2.6.5) by taking an additional derivative:

$$
\partial X (z) \partial X \left(z _ {0}\right) = \frac {\lambda}{\left(z - z _ {0}\right) ^ {2}} + O \left(z - z _ {0}\right) ^ {0} \tag {2.6.7}
$$

The action (2.6.1) is invariant under the transformation  $\delta_{\varepsilon}X = \varepsilon^{i}\partial_{i}X$ . Hence, we have that the energy-momentum tensor eq. (2.2.8) is classically conserved. It is given by:

$$
T _ {i j} = \frac {1}{2 \lambda} \partial_ {i} X (x) \partial_ {j} X (x) - \frac {1}{4 \lambda} g _ {i j} g ^ {k l} \partial_ {k} X (x) \partial_ {l} X (x). \tag {2.6.8}
$$

In the complex basis,  $T_{ij}$  has only two non-vanishing components:

$$
T _ {z z} = \frac {1}{2 \lambda} \partial X \partial X, \quad T _ {\bar {z} \bar {z}} = \frac {1}{2 \bar {\lambda}} \bar {\partial} X \bar {\partial} X. \tag {2.6.9}
$$

However, these expression are not well-defined in the quantum case, due to short-distance singularities. We use point-splitting regularisation to define the quantum operator. In the OPE formalism, this becomes  $T = \frac{1}{2\lambda} [\partial X\partial X]_0$ . We can now use the OPE (2.6.7), and the rules of subsection 2.3.2 to compute the OPE of  $T$  with itself. One ends up with the correct Virasoro OPE (2.3.2), with a central charge  $c = 1$ . Furthermore, it is easy to check that  $X$  is a primary field with respect to  $T$

of dimension zero. Moreover,  $\partial X$  is also primary, having dimension one. In the rest of this subsection we set  $\lambda = 1$ .

Let us now check for the mode algebra eq. (2.4.4) corresponding to the OPE of  $\partial X$  with itself eq. (2.6.7). The modes of  $\partial X$  are traditionally denoted with  $\alpha_{m}$ . We find:

$$
[ \alpha_ {m}, \alpha_ {n} ] = m \delta_ {m + n}, \tag {2.6.10}
$$

which is related to the standard harmonic oscillator commutation relations via a rescaling with  $\sqrt{|m|}$ .

In string theory, computing scattering amplitudes is done by inserting local operators of the correct momentum in the path integral and integrating over the coordinates [127, 105]. These local operators are (composites with) normal ordered exponentials of the scalar field. Because  $X$  cannot be treated in the present OPE scheme, we should resort to different techniques to define these exponentials. This is most conveniently done using a mode expansion of  $X$ . One shows that one can define a chiral operator which we write symbolically as:

$$
V _ {a} (z) =: \exp (a X (z)):, \quad a \in \mathbb {C}, \tag {2.6.11}
$$

for which the following identities hold:

$$
\partial V _ {a} = a: \partial X V _ {a}: \tag {2.6.12}
$$

$$
\partial X (z) V _ {a} (w) = a \frac {V _ {a} (w)}{(z - w)} + O (z - w) ^ {0} \tag {2.6.13}
$$

and

$$
\begin{array}{l} V _ {a} (z) V _ {b} (w) \\ = (z - w) ^ {a b}: \exp (a X (z)) \exp (b X (w)): \\ = (z - w) ^ {a b}: \left(1 + (z - w) a \partial X (w) + \right.. \\ \left. \frac {(z - w) ^ {2}}{2} \left(a \partial^ {2} X (w) + a ^ {2} \partial X (w) \partial X (w)\right) \dots\right) V _ {a + b} (w):, \tag {2.6.14} \\ \end{array}
$$

where the last line is a well-defined expression where normal ordering, as defined in eq. (2.3.6), from right to left is understood. We will use these equations as the definition for the vertex operators.

Using the eqs. (2.6.12) and (2.6.13), it is easy to check that  $V_{a}$  is primary with conformal dimension  $a^2 / 2$  with respect to  $T$  (2.6.9). This differs from the classical dimension which is zero as  $X$  has dimension zero.

The OPE (2.6.14) is distinctly different from the ones considered before. Indeed, noninteger poles are possible. The rules constructed in section 2.3.2 are not valid

# 2.6. A few examples

for such a case. However, every OPE with vertex operators we will need, will follow immediately from the above definitions. Another peculiarity is that when  $ab \in \mathbb{Z}$ , the definition (2.6.14) fixes the regular part of the OPE of two vertex operators in terms of normal ordered products of  $\partial X$  and a vertex operator. An interesting example of this is for  $V_{\pm 1}$ :

$$
V _ {\pm 1} (z) V _ {\pm 1} (w) = (z - w) V _ {2} (w) + \dots
$$

$$
V _ {\pm 1} (z) V _ {\mp 1} (w) = \frac {1}{(z - w)} \pm \partial X (w) + \dots \tag {2.6.15}
$$

We see that we have two operators:

$$
\frac {V _ {1} + V _ {- 1}}{\sqrt {2}} = \sqrt {2} \cosh X \quad \text {a n d} \quad \frac {V _ {1} - V _ {- 1}}{\sqrt {2} i} = - i \sqrt {2} \sinh X \tag {2.6.16}
$$

which satisfy the OPEs of two free fermions (see eq. (2.6.26) below). Vertex operators can also be used to construct realisations of affine Lie algebras, as we will see in section 4.6.

When considering realisations of  $\mathcal{W}$ -algebras using scalars, it is obviously a problem that the energy-momentum tensor of the free scalars is restricted to an integer central charge  $n$  by considering  $n$  free bosons. However, we can modify the eq. (2.6.9) to:

$$
T = \frac {1}{2}: \partial X \partial X: - q \partial^ {2} X, \tag {2.6.17}
$$

which has a central charge  $1 - 12q^2$ . This corresponds to adding a background charge to the free field action [73, 91, 59].

# 2.6.2 The free Majorana fermion

The free Majorana fermion in two dimensions has the following action in the conformal gauge:

$$
S _ {f} [ \psi ] = \frac {1}{2 \pi \lambda} \int d ^ {2} x \psi (x) \bar {\partial} \psi (x), \tag {2.6.18}
$$

where  $\lambda$  is a normalisation constant. The equation of motion for  $\psi$  is:

$$
\bar {\partial} \psi = 0. \tag {2.6.19}
$$

Hence,  $\psi$  is a chiral field and we will write  $\psi (z)$ . Note that under conformal transformations it has to transform as a primary field of dimension  $h = 1 / 2$ ,  $\bar{h} = 0$  to have a conformal invariant action.  $-\frac{1}{\lambda}\psi (x)$  is the Noether current corresponding to the transformation:

$$
\delta_ {\varepsilon} \psi (x) = \varepsilon (x), \tag {2.6.20}
$$

where  $\varepsilon$  is a Grasmann-odd field. Although we can proceed as before, we will derive the OPEs of  $\psi$  using a different technique. Consider the generating functional (see eq. (2.5.1)):

$$
Z [ \mu ] = \frac {1}{\mathcal {N}} \int [ d \psi ] \exp \left(- S _ {f} [ \psi ] - \frac {1}{\pi} \int d ^ {2} x \mu (x) \psi (x)\right), \tag {2.6.21}
$$

where  $\mu$  is Grasmann-odd. Such a pathintegral can be computed by converting it to a Gaussian path integral, as we will now show. We start by rewriting the exponent in the following way:

$$
Z [ \mu ] = \frac {1}{\mathcal {N}} \int [ d \psi ] \exp \left(- \frac {1}{2 \pi \lambda} \int (\psi + \lambda (\bar {\partial} ^ {- 1} \mu)) \bar {\partial} (\psi + \lambda \bar {\partial} ^ {- 1} \mu)\right)
$$

$$
\exp \left(\frac {\lambda}{2 \pi} \int (\bar {\partial} ^ {- 1} \mu) \mu\right), \tag {2.6.22}
$$

where we used the definition of the inverse derivative of appendix A, supposing that  $\mu$  decays sufficiently fast at infinity. One can now shift the integration variables to  $\tilde{\psi} = \psi +\lambda \bar{\partial}^{-1}\mu$ . This shift has a Jacobian 1. One arrives at:

$$
Z [ \mu ] = \frac {1}{\mathcal {N}} \int [ d \tilde {\psi} ] \exp \left(\frac {- 1}{2 \lambda \pi} \int \tilde {\psi} \bar {\partial} \tilde {\psi}\right) \exp \left(\frac {\lambda}{2 \pi} \int (\bar {\partial} ^ {- 1} \mu) \mu\right). \tag {2.6.23}
$$

The remaining path integral cancels  $\mathcal{N}$  (see eq. (2.2.2)) and we have:

$$
Z [ \mu ] = \exp \left(- \frac {\lambda}{2 \pi^ {2}} \int d x ^ {2} d x _ {0} ^ {2} \mu (z, \bar {z}) \frac {1}{z - z _ {0}} \mu \left(z _ {0}, \bar {z} _ {0}\right)\right). \tag {2.6.24}
$$

From this result we immediately see that the one-point function  $< \psi(x) >$  is zero, while the two-point function is given by:

$$
\begin{array}{l} <   \psi (z) \psi (z _ {0}) > = \pi^ {2} \frac {\delta}{\delta \mu (z , \bar {z})} \frac {\delta}{\delta \mu (z _ {0} , \bar {z} _ {0})} Z [ \mu ] \Big | _ {\mu = 0} \\ = \frac {\lambda}{z - z _ {0}}, \tag {2.6.25} \\ \end{array}
$$

where we used left-functional derivatives. The corresponding OPE is<sup>11</sup>:

$$
\psi (z) \psi (w) = \frac {\lambda}{z - w} + O (z - w) ^ {0}. \tag {2.6.26}
$$

To find the energy-momentum tensor for the action (2.6.18), we use the definition (2.2.7) with  $\varepsilon^i$  decaying at infinity, together with the transformation law (2.1.6) for  $\psi$ . We find:

$$
\begin{array}{l} \delta S _ {f} = \frac {1}{\lambda \pi} \int d ^ {2} x \left(\varepsilon^ {z} \partial \psi + 1 / 2 \partial \varepsilon^ {z} \psi + \varepsilon^ {\bar {z}} \bar {\partial} \psi\right) \bar {\partial} \psi \\ = \frac {1}{2 \lambda \pi} \int d ^ {2} x \bar {\partial} \varepsilon^ {z} \psi \partial \psi . \tag {2.6.27} \\ \end{array}
$$

Hence, we find only one component  $T^{\bar{z}\bar{z}}$  non-zero, according to (2.2.12) we have:

$$
T = \frac {1}{2 \lambda}: \partial \psi \psi :, \tag {2.6.28}
$$

which has central charge  $c = 1/2$ .

# 2.6.3 Other first order systems

We will need two other first order systems, the fermionic  $b, c$  and the bosonic  $\beta, \gamma$  system [91]. They will arise as the ghosts in the BRST-quantisation of systems with conformal invariance. We will combine both using a supermatrix notation (see appendix B). The action is:

$$
S _ {\mathrm {b c}} [ b, c ] = \frac {1}{\pi} \int d ^ {2} x b ^ {i} (x) _ {i} A ^ {j} \bar {\partial} (_ {j} c) (x), \tag {2.6.29}
$$

where  ${}_iA^j$  is a constant bosonic supermatrix which is invertible. We assume that it does not mix fermions and bosons. We also take  $b$  and  $c$  to be bosonic matrices<sup>12</sup>. The equations of motion again show that  $b^i, c_i$  are chiral fields. Conformal invariance requires the  $b^i$  to be primary fields with dimension  $h^i$  and  $c_i$  also to be primary such that  $\dim(_i A^j)_j c = 1 - h^i$  (and zero for  $\bar{h}$ ). We proceed now as in the previous subsection:

$$
\begin{array}{l} Z [ \mu , \nu ] = \frac {1}{\mathcal {N}} \int [ d b ^ {i} ] [ d c _ {j} ] \exp - \left(S _ {\mathrm {b c}} [ b, c ] + \frac {1}{\pi} \int \mu_ {k} ^ {k} b + \nu^ {k} _ {k} c\right) \\ = \frac {1}{\mathcal {N}} \int [ d b ^ {i} ] [ d c _ {j} ] \exp \left(\frac {1}{\pi} \int (b + \bar {\partial} ^ {- 1} \nu A ^ {- 1}) A \bar {\partial} (c - (A \bar {\partial}) ^ {- 1} \mu))\right) \\ \exp \left(\frac {1}{\pi} \int (\bar {\partial} ^ {- 1} \nu) A ^ {- 1} \mu\right) \\ = \exp \left(- \frac {1}{\pi} \int d ^ {2} x d ^ {2} x _ {0} \nu^ {i} (z, \bar {z}) \frac {i (A ^ {- 1}) ^ {j}}{z - z _ {0}} j \mu \left(z _ {0}, \bar {z} _ {0}\right)\right). \tag {2.6.30} \\ \end{array}
$$

Carefully keeping track of the signs, we get for the two-point functions:

$$
\begin{array}{l} <   _ {i} c (z) b ^ {j} \left(z _ {0}\right) > = \pi^ {2} \frac {\delta}{\delta \nu^ {i} \left(z , \bar {z}\right)} \frac {\delta}{\delta \mu_ {j} \left(z _ {0} , \bar {z} _ {0}\right)} Z [ \mu , \nu ] \Big | _ {\mu , \nu = 0} \\ = \frac {- i (A ^ {- 1}) ^ {j}}{z - z _ {0}}. \tag {2.6.31} \\ \end{array}
$$

This OPE does not depend on the dimensions of the fields, while the energy-momentum tensor of course does. To find the energy-momentum tensor, we proceed as in the previous subsection. From:

$$
\delta_ {\varepsilon} (i (A c)) = \varepsilon^ {z} \partial (i (A c)) + \left(1 - h ^ {i}\right) \partial \varepsilon^ {z} _ {i} (A c), \tag {2.6.32}
$$

we get:

$$
T _ {\mathrm {b c}} =: b A \partial c - \left(1 - h ^ {i}\right) \partial \left(b ^ {i} {} _ {i} A ^ {j} {} _ {j} c\right):, \tag {2.6.33}
$$

with central charge:

$$
c = \sum_ {i} (- 1) ^ {i} 2 \left(6 \left(h ^ {i}\right) ^ {2} - 6 h ^ {i} + 1\right), \tag {2.6.34}
$$

where the phase factor is  $-1$  for  $b^{i}$  fermionic. Note that the second term of eq. (2.6.33) is the derivative of the ghost current. It is not present in eq. (2.6.28) because :  $\psi \psi \coloneqq 0$

# 2.6.4 Wess-Zumino-Novikov-Witten models

A final example of a conformal field theory is the WZNW-model [156, 209]. It is a nonlinear sigma model with as target space a group manifold. The (super) Lie group  $\mathcal{G}$  is required to be semisimple (see [154, 81] for a weakening of this condition). We denote its (super) Lie algebra with  $\bar{g}$ , see appendix B for conventions.

The WZWN action  $\kappa S^{+}[g]$  is a functional of a  $\mathcal{G}$ -valued field  $g$  and is given by:

$$
\begin{array}{l} \kappa S ^ {+} [ g ] = \frac {\kappa}{4 \pi x} \int_ {\partial \Omega} d ^ {2} x s t r \left\{\partial g ^ {- 1} \bar {\partial} g \right\} \\ + \frac {\kappa}{1 2 \pi x} \int_ {\Omega} d ^ {3} x \varepsilon^ {\alpha \beta \gamma} s t r \left\{g _ {, \alpha} g ^ {- 1} g _ {, \beta} g ^ {- 1} g _ {, \gamma} g ^ {- 1} \right\}, \tag {2.6.35} \\ \end{array}
$$

where  $\Omega$  is a three-manifold with boundary  $\partial \Omega$ . It satisfies the Polyakov-Wiegman identity [161]:

$$
S ^ {+} [ h g ] = S ^ {+} [ h ] + S ^ {+} [ g ] - \frac {1}{2 \pi x} \int s t r (h ^ {- 1} \partial h \bar {\partial} g g ^ {- 1}), \tag {2.6.36}
$$

which is obtained through direct computation. We also introduce a functional  $S^{-}[g]$  which is defined by:

$$
S ^ {-} [ g ] = S ^ {+} \left[ g ^ {- 1} \right]. \tag {2.6.37}
$$

# 2.6. A few examples

It is with this functional that we now continue. Using the Polyakov-Wiegman identity (2.6.36), we can show that the action  $S^{-}[g]$  transforms under:

$$
\delta_ {\eta} g = \eta g \tag {2.6.38}
$$

as

$$
\delta_ {\eta} S ^ {-} [ g ] = - \frac {\kappa}{2 \pi x} \int s t r \left(\left(\partial g g ^ {- 1}\right) \bar {\partial} \eta\right). \tag {2.6.39}
$$

We see that  $S^{-}[g]$  is invariant when  $\bar{\partial}\eta = 0$ . The corresponding conserved current is:

$$
J _ {z} = \frac {\kappa}{2} \partial g g ^ {- 1}, \tag {2.6.40}
$$

which is chiral,  $\bar{\partial} J_{z} = 0$ . Similarly, for  $\delta_{\bar{\eta}}g = g\bar{\eta}$  with  $\partial \bar{\eta} = 0$ , we find a conserved current:

$$
J _ {\bar {z}} = - \frac {\kappa}{2} g ^ {- 1} \bar {\partial} g. \tag {2.6.41}
$$

$J_{z}$  transforms under eq. (2.6.38) as:

$$
\begin{array}{l} \delta_ {\eta} J _ {z} = \frac {\kappa}{2} \partial \eta + [ \eta , J _ {z} ], (2.6.42) \\ = - \int d y \left\{s t r (\eta (y) J _ {z} (y)), J _ {z} (x) \right\} _ {\mathrm {P B}} (2.6.43) \\ \end{array}
$$

where the Poisson bracket is given by<sup>13</sup>:

$$
J _ {z} ^ {a} \times J _ {z} ^ {b} = \ll - \frac {\kappa}{2} g ^ {a b} | ^ {a} g ^ {d} J _ {z c d} ^ {c} f ^ {b} \gg , \tag {2.6.44}
$$

which defines a current algebra of level  $\kappa$ . It can be argued that the relation eq. (2.6.44) does not renormalise when going to the quantum theory, and we will take (2.6.44) as the definition of the OPE. However, the relation (2.6.40) can be renormalised to:

$$
J _ {z} = \frac {\alpha_ {\kappa}}{2} \partial g g ^ {- 1}. \tag {2.6.45}
$$

Using OPE techniques, it is argued in [135] that for a current algebra of level  $\kappa$ , normalising the currents as in (2.6.44) gives:

$$
\alpha_ {\kappa} = \kappa + \tilde {h}, \tag {2.6.46}
$$

where the dual Coxeter number  $\tilde{h}$  is the eigenvalue of the quadratic Casimir in the adjoint representation (see also appendix B). This follows from consistency requirements in the OPA of the currents with  $g$ . We will need this renormalisation in chapter 7.

The OPA generated by the currents  $J_{z}^{a}$  with the OPEs (2.6.44) is known as a Kac-Moody algebra or affine Lie algebra. Kac-Moody algebras were studied in the mathematical literature [126, 151] much earlier than WZNW-models. For a review on the algebraic aspects of Kac-Moody algebras, see [101]. We will denote the Kac-Moody algebra corresponding to  $\bar{g}$  as  $\hat{g}$ .

One can check that the Sugawara tensor:

$$
T _ {S} = \frac {1}{x (\kappa + \tilde {h})} \operatorname {s t r} \left[ J _ {z} J _ {z} \right] _ {0}, \tag {2.6.47}
$$

satisfies the Virasoro algebra with the central extension given by:

$$
c = \frac {k \left(d _ {B} - d _ {F}\right)}{k + \tilde {h}}. \tag {2.6.48}
$$

The currents  $J_{z}^{a}$  have conformal dimension 1 with respect to  $T_{S}$ .

We now briefly review some basic formulas for gauged WZNW-models [156, 209, 161, 3, 58] which we will need in chapter 6. We can gauge the symmetries generated by  $J_{z}$  by introducing gauge fields  $A_{\bar{z}}$ . The action:

$$
\kappa S ^ {-} [ g ] + \frac {1}{\pi x} \int d ^ {2} x s t r \left(J _ {z} (x) A _ {\bar {z}} (x)\right) \tag {2.6.49}
$$

is (classically) invariant when the gauge fields transform as:

$$
\delta_ {\eta} A _ {\bar {z}} = \bar {\partial} \eta + [ \eta , A _ {\bar {z}} ]. \tag {2.6.50}
$$

The induced action (2.5.4)  $\Gamma[A_{\bar{z}}]$  for the gauge fields  $A_{\bar{z}}$  is:

$$
\exp \left(- \Gamma [ A _ {\bar {z}} ]\right) = \langle \exp - \frac {1}{\pi x} \int d ^ {2} x s t r \Big (J _ {z} (x) A _ {\bar {z}} (x) \Big) \rangle . \tag {2.6.51}
$$

By using the OPEs eq. (2.6.44), and the general formula for the functional equation for  $\Gamma$  (2.5.18) we find the following Ward identity:

$$
\bar {\partial} u _ {z} - \left[ A _ {\bar {z}}, u _ {z} \right] = \partial A _ {\bar {z}}, \tag {2.6.52}
$$

where:

$$
u _ {z} ^ {a} (x) = - \frac {2 \pi}{\kappa} g ^ {a b} \frac {\delta \Gamma [ A _ {\bar {z}} ]}{\delta^ {b} A _ {\bar {z}} (x)}. \tag {2.6.53}
$$

The Ward identity is independent of  $\kappa$ , therefore:

$$
\Gamma \left[ A _ {\bar {z}} \right] = \kappa \Gamma^ {(0)} \left[ A _ {\bar {z}} \right], \tag {2.6.54}
$$

where  $\Gamma^{(0)}[A_{\bar{z}}]$  is independent of  $\kappa$ . In [161, 3], it was observed that eq. (2.6.52) states that the curvature for the Yang-Mills fields  $\{A,u\}$  vanishes.

# 2.A Appendix: Combinatorics

This appendix combines some formulas for binomial coefficients and related definitions. We define first the Pochhammer symbol:

$$
(a) _ {n} = \frac {\Gamma (a + n)}{\Gamma (a)} = \prod_ {j = 0} ^ {n - 1} (a + j), \quad a \in \mathbb {R}, n \in \mathbb {N} \tag {2.A.1}
$$

where the last definition is also valid if  $a$  is a negative integer. We then define the binomial as:

$$
\binom {a} {n} = \frac {(a) _ {n}}{n !} \quad a \in \mathbb {R}, n \in \mathbb {N}. \tag {2.A.2}
$$

In particular,

$$
\binom {- n} {m} = (- 1) ^ {m} \binom {m + n - 1} {m}, \quad m, n \in \mathbb {N}. \tag {2.A.3}
$$

We now list some identities for sums of binomial coefficients which are used in section 2.4.

$$
\sum_ {l \in \mathbb {N}} (- 1) ^ {l} \binom {r} {l} \binom {r + s - l} {p - l} = \binom {s} {p}, r, s \in \mathbb {Z}, p \in \mathbb {N} \tag {2.A.4}
$$

which can be proven by multiplying with  $x^p$  and summing over  $p$ .

$$
\sum_ {p \in \mathbb {N}} \binom {n} {p} \binom {m} {s - p} = \binom {m + n} {s}, \quad m, n \in \mathbb {Z}, s \in \mathbb {N} \tag {2.A.5}
$$

which is proven by multiplying with  $x^{s}$  and summing over  $s$ .

# Chapter 3

# Operator Product Expansions in Mathematica

The previous chapter indicates that the method of Operator Product Expansions (OPEs) gives us a very useful algebraic tool to study a conformal field theory. However, the relevant formulas show that calculations can become quite tedious and error-prone, especially when nested normal ordered products are involved. This was the motivation to implement an algorithm to handle OPEs and normal ordered products. By now, the OPErefs package is widely used, showing that the need for automated OPE calculations was (and is) present among conformal field theorists.

This chapter discusses the OPErefs package. In the first section, the desired capabilities of the package are discussed and a brief history of the development of OPErefs is given. Section 3.2 is devoted to design considerations. The implementation we used is explained in section 3.3. We first give the algorithm, and prove that it is finite. Then some crucial points in the actual code of OPErefs are highlighted. Section 3.4 is a user's guide to OPErefs, and the next section presents an explicit example of a calculation. Section 3.6 gives an outlook on possible extensions. We end this chapter with a review of other symbolic manipulation packages used in conformal field theory.

Some of the material found in this chapter was published in [192] and [193].

All registered trademarks used in this chapter are acknowledged.

# 3.1 Intention and history

We want to construct a package which, given the OPEs of a set of generators, is able to compute any OPE in the Operator Product Algebra  $(\mathrm{OPA})^1$ . Moreover, as the definition of normal ordering is noncommutative and nonassociative, the package should bring composite operators in a standard form. This can be done by introducing an order for the generators and reordering composites such that generators are ordered from left to right, using the associativity conditions of the previous chapter.

We will restrict ourselves to meromorphic conformal field theory. As we use only the concept of an OPA (see def. 2.3.2), we will be able to use operators with any (also negative) conformal dimension.

An important step in the calculations with OPEs is to check that the Jacobi identities (2.3.21) discussed in the previous chapter are satisfied. The package should be able to verify this. This will also enable us to fix some free parameters in the OPEs of the generators by requiring associativity.

A first version of my package was used to perform the calculations on the Casimir algebra of  $B_{2}$  [78] (see section 4.6), where OPEs of composites of up to four fields are involved. After some optimisations, OPErefs 2.0 [192] was made available to the public. Several well-known results were checked with this version, e.g. the Sugawara construction for  $\widehat{SU(n)}$ , the free field realisations and homological constructions of  $W_{3}$  based on a  $\widehat{SU(3)}$ -Kač-Moody algebra [56]. The power of having all calculations automated was already clear at this point, as we found a small misprint in  $[56]^{2}$ , while this calculation was highly nontrivial by hand.

Important improvements (version 3.0) were published in [193]. It was necessary to introduce a different syntax, to be able to make extensive changes in the future. Since then, many other small changes have been made to the package. The latest version, which will be discussed in this chapter, is 3.1.

The package is freely available from the author<sup>3</sup>.

# 3.2 Design considerations

In this section, we comment on how OPErefs was designed. Several aspects are important.

# 3.3. Implementation

- How will the program present itself to the user? (interface)

Preferably, concepts and terminology that the user has to handle should be familiar to any conformal field theorist.

- How will the data be represented?

The representation has to be chosen such that storage space and computation time are minimised. Furthermore, it should be possible to incorporate any related problems (for instance Poisson brackets) without changing the data representation.

Whatever representation is chosen, it is important to hide this choice for the user. Ideally, one should be able to change the representation (and thus the internal working of the program) without affecting the interface.

- Make sure the program can be easily extended.

It is important to keep the main ideas of the algorithm reflected in the structure of the program, separated from the actual details of the computation. This implies a modular style of writing. Of course, a high level programming language enhances readability and makes changes and additions easier.

- How can the results be manipulated?

It is of course not enough to be able to compute a certain OPE or normal ordered product. Different results have to be added, coefficients simplified or extracted, equations solved... It would lead the programmer of the package way to far if he/she had to implement all this functionality.

These remarks point towards an implementation in an interactive environment where it is possible to perform symbolic computations, and where a high level language is available, i.e. in a symbolic manipulation program. Indeed, these environments are developed to perform the dull - or algorithmic - calculations, once the necessary rules are implemented. The user can then use the programming power of the environment to solve his problem. It would be perfectly possible to write the proposed program in an "ordinary" programming language like  $C^{++}$  or lisp. After all, the symbolic manipulation programs are written in one of these languages. However, the amount of code required would probably go beyond the scope of a physics PhD.

Many of the design principles we mentioned, are reminiscent of object oriented programming. As this is a fairly new concept, one immediately is directed towards the more recent environments. Of these, the best known and more powerful are Axiom [123], Maple [36] and Mathematica [210]. All three are suitable for implementation of the algorithm for working with OPEs.

From the viewpoint of the programmer, Axiom seems to be the most attractive choice. It is indeed the most recent environment, and its design reflects all points made above. The major disadvantages of this system are that it is not very popular yet, and the hardware requirements are very high (a powerful Unix workstation is

needed).

Maple is being used more and more over the last few years. Its design is concentrated on run-time efficiency and small memory requirements. Presumably, a consequence of this is that the Maple programming language is more like  $C$  or Pascal, although with powerful extensions for symbolic manipulations. However, it does not seem to permit a high abstraction level and modularity compared to the other two environments.

Finally, the Mathematica environment is somewhere in between the two others. It is probably not as efficient as Maple, but it is a lot easier to program in. It runs on a lot of machines, from PC's to supercomputers and all versions are completely compatible. Many of the ideas of object oriented programming can be used in Mathematica and its pattern matching capabilities are very useful for any mathematical programming.

# 3.3 Implementation

In this section, we first discuss the algorithm which we used. In principle, it simply consists of repeated application of the rules given in subsection 2.3.2. In the next subsection, we give some parts of the actual code of OPErefs. Finally, we comment on the performance of the package.

# 3.3.1 Algorithm

We will discuss the algorithm in two steps. First, we show how the singular part of the OPEs can be computed. Then we give the algorithm for simplifying composites. In OPErefs, composites are simplified at each step of the calculation.

OPEdefines 3.1 can calculate OPEs and Poisson brackets. We concentrate on OPEs here and briefly comment which changes are needed for the other case.

# Computing OPEs

The package should compute OPEs of arbitrarily complicated composites when a set of generators and their OPEs is given. Obviously, we will need the rules constructed in subsection 2.3.2. For easy reference, we list the rules we need, aside from linearity. First, there are the rules involving derivatives:

$$
[ \partial A B ] _ {q} = - (q - 1) [ A B ] _ {q - 1} \tag {3.3.1}
$$

$$
[ A \partial B ] _ {q} = (q - 1) [ A B ] _ {q - 1} + \partial [ A B ] _ {q}. \tag {3.3.2}
$$

Next, we need the OPE of  $B$  with  $A$ , given the OPE of  $A$  with  $B$ :

$$
[ B A ] _ {q} = (- 1) ^ {| A | | B |} \sum_ {l \geq q} \frac {(- 1) ^ {l}}{(l - q) !} \partial^ {(l - q)} [ A B ] _ {l}. \tag {3.3.3}
$$

OPEs with composites can be calculated using eq. (2.3.21):

$$
\begin{array}{l} [ A [ B C ] _ {0} ] _ {q} = (- 1) ^ {| A | | B |} [ B [ A C ] _ {q} ] _ {0} + \left[ \left[ A B \right] _ {q} C \right] _ {0} \\ + \sum_ {l = 1} ^ {q - 1} \binom {q - 1} {l} \left[ [ A B ] _ {q - l} C \right] _ {l} \tag {3.3.4} \\ \end{array}
$$

where  $q\geq 1$

These rules are the only ones needed to compute every OPE in the OPA. Indeed, when computing the OPE of  $A$  with  $B$ , we apply the following procedure:

- if  $A$  and  $B$  are generators whose OPE we know, return it as the result.  
- apply linearity if necessary.  
- if  $A$  is an operator with derivatives, use eq. (3.3.1).  
- if  $B$  is an operator with derivatives, use eq. (3.3.2).  
- if  $B$  is a composite, use eq. (3.3.4).  
- if  $A$  is a composite, use eq. (3.3.3).  
- if the OPE  $B(z) A(w)$  is known, compute the OPE  $A(z) B(w)$  using eq. (3.3.3).

This list should be used recursively until none of the rules applies, which means that the OPE has been calculated. The order in which we check the rules is in this case not important, but we will check them in a "top-down" order. Note that to compute an OPE of a composite with a generator, first eq. (3.3.3) is used, and eq. (3.3.4) in the next step. The question now is if this is a finite procedure.

First of all, all sums in the previous equations have a finite number of terms. This follows from the assumption that all OPEs between generators in the OPA have a finite number of poles. The only possible problem is that eq. (3.3.4) would require an infinite number of steps. Indeed, when computing the OPE  $A[BC]_0$ , one needs all OPEs  $[AB]_nC$  for  $n > 0$ . It is now possible that  $[AB]_n$  contains a composite again, leading to a new application of eq. (3.3.4). For graded OPAs (see def. 2.3.7), we can prove the following theorem.

Theorem 3.3.1 If the OPA is graded and there is a lower bound  $D$  on the dimension of all the operators in the OPA, then the application of eqs. (3.3.1-3.3.4) requires only a finite number of steps.

# Proof :

The sum of the dimensions of the fields in the OPE  $[A_1[B_1C_1]_0]$  is  $\Delta = a_{1} + b_{1} + c_{1}$ . The application of eq. (3.3.4) leads to the computation of the OPEs  $[\left[AB\right]_nC]$  for  $n > 0$ .

The sum of the dimensions is in this case  $a_1 + b_1 + c_1 - n$ . As this is strictly smaller than  $\Delta$ , we see that after (the integer part of)  $\Delta - 2D$  steps, we would need the OPE of a field with dimension lower than  $D$ , which is not present in the OPA by assumption.

Note that the assumption of the existence of a lower bound on the dimension is not very restrictive. First of all, all conformal OPAs (see def. 2.3.6) with generators of positive dimension trivially satisfy the criterion. Moreover, a finite number of fermionic generators can have negative dimension. Indeed, to form composites of negative dimension with these fermions, one needs extra derivatives - e.g.  $[\partial A A]_0$  imposed by eq. (3.3.8) below.

An important remark here is that the dimension used in the proof is not necessarily the conformal dimension of the operators. For instance, for a bosonic  $bc$ -system which commutes with all the other generators of the OPA, we can assign a dimension of  $1/2$  to both  $b$  and  $c$ , such that the theorem can be used.

The algorithm used to compute Poisson brackets is essentially the same. Some small changes in the rule eq. (3.3.4), as discussed in subsection 2.3.5, are needed.

# Simplifying composites

We now discuss how we can reduce normal ordered products to a standard form. We define an order on the generators and their derivatives, e.g. lexicographic ordering. Given a composite, we apply the rules of subsection 2.3.2 until all composites are normal ordered from right to left and the operators are ordered, i.e.  $[A[B[C\ldots]_0]_0]$ , and  $A \leq B \leq C$ . The relevant formulas are:

$$
\partial [ A B ] _ {0} = [ A \partial B ] _ {0} + [ \partial A B ] _ {0} \tag {3.3.5}
$$

$$
[ A B ] _ {- q} = \frac {1}{q !} [ (\partial^ {\varphi} A) B ] _ {0}, \quad q \geq 1 \tag {3.3.6}
$$

$$
[ B A ] _ {0} = (- 1) ^ {| A | | B |} [ A B ] _ {0} + (- 1) ^ {| A | | B |} \sum_ {l \geq 1} \frac {(- 1) ^ {l}}{l !} \partial^ {l} [ A B ] _ {l} \tag {3.3.7}
$$

$$
[ A A ] _ {0} = - \sum_ {l > 0} \frac {(- 1) ^ {l}}{2 l !} \partial^ {l} [ A A ] _ {l}, \quad \text {f o r} A \text {f e r m i o n i c}. \tag {3.3.8}
$$

$$
[ A [ B C ] _ {0} ] _ {0} = (- 1) ^ {| A | | B |} [ B [ A C ] _ {0} ] _ {0} + \tag {3.3.9}
$$

$$
[ ([ A B ] _ {0} - (- 1) ^ {| A | | B |} [ B A ] _ {0}) C ] _ {0} \tag {3.3.10}
$$

$$
[ A [ A C ] _ {0} ] _ {0} = \left[ \left[ A A \right] _ {0} C \right] _ {0}, \quad \text {f o r} A \text {f e r m i o n i c}. \tag {3.3.11}
$$

The rule (3.3.10) follows from eq. (2.3.21).

# 3.3. Implementation

We have to wonder if we again require only a finite number of steps to reach the end-result. We can not prove this along the same lines as theorem 3.3.1. Indeed, in eq. (3.3.10) all terms necessarily have the same dimension. We give an example of infinite recursion in the intermezzo.

# Intermezzo 3.3.1

Consider an OPA with bosonic operators  $A, B$  such that:

$$
[ A B ] _ {1} = a B. \tag {3.3.12}
$$

We do not specify the OPE of  $B$  with itself. We want to reorder  $[[AB]_0B]_0$  into standard order, i.e.  $[A[BB]_0]_0 +$  corrections. In the computation which follows, an ellipsis denotes terms containing more than one derivative or containing  $[BB]_n$  with  $n > 0$ . (The equation numbers at the right correspond to the formula used in a certain step.)

$$
\begin{array}{l} [ [ A B ] _ {0} B ] _ {0} = \left. [ B [ A B ] _ {0} ] _ {0} - \partial [ B [ A B ] _ {0} ] _ {1} + \dots \right. (3.3.7) \\ = [ A [ B B ] _ {0} ] _ {0} + [ \partial [ B A ] _ {1} B ] _ {0} - e q. (3.3.10) \\ \partial \left([ A [ B B ] _ {1} ] _ {0} + \left[ [ B A ] _ {1} B \right] _ {0}\right) + \dots \quad e q. (3.3.4) \\ = [ A [ B B ] _ {0} ] _ {0} - a [ \partial B B ] _ {0} + a \partial [ B B ] _ {0} + \dots e q s. (3. 3. 5), (3. 3. 7) \\ = [ A [ B B ] _ {0} ] _ {0} + a [ \partial B B ] _ {0} + \dots \\ \end{array}
$$

Now, suppose that there is a null operator:

$$
\partial B - \frac {1}{a} [ A B ] _ {0} = 0. \tag {3.3.13}
$$

in the OPA and we use this relation to eliminate  $\partial B$ . We find:

$$
[ [ A B ] _ {0} B ] _ {0} = [ A [ B B ] _ {0} ] _ {0} + [ [ A B ] _ {0} B ] _ {0} + \dots ,
$$

which means that we cannot determine  $[[AB]_0B]_0$  in this way. Note that there is no problem if we use eq. (3.3.13) to eliminate  $[AB]_0$ .

This example occurs for a scalar field  $X$  and the vertex operators defined in subsection 2.6.1. Take  $A = \partial X$  and  $B =: \exp(aX)$ . If we normalised  $X$  such that  $[\partial X \partial X]_2 = 1$ , eq. (3.3.12) is satisfied. On the other hand,  $\partial(: \exp(aX)) = a[\partial X : \exp(aX)]_0$ . Clearly, this corresponds to eq. (3.3.13) if  $a = 1$ .

The example in the intermezzo shows that infinite recursion can occur when null operators are present. If  $N$  is a null operator, we call  $N = 0$  a null-relation. As shown in the intermezzo, the algorithm does not work if we use a null-relation to eliminate operators with "less" derivatives in favour of those with "more" derivatives.

We make this more precise. For an operator  $A$ , we define  $d(A)$  as the number of derivatives occurring in  $A$ . When we do not use any null-relations, we have:

$$
\begin{array}{l} d (A) = 0 \quad \text {i f} A \text {i s a g e n e r a t o r} \\ d (\partial A) = d (A) + 1 \\ d ([ A B ] _ {0}) = d (A) + d (B) \\ d (A + B) = \min  (d (A), d (B)). \tag {3.3.14} \\ \end{array}
$$

We say that we "use a null-relation in increasing derivative-number" if we use it to eliminate (one of) the operator(s) with the smallest derivative-number. In this case, two of the relations eq. (3.3.14) are converted to inequalities:  $d(\partial^n A) \geq n$  and  $d([AB]_0) \geq d(A) + d(B)$ .

We are now in a position to prove the following theorem.

Theorem 3.3.2 In an OPA satisfying the conditions of theorem 3.3.1, and where any null-relations are used in increasing derivative-number, the application of eqs. (3.3.5-3.3.11) (together with eqs. (3.3.1-3.3.4)) to order the operators in a normal ordered product requires only a finite number of steps.

# Proof :

We define:

$$
n \equiv d \left(\left[ A [ B C ] _ {0} \right] _ {0}\right).
$$

We wish to show that the derivative number of the correction terms in eq. (3.3.10) is always bigger than  $n$ :

$$
d \left([ A [ B C ] _ {0} ] _ {0} - (- 1) ^ {| A | | B |} [ B [ A C ] _ {0} ] _ {0}\right) > n.
$$

When  $n = 0$  this is immediately clear.

For  $n = 1$  we have to consider the different orderings. Let us look first at the case where  $A = \partial X$ :

$$
\begin{array}{l} [ \partial X [ B C ] _ {0} ] _ {0} = [ B [ \partial X C ] _ {0} ] _ {0} - \sum_ {l \geq 1} \frac {(- 1) ^ {l}}{l !} [ \partial^ {l} [ \partial X B ] _ {l} C ] _ {0} \\ = \left[ B \left[ \partial X C \right] _ {0} \right] _ {0} + \sum_ {l \geq 2} \frac {(- 1) ^ {l} (l - 1)}{l !} \left[ \partial^ {l} [ X B ] _ {l - 1} C \right] _ {0}, \tag {3.3.15} \\ \end{array}
$$

where we used eq. (3.3.1) in the last step. We see that the derivative-number for the second term in the  $\text{rhs}$  is at least 2. The case where  $B = \partial X$  follows from eq. (3.3.15) by moving the correction terms to the  $\text{lhs}$ . Finally, when  $C = \partial X$ , the correction terms will be of the form  $[\partial^l [AB]_l\partial X]_0$ , which again have derivative number greater than 2 because  $l \geq 1$ .

It is now easy to see that the same reasoning can be used for general  $n$ . For this we will need that:

$$
d \left(\partial^ {l} \left[ \partial^ {m} A \partial^ {n} B \right] _ {l}\right) \geq m + n + 1,
$$

(or  $[\partial^m A\partial^n B]_l = 0)$  which follows from eqs. (3.3.1) and (3.3.2).

To conclude the proof, we note that with each application of eq. (3.3.10) the derivative number of the correction terms increases. On the other hand, the dimension of the correction terms remains constant. However, for any operator  $A$ , one has that:

$$
d (A) \leq \dim (A) - D,
$$

where  $D$  is the smallest dimension occurring in the OPA. This shows that after a finite number of applications of eq. (3.3.10), no correction terms can appear.

Theorem 3.3.2 guarantees that the algorithm to reorder composites ends in a finite number of steps. However, we have not yet proven that the end-result of the algorithm gives a "standard form", i.e. is the result unique if two expressions are the same up to null-relations?

# Intermezzo 3.3.2

Consider the following counterexample. Suppose there is a null-relation:

$$
[ A B ] _ {0} - C = 0,
$$

where  $A, B, C$  are bosonic operators. If we use this relation to eliminate the composite  $[AB]_0$ , we find immediately

$$
[ [ A B ] _ {0} D ] _ {0} = [ C D ] _ {0},.
$$

On the other hand, two applications of eq. (3.3.7) and one of eq. (3.3.10) show that

$$
[ [ A B ] _ {0} D ] _ {0} = [ A [ B D ] _ {0} ] _ {0} + \dots
$$

where the ellipsis denotes terms with derivatives. Together, these relations show that  $[A[BD]_0]_0 = [CD]_0 + \ldots$  while both  $lhs$  and  $rhs$  are not changed by the algorithm.

The intermezzo shows that null-relations should be used to eliminate the "simpler" fields in terms of the composites. We need the notion of "composite-number"  $c(A)$ . When we do not use any null-relations, we have:

$$
c (A) = 1 \quad \text {i f} A \text {i s a g e n e r a t o r}
$$

$$
\begin{array}{r c l} c (\partial A) & = & c (A) \end{array}
$$

$$
c ([ A B ] _ {0}) = c (A) + c (B)
$$

$$
c \left(\sum_ {i} A _ {i}\right) = \min  _ {\forall i \text {s u c h t h a t} d \left(A _ {i}\right) = d \left(\sum_ {j} A _ {j}\right)} c \left(A _ {i}\right). \tag {3.3.16}
$$

With these definitions, the composite-number is equal for both sides of all eqs. (3.3.5-3.3.10). We say that we "use a null-relation in increasing derivative- and composite-number" if we use it to eliminate (one of) the operator(s) with the smallest composite-number of those with the smallest derivative-number. Note that we have  $c(\partial A) \geq c(A)$  in this case. We can now assert:

Theorem 3.3.3 In an OPA satisfying the conditions of theorem 3.3.1, and where any null-relations are used in increasing derivative- and composite-number, the application of eqs. (3.3.5-3.3.11) (together with eqs. (3.3.1-3.3.4)) to order the operators in a normal ordered product defines a "standard form" on the elements of the algebra.

For the case of Poisson brackets, the rules (3.3.7) and (3.3.10) drastically simplify, only the first terms remain. Also, eq. (3.3.8) is changed to  $[AA]_0 = 0$  when  $A$  is fermionic.

# Improvements

The rules given in this section up to now are sufficient to compute any OPE, and to reorder any composite into a standard form. However, some shortcuts exist.  $[A[BC]_0]_q$  can be computed using eq. (3.3.4), but an alternative follows from eq. (2.3.22) using eq. (3.3.6):

$$
\begin{array}{l} [ A [ B C ] _ {0} ] _ {q} = (- 1) ^ {| A | | B |} \left([ B [ A C ] _ {q} ] _ {0} + \sum_ {l \geq 0} \frac {(- 1) ^ {l + q}}{l !} [ \partial^ {l} [ B A ] _ {l + q} C ] _ {0} \right. \\ \left. + \sum_ {l = 1} ^ {q - 1} (- 1) ^ {l} \left[ \left[ B A \right] _ {l} C \right] _ {q - l}\right). \tag {3.3.17} \\ \end{array}
$$

This rule is more convenient when we know the OPE  $B(z)A(w)$  while  $A(z)B(w)$  has to be computed using eq. (3.3.3).

Similarly, to compute an OPE where the first operator is a composite, the algorithm as presented above uses eqs. (3.3.3) and (3.3.4). Clearly, eq. (2.3.24) and eq. (3.3.6) implement this in one step:

$$
\begin{array}{l} [ [ A B ] _ {0} C ] _ {q} = \sum_ {l \geq 0} \frac {1}{l !} [ \partial^ {l} A [ B C ] _ {l + q} ] _ {0} + (- 1) ^ {| A | | B |} \sum_ {l \geq 0} \frac {1}{l !} [ \partial^ {l} B [ A C ] _ {l + q} ] _ {0} \\ + (- 1) ^ {| A | | B |} \sum_ {l = 1} ^ {q - 1} [ B [ A C ] _ {q - l} ], \tag {3.3.18} \\ \end{array}
$$

where  $q\geq 1$  and:

$$
\begin{array}{l} \left[ \left[ A B \right] _ {0} C \right] _ {0} = \left[ A \left[ B C \right] _ {0} \right] _ {0} + \\ \sum_ {l > 0} \frac {1}{l !} [ \partial^ {l} A [ B C ] _ {l} ] _ {0} + (- 1) ^ {| A | | B |} \sum_ {l > 0} \frac {1}{l !} [ \partial^ {l} B [ A C ] _ {l} ] _ {0}. \tag {3.3.19} \\ \end{array}
$$

We will see in subsection 3.3.2 where these rules are applied in OPErefs 3.1.

# 3.3.2 The internals of OPErefs

In this subsection, some crucial points in the implementation of the program are discussed. The rest of this chapter does not depend on the material presented here. Familiarity with the Mathematica programming language and packages is assumed. More information concerning these topics can be found in appendix C and [210]. For clarity, some of the rules in this subsection are slightly simplified versions of those appearing in the actual code of OPErefs. Input for Mathematica is written in typeset font.

# 3.3. Implementation

# 3.3.3 Operator handling

Before any calculations can be done, we have to distinguish between operators and scalars. The user of the package has to declare the bosonic and fermionic operators he/she want to use with the functions Bosonic and Fermionic e.g.

Bosonic[A,B]

This declaration simply sets the appropriate values for some internal functions: a function which distinguishes between operators and scalars (OperatorQ), and another for bosons and fermions (BosonQ).

$$
\begin{array}{l} \text {B o s o n i c H e l p [ A _ {-} ] : = (B o s o n Q [ A ] = T r u e ; O p e r a t o r Q [ A ] = T r u e ;} \\ \text {O P E p o s i t i o n [ A ] = O P E p o s i t i o n C o u n t e r + +)} \end{array}
$$

$$
\text {B o s o n i c} [ \mathrm {A} _ {- -} ] := \text {S c a n} [ \text {B o s o n i c H e l p}, \{\mathrm {A} \} ]
$$

$$
O P E \text {e p i s i o n C o u n t e r} = 0
$$

and similar rules for Fermionic. The function OPEposition is used to record the order between the operators, which is determined by the order in which they were declared and the standard Mathematica order for operators declared with the same pattern:

$$
\begin{array}{r l} \text {O P E O r d e r [ a _ {，} , b _ {]} : =} & \\ \text {B l o c k [ \{r e s = O P E p o s i t i o n [ b ] - O P E p o s i t i o n [ a ] \} ,} & \\ \text {I f [ r e s = = 0 , O r d e r [ a , b ] , r e s ]} & \\ ] & \end{array}
$$

We now list all rules that are necessary for testing if an expression is an operator or not.

$$
\operatorname {O p e r a t o r} Q [ A _ {-} + B _ {-} ] := \operatorname {O p e r a t o r} Q [ A ]
$$

$$
\operatorname {O p e r a t o r Q} [ \mathrm {A} _ {\text {T i m e s}} ] :=
$$

$$
\text {A p p l y} [ \text {O r}, \text {M a p} [ \text {O p e r a t o r Q}, \text {A p p l y} [ \text {L i s t}, A ] ]
$$

$$
\operatorname {O p e r a t o r Q} [ \text {D e r i v a t i v e} [ \_ ] [ A _ {\_} ] := \operatorname {O p e r a t o r Q} [ A ]
$$

$$
\text {O p e r a t o r Q [ \_ N O ]} = \text {T r u e}
$$

$$
\operatorname {O p e r a t o r} Q [ 0 ] = \text {T r u e}
$$

$$
\operatorname {O p e r a t o r Q} [ \_ ] = \text {F a l s e} \tag {3.3.20}
$$

The second rule is the most complicated one. It handles the testing of a product by testing if any of the factors is an operator. The fourth rule declares any composites (which we will give the head  $\mathsf{NO}$ ) to be operators. The last rule says that all expressions which were not handled by the previous rules (and the rules set by any declarations) are scalars.

Note that this simple way of defining OperatorQ is only possible because Mathematica reorders the rules such that more general rules are checked last. This means that when a user declares A to be bosonic, a OperatorQ[A] = True will be added on top of the list of rules given in (3.3.20).

The only rules we need for the function BosonQ are given below, together with the definition of SwapSign which gives the sign when interchanging two operators, i.e.  $-1$  for two fermions and  $+1$  otherwise:

$$
B o s o n Q [ \text {D e r i v a t i v e} [ \_ ] [ A _ {\_} ] := B o s o n Q [ A ]
$$

$$
\text {L i t e r a l} [ \text {B o s o n Q} [ \text {N O} [ \text {A} _ {,}, \text {B} _ {,} ] ] ] := \text {N o t} [ \text {X o r} [ \text {B o s o n Q} [ \text {A} ], \text {B o s o n Q} [ \text {B} ] ] ]
$$

$$
\operatorname {S w a p S i g n} [ A _ {-}, B _ {-} ] := \operatorname {I f} [ O r [ B o s o n Q [ A ], B o s o n Q [ B ] ], 1, - 1 ]
$$

The unit operator  $\mathbf{1}$  is declared in OPErefs, and is written as One. Introducing an explicit unit operator considerably simplifies many of the functions in OPErefs.

# Data representation

The representation we used in OPErefs 2.0 was a Laurent series representation. This allowed us to use the built-in rules of Mathematica to add OPEs, take derivatives with respect to one of the arguments, .... However, storage space turned out to be a problem for complicated algebras. Therefore, we now use a representation for the singular part of an OPE which is a list of the operators at the different poles. The highest order pole occurs first in the list, the first order pole is the last. As an example, the OPE of a Virasoro operator with itself, eq. (2.3.2), is represented as:

OPData[{c/2 One, 0, 2T, T}]

This representation has the additional advantage that it is more general than a Laurent series representation. Poisson brackets form an obvious example, see subsection 2.3.5.

To hide the details of the representation to the user, the user-interface for constructing an OPE is the function MakeOPE. This enables us to allow a different syntax for the interface and the representation, see section 3.4.

With the representation we use, it is easy to determine the order of the highest pole:

MaxPole[OPEData[A_List]] := Length[A]

For the Virasoro OPE this obviously gives 4.

The actual definitions for the OPedata structure enable us to multiply an OPE with a scalar, and add two OPEs. Also higher order poles that are zero are removed. The addition turns out to be the most complicated to achieve. One is tempted to use the Listable attribute of Plus, i.e. adding two list of the same length gives a list with the sum of the corresponding elements. To be able to use this feature, we should make the lists of poles the same length by extending the shorter list with extra zeroes (corresponding to zero operators at those poles). As addition of OPEs will occur frequently, the corresponding rule should be as efficient as possible:

```txt
OPEData /: n_ \* OPEData[A_List] := OPEData[n*A]  
OPEData /: A1_OPEData + A2__OPEData := Block[{maxP = Max[Map[MaxPole, {A1,A2}]]}}, OPEData[ Plus @@ Map[Join[Table[0,{maxP-Length[#]}], #]&, Map[First, {A1,A2}] ] ]  
]  
OPEData[{0, A__}]:= OPEData[{A}]
```

# Intermezzo 3.3.3

In this intermezzo we explain the rule for the addition of OPData structures. The lhs of the definition makes sure that  $\{\mathbf{A1},\mathbf{A2}\}$  will be the list of all OPDatas which are added<sup>5</sup>. When adding a number of OPData structures, maxP is set to the order of the highest pole occurring in the sum. Then First extracts the lists of poles, and we add the correct number of zeroes to the left. Finally, the sum of all these lists is made and OPData wrapped around the result. As an example, consider the sum of OPData[\{\mathbf{T},0\} ] and OPData[\{\mathbf{T}'\}], we get:

```txt
{A1,A2} -> {OPEData[{T,0}], OPEData[{T}]]}  
maxP -> 2  
First ... -> {{T,0}, {T'}}  
Join ... -> {{T,0}, {0,T'}}  
Plus @@ ... -> {T,T'}
```

# OPE rules

We use the head OPE to represent an OPE. The bilinearity of the OPE function is defined by:

```txt
Literal[OPE[A_,s_,B_,C_]] := s OPE[A,B,C] /; OperatorQ[B]  
Literal[OPE[a_,b_,Plus,c_]] := Distribute[ Lineartmp[a,b,c], Plus,Lineartmp, Plus,OPE ]
```

# Intermezzo 3.3.4

The rule for an OPE of a sum uses the Distribute[f[_,g,f,gn,fn] built-in function, which implements distributivity of f with respect to g, replacing g with gn and f with fn in the result. A more obvious rule would be:

$$
O P E [ a \_ \_ \_, b \_ + c \_ \_, d \_ \_ ] := O P E [ a, b, d ] + O P E [ a, c, d ] \tag {3.3.22}
$$

However, this rule is much slower in handling sums with more than two terms. Indeed when there is a sum of  $n_1$  ( $n_2$ ) terms in the first (second) argument of OPE, this rule would have to be applied  $n_1 \times n_2$  times, while Distribute handles all cases in one application. The difference in performance is shown in fig. 3.1. During the execution of

![](images/6f6992b5a07144b89951e12b361bc7b909b7ad52aa1f6e454a3080838468daa3.jpg)  
Figure 3.1: Timings for two different rules to implement bilinearity of OPEs. The timings are for OPE[A,B] with A of a sum of 10 terms and B a sum of  $i$  terms. The plain line is for Distribute (3.3.21) and the dashed line is for the sum-rule (3.3.22).

the algorithm, OPEs of sums occur frequently, hence it is important to select the most efficient rule.

The other rules that are needed for the OPE function are:

# 3.3. Implementation

```txt
Literal[OPE[Derivative[i_][A_,B_]]:=  
    OPEDerivativeHelpL[A,B,i]  
Literal[OPE[A_,Derivative[i_][B_]]]:=  
    OPEDerivativeHelpR[A,B,i]  
Literal[OPE[A_,NO[B_,C_]]]:=  
    CallAndSave[OPECompositeHelpR,A,B,C]  
Literal[OPE[NO[A_,B_,C_]]]:=  
    CallAndSave[OPECompositeHelpL,A,B,C] /;  
Not[SameQ[Head[B],NO]]  
Literal[OPE[B_,A_]]:=  
    OPECommuteHelp[B,A] /;  
Or[SameQ[Head[B],NO],OPEOrder[A,B]>0]  
Literal[OPE[_,_]]=OPEData{{}}
```

Here, OPEDerivativeHelpL, OPEDerivativeHelpR, OPECompositeHelpL, OPECommuteHelp implement eqs. (3.3.1), (3.3.2), (3.3.18) and (3.3.3) respectively. OPECompositeHelpR differentiates between (3.3.4) and (3.3.17). CallAndSave[f_,args__] is defined in OPErefs to call the function f with arguments args, saving the result when a global switch is set. The last rule declares all non-defined OPEs to be regular.

Using the Help functions has several advantages compared to inserting their definitions "in-place". First, the structure of the program becomes much clearer. Second, by simply redefining OPECompositeHelpL and OPECompositeHelpR, we can switch between OPEs and Poisson brackets.

The actual definition of the Help functions is of course rather technical, we will give one example in intermezzo 3.3.5. However, more important is to understand how the rules for OPE handle an OPE with arbitrarily complicated operators. First of all it is important to remember the evaluation sequence of Mathematica. The rules (3.3.23) are stored in the order they are defined, preceded by the definitions the user has given for the OPEs of the generators. Presented an OPE to evaluate, Mathematica starts on top of the list of rules and applies the first matching rule. Then the evaluation sequence is restarted. In this way, these rules obviously form an implementation of the algorithm presented in section 3.3.1.

Note that the rules (3.3.23) avoid making use of the more complicated function OPECompositeHelpL when the first argument of the OPE is a nested composite, instead OPECommuteHelp is used on the result of OPECompositeHelpR. This is more efficient in most cases.

Let us as an example follow OPEdefines in the computation of a fairly complicated OPE.

```latex
$\operatorname{In}[1] := \text{Bosonic}[\mathrm{A}, \mathrm{B}, \mathrm{C}, \mathrm{D}, \mathrm{E}]$
```

```txt
In[2] := TracePrint[OPE[NO[A, NO[C', E]], NO[B, D]], (OPEDerivativeHelpL|OPEDerivativeHelpR| OPECompositeHelpL|OPECompositeHelpR| OPECommuteHelp)[_]]
```

```txt
Out[2] = OPECompositeHelpR[NO[A, NO[C', E]], B, D]  
OPECommuteHelp[NO[A, NO[C', E]], B]  
OPECompositeHelpR[B, A, NO[C', E]]  
OPECommuteHelp[B, A]  
OPECompositeHelpR[B, C', E]  
OPEDerivativeHelpR[B, C, 1]  
OPECommuteHelp[NO[A, NO[C', E]], D]  
OPECompositeHelpR[D, A, NO[C', E]]  
OPECommuteHelp[D, A]  
OPECompositeHelpR[D, C', E]  
OPEDerivativeHelpR[D, C, 1]  
OPECommuteHelp[D, C]
```

# Intermezzo 3.3.5

In this intermezzo, we give an example of the implementation of a Help function, OPECompositeHelpRQ which corresponds to eq. (3.3.4).

```txt
CompositeHelpRQ[A_,B_,C_]:=  
Block{{q,l,sign = SwapSign[A,B], ABC, AB, AC, maxAB, maxABC, maxq},  
AB = OPE[A,B];  
AC = If[ SameQ[B,C], AB, OPE[A,C]];  
maxAB = MaxPole[AB];  
ABC = Table[  
OPE[OPEPole[q][AB], C],  
{q,maxAB}];  
maxABC = Map[MaxPole, ABC];  
maxq = Max[maxABC + Range[maxAB, MaxPole[AC]];  
maxABC = Max[maxABC,0];  
OPEData[  
Table[  
OPESimplify[  
sign * NO[B,OPEPole[q][AC]] +  
NO[OPEPole[q][AB], C] +  
Sum[Binomial[q-1,1] *  
OPEPole[l][ABC[[q-1]] ],  
{1,Max[1,q-maxAB, Min[q-1, maxABC]}]  
],  
{q,maxq,1,-1}  
]
```

Although this is rather lengthy, most of the lines are quite straightforward. The most important part of the routine is formed by the statements bracketed with OPData. Due to the internal representation for OPEs that is used, the result of the routine is just a table of all the poles in the OPE - starting with the highest pole - with OPData wrapped around it. Note that ABC is defined such that OPEPole[1][ABC[[q-1]]] is equal to  $[AB]_{q-l}C]_l$ .

The determination of the order  $\max q$  of the highest pole in the OPE requires some explanation. Clearly, the first term in the  $rhs$  of eq. (3.3.4) gives  $\max q \geq \max AC$ . Furthermore, the terms in the sum over 1 are zero unless

$$
1 \leq \operatorname {M a x P o l e} [ A B C [ [ q - 1 ] ] ]
$$

or also

$$
q - 1 \leq \operatorname {M a x P o l e} [ \mathrm {A B C} [ [ 1 ] ] ]
$$

hence

$$
q \leq \operatorname {M a x P o l e} [ \mathrm {A B C} [ [ 1 ] ] ] + 1
$$

which is exactly the other boundary on maxq which is used.

This shows that the boundaries on 1 and  $\mathbf{q}$  can be computed without any notion of the dimension of a field. The rule is thus suitable for any OPA, which is also true for all other Help functions.

# Rules for derivatives

The Mathematica symbol Derivative has almost no rules associated with it, except derivatives of the standard functions like Sin. However, we do not need these standard functions, but we do need linearity of derivatives  $-\partial (A + B) = \partial A + \partial B$  – which is not included in the built-in definitions for Derivative. We chose not to use a new symbol for derivatives because of the convenient notation A' for Derivative[1][A]. The rules are similar to declaring linearity of OPEs:

Derivative[i_][a_].  $\coloneqq$  Map[Derivative[i], a]

Derivative[i_][a_b_]:=

b Derivative[i][a] /;OperatorQ[a]

Derivative[_[0] = 0 (3.3.24)

The first rule handles the  $i$ -th derivative of a sum. As Derivative[i_] expects only one argument, a call to Distribute is not necessary. Instead, we map the  $i$ -th derivative on all terms in the sum.

The second rule handles operators multiplied with scalars. We never need products of operators in our framework.

# Intermezzo 3.3.6

We comment on the efficiency of the rule for a product given in (3.3.24). In Mathematica standard order,  $\mathsf{b}$  follows a. Because Times has the attribute Orderless, this means that during pattern matching the  $\mathsf{a}_{-}$  pattern will be matched sequentially to every factor in

the product, while  $\mathsf{b\_}$  will be the rest of the product. As we know that the argument of Derivative[i] consists of some scalars times only one operator, we see that for each scalar  $\mathsf{OperatorQ}$  is called once, until the operator is found. Hence, it takes maximum  $n$  evaluations of  $\mathsf{OperatorQ}$  for a product of  $n$  factors. The steps in evaluating (a b c), where only one of these is an operator, are:

```lua
ifOperatorQ[a] then return(bca') elseifOperatorQ[b] then return(acb') elseifOperatorQ[c] then return(abc')endif
```

This has to be contrasted with the rule

Derivative[i__][b_a_]:= a Derivative[i][b] /;OperatorQ[b] (3.3.25)

As a is ordered before b the steps in evaluating the same derivative (a b c)’ are:

```txt
ifOperatorQ[bc] then return(Times[a, ifOperatorQ[c] then return(bc') elseifOperatorQ[b] then return(cb')endif ]) elseifOperatorQ[a c] then ... elseifOperatorQ[a b] then ...endif
```

This is clearly more complicated. Moreover, to test  $\mathsf{OperatorQ[b c]}$  may require two evaluations of  $\mathsf{OperatorQ}$ . In general, for  $n$  factors the worst case would be that  $n(n - 1)/2$  tests are needed. The difference in performance is rather drastic, see fig. 3.2 where we include also two other possibilities using Not[OperatorQ[a]].

# Other rules in OPErefs

With the above rules we have already a working package to compute OPEs of arbitrary nested composites. We also need to bring composites in a standard form. This is done by NO. The rules attached to NO are entirely analogous to those for OPE.

One also needs some extra functions like OPESimplify. We do not give their implementation here. One can always consult the source of OPEdfs for further information.

One major ingredient of OPErefs is not yet discussed: OPEPole. For the operation of the above rules to work, only a very simple definition for the OPEPole of an OPEData is needed:

![](images/3a2a72391c573cf00d4aef21f089272d96b88947f887b6436b88d9b0a432ae7f.jpg)  
Figure 3.2: Timings for four different rules to implement extracting of scalars as a function of the number of scalars. The timings are for the evaluation of  $\mathsf{d}[\mathsf{A}\mathsf{B}]$  where  $\mathsf{A}$  is a product of  $i$  scalars and  $\mathsf{B}$  is an operator.

Plain line for d[a_b_]:= bd[a] /; 0p[a], dashed line for d[a_b_]:= bd[a] /; !0p[b], dash-dot line for d[a_b_]:= ad[b] /; !0p[a], dotted line for d[a_b_]:= ad[b] /; 0p[b].

$$
\begin{array}{l} \text {O P E P o l e [ n ] [ O P E D a t a [ A ] ] : =} \\ \text {I f [ 1 <   = n <   = L e n g t h [ A ] , A [ [ - n ] ] , 0 ]} \end{array}
$$

OPErefs allows also a syntax OPEPole[i][A,B]. This can be used to compute only one pole of the complete OPE, avoiding the computation of all the other poles. Currently, the rules for OPEPole are completely independent of those to compute a full OPE. This is because when computing e.g. OPE[A, NO[B, C]] (see OPECompositeHelpR above), it is more efficient to compute OPE[A, B] and OPE[A, C] once and save these results. It is not a good programming practice to keep two separate sets of rules computing essentially the same thing, but efficiency seemed to be more important at this point.

# 3.3.4 Performance

In [192], a Wakimoto [202] realisation for the Kac-Moody algebra  $\widehat{B_2}$  with level  $k$  using 2 free scalars and 4 (bosonic)  $\beta$ ,  $\gamma$  systems was constructed. Describing this realisation would lead us too far here, but to give an idea of the complexity of the calculation, the total number of composites in the realisation is roughly sixty, and composites of up to four free fields are used.

In table 3.1, we tabulate CPU times for computing an OPE of two of the currents, and the Sugawara tensor for this realisation. The first time given in the table is

the time for evaluating the statement after loading the package and defining the realisation. The time between brackets is measured when the statement is repeated. The second execution is much faster because OPEdfs stores some of the OPEs with composites. Note that version 2.0 of Mathematica is roughly 1.4 times slower than version 1.2.

Table 3.1: CPU time for the computation of the OPE of the currents corresponding to the positive simple root of  $\widehat{B_2}$  (statement 9) and the computation of the Sugawara tensor (statement 11) (see Ref. [192]) for Mathematica running on a PC 386 (25 Mhz).  

<table><tr><td>Mathematica-version
OPErefs-version</td><td>1.2
2.0</td><td>1.2
3.1</td><td>2.0
3.1</td></tr><tr><td>In[9]
In[11]</td><td>23.5 (4.5) s
43.2 (11.6) s</td><td>14.9 (2.8) s
31.3 (9.4) s</td><td>19.3 (3.8) s
40.7 (12.1) s</td></tr></table>

# 3.4 User's Guide

This section is intended as a user's guide to the package OPErefs 3.1. Explicit examples are given for most operations. Note that OPErefs 3.1 requires Mathematica 1.2 or later.

We introduce some special notations. Input for and output from Mathematica is written in typeset font. Input lines are preceded by “ $In[n] :=$ ”, and corresponding output statements by “ $Out[n] =$ ”, as in Mathematica.

As OPErefs is implemented as a Mathematica package, it has to be loaded before any of its global symbols is used. Loading the package a second time will clear all previous definitions of operators and OPEs, as well as all stored intermediate results. Assuming that the package is located in the Mathematica-path, e.g. in your current directory, issue:

$$
I n [ 1 ] := \quad <   <   O P E d e f s. m
$$

After loading OPErefs into Mathematica, help for all the global symbols is provided using the standard help-mechanism, e.g. ?OPE.

Now, you need to declare the operators that will be used. If you want to define bosonic operators  $\mathbf{T}$  and  $\mathbf{J}[\mathbf{i}]$  (any index could be used), and fermionic operators  $\mathbf{psi}[\mathbf{i}]$ , the corresponding statements are:

$$
\begin{array}{l l} I n [ 2 ] := & \text {B o s o n i c} [ T, J [ i _ {-} ] ] \\ I n [ 3 ] := & \text {F e r m i o n i c} [ p s i [ i _ {-} ] ] \end{array}
$$

The order of the declarations fixes also the ordering of operators used by the program:

$$
T <   J [ 1 ] ^ {\prime} <   J [ 1 ] <   J [ 2 ] <   J [ i ] <   p s i [ 1 ] <   \dots \tag {3.4.1}
$$

By default, derivatives of an operator are considered "smaller" than the operator itself. This can be reversed using the global options NOOrdering (see below).

Finally, the non-regular OPEs between the basic operators have to be given. An OPE can be specified in two different ways.

The first way is by listing the operators that occur at the poles, the first operator in the list is the one at the highest non-zero pole, the last operator has to be the one at the first order pole, e.g.:

$$
I n [ 4 ] := \quad O P E [ T, T ] = M a k e O P E [ \{c / 2 O n e, 0, 2 T, T ^ {\prime} \} ];
$$

Note the operator One which specifies the unit-operator.

The second way is by giving the OPE as a Laurent series expansion, adding the symbol  $\mathsf{Ord}$  which specifies the (implicit) arguments of the operators for which the OPE is defined<sup>6</sup>. The arguments for the operators can be any Mathematica expression.

Warning: it is important that the operators occurring as arguments of OPE in a definition should be given in standard order (3.4.1), otherwise wrong results will be generated.

The following statements define a  $\widehat{SU(2)}_k$ -Kac-Moody algebra:

$$
\begin{array}{r l} I n [ 5 ] := & O P E [ J [ i _ {-} ], J [ i _ {-} ] ] := \\ & \text {M a k e O P E} [ - k / 2 (z - w) ^ {\wedge} - 2 + 0 r d [ z, w, 0 ] ] \end{array}
$$

$$
\begin{array}{r l} \operatorname {I n} [ 6 ] := & \operatorname {O P E} [ J [ 1 ], J [ 2 ] ] = \\ & \text {M a k e O P E} [ J [ 3 ] [ w ] (z - w) ^ {\wedge} - 1 + O r d [ z, w, 0 ] ]; \end{array}
$$

$$
\begin{array}{r l} I n [ 7 ] := & O P E [ J [ 2 ], J [ 3 ] ] = \\ & \text {M a k e O P E} [ J [ 1 ] [ w ] (z - w) ^ {-} - 1 + O r d [ z, w, 0 ] ]; \end{array}
$$

$$
\begin{array}{r l} I n [ 8 ] := & O P E [ J [ 1 ], J [ 3 ] ] = \\ & \text {M a k e O P E} [ - J [ 2 ] [ w ] (z - w) ^ {\wedge} - 1 + O r d [ z, w, 0 ] ]; \end{array}
$$

In fact, with the above definitions, one has to use always the explicit indices  $1,2,3$  for the currents  $J$ . If we would compute an OPE with current  $\mathbf{J}[\mathbf{i}]$  where the index  $i$  is not 1,2 or 3, wrong results will be given. One can circumvent this peculiarity by reformulating the definitions.

A normal ordered product  $[AB]_0$  is entered in the form  $\mathsf{NO}[\mathsf{A},\mathsf{B}]$ . Multiple composites can be entered using only one NO head, e.g.  $\mathsf{NO}[\mathsf{A},\mathsf{B},\mathsf{C}]$ . This input is effectively translated into  $\mathsf{NO}[\mathsf{A},\mathsf{NO}[\mathsf{B},\mathsf{C}]]$ . All output is normal ordered with the same convention, i.e. from right to left (input can be in any order). Also, the operators in

composites will always be ordered according to the standard order (3.4.1).

As an example, we can define the Sugawara energy-momentum tensor for  $\widehat{SU(2)}_k$ . The Mathematica output of an OPE is a list of the operators at the poles.

$$
\begin{array}{l} I n [ 9 ] := \quad \mathrm {T s} = - 1 / (\mathrm {k} + 2) (\mathrm {N O} [ \mathrm {J} [ 1 ], \mathrm {J} [ 1 ]) + \\ \left. \right.\left. \right.\left. \right.\left. \right.\left. \right.\left. \right.\left. \right.\left. \right.\left. \right.\left. \right.\left. \right.\left. \right.\left. \right.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left.\left. \text {J} [ 2 ], J [ 2 ]\right] + \text {N O} [ J [ 3 ], J [ 3 ] ]\right]\right)\right)\right)\right)\right)\right)\right)\right)\right)\right)\right)\right)\right)\right)\right)\right)\right)\right)\right) \\ I n [ 1 0 ] := \text {O P E S i m p l i f y} [ O P E [ T s, J [ 1 ] ] ] \\ O u t [ 1 0 ] = \ll 2 | | J [ 1 ] | | 1 | | J [ 1 ], > > \\ \end{array}
$$

Warning: when computing OPEs with composites, or when reordering composites, OPEdefines remembers by default some intermediate results. Thus, it is dangerous to change the definition of the basic OPEs after some calculations have been performed. For example, consider a constant  $a$  in an OPE. If calculations are performed after assigning a value to  $a$ , the intermediate results are stored with this value. Changing  $a$  afterwards will give wrong results.

The other globally defined functions available from the package are:

- OPEOperator[operator_, parity_] provides a more general way to declare an operator than Bosonic and Fermionic. The second argument is the parity of the operator such that  $(-1)^{\text{parity}}$  is  $+1$  for a boson, and  $-1$  for a fermion. It can be a symbolic constant. This is mainly useful for declaring a bc-system of unspecified parity, or a Kac-Moody algebra based on a super-Lie algebra. In such cases, the operator can contain a named pattern:

$$
I n [ 1 1 ] := \quad O P E O p e r a t o r [ J [ i _ {-} ], p a r i t y [ i ] ]
$$

If one wants to declare more operators, one can group each operator and its parity in a list:

$$
I n [ 1 2 ] := O P E O p e r a t o r \left[ \{b [ i _ {-} ], \text {p a r i t y} [ i ] \}, \{c [ i _ {-} ], \text {p a r i t y} [ i ] \} \right]
$$

See also SetOPE0ptions[ParityMethod, _].

- OPEPole[n_][ope_] gets a single pole term of an OPE:

$$
I n [ 1 3 ] := \quad O P E P o l e [ 2 ] [ O u t [ 1 0 ] ]
$$

$$
O u t [ 1 3 ] = J [ 1 ]
$$

OPEPole[n_][A_,B_] can also be used to compute only one pole term of an OPE:

$$
I n [ 1 4 ] := \text {F a c t o r} [ O P E P o l e [ 4 ] [ T s, T s ] ]
$$

$$
O u t [ 1 4 ] = (3 k 0 n e) / (2 (2 + k))
$$

OPEPoIe can also give terms in the regular part of the OPE:

$$
I n [ 1 5 ] := O P E P o l e [ - 1 ] [ T, T ]
$$

$$
O u t [ 1 5 ] = N O [ T ^ {\prime}, T ]
$$

- MaxPole[ope_] gives the order of the highest pole in the OPE.  
- OPEParity[A] returns an even (odd) integer of  $A$  is bosonic (fermionic).

# 3.5. Example : The conformal anomaly in superstring theory

- OPESimplify[ope_, function_] "collects" all terms in ope with the same operator and applies function on the coefficients. When no second argument is given, the coefficients are Expanded.

$$
\begin{array}{r l} I n [ 1 6 ] := & O P E S i m p l i f y [ O P E [ J [ 1 ], N O [ J [ 2 ], J [ 1 ] ] ] ] \\ O u t [ 1 6 ] = & \ll 2 | | (1 - k / 2) J [ 2 ] | | 1 | | \\ & N O [ J [ 1 ], J [ 3 ] ] + J [ 2 ] ^ {\prime} \gg \end{array}
$$

OPESimplify[pole_, function_] does the same simplifications on sums of operators.

- OPEMap[function_, ope_] maps function to all poles of ope.  
- GetCoefficients [expr_] returns a list of all coefficients of operators in expr which can be (a list of) OPEs or poles.  
- OPEJacobi[op1_, op2_, op3_] computes the Jacobi-identities (2.3.21) for the singular part of the OPEs of the three arguments. Due to the nature of eq. (2.3.21), the computing time will be smallest (in most cases) when op1 ≤ op2 ≤ op3 in the order (3.4.1). In nonlinear algebras, the computation will use rules for OPEs of composites which assume that the Jacobi identities hold. This means that the result of OPEJacobi gives only necessary conditions. In some cases, different orderings of op1-op3 have to be tried to find all conditions.

The result of OPEJacobi is a double list of operators. It is generated by Table [OPEPole[n] [A, OPEPole[m] [B, C]] + corrections, {m, maxm}, {n, maxn}]

All elements of the list should be zero up to null operators for the OPA to be associative.

- Delta[i_,j_] is the Kronecker delta symbol  $\delta_{ij}$ .

- ClearOPESavedValues[] clears all stored intermediate results, but not the definition of the operators and their OPEs. To clear everything, reload the package.  
- OPETO Series [ope_] converts an OPE to a Laurent series expansion in  $\mathbf{z}$  and  $\mathbf{w}$ . The arguments can be set to  $\mathbf{x}$  and  $\mathbf{y}$  with:

$$
I n [ 1 7 ] := \quad \text {S e t O P E O p t i o n s} [ \text {S e r i e s A r g u m e n t s}, \{x, y \} ]
$$

- TeXForm[ope_] gives TExoutput for an OPE. The same arguments are used as in OPETOseries.  
- OPESave[filename_] (with filename a string between double quotes) saves the intermediate results that OPEdfs remembers to file (see the option OPESaving below).

- SetOPE0options is a function to set the global options of the package. The current options are:

- SetOPEOptions[SeriesArguments, {arg1_, arg2_}]: sets arguments to be used by TeXForm and OPEToSeries. One can use any Mathematical expression for arg1 and arg2.  
- SetOPEOptions[NOOrdering, n_]: if n is negative, order higher derivatives to the left (default), if n is positive, order them to the right.  
- SetOPEOptions[ParityMethod, 0|1]: makes it possible to use operators of an unspecified parity. When the second argument is 0 (default), all operators have to be declared to be bosonic or fermionic. When the argument is 1, OPEOperator can be used with a symbolic parity. Note that in this case, powers of  $-1$  are used to compute signs, which is slightly slower than the boolean function which is used by the first method.

This option is not normally needed as the use of OPEOperator with a non-integer second argument sets this option automatically.

- SetOPEOptions[OPESaving, boolean_]: if boolean evaluates to True (default), OPErefs stores the intermediate results when computing OPEs of composites and when reordering composites. This option is useful if Mathematica runs short of memory in a large calculation, or when computing with dummy indices<sup>7</sup>.  
- SetOPEOptions[OPEMethod, method_]: with the parameter method equal to QuantumOPEs enables normal OPE computations (default setting), while ClassicalOPEs enables Poisson bracket computations. Using this option implicitly calls ClearOPESavedValues[].

# 3.5 Example : The conformal anomaly in superstring theory

We consider only one free boson field  $X$  and one free fermion field  $\psi$  because additional free fields will have exactly the same OPEs and commute with each other. We denote  $\partial X$  with  $\mathbf{J}$  and  $\psi$  with  $\mathbf{psi}$  (we normalise them such that they have a  $+1$  in their OPEs (2.6.7) and (2.6.26)). The ghosts are a fermionic  $b$ ,  $c$  system (operators  $\mathbf{b},\mathbf{c}$ ) and a bosonic  $\beta$ ,  $\gamma$  system (operators  $\mathbf{B},\mathbf{G}$ ) (normalised such that  $iA^{j} = (-1)^{i + 1}{}_{i}\delta^{j}$  in eq. (2.6.31)).  $b$  has conformal dimension 2 and  $\beta$  has  $3/2$ . It is now a trivial task to compute the conformal anomaly:

$$
I n [ 1 ] := \quad <   <   O P E d e f s. m
$$

```latex
$\begin{array}{rl}\text{In[2]}\coloneqq & \text{Bosonic[J,B,G];Fermionic[b,c,psi]};\\ & \text{OPE[J,J] = MakeOPE[\{One,0\}]};\\ & \text{OPE[psi,psi] = MakeOPE[\{One\}]};\\ & \text{OPE[b,c] = MakeOPE[\{One\}]};\\ & \text{OPE[B,G] = MakeOPE[\{One\}]};\\ & \text{Tb = 1 / 2 NO[J,J];Tf = -1 / 2 NO[psi,psi]};\\ & \text{Tbc = -2 NO[b,c]} -\text{NO[b',c]};\\ & \text{TBG = 3 / 2 NO[B,G]} +1 / 2\text{NO[B',G]}; \end{array}$
```

```latex
$\operatorname{In}[3] := \quad \text{OPESimplify}[\text{OPE}[\text{Tb}, \text{Tb}]]$
```

```txt
Out[3] = << 4||One/2 ||3||0 ||2||NO[J, J] ||1|| NO[J', J] >>
```

```latex
$\operatorname{In}[4] := \quad \text{OPESimplify}[\text{OPE[Tf,Tf}]$
```

```txt
Out[4] = << 4|| One/4 ||3|| 0 ||2|| NO[psi', psi] ||1|| NO[psi'], psi]/2>>
```

```latex
$\operatorname{In}[5] := \quad$  OPESimplify[OPE[Tbc,Tbc]]
```

```javascript
Out[5] = << 4 || -13 0ne || 3 || 0 || 2 ||
```

```csv
-4 NO[b, c'] - 2 NO[b', c]||1||
```

```txt
-2 NO[b, c'] - 3 NO[b', c'] - NO[b'', c] >>
```

```txt
$\operatorname{In}[6] := \quad$  OPESimplify[OPE[TBG, TBG] - MakeOPE[{2 TBG, TBG}]
```

```javascript
Out[6] = << 4|| 11 0ne/2 ||3|| 0 ||2|| 0 ||1|| 0 >>
```

We see that each bosonic (fermionic) field will contribute a central charge  $1(1/2)$  to the total central charge of the theory. The  $b$ ,  $c$  system contributes  $-26$ , and the  $\beta$ ,  $\gamma$  system 11. This gives the well-known relation for the critical dimensions of the bosonic string  $D_{b} - 26 = 0$  and the superstring  $3/2D_{s} - 26 + 11 = 0$ . Moreover, we can easily verify that the energy-momentum tensors obey the Virasoro algebra.

The reader without experience in CFT is invited at this point to take out some time and compute the OPE for  $T_{BG}$ , for instance, by hand. Although this computation is rather trivial with OPErefs, the same calculation was attempted in [179] using the mode-algebra. There it proved not to be possible to compute the Virasoro algebra automatically due to difficulties with the infinite sums in the normal ordered products.

# 3.6 Future developments

A first extension would be to add the possibility to specify a range of poles one wants. This would unify the current implementations of OPE and OPEPole. The main difficulty in programming this, is that when saving intermediate results, one has to see which poles are already computed and which not.

The objective in writing OPErefs, was to make a package available which is as general as possible. One could write specialised extensions which would outperform

OPErefs. For instance, the restriction to free fields would be very useful. Also  $\mathcal{W}$ -algebras, a conformal OPA generated by quasiprimaries, form a preferred subclass of the OPAs. In this case it would be advantageous to compute in a basis of quasiprimaries, see section 4.3. A package under development [194] collects some formulas for working with (quasi)primaries in a  $\mathcal{W}$ -algebra, but it relies on OPErefs for computing OPEs.

The main restriction of OPErefs is of course the requirement of poles of integer order. In particular, vertex operators are widely used in conformal field theory. Here the order of the powers in the generalised Laurent expansion remain integer spaced. This should make it possible to extend OPErefs to handle this case. However, the notion of the singular and regular part of an OPE is not so important when using vertex operators as in other cases. Indeed, the normal ordered product of two vertex operators should not be defined as the zeroth order pole in their OPE. This makes it desirable to use a data representation which keeps the information to generate any "pole", but stores already computed results<sup>8</sup>. This would make it possible to work with vertex operators  $V_{a}$  with a symbolic weight  $a$ , and not only with fixed numbers. One could then also define an OPE between fields of (non-numeric) dimension  $h_1, h_2$ , e.g. with the unit operator at the pole of order  $h_1 + h_2$ . Finally, the extension to nonmeromorphic OPEs would allow to treat parafermions [67, 68, 95].

We are currently working at a version for computing with super OPEs in  $N = 2$  superfields [137]. Extension to arbitrary  $N$  will not give great difficulties.

Finally, it would be convenient to be able to use dummy indices in OPErefs. I wrote a separate package Dummies to handle this. However, this package uses only a rudimentary algorithm to simplify expressions with dummy indices. There seems to exist no algorithm (except exhaustive enumeration) to do this simplification for the case at hand. The main difficulty is that correction terms are needed when interchanging operators in a normal ordered product.

# 3.7 Other packages

As shown in section 2.4, an other approach to the same problem would be to use modes. An attempt to compute commutators of the modes of normal ordered operators in REDUCE [179] was not completely successful due to difficulties in the reshuffling of the indices in infinite sums. It should be possible to avoid this by constructing the formulas for commutation with a mode of a composite operator by looking at the corresponding formula for OPEs. Apparently, this is done in a package by L. Romans, which is not published, but acknowledged by a few authors. The modeapproach has probably the advantage that contributions of the tower of derivatives of a quasiprimary operator (see section 4.3) are summed in advance. However, this

# 3.7. Other packages

definition (2.4.1) which makes the formulas for computing with modes convenient is restricted to a conformal OPA.

Related to mode calculations is the approach taken by H. Kausch (using another symbolic manipulation program Maple). He uses the equivalence between conformal fields  $\Phi(x)$  and states  $\Phi(0)|0 >$  discussed in [98]. As this reference restricts itself to bosonic fields of integer dimension, and fermionic fields of half-integer dimension, we suspect that Kausch's program will not be able to handle algebras outside this framework, like [160, 21, 49], while they present no problems for OPErefs. However, Kausch did not publish his work.

Finally, A. Fujitsu recently developed a package, ope.math, in Mathematica for computing OPEs of free fields [93]. It is able to treat vertex operators at non-integer poles (this is not possible in  $OPE_$  defs). However, one cannot compute OPEs in a  $\mathcal{W}$ -algebra or a Kac-Moody algebra, except by working in a realisation with free fields. The current version is much slower than  $OPE_$  defs for calculations without vertex operators. ope.math uses a fixed notation for all fields, e.g. the  $i$ -th derivative of a complex boson (part of a bosonic  $b$ ,  $c$ -system, see section 2.6) has to be called be[i, j, z] where  $j$  is an index and  $z$  the coordinate. The package is currently not able to work with an unspecified number of fields, i.e. ope[be[0,i, z], ga[0,j, w]] returns zero.

# Chapter 4

# $\mathcal{W}$ -algebras

In this chapter, we focus our attention on  $\mathcal{W}$ -algebras. Recall that we defined a  $\mathcal{W}$ -algebra as an Operator Product Algebra (OPA) generated by a Virasoro operator  $T$  and quasiprimary operators  $W^i$ . We first discuss the representations of  $\mathcal{W}$ -algebras, section 4.2. We then show how the representations of the global conformal group restrict the possible Operator Product Expansions (OPEs) of the generators, section 4.3. In the next section, the full conformal group is studied. We demonstrate by using the Jacobi identities how to reconstruct the complete OPEs from the knowledge of the coefficients of the primary operators in these OPEs. The ultimate goal of section 4.3 is to provide the formulas for an algorithmic computation of Virasoro descendants and the coefficients with which they appear in OPEs.

In section 4.5, we discuss a few of the more important methods to construct  $\mathcal{W}$ -algebras and comment on the classification of the  $\mathcal{W}$ -algebras. As an illustration of the ideas in this chapter, we conclude with an example in section 4.6: the  $\mathcal{W}_cB_2$  algebra.

The presentation used in sections 4.3 and 4.4 is slightly more general than the literature, and in subsections 4.3.2, 4.4.2 and 4.4.3 new developments are given. The results on the  $\mathcal{W}_cB_2$  algebra were published in [78].

# 4.1 Introduction

The chiral symmetry generators of a conformal field theory in general form a  $\mathcal{W}$ -algebra. The representation theory of  $\mathcal{W}$ -algebras will thus have its direct consequences for the correlation functions of the theory. For physical applications, highest weight representations are most useful. In certain special cases, the correlation functions are completely fixed by the symmetry algebra. These representations are the minimal models, pioneered for the Virasoro algebra in [13] and extended to other  $\mathcal{W}$ -algebras in [66, 63, 65, 6]. They consist of a finite number of highest weight operators which each give rise to degenerate a representation. Minimal models were first studied as

a simple model for infinite dimensional highest weight representations. Almost all exactly solvable conformal field theories (CFT) are based on minimal models. They also appear in  $\mathcal{W}$ -string theory, see chapter 8.

The study of the  $\mathcal{W}$ -algebras is complicated by the fact that the OPEs can contain composite terms, or equivalently the mode algebra is nonlinear. The prototype of a nonlinear  $\mathcal{W}$ -algebra is the  $\mathcal{W}_3$  algebra [211]. It is generated by the energy-momentum tensor  $T$  and a dimension 3 current  $W$  with OPEs given by:

$$
T (z) T (w) = \frac {c}{2} (z - w) ^ {- 4} + 2 (z - w) ^ {- 2} T (w) + (z - w) ^ {- 1} \partial T (w) + \dots
$$

$$
T (z) W (w) = 3 (z - w) ^ {- 2} W (w) + (z - w) ^ {- 1} \partial W (w) + \dots
$$

$$
\begin{array}{l} {W (z) W (w)} = {\frac {c}{3} (z - w) ^ {- 6} + 2 (z - w) ^ {- 4} T (w) + (z - w) ^ {- 3} \partial T (w)} \\ + (z - w) ^ {- 2} \left[ 2 \beta \Lambda (w) + \frac {3}{1 0} \partial^ {2} T (w) \right] \\ + (z - w) ^ {- 1} \left[ \beta \partial \Lambda (w) + \frac {1}{1 5} \partial^ {3} T (w) \right] + \dots , \tag {4.1.1} \\ \end{array}
$$

where:

$$
\Lambda = [ T T ] _ {0} - \frac {3}{1 0} \partial^ {2} T \tag {4.1.2}
$$

and:

$$
\beta = \frac {1 6}{2 2 + 5 c}. \tag {4.1.3}
$$

The operator  $\Lambda$  is defined such that it is quasiprimary, i.e.  $[T\Lambda]_3 = 0$ . Although this is the simplest nonlinear  $\mathcal{W}$ -algebra for generic  $c$ , checking the Jacobi identities by hand is already a nontrivial task.

We notice that the  $\text{rhs}$  of the OPE  $[WW]$  can be written in terms of quasiprimaries or their derivatives. Using the global conformal transformations, generated by  $L_{\pm 1}, L_0$ , we will see in section 4.3 that the coefficients of the derivatives of the quasiprimaries, are numbers depending only on the dimensions of the operators involved. This result was already contained in the paper of Belavin, Polyakov and Zamolodchikov [13]. We present some new results on how the quasiprimaries in a given OPE can be found.

The full conformal group gives even more information. In case that the OPA is generated by primary operators and their descendants, the Jacobi identities with  $\{T,\Psi_i,\Psi_j\}$  fix the form of the OPE of two primaries  $\Psi_{i},\Psi_{j}$ . In fact, the coefficients of the descendants of the primaries occurring in the OPE now depend on the dimensions of the operators and the central charge. In the example of  $\mathcal{W}_3$ , the only primary in the singular part of the  $W(z)W(w)$  OPE is the unit operator at the sixth order pole. The rest of this OPE is then completely determined by the Jacobi identity  $\{T,W,W\}$ .

It was shown in [13, 74] how the coefficients of the Virasoro descendants can be found. As they can be very complicated, automation of the computation of these coefficients is highly desirable<sup>1</sup>. However, in the basis for the descendants used in [13, 74] this computation is very CPU-intensive. We will therefore construct a basis of quasiprimaries, and give the necessary formulas to compute descendants and coefficients in this basis. A Mathematica package that implements these algorithms is now in testing phase [194].

The conclusion will be that we can reconstruct the OPEs of the primary generators of a  $\mathcal{W}$ -algebra from the list of the coefficients of all primaries. For primary operators  $\Psi_{i}$ , we will write symbolically:

$$
\Psi_ {i} \times \Psi_ {j} \longrightarrow C _ {i j} ^ {k} [ \Psi_ {k} ], \tag {4.1.4}
$$

where  $[\Psi_k]$  denotes the conformal family of the primary operator  $\Psi_{k}$ . We only have to keep track of the primaries in the singular part of the OPE, because any primaries in the regular part can be constructed from the information contained in the singular part of the OPE. The structure constants  $C_{ij}^{k}$  are still restricted by the Jacobi identities discussed in 2.3, but now for triples of primaries. Unfortunately, the equations for the coefficients are too complicated to solve in the general case, preventing to classify the  $\mathcal{W}$ -algebras in this way. We will discuss briefly some other attempts towards the classification.

# 4.2 Highest weight representations and minimal models

In this section we consider representations of a  $\mathcal{W}$ -algebra which are constructed by acting repeatedly with the generators of the  $\mathcal{W}$ -algebra on a highest weight state  $|\phi\rangle$ .

Definition 4.2.1 For a mode algebra with generators  $W_{n}^{i}$  (with  $W^0 \equiv T$ ), a highest weight state (HWS)  $|\phi\rangle$  with weights  $w^{i}$  satisfies:

$$
W ^ {i} _ {0} | \phi > = w ^ {i} | \phi >, \quad W ^ {i} _ {n} | \phi > = 0 n > 0,
$$

In this definition it is understood that a generator with half-integral modes has no weight associated to it. We will use the correspondence between states and fields in meromorphic conformal field theory, as discussed in [98]. So, with every highest weight operator  $\phi$  corresponds a highest weight state  $|\phi\rangle$ .

For the Virasoro algebra it is natural to take  $\phi$  to be a primary operator of dimension  $h$ . This gives an extra condition  $\widehat{L}_{-1}\phi = \partial \phi$ . The notion of a primary operator, and in particular the requirement that the first order pole of the OPE of the energy-momentum tensor with a primary operator is the derivative of the operator, arises from the geometrical interpretation of a conformal transformation eq. (2.1.4). As no geometrical meaning is currently known for the transformations generated by currents with a nonlinear OPE, the concept of a  $\mathcal{W}$ -primary cannot be defined at present.

To define a highest weight representation, we first introduce some notation. We will denote a sequence  $(W^{i})_{-n_{1}}\ldots (W^{i})_{-n_{k}}$  where  $n_j\geq n_{j + 1} > 0$ , as  $(W^{i})_{-\{n\}}$ . This notation is extended to sequences of modes of different generators  $W^0_{-\{n_0\}}\dots W^d_{-\{n_d\}}$  which we write as  $W_{-\{n\}}$ , where  $\{n\}$  forms an ordered partition of  $N$  with a different "colour" for every generator. We use the convention that for a sequence of positive modes, the order is the reverse as for negative modes, i.e.  $W_{\{n\}}$  acts like  $W^{d}_{\{n_{k,d}^{d}\}}\ldots W^{d}_{\{n_{1}^{d}\}}\ldots W^{0}_{\{n_{1}^{0}\}}$ .

A Verma module is then defined as the space of all "descendants"  $W_{-\{n\}}|\phi >$ . The level of the descendant is defined as its  $L_0$  weight minus the  $L_0$  weight of  $|\phi >$ , e.g. the level of  $L_{-n}(W^i)_{-m}\phi$  is  $n + m$ . Although the dimension of the Verma modulo is infinite, the dimension of the space of descendants at a certain level is finite.

The representation can then be computed by using the (graded) commutators of the modes and the definition 4.2.1. The representation will (at least) depend on the weights  $w^i$  and on the central charge  $c$  of the  $\mathcal{W}$ -algebra, we will write  $\mathcal{R}(w^i, c)$ .

The representation is reducible, or degenerate, if at a certain level  $N$  there occurs a new HWS, which will be the starting point for a new tower of states. One can then divide out this new representation  $\mathcal{R}(\tilde{w}^i, c)$ , where  $\tilde{h} = \tilde{w}^0 = h + N$ . If this process is repeated for all descendant HWSs, one finally obtains an irreducible representation.

# 4.3. Consequences of the global conformal group

The descendant state which is also a HWS is often called a singular vector of the representation.

We can define an inner product on the Verma module indirectly via the adjoint operation:

$$
\left(W ^ {i} _ {n}\right) ^ {\dagger} = W ^ {i} _ {- n} \quad <   \phi | \phi > = 1. \tag {4.2.1}
$$

For two states with weights  $w^i, \tilde{w}^i$ , eq. (4.2.1) implies that a nonzero inner product can only occur if  $w^i = \tilde{w}^i$ . In particular, two descendant states of different level have zero inner product. Singular vectors have inner product with all other states.

Singular vectors can be constructed using screening operators  $S$ . These are characterised by the fact that the singular part of the OPE of a generator  $W^i$  with  $S$  can be written as a total derivative:

$$
W ^ {i} (z) S (w) = \sum_ {n} \frac {[ W ^ {i} S ] _ {n} (w)}{(z - w) ^ {n}} = \frac {d}{d w} [ \mathrm {s o m e t h i n g} ] + O (z - w) ^ {0}. \qquad \qquad (4. 2. 2)
$$

Using eq. (2.3.16), this is seen to be equivalent to  $[SW^i]_1 = 0$ . For any screening current  $S$ , we can define an "intertwiner"  $Q_S$  whose action on an operator  $X$  is given by:

$$
Q _ {S} X (w) = \oint_ {\mathcal {C} _ {w}} \frac {d z}{2 \pi i} S (z) X (w) = [ S X ] _ {1} (w). \tag {4.2.3}
$$

Due to eq. (4.2.2),  $Q_{S}$  commutes with the modes  $W^{i}_{m}$ . This means that  $Q_{S}$  on a HWS  $|\phi\rangle$  is either zero, or another HWS. In general,  $Q_{S}|\phi\rangle$  will be a descendant of another HWS. Hence, the intertwiners can be used to construct singular vectors.

The matrix  $S$  of the inproducts of all descendants at level  $N$  is called the Sapovalov form. It depends only on the weights  $w^{i}$  of the HWS and the central charge of the  $\mathcal{W}$ -algebra. The determinant of  $S$  is the Kac-determinant [126]. Its zeroes are related to the singular vectors in the representation.

Because highest weight descendants are null operators, they generate partial differential equations on the correlation functions [13]. Completely degenerate representations have as much independent null vectors as possible [31]. It can then be proven [13] that the HWSs which give rise to a singular vector form a closed (possibly non-meromorphic) OPA. Minimal models are then defined as those cases where this OPA is finitely generated.

# 4.3 Consequences of the global conformal group

In this section, we will study the OPEs of the quasiprimaries in a  $\mathcal{W}$ -algebra. The global conformal group, generated by the modes  $L_{-1}, L_0, L_{+1}$  of the energy-momentum tensor  $T$ , puts strong restrictions on these OPEs.

Recall the definition of a quasiprimary operator  $\Phi_i$  of dimension  $h_i$  (def. 2.3.6):

$$
\begin{array}{l} \mathrm {L} _ {1} \Phi_ {i} = 0 \\ \mathrm {L} _ {0} \Phi_ {i} = h _ {i} \Phi_ {i} \tag {4.3.1} \\ \mathrm {L} _ {- 1} \Phi_ {i} = \partial \Phi_ {i}, \\ \end{array}
$$

where the modes are defined, eq. (2.4.1), as:

$$
L _ {n} \Phi \equiv [ T \Phi ] _ {n + 2}. \tag {4.3.2}
$$

From the definitions in section 4.2, we see that a quasiprimary operator furnishes a highest weight representation for the global conformal algebra.

We will use the following assumption:

Assumption 4.3.1 All elements of the OPA are linear combinations of quasiprimaries  $\Phi_i$  and their global conformal descendants, namely their derivatives.

This is a natural - and commonly used - assumption, but see intermezzo 4.3.1.

# Intermezzo 4.3.1

Consider a free field theory with background charge, see subsection 2.6.1. Using the definition of the energy-momentum tensor, eq. (2.6.17),  $T = 1 / 2[\partial X\partial X]_0 - q\partial^2 X$ , we see that  $\partial X$  is not quasiprimary when  $q$  is not zero:

$$
T \partial X = \ll 2 q | \partial X | \partial^ {2} X \gg .
$$

$X$  is not even a scaling operator:

$$
T X = \ll q | \partial X \gg .
$$

Hence, including  $\partial X$  as a generator of the OPA, would make the assumption 4.3.1 invalid.

Consider now a vertex operator  $V_{a}$ , which is primary with conformal dimension  $h_a = a(a + 2q) / 2$  with respect to  $T$ . Using eq. (2.6.14), we find for the OPE of  $V_{a}$  with  $V_{-a}$ :

$$
V _ {a} (z) V _ {- a} (w) = (z - w) ^ {- a ^ {2}} (1 + (z - w) a \partial X (w) + \dots)
$$

This OPE obeys all associativity conditions. In particular, this means that it satisfies Jacobi identities with  $T$ , given in eq. (4.3.6), for all  $q$ . Clearly, it provides an example of an OPE of two primary operators which cannot be written in a basis of quasiprimaries (if  $q \neq 0$ ), even not in a basis of highest weight operators.

Of course, by taking a different  $T$  as the Virasoro operator, the results of this chapter can be applied.

# 4.3.1 Consequences for OPEs

We will now investigate how this assumption restricts the form of the OPEs. Consider the OPE of two quasiprimary operators. The assumption 4.3.1 implies that it can be written as:

$$
\Phi_ {i} (z) \Phi_ {j} (w) = \sum_ {k} \sum_ {p \geq 0} a _ {i j} ^ {k} (p) \partial^ {p} \Phi_ {k} (w) (z - w) ^ {p - h _ {i j k}}, \tag {4.3.3}
$$

where we introduced the notation:

$$
h _ {i j k} = h _ {i} + h _ {j} - h _ {k}. \tag {4.3.4}
$$

The coefficients  $a_{ij}^{k}(p)$  are (partially) determined by the Jacobi identities, eq. (2.3.21), in the OPA. We now set out to find these coefficients. The method we use consists in acting with  $L_{1}$  on eq. (4.3.3). For  $L_{1}$  acting on the  $lhs$ , we use the Jacobi identities of  $\{T, \Phi_i, \Phi_j\}$ , while for the  $rhs$  of eq. (4.3.3) we use the identities of  $\{T, T, \Phi_k\}$ .

We write the  $lhs$  of eq. (4.3.3) with the notation introduced in eq. (2.3.3):

$$
\Phi_ {i} (z) \Phi_ {j} (w) = \sum_ {n} \left[ \Phi_ {i} \Phi_ {j} \right] _ {n} (z - w) ^ {- n}. \tag {4.3.5}
$$

From the Jacobi identities (2.3.21), or alternatively, using the commutation rules of the modes of  $\Phi_i$  with  $L_{m}$  given in eq. (2.4.7)<sup>3</sup>, gives:

$$
L _ {1} ^ {m} \left[ \Phi_ {i} \Phi_ {j} \right] _ {n} = \left(2 h _ {i} - m - n\right) _ {m} \left[ \Phi_ {i} \Phi_ {j} \right] _ {m + n}, \tag {4.3.6}
$$

where the Pochhammer symbol  $(a)_n$  is defined in eq. (2.A.1).

For the action of  $L_{1}^{m}$  on the rhs of eq. (4.3.3), we use the following formula:

$$
L _ {1} ^ {m} L _ {- 1} ^ {p} \Phi_ {k} = (p - m + 1) _ {m} (2 h _ {k} + p - m) _ {m} L _ {- 1} ^ {p - m} \Phi_ {k}, \tag {4.3.7}
$$

which can be derived using (4.A.1). Note that we dropped the term  $L_{1}\Phi_{k}$  as it is zero.

Comparing the powers of  $(z - w)^{-n}$ , we find:

$$
(2 h _ {i} - m - n) _ {m} \sum_ {k} a _ {i j} ^ {k} (h _ {i j k} - m - n) \partial^ {h _ {i j k} - m - n} \Phi_ {k} =
$$

$$
\sum_ {k} a _ {i j} ^ {k} \left(h _ {i j k} - n\right) \left(h _ {i j k} - n - m + 1\right) _ {m}
$$

$$
\left(h _ {i} + h _ {j} + h _ {k} - n - m\right) _ {m} \partial^ {h _ {i j k} - m - n} \Phi_ {k}. \tag {4.3.8}
$$

We now temporarily assume that the  $\{\Phi_i\}$  form an independent set, i.e. no linear combination of the  $\Phi_{i}$  and their derivatives can be made which is a null operator. This means that the coefficients of all  $\partial^p\Phi_k$  in the  $lhs$  and  $rhs$  have to be equal. Choosing  $m = h_{ijk} - n$ , gives:

$$
a _ {i j} ^ {k} (m) m! \left(2 h _ {k}\right) _ {m} = \left(h _ {i} - h _ {j} + h _ {k}\right) _ {m} a _ {i j} ^ {k} (0). \tag {4.3.9}
$$

Introducing the notation:

$$
a _ {i j} ^ {k} (0) = \mathcal {C} _ {i j} ^ {k} \quad a _ {i j} ^ {k} (n) = \mathcal {C} _ {i j} ^ {k} \alpha \left(h _ {i}, h _ {j}, h _ {k}, n\right), \tag {4.3.10}
$$

we conclude that if all  $h_k > 0$ :

$$
\Phi_ {i} (z) \Phi_ {j} (w) = \sum_ {k} \sum_ {n \geq 0} \mathcal {C} _ {i j} ^ {k} \alpha \left(h _ {i}, h _ {j}, h _ {k}, n\right) \partial^ {n} \Phi_ {k} (w) (z - w) ^ {n - h _ {i j k}}, \tag {4.3.11}
$$

where:

$$
\alpha \left(h _ {i}, h _ {j}, h _ {k}, n\right) = \frac {\left(h _ {i} - h _ {j} + h _ {k}\right) _ {n}}{n ! \left(2 h _ {k}\right) _ {n}}. \tag {4.3.12}
$$

Eqs. (4.3.11,4.3.12) were derived in [13]. They enable us to reconstruct the complete OPE when the structure constants  $\mathcal{C}_{ij}^{k}$  are given, i.e. the coefficients of the quasiprimaries at every pole. It is easily verified that the OPEs of the  $\mathcal{W}_3$ -algebra have the structure given in eqs. (4.3.11,4.3.12).

We now treat some special cases where eqs. (4.3.11, 4.3.12) are not valid.

When  $h_k$  is half-integer and not strictly positive, eq. (4.3.9) indicates that  $a_{ij}^k (1 - 2h_k)$  is not determined by the global conformal transformations. This is because in this case  $\partial^{1 - 2h_k}\Phi_k$  is a quasiprimary (with dimension  $1 - h_k$ ), as is easily checked using eq. (4.3.7). This situation will be mirrored in the next section when considering the full conformal group. It occurs because the highest weight representation generated by  $\Phi_k$  for  $h_k \leq 0$  is reducible.

When null operators occur in the OPA, extra free coefficients in eq. (4.3.8) can appear. However, we can continue to use eqs. (4.3.11,4.3.12) as indeed coefficients if null operators are arbitrary.

# 4.3.2 Finding the quasiprimaries in an OPE

For practical applications, the quasiprimaries at each pole have to be identified when the complete OPE is given. This is the reverse problem of the previous subsection.

We write the quasiprimary at the  $m$ -th order as  $QP^m(\Phi_i, \Phi_j)$ . It can be determined in a recursive way by using the results of the previous subsection. Indeed, the highest pole of the OPE is by assumption 4.3.1 quasiprimary. We can subtract the complete "tower" of derivatives of this quasiprimary from the OPE to end up with

# 4.3. Consequences of the global conformal group

a new Laurent series where the most singular term is again quasiprimary. We find:

$$
\begin{array}{l} Q P ^ {m} \left(\Phi_ {i}, \Phi_ {j}\right) = [ A B ] _ {m} \\ - \sum_ {n \geq 1} \alpha \left(h _ {i}, h _ {j}, h _ {i} + h _ {j} - m - n, n\right) \partial^ {n} \left(Q P ^ {m + n} \left(\Phi_ {i}, \Phi_ {j}\right)\right). \tag {4.3.13} \\ \end{array}
$$

In particular, for  $m = 0$  this formula provides a definition of a quasiprimary normal ordered product of two quasiprimaries, for example:

$$
\mathcal {N O} \left(T \Phi_ {i}\right) = \left[ T \Phi_ {i} \right] _ {0} - \frac {3}{2 \left(2 h _ {i} + 1\right)} \partial^ {2} \Phi_ {i}, \tag {4.3.14}
$$

of which eq. (4.1.2) is a special case. More general, we see that  $QP^{m}(\Phi_{i},\Phi_{j})$  for  $m\leq 0$  are composite quasiprimaries. The formula (4.3.13) appears also in [26].

The recursive definition eq. (4.3.13) of the operators  $QP^{m}$  is a very inefficient, and complicated, way of computing the quasiprimaries. To simplify this definition, we note that it can be rewritten as:

$$
Q P ^ {m} \left(\Phi_ {i}, \Phi_ {j}\right) = \sum_ {n \geq 0} a _ {n} ^ {m} \left(h _ {i}, h _ {j}\right) \partial^ {n} \left[ \Phi_ {i} \Phi_ {j} \right] _ {n + m}, \tag {4.3.15}
$$

where the coefficients  $a_{n}^{m}(h_{i},h_{j})$  are determined by eq. (4.3.13). We find:

$$
\sum_ {k = 0} ^ {n} a _ {n - k} ^ {m + k} \left(h _ {i}, h _ {j}\right) \alpha \left(h _ {i}, h _ {j}, h _ {i} + h _ {j} - m - k, k\right) = 0. \tag {4.3.16}
$$

We will provide a closed form for the solution of eq. (4.3.16) by requiring that eq. (4.3.15) defines a quasiprimary operator:

$$
L _ {1} Q P ^ {m} \left(\Phi_ {i}, \Phi_ {j}\right) = 0. \tag {4.3.17}
$$

To do this, we need the action of  $L_{1}$  on all the terms in eq. (4.3.15):

$$
L _ {1} \left(\partial^ {n} \left[ \Phi_ {i} \Phi_ {j} \right] _ {p}\right) = (2 h _ {i} - p - 1) \partial^ {n} \left[ \Phi_ {i} \Phi_ {j} \right] _ {p + 1} +
$$

$$
n \left(2 h _ {i} + 2 h _ {j} - 2 p + n - 1\right) \partial^ {n - 1} \left[ \Phi_ {i} \Phi_ {j} \right] _ {p}, \tag {4.3.18}
$$

which can be derived by commuting  $L_{1}$  to the right. The result contains a term similar to eq. (4.3.7), and a term coming from the action of the  $L_{1}$  on  $[\Phi_i\Phi_j]_p$ , see eq. (4.3.6). We find for eq. (4.3.17):

$$
\begin{array}{l} \sum_ {n \geq 1} \partial^ {n - 1} \left[ \Phi_ {i} \Phi_ {j} \right] _ {m + n} \left(a _ {n - 1} ^ {m} \left(h _ {i}, h _ {j}\right) (2 h _ {i} - n - m) + \right. \\ \left. a _ {n} ^ {m} \left(h _ {i}, h _ {j}\right) n \left(2 h _ {i} + 2 h _ {j} - 2 m - n - 1\right)\right) = 0. \tag {4.3.19} \\ \end{array}
$$

In the general case, this amounts to a recursive relation<sup>4</sup>:

$$
a _ {n} ^ {m} \left(h _ {i}, h _ {j}\right) = - \frac {2 h _ {i} - n - m}{n \left(2 h _ {i} + 2 h _ {j} - 2 m - n - 1\right)} a _ {n - 1} ^ {m} \left(h _ {i}, h _ {j}\right). \tag {4.3.20}
$$

Choosing  $a_0^m$  equal to 1, gives:

$$
a _ {n} ^ {m} \left(h _ {i}, h _ {j}\right) = (- 1) ^ {n} \frac {\left(2 h _ {i} - n - m\right) _ {n}}{n ! \left(2 h _ {i} + 2 h _ {j} - 2 m - n - 1\right) _ {n}}. \tag {4.3.21}
$$

# Intermezzo 4.3.2

We can directly check that the coefficients given in eq. (4.3.21) indeed satisfy eq. (4.3.16). This involves the summation of terms which are a product of factorials. This summation, and most of the other sums that will be used in this chapter, can be done using the Algebra 'SymbolicSum' package of Mathematica, which provides an implementation of the Gosper algorithm. However, the current version of this package does not handle the Pochhammer function. This can be remedied by converting the Pochhammer functions to a quotient of  $\Gamma$ -functions, eq. (2.A.1). The results returned by SymbolicSum in the case at hand contain a cosecans. This can again be converted to  $\Gamma$ -functions using the identity:

$$
\sin (\pi x) = \frac {\pi}{\Gamma (x) \Gamma (1 - x)}.
$$

The result then contains terms like  $\Gamma(-n)$  which are infinite for positive  $n$ . However, as all sums in this chapter are finite, these infinities necessarily disappear against terms like  $\Gamma(2 - n)^{-1}$ . It is possible to write a set of Mathematica rules which checks this cancellation automatically by converting quotients of  $\Gamma$ -functions back to Pochhammer symbols. All sums in this chapter are computed in this way.

In the case at hand, eq. (4.3.16), the final result contains a factor:

$$
- 2 n \sqrt {\pi} \Gamma (2 n) + 4 ^ {n} \Gamma (\frac {1}{2} + n) \Gamma (1 + n),
$$

which is zero (except at the singularities), proving that eq. (4.3.21) provides the solution of eq. (4.3.16).

As the operator  $QP^{m}(\Phi_{j},\Phi_{i})$  extracts the quasiprimary at the  $m$ -th order pole, we expect that reversing  $i$  and  $j$  does not give a new operator. Indeed, an explicit calculation gives<sup>5</sup>:

$$
Q P ^ {m} \left(\Phi_ {j}, \Phi_ {i}\right) = (- 1) ^ {i j + m} Q P ^ {m} \left(\Phi_ {i}, \Phi_ {j}\right), \tag {4.3.22}
$$

where  $ij$  in the phase factor gives a sign depending on the parity of the quasiprimaries (-1 if both are fermionic, 1 otherwise). This equation can be proven using eq. (2.3.16) and the identity:

$$
\sum_ {m = 0} ^ {n} \frac {(- 1) ^ {m}}{(n - m) !} a _ {m} ^ {q} \left(h _ {j}, h _ {i}\right) = a _ {n} ^ {q} \left(h _ {i}, h _ {j}\right). \tag {4.3.23}
$$

To conclude this subsection, let us discuss a special case where eq. (4.3.21) is not valid, namely when  $n \geq n_c$ , where we define:

$$
n _ {c} = 2 \left(h _ {i} + h _ {j} - m\right) - 1. \tag {4.3.24}
$$

We observe that the coefficient  $a_{n_c}^m (h_i, h_j)$  is not determined by eq. (4.3.19). On the other hand, eq. (4.3.15) shows that for such  $n \geq n_c$ , the  $a_n^m$  are coefficients of the  $n$ -th derivative of a operator with dimension  $h \leq 0$ . As discussed in the previous subsection, for any quasiprimary of dimension  $h_k \leq 0$ ,  $\partial^{1 - 2h_k}\Phi_k$  is also quasiprimary. This is clearly the origin of the freedom in  $a_{n_c}^m (h_i, h_j)$ .

Let us take as an example  $n = n_c$ , and assume that there are no poles of order higher than  $m + n_c$ . In this case  $[\Phi_i\Phi_j]_{m + n_c}$  is a quasiprimary of negative or zero conformal dimension. In this case,  $QP^m (\Phi_i,\Phi_j)$  is indeed quasiprimary for arbitrary  $a_{n_c}^m (h_i,h_j)$ .

# 4.4 Consequences of the full conformal group

In the previous section, the global conformal transformations and their implications on OPEs were studied. Here, we extend the analysis to include all modes  $L_{n}$  of the energy-momentum tensor  $T$ , acting on primary operators. A primary operator  $\Psi_{i}$  of dimension  $h_i$  satisfies (def. 2.3.6):

$$
L _ {n} \Psi_ {i} \quad = 0 \quad n > 0
$$

$$
L _ {0} \Psi_ {i} = h _ {i} \tag {4.4.1}
$$

$$
{L _ {- 1} \Psi_ {i}} {= \partial \Psi_ {i}.}
$$

Similarly to the previous section, we use the following assumption.

Assumption 4.4.1 All elements of the OPA are linear combinations of the primaries  $\Psi_{i}$  and their Virasoro descendants.

Note that we included the unit operator  $\mathbb{1}$  in the list of primaries which generate the OPA. The energy-momentum tensor  $T$  is then a Virasoro descendant of  $\mathbb{1}$ :

$$
L _ {- 2} \mathbb {1} = T. \tag {4.4.2}
$$

# 4.4.1 Restrictions of conformal covariance on the OPEs

As in the previous section, assumption 4.4.1 implies that the OPE of two primary operators can be written as:

$$
\Psi_ {i} (z) \Psi_ {j} (w) = \sum_ {k} \sum_ {\{p \}} C _ {i j} ^ {k} \beta \left(h _ {i}, h _ {j}, h _ {k}, \{p \}\right) L _ {- \{p \}} \Psi_ {k} (w) (z - w) ^ {P - h _ {i j k}}, \tag {4.4.3}
$$

where the sum is over all ordered sequences  $\{p\}$  and  $P$  is the level of  $\{p\}$  (see section 4.2 for notations). To rewrite eq. (4.4.3) in terms of composites with  $T$ , we use eq. (2.3.14):

$$
L _ {- p} \Psi = [ T \Psi ] _ {2 - p} = \frac {1}{(p - 2) !} [ \partial^ {p - 2} T \Psi ] _ {0} \quad p \geq 2. \tag {4.4.4}
$$

The  $\beta$ -coefficients defined in eq. (4.4.3) depend only on the dimensions of the primaries. They were introduced in [13], where the first two levels where computed explicitly. The equations which determine these coefficients where explicitly written down in [74] using the inner product given in section 4.2. We now show how these equations can be derived independently of the inproduct, in a way entirely analogous to the previous section.

We will act with positive modes  $L_{n}$  on eq. (4.4.3). For the lhs, the Jacobi identities (2.3.21) with  $T$  imply:

$$
L _ {p} \left[ \Psi_ {i} \Psi_ {j} \right] _ {n} = h _ {j} (p - 1) \left[ \Psi_ {i} \Psi_ {j} \right] _ {n} + \left[ \left(L _ {- 1} \Psi_ {i}\right) \Psi_ {j} \right] _ {n + 1}. \tag {4.4.5}
$$

Because  $\Psi_{i}$  is a primary operator, we have that  $L_{-1}\Psi_{i} = \partial \Psi_{i}$ , eq. (4.4.1), which gives together with eq. (2.3.13):

$$
L _ {p} \left[ \Psi_ {i} \Psi_ {j} \right] _ {h _ {i j k} - P} = \left(h _ {i} p - h _ {j} + h _ {k} + P - p\right) \left[ \Psi_ {i} \Psi_ {j} \right] _ {h _ {i j k} - P + p}. \tag {4.4.6}
$$

Extending this result to a partition  $\{p\}$  of  $P$ , we find:

$$
L _ {\{p \}} \left(\left[ \Psi_ {i} \Psi_ {j} \right] _ {h _ {i j k} - P}\right) = f \left(h _ {i}, h _ {j}, h _ {k}, \{p \}\right) \left[ \Psi_ {i} \Psi_ {j} \right] _ {h _ {i j k}}, \tag {4.4.7}
$$

where

$$
f \left(h _ {i}, h _ {j}, h _ {k}, \{p \}\right) \equiv \prod_ {l} \left(h _ {i} p _ {l} - h _ {j} + h _ {k} + \sum_ {l ^ {\prime} > l} p _ {l ^ {\prime}}\right). \tag {4.4.8}
$$

On the other hand, the action of  $L_{\{p\}}$  on the rhs of eq. (4.4.3) can be computed by commuting the positive modes to the right, where they annihilate  $\Psi_k$ . The result depends on the dimensions of the fields, but also on the central charge  $c$  appearing in the commutators of the Virasoro algebra, eq. (2.4.5). No easy formula can be given for the resulting expression, we write:

$$
L _ {\{m \}} L _ {- \{n \}} \Psi_ {k} = S ^ {\{m \}, \{n \}} \left(h _ {k}, c\right) \Psi_ {k}, \tag {4.4.9}
$$

# 4.4. Consequences of the full conformal group

for partitions of the same level  $N$ . The matrix  $S$  is the Šapovalov form defined in section 4.2. The dimension of  $S$  is given by the number of (ordered) partitions of  $N$ , which we denote by  $p(N)$ .

# Intermezzo 4.4.1

As an example we give the matrix  $S$  at level 2:

$$
\begin{array}{r c l} L _ {2} L _ {- 2} \Psi & = & (4 h + c / 2) \Psi \\ L _ {1} ^ {2} L _ {- 2} \Psi & = & 6 h \Psi \end{array} \left| \begin{array}{c c c} L _ {2} L _ {- 1} ^ {2} \Psi & = & 6 h \Psi \\ L _ {1} ^ {2} L _ {- 1} ^ {2} \Psi & = & 4 h (2 h + 1) \Psi , \end{array} \right.
$$

where  $\Psi$  is a primary with dimension  $h$ .  $S$  is a symmetric matrix. This follows from the fact that the Virasoro commutators (2.4.5) are invariant under the substitution  $L_{n} \to L_{-n}$ .

Assuming independence of the Virasoro descendants, we get:

$$
\sum_ {\{m \}} S ^ {\{n \}, \{m \}} \beta \left(h _ {i}, h _ {j}, h _ {k}, \{n \}\right) = f \left(h _ {i}, h _ {j}, h _ {k}, \{n \}\right), \tag {4.4.10}
$$

where the sum is over all partitions of  $N$ . This equation determines the  $\beta$ -coefficients at level  $N$  completely, unless the matrix  $S$  is singular. The singular vectors of  $S$  correspond to primary Virasoro descendants. Clearly, their coefficients cannot be determined by the Jacobi identities with  $T$ . This corresponds to the poles in the  $\beta$ -coefficients.

# Intermezzo 4.4.2

From the previous intermezzo, we can give the  $\beta$ -coefficients at level two:

$$
\beta \left(h _ {i}, h _ {j}, h _ {k}, \{1, 1 \}\right) = \left(c \left(h _ {i} - h _ {j} + h _ {k}\right) \left(1 + h _ {i} - h _ {j} + h _ {k}\right) + \right.
$$

$$
4 h _ {k} \left(- 4 h _ {i} + 2 h _ {i} ^ {2} + h _ {j} - 4 h _ {i} h _ {j} + 2 h _ {j} ^ {2} - \right.
$$

$$
\left. \left. h _ {k} + 4 h _ {i} h _ {k} - 4 h _ {j} h _ {k} + 2 h _ {k} ^ {2}\right)\right) /
$$

$$
\left(4 h _ {k} \left(c \left(1 + 2 h _ {k}\right) + 2 h _ {k} (- 5 + 8 h _ {k})\right)\right)
$$

$$
\beta \left(h _ {i}, h _ {j}, h _ {k}, \{2 \}\right) = \left(h _ {i} - 3 h _ {i} ^ {2} + h _ {j} + 6 h _ {i} h _ {j} - 3 h _ {j} ^ {2} - h _ {k} \right.
$$

$$
\left. \left. + 2 h _ {i} h _ {k} + 2 h _ {j} h _ {k} + h _ {k} ^ {2}\right) / \right.
$$

$$
\left(c \left(1 + 2 h _ {k}\right) + 2 h _ {k} \left(- 5 + 8 h _ {k}\right)\right).
$$

These expressions are already complicated. The coefficient of  $L_{-1}^2\Psi_k$  is undetermined when  $h_k = 0$  because  $\partial \Psi_k$  is a primary of dimension 1 in this case. Also, when

$$
c = \frac {2 (5 - 8 h _ {k}) h _ {k}}{1 + 2 h _ {k}}
$$

the quasiprimary combination:

$$
\left[ T \Psi_ {k} \right] _ {0} - \frac {3 \partial^ {2} \Psi_ {k}}{2 \left(1 + 2 h _ {k}\right)}
$$

turns out to be primary, which explains the second pole in the  $\beta$ -coefficients.

When trying to work out the  $\beta$ -coefficients at higher level, the main problem lies in solving eq. (4.4.10). Standard numerical algorithms cannot be used, as the matrix  $S$  contains the constant  $c$ , which we want to keep as a free parameter. Also, the dimension  $p(N)$  of the matrix increases rapidly with the required level, which makes solving eq. (4.4.10) for high levels unfeasible. In the next subsection, we will work in a basis of quasiprimaries. This will increase the complexity of computing  $S$ , but brings it to block-diagonal form, making the solution of equation eq. (4.4.10) much easier. Also, the primary descendants are of course necessarily quasiprimaries, so they are much easier to find in this basis.

We wish to stress that this analysis for Virasoro-primary operators cannot be extended to highest weight operators. Indeed, to determine the coefficients of the Virasoro descendants in an OPE, we explicitly used extra information about  $L_{-1}$  on a primary operator, namely that  $L_{-1}\Psi = \partial \Psi$  (see eq. (4.4.5)). As we do not have such information at present for nonlinear algebras, it is not yet possible to write down  $\mathcal{W}$ -covariant OPEs. See also section 4.1.1 in [177].

# 4.4.2 Virasoro descendants in a quasiprimary basis

Our aim in this subsection is to use the information of section 4.3 to bring the system of equations in eq. (4.4.10) in a block-diagonal form. This will be done by changing to a basis of quasiprimary Virasoro descendants and their derivatives. This basis is defined in terms of operators  $\tilde{L}_m$  and  $L_{\pm 1}$ , where:

$$
\tilde {L} _ {m} \Psi \equiv Q P ^ {m + 2} (T, \Psi), \tag {4.4.11}
$$

for  $\Psi$  quasiprimary (see eqs. (4.3.15) and (4.3.21))<sup>6</sup> The definition of  $\tilde{L}_m$  depends on the dimension of the operator on which it is acting. This can be avoided by using  $L_0$ . However, in the normalisation we chose, one needs to introduce factors  $(2L_0 + i)^{-1}$ . To avoid this, we will write  $\tilde{L}_m(h)$  when confusion can arise.

The action of  $\tilde{L}_m$  is well-defined for any quasiprimary  $\Phi$ . When  $m < -1$ , the sum in eq. (4.3.15) contains only a finite number of terms because  $a_n^{m+2}(2, h_j) = 0$  for  $n \geq 2 - m$ . When  $m > 1$ ,  $L_{n+m}\Phi = [T\Phi]_{n+m+2}$  is zero for  $n$  large enough. In particular, when  $\Phi$  is a level  $N$  descendant of a primary operator,  $L_{n+m}\Phi = 0$  if  $n + m > N$ .

We will denote partitions which do not contain the number 1 as  $\{\tilde{n}\}$ , and write  $p_1(N)$  for the number of such partitions. Note that all partitions of level  $N$  which do contain a 1 can be obtained by adding a 1 to the normal partitions of level  $N - 1$ . This implies  $p_1(N) = p(N) - p(N - 1)$ .

The descendants

$$
L _ {- 1} ^ {N - N _ {k}} \tilde {L} _ {- \left\{\tilde {n} _ {k} \right\}} \Psi_ {j} \tag {4.4.12}
$$

(where the partition  $\{\tilde{n}_k\}$  has level  $N_{k}$ ) for  $N_{k} = 0,2,3,\ldots N$  span the space of all Virasoro descendants at level  $N$ . When no primary descendants exist, they are independent.

Once we know the coefficients of the quasiprimaries at all levels  $M \leq N$ , we can use the results of the previous section to find the complete OPE:

$$
\Psi_ {i} (z) \Psi_ {j} (w) = \sum_ {k} \sum_ {\{\tilde {n} \}} C _ {i j} ^ {k} \tilde {\beta} (h _ {i}, h _ {j}, h _ {k}, \{\tilde {n} \})
$$

$$
\sum_ {m \geq 0} \alpha \left(h _ {i}, h _ {j}, h _ {k} + N, m\right) \partial^ {m} \tilde {L} _ {- \{n \}} \Psi_ {k} (w) (z - w) ^ {N + m - h _ {i j k}}. \tag {4.4.13}
$$

To determine the  $\tilde{\beta}$ -coefficients, we proceed as before, after a projection on the quasiprimaries in eq. (4.4.13). We act with a sequence of  $\tilde{L}_m$  operators,  $m \geq 2$ , on eq. (4.4.13). For the rhs, we need the matrix  $\tilde{S}$  which appears in:

$$
\tilde {L} _ {\{\tilde {m} \}} \tilde {L} _ {- \{\tilde {n} \}} \Psi_ {k} = \tilde {S} ^ {\{\tilde {m} \}, \{\tilde {n} \}} \left(h _ {k}, c\right) \Psi_ {k}, \tag {4.4.14}
$$

To compute  $\tilde{S}$ , we have to know the commutation rules for the  $\tilde{L}_n$  operators. We will need only the case where we commute a positive mode through a negative mode. We find for  $m, n > 1$ :

$$
\begin{array}{l} \tilde {L} _ {m} (h + n) \tilde {L} _ {- n} (h) = \delta_ {m - n} f _ {1} (h, m) + \tilde {L} _ {m - n} f _ {2} (h, m, - n) \\ + \sum_ {p \geq 2 - \min  (n, m)} \tilde {L} _ {- n - p} (h - m - p) \tilde {L} _ {m + p} (h) f _ {3} (h, m, - n, p), \tag {4.4.15} \\ \end{array}
$$

where we take  $\tilde{L}_{\pm 1} = 0$  and  $\tilde{L}_0 = L_0$ . The coefficients  $f_{i}$  appearing in eq. (4.4.15) are quite involved and are given in appendix 4.A, see eqs. (4.A.10, 4.A.16). The sum in eq. (4.4.15) contains an infinite number of terms, but we have to keep only the terms with  $p + m \leq N$  when acting on a Virasoro descendant of level  $N$ .

Applying  $\tilde{L}_{\{\tilde{m}\}}$  to the  $lhs$  of eq. (4.4.13), after projecting on the quasiprimaries, we get the analogue of eq. (4.4.7):

$$
\tilde {L} _ {\{\tilde {n} \}} \left(Q P ^ {h _ {i j k} - \tilde {N}} \left(\Psi_ {i}, \Psi_ {j}\right)\right) = \tilde {f} \left(h _ {i}, h _ {j}, h _ {k}, \{\tilde {n} \}\right) Q P ^ {h _ {i j k}} \left(\Psi_ {i}, \Psi_ {j}\right), \tag {4.4.16}
$$

where the  $\tilde{f} (h_i,h_j,h_k,\{\tilde{n}\}$  are given in the appendix, eq. (4.A.25).

We have the final equation:

$$
\sum_ {\{\tilde {m} \}} \tilde {S} ^ {\{\tilde {n} \}, \{\tilde {m} \}} \tilde {\beta} \left(h _ {i}, h _ {j}, h _ {k}, \{\tilde {n} \}\right) = \tilde {f} \left(h _ {i}, h _ {j}, h _ {k}, \{\tilde {n} \}\right), \tag {4.4.17}
$$

where the sum is over all partitions of  $N$ , not containing the number 1.

# Intermezzo 4.4.3

As an example, at level 2 we find:

$$
\tilde {S} ^ {\{2 \}, \{2 \}} (h, c) = \frac {c - 1 0 h + 2 c h + 1 6 h ^ {2}}{2 (1 + 2 h)}. \tag {4.4.18}
$$

This gives for  $\tilde{\beta}(h_i, h_j, h_k, \{2\})$  immediately the result of intermezzo 4.4.2. In contrast to when using the  $L_{-\{n\}}$  basis, no problem of inverting  $\tilde{S}$  occurs when we consider a dimension 0 operator<sup>7</sup>. We find:

$$
\tilde {\beta} \left(h _ {i}, h _ {j}, 0, \{2 \}\right) = \frac {h _ {i} - 3 h _ {i} ^ {2} + h _ {j} + 6 h _ {i} h _ {j} - 3 h _ {j} ^ {2}}{c}, \tag {4.4.19}
$$

which reduces to  $2h_{i} / c$  for  $h_j = h_i$ .

A particular example of a dimension zero operator is the unit operator  $\mathbb{1}$ . We see that the conventional normalisation of a primary operator of dimension  $h$  with OPE  $\Psi(z)\Psi(w) = c/h(z - w)^{-2h} + \ldots$  gives as level two descendant simply  $2T$ , i.e. with a  $c$ -independent coefficient.

The complexity of the coefficients in eqs. (4.4.15) and (4.4.16) makes these formulas unsuited for pen and paper calculations, but is no problem when using a symbolic manipulation program. As discussed in intermezzo 4.4.4, it turns out that the calculation of  $S$  takes even more time than needed for  $\tilde{S}$ . Actually, we are more interested in the quasiprimary basis, so we do not need  $S$  at all. Even more important for computer implementation [194] is that we reduced the dimension of the system of equations eq. (4.4.17) from  $p(N)$  to  $p_1(N) = p(N) - p(N - 1)$ . That this is a significant simplification is easily seen on a few examples. For levels 2 through 8,  $p(N)$  is  $\{2,3,5,7,11,15,22\}$  while  $p_1(N)$  is  $\{1,1,2,2,4,4,7\}$ .

# Intermezzo 4.4.4

As an illustration of the above arguments, we will present some timings when using an implementation in Mathematica of the formulas in this and the previous section. Timings were done in Mathematica 2.2 running on a 486 (50 MHz).

To compare the computation of the  $S$  and  $\tilde{S}$  matrices, we compute the case of an OPE

where descendants of a primary occur up to level 7, i.e. we need the matrices for level 1 to 7. As the matrices are computed recursively, we use an algorithm that stores and reuses its results. For  $S$ , the total time needed is  $20s$ , while for  $\tilde{S}$  it is only  $4s$ . For level 9, the timings become  $110s$  and  $40s$  respectively. When going to even higher level, the difference in timing decreases. However, for the most complicated algebras explicitly constructed up to now (e.g.  $\mathcal{W}A_5$  in [113]) the maximum level of a descendant of a primary (not equal to the unit operator) is 8.

The advantage of using the quasiprimary basis shows up even more when solving the equations for the  $\beta, \tilde{\beta}$ -coefficients, eqs. (4.4.10) and (4.4.17). We computed the coefficients for  $h_i = 6, h_j = 7, h_k = 10$  and  $c$  arbitrary. Although the equations we need to solve are linear equations, serious problems occurred for moderately high level.

The Mathematica (version 2.2 or lower) function Solve has a serious deficiency in that it does not simplify intermediate results. This means that when arbitrary constants are present in the equations, the solution is usually a quite big expression which simplifies drastically after factorisation. Unfortunately, this factorisation can take a very long time. For even more complicated equations Solve runs into problems because the intermediate results grow too fast. K. Hornfeck and I wrote a separate package to solve linear equations using Gaussian elimination, simplifying the coefficients of the variables at each step.

We find the following CPU-times for the solution of the equations (4.4.10) and (4.4.17).

<table><tr><td>level</td><td>β with Solve</td><td>β with LinSolve</td><td>β with Solve</td><td>β with LinSolve</td></tr><tr><td>4</td><td>6.3s</td><td>4.7s</td><td>1.6s</td><td>1.9s</td></tr><tr><td>5</td><td>* 18.4s</td><td>14.5s</td><td>1.6s</td><td>1.9s</td></tr><tr><td>6</td><td>* &gt; 400s</td><td>69.6s</td><td>* 15.3s</td><td>54.6s</td></tr></table>

In this table, a star means that factorisation of the result of Solve did not succeed in a reasonable amount of time.

As before, eq. (4.4.17) only determines the  $\tilde{\beta}$ -coefficients when  $\tilde{S}$  is nonsingular. Clearly, any singular vectors of  $\tilde{S}$  correspond to primary descendants, whose coefficients remain unfixed by Jacobi identities with  $T$ . These free coefficients will appear automatically when solving the linear equations.

# 4.4.3 Virasoro descendants of the unit operator

We now consider the case of the Virasoro descendants of the unit operator, which is a primary of dimension zero. In this case, it turns out that the matrix  $\tilde{S}$  has a large amount of singular vectors, even for general  $c$ . We usually need the descendants of  $\mathbb{1}$  to a much higher level, because most  $\mathcal{W}$ -algebras studied up to now have primary generators with dimension larger than 2. This means we should reduce the number of (linear) equations in eq. (4.4.17) as much as possible when computing  $\alpha(h_i, h_j, 0, m)$ .

The singular vectors of  $\tilde{S}$  for a dimension zero field are in general primary fields. For the unit operator however, the singular vectors are exactly zero. There is thus

no use in keeping them in the calculation. This provides an additional reason for studying the singular vectors of  $\tilde{S}$  for dimension zero.

Some examples of these singular vectors are easily constructed. Obviously we have that  $L_{-1}\mathbb{1} = \partial \mathbb{1} = 0$ .  $\tilde{L}_{-n}\mathbb{1}$  is zero for  $n > 2$  as well. Indeed, the definition (4.3.2) shows that  $\tilde{L}_{-n}\mathbb{1}$  is proportional to  $\partial^{n - 2}T$ , which is clearly not quasiprimary unless  $n = 2$ , hence it is zero for  $n \neq 2$ . One can also prove that  $\tilde{L}_{-n}T$  is zero for odd  $n$ . More complicated relations at higher level exist.

Before discussing how to find a basis for the quasiprimary descendants of  $\mathbb{1}$ , let us determine the number of independent quasiprimaries at level  $N$ , which we will denote as  $p_2(N)$ . First, the number of independent descendants is simply  $p_1(N)$ . Indeed, when ordering the  $L_{-\{n\}}$  as before, the  $L_{-1}$  will act on  $\mathbb{1}$  first, giving zero. Hence, all partitions containing 1 should be dropped. For general central charge, no further relations between the remaining descendants exist. To compute the number of quasiprimary descendants of  $\mathbb{1}$  at level  $N$ , we take all descendants, and "subtract" those which are derivatives of the descendants at the previous level, i.e.  $p_2(N) = p_1(N) - p_1(N - 1)$ .

We can use this information to find other singular vectors than the ones given already. Because  $p_2(7) = 0$ , there are no descendants at level 7, which means that  $\tilde{L}_{-3}\tilde{L}_{-2}\tilde{L}_{-2}\mathbb{1}$  is zero, which can be checked explicitly. At level 9, there is only one independent quasiprimary. Eliminating the zeroes we found at previous levels, two candidates remain. We find:

$$
\tilde {L} _ {- \{3, 2, 2, 2 \}} \mathbb {1} = - \frac {8}{5} \tilde {L} _ {- \{5, 2, 2 \}} \mathbb {1}. \tag {4.4.20}
$$

Similar relations exist at higher (not necessarily odd) levels.

Remembering the rule for the number  $p_2(N)$  of independent quasiprimary descendants of  $\mathbf{1}$ , we propose<sup>8</sup>:

Conjecture 4.4.1 The partitions which give the independent descendants of the unit operator at level  $N$  can be found as follows:

- Take the partitions (not containing the number 1) at level  $N - 1$  and order them in increasing lexicographic order ( $\{32\} < \{33\}$ ).  
- "Increment" these partitions one by one, such that no partition of level  $N$  is obtained twice.  
- The Virasoro descendants of  $\mathbb{1}$  at level  $N$  correspond to partitions at level  $N$  which are not in the list constructed in the previous step.

With "incrementing" a partition  $\{n_1, n_2, \ldots, n_k\}$ , we mean adding a 1 at a position  $i$  where  $n_{i-1} \geq n_i + 1$ . As an example, incrementing  $\{3, 2, 2\}$  gives  $\{3, 3, 2\}$  or  $\{4, 2, 2\}$ .

The first nontrivial example of this procedure is at level 9. At level 8 the partitions are  $\{22223324453628\}$ . Incrementing the first three partitions gives  $\{322233354\}$ . We now have to increment the partition  $\{53\}$ . The partition  $\{54\}$  is already in our list, so we should take  $\{63\}$ . In the next steps  $\{72\}$  and  $\{9\}$  are found. The only partition of level 9 that remains is  $\{522\}$ , in agreement with eq. (4.4.20).

This conjecture was checked up to level 16.

# 4.4.4 Finding the primaries in an OPE

Entirely analogously to eq. (4.3.13), we can define:

$$
P ^ {m} (\Psi_ {i}, \Psi_ {j}) = Q P ^ {m} (A, B) -
$$

$$
\sum_ {\{\bar {n} \}} \tilde {\beta} \left(h _ {i}, h _ {j}, h _ {i} + h _ {j} - m - N, \{n \}\right) \tilde {L} _ {- \{\bar {n} \}} P ^ {m + N} \left(\Phi_ {i}, \Phi_ {j}\right). \tag {4.4.21}
$$

This can be rewritten as:

$$
P ^ {m} \left(\Psi_ {i}, \Psi_ {j}\right) = \sum_ {\{\tilde {n} \}} a _ {\{\tilde {n} \}} ^ {m} \left(h _ {i}, h _ {j}\right) \tilde {L} _ {- \{\tilde {n} \}} Q P ^ {m + N} \left(\Psi_ {i}, \Psi_ {j}\right), \tag {4.4.22}
$$

where the sum is over all partitions of level not equal to one. We have not found a closed form for this expression. Even determining the recursion relations for the coefficients  $a_{\{\tilde{n}\}}^{m}(h_{i},h_{j})$  seems unfeasible in general.

In many cases, the action of  $P^m$  will give zero, as there simply is no primary at a certain pole.

The operators  $P^m$  are particularly convenient to construct primary composites. Indeed, for  $m \geq 2$ , the operator  $P^{-m}(\Psi_i, \Psi_j)$  contains  $[\Psi_i \Psi_j]_{-m} \sim [\partial^m \Psi_i \Psi_j]_0$ . In [201], a formula is given to count the number of primaries that can be constructed.

# 4.5 An overview of  $\mathcal{W}$ -algebras

In this section, we present a overview on the  $\mathcal{W}$ -algebras which are known. For a more in-depth review, see [31], which we also follow for the nomenclature of the  $\mathcal{W}$ -algebras.

We first distinguish different classes of  $\mathcal{W}$ -algebras. "Deformable" algebras are algebras that exist for generic values of the central charge.  $\mathcal{W}$ -algebras which are only associative for specific values of  $c$  are called non-deformable. In a "freely generated"  $\mathcal{W}$ -algebra, for generic  $c$  no relations exist between the generators, i.e. no null operators appear.

Up to now, all efforts have been concentrated on algebras in the following class:

Assumption 4.5.1 The OPA is generated by a set of primary operators and the Virasoro operator.

- All generators have strictly positive conformal dimension.  
- The unit operator occurs only in OPEs between primaries of the same conformal dimension.  
- Generators which are null operators are discarded.

Not all  $\mathcal{W}$ -algebras of interest fall within this class. For instance, it was recently shown [136] that the (nonlinear)  $\mathcal{W}_3$  algebra and the Bershadsky algebra [160, 21] with bosonic generators of dimensions  $(2,3/2,3/2,2)$  can be viewed as subalgebras of linear algebras with a null operator as generator, and a non-primary dimension 1 operator. Yet even with the restrictions 4.5.1, no complete classification has been found yet, see [62].

The current status in our knowledge of  $\mathcal{W}$ -algebras can be compared to the study of Lie algebras before Cartan presented his classification. Different construction methods of  $\mathcal{W}$ -algebras are known, giving rise to series of  $\mathcal{W}$ -algebras with common features. In a sense, they are analogous to the "classical Lie algebras", which were also known via a realisation. A few isolated cases have also been constructed, although they currently seem to fit in the general pattern.

We will now comment on some of the more important construction methods and end with the (incomplete) classification proposed in [34].

# 4.5.1 Direct construction

One can start with a list of primaries of a certain dimension and attempt to construct a  $\mathcal{W}$ -algebra. One has then to determine the structure constants of the primaries. Several methods are in use to check that the OPA is associative.

A first method is based on the claim in [13] that it is sufficient to check whether all three- and four-point functions are crossing symmetric. In the perturbative conformal bootstrap method [13, 29, 76], one analyses:

$$
G _ {n m} ^ {l k} (x) \equiv \lim  _ {\substack {z \rightarrow \infty\\\varepsilon \rightarrow 0}} z ^ {2 h _ {k}} \langle \Psi_ {j} (z) \Psi_ {l} (1) \Psi_ {n} (x) \Psi_ {m} (\varepsilon) \rangle . \tag{4.5.1}
$$

The structure constants are then constrained by requiring crossing symmetry of the  $G_{nm}^{lk}(x)$ . This can be checked perturbatively around  $x = 0$ . The crossing symmetry constraints can be implemented using a group theoretical method due to Bouwknegt [29], that was generalised in [76].

The equations resulting from this method simplify drastically when considering the  $c \to \infty$  limit. This is used in [201] to check in a very efficient way if an algebra exists in this limit, which can be considered to be a necessary condition for the algebra to

exist at generic values of the central charge  $c$  (see also subsection 4.5.5). Of course, this method gives the structure constants only in the large  $c$  limit.

Another approach is to check the Jacobi identities with the help of the mode algebra. In [33], the necessary formulas are given to compute the commutator of two quasiprimaries in a basis of modes of quasiprimaries. This was used in [26, 131, 129] to explicitly construct the  $\mathcal{W}$ -algebras with a small number of primary operators of low dimensions.

Finally, one can use OPEs to construct the OPA, and use the formulas of subsection 2.3.2 (in particular eq. (2.3.21) for  $q,p > 0$ ) to check the associativity. In [113] OPEdefines was used to construct all algebras with primary generators with dimensions 3, 4, 5 and 3, 4, 5, 6. These are the most complicated  $\mathcal{W}$ -algebras explicitly constructed up to now.

It is clear that none of these methods can give a classification. Direct construction is however the only known way that can give an exhaustive list of algebras with a certain operator-content.

# 4.5.2 Subalgebras of known  $\mathcal{W}$ -algebras

Any subalgebra of an OPA trivially satisfies the associativity conditions. Two main methods exist to construct subalgebras. First, we can factor out any free fermions and bosons, as is discussed in chapter 5. This is a powerful result, as the classification can now be restricted to  $\mathcal{W}$ -algebras where such operators are not present.

A second method to construct a subalgebra is by "orbifolding". If a certain discrete symmetry exists in the  $\mathcal{W}$ -algebra, the elements which are invariant under this symmetry necessarily form a subalgebra. As an example, the  $N = 1$  superconformal algebra is invariant under a sign change of the supersymmetry generator  $G$ . The corresponding subalgebra is purely bosonic and is generated by  $T$  and two more operators of dimension 4 and 6 respectively [111, 130].

Deformable  $\mathcal{W}$ -algebras will contain null operators at a certain value of  $c$ . This gives rise to a quotient algebra. In many cases, this algebra is only associative for this value of  $c$ . It is a common belief that all non-deformable algebras can be obtained in this way. An example of this mechanism will be given in section 4.6. Recently, it has been found that in some cases, a new deformable algebra results [27]. In the quantum case, these new algebras are finitely but non-freely generated, i.e. contain null operators for all values of  $c$ . In the classical case, an infinite number of generators is needed to construct the complete OPA. A particular feature is that the  $\mathcal{W}$ -algebras based on Lie-algebras in the same Cartan series (see next subsection) seem to contain the same subalgebra of this new type (for different values of  $c$ ). The study of these "unifying"  $\mathcal{W}$ -algebras is based on a formula for the structure constants of the operators of the lowest dimensions in the  $\mathcal{W}$ -algebra [114].

# 4.5.3 Constructing a  $\mathcal{W}$ -algebra via a realisation

If an algebra is found that is realised in terms of the generators of an associative algebra, the Jacobi identities are automatically satisfied. One has to prove that the algebra closes on a finite number of operators.

The main advantage of using a representation is that the representation theory of the underlying OPA can be used to study the representations of the  $\mathcal{W}$ -algebra itself. As we will discuss below, a whole series of  $\mathcal{W}$ -algebras has been found via realisations in a Kac-Moody algebra. Using the representation theory of the Kac-Moody algebras, the minimal models for the corresponding  $\mathcal{W}$ -algebras have been completely characterised [31].

A first example is given by the coset construction [99, 6]. For a Kac-Moody algebra  $\hat{g}$  of level  $k$  which has a subalgebra  $\hat{g}'$ , one defines the coset algebra  $\mathcal{W}_c[\hat{g} / \hat{g}', k]$  as the set of elements in  $\hat{g}$  which commute with the elements of  $\hat{g}'$ . The structure of the algebra generically depends on  $k$ , due to the appearance of null operators at specific levels. For certain cases, a representation of the same deformable algebra is found for all  $k$ , e.g. for a diagonal coset  $\mathcal{W}_c[\hat{g} \oplus \hat{g} / \hat{g}_{\mathrm{diag}}, k]$ .

One can define a second type of coset algebras by considering all elements in  $\hat{g}$  which commute with the originating Lie algebra  $g$ . These algebras can be obtained in a certain limit from the diagonal coset algebras [6]. For a simply-laced Lie algebra  $\bar{g}$  and at level  $k = 1$ , the resulting  $\mathcal{W}$ -algebras are the same as the Casimir algebras defined in [5]. This is proven in [31] using character techniques. For non-simply laced algebras, the situation is more complicated. Nevertheless, one still refers to  $\mathcal{W}_c[\hat{g} / g, k]$  as the Casimir algebra of  $\hat{g}$  at level  $k$ .

As an example, the Casimir algebra for  $B_{n}$  at level 1 turns out to be the bosonic projection of an algebra generated by the Casimir operators of  $\hat{B}_n$  and a fermionic operator of dimension  $n + 1 / 2$ . For convenience, we call the underlying  $\mathcal{W}$ -algebra, the Casimir algebra of  $\hat{B}_n$ ,  $\mathcal{W}_c B_n$ . Section 4.6 contains an explicit construction of  $\mathcal{W}_c B_2$ .

A third method to construct realisations of  $\mathcal{W}$ -algebras is by using Drinfeld-Sokolov reduction. This method was developed for classical  $\mathcal{W}$ -algebras in [60, 7], quantisation is discussed in chapter 6. Given a (super) Kac-Moody algebra, one chooses an  $sl(2)$  embedding. One then puts certain constraints on the currents of the Kac-Moody algebra. These constraints are characterised by the  $sl(2)$ . As such, one finds a realisation of a finitely generated  $\mathcal{W}$ -algebra for any  $sl(2)$ -embedding of a (super) Lie algebra. The dimensions of the operators are given by  $j_{\alpha} + 1$ , where the adjoint representation of the Lie algebra contains the  $sl(2)$ -irreducible representations  $\underline{j}_{\alpha}$ . Clearly, no dimension  $1/2$  operators can be present in a  $\mathcal{W}$ -algebra arising from Drinfeld-Sokolov reduction. This is no severe restriction, as it is shown in chapter 5 that any  $\mathcal{W}$ -algebra can be written as the direct product of some free operators with a  $\mathcal{W}$ -algebra without dimension  $1/2$  operators.

All known finitely generated  $\mathcal{W}$ -algebras in the class 4.5.1 can be obtained from

Drinfeld-Sokolov reduction, either directly or by applying the methods of the previous subsection, see [71]. As an example, it is argued in [31] that, for simply laced algebras, the diagonal coset algebras  $\mathcal{W}_c$ $[\hat{g} \oplus \hat{g} / \hat{g}_{\text{diag}}, k]$  are the same as the Drinfeld-Sokolov reductions  $\mathcal{W}_{DS}[\hat{g}, k']$  if one takes the principal  $sl(2)$  embedding of  $g^9$ . For non-simply laced algebras, the relation between coset and DS-algebras is more complicated. As an example, one finds at least classically, that  $\mathcal{W}_c B_n = \mathcal{W}_{DS}[B(0,n), k]$  [121, 205], for generic  $n$ , see also section 4.6.

Before the Drinfeld-Sokolov reduction was quantised, several other attempts have been made to find a method to construct realisations of quantum  $\mathcal{W}$ -algebras. We mention only the work of Fateev and Lukyanov [63, 64, 65], which is based on a quantisation of the Miura transformation. This transformation relates different gauge choices in the constrained phase space, and gives a realisation in terms of "simpler" operators [45], in particular in terms of free operators for the simply laced algebras. In section 4.6, a free operator realisation of  $\mathcal{W}_cB_2$  is constructed which is shown to correspond to the Fateev-Lukyanov construction.

# 4.5.4 Superconformal algebras

For superconformal algebras, the extra structure given by the supersymmetry transformations puts strong restrictions on the number of primary operators in the theory. By studying the Jacobi identities for the linear superconformal algebras with nonzero central extension, we were able to classify all possible linear algebras with generators of positive dimension [109]. The only such algebras that exist have a number of supersymmetry generators  $N \leq 4$ . The  $N = 3$  and the "large"  $N = 4$  algebras appear in section 5.3. The  $N = 1,2,3$  and the "small"  $N = 4$  algebra are all subalgebras of the large  $N = 4$ . The superconformal algebras with quadratic nonlinearity were classified in [85].

# 4.5.5 Attempts towards a classification

The methods which rely on a construction of a  $\mathcal{W}$ -algebra, either directly or via a realisation, can give no clue if all  $\mathcal{W}$ -algebras are obtained. In [34], a first attempt is made to classify all classical positive-definite  $\mathcal{W}$ -algebras, which are algebras in the class 4.5.1 with the additional condition that the central extensions define a positive-definite metric. It is proven that the finite-dimensional algebra defined by the linearised commutators of the "vacuum preserving modes"  $(\Psi_i)_m$  for  $|m| < h_i$ , is the direct sum of a semisimple Lie algebra with an abelian algebra<sup>10</sup>. Moreover, this Lie algebra necessarily contains an  $sl(2)$  formed by the  $-1,0,+1$  modes of  $T$ . The vacuum preserving modes of  $\Psi_i$  form a spin  $h_i - 1$  representation of the  $sl(2)$ . This

means that by classifying all possible  $sl(2)$ -embeddings, the dimensions which can occur in a classical positive-definite  $\mathcal{W}$ -algebra are obtained. Ref. [34] then proceeds by investigating under which conditions a quantum  $\mathcal{W}$ -algebra has a positive-definite classical limit ( $c \to \infty$ ).

The algebras arising from Drinfeld-Sokolov reduction (both quantum and classical) satisfy the criteria of [34], proving that at least one algebra exists for every  $sl(2)$ -embedding. It is not proven that this algebra is unique. Moreover, for many algebras constructed by orbifolding no classical limit can be found, and hence they fall outside this attempt towards classification. For more details, see [70, 43].

# 4.6 An example:  $\mathcal{W}_cB_2$

In this section we will construct  $\mathcal{W}_cB_2$ , the Casimir algebra of  $\hat{B}_2$ . We refer to subsection 4.5.3 for the terminology. This algebra contains, besides  $T$ , an extra dimension 4 operator and a fermionic dimension  $5/2$  operator. The algebra will be written down in terms of quasiprimary families. We will construct the most general realisation of this algebra with two free bosons and one free fermion and show that it is equivalent to the realisation proposed by Fateev and Lukyanov [65]. Finally, using this realisation, we check that the screening operators are related to the long and short root of  $B_2$  leading to the degenerate representations.

# 4.6.1 The  $\mathcal{W}_cB_2$ -algebra

The Lie algebra  $B_{2}$  has two independent Casimir operators of order two and four. The energy-momentum tensor  $T$  of  $\mathcal{W}_cB_2$  corresponds to the second order Casimir. The  $\mathcal{W}$ -algebra contains also a Virasoro primary operator  $W$  of dimension 4, which corresponds to the fourth order Casimir. Although there exists a deformable OPA with only this operator content [29, 108], this is not yet the Casimir algebra of  $B_{2}$ , see subsection 4.5.3. To get  $\mathcal{W}_cB_2$ , one has to introduce a primary weight  $5/2$  operator  $Q$ . This is reminiscent of the fact that the Frenkel-Kac level 1 realisation of the non-simply laced affine Lie algebras  $\widehat{B}_n$ , requires, due to the short root, the introduction of an extra fermion [100].

The OPA of  $\mathcal{W}_cB_2$  can be written schematically as:

$$
Q \times Q \longrightarrow \frac {2 c}{5} [ \mathbb {1} ] + C _ {\frac {5}{2}, \frac {5}{2}} ^ {4} [ W ],
$$

$$
W \times W \longrightarrow \frac {c}{4} [ 0 ] + C _ {4 4} ^ {4} [ W ] + C _ {4 4} ^ {6} [ \Psi ],
$$

$$
Q \times W \longrightarrow C _ {\frac {5}{2} 4} ^ {\frac {5}{2}} [ Q ], \tag {4.6.1}
$$

with  $[\mathbb{1}]$  the conformal family of the identity. From the Jacobi identities we find the symmetry property  $C_{\frac{5}{2}\frac{5}{2}}^{\frac{4}{2}} = \frac{8}{5} C_{\frac{5}{2}4}^{\frac{5}{2}} = \frac{8}{5} C_{4\frac{5}{2}}^{\frac{5}{2}}$ . The algebra generated by  $\{T,W\}$ ,

for  $C_{44}^{6} = 0$ , is well known to be associative for all values of  $c$  [29]. Here we discuss the solution where  $\Psi$  is the dimension 6 Virasoro primary  $\Psi \propto P^{-1}(Q, Q) = \partial Q Q + \text{corrections}$ , see eq. (4.4.22). Using Jacobi identities (2.3.21) or the conformal bootstrap, we find that the algebra is associative for all values of the central charge<sup>11</sup> provided the couplings are given by:

$$
\begin{array}{l} C _ {\frac {5}{2} \frac {5}{2}} ^ {4} = \sqrt {\frac {6 (1 4 c + 1 3)}{5 c + 2 2}} \epsilon_ {1}, \\ {C _ {4 4}} ^ {4} = \frac {3 \sqrt {6} (2 c ^ {2} + 8 3 c - 4 9 0)}{\sqrt {(1 4 c + 1 3) (5 c + 2 2)} (2 c + 2 5)} \epsilon_ {1}, \\ C _ {4 4} ^ {6} = \frac {1 2 \sqrt {5 (6 c + 4 9) (4 c + 1 1 5) (c - 1)} (5 c + 2 2)}{\sqrt {(1 4 c + 1 3) (7 c + 6 8) (2 c - 1) (c + 2 4)} (2 c + 2 5)} \epsilon_ {2}, \tag {4.6.2} \\ \end{array}
$$

where  $\epsilon_{1}$  and  $\epsilon_{2}$  are arbitrary signs.

For  $c = -13/14$  the subalgebra generated by  $\{T, Q\}$  corresponds to the spin  $5/2$  algebra of Zamolodchikov [211]. This is an example of a non-deformable  $\mathcal{W}$ -algebra.

As an application of sections 4.3 and 4.4, we now present the more complicated OPEs:

$$
Q \times Q = \ll \frac {2 c}{5} | 0 | 2 T | \partial T | \frac {3}{1 0} \partial^ {2} T + \frac {2 7}{5 c + 2 2} \Lambda + C _ {\frac {5}{2}, \frac {5}{2}} ^ {4} W \gg \tag {4.6.3}
$$

with  $\Lambda$  defined in eq. (4.1.2). For the OPE  $W(z)W(w)$ , we only list the quasiprimaries appearing in its singular part. The OPE can then be reconstructed using eqs. (4.3.11, 4.3.12). Denoting  $QP^{m}(W,W)$  as  $\phi_{8 - m}$ , we have:

$$
\begin{array}{l} \phi_ {0} = \frac {c}{4} \mathbb {1}, \quad \phi_ {2} = 2 T, \\ \phi_ {4} = \frac {4 2}{5 c + 2 2} \Lambda + C _ {4 4} ^ {4} W, \\ \phi_ {6} = - \frac {9 5 c ^ {2} + 1 2 5 4 c - 1 0 9 0 4}{6 (7 c + 6 8) (5 c + 2 2) (2 c - 1)} \left(\partial T \partial T - \frac {4}{5} \partial^ {2} T T - \frac {1}{4 2} \partial^ {4} T\right) \\ + \frac {2 4 (7 2 c + 1 3)}{(7 c + 6 8) (5 c + 2 2) (2 c - 1)} \left((T (T T)) - \frac {9}{1 0} \partial^ {2} T T - \frac {1}{2 8} \partial^ {4} T\right) \\ - \frac {1 4}{9 (c + 2 4)} C _ {4 4} ^ {4} \left(\partial^ {2} W - 6 T W\right) + C _ {4 4} ^ {6} \Psi . \tag {4.6.4} \\ \end{array}
$$

Finally,  $\Psi$  is the unique dimension 6 Virasoro primary given by:

$$
\begin{array}{l} \Psi = \mathcal {N} \left(Q \partial Q + \alpha_ {1} \partial^ {2} W + \alpha_ {2} T W + \alpha_ {3} \partial T \partial T + \right. \\ \left. \alpha_ {4} T \partial^ {2} T + \alpha_ {5} \partial^ {4} T + \alpha_ {6} (T (T T))\right), \tag {4.6.5} \\ \end{array}
$$

where the normalisation factor<sup>12</sup>:

$$
\mathcal {N} = - \sqrt {\frac {5 (2 c - 1) (c + 2 4) (7 c + 6 8)}{(c - 1) (4 c + 1 1 5) (6 c + 4 9) (1 4 c + 1 3)}} \tag {4.6.6}
$$

is fixed such that  $\Psi (z)\Psi (w) = (c / 6)(z - w)^{-12} + \ldots$  and the coefficients  $\alpha_{1},\ldots ,\alpha_{6}$  are:

$$
\begin{array}{l} \alpha_ {1} = - \frac {1 3 c + 3 5 0}{3 6 (c + 2 4)} C _ {\frac {5}{2}, \frac {5}{2}} ^ {4}, \alpha_ {2} = \frac {1 9}{3 (c + 2 4)} C _ {\frac {5}{2}, \frac {5}{2}} ^ {4}, \\ \alpha_ {3} = - 6 (5 6 6 c ^ {2} + 5 2 9 5 c - 3 9 9 8) / N, \quad \alpha_ {4} = - 6 (5 3 0 c ^ {2} + 6 1 4 1 c - 1 4 8 7) / N, \\ \alpha_ {5} = - (4 6 c ^ {3} - 1 2 5 c ^ {2} - 6 2 3 5 c + 1 7 7 8) / N, \quad \alpha_ {6} = 1 2 (7 3 4 c + 4 9) / N, \\ \end{array}
$$

with  $N = 12(7c + 68)(5c + 22)(2c - 1)$

The fermionic primary  $Q$  can be viewed as some kind of "generalised supersymmetry" generator. The appearance of the dimension 4 operator is then similar to the appearance of the affine  $U(1)$  current in the  $N = 2$  super Virasoro algebra. Similarly, one can view the  $\mathcal{W}_cB_n$  algebras as some "generalised supersymmetry" algebras, the role of the higher spin operators being to ensure associativity of the OPA for generic values of the central charge.

# 4.6.2 Coulomb Gas Realisation

The next step in the analysis of the  $\mathcal{W}_cB_2$  algebra consists in constructing a Coulomb gas realisation, i.e. a realisation of the primary operators of the OPA in terms of free operators. In order to find this realisation, one needs two free bosons and one free fermion. The OPEs for the free operators are defined to be:

$$
\begin{array}{l} \partial \varphi_ {i} \times \partial \varphi_ {j} = \ll \delta_ {i j} | 0 \gg \\ \psi \times \psi = \ll 1 \gg . \tag {4.6.7} \\ \end{array}
$$

The energy-momentum tensor is simply the free energy-momentum tensor of the two free bosons, with a background charge, and the free fermion. Due to rotational invariance, one can always transform the background charge into one direction. Without limiting the generality, we take:

$$
T = \frac {1}{2} \partial \varphi_ {1} \partial \varphi_ {1} + \frac {1}{2} \partial \varphi_ {2} \partial \varphi_ {2} + \alpha_ {0} \partial^ {2} \varphi_ {1} + \frac {1}{2} \partial \psi \psi , \tag {4.6.8}
$$

which satisfies a Virasoro algebra with central charge:

$$
c = \frac {5}{2} - 1 2 \alpha_ {0} ^ {2}. \tag {4.6.9}
$$

It is a long and boring task to find the explicit form of the dimension 4 and dimension  $5/2$  primaries. Using OPErefs, one can try to construct the most general primary dimension four operator and require the  $W(z)W(w)$  OPE to be satisfied. This, however, leads to a system of quadratic equations which is very difficult to solve. A somewhat easier way is to construct the most general primary dimension  $5/2$  operator. This leads to a four parameter family of such primaries. The  $Q(z)Q(w)$  OPE gives (up to some discrete automorphisms) three possible solutions for these parameters, and hence also three candidate dimension four primary operators. Finally, matching the  $W(z)W(w)$  OPE eliminates two of the solutions, and yields a unique construction. We wish to stress that we have checked all these statements explicitly using OPErefs, including the appearance of the dimension 6 primary mentioned earlier.

Let us now present the explicit solution<sup>13</sup>. The dimension 5/2 operator is given by:

$$
\begin{array}{l} Q = \xi \left(\frac {3}{2} \partial \varphi_ {1} \partial \varphi_ {1} \psi - \frac {3}{2} \partial \varphi_ {2} \partial \varphi_ {2} \psi + 4 \partial \varphi_ {1} \partial \varphi_ {2} \psi + \alpha_ {0} \partial^ {2} \varphi_ {1} \psi + \right. \\ \left. 3 \alpha_ {0} \partial^ {2} \varphi_ {2} \psi + 4 \alpha_ {0} \partial \varphi_ {1} \partial \psi + 2 \alpha_ {0} \partial \varphi_ {2} \partial \psi + 2 \alpha_ {0} ^ {2} \partial^ {2} \psi\right), \tag {4.6.10} \\ \end{array}
$$

where  $\xi = 1 / \sqrt{5(5 - 4\alpha_0^2)}$ . The dimension 4 operator is somewhat more complicated:

$$
\begin{array}{l} W = \sigma \Big (N ^ {i j k l} \partial \varphi_ {i} \partial \varphi_ {j} \partial \varphi_ {k} \partial \varphi_ {l} + N ^ {i j k} \partial^ {2} \varphi_ {i} \partial \varphi_ {j} \partial \varphi_ {k} + N ^ {i j} \partial^ {3} \varphi_ {i} \partial \varphi_ {j} \\ + \tilde {N} ^ {i j} \partial^ {2} \varphi_ {i} \partial^ {2} \varphi_ {j} + N ^ {i} \partial^ {4} \varphi_ {i} + M ^ {i j} \partial \varphi_ {i} \partial \varphi_ {j} \partial \psi \psi \\ + M ^ {i} \partial^ {2} \varphi_ {i} \partial \psi \psi + \tilde {M} ^ {i} \partial \varphi_ {i} \partial^ {2} \psi \psi + K _ {1} \partial^ {3} \psi \psi + K _ {2} \partial^ {2} \psi \partial \psi), \tag {4.6.11} \\ \end{array}
$$

with  $N^{ijkl}$ ,  $\tilde{N}^{ij}$  and  $M^{ij}$  completely symmetric and  $N^{ijk}$  symmetric in the last two indices. The coefficients appearing in (4.6.11) are given explicitly by:

$$
\begin{array}{l} N ^ {1 1 1 1} \quad = \quad N ^ {2 2 2 2} = 8 1 / 8 0, \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \left. N ^ {1 1 1 2} \right. = - N ^ {1 2 2 2} = - \mu / 2 0, \\ N ^ {1 1 2 2} = (5 6 0 \alpha_ {0} ^ {2} - 7 9) / 7 2 0, \\ N ^ {1 1 1} \quad = \quad 8 1 \alpha_ {0} / 2 0, \\ {N ^ {1 2 2}} = {\alpha_ {0} (8 0 \alpha_ {0} ^ {2} + 1 9 7) / 6 0,} \\ N ^ {2 1 2} \quad = - \alpha_ {0} \mu / 1 0, \\ N ^ {1 1} \quad = \quad (1 2 8 \alpha_ {0} ^ {2} - 2 5) / 6 0, \\ N ^ {2 1} = - \alpha_ {0} ^ {2} \mu / 5, \\ \tilde {N} ^ {1 1} = (3 4 \alpha_ {0} ^ {2} + 2 5) / 4 0, \\ \tilde {N} ^ {2 2} \quad = \quad (8 0 \alpha_ {0} ^ {4} - 1 7 4 \alpha_ {0} ^ {2} + 2 5) / 4 0, \\ N ^ {1} \quad = \quad \alpha_ {0} (1 2 8 \alpha_ {0} ^ {2} - 2 5) / 3 6 0, \\ M ^ {1 1} \quad = \quad M ^ {2 2} = (8 2 \alpha_ {0} ^ {2} - 3 5) / 6, \\ M ^ {1} \quad = \quad \alpha_ {0} (2 \alpha_ {0} ^ {2} + 1 1) / 3, \\ \tilde {M} ^ {1} \quad = \quad - \alpha_ {0} \mu / 3, \\ K _ {1} \quad = \quad (3 6 \alpha_ {0} ^ {4} - 2 2 \alpha_ {0} ^ {2} + 5) / 1 8, \\ \end{array}
$$

where

$$
\sigma = \sqrt {\frac {3}{2 (2 - 7 \alpha_ {0} ^ {2}) (2 3 - 4 0 \alpha_ {0} ^ {2})}} \frac {1}{4 \alpha_ {0} ^ {2} - 5}, \quad \mu = 2 3 - 4 0 \alpha_ {0} ^ {2}. \qquad (4. 6. 1 2)
$$

This solution corresponds to the sign choices  $\epsilon_{1} = \epsilon_{2} = +1$  in (4.6.2).

Some remarks are in order.

The primary dimension  $5/2$  operator  $Q$  can be rewritten in a more suggestive form. Indeed, rotating the free scalars:

$$
{\bar {\varphi} _ {1}} = {\frac {1}{\sqrt {1 0}} (3 \varphi_ {1} - \varphi_ {2}),}
$$

$$
\bar {\varphi} _ {2} = \frac {1}{\sqrt {1 0}} \left(\varphi_ {1} + 3 \varphi_ {2}\right), \tag {4.6.13}
$$

one can rewrite:

$$
Q = 5 \left(\sqrt {\frac {2}{5}} \alpha_ {0} \partial + \partial \bar {\varphi} _ {1}\right) \left(\sqrt {\frac {2}{5}} \alpha_ {0} \partial + \partial \bar {\varphi} _ {2}\right) \psi , \tag {4.6.14}
$$

which is exactly the starting point of the analysis of Fateev and Lukyanov [65]. Our construction hence proves that the algebra generated by (4.6.14) is indeed finitely generated for all values of the central charge, as was conjectured in [65].

A second remark concerns the fact that, (up to some discrete automorphisms), there is (up to some discrete automorphisms) only one free field realisation of  $\mathcal{W}_c$ $B_{2}$  with two free bosons and one free fermion. In the case of  $\mathcal{W}_3$ , Fateev and Zamolodchikov found two inequivalent free field realisations with two free bosons [66]. The simplest one was related to  $su(3)$  and led Fateev and Lukyanov to generalise this to the  $\mathcal{W}A_{n}$  algebras using  $n$  free bosons [63]. The more complicated one has

# 4.7. Discussion

been [155] shown to be related to parafermions (at least for a specific value of the central charge). In fact, it was argued in [155] that such a realisation exists (for fixed  $c$ ) for all  $\mathcal{W}A_{n}$  algebras, and that in the limit  $n\to \infty$  it corresponds to the  $c = 2$  free field realisation of  $\mathcal{W}_{\infty}$  by Bakas and Kiritsis [8]. For  $\mathcal{W}_cB_2$ , there does not seem to be a similar construction<sup>14</sup>.

# 4.6.3 Highest weight representations

Examples of  $\mathcal{W}_cB_2$ -highest weight operators are easily constructed in the Coulomb gas realisation. They are given by the vertex operators defined in subsection 2.6.1:

$$
V _ {\vec {\beta}} (z) = e ^ {\vec {\beta}. \vec {\varphi}} (z), \tag {4.6.15}
$$

with  $\vec{\beta} \equiv (\beta_{1},\beta_{2})$  and  $\vec{\varphi}(z) \equiv (\varphi_{1}(z),\varphi_{2}(z))$ .  $V_{\vec{\beta}}(z)$  has Virasoro dimension:

$$
\Delta_ {\bar {\beta}} = \frac {1}{2} \beta_ {1} ^ {2} + \frac {1}{2} \beta_ {2} ^ {2} - \beta_ {1} \alpha_ {0}, \tag {4.6.16}
$$

To write down the  $W$ -weight, we introduce the following notation, connected to the root system of  $B_{2}$ .  $\vec{e}_{L} \equiv \sqrt{\frac{2}{5}} (1, -2)$  and  $\vec{e}_{S} \equiv \sqrt{\frac{1}{10}} (1, 3)$  are the positive simple roots.  $\vec{\rho} = \frac{1}{2} (3\vec{e}_{L} + 4\vec{e}_{S})$  is half the sum of the positive roots; note that  $\vec{\alpha}_0 \equiv (\alpha_0, 0)$  is parallel to  $\vec{\rho}$ . Finally,  $W$  denotes the Weyl group of  $B_{2}$ . With this notation, the  $W$ -weight of the vertex operator (4.6.15) can be written as:

$$
\begin{array}{l} w _ {\vec {\beta}} = \frac {\sigma}{4 8 0} \Big ((4 0 \alpha_ {0} ^ {2} - 2 3) \left(\prod_ {w \in W} w (\vec {\beta} - \vec {\alpha} _ {0}). \vec {\rho}\right) ^ {1 / 2} \\ \left. + 8 \left(1 2 8 \alpha_ {0} ^ {2} - 2 5\right) \Delta_ {\vec {\beta}} + 1 9 4 4 \Delta_ {\vec {\beta}} ^ {2}\right), \tag {4.6.17} \\ \end{array}
$$

where  $\sigma$  was defined in eq. (4.6.12). From this formula, it follows immediately that the weights  $\Delta_{\vec{\beta}}$  and  $w_{\vec{\beta}}$  are invariant under  $\vec{\beta} \rightarrow w(\vec{\beta} - \vec{\alpha_0}) + \vec{\alpha_0}$ , with  $w \in \mathcal{W}$  [65].

In the case at hand, one finds at the third order pole  $[W V_{\vec{\beta}}]_3$  a new Virasoro primary operator, proportional to  $(\beta_2\partial \varphi_1 + (2\alpha_0 - \beta_1)\partial \varphi_2))V_{\vec{\beta}}$ . As mentioned in section 4.2, we have no intrinsic geometric way to express this new Virasoro primary as well as the ones appearing in the first and second order poles.

Using the realisation in terms of free fields of the previous section, one can construct two different kinds of screening operators. They can be written in such a way

that their relation with the root system of  $B_{2}$  is manifest:

$$
\begin{array}{l} {V _ {\pm} ^ {L}} = {\exp (\beta_ {\pm} \vec {e} _ {L}. \vec {\varphi})} \\ V _ {\pm} ^ {S} = \psi \exp \left(\beta_ {\pm} \vec {e} _ {S}. \vec {\varphi}\right), \tag {4.6.18} \\ \end{array}
$$

with  $\beta_{+} + \beta_{-} = \sqrt{\frac{2}{5}}\alpha_{0}$  and  $\beta_{+}\beta_{-} = -1$ . They are thus seen to be equal to the screening charges presented in [65].

Given the screening operators, it is a standard construction to derive the degenerate representations, see [128, 66, 65, 203] to which we refer for details.

# 4.7 Discussion

In this chapter we analysed in detail what the consequences are of the global and local conformal group on the OPEs in a  $\mathcal{W}$ -algebras. We provided explicit formulas for working with primaries and quasiprimaries. Although some of the formulas are quite complicated, computer implementation presents no problem [194]. A further step would be to implement the formulas of working with OPEs of (quasi)primaries in a future version of  $OPE$  defs, i.e. without expanding the (quasi)primaries in order to calculate an OPE.

As an example of the power of using the techniques of primary quasiprimary operators, combined with automated OPEs, we have proven the existence, for generic  $c$ , of the Casimir algebra of  $B_{2}$  by explicitly constructing it. Using a Coulomb gas realisation in terms of two free bosons and one free fermion, we have been able to show the equivalence with the results conjectured by Fateev and Lukyanov [65].

# 4.A Appendix

This appendix collects some of the more technical details of this chapter.

# A few formulas with the modes of the energy-momentum tensor

First we give a number of identities which follow from the Virasoro algebra of the modes of the energy-momentum tensor (2.4.5).

For all  $m\in \mathbb{Z}$ $n\in \mathbb{N}$  we have:

$$
L _ {m} L _ {- 1} ^ {n} = \sum_ {k = 0} ^ {m + 1} \binom {n} {k} (m - k + 2) _ {k} L _ {- 1} ^ {n - k} L _ {m - k} \tag {4.A.1}
$$

$$
L _ {- 1} ^ {n} L _ {m} = \sum_ {k = 0} ^ {m + 1} \binom {n} {k} (- m - 1) _ {k} L _ {m - k} L _ {- 1} ^ {n - k}. \tag {4.A.2}
$$

# Derivation of eq. (4.4.15)

We repeat eq. (4.4.15) here for convenience:

$$
\tilde {L} _ {m} (h + n) \tilde {L} _ {- n} (h) = \delta_ {m - n} f _ {1} (h, m) + \tilde {L} _ {m - n} f _ {2} (h, m, - n) +
$$

$$
\sum_ {p \geq 2 - \min  (n, m)} \tilde {L} _ {- n - p} (h - m - p) \tilde {L} _ {m + p} (h) f _ {3} (h, m, - n, p), \tag {4.A.3}
$$

where  $n, m > 1$  for the remainder of this appendix.

We also introduce a shorter notation for the coefficients in eq. (4.4.11):

$$
\tilde {L} _ {m} (\Phi) = \sum_ {n \geq 0} b _ {n} ^ {m} (h) L _ {- 1} ^ {n} L _ {n + m} \Phi , \quad m \in \mathbb {Z} \backslash \{1, 0, - 1 \}, \tag {4.A.4}
$$

where we defined:

$$
b _ {n} ^ {m} (h) \equiv a _ {n} ^ {m + 2} (2, h) = (- 1) ^ {n} \frac {(2 - n - m) _ {n}}{n ! (2 h - 2 m - n - 1) _ {n}}. \tag {4.A.5}
$$

We will need the following lemma.

Lemma 4.A.1 On a quasiprimary  $\Phi$  of dimension  $h$ , the action of  $\tilde{L}_m$  can be written as:

$$
\tilde {L} _ {- m} (\Phi) = \sum_ {n \geq 0} \tilde {b} _ {n} ^ {m} (h) L _ {n - m} L _ {- 1} ^ {n} \Phi , \quad m > 1, \tag {4.A.6}
$$

where:

$$
\tilde {b} _ {n} ^ {m} (h) = (- 1) ^ {n} \binom {1 + m} {n} \frac {(2 h + n) _ {m - n - 2}}{(2 h + m + 1) _ {m - 2}}. \tag {4.A.7}
$$

The proof follows immediately from the definitions eq. (4.A.4) and eq. (4.A.5) of  $\tilde{L}$  and eq. (4.A.2).

We now set out to prove eq. (4.A.3). Our strategy consists of ordering all  $L_{k}$  modes in (4.A.3) such that higher modes are moved to the right. After reordering, we look only at terms which do not contain  $L_{-1}$ . This is sufficient as the other terms are fixed by requiring that both  $lhs$  and  $rhs$  of eq. (4.A.3) are quasiprimary.

For the terms in the  $\text{rhs}$  of eq. (4.A.3), the lemma immediately gives:

$$
\tilde {L} _ {- n} (h - m) \tilde {L} _ {m} (h) \rightarrow \frac {(2 h - 2 m) _ {n - 2}}{(2 h - 2 m + n + 1) _ {n - 2}} L _ {- n} L _ {m}, \tag {4.A.8}
$$

where the rightarrow means that we drop terms containing  $L_{-1}$ .

For the  $lhs$  of eq. (4.A.3), we find:

$$
\begin{array}{l} = \sum_ {i, j \geq 0} b _ {i} ^ {m} (h + n) \tilde {b} _ {j} ^ {n} (h) L _ {- 1} ^ {i} \left(L _ {- n + j} L _ {m + i} + \right. \\ (m + i + n - j) L _ {m + i - n + j} + \frac {c}{2} \left( \begin{array}{c} m + i \\ 3 \end{array} \right) \delta_ {m + i - n + j} \Big) L _ {- 1} ^ {j} \\ \end{array}
$$

$$
\begin{array}{l} \tilde {L} _ {m} (h + n) \tilde {L} _ {- n} (h) \\ = \sum_ {i, j \geq 0} b _ {i} ^ {m} (h + n) \tilde {b} _ {j} ^ {n} (h) L _ {- 1} ^ {i} \left(L _ {- n + j} L _ {m + i} + \right. \\ (m + i + n - j) L _ {m + i - n + j} + \frac {c}{2} \left( \begin{array}{c} m + i \\ 3 \end{array} \right) \delta_ {m + i - n + j} \Big) L _ {- 1} ^ {j} \\ \rightarrow \sum_ {i, j \geq 0} b _ {i} ^ {m} (h + n) \tilde {b} _ {j} ^ {n} (h) \\ \Big (  (n - j - 1) _ {i}   (m + i - j + 2) _ {j}     L _ {- n + j - i} L _ {m + i - j} + \frac {c}{2} \left( \begin{array}{c} m \\ 3 \end{array} \right) \delta_ {m - n} \delta_ {i} \delta_ {j} \\ + (m + i + n - j) L _ {m - n} \times \\ \left. \left(\left(m - n - 2\right) _ {j} \delta_ {i} \delta_ {m - n + j \geq 0} + \left(n - m - i - 1\right) _ {i} \delta_ {j} \delta_ {m - n + i <   0}\right)\right). \tag {4.A.9} \\ \end{array}
$$

Let us look at the term with  $L_{-n + j - i}L_{m + i - j}$ . An additional reordering of the modes is necessary when  $-n + j - i > m + i - j$ , which can only happen if  $n > m$ . This reordering will give an additional contribution to the term proportional to  $L_{m - n}$  in eq. (4.A.9).

We will now find the lower boundary for  $p$  in the sum  $\tilde{L}_{-n-p}\tilde{L}_{m+p}$  in eq. (4.A.3). We notice that the coefficient of  $L_{-n+j-i}L_{m+i-j}$  is zero unless  $m+i-j+2 > 0$  or  $j = 0$ . This means that the lowest  $p$  that occurs is certainly larger than  $-2 - m$ . Furthermore,  $p = -1 - m$  gives a  $L_{-1}$  contribution which we dropped,  $p = 1 - m$  gives  $L_1$  which is zero on a quasiprimary and so does not contribute either. Finally, the term  $p = -m$  simply gives  $L_0$ , and hence  $h$ , and should be added to the linear term in eq. (4.A.3). We can conclude that the lowest  $p$  which gives a quadratic contribution is larger than or equal to  $2 - m$ , such that the rightmost  $\tilde{L}$  mode of the quadratic term is always a positive mode. On the other hand, from the definition of  $\tilde{b}_j^n$  eq. (4.A.7), we have that  $j \leq n + 1$ . Together with the factor  $(n-j-1)_i$ , this means that for the leftmost  $\tilde{L}$  the highest possible mode has  $p = -n$ . For the same reasons as before, we can conclude that the minimal  $p$  should also be larger than or equal to  $2 - n$  (making the leftmost  $\tilde{L}$  always a negative mode). In this way, we see that eq. (4.A.3) has the correct form.

To determine the coefficients in eq. (4.A.3), we compare eq. (4.A.9) to eq. (4.A.8). Some of the sums in these coefficients can be found using Mathematica. We find for  $f_{1}$  and  $f_{2}$ :

$$
{f _ {1} (h, m)} = {\frac {c}{1 2} (m - 1) m (m + 1) \frac {(2 h) _ {m - 2}}{(2 h + m + 1) _ {m - 2}}}
$$

$$
f _ {2} (h, m, - n) = \left((- 1) ^ {m} (- 2 - 2 h + 2 h ^ {2} + m - 3 h m + 2 h ^ {2} m + m ^ {2} - h m ^ {2} - n \right.
$$

# 4.A. Appendix

$$
+ 2 h n - 2 m n + 2 h m n + n ^ {2}) (2 - m + n) _ {m} +
$$

$$
(2 - m + 2 h m - m ^ {2} - n + 2 h n + 2 m n + n ^ {2}) (- 2 + 2 h - m + n) _ {m}
$$

$$
\left. \right) / (- 2 + 2 h - m + 2 n) _ {1 + m}, \quad \text {f o r} n \geq m \tag {4.A.10}
$$

In the case  $n < m$ ,  $f_{2}$  can be determined from the following lemma.

# Lemma 4.A.2

$$
f _ {2} (h, m, - n) = f _ {2} (h - m + n, n, - m). \quad \text {f o r} n <   m \tag {4.A.11}
$$

Proof :

This relation is most easily proven by considering the inproduct (for  $n > m$ ):

$$
<   \Psi | \tilde {L} _ {n - m} (h) \tilde {L} _ {m} (h + n) \tilde {L} _ {- n} (h) \Psi > =
$$

$$
f _ {2} (h, m, - n) <   \Psi | \tilde {L} _ {n - m} (h - m + n) \tilde {L} _ {m - n} (h) \Psi >, \tag {4.A.12}
$$

where  $\Psi$  is a primary of dimension  $h$  (with  $<  \Psi |\Psi >$  nonzero) and we used eq. (4.A.3) on the two last operators of the  $lhs$ . The inproduct is defined in eq. (4.2.1). We can compute the inproduct in the  $lhs$  of eq. (4.A.12) as:

$$
<   \Psi | (\tilde {L} _ {- n} (h)) ^ {+} (\tilde {L} _ {m} (h + n)) ^ {+} (\tilde {L} _ {n - m} (h)) ^ {+} \Psi >. \tag {4.A.13}
$$

Now, substituting the definition eq. (4.4.11) for the rightmost operator, only the term  $L_{m - n}$  of the sum remains, as  $L_{1}\Psi$  is zero. Also, the first two operators acting on the left state create a quasiprimary state, which is annihilated (from the right) by  $L_{-1}$ . Hence, we can effectively replace  $\left(\tilde{L}_{n - m}(h)\right)^{+}$  by  $\tilde{L}_{m - n}(h)$  in eq. (4.A.13). The same reasoning can be followed for the other operators, but we have to shift the dimensions<sup>15</sup>. We get:

$$
<   \Psi | \tilde {L} _ {n} (h + n) \tilde {L} _ {- m} (h - m + n) \tilde {L} _ {m - n} (h) \Psi >. \tag {4.A.14}
$$

However, using eq. (4.A.3) on the two first operators, this is also equal to:

$$
f _ {2} (h - m + n, n, - m) <   \Psi | \tilde {L} _ {n - m} (h - m + n) \tilde {L} _ {m - n} (h) \Psi >, \tag {4.A.15}
$$

which proves eq. (4.A.11).

To conclude the computation of eq. (4.A.3), we give the expression for  $f_{3}$ :

$$
f _ {3} (h, n, m, p) = f _ {4} (h, n, m, p) \frac {\left(1 + 2 h - 2 m + n - p\right) _ {- 2 + n + p}}{\left(2 h - 2 m - 2 p\right) _ {- 2 + n + p}}, \tag {4.A.16}
$$

where

![](images/b571d106902afa106c9d2e3f93f6e27ae183c4fcf378f0e05ca750f5ee7a1b93.jpg)

We did not find a simple expression for the sums that are involved. We rewrote them in terms of generalised hypergeometric functions:

$$
{ } _ { p } \mathcal { F } _ { q } ( n _ { i } ; m _ { j } ; z ) = \sum _ { k \geq 0 } \frac { \prod _ { i } ( n _ { i } ) _ { k } } { k ! \prod _ { j } ( m _ { j } ) _ { k } } z ^ { k } , \tag {4.A.18}
$$

where  $i$  runs from 1 to  $p$  and  $j$  from 1 to  $q$ . In  $f_4$ , the infinite sum always reduces to a finite number of terms because one of the  $n_i$  is negative. One can now use identities for the generalised hypergeometric functions [188] to prove that:

$$
f _ {3} (h, m, - n, p) = f _ {3} (h + n - m, n, - m, p). \tag {4.A.19}
$$

Alternatively, this can be checked using an inproduct with four  $\tilde{L}$  operators.

# Determination of the coefficient in eq. (4.4.16)

We first prove:

$$
\tilde {L} _ {n} \left(h _ {j} - m\right) Q P ^ {h _ {i} + m} \left(\Psi_ {i}, \Psi_ {j}\right) = f _ {5} \left(h _ {i}, h _ {j}, m, n\right) Q P ^ {h _ {i} + m + n} \left(\Psi_ {i}, \Psi_ {j}\right), \quad n \geq 2. \tag {4.A.20}
$$

along the same lines as eq. (4.A.3), i.e. we will order the modes, and drop  $L_{-1}$  contributions. We write down the definition of the lhs using modes:

$$
\sum_ {k, l \geq 0} a _ {k} ^ {n + 2} (2, h _ {j} - m) L _ {- 1} ^ {k} L _ {n + k} a _ {l} ^ {m + h _ {i}} \left(h _ {i}, h _ {j}\right) L _ {- 1} ^ {l} \left(\widehat {\Psi_ {i}}\right) _ {m + l} \Psi_ {j}, \tag {4.A.21}
$$

where the  $a_{n}^{m}$  are given in eq. (4.3.21). We move the  $L_{-1}$  to the left using eq. (4.A.1), and drop terms containing  $L_{-1}$ :

$$
\sum_ {l \geq 0} a _ {l} ^ {m + h _ {i}} \left(h _ {i}, h _ {j}\right) (n - l + 2) _ {l} L _ {n - l} \widehat {\Psi_ {i}} _ {m + l} \Psi_ {j}. \tag {4.A.22}
$$

Now, the factor  $(n - l + 2)_l$  restricts  $l$  to be smaller than  $n + 2$ . In fact, for  $l = n + 1$  we get a  $L_{-1}$  term which we should drop. For  $l$  less than  $n$ , we can commute  $L_{n - l}$  through the  $\Psi_i$  mode. Only the commutator eq. (2.4.7) remains, as  $L_{n - l}$  annihilates the primary operator  $\Psi_j$ . We get for the  $lh_s$  of eq. (4.A.20), dropping  $L_{-1}$  terms:

$$
\begin{array}{l} \Big (\sum_ {l = 0} ^ {n - 1} a _ {l} ^ {m + h _ {i}} \left(h _ {i}, h _ {j}\right) (n - l + 2) _ {l} (n \left(h _ {i} - 1\right) - m - l) \\ \left. + a _ {n} ^ {m + h _ {i}} \left(h _ {i}, h _ {j}\right) (2) _ {n} \left(h _ {j} - m - n\right)\right) \widehat {\Psi_ {i}} _ {m + n} \Psi_ {j}. \tag {4.A.23} \\ \end{array}
$$

The coefficient of  $\widehat{\Psi}_{im + n}\Psi_{j}$  in this equation is  $f_{5}$ , as the only term in the rhs of eq. (4.A.20) without  $L_{-1}$  is simply  $\widehat{\Psi}_{im + n}\Psi_{j}$ . After summation, we find:

$$
\begin{array}{l} f _ {5} \left(h _ {i}, h _ {j}, m, n\right) = \\ \left((- 1) ^ {n} \left(- h _ {i} \left(h _ {i} - 1\right) + h _ {j} \left(h _ {j} - 1\right) + M (M - 1) + \right. \right. \\ \left. h _ {j} n (2 M + n - 1)\right) \left(h _ {i} - h _ {j} + M\right) _ {n} \\ + \left(h _ {i} \left(h _ {i} - 1\right) - h _ {j} \left(h _ {j} - 1\right) + M (M - 1) + \right. \\ h _ {i} n (2 M + n - 1)) (- h _ {i} + h _ {j} + M) _ {n} \\ \bigg) / (2 M + n - 2) _ {n + 1} \tag {4.A.24} \\ \end{array}
$$

where  $M = h_{j} - m - n$ .

Although we only looked at the term free of  $L_{-1}$ , the others have to be such that the rhs of eq. (4.A.20) is quasiprimary. Moreover, we see that the lhs of eq. (4.A.20) is of the form  $\sum x_l L_{-1}^l (\widehat{\Psi_i})_{m+n+l} \Psi_j$ . We found in subsection 4.3.2 that requiring this form to be quasiprimary fixed all  $x_j$  in terms of  $x_0$ .

It is now clear that the proportionality constant in eq. (4.4.16) is given by:

$$
\tilde {f} \left(h _ {i}, h _ {j}, h _ {k}, \{\tilde {n} \}\right) \equiv \prod_ {l} f _ {5} \left(h _ {i}, h _ {j}, h _ {j} - h _ {k} - \sum_ {k \geq l} n _ {k}, n _ {l}\right). \tag {4.A.25}
$$
