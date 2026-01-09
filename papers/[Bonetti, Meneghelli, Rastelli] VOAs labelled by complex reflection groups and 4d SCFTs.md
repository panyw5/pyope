# VOAs labelled by complex reflection groups and 4d SCFTs

# Federico Bonetti, $^{a}$  Carlo Meneghelli, $^{b}$  and Leonardo Rastellic

$^{a}$ Department of Physics and Astronomy, Johns Hopkins University, 3400 North Charles Street, Baltimore, MD 21218, USA  
$^{b}$  Mathematical Institute, University of Oxford, Andrew Wiles Building, Radcliffe Observatory Quarter, Woodstock Road, Oxford, OX2 6GG, United Kingdom  
$^{c}$ C. N. Yang Institute for Theoretical Physics, Stony Brook University, Stony Brook, NY 11794, USA

ABSTRACT: We define and study a class of  $\mathcal{N} = 2$  vertex operator algebras  $\mathcal{W}_{\mathsf{G}}$  labelled by complex reflection groups. They are extensions of the  $\mathcal{N} = 2$  super Virasoro algebra obtained by introducing additional generators, in correspondence with the invariants of the complex reflection group  $\mathsf{G}$ . If  $\mathsf{G}$  is a Coxeter group, the  $\mathcal{N} = 2$  super Virasoro algebra enhances to the (small)  $\mathcal{N} = 4$  superconformal algebra. With the exception of  $\mathsf{G} = \mathbb{Z}_2$ , which corresponds to just the  $\mathcal{N} = 4$  algebra, these are non-deformable VOAs that exist only for a specific negative value of the central charge. We describe a free-field realization of  $\mathcal{W}_{\mathsf{G}}$  in terms of rank  $(\mathsf{G})$ $\beta \gamma bc$  ghost systems, generalizing a construction of Adamovic for the  $\mathcal{N} = 4$  algebra at  $c = -9$ . If  $\mathsf{G}$  is a Weyl group,  $\mathcal{W}_{\mathsf{G}}$  is believed to coincide with the  $\mathcal{N} = 4$  VOA that arises from the four-dimensional super Yang-Mills theory whose gauge algebra has Weyl group  $\mathsf{G}$ . More generally, if  $\mathsf{G}$  is a crystallographic complex reflection group,  $\mathcal{W}_{\mathsf{G}}$  is conjecturally associated to an  $\mathcal{N} = 3$  4d superconformal field theory. The free-field realization allows to determine the elusive “ $R$ -filtration” of  $\mathcal{W}_{\mathsf{G}}$ , and thus to recover the full Macdonald index of the parent  $4d$  theory.

# Contents

# 1 Introduction and summary 1

1.1 Main results 4  
1.2 Connection with 4d physics 9  
1.3 Outlook 12

# 2 Preliminaries 13

2.1  $\mathcal{N} = 2$  and small  $\mathcal{N} = 4$  SCAs 13  
2.2 Coxeter groups, complex reflection groups, rings of invariants 15

# 3 Free-field realizations 18

3.1 Realization of the  $\mathcal{N} = 2$  SCA in terms of  $\beta \gamma bc$  systems 18  
3.2 Realization of the remaining generators in the  $\mathcal{N} = 4$  case 19  
3.3 Realization of the remaining generators in the  $\mathcal{N} = 2$  case 21  
3.4 Free-field realization and classical Poisson structure 22

# 4 Examples of  $\mathcal{N} = 4$  VOA  $\mathcal{W}_{\Gamma}$  23

4.1 Rank 1:  $\Gamma = A_{1}$  23  
4.2 Rank 2:  $\Gamma = I_2(p)$  25

4.2.1 The Coxeter groups  $I_2(p)$  and associated rings 25  
4.2.2 Free-field realization of all generators 26  
4.2.3 Null states 28  
4.2.4 Classical limit: relation between  $\beta \gamma$  and  $z^{\pm}$  29  
4.2.5 Comments on the screening operator 29

4.3 Rank 3:  $\Gamma = A_3, B_3, H_3$  30

4.3.1 Example:  $\Gamma = A_{3}$  30  
4.3.2 Example:  $\Gamma = B_{3}$  31  
4.3.3 Example:  $\Gamma = H_{3}$  32

4.4 A rank 4 example:  $\Gamma = D_4$  32

# 5 Examples of  $\mathcal{N} = 2$  VOA  $\mathcal{W}_{\mathrm{G}}$  35

5.1 Rank 1:  $\mathsf{G} = \mathbb{Z}_p$  35

5.1.1 Construction of the generators and OPEs 35  
5.1.2 Null states 37  
5.1.3 Screening operator 38

5.2 A rank 2 example:  $\mathsf{G} = G(3,1,2)$  39

# 6 The  $R$  -filtration from free fields 41

6.1 The  $R$  -filtration 43

6.1.1 The  $4d / 2d$  map and Schur operators 43  
6.1.2 Gradings and filtration 44  
6.1.3 The superconformal index 45  
6.1.4 The  $\mathcal{N} = 4$  case 46

6.2 The  $\mathcal{R}$  -filtration 48  
6.3 Identification of  $\mathcal{R}$  -filtration and  $R$  -filtration 50  
6.4 The Hall-Littlewood limit 53

6.4.1 Examples 55

# A Some basic facts on OPEs and VOAs 57

# B Hilbert series and indices 62

B.1 Molien series 62  
B.2 The index in some limits 64

# 1 Introduction and summary

Any four-dimensional  $\mathcal{N} = 2$  superconformal field theory (SCFT) contains a subsector isomorphic to a vertex operator algebra (VOA) [1]. This  $4d / 2d$  correspondence (see [2-20] for some further developments) promises to become an organizing principle for the whole landscape of  $\mathcal{N} = 2$  SCFTs. The aspiration is to combine the rigid associativity constraints of VOAs with additional physical requirements, such as unitary of the  $4d$  theory, in order to constrain and ideally classify the set of  $\mathcal{N} = 2$  SCFTs. $^{1}$  To make inroads into this program, it may be fruitful to start with theories that enjoy enhanced supersymmetry. The map of [1] associates to a generic  $\mathcal{N} = 2$  SCFTs a conformal VOA whose Virasoro subalgebra has central charge  $c = -12c_{4d}$ , where  $c_{4d}$  is the (Weyl) $^{2}$  conformal anomaly coefficient. If the  $4d$  SCFT has  $\mathcal{N} = 3$  or  $\mathcal{N} = 4$  supersymmetry, the associated VOA is supersymmetric, containing an  $\mathcal{N} = 2$  or (small)  $\mathcal{N} = 4$  super Virasoro subalgebra, respectively.

With this broad motivation in mind, we have undertaken a systematic study of the classes of  $\mathcal{N} = 2$  and  $\mathcal{N} = 4$  VOAs that arise from 4d SCFTs. In this paper we describe a uniform construction of all the known examples (and, as we shall see, of additional examples, some of which are unlikely to have a four-dimensional interpretation). Our aim here is descriptive rather than taxonomic. The general classification program outlined above will require different methods, and is left for future work.

A central class of examples are the VOAs associated to the  $\mathcal{N} = 4$  super Yang-Mills (SYM) theories. The  $\mathcal{N} = 4$  4d theory  $\mathrm{SYM}_{\mathfrak{g}}$  with gauge algebra  $\mathfrak{g}$  descends to an  $\mathcal{N} = 4$  VOA  $\chi[\mathrm{SYM}_{\mathfrak{g}}]$  with central charge  $c = -3\dim(\mathfrak{g})$ . As the 4d theory has a Lagrangian description, one can apply the methods of [1] to give a description of the associated VOA, as a subalgebra of  $\dim(\mathfrak{g})$  copies of the  $\beta \gamma bc$  ghost system. The subalgebra is defined by passing to the cohomology of a certain nilpotent operator  $Q_{\mathrm{BRST}}$  built in terms of the  $\beta \gamma bc$  ghosts. While giving in principle a complete definition of the VOA, this description is cumbersome and redundant. The calculation of the requisite BRST cohomology is a difficult problem that so far has only been solved in examples, by brute force level-by-level calculation up to some maximum conformal weight. The simplest example is the VOA associated to  $\mathcal{N} = 4$  SYM theory with gauge algebra  $\mathfrak{sl}(2)$ . There is strong evidence [1] that  $\chi[\mathrm{SYM}_{\mathfrak{sl}(2)}]$  coincides with the small  $\mathcal{N} = 4$  superVirasoro algebra with central charge  $c = -9$ . More generally, for a simple Lie algebra  $\mathfrak{g} \neq \mathfrak{sl}(2)$ ,  $\chi[\mathrm{SYM}_{\mathfrak{g}}]$  is an extension of  $\mathrm{Vir}_{\mathcal{N} = 4}$  obtained by introducing additional strong generators. The list of generators includes one short superprimary of  $\mathrm{Vir}_{\mathcal{N} = 4}$  for each Casimir invariant of  $\mathfrak{g}$ .

A concrete description  $\chi[\mathrm{SYM}_{\mathfrak{g}}]$  as a W-algebra (i.e., in terms of the singular OPE of its strong generators), becomes more and more involved as the rank of  $\mathfrak{g}$  increases. We have found such a W-algebra description for a few low rank cases, see section 4 and 5, but it seems very difficult to find such an explicit presentation in the general case. Instead, the main subject of this paper is a proposal for a novel free-field realization of these VOAs, much simpler and more explicit than the cohomological description of [1]. Our new free-field realization is in terms  $\mathrm{rank}(\mathfrak{g})$  copies of the  $\beta \gamma bc$  ghost system. There is a heuristic understanding of these free fields as corresponding to the low-energy degrees of freedom of the  $4d$  theory at a generic point on its Higgs branch of vacua. This physical picture will be discussed elsewhere [25]: it appears to be much more general, possibly valid for all VOAs that arise from  $\mathcal{N} = 2$  SCFTs.

Our proposal generalizes to all  $\mathfrak{g}$  a construction of Adamovic [26], who exhibited a free-field realization of  $\mathrm{Vir}_{\mathcal{N} = 4}$  with  $c = -9$  (in our framework, the  $\mathfrak{g} = \mathfrak{sl}(2)$  case) in terms of a single  $\beta \gamma bc$  system. The  $\mathrm{Vir}_{\mathcal{N} = 4}$  VOA with  $c = -9$  contains a large set of null vectors, and a remarkable feature of Adamovic's construction is that they identically vanish when expressed in terms of the free fields. In other terms, this is a construction of the simple quotient of the VOA. We find strong evidence that our generalization to  $\chi[\mathrm{SYM}_{\mathfrak{g}}]$  shares the same property.

A notable corollary of our proposal is a compelling conjecture for the “ $R$ -filtration” of this class of VOAs. As we review in detail below, any VOA that descends from a  $4d\mathcal{N} = 2$  SCFT inherits a filtration associated to the  $4dR$ -symmetry quantum number - the details of the cohomological construction of [1] imply that while  $R$  is in general not preserved by the OPE, it can at most decrease. The  $R$ -filtration is of paramount importance in extracting four-dimensional physical information from the VOA, but it is completely hidden in an abstract W-algebra presentation. By contrast, our free-field realizations come equipped with a natural

filtration, which coincides with the  $R$ -filtration in all examples that we have been able to check.

Curiously, we overshoot our initial target, finding a larger set of  $\mathcal{N} = 4$  VOAs than originally expected. We found free-field constructions for  $\mathcal{N} = 4$  VOAs  $\mathcal{W}_{\Gamma}$  labelled by a general Coxeter group  $\Gamma$ . So far, we have described the situation when  $\Gamma$  is the Weyl group  $\mathrm{Weyl}(\mathfrak{g})$  of a simple Lie algebra  $\mathfrak{g}$ , hence both Coxeter and crystallographic. Our basic contention is that  $\mathcal{W}_{\mathrm{Weyl}(\mathfrak{g})} = \chi[\mathrm{SYM}_{\mathfrak{g}}]$ . However, our construction goes through even if  $\Gamma$  is not crystallographic, in which case  $\mathcal{W}_{\Gamma}$  does not have any obvious four-dimensional interpretation. Clearly, it cannot descend from a  $4d$  SYM theory. Even if one is willing to entertain the possibility that the SYM theories do not exhaust the set of  $4d\mathcal{N} = 4$  SCFTs, the conventional wisdom is that none of them can give rise to  $\mathcal{W}_{\Gamma}$  if  $\Gamma$  is a non-crystallographic Coxeter group. Indeed, as we explain below, the moduli space of vacua of the putative parent  $4d$  theory would be the orbifold  $\mathbb{R}^{6n}/\Gamma$ , where  $n$  is the rank of  $\Gamma$ , but general consistency conditions on the low-energy effective theory restrict  $\Gamma$  to be crystallographic [23, 29].

This whole circle of ideas admits a natural extension to a class of VOAs  $\mathcal{W}_{\mathsf{G}} \supset \mathcal{W}_{\Gamma}$ , labelled by a general complex reflection group  $\mathsf{G}$ . These vertex algebras are extensions of the  $\mathcal{N} = 2$  superVirasoro algebra by additional generators, including one short superprimary for each of the fundamental invariants of  $\mathsf{G}$ . Their central charge is fixed in terms of the degrees of the primitive invariants of  $\mathsf{G}$ , see (1.4). We propose a free-field construction of the simple quotient of these algebras in terms of  $\mathrm{rank}(\mathsf{G})$ $\beta \gamma bc$  ghost systems. If (and only if)  $\mathsf{G} = \Gamma$  is a Coxeter group, its lowest fundamental invariant has degree two, and the corresponding short generator of the VOA induces an enhancement of the superconformal algebra from  $\mathcal{N} = 2$  to  $\mathcal{N} = 4$ , so that we recover the construction of  $\mathcal{W}_{\Gamma}$  discussed above. If  $\mathsf{G}$  is a crystallographic complex reflection group (which is not a Coxeter group),  $\mathcal{W}_{\mathsf{G}}$  may descend via the map of [1] from an  $\mathcal{N} = 3$  4d SCFT (which is not an  $\mathcal{N} = 4$  SCFT). Many examples of  $\mathcal{N} = 2$  VOAs associated to  $\mathcal{N} = 3$  SCFTs have been described in the literature [30, 31], and we are able to identify each of them with  $\mathcal{W}_{\mathsf{G}}$  for some choice of  $\mathsf{G}$ . An Euler-Venn diagram of complex reflection groups is presented in figure 1.

In summary, we have found a uniform construction for:

(i) All the VOAs that descend from currently known  $\mathcal{N} = 3$  and  $\mathcal{N} = 4$  4d SCFTs. They are labelled by Weyl groups in the  $\mathcal{N} = 4$  case and by a subset of the crystallographic complex reflection groups in the  $\mathcal{N} = 3$  case.  
(ii) Additional VOAs, labelled by the remaining crystallographic complex reflection groups, which are candidates for  $4d$  uplifts to  $\mathcal{N} = 3$  SCFTs, but whose  $4d$  counterparts are currently unknown.  
(iii) VOAs labelled by non-crystallographic Coxeter and complex reflection groups, which presumably do not correspond to standard  $4d$  theories.

![](images/e322ccdbb42045bebd16592a7d74086542da125bb4619839319c0e1f1cb40326.jpg)  
Figure 1: Euler-Venn diagram depicting the relations among complex reflection groups, crystallographic complex reflection groups, Coxeter groups, and Weyl groups.

Perhaps the most interesting open question is whether abstract bootstrap methods can lead to a rigorous classification of all  $\mathcal{N} \geq 2$  VOAs that descend from  $\mathcal{N} \geq 3$  4d SCFTs. The bootstrap constraints would have to incorporate intrinsic 4d conditions, such as 4d unitarity. The simplest conjecture generalizing all available data is that the complete list of such VOAs coincides with  $\mathcal{W}_{\mathsf{G}}$  with  $\mathsf{G}$  crystallographic, as well as of the additional VOAs obtained from the above set by performing discrete quotients (see, e.g., [29, 32]). We have made some preliminary progress in this direction, which we will report in an upcoming publication [33]. We preview some our findings in the outlook subsection of this introduction.

# 1.1 Main results

In this work we study a class of  $\mathcal{N} = 2$  VOAs labelled by complex reflection groups. By an  $\mathcal{N} = 2$  VOA we mean a VOA that contains the  $\mathcal{N} = 2$  superconformal algebra (SCA) as a subalgebra. The  $\mathcal{N} = 2$  SCA is generated by the stress tensor, together with an affine  $\mathfrak{gl}(1)$  current and two supercurrents. The global part of the  $\mathcal{N} = 2$  SCA is  $\mathfrak{osp}(2|2)$ , and it is natural to organize the operator content of our VOAs into representations of  $\mathfrak{osp}(2|2)$ . We will discuss  $\mathfrak{osp}(2|2)$  representations in more detail in section 2.1 and in appendix A, but for the purposes of the present introduction, it suffices to recall they are labelled by two quantum numbers  $(h,m)$ , which are the conformal dimension and the  $\mathfrak{gl}(1)$  charge of the highest weight state in the multiplet. Representations with  $h = \pm m$  obey a shortening condition, and are referred to as chiral, anti-chiral, respectively. A representation with  $h \neq |m|$  will be referred to as non-chiral.

The class of VOAs that we analyze in this work are extensions of the  $\mathcal{N} = 2$  SCA obtained by introducing additional strong generators. The additional strong generators, as well as several other interesting properties of the VOAs under examination, are intimately related to the theory of invariants of the complex reflection group  $\mathsf{G}$ . Therefore, before we proceed, we need to introduce some algebraic structures related to complex reflection groups.

Let  $\mathsf{G}$  be a complex reflection group, regarded as a subgroup of  $GL(V_{\mathsf{G}})$  with  $V_{\mathsf{G}} \cong \mathbb{C}^{\mathrm{rank}(\mathsf{G})}$ . According to the Chevalley-Shephard-Todd theorem, the ring of invariants  $\mathbb{C}[V_{\mathsf{G}}]^{\mathsf{G}}$  is a freely generated polynomial ring. The generators of the ring  $\mathbb{C}[V_{\mathsf{G}}]^{\mathsf{G}}$  are usually referred to as the

fundamental invariants of  $\mathsf{G}$ . Their number equals the rank of  $\mathsf{G}$ , and their degrees are denoted  $p_{\ell}$ , with  $\ell = 1, \ldots, \mathrm{rank}(\mathsf{G})$ . For example, the rank-one complex reflection group  $\mathsf{G} = \mathbb{Z}_p$  has a unique fundamental invariant of degree  $p_1 = p$ . If we introduce the coordinate  $z$  on  $V_{\mathbb{Z}_p} \cong \mathbb{C}$ , the action of  $\mathbb{Z}_p$  on  $V_{\mathbb{Z}_p}$  is simply  $z \mapsto e^{2\pi i / p} z$ , and the ring of invariants  $\mathbb{C}[V_{\mathbb{Z}_p}]^{\mathbb{Z}_p}$  is freely generated by  $z^p$ .

For our purposes, we need to consider the canonical symplectic variety associated to the action of  $\mathsf{G}$  on  $V_{\mathsf{G}}$ . More precisely, let us define the variety

$$
\mathcal {M} _ {\mathrm {G}} = \frac {V _ {\mathrm {G}} \oplus V _ {\mathrm {G}} ^ {*}}{\mathrm {G}}, \tag {1.1}
$$

where  $V_{\mathsf{G}}^{*}$  denotes the dual of the vector space  $V_{\mathsf{G}}$ . The space  $V_{\mathsf{G}} \oplus V_{\mathsf{G}}^{*}$  admits a canonical symplectic structure, which is preserved by the action of  $\mathsf{G}$ . We also define the ring associated to  $\mathcal{M}_{\mathsf{G}}$ ,

$$
\mathcal {R} _ {\mathsf {G}} = \mathbb {C} [ \mathcal {M} _ {\mathsf {G}} ] = \mathbb {C} \left[ V _ {\mathsf {G}} \oplus V _ {\mathsf {G}} ^ {*} \right] ^ {\mathsf {G}}. \tag {1.3}
$$

As opposed to the ring  $\mathbb{C}[V_{\mathsf{G}}]^{\mathsf{G}}$ , the ring  $\mathcal{R}_{\mathsf{G}}$  is not freely generated. For instance, in the case  $\mathsf{G} = \mathbb{Z}_p$ , we may parametrize the  $\mathsf{G}$  action on  $V_{\mathsf{G}} \oplus V_{\mathsf{G}}^{*}$  as  $(z,\bar{z}) \mapsto (e^{2\pi i / p}z,e^{-2\pi i / p}\bar{z})$ . The ring  $\mathcal{R}_{\mathsf{G}}$  is then generated by the monomials  $\mathfrak{j} = z\bar{z}$ ,  $\mathsf{w} = z^p$ ,  $\bar{\mathsf{w}} = \bar{z}^p$ , subject to the relation  $\mathsf{w}\bar{\mathsf{w}} = \mathfrak{j}^p$ .

As  $\mathcal{M}_{\mathsf{G}}$  admits a natural  $GL(1)\times GL(1)$  action, the ring  $\mathcal{R}_{\mathsf{G}}$  can always be given a presentation in terms of generators and relations, each possessing definite quantum numbers  $(h,m)$  under the  $GL(1)\times GL(1)$  action. The generators with  $h = m$  are in 1-to-1 correspondence with the fundamental invariants of  $\mathsf{G}$ . Our normalization of the quantum numbers  $h$ ,  $m$  is such that  $h = m = p_{\ell} / 2$  for the  $\ell$ -th fundamental invariant of  $\mathsf{G}$ . For any complex reflection group  $\mathsf{G}$ , there is a unique generator with  $(h,m) = (1,0)$ . In the example of  $\mathsf{G} = \mathbb{Z}_p$ , the quantum numbers of the generators  $\mathbf{j},\mathbf{w},\bar{\mathbf{w}}$  are  $(1,0)$ ,  $(\frac{p}{2},\frac{p}{2})$ ,  $(\frac{p}{2}, - \frac{p}{2})$ , respectively. We refer the reader to section 2.2 for more information on the general case.

We are now in a position to articulate our main proposal. Given a complex reflection group  $\mathsf{G}$ , we claim that there exists an  $\mathcal{N} = 2$  VOA, denoted  $\mathcal{W}_{\mathsf{G}}$ , satisfying the properties (i)-(vii) listed below. Before giving the list of properties, we would like to caution the reader that we do not have a general existence proof for  $\mathcal{W}_{\mathsf{G}}$ , but we have nonetheless gathered a considerable amount of evidence in favor of our proposal. In particular, in the case in which the complex reflection group is the Weyl group of a semisimple Lie algebra, one can resort to the BRST construction of [1] to demonstrate the existence of  $\mathcal{W}_{\mathsf{G}}$ . Beyond Weyl groups, we have a fully explicit construction of  $\mathcal{W}_{\mathsf{G}}$  for the infinite series  $\mathsf{G} = I_2(p)$  and  $\mathsf{G} = \mathbb{Z}_p$ , as well as for all rank-three Coxeter groups, and the rank-two complex reflection group  $G(3,1,2)$  (whose definition is recalled in section 5.2). We now give the list of properties enjoyed by  $\mathcal{W}_{\mathsf{G}}$ :

$$
\omega \left(x _ {1} \oplus \xi_ {1}, x _ {2} \oplus \xi_ {2}\right) = \xi_ {2} \left(x _ {1}\right) - \xi_ {1} \left(x _ {2}\right), \tag {1.2}
$$

(i)  $\mathcal{W}_{\mathsf{G}}$  is a simple  $\mathcal{N} = 2$  VOA that is strongly generated by a finite number of operators.7 All strong generators of  $\mathcal{W}_{\mathsf{G}}$  are organized in multiplets of  $\mathfrak{osp}(2|2)$ .  
(ii) Chiral and anti-chiral multiplets of strong generators of  $\mathcal{W}_{\mathsf{G}}$  come in pairs, which are in 1-to-1 correspondence with the fundamental invariants of  $\mathsf{G}$ . The quantum numbers of the chiral/anti-chiral pairs are  $(h,m) = (p_{\ell} / 2,\pm p_{\ell} / 2)$ , with  $\ell = 1,\ldots ,\mathrm{rank}(\mathsf{G})$ , where  $p_{\ell}$  are the degrees of the fundamental invariants of  $\mathsf{G}$ .  
(iii) For each generator of the ring  $\mathcal{R}_{\mathsf{G}}$  with  $GL(1)\times GL(1)$  quantum numbers  $(h,m)$ , there is an  $\mathfrak{osp}(2|2)$  multiplet of strong generators of  $\mathcal{W}_{\mathsf{G}}$ , labelled by the same quantum numbers. Because of point (ii), this observation is trivial if  $h = \pm m$ , but it is non-trivial for generators of  $\mathcal{R}_{\mathsf{G}}$  with  $h\neq |m|$ , which are mapped to non-chiral multiplets of strong generators of  $\mathcal{W}_{\mathsf{G}}$ . The unique generator of  $\mathcal{R}_{\mathsf{G}}$  with  $(h,m) = (1,0)$  is mapped to the set of generators of the  $\mathcal{N} = 2$  SCA.

(iv) The central charge of  $\mathcal{W}_{\mathbf{G}}$  is given by

$$
c = - 3 \sum_ {\ell = 1} ^ {\operatorname {r a n k} (\mathbf {G})} (2 p _ {\ell} - 1), \tag {1.4}
$$

where  $p_{\ell}$  are the degrees of the fundamental invariants of  $\mathsf{G}$ .

(v)  $\mathcal{W}_{\mathsf{G}}$  can be realized as a subalgebra of  $\mathrm{rank}(\mathsf{G})$  copies of a standard  $\beta \gamma bc$  system. We may write

$$
\mathcal {W} _ {\mathbf {G}} \subset \mathbb {M} _ {\beta \gamma b c} ^ {(\mathbf {G})} := \bigotimes_ {\ell = 1} ^ {\operatorname {r a n k} (\mathbf {G})} M _ {\beta \gamma b c} ^ {(p _ {\ell})}. \tag {1.5}
$$

The symbol  $M_{\beta \gamma bc}^{(p_\ell)}$  denotes a standard  $\beta \gamma bc$  system, endowed with its canonical  $\mathcal{N} = 2$  superconformal structure labelled by the quantum number  $p_\ell$ . More explicitly, each tensor factor  $M_{\beta \gamma bc}^{(p_\ell)}$  comes with a free-field realization of the  $\mathcal{N} = 2$  SCA, in which the quantum numbers of the free fields in  $M_{\beta \gamma bc}^{(p_\ell)}$  are determined by  $p_\ell$ , see table (3.4). For  $\mathsf{G} = \mathbb{Z}_2$ , the free-field realization (1.5) of  $\mathcal{W}_{\mathsf{G}}$  coincides with the one given in [26].

Crucially, we claim that the above free-field realization is a realization of  $\mathcal{W}_{\mathsf{G}}$  regarded as simple quotient of the span of the strong generators. More explicitly, we propose that all null states built from the strong generators are automatically identically zero in the free-field realization. This has been proven for  $\mathsf{G} = \mathbb{Z}_2$  in [26].

(vi) Thanks to the free-field realization of point (v) above, we have a natural way to associate a Poisson algebra to  $\mathcal{W}_{\mathsf{G}}$ . This is achieved by means of the map

$$
\mathcal {P}: \mathcal {W} _ {\mathsf {G}} \rightarrow \mathbb {M} _ {\beta \gamma} ^ {(\mathsf {G}) \mathrm {c l}}, \tag {1.6}
$$

where  $\mathbb{M}_{\beta \gamma}^{(\mathsf{G})\mathrm{cl}}$  denotes the Poisson algebra of polynomials in the indeterminates  $\beta_{\ell}$ ,  $\gamma_{\ell}$  with  $\ell = 1, \ldots, \operatorname{rank}(\mathsf{G})$ , with each pair  $(\beta_{\ell}, \gamma_{\ell})$  regarded as a pair of canonically conjugate variables. The image of an operator  $\mathcal{O}$  in  $\mathcal{W}_{\mathsf{G}}$  under the map  $\mathcal{P}$  is obtained starting from the free-field realization of  $\mathcal{O}$  and setting to zero all derivatives and all factors of  $b_{\ell}$ ,  $c_{\ell}$ ,  $\ell = 1, \ldots, \operatorname{rank}(\mathsf{G})$ . More details on the map  $\mathcal{P}$  are found in section 3.4. The image of  $\mathcal{W}_{\mathsf{G}}$  under the map  $\mathcal{P}$  is a Poisson subalgebra  $\mathcal{P}(\mathcal{W}_{\mathsf{G}})$  of  $\mathbb{M}_{\beta \gamma}^{(\mathsf{G})\mathrm{cl}}$ . We propose that  $\mathcal{P}(\mathcal{W}_{\mathsf{G}})$  is isomorphic to the ring of invariants  $\mathcal{R}_{\mathsf{G}}$ , regarded as a Poisson algebra,

$$
\mathcal {P} \left(\mathcal {W} _ {\mathrm {G}}\right) \cong \mathcal {R} _ {\mathrm {G}}. \tag {1.7}
$$

(vii) Following the general construction of [34], one can define the vector subspace  $C_2(\mathcal{W}_{\mathsf{G}}) \subset \mathcal{W}_{\mathsf{G}}$  and the commutative Zhu algebra  $\mathcal{W}_{\mathsf{G}} / C_2(\mathcal{W}_{\mathsf{G}})$ . The latter is naturally endowed with the structure of a Poisson algebra. Crucially,  $\mathcal{W}_{\mathsf{G}} / C_2(\mathcal{W}_{\mathsf{G}})$  needs not be a reduced ring, i.e.  $\mathcal{W}_{\mathsf{G}} / C_2(\mathcal{W}_{\mathsf{G}})$  may contain nilpotent elements. One can mod out nilpotent elements of  $\mathcal{W}_{\mathsf{G}} / C_2(\mathcal{W}_{\mathsf{G}})$  and obtain a reduced ring, denoted  $(\mathcal{W}_{\mathsf{G}} / C_2(\mathcal{W}_{\mathsf{G}}))_{\mathrm{red}}$ . The latter is still a Poisson algebra. We propose the following Poisson algebra isomorphism,

$$
\left(\mathcal {W} _ {\mathrm {G}} / C _ {2} \left(\mathcal {W} _ {\mathrm {G}}\right)\right) _ {\text {r e d}} \cong \mathscr {R} _ {\mathrm {G}}. \tag {1.8}
$$

In particular, this means that the associated variety to  $\mathcal{W}_{\mathsf{G}}$ , in the sense of [34], is precisely the variety  $\mathcal{M}_{\mathsf{G}}$  defined in (1.1).

In the case of a complex reflection group that is also a Coxeter group, the structure of the corresponding VOA is richer. In what follows, we reserve the symbol  $\Gamma$  for a Coxeter group.

First of all, supersymmetry is enhanced from  $\mathcal{N} = 2$  to small  $\mathcal{N} = 4$ . The small  $\mathcal{N} = 4$  SCA is generated by the stress tensor, an  $\mathfrak{sl}(2)$  triplet of affine currents, and two  $\mathfrak{sl}(2)$  doublets of supercurrents. The explicit embedding of the  $\mathcal{N} = 2$  SCA into the small  $\mathcal{N} = 4$  SCA is given in (2.8). The global part of the small  $\mathcal{N} = 4$  SCA is  $\mathfrak{psl}(2|2)$ . Multiplets of  $\mathfrak{psl}(2|2)$  are labelled by a pair  $(h,j)$ , where  $h$  is the conformal dimension and  $j$  is the  $\mathfrak{sl}(2)$  spin of the highest weight state. Multiplets with  $h = j$  obey a shortening condition, and are referred to as short in what follows. Generic multiplets with  $h > j$  are referred to as long multiplets. We refer the reader to section 2.1 and appendix A for more information.

Recall that Coxeter groups are characterized within complex reflection groups by the property of possessing exactly one fundamental invariant of degree two. By virtue of point (ii) above, it follows that the VOA  $\mathcal{W}_{\Gamma}$  associated to a Coxeter group  $\Gamma$  has exactly one chiral/anti-chiral pair of strong generators with  $h = \pm m = 1$ . This chiral/anti-chiral pair combines with the set of generators of the  $\mathcal{N} = 2$  SCA to give a short multiplet of  $\mathfrak{psl}(2|2)$  that encompasses all generators of the small  $\mathcal{N} = 4$  SCA.

This phenomenon extends to all chiral/anti-chiral pair of multiplets of strong generators of  $\mathcal{W}_{\Gamma}$ . More precisely, each chiral/anti-chiral pair in  $\mathcal{N} = 2$  language is paired with non-chiral  $\mathcal{N} = 2$  multiplets to give a short multiplet of  $\mathfrak{psl}(2|2)$ . Crucially, not all  $\mathcal{N} = 2$  non-chiral multiplets are necessarily paired with chiral/anti-chiral pairs in the enhancement to small  $\mathcal{N} = 4$  supersymmetry. As a result, the VOA  $\mathcal{W}_{\Gamma}$  generically admits long multiplets of strong generators.

The supersymmetry enhancement at the level of the VOA is mirrored by an enhancement of the isometry group of the variety  $\mathcal{M}_{\mathsf{G}}$ . Indeed, if  $\mathsf{G}$  is a Coxeter group  $\Gamma$ , the variety (1.1) admits an alternative presentation,

$$
\mathcal {M} _ {\Gamma} = \frac {\mathbb {C} ^ {2} \otimes V _ {\Gamma} ^ {\mathbb {R}}}{\Gamma}, \tag {1.9}
$$

where  $\Gamma$  acts trivially on the  $\mathbb{C}^2$  factor. In the above expression, we have implicitly identified the Coxeter group  $\Gamma$  with a subgroup of  $O(V_{\Gamma}^{\mathbb{R}})$ , where  $V_{\Gamma}^{\mathbb{R}} \cong \mathbb{R}^{\mathrm{rank}(\Gamma)}$  is a real vector space. $^9$  The presentation (1.9) makes it manifest that the  $GL(1) \times GL(1)$  action is enhanced to a  $GL(2) \cong GL(1) \times SL(2)$  action. It follows that the ring  $\mathcal{R}_{\Gamma}$  associated to  $\mathcal{M}_{\Gamma}$  admits a presentation in terms of generators and relations with a definite  $GL(1)$  weight  $h$  and  $SL(2)$  spin  $j$ . The simplest example of Coxeter group is  $\Gamma = \mathbb{Z}_2$ . In this case, the generators  $\mathbf{j}, \mathbf{w}, \bar{\mathbf{w}}$  (introduced above for  $\mathbb{Z}_p$ ) all have the same weight and form a triplet of  $SL(2)$ . The relation  $\mathsf{w}\bar{\mathsf{w}} = \mathsf{j}^2$  is a singlet of  $SL(2)$ .

Let us summarize how points (i)-(iii) above can be rephrased in the Coxeter case.

(i)'  $\mathcal{W}_{\Gamma}$  is a simple VOA with small  $\mathcal{N} = 4$  supersymmetry that is strongly generated by a finite number of operators. All strong generators of  $\mathcal{W}_{\Gamma}$  are organized in multiplets of  $\mathfrak{psl}(2|2)$ .  
(ii)' Short multiplets of strong generators of  $\mathcal{W}_{\Gamma}$  are in 1-to-1 correspondence with the fundamental invariants of  $\Gamma$ . The quantum numbers of the short multiplets of strong generators are  $h = j = p_{\ell} / 2$ , with  $\ell = 1, \ldots, \mathrm{rank}(\mathsf{G})$ , where  $p_{\ell}$  are the degrees of the fundamental invariants of  $\Gamma$ . The short multiplet containing the generators of the small  $\mathcal{N} = 4$  SCA is in correspondence with the fundamental invariant of  $\Gamma$  of degree two.  
(iii)' Upon expressing generators and relations of the ring  $\mathcal{R}_{\Gamma}$  in a  $GL(2)$  covariant way, for each  $GL(2)$  multiplet of generators of  $\mathcal{R}_{\Gamma}$  with quantum numbers  $(h,j)$ , there is a  $\mathfrak{psl}(2|2)$  multiplet of strong generators of  $\mathcal{W}_{\Gamma}$ , labelled by the same quantum numbers.

Because of point (ii)' above, this observation is trivial if  $h = j$ , but it is non-trivial for generators of  $\mathcal{R}_{\Gamma}$  with  $h \neq j$ , which are mapped to long multiplets of strong generators of  $\mathcal{W}_{\Gamma}$ .

The points (iv) to (vii) in the list of properties in the complex reflection case are specialized to the Coxeter case with obvious modifications. In particular, we still have a free-field realization of  $\mathcal{W}_{\Gamma}$ , in which all null states are automatically zero.

# 1.2 Connection with 4d physics

The class of VOAs labelled by complex reflection groups that we have described in the previous section constitutes a rich and novel set of VOAs, worth studying in its own right. We now come back to our original motivation for the study of these  $2d$  algebraic structures, which is rooted in the analysis of 4d SCFTs. We propose that the class of VOAs  $\mathcal{W}_{\mathsf{G}}$  associated to complex reflection groups  $\mathsf{G}$  provides a unified framework to describe VOAs originating from 4d SCFT with  $\mathcal{N} \geq 3$  via the map of [1].

In a generic  $4d\mathcal{N} = 2$  SCFT, there is no obvious relation between the geometries of the Coulomb and Higgs branches. If the theory has  $\mathcal{N}\geq 3$  supersymmetry, however, the geometry of the Higgs branch is modeled after the geometry of the Coulomb branch. This is a consequence of the fact that both branches are subspaces of the full moduli space of the theory, which is highly constrained by  $\mathcal{N}\geq 3$  supersymmetry.

Let us explain this point in more detail. The moduli space of a  $4d\mathcal{N} = 3$  theory is parametrized by the expectation value of scalars in  $\mathcal{N} = 3$  vector multiplets, which are identical to  $\mathcal{N} = 4$  vector multiplets. Therefore, the moduli space of a  $4d\mathcal{N}\geq 3$  has real dimension  $6r$ , where  $r$  is the rank of the SCFT, and is locally flat. Upon selecting an  $\mathcal{N} = 2$  subalgebra of the  $\mathcal{N}\geq 3$  superconformal symmetry, the moduli space geometry is reformulated in terms of a Coulomb branch and a Higgs branch. Let us assume that the Coulomb branch can be written globally as a quotient of  $\mathbb{C}^r$  by a discrete group,  $\mathbb{C}^r /G$ . The analysis of [23, 29] then reveals that  $G$  must be a crystallographic complex reflection group, whose rank equals the rank  $r$  of the SCFT. Moreover, if the Coulomb branch is  $\mathbb{C}^r /G$ ,  $\mathcal{N}\geq 3$  supersymmetry implies that the Higgs branch must coincide with the variety (1.1) (recall  $V_{\mathsf{G}}\cong \mathbb{C}^{\mathsf{r}}$ ).

Let us now consider the VOA associated to the SCFT. According to the conjecture of [13], the associated variety to the VOA must coincide with the Higgs branch of the  $4d$  theory, which is the variety (1.1). Moreover, we know from [1] that the  $2d$  central charge  $c_{2d}$  must be given by  $c_{2d} = -12c_{4d}$ , where  $c_{4d}$  is the  $(\mathrm{Weyl})^2$  trace anomaly coefficient. On the other hand, in any  $4d\mathcal{N}\geq 3$  SCFT, the two trace anomaly coefficients are equal,  $a_{4d} = c_{4d}$  [35]. Furthermore, we can use the Shapere-Tachikawa formula [36, 37] to relate them to the dimensions of the Coulomb branch generators  $D_{\ell}$ ,  $\ell = 1,\ldots ,\mathbf{r}$ ,

$$
c _ {4 \mathrm {d}} = \frac {1}{4} \sum_ {\ell = 1} ^ {\mathrm {r}} \left(2 D _ {\ell} - 1\right). \tag {1.11}
$$

By assumption, the Coulomb branch is  $\mathbb{C}^{\mathbf{r}} / \mathsf{G}$ , which implies that the dimensions of the Coulomb generators coincide with the degrees of the fundamental invariants of  $\mathsf{G}$ ,  $D_{\ell} = p_{\ell}$ ,  $\ell = 1,\ldots ,\mathsf{r}$ . It follows that the central charge of the VOA is given by the formula (1.4).

The above considerations provide strong evidence that, if we start with a  $4d\mathcal{N}\geq 3$  SCFT with Coulomb branch  $\mathbb{C}^r /\mathbb{G}$ , with  $\mathsf{G}$  a crystallographic complex reflection group of rank  $\mathbf{r}$ , the associated VOA is precisely the VOA  $\mathcal{W}_{\mathsf{G}}$  described in section 1.1. Conversely, it seems natural to expect that, for any crystallographic complex reflection group  $\mathsf{G}$ , the VOA  $\mathcal{W}_{\mathsf{G}}$  should originate from a  $4d\mathcal{N}\geq 3$  SCFT.

In table 1 we list all irreducible crystallographic complex reflection groups. $^{10}$  If  $\mathsf{G}$  is a crystallographic complex reflection group that is also Coxeter, i.e. a Weyl group, the identification of the parent  $4d$  theory for the VOA  $\mathcal{W}_{\mathsf{G}}$  is straightforward: it is simply  $4d\mathcal{N} = 4$  SYM with the appropriate gauge algebra  $\mathfrak{g}$ , such that  $\mathsf{G} = \mathrm{Weyl}(\mathfrak{g})$ . If we consider a crystallographic complex reflection group that is not Coxeter, the putative  $4d$  parent theory in the sense of [1] should be a genuine  $\mathcal{N} = 3$  SCFT. For some of the entries in table 1 the parent  $\mathcal{N} = 3$  theory has been identified in [30, 38, 39]. $^{11}$  It would be interesting to confirm or rule out the existence of parent  $\mathcal{N} = 3$  theories for all other entries of table 1.

Any VOA arising from a  $4d$  SCFT via the construction of [1] is equipped with the so-called "R-filtration" [13], which originates from the  $\mathfrak{sl}(2)_R$  symmetry of the parent  $4d\mathcal{N}\geq 2$  SCFT. While natural from a  $4d$  perspective, the  $R$ -filtration does not seem to be intrinsic to the VOA in any obvious way. For instance, given an abstract presentation of the VOA in terms of its strong generators and the singular OPEs among them, it is not clear in general how to recover the  $R$ -filtration. This is a pressing problem, since the  $R$ -filtration is instrumental in achieving a detailed understanding of the map from  $4d$  to  $2d$  operators, which in turn is pivotal in many applications of the VOA technology to  $4d$  physics.

In this work, we propose a simple solution to the problem of the  $R$ -filtration for all VOAs that admit a  $4d$  origin via the map of [1], and at the same time can be identified with one of the VOAs  $\mathcal{W}_{\mathsf{G}}$  for some complex reflection group  $\mathsf{G}$ . In this case we can utilize the free-field construction to define a filtration of  $\mathcal{W}_{\mathsf{G}}$ , which will be referred to as  $\mathcal{R}$ -filtration. This filtration is specified in a simple way by assigning weights to the free field "letters"  $\beta_{\ell}, \gamma_{\ell}, b_{\ell}, c_{\ell}, \ell = 1, \ldots, \mathrm{rank}(\mathsf{G})$ , see (6.22). We propose the identification of this novel  $\mathcal{R}$ -filtration with the sought-for  $R$ -filtration, and we perform several tests of this proposal.

The identification of  $R$ -filtration and  $\mathcal{R}$ -filtration allows us to recover the Macdonald limit of the superconformal index of the parent  $4d$  theory from the corresponding VOA. This is achieved with the refinement the vacuum character of the VOA by the  $\mathcal{R}$ -filtration.

<table><tr><td colspan="3">Non-Coxeter groups</td></tr><tr><td>Group</td><td>Rank</td><td>Degrees</td></tr><tr><td>G(3,1,n)</td><td>n</td><td>3,6,...,3n</td></tr><tr><td>G(3,3,n), n ≥ 3</td><td>n</td><td>3,6,...,3(n-1); n</td></tr><tr><td>G(4,1,n)</td><td>n</td><td>4,8,...,4n</td></tr><tr><td>G(4,2,n)</td><td>n</td><td>4,8,...,4(n-1); 2n</td></tr><tr><td>G(4,4,n), n ≥ 3</td><td>n</td><td>4,8,...,4(n-1); n</td></tr><tr><td>G(6,1,n)</td><td>n</td><td>6,12,...,6n</td></tr><tr><td>G(6,2,n)</td><td>n</td><td>6,12,...,6(n-1); 3n</td></tr><tr><td>G(6,3,n)</td><td>n</td><td>6,12,...,6(n-1); 2n</td></tr><tr><td>G(6,6,n), n ≥ 3</td><td>n</td><td>6,12,...,6(n-1); n</td></tr><tr><td>Zk, k ∈ {3,4,6}</td><td>1</td><td>k</td></tr><tr><td>G4</td><td>2</td><td>4,6</td></tr><tr><td>G5</td><td>2</td><td>6,12</td></tr><tr><td>G8</td><td>2</td><td>8,12</td></tr><tr><td>G12</td><td>2</td><td>6,8</td></tr><tr><td>G24</td><td>3</td><td>4,6,14</td></tr><tr><td>G25</td><td>3</td><td>6,9,12</td></tr><tr><td>G26</td><td>3</td><td>6,12,18</td></tr><tr><td>G29</td><td>4</td><td>4,8,12,20</td></tr><tr><td>G31</td><td>4</td><td>8,12,20,24</td></tr><tr><td>G32</td><td>4</td><td>12,18,24,30</td></tr><tr><td>G33</td><td>5</td><td>4,6,10,12,18</td></tr><tr><td>G34</td><td>6</td><td>6,12,18,24,30,42</td></tr></table>

Table 1: Irreducible crystallographic complex reflection groups, partitioned into non-Coxeter and Coxeter groups. For each group, we give the rank and the degrees of the fundamental invariants. A crystallographic Coxeter group is a Weyl group: in this case we also include the corresponding Lie algebra(s). The notations  $G(m,p,n)$  and  $G_{n}$  refer to the original list by Shephard and Todd. The symbol  $S_{n}$  denotes the symmetric group of permutations of  $n$  objects. Unless otherwise stated, it is understood that  $n \geq 2$ . The specifications  $n \geq 3$  exclude  $G(3,3,2) \cong \mathrm{Weyl}(\mathfrak{a}_2)$ ,  $G(6,6,2) \cong \mathrm{Weyl}(\mathfrak{g}_2)$ , and  $G(4,4,2)$ , which is conjugate in  $U(2)$  to  $\mathrm{Weyl}(\mathfrak{b}_2)$ . The shaded entries correspond to crystallographic complex reflection groups that govern the Coulomb branch geometry of known  $4d\mathcal{N} \geq 3$  SCFTs.

<table><tr><td colspan="4">Coxeter groups</td></tr><tr><td>Group</td><td>g</td><td>Rank</td><td>Degrees</td></tr><tr><td>Sn</td><td>an-1</td><td>n-1</td><td>2,3,...,n</td></tr><tr><td>G(2,1,n)</td><td>bn,cn</td><td>n</td><td>2,4,...,2n</td></tr><tr><td>G(2,2,n)</td><td>dn</td><td>n</td><td>2,4,...,2(n-1);n</td></tr><tr><td>G(6,6,2)</td><td>g2</td><td>2</td><td>2,6</td></tr><tr><td>G28</td><td>f4</td><td>4</td><td>2,6,8,12</td></tr><tr><td>G35</td><td>e6</td><td>6</td><td>2,5,6,8,9,12</td></tr><tr><td>G36</td><td>e7</td><td>7</td><td>2,6,8,10,12,14,18</td></tr><tr><td>G37</td><td>e8</td><td>8</td><td>2,8,12,14,18,20,24,30</td></tr></table>

# 1.3 Outlook

The main goal of this work is to provide a unified description of a large class of supersymmetric VOAs, which includes all known VOAs originating from 4d SCFTs with  $\mathcal{N} \geq 3$  via the map of [1]. A full classification of all VOAs with  $\mathcal{N} = 2$  or small  $\mathcal{N} = 4$  supersymmetry is a more ambitious task, which would require different tools. One way to address the classification problem consists in posing and trying to solve suitable supersymmetric VOA bootstrap problems. We plan to report progress in this direction in an upcoming paper [33].

A simple example of VOA bootstrap problem is the following. Consider a VOA with small  $\mathcal{N} = 4$  supersymmetry, and suppose it is strongly generated by the generators of the small  $\mathcal{N} = 4$  SCA, together with additional strong generators, organized in a single short  $\mathfrak{psl}(2|2)$  multiplet with prescribed quantum numbers  $h = j = p / 2$ .<sup>12</sup> We may then ask: for a given  $p \geq 3$ , for which values of the central charge  $c$  does the VOA exist? To address this question, one tries to determine the singular OPEs of the strong generators in such a way that all axioms of a bona fide VOA are satisfied. The outcome of this analysis is that the VOA exists in three cases:

(A):  $c = -6(p + 1)$ , for any integer  $p \geq 3$ ;  
(B):  $c = -6(\frac{1}{2} p + 1)$ , for even  $p \geq 4$ ; (1.12)  
(C):  $c = 3(p - 2)$ , for odd  $p \geq 3$ .

The interpretation of case (A) is clear: the VOA is the VOA  $\mathcal{W}_{I_2(p)}$  associated to the Coxeter group  $I_{2}(p)$ , which is discussed in detail in section 4.2 below. The VOAs of cases (B) and (C) cannot be identified with any VOA  $\mathcal{W}_{\Gamma}$  with  $\Gamma$  Coxeter group. We do have, however, an interpretation for case (B) in terms of a quotient of the VOA  $\mathcal{W}_{I_2(p/2)}$  associated to the Coxeter group  $I_{2}(p/2)$ . More precisely, the VOA  $\mathcal{W}_{I_2(p/2)}$  admits a  $\mathbb{Z}_2$  automorphism, and we propose the identification of the VOA of case (B) with  $\mathcal{W}_{I_2(p/2)}/\mathbb{Z}_2$ . A first trivial check of this proposal is the value of the central charge; more non-trivial checks can be performed by analyzing null states in the VOA. As far as the VOA of case (C) is concerned, we restrict ourselves to the simple observation that, since its central charge is positive, it cannot originate from any unitary  $4d$  SCFT via the map of [1].

Another simple VOA bootstrap problem we can address is the following. Consider a small  $\mathcal{N} = 4$  VOA, which is by assumption strongly generated by the generators of the small  $\mathcal{N} = 4$  SCA, together with additional generators, organized in two short  $\mathfrak{psl}(2|2)$  multiplets with given quantum numbers  $h = j = p_1 / 2$  and  $h = j = p_2 / 2$ .<sup>13</sup> Once again, we imagine to fix  $p_1 \leq p_2$ , and we investigate if there is any value of  $c$  for which the VOA exists. We find

that the VOA exists only in four cases:

(a):  $(p_1, p_2) = (3, 4)$ ,  $c = -36$ ;  
(b):  $(p_1, p_2) = (4, 6)$ ,  $c = -54$ ; (1.13)  
(c):  $(p_1, p_2) = (6, 10)$ ,  $c = -90$ ;  
(d):  $(p_1, p_2) = (4, 6)$ ,  $c = -36$ .

We can interpret all these four cases in terms of VOAs associated to a Coxeter group. The VOAs in the cases (a), (b), (c) are the VOAs  $\mathcal{W}_{A_3}$ ,  $\mathcal{W}_{B_3}$ ,  $\mathcal{W}_{H_3}$ , respectively. The VOA of case (d) is the quotient  $\mathcal{W}_{A_3} / \mathbb{Z}_2$ .

As a final comment, we have studied similar VOA bootstrap problems involving more strong generators. The picture emerging from the bootstrap analysis is compatible with the expectation that the VOAs labelled by Coxeter groups introduced in this work exhaust the complete list of (small)  $\mathcal{N} = 4$  W-algebras with certain "good" properties. The fundamental question of identifying such properties will be addressed in [33].

# 2 Preliminaries

In the first part of this section, we briefly review some basic features of the  $\mathcal{N} = 2$  and small  $\mathcal{N} = 4$  SCAs and of their representation theory. In the second part of this section, we collect some standard material on Coxeter and complex reflection groups. In particular, we describe their ring of invariants in the cases of interest for applications in the rest of the paper.

Unless otherwise stated, all Lie (super)algebras in this work are understood to be Lie (super)algebras over the complex numbers.

# 2.1  $\mathcal{N} = 2$  and small  $\mathcal{N} = 4$  SCAs

The small  $\mathcal{N} = 4$  super-Virasoro algebra is generated by affine  $\mathfrak{sl}(2)$  currents  $J^{0,\pm}$ , a stress tensor  $T$ , and four fermionic operators  $\widetilde{G}^{\pm}$  and  $G^{\pm}$ . The bosonic sub-VOA has OPEs

$$
J ^ {0} \left(z _ {1}\right) J ^ {0} \left(z _ {2}\right) \sim \frac {2 k}{\left(z _ {1} - z _ {2}\right) ^ {2}}, \tag {2.1}
$$

$$
J ^ {0} \left(z _ {1}\right) J ^ {\pm} \left(z _ {2}\right) \sim \frac {\pm 2 J ^ {\pm}}{\left(z _ {1} - z _ {2}\right)}, \tag {2.2}
$$

$$
J ^ {+} \left(z _ {1}\right) J ^ {-} \left(z _ {2}\right) \sim \frac {- k}{\left(z _ {1} - z _ {2}\right) ^ {2}} + \frac {- J ^ {0}}{\left(z _ {1} - z _ {2}\right)}, \tag {2.3}
$$

$$
T \left(z _ {1}\right) T \left(z _ {2}\right) \sim \frac {c / 2}{\left(z _ {1} - z _ {2}\right) ^ {4}} + \frac {2 T \left(z _ {2}\right)}{\left(z _ {1} - z _ {2}\right) ^ {2}} + \frac {\partial T \left(z _ {2}\right)}{\left(z _ {1} - z _ {2}\right)}, \tag {2.4}
$$

the OPE between  $T$  and  $J$  express the fact that  $J^{0,\pm}$  are Virasoro primaries of conformal dimension  $h = 1$ , see appendix A. The level and the central charge are related as

$$
c = 6 k. \tag {2.5}
$$

The fermionic generators are Virasoro and affine Kac-Moody (AKM) primaries with weights  $h = \frac{3}{2}$  and  $j = \frac{1}{2}$ . Their OPE takes the form

$$
G ^ {I} (z _ {1}) \widetilde {G} ^ {J} (z _ {2}) \sim \frac {2 k \epsilon^ {I J}}{(z _ {1} - z _ {2}) ^ {3}} + \frac {2 J ^ {I J} (z _ {2})}{(z _ {1} - z _ {2}) ^ {2}} + \frac {\epsilon^ {I J} T (z _ {2}) + \partial J ^ {I J} (z _ {2})}{z _ {1} - z _ {2}}, \tag {2.6}
$$

where  $\epsilon^{+ - } = 1$ ,  $J^{\pm \pm} = J^{\pm}$  and  $J^{+ - } = J^{- + } = \frac{1}{2} J^{0}$ . The remaining OPE among fermionic generators are regular. The small  $\mathcal{N} = 4$  SCA possesses an  $SL(2)$  outer automorphism that rotates  $G$  and  $\widetilde{G}$  as a doublet and acts trivially on the bosonic generators. We denote by  $GL(1)_r$  the corresponding Cartan generator, normalized as  $r[G] = \frac{1}{2}$ ,  $r[\widetilde{G}] = -\frac{1}{2}$ . In order to avoid keeping track of  $\mathfrak{sl}(2)$  indices we use a standard index free notation by introducing the auxiliary variable  $y$  and write

$$
J (y) = J ^ {+} + J ^ {0} y + J ^ {-} y ^ {2}, \quad G (y) = G ^ {+} + G ^ {-} y, \quad \widetilde {G} (y) = \widetilde {G} ^ {+} + \widetilde {G} ^ {-} y, \tag {2.7}
$$

where we omitted the explicit  $z$ -dependence of the operators. Because of the auxiliary variable  $y$ , we often refer to the  $\mathfrak{sl}(2)$  R-symmetry of the small  $\mathcal{N} = 4$  SCA as  $\mathfrak{sl}(2)_y$ . In contrast, the conformal algebra  $\mathfrak{sl}(2)$  on the Riemann sphere with coordinate  $z$  will be referred to as  $\mathfrak{sl}(2)_z$ . More details are given in appendix A.

The  $\mathcal{N} = 2$  SCA can be defined as the sub-VOA of the small  $\mathcal{N} = 4$  SCA generated by

$$
\mathcal {J} = J ^ {0}, \quad \mathcal {G} = G ^ {-}, \quad \widetilde {\mathcal {G}} = \widetilde {G} ^ {+}, \quad \mathcal {T} = T. \tag {2.8}
$$

The remaining generators of the small  $\mathcal{N} = 4$  SCA, namely  $(J^{+},G^{+})$  and  $(J^{-},\widetilde{G}^{-})$ , are  $\mathcal{N} = 2$  chiral and anti-chiral multiplets respectively.

The global part of the small  $\mathcal{N} = 4$  and  $\mathcal{N} = 2$  SCA are  $\mathfrak{psl}(2|2)$  and  $\mathfrak{osp}(2|2)$  respectively. The representations of  $\mathfrak{psl}(2|2)$  and  $\mathfrak{osp}(2|2)$  relevant for this work are summarized in table 2. The bosonic subalgebra of  $\mathfrak{psl}(2|2)$  is  $\mathfrak{sl}(2)_z \oplus \mathfrak{sl}(2)_y$ , and representations are labelled by the conformal dimension  $h$  and the (half-integer) spin  $j$  of the superconformal primary. The bosonic subalgebra of  $\mathfrak{osp}(2|2)$  is  $\mathfrak{sl}(2)_z \oplus \mathfrak{gl}(1)$ , and representations are labelled by the conformal dimension  $h$  and the (half-integer)  $\mathfrak{gl}(1)$  charge  $m$  of the superconformal primary.

Finally we recall the definition of super-Virasoro primary. The operator  $\mathcal{O}$  is a super-Virasoro primary if the OPE of any super-Virasoro generator with  $\mathcal{O}$  does not contain any pole of order higher than one, with the obvious exception of the order-two pole in the  $T\mathcal{O}$  OPE, which encodes the conformal weight. The first order poles in the OPEs of the super-Virasoro generators with  $\mathcal{O}$  encode the action of the global part of the SCA. See appendix A for more details.

Notation: Occasionally we will use the following notation for poles in the OPE:

$$
A \left(z _ {1}\right) B \left(z _ {2}\right) = \sum_ {n} \frac {\left\{A B \right\} _ {n} \left(z _ {2}\right)}{\left(z _ {1} - z _ {2}\right) ^ {n}}. \tag {2.9}
$$

<table><tr><td>algebra</td><td>quantum numbers of s.c.p. X</td><td>shortening conditions</td><td>notation</td></tr><tr><td rowspan="2">psl(2|2)</td><td>h &gt; j</td><td>-</td><td>L(h,j)</td></tr><tr><td>h = j</td><td>G↑X = G↑X = 0</td><td>Sh</td></tr><tr><td rowspan="3">osp(2|2)</td><td>h &gt; |m|</td><td>-</td><td>X(h,m)</td></tr><tr><td>h = +m</td><td>G·X = 0</td><td>Ch</td></tr><tr><td>h = -m</td><td>G·X = 0</td><td>Ch</td></tr></table>

Table 2: Representations of  $\mathfrak{psl}(2|2)$  and  $\mathfrak{osp}(2|2)$  that appear in this work. We use  $X$  to the denote the superconformal primary, or s.c.p. for short. The symbol  $G^{\uparrow}X$  denotes the descendant of  $X$  with weight  $h + \frac{1}{2}$  and spin  $j + \frac{1}{2}$  obtained by acting once with the supercharge  $G$ . The notation  $\mathcal{G} \cdot X$  stands for the descendant of  $X$  with weight  $h + \frac{1}{2}$  and charge  $m + \frac{1}{2}$  obtained by acting once with the supercharge  $\mathcal{G}$ . Similar remarks apply to  $\widetilde{G}^{\uparrow}X$ ,  $\widetilde{\mathcal{G}} \cdot X$ . More details on our notation can be found in appendix A.

We will also use  $(AB)_n$  to represent the completion of  $\{AB\}_n$  to a quasiprimary, i.e. an  $\mathfrak{sl}(2)_z$  primary. The object  $(AB)_n$  is obtained from  $\{AB\}_n$  by adding  $z$ -derivatives of higher-order poles in the  $AB$  OPE, the explicit formula is given in (A.7). See e.g. [42] for details. Concerning the  $\mathfrak{sl}(2)_y$  structures, we will use the notation  $(AB)^j$  for the spin  $j$  projection of the relevant product of  $A$  with  $B$ , see appendix A for more details. Finally, given a quasiprimary  $X$  of  $\mathfrak{sl}(2)_y$  with spin  $j$ , we introduce the shorthand notation for its supersymmetric descendant

$$
G ^ {\downarrow} X = \left(G X\right) _ {1} ^ {j - \frac {1}{2}}. \tag {2.10}
$$

This is the quasiprimary in the order-one pole of the OPE of  $G$  with  $X$ , projected onto the component with spin  $j - \frac{1}{2}$ . Similar remarks apply to the operations  $G^{\uparrow}$  and  $\widetilde{G}^{\uparrow, \downarrow}$ .

# 2.2 Coxeter groups, complex reflection groups, rings of invariants

We will now describe the symplectic varieties (1.9) and (1.1) and the associated ring of functions in more detail. The case of Coxeter groups is described first, the generalization to complex reflection groups is presented in the end of this section. Let us start by fixing a basis of  $V_{\Gamma}$  to be  $z_{1},\ldots ,z_{r}$ , where  $r = \mathrm{rank}(\Gamma)$ . The action of the (finite) Coxeter group  $\Gamma$  on the real vector space  $V_{\Gamma}$  is generated by reflections. This can be taken to be the definition of Coxeter group. Their full list is given by the Weyl groups of finite dimensional semisimple Lie algebras  $ABCDEFG$  together with  $I_2(p)$ ,  $H_{3}$ ,  $H_{4}$ , see table 3 for more details. The expression for the central charge given in (1.4) can be rewritten as

$$
- \frac {c}{3} = \sum_ {\ell = 1} ^ {\operatorname {r a n k} (\Gamma)} (2 p _ {\ell} - 1) = | \Phi_ {\Gamma} | + \operatorname {r a n k} (\Gamma) \stackrel {\text {W e y l}} {=} \dim (\mathfrak {g} _ {\Gamma}), \tag {2.11}
$$

<table><tr><td>Γ</td><td>p1, ..., pr</td><td>(h, j) quantum numbers of long generators*</td><td>-c/3</td></tr><tr><td>An-1</td><td>2, 3, ..., n</td><td>-</td><td>n2 - 1</td></tr><tr><td>Bn</td><td>2, 4, ..., 2n</td><td>-</td><td>n(2n + 1)</td></tr><tr><td>Dn</td><td>2, 4, ..., 2(n - 1); n</td><td>(n+2/2, n-4/2),...</td><td>n(2n - 1)</td></tr><tr><td>E6</td><td>2, 5, 6, 8, 9, 12</td><td>(4, 0), (9/2, 3/2), (6, 3)</td><td>78</td></tr><tr><td>E7</td><td>2, 6, 8, 10, 12, 14, 18</td><td>(5, 1), (6, 3), (7, 3), (8, 5)</td><td>133</td></tr><tr><td>E8</td><td>2, 8, 12, 14, 18, 20, 24, 30</td><td>(6, 0), (7, 3), (9, 6), (9, 4), (9, 3), (10, 6), (10, 0)</td><td>248</td></tr><tr><td>F4</td><td>2, 6, 8, 12</td><td>(4, 0), (6, 3)</td><td>52</td></tr><tr><td>H3</td><td>2, 6, 10</td><td>-</td><td>33</td></tr><tr><td>H4</td><td>2, 12, 20, 30</td><td>(6, 0), (10, 6), (10, 0)</td><td>124</td></tr><tr><td>I2(p)</td><td>2, p</td><td>-</td><td>2(p + 1)</td></tr></table>

Table 3: The notation for Coxeter groups is such that when restricting to Weyl groups one has the identification  $A = \mathrm{Weyl}(\mathfrak{a})$ ,  $B = \mathrm{Weyl}(\mathfrak{b})$  and so on. Moreover, the Weyl group  $G_{2}$  appears as  $G_{2} = I_{2}(6)$ . Recall that  $\mathrm{Weyl}(\mathfrak{b}_n) = \mathrm{Weyl}(\mathfrak{c}_n)$ . The asterisk * above indicates that only long generators whose conformal weight  $h$  is smaller than the lightest relation are listed. These are the one that can be extracted unambiguously from the associated Hilbert series but there might be more generators. The ... in the  $D_{n}$  series are given explicitly for  $4\leq n\leq 9$  in equation (B.4). See appendix B.1 for more details.

where  $|\Phi_{\Gamma}|$  is the cardinality of the root system associated to  $\Gamma$

The main ring of interest in the following is the ring of invariant polynomials in two sets of variables  $z_i^{\pm}$

$$
\mathcal {R} _ {\Gamma} = \mathbb {C} \left[ z _ {1} ^ {+}, \dots , z _ {\mathrm {r}} ^ {+}, z _ {1} ^ {-}, \dots , z _ {\mathrm {r}} ^ {-} \right] ^ {\Gamma} = \mathbb {C} \left[ \mathcal {M} _ {\Gamma} \right], \tag {2.12}
$$

where the Coxeter group acts independently on  $z_{i}^{+}$  and  $z_{i}^{-}$ . This ring carries the action of  $GL(2) = SL(2) \times GL(1)$ , where  $z^{\pm}$  transform as a doublet under  $SL(2)$  and have the same  $GL(1)$  weight. This ring is thus graded by  $m$ , the eigenvalues of the Cartan of  $SL(2)$ , and  $h$  with  $m[z_{i}^{\pm}] = \pm \frac{1}{2}$ ,  $h[z_{i}^{\pm}] = \frac{1}{2}$ .

The ring  $\mathcal{R}_{\Gamma}$  has an alternative description in terms of generators and relations. Let us present the simplest example of  $\Gamma = \mathrm{Weyl}(\mathfrak{a}_1) = \mathbb{Z}_2$ , whose action is generated by  $\sigma \cdot z_1^{\pm} = -z_1^{\pm}$ . In this case the invariants are  $j^{\pm} = z_1^{\pm}z_1^{\pm}$  and  $j^0 = 2z_1^+ z_1^-$  with the obvious relation  $j^{+}j^{-} = \frac{1}{4} j^{0}j^{0}$ . For groups  $\Gamma$  of higher rank giving an explicit description of this type is more involved but in principle straightforward. We do so in the low rank examples of  $I_{2}(p)$ ,  $A_{3}$ ,  $B_{3}$ ,  $H_{3}$  and comment on the general case in appendix B.1. One can develop an opinion about the set of generators and relations by considering the Hilbert series of  $\mathcal{R}_{\Gamma}$  that we will now review.

Hilbert series. The (refined) Hilbert series of  $\mathcal{R}_{\Gamma}$  is defined as

$$
\mathrm {H S} _ {\Gamma} (\tau , x) = \operatorname {T r} _ {\mathscr {R} _ {\Gamma}} \left(\tau^ {2 h} x ^ {2 m}\right). \tag {2.13}
$$

In the case of ring of invariants as (2.12) the Hilbert series can be computed by averaging over  $\Gamma$  the Hilbert series of the freely generated ring  $\mathbb{C}[z^{+},z^{-}]$ . This is the content of the Molien formula

$$
\mathsf {H S} _ {\Gamma} (\tau , x) = \mathsf {M o l i e n} _ {\Gamma} (\tau , x) := \frac {1}{| \Gamma |} \sum_ {g \in \Gamma} \frac {1}{\det  _ {\mathbb {C} ^ {2} \otimes V _ {\Gamma}} (1 - h \otimes g)}, \quad h = \tau \left( \begin{array}{c c} x & 0 \\ 0 & x ^ {- 1} \end{array} \right), \tag {2.14}
$$

where  $|\Gamma|$  is the order of  $\Gamma$ . The simplest example is given by

$$
\operatorname {M o l i e n} _ {\mathbb {Z} _ {2}} (\tau , x) = \frac {1 - \tau^ {4}}{(1 - \tau^ {2} x ^ {- 2}) (1 - \tau^ {2}) (1 - \tau^ {2} x ^ {+ 2})}. \tag {2.15}
$$

In this formula the denominator can be interpreted as the contribution of the generators  $j^{+}, j^{-}, j^{0}$  defined above and the subtraction of the terms  $\tau^4$  in the numerator corresponds to the relation  $j^{+}j^{-} = \frac{1}{4} j^{0}j^{0}$ . While it is not possible in general to extract the set of generators and relations, together with their  $h, m$  quantum numbers, from the series (2.13) alone, one can show that certain generators and relations must be present, see e.g. [43]. This can be done efficiently by using the so-called plethystic logarithm. The expressions are collected in appendix B.1 and the resulting set of generators is given in table 3. For convenience of the reader we collect the expression for Molien generating functions for all Coxeter groups in an ancillary Mathematica file. It should be noticed that the generators of  $\mathcal{R}_{\Gamma}$  can be divided in two groups: one consisting of elements with quantum numbers  $h = j$ , which will be referred to as short generators, and the other with quantum numbers  $h > j$ , which will be referred to as long generators. We do not give a complete list of long generators for general  $\Gamma$  but list generators whose existence can be shown unambiguously using the Molien series (2.14). The result is collected in table 3.

Symplectic structure. The ring  $\mathcal{R}_{\Gamma}$  possesses the important property of being a Poisson algebra with Poisson bracket

$$
\left\{z _ {i} ^ {I}, z _ {j} ^ {J} \right\} _ {\mathrm {P B}} = \eta_ {i j} \epsilon^ {I J}, \quad i, j = 1, 2, \dots , r, \quad I, J = \pm , \tag {2.16}
$$

where  $\epsilon^{IJ}$  is antisymmetric and normalized as  $\epsilon^{+ - } = 1$  while  $\eta_{ij}$  is symmetric and nondegenerate. This implies that (1.9) is a symplectic variety. It is straightforward to compute the Poisson bracket of the generators of  $\mathcal{R}_{\Gamma}$  using (2.16). In the simple example of  $\Gamma = \mathbb{Z}_2$  this gives the Lie algebra of  $SL(2)$ . For higher rank this procedure will produce an extension of it. A priori, there is no canonical choice for the set of generators, but there is a distinguished one originating from the associated VOA. We remark that when long generators are present they can be generated by Poisson brackets of short generators, as the example of  $D_4$  given in section 4.4 illustrates.

Complex reflection groups. The case of complex reflection groups is very similar so we will be brief. The first important difference is that the  $GL(1)\times SL(2)$  symmetry of the Coxeter case is reduced to  $GL(1)\times GL(1)$ . The Molien formula is a slight generalization of (2.14) to

$$
\operatorname {M o l i e n} _ {\mathsf {G}} (\tau , x) = \frac {1}{| \mathsf {G} |} \sum_ {g \in \mathsf {G}} \frac {1}{\det  _ {V _ {\mathsf {G}}} (1 - \tau x g) \det  _ {V _ {\mathsf {G}} ^ {*}} (1 - \tau x ^ {- 1} g)}. \tag {2.17}
$$

Finally the generalization of Poisson brackets (2.16) uses the canonical pairing between  $V_{\mathsf{G}}$  and its dual  $V_{\mathsf{G}}^{*}$ .

Remark: One might refer to (2.12) as the Higgs branch chiral rings. It contains a subring defined as the graded component with  $h = m$ . The latter can be referred to as Coulomb branch chiral ring and is a freely generated polynomial ring. This is the content of a famous theorem of Chevalley, Shephard and Todd states that the ring of invariants  $\mathbb{C}[V]^{\mathsf{G}}$  is a freely generated polynomial ring if and only if  $\mathsf{G}$  acts as a complex reflection group on  $V$ . These finite groups have been classified by Shephard and Todd, see e.g. [44].

# 3 Free-field realizations

This section is devoted to our proposal for a free-field realization of the VOAs associated to Coxeter and complex reflection groups.

# 3.1 Realization of the  $\mathcal{N} = 2$  SCA in terms of  $\beta \gamma bc$  systems

To begin with, we present a free-field realization of the  $\mathcal{N} = 2$  SCA. According to our proposal, the relevant free fields consist of  $\mathsf{r} = \mathrm{rank}(\mathsf{G})$  copies of a free  $\beta \gamma bc$  system. The relevant nontrivial OPEs are simply

$$
\beta_ {\ell_ {1}} (z _ {1}) \gamma_ {\ell_ {2}} (z _ {2}) = - \frac {\delta_ {\ell_ {1} \ell_ {2}}}{z _ {1 2}} + \operatorname {r e g .}, \quad b _ {\ell_ {1}} (z _ {1}) c _ {\ell_ {2}} (z _ {2}) = \frac {\delta_ {\ell_ {1} \ell_ {2}}}{z _ {1 2}} + \operatorname {r e g .}, \tag {3.1}
$$

where  $\ell_1, \ell_2$  run from 1 to  $r$  and  $\delta_{\ell_1 \ell_2}$  denotes the Kronecker delta. The generators of the  $\mathcal{N} = 2$  SCA take a simple expression in terms of the free fields, consisting of a direct sum of terms, one for each copy of the  $\beta \gamma bc$  system. More precisely,

$$
\mathcal {J} = \sum_ {\ell = 1} ^ {r} \left[ p _ {\ell} \beta_ {\ell} \gamma_ {\ell} + (p _ {\ell} - 1) b _ {\ell} c _ {\ell} \right],
$$

$$
\mathcal {G} = \sum_ {\ell = 1} ^ {r} b _ {\ell} \gamma_ {\ell}, \qquad \widetilde {\mathcal {G}} = \sum_ {\ell = 1} ^ {r} \left[ p _ {\ell} \beta_ {\ell} \partial c _ {\ell} + (p _ {\ell} - 1) \partial \beta_ {\ell} c _ {\ell} \right],
$$

$$
\mathcal {T} = \sum_ {\ell = 1} ^ {r} \left[ - \frac {1}{2} p _ {\ell} \beta_ {\ell} \partial \gamma_ {\ell} + \left(1 - \frac {1}{2} p _ {\ell}\right) \partial \beta_ {\ell} \gamma_ {\ell} - \frac {1}{2} (p _ {\ell} + 1) b _ {\ell} \partial c _ {\ell} + \frac {1}{2} (1 - p _ {\ell}) \partial b _ {\ell} c _ {\ell} \right]. \quad (3. 2)
$$

The quantities  $p_{\ell}$  are the degrees of the invariants of the Coxeter group. The central charge of the  $\mathcal{N} = 2$  SCA is given by

$$
c = - 3 \sum_ {\ell = 1} ^ {r} \left(2 p _ {\ell} - 1\right), \tag {3.3}
$$

compare to (1.4). The conformal weights  $h$ , the  $\mathfrak{gl}(1)$  charges  $m$  and the  $\mathfrak{gl}(1)_r$  charge  $r$  defined below (2.6) of the free fields are summarized in table (3.4). The charge  $m$  is normalized in such a way that  $h = m$  for chiral primary operators. Notice that, although the combined  $\beta \gamma bc$  system contains states of negative conformal dimension, all states have a non-negative "twist"  $h - m$ . Furthermore, the space of states with given twist and charge is finite-dimensional.

<table><tr><td></td><td>h</td><td>m</td><td>h-m</td><td>h+m</td><td>r</td></tr><tr><td>βl</td><td>1/2p l</td><td>1/2p l</td><td>0</td><td>p l</td><td>0</td></tr><tr><td>bl</td><td>1/2(p l+1)</td><td>1/2(p l-1)</td><td>1</td><td>p l</td><td>+1/2</td></tr><tr><td>cl</td><td>-1/2(p l-1)</td><td>-1/2(p l-1)</td><td>0</td><td>1-p l</td><td>-1/2</td></tr><tr><td>γl</td><td>1-1/2p l</td><td>-1/2p l</td><td>1</td><td>1-p l</td><td>0</td></tr><tr><td>δ</td><td>1</td><td>0</td><td>1</td><td>1</td><td>0</td></tr></table>

(3.4)

A comment on normal-ordered products in the  $\beta \gamma bc$  system is in order. In all equations in (3.2), the juxtaposition of free fields is understood as their  $\{\cdot \cdot \}_{0}$  normal-ordered product, see (2.9). More generally, let  $X_{1},\ldots ,X_{n}$  stand for any of the free fields  $\beta ,\gamma ,b,c,$  or any  $z$  derivative thereof. A natural object is the nested normal-ordered product

$$
: X _ {1} X _ {2} \dots X _ {n}: = \left\{X _ {1} \left\{X _ {2} \left\{\dots \left\{X _ {n - 1} X _ {n} \right\} _ {0} \dots \right\} _ {0} \right\} _ {0} \right.. \tag {3.5}
$$

Since we are considering a free theory, and since the  $X_{i}$  are (derivatives of) free fields, in the normal-ordered product:  $X_{1}X_{2}\ldots X_{n}$ : we are free to permute the factors  $X_{i}$ , up to Grassmann signs, exactly as we would do in a supercommutative algebra. For instance

$$
: \beta \partial \beta \gamma \gamma \gamma : = \left\{\beta \left\{\partial \beta \left\{\gamma \left\{\gamma \gamma \right\} _ {0} \right\} _ {0} \right\} _ {0} \right\} _ {0} = \left\{\gamma \left\{\partial \beta \left\{\gamma \left\{\gamma \beta \right\} _ {0} \right\} _ {0} \right\} _ {0} \right\} _ {0} =: \gamma \partial \beta \gamma \gamma \beta :. \tag {3.6}
$$

This allows to write compactly:  $\beta \partial \beta \gamma^3$ : without ambiguities. Notice that such manipulations are not allowed in a generic VOA. For the sake of brevity, we omit the colons from the normal-ordered products in the rest of this work.

# 3.2 Realization of the remaining generators in the  $\mathcal{N} = 4$  case

Let us first discuss the case of the  $\mathcal{N} = 4$  VOA  $\mathcal{W}_{\Gamma}$  associated to a Coxeter group  $\Gamma$ . As outlined in section 1.1, the set of strong generators of this VOA includes elements transforming in short representations  $\mathfrak{S}_{\frac{p_{\ell}}{2}}$  of the global conformal algebra  $\mathfrak{psl}(2|2)$ , see table 2. Among these, there is a distinguished invariant of degree two which corresponds to the generators of the small  $\mathcal{N} = 4$  super-Virasoro VOA. This will be labelled by  $\ell = 1$  so that  $p_1 = 2$ . The remaining short

generators are denoted as  $^{16}W_{\ell}$ . A notable feature of the proposed free-field realization is that the highest weight states of the short generators are identified with  $\beta$ s, see (3.1), more precisely

$$
J ^ {+} = \beta_ {1}, \quad G ^ {+} := \left\{G ^ {-} J ^ {+} \right\} _ {1} \quad = b _ {1}, \tag {3.7a}
$$

$$
W _ {\ell} ^ {\mathrm {h . w .}} = \beta_ {\ell}, \quad G _ {W _ {\ell}} ^ {\mathrm {h . w .}} := \left\{G ^ {-} W _ {\ell} ^ {\mathrm {h . w .}} \right\} _ {1} = b _ {\ell}, \tag {3.7b}
$$

$\ell = 2, \ldots, \mathrm{rank}(\Gamma)$ . Two remarks are in order. First notice that the  $(h, m)$  weights assignment of these object is consistent by construction, see table (3.4). The  $\mathcal{N} = 2$  subalgebra, embedded as specified by (2.8), is realized as in (3.2). The second remark is that each pair  $(\beta_{\ell}, b_{\ell})$  forms an  $\mathcal{N} = 2$  chiral multiplet as (3.7) indicates.

The next generator that needs to be constructed is  $J^{-}$ , the  $\widehat{\mathfrak{sl}(2)}$  affine Kac-Moody current with weight  $m = -1$ , see section 2.1. Once this operator is constructed one can build the whole  $\mathcal{N} = 4$  super-conformal multiplets to which (3.7) belong by taking appropriate poles in the OPE with  $J^{-}$  and its  $\mathcal{N} = 2$  partner  $\widetilde{G}^{-}$ . Next, one takes the OPE of the generators obtained in this way. If  $\mathcal{W}_{\Gamma}$  does not have long generators, see table 3, the latter OPE needs to close on the generators that have already been constructed. If long generators are present, they will be defined by the failure of these OPE to close on the short generators. An example of this mechanism, which is already at play at the classical level when closing the Poisson brackets of the short generators of the ring  $\mathcal{R}_{\Gamma}$ , is given for the example of  $D_{4}$  in section 4.4.

In order to construct  $J^{-}$ , we build an Ansatz and impose necessary conditions that this operator must obey. The construction goes as follows:

1. Construct the most general Ansatz for  $J^{-}$  in terms of the  $\beta \gamma bc$  free fields. There is always a finite number of terms in the Ansatz. This is easy to verify by recalling that  $J^{-}$  has weights  $h = 1$ ,  $h - m = 2$  and by staring at weight assignments of the constituent free fields given in table 3.4.  
2. Impose that the small  $\mathcal{N} = 4$  algebra closes, this is equivalent to:

a. Linear constraints:  $J^{-}$  is an anti-chiral  $\mathcal{N} = 2$  super-Virasoro primary of weight  $h = 1$  and its OPE with  $J^{+}$  closes on  $J^{0}$ .  
b. Non-linear constraints:  $J^{-}$  has a regular OPE with itself.

3. Impose that the short generators  $W_{\ell}^{\mathrm{h.w.}}$  given in (3.7) must be super-Virasoro primary, this implies:

$$
\left\{J ^ {-} \beta_ {\ell} \right\} _ {n \geq 2} = 0, \quad \ell = 2, \dots , \operatorname {r a n k} (\Gamma). \tag {3.8}
$$

4. Impose that the short generators  $W_{\ell}$  have non-zero norms.  
5. Impose that the VOA closes on the strong generators.

It is convenient to write the Ansatz for  $J^{-}$  in the form

$$
J ^ {-} = J _ {\mathrm {m i n}} ^ {-} + J _ {\mathrm {n o r m s}} ^ {-}, \tag {3.9}
$$

where

$$
J _ {\min } ^ {-} := k \partial \gamma_ {1} + \beta_ {1} (\gamma_ {1}) ^ {2} + \gamma_ {1} b _ {1} c _ {1} + \gamma_ {1} \widehat {\mathcal {I}} - c _ {1} \widehat {\mathcal {G}}. \tag {3.10}
$$

The hat on  $\mathcal{J}$  and  $\mathcal{G}$  signals the omission of all terms built with  $\beta_{1},\gamma_{1},b_{1},c_{1}$  in (3.2). For all  $\Gamma$  different from  $A_{1}$  the factor  $J_{\mathrm{norms}}^{-}$  has to be non-zero. As the name suggests, it must be there in order for point 4. above to be satisfied. To illustrate this point, let us explain what happens if we set  $J_{\mathrm{norms}}^{-} = 0$  in (3.9). It is easy to verify that in this case the  $\mathcal{N} = 4$  superVirasoro subalgebra is correctly reproduced. Next, let us construct  $W_{\ell}^{\mathrm{min}}(z,y)$  by summing up the  $\mathfrak{sl}(2)_y$  descendants of (3.7b) defined by the action of  $J_{\mathrm{min}}^{-}$ . A little computation shows that

$$
W _ {\ell} ^ {\min } (y) = \left(1 + y \gamma_ {1}\right) ^ {p _ {\ell}} \beta_ {\ell} - \left(1 + y \gamma_ {1}\right) ^ {p _ {\ell} - 1} c _ {1} b _ {\ell}, \quad \ell = 2, \dots , \operatorname {r a n k} (\Gamma). \tag {3.11}
$$

These operators have obviously regular OPE among themselves, in particular they have zero norm. This explains the necessity of adding  $J_{\mathrm{norms}}^{-}$ . By using an Ansatz of the form (3.9) in the steps above one quickly verifies that  $J_{\mathrm{norms}}^{-}$  does not include  $\gamma_{1}$  and  $c_{1}$ .

Before describing various examples of this construction in the next section, let us make a few remarks:

Remark 1: The generators of the  $\mathcal{N} = 2$  subalgebra in (3.2) are invariant under the transformation<sup>17</sup>

$$
\left(\beta_ {\ell}, \gamma_ {\ell}, b _ {\ell}, c _ {\ell}\right) \mapsto \left(\lambda_ {\ell} \beta_ {\ell}, \lambda_ {\ell} ^ {- 1} \gamma_ {\ell}, \lambda_ {\ell} b _ {\ell}, \lambda_ {\ell} ^ {- 1} c _ {\ell}\right), \quad \lambda_ {\ell} \in G L (1), \tag {3.12}
$$

for  $\ell = 2,\ldots ,\mathrm{rank}(\Gamma)$ . We claim that  $J^{-}$  is uniquely determined by the steps above up to this ambiguity.

Remark 2: For low ranks  $\mathrm{rank}(\Gamma) = 1,2,3$  the steps 1.-4. are sufficient to determine  $J^{-}$  up to the action (3.12) and condition 5. holds automatically. For higher ranks condition 5. needs to be used as well. As the example of  $D_{4}$  illustrates, see section 4.4, there is a subset of these conditions that is easy to implement and is sufficient to fully determine  $J^{-}$ .

Remark 3. It is natural to ask what happens if one applies the procedure presented above to a set of weights  $p_1 = 2, p_2, \ldots, p_r$  that do not correspond to a Coxeter group. We observed experimentally, by looking at rank 3 examples with  $2 \leq p_2 \leq p_3 \leq 10$ , that the norms mentioned above are non-zero if and only if the weights are the one associated to a Coxeter group, which in rank 3 are  $A_3, B_3, H_3$ , see table 3.

# 3.3 Realization of the remaining generators in the  $\mathcal{N} = 2$  case

The realization of the  $\mathcal{N} = 2$  VOA associated to a complex reflection group  $\mathsf{G}$  is qualitatively similar. The generators of the  $\mathcal{N} = 2$  SCA algebra are given in (3.2) and the additional chiral

generators have the form

$$
W _ {\ell} = \beta_ {\ell}, \quad \mathcal {G} _ {W _ {\ell}} := \left\{\mathcal {G} W _ {\ell} \right\} _ {1} = b _ {\ell}, \quad \ell = 1, \dots , \operatorname {r a n k} (\mathrm {G}). \tag {3.13}
$$

As outlined in section 1.1 the complete set of generators can be found from the set of generators of the ring  $\mathcal{R}_{\mathsf{G}}$ . This include in particular anti-chiral operators with conformal weights  $\frac{p_{\ell}}{2}$  which will be denoted as  $\overline{W}_{\ell}$ . Concerning the construction of the remaining generators we propose the following strategy:

1. Make an Ansatz for the remaining generators in terms of free fields and impose that they are  $\mathcal{N} = 2$  super-Virasoro primary with the correct  $(h,m)$  weights. As in the  $\mathcal{N} = 4$  case there is a finite number of terms in the Ansatz.  
2. Impose that the VOA closes on the strong generators. This gives both linear and nonlinear conditions on the coefficients of the Ansatz. It is practically convenient to first solve the linear constraints.

We will present all rank 1 examples as well as a rank 2 example of this procedure in section 5.

# 3.4 Free-field realization and classical Poisson structure

We will now show how the Poisson algebra (2.12), (2.16) can be obtained starting from the free-field realization of  $\mathcal{W}_{\Gamma}$ . Let  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)}$  denote the free  $\beta \gamma bc$  system associated to  $\Gamma$ , consisting of  $r$  copies of a single  $\beta \gamma bc$  system. Let  $\mathbb{M}_{\beta \gamma}^{(\Gamma)\mathrm{cl}}$  denote the classical Poisson algebra comprised by all polynomials in the variables  $\beta_{\ell}$ ,  $\gamma_{\ell}$ ,  $\ell = 1, \ldots, r$ , with Poisson bracket

$$
\{f, g \} _ {\mathrm {P B}} = \sum_ {\ell = 1} ^ {r} \left[ \partial_ {\gamma_ {\ell}} f \partial_ {\beta_ {\ell}} g - \partial_ {\beta_ {\ell}} f \partial_ {\gamma_ {\ell}} g \right]. \tag {3.14}
$$

We can define a linear map

$$
\mathcal {P}: \mathbb {M} _ {\beta \gamma b c} ^ {(\Gamma)} \rightarrow \mathbb {M} _ {\beta \gamma} ^ {(\Gamma) \mathrm {c l}}, \tag {3.15}
$$

according to the following prescription. Any element of  $\mathbb{M}_{\beta \gamma bc}^{\mathrm{(F)}}$  can be cast as a linear combination of nested normal ordered products of derivatives of free fields. The image under  $\mathcal{P}$  of such an object is simply obtained by dropping all terms with derivatives and/or fermionic free fields  $b_{\ell}, c_{\ell}$ , and by replacing all normal-ordered products with regular products in the algebra of polynomials in the variables  $\beta_{\ell}, \gamma_{\ell}$ . One may verify that  $\mathcal{P}$  is well-defined. Furthermore, the map  $\mathcal{P}$  satisfies

$$
\mathcal {P} \left(\left\{X _ {1} X _ {2} \right\} _ {0}\right) = \mathcal {P} \left(X _ {1}\right) \mathcal {P} \left(X _ {2}\right),
$$

$$
\mathcal {P} \left(\left\{X _ {1} X _ {2} \right\} _ {1}\right) = \left\{\mathcal {P} \left(X _ {1}\right), \mathcal {P} \left(X _ {2}\right) \right\} _ {\mathrm {P B}}. \tag {3.16}
$$

On the right hand side of the first relation, the product is the commutative and associative product in the algebra of polynomials in  $\beta_{\ell}$ ,  $\gamma_{\ell}$ . Since the VOA  $\mathcal{W}_{\Gamma}$  associated to a given Coxeter group  $\Gamma$  is realized as a subalgebra of  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)}$ , we can apply the map  $\mathcal{P}$  to any element

of  $\mathcal{W}_{\Gamma}$ , thus defining the Poisson algebra  $\mathcal{P}(W_{\Gamma})$ . We claim the following isomorphism of Poisson algebras,

$$
\mathcal {P} \left(\mathcal {W} _ {\Gamma}\right) \simeq \mathcal {R} _ {\Gamma}, \tag {3.17}
$$

where  $\mathcal{R}_{\Gamma}$  is defined in (2.12) and its symplectic structure is given in (2.16). The explicit form of the isomorphism (3.17) will be given in some examples in section 4. The case of complex reflection groups is identical.

In close analogy to the above discussion about the map  $\mathcal{P}$ , we can also define the map

$$
\mathcal {P} ^ {\prime}: \mathbb {M} _ {\beta \gamma b c} ^ {(\Gamma)} \rightarrow \mathbb {M} _ {\beta \gamma b c} ^ {(\Gamma) \mathrm {c l}}. \tag {3.18}
$$

In the above expression,  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)\mathrm{cl}}$  denotes the Poisson superalgebra of functions of the classical Grassmann even variables  $\beta_{\ell}, \gamma_{\ell}$  and Grassmann odd variables  $b_{\ell}, c_{\ell}$ . The usual properties of the Poisson bracket in a Poisson algebra hold, up to the obvious modifications due to Grassmann signs. In particular, the Poisson bracket on  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)\mathrm{cl}}$  is entirely specified by

$$
\left\{\beta_ {\ell_ {1}}, \gamma_ {\ell_ {2}} \right\} _ {\mathrm {P B}} = - \delta_ {\ell_ {1} \ell_ {2}}, \quad \left\{b _ {\ell_ {1}}, c _ {\ell_ {2}} \right\} _ {\mathrm {P B}} = \delta_ {\ell_ {1} \ell_ {2}}. \tag {3.19}
$$

The map  $\mathcal{P}'$  is defined as follows. Given any object in the VOA  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)}$ , presented as a polynomial of normal orders of derivatives of the free fields, its image under  $\mathcal{P}'$  is obtained by setting to zero all  $z$ -derivatives and by replacing the normal ordered product of the VOA with the associative, supercommutative product in the Poisson superalgebra  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)\mathrm{cl}}$ . The map  $\mathcal{P}'$  will be useful in section 6.4 in relation to the discussion of the Hall-Littlewood ring.

# 4 Examples of  $\mathcal{N} = 4$  VOA  $\mathcal{W}_{\Gamma}$

In this section we present the proposed free-field construction of  $\mathcal{W}_{\Gamma}$  in some examples. We start by reviewing the rank one case  $\Gamma = A_{1}$  following [26]. Next we present all rank two and three cases, namely  $I_2(p),A_3,B_3,H_3$  , and some aspects of the interesting example of  $D_{4}$  . All algebraic manipulations were performed on a laptop using the Mathematica package introduced in [42]. The analysis of higher rank Coxeter groups along the lines of this paper will require the use of more computational power and/or packages that deal more efficiently wit  $\beta \gamma bc$  free fields.

# 4.1 Rank 1:  $\Gamma = A_{1}$

The VOA associated to the Coxeter group  $\Gamma = A_{1}$  is simply the small  $\mathcal{N} = 4$  SCA with central charge  $c = -9$ . Our proposed free-field realization reduces in this case to the free-field realization studied in [26]. All generators of the small  $\mathcal{N} = 4$  SCA are expressed in terms of a single  $\beta \gamma bc$  system. In this simple case, the object  $J_{\mathrm{min}}^{-}$  introduced in (3.10) is actually sufficient to obtain the desired free-field realization of all the generators of the small  $\mathcal{N} = 4$

SCA. For the convenience of the reader, we summarize here all the relevant formulae,

$$
J ^ {+} = \beta ,
$$

$$
J ^ {0} = b c + 2 \beta \gamma ,
$$

$$
J ^ {-} = \beta \gamma \gamma + \gamma b c - \frac {3}{2} \partial \gamma ,
$$

$$
G ^ {+} = b,
$$

$$
G ^ {-} = b \gamma ,
$$

$$
\widetilde {G} ^ {+} = c \partial \beta + 2 \partial c \beta ,
$$

$$
\tilde {G} ^ {-} = - b \partial c c + 2 \beta \gamma \partial c + \partial \beta \gamma c - \frac {3}{2} \partial^ {2} c,
$$

$$
T = - \frac {3}{2} b \partial c - \beta \partial \gamma - \frac {1}{2} \partial b c. \tag {4.1}
$$

The ring  $\mathcal{R}_{A_1}$  defined in (2.12) is generated by the  $SL(2)$  triplet

$$
j (y) = j ^ {+} + y j ^ {0} + y ^ {2} j ^ {-} = z _ {1} (y) z _ {1} (y), \quad z _ {1} (y) = z _ {1} ^ {+} + y z _ {1} ^ {-}, \tag {4.2}
$$

subject to the relation

$$
j ^ {+} j ^ {-} - \frac {1}{4} \left(j ^ {0}\right) ^ {2} = 0. \tag {4.3}
$$

The corresponding composite operator in the small  $\mathcal{N} = 4$  SCA is the  $\mathfrak{psl}(2|2)$  primary operator

$$
\mathfrak {L} _ {2, 0} = (J J) _ {0} ^ {0} + \frac {1}{3} T = \frac {2}{3} \left\{J ^ {+} J ^ {-} \right\} _ {0} - \frac {1}{6} \left\{J ^ {0} J ^ {0} \right\} _ {0} + \frac {1}{3} \partial_ {z} J ^ {0} + \frac {1}{3} T, \tag {4.4}
$$

where  $(JJ)_0^0$  denotes the quasiprimary completion of the normal ordered product of two  $J$ 's, projected onto the spin-0 component, see appendix A for more details on the notation. It is straightforward to check that the composite operator  $\mathfrak{L}_{2,0}$  is identically zero in the free-field realization (4.1). This comes as no surprise, since it has been proven in [26] that (4.1) is a free-field realization of the simple quotient of the small  $\mathcal{N} = 4$  SCA at  $c = -9$ , implying that all super-Virasoro descendants of the identity operator that are null for  $c = -9$  are automatically zero in the free-field realization.

It is worth recalling that in [26] it is also proven that the simple quotient of the small  $\mathcal{N} = 4$  SCA at  $c = -9$  can be characterized as the kernel of a suitable screening operator  $\mathbb{S}$  acting on the free  $\beta \gamma bc$  system. In order to write down the screening operator, we first have to express  $\beta$  and  $\gamma$  in terms of chiral bosons  $\chi$ ,  $\phi$ ,

$$
\beta = e ^ {\chi + \phi}, \quad \gamma = \partial \chi e ^ {- \chi - \phi}. \tag {4.5}
$$

The chiral bosons non-trivial OPEs are

$$
\chi \left(z _ {1}\right) \chi \left(z _ {2}\right) = + \log z _ {1 2} + \operatorname {r e g}. , \quad \phi \left(z _ {1}\right) \phi \left(z _ {2}\right) = - \log z _ {1 2} + \operatorname {r e g}. . \tag {4.6}
$$

Using this notation, the screening operator and the screening current read

$$
\mathbb {S} = \int d z \mathrm {J} (z), \quad \mathrm {J} = b e ^ {- \frac {1}{2} (\chi + \phi)}. \tag {4.7}
$$

The screening current  $J$  has conformal dimension 1 and  $J^0$ -eigenvalue 0. It acts on the  $\beta \gamma bc$  system via the order-one pole in the OPE,

$$
X \mapsto \mathbb {S} \cdot X = \{\mathrm {J} X \} _ {1}, \tag {4.8}
$$

where  $X$  is any operator in the  $\beta \gamma bc$  system. It is also worth pointing out that  $J$  can be written as a supersymmetry descendant of a chiral operator  $K$  with dimension  $1/2$  and  $J^0$ -eigenvalue  $1/2$ ,

$$
J = \left\{G ^ {-} K \right\} _ {1}, \quad K = e ^ {\frac {1}{2} (\chi + \phi)}. \tag {4.9}
$$

# 4.2 Rank 2:  $\Gamma = I_2(p)$

This section is devoted to a detailed description of the free-field realization of the  $\mathcal{N} = 4$  VOA associated to the Coxeter group  $I_2(p)$ ,  $p \geq 3$ .

# 4.2.1 The Coxeter groups  $I_2(p)$  and associated rings

The group  $I_2(p)$  is the symmetry group of the regular  $p$ -gon on the plane. As a Coxeter group, it is the subgroup of  $O(2,\mathbb{R})$  generated by two reflections with respect to two lines forming an angle  $\pi / p$ . Equivalently, we may regard it as generated by a reflection  $\sigma$  and a rotation  $\rho$  by an angle  $2\pi / p$ ,

$$
\sigma = \left( \begin{array}{c c} 1 & 0 \\ 0 & - 1 \end{array} \right), \quad \rho = \left( \begin{array}{c c} \cos \frac {2 \pi}{p} & - \sin \frac {2 \pi}{p} \\ \sin \frac {2 \pi}{p} & \cos \frac {2 \pi}{p} \end{array} \right). \tag {4.10}
$$

The action on  $I_2(p)$  on  $\mathbb{R}^2$  is extended naturally to  $\mathbb{C}^2$ . It is then convenient to perform a change of basis, and introduce coordinates  $z^{1,2}$  in  $\mathbb{C}^2$  such that

$$
\sigma : \left( \begin{array}{c} z _ {1} \\ z _ {2} \end{array} \right) \mapsto \left( \begin{array}{c} z _ {2} \\ z _ {1} \end{array} \right), \quad \rho : \left( \begin{array}{c} z _ {1} \\ z _ {2} \end{array} \right) \mapsto \left( \begin{array}{c} e ^ {+ 2 \pi i / p} z _ {1} \\ e ^ {- 2 \pi i / p} z _ {2} \end{array} \right). \tag {4.11}
$$

The original space  $\mathbb{R}^2\subset \mathbb{C}^2$  is recovered via the reality condition  $z_{2} = (z_{1})^{*}$ . In terms of  $z_{1,2}$  the invariants of  $I_{2}(p)$  take the simple form

$$
\mathcal {I} _ {2} = z _ {1} z _ {2}, \quad \mathcal {I} _ {p} = z _ {1} ^ {p} + z _ {2} ^ {p}. \tag {4.12}
$$

The Coxeter group  $I_2(p)$  has the property of being crystallographic only for  $p = 3,4,6$ . This is equivalent to the fact that the plane  $\mathbb{R}^2$  admits a tessellations by regular  $p$ -gons only for  $p = 3,4,6$ , i.e. triangles, squares and hexagons. These are also the value for which it coincides with a Weyl group, namely

$$
I _ {2} (3) \cong \operatorname {W e y l} (\mathfrak {a} _ {2}), \quad I _ {2} (4) \cong \operatorname {W e y l} (\mathfrak {b} _ {2}) \cong \operatorname {W e y l} (\mathfrak {c} _ {2}), \quad I _ {2} (6) \cong \operatorname {W e y l} (\mathfrak {g} _ {2}). \tag {4.13}
$$

Let us describe the ring  $\mathcal{R}_{I_2(p)}$  defined in (2.12). The set of generators in this case is obtained by promoting the generators (4.12) to  $SL(2)$  multiplets

$$
j (y) = z _ {1} (y) z _ {2} (y), \quad w (y) = z _ {1} (y) ^ {p} + z _ {2} (y) ^ {p}, \tag {4.14}
$$

where  $z_{i}(y) = z_{i}^{+} + yz_{i}^{-}$ . While (4.12) are algebraically independent the generators  $j(y)$  and  $w(y)$  satisfy the following relations

$$
\left. (j w) \right| _ {\frac {p}{2} - 1} = 0, \qquad \left. (w w) \right| _ {p - 2 m} + c _ {p, m} \left. (j ^ {2} \big | _ {0}\right) ^ {m} j ^ {p - 2 m} \big | _ {p - 2 m} = 0 \qquad \qquad (4. 1 5)
$$

$m = 1,\ldots ,\left\lfloor \frac{p}{2}\right\rfloor$  and the notation  $\big|_{*}$  denotes the projection onto the  $SL(2)$  spin  $*$  component. Using the same normalization of  $SL(2)$  projections as in appendix A, the coefficients  $c_{p,m}$  take the form

$$
c _ {p, m} = - \frac {2 \left(\frac {3}{2}\right) ^ {m} (\Gamma (p + 1)) ^ {2} \Gamma (m - p - \frac {1}{2})}{\Gamma (2 m + 1) \Gamma (1 - 2 m + p) \Gamma (1 - m + p) \Gamma (2 m - p - \frac {1}{2})}. \tag {4.16}
$$

We have obtained an alternative description of the ring  $\mathcal{R}_{I_2(p)}$  given in (2.12) as the ring generates by the  $y$ -components of (4.14) subject to the relations (4.15). As we will show shortly, these relations are the image of null states in the chiral algebra  $\mathcal{W}_{I_2(p)}$ .

Finally, let us comment on the Poisson algebra structure of  $\mathcal{R}_{I_2(p)}$ . It is given by (2.16) with  $\eta = \left( \begin{array}{cc}0 & -1\\ -1 & 0 \end{array} \right)$ . This choice of  $\eta$  ensures

$$
\left\{j \left(y _ {1}\right), j \left(y _ {2}\right) \right\} _ {\mathrm {P B}} = 2 y _ {1 2} \left(1 + \frac {1}{2} y _ {1 2} \partial_ {y _ {2}}\right) j \left(y _ {2}\right), \tag {4.17}
$$

which is the Poisson counterpart of the  $JJ$  OPE, with the same normalization conventions.

# 4.2.2 Free-field realization of all generators

We will now apply the procedure outlined in section 3 to obtain a free-field realization of the VOA  $\mathcal{W}_{I_2(p)}$ . Since  $I_{2}(p)$  has rank 2, we need two copies of the  $\beta \gamma bc$  system, denoted  $\beta_{1}$ ,  $\beta_{2}$ , and so on. We associate the label 1 to the invariant of degree 2,

$$
p _ {1} = 2, \quad p _ {2} = p. \tag {4.18}
$$

As explained in section 3.2 the only quantity that has to be constructed is  $J_{\mathrm{norms}}^{-}$  in (3.9). Following the steps given below (3.7a) one obtains the unique solution

$$
J _ {\text {n o r m s}} ^ {-} = \Lambda (\beta_ {1}) ^ {p - 2} \gamma_ {2} (\beta_ {1} \gamma_ {2} + (p - 1) b _ {1} c _ {2}). \tag {4.19}
$$

The normalization  $\Lambda$  could be scaled to one by the transformation (3.12) but it is instructive to keep it as a parameter. The remaining generators of the small  $\mathcal{N} = 4$  SCA take the form given in (2.8), (3.2) and (3.7a) together with  $\widetilde{G}^{-} = \{G^{-}J^{-}\}_{1}$ . For convenience of the reader, we summarize here the expression of all generators of the small  $\mathcal{N} = 4$  SCA in terms of the

two  $\beta \gamma bc$  systems:

$$
J ^ {+} = \beta_ {1},
$$

$$
J ^ {0} = b _ {1} c _ {1} + 2 \beta_ {1} \gamma_ {1} + (p - 1) b _ {2} c _ {2} + p \beta_ {2} \gamma_ {2},
$$

$$
\begin{array}{l} J ^ {-} = b _ {1} c _ {1} \gamma_ {1} + \beta_ {1} \gamma_ {1} \gamma_ {1} + (p - 1) \gamma_ {1} b _ {2} c _ {2} + p \gamma_ {1} \beta_ {2} \gamma_ {2} - c _ {1} b _ {2} \gamma_ {2} - (p + 1) \partial \gamma_ {1} \\ + \Lambda (\beta_ {1}) ^ {p - 1} (\gamma_ {2}) ^ {2} + (p - 1) \Lambda b _ {1} (\beta_ {1}) ^ {p - 2} c _ {2} \gamma_ {2}, \\ \end{array}
$$

$$
G ^ {+} = b _ {1},
$$

$$
G ^ {-} = b _ {1} \gamma_ {1} + b _ {2} \gamma_ {2},
$$

$$
\widetilde {G} ^ {+} = c _ {1} \partial \beta_ {1} + 2 \partial c _ {1} \beta_ {1} + (p - 1) c _ {2} \partial \beta_ {2} + p \partial c _ {2} \beta_ {2},
$$

$$
\begin{array}{l} \widetilde {G} ^ {-} = - b _ {1} \partial c _ {1} c _ {1} - c _ {1} b _ {2} \partial c _ {2} + c _ {1} \partial \beta_ {1} \gamma_ {1} + c _ {1} \partial \beta_ {2} \gamma_ {2} + (p - 1) \partial c _ {1} b _ {2} c _ {2} + 2 \partial c _ {1} \beta_ {1} \gamma_ {1} + p \partial c _ {1} \beta_ {2} \gamma_ {2} \\ + (p - 1) \gamma_ {1} c _ {2} \partial \beta_ {2} + p \gamma_ {1} \partial c _ {2} \beta_ {2} - (p + 1) \partial^ {2} c _ {1} \\ - (p - 1) \Lambda b _ {1} (\beta_ {1}) ^ {p - 2} \partial c _ {2} c _ {2} + 2 \Lambda (\beta_ {1}) ^ {p - 1} \partial c _ {2} \gamma_ {2} + (p - 1) \Lambda \partial \beta_ {1} (\beta_ {1}) ^ {p - 2} c _ {2} \gamma_ {2}, \\ \end{array}
$$

$$
T = - \frac {3}{2} b _ {1} \partial c _ {1} - \frac {1}{2} \partial b _ {1} c _ {1} - \beta_ {1} \partial \gamma_ {1} - \frac {p + 1}{2} b _ {2} \partial c _ {2} - \frac {p - 1}{2} \partial b _ {2} c _ {2} - \frac {p}{2} \beta_ {2} \partial \gamma_ {2} - \left(\frac {p}{2} - 1\right) \partial \beta_ {2} \gamma_ {2}. \tag {4.20}
$$

Given the small  $\mathcal{N} = 4$  SCA in terms of free fields as above, we proceed building the additional short generator  $W$  of the VOA associated to  $I_{2}(p)$ . According to our general prescription, we simply set

$$
W ^ {\mathrm {h . w .}} = \beta_ {2}, \tag {4.21}
$$

where h.w. stand for highest weight and refers to the component of  $W$  with charges  $h = m = p / 2$ . The whole  $\mathfrak{psl}(2|2)$  short supermultiplet  $\mathbb{W} = \{W,G_W,\widetilde{G}_W,T_W\}$ , see appendix A for our notation, is generated from  $W^{\mathrm{h.w.}}$ . We have now entirely specified our free-field realization. We refrain from giving the expressions for the other components of  $W$ , since their complexity grows quickly.

The next step is to verify that the set of strong generators of the VOA  $\mathcal{W}_{I_2(p)}$  is given by  $\mathbb{W} = \{W,G_W,\widetilde{G}_W,T_W\}$  together with the generators of the small  $\mathcal{N} = 4$  super-Virasoro algebra  $\mathbb{J} = \{J,G,\widetilde{G},T\}$ . To do so we need to close the OPE on this set of generators. The  $\mathbb{J}-\mathbb{J}$  and  $\mathbb{J}-\mathbb{W}$  OPEs take the required form by construction. The  $\mathbb{W}-\mathbb{W}$  OPEs are fixed in terms of the  $WW$  OPE. By means of a direct computation we verified that

$$
W \times W \sim g _ {W W} [ \mathrm {i d} ], \quad g _ {W W} = \frac {(2 p) !}{p ! p ^ {2}} \Lambda , \tag {4.22}
$$

where the notation [id] stands for the small  $\mathcal{N} = 4$  super-Virasoro family of the identity operator. The first few terms are given by

$$
[ \mathrm {i d} ] = \mathrm {i d} - \frac {6 p}{c} (J - \frac {1}{6} T) + \frac {1 8 p (p - 1)}{c (c - 6)} (J J) _ {0} ^ {2} - \frac {9 p (p + 1)}{c (c + 9)} ((J J) _ {0} ^ {0} + \frac {1}{3} T) + \dots \tag {4.23}
$$

with  $c = -6(p + 1)$ . It is worth remarking that the stress tensor  $T$  appears as  $\mathfrak{psl}(2|2)$  descendants of  $J$  and as completion of  $(JJ)_0^0$  to a  $\mathfrak{psl}(2|2)$  primary. Note that all factors  $z_{12}$ ,

$y_{12}$ , as well as all  $\partial_z$ ,  $\partial_y$  operators, are implicit, since they can be unambiguously restored exploiting  $\mathfrak{sl}(2)$  covariance. This compact notation for OPEs is described in more detail in appendix A, and is also utilized below in other examples. The formula (4.22) has been tested explicitly up to  $p = 7$ . As anticipated, the free parameter  $\Lambda$  entering (4.19) is crucial in order to obtain a viable realization of the full VOA.

# 4.2.3 Null states

Our proposed free-field realization is conjectured to enjoy the highly non-trivial property that all null states in the abstract VOA are realized manifestly as zero. In this section we check that the null states in the VOA  $\mathcal{W}_{I_2(p)}$  associated to the "Higgs Branch" relations (4.15) are identically zero in the free field realization. We expect that these nulls generated the maximal ideal of  $\mathcal{W}_{I_2(p)}$ , but at the moment we do not have a complete proof of this fact.

To begin with, let us consider the long composite operator linear in the extra generator  $W$ ,

$$
\mathfrak {L} _ {\frac {p}{2} + 1, \frac {p}{2} - 1} = (J W) _ {0} ^ {\frac {p}{2} - 1} + \frac {1}{p + 1} T _ {W}. \tag {4.24}
$$

The first term in the expression above is the quasiprimary completion of the normal ordered product  $JW$  projected on the  $\mathfrak{sl}(2)_y$  spin  $\frac{p}{2} - 1$  component. The second term contains the supersymmetry descendant of  $W$ , namely  $T_W = -G^\downarrow \widetilde{G}^\downarrow W$ , see (2.10) for the notation. This term is necessary in order for (4.24) to be a  $\mathfrak{psl}(2|2)$  primary. Being linear in the new generator  $W$ , the operator  $\mathfrak{L}_{\frac{p}{2} + 1, \frac{p}{2} - 1}$  can be defined for any value of the central charge without any reference to the free field realization and by construction is a  $\mathcal{N} = 4$  super-Virasoro descendant of  $W$  itself. If the central charge takes the special value  $c = -6(p + 1)$  the operator  $\mathfrak{L}_{\frac{p}{2} + 1, \frac{p}{2} - 1}$  becomes an  $\mathcal{N} = 4$  super-Virasoro primary operator. Being a primary and a descendant at the same time it must be null. In our free-field realization the operator  $\mathfrak{L}_{\frac{p}{2} + 1, \frac{p}{2} - 1}$  vanishes identically. We have thus recovered the VOA counterpart of the first "Higgs Branch" relation in (4.15). Let us turn to the second set of relations in (4.15). Let us consider the composite operators

$$
\mathfrak {L} _ {p, j} ^ {W W} = (W W) _ {0} ^ {j} + (\text {s u p e r - V i r a s o r o p r i m a r y c o m p l e t i o n}), \quad j <   p, \quad p - j \text {e v e n}. \tag {4.25}
$$

In other words,  $\mathfrak{L}_{p,j}^{WW}$  is defined to be the  $\mathcal{N} = 4$  super-Virasoro primary completion of the normal ordered product  $WW$  projected onto the component with  $\mathfrak{sl}(2)_y$  spin  $j$ . The requirement that  $p - j$  be even stems from Bose symmetry. One can verify that this definition is well-posed, in the sense that, making only use of the OPEs of the abstract VOA, one can check that there exists a unique super-Virasoro primary operator starting with  $(WW)_0^j$ . Having unambiguously defined the composite  $\mathfrak{L}_{p,j}^{WW}$  in the abstract VOA, we can resort to our free-field realization and verify, in a few examples, that this object is indeed identically vanishing. This finding is in perfect agreement with the bootstrap analysis of [33].

# 4.2.4 Classical limit: relation between  $\beta \gamma$  and  $z^{\pm}$

Using the map  $\mathcal{P}$  defined in section 3.4, we can define the classical objects associated to the generators  $J$  and  $W$ ,

$$
J _ {\mathrm {c l}} := \mathcal {P} (J) , \quad W _ {\mathrm {c l}} := \mathcal {P} (W) , \tag {4.26}
$$

more explicitly

$$
J _ {\mathrm {c l}} ^ {+} = \beta_ {1}, J _ {\mathrm {c l}} ^ {0} = 2 \beta_ {1} \gamma_ {1} + p \beta_ {2} \gamma_ {2}, J _ {\mathrm {c l}} ^ {-} = (\beta_ {1} \gamma_ {1} + p \beta_ {2} \gamma_ {2}) \gamma_ {1} + \Lambda (\beta_ {1}) ^ {p - 1} (\gamma_ {2}) ^ {2}, (4. 2 7)
$$

and  $W_{\mathrm{cl}} = \beta_2 +$  descendants. Recall that, as explained in section 3.4,  $\beta \gamma$  are now commuting variables. The combinations (4.26) satisfy the same relations (4.15) as (4.14). This implies that  $J_{\mathrm{cl}}, W_{\mathrm{cl}}$  provide a realization of the ring  $\mathcal{R}_{I_2(p)}$  as a subring of  $\mathbb{C}[\beta_1, \gamma_1, \beta_2, \gamma_2]$ . It is instructive to determine  $\beta_{1,2}, \gamma_{1,2}$  in terms of the quotient variables  $z_{1,2}^{\pm}$  by equating (4.14) with (4.26). This gives the remarkably simple expressions

$$
\beta_ {1} = z _ {1} ^ {+} z _ {2} ^ {+}, \gamma_ {1} = \frac {(z _ {1} ^ {+}) ^ {p - 1} z _ {1} ^ {-} - (z _ {2} ^ {+}) ^ {p - 1} z _ {2} ^ {-}}{(z _ {1} ^ {+}) ^ {p} - (z _ {2} ^ {+}) ^ {p}},
$$

$$
\beta_ {2} = \left(z _ {1} ^ {+}\right) ^ {p} + \left(z _ {2} ^ {+}\right) ^ {p}, \quad \gamma_ {2} = \frac {1}{p} \frac {z _ {1} ^ {+} z _ {2} ^ {-} - z _ {2} ^ {+} z _ {1} ^ {-}}{\left(z _ {1} ^ {+}\right) ^ {p} - \left(z _ {2} ^ {+}\right) ^ {p}}, \quad \Lambda = p ^ {2}. \tag {4.28}
$$

Notice that  $\beta_{1,2}$ ,  $\gamma_{1,2}$  are rational functions that are invariant under the action (4.11) of  $I_2(p)$ . The Poisson brackets (2.16) with  $\eta = \left( \begin{array}{cc}0 & -1\\ -1 & 0 \end{array} \right)$  imply the expected Poisson brackets

$$
\left\{\beta_ {\ell_ {1}}, \gamma_ {\ell_ {2}} \right\} _ {\mathrm {P B}} = - \delta_ {\ell_ {1}, \ell_ {2}}, \quad \left\{\beta_ {\ell_ {1}}, \beta_ {\ell_ {2}} \right\} _ {\mathrm {P B}} = \left\{\gamma_ {\ell_ {1}}, \gamma_ {\ell_ {2}} \right\} _ {\mathrm {P B}} = 0, \quad \ell_ {1}, \ell_ {2} = 1, 2. \tag {4.29}
$$

The minus sign in the first equation is a consequence of our conventions for the  $\beta \gamma$  OPEs, see (3.1).

# 4.2.5 Comments on the screening operator

It is natural to ask if  $\mathcal{W}_{I_2(p)}$  can be identified with the kernel of a suitable screening operator acting on the free field VOA  $\mathbb{M}_{\beta \gamma bc}^{(I_2(p))}$ . A simpler version of this problem is obtained using the map  $\mathcal{P}'$  of section 3.4. More precisely, we aim at identifying  $\mathcal{P}'(\mathcal{W}_{I_2(p)})$  with the kernel of a suitable object  $\mathsf{J}_{\mathrm{cl}}$  in the classical Poisson superalgebra  $\mathbb{M}_{\beta \gamma bc}^{(I_2(p))\mathrm{cl}}$ . The object  $\mathsf{J}_{\mathrm{cl}}$  acts via Poisson bracket.

The object  $\mathsf{J}_{\mathrm{cl}}$  can be presented as

$$
\mathsf {J} _ {\mathrm {c l}} = b _ {1} \left(\beta_ {1}\right) ^ {- \frac {1}{2}} \left[ \frac {1}{2} F (x) - \frac {1}{2} p x F ^ {\prime} (x) \right] + b _ {2} \left(\beta_ {1}\right) ^ {\frac {1}{2} - \frac {p}{2}} F ^ {\prime} (x) = \{\mathcal {P} ^ {\prime} (G ^ {-}), \mathsf {K} _ {\mathrm {c l}} \} _ {\mathrm {P B}}, \tag {4.30}
$$

where in the last step we introduced the auxiliary object

$$
\mathrm {K} _ {\mathrm {c l}} = \left(\beta_ {1}\right) ^ {\frac {1}{2}} F (x), \quad x := \left(\beta_ {1}\right) ^ {- \frac {p}{2}} \beta_ {2}. \tag {4.31}
$$

The function  $F(x)$  is required to be a solution to the differential equation

$$
\left(p ^ {2} x ^ {2} - 4 \Lambda\right) F ^ {\prime \prime} (x) + p ^ {2} x F ^ {\prime} (x) - F (x) = 0. \tag {4.32}
$$

We have checked that, by virtue of the above equation, one has

$$
\{\mathrm {J} _ {\mathrm {c l}}, \mathcal {P} ^ {\prime} (X) \} _ {\mathrm {P B}} = 0 \quad \text {f o r} \quad X \in \{J ^ {+}, J ^ {0}, J ^ {-}, G ^ {+}, G ^ {-}, \widetilde {G} ^ {+}, \widetilde {G} ^ {-}, T, W ^ {\mathrm {h . w .}} \}. \tag {4.33}
$$

This is enough to guarantee that  $\mathcal{P}'(\mathcal{W}_{I_2(p)})$  lies inside the kernel of  $\mathsf{J}_{\mathrm{cl}}$  acting of the Poisson superalgebra  $\mathbb{M}_{\beta \gamma bc}^{(I_2(p))\mathrm{cl}}$ . It seems natural to conjecture that  $\mathcal{P}'(\mathcal{W}_{I_2(p)})$  is actually the entirety of the kernel of  $\mathsf{J}_{\mathrm{cl}}$ , but we do not have a proof of this fact.

In order to promote the results of the previous paragraphs from the level of the Poisson algebra to the level of the full VOA, we have to be able to make sense of expressions like (4.31) in the context of the VOA. This is possible by expressing  $\beta_{1}$ ,  $\beta_{2}$  in terms of chiral bosons and making use of vertex operators. We refrain, however, from pursuing this direction further.

# 4.3 Rank 3:  $\Gamma = A_3, B_3, H_3$

We will now present the free field construction of  $\mathcal{W}_{\Gamma}$  for  $\mathrm{rank}(\Gamma) = 3$ . There are only three examples in this case, namely  $\Gamma = A_3, B_3, H_3$ . The structure of these VOA is analyzed in less details compared to the rank 2 series discussed in the previous section. A few remarks are in order. Notice that in all three cases the degrees are  $(2,3,4)$ ,  $(2,4,6)$ ,  $(2,6,10)$  so that  $p_3 = 2p_2 - 2$ . Moreover all super-Virasoro primary operators of the form  $(W_{\ell_1}W_{\ell_2} + \ldots)_{\mathfrak{L}}$  are null except for  $(W_{2}W_{3} + \ldots)_{\mathfrak{L}(\frac{p_2 + p_3}{2},\frac{p_2 + p_3}{2} -1)}$ . Finally, as in the rank two series, the classical limit of  $J, W_2, W_3$  obtained by applying the map  $\mathcal{P}$  defined in section 3.4, gives a realization of the rings  $\mathcal{R}_{A_3}$ ,  $\mathcal{R}_{B_3}$  and  $\mathcal{R}_{H_3}$  as subrings of  $\mathbb{C}[\beta_1,\gamma_1,\beta_2,\gamma_2,\beta_3,\gamma_3]$ . This properties is not manifest but has been checked by verifying that all the relations are satisfied<sup>18</sup>.

Notation: In order to make the equations easier to read we will label  $W$  generators as well as  $\beta \gamma bc$  by their weight with the gothic suffix  $\mathfrak{p} \in \{3,4,5,6,7,\ldots\}$ .

# 4.3.1 Example:  $\Gamma = A_{3}$

Free field realization By following the steps described in section 3 we find a unique solution for  $J_{\mathrm{norms}}^-$  up to the rescaling (3.12). Its explicit form is rather long so we present only its classical limit:

$$
\mathcal {P} \left(J _ {\text {n o r m s}} ^ {-}\right) = \Lambda_ {2} \left(\Lambda_ {1} \beta_ {2} ^ {2} - \beta_ {4}\right) \gamma_ {3} ^ {2} - \frac {1 6 \Lambda_ {1}}{3} \beta_ {2} \beta_ {3} \gamma_ {3} \gamma_ {4} + \frac {\Lambda_ {1}}{1 2 \Lambda_ {2}} \left(2 0 \Lambda_ {1} \Lambda_ {2} \beta_ {2} ^ {3} + 5 1 \beta_ {3} ^ {2} + 2 8 \Lambda_ {2} \beta_ {2} \beta_ {4}\right) \gamma_ {4} ^ {2}, \tag {4.34}
$$

where  $\mathcal{P}$  is defined in section 3.4. Given  $J_{\mathrm{norms}}^-$  we can construct the remaining strong generators as descendants of  $\beta_{2}$  and  $\beta_{3}$ . The parameters  $\Lambda_1, \Lambda_2$  are related to the normalization appearing below as  $g_{33} = \frac{85}{2}\Lambda_1\Lambda_2$ ,  $g_{44} = 595\Lambda_1^2$ .

Closing the OPE. Let us present the OPE of strong generators in this case:

$$
W _ {3} \times W _ {3} \sim g _ {3 3} [ \mathrm {i d} ] + \lambda_ {3 3} ^ {4} [ W _ {4} ] \tag {4.35a}
$$

$$
W _ {3} \times W _ {4} \sim \lambda_ {3 4} ^ {3} [ W _ {3} ] \tag {4.35b}
$$

$$
W _ {4} \times W _ {4} \sim g _ {4 4} [ \mathrm {i d} ] + \lambda_ {4 4} ^ {4} [ W _ {4} ] + \lambda_ {4 4} ^ {(3 3)} [ (W _ {3}) _ {\mathfrak {S}} ^ {2} ] \tag {4.35c}
$$

where

$$
\frac {\sqrt {g _ {4 4}}}{g _ {3 3}} \lambda_ {3 3} ^ {4} = \frac {1}{\sqrt {g _ {4 4}}} \lambda_ {3 4} ^ {3} = - 4 \sqrt {\frac {7}{8 5}}, \quad \frac {1}{\sqrt {g _ {4 4}}} \lambda_ {4 4} ^ {4} = \frac {1 1}{3} \sqrt {\frac {5}{1 1 9}}, \quad \frac {g _ {3 3}}{g _ {4 4}} \lambda_ {4 4} ^ {(3 3)} = \frac {1 7}{2 8}, \tag {4.36}
$$

and

$$
\left(W _ {3}\right) _ {\mathfrak {S}} ^ {2} = \left(\beta_ {3}\right) ^ {2} - \frac {8 g _ {3 3}}{2 3 \sqrt {g _ {4 4}}} \sqrt {\frac {7}{8 5}} \beta_ {2} \beta_ {4} + \frac {8 g _ {3 3}}{4 8 4 5} \left(\beta_ {2}\right) ^ {3} + \text {d e s c e n d a n t s}, \tag {4.37}
$$

is the completion of  $(\beta_{3})^{2}$  to a super-Virasoro primary. Its norm is given by  $\frac{1344}{437} g_{33}^{2}$ . As usual,  $[X]$  denotes the contribution from the  $\mathcal{N} = 4$  super-Virasoro family of the primary  $X$ . Notice that there is a null state of type  $\mathfrak{L}_{(3,1)}$  of the schematic form  $W_{3}W_{3} + JW_{4} + J^{3} + \ldots$ . This is a primary operator that, if not null, could appear in the right hand side of (4.35c). The quantum numbers of the relations among the generators of  $\mathcal{R}_{A_3}$  are

$$
\mathfrak {L} _ {(3, 1)}, \quad \mathfrak {L} _ {(\frac {7}{2}, \frac {3}{2})}, \quad \mathfrak {L} _ {(\frac {7}{2}, \frac {1}{2})}, \quad \mathfrak {L} _ {(4, 2)}, \quad \mathfrak {L} _ {(4, 0)}. \tag {4.38}
$$

They all correspond to null operators in the VOA. Notice that all the operators of the type  $(W_{\mathfrak{p}_1}W_{\mathfrak{p}_2} + \ldots)_{\mathfrak{L}}$  are null except for  $(W_3W_4 + \ldots)_{\mathfrak{L}_{(\frac{7}{2},\frac{5}{2})}}$ .

# 4.3.2 Example:  $\Gamma = B_{3}$

Free field realization As before, following the recipe given in section 3 we find a unique solution for  $J_{\mathrm{norms}}^{-}$ . In the classical limit it reads

$$
\begin{array}{l} \mathcal {P} \left(J _ {\text {n o r m s}} ^ {-}\right) = \left(u _ {1} \beta_ {2} ^ {3} + u _ {2} \beta_ {2} \beta_ {4} + u _ {3} \beta_ {6}\right) \gamma_ {4} ^ {2} + (4.39a) \\ + \left(u _ {4} \beta_ {2} ^ {4} + u _ {5} \beta_ {2} ^ {2} \beta_ {4} + u _ {6} \beta_ {4} ^ {2} + u _ {7} \beta_ {2} \beta_ {6}\right) \gamma_ {4} \gamma_ {6} + (4.39b) \\ + \left(u _ {8} \beta_ {2} ^ {5} + u _ {9} \beta_ {2} ^ {3} \beta_ {4} + u _ {1 0} \beta_ {2} \beta_ {4} ^ {2} + u _ {1 1} \beta_ {2} ^ {2} \beta_ {6} + u _ {1 2} \beta_ {4} \beta_ {6}\right) \gamma_ {6} ^ {2} (4.39c) \\ \end{array}
$$

The explicit form of the coefficients  $u_{k}$  is not very illuminating so we omit it. With this ingredient we can produce all the strong generators.

Closing the OPE. Let us present the OPE of strong generators in this case, setting the normalizations to one,

$$
W _ {4} \times W _ {4} \sim [ \mathrm {i d} ] + \lambda_ {4 4} ^ {4} \left[ W _ {4} \right] + \lambda_ {4 4} ^ {6} \left[ W _ {6} \right] \tag {4.40a}
$$

$$
W _ {4} \times W _ {6} \sim \lambda_ {4 6} ^ {4} \left[ W _ {4} \right] + \lambda_ {4 6} ^ {6} \left[ W _ {6} \right] + \lambda_ {4 6} ^ {(4 4)} \left[ \left(W _ {4}\right) _ {\mathfrak {S}} ^ {2} \right] \tag {4.40b}
$$

$$
W _ {6} \times W _ {6} \sim [ \mathrm {i d} ] + \lambda_ {6 6} ^ {4} [ W _ {4} ] + \lambda_ {6 6} ^ {6} [ W _ {6} ] + \lambda_ {6 6} ^ {(4 4)} [ (W _ {4}) _ {\mathfrak {S}} ^ {2} ] + \lambda_ {6 6} ^ {(6 4)} [ (W _ {6} W _ {4}) _ {\mathfrak {S}} ] \tag {4.40c}
$$

where

$$
\lambda_ {4 4} ^ {4} = - \frac {1 4 3}{3 \sqrt {2 4 1 5}}, \quad \lambda_ {4 4} ^ {6} = \lambda_ {4 6} ^ {4} = \frac {4 6}{5} \sqrt {\frac {1 1}{6 0 9}}, \quad \lambda_ {4 6} ^ {(4 4)} = - 1 2 \sqrt {\frac {5}{7 3 3 7}}, \tag {4.41}
$$

$$
\lambda_ {4 6} ^ {6} = \lambda_ {6 6} ^ {4} = \frac {6 5}{5 8} \sqrt {\frac {3 5}{6 9}}, \quad \lambda_ {6 6} ^ {6} = - \frac {2 7}{2 9} \sqrt {\frac {2 1}{3 1 9}}, \lambda_ {6 6} ^ {(4 4)} = \frac {2 8 5}{3 1 9}, \quad \lambda_ {4 4} ^ {(6 4)} = 2 \sqrt {\frac {5 5}{6 6 7}}, \tag {4.42}
$$

and  $(W_4)_{\mathfrak{S}}^2$ ,  $(W_6W_4)_{\mathfrak{S}}$  are defined in a similar way to (4.37) so that they are  $\mathcal{N} = 4$  super-Virasoro primaries. Notice that there are nulls of type  $\mathfrak{L}_{(4,2)}$  and  $\mathfrak{L}_{(4,0)}$  which are the unique super-Virasoro primary completions of  $(W_4)_{\mathfrak{L}}^2$ . These are primary operators that, if not null,

could appear in the right hand side of (4.40b). Similarly, there are nulls of type  $\mathfrak{L}_{(5,3)}$ ,  $\mathfrak{L}_{(5,2)}$ ,  $\mathfrak{L}_{(5,1)}$  relevant for the OPE (4.40c) of the schematic form  $(W_4W_6)_{\mathfrak{L}}$ . The quantum numbers of the relations are

$$
\mathfrak {L} _ {(4, 2)}, \quad \mathfrak {L} _ {(4, 0)}, \quad \mathfrak {L} _ {(5, 3)}, \quad \mathfrak {L} _ {(5, 2)}, \quad \mathfrak {L} _ {(5, 1)}, \quad \mathfrak {L} _ {(6, 4)}, \quad \mathfrak {L} _ {(6, 2)}, \quad \mathfrak {L} _ {(6, 0)}. \tag {4.43}
$$

As in the case of  $A_{3}$  all the super-Virasoro primary operators of the form  $(W_{\mathfrak{p}_1}W_{\mathfrak{p}_2} + \ldots)_{\mathfrak{L}}$  are null except for  $(W_6W_4 + \dots)_{\mathfrak{L}_{(5,3)}}$ .

# 4.3.3 Example:  $\Gamma = H_{3}$

Free field realization Also in this case the procedure outlined in section 3 gives a unique solution for  $J_{\mathrm{norms}}^{-}$ . In the classical limit it reads

$$
\begin{array}{l} \mathcal {P} \left(J _ {\text {n o r m s}} ^ {-}\right) = \left(u _ {1} \beta_ {2} ^ {5} + u _ {2} \beta_ {2} ^ {2} \beta_ {6} + u _ {3} \beta_ {1 0}\right) \gamma_ {6} ^ {2} + (4.44a) \\ + \left(u _ {4} \beta_ {2} ^ {7} + u _ {5} \beta_ {2} ^ {4} \beta_ {6} + u _ {6} \beta_ {2} \beta_ {6} ^ {2} + u _ {7} \beta_ {2} ^ {2} \beta_ {1 0}\right) \gamma_ {6} \gamma_ {1 0} + (4.44b) \\ + \left(u _ {8} \beta_ {2} ^ {9} + u _ {9} \beta_ {2} ^ {6} \beta_ {6} + u _ {1 0} \beta_ {2} ^ {3} \beta_ {6} ^ {2} + u _ {1 1} \beta_ {6} ^ {3} + u _ {1 2} \beta_ {2} ^ {4} \beta_ {1 0} + u _ {1 3} \beta_ {2} \beta_ {6} \beta_ {1 0}\right) \gamma_ {1 0} ^ {2}, (4.44c) \\ \end{array}
$$

where we omit the explicit form of the coefficients  $u_{k}$ .

Closing the OPE. Let us present the OPE of strong generators in this case:

$$
W _ {6} \times W _ {6} \sim [ \mathrm {i d} ] + \lambda_ {6 6} ^ {6} \left[ W _ {6} \right] + \lambda_ {6 6} ^ {1 0} \left[ W _ {1 0} \right] \tag {4.45a}
$$

$$
W _ {6} \times W _ {1 0} \sim \lambda_ {6, 1 0} ^ {6} \left[ W _ {6} \right] + \lambda_ {6, 1 0} ^ {1 0} \left[ W _ {1 0} \right] + \lambda_ {6, 1 0} ^ {(6 6)} \left[ \left(W _ {6}\right) _ {\mathfrak {S}} ^ {2} \right] \tag {4.45b}
$$

$$
\begin{array}{l} W _ {1 0} \times W _ {1 0} \sim [ \mathrm {i d} ] + \lambda_ {1 0, 1 0} ^ {6} \left[ W _ {6} \right] + \lambda_ {1 0, 1 0} ^ {1 0} \left[ W _ {1 0} \right] + \lambda_ {1 0, 1 0} ^ {(6 6)} \left[ \left(W _ {6}\right) _ {\mathfrak {S}} ^ {2} \right] + \tag {4.45c} \\ + \lambda_ {1 0, 1 0} ^ {(6, 1 0)} \left[ \left(W _ {6} W _ {1 0}\right) _ {\mathfrak {S}} \right] + \lambda_ {1 0, 1 0} ^ {(6 6 6)} \left[ \left(W _ {6}\right) _ {\mathfrak {S}} ^ {3} \right], \\ \end{array}
$$

where

$$
\lambda_ {6 6} ^ {6} = \frac {5 7}{2} \sqrt {\frac {3 5}{1 0 5 8 2}}, \quad \lambda_ {6 6} ^ {1 0} = \lambda_ {6 1 0} ^ {6} = \frac {3 7}{2} \sqrt {\frac {4 8 4 5}{5 5 1 1 2 2}}, \quad \lambda_ {6, 1 0} ^ {1 0} = \lambda_ {1 0, 1 0} ^ {6} = - \frac {1 1 8 9}{9 4} \sqrt {\frac {5 5}{6 7 3 4}}, \tag {4.46}
$$

$$
\lambda_ {6, 1 0} ^ {(6 6)} = \frac {8 1 7}{2} \sqrt {\frac {4 1}{1 1 7 9 5 6 3 7}}, \quad \lambda_ {1 0, 1 0} ^ {1 0} = \frac {2 7 4 7 3 5 7}{1 9 7 4} \sqrt {\frac {5 5}{4 8 5 4 8 8 3 8}}, \quad \lambda_ {1 0, 1 0} ^ {(6 6)} = \frac {4 7 9 1 6 7}{2 9 5 6 3 0}, \tag {4.47}
$$

$$
\lambda_ {1 0, 1 0} ^ {(6, 1 0)} = - \frac {2 3 2}{3} \sqrt {\frac {2 8 7}{1 6 8 5 0 9 1}}, \quad \lambda_ {1 0, 1 0} ^ {(6 6 6)} = \frac {1 5 4 1 6}{3 5 8 5 3} \sqrt {\frac {2 8 6}{1 2 9 5}}. \tag {4.48}
$$

In this case the defining relations for the ring  $\mathcal{R}_{H_3}$  are of type

$$
\left\{\mathfrak {L} _ {(6, j)} \right\} _ {j \in \{0, 2, 4 \}}, \quad \left\{\mathfrak {L} _ {(8, j)} \right\} _ {j \in \{2, 3, 4, 5, 6 \}}, \quad \left\{\mathfrak {L} _ {(1 0, j)} \right\} _ {j \in \{0, 2, 4, 6, 8 \}}. \tag {4.49}
$$

Notice that the super-Virasoro primary operators of the form  $(W_6)_{\mathfrak{L}}^3$  are null as a consequence of  $(W_6)_{\mathfrak{L}}^2$  being null.

# 4.4 A rank 4 example:  $\Gamma = D_4$

In this example we will encounter two new features: (1) the ring of invariants  $\mathcal{R}_{D_4}$  and, according to our proposal, the VOA  $\mathcal{W}_{\Gamma}$  has long generators, (2) the form of  $J^{-}$  is not uniquely determined by the first four steps given in section (3.2) and the fifth condition needs to be incorporated.

Generators of the ring  $\mathcal{R}_{D_4}$ . The action of the Weyl group of type  $D_4$  on  $\mathbb{R}^4$  is generated by

$$
s _ {1}: \quad \left(z _ {1}, z _ {2}, z _ {3}, z _ {4}\right) \mapsto \left(- z _ {1}, z _ {2}, z _ {3}, z _ {4}\right) \tag {4.50}
$$

$$
s _ {2}: \quad \left(z _ {1}, z _ {2}, z _ {3}, z _ {4}\right) \mapsto \frac {1}{2} \left(z _ {1 2} ^ {+} - z _ {3 4} ^ {+}, z _ {1 2} ^ {+} + z _ {3 4} ^ {+}, - z _ {1 2} ^ {-} + z _ {3 4} ^ {-}, - z _ {1 2} ^ {-} - z _ {3 4} ^ {-}\right) \tag {4.51}
$$

$$
s _ {3}: \quad \left(z _ {1}, z _ {2}, z _ {3}, z _ {4}\right) \mapsto \left(z _ {1}, z _ {2}, - z _ {3}, z _ {4}\right) \tag {4.52}
$$

$$
s _ {4}: \quad \left(z _ {1}, z _ {2}, z _ {3}, z _ {4}\right) \mapsto \left(z _ {1}, z _ {2}, z _ {3}, - z _ {4}\right) \tag {4.53}
$$

where  $z_{ij}^{\pm} = z_i \pm z_j$ . In this basis, a choice of generators for the ring of invariants  $\mathbb{C}[z_1, \ldots, z_4]^{D_4}$  is given by

$$
\mathcal {I} _ {2} = \frac {1}{2} \left(z _ {1} ^ {2} + z _ {2} ^ {2} + z _ {3} ^ {2} + z _ {4} ^ {2}\right), \quad \mathcal {I} _ {6} = \left(z _ {1} ^ {2} z _ {2} ^ {2} - z _ {3} ^ {2} z _ {4} ^ {2}\right) \left(z _ {1} ^ {2} + z _ {2} ^ {2} - z _ {3} ^ {2} - z _ {4} ^ {2}\right) \tag {4.54a}
$$

$$
\mathcal {I} _ {4} ^ {(1)} = \left(z _ {2} ^ {2} - z _ {3} ^ {2}\right) \left(z _ {1} ^ {2} - z _ {4} ^ {2}\right), \quad \mathcal {I} _ {4} ^ {(2)} = \left(z _ {1} ^ {2} - z _ {3} ^ {2}\right) \left(z _ {2} ^ {2} - z _ {4} ^ {2}\right). \tag {4.54b}
$$

These generator are algebraically independent. As discussed in section (2.2) we need to analyze the ring (2.12) in which the Coxeter group acts on two copies of  $z$ , called  $z^{\pm}$ . Some of the generators of this ring are immediately identified starting from (4.54) and promoting each  $z_i$  to  $z_i(y) = z_i^+ + yz_i^-$ . These are the so-called short generators. As opposed to all the examples encountered so far, namely  $A_1, I_2(p), A_3, B_3, H_3$ , this is the first example in which these are not all the generators of  $\mathcal{R}_{D_4}$ . The missing generator is given by

$$
w _ {(3, 0)} := X _ {1 2 3} + X _ {1 3 4} - X _ {1 2 4} - X _ {2 3 4}, \quad X _ {i j k} := \langle i j \rangle \langle i k \rangle \langle j k \rangle , \tag {4.55}
$$

where  $\langle ij\rangle = \epsilon_{IJ}z_i^I z_j^J$ . It is rather clear that this invariant cannot be written as a composite of the short generators. What is less obvious is that there is no additional generator. We claim that it is the case. This fact can be checked by matching the Hilbert series computed from the proposed description of  $\mathcal{R}_{D_4}$  in terms of generators and relations with the Molien series obtained from the quotient description. Alternatively one can verify that the set of generators we propose closes under the Poisson bracket (2.16). We followed the second strategy.

Finally let us explain how the triality automorphism of  $D_4$  acts on the space of invariants. Its action in terms of the  $z_{i=1,\ldots,4}$  variables can be defined as the group  $S_3$  of permutations of  $z_1, z_2, z_3$ . It is a simple exercise to verify that, up to a simple redefinition of  $\mathcal{I}_6$ , the invariants of degree 2 and 6 transform trivially under  $S_3$ , the invariants of degree 4 transform in a two-dimensional representation<sup>19</sup> and the long invariant  $w_{(3,0)}$  transform in the non-trivial one-dimensional representation corresponding to the sign of the permutation.

Free field realization. This is the first example in which the first four steps described in section (3.2) to construct  $J^{-}$  do not give a unique result. There is a number or relatively simple conditions that the free field realization must satisfy that are easy to add. The first one is the following. Consider the OPE  $W_{\ell_1}W_{\ell_2}$ . The term with next to extremal  $\mathfrak{sl}(2)_y$  spin  $j = \frac{1}{2} (p_{\ell_1} + p_{\ell_2}) - 1$  in the singular part of these OPE can appear only in the first order pole.

By a simple quantum number analysis this term must be a short operator. We require that this object is indeed a composite operator of the short generators that have been postulated in the free field realization. The second condition is given by focusing on the first two most singular terms in the WW OPE. Their form is fixed by super-Virasoro symmetry to be

$$
W _ {\ell_ {1}} (z _ {1}) W _ {\ell_ {2}} (z _ {1}) = \delta_ {\ell_ {1}, \ell_ {2}} \left(\mathrm {i d} - \frac {6 p _ {\ell_ {1}}}{c} J\right) + \dots \tag {4.56}
$$

This form of the OPE is added as an extra requirement. We found experimentally that these two conditions allow to completely fix the free field realization. It is possible that in more complicated examples more conditions need to be added, but we expect that the Ansatz we propose is sufficient.

In this case we obtain<sup>20</sup>

$$
\begin{array}{l} \mathcal {P} \left(J _ {\text {n o r m s}} ^ {-}\right) = \Lambda_ {1} \beta_ {2} ^ {3} \left(\left(\gamma_ {4} ^ {+}\right) ^ {2} + \left(\gamma_ {4} ^ {-}\right) ^ {2}\right) + 4 \sqrt {\frac {\Lambda_ {1}}{3}} \beta_ {2} \left(\beta_ {4} ^ {+} \left(\left(\gamma_ {4} ^ {-}\right) ^ {2} - \left(\gamma_ {4} ^ {+}\right) ^ {2}\right) + 2 \beta_ {4} ^ {-} \gamma_ {4} ^ {+} \gamma_ {4} ^ {-}\right) + (4. 5 7 a) \\ + \left(\frac {1 6 \Lambda_ {2}}{\Lambda_ {1}} \beta_ {2} \left((\beta_ {4} ^ {+}) ^ {2} + (\beta_ {4} ^ {-}) ^ {2}\right) + \Lambda_ {2} \beta_ {2} ^ {5}\right) \gamma_ {6} ^ {2} + 2 8 \sqrt {\frac {\Lambda_ {2}}{1 5}} \beta_ {2} ^ {2} \left(\beta_ {4} ^ {+} \gamma_ {4} ^ {+} + \beta_ {4} ^ {-} \gamma_ {4} ^ {-}\right) \gamma_ {6} + (4. 5 7 b) \\ + \frac {1}{\sqrt {1 5 \Lambda_ {2}}} \beta_ {6} \left(5 \Lambda_ {1} \left(\left(\gamma_ {4} ^ {+}\right) ^ {2} + \left(\gamma_ {4} ^ {-}\right) ^ {2}\right) - 2 2 \Lambda_ {2} \beta_ {2} ^ {2} \gamma_ {6} ^ {2}\right) + (4.57c) \\ + 1 6 \sqrt {\frac {\Lambda_ {2}}{5 \Lambda_ {1}}} \left(\left(\left(\beta_ {4} ^ {-}\right) ^ {2} - \left(\beta_ {4} ^ {+}\right) ^ {2}\right) \gamma_ {4} ^ {+} + 2 \beta_ {4} ^ {+} \beta_ {4} ^ {-} \gamma_ {4} ^ {-}\right) \gamma_ {6}. (4.57d) \\ \end{array}
$$

The parameters  $\Lambda_1, \Lambda_2$  are related to the norms of the short generators as  $g_{44} = 1680\Lambda_1$ ,  $g_{66} = 665280\Lambda_2$ .

Closing the OPE The extra generators in this case are  $W_4^{\pm}, W_6$  and  $W_{(3,0)}$  and their OPEs take the form

$$
W _ {4} ^ {\pm} \times W _ {4} ^ {\pm} \sim g _ {4 4} [ \mathrm {i d} ] \pm \lambda_ {4 4} ^ {4 ^ {+}} \left[ W _ {4} ^ {+} \right] + \lambda_ {4 4} ^ {6} \left[ W _ {6} \right] \tag {4.58a}
$$

$$
W _ {4} ^ {+} \times W _ {4} ^ {-} \sim \lambda_ {4 ^ {+} 4 ^ {-}} ^ {4 ^ {-}} \left[ W _ {4} ^ {-} \right] + \left[ W _ {(3, 0)} \right] \tag {4.58b}
$$

$$
W _ {4} ^ {+} \times W _ {6} \sim \lambda_ {4 6} ^ {4} \left[ W _ {4} ^ {+} \right] + \lambda_ {4 6} ^ {(4 4)} \left[ \left(W _ {4}\right) _ {\mathfrak {S}} ^ {2, +} \right] + \dots \tag {4.58c}
$$

$$
W _ {4} ^ {-} \times W _ {6} \sim \lambda_ {4 6} ^ {4} \left[ W _ {4} ^ {-} \right] + \lambda_ {4 6} ^ {(4 4)} \left[ \left(W _ {4}\right) _ {\mathfrak {S}} ^ {2, -} \right] + \dots \tag {4.58d}
$$

$$
W _ {6} \times W _ {6} \sim g _ {6 6} [ \mathrm {i d} ] + \lambda_ {6 6} ^ {6} \left[ W _ {6} \right] + \lambda_ {6 6} ^ {(4 4) C} \left[ \left(W _ {4}\right) _ {\mathfrak {S}} ^ {2, 0} \right] + \dots \tag {4.58e}
$$

where

$$
- \frac {1}{\sqrt {g _ {4 4}}} \lambda_ {4 4} ^ {4 ^ {+}} = \frac {1}{\sqrt {g _ {4 4}}} \lambda_ {4 ^ {+} 4 ^ {-}} ^ {4 ^ {-}} = \frac {6}{\sqrt {3 5}}, \quad \frac {\sqrt {g _ {6 6}}}{g _ {4 4}} \lambda_ {4 4} ^ {6} = \sqrt {\frac {1 1}{7}}, \quad \frac {1}{\sqrt {g _ {6 6}}} \lambda_ {4 6} ^ {4} = \sqrt {\frac {1 1}{7}}, \tag {4.59}
$$

$$
\sqrt {\frac {g _ {4 4}}{g _ {6 6}}} \lambda_ {4 6} ^ {(4 4)} = \frac {8}{3 \sqrt {5 5}}, \quad \frac {1}{\sqrt {g _ {6 6}}} \lambda_ {6 6} ^ {6} = - \frac {2}{3} \sqrt {\frac {7}{1 1}}, \quad \frac {g _ {4 4}}{g _ {6 6}} \lambda_ {6 6} ^ {(4 4)} = \frac {8}{9} \tag {4.60}
$$

and

$$
\left(W _ {4}\right) _ {\mathfrak {S}} ^ {2, +} = \left(\beta_ {4} ^ {-}\right) ^ {2} - \left(\beta_ {4} ^ {+}\right) ^ {2} - \frac {2 \sqrt {g _ {4 4}}}{5 7 \sqrt {3 5}} \beta_ {2} ^ {2} \beta_ {4} ^ {+} + \operatorname {d e s c}. \tag {4.61}
$$

$$
\left(W _ {4}\right) _ {\mathfrak {S}} ^ {2, -} = 2 \beta_ {4} ^ {+} \beta_ {4} ^ {-} - \frac {2 \sqrt {g _ {4 4}}}{5 7 \sqrt {3 5}} \beta_ {2} ^ {2} \beta_ {4} ^ {-} + \operatorname {d e s c}. \tag {4.62}
$$

$$
\left(W _ {4}\right) _ {\mathfrak {S}} ^ {2, 0} = \left(\beta_ {4} ^ {+}\right) ^ {2} + \left(\beta_ {4} ^ {-}\right) ^ {2} + \frac {1}{1 0} \sqrt {\frac {1 1}{7}} \frac {g _ {4 4}}{\sqrt {g _ {6 6}}} \beta_ {2} \beta_ {6} - \frac {g _ {4 4}}{2 8 5 6 0} \beta_ {2} ^ {4} + \operatorname {d e s c}. \tag {4.63}
$$

Notice that we omitted the OPEs of  $W_{(3,0)}$  with the remaining generators as well as its explicit form in terms of free fields. The ... in the last three OPEs indicate additional contributions. For example in the OPEs (4.58c) and (4.58d) operators of type  $\mathfrak{L}_{(4,3)}$ ,  $\mathfrak{L}_{(4,2)}$  and  $\mathfrak{L}_{(4,1)}$  of the schematic form  $W_{4}W_{4}$  could appear. Let us also collect the quantum numbers of the null operators with small conformal weight as read off from the Hilbert series of  $\mathcal{R}_{D_4}$ :

$$
\mathfrak {L} _ {(4, 2)}, \quad \mathfrak {L} _ {(4, 1)}, \quad 2 \mathfrak {L} _ {(4, 0)}, \quad 2 \mathfrak {L} _ {(5, 3)}, \quad 2 \mathfrak {L} _ {(5, 2)}, \quad 2 \mathfrak {L} _ {(5, 1)}. \tag {4.64}
$$

The corresponding null operators have the schematic form  $W_4W_4$  and  $W_4W_6$ .

Remark: The OPEs above have an  $S_{3} \subset O(2)$  symmetry that acts non-trivially only on the generators  $W_{4}^{\pm}$  and is generated by the reflections  $s_{1}, s_{2}$  as

$$
\binom {\beta_ {4} ^ {+}} {\beta_ {4} ^ {-}} \mapsto s _ {k} \binom {\beta_ {4} ^ {+}} {\beta_ {4} ^ {-}}, \quad s _ {1} = \left( \begin{array}{l l} 1 & 0 \\ 0 & - 1 \end{array} \right), \quad s _ {2} = \frac {1}{2} \binom {- 1 \sqrt {3}} {\sqrt {3} 1}. \tag {4.65}
$$

Notice that the product  $s_1 s_2$  generates a  $\mathbb{Z}_3$  subgroup. In order to check this claim it is convenient to observe that  $(W_4)_{\mathfrak{S}}^{2,\pm}$  transform as  $\beta_4^{\pm}$  under (4.65). The relations (4.64) have definite transformation properties under  $S_3$ :

# 5 Examples of  $\mathcal{N} = 2$  VOA  $\mathcal{W}_{\mathrm{G}}$

# 5.1 Rank 1:  $G = \mathbb{Z}_p$

In this section we examine in detail the proposed free-field realization for the  $\mathcal{N} = 2$  VOA  $\mathcal{W}_{\mathbb{Z}_p}$  associated to the rank-1 complex reflection group  $\mathbb{Z}_p$ ,  $p \geq 2$ . These algebras have been first analyzed by means of direct bootstrap techniques in [30].

# 5.1.1 Construction of the generators and OPEs

According to our prescription, we need one copy of a  $\beta \gamma bc$  system. For the sake of simplicity, we omit the subscript 1 from the free fields in this section. The  $\mathcal{N} = 2$  SCA algebra is realized according to the formulae of section 3.1. In particular, the central charge is  $c = -3(2p - 1)$  and the conformal weight of  $\beta$  is  $p / 2$ .

Let  $W$ ,  $\overline{W}$  denote the chiral, antichiral extra generators of the VOA. The only non-trivial task at hand is the construction of the antichiral generator  $\overline{W}$ , of conformal weight  $h = p / 2$  and charge  $q = -p / 2$  using the strategy presented in section 3.3.

Before proceeding, let us stress that, for  $p = 2$ , supersymmetry enhances from  $\mathcal{N} = 2$  to small  $\mathcal{N} = 4$ , and the sought-for VOA is nothing but the small  $\mathcal{N} = 4$  SCA at central charge  $c = -9$ , with the identification  $\overline{W} = J^{-}$ . In this case the free-field realization coincides with the one of [26], which has been reviewed in section 4.1.

Let us now discuss the case of generic  $p$ . The direct analysis of a few examples reveals that there exists a unique, up to scaling, chiral  $\mathcal{N} = 2$  superVirasoro primary of conformal weight  $\frac{p}{2}$  that can be constructed using a single  $\beta \gamma bc$  of the type given in table (3.4). For definiteness, we fix the normalization of the generators  $W$  and  $\overline{W}$  as follows,

$$
W = \beta , \quad \bar {W} = \beta^ {p - 1} \gamma^ {p} + \dots . \tag {5.1}
$$

All omitted terms in  $\overline{W}$  contain at least one derivative or a pair of fermionic free fields.

In terms of the map  $\mathcal{P}$  of section 3.4,  $\mathcal{P}(\overline{W}) = \beta^{p - 1}\gamma^p$ . We may regard  $\overline{W}$  as the unique  $\mathcal{N} = 2$  super-Virasoro primary completion of the monomial  $\beta^{p - 1}\gamma^p$ .

Once the form of the antichiral generator  $\overline{W}$  is fixed, we may check that the OPE  $\overline{W}\overline{W}$  is regular, and that the singular part of the  $W\overline{W}$  OPE is expressed entirely in terms of  $\mathcal{N} = 2$  super-Virasoro descendants of the identity, as it must be. In particular, we can verify that

$$
W \times \bar {W} = g _ {W \bar {W}} [ \mathrm {i d} ], \quad g _ {W \bar {W}} = (-) ^ {p} \frac {(2 p - 1) !}{p ! p ^ {p - 1}}, \tag {5.2}
$$

where [id] denotes the contribution for the  $\mathcal{N} = 2$  super-Virasoro family of the identity operator. More explicitly

$$
[ \mathrm {i d} ] = \mathrm {i d} + \frac {3 p}{c} \mathcal {J} + \frac {p}{c} \mathcal {T} + \frac {3 p (3 p - 1)}{2 c (c - 1)} \left(\left(\mathcal {J} \mathcal {J}\right) _ {0} - \frac {2}{3} \mathcal {T}\right) + \dots \tag {5.3}
$$

Similarly to the  $\mathcal{N} = 4$  case given in (4.23), the contribution of the stress tensor  $\mathcal{T}$  to this OPE is split in two parts: the first correspond to the  $\mathfrak{osp}(2|2)$  descendant of  $\mathcal{I}$ , the second to the completion of  $(\mathcal{J}\mathcal{J})_0$  to a  $\mathfrak{osp}(2|2)$  primary whose norm is given by  $\frac{2}{9} c(c - 1)$ . The relation (5.2) has been checked explicitly for  $p = 2,3,\ldots,7$ . For  $p = 2,3,4,5$  we can record the entire content of the singular part of the  $W\overline{W}$  OPE in terms of quasiprimary fields

$$
p = 2 : \quad [ \mathrm {i d} ] = \mathrm {i d} - \frac {2}{3} \mathcal {J} ,
$$

$$
p = 3: \quad [ \mathrm {i d} ] = \mathrm {i d} - \frac {3}{5} \mathcal {J} - \frac {3}{1 0} \mathcal {T} + \frac {3}{2 0} (\mathcal {J} \mathcal {J}) _ {0},
$$

$$
p = 4: \quad [ \mathrm {i d} ] = \mathrm {i d} - \frac {4}{7} \mathcal {J} - \frac {2}{7} \mathcal {T} + \frac {1}{7} (\mathcal {J} \mathcal {J}) _ {0} - \frac {2}{1 0 5} (\mathcal {J} (\mathcal {J} \mathcal {J}) _ {0}) _ {0} + \frac {2}{3 5} (\mathcal {G} \bar {\mathcal {G}}) _ {0} + \frac {6}{3 5} (\mathcal {J} \mathcal {T}) _ {0},
$$

$$
\begin{array}{l} p = 5: \quad [ \mathrm {i d} ] = \mathrm {i d} - \frac {5}{9} \mathcal {J} - \frac {5}{1 8} \mathcal {T} + \frac {5}{3 6} (\mathcal {J} \mathcal {J}) _ {0} - \frac {5}{2 5 2} (\mathcal {J} (\mathcal {J} \mathcal {J}) _ {0}) _ {0} \\ + \frac {5}{1 2 6} (\mathcal {G} \widetilde {\mathcal {G}}) _ {0} + \frac {1 0}{6 3} (\mathcal {J T}) _ {0} + \frac {2 5}{5 0 4} (\mathcal {T T}) _ {0} + \frac {5}{3 0 2 4} (\mathcal {J} (\mathcal {J} (\mathcal {J T}) _ {0}) _ {0}) _ {0} \\ - \frac {5}{1 2 6} (\mathcal {J} (\mathcal {J T}) _ {0}) _ {0} - \frac {5}{2 5 2} (\mathcal {J} (\mathcal {G} \widetilde {\mathcal {G}}) _ {0}) _ {0} - \frac {5}{2 5 2} (\mathcal {G} \widetilde {\mathcal {G}}) _ {- 1} + \frac {9 5}{3 0 2 4} (\mathcal {J J}) _ {- 2}. \tag {5.4} \\ \end{array}
$$

We are using a compact notation in which all  $z_{12}$  factors and  $z$  derivatives are omitted, since they can be straightforwardly recovered from  $\mathfrak{sl}(2)$  covariance. For further details on the notation, the reader is referred to appendix A. We are only recording quasiprimary operators

that enter the singular part of the OPEs. Clearly, all OPE coefficients on the RHS of the OPE  $W\overline{W}$  can in principle be recovered from the two-point function coefficient  $g_{W\overline{W}}$  by exploiting the  $\mathcal{N} = 2$  superVirasoro symmetry.

# 5.1.2 Null states

An essential feature of our proposed free-field realization is that all null states of the VOA are realized manifestly as zero. First of all, let us verify this claim for the null state corresponding to the "Higgs branch" relation of the associated variety  $\mathcal{M}_{\mathbb{Z}_p} = (\mathbb{C}\oplus \mathbb{C}^*) / \mathbb{Z}_p$ . To describe this variety we introduce a complex coordinate  $z$  and its conjugate  $\bar{z}$ , and let the generator of  $\mathbb{Z}_p$  act on  $z$ ,  $\bar{z}$  as

$$
z \mapsto e ^ {2 \pi i / p} z, \quad \bar {z} \mapsto e ^ {- 2 \pi i / p} \bar {z}. \tag {5.5}
$$

The invariants are clearly

$$
w = z ^ {p}, \quad \bar {w} = \bar {z} ^ {p}, \quad j = z \bar {z}, \tag {5.6}
$$

satisfying one relation,

$$
w \bar {w} = j ^ {p}. \tag {5.7}
$$

At the level of the VOA, this invariant motivates us to consider the  $\mathcal{N} = 2$  superVirasoro primary composite operator

$$
\mathfrak {X} _ {p, 0} ^ {W \bar {W}} = (W \bar {W}) _ {0} + \dots . \tag {5.8}
$$

The dots represent the terms needed to obtain a superVirasoro primary, and are uniquely determined by the OPEs of the abstract VOA. For example, the full expressions of this composite for  $p = 2,3,4$  are

$$
p = 2: \mathfrak {X} _ {2, 0} ^ {W \overline {{W}}} = (W \overline {{W}}) _ {0} - \frac {1}{4} (\mathcal {J} \mathcal {J}) _ {0} + \frac {1}{2} T,
$$

$$
p = 3: \mathfrak {X} _ {3, 0} ^ {W \overline {{W}}} = (W \overline {{W}}) _ {0} - \frac {1}{2 \overline {{7}}} (\mathcal {J} (\mathcal {J} \mathcal {J}) _ {0}) _ {0} + \frac {2}{9} (\mathcal {G} \widetilde {\mathcal {G}}) _ {0} + \frac {4}{9} (\mathcal {J} T) _ {0},
$$

$$
\begin{array}{l} p = 4: \quad \mathfrak {X} _ {4, 0} ^ {W \overline {{W}}} = (W \overline {{W}}) _ {0} - \frac {1}{2 5 6} (\mathcal {J} (\mathcal {J} (\mathcal {J}) _ {0}) _ {0}) _ {0} - \frac {3}{1 6} (T T) _ {0} + \frac {9}{6 4} (\mathcal {J} (\mathcal {J} T) _ {0}) _ {0} \\ + \frac {3}{3 2} (\mathcal {J} (\mathcal {G} \widetilde {\mathcal {G}}) _ {0}) _ {0} + \frac {3}{3 2} (\mathcal {G} \widetilde {\mathcal {G}}) _ {- 1} - \frac {1 5}{1 2 8} (\mathcal {J} \mathcal {J}) _ {- 2}. \tag {5.9} \\ \end{array}
$$

We checked that, in our free-field realization, these composites are indeed identically zero. Even though the full expression of the operator  $\mathfrak{X}_{p,0}^{W\overline{W}}$  in the VOA becomes increasingly complex as we increase  $p$ , the classical counterpart of  $\mathfrak{X}_{p,0}^{W\overline{W}}$  via the map  $\mathcal{P}$  of section 3.4 has a very simple structure. Indeed, one verifies that

$$
W _ {\mathrm {c l}} := \mathcal {P} (W) = \beta , \quad \overline {{W}} _ {\mathrm {c l}} := \mathcal {P} (\overline {{W}}) = \beta^ {p - 1} \gamma^ {p}, \quad \mathcal {J} _ {\mathrm {c l}} := \mathcal {P} (\mathcal {J}) = p \beta \gamma ,
$$

$$
\left(\mathfrak {X} _ {p, 0} ^ {W \overline {{W}}}\right) _ {\mathrm {c l}} := \mathcal {P} \left(\mathfrak {X} _ {p, 0} ^ {W \overline {{W}}}\right) = W _ {\mathrm {c l}} \overline {{W}} _ {\mathrm {c l}} - \frac {1}{p ^ {p}} \left(\mathcal {J} _ {\mathrm {c l}}\right) ^ {p}. \tag {5.10}
$$

These expressions show that the operator  $\mathfrak{X}_{p,0}^{WW}$  is indeed the null operator associated to the "Higgs branch" relation (5.7). It is also straightforward to find the map between the classical Poisson variables  $z$ ,  $\bar{z}$  and  $\beta$ ,  $\gamma$ ,

$$
\beta = z ^ {p}, \quad \gamma = z ^ {1 - p} \bar {z}. \tag {5.11}
$$

Let us now discuss another pair of null states that are expected in the VOA  $\mathcal{W}_{\mathbb{Z}_p}$ . They are a pair of non-chiral,  $\mathcal{N} = 2$  superVirasoro primary operators linear in  $W_{p}$  and  $\overline{W}_p$  respectively, given by

$$
\mathfrak {X} _ {\frac {p}{2} + \frac {3}{2}, - 1} ^ {\mathcal {G} W} = (\mathcal {G} W) _ {0} - \frac {1}{p} \left(\mathcal {J} \mathcal {G} _ {W}\right) _ {0}, \quad \mathfrak {X} _ {\frac {p}{2} + \frac {3}{2}, + 1} ^ {\widetilde {\mathcal {G}} \overline {{W}}} = \left(\widetilde {\mathcal {G}} \overline {{W}}\right) _ {0} + \frac {1}{p} \left(\mathcal {J} \widetilde {\mathcal {G}} _ {\overline {{W}}}\right) _ {0}. \tag {5.12}
$$

The operator  $\mathcal{G}_W$  is the supersymmetry descendant of  $W_{p}$ , defined as  $\{\mathcal{G}W\}_{1}$ . A similar remark applies to  $\widetilde{\mathcal{G}}_{\overline{W}}$ . Let us stress that, in order to verify that  $\mathfrak{X}_{\frac{p}{2} +\frac{3}{2}, - 1}^{\mathcal{G}W}$ ,  $\mathfrak{X}_{\frac{p}{2} +\frac{3}{2}, + 1}^{\widetilde{\mathcal{G}}\overline{W}}$  are superVirasoro primaries, we only need to use the OPEs of the abstract VOA, and not our specific free-field realization. As already argued in [30], these composite operators must be null in order for the VOA to exist. In our free-field construction one can indeed verify that both  $\mathfrak{X}_{\frac{p}{2} +\frac{3}{2}, - 1}^{\mathcal{G}W}$  and  $\mathfrak{X}_{\frac{p}{2} +\frac{3}{2}, + 1}^{\widetilde{\mathcal{G}}\overline{W}}$  vanish identically.

Finally, let us comment on the relation between the "Higgs branch" null state  $\mathfrak{X}_{p,0}^{W\overline{W}}$  and the fermionic null states  $\mathfrak{X}_{\frac{p}{2} +\frac{3}{2}, - 1}^{\mathcal{G}W}$ ,  $\mathfrak{X}_{\frac{p}{2} +\frac{3}{2}, + 1}^{\widetilde{\mathcal{G}}\overline{W}}$ . The latter can be obtained from  $\mathfrak{X}_{p,0}^{W\overline{W}}$  by taking a singular OPE with  $\mathcal{G}_W$ ,  $\widetilde{\mathcal{G}}_{\overline{W}}$ . For example,

$$
p = 3: \qquad \{\mathcal {G} _ {W}   \mathfrak {X} _ {3, 0} ^ {W \overline {{W}}} \} _ {2} = \frac {2}{3}   \mathfrak {X} _ {3, - 1} ^ {\mathcal {G} W}, \qquad \{\widetilde {\mathcal {G}} _ {\overline {{W}}}   \mathfrak {X} _ {3, 0} ^ {W \overline {{W}}} \} _ {2} = \frac {4}{3}   \mathfrak {X} _ {3, + 1} ^ {\widetilde {\mathcal {G}} \overline {{W}}},
$$

$$
p = 4: \quad \left\{\mathcal {G} _ {W} \mathfrak {X} _ {4, 0} ^ {W \overline {{W}}} \right\} _ {3} = \frac {3}{8} \mathfrak {X} _ {\frac {7}{2}, - 1} ^ {\mathcal {G} W}, \quad \left\{\widetilde {\mathcal {G}} _ {\overline {{W}}} \mathfrak {X} _ {4, 0} ^ {W \overline {{W}}} \right\} _ {3} = \frac {2 7}{8} \mathfrak {X} _ {\frac {7}{2}, + 1} ^ {\widetilde {\mathcal {G}} \overline {{W}}}. \tag {5.13}
$$

For general  $p$ ,  $\mathfrak{X}_{\frac{p}{2} + \frac{3}{2}, -1}^{\mathcal{G}W}$  enters the order  $(p - 1)$  pole of the  $\mathcal{G}_W \mathfrak{X}_{p,0}^{W\overline{W}}$  OPE, and similarly  $\mathfrak{X}_{\frac{p}{2} + \frac{3}{2}, -1}^{\widetilde{\mathcal{G}}\overline{W}}$  enters the order  $(p - 1)$  pole of the  $\widetilde{\mathcal{G}}_{\overline{W}} \mathfrak{X}_{p,0}^{W\overline{W}}$  OPE. In other words, the operators  $\mathfrak{X}_{\frac{p}{2} + \frac{3}{2}, -1}^{\mathcal{G}W}$ ,  $\mathfrak{X}_{\frac{p}{2} + \frac{3}{2}, -1}^{\widetilde{\mathcal{G}}\overline{W}}$  belong to the ideal of  $\mathcal{W}_{\mathbb{Z}_p}$  generated by  $\mathfrak{X}_{p,0}^{W\overline{W}}$ . This observation is consistent with the expectation that the "Higgs branch" null state  $\mathfrak{X}_{p,0}^{W\overline{W}}$  generates all null states in  $\mathcal{W}_{\mathbb{Z}_p}$ .

# 5.1.3 Screening operator

We have a conjectural characterization of the VOA  $\mathcal{W}_{\mathbb{Z}_p}$  as a subalgebra of the  $\beta \gamma bc$  system in terms of the kernel of a screening operator  $\mathbb{S}$ . Our proposal is a natural generalization of the screening operator discussed in [26] in the case  $p = 2$ , i.e. the small  $\mathcal{N} = 4$  SCA at  $c = -9$ .

In order to define  $\mathbb{S}$ , we first have to express  $\beta$  and  $\gamma$  in terms of chiral bosons  $\chi$ ,  $\phi$ . The OPEs of the chiral bosons are recorded in (4.6), and the expressions of  $\beta$ ,  $\gamma$  in terms of the chiral bosons are given in (4.5). We may now define the screening current

$$
\mathbf {J} = b e ^ {(p ^ {- 1} - 1) (\chi + \phi)}. \tag {5.14}
$$

This operator has conformal dimension 1 and charge 0 under  $\mathcal{I}$ . We conjecture that the VOA  $\mathcal{W}_{\mathbb{Z}_p}$  coincides with the kernel of  $\mathbb{S} = \int \mathsf{J}$  acting on the  $\beta \gamma bc$  system. Its action on any object

$X$  is defined as

$$
\mathbb {S} \cdot X = \{\mathrm {J} X \} _ {1}. \tag {5.15}
$$

In the case  $p = 2$ , the conjecture is proven in [26]. It is worth pointing out that the operator  $\mathbf{J}$  can be expressed as a supersymmetry descendant of an operator  $\mathsf{K}$  with dimension and charge  $1/2$ ,

$$
J = \left\{\mathcal {G} K \right\} _ {1}, \quad K = e ^ {p ^ {- 1} (\chi + \phi)}. \tag {5.16}
$$

The operator  $\mathsf{K}$  is a chiral  $\mathcal{N} = 2$  super-Virasoro primary.

It is a matter of straightforward computation to verify that all generators of the VOA  $\mathcal{W}_{\mathbb{Z}_p}$  are annihilated by the action of  $\mathbb{S}$ . It is then immediate that  $\mathcal{W}_{\mathbb{Z}_p}$  is contained in the kernel of  $\mathbb{S}$ . We do not have a general argument for the reversed inclusion for  $p \geq 3$ , but we have checked our claim in a few examples, for some states with low twist  $h - m$ . More precisely, we worked at generic  $p$  and considered operators of twist up to 2. A direct computation of the kernel of  $\mathbb{S}$  in the  $\beta \gamma bc$  system shows a perfect agreement with a counting of states in  $\mathcal{W}_{\mathbb{Z}_p}$ . In order for this match to work, it is essential that in our free-field construction all null states are (conjecturally) identically zero.

Let us close this section by analyzing a few properties of the classical counterpart of the screening current  $\mathbf{J}$  of (5.14). The map  $\mathcal{P}'$  of section 3.4 maps the VOA  $\mathcal{W}_{\mathbb{Z}_p}$  to the Poisson superalgebra of functions of the variables  $\beta, \gamma, b, c$ . The classical counterpart of (5.14) is simply

$$
J _ {\mathrm {c l}} = b \beta^ {\frac {1}{p} - 1}. \tag {5.17}
$$

Its action on the classical variables  $\beta, \gamma, b, c$  is the following,

$$
\{\mathrm {J} _ {\mathrm {c l}}, b \} _ {\mathrm {P B}} = 0, \quad \{\mathrm {J} _ {\mathrm {c l}}, \beta \} _ {\mathrm {P B}} = 0,
$$

$$
\left\{\mathrm {J} _ {\mathrm {c l}}, c \right\} _ {\mathrm {P B}} = \beta^ {\frac {1}{p} - 1}, \quad \left\{\mathrm {J} _ {\mathrm {c l}}, \gamma \right\} _ {\mathrm {P B}} = \left(1 - \frac {1}{p}\right) b \beta^ {\frac {1}{p} - 2}. \tag {5.18}
$$

# 5.2 A rank 2 example:  $G = G(3,1,2)$

In this section we discuss the proposed free-field realization for the  $\mathcal{N} = 2$  VOA associated to the rank-2 complex reflection group  $G(3,1,2)$ .

The ring  $\mathcal{R}_{\mathsf{G}}$  for  $\mathsf{G} = G(3,1,2)$ . To begin with let us describe the family of complex reflection groups  $G(k,1,\mathfrak{r})$ . Recall that  $\mathfrak{r}$  is the rank of  $G(k,1,\mathfrak{r})$  and the invariants have degrees

$$
k, 2 k, \dots , r k. \tag {5.19}
$$

The action of  $G(k,1,\mathfrak{r})$  on  $(z_{1},\ldots ,z_{\mathfrak{r}})\in V_{\mathsf{G}}\simeq \mathbb{C}^{\mathfrak{r}}$  is generated by permutations  $\sigma_{i} = p_{i,i + 1}$  (with an obvious action on the coordinates) together with

$$
\tau : \left(z _ {1}, \dots , z _ {\mathrm {r}}\right) \mapsto \left(\omega z _ {1}, z _ {2} \dots , z _ {\mathrm {r}}\right), \quad \omega^ {k} = 1. \tag {5.20}
$$

The ring of invariants is freely generated by the elementary symmetric polynomials in  $z_1^k, \ldots, z_r^k$ . Notice that for  $k = 2$  this group is the Coxeter group  $B_r$ .

In the following we consider in more details the example of  $G(3,1,2)$ . In this case the ring  $\mathcal{R}_{\mathsf{G}}$  is generated by

$$
w _ {3} = z _ {1} ^ {3} + z _ {2} ^ {3}, \qquad w _ {6} = z _ {1} ^ {6} + z _ {2} ^ {6}, \qquad O = z _ {1} ^ {4} \bar {z} _ {1} + z _ {2} ^ {4} \bar {z} _ {2},
$$

$$
\bar {w} _ {3} = \bar {z} _ {1} ^ {3} + \bar {z} _ {2} ^ {3}, \quad \bar {w} _ {6} = \bar {z} _ {1} ^ {6} + \bar {z} _ {2} ^ {6}, \quad \bar {O} = z _ {1} \bar {z} _ {1} ^ {4} + z _ {2} \bar {z} _ {2} ^ {4}, \tag {5.21}
$$

$$
j = z _ {1} \bar {z} _ {1} + z _ {2} \bar {z} _ {2}, U = z _ {1} ^ {2} \bar {z} _ {1} ^ {2} + z _ {2} ^ {2} \bar {z} _ {2} ^ {2}.
$$

The relations of lowest conformal weight, namely  $h = 4$  and  $h = 4 + \frac{1}{2}$ , take the explicit form

$$
U ^ {2} - \bar {O} w _ {3} - O \bar {w} _ {3} + j w _ {3} \bar {w} _ {3} + + \frac {1}{2} j ^ {2} U - \frac {1}{2} j ^ {4} = 0, \tag {5.22}
$$

$$
\left(w _ {6} + w _ {3} w _ {3}\right) \bar {w} _ {3} - O U - j ^ {2} O + 2 j U w _ {3} = 0. \tag {5.23}
$$

The second equation is accompanied by it conjugate.

Free-field realization. As explained in section 3 the free fields in this case are two copies of the  $\beta \gamma bc$  system, with weights determined by the degree of the invariants of  $\mathbb{C}^2 / G(3,1,2)$ ,

$$
(p _ {1}, p _ {2}) = (3, 6). \tag {5.24}
$$

The small  $\mathcal{N} = 2$  algebra with central charge  $c = -48$  is realized according to the formulae in 3.1, and the chiral generators  $\mathcal{W}_3$ ,  $\mathcal{W}_6$  are simply realized as

$$
\mathcal {W} _ {3} = \beta_ {1}, \quad \mathcal {W} _ {6} = \beta_ {2}. \tag {5.25}
$$

The remaining generators $^{22}$  with their quantum numbers are summarized in the following table

<table><tr><td></td><td>W3</td><td>W6</td><td>J</td><td>O</td><td>U</td><td>W3</td><td>O</td><td>W6</td></tr><tr><td>h</td><td>3/2</td><td>3</td><td>1</td><td>5/2</td><td>2</td><td>3/2</td><td>5/2</td><td>3</td></tr><tr><td>m</td><td>3/2</td><td>3</td><td>0</td><td>3/2</td><td>0</td><td>-3/2</td><td>-3/2</td><td>-3</td></tr><tr><td>h-m</td><td>0</td><td>0</td><td>1</td><td>1</td><td>2</td><td>3</td><td>4</td><td>6</td></tr><tr><td>#</td><td>-</td><td>-</td><td>-</td><td>4</td><td>8</td><td>13</td><td>36</td><td>104</td></tr></table>

According to the prescription given in section 3.3, the first step is to construct all  $\mathcal{N} = 2$  superVirasoro primary with the appropriate weights. The entry  $\sharp$  denotes the number of such objects. The symbol - indicates that the corresponding entry has already been constructed so we do not need an Ansatz for it. After imposing that the algebra closes one finds the

complete list of OPEs to be

$$
\mathcal {W} _ {3} \times \overline {{\mathcal {W}}} _ {3} \sim [ \mathrm {i d} ] + [ \mathcal {U} ], \tag {5.27a}
$$

$$
\mathcal {W} _ {6} \times \overline {{\mathcal {W}}} _ {6} \sim [ \mathrm {i d} ] + \frac {7}{1 4 3} [ \mathcal {U} ] - \frac {2 6 8}{4 2 9} [ \mathcal {W} _ {3} \overline {{\mathcal {W}}} _ {3} ] + \dots , \tag {5.27b}
$$

$$
\mathcal {O} \times \overline {{\mathcal {O}}} \sim \frac {1 3 3 6 5}{5 9 5 8 4} \left([ \mathrm {i d} ] - \frac {2 9 5}{6 2 7} [ \mathcal {U} ] - \frac {3 4 7}{2 8 5} [ \mathcal {W} _ {3} \overline {{\mathcal {W}}} _ {3} ]\right), \qquad \mathcal {U} \times \mathcal {U} \sim \frac {9 9}{3 9 2} [ \mathrm {i d} ] - \frac {4 5}{9 8} [ \mathcal {U} ], \qquad (5. 2 7 \mathrm {c})
$$

$$
\mathcal {W} _ {3} \times \mathcal {U} \sim \frac {9 9}{3 9 2} [ \mathcal {W} _ {3} ] + [ \mathcal {O} ], \quad \mathcal {W} _ {3} \times \overline {{\mathcal {O}}} \sim - \frac {1 3 5}{1 5 2} [ \mathcal {U} ] + \frac {2 7}{7 6} [ \mathcal {W} _ {3} \overline {{\mathcal {W}}} _ {3} ], \tag {5.27d}
$$

$$
\mathcal {W} _ {3} \times \mathcal {O} \sim \frac {2 7}{2 1 2 8} \sqrt {\frac {2 1 4 5}{2}} [ \mathcal {W} _ {6} ], \quad \mathcal {W} _ {6} \times \mathcal {O} \sim \frac {2 7}{1 0 6 4} [ \mathcal {W} _ {3} \mathcal {W} _ {6} ] + \frac {8 1}{2 6 6} \sqrt {\frac {1 5}{2 8 6}} [ \mathcal {W} _ {3} \mathcal {W} _ {3} \mathcal {W} _ {3} ], \quad (5. 2 7 \mathrm {e})
$$

$$
\mathcal {W} _ {6} \times \overline {{\mathcal {O}}} \sim - \frac {2 7}{2 1 2 8} \sqrt {\frac {2 1 4 5}{2}} \left[ \mathcal {W} _ {3} \right] + \frac {4 3}{1 3 3} \sqrt {\frac {5}{8 5 8}} \left[ \mathcal {O} \right] + - \frac {4 3 5}{2 6 6} \sqrt {\frac {1 5}{2 8 6}} \left[ \mathcal {W} _ {3} \mathcal {U} \right] + \dots , \tag {5.27f}
$$

$$
\mathcal {W} _ {3} \times \overline {{\mathcal {W}}} _ {6} \sim \frac {3}{1 4} \sqrt {\frac {1 6 5}{2 6}} [ \overline {{\mathcal {W}}} _ {3} ] - \frac {1 4}{3} \sqrt {\frac {2 6}{1 6 5}} [ \overline {{\mathcal {O}}} ] - 6 \sqrt {\frac {6}{7 1 5}} [ \mathcal {U} \overline {{\mathcal {W}}} _ {3} ], \tag {5.27g}
$$

$$
\mathcal {W} _ {6} \times \mathcal {U} \sim - \frac {1 1 7}{9 3 1} \left[ \mathcal {W} _ {6} \right] + \frac {2 7}{2 6 6} \sqrt {\frac {1 6 5}{2 6}} \left[ \mathcal {W} _ {3} \mathcal {W} _ {3} \right] + 4 \sqrt {\frac {1 0}{4 2 9}} \left[ \mathcal {W} _ {3} \mathcal {O} \right], \tag {5.27h}
$$

$$
\mathcal {O} \times \mathcal {O} \sim \frac {4 0 5}{1 6 1 7 2 8} \sqrt {\frac {2 1 4 5}{2}} \left[ \mathcal {W} _ {6} \right] - \frac {4 0 0 9 5}{2 8 3 0 2 4} \left[ \mathcal {W} _ {3} \mathcal {W} _ {3} \right], \tag {5.27i}
$$

$$
\mathcal {U} \times \mathcal {O} \sim \frac {1 3 3 6 5}{5 9 5 8 4} [ \mathcal {W} _ {3} ] - \frac {8 8 5}{7 4 4 8} [ \mathcal {O} ] - \frac {2 7}{5 6} [ \mathcal {U W} _ {3} ]. \tag {5.27j}
$$

while the  $\mathcal{W}$ - $\mathcal{W}$ ,  $\overline{\mathcal{W}}$ - $\overline{\mathcal{W}}$  OPEs are regular. The ... above indicate contributions from super-Virasoro primaries with higher conformal dimensions. They have the same quantum numbers as the relations (5.22). It is easy to check that the expression of the corresponding operators in the free-field realization is zero. Notice that the first four OPEs are real, the remaining nine are complex so they are accompanied by their conjugate.

A few additional remarks are in order. First of all, in order to constrain the form of the generators in terms of free fields, it is convenient to start by imposing all the linear constraints originating form having the right structure of the OPE of the chiral generators  $\mathcal{W}_3, \mathcal{W}_6$  with everything else. Second of all, observe that once  $\overline{\mathcal{W}}_3$  is constructed, all the remaining generators, namely  $\mathcal{U}, \mathcal{O}, \overline{\mathcal{O}}$  and  $\overline{\mathcal{W}}_6$  are generated $^{23}$  in the OPE. Finally, let us point out that the procedure given in section 3.3 lives the freedom of redefining  $\overline{\mathcal{W}}_6 \to x(\alpha) \overline{\mathcal{W}}_6 + \alpha \overline{\mathcal{W}}_3 \overline{\mathcal{W}}_3$ . We fix this freedom by imposing that the OPEs are manifestly  $\mathbb{Z}_2$  conjugation symmetric.

# 6 The  $R$ -filtration from free fields

This work was largely motivated by the desire to achieve a better understanding of the VOAs that arise from four-dimensional  $\mathcal{N} = 4$  superconformal field theories via the map introduced in [1]. Applying the construction of [1] to the  $\mathcal{N} = 4$  super Yang-Mills theories, one finds a class of  $\mathcal{N} = 4$  VOAs labelled by Weyl groups. In this paper, we have offered an alternative construction of  $\mathcal{N} = 4$  VOAs labelled by Weyl groups, and conjectured that they are in fact the same as the ones that arise from SYM theories.[24]

So far, we have studied these algebras in their own right, as novel and interesting examples of VOAs. In this section we wish to go back to their four-dimensional interpretation. We propose that the free-field representations that we have introduced allow to solve a longstanding open problem, the assignment of the “ $R$ -filtration” [13]. Let us briefly review the terms of the problem, for the general case of VOAs that descend from arbitrary  $\mathcal{N} = 2$  4d SCFTs.

According to the map of [1], the operators of the VOA are in one-to-one correspondence with the so-called "Schur operators" of the parent  $\mathcal{N} = 2$  SCFT. Schur operators belong to certain (semi)short representations of the  $\mathfrak{sl}(4|2)$  superconformal algebra - as such, due to the shortening conditions, they are labelled by three out of the five Cartan quantum numbers[25] of  $\mathfrak{sl}(4|2)$ . Of these quantum numbers, all except one survive as good quantum number of the VOA, defining gradings respected by the operator product expansion. The exception is the  $\mathfrak{sl}(2)_R$  symmetry Cartan, denoted by  $R$ . While Schur operators are all highest weights of  $\mathfrak{sl}(2)_R$ , the cohomological construction of [1] involves the lower components of the  $\mathfrak{sl}(2)_R$  multiplet to which a Schur operator belongs, and as a result  $R$  does not descend to a grading of the VOA. It is however easy to argue [13] that  $R$  can only decrease or remain constant in the OPE, and as such it defines a filtration of the VOA. It follows that any VOA associated to a  $4d$  SCFT by the map of [1] must admit such an  $R$ -filtration. However, the  $R$ -filtration does not appear to be intrinsic to the VOA - at least, not in any obvious way. Given an abstract presentation of the VOA (say in terms of strong generators and their singular OPEs) it is a priori unclear how to determine its  $R$ -filtration. This is a severe limitation if one's goal is to use the VOA as a tool to study the parent  $4d$  theory, because without knowledge of the  $R$  quantum number the identification of  $4d$  Schur representations is ambiguous.[26]

Our main new observation is that the free-field constructions analyzed in this paper come equipped with a natural filtration, which we call the “ $\mathcal{R}$ -filtration”. We conjecture that for the  $\mathcal{N} = 4$  VOAs that arise from  $\mathcal{N} = 4$  SYM theories, the  $\mathcal{R}$ - and  $R$ -filtrations coincide. We have performed several successful checks of this conjecture. We expect this statement to generalize to the VOAs labelled by complex reflection groups that descend from  $4d\mathcal{N} = 3$  SCFTs, but we did not perform any check of that more general conjecture.

We begin in the next subsection with a review of the  $R$ -filtration for VOAs that arise from general  $\mathcal{N} = 2$  SCFTs, and indicate its obvious extension to the  $\mathcal{N} = 4$  case. In section 6.2 we show that the VOAs studied in this paper admit a natural “ $\mathcal{R}$ -filtration”, defined in terms of their free-field realizations. According to the basic dictionary of [1], the vacuum character of the VOA reproduces the Schur limit of the superconformal index of its parent  $4d$  SCFT, which is a function of a single superconformal fugacity  $q$ . Knowledge of the  $R$ -filtration allows to refine the vacuum character, so that it yields the full Macdonald index of the  $4d$  SCFT, a function of two superconformal fugacities  $q$  and  $t$ . The basic check that we have correctly

identified the  $R$ -filtration consists in matching the refined vacuum character computed from our free-field realization with the well-known Madconald index of an  $\mathcal{N} = 4$  SYM theory. In section 6.3, we perform this check in several examples, up to some order in an expansion in the conformal weight. In the limit  $q \to 0$ , the Macdonald index reduces to the so-called Hall-Littlewood index, which is much simpler to compute. In section 6.4 we collect some further evidence for our proposal, showing that it correctly reproduces the full Hall-Littlewood index in rank-one and rank-two examples.

# 6.1 The  $R$ -filtration

For the reader's convenience, we start by reviewing the salient facts of the  $4d / 2d$  correspondence introduced in [1].

# 6.1.1 The  $4d / 2d$  map and Schur operators

Given an  $\mathcal{N} = 2$  SCFT, the associated VOA is obtained by passing to the cohomology of a certain nilpotent fermionic generator  $\mathfrak{Q}$  of the  $\mathcal{N} = 2$  superconformal algebra  $\mathfrak{sl}(4|2)$ . We denote the Cartan quantum numbers of  $\mathfrak{sl}(4|2)$  by  $(E,j_1,j_2,R,r)$ , where  $E$  is the conformal dimension and  $j_{1}, j_{2}, R, r$  are eigenvalues with respect to the Cartan generators of  $\mathfrak{sl}(2)_1 \oplus \mathfrak{sl}(2)_2 \oplus \mathfrak{sl}(2)_R \oplus \mathfrak{gl}(1)_r$ , respectively. The nontrivial cohomology classes of local operators inserted at the origin in  $\mathbb{R}^4 \cong \mathbb{C}^2$  ( $z = \bar{z} = w = \bar{w} = 0$ ) have canonical representatives which are the Schur operators [45]. These are local operators whose quantum numbers satisfy the relations<sup>27</sup>

$$
E - \left(j _ {1} + j _ {2}\right) - 2 R = 0, \tag {6.1}
$$

$$
r + j _ {1} - j _ {2} = 0.
$$

Schur operators are always the highest weight states of their respective  $\mathfrak{sl}(2)_1 \oplus \mathfrak{sl}(2)_2 \oplus \mathfrak{sl}(2)_R$  modules. The various (unitary) supermultiplets that contain Schur operators and the positioning of Schur operators within those multiplets is summarized in table 4.

Finite linear combinations of local operators inserted away from the origin cannot define nontrivial  $\mathbb{Q}$ -cohomology classes unless  $w = \bar{w} = 0$ . A canonical choice of representatives for local operators inserted on the  $w = \bar{w} = 0$  plane,  $\mathbb{C}_{[z,\bar{z} ]}$ , is given by twisted translated Schur operators,

$$
\mathcal {O} (z) \equiv e ^ {z L _ {- 1} + \bar {z} (\bar {L} _ {- 1} + R ^ {-})} \mathcal {O} _ {\mathrm {S c h}} (0) e ^ {- z L _ {- 1} - \bar {z} (\bar {L} _ {- 1} + R ^ {-})}, \tag {6.2}
$$

where  $L_{-1}$  and  $\bar{L}_{-1}$  are the generators of holomorphic and antiholomorphic translations in  $\mathbb{C}_{[z,\bar{z} ]}$ ,  $R^{-}$  is the lowering operator of  $\mathfrak{su}(2)_R$ , and  $\mathcal{O}_{\mathrm{Sch}}(z,\bar{z})$  is a Schur operator. The OPE of twisted-translated Schur operators, taken at the level of  $\mathbb{Q}$ -cohomology, is  $\mathfrak{sl}(2)_z$  covariant and  $\mathfrak{sl}(2)_{\bar{z}}$  invariant, with the holomorphic dimension of the twisted-translated operator  $\mathcal{O}(z)$  being determined in terms of the quantum numbers of the corresponding Schur operator according to

$$
\left[ L _ {0}, \mathcal {O} (z) \right] = h \mathcal {O} (z), \quad h = \frac {E + j _ {1} + j _ {2}}{2} = E - R. \tag {6.3}
$$

<table><tr><td>Multiplet</td><td>OSchur</td><td>h</td><td>r</td></tr><tr><td>B_R</td><td>Ψ11...1</td><td>R</td><td>0</td></tr><tr><td>DR(0,j2)</td><td>Q1+ Ψ11...1+...+</td><td>R+j2+1</td><td>j2+1/2</td></tr><tr><td>DR(j1,0)</td><td>Q1+ Ψ11...1+...+</td><td>R+j1+1</td><td>-j1-1/2</td></tr><tr><td>C_R(j1,j2)</td><td>Q1+ Q1+ Ψ11...1+...+</td><td>R+j1+j2+2</td><td>j2-j1</td></tr></table>

Table 4: Summary of the appearance of Schur operators in short multiplets of the  $\mathcal{N} = 2$  superconformal algebra,  $\mathfrak{sl}(4|2)$ . The superconformal primary in a supermultiplet is denoted by  $\Psi$ . There is a single conformal primary Schur operator  $\mathcal{O}_{\mathrm{Schur}}$  in each listed superconformal multiplet. The holomorphic dimension  $h$  and  $\mathfrak{gl}(1)_r$  charge  $r$  of  $\mathcal{O}_{\mathrm{Schur}}$  are given in terms of the quantum numbers  $(R,j_1,j_2)$  that label the shortened multiplet (left-most column).

This holomorphic OPE endows the vector space of Schur operators with the structure of a vertex operator algebra.

In this work, we use the notation  $\{AB\}_{n}$  to denote the coefficient of  $z_{12}^{-n}$  in the holomorphic OPE  $A(z_1)B(z_2)$ , see (2.9). We find it convenient to introduce alternative notations for the special cases  $n = 0,1$ . More precisely, we define

$$
\mathrm {N O} [ A, B ] := \{A B \} _ {0}, \quad \{\{A, B \} \} := \{A B \} _ {1}. \tag {6.4}
$$

We refer to these operations as the normal order product of  $A$  and  $B$ , and the bracket of  $A$  and  $B$ , respectively.

# 6.1.2 Gradings and filtration

The vector space  $\mathcal{V}$  of Schur operators has a triple grading by  $(h,R,r)\in \frac{1}{2}\mathbb{Z}_{+}\times \frac{1}{2}\mathbb{Z}_{+}\times \frac{1}{2}\mathbb{Z}$ ,

$$
\mathcal {V} = \bigoplus_ {h, R, r} \mathcal {V} _ {h, R, r}. \tag {6.5}
$$

The normally-ordered product preserves  $h$  and  $r$  but not  $R$ , making the  $R$  grading unnatural from the point of view of the VOA structure. However, the specifics of the twisted translation construction implies that  $R$ -charge violation occurs with a definite sign,

$$
\mathrm {N O} \left(\mathcal {V} _ {h _ {1}, R _ {1}, r _ {1}}, \mathcal {V} _ {h _ {2}, R _ {2}, r _ {2}}\right) \subseteq \bigoplus_ {k \geqslant 0} \mathcal {V} _ {h _ {1} + h _ {2}, R _ {1} + R _ {2} - k, r _ {1} + r _ {2}}. \tag {6.6}
$$

Consequently, there is a filtration by  $R$  that is preserved by the normally-ordered product. That is, if we define,

$$
\mathcal {F} _ {h, R, r} = \bigoplus_ {k \geqslant 0} \mathcal {V} _ {h, R - k, r}, \tag {6.7}
$$

then we have the following filtered property for normally-ordered multiplication,

$$
\mathrm {N O} \left(\mathcal {F} _ {h _ {1}, R _ {1}, r _ {1}}, \mathcal {F} _ {h _ {2}, R _ {2}, r _ {2}}\right) \subseteq \mathcal {F} _ {h _ {1} + h _ {2}, R _ {1} + R _ {2}, r _ {1} + r _ {2}}. \tag {6.8}
$$

In addition, the bracket operation  $\{\{\cdot ,\cdot \} \}$  defined in (6.4) obeys

$$
\{\left\{\mathcal {F} _ {h, R, r}, \mathcal {F} _ {h ^ {\prime}, R ^ {\prime}, r ^ {\prime}} \right\} \} \subseteq \mathcal {F} _ {h + h ^ {\prime} - 1, R + R ^ {\prime} - 1, r + r ^ {\prime}}, \tag {6.9}
$$

so, the bracket is filtered of tri-degree  $(-1, -1, 0)$ . The associated graded of our filtered VOA is defined in the usual fashion,

$$
\operatorname {g r} _ {\mathcal {F}} \mathcal {V} = \bigoplus_ {h, R, r} \mathcal {G} _ {h, R, r}, \quad \mathcal {G} _ {h, R, r} = \mathcal {F} _ {h, R, r} / \mathcal {F} _ {h, R - 1, r}. \tag {6.10}
$$

On this space, which is isomorphic as a vector space to  $\mathcal{V}$ , one can show that the normally-ordered product induces a grade-preserving (with respect to all the gradings) commutative, associative product, and the bracket induces an anti-symmetric bracket of tri-degree  $(-1, -1, 0)$  that obeys the Jacobi identity.

Table 4 makes it clear that knowledge of only the  $h$  and  $r$  quantum numbers state of the VOA is not sufficient to characterize uniquely the  $4d$  protected operator from which it descends—the  $R$  quantum number is also needed for a precise uplift.

# 6.1.3 The superconformal index

The superconformal index of a  $4d\mathcal{N} = 2$  SCFT  $\mathcal{T}$  (see, e.g., [45, 46]) is defined as

$$
\mathcal {I} ^ {\mathcal {T}} (p, q, t) := \operatorname {S T r} _ {\mathcal {H}} \left(p ^ {\frac {1}{2} (E + 2 j _ {1} - 2 R - r)} q ^ {\frac {1}{2} (E - 2 j _ {1} - 2 R - r ]} t ^ {R + r}\right), \tag {6.11}
$$

where  $\mathrm{STr}$  denotes the supertrace, and  $\mathcal{H}$  is the Hilbert space of local operators of the SCFT. The index receives contributions only from the states that lie in short representations of the superconformal algebra, with the contributions being such that the index is insensitive to recombinations. Let us also recall the definition of two special limits of the full superconformal index: the Macdonald index and the Schur index. The Macdonald index depends on two fugacities only, and is given by

$$
\mathcal {I} _ {\text {M a c d o n a l d}} ^ {\mathcal {T}} (q, t) := \operatorname {S T r} _ {\mathcal {H} _ {M}} \left(q ^ {E - 2 R - r} t ^ {R + r}\right), \tag {6.12}
$$

where  $\mathcal{H}_M$  denotes the subspace of states in  $\mathcal{H}$  satisfying  $E + 2j_{1} - 2R - r = 0$ . The Schur index depends on only one fugacity, and reads

$$
\mathcal {I} _ {\text {S c h u r}} ^ {\mathcal {T}} (q) := \operatorname {S T r} _ {\mathcal {H}} \left(q ^ {E - R}\right). \tag {6.13}
$$

Notice that, even if the trace is taken over the entire Hilbert space  $\mathcal{H}$ , only Schur operators contribute to the Schur limit of the superconformal index.

Under the map of [1], the Schur index of a  $4d\mathcal{N} = 2$  SCFT is mapped to the vacuum character of the associated VOA, $^{28}$

$$
\mathcal {V} = \chi [ \mathcal {T} ], \qquad \chi_ {\mathcal {V}} (q) := \mathrm {S T r} _ {\mathcal {V}} (q ^ {L _ {0}}), \qquad \chi_ {\mathcal {V}} (q) = \mathcal {I} _ {\mathrm {S c h u r}} ^ {\mathcal {T}} (q). \tag {6.14}
$$

In order to reconstruct the Macdonald index (6.12) of the parent  $4d$  theory from VOA data, we need control over the  $r$  grading and the  $R$  filtration discussed in section 6.1.2. If these pieces of data are successfully identified, the Macdonald index can be recovered via a refinement of the VOA vacuum character,

$$
\mathcal {V} = \chi [ \mathcal {T} ], \qquad \mathcal {I} _ {\mathrm {M a c d o n a l d}} ^ {\mathcal {T}} (q, t) = \chi \nu (q, t) := \mathrm {S T r} \nu (q ^ {L _ {0} - R - r} t ^ {R + r}). \qquad (6. 1 5)
$$

In the definition of  $\chi_{\mathcal{V}}(q,t)$ , the supertrace is implicitly understood to be taken on the associated graded algebra (6.10), in which both quantum numbers  $R$  and  $r$  define gradings.

No simple recipe is known to extract the  $R$ -filtration from a presentation of the VOA in terms of strong generators and their singular OPEs. In this work, we propose an efficient way to identify the  $R$ -filtration in the case in which the VOA coincides with  $\mathcal{W}_{\mathsf{G}}$  for some (crystallographic) complex reflection group  $\mathsf{G}$ . Our proposal relies on the free-field realization of the VOA. More precisely, the free-field realization of  $\mathcal{W}_{\mathsf{G}}$  allows us to introduce a filtration on  $\mathcal{W}_{\mathsf{G}}$ , denoted  $\mathcal{R}$ -filtration and defined in section 6.2. We propose the identification of this new  $\mathcal{R}$ -filtration with the sought-for  $R$ -filtration. Several checks of this proposal are discussed in sections 6.3 and 6.4.

# 6.1.4 The  $\mathcal{N} = 4$  case

Unitary irreducible highest weight representations of the four dimensional  $\mathcal{N} = 4$  superconformal algebra  $\mathfrak{psu}(2,2|4)$  are labelled by six quantum numbers  $\{E,(j_2,j_2),[q_1,p,q_2]\}$  where  $E,(j_{1},j_{2})$  are as in the previous discussion while  $[q_1,p,q_2]$  are  $\mathfrak{su}(4)_R$  Dynkin labels, see [47]. When these quantum numbers satisfy certain relations the corresponding supermultiplet obeys shortening conditions. One of the properties of the chiral algebra map  $\chi$  introduced in [1] specialized to four dimensional  $\mathcal{N} = 4$  SCFT is that it acts on irreducible representations as

$$
\chi : \operatorname {R e p s} \mathfrak {p s l} (4 | 4) \rightarrow \operatorname {R e p s} \mathfrak {p s l} (2 | 2). \tag {6.16}
$$

A little inspection shows that, using the notation of [47] for four dimensional multiplets

$$
\chi (\mathcal {B} _ {[ 0, p, 0 ]}) = \mathfrak {S} _ {h = \frac {1}{2} p}
$$

$$
\chi (\mathcal {B} _ {[ q, p, q ]}) = \mathfrak {L} _ {(h, j) = (q + \frac {1}{2} p, \frac {1}{2} p)}
$$

$$
\chi \left(\mathcal {C} _ {\left[ q _ {1}, p, q _ {2} \right], \left(j _ {1}, j _ {2}\right)}\right) = \mathfrak {L} _ {(h, j) = (\frac {1}{2} (p + q _ {1} + q _ {2}) + j _ {1} + j _ {2} + 2, \frac {1}{2} p)} \tag {6.17}
$$

$$
\chi (\mathcal {D} _ {[ q _ {1}, p, q _ {2} ], (0, j _ {2})}) = \mathfrak {L} _ {(h, j) = (\frac {1}{2} (p + q _ {1} + q _ {2}) + j _ {2} + 1, \frac {1}{2} p)}
$$

$$
\chi (\bar {\mathcal {D}} _ {[ q _ {1}, p, q _ {2} ], (j _ {1}, 0)}) = \mathfrak {L} _ {(h, j) = (\frac {1}{2} (p + q _ {1} + q _ {2}) + j _ {1} + 1, \frac {1}{2} p)}
$$

$$
\chi (\mathcal {A}) = 0.
$$

The action (6.17) can be determined by considering the decomposition of  $\mathcal{N} = 4$  superconformal multiplets in  $\mathcal{N} = 2$  superconformal multiplets and the results of table 4. The  $\mathcal{N} = 2$  superconformal algebra is embedded as

$$
\mathfrak {p s l} (4 | 4) \supset \mathfrak {s l} (4 | 2) \oplus \mathfrak {s l} (2) _ {y}, \tag {6.18}
$$

with

$$
\mathfrak {s l} (4) _ {R} \supset \mathfrak {s l} (2) _ {R} \oplus \mathfrak {s l} (2) _ {y} \oplus \mathfrak {g l} (1) _ {r}, \quad [ 1, 0, 0 ] \mapsto (\frac {1}{2}, 0) _ {+ \frac {1}{2}} \oplus (0, \frac {1}{2}) _ {- \frac {1}{2}}. \tag {6.19}
$$

The decomposition of the short supermultiplets with respect to the embedding (6.18) is given by<sup>29</sup>

$$
\mathcal {B} _ {[ q, p, q ]} \rightarrow \hat {\mathcal {B}} _ {R = \frac {p}{2} + q} \otimes \left[ \frac {p}{2} \right] + \dots \tag {6.20a}
$$

$$
\mathcal {C} _ {[ q _ {1}, p, q _ {2} ], (j _ {1}, j _ {2})} \rightarrow \hat {\mathcal {C}} _ {R = \frac {p + q _ {1} + q _ {2}}{2}, (j _ {1}, j _ {2})} \otimes \left[ \frac {p}{2} \right] + \dots \tag {6.20b}
$$

$$
\mathcal {D} _ {[ q _ {1}, p, q _ {2} ], (0, j _ {2})} \rightarrow \mathcal {D} _ {R = \frac {p + q _ {1} + q _ {2}}{2}, (0, j _ {2})} \otimes \left[ \frac {p}{2} \right] + \dots \tag {6.20c}
$$

$$
\bar {\mathcal {D}} _ {[ q _ {1}, p, q _ {2} ], (j _ {1}, 0)} \rightarrow \bar {\mathcal {D}} _ {R = \frac {p + q _ {1} + q _ {2}}{2}, (j _ {1}, 0)} \otimes \left[ \frac {p}{2} \right] + \dots \tag {6.20d}
$$

where the ... indicate terms that are either obtained by acting with the  $Q, \bar{Q}$  generators in  $\mathfrak{psl}(2|2)$  or that do not contain Schur operators<sup>30</sup>. The chiral algebra map (6.17) follows:  $h$  is the same as in the  $\mathcal{N} = 2$  case, see table 4 and  $j = \frac{p}{2}$ . To be precise, the image of the map  $\chi$  is a representation of  $\mathfrak{pl}(2|2)$ , which is the extension of  $\mathfrak{psl}(2|2)$  by the outer automorphism  $GL(1)_r \subset SL(2)$ . The  $GL(1)_r$  quantum number of the  $\mathfrak{pl}(2|2)$  primary is obtained combining (6.20) with the content of table 4. Similarly the  $R$  quantum number is found to be

$$
R \left[ \mathcal {B} _ {[ q, p, q ]} \right] = q + \frac {1}{2} p, \quad R \left[ \mathcal {C} _ {[ q _ {1}, p, q _ {2} ], (j _ {1}, j _ {2})} \right] = \frac {1}{2} (p + q _ {1} + q _ {2}) + 1, \tag {6.21a}
$$

$$
R \left[ \mathcal {D} _ {\left[ q _ {1}, p, q _ {2} \right], (0, j _ {2})} \right] = \frac {1}{2} \left(p + q _ {1} + q _ {2}\right) + \frac {1}{2}, \quad R \left[ \bar {\mathcal {D}} _ {\left[ q _ {1}, p, q _ {2} \right], (j _ {1}, 0)} \right] = \frac {1}{2} \left(p + q _ {1} + q _ {2}\right) + \frac {1}{2}, \tag {6.21b}
$$

where the first formula can be applied to all values of  $q$  including  $q = 0$ . An important feature of the map  $\chi$  is that in general it cannot be inverted since different types of four dimensional multiplets correspond to the same  $\mathfrak{pl}(2|2)$  multiplet. The pair of maps  $(\chi, R)$ , on the other hand, can be inverted.

# 6.2 The  $\mathcal{R}$ -filtration

According to our free-field construction,  $\mathcal{W}_{\Gamma}$  is realized as a subalgebra of the total  $\beta \gamma bc$  system  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)}$  generated by  $\beta_{\ell}, \gamma_{\ell}, b_{\ell}, c_{\ell}$ , with  $\ell = 1, \ldots, r$ , where  $r = \mathrm{rank}(\Gamma)$ . We first define a filtration  $\widetilde{\mathcal{R}}$  at the level of the  $\beta \gamma bc$  system. To this end, it is sufficient to assign an additive weight, referred to as  $\mathcal{R}$ -weight in the following, to the fundamental 'letters'  $\beta_{\ell}, \gamma_{\ell}, b_{\ell}, c_{\ell}$ , according to the following table,

<table><tr><td></td><td>βl</td><td>γl</td><td>bl</td><td>cl</td><td>δ</td></tr><tr><td>R</td><td>1/2pl</td><td>1-1/2pl</td><td>1/2pl</td><td>1-1/2pl</td><td>0</td></tr><tr><td>h-R</td><td>0</td><td>0</td><td>+1/2</td><td>-1/2</td><td>1</td></tr></table>

(6.22)

For convenience we reported here also the combination  $h - \mathcal{R}$ . Together with (3.4), we have thus assigned quantum numbers  $h, m, r, \mathcal{R}$  to the free fields. Given any monomial in derivatives of  $\beta_{\ell}, \gamma_{\ell}, b_{\ell}, c_{\ell}$ , its  $\mathcal{R}$ -weight is simply the sum of the  $\mathcal{R}$ -weights of its letters. Given a polynomial, we define its  $\mathcal{R}$ -weight as the maximal  $\mathcal{R}$ -weight of its monomials. Given  $\mathcal{R} \in \frac{1}{2}\mathbb{Z}$ , we may then introduce the vector space

$$
\widetilde {\mathbb {V}} _ {\mathcal {R}} := \text {l i n e a r s p a n o f s t a t e s i n} \mathbb {M} _ {\beta \gamma b c} ^ {(\Gamma)} \text {w i t h} \mathcal {R} \text {- w e i g h t} \mathcal {R} - k, k \in \mathbb {Z} _ {\geq 0}. \tag {6.23}
$$

The collection of vector spaces  $\{\widetilde{\mathbb{V}}_{\mathcal{R}}\}_{\mathcal{R}\in \frac{1}{2}\mathbb{Z}}$  defines the filtration  $\widetilde{\mathcal{R}}$  of  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)}$ ,

$$
\mathbb {M} _ {\beta \gamma b c} ^ {(\Gamma)} = \bigcup_ {R \in \frac {1}{2} \mathbb {Z}} \widetilde {\mathbb {V}} _ {\mathcal {R}}. \tag {6.24}
$$

To verify the compatibility of the filtration with the VOA structure of  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)}$ , we have to ensure that, in the OPE of two operators of  $\mathcal{R}$ -weights  $\mathcal{R}_1$  and  $\mathcal{R}_2$ , only operators with  $\mathcal{R}$ -weight  $\leq \mathcal{R}_1 + \mathcal{R}_2$  appear, both in the singular and in the regular part. This property is easily verified by noticing that OPEs in  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)}$  are computed via Wick contractions. Since the  $\mathcal{R}$ -weights of a pair  $\beta_{\ell} \gamma_{\ell}$  or  $b_{\ell} c_{\ell}$  is one, every Wick contraction decreases the  $\mathcal{R}$ -weight by one unit. The same observation implies that in the singular OPEs the  $\mathcal{R}$ -weights are strictly smaller than  $\mathcal{R}_1 + \mathcal{R}_2$ .

The filtration  $\widetilde{\mathcal{R}}$  descends to a filtration  $\mathcal{R}$  of the VOA  $\mathcal{W}_{\Gamma}$ , since the latter is a subalgebra of  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)}$ . The filtration  $\mathcal{R}$  is specified by the collection of vector spaces

$$
\mathbb {V} _ {\mathcal {R}} := \widetilde {\mathbb {V}} _ {\mathcal {R}} \cap \mathcal {W} _ {\Gamma}, \quad \mathcal {R} \in \frac {1}{2} \mathbb {Z}. \tag {6.25}
$$

Notice that in  $\mathbb{M}_{\beta \gamma bc}^{(\Gamma)}$  we can construct operators with arbitrarily low  $\mathcal{R}$ -weight. In contrast, all operators in  $\mathcal{W}_{\Gamma}$  have  $\mathcal{R} \geq 0$ . As a result, we can write

$$
\mathcal {W} _ {\Gamma} = \bigcup_ {\mathcal {R} \in \frac {1}{2} \mathbb {Z} _ {\geq 0}} \mathbb {V} _ {\mathcal {R}}. \tag {6.26}
$$

Given any filtered algebra such as the pair  $(\mathcal{W}_{\Gamma},\mathcal{R})$ , there is a standard coset construction that yields a graded algebra, the associated graded algebra  $\mathcal{G}(\mathcal{W}_{\Gamma})$ ,

$$
\mathcal {G} \left(\mathcal {W} _ {\Gamma}\right) = \bigoplus_ {\mathcal {R} \in \frac {1}{2} \mathbb {Z} _ {\geq 0}} \mathcal {H} _ {\mathcal {R}}, \quad \mathcal {H} _ {0} = \mathbb {V} _ {0} \simeq \mathbb {C}, \quad \mathbb {V} _ {\frac {1}{2}} = 0, \quad \mathcal {H} _ {\mathcal {R}} = \mathbb {V} _ {\mathcal {R}} / \mathbb {V} _ {\mathcal {R} - 1}, \quad \mathcal {R} \geq 1. \tag {6.27}
$$

Notice that  $\mathcal{G}(\mathcal{W}_{\Gamma})$  and  $\mathcal{W}_{\Gamma}$  are isomorphic as vector spaces. The graded algebra  $\mathcal{G}(\mathcal{W}_{\Gamma})$  will be useful in the next paragraph for the definition of the refined vacuum character.

It should be noticed that all generators of the small  $\mathcal{N} = 4$  algebra have  $\mathcal{R} = 1$  and they actually exhaust the  $\mathcal{R} = 1$  component, i.e.,  $\mathbb{V}_1\simeq \mathfrak{S}_1$  as  $\mathfrak{psl}(2|2)$  modules. In more generality one can argue that the action of  $\mathfrak{psl}(2|2)$  does not change the  $\mathcal{R}$ -weight. A particular example is given by operators transforming in short representations of  $\mathfrak{psl}(2|2)$ . In this case the weight of the super-multiplet can be determined unambiguously since the corresponding quasi-primary, due to the  $h = m$  condition, is necessarily a function of  $\{\beta_1,\ldots ,\beta_{\mathrm{r}}\}$  so its  $\mathcal{R}$ -weight is equal to its conformal dimension,  $\mathcal{R} = h$ . As we will see in the next subsection, this fact has a clear four-dimensional interpretation.

Refined vacuum character. By means of the filtration  $\mathcal{R}$ , we can refine the vacuum character of the VOA  $\mathcal{W}_{\Gamma}$ . This is achieved in a standard way via the associated graded algebra  $\mathcal{G}(\mathcal{W}_{\Gamma})$  defined in (6.27):

$$
\chi_ {\mathcal {W} _ {\Gamma}} (q, \xi , z) = \sum_ {h, r, m, \mathcal {R}} \operatorname {s d i m} \left(\mathcal {H} _ {h, r, m, \mathcal {R}}\right) q ^ {h} \xi^ {\mathcal {R} + r} z ^ {m}, \quad \mathcal {H} _ {\mathcal {R}} = \bigoplus_ {h, r, m} \mathcal {H} _ {h, r, m, \mathcal {R}}. \tag {6.28}
$$

Some comments are in order. The variables  $q, \xi, z$  are fugacities. The quantum numbers  $h$ ,  $m$ ,  $r$  were defined in section 2, their characterization is reported here for convenience:  $h$  is the conformal dimension,  $r$  is associated to the outer automorphisms of the small  $\mathcal{N} = 4$  and  $\mathcal{N} = 2$  SCA normalized as  $r[G] = 1/2$  and  $r[\mathcal{G}] = 1/2$  respectively, finally  $m \in \frac{1}{2}\mathbb{Z}$  is the Cartan of the  $\mathfrak{sl}(2)_y$  of the small  $\mathcal{N} = 4$  SCA or the  $\mathfrak{gl}(1)$  of the  $\mathcal{N} = 2$  SCA. The superdimension  $\mathrm{sdim}(\mathcal{H})$  is defined as  $\dim(\mathcal{H}^{(\mathrm{bos})}) - \dim(\mathcal{H}^{(\mathrm{ferm})})$ . Notice that the character (6.28) can be further refined by applying the substitution  $\xi^{\mathcal{R} + r} \mapsto \hat{\xi}_1^\mathcal{R} y^{2r}$  in (6.28).

On the computation of  $\mathsf{sdim}(\mathcal{H}_{h,r,m,\mathcal{R}})$ . Let us comment briefly on the computation of the refined character (6.28) level by level. The VOA  $\mathcal{W}_{\Gamma}$  is triply graded as

$$
\mathcal {W} _ {\Gamma} = \bigoplus_ {h \geq 0, r, m} V _ {h, r, m}, \quad V _ {h, r, m} = V _ {h, r, m} ^ {(\mathrm {b o s})} \oplus V _ {h, r, m} ^ {(\mathrm {f e r m})}, \tag {6.29}
$$

where each  $V_{h,r,m}^{(\mathrm{bos})}$  and  $V_{h,r,m}^{(\mathrm{ferm})}$  are finite-dimensional. Given any base of  $V_{h,r,m}^{(\mathrm{bos})}$ , we can determine the  $\mathcal{R}$ -weights of each element of the base using the free-field construction. Let us select a basis of  $V_{h,r,m}^{(\mathrm{bos})}$  that is minimal, in the sense that the  $\mathcal{R}$ -weights of its elements are the lowest possible. Let  $\mathcal{R}_1 < \dots < \mathcal{R}_k$  be the distinct  $\mathcal{R}$ -weights of the elements of a minimal basis, occurring with multiplicities  $n_1, \ldots, n_k$ . In an analogous fashion, we can determine a minimal basis of  $V_{h,r,m}^{(\mathrm{ferm})}$ , in which the distinct  $\mathcal{R}$ -weights are  $\mathcal{R}_1' < \dots < \mathcal{R}_{k'}'$ , occurring with multiplicities  $n_1', \ldots, n_{k'}'$ . The contribution of  $V_{h,r,m}$  to the refined character (6.28) is then

$$
\sum_ {i = 1} ^ {k} n _ {i} q ^ {h} \xi^ {\mathcal {R} _ {i} + r} z ^ {m} - \sum_ {i = 1} ^ {k ^ {\prime}} n _ {i} ^ {\prime} q ^ {h} \xi^ {\mathcal {R} _ {i} ^ {\prime} + r} z ^ {m}. \tag {6.30}
$$

# 6.3 Identification of  $\mathcal{R}$ -filtration and  $R$ -filtration

We propose the following identification:31

Suppose  $\mathcal{W}_{\Gamma}$  is the VOA associated to a four-dimensional  $\mathcal{N} = 4$  SCFT  $\mathcal{T}$ ,  $\mathcal{W}_{\Gamma} = \chi[\mathcal{T}]$ . Then the  $\mathcal{R}$ -filtration of  $\mathcal{W}_{\Gamma}$ , defined in (6.22), (6.23), (6.25), (6.26), coincides with the  $R$ -filtration canonically attached to  $\mathcal{W}_{\Gamma}$  by its definition in terms of four-dimensional CFT data of  $\mathcal{T}$ .

This proposed identification implies in particular the equality of the Macdonald index of  $\mathcal{T}$  with the refined vacuum character (6.28) of  $\mathcal{W}_{\Gamma}$ ,

$$
\chi \mathcal {W} _ {\Gamma} (q, \xi , z) = \mathcal {I} _ {\text {M a c d o n a l d}} ^ {\mathcal {T}} (q, \xi , z), \quad \mathcal {W} _ {\Gamma} = \chi [ \mathcal {T} ], \tag {6.31}
$$

where we have refined the  $4d$  Macdonald index using the  $\mathfrak{sl}(2)_y$  flavor fugacity  $z$ . To provide evidence for (6.31), we followed two approaches. Firstly, we performed a direct match of the Macdonald index of  $4d\mathcal{N} = 4$  SYM with gauge algebra  $\mathfrak{a}_1$ ,  $\mathfrak{a}_2$ , and the refined vacuum character of  $\mathcal{W}_{\Gamma}$  for  $\Gamma = A_1, A_2$ , working up to and including terms  $q^5$ ,  $q^3$  respectively. Secondly, we analyzed the Hall-Littlewood limits of the Macdonald index of  $4d\mathcal{N} = 4$  SYM with gauge algebra  $\mathfrak{g}$  and of the refined vacuum character of  $\mathcal{W}_{\mathrm{Weyl}(\mathfrak{g})}$ , and we provided general arguments for their equality, for any simple Lie algebra  $\mathfrak{g}$ .

The Macdonald index of  $4d\mathcal{N} = 4$  SYM. The definition of the Macdonald index was recalled in (6.12). For  $4d\mathcal{N} = 4$  SYM with gauge group  $G$ , the Macdonald index can be written as an integral over the maximal torus of  $G$  as follows

$$
\mathcal {I} _ {\text {M a c d o n a l d}} ^ {\text {S Y M} (\mathfrak {g})} (q, \xi , z) = \int [ d u ] \text {P . E .} \left[ \xi^ {\frac {1}{2}} \chi_ {\mathfrak {S} _ {\frac {1}{2}}} ^ {\mathfrak {p l} (2 | 2)} (q, \xi , z) \chi_ {\text {A d j}} ^ {\mathfrak {g}} (u) \right]. \tag {6.32}
$$

The notation deserves some explanations:  $\mathfrak{g} = \operatorname{Lie}(G)$  and  $[du]$  is the corresponding normalized Haar measure, the  $\mathfrak{pl}(2|2)$  character of the extra-short representation  $\mathfrak{S}_{\frac{1}{2}}$  is given by

$$
\chi_ {\mathfrak {S} _ {\frac {1}{2}}} ^ {\mathrm {p l} (2 | 2)} (q, \xi , z) = \frac {q ^ {\frac {1}{2}}}{1 - q} \left(z ^ {\frac {1}{2}} - z ^ {- \frac {1}{2}}\right) - \frac {q}{1 - q} \left(\xi^ {\frac {1}{2}} - \xi^ {- \frac {1}{2}}\right). \tag {6.33}
$$

Finally the plethystic exponential is defined as

$$
\mathrm {P . E .} [ f (q, \xi , z, u) ] = \exp \left[ \sum_ {m = 1} ^ {\infty} \frac {1}{m} f \left(q ^ {m}, \xi^ {m}, z ^ {m}, u ^ {m}\right) \right]. \tag {6.34}
$$

Match between index and character for  $\Gamma = A_{1}$ . In the case  $\Gamma = A_{1}$ , the VOA  $\mathcal{W}_{\Gamma}$  is (the simple quotient of) the small  $\mathcal{N} = 4$  SCA with central charge  $c = -9$ . To compute the refined vacuum character (6.28) up to and including  $q^{5}$  terms, we performed a counting of states with definite quantum numbers  $h$ ,  $j$ ,  $r$ ,  $\mathcal{R}$ , up to  $h = 5$ . The results are collected in table 5. Clearly, all states with integer  $h$  are bosons, and all states with half-integer  $h$  are fermions. Let us stress that this counting differs from the analogous counting at generic central charge  $c$ , because of the states that become null upon setting  $c = -9$ . Fortunately, the counting of states is facilitated by the fact that, in the free-field realization, all null states are automatically zero, as proven in [26].

<table><tr><td>h</td><td>su(2) multiplets with quantum numbers (j)rR</td></tr><tr><td>0</td><td>(0)00</td></tr><tr><td>1/2</td><td>—</td></tr><tr><td>1</td><td>(1)10</td></tr><tr><td>3/2</td><td>(1/2)1±1/2</td></tr><tr><td>2</td><td>(2)20+(1)10+(0)10</td></tr><tr><td>5/2</td><td>(3/2)2±1/2+(1/2)1±1/2</td></tr><tr><td>3</td><td>(3)30+(2)20+2(1)02+(1)10+(0)10</td></tr><tr><td>7/2</td><td>(5/2)3±1/2+2(3/2)2±1/2+(1/2)2±1/2+(1/2)1±1/2</td></tr><tr><td>4</td><td>(4)40+(3)30+2(2)03+2(2)02+4(1)02+(1)01+2(0)02+(0)01+(1)2±1</td></tr><tr><td>9/2</td><td>(7/2)4±1/2+2(5/2)3±1/2+(3/2)3±1/2+3(3/2)2±1/2+(1/2)1±1/2</td></tr><tr><td>5</td><td>(5)50+(4)40+2(3)04+2(3)03+5(2)03+2(2)02+2(1)03+7(1)02+(1)01+3(0)02+(0)1+(2)3±1+(1)2±1+(0)2±1</td></tr></table>

Table 5: States in the VOA  $\mathcal{W}_{\Gamma}$  with  $\Gamma = A_{1}$  with  $h \leq 5$ . The notation  $(j)_{r}^{\mathcal{R}}$  encodes the spin  $j$  under the  $\mathfrak{sl}(2)_y$  symmetry of the small  $\mathcal{N} = 4$  SCA, the outer automorphism quantum number  $r$ , and the  $\mathcal{R}$ -weight. When we write, for instance,  $(\frac{1}{2})_{\pm 1/2}^{1}$ , we mean that multiplets with  $r = 1/2$  and  $r = -1/2$  are found with the same multiplicity.

Upon assembling the refined vacuum character from the data of table 5, we found a perfect match with the Macdonald index of  $4d\mathcal{N} = 4$  SYM with gauge algebra  $\mathfrak{a}_1$ , up to  $q^5$  terms. The Macdonald index is computed using (6.32).

The content of table 5 can also be encoded more compactly in the expression

$$
\chi_ {A _ {1}} (q, \hat {\xi}, y, z) = 1 + \chi_ {\mathrm {S}} ^ {(A _ {1})} + \hat {\xi} ^ {2} \mathfrak {L} _ {(3, 1)} + \hat {\xi} ^ {3} \mathfrak {L} _ {(4, 2)} + \hat {\xi} ^ {2} \mathfrak {L} _ {(4, 0)} + \hat {\xi} ^ {4} \mathfrak {L} _ {(5, 3)} + \hat {\xi} ^ {3} \mathfrak {L} _ {(5, 2)} + (\hat {\xi} ^ {3} + \hat {\xi} ^ {2}) \mathfrak {L} _ {(5, 1)} + \dots \tag {6.35}
$$

The RHS is obtained from (6.28) by means of the further refinement  $\xi^{R + r} \mapsto \hat{\xi}_1^R y^{2r}$ . On the LHS,  $\chi_{\mathrm{S}}^{(A_1)} = \sum_{R=1}^{\infty} \hat{\xi}^R \mathfrak{S}_R$  and  $\mathfrak{S}_R(q,y,z)$ ,  $\mathfrak{L}_{(h,j)}(q,y,z)$  denote the  $\mathfrak{pl}(2|2)$  character of the corresponding representations introduced in table 2.

Match between index and character for  $\Gamma = A_2$ . The counting of states up to  $h = 3$  in the case  $\Gamma = A_2$  is summarized in table 6. In contrast to the case  $\Gamma = A_1$ , this time we have a non-trivial interplay between bosons and fermions with the same weight  $h$ . Once again, it is essential to take into account the states that become null for  $c = -24$ , but luckily this task is efficiently performed with the help of the free-field realization.

Constructing the refined vacuum character (6.28) up to  $q^3$  from the data in table (6), we find a perfect match with the Macdonald index of  $4d\mathcal{N} = 4$  SYM with gauge algebra  $\mathfrak{a}_2$ , as given by (6.32).

We can repackage the content of table 6 in the compact expression

$$
\chi_ {A _ {2}} (q, \hat {\xi}, y, z) = 1 + \chi_ {S} ^ {(A _ {2})} + \hat {\xi} ^ {2} \mathfrak {L} _ {(2, 0)} + (\hat {\xi} ^ {2} + \hat {\xi} ^ {3}) \mathfrak {L} _ {(3, 1)} + \hat {\xi} ^ {5 / 2} \mathfrak {L} _ {(\frac {5}{2}, \frac {3}{2})} + \dots \tag {6.36}
$$

where we are using the same notation as is (6.35), and on the RHS  $\chi_{\mathrm{S}}^{(A_2)} = \chi_{\mathrm{S}}^{(A_1)} + \hat{\xi}^{3 / 2}\mathfrak{S}_{\frac{3}{2}}+$ $\hat{\xi}^{5 / 2}\mathfrak{S}_{\frac{5}{2}} + \hat{\xi}^{3}\mathfrak{S}_{3} + \dots$

<table><tr><td>h</td><td>bosons</td><td>fermions</td></tr><tr><td>0</td><td>(0)00</td><td>—</td></tr><tr><td>1/2</td><td>—</td><td>—</td></tr><tr><td>1</td><td>(1)10</td><td>—</td></tr><tr><td>3/2</td><td>(3/2)3/20</td><td>(1/2)1±1/2</td></tr><tr><td>2</td><td>(2)20+(1)10+(0)20+(0)10</td><td>(1)3/2±1/2</td></tr><tr><td>5/2</td><td>(5/2)5/2+(3/2)05/2+(3/2)03/2+(1/2)03/2</td><td>(3/2)2±1/2+(1/2)2±1/2+(1/2)1±1/2</td></tr><tr><td>3</td><td>2(3)03+(2)02+(1)03+3(1)02+(1)01+2(0)02+(0)01+(0)2±1</td><td>2(2)5/2±1/2+(1)5/2±1/2+(1)3/2±1/2</td></tr></table>

Table 6: States in the VOA  $\mathcal{W}_{\Gamma}$  with  $\Gamma = A_{2}$  with  $h \leq 3$ . The notation  $(j)_{r}^{\mathcal{R}}$  encodes the spin  $j$  under the  $\mathfrak{sl}(2)_y$  symmetry of the small  $\mathcal{N} = 4$  SCA, the outer automorphism quantum number  $r$ , and the  $\mathcal{R}$ -weight.

# 6.4 The Hall-Littlewood limit

In this section we give additional evidence that the  $R$ -filtration, defined in terms of the 4d parent theory SCFT data in section 6, and the  $\mathcal{R}$ -filtration, defined in terms of the free-field realization in section 6.2, coincide by focusing on a subsector of operators known as Hall-Littlewood (HL) chiral ring [45]. This ring is obtained by restricting to operators satisfying the condition  $h = R + r = \mathcal{R} + r$ . As the Higgs branch chiral ring, the HL chiral ring carries the structure of a Poisson algebra, see [13]. At the level of the character and of the Macdonald index, the HL limit corresponds to

$$
\chi_ {\mathcal {W} _ {\Gamma}} ^ {\mathrm {H L}} (\tau , x) = \lim _ {q \to 0} \chi_ {\mathcal {W} _ {\Gamma}} (q, q ^ {- 1} \tau^ {2}, x ^ {2}),
$$

$$
\mathcal {I} _ {\mathrm {H L}} ^ {\mathrm {S Y M} (\mathfrak {g})} (\tau , x) = \lim  _ {q \rightarrow 0} \mathcal {I} _ {\text {M a c d o n a l d}} ^ {\mathrm {S Y M} (\mathfrak {g})} \left(q, q ^ {- 1} \tau^ {2}, x ^ {2}\right), \tag {6.37}
$$

where the character  $\chi_{\mathcal{W}_{\Gamma}}$  is defined in (6.28) and the index in (6.32).

In the following we describe more explicitly the HL chiral ring as obtained from the free-field description by applying the map  $\mathcal{P}'$  introduced in (3.18) and the four dimensional HL chiral ring denoted by  $\mathcal{R}_{HL}[\mathrm{SYM}_{\mathfrak{g}}]$ . As a consistency check of our proposed equality  $\mathcal{W}_{\mathrm{Weyl}(\mathfrak{g})} = \chi[\mathrm{SYM}_{\mathfrak{g}}]$  and equivalence of  $\mathcal{R}$ - and  $R$ -filtrations we have

$$
\mathcal {P} ^ {\prime} \left(\mathcal {W} _ {\Gamma}\right) \simeq \mathscr {R} _ {H L} \left[ \mathrm {S Y M} _ {\mathfrak {g}} \right], \quad \Gamma = \operatorname {W e y l} (\mathfrak {g}), \tag {6.38}
$$

where the right hand side is given in (6.42). We verify this isomorphism explicitly for  $\Gamma = A_1, I_2(p)$  with  $p = 3, 4, 6$ .

The HL chiral ring from the free-field description. Restricting to HL operators, defined by the condition  $h = \mathcal{R} + r$ , in the free-field description is straightforward. According to the weight assignments given in (3.4), (6.22) all the constituent  $\{\beta_{\ell}, \gamma_{\ell}, b_{\ell}, c_{\ell}\}$  satisfy the condition  $h = \mathcal{R} + r$ , while adding a derivative will violate this condition<sup>33</sup>. This implies that the (candidate) HL chiral ring is obtained by setting to zero derivatives in the free-field realization of the strong generators of  $\mathcal{W}_{\Gamma}$ . This is the definition of the map  $\mathcal{P}'$  introduced in (3.18) so that the candidate HL chiral ring is  $\mathcal{P}'(\mathcal{W}_{\Gamma})$  with the Poisson ring structure defined in section 3.4.

<table><tr><td></td><td>βl</td><td>γl</td><td>b1</td><td>c1</td><td>θ</td></tr><tr><td>h-(R+r)</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td></tr><tr><td>h-(R-r)</td><td>0</td><td>0</td><td>+1</td><td>-1</td><td>1</td></tr></table>

as follows from (3.4) and (6.22).

Interestingly, this ring (conjecturally) admits an alternative description as the subring of  $\mathbb{C}[\beta, \gamma, b, c]$  annihilated by a certain set of nilpotent operators which are interpreted as the image of the screening charges of the VOA under the map  $\mathcal{P}'$ . In equation

$$
\mathcal {P} ^ {\prime} \left(\mathcal {W} _ {\Gamma}\right) = \operatorname {K e r n e l} \left(J _ {\mathrm {c l}} ^ {(1)}\right) \cap \dots \cap \operatorname {K e r n e l} \left(J _ {\mathrm {c l}} ^ {(\mathrm {r})}\right), \tag {6.40}
$$

where  $\mathsf{J}_{\mathrm{cl}}^{(\ell)}$  acts by Poisson brackets  $\{\mathsf{J}_{\mathrm{cl}}^{(\ell)},\cdot \}_{\mathrm{PB}}$ . Some examples are given in the end of this section.

The HL chiral ring  $\mathcal{R}_{HL}[\mathrm{SYM}_{\mathfrak{g}}]$ . Let us recall that the HL chiral ring is a consistent truncation of the usual  $\mathcal{N} = 1$  chiral ring obtained by restricting to Schur operators. The  $\mathcal{N} = 1$  chiral ring  $\mathcal{R}_{\mathcal{N} = 1}[\mathrm{SYM}_{\mathfrak{g}}]$  of  $\mathcal{N} = 4$  SYM is given by

$$
\mathcal {R} _ {\mathcal {N} = 1} \left[ \mathrm {S Y M} _ {\mathfrak {g}} \right] = \mathbb {C} \left[ \mathcal {M} _ {\Gamma} ^ {\prime \prime} \right], \quad \mathcal {M} _ {\Gamma} ^ {\prime \prime} = \frac {\mathbb {C} ^ {(3 | 2)} \otimes V _ {\Gamma}}{\Gamma}, \tag {6.41}
$$

with  $V_{\Gamma} \simeq \mathbb{R}^{\mathrm{rank}(\Gamma)}$ , generalizing (1.9) and (2.12), see [48], [49]. By restricting to Schur operators we conclude that

$$
\mathcal {R} _ {H L} \left[ \mathrm {S Y M} _ {\mathfrak {g}} \right] = \mathcal {R} _ {\Gamma} ^ {\prime} := \mathbb {C} \left[ \mathcal {M} _ {\Gamma} ^ {\prime} \right], \tag {6.42}
$$

$$
\mathcal {R} _ {\Gamma} ^ {\prime} := \mathbb {C} \left[ z _ {1} ^ {+}, \dots , z _ {r} ^ {+}, z _ {1} ^ {-}, \dots , z _ {r} ^ {-}, \theta_ {1}, \dots , \theta_ {r} \right] ^ {\Gamma}, \quad \mathcal {M} _ {\Gamma} ^ {\prime} = \frac {\mathbb {C} ^ {(2 | 1)} \otimes V _ {\Gamma}}{\Gamma}. \tag {6.43}
$$

As  $\mathcal{R}_{\Gamma}$  the Hall-Littlewood chiral ring  $\mathcal{R}_{\Gamma}^{\prime}$  carries the action of  $GL(2)$ . This is actually extended to the action of  $SL(2|1)$ .

It should be remarked that the HL chiral ring admits an alternative description which involves solving a certain BRST cohomology problem. This is the truncation of the BRST definition of the VOA  $\chi[\mathrm{SYM}_{\mathfrak{g}}]$ , see [1]. It is a non-trivial fact<sup>34</sup> that this definition reproduces the HL chiral ring (6.42). Further evidence of this equivalence can be obtained by matching the corresponding Hilbert series. The Hilbert series of (6.42) can be computed using the Molien formula<sup>35</sup>

$$
\mathsf {H S} _ {\Gamma} ^ {\prime} (\tau , x) = \frac {1}{| \Gamma |} \sum_ {g \in \Gamma} \frac {1}{\operatorname {s d e t} _ {\mathbb {C} ^ {(2 | 1)} \otimes V _ {\Gamma}} (1 - h \otimes g)}, \qquad h = \left( \begin{array}{c c c} \tau x & 0 & 0 \\ 0 & \tau x ^ {- 1} & 0 \\ 0 & 0 & \tau^ {2} \end{array} \right). \qquad \tag {6.45}
$$

The BRST definition of the HL chiral ring on the other hand implies that its Hilbert series is obtained by taking the HL limit given in (6.37) of the integral (6.32) and gives

$$
\mathcal {I} _ {\mathrm {H L}} ^ {\mathrm {S Y M} (\mathfrak {g})} (\tau , x) = \int [ d u ] \mathrm {P . E .} \left[ \left((x + x ^ {- 1}) \tau - \tau^ {2}\right) \chi_ {\mathrm {A d j}} ^ {\mathfrak {g}} (u) \right]. \tag {6.46}
$$

It is possible to further refine (6.45) by keeping track of the  $\mathfrak{gl}(1)_r$  quantum number. This is achieved by taking  $h = \mathrm{diag}(\tau x, \tau x^{-1}, z\tau^2)$ , where  $z$  is a  $\mathfrak{gl}(1)_r$  fugacity. This refinement of the Hilbert series cannot be obtained as a specialization of the four dimensional index.

Since the argument of the plethystic exponential is a Laurant polynomial, the integrand is a rational function, see (B.14). The equivalence of the two descriptions implies that

$$
\mathcal {I} _ {\mathrm {H L}} ^ {\mathrm {S Y M} (\mathfrak {g})} (\tau , x) = \mathsf {H S} _ {\Gamma} ^ {\prime} (\tau , x), \quad \Gamma = \operatorname {W e y l} (\mathfrak {g}), \tag {6.47}
$$

See appendix B.2 for more details.

# 6.4.1 Examples

Rank 1 example:  $\Gamma = A_{1}$ . Let us describe the entries of (6.38) in this example. In this case the VOA  $\mathcal{W}_{A_1}$  is isomorphic to the simple quotient of the small  $\mathcal{N} = 4$  super-Virasoro algebra at  $c = -9$ . Applying the map  $\mathcal{P}'$  to the free-field realization of the generators of  $\mathcal{W}_{A_1}$  gives

$$
j ^ {+} = \beta , \quad j ^ {0} = 2 \beta \gamma + b c, \quad j ^ {-} = \gamma (\beta \gamma + b c), \quad g ^ {+} = b, \quad g ^ {-} = b \gamma , \tag {6.48}
$$

where we used the notation  $\mathcal{P}'(X) = x$  when acting on the generators. Above  $(\beta, \gamma, b, c)$  are super-commuting variables, they can be thought of as coordinates of  $\mathbb{C}^{2|2}$ . It should be noticed that the image of  $T$  and  $\widetilde{G}^{\pm}$  is zero since derivatives are set to zero. It is easy to verify that the following combinations

$$
\operatorname {N u l l s} _ {A _ {1}} = \left\{j ^ {+} j ^ {-} - \frac {1}{4} j ^ {0} j ^ {0}, j ^ {\pm} g ^ {\mp} - \frac {1}{2} j ^ {0} g ^ {\pm}, g ^ {+} g ^ {-} \right\} \tag {6.49}
$$

are identically zero. Recall that the relations  $g^{\pm}g^{\pm} = 0$  and  $g^{+}g^{-} + g^{-}g^{+} = 0$  hold by definition in a supercommutative ring. We conclude that

$$
\mathcal {P} ^ {\prime} \left(\mathcal {W} _ {A _ {1}}\right) \simeq \frac {\mathbb {C} \left[ j ^ {+} , j ^ {0} , j ^ {-} , g ^ {+} , g ^ {-} \right]}{\operatorname {N u l l s} _ {A _ {1}} = 0}. \tag {6.50}
$$

The Hilbert series of this ring can be computed by standard methods, the result is given by (6.52) below. Next let us consider the ring (6.43) with  $\Gamma = A_{1}$ . It is rather clear the it is generated by

$$
j ^ {+} = z ^ {+} z ^ {+}, \quad j ^ {0} = 2 z ^ {+} z ^ {-}, \quad j ^ {-} = z ^ {-} z ^ {-}, \quad g ^ {+} = z ^ {+} \theta , \quad g ^ {-} = z ^ {-} \theta . \tag {6.51}
$$

The character of this ring is given by the Molien formula (6.45) and in this case gives

$$
\mathsf {H S} _ {A _ {1}} ^ {\prime} (\tau , x) = \frac {1 - (x + x ^ {- 1}) \tau^ {3} - \tau^ {4} + (x + x ^ {- 1}) \tau^ {5}}{(1 - \tau^ {2} x ^ {- 2}) (1 - \tau^ {2}) (1 - \tau^ {2} x ^ {+ 2})}. \tag {6.52}
$$

This expression is a generalization of (2.15). Similarly to the Higgs branch chiral ring, the most efficient method to show the equivalence (6.38) is to establish a relation between the  $\beta \gamma bc$  coordinates and the quotient coordinates by equating the genetators (6.48) with (6.51). This gives

$$
\beta = z ^ {+} z ^ {+}, \quad \gamma - \frac {1}{2} \beta^ {- 1} b c = \frac {z ^ {-}}{z ^ {+}}, \quad b = z ^ {+} \theta . \tag {6.53}
$$

Notice that the expressions (6.48) are invariant under the transformation  $\gamma \mapsto \gamma +\eta b,c\mapsto c + 2\eta \beta$  where  $\eta$  is a fermionic parameter. For this reason only the invariant combinations under this transformation, given in (6.53), can be determined. As already anticipated in (6.40),  $\mathcal{P}'(\mathcal{W}_{A_1})$  can be identified with the kernel of  $\mathsf{J}_{\mathrm{cl}} = b\beta^{-1 / 2}$ , which generates the fermionic symmetry described above, in  $\mathbb{C}[\beta ,\gamma ,b,c]$ .

Rank 2 example:  $\Gamma = I_2(p)$ . Let us start from describing the generators of  $\mathcal{R}_{I_2(p)}'$ . They are given by  $j(y), w(y)$  defined in (4.14) together with their  $\mathfrak{sl}(2|1)$  fermionic partners

$$
g (y) = z _ {1} (y) \theta_ {2} + z _ {2} (y) \theta_ {1}, \quad g _ {w} (y) = p \left(z _ {1} (y) ^ {p} \theta_ {1} + z _ {2} (y) ^ {p} \theta_ {2}\right). \tag {6.54}
$$

These generators should be compared to the image of the  $\mathcal{W}_{I_2(p)}$  generators, namely (4.20) and the  $\mathfrak{psl}(2|2)$  descendants of  $W^{\mathrm{h.w.}} = \beta_{2}$ , under  $\mathcal{P}'$ . The resulting expressions are rather involved but one can show that they coincide with (4.14), (6.54) upon using the identification (4.28) and

$$
b _ {1} = z _ {1} ^ {+} \theta_ {2} + z _ {2} ^ {+} \theta_ {1}, \quad b _ {2} = p \left(\left(z _ {1} ^ {+}\right) ^ {p - 1} \theta_ {1} + \left(z _ {2} ^ {+}\right) ^ {p - 1} \theta_ {2}\right). \tag {6.55}
$$

This concludes the proof of the isomorphism (6.38) for  $\Gamma = I_2(p)$ .

As in the case  $\Gamma = A_1$ , see (6.53), also in this case the expression of  $\gamma_{1},\gamma_{2}$  can be determined only up to certain nilpotent quantities. This is due to the fact that  $\mathcal{P}'(x)$  with  $x\in \mathcal{W}_{I_2(p)}$  are invariant under the transformations

$$
\begin{array}{l} \gamma_ {1} \mapsto \gamma_ {1} + \Lambda \eta \beta_ {1} ^ {p - 2} b _ {1} + \tilde {\eta} b _ {2}, \\ \gamma_ {2} \mapsto \gamma_ {2} + \eta b _ {2} + \tilde {\eta} b _ {1}, \\ c _ {1} \mapsto c _ {1} + 2 \Lambda \eta \beta_ {1} ^ {p - 1} + p \tilde {\eta} \beta_ {2}, \\ c _ {2} \mapsto c _ {2} + \frac {p}{p - 1} \eta \beta_ {2} + \frac {2}{p - 1} \tilde {\eta} \beta_ {1}, \\ \end{array}
$$

where  $\eta, \tilde{\eta}$  are fermionic parameters. These transformations are generated by  $^36$ $\mathsf{J}_{\mathrm{cl}}^{\pm}$  given in (4.30) with the suffix  $\pm$  referring to the two solutions of the differential equation (4.32).

# Acknowledgments

The authors are grateful to Philip Argyres, Chris Beem, Mario Martone for useful conversations and correspondence. The work of CM is supported in part by grant #494786 from the Simons Foundation. The work of LR is supported in part by NSF Grant PHY-1620628. We thank the Galileo Galilei Institute for Theoretical Physics for the hospitality and the INFN for partial support.

# A Some basic facts on OPEs and VOAs

For the convenience of the reader, in this appendix we collect some standard facts about OPEs in VOAs, with emphasis on the cases with  $\mathcal{N} = 2$  and small  $\mathcal{N} = 4$  supersymmetry.

Covariance under  $\mathfrak{sl}(2)_z$ . Given the operators  $A(z)$ ,  $B(z)$ , we adopt the notation

$$
A \left(z _ {1}\right) B \left(z _ {2}\right) = \sum_ {n \in \mathbb {Z}} \frac {\left\{A B \right\} _ {n} \left(z _ {2}\right)}{z _ {1 2} ^ {n}}, \quad z _ {1 2} := z _ {1} - z _ {2}. \tag {A.1}
$$

An operator  $\mathcal{O}(z)$  is an  $\mathfrak{sl}(2)_z$  primary (or quasi-primary) if it transforms tensorially under the global part of the conformal group on the Riemann sphere, denoted  $SL(2)_z$ ,

$$
\mathcal {O} ^ {\prime} \left(z ^ {\prime}\right) = \left(\frac {\partial z ^ {\prime}}{\partial z}\right) ^ {- h _ {\mathcal {O}}} \mathcal {O} (z), \tag {A.2}
$$

where  $h_{\mathcal{O}}$  is the conformal dimension of  $\mathcal{O}$  and

$$
z ^ {\prime} = \frac {a z + b}{c z + d}, \quad \left( \begin{array}{l l} a & b \\ c & d \end{array} \right) \in S L (2) _ {z}. \tag {A.3}
$$

An  $\mathfrak{sl}(2)_z$  primary operator  $\mathcal{O}$  satisfies

$$
\{T \mathcal {O} \} _ {3} = 0, \quad \{T \mathcal {O} \} _ {2} = h _ {\mathcal {O}} \mathcal {O}, \quad \{T \mathcal {O} \} _ {1} = \partial \mathcal {O}, \tag {A.4}
$$

where  $T$  is the stress tensor.

The OPE of  $\mathfrak{sl}(2)_z$  primary operators  $\mathcal{O}_1$ ,  $\mathcal{O}_2$  with dimensions  $h_1$ ,  $h_2$  is constrained by  $\mathfrak{sl}(2)_z$  covariance to take the form

$$
\mathcal {O} _ {1} \left(z _ {1}\right) \mathcal {O} _ {2} \left(z _ {2}\right) = \sum_ {\mathcal {O} \in B _ {h}} \lambda_ {\mathcal {O} _ {1} \mathcal {O} _ {2}} ^ {\mathcal {O}} \frac {1}{z _ {1 2} ^ {h _ {1} + h _ {2} - h}} \mathcal {D} _ {h _ {1}, h _ {2}; h} \left(z _ {1 2}, \partial_ {z _ {2}}\right) \mathcal {O} \left(z _ {2}\right), \tag {A.5}
$$

where  $B_{h}$  is a basis in the space of  $\mathfrak{sl}(2)_z$  primary operators with dimension  $h$ ,  $\lambda_{\mathcal{O}_1\mathcal{O}_2}^{\mathcal{O}}$  are OPE coefficients, and the differential operator  $\mathcal{D}_{h_1,h_2;h}(z_{12},\partial_{z_2})$  is given by

$$
\mathcal {D} _ {h _ {1}, h _ {2}; h} \left(z _ {1 2}, \partial_ {z _ {2}}\right) = \sum_ {k = 0} ^ {\infty} \frac {\left(h + h _ {1} - h _ {2}\right) _ {k}}{k ! \left(2 h\right) _ {k}} z _ {1 2} ^ {k} \partial_ {z _ {2}} ^ {k}. \tag {A.6}
$$

The quantity  $(x)_k$  is the ascending Pochhammer symbol,  $(x)_k = \prod_{i=0}^{k-1}(x + i)$ .

Notice that, if  $\mathcal{O}_1$ ,  $\mathcal{O}_2$  are  $\mathfrak{sl}(2)_z$  primaries, the operators  $\{\mathcal{O}_1\mathcal{O}_2\}_n$ ,  $n \in \mathbb{Z}$  are not necessarily  $\mathfrak{sl}(2)_z$  primaries. As a consequence of (A.5), however,  $\{\mathcal{O}_1\mathcal{O}_2\}_n$  is generically given as the sum of an  $\mathfrak{sl}(2)_z$  primary of weight  $h_1 + h_2 - n$  and derivatives of other  $\mathfrak{sl}(2)_z$  primaries of lowest weight. There is a standard formula for extracting the  $\mathfrak{sl}(2)_z$  primary with  $h = h_1 + h_2 - n$  from  $\{\mathcal{O}_1\mathcal{O}_2\}_n$ . In our normalization conventions, the formula reads

$$
\left(\mathcal {O} _ {1} \mathcal {O} _ {2}\right) _ {n} (z) = \sum_ {p \geq 0} \mathcal {K} _ {h _ {1}, h _ {2}, n, p} \partial_ {z} ^ {p} \left\{\mathcal {O} _ {1} \mathcal {O} _ {2} \right\} _ {n + p} (z), \tag {A.7}
$$

where

$$
\mathcal {K} _ {h _ {1}, h _ {2}, n, p} = \frac {(-) ^ {p} (2 h _ {1} - n - p) _ {p}}{p ! (2 h _ {1} + 2 h _ {2} - 2 n - p - 1) _ {p}}. \tag {A.8}
$$

Our normalization is chosen in such a way that

$$
\mathcal {O} _ {1} \left(z _ {1}\right) \mathcal {O} _ {2} \left(z _ {2}\right) = \sum_ {n \in \mathbb {Z}} \frac {1}{z _ {1 2} ^ {h _ {1} + h _ {2} - h}} \mathcal {D} _ {h _ {1}, h _ {2}; h} \left(z _ {1 2}, \partial_ {z _ {2}}\right) \left(\mathcal {O} _ {1} \mathcal {O} _ {2}\right) _ {n} \left(z _ {2}\right). \tag {A.9}
$$

Covariance under  $\mathfrak{sl}(2)_y$  The operator content of VOAs with small  $\mathcal{N} = 4$  superconformal symmetry falls into representations of  $\mathfrak{sl}(2)$  R-symmetry. We find it convenient to study irreps of  $\mathfrak{sl}(2)$  by means of a standard index-free formalism, based on the introduction of an auxiliary variable  $y$ . The group  $SL(2)$  acts on  $y$  via Möbius transformations. We often denote this  $SL(2)$  group as  $SL(2)_y$ , in order to distinguish it from the global part of the conformal group on the Riemann sphere, which is denoted  $SL(2)_z$ . An object  $\mathcal{O}(y)$  transforms in the irrep of  $\mathfrak{sl}(2)_y$  with spin  $j$  if it satisfies

$$
\mathcal {O} ^ {\prime} \left(y ^ {\prime}\right) = \left(\frac {\partial y ^ {\prime}}{\partial y}\right) ^ {j} \mathcal {O} (y), \tag {A.10}
$$

where we suppressed the  $z$  dependence, and

$$
y ^ {\prime} = \frac {\hat {a} y + \hat {b}}{\hat {c} y + \hat {d}}, \quad \left( \begin{array}{l l} \hat {a} & \hat {b} \\ \hat {c} & \hat {d} \end{array} \right) \in S L (2) _ {y}. \tag {A.11}
$$

As a function of  $y$ ,  $\mathcal{O}$  is a polynomial of degree  $2j$ .

Consider any objects  $\mathcal{O}_1(y),\mathcal{O}_2(y)$  transforming under  $SL(2)_y$  according to (A.10) with spins  $j_{1},j_{2}$ . Let  $\mathcal{B}$  denote any bilinear operation. We are mainly interested in the cases  $\mathcal{B}(\cdot ,\cdot) = \{\cdot \cdot \}_{n}$  or  $\mathcal{B}(\cdot ,\cdot) = (\cdot \cdot)_n$  in a VOA, but the following considerations also apply if  $\mathcal{B}$  is, for instance, the commutative product in a polynomial ring, or the Poisson bracket in a Poisson algebra. It is useful to have a formula to decompose the product  $\mathcal{B}(\mathcal{O}_1(y_1),\mathcal{O}_2(y_2))$  into contributions of definite spin  $j$ . Such a formula reads

$$
\mathcal {B} \left(\mathcal {O} _ {1} \left(y _ {1}\right), \mathcal {O} _ {2} \left(y _ {2}\right)\right) = \sum_ {j} y _ {1 2} ^ {j _ {1} + j _ {2} - j} \widehat {\mathcal {D}} _ {j _ {1}, j _ {2}; j} \left(y _ {1 2}, \partial_ {y _ {2}}\right) \mathcal {B} \left(\mathcal {O} _ {1}, \mathcal {O} _ {2}\right) ^ {j} \left(y _ {2}\right). \tag {A.12}
$$

Some comments are in order. The range of the summation over  $j$  is determined by the usual rules for composing angular momenta,

$$
j \in \left\{\left| j _ {1} - j _ {2} \right|, \left| j _ {1} - j _ {2} \right| + 1, \dots , j _ {1} + j _ {2} \right\}. \tag {A.13}
$$

The differential operator  $\widehat{\mathcal{D}}_{j_1,j_2;j}(y_{12},\partial_{y_2})$  is given by

$$
\widehat {\mathcal {D}} _ {j _ {1}, j _ {2}; j} \left(y _ {1 2}, \partial_ {y _ {2}}\right) = \sum_ {k = 0} ^ {\infty} \frac {\left(- j - j _ {1} + j _ {2}\right) _ {k}}{k ! (- 2 j) _ {k}} y _ {1 2} ^ {k} \partial_ {y _ {2}} ^ {k}, \tag {A.14}
$$

and can be thought of as the continuation of the operator (A.6) to  $h = -j$ . As usual,  $y_{12} = y_1 - y_2$ . Notice that the sum over  $k$  always truncates to a finite sum. The object  $\mathcal{B}(\mathcal{O}_1,\mathcal{O}_2)^j$  is the projection onto the part with definite spin  $j$ . In order to define it more precisely, we introduce the notation

$$
\mathcal {O} _ {1} (y) = \sum_ {k _ {1} = 0} ^ {2 j _ {1}} \mathcal {O} _ {1, k _ {1}} y ^ {k _ {1}}, \quad \mathcal {O} _ {2} (y) = \sum_ {k _ {2} = 0} ^ {2 j _ {2}} \mathcal {O} _ {1, k _ {2}} y ^ {k _ {2}}. \tag {A.15}
$$

We may then write

$$
\mathcal {B} \left(\mathcal {O} _ {1}, \mathcal {O} _ {2}\right) ^ {j} \left(y _ {2}\right) = \sum_ {k _ {1} = 0} ^ {2 j _ {1}} \sum_ {k _ {2} = 0} ^ {2 j _ {2}} \mathcal {C} _ {j _ {1}, j _ {2}, j, k _ {1}, k _ {2}} \mathcal {B} \left(\mathcal {O} _ {1, k _ {1}}, \mathcal {O} _ {2, k _ {2}}\right), \tag {A.16}
$$

where the coefficients  $\mathcal{C}$  are given by

$$
\mathcal {C} _ {j _ {1}, j _ {2}, j, k _ {1}, k _ {2}} = \sum_ {r = 0} ^ {j _ {1} + j _ {1} - j} \frac {(-) ^ {r} (j _ {1} - j _ {2} + j + r) _ {r} ^ {\downarrow}}{r ! (2 j + r + 1) _ {r} ^ {\downarrow}} \frac {1}{(j _ {1} + j _ {2} - j - r) !} \sum_ {s = 0} ^ {r} \binom {p} {s} \left(k _ {1}\right) _ {j _ {1} + j _ {2} - j - s} ^ {\downarrow} \left(k _ {2}\right) _ {s} ^ {\downarrow}. \tag {A.17}
$$

We have used the descending Pochhammer symbol,  $(x)_{k}^{\downarrow} = \prod_{i = 0}^{k - 1}(x - i)$ .

A compact notation for OPEs. As discussed above, covariance under  $\mathfrak{sl}(2)_z$  and  $\mathfrak{sl}(2)_y$  completely fixes the way the  $\partial_z$  and  $\partial_y$  derivative of a  $\mathfrak{sl}(2)_z$ ,  $\mathfrak{sl}(2)_y$  primary operator enter the OPE of two  $\mathfrak{sl}(2)_z$ ,  $\mathfrak{sl}(2)_y$  primary operators. This allows us to use a compact notation in which all factors  $z_{12}$ ,  $y_{12}$  and all terms with  $\partial_z$  and/or  $\partial_y$  derivatives are omitted. If needed, such elements can be easily reconstructed unambiguously with the formulae recorded above. For example, the non-trivial OPEs of the small  $\mathcal{N} = 4$  algebra at level  $k$  in compact notation read

$$
J J = - k \operatorname {i d} + 2 J, \quad J G = G, \quad J \widetilde {G} = \widetilde {G},
$$

$$
T J = J, \qquad T G = \frac {3}{2} G, \qquad T \widetilde {G} = \frac {3}{2} \widetilde {G},
$$

$$
T T = 3 k \mathrm {i d} + 2 T, \quad G \widetilde {G} = - 2 k \mathrm {i d} + 2 J - T. \tag {A.18}
$$

To make contact with the OPEs is section 2.1 in the main text, one uses the parametrization

$$
J (y) = J ^ {+} + J ^ {0} y + J ^ {-} y ^ {2}, \quad G (y) = G ^ {+} + G ^ {-} y, \quad \widetilde {G} (y) = \widetilde {G} ^ {+} + \widetilde {G} ^ {-} y. \tag {A.19}
$$

In a completely analogous fashion, the full  $\mathcal{N} = 2$  SCA at level  $k$  is encoded in the non-trivial OPEs

$$
\mathcal {J} \mathcal {J} = 2 k \text {i d}, \qquad \qquad \mathcal {J} \mathcal {G} = - \mathcal {G}, \qquad \qquad \mathcal {J} \widetilde {\mathcal {G}} = \widetilde {\mathcal {G}},
$$

$$
\mathcal {T} \mathcal {J} = \mathcal {J}, \qquad \qquad \qquad \mathcal {T} \mathcal {G} = \frac {3}{2} \mathcal {G}, \qquad \qquad \qquad \mathcal {T} \widetilde {\mathcal {G}} = \frac {3}{2} \widetilde {\mathcal {G}},
$$

$$
\mathcal {T} \mathcal {T} = 3 k \mathrm {i d} + 2 \mathcal {T}, \quad \mathcal {G} \widetilde {\mathcal {G}} = - 2 k \mathrm {i d} + \mathcal {J} - \mathcal {T}. \tag {A.20}
$$

Primary conditions. Let us summarize the different notions of primary operators encountered in this work in the case of VOAs with small  $\mathcal{N} = 4$  supersymmetry.

-  $\mathfrak{sl}(2)_z$  primary of dimension  $h$ :

$$
\{T \mathcal {O} (y) \} _ {3} = 0, \quad \{T \mathcal {O} (y) \} _ {2} = h \mathcal {O} (y), \quad \{T \mathcal {O} (y) \} _ {1} = \partial_ {z} \mathcal {O} (y). \tag {A.21}
$$

- Virasoro primary of dimension  $h$ :

$$
\left\{T \mathcal {O} (y) \right\} _ {n \geq 3} = 0, \quad \left\{T \mathcal {O} (y) \right\} _ {2} = h \mathcal {O} (y), \quad \left\{T \mathcal {O} (y) \right\} _ {1} = \partial_ {z} \mathcal {O} (y). \tag {A.22}
$$

Operator with definite  $\mathfrak{sl}(2)_y$  spin  $j$ :

$$
\left\{J \left(y _ {1}\right) \mathcal {O} \left(y _ {2}\right) \right\} _ {1} = 2 j y _ {1 2} \widehat {D} _ {1, j; j} \mathcal {O} \left(y _ {2}\right). \tag {A.23}
$$

- AKM primary with definite  $\mathfrak{sl}(2)_y$  spin  $j$ :

$$
\{J \left(y _ {1}\right) \mathcal {O} \left(y _ {2}\right) \} _ {n \geq 2} = 0, \quad \{J \left(y _ {1}\right) \mathcal {O} \left(y _ {2}\right) \} _ {1} = 2 j y _ {1 2} \hat {D} _ {1, j; j} \mathcal {O} \left(y _ {2}\right). \tag {A.24}
$$

-  $\mathfrak{psl}(2|2)$  primary with quantum numbers  $(h,j)$ :

$$
\{T \mathcal {O} (y) \} _ {3} = 0, \quad \{T \mathcal {O} (y) \} _ {2} = h \mathcal {O} (y), \quad \{T \mathcal {O} (y) \} _ {1} = \partial_ {z} \mathcal {O} (y), \tag {A.25}
$$

$$
\{J \left(y _ {1}\right) \mathcal {O} \left(y _ {2}\right) \} _ {1} = 2 j y _ {1 2} \widehat {D} _ {1, j; j} \mathcal {O} \left(y _ {2}\right), \tag {A.26}
$$

$$
\left\{G \left(y _ {1}\right) \mathcal {O} \left(y _ {2}\right) \right\} _ {2} = 0, \quad \left\{\widetilde {G} \left(y _ {1}\right) \mathcal {O} \left(y _ {2}\right) \right\} _ {2} = 0. \tag {A.27}
$$

- small  $\mathcal{N} = 4$  super-Virasoro primary with quantum numbers  $(h,j)$ :

$$
\{T \mathcal {O} (y) \} _ {n \geq 3} = 0, \quad \{T \mathcal {O} (y) \} _ {2} = h \mathcal {O} (y), \quad \{T \mathcal {O} (y) \} _ {1} = \partial_ {z} \mathcal {O} (y), \tag {A.28}
$$

$$
\{J \left(y _ {1}\right) \mathcal {O} \left(y _ {2}\right) \} _ {n \geq 2} = 0, \quad \{J \left(y _ {1}\right) \mathcal {O} \left(y _ {2}\right) \} _ {1} = 2 j y _ {1 2} \widehat {D} _ {1, j; j} \mathcal {O} \left(y _ {2}\right), \tag {A.29}
$$

$$
\left\{G \left(y _ {1}\right) \mathcal {O} \left(y _ {2}\right) \right\} _ {n \geq 2} = 0, \quad \left\{\widetilde {G} \left(y _ {1}\right) \mathcal {O} \left(y _ {2}\right) \right\} _ {n \geq 2} = 0. \tag {A.30}
$$

The differential operator  $\widehat{D}_{1,j;j}$  was defined in (A.14) and takes the simple form

$$
\widehat {D} _ {1, j; j} \left(y _ {1 2}, \partial_ {y _ {2}}\right) = 1 + \frac {1}{2 j} y _ {1 2} \partial_ {y _ {2}}. \tag {A.31}
$$

It is useful to notice that the AKM primary condition, combined with the  $\mathfrak{psl}(2|2)$  primary condition, is equivalent to the super-Virasoro primary condition.

The analogous notions of primary operators in the case with  $\mathcal{N} = 2$  supersymmetry are obtained with minimal modifications. The R-symmetry of the small  $\mathcal{N} = 2$  SCA is  $\mathfrak{gl}(1)$ , and therefore we do not need the auxiliary variable  $y$ .

Operator with definite  $\mathfrak{gl}(1)$  charge  $m$ :

$$
\{\mathcal {I} \mathcal {O} \} _ {1} = 2 m \mathcal {O}. \tag {A.32}
$$

- AKM primary with definite  $\mathfrak{gl}(1)$  charge  $m$ :

$$
\{\mathcal {J} \mathcal {O} \} _ {n \geq 2} = 0, \quad \{\mathcal {J} \mathcal {O} \} _ {1} = 2 m \mathcal {O}. \tag {A.33}
$$

-  $\mathfrak{osp}(2|2)$  primary with quantum numbers  $(h,m)$ :

$$
\{T \mathcal {O} \} _ {3} = 0, \quad \{T \mathcal {O} \} _ {2} = h \mathcal {O}, \quad \{T \mathcal {O} \} _ {1} = \partial_ {z} \mathcal {O}, \tag {A.34}
$$

$$
\{\mathcal {J} \mathcal {O} \} _ {1} = 2 m \mathcal {O}, \tag {A.35}
$$

$$
\left\{\mathcal {G} \mathcal {O} \right\} _ {2} = 0, \quad \left\{\widetilde {\mathcal {G}} \mathcal {O} \right\} _ {2} = 0. \tag {A.36}
$$

-  $\mathcal{N} = 2$  super-Virasoro primary with quantum numbers  $(h,m)$ :

$$
\{T \mathcal {O} \} _ {n \geq 3} = 0, \quad \{T \mathcal {O} \} _ {2} = h \mathcal {O}, \quad \{T \mathcal {O} \} _ {1} = \partial_ {z} \mathcal {O}, \tag {A.37}
$$

$$
\left\{\mathcal {J} \mathcal {O} \right\} _ {n \geq 2} = 0, \quad \left\{\mathcal {J} \mathcal {O} \right\} _ {1} = 2 m \mathcal {O}, \tag {A.38}
$$

$$
\left\{\mathcal {G} \mathcal {O} \right\} _ {n \geq 2} = 0, \quad \left\{\widetilde {\mathcal {G}} \mathcal {O} \right\} _ {n \geq 2} = 0. \tag {A.39}
$$

In analogy with the previous case, the AKM primary condition, combined with the  $\mathfrak{osp}(2|2)$  primary condition, is equivalent to the  $\mathcal{N} = 2$  super-Virasoro primary condition.

Superconformal multiplets. The action of the fermionic generators of  $\mathfrak{psl}(2|2)$  is encoded in the OPE of an operator  $\mathcal{O}$  with the supersymmetry currents  $G$ ,  $\widetilde{G}$ . The fermionic generators of  $Q$  type are encoded in the order-one pole of the OPE, while the fermionic generators of  $S$  type are encoded in the order-two pole. We are mainly interested in the action of generators of  $Q$  type. To describe it efficiently, we introduce the notation

$$
G ^ {\uparrow} \mathcal {O} := (G \mathcal {O}) _ {1} ^ {j + \frac {1}{2}}, \quad G ^ {\downarrow} \mathcal {O} := (G \mathcal {O}) _ {1} ^ {j - \frac {1}{2}}, \quad \widetilde {G} ^ {\uparrow} \mathcal {O} := (\widetilde {G} \mathcal {O}) _ {1} ^ {j + \frac {1}{2}}, \quad \widetilde {G} ^ {\downarrow} \mathcal {O} := (\widetilde {G} \mathcal {O}) _ {1} ^ {j - \frac {1}{2}}, \tag {A.40}
$$

where  $\mathcal{O}$  has spin  $j$ . Suppose  $W$  is a  $\mathfrak{psl}(2|2)$  primary operator with  $h = j$ . The corresponding supersymmetry multiplet is a short multiplet, denoted  $\mathfrak{S}_h$ . The content of  $\mathfrak{S}_h$  is the following

W

$$
G _ {W} := G ^ {\downarrow} W \quad \widetilde {G} _ {W} := \widetilde {G} ^ {\downarrow} W \tag {A.41}
$$

$$
T _ {W} := - G ^ {\downarrow} \widetilde {G} ^ {\downarrow} W.
$$

Let us now consider a  $\mathfrak{psl}(2|2)$  primary operator  $X$  with generic  $h > j$ . In this case, the relevant supersymmetry multiplet is a long multiplet denoted  $\mathfrak{L}_{h,j}$ . Its content can be presented in the following way,

$$
j - 1 \qquad j - \frac {1}{2} \qquad \qquad j \qquad \qquad j + \frac {1}{2} \qquad \qquad j + 1
$$

$h$

$$
h + \frac {1}{2} \qquad \qquad G ^ {\downarrow} X  ,   \tilde {G} ^ {\downarrow} X \qquad \qquad \qquad \qquad G ^ {\uparrow} X  ,   \tilde {G} ^ {\uparrow} X
$$

$$
h + 1   G ^ {\downarrow} \tilde {G} ^ {\downarrow} X \qquad \qquad \qquad G ^ {\downarrow} G ^ {\uparrow} X  ,     \tilde {G} ^ {\downarrow} \tilde {G} ^ {\uparrow} X  ,     G ^ {\downarrow} \tilde {G} ^ {\uparrow} X  ,     \tilde {G} ^ {\downarrow} G ^ {\uparrow} X \qquad \qquad G ^ {\uparrow} \tilde {G} ^ {\uparrow} X
$$

$$
h + \frac {3}{2} \qquad \qquad G ^ {\downarrow} \tilde {G} ^ {\downarrow} \tilde {G} ^ {\uparrow} X, \tilde {G} ^ {\downarrow} G ^ {\downarrow} G ^ {\uparrow} X \qquad \qquad G ^ {\uparrow} \tilde {G} ^ {\downarrow} \tilde {G} ^ {\uparrow} X, \tilde {G} ^ {\uparrow} G ^ {\downarrow} G ^ {\uparrow} X
$$

$$
h + 2 \quad G ^ {\downarrow} \tilde {G} ^ {\downarrow} G ^ {\uparrow} \tilde {G} ^ {\uparrow} X
$$

The cases  $j = 1/2$  and  $j = 0$  deserve special attention. If  $j = 1/2$ , the state  $G^{\downarrow} \widetilde{G}^{\downarrow} X$ , which would have spin  $-1/2$ , is identically zero. In the case  $j = 0$ , all states that would have negative spin in the above table are identically zero. Moreover, the spin-0 states  $G^{\downarrow} \widetilde{G}^{\uparrow} X$  and  $\widetilde{G}^{\downarrow} G^{\uparrow} X$  become linearly dependent because of the identity

$$
G ^ {\downarrow} \tilde {G} ^ {\uparrow} X = \tilde {G} ^ {\downarrow} G ^ {\uparrow} X. \tag {A.42}
$$

Similar considerations apply to the case with  $\mathcal{N} = 2$  supersymmetry. In that case, we introduce the notation

$$
\mathcal {G} \cdot \mathcal {O} = (\mathcal {G O}) _ {1}, \quad \widetilde {\mathcal {G}} \cdot \mathcal {O} = (\widetilde {\mathcal {G O}}) _ {1}. \tag {A.43}
$$

If we start with an  $\mathfrak{osp}(2|2)$  primary  $X$  with generic  $h \neq m$ , the relevant supersymmetry multiplet is non-chiral and denoted  $\mathfrak{X}_{h,m}$ . Its content is

X

$$
\mathcal {G} _ {X} := \mathcal {G} \cdot X \quad \widetilde {\mathcal {G}} _ {X} := \widetilde {\mathcal {G}} \cdot X \tag {A.44}
$$

$$
\mathcal {T} _ {X} := - \mathcal {G} \cdot (\tilde {\mathcal {G}} \cdot X).
$$

In the special cases  $h = \pm m$ , the primary  $X$  is annihilated by the action of  $\mathcal{G}$  or  $\widetilde{\mathcal{G}}$ , and therefore we obtain a chiral or antichiral multiplet  $\mathfrak{C}_h$ ,  $\bar{\mathfrak{C}}_h$ , which only contains two states.

# B Hilbert series and indices

# B.1 Molien series

Recall the definition of the plethystic logarithm applied to the Hilbert series (2.14):

$$
\mathrm {P L} _ {\Gamma} (\tau , x) = \sum_ {k = 1} ^ {\infty} \frac {\mu (k)}{k} \log \left(\operatorname {M o l i e n} _ {\Gamma} \left(\tau^ {k}, x ^ {k}\right)\right), \tag {B.1}
$$

where  $\mu(k)$  is the Möbius function. It is convenient to remove the contribution of short generators from  $\mathsf{PL}_{\Gamma}$  and define

$$
\mathrm {X} _ {\Gamma} (\tau , x) := \mathrm {P L} _ {\Gamma} (\tau , x) - \sum_ {\ell = 1} ^ {r} \chi_ {\frac {p _ {\ell}}{2}} (x) \tau^ {p _ {\ell}}, \tag {B.2}
$$

where  $r = \mathrm{rank}(\Gamma)$  and the degrees of the invariants  $\{p_1, \ldots, p_r\}$  given in table 3 and  $\chi_j$  are  $SL(2)$  characters

$$
\chi_ {j} := \chi_ {j} (x) = \frac {x ^ {2 j + 1} - x ^ {- 2 j - 1}}{x - x ^ {- 1}}. \tag {B.3}
$$

Lets collect the  $\tau$  expansion of  $\mathsf{X}_{\Gamma}(\tau ,x)$  for all Coxeter groups up to the first relation, i.e., the first negative sign in the expansion in  $GL(1)\times SL(2)$  characters:

$$
A _ {1}: \quad - \tau^ {4}. \tag {B.4a}
$$

$$
A _ {N - 1}: \quad - \chi_ {\frac {N}{2} - 1} \tau^ {N + 2} + \dots , \tag {B.4b}
$$

$$
B _ {N} / C _ {N}: \quad - \tau^ {2 N + 2} \sum_ {\ell = 1} ^ {\left\lfloor \frac {N + 1}{2} \right\rfloor} \chi_ {N - 2 \ell + 1} + \dots , \tag {B.4c}
$$

$$
\left\{ \begin{array}{l l} D _ {4}: & \tau^ {6} - \tau^ {8} (\chi_ {2} + \chi_ {1} + 2) + \dots \\ D _ {5}: & \tau^ {7} \chi_ {1 / 2} - \tau^ {9} \chi_ {1 / 2} + \dots \\ D _ {6}: & \tau^ {8} \chi_ {1} - \tau^ {1 2} (\chi_ {4} + 2 \chi_ {2} + \chi_ {1} + 2) + \dots \\ D _ {7}: & \tau^ {9} \chi_ {3 / 2} + \tau^ {1 1} \chi_ {1 / 2} - \tau^ {1 3} (\chi_ {3 / 2} + \chi_ {1 / 2}) + \dots \\ D _ {8}: & \tau^ {1 0} \chi_ {2} + \tau^ {1 2} (\chi_ {1} + 1) - \tau^ {1 4} \chi_ {1} + \dots \\ D _ {9}: & \tau^ {1 1} \chi_ {5 / 2} + \tau^ {1 3} (\chi_ {3 / 2} + \chi_ {1 / 2}) - \tau^ {1 7} (\chi_ {5 / 2} + \chi_ {3 / 2} + 2 \chi_ {1 / 2}) + \dots \end{array} \right. \tag {B.4d}
$$

$$
E _ {6}: \quad \tau^ {8} + \tau^ {9} \chi_ {\frac {3}{2}} + \tau^ {1 2} \chi_ {3} - \tau^ {1 1} \chi_ {\frac {1}{2}} + \dots , \tag {B.4e}
$$

$$
E _ {7}: \quad \tau^ {1 0} \chi_ {1} + \tau^ {1 2} \chi_ {3} + \tau^ {1 4} \chi_ {3} + \tau^ {1 6} \chi_ {5} - \tau^ {1 6} \chi_ {2} + \dots , \tag {B.4f}
$$

$$
E _ {8}: \quad \tau^ {1 2} + \tau^ {1 4} \chi_ {3} + \tau^ {1 8} \left(\chi_ {6} + \chi_ {4} + \chi_ {3}\right) + \tau^ {2 0} \left(\chi_ {6} + 1\right) - \tau^ {2 2} \chi_ {4} + \dots \tag {B.4g}
$$

$$
F _ {4}: \quad \tau^ {8} + \tau^ {1 2} \chi_ {3} - \tau^ {1 2} \chi_ {2} - \tau^ {1 4} \left(\chi_ {5} + \chi_ {4} + \chi_ {3} + \chi_ {2} + \chi_ {1}\right) + \dots \tag {B.4h}
$$

$$
H _ {3}: \quad - \tau^ {1 2} (\chi_ {4} + \chi_ {2} + 1) + \dots \tag {B.4i}
$$

$$
H _ {4}: \quad \tau^ {1 2} + \tau^ {2 0} (\chi_ {6} + 1) - \tau^ {2 4} (\chi_ {8} + \chi_ {6} + \chi_ {4} + \chi_ {2} + 1) + \dots \tag {B.4j}
$$

$$
I _ {2} (p): \quad - \chi_ {\frac {p}{2} - 1} \tau^ {p + 2} + \dots , \tag {B.4k}
$$

The Weyl group  $G_{2} = I_{2}(6)$ . From the expressions in (B.4) one immediately reads off the long generators from table 3 as well as the quantum number of the lightest relation by identifying terms of the form  $\tau^{2h}\chi_j(x)$  with quantities with quantum numbers  $(h,j)$ . Notice that for  $E_{6,7,8}, H_{4}$  there are short generators with smaller  $h$ -weight than the lightest relation.

The Molien series for all Coxeter groups are recorded in an ancillary Mathematica file. They are obtained either by direct summation over the elements of the Coxeter group $^{37}$  or by summing over conjugacy classes and using the results of [51]. The second approach is particularly convenient in type  $E_{6,7,8}$  for which the order of the Weyl group  $|\Gamma| = \prod_{\ell} p_{\ell}$  is very large.

For the classical cases  $A$  and  $B / C$  the Molien series can be extracted from a simple

generating function in the following way. Let

$$
\mathcal {Z} _ {A} \left(p, v _ {1}, v _ {2}\right) := \frac {1}{\prod_ {n , m = 0} ^ {\infty} \left(1 - p v _ {1} ^ {n} v _ {2} ^ {m}\right)} = \exp \left(\sum_ {k = 1} ^ {\infty} \frac {p ^ {k}}{k} \frac {1}{\left(1 - v _ {1} ^ {k}\right) \left(1 - v _ {2} ^ {k}\right)}\right), \tag {B.5a}
$$

$$
\mathcal {Z} _ {B / C} \left(p, v _ {1}, v _ {2}\right) := \mathcal {Z} _ {A} \left(p, v _ {1} ^ {2}, v _ {2} ^ {2}\right) \mathcal {Z} _ {A} \left(p v _ {1} v _ {2},, v _ {1} ^ {2}, v _ {2} ^ {2}\right). \tag {B.5b}
$$

The following relations hold

$$
\mathcal {Z} _ {A} (p, \tau x, \tau x ^ {- 1}) = 1 + z _ {u (1)} \sum_ {k = 1} ^ {\infty} p ^ {k} \operatorname {M o l i e n} _ {A _ {k - 1}} (\tau , x), \tag {B.6a}
$$

$$
\mathcal {Z} _ {B / C} (p, \tau x, \tau x ^ {- 1}) = 1 + \sum_ {k = 1} ^ {\infty} p ^ {k} \operatorname {M o l i e n} _ {B / C _ {k}} (\tau , x), \tag {B.6b}
$$

where  $z_{u(1)} = (1 - x^{-1}\tau)^{-1}(1 - x\tau)^{-1}$  and  $\mathsf{Molien}_{A_0} = 1$

The ring  $\mathcal{R}_{\Gamma}$  as a finitely generated  $\mathcal{I}_{\Gamma}$ -module. It is interesting to observe that the Hilbert series of  $\mathcal{R}_{\Gamma}$  defined in (2.13) can be rewritten as

$$
\operatorname {M o l i e n} _ {\Gamma} (\tau , x) = \mathcal {Z} _ {\mathrm {C B}} (\tau x) \mathcal {Z} _ {\mathrm {C B}} \left(\tau x ^ {- 1}\right) \left(\mathrm {P} _ {\Gamma} (\tau , x) + \tau^ {| \Phi_ {\Gamma} |} \mathrm {P} _ {\Gamma} \left(\tau^ {- 1}, x\right)\right), \tag {B.7}
$$

where

$$
\mathcal {Z} _ {\mathrm {C B}} (y) = \prod_ {\ell = 1} ^ {r} \frac {1}{1 - y ^ {p _ {\ell}}} \tag {B.8}
$$

and  $|\Phi_{\Gamma}|$  is the cardinality of the root system associated to  $\Gamma$ .  $\mathsf{P}_{\Gamma}(\tau, x)$  is a polynomial in  $\tau$  and satisfies  $2\mathsf{P}_{\Gamma}(1,1) = |\Gamma|$ . For example, for  $r = 1$  one has  $\mathsf{P}_{A_1}(\tau, x) = 1$ . In rank two

$$
\mathsf {P} _ {I _ {2} (p)} (\tau , x) + \tau^ {2 p} \mathsf {P} _ {I _ {2} (p)} (\tau^ {- 1}, x) = \frac {1 - \tau^ {2 p + 2}}{1 - \tau^ {2}} + \frac {x ^ {p - 1} - x ^ {1 - p}}{x - x ^ {- 1}} \tau^ {p}. \tag {B.9}
$$

The presentation (B.7) can be interpreted as follows. Let

$$
\mathscr {I} _ {\Gamma} := \mathbb {C} \left[ z _ {1} ^ {+}, \dots , z _ {r} ^ {+} \right] ^ {\Gamma} \otimes \mathbb {C} \left[ z _ {1} ^ {-}, \dots , z _ {r} ^ {-} \right] ^ {\Gamma} = \mathbb {C} \left[ \mathcal {E} _ {1} ^ {+}, \dots , \mathcal {E} _ {r} ^ {+}, \mathcal {E} _ {1} ^ {-}, \dots , \mathcal {E} _ {r} ^ {-} \right], \tag {B.10}
$$

where  $\mathcal{E}_{\ell}^{\pm}$  are the  $\Gamma$ -invariants and are algebraically independent. The ring of invariants (2.12) is a finitely generated free  $\mathcal{I}_{\Gamma}$ -module. The Hilbert series (B.7) makes this fact manifest. Let us finally observe that

$$
\operatorname {M o l i e n} _ {\Gamma} \left(\tau^ {- 1}, x\right) = \tau^ {2 \operatorname {r a n k} (\Gamma)} \operatorname {M o l i e n} _ {\Gamma} (\tau , x). \tag {B.11}
$$

This can be easily verified from the form (B.7) and (2.11).

# B.2 The index in some limits

In this appendix we include some details about the integral representations of the index (6.32) in two limits.

Coulomb branch index. As a warm up lets consider a limit of the index (6.32) which reproduced the Hilbert series of  $\mathbb{C}[z_1,\ldots ,z_r]^\Gamma$ . It is given by

$$
\mathcal {I} _ {\mathrm {C B}} ^ {\mathrm {S Y M} (\mathfrak {g})} (t) = \int [ d u ] \text {P . E .} \left[ t \chi_ {\mathrm {A d j}} ^ {\mathfrak {g}} (u) \right] = \frac {1}{| \Gamma |} \frac {1}{(1 - t) ^ {\mathfrak {r}}} \int \frac {d ^ {\mathfrak {r}} u}{(2 \pi i u) ^ {\mathfrak {r}}} \prod_ {\alpha \in \Phi_ {\Gamma}} \left(\frac {1 - u ^ {\alpha}}{1 - t u ^ {\alpha}}\right) = \prod_ {\ell = 1} ^ {\mathfrak {r}} \frac {1}{1 - t ^ {p _ {\ell}}} \tag {B.12}
$$

where  $\mathsf{r} = \mathrm{rank}(\Gamma)$ ,  $\Gamma = \operatorname{Weyl}(\mathfrak{g})$  and  $p_{\ell}$  are the degrees of the invariants, see table 3. The explicit evaluation of the integral (B.12) is non-trivial. The simplest example is

$$
\mathcal {I} _ {\mathrm {C B}} ^ {\mathrm {S Y M} (\mathfrak {s u} (2))} (t) = \frac {1}{2} \frac {1}{(1 - t)} \int_ {0} ^ {2 \pi} \frac {d \theta}{2 \pi} \frac {(1 - e ^ {i \theta}) (1 - e ^ {- i \theta})}{(1 - t e ^ {i \theta}) (1 - t e ^ {- i \theta})} = \frac {1}{1 - t ^ {2}}. \qquad \mathrm {(B . 1 3)}
$$

This expression is valid for  $|t| < 1$ . The integral above can be computed by summing the two residues.

Hall-Littlewood index. The HL index given in (6.46) can be massaged to the form

$$
\mathcal {I} _ {\mathrm {H L}} ^ {\mathrm {S Y M} (\mathfrak {g})} (\tau , x) = \frac {1}{| \Gamma |} \left(\frac {(1 - \tau^ {2})}{(1 - x \tau) (1 - x ^ {- 1} \tau)}\right) ^ {r} \int \frac {d ^ {r} u}{(2 \pi i u) ^ {r}} \prod_ {\alpha \in \Phi_ {\Gamma}} \frac {(1 - u ^ {\alpha}) (1 - \tau^ {2} u ^ {\alpha})}{(1 - x \tau u ^ {\alpha}) (1 - x ^ {- 1} \tau u ^ {\alpha})}. \tag {B.14}
$$

In this example one should take  $|\tau x|, |\tau x^{-1}| < 1$ . We checked in various examples that the integral (B.14) coincides with the Molien series of  $\mathcal{R}_{\Gamma}^{\prime}$  defined in (6.45). It is likely that the equality could be proven by showing that the residues in (B.14) are in one to one correspondence with elements of the Weyl group.

A trained eye might recognize that for  $\mathfrak{g} = \mathfrak{su}(N)$  the index (B.14) is related to the Nekrasov instanton partition function, see e.g. equation (14) in [52]. Introduce  $a_{N}$  via the generating function

$$
\sum_ {N = 0} ^ {\infty} a _ {N} \left(q _ {1}, q _ {2}\right) z ^ {N} = \exp \left(\sum_ {N = 1} ^ {\infty} \frac {1 - q _ {1} ^ {N} q _ {2} ^ {N}}{\left(1 - q _ {1} ^ {N}\right) \left(1 - q _ {2} ^ {N}\right)} \frac {z ^ {N}}{N}\right). \tag {B.15}
$$

The expressions (B.14) can be integrated to

$$
\mathcal {I} _ {\mathrm {H L}} ^ {\operatorname {S Y M} (\mathfrak {s u} (N))} (\tau , x) = \frac {(1 - x \tau) (1 - x ^ {- 1} \tau)}{(1 - \tau^ {2})} a _ {N} (x \tau , x ^ {- 1} \tau). \tag {B.16}
$$

The relative factor in this equations is interpreted as a  $\mathcal{I}_{\mathrm{HL}}^{\mathrm{SYM}(\mathfrak{u}(1))}(\tau, x)$ . Similarly, the generating series for the  $B_{r}$  family reads

$$
1 + \sum_ {r = 1} ^ {\infty} \mathcal {I} _ {\mathrm {H L}} ^ {\mathrm {S Y M} \left(\mathfrak {b} _ {r}\right)} (\tau , x) = \exp \left(\sum_ {k = 1} ^ {\infty} \frac {1 + q _ {1} ^ {k} q _ {2} ^ {k} - q _ {1} ^ {k} q _ {2} ^ {k} \left(q _ {1} ^ {k} + q _ {2} ^ {k}\right)}{\left(1 - q _ {1} ^ {2 k}\right) \left(1 - q _ {2} ^ {2 k}\right)} \frac {z ^ {k}}{k}\right), \tag {B.17}
$$

with  $q_{1} = x^{+1}\tau$ ,  $q_{2} = x^{-1}\tau$ .

# References

[1] C. Beem, M. Lemos, P. Liendo, W. Peelaers, L. Rastelli, and B. C. van Rees, “Infinite Chiral Symmetry in Four Dimensions,” Commun. Math. Phys. 336 no. 3, (2015) 1359–1433, arXiv:1312.5344 [hep-th].  
[2] C. Beem, W. Peelaers, L. Rastelli, and B. C. van Rees, “Chiral algebras of class S,” JHEP 05 (2015) 020, arXiv:1408.6522 [hep-th].  
[3] M. Lemos and W. Peelaers, “Chiral Algebras for Trinion Theories,” JHEP 02 (2015) 113, arXiv:1411.3252 [hep-th].  
[4] M. Lemos and P. Liendo, “ $\mathcal{N} = 2$  central charge bounds from  $2d$  chiral algebras,” JHEP 04 (2016) 004, arXiv:1511.07449 [hep-th].  
[5] S. Cecotti, J. Song, C. Vafa, and W. Yan, "Superconformal Index, BPS Monodromy and Chiral Algebras," JHEP 11 (2017) 013, arXiv:1511.01516 [hep-th].  
[6] T. Arakawa and K. Kawasetsu, “Quasi-lisse vertex algebras and modular linear differential equations,” arXiv:1610.05865 [math.QA].  
[7] F. Bonetti and L. Rastelli, "Supersymmetric localization in  $\mathrm{AdS}_5$  and the protected chiral algebra," JHEP 08 (2018) 098, arXiv:1612.06514 [hep-th].  
[8] J. Song, “Macdonald Index and Chiral Algebra,” JHEP 08 (2017) 044, arXiv:1612.08956 [hep-th].  
[9] L. Fredrickson, D. Pei, W. Yan, and K. Ye, “Argyres-Douglas Theories, Chiral Algebras and Wild Hitchin Characters,” JHEP 01 (2018) 150, arXiv:1701.08782 [hep-th].  
[10] C. Cordova, D. Gaiotto, and S.-H. Shao, "Surface Defects and Chiral Algebras," JHEP 05 (2017) 140, arXiv:1704.01955 [hep-th].  
[11] J. Song, D. Xie, and W. Yan, “Vertex operator algebras of Argyres-Douglas theories from M5-branes,” JHEP 12 (2017) 123, arXiv:1706.01607 [hep-th].  
[12] M. Buican, Z. Laczko, and T. Nishinaka, “ $\mathcal{N} = 2$  S-duality revisited,” JHEP 09 (2017) 087, arXiv:1706.03797 [hep-th].  
[13] C. Beem and L. Rastelli, “Vertex operator algebras, Higgs branches, and modular differential equations,” JHEP 08 (2018) 114, arXiv:1707.07679 [hep-th].  
[14] Y. Pan and W. Peelaers, “Chiral Algebras, Localization and Surface Defects,” JHEP 02 (2018) 138, arXiv:1710.04306 [hep-th].  
[15] M. Fluder and J. Song, “Four-dimensional Lens Space Index from Two-dimensional Chiral Algebra,” JHEP 07 (2018) 073, arXiv:1710.06029 [hep-th].  
[16] J. Choi and T. Nishinaka, “On the chiral algebra of Argyres-Douglas theories and S-duality,” JHEP 04 (2018) 004, arXiv:1711.07941 [hep-th].  
[17] T. Arakawa, “Representation theory of W-algebras and Higgs branch conjecture,” in International Congress of Mathematicians (ICM 2018) Rio de Janeiro, Brazil, August 1-9, 2018. 2017. arXiv:1712.07331 [math.RT].  
[18] V. Niarchos, “Geometry of Higgs-branch superconformal primary bundles,” Phys. Rev. D98 no. 6, (2018) 065012, arXiv:1807.04296 [hep-th].

[19] B. Feigin and S. Gukov, “VOA[M4],” arXiv:1806.02470 [hep-th].  
[20] T. Creutzig, “Logarithmic W-algebras and Argyres-Douglas theories at higher rank,” arXiv:1809.01725 [hep-th].  
[21] P. C. Argyres, C. Long, and M. Martone, “The Singularity Structure of Scale-Invariant Rank-2 Coulomb Branches,” JHEP 05 (2018) 086, arXiv:1801.01122 [hep-th].  
[22] P. C. Argyres and M. Martone, “Scaling dimensions of Coulomb branch operators of 4d  $\mathbf{N} = 2$  superconformal field theories,” arXiv:1801.06554 [hep-th].  
[23] M. Caorsi and S. Cecotti, “Geometric classification of 4d  $\mathcal{N} = 2$  SCFTs,” JHEP 07 (2018) 138, arXiv:1801.04542 [hep-th].  
[24] M. Caorsi and S. Cecotti, “Special Arithmetic of Flavor,” JHEP 08 (2018) 057, arXiv:1803.00531 [hep-th].  
[25] C. Beem, C. Meneghelli, and L. Rastelli, “Free Field Realizations from the Higgs Branch,” arXiv:1903.07624 [hep-th].  
[26] D. Adamovic, “A realization of certain modules for the  $N = 4$  superconformal algebra and the affine Lie algebra  $A_2^{(1)}$ ,” arXiv:1407.1527 [math.QA].  
[27] M. Geck and G. Malle, “Reflection Groups. A Contribution to the Handbook of Algebra,” arXiv Mathematics e-prints (Nov, 2003) math/0311012, arXiv:math/0311012 [math.RT].  
[28] I. V. Dolgachev, “Reflection groups in algebraic geometry,” Bull. Amer. Math. Soc. (N.S.) 45 no. 1, (2008) 1-60.  
[29] P. C. Argyres and M. Martone, “Coulomb branches with complex singularities,” JHEP 06 (2018) 045, arXiv:1804.03152 [hep-th].  
[30] T. Nishinaka and Y. Tachikawa, “On 4d rank-one  $\mathcal{N} = 3$  superconformal field theories,” JHEP 09 (2016) 116, arXiv:1602.01503 [hep-th].  
[31] M. Lemos, P. Liendo, C. Meneghelli, and V. Mitev, “Bootstrapping  $\mathcal{N} = 3$  superconformal theories,” JHEP 04 (2017) 032, arXiv:1612.01536 [hep-th].  
[32] T. Bourton, A. Pini, and E. Pomoni, “4d  $\mathcal{N} = 3$  indices via discrete gauging,” arXiv:1804.05396 [hep-th].  
[33] F. Bonetti, C. Meneghelli, and L. Rastelli, "Bootstrapping  $\mathcal{N} = 4$  VOA," to appear.  
[34] T. Arakawa, “A remark on the c 2-cofiniteness condition on vertex algebras,” Mathematische Zeitschrift 270 no. 1, (Feb, 2012) 559-575, arXiv:1004.1492 [math].  
[35] O. Aharony and M. Evtikhiev, “On four dimensional  $\mathrm{N} = 3$  superconformal theories,” JHEP 04 (2016) 040, arXiv:1512.03524 [hep-th].  
[36] P. C. Argyres and J. R. Wittig, “Infinite coupling duals of  $N = 2$  gauge theories and new rank 1 superconformal field theories,” JHEP 01 (2008) 074, arXiv:0712.2028 [hep-th].  
[37] A. D. Shapere and Y. Tachikawa, “Central charges of  $\mathrm{N} = 2$  superconformal field theories in four dimensions,” JHEP 09 (2008) 109, arXiv:0804.1957 [hep-th].  
[38] I. García-Etxebarria and D. Regalado, “ $\mathcal{N} = 3$  four dimensional field theories,” JHEP 03 (2016) 083, arXiv:1512.06434 [hep-th].

[39] O. Aharony and Y. Tachikawa, “S-folds and 4d N=3 superconformal field theories,” JHEP 06 (2016) 044, arXiv:1602.08638 [hep-th].  
[40] I. García-Etxebarria and D. Regalado, “Exceptional  $\mathcal{N} = 3$  theories,” JHEP 12 (2017) 042, arXiv:1611.05769 [hep-th].  
[41] V. G. Kac and M. Wakimoto, “Quantum reduction and representation theory of superconformal algebras,” Advances in Mathematics 185 no. 2, (2004) 400 - 458.  
[42] K. Thielemans, An Algorithmic approach to operator product expansions, W algebras and W strings. PhD thesis, Leuven U., 1994. arXiv:hep-th/9506159 [hep-th].  
[43] R. P. Stanley, “Invariants of finite groups and their applications to combinatorics,” Bull. Amer. Math. Soc. (N.S.) 1 no. 3, (05, 1979) 475-511.  
[44] G. I. Lehrer and D. E. Taylor, Unitary Reflection Groups. Cambridge University Press, 2009.  
[45] A. Gadde, L. Rastelli, S. S. Razamat, and W. Yan, “Gauge Theories and Macdonald Polynomials,” Commun. Math. Phys. 319 (2013) 147–193, arXiv:1110.3740 [hep-th].  
[46] L. Rastelli and S. S. Razamat, “The Superconformal Index of Theories of Class  $S$ ," in New Dualities of Supersymmetric Gauge Theories, J. Teschner, ed., pp. 261-305. 2016. arXiv:1412.7131 [hep-th]. https://inspirehep.net/record/1335343/files/arXiv:1412.7131.pdf.  
[47] F. A. Dolan and H. Osborn, “On short and semi-short representations for four-dimensional superconformal symmetry,” Annals Phys. 307 (2003) 41–89, arXiv:hep-th/0209056 [hep-th].  
[48] F. Cachazo, M. R. Douglas, N. Seiberg, and E. Witten, “Chiral rings and anomalies in supersymmetric gauge theory,” JHEP 12 (2002) 071, arXiv:hep-th/0211170 [hep-th].  
[49] J. Kinney, J. M. Maldacena, S. Minwalla, and S. Raju, “An Index for 4 dimensional super conformal theories,” Commun. Math. Phys. 275 (2007) 209–254, arXiv:hep-th/0510251 [hep-th].  
[50] Y. Berest, G. Felder, S. Patotski, A. C. Ramadoss, and T. Willwacher, “Representation Homology, Lie Algebra Cohomology and Derived Harish-Chandra Homomorphism,” arXiv:1410.0043 [math].  
[51] R. W. Carter, “Conjugacy classes in the weyl group,” Compositio Mathematica 25 no. 1, (1972) 1-59.  
[52] G. Felder and M. Mueller-Lennert, “Analyticity of Nekrasov Partition Functions,” arXiv:1709.05232 [math-ph].