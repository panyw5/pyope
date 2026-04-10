# Testing material for $\mathcal{W}_{\mathbb{Z}_3}$ algebra

## 参考论文: VOAs labelled by complex reflection groups and 4d SCFTs

## $\mathcal{W}_{\mathbb{Z}_3}$ VOA 与 free field realization

Bosonic strong generators: $T, J, W, \bar W$

Fermionic strong generators: $G, \bar G, G_W, \bar G_{\bar W}$

> 在 `OPEdefs.m` 中声明算符时，请**按照**上述顺序，因为下面的 `OPE` 主要是按照这个顺序写的。

> $J, T, G, \bar G$ 形成 2d $\mathcal{N} = 2$ superconformal algebra

### Free field

这些 strong generators 可以用 free field realization 写成 $bc\beta\gamma$ 系统的正规乘积，具体如下。

$$
b(z)c(w)\sim \frac{1}{z-w},\qquad \beta(z)\gamma(w)\sim -\frac{1}{z-w},
$$

并且采用论文 5.1 的 $p=3$ 记号，那么生成元可写为

$$
\mathcal{J}=2:bc:+3:\beta\gamma:,
$$

$$
G=:\gamma b:,
$$

$$
\bar{G}=2:(\partial\beta)c:+3:\beta\,\partial c:,
$$

$$
\mathcal{T}
=-2:b\,\partial c:
-\frac{3}{2}:\beta\,\partial\gamma:
-:(\partial b)c:
-\frac{1}{2}:(\partial\beta)\gamma:,
$$

$$
W=\beta,
$$

$$
\begin{aligned}
\bar{W}
={}&:\beta^2\gamma^3:
+2:\beta\gamma^2bc:
-4:\beta(\partial\gamma)\gamma:
-\frac{4}{3}:\gamma\,b\,\partial c: \\
&+\frac{2}{3}:\gamma\,(\partial b)c:
+\frac{2}{3}:(\partial\beta)\gamma^2:
-\frac{8}{3}:(\partial\gamma)bc:
+\frac{10}{9}\partial^2\gamma.
\end{aligned}
$$

$G_W$ 和 $\bar G_{\bar W}$ 可以写成 OPE 极点

$$
G_W \coloneqq \{GW\}_1 = b,
$$


$$
\begin{aligned}
\bar{G}_{\bar W} \coloneqq \{\bar{G}{\bar W}\}_1 = & -\frac{10}{9}\partial^3 c
-\frac{8}{3}:b\,(\partial^2 c)\,c:
-3:\beta^2\gamma^2\partial c:
+4:\beta\gamma\,b\,(\partial c)\,c: \\
&+4:\beta\gamma\,\partial^2 c:
+4:\beta\,(\partial\gamma)\,\partial c:
-\frac{2}{3}:(\partial^2\beta)\gamma c:
+\frac{2}{3}:(\partial b)(\partial c)c: \\
&-2:(\partial\beta)\beta\gamma^2 c:
+\frac{8}{3}:(\partial\beta)(\partial\gamma)c: .
\end{aligned}
$$

```Mathematica
(* h[\[Beta]]=3/2, h[\[Gamma]]=-(1/2), h[b]=2, h[c]=-1 *)
p = 3;
J = p NO[\[Beta], \[Gamma]] + (p - 1) NO[b, c];
(* h[\[ScriptCapitalG]]=h[\[ScriptCapitalG]bar]=3/2 *)
G = \[ScriptCapitalG] = NO[b, \[Gamma]];
Gbar = \[ScriptCapitalG]bar = 
   p NO[\[Beta], Derivative[1][c]] + (p - 1) NO[Derivative[
      1][\[Beta]], c];
T = -(1/2) p NO[\[Beta], Derivative[1][\[Gamma]]] + (1 - 1/2 p) NO[
     Derivative[1][\[Beta]], \[Gamma]] - 
   1/2 (p + 1) NO[b, Derivative[1][c]] + 
   1/2 (1 - p) NO[Derivative[1][b], c];
(* h[W]=3/2 *)
W = \[Beta];
(* h[Wbar] = 3/2 *)
Wbar = NO[\[Beta], 
    NO[\[Beta], NO[\[Gamma], NO[\[Gamma], \[Gamma]]]]] + 
   2 NO[\[Beta], NO[\[Gamma], NO[\[Gamma], NO[b, c]]]] - 
   4 NO[\[Beta], NO[Derivative[1][\[Gamma]], \[Gamma]]] - 
   4/3 NO[\[Gamma], NO[b, Derivative[1][c]]] + 
   2/3 NO[\[Gamma], NO[Derivative[1][b], c]] + 
   2/3 NO[Derivative[1][\[Beta]], NO[\[Gamma], \[Gamma]]] - 
   8/3 NO[Derivative[1][\[Gamma]], NO[b, c]] + (
   10 \[Gamma]^\[Prime]\[Prime])/9;

(* 计算额外的生成元 \[ScriptCapitalG]W, h[\[ScriptCapitalG]W] = 2 *)
GW = OPEPole[1][G, W];
GbarWbar = OPEPole[1][Gbar, Wbar];
```

```python

J = 2 * NO(b, c) + 3 * NO(β, γ)

G = NO(γ, b)
Gbar = 2 * NO(d(β), c) + 3 * NO(β, d(c))

T = - 2 * NO(b, d(c)) - Fraction(3,2) * NO(β, d(γ)) - NO(d(b), c) - Fraction(1,2) * NO(d(β), γ)


W = β     # 权重 3/2
GW = b    # 权重 2

Wbar = (NO(β, NO(β, NO(γ, NO(γ, γ))))
      + 2 * NO(β, NO(γ, NO(γ, NO(b, c))))
      - 4 * NO(β, NO(d(γ), γ))
      - Fraction(4, 3) * NO(γ, NO(b, d(c)))
      + Fraction(2, 3) * NO(γ, NO(d(b), c))
      + Fraction(2, 3) * NO(d(β), NO(γ, γ))
      - Fraction(8, 3) * NO(d(γ), NO(b, c))
      + Fraction(10, 9) * d(d(γ)))

GbarWbar = (Fraction(8, 3) * NO(b, NO(d(d(c)), c))
       + 3 * NO(β, NO(β, NO(γ, NO(γ, d(c)))))
       - 4 * NO(β, NO(γ, NO(b, NO(d(c), c))))
       - 4 * NO(β, NO(γ, d(d(c))))
       - 4 * NO(β, NO(d(γ), d(c)))
       - Fraction(2, 3) * NO(d(b), NO(d(c), c))
       + 2 * NO(d(β), NO(β, NO(γ, NO(γ, c))))
       - Fraction(8, 3) * NO(d(β), NO(d(γ), c))
       + Fraction(2, 3) * NO(d(d(β)), NO(γ, c))
       + Fraction(10, 9) * d(d(d(c))))
```


### OPE

Nontrivial OPE 如下。

首先，各 strong generator 的conformal weight/$J$-荷：

- $h(J)=1$
- $h(G)=h(\bar G)=h(W)=h(\bar W)=\frac32$
- $h(G_W)=h(\bar G_{\bar W})=2$

所以有：

- $T(z)G(w)\sim \frac{3G/2}{(z-w)^2}+\frac{\partial G}{z-w}$
- $T(z)\bar G(w)\sim \frac{3\bar G/2}{(z-w)^2}+\frac{\partial \bar G}{z-w}$
- $T(z)W(w)\sim \frac{3W/2}{(z-w)^2}+\frac{\partial W}{z-w}$
- $T(z)\bar W(w)\sim \frac{3\bar W/2}{(z-w)^2}+\frac{\partial \bar W}{z-w}$
- $T(z)G_W(w)\sim \frac{2G_W}{(z-w)^2}+\frac{\partial G_W}{z-w}$
- $T(z)\bar G_{\bar W}(w)\sim \frac{2\bar G_{\bar W}}{(z-w)^2}+\frac{\partial \bar G_{\bar W}}{z-w}$

以及

- $J(z)G(w)\sim -\frac{G}{z-w}$
- $J(z)\bar G(w)\sim \frac{\bar G}{z-w}$
- $J(z)W(w)\sim \frac{3W}{z-w}$
- $J(z)\bar W(w)\sim -\frac{3\bar W}{z-w}$
- $J(z)G_W(w)\sim \frac{2G_W}{z-w}$
- $J(z)\bar G_{\bar W}(w)\sim -\frac{2\bar G_{\bar W}}{z-w}$

$W\bar W$ OPE：

$$
W(z)\bar W(w) = \frac{-20/9}{(z - w)^3} + \frac{\frac{4}{3} J(w)}{(z - w)^2} + \frac{\frac{2}{3}T(w) - \frac{1}{3} (JJ)(w) + \frac{2}{3}J'(w)}{(z - w)} .
$$

shortening / descendant 关系：

- $G(z)W(w)\sim \frac{G_W(w)}{z-w}$
- $\bar G(z) W(w) \sim 0$ ($W$ 是 chiral primary)
- $\bar G(z)\bar W(w)\sim \frac{\bar G_{\bar W}(w)}{z-w}$
- $G(z)\bar W(w)\sim 0$ ($\bar W$ 是 anti-chiral primary)

$$
G(z)\bar G(w)
= \frac{5}{(z-w)^3}
+ \frac{J(w)}{(z-w)^2}
+ \frac{-T(w) + \frac12 \partial J(w)}{z-w} .
$$

$$
\bar G(z) G_W(w)
= - \frac{3 W(w)}{(z-w)^2}
- \frac{\partial W(w)}{z-w} .
$$

$$
G(z) \bar G_{\bar W}(w)
= - \frac{3 \bar W(w)}{(z-w)^2}
- \frac{\partial \bar W(w)}{z-w} .
$$

$$
\bar W(z) G_W(w)
= \frac{\frac{4}{3} G(w)}{(z-w)^2}
+ \frac{\frac23 (JG)(w) + \frac23 \partial G(w)}{z-w} .
$$

$$
W(z) \bar G_{\bar W}(w)
= - \frac{\frac{4}{3}\bar G(w)}{(z-w)^2}
+ \frac{\frac23 (J\bar G)(w) - \frac23 \partial \bar G(w)}{z-w} .
$$

$$
G_W(z)\bar G_{\bar W}(w)
= -\frac{20/3}{(z-w)^4}
+ \frac{\frac83 J(w)}{(z-w)^3}
+ \frac{2T(w) - \frac13 (JJ)(w) + \frac{4}{3}\partial J(w)}{(z-w)^2}
+ \frac{\frac23 (G\bar G)(w) - \frac23 (JT)(w) - \frac13 (J\partial J)(w) + \frac43 \partial T(w) + \frac13 \partial^2 J(w)}{z-w} .
$$

其中二级、一级极点分别满足

$$
\{G_W\bar G_{\bar W}\}_2
= 2T - \frac13 (JJ) + \frac43 \partial J
$$

以及

$$
\{G_W\bar G_{\bar W}\}_1
= \frac23 (G\bar G)
- \frac23 (JT)
- \frac13 (J\partial J)
+ \frac43 \partial T
+ \frac13 \partial^2 J .
$$


## Jacobi identity


## null states

# Weight-8 T4 Null Family In Plain Operator Names

There are `19` basis directions containing `NO(T,NO(T,NO(T,T)))`.

## Basis 0

- term count: `11`


```wl
(
  -8/9 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + 4/3 * NO[T,NO[T,NO[T,T]]]
  + 4/9 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -8/9 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + 2/3 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + -3 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -1/9 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -4/9 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + -2/9 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + 5 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 1 * NO[T,NO[T,NO[Derivative[2][J],J]]]
) == 0
```

## Basis 1

- term count: `14`

```wl
(
  -13/72 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + -17/12 * NO[T,NO[T,NO[T,T]]]
  + -115/72 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + 1/8 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 17/18 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + -27/8 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + -11/6 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + 3/2 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + 9/8 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + -91/144 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -5/18 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + -8/9 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + -83/8 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 1 * NO[T,NO[T,Derivative[3][J]]]
) == 0
```

## Basis 2

- term count: `33`

```wl
(
  -24412/775929 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + -97648/258643 * NO[T,NO[T,NO[T,T]]]
  + -97648/775929 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + 919/3359 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 195296/775929 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + -48824/258643 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + 219708/258643 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -23773/70539 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + 97648/775929 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + 48824/775929 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + -366180/258643 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 56653/517286 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + 453463/517286 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + -1529631/517286 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -34494/23513 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + 1338642/258643 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + -302055/258643 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 189557/258643 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + -214317/23513 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + 1399581/517286 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + 302175/517286 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + 80460/258643 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + -350550/258643 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + 206880/258643 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + -164352/258643 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + 1441557/258643 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + -660579/258643 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + -7908/36949 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + -289332/258643 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -437571/517286 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + -67137/73898 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + 1104273/517286 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 1 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],NO[J,J]]]]
) == 0
```

## Basis 3

- term count: `49`

```wl
(
  1727237781242/15582441284487 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + -2735163527512/5194147094829 * NO[T,NO[T,NO[T,T]]]
  + -2745610635568/15582441284487 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + 427883065/22485485259 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 859465187720/15582441284487 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + -603636/743789 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + -1357134655700/5194147094829 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + 682459387714/577127454981 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + 136708/743789 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + -327227706283/1416585571317 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + 2837646677974/15582441284487 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + -3838765262740/15582441284487 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + -4897527467198/1731382364943 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + -1911545190/192375818327 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + -213595793165/1731382364943 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + 51628790568/192375818327 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + 48872153140/157398396813 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + -189947456548/577127454981 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + 131684723141/1731382364943 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + -382204416106/1731382364943 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + 44141070367/52466132271 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + -224120893559/577127454981 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + -102358717871/577127454981 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + 3992757100/192375818327 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + -15683827436/577127454981 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + 29973055132/1731382364943 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + 712329374116/1731382364943 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + -80844929896/192375818327 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + 445465614452/1731382364943 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + 10827805310/247340337849 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + 46774809041/192375818327 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + 748044/743789 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + 210540/743789 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + -1158696/743789 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + 46464/743789 * NO[T,NO[G,Derivative[3][Gbar]]]
  + -250008/743789 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + 138204/743789 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + -1782/743789 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + 50820/743789 * NO[T,NO[W,Derivative[3][Wbar]]]
  + -37620/743789 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + 201732966187/3462764729886 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + 24644884247/247340337849 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + 13642971319/192375818327 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + -56894/743789 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + -368107/743789 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + -358132/2231367 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + -444236/6694101 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + 83626/6694101 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,NO[Derivative[1][J],NO[Derivative[1][W],Wbar]]]
) == 0
```

## Basis 6

- term count: `49`

```wl
(
  13351306326725/41553176758632 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + -9690709968437/6925529459772 * NO[T,NO[T,NO[T,T]]]
  + -27833346838369/41553176758632 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + 8516950673/119922588048 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + -1159631320214/5194147094829 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + -11326497/5950312 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + -1175386949207/3462764729886 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + 717146462015/192375818327 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + 1622127/5950312 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + -1820575257283/3777561523512 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + 17790305360977/41553176758632 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + -9696179099521/10388294189658 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + -28977335433521/4617019639848 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 33535889259/6156026186464 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + -6615328461391/55404235678176 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + -909488322213/6156026186464 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + 812372863363/1259187174504 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + 19214234058/192375818327 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + 334595352001/9234039279696 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + -13550965950649/27702117839088 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + -133972788217/279819372112 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + -5050119440661/6156026186464 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + -31394690698663/55404235678176 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + -97873473443/769503273308 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + 5722322755/1539006546616 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + -5840540831/1154254909962 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + 5052445115513/3462764729886 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + 1055923555227/3078013093232 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + 2241702792037/9234039279696 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + 50053093201/329787117132 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + -260590313679/1539006546616 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + 4760109/2975156 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + 5704049/2975156 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + -5262405/1487578 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + 461938/2231367 * NO[T,NO[G,Derivative[3][Gbar]]]
  + -2271217/1487578 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + 4857873/2975156 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + 104895/5950312 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + -8147/2975156 * NO[T,NO[W,Derivative[3][Wbar]]]
  + -380353/2975156 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + -590341596919/18468078559392 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + 3828162157057/7914890811168 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + 771193488651/6156026186464 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + -1807791/2975156 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + 101027/743789 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + -1309867/2975156 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + 120881/8925468 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + -1845949/8925468 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,NO[Derivative[1][Gbar],NO[Wbar,GW]]]
) == 0
```

## Basis 7

- term count: `49`

```wl
(
  1073538361195/15582441284487 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + -5090622458030/5194147094829 * NO[T,NO[T,NO[T,T]]]
  + 524733038875/15582441284487 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + 20834755/22485485259 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 9093516606700/15582441284487 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + -207531/743789 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + -8160666725920/5194147094829 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + 460647855830/577127454981 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -305295/743789 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + -402375259555/2833171142634 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + 5523241671980/15582441284487 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + 9884020553500/15582441284487 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + -6482592286903/1731382364943 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 25182105065/1154254909962 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + 225513962965/1154254909962 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + -228015443631/384751636654 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -4849505530/52466132271 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + 24008228161/192375818327 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + 47259950770/192375818327 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + -18228220275/192375818327 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + -36167301257/17488710757 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + -284271002111/384751636654 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + 224466337565/1154254909962 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + -162042972640/192375818327 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + 3901476641/192375818327 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + 132251174864/192375818327 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + -173742915260/192375818327 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + -215879296707/192375818327 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + -206254009723/192375818327 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + 2074389550/27482259761 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + -292049870506/192375818327 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -3383370/743789 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + 436878/743789 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + 5064332/743789 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + 330280/2231367 * NO[T,NO[G,Derivative[3][Gbar]]]
  + 1154400/743789 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + -1109400/743789 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + 71847/743789 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + 182397/743789 * NO[T,NO[W,Derivative[3][Wbar]]]
  + 29192/743789 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + -15817643360/192375818327 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + -44639648275/164893558566 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + -26936064585/384751636654 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 197730/743789 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + -203418/743789 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + 710968/743789 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + -68480/743789 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + -156750/743789 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,NO[Derivative[1][W],Derivative[2][Wbar]]]
) == 0
```

## Basis 8

- term count: `49`

```wl
(
  7008328421951/62329765137948 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + -1050598621745/10388294189658 * NO[T,NO[T,NO[T,T]]]
  + -5643751528375/62329765137948 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -170325478/22485485259 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 17071727660891/15582441284487 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + -519417/2975156 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + 622988915785/5194147094829 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + 400371904045/577127454981 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -4505145/2975156 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + -962827479563/11332684570536 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + 6958228726583/31164882568974 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + 22937998178090/15582441284487 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + -10405377597719/6925529459772 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 69539785111/4617019639848 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + 5305842884887/13851058919544 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + -631255428123/1539006546616 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -289863119551/314796793626 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + -4037052508/192375818327 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + 179494775497/769503273308 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 3552053918449/6925529459772 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + -101022901915/69954843028 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + 2199678895149/1539006546616 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + 8975549323711/13851058919544 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + -117549842609/192375818327 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + 1092378845/384751636654 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + 56467142407/192375818327 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + -2566673454250/1731382364943 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + -288964040003/769503273308 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + -319467347055/769503273308 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + -13565541497/27482259761 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + -442354279987/384751636654 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -2788335/743789 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + -1627014/743789 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + 5630220/743789 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + -348806/743789 * NO[T,NO[G,Derivative[3][Gbar]]]
  + 955530/743789 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + -999465/743789 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + 140859/1487578 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + 222822/743789 * NO[T,NO[W,Derivative[3][Wbar]]]
  + -733/743789 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + -87097848811/1539006546616 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + -1605145917769/1978722702792 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + 246848494545/1539006546616 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 930071/743789 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + -73623/743789 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + 706165/743789 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + -2561/743789 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + -32907/743789 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,NO[Derivative[1][GW],Derivative[1][GbarWbar]]]
) == 0
```

## Basis 9

- term count: `33`


```wl
(
  5836/258643 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + 70032/258643 * NO[T,NO[T,NO[T,T]]]
  + 23344/258643 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -1308/3359 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + -46688/258643 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + 35016/258643 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + -157572/258643 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -3020/23513 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -23344/258643 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + -11672/258643 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + 262620/258643 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 5854/258643 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + -121734/258643 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + -158058/258643 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + 13458/23513 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + -55719/258643 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + -15144/258643 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + -193542/258643 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + 19863/23513 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + -333486/258643 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + -44598/258643 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + -126360/258643 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + -136872/258643 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + 133368/258643 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + 133128/258643 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + 79488/258643 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + -201984/258643 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + 3492/36949 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + 173178/258643 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -62595/258643 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + 7338/36949 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + 23382/258643 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 1 * NO[T,NO[Derivative[2][J],NO[J,NO[J,J]]]]
) == 0
```

## Basis 12

- term count: `33`

```wl
(
  -66/3359 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + -792/3359 * NO[T,NO[T,NO[T,T]]]
  + -264/3359 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + 621/3359 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 528/3359 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + -396/3359 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + 1782/3359 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -2463/3359 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + 264/3359 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + 132/3359 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + -2970/3359 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 141/3359 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + 762/3359 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + -3807/3359 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -576/3359 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + 9198/3359 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + -2847/3359 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 1284/3359 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + -11619/3359 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + 4248/3359 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + 1050/3359 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + -2052/3359 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + -2016/3359 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + 2028/3359 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + -732/3359 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + 7182/3359 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + -5106/3359 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + -1050/3359 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + -3699/3359 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -2033/6718 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + -3834/3359 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + -3186/3359 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 1 * NO[T,NO[Derivative[2][J],NO[Derivative[1][J],J]]]
) == 0
```

## Basis 13

- term count: `33`

```wl
(
  -3611176/2327787 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + 4177592/775929 * NO[T,NO[T,NO[T,T]]]
  + 4177592/2327787 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -181/3359 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + -8355184/2327787 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + 2088796/775929 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + -3133194/258643 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + 221939/211617 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -4177592/2327787 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + -2088796/2327787 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + 5221990/258643 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + -7481/517286 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + 125877/517286 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + 201987/517286 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -18492/23513 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + -136377/258643 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + -82311/258643 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 205713/258643 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + 35730/23513 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + 2718171/517286 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + 974217/517286 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + -218286/258643 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + 477144/258643 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + -245334/258643 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + 68160/258643 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + 1852317/258643 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + -444867/258643 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + -5526/36949 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + -64926/258643 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -729/517286 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + -182067/73898 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + 1067607/517286 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 1 * NO[T,NO[Derivative[2][J],Derivative[2][J]]]
) == 0
```

## Basis 14

- term count: `49`

```wl
(
  -118436536077119/166212707034528 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + 52667617453085/27702117839088 * NO[T,NO[T,NO[T,T]]]
  + 273089400164311/166212707034528 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -34499241989/239845176096 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 29088200602519/41553176758632 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + 83109519/23801248 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + -13189195666777/13851058919544 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + -11244715584037/1539006546616 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -23958261/23801248 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + 21711417774533/30220492188096 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -31777803019549/41553176758632 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + 18310217432347/10388294189658 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + 199047684498287/18468078559392 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + -27848131649/1539006546616 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + 2013546273257/13851058919544 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + 189093972231/384751636654 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -273530304245/314796793626 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + -775752811653/1539006546616 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + 434623594107/1539006546616 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 995782335518/1731382364943 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + 106108605951/69954843028 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + 1662046371047/1539006546616 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + 14452461540617/13851058919544 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + -101623543872/192375818327 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + -47211683721/1539006546616 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + -79331498661/769503273308 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + -3858511795814/1731382364943 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + -1385200306131/769503273308 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + 161286579975/769503273308 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + -84594114263/109929039044 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + 695627360121/769503273308 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -9026847/2975156 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + -12559941/2975156 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + 8831523/1487578 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + -637475/2231367 * NO[T,NO[G,Derivative[3][Gbar]]]
  + 1846287/743789 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + -8655975/2975156 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + -467289/5950312 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + -15453/1487578 * NO[T,NO[W,Derivative[3][Wbar]]]
  + 880803/1487578 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + 395394081903/3078013093232 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + -1786269956177/1978722702792 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + -514769961357/769503273308 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 3285493/1487578 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + -9808659/5950312 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + 13526199/5950312 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + -1921821/5950312 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + -650523/743789 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,NO[Derivative[2][G],Derivative[1][Gbar]]]
) == 0
```

## Basis 15

- term count: `49`

```wl
(
  -121991554606237/373978590827688 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + 126718290392755/62329765137948 * NO[T,NO[T,NO[T,T]]]
  + 175301255979125/373978590827688 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -4648543547/67456455777 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + -30190936604197/93494647706922 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + 13553661/5950312 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + 51213403799785/31164882568974 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + -13364398690025/3462764729886 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + 1410745/5950312 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + 42758251486801/67996107423216 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -119214571231201/186989295413844 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + 2915631474835/46747323853461 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + 287509447830565/41553176758632 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + -106880814077/27702117839088 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + 8638343477251/83106353517264 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + 332201715187/3078013093232 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -1115088571663/1888780761756 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + 52582786968/192375818327 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + -1787126024357/13851058919544 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 22411289923717/41553176758632 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + 196575346909/419729058168 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + 11735298923329/9234039279696 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + 18682915762843/83106353517264 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + 347991493963/1154254909962 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + 9912471913/769503273308 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + -1087573356095/1731382364943 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + -1752069048845/5194147094829 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + 850658652385/4617019639848 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + 10169635232935/13851058919544 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + -69384414653/494680675698 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + 2563123543091/2308509819924 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + 839115/1487578 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + -1517579/1487578 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + 51777/743789 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + -1640482/6694101 * NO[T,NO[G,Derivative[3][Gbar]]]
  + -167785/743789 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + 1285085/1487578 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + -301671/2975156 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + -161129/1487578 * NO[T,NO[W,Derivative[3][Wbar]]]
  + -209149/1487578 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + -145249668829/27702117839088 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + -2503169379277/11872336216752 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + -414885086477/3078013093232 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 191923/4462734 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + 477250/743789 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + -1193237/4462734 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + 1841471/13388202 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + 3889157/13388202 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,NO[Derivative[2][W],Derivative[1][Wbar]]]
) == 0
```

## Basis 16

- term count: `49`


```wl
(
  -1005733874243/1978722702792 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + 843989657117/329787117132 * NO[T,NO[T,NO[T,T]]]
  + 1869878243611/1978722702792 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -1107741614/7495161753 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 179998745221/494680675698 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + 22853223/5950312 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + 165522681935/164893558566 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + -353669092401/54964519522 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + 2628507/5950312 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + 327195514799/359767764144 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -650426630471/989361351396 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + 177259949105/247340337849 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + 2370693103883/219858078088 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 989064743/439716156176 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + 1315773661093/3957445405584 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + -24780000753/439716156176 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -92144409529/89941941036 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + 2259622954/27482259761 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + 219261342533/659574234264 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 1516217173939/1978722702792 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + -4138987181/19987098008 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + 514432656087/439716156176 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + 2102520656077/3957445405584 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + -39460550571/54964519522 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + -255154411/109929039044 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + -37727607397/82446779283 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + -391779057311/247340337849 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + -209382188457/219858078088 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + 547595383121/659574234264 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + -22836129433/164893558566 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + 112574612505/109929039044 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -6527871/1487578 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + -1541985/1487578 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + 5320329/743789 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + -561746/2231367 * NO[T,NO[G,Derivative[3][Gbar]]]
  + 1055421/743789 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + 246879/1487578 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + -351621/2975156 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + 551121/1487578 * NO[T,NO[W,Derivative[3][Wbar]]]
  + -736399/1487578 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + -7125671699/1319148468528 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + -3462918052621/3957445405584 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + -145824238497/439716156176 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 641479/1487578 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + -79407/743789 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + 934371/1487578 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + 131239/4462734 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + 2175409/4462734 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,NO[Derivative[2][GW],GbarWbar]]
) == 0
```

## Basis 17

- term count: `33`

```wl
(
  -13152/258643 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + -157824/258643 * NO[T,NO[T,NO[T,T]]]
  + -52608/258643 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -927/3359 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 105216/258643 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + -78912/258643 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + 355104/258643 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -6393/23513 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + 52608/258643 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + 26304/258643 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + -591840/258643 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 57111/517286 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + 740667/517286 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + -1541997/517286 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -77226/23513 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + -107370/258643 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + 344181/258643 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 163341/258643 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + -231723/23513 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + 1772721/517286 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + 239835/517286 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + 514512/258643 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + 211131/258643 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + -312258/258643 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + -351072/258643 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + 285147/258643 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + 388713/258643 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + -576/36949 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + -275400/258643 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -266181/517286 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + 12435/73898 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + 1277883/517286 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 1 * NO[T,NO[Derivative[3][J],NO[J,J]]]
) == 0
```

## Basis 18

- term count: `49`


```wl
(
  80798229882/192375818327 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + -345527784936/192375818327 * NO[T,NO[T,NO[T,T]]]
  + -115650796860/192375818327 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + 379010301/2498387251 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 20767226292/192375818327 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + -2222478/743789 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + -171339286824/192375818327 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + 775803409632/192375818327 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + 503334/743789 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + -21191831841/17488710757 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -37564143480/192375818327 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + -179064173412/192375818327 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + -1900599989886/192375818327 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + -7578358353/192375818327 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + -104239799127/192375818327 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + 204678525780/192375818327 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + 23981620644/17488710757 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + -75307553289/192375818327 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + -12496018449/192375818327 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + -169410457548/192375818327 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + 60975251262/17488710757 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + -358947021177/192375818327 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + -135119065071/192375818327 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + 33786096192/192375818327 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + -115076823705/192375818327 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + 69785534226/192375818327 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + 328339840704/192375818327 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + -88497321390/192375818327 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + 61354263030/192375818327 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + 7976415990/27482259761 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + 274399664202/192375818327 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + 2754162/743789 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + 775170/743789 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + -4266108/743789 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + 171072/743789 * NO[T,NO[G,Derivative[3][Gbar]]]
  + -920484/743789 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + 508842/743789 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + -6561/743789 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + 187110/743789 * NO[T,NO[W,Derivative[3][Wbar]]]
  + -138510/743789 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + 83947184775/384751636654 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + 19441264131/27482259761 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + -186573187296/192375818327 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + -818028/743789 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + -2268135/743789 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + 169029/743789 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + -587436/743789 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + -574344/743789 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,NO[Derivative[3][J],Derivative[1][J]]]
) == 0
```

## Basis 19

- term count: `49`

```wl
(
  -118436536077119/166212707034528 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + 52667617453085/27702117839088 * NO[T,NO[T,NO[T,T]]]
  + 273089400164311/166212707034528 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -34499241989/239845176096 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + -81720270753833/41553176758632 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + 83109519/23801248 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + -13189195666777/13851058919544 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + -11244715584037/1539006546616 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -23958261/23801248 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + 21711417774533/30220492188096 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -31777803019549/41553176758632 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + -44019547705601/10388294189658 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + 199047684498287/18468078559392 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + -27848131649/1539006546616 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + -1144278455519/1539006546616 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + 189093972231/384751636654 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + 78426388683/34977421514 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + -775752811653/1539006546616 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + 434623594107/1539006546616 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + -295484245855/192375818327 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + 106108605951/69954843028 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + -6032986362033/1539006546616 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + -1814185487967/1539006546616 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + -101623543872/192375818327 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + -47211683721/1539006546616 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + -79331498661/769503273308 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + 554530649692/192375818327 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + -1385200306131/769503273308 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + 161286579975/769503273308 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + 25334924781/109929039044 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + 695627360121/769503273308 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + 17749557/2975156 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + 14216463/2975156 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + -22407615/1487578 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + 779227/743789 * NO[T,NO[G,Derivative[3][Gbar]]]
  + 1846287/743789 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + -8655975/2975156 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + -467289/5950312 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + -15453/1487578 * NO[T,NO[W,Derivative[3][Wbar]]]
  + 880803/1487578 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + 395394081903/3078013093232 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + 290099067287/219858078088 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + -514769961357/769503273308 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + -2664819/1487578 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + -9808659/5950312 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + 13526199/5950312 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + -1921821/5950312 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + -650523/743789 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,NO[Derivative[3][G],Gbar]]
) == 0
```

## Basis 20

- term count: `49`


```wl
(
  -8177164555237/41553176758632 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + 8429300171707/6925529459772 * NO[T,NO[T,NO[T,T]]]
  + 11475651100589/41553176758632 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -319316338/7495161753 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + -3086968484245/10388294189658 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + 2231005/5950312 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + 3453062353633/3462764729886 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + -673512111449/384751636654 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + 2817849/5950312 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + 2640855435865/7555123047024 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -8035629489841/20776588379316 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + -1242671551529/5194147094829 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + 13753860572989/4617019639848 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + -31090516917/3078013093232 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + -822458774495/27702117839088 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + 818662508995/3078013093232 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -102442539253/629593587252 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + -252588237064/577127454981 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + -241495808863/4617019639848 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 2985074236999/13851058919544 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + 148333380967/139909686056 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + 8773812882145/9234039279696 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + -15238136279/27702117839088 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + 95989262205/384751636654 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + 898000278547/2308509819924 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + -273712559485/577127454981 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + 179663025313/1731382364943 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + 2002562513579/1539006546616 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + 638569527551/1539006546616 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + -2877446467/164893558566 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + -3096464543/769503273308 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + 2198055/1487578 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + -799627/1487578 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + -1092837/743789 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + -192890/2231367 * NO[T,NO[G,Derivative[3][Gbar]]]
  + -316125/743789 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + 1282581/1487578 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + 542349/2975156 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + -295605/1487578 * NO[T,NO[W,Derivative[3][Wbar]]]
  + -225517/1487578 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + 341989319305/9234039279696 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + -137762466415/3957445405584 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + 2163895151315/3078013093232 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + -211319/1487578 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + 282340/743789 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + -526239/1487578 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + 250141/4462734 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + 1051015/4462734 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,NO[Derivative[3][W],Wbar]]
) == 0
```

## Basis 21

- term count: `33`


```wl
(
  -52880/258643 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + -634560/258643 * NO[T,NO[T,NO[T,T]]]
  + -211520/258643 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + 870/3359 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 423040/258643 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + -317280/258643 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + 1427760/258643 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -55310/23513 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + 211520/258643 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + 105760/258643 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + -2379600/258643 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 26553/258643 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + 254775/258643 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + -716931/258643 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -43620/23513 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + 989892/258643 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + -444594/258643 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 897450/258643 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + -225522/23513 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + 540765/258643 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + -412245/258643 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + 1049220/258643 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + -853056/258643 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + 218964/258643 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + -1185000/258643 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + -1495638/258643 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + 629322/258643 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + -34680/36949 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + -1617030/258643 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -126855/258643 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + 39495/36949 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + -791019/258643 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 1 * NO[T,NO[Derivative[4][J],J]]
) == 0
```

## Basis 22

- term count: `49`

```wl
(
  -5378343119595/769503273308 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + 6415701130935/384751636654 * NO[T,NO[T,NO[T,T]]]
  + 7188625316475/769503273308 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -963161535/2498387251 * NO[T,NO[T,NO[J,NO[J,NO[J,J]]]]]
  + 714112646985/192375818327 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + 20590065/2975156 * NO[T,NO[T,NO[J,NO[W,Wbar]]]]
  + -579693139155/192375818327 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + -10441763701785/192375818327 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + 6066765/2975156 * NO[T,NO[T,NO[GW,GbarWbar]]]
  + 474116377335/139909686056 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -3751917416895/384751636654 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + 522228822030/192375818327 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + 42093509740155/769503273308 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + -323319044745/1539006546616 * NO[T,NO[J,NO[J,NO[J,NO[J,NO[J,J]]]]]]
  + -296700091155/1539006546616 * NO[T,NO[J,NO[J,NO[J,NO[G,Gbar]]]]]
  + 8722770514335/1539006546616 * NO[T,NO[J,NO[J,NO[J,NO[W,Wbar]]]]]
  + -144694956345/34977421514 * NO[T,NO[J,NO[J,NO[G,Derivative[1][Gbar]]]]]
  + -979549523325/192375818327 * NO[T,NO[J,NO[J,NO[W,Derivative[1][Wbar]]]]]
  + -574687388685/769503273308 * NO[T,NO[J,NO[J,NO[GW,GbarWbar]]]]
  + 4201614334035/769503273308 * NO[T,NO[J,NO[J,NO[Derivative[1][G],Gbar]]]]
  + 1319185091115/69954843028 * NO[T,NO[J,NO[J,NO[Derivative[1][W],Wbar]]]]
  + -5439403158345/1539006546616 * NO[T,NO[J,NO[G,NO[W,GbarWbar]]]]
  + 6885323229045/1539006546616 * NO[T,NO[J,NO[G,Derivative[2][Gbar]]]]
  + -46985316435/192375818327 * NO[T,NO[J,NO[Gbar,NO[Wbar,GW]]]]
  + 1103649978645/384751636654 * NO[T,NO[J,NO[W,Derivative[2][Wbar]]]]
  + -411153361620/192375818327 * NO[T,NO[J,NO[GW,Derivative[1][GbarWbar]]]]
  + -2715585007230/192375818327 * NO[T,NO[J,NO[Derivative[1][G],Derivative[1][Gbar]]]]
  + 1623684402855/769503273308 * NO[T,NO[J,NO[Derivative[1][W],Derivative[1][Wbar]]]]
  + -3879877604205/769503273308 * NO[T,NO[J,NO[Derivative[1][GW],GbarWbar]]]
  + 36342047265/27482259761 * NO[T,NO[J,NO[Derivative[2][G],Gbar]]]
  + 1246676059665/384751636654 * NO[T,NO[J,NO[Derivative[2][W],Wbar]]]
  + -37487205/743789 * NO[T,NO[G,NO[Gbar,NO[W,Wbar]]]]
  + -10550925/743789 * NO[T,NO[G,NO[W,Derivative[1][GbarWbar]]]]
  + 58066470/743789 * NO[T,NO[G,NO[Derivative[1][W],GbarWbar]]]
  + -2328480/743789 * NO[T,NO[G,Derivative[3][Gbar]]]
  + 12528810/743789 * NO[T,NO[Gbar,NO[Wbar,Derivative[1][GW]]]]
  + -6925905/743789 * NO[T,NO[Gbar,NO[Derivative[1][Wbar],GW]]]
  + 178605/1487578 * NO[T,NO[W,NO[W,NO[Wbar,Wbar]]]]
  + -2546775/743789 * NO[T,NO[W,Derivative[3][Wbar]]]
  + 1885275/743789 * NO[T,NO[GW,Derivative[2][GbarWbar]]]
  + 1768294248555/1539006546616 * NO[T,NO[Derivative[1][J],NO[J,NO[J,NO[J,J]]]]]
  + -721437009075/219858078088 * NO[T,NO[Derivative[1][J],NO[J,NO[G,Gbar]]]]
  + -771150832785/1539006546616 * NO[T,NO[Derivative[1][J],NO[J,NO[W,Wbar]]]]
  + 7415325/743789 * NO[T,NO[Derivative[1][J],NO[G,Derivative[1][Gbar]]]]
  + -8177085/743789 * NO[T,NO[Derivative[1][J],NO[W,Derivative[1][Wbar]]]]
  + 6996690/743789 * NO[T,NO[Derivative[1][J],NO[GW,GbarWbar]]]
  + -1921530/743789 * NO[T,NO[Derivative[1][J],NO[Derivative[1][J],Derivative[1][J]]]]
  + -3339375/743789 * NO[T,NO[Derivative[1][J],NO[Derivative[1][G],Gbar]]]
  + 1 * NO[T,Derivative[5][J]]
) == 0
```

## Particular T4 Relation


```wl
(
  -2/3 * NO[T,NO[T,NO[T,NO[J,J]]]]
  + 1 * NO[T,NO[T,NO[T,T]]]
  + 1/3 * NO[T,NO[T,NO[T,Derivative[1][J]]]]
  + -2/3 * NO[T,NO[T,NO[J,NO[G,Gbar]]]]
  + 1/2 * NO[T,NO[T,NO[G,Derivative[1][Gbar]]]]
  + -9/4 * NO[T,NO[T,NO[W,Derivative[1][Wbar]]]]
  + -1/12 * NO[T,NO[T,NO[Derivative[1][J],NO[J,J]]]]
  + -1/3 * NO[T,NO[T,NO[Derivative[1][J],Derivative[1][J]]]]
  + -1/6 * NO[T,NO[T,NO[Derivative[1][G],Gbar]]]
  + 15/4 * NO[T,NO[T,NO[Derivative[1][W],Wbar]]]
  + 3/4 * NO[T,NO[T,NO[Derivative[2][J],J]]]
) == 0
```

