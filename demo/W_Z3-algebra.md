##
Ref: VOAs labelled by complex reflection groups and 4d SCFTs

## $\mathcal{W}_{\mathbb{Z}_3}$ VOA 与 free field realization

strong generators: $J, G, \bar G, T, W, \bar W, G_W, \bar G_{\bar W}$.

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

以及反手征生成元

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
G_W=\{GW\}_1 = b,
$$

而

$$
\bar{G}_{\bar W}=\{\bar{G}{\bar W}\}_1
$$

是由上面的自由场表达式直接取一级极点得到，具体为

$$
\begin{aligned}
\bar G_{\bar W} ={}& -\frac{10}{9}\partial^3 c
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
- $J(z)\bar G_{\bar W}(w)\sim -\frac{\bar G_{\bar W}}{z-w}$

shortening / descendant 关系：

- $G(z)W(w)\sim \frac{G_W(w)}{z-w}$
- $\bar G(z)\bar W(w)\sim \frac{\bar G_{\bar W}(w)}{z-w}$

$W \bar W$

$$
W(z)\bar W(w) = \frac{-20/9}{(z - w)^3} + \frac{\frac{4}{3} J(w)}{(z - w)^2} + \frac{\frac{2}{3}T(w) - \frac{1}{3} (JJ)(w) + \frac{2}{3}J'(w)}{(z - w)} .
$$

并且

- $G(z)\bar W(w)$ non-singular
- $\bar G(z)W(w)$ non-singular
- $\bar W(z)\bar W(w)$ non-singular

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

这里的系数对应于本文档前面采用的归一化

$$
J = 2:bc: + 3:\beta\gamma: .
$$