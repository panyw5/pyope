# convention

- Normal order product: $(ab)(w)$ means the regular part of the OPE $a(z)b(w)$ as $z \to w$, or equivalently
  $$
  (ab)(w) = \oint \frac{dz}{2\pi i} \frac{1}{z-w} a(z) b(w)
  $$

# $bc$ system

Fermionic $b, c$ fields, with OPE
$$
b(z) c(w) \sim \frac{1}{z - w}
$$
$U(1)$ current
$$
J(z) = (cb)(z)
$$
OPE
$$
J(z)b(w) = - \frac{b(w)}{z - w}, \qquad
J(z)c(w) = \frac{c(w)}{z - w}, \qquad
J(z)J(w) = + \frac{1}{z - w} \ .
$$
The stress tensor
$$
T_\lambda = (\lambda - 1)(b \partial c) - \lambda c\partial b
= - b \partial c - \lambda \partial J
= T_{\lambda = 0} - \lambda \partial J \ .
$$
Under $T_\lambda$ the conformal weights of $b, c$ are
$$
h_b = 1 - \lambda, \qquad h_c = \lambda
$$
The central charge is $c_\lambda = -2(6 \lambda(\lambda - 1) + 1)$


# $\beta \gamma$ system

Bosonic $\beta, \gamma$ fields, with OPE
$$
\beta(z) \gamma(w) \sim \frac{1}{z - w}
$$
> Some literature uses the opposite sign convention for the OPE.
The $U(1)$ current
$$
J(z) = (\gamma \beta)(z)
$$
with OPE
$$
J(z) \beta(w) = - \frac{\beta(w)}{z - w}, \qquad
J(z) \gamma(w) = \frac{\gamma(w)}{z - w}, \qquad
J(z)J(w) = - \frac{1}{z - w} \ .
$$
The stress tensor
$$
T_\lambda = (1 - \lambda) (\beta \partial \gamma) - \lambda (\gamma \partial \beta)
= + \beta \partial \gamma - \lambda \partial J
$$
Under $T_\lambda$ the conformal weights of $\beta, \gamma$ are
$$
h_\beta = 1 - \lambda, \qquad h_\gamma = \lambda
$$
