# Fourier decay on odd superellipsoids

This repository contains a corrected proof for the hypersurface

\[
\Sigma_k=\left\{x\in\mathbb R^d:\sum_{j=1}^d|x_j|^k=1\right\},
\qquad k\ge3\ \text{odd}.
\]

The measure is the **Gelfand–Leray measure**

\[
d\nu=\frac{d\mathcal H^{d-1}}{|\nabla P|},
\qquad P(x)=\sum_{j=1}^d|x_j|^k,
\]

and the Fourier transform convention is

\[
\widehat\nu(\xi)=\int_{\Sigma_k}e^{-i\langle\xi,x\rangle}\,d\nu(x).
\]

The proved estimate is

\[
|\widehat\nu(\xi)|
\le C_{d,k}(1+|\xi|)^{-\min\{k,(d-1)/k\}}.
\]

For \(|\xi|\ge1\), this is equivalent to

\[
|\widehat\nu(\xi)|
\le C_{d,k}\max\left\{|\xi|^{-k},|\xi|^{-(d-1)/k}\right\}.
\]

In particular, when \(d\le k^2+1\), the decay is
\((1+|\xi|)^{-(d-1)/k}\).

## Corrections made in the proof

- The measure is defined explicitly. In a graph chart
  \(x_1=F(x')\), the Gelfand–Leray density is
  \(dx'/(k|F(x')|^{k-1})\); it is not the Euclidean graph Jacobian.
- The incomplete \(I_{20}\) lemma and the duplicated \(I_2\) sections have
  been removed.
- The vertex estimate is based on a stated multivariable finite-type lemma.
  A one-dimensional radial van der Corput estimate alone does not imply one
  factor of \(|\xi|^{-1/k}\) for every tangential direction.
- Integration by parts across coordinate hyperplanes is performed one-sidedly.
  Jets of orders below \(k\) cancel, while the first possible jump contributes
  \(O(|\xi|^{-k})\).
- The estimate is written with \(1+|\xi|\), so it is meaningful also at
  \(\xi=0\).

The full statement and proof are in [`estimates_for_nu.tex`](estimates_for_nu.tex).
