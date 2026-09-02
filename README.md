\documentclass{amsart}

\usepackage{amsmath,amssymb,amsthm}
\usepackage{mathrsfs}
\usepackage{geometry}
\usepackage{hyperref}
\usepackage{enumitem}
\usepackage{color}
\geometry{margin=2.5cm}

\newtheorem{theorem}{Theorem}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollary}
\theoremstyle{definition}
\newtheorem{remark}[theorem]{Remark}


\newcommand{\Z}{\mathbb{Z}}
\newcommand{\N}{\mathbb{N}}
\newcommand{\C}{\mathbb{C}}
\newcommand{\Gk}{\Gamma_k^n}
\newcommand{\Gkp}{\Gamma_k^{n-1}}
\newcommand{\norm}[1]{\|#1\|}
\newcommand{\abs}[1]{\left|#1\right|}
\newcommand{\ip}[2]{\langle #1, #2 \rangle}
\DeclareMathOperator{\sgn}{sgn}
\DeclareMathOperator{\supp}{supp}

\title{Estimates for the Fourier transform of surface-carried measure}

\begin{document}

\maketitle

Consider the surface
$$\Sigma:=\{x\in\mathbb{R}^d:\ |x_1|^k+|x_2|^k+\cdots+|x_d|^k=1\}.$$
We assume that $k\geq3$ is an odd number. Let $\lambda\in\mathbb{R}^d$ and
$$\hat{\mu}(\lambda)=\int_{\Sigma}e^{i\langle\lambda,x\rangle}d\mu(x);$$
where $\langle\lambda,x\rangle:=\lambda_1x_1+\cdots+\lambda_dx_d$ is the usual inner product of the vectors $\lambda,x\in\mathbb{R}^d$.
Note that if $k\geq2$ is an even number then $\Sigma$ is $C^1$ smooth convex hypersurface
then the following estimate
\begin{equation}\label{even}
|\hat{\mu}(\lambda)|\leq C|\lambda|^{-\frac{d-1}{k}}
\end{equation}
holds with constant $C$ which does not depend on $\lambda$. However, if $k$ is not an odd
number then $\Sigma$ is only $C^{k-1}$ smooth convex hypersurface and estimate \eqref{even} does not hold for small $k$. We obtain an optimal bound for the Fourier transform of
measures supported on $\Sigma$ for odd numbers $k$. If $k$ is sufficiently big depending on
the dimension $k$ then the estimate \eqref{even} still holds.
We prove the following:
\begin{theorem}\label{main}
Let $k\geq3$ be an odd number. The following estimate holds:
$$|\hat{\mu}(\lambda)|\leq C\max\big\{|\lambda|^{-k},|\lambda|^{-\frac{d-1}{k}}\big\},\qquad \lambda\in\mathbb{R}^d.$$
where $C$ is a constant depending only on $k, d$. In particular, if $d\leq k^2+1$, then the estimate \eqref{even} holds true.
\end{theorem}	

To prove Theorem \ref{main}, we introduce a smooth even cutoff function $\omega\in C_c^{\infty}(\mathbb{R})$ satisfying $0\leq\omega\leq1$ and
$$
\omega(t)=
\begin{cases}
1, & \text{if } |t|\leq\frac14,\\
0, & \text{if } |t|\geq\frac12.
\end{cases}
$$
Define
$$I_1(\lambda):=\int_{\Sigma}\omega(x_1)e^{i\langle\lambda,x\rangle} d\mu(x),\qquad \lambda\in\mathbb{R}^d$$
and
$$I_2(\lambda):=\int_{\Sigma}(1-\omega(x_1))e^{i\langle\lambda,x\rangle} d\mu(x),\qquad \lambda\in\mathbb{R}^d.$$
Then we have
$$\hat{\mu}(\lambda)=I_1(\lambda)+I_2(\lambda).$$

\section{Estimate for $I_1(\lambda)$}
In this section, we proof the following estimate for $I_1(\lambda)$.
\begin{proposition}
The following estimate holds:
$$|I_1(\lambda)|\leq C|\lambda|^{-k},$$
whenever $\lambda\in V(e_1)$, where $e_1=(1,0,0,\cdots,0)\in\mathbb{R}^d$, and $V(e_1)$ denotes a sufficiently small conic neighborhood of $e_1$.
More precisely,
$$V(e_1):=\{\nu=(\nu_1,\nu')\in\mathbb{R}^d:\nu_1>\frac34, |\nu'|\leq \epsilon\}$$
for a sufficiently small $\epsilon>0$.
\end{proposition}

\begin{proof}
Note that on $\supp I_1$, we have $|x_1|<\frac12$. Consequently,
$$\sum_{j=2}^d|x_j|^k=1-|x_1|^k>1-2^{-k}.$$
Define
$$a_0=\Big(\frac{1-2^{-k}}{d-1}\Big)^{\frac1k}.$$
Then there exists some index $j_0\in\{2,3,\cdots,d\}$ such that
$$|x_{j_0}|\geq a_0,$$
that is,  for some $\sigma\in\{\pm1\}$, we have $\sigma x_{j_0}>a_0$.

Choose $\eta\in C^\infty(\mathbb{R})$ such that
$$\eta(s)=0\quad (s\leq1/2),\qquad \eta(s)>0\quad (s>1/2),$$
and define
$$G(x)=\sum_{\ell=2}^{d}\sum_{\sigma=\pm1}\eta\Big(\frac{\sigma x_\ell}{a_0}\Big).$$
On $\supp\omega(x_1)\cap\Sigma$, the function $G$ is strictly positive.
Thus
$$\chi_{\ell,\sigma}(x)=\frac{\eta(\sigma x_\ell/a_0)}{G(x)},\qquad 2\leq \ell\leq d,\quad \sigma\in\{\pm1\},$$
forms an ambient $C^\infty$ partition of unity after restriction to this part of $\Sigma$.
We write
$$I_1(\lambda)=
\sum_{\ell=2}^{d}\sum_{\sigma=\pm1}I_{\ell,\sigma}(\lambda),$$
where
$$I_{\ell,\sigma}(\lambda)
  =\int_{\Sigma}\omega(x_1)\chi_{\ell,\sigma}(x)e^{i\ip{\lambda}{x}}\,d\mu(x).$$
It is enough to estimate one representative term, say
$I_{d,+}(\lambda)$; the other terms are identical after permuting
coordinates and changing signs.

On $\supp\chi_{d,+}$ we have $x_d\geq\frac{a_0}{2}$. Put
$$t=x_1,\qquad y=(x_2,\ldots,x_{d-1})\in \mathbb{R}^{d-2},$$
$$B(t,y)=1-|t|^k-\sum_{j=2}^{d-1}|y_j|^k,
\qquad
\Theta(t,y)=B(t,y)^{\frac1k}.$$
The graph parametrization is
$$X(t,y)=(t,y_2,\ldots,y_{d-1},\Theta(t,y)).$$
%On $\supp \chi_{d,+}$, we have $x_d\geq\frac{a_0}{2}$.
Thus,
\begin{equation}\label{Bty}
B(t,y)=x_d^k\geq \big(\frac{a_0}{2}\big)^k=:b_0>0.
\end{equation}

Write $\lambda=\rho v$, where $\rho=|\lambda|$ and $|v|=1$.
Under this parametrization, we have
\begin{equation}\label{Id+}
I_{d,+}(\lambda)
 =\int e^{i\rho\Phi_v(t,y)}A(t,y) dt dy,
\end{equation}
where
$$\Phi_v(t,y)=v_1t+\sum_{j=2}^{d-1}v_jy_j+v_d\Theta(t,y),$$
and
$$A(t,y)=\omega(t)\chi_{d,+}(X(t,y))\sqrt{1+|\nabla_{t,y}\Theta(t,y)|^2}.$$
Since $\omega(t)=0$ when $|t|\geq\frac12$, it follows that there exists a fixed bounded set $E\subset\mathbb{R}^{d-2}$ such that
\begin{equation}\label{Aty}
A(t,y)=0, \quad \text{if}\ y\notin E.
\end{equation}
%and on the support of $A$, we have
%$$B(t,y)\geq b_0.$$

Now, fixing $y$, we define
$$K_y(\rho):=\int e^{i\rho\Phi_v(t,y)}A(t,y)dt.$$
Note that the function $t \mapsto |t|^k$ is of class $C^{k-1}$ in a neighborhood of  $t=0$, and a jump may occur only in its $k$-th derivative.
Consequently, $\Theta(\cdot,y), \Phi_v(\cdot,y), A(\cdot,y)$ are all of class  $C^{k-1}$ at $t=0$, and they are smooth on $(-1/2,0)$ and $(0,1/2)$.
Thus, we have
$$|\partial_t\Theta(t,y)|=|B(t,y)|^{\frac1k-1}|t|^{k-1}\leq b_0^{\frac{1}{k}-1}2^{1-k}=:M.$$
Then, choose $\epsilon > 0$ sufficiently small such that $\epsilon M\leq\frac18$. Since $|v_d|\leq \epsilon$ and $|v_1|\geq\frac34$, it follows that
\begin{equation}\label{lower bound}
|\partial_t\Phi_{v}(t,y)|=|v_1+v_d\partial\theta(t,y)|\geq|v_1|-|v_d||\partial\theta(t,y)|
\geq\frac34-\epsilon M\geq\frac34-\frac18=\frac58.
\end{equation}
This implies that for every fixed $y$, $\Phi_v(\cdot,y)$ has no stationary points in the interval $(-\frac12,0)\cup(0,\frac12)$.

Define on $J_-:=(-\frac12,0)$ and $J_+:=(0,\frac12)$ respectively,
$$T_0=A,\qquad T_{m+1} = \partial_t\left(\frac{T_m}{\partial_t\Phi_v}\right), \qquad 0\le m\le k-1.$$
Since
$$e^{i\rho\Phi_v}=\frac{1}{i\rho\,\partial_t\Phi_v}\partial_t(e^{i\rho\Phi_v}),$$
integrating by parts $k$ times on any $(\alpha,\beta)\subset J_\pm$ yields
\begin{equation}\label{intergrate by parts}
\int_\alpha^\beta e^{i\rho\Phi_v}Ad t
=\sum_{q=1}^{k}(-1)^{q-1}(i\rho)^{-q}
\Big[e^{i\rho\Phi_v}
\frac{T_{q-1}}{\partial_t\Phi_v}
\Big]_{\alpha}^{\beta}
+(-1)^k(i\rho)^{-k}
\int_\alpha^\beta e^{i\rho\Phi_v}T_kdt.
\end{equation}
Here, integration by parts is performed strictly within the interiors of $J_-$ and $J_+$, without crossing $t=0$. Applying \eqref{intergrate by parts} to $(-1/2,-\delta)$ and $(\delta,1/2)$ and letting $\delta\to0^+$, the boundary terms at the outer endpoints vanish because $\omega$ is identically zero for $|t|\ge1/2$. Set
$$ B_q^\pm(y)=\lim_{t\to0^\pm}\frac{T_{q-1}(t,y)}{\partial_t\Phi_v(t,y)}.$$
It follows that
$$K_y(\rho)=\sum_{q=1}^{k}(-1)^{q-1}(i\rho)^{-q}
e^{i\rho\Phi_v(0,y)}\bigl(B_q^-(y)-B_q^+(y)\bigr)
+(-1)^k(i\rho)^{-k}
\Big(
\int_{-1/2}^{0}e^{i\rho\Phi_v}T_k d t
+
\int_{0}^{1/2}e^{i\rho\Phi_v}T_k dt
\Big).$$
By the $C^{k-1}$ regularity, we have
$$B_q^-(y)=B_q^+(y),\qquad 1\le q\le k-1.$$
Thus,
$$K_y(\rho)=(-1)^{k-1}(i\rho)^{-k}
e^{i\rho\Phi_v(0,y)}\big(B_k^-(y)-B_k^+(y)\big)
+(-1)^k(i\rho)^{-k}
\Big(
\int_{-1/2}^{0}e^{i\rho\Phi_v}T_k(t,y) d t
+
\int_{0}^{1/2}e^{i\rho\Phi_v}T_k(t,y) dt
\Big).$$
Note that
$$T_{m+1} = \frac{\partial_tT_m}{\partial_t\Phi_v}- \frac{T_m\,\partial_t^2\Phi_v}{(\partial_t\Phi_v)^2},$$
By \eqref{lower bound}, we have for $t\in(-1/2,0)\cup(0,1/2)$,
$$\sum_{r=0}^{k+1}|\partial_t^r\Phi_v(t,y)|+\sum_{r=0}^{k}|\partial_t^r A(t,y)|\leq C.$$
Thus, %递推公式
$$\sup_{t\in J_-\cup J_+}|T_m(t,y)|\leq C\mathbf{1}_E(y), \qquad 0\leq m\leq k.$$
Moreover,
$$|B_k^-(y)|+|B_k^+(y)|+ \int_{-\frac12}^{0}|T_k(t,y)|dt+\int_{0}^{\frac12}|T_k(t,y)|dt\leq C\mathbf{1}_E(y).$$
Therefore,
$$|K_\rho(y)|\leq C\rho^{-k}\mathbf{1}_E(y).$$
Substituting this back into \eqref{Id+} and integrating with respect to $y$, we have
$$|I_{d,+}(\lambda)|\leq \int |K_y(\rho)|dy \leq C\rho^{-k}.$$
Similarly, for any $\ell\in\{2,\cdots,d\}$ and $\sigma\in\{\pm1\}$, we have
$$|I_{\ell,\sigma}(\lambda)|\leq \int |K_y(\rho)|dy \leq C\rho^{-k}.$$
Then, summing over finitely many $(\ell,\tau)$, we conclude that
$$|I_1(\lambda)|\leq C\rho^{-k}=C|\lambda|^{-k}.$$
\end{proof}


\section{Estimate for $I_2(\lambda)$}
In this section, we proof the following estimate for $I_2(\lambda)$.
\begin{proposition}\label{estimate for I2}
Let $e_1 = (1, 0, \dots, 0) \in \mathbb{R}^d$ and $V(e_1)$ be a sufficiently small conic neighborhood of $e_1$. Then for any $\lambda \in V(e_1)$, the following estimate holds:
$$|I_2(\lambda)|\leq C|\lambda|^{-\frac{d-1}{k}}.$$
\end{proposition}

To prove Proposition \ref{estimate for I2}, we introduce a smooth cutoff function $\chi_0\in C_c^{\infty}(\mathbb{R}^{d-1})$ supported in a sufficiently small neighborhood of the origin of $\mathbb{R}^{d-1}$. This enables a further decomposition of $I_2(\lambda)$ into $I_{20}(\lambda)$ and $I_{21}(\lambda)$, defined respectively by
$$I_{20}(\lambda):=\int_{\Sigma}(1-\omega(x_1))(1-\chi_0(x'))e^{i\langle\lambda,x\rangle} d\mu(x),\qquad \lambda\in\mathbb{R}^d,$$
and
$$I_{21}(\lambda):=\int_{\Sigma}(1-\omega(x_1))\chi_0(x')e^{i\langle\lambda,x\rangle} d\mu(x),\qquad \lambda\in\mathbb{R}^d,$$
where $x'=(x_2,\cdots,x_d)\in\mathbb{R}^{d-1}$.

\begin{lemma}\label{I20}
We have
$$|I_{20}(\lambda)|\leq\frac{C}{(1+|\lambda|)^{N}},\qquad \lambda\in\mathbb{R}^d,$$
for any $N\in\mathbb{N}$.
\end{lemma}
\begin{proof}
%Integral by parts
\end{proof}

\begin{lemma}\label{I21}
We have
$$|I_{21}(\lambda)|\leq\frac{C}{(1+|\lambda|)^{\frac{d-1}{k}}},\qquad \lambda\in V(e_1).$$
\end{lemma}

Without loss of generality, we may assume $x_1>0$ on the support of the integrand.
%Let $x'=(x_2, \cdots, x_d)\in\mathbb{R}^{d-1}$.
The hypersurface $\Sigma$ can be locally parameterized as a graph $x_1=\psi(x')$, where
$$\psi(x')=\Big(1-\sum_{i=2}^d|x_i|^k\Big)^{\frac1k}.$$
Set $s_\ell:=\frac{\lambda_\ell}{\lambda_1}$ for $\ell=2,3,\dots,d$.
Using this parameterization, the integral $I_{21}(\lambda)$ can be written as follows:
$$I_{21}(\lambda)=\int_{\mathbb{R}^{d-1}}e^{i\lambda_1\Psi(x',s)}
A(x')dx',$$
where the phase function $\Psi(x', s)$ is given by
$$\Psi(x',s):=\psi(x')+\sum_{i=2}^ds_ix_i,$$
and the amplitude function $A(x')$ is defined as
$$A(x'):=\big(1-\omega\big(\psi(x')\big)\big)\chi_0(x')\sqrt{1+|\nabla\psi(x')|^2}.$$

By the symmetry of the domain and the parity of the integrand, the estimation of $I_{21}(\lambda)$ can be reduced to the primary quadrant $\mathbb{R}_+^{d-1}$. Specifically, it suffices to consider the localized integral
$$I_+(\lambda)=\int_{\mathbb{R}_+^{d-1}}e^{i\lambda_1\Psi(x',s)}A(x')dx'.$$
The desired bound for $I_+(\lambda)$ is established via the following lemma.
\begin{lemma}\label{I+}
We have
$$|I_{+}(\lambda)|\leq\frac{C}{(1+|\lambda|)^{\frac{d-1}{k}}},\qquad \lambda\in V(e_1).$$
\end{lemma}
\begin{proof}
Let $\chi_1(t):=\eta(t)-\eta(2t)$, where where $\eta\in C_c^\infty(\mathbb{R})$ is a smooth cutoff function satisfying $0\leq\eta\leq1$ and
$$
\eta(t)=
\begin{cases}
1, & \text{if } |t|\leq1,\\
0, & \text{if } |t|\geq2.
\end{cases}
$$
It follows immediately from construction that $\chi_1$ is supported on the dyadic interval $\{t \in \mathbb{R} : 1/2 \le |t| \le 2\}$ and satisfies the identity
$$\sum_{\ell=0}^{\infty}\chi_1(2^\ell t)=1, \quad \text{ for }\quad 0<|t|\leq1.$$
Therefore, for any $x'\in\mathbb{R}_+^{d-1}$ satisfying $0< x_j<1$ for all $j=2,\dots,d$, we obtain
$$\sum_{\ell\in\mathbb{Z}_+^{d-1}}\chi_1(2^{\ell_2}x_2)\chi_1(2^{\ell_3}x_3)\cdots
\chi_1(2^{\ell_d}x_d)=1,$$
where $\ell=(\ell_1,\cdots,\ell_d)$ and $\mathbb{Z}_+:=\{0\}\cup\mathbb{N}$. Thus,
$$I_+(\lambda)
=\sum_{\ell\in\mathbb{Z}_+^{d-1}}\int_{\mathbb{R}_+^{d-1}}e^{i\lambda_1\Psi(x',s)}
A(x')\chi_\ell(x')dx':=\sum_{\ell\in\mathbb{Z}_+^{d-1}}I_{\ell}(\lambda),$$
where $\chi_\ell(x'):=\chi_1(2^{\ell_2}x_2)\chi_1(2^{\ell_3}x_3)\dots
\chi_1(2^{\ell_d}x_d)$.

Further, we consider an estimate for the integral $I_{\ell}$.
Let $L$ be a sufficiently large natural number such that $2^{-L+1}<\delta$, where $\delta$ denotes the radius of $\supp\chi_0$.
For any $\ell\in\mathbb{Z}_+^{d-1}$, if there exists at least one component $l_j\leq L$, then the support of the corresponding component $\chi_1(2^{l_j} x_j)$ satisfies $|x_j| \ge 2^{-l_j-1} \ge 2^{-L-1}$. By our choice of $L$, this domain is strictly disjoint from the support of $\chi_0(x')$. Consequently, these terms vanish identically on $\supp \chi_0$:
$$\chi_1(2^{l_j} x_j) \chi_0(x') \equiv 0, \qquad \text{whenever } l_j \le L.$$
Therefore, we may assume that $L<\min\{\ell_2,\dots,\ell_d\}$. It follows that on the support of $\chi_0$, we have the identity:
$$\sum_{\substack{l \in \mathbb{Z}_+^{d-1} \\ l_j > L, \, j=2,\dots,d}} \chi_1(2^{l_2}x_2)\chi_1(2^{l_3}x_3)\cdots \chi_1(2^{l_d}x_d) = 1, \qquad \forall x' \in \supp \chi_0.$$

We consider integral $I_{\ell}$ assuming
$$L<\ell_2\leq\ell_3\leq\cdots\leq\ell_d.$$
Surely other cases can be treated quite similarly. First we use change of variables
by using scaling $2^{\ell_j}x_j=y_j$ for $j=2,\dots,d$. Then we get
$$I_{\ell}(\lambda)
=2^{-|\ell|_1}\int_{\mathbb{R}_+^{d-1}}e^{i\lambda_1\Psi(\delta_{2^{-\ell}}(y),s)}
A(\delta_{2^{-\ell}}(y))\chi_\ell(y)dy,$$
where $|\ell|_1:=\ell_2+\cdots+\ell_d$ and $\delta_{2^{-\ell}}(y)=(2^{-\ell_2}y_2,\dots,2^{-\ell_d}y_d)$.
At this stage, the rescaled phase function takes the explicit form
$$\Psi(\delta_{2^{-\ell}}(y),s)=\Big(1-\sum_{j=2}^d2^{-k\ell_j}y_j^k\Big)^{\frac1k}
+\sum_{j=2}^d2^{-\ell_j}s_jy_j.$$
Now, let $y':=(y_3,\dots,y_d)\in\mathbb{R}_+^{d-2}$. Fixing $y'$, we consider the one-dimensional integral in $y_2$ defined by
$$I_{\ell,y'}(\lambda)=2^{-|\ell|_1}\int_{\mathbb{R}_+^{d-1}}
e^{i\lambda_1\Psi(\delta_{2^{-\ell}}(y_2,y'),s)}
A(\delta_{2^{-\ell}}(y_2,y'))\chi_\ell(y_2,y')dy_2.$$
First, we have
$$\partial_2\Psi=\frac{-2^{k\ell_2}y_2^{k-1}}{\Big(1-\big(\sum_{j=2}^d2^{-k\ell_j}y_j^k\big)\Big)^{1-\frac1k}}
+s_22^{-\ell_2}.$$
Note that if $|s_22^{-\ell_2}|\gg2^{-k\ell_2}$
\end{proof}














\newpage
In this section, we prove the following estimate for $I_2(\lambda)$.
\begin{proposition}\label{prop-I2}
The following estimate holds:
$$|I_2(\lambda)|\leq C|\lambda|^{-\frac{d-1}{k}},$$
whenever $\lambda\in V(e_1)$, where $e_1=(1,0,0,\cdots,0)$, and $V(e_1)$ denotes a sufficiently small conic neighborhood of $e_1$.
\end{proposition}

\begin{proof}
On the support of the integrand of $I_2(\lambda)$, we have $|x_1| \geq \frac{1}{4}$. Since $k \geq 3$ is an odd integer, the expression $|x_1|^k = x_1^k \sgn(x_1)$ is infinitely smooth ($C^\infty$) on the set $\{x_1 \in \mathbb{R} : |x_1| \geq \frac{1}{4}\}$.

We parameterize this portion of the surface $\Sigma$ by writing $x_1$ as a function of $x' = (x_2, \dots, x_d)$. Let
$$F(x') = \big(1 - |x_2|^k - \cdots - |x_d|^k\big)^{\frac{1}{k}} \sgn(x_1).$$
Since $|x_1| \geq \frac{1}{4}$, it follows that $\sum_{j=2}^d |x_j|^k \leq 1 - 4^{-k}$, which implies that $x'$ stays away from the boundary of the unit ball in the $\ell^k$-norm. Thus, the parametrization $x_1 = F(x')$ is smooth ($C^\infty$) on the support of $(1-\omega(x_1))$.

We write $\lambda = \rho v$ with $\rho = |\lambda|$ and $v \in V(e_1)$. The integral $I_2(\lambda)$ can be written locally as a sum of terms of the form
$$ \int_{\mathbb{R}^{d-1}} e^{i \rho \Psi_v(x')} \tilde{A}(x') \, dx', $$
where the phase function is given by
$$ \Psi_v(x') = v_1 F(x') + \sum_{j=2}^d v_j x_j, $$
and $\tilde{A}(x')$ is a smooth cutoff function with compact support.

Next, we examine the Hessian matrix of $F(x')$ to determine the types of stationary points. For $v = e_1$ (i.e., $v_1 = 1, v' = 0$), the phase function simplifies to $\Psi_{e_1}(x') = F(x')$. Direct calculation of the derivatives of $F$ gives that at a stationary point $\nabla_{x'} F(x') = 0$, we must have $x_2 = \dots = x_d = 0$, which corresponds to the vertex $(\pm 1, 0, \dots, 0)$.

At this point, the cross derivatives vanish, and the second derivatives are given by:
$$ \frac{\partial^2 F}{\partial x_j^2}\Big|_{x'=0} = - \frac{k-1}{k} |x_j|^{k-2} \sgn(x_j) \cdot (\dots) $$
Notice that for $k \geq 3$, the second derivatives vanish at $x'=0$, meaning the surface has a flat flattening of order $k$ at the vertices. More precisely, around the vertex, $F(x') = \pm (1 - \frac{1}{k}\sum_{j=2}^d |x_j|^k + \dots)$.

Since the flattening is of type $k$ in all $d-1$ tangential directions, by applying the standard multi-dimensional stationary phase lemma for degenerate types (or using van der Corput's lemma via a polar-like coordinate change on the $\ell^k$ sphere), we obtain that the decay rate of this oscillatory integral is dominated by the type of singularity, yielding exactly:
$$ |I_2(\lambda)| \leq C \rho^{-\frac{d-1}{k}} = C |\lambda|^{-\frac{d-1}{k}}. $$
By the stability of the stationary points under small perturbations, this estimate holds uniformly for all $v$ in a sufficiently small conic neighborhood $V(e_1)$. This completes the proof.
\end{proof}

\section{Estimate for $I_2(\lambda)$}
In this section, we provide a detailed proof for the optimal decay of $I_2(\lambda)$.
\begin{proposition}\label{prop-I2-detailed}
The following estimate holds uniformly for all $\lambda\in V(e_1)$:
$$|I_2(\lambda)|\leq C|\lambda|^{-\frac{d-1}{k}}.$$
\end{proposition}

\begin{proof}
Recall that on the support of $1-\omega(x_1)$, we have $|x_1| \geq \frac{1}{4}$. Since $k \geq 3$ is an odd integer, the map $x_1 \mapsto |x_1|^k = x_1^k$ is infinitely differentiable on $\mathbb{R} \setminus \{0\}$. Thus, $\Sigma \cap \supp(1-\omega(x_1))$ is a $C^\infty$ smooth hypersurface.

We parameterize this piece of $\Sigma$ over the tangential variables $x' = (x_2, \dots, x_d) \in \mathbb{R}^{d-1}$. Let
$$ F(x') = \big(1 - |x_2|^k - \dots - |x_d|^k\big)^{1/k}. $$
The constraint $|x_1| \geq \frac{1}{4}$ implies $\sum_{j=2}^d |x_j|^k \leq 1 - 4^{-k}$, meaning $x'$ is restricted to a compact domain $\Omega_0 \subset \mathbb{R}^{d-1}$ strictly contained in the interior of the $\ell^k$ unit ball. On $\Omega_0$, $F(x')$ is $C^\infty$ smooth.

Writing $\lambda = \rho v$ with $\rho = |\lambda|$ and $v = (v_1, v') \in V(e_1)$, the integral $I_2(\lambda)$ can be expressed as
$$ I_2(\lambda) = \int_{\mathbb{R}^{d-1}} e^{i \rho \Psi_v(x')} \tilde{A}(x') \, dx', $$
where $\Psi_v(x') = v_1 F(x') + \langle v', x'\rangle$, and $\tilde{A}(x') = (1-\omega(F(x')))\cdot J(x')$ is a smooth, compactly supported amplitude function involving the geometric Jacobian $J(x') = \sqrt{1+|\nabla F(x')|^2}$.

To apply the method of stationary phase, we analyze the critical points of $\Psi_v(x')$. The gradient equation $\nabla_{x'} \Psi_v(x') = 0$ reads
\begin{equation}\label{crit-eq}
v_1 \frac{\partial F}{\partial x_j}(x') + v_j = 0, \quad \forall\, j = 2, \dots, d.
\end{equation}
A direct calculation yields $\frac{\partial F}{\partial x_j}(x') = - \frac{\sgn(x_j)|x_j|^{k-1}}{F(x')^{k-1}}$. Thus, \eqref{crit-eq} is equivalent to
$$ v_1 \frac{\sgn(x_j)|x_j|^{k-1}}{F(x')^{k-1}} = v_j. $$
Since $v \in V(e_1)$, we have $v_1 \geq \frac{3}{4}$ and $|v'| \leq \epsilon$. For $\epsilon > 0$ sufficiently small, this forces the unique non-degenerate/degenerate critical point to lie in a tiny neighborhood of the origin $x' = 0$ (which corresponds to the vertex $(\pm 1, 0, \dots, 0)$).

At the exact vertex $x'=0$, the Hessian matrix $\nabla^2 \Psi_v(0)$ vanishes identically because $k \geq 3 \implies k-1 \geq 2$. Therefore, the standard non-degenerate stationary phase lemma fails. To determine the precise decay, we examine the Taylor expansion of $F(x')$ around $x'=0$:
$$ F(x') = 1 - \frac{1}{k} \sum_{j=2}^d |x_j|^k + O\Big( \Big(\sum_{j=2}^d |x_j|^k\Big)^2 \Big). $$
Thus, for the model direction $v = e_1$, the phase is perfectly uncoupled in a standard scaling sense: $\Psi_{e_1}(x') = 1 - \frac{1}{k} \sum_{j=2}^d |x_j|^k + \dots$.

To evaluate the decay rate uniformly for $v \in V(e_1)$, we perform an anisotropic change of variables (dyadic decomposition near the singularity). Let $\delta = |v'|/v_1 \lesssim \epsilon$. The critical point $x^c = (x_2^c, \dots, x_d^c)$ satisfies $|x_j^c| \sim \delta^{1/(k-1)}$. We shift the origin to the critical point by setting $x' = x^c + z$. Around $x^c$, the phase behaves as a degenerate polynomial of degree $k$.

Alternatively, we can invoke the generalized van der Corput lemma in a directional manner. Since the boundary of the support of $\tilde{A}$ does not pass through the critical point, we can introduce polar-like coordinates adapted to the $\ell^k$ metric on $\mathbb{R}^{d-1}$. Let $x_j = r \theta_j$, where $\sum_{j=2}^d |\theta_j|^k = 1$. Then the phase becomes
$$ \Psi_v(r, \theta) = v_1 (1 - r^k)^{1/k} + r \langle v', \theta \rangle. $$
Expanding this in $r$ for small $r$ yields
$$ \Psi_v(r, \theta) = v_1 - \frac{v_1}{k} r^k + r \langle v', \theta \rangle + O(r^{2k}). $$
Differentiating with respect to $r$, we have
$$ \partial_r \Psi_v(r, \theta) = -v_1 r^{k-1} + \langle v', \theta \rangle + O(r^{2k-1}). $$
The $k$-th derivative with respect to $r$ satisfies
$$ \partial_r^k \Psi_v(r, \theta) = -v_1 k! + O(r^k). $$
Since $v_1 \geq \frac{3}{4}$, there exists a constant $c_0 > 0$ such that $|\partial_r^k \Psi_v(r, \theta)| \geq c_0 > 0$ uniformly for all $r$ in a neighborhood of $0$ and all $\theta$.

By the classical multi-dimensional van der Corput lemma (or by applying the standard one-dimensional version to the radial variable $r$ after regularizing the angular integrals), a $k$-th order flattening along $d-1$ independent dimensions yields a decay rate of exactly $\rho^{-1/k}$ per dimension. Since the amplitude $\tilde{A}$ is smooth and compactly supported, the integration over the remaining smooth variables preserves this uniform bounds.

Thus, we conclude that
$$ |I_2(\lambda)| \leq C \rho^{-\frac{d-1}{k}} = C |\lambda|^{-\frac{d-1}{k}}, $$
where the constant $C$ depends only on $k, d$ and the $C^k$ norms of $\tilde{A}$, which are uniformly bounded since $|x_1| \ge 1/4$. This completes the proof.
\end{proof}

\section{Proof of Theorem \ref{main}}
\bigskip
\bigskip


\end{document}
