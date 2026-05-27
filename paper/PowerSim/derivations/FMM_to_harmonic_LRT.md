# From the FMM Model to the K-Harmonic LRT

**Purpose:** complete derivation showing that the K-harmonic LRT used in
`detect_FMM()` is the linearized rhythmicity test motivated by the
truncated Fourier expansion of the FMM signal.

This document is meant to be a self-contained reference: it derives the
Möbius form, the Fourier expansion, the FMM-to-harmonic identification,
and the hypothesis test in one place.

---

## 1. The FMM signal model

The Frequency Modulated Möbius (FMM) model (Rueda, Larriba & Peddada, 2019)
represents a single rhythmic gene's expression profile as

$$y(t) = M + A\cos\{\beta + u(t)\} + \varepsilon, \qquad \varepsilon \sim \mathcal{N}(0,\sigma^2)$$

where the **phase transformation** is

$$u(t) = 2\arctan\!\left(\omega\tan\!\frac{t-\alpha}{2}\right), \qquad \omega \in (0, 1].$$

Parameter interpretations:
- $M$ — mesor (baseline expression)
- $A \ge 0$ — amplitude
- $\alpha \in [0, 2\pi)$ — location/peak time (phase)
- $\beta \in [0, 2\pi)$ — orientation
- $\omega \in (0, 1]$ — asymmetry / peak sharpness; $\omega = 1$ recovers a pure sinusoid (cosinor)

The FMM permits sharp, asymmetric peaks while remaining periodic. The
question we address is: **how do we test whether $A > 0$ (the signal is
rhythmic)?**

The natural likelihood-ratio test is intractable because $\{\beta, \alpha,
\omega\}$ are unidentifiable when $A = 0$ (the **Davies problem** —
nuisance parameters present only under the alternative). We sidestep this
by reparameterizing the FMM signal as a Fourier series and testing in that
linear basis.

---

## 2. Möbius form of the FMM phase

Let $s = t - \alpha$. The FMM signal is

$$y(t) = M + A\cos\{\beta + u\} + \varepsilon, \qquad u = 2\arctan(\omega\tan(s/2)).$$

We compute $e^{iu}$ explicitly. Using the identity
$e^{2i\arctan x} = (1 + ix)/(1 - ix)$ with $x = \omega\tan(s/2)$,

$$e^{iu} = \frac{1 + i\omega\tan(s/2)}{1 - i\omega\tan(s/2)}.$$

The half-angle-to-exponential identity is

$$\tan(s/2) \;=\; \frac{e^{is} - 1}{i(e^{is} + 1)} \;\;\Longrightarrow\;\; i\tan(s/2) = \frac{e^{is}-1}{e^{is}+1}.$$

Substituting into the numerator and denominator separately:

$$\begin{aligned}
1 + i\omega\tan(s/2)
&= 1 + \frac{\omega(e^{is}-1)}{e^{is}+1}
= \frac{(e^{is}+1) + \omega(e^{is}-1)}{e^{is}+1} \\[4pt]
1 - i\omega\tan(s/2)
&= 1 - \frac{\omega(e^{is}-1)}{e^{is}+1}
= \frac{(e^{is}+1) - \omega(e^{is}-1)}{e^{is}+1}.
\end{aligned}$$

Therefore

$$e^{iu} = \frac{(e^{is}+1) + \omega(e^{is}-1)}{(e^{is}+1) - \omega(e^{is}-1)}
       = \frac{(1+\omega)e^{is} + (1-\omega)}{(1-\omega)e^{is} + (1+\omega)}.$$

Dividing numerator and denominator by $(1+\omega)$ and defining

$$\tilde r \equiv \frac{1-\omega}{1+\omega}, \qquad \tilde r \in [0, 1),$$

we obtain the **Möbius / Blaschke form**:

$$\boxed{\;e^{iu(t)} \;=\; \frac{e^{is} + \tilde r}{1 + \tilde r\, e^{is}} \;=\; \frac{e^{is} - r}{1 - r\, e^{is}},\quad r \equiv -\tilde r = \frac{\omega - 1}{\omega + 1}.\;}$$

For $\omega \in (0, 1]$, $r \in (-1, 0]$. At $\omega = 1$ (pure cosinor),
$r = 0$ and the Möbius transformation reduces to the identity.

**Sanity check** (numerical): at $\omega = 0.5$, $s = \pi/2$:
- Direct: $u = 2\arctan(0.5\cdot 1) = 0.927$, so $e^{iu} = \cos 0.927 + i\sin 0.927 = 0.6 + 0.8i$.
- Möbius form: $r = -1/3$, $e^{is} = i$, so $(i - (-1/3))/(1 - (-1/3)\cdot i) = (i + 1/3)/(1 + i/3) = 0.6 + 0.8i.$ ✓

---

## 3. Geometric Fourier expansion

Since $|r| < 1$ for $\omega \in (0, 1]$, the geometric series

$$\frac{1}{1 - r\,e^{is}} = \sum_{j=0}^{\infty} r^j e^{ijs}$$

converges. Therefore

$$\begin{aligned}
e^{iu(t)}
&= (e^{is} - r) \sum_{j=0}^{\infty} r^j e^{ijs} \\[4pt]
&= \sum_{j=0}^{\infty} r^j e^{i(j+1)s} \;-\; r\sum_{j=0}^{\infty} r^j e^{ijs} \\[4pt]
&= \sum_{k=1}^{\infty} r^{k-1} e^{iks} \;-\; r \;-\; \sum_{k=1}^{\infty} r^{k+1} e^{iks} \\[4pt]
&= -r \;+\; \sum_{k=1}^{\infty} (r^{k-1} - r^{k+1}) e^{iks} \\[4pt]
&= \boxed{\;-r \;+\; (1 - r^2) \sum_{k=1}^{\infty} r^{k-1} e^{iks}.\;}
\end{aligned}$$

Taking real parts and incorporating $\beta$:

$$\begin{aligned}
\cos\{\beta + u(t)\}
&= \mathrm{Re}\{e^{i\beta} e^{iu(t)}\} \\
&= -r \cos\beta \;+\; (1 - r^2) \sum_{k=1}^{\infty} r^{k-1} \cos\{ks + \beta\} \\
&= -r \cos\beta \;+\; (1 - r^2) \sum_{k=1}^{\infty} r^{k-1} \cos\{k(t - \alpha) + \beta\}.
\end{aligned}$$

The full FMM signal becomes

$$y(t) = M' \;+\; A(1 - r^2) \sum_{k=1}^{\infty} r^{k-1} \cos\{k(t - \alpha) + \beta\} \;+\; \varepsilon, \qquad M' \equiv M - Ar\cos\beta.$$

So **the FMM signal is exactly an infinite Fourier series in $t$**, with
harmonic amplitudes that decay geometrically:

$$|c_k| = A(1 - r^2)\, |r|^{k-1}, \qquad k = 1, 2, 3, \ldots$$

The **decay rate** is $|r| = (1 - \omega)/(1 + \omega)$, increasing from
$0$ (cosinor, $\omega = 1$) toward $1$ (very sharp peak, $\omega \to 0$).

---

## 4. Variance captured by truncation at K harmonics

By Parseval, the contribution of harmonic $k$ to total signal variance is
$\propto |c_k|^2 = A^2(1-r^2)^2 r^{2(k-1)}$. Total variance:

$$V_\infty = A^2(1-r^2)^2 \sum_{k=1}^{\infty} r^{2(k-1)} = A^2(1-r^2)^2 \cdot \frac{1}{1-r^2} = A^2(1-r^2).$$

Variance in the first $K$ harmonics:

$$V_K = A^2(1-r^2)^2 \cdot \frac{1 - r^{2K}}{1-r^2} = A^2(1-r^2)(1 - r^{2K}).$$

Therefore the **variance fraction captured by truncating at K** is

$$\boxed{\;\frac{V_K}{V_\infty} \;=\; 1 - r^{2K} \;=\; 1 - \left(\frac{1-\omega}{1+\omega}\right)^{2K}.\;}$$

| $\omega$ | $\|r\|$ | $K=1$ (DCP) | $K=2$ | $K=3$ | $K=4$ |
|---|---|---|---|---|---|
| 0.50 | 0.333 | 88.9% | 98.8% | 99.9% | 99.98% |
| 0.40 | 0.429 | 81.6% | 96.6% | 99.4% | 99.9% |
| 0.30 | 0.538 | 71.0% | 91.6% | 97.6% | 99.3% |
| 0.20 | 0.667 | 55.6% | 80.2% | 91.2% | 96.1% |
| 0.10 | 0.818 | 33.0% | 55.2% | 70.1% | 80.1% |

For empirically observed $\omega \in [0.20, 0.40]$, $K = 2$ captures
80–97% of FMM variance; $K = 3$ captures 91–99%.

---

## 5. From the FMM signal to the K-harmonic regression model

Truncating the Fourier series at $K$ and re-expressing each harmonic
$\cos\{k(t-\alpha) + \beta\}$ as a sum of $\cos$ and $\sin$:

$$\cos\{k(t-\alpha) + \beta\} = \cos(k\alpha - \beta)\cos(k\omega_0 t) + \sin(k\alpha - \beta)\sin(k\omega_0 t),$$

where $\omega_0 = 2\pi/T$ is the fundamental angular frequency. Defining

$$a_k = A(1-r^2) r^{k-1} \cos(\beta - k\alpha), \qquad b_k = A(1-r^2) r^{k-1} \sin(\beta - k\alpha),$$

the truncated FMM signal is

$$y(t) \approx M + \sum_{k=1}^{K} \big[a_k \cos(k\omega_0 t) + b_k \sin(k\omega_0 t)\big] + \varepsilon.$$

The FMM model imposes a **constraint** on the harmonic amplitudes:

$$\sqrt{a_k^2 + b_k^2} = A(1-r^2) |r|^{k-1}\quad\text{(geometric in }k\text{).}$$

**Releasing this constraint** — letting $(a_k, b_k)$ be free for each
$k = 1, \ldots, K$ — yields the **K-harmonic linear regression model**:

$$\boxed{\;y(t) \;=\; M \;+\; \sum_{k=1}^{K} \big[a_k \cos(k\omega_0 t) + b_k \sin(k\omega_0 t)\big] \;+\; \varepsilon.\;}$$

This model has $1 + 2K$ parameters: an intercept and $2K$ harmonic
coefficients. The FMM model has 5 parameters $(M, A, \beta, \alpha, \omega)$,
so at $K = 2$ both models have the **same parameter count**. The
K-harmonic regression is a strict generalization in the sense that it
spans a larger family of waveforms (any 2-harmonic signal) at the same
df cost, while sacrificing the geometric-decay constraint.

---

## 6. The hypothesis test

We test the null hypothesis of **no rhythm**:

$$H_0: A = 0 \quad \Longleftrightarrow \quad a_1 = b_1 = a_2 = b_2 = \cdots = a_K = b_K = 0$$

against the alternative

$$H_1: \exists k \in \{1, \ldots, K\} \text{ with } (a_k, b_k) \ne (0, 0).$$

### 6.1. Models under H₀ and H₁

Stack the $n$ time points into a column vector $\mathbf{y}\in\mathbb{R}^n$. Build the design matrix

$$\mathbf{X}_K = \big[\,\mathbf{1},\; \cos(\omega_0 \mathbf{t}),\; \sin(\omega_0 \mathbf{t}),\; \ldots,\; \cos(K\omega_0 \mathbf{t}),\; \sin(K\omega_0 \mathbf{t})\,\big] \in \mathbb{R}^{n \times (2K+1)}.$$

Under $H_0$ the design is just the intercept column $\mathbf{X}_0 = \mathbf{1}$.

### 6.2. Sum-of-squared errors

$$\mathrm{SSE}_0 = \min_M \|\mathbf{y} - M\mathbf{1}\|^2 = \sum_{i=1}^n (y_i - \bar y)^2$$

$$\mathrm{SSE}_K = \min_{(M, a_k, b_k)} \|\mathbf{y} - \mathbf{X}_K \boldsymbol{\theta}\|^2 = \mathbf{y}^\top(\mathbf{I} - \mathbf{H}_K)\mathbf{y},$$

where $\mathbf{H}_K = \mathbf{X}_K(\mathbf{X}_K^\top \mathbf{X}_K)^{-1}\mathbf{X}_K^\top$ is the projection onto the column space of $\mathbf{X}_K$.

### 6.3. Test statistic

The standard nested F-test:

$$\boxed{\;F = \frac{(\mathrm{SSE}_0 - \mathrm{SSE}_K) / (2K)}{\mathrm{SSE}_K / (n - 2K - 1)}.\;}$$

Equivalently, the LRT statistic in Gaussian linear regression is

$$\Lambda = -n \log\!\frac{\mathrm{SSE}_K}{\mathrm{SSE}_0} = -n \log(1 - R_K^2),$$

which is a strictly monotone function of $F$, so the two tests give
identical p-value rankings.

### 6.4. Null distribution — exact

Under Gaussian errors, the F statistic has **exactly** the F-distribution
under $H_0$:

$$\boxed{\;F \mid H_0 \;\sim\; F\!\left(2K,\; n - 2K - 1\right) \quad\text{(exact, finite-sample)}.\;}$$

The p-value is

$$p = P\!\left[F(2K, n-2K-1) > F_{\text{obs}}\right] = 1 - F_{\text{cdf}}(F_{\text{obs}};\, 2K,\, n-2K-1).$$

No asymptotics, no boundary correction, no Davies supremum, no empirical
calibration table — the null distribution is exact because the K-harmonic
model is **fully linear** in $(M, a_1, b_1, \ldots, a_K, b_K)$, and all
parameters are identifiable under both $H_0$ (where they equal zero) and
$H_1$.

### 6.5. Why this avoids the FMM-LRT's boundary problem

The FMM model directly tests $H_0: A = 0$ in the nonlinear
parameterization. Under $H_0$ with $A = 0$, the parameters $\{\beta,
\alpha, \omega\}$ become **unidentifiable** — they multiply zero. Davies
(1977, 1987) showed that the LRT statistic in this setting follows the
supremum of a Gaussian process over the unidentifiable nuisance space,
not a chi-squared distribution. Empirical calibration is required, with
no closed-form null. The "Davies tax" inflates critical values
substantially relative to nominal F-test thresholds.

The K-harmonic regression escapes this entirely: every parameter $a_k, b_k$
is identifiable as zero under $H_0$, so the F-distribution applies
exactly. The trade-off is the geometric-decay constraint, which is
abandoned (slightly larger alternative space at the same df).

---

## 7. Summary of the construction

| Step | Result |
|---|---|
| FMM model | $y(t) = M + A\cos(\beta + u(t)) + \varepsilon,\; u = 2\arctan(\omega\tan((t-\alpha)/2))$ |
| Möbius form | $e^{iu} = (e^{is} - r)/(1 - re^{is}),\; r = (\omega-1)/(\omega+1)$ |
| Fourier expansion | $y(t) = M' + A(1-r^2)\sum_{k\ge1} r^{k-1}\cos(k(t-\alpha)+\beta) + \varepsilon$ |
| Truncation at K | $y(t) \approx M + \sum_{k=1}^K [a_k\cos(k\omega_0 t) + b_k\sin(k\omega_0 t)] + \varepsilon$ |
| Test | $H_0: a_1 = b_1 = \cdots = a_K = b_K = 0$ |
| Statistic | $F = \frac{(\mathrm{SSE}_0 - \mathrm{SSE}_K)/(2K)}{\mathrm{SSE}_K/(n - 2K - 1)}$ |
| Null distribution | $F(2K, n - 2K - 1)$ exact under Gaussianity |
| Variance retained | $1 - r^{2K} = 1 - ((1-\omega)/(1+\omega))^{2K}$ |

The K-harmonic LRT (`detect_FMM` in this codebase) is therefore the
**linearized, untied rhythmicity test motivated by the truncated Fourier
expansion of the FMM model**. At $K = 1$ it reduces to standard cosinor
(DCP). $K = 2$ is the recommended default — it captures $\geq 80\%$ of
FMM variance for empirical $\omega \geq 0.2$ and matches the FMM
parameter budget at exact F-distribution calibration.

---

## 8. References

- Rueda, C., Larriba, Y., & Peddada, S. D. (2019). Frequency Modulated
  Möbius Model for the Estimation of Rhythmic Signals. *Scientific
  Reports*, 9, 18138.
- Davies, R. B. (1977). Hypothesis testing when a nuisance parameter is
  present only under the alternative. *Biometrika*, 64(2), 247–254.
- Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is
  present only under the alternative. *Biometrika*, 74(1), 33–43.
- Hughes, M. E., et al. (2017). Guidelines for genome-scale analysis of
  biological rhythms. *Journal of Biological Rhythms*, 32(5), 380–393.
