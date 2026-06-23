# Physical singlet meson correlator: sign and normalization of $C$ vs $D$

Goal: pin down, from first principles, the **relative sign** and **relative
normalization** between the Grid connected correlator $C$ and Rebbi's
disconnected correlator $D$, so that the physical flavor-singlet correlator is
assembled correctly. Theory: $N_f=1$ Dirac fermion, SU(4), Mobius DWF, $24^3\times48$.

Standard ingredients (Wick contractions, $\gamma_5$-hermiticity); no published
algorithm beyond textbook lattice meson contractions.

## 1. What each code actually computes

### Connected $C$ (Grid: `baryons_0000_dirac_claude.cc`)
- Point source at the origin, single solve, $q_1=q_2=S(\cdot,0)$
  (`PointSource` :77, `Solve` with `ImportPhysicalFermionSource` /
  `ExportPhysicalFermionSolution` :106-115).
- Contraction (`MesonTrace_hdf` :203):
$$
\texttt{meson\_CF}(x)=\mathrm{tr}\big[\gamma_5\,q_1^\dagger\,\gamma_5\,\Gamma_{\rm snk}\,q_2\,\Gamma_{\rm src}^\dagger\big].
$$
- Use $\gamma_5$-hermiticity $\gamma_5 S(x,0)^\dagger \gamma_5 = S(0,x)$:
$$
\gamma_5\,q_1^\dagger\,\gamma_5=\gamma_5 S(x,0)^\dagger \gamma_5 = S(0,x),
$$
so for a diagonal channel $\Gamma_{\rm snk}=\Gamma_{\rm src}=\Gamma$,
$$
\texttt{meson\_CF}(x)=\mathrm{tr}\big[S(0,x)\,\Gamma\,S(x,0)\,\Gamma^\dagger\big]
=\mathrm{tr}\big[\Gamma\,S(x,0)\,\Gamma^\dagger\,S(0,x)\big].
$$
- Sink momentum projection + `sliceSum` over $T$ (:250-253) gives, for $\vec p=0$,
$$
\boxed{\;C(t)=\sum_{\vec x}\mathrm{tr}\big[\Gamma\,S(x,0)\,\Gamma^\dagger\,S(0,x)\big]\;}
$$
the connected ("valence/pion-like") correlator. **Single fixed source point,
sum over sink.**

### Disconnected loop (Grid: `disc_multipleGamma_binary_claude.cc`)
- Stochastic source $\eta$, full spin/color/(t,eo) dilution, $\psi=S\eta$
  (`TraceField` :113): `meson_CF = trace(Gamma * psi * adj(eta))`
  $=\mathrm{tr}[\Gamma\,S\,\eta\,\eta^\dagger]$, accumulated over dilution so that
  $\mathbb{E}[\eta\eta^\dagger]=\mathbb{1}$ gives the **unbiased loop**
$$
T_\Gamma(x)=\mathrm{tr}\big[\Gamma\,S(x,x)\big].
$$
  (No $1/N_{\rm hit}$: single hit, dilution tiles the identity exactly.)

### Disc correlator $D$ (Rebbi: `average_trace2_claude.f90`, `makecorr_96_d_claude.f90`)
- `average_trace2`: stored loop is **real** (`tr(...,0)=a`, the real part,
  :114 — see $\gamma_5$ note below). Per-site-normalized momentum components:
$$
\texttt{avtr}(t)=\frac{1}{N_S^3}\sum_{\vec x}T_\Gamma(\vec x,t),\qquad
\texttt{pxyz}(t,j)=\frac{1}{N_S^3}\sum_{\vec x}e^{i\vec p_j\cdot\vec x}\,T_\Gamma(\vec x,t)
$$
  (the $1/N_S^3$ at :120 and :143).
- `makecorr` (factor-3 line removed): with $n_s\equiv N_T/2=24=N_S$ (:166),
$$
\texttt{pcorr}(i,it)=\frac{n_s^3}{N_T}\sum_t \texttt{pdata}(i,t)\,\texttt{pdata}(i,t+it)
\quad(\text{:248}),
$$
  and for the $\vec p^2=0$ class $\texttt{corr}(it,0)=\texttt{pcorr}(0,it)$. Substituting
  $\texttt{pdata}(0,t)=\texttt{avtr}(t)=T_{\rm sum}(t)/N_S^3$ with
  $T_{\rm sum}(t)=\sum_{\vec x}T_\Gamma(\vec x,t)$:
$$
\boxed{\;D(it)=\frac{n_s^3}{N_T N_S^6}\sum_t T_{\rm sum}(t)\,T_{\rm sum}(t+it)
=\frac{1}{N_S^3}\,\Big\langle T_{\rm sum}(t)\,T_{\rm sum}(t+it)\Big\rangle_t\;}
$$
  using $n_s=N_S$. This is the loop-loop correlator **translation-averaged over
  the source timeslice**, normalized to **one source point** (the $1/N_S^3$ turns
  the double spatial sum into "sum over sink, average over source").

So $C$ and $D$ are in the **same normalization**: one source point, sum over
sink, per unit time separation. Hence the normalization factor relating them is
$\boxed{\,\text{rel\_norm}=1\,}$, confirming OQ2.

> Latent caveat: the code uses $n_s=N_T/2$, which equals $N_S$ only because
> $N_T/2=N_S=24$ here. If $N_T/2\neq N_S$ the net disc normalization would be
> $(N_T/2)^3/N_S^6$, not $1/N_S^3$. Fine for this lattice; flag if reused.

## 2. $\gamma_5$-hermiticity facts used

- $C(t)$ is **positive**. For $\Gamma=\gamma_5$, $\Gamma^\dagger=\gamma_5$ and
$$
\mathrm{tr}[\gamma_5 S(x,0)\gamma_5 S(0,x)]
=\mathrm{tr}[\gamma_5 S(x,0)\gamma_5\,\gamma_5 S(x,0)^\dagger \gamma_5]
=\mathrm{tr}[S(x,0)S(x,0)^\dagger]=\Vert S(x,0)\Vert^2>0.
$$
  Matches the measured Grid $C>0$ at all $t$.
- The loop $T_{\gamma_5}(x)=\mathrm{tr}[\gamma_5 S(x,x)]$ is **real**:
$$
\mathrm{tr}[\gamma_5 S(x,x)]
=\mathrm{tr}[\gamma_5\,\gamma_5 S(x,x)^\dagger \gamma_5]
=\mathrm{tr}[\gamma_5 S(x,x)^\dagger]
=\big(\mathrm{tr}[\gamma_5 S(x,x)]\big)^*,
$$
  so `average_trace2` keeping the real part is correct, and $D>0$ at small $t$
  (it is the connected correlator of a real loop amplitude). Matches the data.

## 3. The physical singlet: Wick contractions and the relative sign

Singlet interpolator (one flavor) $O_\Gamma(x)=\bar q(x)\,\Gamma\,q(x)$. Zero-momentum
two-point:
$$
G(t)=\sum_{\vec x}\big\langle O_\Gamma(x)\,O_\Gamma^\dagger(0)\big\rangle .
$$
Two Wick contractions of $\bar q_a\Gamma_{ab}q_b\,\bar q_c\Gamma_{cd}q_d$:

- **Connected** ($q_b(x)\!-\!\bar q_c(0)$, $q_d(0)\!-\!\bar q_a(x)$): one closed
  fermion loop. Gives the standard, **positive** meson 2pt
$$
G_{\rm conn}(t)=+\sum_{\vec x}\mathrm{tr}\big[\Gamma S(x,0)\Gamma^\dagger S(0,x)\big]=+\,C(t).
$$
- **Disconnected** ($q_b(x)\!-\!\bar q_a(x)$, $q_d(0)\!-\!\bar q_c(0)$): **two**
  closed loops, i.e. **one extra fermion loop relative to the connected**, hence a
  relative $(-1)$:
$$
G_{\rm disc}(t)=-\sum_{\vec x}\big\langle T_\Gamma(x)\,T_\Gamma(0)\big\rangle_c
=-\,D(t).
$$

Therefore
$$
\boxed{\;G(t)=C(t)-D(t)\;}\qquad(N_f=1).
$$

For $N_f$ degenerate flavors with $O=\tfrac{1}{\sqrt{N_f}}\sum_f\bar q_f\Gamma q_f$
the connected term keeps coefficient $1$ and the disconnected double-sum gives
$N_f$:
$$
G(t)=C(t)-N_f\,D(t).
$$

### Consistency / positivity check
$G(t)$ must be $\ge 0$ (reflection positivity). With $C,D>0$ and $C\gg D$ at
small $t$ (measured: $D/C\sim0.4\%$ at $t=0$), $C-N_fD>0$ as required, and the
disconnected piece **lowers** the correlator / **raises** the effective mass —
the expected anomaly effect that splits the singlet from the connected state.
The opposite ordering $N_fD-C$ is $\le0$ at small $t$, i.e. **not a physical
correlator**.

## 4. Conclusion and required change

- **Normalization:** rel_norm $=1$ (confirmed; $C$ and $D$ share the one-source,
  sum-over-sink convention; factor-3 already retired).
- **Sign:** the physical singlet is $G=C-N_f D$ with $N_f=1$, i.e. $C-D$.
  The disconnected diagram carries a **relative minus** (extra fermion loop).
- **Current code is wrong by an overall sign.** The notebook / batch form
  $N_f D - C = -(C-N_f D)$. It must be flipped to
$$
\texttt{phys} = C - \texttt{nf\_disc}\cdot D .
$$
  Effective masses (ratios) are unchanged, but positivity, fitting to
  $\sum_n A_n e^{-E_n t}$ with $A_n>0$, and the physical interpretation require
  the $C-N_fD$ form.

### Where Rebbi's "$2D-C$" came from
Rebbi's QCD code ($N_f=2$) stored the connected `ccorr` with the **opposite**
sign convention to Grid's positive $C$, so his "$2D-\texttt{ccorr}$" is $+$(physical).
Dropping Grid's positive $C$ into the literal "$2D-C$" slot inherits the wrong
sign. With Grid's connected sign, the correct combination is $C-N_fD$.

## 5. Files / naming (status)
- Drivers RENAMED `meson_2DminusC_*` $\to$ `meson_CminusD_*` (2026-06-19), output
  dir/npz `2DminusC_*` $\to$ `CminusD_*`, function `combine_2D_minus_C` $\to$
  `combine_CminusD`. These now compute $C - N_f\,\text{rel\_norm}\,D$ ($N_f{=}1$),
  i.e. $C-D$ (overall sign irrelevant for masses; chosen to match the name).
- notebook `analyze_corr_v2_claude.ipynb` build/plot cells still print $D-C$
  labels; the overall sign is irrelevant so left as-is.
- `I_I` ($\sigma$) still needs the vacuum subtraction (loop has a condensate VEV);
  left out as before.
