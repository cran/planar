Decay rates of a dipole near a planar multilayer stack
========================================================
<!-- 
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{decay_rates}
-->

Following (Le Ru, Etchegoin, p. 571), and (Novotny, Hecht pp. 335--360), the enhancement factor for the total decay rate for a dipole perpendicular to the interface is

$$
  M^\perp_\text{tot} = 1 + \frac 3 2 \int_0^\infty \Re\left\{ \frac{q^3}{\sqrt{1 - q^2}} r^p (q) \exp\left(2 i k_1 d\sqrt{1 - q^2} \right)\right\}\text{d}q 
$$

The integrand diverges as $q\to 1$, it is therefore advantageous to perform the substitution $u:=\sqrt{1 - q^2}$. In order to maintain a real path of integration, the integral is first split into a radiative region ($0\leq q\leq 1$, $u:=\sqrt{1 - q^2}\geq 0$), and an evanescent region ($1\leq q\leq \infty$, $-i u:=\sqrt{q^2 - 1}\geq 0$). After some algebraic manipulation, we obtain,
$$
  M^\perp_\text{tot} = 1 + \frac 3 2 \left(I_1 + I_2\right)
$$
where
$$
  \begin{split}
I_1 + I_2 =& \int_0^1 \left[1 - u^2\right]\cdot\Re\left\{ r^p (\sqrt{1 - u^2}) \exp\left(2 id k_1 u\right)\right\}\text{d}u \\ & +\int_0^\infty \left[1 + u^2\right]\cdot\exp\left(-2d k_1  u\right)\cdot\Im\left\{ r^p (\sqrt{1 + u^2}) \right\}\text{d}u 
  \end{split}
$$
Similarly, for the parallel dipole
$$
  M^\parallel_\text{tot} = 1 + \frac 3 4 \int_0^\infty \Re\left\{ \left[ \frac{r^s(q)}{\sqrt{1 - q^2}} - r^p (q) \sqrt{1 - q^2}\right] \cdot q\cdot\exp\left(2 i k_1 d\sqrt{1 - q^2} \right)\right\}\text{d}q  \label{eq:Mstot1}
$$
which can be rewritten as,
$$
  M^\parallel_\text{tot} = 1 + \frac 3 4 \left(I^\parallel_1 + I^\parallel_2\right)
$$
where
$$
  \begin{split}
I^\parallel_1 + I^\parallel_2 =& \int_0^1 \Re\left\{\left[r^s (\sqrt{1 - u^2})  - u^2\cdot r^p (\sqrt{1 - u^2}) \right] \exp\left(2 id k_1 u\right)\right\}\text{d}u \\ & 
+\int_0^\infty \exp\left(-2d k_1  u\right)\cdot \Im\left\{r^s (\sqrt{1 + u^2}) + u^2\cdot r^p (\sqrt{1 + u^2}) \right\}\text{d}u 
  \end{split}
$$

