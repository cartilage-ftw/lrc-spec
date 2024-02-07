## Possibilities with Yb II
Yb II has a simpler, less crowded electronic structure and after Lu, it might be the simplest lanthanide element we can reach with LRC. With little literature reporting Yb II and no reported hyperfine splitting constants $^\dagger$ we will also be able to fill gaps in atomic data.

$^\dagger$ To my knowledge


The following ground state E1 transitions of Yb II has branching fractions (and absolute transition probabilities) reported by [Kedzierski _et al._ (2010)](https://www.sciencedirect.com/science/article/pii/S0584854710000522). Lifetime measurements come from several studies.

Yb II has a ground state configuration of $4f^{14}6s$ (term symbol $^2S_{1/2}$). I have listed below which Sirah dye would be needed to probe these with our existing laser and frequency doubling arrangement.

| $E_u$ [cm-1] | configuration | Term | $\tau$ [ns] | $\lambda$ [nm] | Sirah dye(s) |
| ---- | ---- | ---- | ---- | ---- | ---- |
| 27,062 | $4f^{14}6p$ | $^2P^o_{1/2}$ | 8.07±0.09 | 369.42 | Pyridine 2, Styryl 8 |
| 30,392 | $4f^{13}6p$ | $^2P^o_{3/2}$ | 6.15±0.09 | 328.94 | DCM (in DMSO) |
| Transitions listed below didn't work |  |  |  |  |  |
| 32,982 | $4f^{13}5d6s$ | $[^2F^o_{7/2}\ ^3D]_{3/2}$ | 196±20 | 303.11 | Rhodamine B (or mixture with Rhod 101) |
| 33,654 | $4f^{13}5d6s$ | $[^2F^o_{7/2}\ ^1D]_{1/2}$ | 39±3 | 297.05 | Rhodamine B (or mix with Rhod 101) |
| 34,575 | $4f^{13}5d6s$ | $[^2F^o_{7/2}\ ^3D]_{3/2}$ | 28.6±0.4 | 289.14 | Pyrromethene 597, Pyrromethene 580 |

**Caution**: Some of the configurations, wavenumbers, and even $J$ values for the transitions involving $4f^{13}5d6s$ are either not on NIST, or seem to disagree with their values. The $J=1/2$ level at $34,575$ cm-1 has been assigned a $[^2F^o_{7/2}\ ^1D]$ instead of $^3D$

In the [DREAM database](https://agif.umons.ac.be/databases/dream.html), maintained by the group at Mons University in Belgium, there are measurements of Yb II transition strengths, which we can use to estimate the degree of branching to the metastable state(s)

| Transition term symbol | $E_L$ [cm-1] | $E_U$ [cm-1] | $gA$ [s-1] |
| --- | --- | --- | --- |
| $^2S_{1/2} \to$ $^2P^o_{1/2}$ | $0$ | $27062$ | $2.31\times 10^8$ |
| $^2P^o_{1/2} \to$ $^2D_{3/2}$ | $22961$ | $27062$ | $1.13\times 10^6$ |
| $^2P^o_{3/2} \to$ $^2D_{3/2}$ | $22961$ | $30392$ | $1.15\times 10^6$ |
| $^2P^o_{3/2} \to$ $^2D_{5/2} $| $24333$ | $30392$ | $5.64\times 10^6$ |

They do not report any measurement of the potential $^2D_{5/2}\to$ $^2F^o_{7/2}$ transition. If the term symbols and $J$ values are correct, this transition shouldn't be forbidden due to angular momentum or parity.
Perhaps what is true is that, this transition would require two electrons to jump, instead of just one, and that's not allowed.


### Charge Radii and Nuclear Shape Deformations from Isotope Chain
Yb has 8 stable isotopes: two odd $^{171}$ Yb (14.2\%) and $^{173}$ Yb (16.1\%), and six even $^{168,170,172,174,176}$ Yb. The odd isotopes have ground state nuclear spins $I_{^{171}\textrm{Yb}}=1/2^-$ and $I_{^{173}\textrm{Yb}} = 5/2^-$. This spin, along with the non-zero $J$ of the ground and excited states will give rise to a rich hyperfine spectrum. Our ability to separate out different isotopes in-flight with the QMS will allow us to study the different isotopes separately. By combining with atomic theory calculations, we will be able to extract the charge radii and nuclear shape deformations of the different isotopes along a broad sequence without the need of going to an online facility.