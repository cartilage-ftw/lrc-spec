# Element Species of Interest


| Species | $E_u$ [cm-1] | Ground state config. | Transition Term Labels | Metastable Level [cm-1] | Metastable config | $\lambda_0$ (nm) | Sirah dye | Comments |
| --- | --- | --- | --- | --- | ---| --- | --- | --- | 
Sr II | $23 715.19$ | $4p^6 4s$ | $^2S_{1/2}\to\ ^2P^o_{1/2}$ | $14 555.90$ | $4p^6 4d$ | $421.6$ $nm$ | Styryl 9 |
| | $24 516.65$ | $4p^6 4s$ | $^2S_{1/2}\to\ ^2P^o_{3/2}$ | $14 555.90$, $14 836.24$ | $4p^6 4d$ | $407.9$ $nm$ | Styryl 9, Styryl 11 |
Ba II | $20 261.56$ | $5p^6 5s$ | $^2S_{1/2}\to\ ^2P^o_{1/2}$ | $4873.852$ | $5p^6 5d$| $493.6$ $nm$ (!!)| not covered | Hyperfine constants have been studied and reported before
| | $21 952.40$ |  | $^2S_{1/2}\to\ ^2P^o_{3/2}$ | $4 873.852$, $5 674.807$ |  | $455.5$ $nm$ (!!) | not covered|same as above |
Yb II | $27062$ | $4f^{14}6s$ | $^2S_{1/2} \to$ $^2P^o_{1/2}$ | $22961$ | $4f^{14}5d$ | $369.5$ $nm$ |  Styryl 8 | Hyperfine constants have been studied and reported before
| | $30392$ | $4f^{14}6s$ | $^2S_{1/2} \to\ ^2P^o_{3/2}$ | $22961$, $24332.69$ | $4f^{14}5d$ | $329.0$ $nm$ |  DCM| 
Hf II | $29160.04$ | $5d 6s^2$  |  $^2D_{3/2}\to\ ^4F^o_{5/2}$ |  $3,644$ | $5d^2 6s$  | $339.9$ $nm$ | Pyridine 1 |
Ta II | $33,706.464$ | $5d^3 6s$ | designation differs between NIST and [Messnarz & Guthöhrlein](https://iopscience.iop.org/article/10.1238/Physica.Regular.067a00059/meta) | $3,180-9,690$ [multiple] | $5d^2 6s^2$ | $296.7$ nm | Rhodamine B | 
Er II | They could not resolve all isotopes; hyperfine constants of the odd isotopes couldn't be measured
Ho II |
Ca II |
Sc II |

A few other interesting targets are a bit outside our wavelength coverage
include Ru II and Nb II ($37,298$ cm-1 ; $4d^4$ ground state ; $^5D_0 \to 5F^o_1$ ; $4d^3 5p$  metastable)

## Remarks on Ta II
The branching fraction for the Ta II transition looks promising (see below). What we may have to worry about are the lower metastable levels. Would they survive the collisions in the drift tube? We do not know until we try.
![](Tantalum_BF.png)

Table from [Quienet _et al_. (2009)](https://www.aanda.org/articles/aa/abs/2009/02/aa11035-08/aa11035-08.html)

[https://www.sciencedirect.com/science/article/pii/S0022407315302764](https://www.sciencedirect.com/science/article/pii/S0022407315302764) provide HFS constants for most Ta II transitions, including this one that I mentioned. [Messnarz & Guthöhrlein (2003)](https://iopscience.iop.org/article/10.1238/Physica.Regular.067a00059/meta) reported HFS constants but didn't study this transition 

## Remarks on Yb II
In the [DREAM database](https://agif.umons.ac.be/databases/dream.html), maintained by the group at Mons University in Belgium, there are measurements of Yb II transition strengths, which we can use to estimate the degree of branching to the metastable state(s)

| Transition term symbol | $E_L$ [cm-1] | $E_U$ [cm-1] | $gA$ [s-1] |
| --- | --- | --- | --- |
| $^2S_{1/2} \to$ $^2P^o_{1/2}$ | $0$ | $27062$ | $2.31\times 10^8$ |
| $^2P^o_{1/2} \to$ $^2D_{3/2}$ | $22961$ | $27062$ | $1.13\times 10^6$ |
| $^2P^o_{3/2} \to$ $^2D_{3/2}$ | $22961$ | $30392$ | $1.15\times 10^6$ |
| $^2P^o_{3/2} \to\ ^2D_{5/2}$ | $24333$ | $30392$ | $5.64\times 10^6$ |

## Remarks on Er II
The prospects of doing Er II are similar to that for Yb II. The branching to metastable states, as far as I have been able to figure out, would be similar. However, Er II has less atomic data known from prior measurements.

## Remarks on Hf II
The branching to the metastable state for the transitions mentioned is substantial. But a big concern is that the lowest metastable level at $3,644$ cm-1 is only a few hundred $cm^{-1}$ above the ground state fine structure sublevel at $3,050$. Refer to [Lawler et al (2009)](https://iopscience.iop.org/article/10.1086/510368/pdf) if searching for measured $A_{ki}$ values.
