
### Generalized Wick's Theorem

$$ \ket{\phi_0} = e^T \ket{\phi_0}$$
$$ \{A\}\{B\} = \{AB\} + \{ \overline{AB} \} $$

$$ T = T_1 + T_2 + ... $$ where  $$T_i = \frac{1}{(i!)^2} \sum_{ij ab}t_{ab}^{ij} \{ a_a^\dagger a_b^\dagger .. a_i a_j\}$$ 

### CCD approximation

$$ T = \frac{1}{4} \sum_{ijab} t_{ij}^{ab}\{a_a^\dagger a_b^\dagger a_i a_j\}$$

$$\bar{H}_N \ket{\phi_0}  = E\ket{\phi_0}$$

$$ \bra{\phi_0} \bar{H}_N \ket{\phi_0} = E$$
$$ \bra{\phi_{ij}^{ab}} \bar{H}_N \ket{\phi_0} = 0$$
where,
$$ \bar{H}_N = H_N + \overline{H_NT_z} + \frac{1}{2} \overline{H_N T_2} \overline{T_2}$$

$$ H_N =  \underbrace{\sum_{p,q} f_q^p \{p^\dagger q\}}_{F} +  \underbrace{\frac{1}{4}\sum_{pqrs} \{pq|v|rs\} \{p^\dagger q^\dagger r s\} }_{V}$$

$$ \bra{\phi_0} V T_2 \ket{\phi_0} = \frac{1}{4} \frac{1}{4} \sum_{pqrs\\ abij} t_{ij}^{ab} \{pq|v|rs\} \{ p^\dagger q^\dagger rs\} \{ a^\dagger b^\dagger ij \} = \frac{1}{4} \sum_{a,b,i,j} t^{ab}_{ij} \{ ij | v| ab\}$$



Let's do some finger exercise to get acquainted ourselves, by calculating a part of the  CCD equation

$$ \bra{\phi_{ij}^{ab}} FT_2 \ket{\phi_0} = \frac{1}{4} \sum_{p,q,c,d,k,l} f_q^p t_{kl}^{cd} \{ i^\dagger j^\dagger ba\} \{p^\dagger q\} \{ c^\dagger d^\dagger l k\} $$
 $$ \{ i^\dagger j^\dagger ba\} \{ p^\dagger q\} \{c^\dagger d^\dagger l k\}$$ which is a contraction with particle states, and $$ .. $$ is a contraction with hole states


and the final result is $$ \frac{1}{4} \sum_c \left( f_c^a t_{ij}^{cb} + f_c^b t_{ij}^{ac} \right) - \frac{1}{4} \left( \sum_k \left( f_j^k t_{ik}^{ab} + f_i^k t_{kj}^{ab}\right)\right)$$
 
 $$ P(ab ) = 1- \underbrace{\hat{P}(ab)}_{\textrm{permutation}}$$

$$ \bra{\phi_{ij}^{ab}} FT_2 \ket{\phi_0} =  P(ab) \sum_c f_c^b t^{ac}_{ij} - P(ij) \sum_k f_j^k t_{ik}^{ab}$$
$$ 0 = \bra{ab} v\ket{ij} + P(ab) \sum_c f_c^b t_{ij}^{ac} - P(ij) \sum_k f_j^k t_{ik}^{ab} + \frac{1}{2}\sum_{cd} \bra{ab}V\ket{cd} t_{ij}^{cd} + \frac{1}{2}\sum_{kl} \bra{kl} v\ket{ij} t_{kl}^{ab} + P(ab)P(ij)\sum_{k,c} \bra{kb}V\ket{cj}t^{ac}_{ik}  + .. + \frac{1}{4}\bra{kl}V\ket{cd} t_{ij}^{cd} t_{kl}^{ab}$$


