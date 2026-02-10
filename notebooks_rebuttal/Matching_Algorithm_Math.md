# Gene Matching Algorithm: Mathematical Description

## Overview

We developed a percentile-based kernel-weighted matching algorithm to match genes based on multiple biological variables simultaneously. The algorithm matches each target gene to a control gene from a candidate pool based on similarity in CDS length, BrainSpan whole-brain expression (WB), and LOEUF (Loss-of-Function Observed/Expected Upper Fraction) score.

## Notation

Let $\mathcal{G} = \{g_1, g_2, \ldots, g_N\}$ denote the set of all genes with complete data for all matching variables.

For each gene $g_i \in \mathcal{G}$, we have three matching variables:
- $x_{i,1}$: CDS length
- $x_{i,2}$: BrainSpan whole-brain expression (WB)
- $x_{i,3}$: LOEUF score

Let $\mathcal{T} = \{t_1, t_2, \ldots, t_M\}$ denote the set of target genes to be matched, and $\mathcal{C} = \mathcal{G} \setminus \mathcal{T}$ denote the candidate pool of control genes.

## Step 1: Percentile Conversion

To normalize across variables with different scales and distributions, we convert each variable to percentile ranks (0-100 scale). For variable $j$ and gene $g_i$:

$$P_{i,j} = \frac{\text{rank}(x_{i,j})}{|\mathcal{G}|} \times 100$$

where $\text{rank}(x_{i,j})$ is the rank of $x_{i,j}$ among all values of variable $j$ in $\mathcal{G}$ (ties are handled by averaging ranks).

The percentile vector for gene $g_i$ is:
$$\mathbf{P}_i = [P_{i,1}, P_{i,2}, P_{i,3}]$$

## Step 2: Distance Calculation

For each target gene $t \in \mathcal{T}$ and candidate gene $c \in \mathcal{C}$, we calculate the Euclidean distance in the 3-dimensional percentile space:

$$d(t, c) = \sqrt{\sum_{j=1}^{3} (P_{t,j} - P_{c,j})^2}$$

This metric simultaneously considers all three variables, ensuring that matched genes are similar across CDS length, WB expression, and LOEUF score.

## Step 3: Kernel Weighting

To assign weights to candidate genes based on their distance from the target gene, we apply a kernel function. Two kernel types are implemented:

### Uniform Kernel

The uniform kernel assigns equal weight to all candidates within the bandwidth $h$:

$$w_{\text{uniform}}(d) = \begin{cases}
1 & \text{if } d \leq h \\
0 & \text{if } d > h
\end{cases}$$

### Tricubic Kernel

The tricubic kernel provides smooth, distance-weighted probabilities that decrease with increasing distance:

$$w_{\text{tricubic}}(d) = \begin{cases}
\left(1 - \left(\frac{d}{h}\right)^3\right)^3 & \text{if } d < h \\
0 & \text{if } d \geq h
\end{cases}$$

The tricubic kernel assigns higher weights to closer matches while smoothly decreasing weights for more distant candidates, ensuring the most similar genes have the highest probability of being selected.

## Step 4: Probability Normalization and Sampling

For target gene $t$, we compute weights for all candidates:

$$w_t(c) = w(d(t, c)) \quad \forall c \in \mathcal{C}$$

where $w$ is either $w_{\text{uniform}}$ or $w_{\text{tricubic}}$ depending on the chosen kernel.

We filter to candidates with non-zero weight:
$$\mathcal{C}_t^* = \{c \in \mathcal{C} : w_t(c) > 0\}$$

If $|\mathcal{C}_t^*| < k_{\min}$ (minimum required candidates, default $k_{\min} = 10$), we expand the bandwidth to $h' = 2h$ and recompute weights.

The weights are normalized to probabilities:

$$p_t(c) = \frac{w_t(c)}{\sum_{c' \in \mathcal{C}_t^*} w_t(c')} \quad \forall c \in \mathcal{C}_t^*$$

Finally, we sample one matched gene $m(t)$ from $\mathcal{C}_t^*$ according to these probabilities:

$$m(t) \sim \text{Categorical}(\mathbf{p}_t)$$

where $\mathbf{p}_t = [p_t(c_1), p_t(c_2), \ldots, p_t(c_{|\mathcal{C}_t^*|})]$ is the probability vector over $\mathcal{C}_t^*$.

## Step 5: Matching with Replacement Prevention

To prevent duplicate matches, once a candidate gene $c$ is matched to a target gene, it is removed from the candidate pool for subsequent matches:

$$\mathcal{C}^{(k+1)} = \mathcal{C}^{(k)} \setminus \{m(t_k)\}$$

where $\mathcal{C}^{(k)}$ is the candidate pool at iteration $k$ and $t_k$ is the $k$-th target gene to be matched.

## Algorithm Parameters

- **Bandwidth ($h$)**: Controls the similarity threshold. Default: $h = 10$ percentile units
- **Kernel type**: `uniform` or `tricubic`. Default: `tricubic`
- **Minimum candidates ($k_{\min}$)**: Minimum number of candidates required. Default: $k_{\min} = 10$

## Advantages of This Approach

1. **Simultaneous multi-variable matching**: The Euclidean distance in percentile space ensures matches are similar across all three variables simultaneously, preventing cases where one variable matches perfectly while others diverge.

2. **Percentile normalization**: By converting to percentiles, we account for different variable scales and distributions, ensuring each variable contributes equally to the distance metric.

3. **Probabilistic sampling**: Rather than selecting the closest match deterministically, weighted sampling allows for variability while still favoring closer matches.

4. **Smooth weighting (tricubic kernel)**: The tricubic kernel provides a smooth decay in matching probability with distance, giving preference to very similar genes while still allowing some flexibility for slightly more distant matches.

## Implementation Notes

- If a target gene has insufficient candidates after bandwidth expansion, it is excluded from the final matched set
- The random seed can be set for reproducibility
- Percentile tables are pre-computed and cached for computational efficiency when matching multiple gene sets








