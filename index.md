---
layout: dis_template
---

The xyz package implements the [xyz](https://arxiv.org/abs/1610.05108) algorithm for fast interaction search in high-dimensional data.

# Simple problem

Given a data matrix \\( X \in \mathbb{R}^{n\times p} \\) and response vector \\( \mathbf{Y} \in \mathbb{R}^n \\), we want to find strong product interaction. That is

$$argmin_{(j,k) \in \{1,...,p\}^2} |\sum_{i=1}^n Y_i X_{ij}X_{ik}|$$
