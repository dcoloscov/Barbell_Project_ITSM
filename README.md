# Random walks on barbell graphs: mixing and "echo chambers"

Individual simulation project for *Introduction to Stochastic Methods* (USI, Autumn 2026).

**Model.** Two communities, each a complete graph $K_k$, joined by a single bridge edge.
A message starts at an interior node $v_0$ and is passed at each step to a uniformly
random neighbor (simple random walk). [Picture](images/barbell.pdf)


**Questions.**
1. How long until the message reaches the other community (hitting time of the bridge)?
2. How fast does the walk settle to its stationary distribution (mixing time)?

**Results.**
- [RESULT 1, e.g. how the crossing time scales with k, with simulation vs. theory](poster/images/cross_plot.pdf)
- [RESULT 2, e.g. mixing time / spectral gap behaviour](poster/images/dist_plot.pdf)


Full poster: [poster/poster.pdf](poster/poster.pdf)

## Reproduce

All parameters are set at the top of `code/main.py`.

```bash
pip install -r code/requirements.txt
python code/main.py
```

This regenerates all data and figures used in the poster.
