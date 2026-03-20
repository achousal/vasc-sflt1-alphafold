# ADR-001: A100 GPU constraint for AF2 2.3.2

Date: 2026-03-16
Status: accepted

## Context

AF2 2.3.2 on Minerva uses CUDA kernels incompatible with H100 (compute capability 9.0). Jobs fail silently or produce garbage on H100 nodes.

## Decision

All AF2 LSF scripts constrain to A100 via `-R "rusage[ngpus_physical=1] select[ngpus_physical>0]" -q gpu -m "gpu-a100"`.

## Consequence

Queue wait times are longer (A100 pool is smaller than general GPU pool). Large jobs (>2000 residues) need 72-96h walltime. Acceptable tradeoff for correctness.
