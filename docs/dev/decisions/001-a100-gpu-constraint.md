# ADR-001: A100 GPU constraint for AF2 2.3.2

Date: 2026-03-16
Status: superseded (2026-04-02)

## Context

AF2 2.3.2 on Minerva uses CUDA kernels incompatible with H100 (compute capability 9.0). Jobs fail silently or produce garbage on H100 nodes.

## Original Decision (2026-03-16)

All AF2 LSF scripts constrain to A100 via `-R "rusage[ngpus_physical=1] select[ngpus_physical>0]" -q gpu -m "gpu-a100"`.

## Updated Decision (2026-04-02)

Allow both V100 and A100. Exclude H100 only. Use `-R "rusage[ngpus_physical=1] select[v100 || a100]"`.

**Evidence:** 181 d1d3 jobs completed successfully on V100 nodes (lg03a* pool) with `-R v100`. Results validated against VEGFA positive control (ipTM 0.816). V100 pool is larger, reducing queue wait times.

## Consequence

Wider GPU pool, shorter queue times. H100 still excluded.
