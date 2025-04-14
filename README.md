# Agent-Based Model of Protein Aggregation under Macromolecular Crowding

## Species Key
- **N**: Native monomer (non-aggregating)
- **A**: Amyloid-prone monomer (aggregation-competent)
- **O**: Oligomer (dimer of A monomers, A₂)
- **F**: Fibril (initiated by O + A, elongated by A + F)

---

## Project Overview
This Agent-Based Model (ABM) simulates early-stage protein aggregation under macromolecular crowding, focusing on the nucleation process. The model includes transitions between native (N) and amyloid-prone (A) monomers and tracks the formation of small oligomers and initial fibrils in a spatially explicit environment.

## Key Features
- N ⇄ A monomer conformational transitions
- A + A ⇄ Oligomer formation (dimers only)
- Oligomers (A₂) can dissociate back into A monomers
- Fibril nucleation: A₂ + A → Fibril (F) initiation
- Fibrils can grow through A + F → F (elongation)
- Spherical crowders influence spatial dynamics and diffusion
- Output: number of monomers, oligomers, and fibrils over time

##Reaction Summary
