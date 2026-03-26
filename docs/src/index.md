# FAIR_Implementation_ABM

Welcome to the documentation for `FAIR_Implementation_ABM`, a Julia-based agent-based model designed to simulate protein aggregation dynamics under intracellular-like conditions.

This model investigates how different environmental and molecular parameters—such as macromolecular crowding and conformational changes—influence the transition from native monomers to pathological aggregate-prone structures.

## Key Features

- Supports pathological aggregation scenarios
- Includes the following states:
  - **N**: Native monomer
  - **A**: Aggregate-prone monomer
  - **O**: Oligomer
  - **F**: Fibril
- Integrates probabilistic reactions, movement logic, and crowding effects
- Built with an emphasis on **FAIR** (Findable, Accessible, Interoperable, Reusable) model design

## API Reference

```@autodocs
Modules = [FAIR_Implementation_ABM]

```

## Model Parameters

This section summarizes the main parameters used to control the agent-based model of protein aggregation.

| Parameter                          | Description                                                    | Example Value |
|-----------------------------------|----------------------------------------------------------------|---------------|
| `Lattice_Size`                    | Size of the 3D FCC lattice                                     | 30            |
| `MAX_NumberMovements`             | Number of time steps                                           | 1000          |
| `Max_NumberMonomers_Native`       | Initial count of native monomers (N)                           | 50            |
| `Max_NumberMonomers_AggregateProne` | Initial count of aggregate-prone monomers (A)                | 20            |
| `Native_to_AggregateProne`        | Probability for N → A transition per timestep                  | 0.1           |
| `AggregateProne_to_Native`        | Probability for A → N transition per timestep                  | 0.05          |
| `Oligomer_Formation`              | Probability of A + A → O                                       | 0.5           |
| `Oligomer_Dissociation_rate`      | Probability of O → A + A                                       | 0.1           |
| `Fibril_Formation`                | Probability of O + A → F                                       | 0.8           |
| `Fibril_Growth`                   | Probability of F + A → F(n) (elongation)                       | 0.9           |
| `Crowder_Concentration_Spheres`   | Volume fraction of crowders in the lattice                     | 0.4           |
| `Obstacle_Radius`                 | Radius of each spherical crowder; valid integer values are 0–6 | 1             |

---

## Reaction Mechanisms

This model captures a simplified kinetic representation of protein aggregation pathways, driven by stochastic rules encoded in the agent-based framework.

### 1. Conformational Change
**Reaction:** `N ⇌ A`  
Native monomers can spontaneously convert into aggregate-prone conformations and vice versa, controlled by the probabilities `Native_to_AggregateProne` and `AggregateProne_to_Native`.

### 2. Oligomerization
**Reaction:** `A + A ⇌ O`  
Two aggregate-prone monomers can associate to form a small oligomer. Oligomers can also dissociate into monomers.

### 3. Fibril Nucleation
**Reaction:** `O + A → F`  
Oligomers can convert into fibrils when interacting with an additional aggregate-prone monomer, initiating fibril formation.

### 4. Fibril Elongation
**Reaction:** `F + A → F(n)`  
Fibrils grow through sequential addition of aggregate-prone monomers, representing elongation.

Each of these reactions is probabilistic and depends on local lattice arrangement and spatial proximity of monomers.