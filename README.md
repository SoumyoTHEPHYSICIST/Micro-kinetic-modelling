# Micro-kinetic-modelling
# Microkinetic Modeling of Methane Pyrolysis (with ZnCl₂ Catalyst)

Welcome!   
This project dives into the **reaction kinetics of methane (CH₄) pyrolysis**, where methane breaks down into **hydrogen and carbon** — with the help of a **gas-phase ZnCl₂ catalyst** at high temperatures.

It’s part of my computational assignment for **CHE631 (Chemical Reaction Engineering)** at IIT Kanpur.

---

## What's This About?

Methane pyrolysis is a cleaner way to make hydrogen without CO₂ emissions — unlike steam reforming. But, it has its own challenges like slow reaction rates and carbon deposition. So we decided to explore how things behave when we throw in **ZnCl₂ as a gas-phase catalyst** and simulate the whole system.

---

##  What I Did

Here’s the breakdown of the project:

- Modeled **24 chemical reactions** (both unimolecular and bimolecular).
- Calculated rate constants using **Transition State Theory**.
- Simulated how **methane concentration changes over time** during pyrolysis.
- Used a sensitivity tool called **Degree of Reaction Control (DRC)** to figure out which reactions are actually important (aka the rate-limiting steps).
- Plotted everything to get a better visual understanding.

---

## The Setup

- **Temperature**: 1273 K  
- **Pressure**: 1 atm  
- **Catalyst**: Zinc Chloride (ZnCl₂)  
- **Gas Mix**: CH₄ + Argon (with ZnCl₂ vapor)

We assume an ideal batch reactor system and model the kinetics using rate equations and stoichiometry.

---

##  Key Results

### CH₄ Degradation Curve
- Methane breaks down exponentially with time.
- Expected behavior: as methane is the primary reactant, it decreases fast.

###  Degree of Reaction Control Plot
- Identified which steps are actually slowing the process (rate-limiting).
- Found that some early reactions dominate initially, but others (like steps 11, 13, and 18) gain importance later as stable products form.

###  Time-Based Behavior
- Some intermediate species disappear quickly — they’re unstable and don’t hang around.
- Final products form gradually and stabilize.

---

##  Takeaways

- ZnCl₂ is an effective gas-phase catalyst for methane pyrolysis.
- DRC is super useful — it highlights which steps are worth focusing on and saves computation time.
- This kind of modeling helps in designing cleaner hydrogen production methods.

---


