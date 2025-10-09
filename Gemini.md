# An Electrochemical Review: Kinetic Analysis of the [Fe(CN)₆]³⁻/⁴⁻ Redox Couple and Performance Benchmarking of Graphite|NMC-811 Lithium-Ion Batteries

## Executive Summary

This report provides a comprehensive, expert-level analysis of two distinct yet fundamentally important electrochemical systems. The first part of the report focuses on the aqueous hexacyanoferrate(III)/(II) redox couple, [Fe(CN)₆]³⁻/⁴⁻, a model system used to benchmark electrochemical kinetics. A detailed review of the literature is presented to establish consensus values for its key kinetic parameters, including the diffusion coefficient (D) and the standard heterogeneous rate constant (k₀). The analysis confirms the system's classification as electrochemically reversible, while highlighting its profound sensitivity to electrode surface conditions, making it an invaluable diagnostic tool. Key literature values for the diffusion coefficient typically range from 6.3 to 9.1 × 10⁻⁶ cm²/s.

The second part of the report transitions to the complex, applied system of a high-energy lithium-ion battery, specifically the Graphite|LiNi₀.₈Co₀.₁Mn₀.₁O₂ (NMC-811) cell chemistry. The fundamental properties of the graphite anode, the nickel-rich NMC-811 cathode, and the standard LP30-type carbonate electrolyte are reviewed in detail. A synthesis of recent literature (2020–2024) is conducted to benchmark the cell's performance, focusing on Coulombic efficiency, capacity retention, and cycling stability under various operating conditions. The analysis reveals a critical trade-off between accessing the high energy density of NMC-811 through high-voltage cycling (up to 4.5 V) and maintaining long-term stability. Advanced strategies, such as cathode coatings and the use of single-crystal materials, demonstrate significant improvements, with capacity retentions reaching 88% after 500 cycles. A central theme that emerges is the role of electrode "crosstalk"—whereby transition metal dissolution from the cathode leads to accelerated degradation of the anode's solid electrolyte interphase—as a primary driver of cell failure.

Finally, the report provides a critical evaluation of the electrochemical characterization techniques central to these analyses. The complementary roles of voltammetric (e.g., Cyclic Voltammetry), galvanostatic (e.g., Charge-Discharge Cycling), and impedimetric (e.g., Electrochemical Impedance Spectroscopy) methods are elucidated. A synergistic workflow is proposed, wherein galvanostatic cycling is used for long-term performance monitoring, its derivatives (e.g., dQ/dV) serve as sensitive state-of-health indicators, and impedance spectroscopy provides the in-depth diagnostics required to deconstruct degradation mechanisms. This multi-technique approach is presented as essential for advancing the understanding and development of next-generation energy storage systems.

## Section 1: Theoretical Foundations of Electrochemical Characterization

### 1.1. Principles of Electrochemical Impedance Spectroscopy (EIS)

Electrochemical Impedance Spectroscopy (EIS) is a powerful, non-destructive technique used to investigate the properties of electrochemical systems by examining their response to a small-amplitude alternating current (AC) signal.¹ The core principle involves applying a sinusoidal potential (or current) perturbation, E(t) = E₀ sin(ωt), to an electrochemical cell and measuring the resulting sinusoidal current (or potential) response, I(t) = I₀ sin(ωt + φ).³ The relationship between the applied voltage and the resulting current defines the system's impedance, Z, which is a complex quantity. The impedance is a function of the angular frequency (ω) and is expressed as Z(ω) = |Z|e^(iφ), where |Z| is the magnitude of the impedance and φ is the phase shift between the voltage and current signals.³

In Cartesian coordinates, the impedance is represented as Z(ω) = Z' + iZ'', where Z' is the real (resistive) component and Z'' is the imaginary (reactive) component.⁵ By systematically varying the frequency of the AC signal over a wide range (typically from MHz to mHz), a complete impedance spectrum is generated, which contains detailed information about the various physical and chemical processes occurring within the cell, each with its own characteristic time constant.¹

Two critical conditions must be met for a valid EIS measurement: linearity and steady-state. The system's response must be pseudo-linear, meaning the current response is directly proportional to the voltage perturbation. This is achieved by using a small excitation signal, typically 5-10 mV, so that the system's behavior can be approximated by a linear transfer function.³ The system must also be at a steady state throughout the measurement, which can often take hours. Any drift in the system—caused by factors such as temperature fluctuations, adsorption of impurities, or ongoing degradation—can distort the spectrum and lead to erroneous results.³ The validity of EIS data is often checked using the Kramers-Kronig relations, which provide a mathematical test for linearity and time-invariance.⁶

The complex impedance data are commonly visualized in two formats: the Nyquist plot and Bode plots. The Nyquist plot graphs the negative of the imaginary component (-Z'') against the real component (Z') for each frequency. This format provides an intuitive visual representation of the electrochemical processes; for instance, a semicircle often corresponds to a parallel resistor-capacitor combination, while a straight line can indicate diffusion.² Bode plots display the impedance magnitude (|Z|) and the phase angle (φ) as a function of log(frequency). These plots are particularly useful for analyzing the frequency dependence of the system and identifying the characteristic time constants of different processes.³

### 1.2. Equivalent Circuit Modeling in EIS

The most common method for interpreting EIS data is equivalent circuit modeling. This approach involves representing the complex electrochemical processes occurring at the electrode-electrolyte interface with an idealized electrical circuit composed of discrete elements.¹ The goal is to construct a model that is both physically meaningful and provides a good mathematical fit to the experimental data. The total impedance of the equivalent circuit is calculated based on the series and parallel combinations of its constituent elements.⁸

The fundamental elements used in these models and their physical analogues are:

**Resistor (R)**: A resistor represents processes that dissipate energy, typically by impeding the flow of charge (either ions or electrons). Its impedance is real, frequency-independent (Z_R = R), and has a phase angle of 0°. Common physical representations include the resistance of the electrolyte (R_s or R_e), the charge-transfer resistance (R_ct) associated with the kinetic barrier of a Faradaic reaction, and the resistance of solid-state layers like the solid electrolyte interphase (R_SEI).²

**Capacitor (C)**: A capacitor represents the ability of an interface to store charge. Its impedance is purely imaginary, frequency-dependent (Z_C = 1/(iωC)), and has a phase angle of -90°. The most significant capacitive element in an electrochemical system is the electrical double-layer capacitance (C_dl), which arises from the charge separation at the electrode-electrolyte interface.⁸

**Constant Phase Element (CPE)**: Real-world interfaces are rarely perfectly smooth or homogeneous. Surface roughness, porosity, and non-uniform current distribution can lead to a distribution of relaxation times, causing the interfacial capacitance to behave non-ideally. The CPE is an empirical element used to model this behavior. Its impedance is given by Z_CPE = 1/(Q(iω)^n), where Q is a capacitive parameter and the exponent n (or α) ranges from 0 to 1. When n = 1, the CPE behaves as an ideal capacitor; when n = 0, it is a resistor. For modeling a non-ideal double layer, n is typically between 0.85 and 1.⁸

**Warburg Impedance (W)**: This element specifically models the impedance arising from the diffusion of electroactive species to or from the electrode surface. It is a key indicator of mass-transport limitations in a reaction.⁸

Once a physically plausible equivalent circuit is proposed, a non-linear least squares regression algorithm, such as the Levenberg-Marquardt method, is used to fit the model's theoretical impedance response to the experimental data.⁸ This process involves adjusting the parameter values of each circuit element to minimize the difference between the model and the data. The quality of the fit is typically assessed using the chi-squared (χ²) value, with smaller values indicating a better fit.⁸ A critical step in this process is the selection of appropriate initial "seed" values for the parameters, as poor starting values can cause the fitting algorithm to fail or converge on a physically meaningless local minimum.⁸

### 1.3. The Randles Circuit: A Model for Interfacial Kinetics

The Randles circuit is one of the most fundamental and widely used equivalent circuit models in electrochemistry. It provides a simplified but powerful representation of a simple Faradaic reaction at a planar electrode, where the reaction rate is controlled by a combination of electron-transfer kinetics and semi-infinite linear diffusion.⁷

The circuit, first proposed by John Randles, consists of four key components arranged in a specific topology⁷:

1. **Solution Resistance (R_s)**: This element represents the resistance of the electrolyte between the working and reference electrodes. It is placed in series with all other components, as it is an uncompensated resistance that affects the entire system.¹¹

2. **Double-Layer Capacitance (C_dl)**: This capacitor is placed in parallel with the Faradaic reaction pathway. It models the non-Faradaic process of charging and discharging the electrical double layer at the electrode surface.⁷

3. **Charge-Transfer Resistance (R_ct)**: This resistor is in series with the Warburg impedance and represents the kinetic barrier to the electron transfer reaction. Its magnitude is inversely proportional to the reaction's exchange current density and thus reflects the intrinsic speed of the reaction.⁷

4. **Warburg Impedance (W)**: This element is in series with R_ct and models the impedance associated with the diffusion of electroactive species to and from the electrode surface. Its presence signifies that at certain frequencies, the overall reaction rate is limited by mass transport.⁷

The characteristic Nyquist plot of a Randles circuit can be deconstructed into three distinct frequency regions:

- **High-Frequency Region**: At very high frequencies (ω → ∞), the impedance of the capacitor (1/(iωC_dl)) approaches zero, effectively short-circuiting the parallel branch. The total impedance is therefore dominated by the solution resistance, R_s. This corresponds to the intercept of the plot on the real axis at the high-frequency end.¹³

- **Mid-Frequency Region**: In this region, the interplay between the charge-transfer resistance and the double-layer capacitance becomes dominant. This parallel R_ctC_dl combination produces a distinct semicircle in the Nyquist plot. The diameter of this semicircle is equal to the charge-transfer resistance, R_ct. The frequency at the apex of the semicircle, ω_max, corresponds to the time constant of the interfacial process, given by τ = R_ctC_dl.¹³

- **Low-Frequency Region**: At very low frequencies (ω → 0), the reaction becomes limited by the rate at which reactants can diffuse to the electrode surface. The Warburg impedance, which is proportional to ω^(-1/2), dominates the system's response. This manifests as a straight line with a slope of 1 (an angle of 45°) extending from the low-frequency end of the semicircle.⁷

The Randles circuit serves as a foundational tool for extracting key kinetic and mass transport parameters from EIS data for simple electrochemical systems.

### 1.4. Characterizing Diffusion-Limited Electron Transfer

When an electrochemical reaction is limited by the rate of mass transport of species to the electrode surface, a concentration gradient develops. The impedance associated with this diffusion process is known as the Warburg impedance (W).⁵ While the classic Randles circuit incorporates a simple Warburg element, real-world systems, particularly those with complex geometries like porous battery electrodes, exhibit different types of diffusion that require more sophisticated models.¹⁴

The three primary types of Warburg impedance are:

**Semi-Infinite Warburg Impedance (W)**: This is the classical model used in the Randles circuit. It assumes one-dimensional diffusion to a planar electrode in a solution of infinite depth, where the diffusion layer can grow without restriction.¹⁴ Its impedance is given by Z_W = σω^(-1/2)(1-i), where σ is the Warburg coefficient (in units of Ω·s^(-1/2)).⁵ The Warburg coefficient is inversely related to the concentration and diffusion coefficient of the electroactive species.¹⁴ In a Nyquist plot, this behavior is characterized by a straight line with a constant phase angle of 45°.¹⁴

**Finite-Length Warburg (FLW) or "Short Warburg" (W_s)**: This model is applicable when diffusion occurs through a finite-thickness layer that is terminated by a conductive or "absorbing" boundary condition.¹⁴ At high frequencies, the response is identical to the semi-infinite case (a 45° line). However, at low frequencies, the diffusing species begins to "feel" the boundary. The Nyquist plot deviates from the 45° line and transitions into a second semicircle.¹⁴ This behavior is often observed in systems like porous electrodes where ions diffuse through a finite layer and then react at a current collector interface.¹⁶

**Finite-Space Warburg (FSW) or "Open Warburg" (W_o)**: This model describes diffusion within a finite space that is terminated by a reflective or "blocking" boundary, where no charge transfer can occur.¹⁴ Again, at high frequencies, the response is a 45° line. At low frequencies, as the diffusing species accumulates at the blocking boundary, the impedance becomes purely capacitive. This is represented in the Nyquist plot by the 45° line transitioning to a vertical line (90° phase angle).¹⁴ This model is particularly relevant for describing the diffusion of lithium ions within the solid-state particles of a battery electrode, where the particle boundary is effectively blocking.¹⁴

The choice of an equivalent circuit model is not merely a mathematical fitting exercise but a hypothesis about the physical processes governing the system. While a simple Randles circuit might provide a reasonable fit to the impedance data of a porous battery electrode over a limited frequency range, it is physically inaccurate. The diffusion of lithium ions within the finite, tortuous pore structure of a battery electrode is fundamentally a finite-space process. Applying a semi-infinite Warburg model can lead to significant errors in the interpretation of diffusion coefficients and a misunderstanding of the rate-limiting steps. The correct approach is to select a model based on the known physical structure of the system—such as a transmission line model or a circuit incorporating finite-space Warburg elements—before attempting to fit the data.¹ This underscores a critical challenge in EIS analysis: model ambiguity. Different physical models can produce similar-looking spectra, necessitating the use of complementary techniques or advanced analysis methods like the Distribution of Relaxation Times (DRT) to validate the chosen equivalent circuit and ensure a physically meaningful interpretation.¹⁶

## Section 2: Kinetic and Reversibility Analysis of the [Fe(CN)₆]³⁻/⁴⁻ Redox System

The hexacyanoferrate(III)/hexacyanoferrate(II) redox couple, [Fe(CN)₆]³⁻/⁴⁻, is a cornerstone of electrochemical research. Its well-behaved, rapid electron transfer kinetics make it an ideal model system for validating electrochemical techniques and for probing the properties of electrode surfaces. This section synthesizes literature data to establish benchmark values for its key kinetic parameters and to formally classify its electrochemical reversibility.

### 2.1. Diffusion Coefficient (D) Determination

The diffusion coefficient, D, quantifies the rate of mass transport of an electroactive species due to a concentration gradient. For the [Fe(CN)₆]³⁻/⁴⁻ couple, accurate knowledge of D is essential for interpreting kinetic data. A survey of the literature reveals a general consensus, with most values falling within a relatively narrow range, as summarized in Table 2.1. The slight variations observed are typically attributable to differences in experimental conditions such as temperature, supporting electrolyte identity, and concentration.

**Table 2.1: Literature Values for the Diffusion Coefficient (D) of [Fe(CN)₆]³⁻/⁴⁻**

| Species | D (× 10⁻⁶ cm²/s) | Supporting Electrolyte & Conc. | Method | Reference |
|---------|------------------|--------------------------------|--------|-----------|
| [Fe(CN)₆]³⁻ (ferricyanide) | 7.63 | Not Specified | Not Specified | 17 |
| [Fe(CN)₆]³⁻ (ferricyanide) | 7.3 | Aqueous | Not Specified | 18 |
| [Fe(CN)₆]⁴⁻ (ferrocyanide) | 6.32 | Not Specified | Not Specified | 19 |
| [Fe(CN)₆]³⁻ (ferricyanide) | 7.2 | Aqueous | Not Specified | 18 |
| [Fe(CN)₆]⁴⁻ (ferrocyanide) | 6.8 | Aqueous | Not Specified | 18 |
| [Fe(CN)₆]³⁻/⁴⁻ couple | 3.09 | 50 mM | Steady-state limiting current | 20 |
| [Fe(CN)₆]³⁻/⁴⁻ couple | 9.1 ± 0.9 | Various | EIS | 19 |

Several electrochemical techniques are commonly employed to determine the diffusion coefficient, each relying on a different governing equation:

**Cyclic Voltammetry (CV)**: The diffusion coefficient can be calculated from the peak current, i_p, observed in a cyclic voltammogram using the Randles-Sevcik equation¹⁷:

```
i_p = 0.4463nFAC√(nFνD/RT)
```

Here, n is the number of electrons transferred (1 for this system), A is the electrode area (cm²), ν is the potential scan rate (V/s), and C is the bulk concentration of the species (mol/cm³). A plot of i_p versus ν^(1/2) should be linear for a diffusion-controlled process, and D can be extracted from the slope.

**Electrochemical Impedance Spectroscopy (EIS)**: In an EIS experiment, D is determined from the Warburg coefficient (σ), which is obtained from the low-frequency diffusion tail of the Nyquist plot. For a reversible system with equal bulk concentrations of the oxidized and reduced species (C_O* = C_R* = C*), the relationship is⁵:

```
σ = RT/(n²F²A√2) × (1/C_O*√D_O + 1/C_R*√D_R)
```

where R is the ideal gas constant, T is the absolute temperature, and F is the Faraday constant. The Warburg coefficient can be extracted from the slope of a plot of Z' or -Z'' versus ω^(-1/2) in the diffusion-limited region.

**Chronopotentiometry**: This galvanostatic technique involves applying a constant current, I, and measuring the time required for the concentration of the electroactive species at the electrode surface to drop to zero. This transition time, τ, is related to the diffusion coefficient by the Sand equation:

```
I√τ = nFAC√(πD/2)
```

By measuring τ at a known current, D can be calculated.

### 2.2. Heterogeneous Rate Constant (k₀) Analysis

The standard heterogeneous rate constant, k₀, is the fundamental measure of the kinetic facility of an electron transfer reaction at the equilibrium potential. It represents the intrinsic rate of electron transfer, independent of mass transport effects. The [Fe(CN)₆]³⁻/⁴⁻ couple is renowned for its fast kinetics, making it a "kinetically facile" or "fast" redox couple.

**Table 2.2: Literature Values for the Standard Heterogeneous Rate Constant (k₀) of [Fe(CN)₆]³⁻/⁴⁻**

| Electrode Material | k₀ (cm/s) | Method | Reference |
|-------------------|-----------|--------|-----------|
| Screen-Printed Carbon | 0.0001545 - 0.01590 | EIS | 22 |
| Hydroxylated Carbon Nanotubes | 0.00259 | CV | 18 |
| Carbon Nanotubes | 0.00063 | CV | 18 |

*Note: The values in Table 2.2 show significant variation, which underscores the profound influence of the electrode material and its surface state on the measured kinetic rate. The values from 22 are calculated from R_ct data across a range of concentrations.*

The primary methods for determining k₀ are:

**From Charge-Transfer Resistance (R_ct)**: This is one of the most direct and reliable methods, utilizing EIS. The charge-transfer resistance, obtained from the diameter of the Nyquist semicircle, is inversely proportional to the exchange current, i₀. The exchange current is related to k₀ by the following expression²³:

```
i₀ = nFAk₀C_O*^(1-α)C_R*^α
```

where α is the transfer coefficient (typically assumed to be 0.5). The charge-transfer resistance is then given by:

```
R_ct = RT/(nFi₀)
```

By measuring R_ct from the EIS spectrum, i₀ can be calculated, and subsequently, k₀ can be determined.¹⁹

**From Tafel Analysis**: In cases where a reaction is kinetically slow (irreversible), the relationship between current and overpotential (η) becomes exponential. A Tafel plot of log|i| versus η yields a linear region. Extrapolating this linear portion back to zero overpotential (η = 0) gives the logarithm of the exchange current, log i₀. While the [Fe(CN)₆]³⁻/⁴⁻ couple is generally too fast for this DC method to be easily applied (as the response is dominated by mass transport before a clear Tafel region develops), the underlying principle is fundamental to electrochemical kinetics.

### 2.3. Classification of System Reversibility

The reversibility of an electrochemical system is a kinetic concept that describes the rate of electron transfer relative to the rate of mass transport.²⁵

**Defining Criteria:**

- **Reversible**: The electron transfer kinetics are so fast (k₀ >> √(Dω)) that the concentrations of the oxidized and reduced species at the electrode surface are always in Nernstian equilibrium with the applied potential. For CV, this is characterized by a peak potential separation (ΔE_p) of approximately 59 mV at 298 K, which is independent of the scan rate. The ratio of anodic to cathodic peak currents (i_pa/i_pc) is unity.²¹

- **Quasi-reversible**: The rates of electron transfer and mass transport are comparable. The system deviates from Nernstian equilibrium. In CV, ΔE_p is greater than 59 mV and increases as the scan rate increases. The peak currents are depressed and broadened compared to the reversible case.²⁵

- **Irreversible**: The electron transfer kinetics are very slow (k₀ << √(Dω)). Significant overpotential is required to drive the reaction. In CV, ΔE_p becomes very large, and one of the peaks may be completely absent on the return scan if the reaction is chemically irreversible as well.²⁵

**Justification for [Fe(CN)₆]³⁻/⁴⁻**: The [Fe(CN)₆]³⁻/⁴⁻ redox couple is universally classified as electrochemically reversible under ideal conditions.³⁰ This classification is justified by its characteristically fast electron transfer kinetics, with a high k₀ value on noble metal electrodes. Cyclic voltammograms of the couple on clean platinum, gold, or glassy carbon electrodes typically exhibit a ΔE_p close to 59 mV at moderate scan rates, and the peak currents follow the Randles-Sevcik relationship, confirming diffusion control.²¹ The redox process involves a simple one-electron transfer with minimal structural rearrangement, as no Fe-C bonds are broken or formed, contributing to its kinetic facility.³¹

However, the observed reversibility of this couple is not an immutable property but rather a sensitive indicator of the electrode's interfacial condition. The system's well-behaved nature makes it an excellent probe for surface phenomena. For instance, the modification of an electrode surface with a blocking layer will impede electron transfer, causing the measured response of the [Fe(CN)₆]³⁻/⁴⁻ couple to shift from reversible to quasi-reversible or even irreversible. This is observed as an increase in ΔE_p in CV or an increase in R_ct in EIS.³² Therefore, while the intrinsic kinetics of the couple are fast, the experimentally observed behavior is a convolution of these intrinsic kinetics and the properties of the specific electrode interface being studied. This makes the [Fe(CN)₆]³⁻/⁴⁻ system an indispensable tool for evaluating electrode surface cleanliness, activity, and the efficacy of surface modifications.

## Section 3: Materials and Performance of Graphite|NMC-811 Lithium-Ion Cells

Transitioning from a model aqueous system to a high-performance, non-aqueous energy storage device, this section examines the Graphite|NMC-811 lithium-ion cell. This chemistry is at the forefront of efforts to increase the energy density of commercial batteries, particularly for electric vehicle applications. The analysis begins with a review of the fundamental properties of the constituent materials before synthesizing recent literature to benchmark the performance and degradation of the full cell.

### 3.1. Electrode and Electrolyte Material Properties

The performance of a lithium-ion cell is dictated by the intrinsic properties of its anode, cathode, and electrolyte, as well as the complex interactions at their interfaces.

**Table 3.1: Summary of Key Properties of Graphite, NCM-811, and LP30 Electrolyte**

| Component | Key Structural/Chemical Property | Theoretical Capacity | Typical Voltage Window vs. Li/Li⁺ | Key Advantages | Key Challenges/Degradation Modes |
|-----------|--------------------------------|---------------------|----------------------------------|----------------|----------------------------------|
| Graphite Anode | Layered sp² carbon structure; Li⁺ intercalation to form LiC₆ | 372 mAh/g³⁵ | 0.01 - 0.2 V (fully lithiated)³⁸ | Low cost, high abundance, structural stability, high conductivity³⁸ | SEI formation (irreversible Li loss), Li plating risk, slow diffusion, volume expansion³⁸ |
| NCM-811 Cathode | LiNi₀.₈Co₀.₁Mn₀.₁O₂; Ni-rich layered oxide structure | ~278 mAh/g⁴³ | 3.0 - 4.6 V⁴⁵ | High specific capacity and energy density, reduced cobalt content⁴³ | Particle cracking, TM dissolution, phase transition (H2→H3), cation mixing, poor thermal stability⁴³ |
| LP30-type Electrolyte | 1M LiPF₆ in Ethylene Carbonate (EC) and Dimethyl Carbonate (DMC) | N/A | ~1.0 - 4.5 V (limited stability window) | Good ionic conductivity, enables stable SEI formation (due to EC) | LiPF₆ thermal/hydrolytic instability (HF generation), decomposition at both electrodes, flammability⁴⁹ |

#### 3.1.1. Graphite Anode

Graphite has been the dominant anode material in commercial lithium-ion batteries since their inception due to its compelling combination of low cost, natural abundance, high electrical conductivity, and excellent structural stability during lithium intercalation/deintercalation cycles.³⁸ It possesses a layered crystal structure that can reversibly host lithium ions between its graphene sheets, forming a series of lithium-graphite intercalation compounds (LiₓC₆). The final, fully lithiated stage is LiC₆, which corresponds to a theoretical specific capacity of 372 mAh/g.³⁵

A key electrochemical characteristic of graphite is its very low operating potential, with the main lithiation/delithiation plateaus occurring below 0.2 V versus Li/Li⁺.³⁸ This low potential is highly advantageous as it maximizes the overall cell voltage when paired with a high-potential cathode, thereby increasing energy density. However, this potential is below the reduction potential of common organic carbonate electrolytes. Consequently, during the first charging cycle, the electrolyte decomposes at the graphite surface to form a passivating layer known as the solid electrolyte interphase (SEI).³⁶ A stable SEI is crucial for long-term battery performance as it is electronically insulating but ionically conducting, preventing further electrolyte decomposition while allowing Li⁺ transport. The formation process, however, consumes active lithium, leading to an irreversible capacity loss in the first cycle.³⁶ A major challenge is the risk of lithium plating, where metallic lithium deposits on the anode surface, particularly during fast charging or at low temperatures. This occurs because the graphite potential is very close to that of pure lithium (0 V vs. Li/Li⁺).³⁸ Lithium plating reduces capacity and can lead to the growth of dendrites, which pose a severe safety hazard.

#### 3.1.2. NCM-811 Cathode (LiNi₀.₈Co₀.₁Mn₀.₁O₂)

LiNi₀.₈Co₀.₁Mn₀.₁O₂ (NCM-811) is a member of the nickel-rich family of layered oxide cathodes, developed to meet the demand for higher energy density batteries.⁴³ The high nickel content (80%) is the primary reason for its high specific capacity, as the Ni²⁺/Ni³⁺ and Ni³⁺/Ni⁴⁺ redox couples are the main contributors to charge storage.⁴⁶ This results in a high theoretical capacity of approximately 278 mAh/g and practical capacities often exceeding 200 mAh/g.⁴³

Despite its high capacity, NCM-811 suffers from several significant stability issues, particularly when cycled to high voltages (>4.3 V vs. Li/Li⁺), which is necessary to extract its full capacity.⁴⁶ The primary degradation mechanisms include:

- **Particle Cracking**: Large anisotropic volume changes during deep delithiation induce mechanical stress, leading to inter- and intra-granular cracking of the cathode particles. This exposes new surfaces to the electrolyte, accelerating parasitic reactions and causing loss of electrical contact.⁴³

- **Irreversible Phase Transition**: At high states of charge, the layered structure can undergo an irreversible transition from the active hexagonal phases (H2, H3) to an inactive, rock-salt-like structure on the particle surface, which impedes Li⁺ diffusion.⁴⁶

- **Transition Metal (TM) Dissolution**: The highly oxidative environment at high potentials can lead to the dissolution of transition metal ions (especially Mn and Ni) into the electrolyte.⁴⁶

- **Cation Mixing**: The similar ionic radii of Ni²⁺ (0.69 Å) and Li⁺ (0.76 Å) can lead to Ni²⁺ ions migrating into the lithium layers, a phenomenon known as cation mixing, which blocks Li⁺ diffusion pathways and reduces capacity.⁴³

#### 3.1.3. LP30-Type Electrolyte (1M LiPF₆ in EC:DMC)

The electrolyte serves as the medium for ion transport between the anode and cathode. A standard formulation, often referred to as LP30, consists of 1.0 M lithium hexafluorophosphate (LiPF₆) salt dissolved in a 1:1 (by volume or weight) mixture of ethylene carbonate (EC) and dimethyl carbonate (DMC).⁵¹ This composition represents a balance of properties:

- **LiPF₆ Salt**: Provides a good combination of ionic conductivity, electrochemical stability, and the ability to passivate the aluminum current collector used for the cathode.⁵⁰ However, it is thermally unstable and highly susceptible to hydrolysis in the presence of trace water, which produces highly corrosive hydrofluoric acid (HF) that attacks both the SEI and the cathode material.⁴⁹

- **EC Solvent**: A cyclic carbonate with a high dielectric constant, which is excellent for dissolving the LiPF₆ salt. Crucially, its reduction on the graphite surface forms a stable, passivating SEI layer.⁵⁰ Its main drawback is its high melting point (~36 °C), which necessitates its use in a mixture.

- **DMC Solvent**: A linear carbonate with low viscosity and a low melting point. It is added to EC to improve the electrolyte's ionic conductivity, especially at low temperatures, and to lower the overall freezing point of the mixture.⁵⁰

The electrolyte's electrochemical stability window is limited. It is thermodynamically unstable at the low potential of the graphite anode and at the high potential of the charged NCM-811 cathode, leading to the formation of the SEI and the cathode-electrolyte interphase (CEI), respectively.⁴²

### 3.2. Performance Comparison of Graphite|NMC-811 Cells (2020–2024 Literature)

Recent research has focused intensely on mitigating the degradation of Graphite|NMC-811 cells to unlock their high-energy potential for practical applications. The performance of these cells is highly dependent on the operating conditions (voltage window, C-rate, temperature) and any material-level modifications, such as surface coatings or using single-crystal morphologies.

A critical finding from the literature is the direct link between the instability of the NCM-811 cathode at high potentials and the overall cell degradation. This occurs through a mechanism of electrode "crosstalk." At high voltages, transition metals dissolve from the cathode, migrate through the electrolyte, and deposit on the anode surface.⁵⁵ These deposited metals are catalytically active, promoting continuous decomposition of the electrolyte at the anode. This disrupts and thickens the SEI, consuming cyclable lithium and increasing the cell's impedance.⁵⁵ This vicious cycle, where the cathode's degradation accelerates the anode's degradation, is a primary driver of capacity fade in high-voltage NMC-811 cells. Therefore, strategies that stabilize the cathode surface, such as the LiF coating reported by Llanos et al., are effective because they not only protect the cathode but also prevent the poisoning of the anode, thereby improving full-cell cycle life.⁴⁶

**Table 3.2: Performance Metrics of Graphite|NMC-811 Full Cells from Recent Literature (2020-2024)**

| Cathode Modification | Capacity Retention (% after N cycles) | Avg. Coulombic Efficiency (%) | C-Rate | Voltage Window (V) | Temp. (°C) | Reference |
|---------------------|--------------------------------------|----------------------------|--------|-------------------|-----------|-----------|
| LiF-coated (ALD) | 88% after 500 | >99.5 (inferred) | Not specified | 2.8 - 4.5 | Not specified | 45 |
| Single-Crystal | ~90% after 100 | Not specified | C/3 | 2.5 - 4.4 | 40 | 61 |
| Hydroborate SSE | 97% after 100 | Not specified | C/5 | Not specified | RT | 62 |
| Hydroborate SSE | 75% after 350 | Not specified | C/2 | Not specified | RT | 62 |
| Polycrystalline | Significant fade at 4.3 V | Lower with EC-free electrolyte | Not specified | Up to 4.3 | Not specified | 55 |
| Si-Graphite Anode | Fade accelerated at high temp. | Not specified | Various | Not specified | 5, 23, 45 | 47 |

The data compiled in Table 3.2 highlight several key trends:

- **Cycling Stability and Voltage**: There is a clear trade-off between capacity utilization and cycle life. While charging to 4.5 V or higher unlocks more capacity from the NMC-811, it drastically accelerates degradation mechanisms, leading to poor capacity retention.⁵⁵

- **Impact of Coatings**: Surface modifications, such as the atomic layer deposition (ALD) of a thin, stable LiF layer on the NMC-811 electrode, are highly effective. This artificial CEI physically separates the active material from the electrolyte, suppressing parasitic side reactions, mitigating TM dissolution, and thereby enabling stable cycling at high voltages. The retention of 88% capacity after 500 cycles at a demanding 4.5 V cutoff is a significant improvement over uncoated materials.⁴⁵

- **Morphology (Single-Crystal vs. Polycrystalline)**: Using single-crystal NMC-811 particles instead of polycrystalline agglomerates is another promising strategy. Single crystals are less prone to the inter-granular cracking that plagues polycrystalline materials, maintaining their structural integrity better during cycling. This leads to improved stability, as demonstrated by the ~90% capacity retention after 100 cycles under harsh conditions (4.4 V, 40 °C).⁵⁵

- **Capacity Fade Rates**: The fade rate is highly non-linear and depends on the dominant degradation mechanism. Initial fade is often dominated by irreversible lithium loss to SEI formation. Subsequent, more gradual fade is caused by a combination of ongoing interfacial reactions, active material loss, and impedance growth. In the study by Hasa et al., the capacity loss was approximately 10% over 100 cycles, translating to an average fade rate of ~0.1% per cycle under those specific aggressive conditions.⁶¹

## Section 4: Comparative Evaluation of Electrochemical Characterization Techniques

The comprehensive analysis of electrochemical systems, from model redox couples to complex batteries, relies on a suite of characterization techniques. Each method provides a unique window into the system's behavior, and their true power lies in their synergistic application. This section compares and contrasts the primary techniques discussed in this report: potentiostatic/voltammetric methods, galvanostatic methods, and impedimetric methods.

### 4.1. Potentiostatic/Voltammetric vs. Galvanostatic Methods for Battery Analysis

The two most fundamental modes of operating an electrochemical cell for battery analysis are potentiostatic (controlling voltage) and galvanostatic (controlling current).

**Fundamental Difference**: In a voltammetric technique like Cyclic Voltammetry (CV), the potential of the working electrode is controlled and swept over a defined range, while the resulting current is measured. In a galvanostatic technique like Galvanostatic Charge-Discharge (GCD) cycling, a constant current is applied to the cell, and the resulting potential is measured as a function of time or charge.⁶³

**Cyclic Voltammetry (CV)**:
- **Pros**: CV is an excellent tool for fundamental mechanistic studies. It rapidly identifies the potentials at which redox reactions (e.g., lithium intercalation/deintercalation stages) occur, visible as peaks in the voltammogram. The shape and separation of these peaks provide qualitative and quantitative information about the reaction's kinetics and electrochemical reversibility. It can also be used to estimate diffusion coefficients.²⁹
- **Cons**: CV is not representative of how a battery operates in a real-world application. The linear potential sweep is an artificial driving force, and the measured current is a convolution of both the desired Faradaic reaction current and the non-Faradaic double-layer charging current, which can complicate quantitative analysis.²⁹
- **Application**: Best suited for initial screening of new electrode materials, determining their stable operating potential window, and gaining initial insights into their reaction mechanisms.

**Galvanostatic Charge-Discharge (GCD)**:
- **Pros**: GCD is the industry-standard method for evaluating battery performance because it directly mimics the constant-current charging and discharging conditions of many applications. It provides the most critical performance metrics: practical specific capacity (mAh/g), Coulombic efficiency (charge out / charge in), energy efficiency (energy out / energy in), average voltage, and cycle life.⁶³ Furthermore, the GCD data can be differentiated to generate differential capacity (dQ/dV vs. V) or differential voltage (dV/dQ vs. Q) plots. These plots transform the flat voltage plateaus of the GCD curve into sharp, well-defined peaks, creating a "fingerprint" of the cell's electrochemical state that is highly sensitive to degradation mechanisms like loss of active material or impedance growth.⁶⁶
- **Cons**: A standard GCD curve provides an overall performance metric but offers limited direct insight into the specific kinetic or mass transport limitations. The measured cell voltage includes the thermodynamic potential plus all overpotentials (ohmic, kinetic, and concentration), and the technique does not inherently separate these contributions.⁶⁹
- **Application**: Essential for all performance validation, including determining capacity, cycle life, rate capability, and for state-of-health diagnostics via differential analysis.

### 4.2. Voltammetric vs. Impedimetric Methods for Interfacial Analysis

While both voltammetry and EIS probe the electrode-electrolyte interface, they do so in fundamentally different ways, providing complementary information.

**Fundamental Difference**: Voltammetry applies a large-amplitude DC potential sweep and measures the total current response. EIS applies a small-amplitude AC potential perturbation over a wide range of frequencies and measures the frequency-dependent complex impedance.⁷⁰

**Voltammetry**:
- **Pros**: The technique is relatively fast and provides a direct, qualitative "fingerprint" of the redox processes occurring at the interface. The interpretation of peak potentials and currents is often straightforward for identifying the presence and potential of electrochemical reactions.⁷¹
- **Cons**: The information obtained is convoluted. The measured DC current is the sum of all processes occurring at a given potential, making it difficult to separate contributions from charge transfer, mass transport, and capacitive charging. It has poor time resolution, meaning it cannot easily distinguish between fast and slow processes.⁷⁰

**Electrochemical Impedance Spectroscopy (EIS)**:
- **Pros**: The primary strength of EIS is its ability to deconstruct complex interfacial processes. By probing the system across a wide frequency range, it can separate phenomena based on their characteristic time constants (e.g., double-layer charging is fast and appears at high frequencies, charge transfer is intermediate, and diffusion is slow and appears at low frequencies). Through equivalent circuit modeling, EIS can provide quantitative values for specific physical parameters like solution resistance (R_s), charge-transfer resistance (R_ct), double-layer capacitance (C_dl), and diffusion coefficients. This makes it an exceptionally powerful diagnostic tool for identifying rate-limiting steps and tracking changes in specific interfacial properties during battery aging.⁶⁵
- **Cons**: Data interpretation can be complex and is subject to model ambiguity—a good mathematical fit does not guarantee a physically correct model. The requirement for the system to be at a steady state can be difficult to achieve in a battery, which is an inherently dynamic system. Measurements, especially at low frequencies, can be very time-consuming.³
- **Application**: Ideal for in-depth mechanistic investigations, diagnosing the root causes of performance limitations (e.g., distinguishing between kinetic and transport losses), and quantifying the evolution of interfacial resistances (e.g., SEI and CEI growth) during cycling.

The most effective approach to battery characterization involves a synergistic workflow. Long-term performance and degradation trends are monitored using GCD cycling. The resulting data are processed into dQ/dV plots, which serve as sensitive indicators of changes in the cell's state of health. For example, a shift in a dQ/dV peak to a higher potential on charging indicates an increase in cell polarization. However, the dQ/dV plot alone cannot identify the source of this increased polarization. At this point, EIS can be employed as a diagnostic tool. By performing an EIS measurement on the aged cell and fitting the data, one can determine whether the increased polarization is due to a growth in the SEI resistance (affecting the high-frequency semicircle), a slowdown in the cathode's charge-transfer kinetics (affecting the mid-frequency semicircle), or an increase in mass transport limitations (affecting the low-frequency Warburg region). This integrated approach transforms battery analysis from simple performance testing into a powerful diagnostic cycle: Monitor with GCD/dQdV, then Diagnose with EIS. This workflow is central to modern battery research and development, enabling a deeper understanding of failure mechanisms and guiding the development of more durable materials.

**Table 4.1: Comparison of Electrochemical Characterization Techniques**

| Technique | Controlled Parameter | Measured Response | Primary Application in Battery R&D | Key Information Provided | Pros | Cons |
|-----------|---------------------|------------------|-----------------------------------|-------------------------|------|------|
| Cyclic Voltammetry (CV) | Potential (swept linearly) | Current | Initial material screening, mechanistic study | Redox potentials, kinetic reversibility, qualitative reaction information | Fast, provides a direct "fingerprint" of redox activity | Not representative of battery operation, convolutes Faradaic and capacitive currents |
| Galvanostatic Charge-Discharge (GCD) | Current (constant) | Potential | Performance evaluation, cycle life testing, state-of-health diagnostics | Practical capacity, efficiency, cycle life, voltage profile, degradation trends (via dQ/dV) | Directly simulates real-world use, provides key performance metrics | Convolutes all sources of overpotential, limited direct kinetic information |
| Electrochemical Impedance Spectroscopy (EIS) | Potential (small AC perturbation) | Complex Impedance (AC current) | Mechanistic diagnosis, deconstruction of resistance sources | R_s, R_ct, C_dl, diffusion parameters | Deconvolutes processes by time constant, highly sensitive to interfacial changes, quantitative | Complex interpretation, requires steady-state, can be time-consuming |

## Section 5: Synthesis and Concluding Remarks

This report has conducted a detailed, comparative analysis of two fundamentally different electrochemical systems: the ideal aqueous [Fe(CN)₆]³⁻/⁴⁻ redox couple and the complex, high-energy Graphite|NMC-811 lithium-ion battery. The examination of the hexacyanoferrate system reaffirms its status as a benchmark for fast, electrochemically reversible electron transfer. Its well-defined kinetic parameters and predictable behavior on ideal electrode surfaces make it an invaluable tool not for its own sake, but as a sensitive probe for quantifying the activity, cleanliness, and modification of electrode interfaces. The deviation of its response from ideal reversibility serves as a direct measure of interfacial impedance, a concept that is central to understanding the performance limitations of more complex systems.

In stark contrast, the analysis of the Graphite|NMC-811 cell highlights the immense challenges involved in translating high theoretical material capacities into practical, long-lasting energy storage devices. The central conflict for this chemistry is the inherent trade-off between energy density and cycle life. Accessing the high specific capacity of the NMC-811 cathode requires charging to high cell voltages (e.g., 4.4-4.5 V), but this very condition accelerates a cascade of deleterious degradation mechanisms. These include mechanical failure via particle cracking, irreversible structural changes at the cathode surface, and parasitic reactions with the electrolyte.

A critical conclusion of this review is that the performance and longevity of the Graphite|NMC-811 cell are not merely limited by the individual stabilities of the anode and cathode but are profoundly dictated by the negative "crosstalk" between them. The dissolution of transition metals from the unstable high-voltage cathode and their subsequent migration and deposition on the anode creates a catalytic cycle of degradation at the SEI. This mechanism effectively links the fate of the two electrodes, meaning that stabilizing the full cell requires a holistic approach. Recent advancements, such as the application of protective surface coatings (e.g., LiF) and the engineering of more robust particle morphologies (e.g., single crystals), have shown significant promise in mitigating these issues by stabilizing the cathode interface, thereby breaking this cycle of degradation and enabling improved performance at higher voltages.

Finally, the evaluation of characterization techniques underscores that no single method is sufficient for understanding and solving these complex challenges. A synergistic and methodical approach is required. Galvanostatic cycling remains the ultimate arbiter of practical performance, providing the essential metrics of capacity and cycle life. The analysis of its derivatives, particularly dQ/dV, offers a highly sensitive, non-invasive diagnostic for tracking the evolution of the cell's state of health over its lifetime. When these monitoring techniques signal the onset of degradation, Electrochemical Impedance Spectroscopy serves as the indispensable diagnostic tool, capable of deconstructing the total cell impedance into its constituent parts. This allows researchers to pinpoint the root cause of failure—be it a growing SEI, slowing charge transfer kinetics, or impeded mass transport—and to rationally design targeted material and electrolyte interventions. The continued advancement of high-energy batteries will therefore depend not only on novel material discovery but also on the sophisticated application of this integrated electrochemical toolkit to unravel and ultimately control the complex interfacial phenomena that govern battery performance and failure.

## References

1. Electrochemical Impedance Spectroscopy Part 1: Fundamentals - J-Stage, accessed October 7, 2025, https://www.jstage.jst.go.jp/article/electrochemistry/90/10/90_22-66071/_html/-char/en

2. Modeling and Applications of Electrochemical Impedance Spectroscopy (EIS) for Lithium-ion Batteries, accessed October 7, 2025, https://www.jecst.org/upload/pdf/jecst-2019-00528.pdf

3. Basics of Electrochemical Impedance Spectroscopy - Gamry ..., accessed October 7, 2025, https://www.gamry.com/assets/Application-Notes/basics-of-electrochemical-impedance-spectroscopy.pdf

4. Basics of Electrochemical Impedance Spectroscopy - Gamry Instruments, accessed October 7, 2025, https://www.gamry.com/application-notes/EIS/basics-of-electrochemical-impedance-spectroscopy/

5. Handbook of Electrochemical Impedance Spectroscopy. DIFFUSION IMPEDANCES - BioLogic, accessed October 7, 2025, https://biologic.net/wp-content/uploads/2022/03/2020_zdiffusion.pdf

6. Electrochemical Impedance Spectroscopy and Related Techniques | Advanced Textbooks in Chemistry, accessed October 7, 2025, https://www.worldscientific.com/worldscibooks/10.1142/q0428

7. Randles circuit - Wikipedia, accessed October 7, 2025, https://en.wikipedia.org/wiki/Randles_circuit

8. EIS Data fitting – How to obtain good starting values of ... - Metrohm, accessed October 7, 2025, https://www.metrohm.com/content/dam/metrohm/shared/documents/application-notes/an-e/AN-EIS-007.pdf

9. Equivalent Circuit Modeling in EIS - Gamry Instruments, accessed October 7, 2025, https://www.gamry.com/assets/Application-Notes/Equivalent-Circuit-Modeling-in-EIS.pdf

10. (PDF) Equivalent Circuit Models and Analysis of Electrochemical Impedance Spectra of Caffeine Solutions and Beverages - ResearchGate, accessed October 7, 2025, https://www.researchgate.net/publication/304558679_Equivalent_Circuit_Models_and_Analysis_of_Electrochemical_Impedance_Spectra_of_Caffeine_Solutions_and_Beverages

11. The Randles circuit and electrochemical impedance spectroscopy | Eying, accessed October 7, 2025, https://viadean.notion.site/The-Randles-circuit-and-electrochemical-impedance-spectroscopy-1bb1ae7b9a3280ecbaceeb6ff63cc607

12. Randles Circuit in Electrochemical Systems - VALIPOD, accessed October 7, 2025, https://www.valipod.com/post/randles-circuit-in-electrochemical-systems

13. Electrochemical Impedance Spectroscopy—A Tutorial - PMC, accessed October 7, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC10288619/

14. Diffusion Impedance - Lithium Inventory, accessed October 7, 2025, https://lithiuminventory.com/experimental-electrochemistry/eis/diffusion-impedance/

15. Understanding EIS Modeling - Energsoft, accessed October 7, 2025, https://energsoft.com/blog/f/understanding-eis-modeling

16. Data-driven analysis of electrochemical impedance spectroscopy using the Loewner framework - PMC, accessed October 7, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC11907484/

17. Study of Ferricyanide in Cyclic Voltammetry Using the BAS 100B - BASi, accessed October 7, 2025, https://basinc.com/assets/pdfs/application_capsules/CAP254lett.pdf

18. The determination of the diffusion coefficients of ferricyanide ..., accessed October 7, 2025, https://www.researchgate.net/figure/The-determination-of-the-diffusion-coefficients-of-ferricyanide-ferrocyanide_fig7_303600838

19. (PDF) Kinetics of electron transfer between Fe(CN) 6 3−/4− and ..., accessed October 7, 2025, https://www.researchgate.net/publication/244151851_Kinetics_of_electron_transfer_between_FeCN_6_3-4-_and_poly34-ethylenedioxythiophene_studied_by_electrochemical_impedance_spectroscopy

20. www.rsc.org, accessed October 7, 2025, https://www.rsc.org/suppdata/cc/b8/b819876d/b819876d.pdf

21. Cyclic Voltammetry of Fe(CN)63-/Fe(CN)64- Couple: Evaluation of CV as an Analytical Method - Bloomsburg, accessed October 7, 2025, https://facstaff.bloomu.edu/dmccurry/instrumentation/chem442/manual/experiments/CV.html

22. Investigation of electrochemical behavior of potassium ferricyanide/ferrocyanide redox probes on screen printed carbon electrode through cyclic voltammetry and electrochemical impedance spectroscopy - PMC, accessed October 7, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC10734727/

23. Solved The exchange current (i0) for Fe(CN)6 3− + e = | Chegg.com, accessed October 7, 2025, https://www.chegg.com/homework-help/questions-and-answers/exchange-current-i0-fe-cn-6-3-fe-cn-6-4-20-ma-pt-fe-cn-6-3-20-fe-cn-6-4-20-nacl-10-cell-ca-q88097988

24. EIS study of the redox reaction of Fe(CN)(6)(3-/4-) poly(3,4-ethylenedioxythiophene) electrodes: influence of dc potential and c(Ox): c(Red) ratio - Åbo Akademi University Research Portal, accessed October 7, 2025, https://research.abo.fi/en/publications/eis-study-of-the-redox-reaction-of-fecn63-4-poly34-ethylenedioxyt

25. Electrochemical Reversibility - SOP4CV, accessed October 7, 2025, https://sop4cv.com/chapters/ElectrochemicalReversibility.html

26. 2. Reversibility – Chemical vs. Electrochemical - Chemistry LibreTexts, accessed October 7, 2025, https://chem.libretexts.org/Bookshelves/Analytical_Chemistry/Supplemental_Modules_(Analytical_Chemistry)/Analytical_Sciences_Digital_Library/Courseware/Analytical_Electrochemistry%3A_The_Basic_Concepts/03_Fundamentals_of_Electrochemistry/B%3A_The_Electrode_Process/02_Reversibility__Chemical_vs._Electrochemical

27. Evaluating Electrochemical Reversibility - SOP4CV, accessed October 7, 2025, https://sop4cv.com/chapters/EvaluatingElectrochemicalReversibility.html

28. Cyclic Voltammetry - Data Analysis - BASi, accessed October 7, 2025, https://www.basinc.com/manuals/EC_epsilon/Techniques/CycVolt/cv_analysis

29. Applications of Voltammetry in Lithium Ion Battery Research - Journal of Electrochemical Science and Technology, accessed October 7, 2025, https://www.jecst.org/upload/pdf/jecst-2019-00619.pdf

30. A fundamental study of the thermoelectrochemistry of ferricyanide/ferrocyanide: cation, concentration, ratio, and heterogeneous - RSC Publishing, accessed October 7, 2025, https://pubs.rsc.org/en/content/articlepdf/2020/se/d0se00440e

31. Ferricyanide - Wikipedia, accessed October 7, 2025, https://en.wikipedia.org/wiki/Ferricyanide

32. Achievement of Near-Reversible Behavior for the [Fe(CN)6]3-/4- Redox Couple Using Cyclic Voltammetry at Glassy Carbon, Gold, and Platinum Macrodisk Electrodes in the Absence of Added Supporting Electrolyte - ResearchGate, accessed October 7, 2025, https://www.researchgate.net/publication/231206872_Achievement_of_Near-Reversible_Behavior_for_the_FeCN63-4-_Redox_Couple_Using_Cyclic_Voltammetry_at_Glassy_Carbon_Gold_and_Platinum_Macrodisk_Electrodes_in_the_Absence_of_Added_Supporting_Electrolyte

33. Use of Inner/Outer Sphere Terminology in Electrochemistry—A Hexacyanoferrate II/III Case Study - MDPI, accessed October 7, 2025, https://www.mdpi.com/2673-3293/4/3/22

34. COMPARATIVE STUDIES ON THE REDOX REACTION OF Fe(CN)6 4-/3- AT MODIFIED GLASSY CARBON ELECTRODES VIA DIAZONIUM SALTS ELECTRORED, accessed October 7, 2025, https://revroum.lew.ro/wp-content/uploads/2012/RRCh_9-10_2012/Art%2004.pdf

35. pubmed.ncbi.nlm.nih.gov, accessed October 7, 2025, https://pubmed.ncbi.nlm.nih.gov/35032965/#:~:text=Abstract,%2Dion%20batteries%20(LIBs).

36. A Brief Introduction to Graphite. This article is contributed by AmirReza… | by BatteryBits Editors - Medium, accessed October 7, 2025, https://medium.com/batterybits/a-brief-introduction-to-graphite-7901d4ed7f19

37. Determination of Si/graphite anode composition for new generation ..., accessed October 7, 2025, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10446933/

38. Recent developments in advanced anode materials for lithium-ion batteries, accessed October 7, 2025, https://www.oaepublish.com/articles/energymater.2021.02

39. What determines the 'Potential Window' of Li-ion battery? - Reddit, accessed October 7, 2025, https://www.reddit.com/r/batteries/comments/j76wi0/what_determines_the_potential_window_of_liion/

40. What is Graphite, and Why is it so Important in Batteries? - AquaMetals, accessed October 7, 2025, https://aquametals.com/recyclopedia/why-is-graphite-so-important-to-lib/

41. Modification of graphite anode for lithium ion battery - ResearchGate, accessed October 7, 2025, https://www.researchgate.net/publication/382369683_Modification_of_graphite_anode_for_lithium_ion_battery

42. Voltage vs Li + /Li+ of the electrode materials of five types of batteries in comparison with the stability window of common liquid organic electrolytes. - ResearchGate, accessed October 7, 2025, https://www.researchgate.net/figure/oltage-vs-Li-Li-of-the-electrode-materials-of-five-types-of-batteries-in-comparison_fig1_327867053

43. Exploring the Properties and Potential of Single-crystal NCM 811 for ..., accessed October 7, 2025, https://www.j-cst.org/data/issue/CST/C002201/C00220100036.pdf

44. How to calculate the theorectic capacity of NMC 622 and NMC 811 ..., accessed October 7, 2025, https://www.researchgate.net/post/How-to-calculate-the-theorectic-capacity-of-NMC-622-and-NMC-811

45. High Voltage Cycling Stability of LiF-Coated NMC811 Electrode ..., accessed October 7, 2025, https://pubs.acs.org/doi/10.1021/acsami.3c14394

46. High Voltage Cycling Stability of LiF-Coated NMC811 Electrode - PMC, accessed October 7, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC10797589/

47. Aging mechanisms of NMC811/Si-Graphite Li-ion batteries | Request PDF - ResearchGate, accessed October 7, 2025, https://www.researchgate.net/publication/379466386_Aging_mechanisms_of_NMC811Si-Graphite_Li-ion_batteries

48. Surface-Coated LiNi0.8Co0.1Mn0.1O2 (NCM811) Cathode Materials by Al2O3, ZrO2, and Li2O-2B2O3 Thin-Layers for Improving the Performance of Lithium Ion Batteries - Frontiers, accessed October 7, 2025, https://www.frontiersin.org/journals/materials/articles/10.3389/fmats.2019.00309/full

49. Fire properties of electrolytes for lithium-ion batteries - DiVA portal, accessed October 7, 2025, https://www.diva-portal.org/smash/get/diva2:1276737/FULLTEXT01.pdf

50. Lithium-ion battery - Wikipedia, accessed October 7, 2025, https://en.wikipedia.org/wiki/Lithium-ion_battery

51. Lithium hexafluorophosphate solution in ethylene carbonate and ..., accessed October 7, 2025, https://www.sigmaaldrich.com/US/en/product/aldrich/746711

52. NIPPON STEEL TECHNICAL REPORT No. 84 JULY 2001 - Development of Anode Materials for Lithium Secondary Battery, accessed October 7, 2025, https://www.nipponsteel.com/en/tech/report/nsc/pdf/8404.pdf

53. Graphene/Li-Ion battery - arXiv, accessed October 7, 2025, https://arxiv.org/pdf/1201.2249

54. Recovery of Degraded Ni-Rich NMC811 Particles for Lithium-Ion Batteries - Scholars' Mine, accessed October 7, 2025, https://scholarsmine.mst.edu/cgi/viewcontent.cgi?article=2043&context=che_bioeng_facwork

55. Capacity Fade of Graphite/NMC811: Influence of Particle Morphology, Electrolyte, and Charge Voltage (Journal Article) | OSTI.GOV, accessed October 7, 2025, https://www.osti.gov/pages/biblio/2429478

56. Aging mechanisms of NMC811/Si-Graphite Li-ion batteries - JRC Publications Repository, accessed October 7, 2025, https://publications.jrc.ec.europa.eu/repository/handle/JRC137168

57. Aging mechanisms of NMC811/Si-Graphite Li-ion batteries - DiVA portal, accessed October 7, 2025, https://www.diva-portal.org/smash/get/diva2:1855872/FULLTEXT01.pdf

58. Lithium hexafluorophosphate solution, in ethylene carbonate and dimethyl carbonate, 1.0 M LiPF6 in EC/DMC=50/50 (v/v), battery grade - Scientific Laboratory Supplies, accessed October 7, 2025, https://www.scientificlabs.co.uk/product/alternative-energy/746711-100ML

59. BATTERY ELECTROLYTE LP 30 - Chemical Cloud Database, accessed October 7, 2025, http://www.chemcd.com/product/CCD00287379.html

60. High Voltage Cycling Stability of LiF-Coated NMC811 Electrode - Aalto Research Portal, accessed October 7, 2025, https://research.aalto.fi/files/134728779/CHEM_Llanos_et_al_High_Voltage_2024_ACS_Applied_Materials_and_Interfaces.pdf

61. Quantifying Electrochemical Degradation in Single-Crystalline ..., accessed October 7, 2025, https://link.aps.org/doi/10.1103/PRXEnergy.3.013004

62. (PDF) Hydroborate Solid-State Lithium Battery with High-Voltage ..., accessed October 7, 2025, https://www.researchgate.net/publication/377814623_Hydroborate_Solid-State_Lithium_Battery_with_High-Voltage_NMC811_Cathode

63. What is the difference between Cyclic Voltammetry(CV) and Galvanostatic Charge-Discharge (GCD)?? | ResearchGate, accessed October 7, 2025, https://www.researchgate.net/post/What-is-the-difference-between-Cyclic-VoltammetryCV-and-Galvanostatic-Charge-Discharge-GCD

64. Basics of Potentiostats and Galvanostats: Principles and Essential Modes for Battery Testing, accessed October 7, 2025, https://www.labx.com/resources/basics-of-potentiostats-and-galvanostats-principles-and-essential-modes-for-battery-testi/5471

65. Analysis of the Charging Current in Cyclic Voltammetry and Supercapacitor's Galvanostatic Charging Profile Based on a Constant-Phase Element | ACS Omega - ACS Publications, accessed October 7, 2025, https://pubs.acs.org/doi/10.1021/acsomega.0c04702

66. Differential Analysis of Galvanostatic Cycle Data from Li-Ion Batteries: Interpretative Insights and Graphical Heuristics | Chemistry of Materials - ACS Publications, accessed October 7, 2025, https://pubs.acs.org/doi/10.1021/acs.chemmater.2c01976

67. Comparison between cyclic voltammetry and differential charge ..., accessed October 7, 2025, https://www.researchgate.net/publication/333623225_Comparison_between_cyclic_voltammetry_and_differential_charge_plots_from_galvanostatic_cycling

68. Differential analysis of galvanostatic cycle data from Li-ion batteries: interpretative insights and graphical heuristics | NPL Publications, accessed October 7, 2025, https://eprintspublications.npl.co.uk/9762/

69. Different Efforts but Similar Insights in Battery R&D: Electrochemical Impedance Spectroscopy vs Galvanostatic (Constant Current) Technique - JuSER, accessed October 7, 2025, https://juser.fz-juelich.de/record/911755/files/acs.chemmater.2c02376.pdf

70. Electrochemical Impedance Spectroscopy in the Characterisation and Application of Modified Electrodes for Electrochemical Sensors and Biosensors - PMC, accessed October 7, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC8911593/

71. 11.4: Voltammetric and Amperometric Methods - Chemistry LibreTexts, accessed October 7, 2025, https://chem.libretexts.org/Bookshelves/Analytical_Chemistry/Analytical_Chemistry_2.1_(Harvey)/11%3A_Electrochemical_Methods/11.04%3A_Voltammetric_and_Amperometric_Methods

72. Interfacing techniques for electrochemical sensors: (A) voltometric; (B) impedometric - ResearchGate, accessed October 7, 2025, https://www.researchgate.net/figure/nterfacing-techniques-for-electrochemical-sensors-A-voltometric-B-impedometric-C_fig5_273954901

73. Electrochemical Impedance Spectroscopy in the Characterisation ..., accessed October 7, 2025, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8911593/