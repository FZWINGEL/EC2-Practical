Electrochemical Methods for Interface Characterization
Practical Course for Electrochemistry 2 Module
part of the M.Sc. Battery Materials and Technology

In this practical course, you will learn the basics of electrochemical impedance spectroscopy and how it can be used to better understand electrode-electrolyte interfaces. Specifically, you will study three different electrode-electrolyte systems and correlate various equivalent electrical circuit elements to the processes that are occurring at the interface. With this, you will familiarize yourself with the various forms of presenting and visually analyzing impedance data depending on the dominant processes at the interface. Finally, you will assemble and characterize a lithium-ion coin-cell battery using various electrochemical techniques.

# 1. Theoretical Background

## 1.1 Lithium-ion batteries¹⁻³

Although lithium-ion battery–related research dates back to the late 1960’s, these types of batteries were commercially introduced only 25 years ago.⁴ Since then, however, the field has exploded and with improved cell engineering the so-called “wireless revolution” has occurred. Lithium-ion batteries (LIBs) are a specific type of secondary (rechargeable) batteries that show high energy density, low self-discharge rates and high cycling stability.

### Working principles

Lithium-ion batteries are based on high-surface area, electrically conducting intercalation materials for the electrodes. The liquid electrolyte ensures ionic connectivity between the electrodes, and an electrolyte-permeable separator prevents short-circuiting of the electrodes. During battery cycling, both electrons and ions are transferred at the electrodes. The Li-ions intercalate/de-intercalate from the electrodes, and migrate and diffuse through the electrolyte; whereas the electrons must travel from the electrodes to the current collectors and perform work through an external circuit (Figure 1.1.1). Typically, batteries are assembled in the discharged state, i.e. that during the first charging cycle Li⁺ is de-intercalated from the positive electrode and is intercalated on the negative electrode.

**Figure 1.1.1** Schematic of a Li-ion battery. A graphite negative electrode (anode) intercalates Li-ions during charging. During the discharging process, the Li-ions and electrons move towards the positive electrode (cathode). From Goodenough et al. *J. Am. Chem. Soc.* **2013**, *135*, 1167–1176.

Carbon-based materials are the preferred choice for negative electrodes, and graphite is the most common of them. The structure of graphite allows the reversible expansion of the layers to accommodate 1 Li⁺ for every 6 C atoms (LiC₆). In order to maintain electroneutrality, the π-network of graphite is reduced for every Li-ion intercalated. A common problem with the layered anodes is the co-intercalation of solvated Li and/or solvent molecules. With increasing number of cycles, the co-intercalated ions and molecules exfoliate the graphite

---

layers and slowly decompose the electrode. Despite these issues, many cycles can be obtained from a graphite anode-based battery (ca. 30,000). Alternate, insertion-based negative electrodes have been developed.⁵ These anodes are usually made out of Si or Sn metal to form Li alloys during charging. Their main drawback is the detrimental effect of large volume changes upon Li-insertion (over 250% volume expansion in Si to form Li₁₅Si₄).

For the cathode electrode, titanium disulfide ($TiS_2$) was one of the first proposed materials and Whittingham successfully demonstrated its viability in 1976.⁶ However, the sulfide cathode leads to low attainable voltages and it encouraged Goodenough and coworkers to explore $LiMO_2$ (M = Co, Ni, Cr) as possible cathode materials with higher voltages due to the oxide environment.⁷ $LiCoO_2$ stood out as the best performing material (in terms of long-term stability and potential window) and with an open-circuit potential of up ca. 4 V vs. $Li/Li^+$. Since then, the study of metal oxide-based materials has dominated the research efforts in the battery cathode community.

Aqueous electrolytes are not used in LIBs because lithium is a highly electropositive material and reacts rapidly with water to evolve hydrogen gas. Therefore, organic solvents and/or ceramic electrolytes are used in LIBs. Typical solvents are polar, aprotic carbonate-based (e.g. ethylene carbonate, propylene carbonate, dimethyl carbonate) in which an organic or inorganic lithium salt is dissolved (e.g. $LiClO_4$, $LiBF_4$, $LiPF_6$, $LiSO_3F$). The resulting conductivities of these electrolytes are in the range of $10^{-3}-10^{-2}$ S/cm. Although these electrolytes are thermodynamically unstable, during the first charging cycle the electrolyte decomposes irreversibly as a new solid phase on the negative electrode. This phase, known as the solid-electrolyte interphase (SEI), is permeable to $Li^+$ but electronically insulating. The SEI prevents further electrolyte decomposition and allows long-term, safe battery cycling. The separator material used to prevent contact between the electrodes is usually a microporous film of polyolefin or fiberglass materials.

## Characterization of batteries

Over the years, certain electrochemical techniques and figures of merit have been established for benchmarking the performance of a battery. The capacity of a battery (Q, Ah/kg or Ah/L) is among the most important figures of merit. The capacity is the charge per mass of active material. The theoretical capacity ($Q_{theo}$) of a material can be calculated using Faraday's laws, by taking into account the molecular weight of the host material and the charge of the electrochemically active component. In practice, Q deviates from $Q_{theo}$ and in the battery research community it is known as *capacity fade*. The causes for capacity fade can be reversible or irreversible. Reversible capacity losses are due to slow diffusion of ions at high current densities. Electrolyte side-reactions, electrode decomposition, and/or SEI layer formation are examples of irreversible capacity losses.

For a specific charge-discharge battery cycle, the Coulombic efficiency (CE) can be determined using Equation 1.1.1. This ratio tells the percent of usable charge in a battery per cycle – when it falls below 80% the battery has reached its total practical lifetime. The values for the discharge capacity ($Q_{discharge}$) and the charge capacity ($Q_{charge}$) are extracted from charge-discharge curves.

$$
CE = \frac{Q_{discharge}}{Q_{charge}} \cdot 100 \quad [1.1.1]
$$

---

LIBs can be charged under controlled-potential or controlled-current conditions. At constant potential conditions, a charging potential ($E_{ch}$) is applied to the cell for a specific time ($t$) until the current reaches a constant value, referred to as the residual current ($i_{res}$). The constant-current charging method relies on providing charge to the battery at a specific current ($i$, A/kg) until a cut-off potential is reached, or the charging period has ended. Figure 1.1.2 shows exemplary how an experiment under constant current conditions is conducted. Constant potential conditions will give a similar response, only voltage and current are exchanged.

**Figure 1.1.2** (left) A constant current is applied over time $t$, which induces either the oxidation or reduction reaction in the battery. (right) As a consequence, the potential of the battery will change, the time during which the reaction is happening is called transition time $\tau$.

The time a charging/discharging step needs is called transition time $\tau$. Depending on the number of electrochemically active species, the charge/discharge step can have several plateaus with distinct transition times (Figure 1.1.3). An experiment where the sign of the applied current is flipped in a regular pattern is called cyclic chronopotentiometry, which is typically used for batteries. Since different batteries have different $Q_{theo}$, the C-rate parameter was introduced (Equation 1.1.2). Here, the charging/discharging rate is normalized to $Q_{theo}$ and allows for a more direct comparison between various systems. A C-rate of 1 means that the battery (no matter what its components are) is (dis)charged completely in exactly 1 hour. Likewise, a 4C-rate indicates that the battery is (dis)charged in 15 minutes.

$$
1C = \frac{i_{applied}}{m_{electrode}} \cdot Q_{theo} \quad [1.1.2]
$$

**Figure 1.1.3** (top) Galvanostatic (controlled-current) charge-discharge curve of graphite during at C/7 vs. Li. The various stages and phase transitions during the (de)intercalation of Li are also identified and correlated to the state-of-charge of the battery (bottom). The phase transition regions are derived from the phase diagram at ca. 300 K by Woo et al.⁸,⁹ Figure taken from *Kinetics and Stage Transitions of Graphite for Lithium-Ion Batteries*, a dissertation by Michael Heß (2013).¹⁰ For more information about the different phases and phase transitions in graphite during (de)intercalation refer to Ref. 10.

---

An advantage of the controlled-current method is the possibility of calculating differential capacity plots (also known as incremental capacity plots). These plots are the derivative of the charge-discharge curves with respect to the potential (dQ/dE, Equation 1.1.3). These plots can provide information about the multiple redox processes, phase transitions and cell degradation processes during battery cycling.<sup>11-13</sup>

$$
\frac{dQ}{dE} \approx \frac{Q_{k+1} - Q_k}{E_{k+1} - E_k} \qquad [1.1.3]
$$

However, there are some drawbacks during constant current cycling. Firstly, the applied current leads to changing potentials, simultaneously a non-faradaic current is present which charges the double-layer capacitance.

$$
i_f = i - i_c \qquad [1.1.4]
$$

Therefore, the faradaic current $i_f$ is given by the total current $i$ minus the charging current $i_c$, when the area of the electrode is not changing over time. As the potential changes during charging/discharging, the contributions of currents to faradaic and non-faradaic processes will adapt.

**Figure 1.1.4** Display of the faradaic current fraction against the transition time for different K values ($K = (\frac{RT}{nF}) \frac{C_d}{nF C_O^*(\pi D_O \tau)^{\frac{1}{2}}}$ represents charging contributions).

Figure 1.1.4 shows clearly the longer the transition time, the higher the contribution of non-faradaic processes. Furthermore, the formation of the double layer during the electrochemical processes also makes the determination of the transition time $\tau$ inaccurate as the faradaic current varies over time. Additionally, problems for longer transition times can occur due to the onset of convection or non-linear diffusion. The convective effects, movement of solution with respect to the electrode, can be counteracted by shielding the electrodes with a glass mantle. For constant current electrolysis reactions, the Sand equation (Equation 1.1.5) gives the relationship between current $i$, transition time $\tau$, concentration $C$ and diffusion coefficient $D$.

$$
\frac{i\tau^{1/2}}{C_O^*} = \frac{nFAD_O^{1/2}\pi^{1/2}}{2} = 85.5 nD_O^{1/2} A \left( \frac{mA - s^2}{mM} \right) \qquad [1.1.5]
$$

The transition time constant $\frac{i\tau^{1/2}}{C_O^*}$ should be independent of $i$ or $C$, if this is not the case, the electrode reactions suffer from complications due to adsorption, double layer charging or convection.

The electrode voltages for fast reaction are typically described by the Nernst equation, which connects the potential with the concentration. This can be substituted to find an expression which links the transition time with the electrode potential (Equation 1.1.6).

---

$$
E = E_{\tau/4} + \frac{RT}{nF} \ln \left( \frac{\tau \frac{1}{2} - t \frac{1}{2}}{t \frac{1}{2}} \right) \qquad [1.1.6]
$$

$$
E_{\tau/4} = E^{0'} - \frac{RT}{2nF} \ln \frac{D_o}{D_R} \qquad [1.1.7]
$$

Within the equation is the quarter-wave potential (Equation 1.1.7), which is the constant current equivalent of the half-wave potential.

## 1.2 Electrochemical Impedance Spectroscopy¹⁴,¹⁵

Potential-step, current-step and potential-sweep electrochemical techniques rely on large perturbations of the electrode that drive the reaction under study away for equilibrium and monitoring its response (usually a transient signal). A different approach to study an electrochemical reaction relies on small alternating signal perturbations and measuring the cell response in the steady-state regime. Electrochemical impedance spectroscopy (EIS) relies on the latter approach and measures the impedance of an electrochemical cell under small AC perturbations over a large range of frequencies (typically 10⁻⁴ to 10⁶ Hz). The main advantage to EIS the possibility to execute high-precision, long-term steady-state measurements of the perturbed system.

The ratio of the potential (E) and the current (j) that defines the resistance (R) of a circuit (Ohm's Law) is limited to DC systems. The resistance of an electrical circuit under an AC perturbation (either potential or current) is called the impedance (Z, Equation 1.2.1). Unlike a typical resistor, Z depends on the frequency (ω, rad/s) of the AC perturbation.

$$
Z(\omega) = \frac{E(t)}{j(t)} = \frac{E_{max} \cdot e^{i\omega t}}{j_{max} \cdot e^{i(\omega t - \phi)}} \qquad [1.2.1]
$$

$$
\omega = 2 \cdot \pi \cdot f
$$

Here, *j(t)* is the current at a specific time (*t*), *j*<sub>max</sub> is the amplitude of the current signal and likewise, *E*<sub>max</sub> is the amplitude of the potential. Here, the variable *i* is the imaginary number (-1)<sup>1/2</sup> and *f* is the frequency in s<sup>-1</sup>. The variable *φ* is the phase shift (also known as the phase angle) between the potential and the current signals. During a typical EIS experiment, amplitudes *E*<sub>max</sub> and *j*<sub>max</sub>, as well as *φ*, are measured at different applied frequencies.

### A quick review on AC circuits

The interpretation of EIS results is based on establishing analogies between electrical circuit elements and electrochemical processes. Therefore, it is important to have some idea about how different circuit elements behave under AC conditions.

For example, consider a pure resistor **R** to which a sinusoidal voltage is applied

$$
E_{appl} = E_{max} \cdot \sin(\omega \cdot t) \qquad [1.2.3]
$$

Since Ohm's law always holds

$$
j_{\text{meas}} = \frac{E_{\text{appl}}}{\mathbf{R}} \Leftrightarrow \frac{E_{\text{max}}\sin(\omega t)}{\mathbf{R}} \Leftrightarrow j_{\text{meas}} = \left(\frac{E_{\text{max}}}{R}\right)\sin(\omega t) \qquad [1.2.4]
$$

---

and the phase angle between the signals is zero (Figure 1.2.1A). Now, consider applying the same sinusoidal voltage across a pure capacitor **C**. Recalling the typical definition of a capacitor

$$
j_{\text{meas}} = C \left(\frac{dE_{\text{appl}}}{dt}\right), \text{ therefore}
$$

$$
j_{\text{meas}} = C \left(\frac{d}{dt} E_{\text{max}} \sin(\omega t)\right) \Leftrightarrow j_{\text{meas}} = C[\omega E_{\text{max}} \cos(\omega t)] \quad [1.2.5]
$$

which is then re-expressed as a sine function as

$$
j_{\text{meas}} = \omega C \left[ E_{\text{max}} \sin\left(\omega t + \frac{\pi}{2}\right) \right] \quad [1.2.6]
$$

and a phase shift of 90° between the signals is measured (Figure 1.2.1B). The impedance is, therefore, a generalized version of Ohm's law. The phase angle measures the balance between resistive and capacitive components in a system. For pure resistances, $\phi = 0$; for pure capacitive behavior, $\phi = 90$; and for a mixture of both, intermediate phase angles are observed.

**Figure 1.2.1.** Relationship between the alternating potential and current across (A) a pure resistor and (B) a pure capacitor.

An electrochemical interface requires the combination of more than one circuit element to describe it, i.e. **equivalent electrical circuits** are needed to study such interfaces. The fundamental physics of electrical circuits also applies to AC elements. Therefore, when two elements are connected in series, the total impedance is simply a sum of both impedances, whereas for elements in parallel it is the inverse of the sum of the reciprocal impedances. Note that the inverse of the impedance is also known as the admittance $Y$. Table 1 summarizes the impedance response of various elements that are often used to build equivalent electrical circuits in the analysis of impedance spectra.

**Table 1.** Commonly used electric elements (and their impedances) for electrochemical impedance spectroscopy analysis

<table><thead><tr><td>Element</td><td>Symbol</td><td>Impedance</td></tr></thead><tbody><tr><td>Resistor</td><td>R</td><td>Z<sub>R</sub> = R</td></tr><tr><td>Capacitor</td><td>C</td><td>Z<sub>C</sub> = <sup>1</sup>&frasl;<sub>i&omega;C</sub></td></tr><tr><td>Inductor</td><td>L</td><td>Z<sub>L</sub> = i&omega;L</td></tr></tbody></table>

---

<table><tr><td>Warburg<br>Used to describe semi-infinite linear diffusion of electroactive species to an electrode</td><td>W</td><td>Z<sub>Warburg</sub> = <sup>A<sub>W</sub></sup>&frasl;<sub>&omega;<sup>1/2</sub></sub> - i <sup>A<sub>W</sub></sup>&frasl;<sub>&omega;<sup>1/2</sub></sub><br>A<sub>W</sub> = <sup>RT</sup>&frasl;<sub>n<sup>2</sup>F<sup>2</sup>A&radic;2</sub> &nbsp; ( <sup>1</sup>&frasl;<sub>c<sub>ox</sub>D<sub>ox</sub><sup>2</sup></sub> + <sup>1</sup>&frasl;<sub>c<sub>red</sub>D<sub>red</sub><sup>2</sup></sub> )<br>R: gas constant, T: temperature, n: number of transferred electrons, A: electrode area, ci: concentration of electroactive species i at surface of the electrode, Di: diffusion coefficient of electroactive species i.</td></tr><tr><td>Constant-phase element (CPE)<br>Often used to describe non-ideal capacitive behavior</td><td>P</td><td>Z<sub>P</sub> = <sup>1</sup>&frasl;<sub>Q(i&omega;)<sup>a</sup></sub><br>Q: admittance, a: ideality factor (ranges from 0 to 1)</td></tr><tr><td colspan="2">Reminder:<br>Z<sub>Series</sub> = Z<sub>1</sub> + Z<sub>2</sub> + Z<sub>3</sub>... and<br>Z<sub>Parallel</sub> = [<sup>1</sup>&frasl;<sub>Z<sub>1</sub></sub> + <sup>1</sup>&frasl;<sub>Z<sub>2</sub></sub> + <sup>1</sup>&frasl;<sub>Z<sub>3</sub></sub> + ...]<sup>-1</sup><br>(1, 2, 3 ... stand for different circuit elements)</td><td></td></tr></table>

**How do electrochemical processes relate to electrical circuit elements and how to choose (or build) the appropriate equivalent circuit?**

With this understanding of AC circuit elements and their impedances, we can start establishing relationships between electrochemical processes and circuit elements that will aid in the analysis of impedance spectra. For instance, the electrochemical double-layer can be described by a capacitor element. This is based on the Gouy-Chapman-Stern double-layer model at electrode-electrolyte interfaces. Another process that affects the impedance response of an electrochemical cell is the conductivity of the electrolyte. By definition, the conductivity is the inverse of the resistivity and a resistor is used to describe bulk ionic conduction in the electrolyte. The rate of electron-transfer to and from and electroactive species and the electrode is also described with a resistor (often called charge-transfer resistance $R_{CT}$). To describe the linear diffusion of redox species close to the electrode surface, a Warburg element ($Z_{Warburg}$) is used.

Now, to build the equivalent electrical circuit, it is necessary to look more closely at the overall current and potential at the electrode of interest. In the simplest case of a planar working electrode in contact with a liquid electrolyte, the overall potential is the sum of the interfacial potential and the bulk electrolyte Ohmic drop (Figure 1.2.3).

---

$Z_{\text{interface}}$

**Figure 1.2.3 Equivalent electrical circuit corresponding to a single reaction on a uniformly accessible (planar) electrode.** Here, a series combination of the electrolyte resistance and the interfacial impedance are used to describe the planar electrode. A parallel combination of the faradaic impedance and the double-layer capacitance comprise the interfacial impedance. Adapted from Figure 9.1 in *Electrochemical Impedance spectroscopy* by Orazem and Tribollet, 2008, pp. 156.

At the interface itself, the overall current includes the faradaic current and the charging current across the double-layer capacitor. The faradaic current is the current related to charge transfer processes, e.g. electron-transfer. Thus, the total impedance of the system includes the Ohmic drop of the bulk electrolyte, as well as the faradaic and capacitive interfacial impedances. Note that the $Z_{\text{faradaic}}$ and $Z_{\text{capacitive}}$ are in parallel because both of these processes occur at the same interface even if one process dominates the measured impedance at certain frequencies. In the absence of faradaic (charge-transfer) processes, e.g. a blocking electrode immersed in an inert electrolyte, the circuit reduces to a simple RC-series circuit. The elements used to describe the faradaic impedance will depend on the electrochemical cell.

In this laboratory experience, you will study three different cells and the equivalent electrical circuits that describe them are shown in Figure 1.2.4. Specifically, circuit (a) describes the interfacial impedance of an ideally polarizable electrode in contact with an electrolyte, where the resistance and the capacitor resemble the bulk electrolyte conductivity and the electrochemical double-layer capacitance, respectively. And circuit (c) describes the impedance of a system in which the charge-transfer resistance ($R_2$) is influenced by the mass-transport ($Z_{\text{MT}}$) of electroactive species to and from the electrode surface.

**Figure 1.2.4 Equivalent electrical circuits of relevance for the present experiment**

## How does the impedance data look like and how to analyze it?

One of the most common ways to present EIS data is in a *Nyquist plot* (–$Z_{\text{Im}}$ vs. $Z_{\text{Re}}$ for different values of $\omega$). Figure 1.2.5 shows the Nyquist plots for two simple RC circuits, one in which the elements are in series and one in which the elements are in parallel.

---

**Figure 1.2.5** Nyquist plots for simple RC circuits connected in series (left) and in parallel (right); R₁ = 50 Ω and C₁ = 1 μF. The arrow indicates the direction in which the frequency is increasing.

The spectra shown in Figure 1.2.5 are not only significantly different in shape, but also in the magnitude of the impedance in Z<sub>Im</sub>, even though the values of R<sub>1</sub> and C<sub>1</sub> are exactly the same. Note, however, that the resistance can be obtained directly from the x-axis on both plots and it is exactly 50 Ω. The difference between the spectra can be easily identified when we calculate the total impedance of each circuit independently:

For the RC-series circuit, based on the basics of electrical circuits, the impedance of each element is simply added:

$$
Z_{\text{Series}} = Z_{R1} + Z_{C1} \qquad [1.2.7]
$$

$$
Z_{\text{Series}} = R_1 + \frac{1}{i\omega C_1} \qquad [1.2.8]
$$

which can be re-expressed as

$$
Z_{\text{Series}} = R_1 - i \frac{1}{\omega C_1} \qquad [1.2.9]
$$

Now, it is evident that the frequency-dependence lies only on the capacitive element and that at low frequencies the imaginary impedance Z<sub>Im</sub> approaches infinity. (Hence the observed straight line with an intercept at 50 Ω at Z<sub>Re</sub> in the Nyquist plot).

In the case of the RC-parallel circuit,

$$
Z_{\text{Parallel}} = \left[ \frac{1}{Z_{R_1}} + \frac{1}{Z_{C_1}} \right]^{-1} \qquad [1.2.10]
$$

$$
Z_{\text{Parallel}} = \left[ \frac{1}{R_1} + (i\omega C_1) \right]^{-1} = \left[ \frac{1 + R_1(i\omega C_1)}{R_1} \right]^{-1} = \frac{R_1}{1 + (i\omega C_1 R_1)} \quad [1.2.11]
$$

where the product R₁C₁ is often referred to as the time constant τ. The impedance of the parallel circuit results in a semicircle in the Nyquist plot with a diameter equal to R₁. The value of the time constant can be obtained directly from the inverse of the frequency at the maximum point in the semicircle in Figure 1.2.3, i.e. at ω ~ 3.26 kHz.

$$
f = \frac{\omega}{2\pi} = 3.26 \text{ kHz} \quad \therefore \quad \omega = 20,000 \text{ Hz} \Rightarrow \omega^{-1} = 50 \text{ } \mu\text{s} \qquad [1.2.12\text{a}]
$$

---

$$
\tau = R_1 C_1 = 50 \ \mu s \qquad [1.2.12b]
$$

The impedance spectra presented so far are the simplest of equivalent electrical circuits. The complexity of measured impedance spectra, as well as its analysis, often requires the use of sophisticated software for correct interpretation based on regression analysis of the fits of the equivalent electrical circuit used and the measured data. Moreover, it is possible that different equivalent electrical circuits show the same impedance response. The circuits shown in Fig 1.2.6, although developed from different physical models, exhibit mathematically equivalent frequency responses because they share time constants of the same magnitude.


**Figure 1.2.6** Two mathematically equivalent electrical circuits and their impedance responses. $R_1 = 25 \ \Omega$, $R_2 = 50 \ \Omega$, $R_3 = 100 \ \Omega$, $C_1 = 1 \ \mu F$, $C_2 = 1 \ mF$

Therefore, a good fit alone does not validate the model used to describe the data. Impedance spectroscopy results must be coupled with complementary techniques that can help verify the proposed model and the values extracted from it. Moreover, a prior knowledge of the system under study is critical to help develop the models and circuits to better describe it. In this laboratory experience you will learn how to analyze and interpret impedance results visually and without the aid of analysis software, rather you will use the theory discussed so far (and later on) to qualitatively obtain information about three different electrochemical interfaces.

## The constant-phase element (CPE, P)

This element was designed to explain deviations from ideal circuit elements in the measured data of real systems.

$$
Z_{\mathrm{P}}=\frac{1}{Q(i \omega)^{\alpha}} \qquad [1.2.13]
$$

Here, Q is the admittance (inverse of the impedance) and $\alpha$ is the ideality factor (its value ranges from 0 to 1). Note that when $\alpha = 1$, Equation 1.2.13 reduces to the definition of an ideal capacitor element. The CPE does not have a frequency-dependent phase angle and hence the name “constant-phase” element. Moreover, the phase angle of the CPE is defined as $-(\alpha \times 90^{\circ})$. The CPE has units of $\Omega^{-1}s^{\alpha}$. The effects of the CPE in parallel with a resistor are depicted in Figure 1.2.7. The CPE causes a depression of the semicircles that depends on the ideality factor.

---

**Figure 1.2.7** Impedance Nyquist plot of a CPE in parallel with a resistor. When $\alpha = 1$, the CPE reduces to a capacitor and the spectrum obtained is that of an RC-parallel circuit. The angle of depression of the semicircles is $(1-\alpha) \times 90^{\circ}$. The value of the ideality factor of each spectrum is identified in the plot. Here, $R_1 = 1 \Omega$ and $Q_1 = 1 \times 10^{-3} \Omega^{-1} s^{\alpha}$.

### A special constant-phase element: The Warburg element

The Warburg element is type of constant-phase element with a characteristic angle of 45°, i.e. $\alpha = 0.5$. This element describes the semi-infinite linear diffusion of electroactive species to the electrode surface (Equation 1.2.14). The semi-infinite linear diffusion condition implies that the diffusion layer grows at a rate that is dependent on the square root of the time ($t^{1/2}$) and approaches the bulk concentration of the electroactive species.

$$
Z_{\text{Warburg}} = \frac{A_W}{\omega^{1/2}} - i \frac{A_W}{\omega^{1/2}} \quad [1.2.14]
$$

Here, $A_W$ is the Warburg coefficient and is defined as

$$
A_w = \frac{RT}{n^2 F^2 A \sqrt{2}} \left( \frac{1}{c_{\text{ox}} D_{\text{ox}}^2} + \frac{1}{c_{\text{red}} D_{\text{red}}^2} \right) \quad [1.2.15]
$$

where $n$ is the number of transferred electrons, $R$ is the ideal gas constant, $T$ is the temperature in K, $F$ is Faraday's constant, $A$ is the electrode area, $c$ is the concentration at the surface of the electrode in $mol/cm^3$, $D$ is the diffusion coefficient in $cm^2/s$ and the subscripts ox and red indicate the oxidized and reduced species, respectively.

The presence of Warburg behavior can be easily identified in the Nyquist plot as a straight line at a 45° angle. It is confirmed by plotting the magnitude of the impedance $|Z|$ versus the frequency and looking for a linear regime with a slope of $-1/2$ (in a logarithmic scale for both axes). Following Equation 1.2.14, to determine the Warburg coefficient, the real and imaginary parts of the impedance are plotted versus the square root of the frequency ($\omega^{1/2}$). The slope of both parts of the impedance should be equal and from it one directly obtains the value of $A_W$.

### Which electrode interface is measured during an EIS experiment?

In a three-electrode cell, i.e. where a reference electrode is used, along with a counter and working electrode, the impedance data is related to processes occurring at the **working electrode** (WE). The counter electrode (CE) is needed to maintain the electroneutrality of the system. The reference electrode (RE) is used to isolate the impedance response of the WE from the CE. Moreover, the RE provides a potential scale to study the electrochemical reactions at the working electrode against (or versus) a stable redox reaction. The use of a reference electrode in electrochemical measurements is analogous to the use of the Kelvin and Celsius scales when measuring changes in temperature. A common reference electrode is the silver-silver chloride (Ag/AgCl) electrode. The standard reduction potential of this electrode is 0.197 V vs. SHE (when immersed in a saturated KCl solution and prepared with

---

a 3 M KCl liquid electrolyte.¹⁴

## 2. Experimental Section

### 2.1 Lithium-ion battery: Assembly and Characterization

In this experiment, you will characterize a coin cell with the configuration Cu|Graphite|LP30|NCM-811|Al. The liquid electrolyte LP30 is a 1 M LiPF₆ salt dissolved in 1:1 volume ratio of ethylene carbonate and dimethyl carbonate solution. The anode and cathode are graphite ($ρ_{area} = 6.3 \text{ mg·cm}^{-2}$) and NCM-811 ($ρ_{area} = 10.5 \text{ mg·cm}^{-2}$) sheets, respectively. The separator is a microporous monolayer membrane from polypropylene (Celgard 2500). For assembly, a disc of a defined diameter will be punched out for the anode (ø14 mm), cathode (ø12 mm) and separator (ø16 mm) material. (*The diameter of the coin cell is 20 mm*) Afterwards, all components are stacked, in between each layer some LP30 electrolyte (total 100 µl) is added. To finalize the coin cell, the coin cell is crimped.

#### Battery efficiency and capacity

a. Galvanostatic charging-discharging cycles at 0.1 C, 1C and 2C with maximum and minimum voltages of 2.4 and 4.2 V vs. Li, respectively.

### 2.2 Electrochemical characterization of two different interfaces

For your electrochemical experiments, you will be using a glassy carbon (GC) working electrode with a diameter of 5.5 mm. The counter electrode will be a Pt wire and the reference electrode is a Ag/AgCl electrode in saturated KCl solution. It is important that you rinse all the electrodes thoroughly with Nanopure® water prior to immersing them in the electrolyte solution. It is also very important that you return the reference electrode to the storage solution (sat'd KCl) when not in use during the experiments, as this ensures that the reference potential is maintained.

#### 2.2.1 Electrode|Inert electrolyte + Electroactive species in solution – The Randles Circuit (Interface No. 1)

The composition of this electrolyte is simply a 2 M potassium chloride with 2.3 mM of $K_3Fe(CN)_6$ in nanopure water.

1. Measure cyclic voltammograms in the potential range of +0.9 V to -0.1 V vs. Ag/AgCl at different scan rates (500, 200, 100, 50, 25, 10 mV/s) → [3 cycles at each scan rate]

2. Measure a single impedance spectrum at the formal potential (equilibrate for 30 s) of the $[Fe(CN)_6]^{3-/4-}$ redox couple in the frequency range of 0.1 MHz to 500 mHz and with a potential amplitude of 10 mV.

3. Measure how the potential changes as a function of time at the following applied currents: -1 µA, -10 µA and -20 µA.

#### 2.2.2 Rotating Electrode|Inert electrolyte + Electroactive species in solution (Interface No. 2)

The composition of this electrolyte is also 2 M potassium chloride with 2.3 mM of $K_3Fe(CN)_6$ in nanopure water.

1. Measure linear sweep voltammograms in the potential range of +0.9 V to -0.1 V vs. Ag/AgCl at a scan rate of 50 mV/s as a function of electrode rotation rate (250, 500, 1000, 1500, 2000, 2500 rpm)

2. Measure a single impedance spectra at the formal potential (equilibrate for 30 s) of the $[Fe(CN)_6]^{3-/4-}$ redox couple in the frequency range of 0.1 MHz to 50 mHz and with a

---

potential amplitude of 10 mV at each rotation rate.

# 3. Tasks and Analysis

## 3.1 Lithium-ion battery characterization

### 1. Cyclic voltammogram

* Plot all three cycles overlaid in the same graph

* Describe the obtained curves and important points like oxidation, reduction and formal potential. Do you observe any differences? What are the possible causes for these changes in the voltammograms? Do you observe any peaks? If so, what electrochemical process(es) could these peaks be related to?

### 2. Galvanostatic charge-discharge cycles

* Plot all the cycles at different C-rates in a single graph of Potential vs. Capacity (mAh/g): (Use the material's specifications to correctly normalize your data)

* Describe the graph. Indicate the charging and discharging regimes. Comment on the chosen voltage window, does it allow to use NMCs full capacity? How does the experimental capacity compare to the theoretical capacity of the graphite material? Does the capacity depend on the C-rate?

* Calculate the Coulombic efficiency (CE) of your battery at all C-rates

* Does the CE depend on the C-rate?

* Generate the derivative plots ($dQ/dE$ vs. $E$) of all cycles and discuss your results. Try to match the charging and discharging processes to an electrochemical process, e.g. extraction of $Li^+$ from NMC. Pay attention that you only monitor cathodic processes.

## 3.2 Diffusion limited electron-transfer – The Randles Circuit (Interface No. 1)

### 1. Cyclic voltammograms (CVs)

* Plot the CVs [only the last cycle at each scan rate] superimposed in a single graph. (Make sure you normalize the measured current by the area of the electrode)

* Describe the plot: what trends do you observe?

* Is the redox couple a reversible? How did you determine this?

* Identify the potentials where the current reaches a maximum (both in the anodic and the cathodic regime) as a function of scan rate.

* Does the peak potential depend on the scan rate? If so, give possible explanations.

* Plot the peak currents (both anodic and cathodic) as a function of the square-root of the scan rate ($i_{peak} \text{ vs. } v^{1/2}$) → [Produce a Randles-Sevcik plot]

* Using the slope of the Randles-Sevcik plots, determine the diffusion coefficient of the redox species.

### 2. Impedance Spectra

* Plot your collected impedance spectra in a Nyquist plot that is normalized to the area of the electrode.

* Determine the electrolyte resistance, the charge-transfer resistance and the double-layer capacitance.

* From the value of the charge-transfer resistance, estimate the $k^0$ value of the redox couple. → How does this value compare to literature values?

---

* Using the data points from the “Warburg tail”, plot both the real and imaginary impedances as a function of $\omega^{1/2}$ and curves use the slope of the curves (which is $A_W$) to determine the diffusion coefficient of the redox species'.

* To determine the diffusion coefficient, use the approximate equation below

$$
A_W = \frac{4RT}{n^2F^2C_O^*D^{1/2}} \quad [3.3.2.1]
$$

* Note that this equation only applies under the following conditions: 1) only the oxidized electroactive species is present in the bulk, 2) the impedance measurement is done at the half-wave potential ($E_{1/2}$) and 3) the surface concentrations of both oxidized and reduced species are estimated based on Langmuir-like behavior.

### 3. Constant-current (Galvanostatic) chronopotentiometry

* Plot all potential transients in a single plot.

* Describe the plot: what trends, if any, do you observe?

* Determine the transition time ($\tau$) of all the transients and from them the quarter-wave potential.

* Using Sand's equation, determine the diffusion coefficient of the ferricyanide.

## 3.3 Diffusion limited electron-transfer with fixed-length diffusion layer (Interface No. 2)

### 1. Linear sweep voltammograms (LSVs)

* Plot all LSVs superimposed in a single graph.

* Describe the plot: what trends do you observe?

* Perform a Koutecky-Levich analysis using the rotation rate-dependent currents at various potentials and determine the kinetic current of this redox couple at those potentials.

* From the obtained kinetic currents, generate a Tafel plot to extract the exchange current density ($j_0$) for this redox couple.

### 2. Impedance Spectra

* Plot all spectra superimposed in a single graph.

* Describe the plot: what trends do you observe?

* By doing a qualitative analysis of the spectra, report the charge-transfer resistance, the double-layer capacitance and the diffusion layer thickness at each rotation rate.

### 4. Discussion and Summary/Conclusions

Use this section to summarize your results from all the measurements performed. Also, you should compare and contrast voltammetric techniques and galvanostatic techniques (in the case of the battery half-cell) and voltammetric vs. impedimetric techniques, in terms of the pros and cons depending on what kind of information one wishes to obtain. You should also address the specific questions below:

* For the LIB cell: How does the graphite material that you measured compare to other graphite materials under similar conditions? (e.g. find a recent literature reference of a graphite||NMC cell and compare its performance to yours). Discuss

---

the overall coulombic efficiency and stability of the material.

* **For Interface No. 1:**

* Were you able to determine the diffusion coefficient of the redox couple?

* How does the value from the CV data compares to the one obtained from the Warburg coefficient to the one obtained from Sand's equation?

* How do these compare to literature values?

* From your charge-transfer resistance value, extract the apparent exchange rate constant ($k^0$) of the redox couple. How does this value compare to those reported in the literature? For deviations from literature values, provide possible explanations.

* Discuss if the measured system is reversible, quasireversible or irreversible.

* **For Interface No. 2:**

* Were you able to determine the diffusion coefficient and $k^0$ of the redox couple?

* How does these values compare to those you obtained for Interface No. 1?

**References**

(1) Goodenough, J. B.; Park, K.-S. *J. Am. Chem. Soc.* **2013**, *135* (4), 1167.
(2) Goodenough, J. B. *J. Solid State Electr.* **2012**, *16* (6), 2019.
(3) Goodenough, J. B. *Acc. Chem. Res.* **2013**, *46* (5), 1053.
(4) Van Noorden, R. *Nature* **2014**, *507* (7490), 26.
(5) Obrovac, M. N.; Chevrier, V. L. *Chem. Rev.* **2014**, *114* (23), 11444.
(6) Whittingham, M. S. *Science* **1976**, *192* (4244), 1126.
(7) Mizushima, K.; Jones, P. C.; Wiseman, P. J.; Goodenough, J. B. *Materials Research Bulletin* **1980**, *15* (6), 783.
(8) Woo, K. C.; Mertwoy, H.; Fischer, J. E.; Kamitakahara, W. A.; Robinson, D. S. *Phys. Rev. B* **1983**, *27* (12), 7831.
(9) Dahn, J. R. *Phys. Rev. B* **1991**, *44* (17), 9170.
(10) Hess, M. Kinetics and stage transitions of graphite for lithium-ion batteries, 2013, pp 1–264.
(11) Dubarry, M.; Liaw, B. Y. *J. Power Sources* **2009**, *194* (1), 541.
(12) Dubarry, M.; Liaw, B. Y.; Chen, M.-S.; Chyan, S.-S.; Han, K.-C.; Sie, W.-T.; Wu, S.-H. *J. Power Sources* **2011**, *196* (7), 3420.
(13) Smith, A. A High Precision Study of Li-Ion Batteries, Dalhousie University: Halifax, Nova Scotia, 2012, pp 1–193.
(14) Bard, A. J.; Faulkner, L. R. *Electrochemical Methods: Fundamentals and Applications*, 2nd ed.; Wiley New York, 2001.
(15) Orazem, M. E.; Tribollet, B. *Electrochemical impedance spectroscopy*; John Wiley & Sons: Hoboken, New Jersey, 2008.
(16) Brown, A. P.; Koval, C.; Anson, F. C. *J. Electroanal. Chem.* **1976**, *72* (3), 379.

**Appendix A.** Guidelines for written report on *Electrochemical Methods for Interface Characterization*

Sections from the text below have been adapted from: Scientific Papers. In *The ACS Style Guide*, 3rd Edition; Coghill, A.M., Garson, L.R., Eds.; Oxford University Press: New York, NY, 2006, pp. 17–26.

---

In general, scientific papers are organized into a standard format that includes (but is not limited to): abstract, introduction, experimental details or theoretical basis, results, discussions and conclusions. This document will briefly describe what is included, as well as suggestions on writing styles and word usage in each section.

1. **Abstract:** Most scientific publications require an informative summary (abstract) that contain between 80 and 250 words in total. Avoid acronyms and abbreviations, but if needed define at first use within the abstract and at first use within the text. It is recommended to write the abstract last (once the full article is written) so that it fully portrays the contents of the manuscript. Moreover, the **verb tenses used per sentence depend on which section of the manuscript it refers to**

2. **Introduction:** This section should start with presenting broad statements regarding the scientific problem and the research motivation. Related and relevant previous work should be presented, accurately cited and discussed within the context of the current project. It is important to highlight how the current work is different and contributes new insights from the previous work. Note that the verb tenses used in this section depend on the type of sentence:

a. When stating a widely accepted fact, the **present tense** is appropriate:

i. "DNA **is** composed of four nucleotides."

b. When referring to a previous study, the **present perfect tense** is used. This tense demonstrates that the work was done in the past but its results are still relevant and applicable in the present (i). It also applies to describe an event that occurred in the past but continues in the present (ii).

i. "Brown and coworkers **have shown** that [...]."

ii. "Patients with XYZ syndrome **have been surveyed** for the last 10 years.

c. When referring to a specific figure (data set), result or paper in a previously published work, the **present tense** is used:

i. "The results of their study **indicate** that ion-dipole interactions **are** dominant".

d. When referring specifically to previously published methodology, the **past tense** is appropriate:

i. Li and coworkers **used** photovoltage decay measurements and **distinguished** various simultaneous interfacial phenomena.

3. **Experimental Details and/or Theoretical Basis** (also known as "*Experimental Methods*", "*Methodology*", "*Materials and Methods*") : This section should provide all the details necessary so that another experienced scientist may repeat the work and obtain similar results. Identify the materials used with correct chemical formulas, chemical names and provide information on purity. A description of non-standard and/or home-built apparatus should be provided. If standard equipment is used, provide the company name and model number in parenthesis. Describe all procedures used, unless standard. In the case of the latter, provide appropriate citations. Note and emphasize any hazards during the procedures, e.g. explosive, pyrophoric, toxicity properties. This section should be **written entirely in past tense** because it reports what has been done during the course of the study.

4. **Results:** Show and summarize all collected data and statistical treatment applied. Use equations, figures and tables when necessary for clarity and brevity. Include only

---

relevant data but with enough detail to justify your conclusions. Since the experiments
done to obtain the data in this section have been completed prior to writing the
manuscript, this section should be written (for the most part) in the **past tense** (a & b).
However, similar to the Introduction section, when referring to a specific figure, table,
section, etc. the present tense is necessary (c).

---

a. "No fluorescence **was detected** in the control sample".
b. "Titration analysis **showed** a final water concentration of (10 ± 1) ppm for the sample."
c. "Figure 1.2 **shows** the particle size distribution of the as-synthesized IrO₂ nanoparticles."

**5. Discussion:** In this section, the results from the study are interpreted and compared. Be objective when establishing connections between your data and describe the logical implications. Consider these guide questions when writing this section: 1) Was the problem resolved?; 2) What has been contributed?; 3) What new insights are these results providing?; 4) How do these results compare to previously published work and the current knowledge in the field? Similar to the "Introduction" section, the verb tense depends on the purpose of the sentence.

a. When referring to your results use the **past tense**, but when establishing a conclusion use the **present tense**:

i. "The extinction coefficient of Complex X at 310 nm **was determined** earlier to be 3300 M⁻¹ cm⁻¹, which **suggests** a ligand-to-metal charge transfer absorption."

ii. "Taken together the potential-dependent spectra **suggest** that the absorbance at 220 nm **is** predominantly from an IrIV species and that the 318 nm band **is** a combination of both IrIII and IrIV species."

**6. Conclusions:** The purpose of this section is to emphasize the main conclusions of the work and how they correlate to the original problem. Do not repeat the content in the "Discussion" section. This section should be written in the **present tense**:

a. "Monomeric Ir(III/IV) anions strongly **adsorb** to the surface of metal oxide photoelectrodes and **promote** IrO₂ nanoparticle deposition. However, the monomer **is** detrimental to the performance of photoanodes in IrO₂-catalyzed water-splitting systems."

**7. References:** This section can either appear as a whole at the end of a manuscript or as footnotes throughout the text. The formatting required for the references depends on the publishing company. However, regardless of the formatting used, the information provided per citation and/or reference must allow any reader to find its content.

For your report, please include the following sections: **Introduction**, **Experimental**, **Results**, **Discussion** and **References**. Format all your references as per The ACS Style Guide, 3rd Edition, Chapter 14 (see Table attached)