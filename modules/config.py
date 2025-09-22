import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "Data"
FIG_DIR = ROOT / "figures"
RES_DIR = ROOT / "results"

# Glassy carbon working electrode (Tasks 3.2/3.3)
GC_DIAMETER_MM = 5.5
GC_RADIUS_CM = (GC_DIAMETER_MM / 10.0) / 2.0
GC_AREA_CM2 = math.pi * (GC_RADIUS_CM ** 2)

# Electrolyte / chemistry (Tasks 3.2 & 3.3)
N_ELECTRONS = 1  # Ferri/ferrocyanide is 1e-
T_K = 298.15
R_GAS = 8.314462618  # J/mol/K
F_CONST = 96485.33212  # C/mol

# 2.3 mM ferricyanide, used in Tasks 3.2 and 3.3
C_BULK_MOL_PER_L = 2.3e-3
C_BULK_MOL_PER_CM3 = C_BULK_MOL_PER_L / 1000.0  # mol/cm^3
C_BULK_mM = 2.3  # for Sand's eq. numeric form in the script (mM)

# Rotating-disk hydrodynamics (Task 3.3)
KCL_MOLAR_MASS_GPMOL = 74.5513      # g/mol
KCL_DENSITY_2M_25C_G_PER_ML = 1.10  # g/mL
M_KCL_MOLARITY = 2.0                 # mol/L

# 1) M -> m
m_KCL_molal = (1000.0 * M_KCL_MOLARITY) / (
    1000.0 * KCL_DENSITY_2M_25C_G_PER_ML - M_KCL_MOLARITY * KCL_MOLAR_MASS_GPMOL
)

# 2) Interpolate viscosity between tabulated molalities (2.0 and 2.5 mol/kg at 25 C)
m1, nu1_mm2s = 2.0, 0.8425  # 2.0 mol/kg
m2, nu2_mm2s = 2.5, 0.8324  # 2.5 mol/kg
nu_mm2s_interp = nu1_mm2s + (m_KCL_molal - m1) * (nu2_mm2s - nu1_mm2s) / (m2 - m1)

# 3) Convert to cm^2/s
NU_CMS2 = nu_mm2s_interp * 0.01


def ensure_dirs() -> None:
    """Create output directories if needed."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    RES_DIR.mkdir(parents=True, exist_ok=True)


def hydrodynamics_summary() -> str:
    """Return the human-readable summary that was previously printed on import."""
    return (
        "[nu] 2.0 M KCl @25C: rho="
        f"{KCL_DENSITY_2M_25C_G_PER_ML:.3f} g/mL -> m={m_KCL_molal:.3f} mol/kg -> "
        f"nu~{nu_mm2s_interp:.4f} mm^2/s -> NU_CMS2={NU_CMS2:.6f} cm^2/s"
    )

