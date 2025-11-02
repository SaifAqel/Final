import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

def main(csv_path):
    df = pd.read_csv(csv_path)

    outdir = os.path.join(os.path.dirname(csv_path), "fig")
    os.makedirs(outdir, exist_ok=True)

    x = "x[m]"
    def save(fig, name):
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, name), dpi=200)
        plt.close(fig)

    # 1) gas T and water T
    fig = plt.figure()
    for col in ["gas_T[K]", "water_T[K]"]:
        if col in df: plt.plot(df[x], df[col], label=col)
    plt.xlabel(x); plt.ylabel("Temperature [K]"); plt.legend()
    save(fig, "01_temp_gas_water.png")

    # 2) qprime and UA_prime
    fig = plt.figure()
    for col in ["qprime[W/m]", "UA_prime[W/K/m]"]:
        if col in df: plt.plot(df[x], df[col], label=col)
    plt.xlabel(x); plt.ylabel("q′ [W/m], UA′ [W/K/m]"); plt.legend()
    save(fig, "02_qprime_UAprime.png")

    # 3) enthalpy: gas_h and water_h
    fig = plt.figure()
    for col in ["gas_h[J/kg]", "water_h[J/kg]"]:
        if col in df: plt.plot(df[x], df[col], label=col)
    plt.xlabel(x); plt.ylabel("Specific enthalpy [J/kg]"); plt.legend()
    save(fig, "03_h_gas_water.png")

    # 4) HTC: h_gas and h_water
    fig = plt.figure()
    for col in ["h_gas[W/m^2/K]", "h_water[W/m^2/K]"]:
        if col in df: plt.plot(df[x], df[col], label=col)
    plt.xlabel(x); plt.ylabel("Heat transfer coefficient [W/m²·K]"); plt.legend()
    save(fig, "04_htc_gas_water.png")

    # 5) gas cp, mu, k, rho
    fig = plt.figure()
    for col in ["gas_cp[J/kg/K]", "gas_mu[Pa*s]", "gas_k[W/m/K]", "gas_rho[kg/m^3]"]:
        if col in df: plt.plot(df[x], df[col], label=col)
    plt.xlabel(x); plt.ylabel("Gas properties"); plt.legend()
    save(fig, "05_gas_props.png")

    # 6) water cp, mu, k, rho
    fig = plt.figure()
    for col in ["water_cp[J/kg/K]", "water_mu[Pa*s]", "water_k[W/m/K]", "water_rho[kg/m^3]"]:
        if col in df: plt.plot(df[x], df[col], label=col)
    plt.xlabel(x); plt.ylabel("Water properties"); plt.legend()
    save(fig, "06_water_props.png")

    # 7) boiling vs x (as 0/1)
    fig = plt.figure()
    if "boiling" in df:
        y = df["boiling"].astype(str).str.upper().isin(["TRUE","1","T","YES"]).astype(int)
        plt.step(df[x], y, where="post", label="boiling (1=true)")
    plt.ylim(-0.1, 1.1)
    plt.xlabel(x); plt.ylabel("Boiling flag")
    save(fig, "07_boiling_flag.png")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python plot_hx.py <data.csv>")
        sys.exit(1)
    main(sys.argv[1])
