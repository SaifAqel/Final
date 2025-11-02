from new_loader import load_all
from combustion.combustor import Combustor
from heat.runner import run_hx

stages, air, fuel, water, drum, operation = load_all(
    stages_path="config/stages.yaml",
    air_path="config/air.yaml",
    fuel_path="config/fuel.yaml",
    water_path="config/water.yaml",
    drum_path="config/drum.yaml",
    operation_path="config/operation.yaml",
)

svc = Combustor(air, fuel, operation["excess_air_ratio"])
combustion_results = svc.run()
print(combustion_results)

result = run_hx(
    stages_yaml="config/stages.yaml",
    streams_yaml="config/streams.yaml",
    drum_yaml="config/drum.yaml",
    target_dx="0.1 m",
    min_steps=20,
    max_steps=400,
    max_passes=20,
    tol_Q="1e-3 W",
    tol_end="1e-3 J/kg",
    write_csv=True,
    outdir="results"
)
