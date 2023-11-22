import ipsuite as ips
import zntrack
from src import LoadModel

project = ips.Project(automatic_node_names=True)

thermostat = ips.calculators.LangevinThermostat(
    temperature=298.15, friction=0.01, time_step=0.5
)

with project:
    model = LoadModel(model_path="2023-08-14-mace-universal.model")

ramp_density = ips.calculators.RescaleBoxModifier(
    density=1086
)

ramp_temperature = ips.calculators.TemperatureRampModifier(
    start_temperature=250, temperature=500,
)

thermostat = ips.calculators.LangevinThermostat(
    temperature=298.15, friction=0.01, time_step=0.5
)

with project.group("BMIM_Cl"):
    cation = ips.configuration_generation.SmilesToAtoms("CCCCN1C=C[N+](=C1)C")
    anion = ips.configuration_generation.SmilesToAtoms("[Cl-]")

    single_structure = ips.configuration_generation.Packmol(
            data=[cation.atoms, anion.atoms],
            count=[1, 1],
            density=1000,
            pbc=False
    )

    structure = ips.configuration_generation.Packmol(
        data=[single_structure.atoms],
        count=[10],
        density=700,
    )
    
    geo_opt = ips.calculators.ASEGeoOpt(
        model=model,
        data=structure.atoms,
        data_id=-1,
        optimizer="FIRE",
        run_kwargs={"fmax": 0.5},
    )

    density_md = ips.calculators.ASEMD(
        data=geo_opt.atoms,
        data_id=-1,
        model=model,
        modifier=[ramp_density],
        thermostat=thermostat,
        steps=1000,
        sampling_rate=10,
    )

    md = ips.calculators.ASEMD(
        data=geo_opt.atoms,
        data_id=-1,
        model=model,
        modifier=[ramp_temperature],
        thermostat=thermostat,
        steps=500_000,
        sampling_rate=100,
    )

    selection = ips.configuration_selection.UniformTemporalSelection(
        data=md.atoms,
        n_configurations=50,
    )

    cp2k = ips.calculators.CP2KSinglePoint(
        data=selection.atoms,
        cp2k_params="config/cp2k.yaml",
        cp2k_files=["GTH_BASIS_SETS", "GTH_POTENTIALS", "dftd3.dat"],
    )


project.build()
