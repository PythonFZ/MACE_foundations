import ipsuite as ips
import zntrack
from src import LoadModel, LoadModelDispersion

project = ips.Project(automatic_node_names=True)

thermostat = ips.calculators.LangevinThermostat(
    temperature=298.15, friction=0.01, time_step=0.5
)

with project:
    model = LoadModel(model_path="2023-08-14-mace-universal.model")
    model_dispersion = LoadModelDispersion(model_path="2023-08-14-mace-universal.model")

ramp_density = ips.calculators.RescaleBoxModifier(
    density=1196
)

thermostat = ips.calculators.LangevinThermostat(
    temperature=303.15, friction=0.01, time_step=0.5
)

mapping = ips.geometry.BarycenterMapping(data=None)

with project.group("BMIM_BF4"):
    cation = ips.configuration_generation.SmilesToAtoms("CCCCN1C=C[N+](=C1)C")
    anion = ips.configuration_generation.SmilesToAtoms("[B-](F)(F)(F)F")

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
        model=model_dispersion,
        data=structure.atoms,
        data_id=-1,
        optimizer="FIRE",
        run_kwargs={"fmax": 0.5},
    )

    density_md = ips.calculators.ASEMD(
        data=geo_opt.atoms,
        data_id=-1,
        model=model_dispersion,
        modifier=[ramp_density],
        thermostat=thermostat,
        steps=1000,
        sampling_rate=10,
    )

    geo_opt = ips.calculators.ASEGeoOpt(
        model=model_dispersion,
        data=structure.atoms,
        data_id=-1,
        optimizer="FIRE",
        run_kwargs={"fmax": 0.5},
    )

    volume_scan = ips.analysis.BoxScale(
            data=geo_opt.atoms,
            mapping=mapping,
            model=model,
            start=0.9,
            stop=2.0,
            num=50,
            data_id=-1,
        )

    cp2k = ips.calculators.CP2KSinglePoint(
            data=volume_scan.atoms,
            cp2k_params="config/cp2k.yaml",
            cp2k_files=["GTH_BASIS_SETS", "GTH_POTENTIALS", "dftd3.dat"],
        )

    prediction = ips.analysis.Prediction(cp2k.atoms, model)
    ips.analysis.ForceDecomposition(prediction)

project.build()