import ipsuite as ips
import zntrack
from src import LoadModel

project = ips.Project(automatic_node_names=True)

thermostat = ips.calculators.LangevinThermostat(
    temperature=298.15, friction=0.01, time_step=0.5
)

with project:
    model = LoadModel(model_path="2023-08-14-mace-universal.model")

with project.group("H2O"):
    h2o = ips.configuration_generation.SmilesToAtoms("O")
    box = ips.configuration_generation.Packmol(data=[h2o.atoms], count=[10], density=997)

    geo_opt = ips.calculators.ASEGeoOpt(data=box.atoms, model=model)

    md = ips.calculators.ASEMD(
        data=geo_opt.atoms,
        data_id=-1,
        model=model,
        thermostat=thermostat,
        steps=50000,
        sampling_rate=10,
    )

with project.group("H2O_H2SO4"):
    h2o = ips.configuration_generation.SmilesToAtoms("O")
    h2so4 = ips.configuration_generation.SmilesToAtoms("OS(=O)(=O)O")
    box = ips.configuration_generation.Packmol(
        data=[h2o.atoms, h2so4.atoms], count=[100, 2], density=997
    )
    geo_opt = ips.calculators.ASEGeoOpt(data=box.atoms, model=model, run_kwargs={"fmax": 0.5})

    md = ips.calculators.ASEMD(
        data=geo_opt.atoms,
        data_id=-1,
        model=model,
        thermostat=thermostat,
        steps=50000,
        sampling_rate=1,
    )

with project.group("DCM_H2SO4"):
    dcm = ips.configuration_generation.SmilesToAtoms("C(Cl)Cl")
    h2so4 = ips.configuration_generation.SmilesToAtoms("OS(=O)(=O)O")
    box = ips.configuration_generation.Packmol(
        data=[dcm.atoms, h2so4.atoms], count=[100, 2], density=1330
    )
    geo_opt = ips.calculators.ASEGeoOpt(data=box.atoms, model=model, run_kwargs={"fmax": 0.5})

    md = ips.calculators.ASEMD(
        data=geo_opt.atoms,
        data_id=-1,
        model=model,
        thermostat=thermostat,
        steps=50000,
        sampling_rate=1,
    )

with project.group("H2O_AcAc"):
    h2o = ips.configuration_generation.SmilesToAtoms("O")
    acac = ips.configuration_generation.SmilesToAtoms("CC(=O)OC(=O)C")

    box = ips.configuration_generation.Packmol(
        data=[h2o.atoms, acac.atoms], count=[10, 1], density=997
    )

    geo_opt = ips.calculators.ASEGeoOpt(data=box.atoms, model=model)

    md = ips.calculators.ASEMD(
        data=geo_opt.atoms,
        data_id=-1,
        model=model,
        thermostat=thermostat,
        steps=200000,
        sampling_rate=1,
    )

with project.group("BMIM_PF6"):
    cation = ips.configuration_generation.SmilesToConformers(smiles="CCCCN1C=C[N+](=C1)C", numConfs=200)
    anion = ips.configuration_generation.SmilesToConformers(smiles="F[P-](F)(F)(F)(F)F", numConfs=200)
    single_structure = ips.configuration_generation.MultiPackmol(
        data=[cation.atoms, anion.atoms],
        count=[1, 1],
        density=1380,
        pbc=False,
        n_configurations=200,
    )

    structure = ips.configuration_generation.MultiPackmol(
        data=[single_structure.atoms],
        count=[10],
        density=1380,
        n_configurations=20,
    )

with project.group("ASA"):
    salicylic_acid = ips.configuration_generation.SmilesToAtoms("C1=CC=C(C(=C1)C(=O)O)O")
    acac = ips.configuration_generation.SmilesToAtoms("CC(=O)OC(=O)C")
    h2so4 = ips.configuration_generation.SmilesToAtoms("OS(=O)(=O)O")

    box = ips.configuration_generation.Packmol(
        data=[salicylic_acid.atoms, acac.atoms, h2so4.atoms],
        count=[20, 20, 2],
        density=1000, # guess
    )

    geo_opt = ips.calculators.ASEGeoOpt(data=box.atoms, model=model)

    md = ips.calculators.ASEMD(
        data=geo_opt.atoms,
        data_id=-1,
        model=model,
        thermostat=thermostat,
        steps=200000,
        sampling_rate=1,
    )


project.build()
