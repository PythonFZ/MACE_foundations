import zntrack
import torch 
import functools
import ase
from mace.calculators.foundations_models import mace_mp
from mace.calculators import MACECalculator
from ase.calculators.singlepoint import SinglePointCalculator
import tqdm

class LoadModel(zntrack.Node):
    model_path: str = zntrack.deps_path()

    def run(self) -> None:
        pass
    
    @functools.lru_cache()
    def get_model(self, **kwargs) -> torch.nn.Module:
        with self.state.fs.open(self.model_path, 'rb') as f:
            if kwargs.get("device", None) is None:
                return torch.load(f)
            else:
                model = torch.load(f, map_location=kwargs["device"])
                model.to(kwargs["device"])
                return model
    
    def get_calculator(self, device=None, default_dtype="float32", **kwargs) -> MACECalculator:
        return mace_mp(model=self.model_path, dispersion=False)

    def predict(self, atoms_list: list[ase.Atoms]) -> list[ase.Atoms]:
        calc = self.get_calculator()
        for atoms in tqdm.tqdm(atoms_list):
            atoms.set_calculator(calc)
            energy = atoms.get_potential_energy()
            forces = atoms.get_forces()
            atoms.set_calculator(
                SinglePointCalculator(atoms, energy=energy, forces=forces)
            )
        return atoms_list


class LoadModelDispersion(LoadModel):
    def get_calculator(self, device=None, default_dtype="float32", **kwargs) -> MACECalculator:
        return mace_mp(model=self.model_path, dispersion=True)