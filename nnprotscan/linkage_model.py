import os

VALID_LINKAGE_MODELS = ["amide"]
MATCH_CST_PATH = os.path.join(os.path.dirname(__file__), "data", "rosetta")

class LinkageModel:
    def __init__(self, linkage_model_name):
        """
        Utility class to handle constraints and linking geometry definitions
        for the rosetta matcher and MD simulation protocols.
        
        Args:
            linkage_model_name (str): A valid linkage model name. Currently
            supported types are amide and sulfide.
        """
        self.linkage_model_name = linkage_model_name
        self._validate()
        
        self._extra_atoms = {
            "amide": {
                "protein": ["HZ1", "HZ2"],   
                "protractor": ["H1"]
            }
        }
        
        self._bond_atoms = {
            "amide": [
                {"protein": "NZ", "protractor": "C1"}
            ]
        }
        
        self._match_res = {
            "amide": "K"
        }
    
    def _validate(self):
        """
        Check if this is a valid linkage model. 
        Returns:
            bool: True if the linkage model is implemented and has associated
            files present.
        """
        if self.linkage_model_name not in VALID_LINKAGE_MODELS:
            raise NotImplementedError(self.linkage_model_name, "not yet!")

    @property
    def match_constraint_file(self):
        """
        Returns the absolute path for the (rosetta enzdes style) match 
        constraint file corresponding to this linkage model
        """
        return os.path.join(MATCH_CST_PATH, self.linkage_model_name + ".cst")
    
    @property
    def extra_atoms(self):
        return self._extra_atoms[self.linkage_model_name]
    
    @property
    def bond_atoms(self):
        return self._bond_atoms[self.linkage_model_name]
    
    @property
    def match_res(self):
        return self._match_res[self.linkage_model_name]

    @property
    def available_models(self):
        return VALID_LINKAGE_MODELS
    