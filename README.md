# ARSCA: Atmospheric Radiation Simulation Computation Application
Correspondence: antti.mikkonen@fmi.fi

## Setting up ARSCA

Clone the repository wherever you wish.
Then add the following lines (you probably need to modify the path accordingly) to the start of your script:
```
import sys
arsca_path = '/home/amikko/Projects/arsca'
sys.path.insert(0,arsca_path)
```

## Examples

The files `model_validation.py` and `snow_covers.py` contain the scripts with which you may recreate the RT results presented in the paper titled "Non-Lambertian snow surface reflection models for simulated top-of-atmosphere radiances in the NIR and SWIR wavelengths".

TODO: Rest of the README.
