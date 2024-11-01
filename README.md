## Magnetic field evolution in hot Exoplanets

The goal of this project is to use the [Christensen](https://ui.adsabs.harvard.edu/abs/2009Natur.457..167C/abstract) relations in order to quantify the planetary dynamo of hot Exoplanets. 
### Prerequisites
1. `MESA` version 12115, found [here](https://sourceforge.net/projects/mesa/files/releases/mesa-r12115.zip/download). Do **not** try to use any other `MESA` version. The associated SDK is [20190830](user.astro.wisc.edu/~townsend/resource/download/mesasdk/old/mesasdk-x86_64-linux-20190830.tar.gz) (clicking this will begin the download). You may try to use a different SDK, but we do not recommend this.
2. Daria Kubyshkina's `MESA` inlists. Found [here](https://zenodo.org/records/4022393)
3. `mesa_reader` to interface with `MESA` data in `python`. Found [here](https://github.com/wmwolf/py_mesa_reader). 
4. `Mors` to account for stellar evolution. Found [here](https://github.com/ColinPhilipJohnstone/Mors).

### Workflow
Follow the inlist instructions to generate `MESA` files.  Consult `The_MESA_log.pdf` for detailed instructions regarding installation and working with the inlists.

After getting them to run, edit `run_star_extras_evol_Mors0.f` with the contents of `Utilities/customsaving.f`. This is necessary in order to prompt `MESA` to adequately sample every period in the planetary lifetime.  To get all the neccecery outputs from `MESA`, look at `Utilities/profile_columns_backup.f`.

The code expects to find the `MESA` output in folder named `data`, in the same directory as the repository, like so:

```
.
├── src # contents of the repo
├── data 
```
