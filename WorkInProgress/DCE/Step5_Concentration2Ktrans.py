import numpy as np
import osipy
import nibabel as nib
from osipy.common.types import AIFType

def concentration_to_ktrans(pathConcentration, pathOut_Ktrans, pathOut_vp, AIF):
    """Fit a Patlak model to concentration values and return Ktrans and vp maps.

    Parameters
    ----------
    pathConcentration : path to NIfTI file
        Path to the NIfTI file containing contrast agent concentration values.

    pathOut_Ktrans : path to NIfTI file
        Path to the output NIfTI file where the fitted Ktrans map will be saved.

    pathOut_vp : path to NIfTI file
        Path to the output NIfTI file where the fitted vp map will be saved.

    AIF : osipy.AIF object
        The arterial input function to be used for fitting the Patlak model.

    Returns
    -------
    dict
        The fitted ``Ktrans`` and ``vp`` values, fit-quality outputs, the AIF,
        and the original OSIPI fit result.
    """    

    # Load the concentration NIfTI file
    dce_dataset = osipy.load_nifti(pathConcentration)


    # DCE timing
    n_timepoints = dce_dataset.shape[-1]
    temporal_resolution = 39.6279  # seconds between frames
    time = np.arange(n_timepoints) * temporal_resolution


    # aif = osipy.ParkerAIF()(time)

    aif_concentration = np.asarray(AIF, dtype=float).reshape(-1)

    aif = osipy.ArterialInputFunction(
        time=time,
        concentration=aif_concentration,
        aif_type=AIFType.MEASURED,
    )


    # print(f"AIF peak: {aif.concentration.max():.2f} mM at {time[aif.concentration.argmax()]:.1f} s")

    # Fit Patlak model
    result = osipy.fit_model(
        "patlak",
        concentration=dce_dataset.data,
        aif=aif,
        time=time,
        fit_delay=True
    )


    # Extract parameter maps (result is a DCEFitResult object)
    ktrans = result.parameter_maps["Ktrans"].values
    vp = result.parameter_maps["vp"].values
    r_squared = result.r_squared_map
    quality_mask = result.quality_mask

    # Save Ktrans and vp maps as NIfTI files
    img = nib.Nifti1Image(ktrans.astype(np.float32),
                        affine=dce_dataset.affine)

    nib.save(img, pathOut_Ktrans)

    img = nib.Nifti1Image(vp.astype(np.float32),
                        affine=dce_dataset.affine)

    nib.save(img, pathOut_vp)    