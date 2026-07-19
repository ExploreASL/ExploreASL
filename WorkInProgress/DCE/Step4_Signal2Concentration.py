from osipy.dce import DCEAcquisitionParams
import osipy
import numpy as np
import nibabel as nib

def signal_to_concentration(signal, pathT1_map, pathVFA):
    """
    Convert DCE-MRI signal to contrast agent concentration using T1 mapping.

    Parameters:
    signal : array-like
        Signal values. May be a 1D time series or a multidimensional array
        whose final dimension represents time.
    - t1_map: ParameterMap containing T1 relaxation times.

    Returns:
    - concentration: Numpy array of contrast agent concentrations.
    """



    # >>>>>>>>>>>> VFA/T1-MAPPING BIT

    # Load VFA images for T1 mapping
    # These should be acquired at different flip angles
    vfa_data = osipy.load_nifti(pathVFA)

    # Define acquisition parameters
    flip_angles = [2, 12]  # degrees, this should not be a numpy array; which would crash "if not params.flip_angles"
    tr = 415  # Repetition time in ms; but this doesn't seem right; there are two different TRs?

    vfa_data.acquisition_params = DCEAcquisitionParams(
        tr=tr, flip_angles=flip_angles
    )

    # Compute T1 map using Variable Flip Angle method
    t1_result = osipy.compute_t1_map(vfa_data, method="vfa")

    # >>>>>>>>>> NOW WE HACK THIS USING A DUMMY T1 value, the DESPOT HIFI method still needs implementing

    # t1_result.t1_map.values = t1_map.data
    # t1_result.t1_map.values[~np.isfinite(t1_result.t1_map.values)] = 1100  # Replace NaN or Inf with 1100 ms

    t1_result.t1_map.values[:] = 1100 # average WM & GM T1 value at 3T


    # >>>>>>>>>>>>>>>>>> Handle MATLAB vectors shaped [1, T], [T, 1], or Python arrays [T].
    signal = np.asarray(signal, dtype=float)
    input_was_1d = False

    
    if signal.ndim == 1 or (
        signal.ndim == 2 and 1 in signal.shape
    ):
        signal = signal.reshape(-1)
        signal = signal.reshape(1, 1, 1, -1)
        input_was_1d = True    

        signal = np.broadcast_to(
            signal.reshape(1, 1, 1, -1),
            t1_result.t1_map.shape + (signal.size,),
        )




    # Define acquisition parameters
    acq_params = DCEAcquisitionParams(
        tr=3.44,            # ms
        te=1.68,            # ms
        flip_angles=[15],  # degrees (DCE flip angle)
        baseline_frames=3,
        relaxivity=4.5,    # mM^-1 s^-1 — Gd-DTPA in water at 3T
                           # (Sasaki 2005, PMID 16462135: R1_3T_water=4.50).
                           # In vivo plasma r1 at 3T is lower; override for
                           # agent-specific in-vivo accuracy.
    )




    # Convert signal to concentration
    concentration = osipy.dce_signal_to_concentration(
        signal=signal,
        t1_map=t1_result.t1_map,
        acquisition_params=acq_params,
    )

    print(f"Concentration shape: {concentration.shape}")

    ## # Check concentration range
    ## valid_conc = concentration[mask]
    ## print(f"Concentration range: {valid_conc.min():.3f} - {valid_conc.max():.3f} mM")

    # Return a normal 1D vector for a single input time series.
    if input_was_1d:
        concentration = np.asarray(concentration)[0, 0, 0, :] # assuming that the T1 map is the same for every voxel

    return concentration


def nifti_signal_to_concentration(pathNifti_In, pathT1_map, pathVFA, pathNifti_Out):
    """
    Load DCE signals from NIfTI, convert them, and save a NIfTI result.
    """
    dce_dataset = osipy.load_nifti(pathNifti_In)

    concentration = signal_to_concentration(
        signal=dce_dataset.data,
        pathT1_map=pathT1_map,
        pathVFA=pathVFA,
    )

    output_image = nib.Nifti1Image(
        concentration.astype(np.float32),
        affine=dce_dataset.affine,
    )
    nib.save(output_image, pathNifti_Out)

    return concentration