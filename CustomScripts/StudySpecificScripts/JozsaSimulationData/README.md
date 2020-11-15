These files were created to process the EPAD data set on a machine with matlab license and extract statistics for tuning a porous finite element model of the human brain.
The process includes the following steps in matlab:

1; Copy and reorganise the EPAD data set:
```
prepare_folders(EPAD_data_dir,copied_data_dir,data_par_template_file,quality,xASL_path,outfile)
```

2; Compute CBF:
```
process_data(start_num,end_num,data_dir,xASL_path,patient_info)
```

3; Obtain statistics for model tuning:
```
extract_stats_for_virtual_brain_model_fct(data_dir,kernel_size,PVC_mask_treshold,quantile_lim)
```

For further details, see the corresponding functions.

The quantile and percentile functions (see quantile_fct.m and prctile_fct.m) are needed only if an older matlab version is used and they are not available by default.

### ToDo
- [x] Test implementation
- [ ] Implement affine transformation matrix computation approximating patient brain with virtual brain
