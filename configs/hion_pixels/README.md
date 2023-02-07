Config to run the pixel track reconstruction for HIon events. Usage: `cmsRun step3_hion_patatrack.py` with options
- `n=N` for number of events;
- `threads=T` for the number of threads and streams;
- `timing=True` not to run the validation;
- `gpu=True` to run patatrack pixel tracks on GPU (if `False` it is forced to run patatrack pixel tracks on CPU);

or `cmsRun step3_hion_leg.py`to run the legacy pixel tracks reco (no `gpu` option of course). On the `*inDQM*` file you get as output run `harvestTrackValidationPlots.py DQMFILE.root -o HARVESTNAME.root` to harvest the DQM plots and on top of  
`HARVESTNAME.root` run `makeTrackValidationPlots.py -o folder_youwant` to get the validation plots and the html to navigate them under `folder_youwant/`.
