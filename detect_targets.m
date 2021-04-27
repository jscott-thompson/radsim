function result = detect_targets(range_doppler_map,detector_params)
%DETECT_TARGETS Detect targets in range Doppler map using CFAR with
%clustering

exceedances = detect_cfar_exceedances(range_doppler_map,detector_params);
result = cluster_exceedances(exceedances,detector_params);
