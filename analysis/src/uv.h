#ifndef UV_H
#define UV_H

float convert_UVperFreq_to_UVmag(float UVperFreq);
float get_UVperFreq_after_time(float timeInYr);
float get_UVcorrectionFactor(float timeOffset, float timeBegin, float timeEnd, float timeBreak);
float calc_UV_lum_SMH(int numSnaps, float *starmassHistoryAll, int thisGal, int currSnap, float *corrFactor);
float *calc_UV_lum_at_snap(dconfObj_t simParam, int numGal, int numSnaps, float *starmassHistoryAll, int currSnap);
float *calc_UV_mag(dconfObj_t simParam, int numGal, int numSnaps, float *starmassHistoryAll, int currSnap);

#endif