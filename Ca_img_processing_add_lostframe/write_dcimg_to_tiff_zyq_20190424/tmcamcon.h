#include <stdint.h>

uint32_t TMCC_OPENDCIMGFILE(uint32_t *filehandle, char path[], int32_t *totalframes, int32_t *Erval);
uint32_t TMCC_CLOSEDCIMGFILE(uint32_t filehandle, int32_t *Erval);
uint32_t TMCC_GETDCIMGFRAMEINFO(uint32_t index, int32_t frameindex, int32_t *width, int32_t *height, int32_t *Erval);
uint32_t TMCC_GETDCIMGFRAMEDATA_A(uint32_t index, int32_t frame, uint16_t imgarray[], double *timestamp, int32_t *Erval);
uint32_t TMCC_REPORTERROR(int32_t error, char Erval[]);
