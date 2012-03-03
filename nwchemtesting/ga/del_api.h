#ifndef DEL_API_H
#define DEL_API_H

struct del_t {
    /* These are the GA handles to the respective GAs for subtrahend,
     * minuend, and deltas - provided as needed. */
    int subHandle;
    int minHandle;
    int delHandle;
    
    /* The size, in bytes, of the input arrays. */
    unsigned long int minSize;
    unsigned long int delSize;
    unsigned long int subSize;

    /* These are provided with function calls that write their
     * output to disk, corresponding to an output delta, or an output
     * recovered minuend. */
    char * subName;
    char * delName;
    char * minName;

    FILE * subFP;
    FILE * minFP;
    FILE * delFP;

    float thresholdExp;
    int doFilter;
    
} ;

/* Delta compression uses in-memory zlib compression. */

/* To-be-added functionality: Support for FPC in-memory compression. */

/* Computes the deltas between these two global arrays of 
 * doubles  and writes the compressed output to  "delName". If
 * "doFilter" is set to 1, then the delta compressor will use the 
 * specified "thresholdExp". 
 *
 * Returns 0 on success, -1 on failure. */
int computeDeltaFromMemToDisk(int subHandle, int minHandle, unsigned long int size, char * delName, int doFilter, float thresholdExp);

/* Takes the GA subHandle and recovers the checkpoint ahead of it, 
 * recalculted using the on-disk delta file whose name is "delName".
 * 
 * Returns a handle to the recovered subtrahend checkpoint  */
int recoverDeltaFromDiskToMem(int subHandle, unsigned long int size, char * delName); 
#endif
