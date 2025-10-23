/* From https://github.com/mzucker/flow_solver */

#ifndef _WIN32
#include <unistd.h>
#include <sys/time.h>
#else
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <stdio.h>
#include <stdint.h>
#include <locale.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <time.h>

double now() {
  
#ifdef _WIN32
	union {
		LONG_LONG ns100; /*time since 1 Jan 1601 in 100ns units */
		FILETIME ft;
	} now;
	GetSystemTimeAsFileTime (&now.ft);
	return (double)now.ns100 * 1e-7; // 100 nanoseconds = 0.1 microsecond
#else
	struct timeval tv;
	gettimeofday(&tv, 0);
	return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
#endif

}
