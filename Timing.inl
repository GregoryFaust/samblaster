/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" -*-

    Project: Genome Alignment
    Author:  Greg Faust

    File:    Timing.inl       Contains inline functions and macros for timing of code sections.
*/

#include <sys/times.h>
#include <sys/resource.h>
#include <time.h>

// Define booleans if they are not already.
#ifndef BOOL
typedef uint8_t BOOL;
#define TRUE (1)
#define FALSE (0)
#endif

// This outputs timing information    
// #define TIMING

// We need this to allow this file to be included in both C and C++ code.

#ifdef __cplusplus
#define STATIC_DECL
#else
#define STATIC_DECL static
#endif

STATIC_DECL inline UINT64 diffTVs (struct timeval * startTV, struct timeval * endTV)
{
    return (((endTV->tv_sec - startTV->tv_sec) * 1000000) + (endTV->tv_usec - startTV->tv_usec));
}

STATIC_DECL inline UINT64 tv2usec (struct timeval * TV)
{
    return (TV->tv_sec * 1000000) + TV->tv_usec;
}

#ifdef TIMING
STATIC_DECL inline void fprintTimer(FILE * out, const char * string, UINT64 timer)
{
    fprintf(out, "%s", string);
    fprintTimeMicroSeconds(out, timer, 3);
    fprintf(out, " User CPU time.\n");
}

STATIC_DECL inline void fprintTimerWithTotalPercent(FILE * out, const char * string, UINT64 timer, UINT64 fullTimer)
{
    fprintf(out, "%s", string);
    fprintTimeMicroSeconds(out, timer, 3);
    fprintf(out, " User CPU time (%.2f%%).\n", ((double)100.0 * timer / fullTimer));
}

#define defineTimer(timer) static UINT64 timer = 0

#define resetTimer(timer)                       \
    timer = 0;

#define USE_RUSAGE_FOR_TIMERS
#ifdef USE_RUSAGE_FOR_TIMERS
#define setupTimers()                                                   \
    static struct rusage __usagebuf;                                    \
    static struct timeval __firstTV = {0, 0};                           \
    static struct timeval __lastTV = {0,0};                             \
    static struct timeval __startTV = {0, 0};                           \
    static struct timeval __endTV = {0, 0};                             \
    static BOOL __started = FALSE;                                      \
    // static UINT64 timerCount = 0;                                     

#define getUTime(timebuf) timebuf.ru_utime

#define startTime()                                    \
    if (! __started)                                   \
    {                                                  \
        getrusage(RUSAGE_SELF, &__usagebuf);           \
        __firstTV = getUTime(__usagebuf);              \
        __startTV = __firstTV;                         \
    }                                                  \
    __started = TRUE;

#define endTime(timer)                                                  \
    getrusage(RUSAGE_SELF, &__usagebuf);                                \
    __lastTV = getUTime(__usagebuf);                                    \
    timer = diffTVs(&__firstTV, &__lastTV);                             \
    // fprintf(stderr, "AddToTimer called %zd times.\n", timerCount);

#define addToTimer(timer)                                               \
    getrusage(RUSAGE_SELF, &__usagebuf);                                \
    __endTV = getUTime(__usagebuf);                                     \
    timer += diffTVs(&__startTV, &__endTV);                             \
    __startTV = __endTV;                                                \
    // timerCount += 1;                                                  

#else // USE_RUSAGE_FOR_TIMERS
// We'll be using clock() to keep time.
#define setupTimers()                                                   \
    static clock_t __firstTV = 0;                                       \
    static clock_t __lastTV = 0;                                        \
    static clock_t __startTV = 0;                                       \
    static clock_t __endTV = 0;                                         \
    static BOOL __started = FALSE;                                      \
    // static UINT64 timerCount = 0;                                     

#define startTime()                                    \
    if (! __started)                                   \
    {                                                  \
        __firstTV = clock();                           \
        __startTV = __firstTV;                         \
    }                                                  \
    __started = TRUE;

#define endTime(timer)                                                  \
    __lastTV = clock();                                                 \
    timer = __lastTV -  __firstTV;                                      \
    // fprintf(stderr, "AddToTimer called %zd times.\n", timerCount);

#define addToTimer(timer)                                               \
    __endTV = clock();                                                  \
    timer += __endTV - __startTV;                                       \
    __startTV = __endTV;                                                \
    // timerCount += 1;                                                  

#endif //USE_RUSAGE_FOR_TIMERS
#else // TIMING
// We want to make this null to avoid overhead when not timing.

#define fprintTimer(a, b, c) ;
#define fprintTimerWithTotalPercent(a, b, c, d) ;

#define setupTimers() ;
#define defineTimer(timer) ;
#define startTime() ;
#define endTime(timer) ;
#define addToTimer(timer) ;
#define resetTimer(timer) ;

#endif //TIMING

#ifdef NOTNOW
static void timeTimers()
{
    struct rusage usagebuf;
    getrusage(RUSAGE_SELF, &usagebuf);
    clock_t clocktime = 0;
    struct timeval startuserTV = usagebuf.ru_utime;
    struct timeval startsysTV = usagebuf.ru_stime;
    int iterCount = 10000000;
    for (int i=0; i<iterCount; i++) clocktime += clock();
    getrusage(RUSAGE_SELF, &usagebuf);
    struct timeval enduserTV = usagebuf.ru_utime;
    struct timeval endsysTV = usagebuf.ru_stime;
    fprintf(stderr, "%d calls to clock took ", iterCount);
    fprintTimeMicroSeconds(stderr, diffTVs(&startuserTV, &enduserTV), 4);
    fprintf(stderr, " user time, and ");
    fprintTimeMicroSeconds(stderr, diffTVs(&startsysTV, &endsysTV), 4);
    fprintf(stderr, " system time, and clock sum of ");
    fprintf(stderr, "%.4lf\n", ((double)(clocktime))/CLOCKS_PER_SEC);

    getrusage(RUSAGE_SELF, &usagebuf);
    startuserTV = usagebuf.ru_utime;
    startsysTV = usagebuf.ru_stime;
    for (int i=0; i<iterCount; i++) getrusage(RUSAGE_SELF, &usagebuf);
    getrusage(RUSAGE_SELF, &usagebuf);
    enduserTV = usagebuf.ru_utime;
    endsysTV = usagebuf.ru_stime;
    fprintf(stderr, "%d calls to rusage took ", iterCount);
    fprintTimeMicroSeconds(stderr, diffTVs(&startuserTV, &enduserTV), 4);
    fprintf(stderr, " user time, and ");
    fprintTimeMicroSeconds(stderr, diffTVs(&startsysTV, &endsysTV), 4);
    fprintf(stderr, " system time.\n");
}
#endif
