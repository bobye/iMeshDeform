#ifndef _UTILITY_H_
#define _UTILITY_H_

// Timing, count in seconds.
static struct timespec start, end;
#define BILLION  1000000000L
inline  void clock_start(std::string description) {
  std::cout << description <<" ... " << std::flush; 
  clock_gettime(CLOCK_MONOTONIC, &start);
}
inline void clock_end() {
  clock_gettime(CLOCK_MONOTONIC, &end);
  printf("\t[done] %.3f seconds\n", (( end.tv_sec - start.tv_sec )+ (double)( end.tv_nsec - start.tv_nsec ) / (double)BILLION ));
}

// Timing, count in nano seconds.
#ifdef _SS_SHOW_DEBUG
static struct timespec nstart, nend;
inline void nclock_start() {clock_gettime(CLOCK_MONOTONIC, &nstart);}
inline void nclock_end() {clock_gettime(CLOCK_MONOTONIC, &nend);     
  printf("\t%ld nanosec", ((long) ( nend.tv_sec - nstart.tv_sec ) * BILLION + ( nend.tv_nsec - nstart.tv_nsec ) ));}
#endif 

#ifdef _SS_SHOW_DEBUG
#define _SS_PROFILE(x) nclock_start(); x nclock_end();
#else
#define _SS_PROFILE(x) x
#endif



#endif /* _UTILITY_H_ */
