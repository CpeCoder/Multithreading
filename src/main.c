// MIT License
// 
// Copyright (c) 2023 Trevor Bakker 
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#include <sys/time.h> //included for getting time
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include <pthread.h>  //allows to use pthread and mutex functions
#include "utility.h"
#include "star.h"
#include "float.h"

#define NUM_STARS 30000 
#define MAX_LINE 1024
#define DELIMITER " \t\n"

struct Star star_array[ NUM_STARS ];
uint8_t   (*distance_calculated)[NUM_STARS];

pthread_mutex_t lock;       // for locking gobals in the function
double  min_of_t = FLT_MAX; // holds min of all threads
double  max_of_t = FLT_MIN; // holds max of all threads
double mean_sum=0;          // sum of average distance from threads parameter
uint64_t  count_sum=0;      // sum of count from function
int num_threads = 1;        // sets threads user want, else default to 1 thread

void showHelp()
{
  printf("Use: findAngular [options]\n");
  printf("Where options are:\n");
  printf("-t          Number of threads to use\n");
  printf("-h          Show this help\n");
}

typedef struct start_end_of_threads
// assigns start and end to each thread, dividing workload
{
    int start;
    int end;
} starts_ends;


// 
// Embarassingly inefficient, intentionally bad method
// to calculate all entries one another to determine the
// average angular separation between any two stars 
void *determineAverageAngularDistance( void *arg )
{
    starts_ends *thread_info = (starts_ends *)arg; // holds start and end which is passed in
    uint32_t i, j;
    uint64_t count = 0;
    double  min  = FLT_MAX;
    double  max  = FLT_MIN;
    double mean = 0;
    for (i = thread_info->start; i < thread_info->end; i++)
    {
      for (j = i; j < NUM_STARS; j++)
      {
        if( i!=j && distance_calculated[i][j] == 0 )
        {
          double distance = calculateAngularDistance( star_array[i].RightAscension, star_array[i].Declination,
                                                      star_array[j].RightAscension, star_array[j].Declination ) ;
          distance_calculated[i][j] = 1;
          distance_calculated[j][i] = 1;
          count++;

          if( min > distance )   // compare and replace
          {
            min = distance;
          }

          if( max < distance )   // compare and replace
          {
            max = distance;
          }
          mean += distance;     // getting sum from all threads so in end it can be
        }                       // combine and calculate by: total_mean / total_distance
      }
    }
    
    pthread_mutex_lock(&lock);
    count_sum += count;
    mean_sum += mean;
    if(min < min_of_t)
    {
      min_of_t = min;
    }
    if(max_of_t < max)
    {
      max_of_t = max;
    }
    pthread_mutex_unlock(&lock);
}


int main( int argc, char * argv[] )
{
  //clock_t start_c, end_c;
  //double execution_time;
  //start_c = clock();
  struct timeval begin;
  struct timeval end;
  FILE *fp;
  uint32_t star_count = 0;

  uint32_t n;

  distance_calculated = malloc(sizeof(uint8_t[NUM_STARS][NUM_STARS]));

  if( distance_calculated == NULL )
  {
    uint64_t num_stars = NUM_STARS;
    uint64_t size = num_stars * num_stars * sizeof(uint8_t);
    printf("Could not allocate %ld bytes\n", size);
    exit( EXIT_FAILURE );
  }

  int i, j;
  // default every thing to 0 so we calculated the distance.
  // This is really inefficient and should be replace by a memset
  for (i = 0; i < NUM_STARS; i++)
  {
    for (j = 0; j < NUM_STARS; j++)
    {
      distance_calculated[i][j] = 0;
    }
  }

  for( n = 1; n < argc; n++ )          
  {
    if( strcmp(argv[n], "-help" ) == 0 )
    {
      showHelp();
      exit(0);
    }
    else if( strcmp(argv[n], "-t" ) == 0 )  // taking in number of thread wanted
    {
      char *x = argv[n+1];
      num_threads = atoi(x);
    }
  }
   
  fp = fopen( "data/tycho-trimmed.csv", "r" );

  if( fp == NULL )
  {
    printf("ERROR: Unable to open the file data/tycho-trimmed.csv\n");
    exit(1);
  }

  char line[MAX_LINE];
  while (fgets(line, 1024, fp))
  {
    uint32_t column = 0;

    char* tok;
    for (tok = strtok(line, " ");
            tok && *tok;
            tok = strtok(NULL, " "))
    {
       switch( column )
       {
          case 0:
              star_array[star_count].ID = atoi(tok);
              break;
       
          case 1:
              star_array[star_count].RightAscension = atof(tok);
              break;
       
          case 2:
              star_array[star_count].Declination = atof(tok);
              break;

          default: 
             printf("ERROR: line %d had more than 3 columns\n", star_count );
             exit(1);
             break;
       }
       column++;
    }
    star_count++;
  }
  printf("%d records read\n", star_count );
  starts_ends thread_info[num_threads];
  pthread_t threads[num_threads];
  int start = 0;
  int each_thread_part = NUM_STARS/num_threads;   // work per thread
  pthread_mutex_init(&lock, NULL);
  gettimeofday( &begin, NULL );
  // Find the average angular distance in the most inefficient way possible
  for(int i=0; i < num_threads; i++)
  {
    thread_info[i].start = start * each_thread_part;       //dividing indexes
    thread_info[i].end = (start+1) * each_thread_part;
    // passing the struct with start and end values to thread individually
    pthread_create(&threads[i], NULL, determineAverageAngularDistance, (void *)&thread_info[i]);
    start++;
  }
  for (int i = 0; i < num_threads; i++)  // join threads
  {
    pthread_join(threads[i], NULL);
  }
  gettimeofday( &end, NULL );
  printf("Average distance found is %lf\n", mean_sum/count_sum ); // getting average distance
  printf("Minimum distance found is %lf\n", min_of_t );
  printf("Maximum distance found is %lf\n", max_of_t );
  
  // calculating time taken which is given by professor
  double time_duration = ( (double)( end.tv_sec - begin.tv_sec ) * 1000000 + 
                         (end.tv_usec - begin.tv_usec)) / 1000000;

  printf("Duration in seconds: %f\n", time_duration );
  return 0;
}
