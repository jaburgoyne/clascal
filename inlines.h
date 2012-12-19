// CLASCAL (inlines.h)
//
// Copyright (c) 2010 by John Ashley Burgoyne and the Royal Institute for the 
// Advancement of Learning (McGill University). All rights reserved.
//
// This source is adapted from Suzanne Winsberg's CLASCAL, version 7.01 (May
// 1993), written in FORTRAN 77. 
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
//   1. Redistributions of source code must retain the above copyright notice, 
//      this list of conditions, and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions, and the following disclaimer in the 
//      documentation and/or other materials provided with the distribution.
//
//   3. Neither the name of McGill University nor the names of its contributors
//      may be used to endorse or promote products derived from this software 
//      without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.

// TODO: Find any dangling pointers not initialised to NULL with ".*\*[^=]*;".

// TODO: Find '+ 1's floating around and convert to SizeSum.

// NOTE: Any arguments passed as non-const in constructors become the 
//       responsibility of the constructed object with respect to memory
//       management. The flip side of this agreement is that any function
//       calling a constructor must ensure that the const arguments to that
//       constructor do not disappear during execution.

/*
 * Requires <limits.h>, <math.h>, <stdbool.h>, <stdint.h>, <stdio.h>, and 
 * <stdlib.h>.
 */

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

/* Check whether a size_t parameter can be converted safely to int. */
static inline bool IsConvertibleToInt(size_t size)
{
        return ((uintmax_t)size <= (uintmax_t)INT_MAX);
}

/* Standard error printing function. */
#define ExitWithError(s) _ExitWithError(s, __FILE__, __LINE__, __func__)
static inline void _ExitWithError(const char * error,
                                  const char * fileName,
                                  int lineNumber,
                                  const char * functionName)
{
#ifdef VERBOSE_EXIT
        fprintf(stderr,
                "ERROR: %s in line %d of %s (%s).\n",
                error,
                lineNumber,
                fileName,
                functionName);
#else
        fprintf(stderr,
                "ERROR: %s.\n",
                error);
#endif
        exit(EXIT_FAILURE);

}

/* 
 * Add two size_t variables safely. An assertion will fail if the result would
 * wrap. It is unnecessary to validate the size of passed arrays: this must be 
 * ensured by the caller. For the same reason, it is unnecessary to call this
 * macro before indexing an array.
 */
#define SizeSum(x, y) _SizeSum(x, y, __FILE__, __LINE__, __func__)
static inline size_t _SizeSum(size_t x, 
                              size_t y, 
                              const char * fileName,
                              int lineNumber,
                              const char * functionName)
{
        if (SIZE_MAX - x < y)
                _ExitWithError("Sum exceeds SIZE_MAX",
                               fileName,
                               lineNumber,
                               functionName);
        return x + y;
}

/* 
 * Multiply two size_t variables safely. An assertion will fail if the result
 * would wrap. It is unnecessary to validate the size of passed arrays: this 
 * must be ensured by the caller. For the same reason, it is unnecessary to call 
 * this macro before indexing an array.
 */
#define SizeProduct(x, y) _SizeProduct(x, y, __FILE__, __LINE__, __func__)
static inline size_t _SizeProduct(size_t x, 
                                  size_t y,
                                  const char * fileName,
                                  int lineNumber,
                                  const char * functionName)
{
        if (x > 0 && SIZE_MAX / x < y) 
                _ExitWithError("Product exceeds SIZE_MAX",
                               fileName,
                               lineNumber,
                               functionName);
        return x * y;
}

/* Convert a string to size_t safely (returns zero if oversize). */
static inline size_t SizeRead(const char * restrict s)
{
        unsigned long long int i = strtoull(s, NULL, 10);
        return ((uintmax_t)i <= (uintmax_t)SIZE_MAX) ? (size_t)i : 0;
}


/*
 * Allocate memory with bounds checking. Exit with error if allocation fails.
 */
#define SafeMalloc(x, y) _SafeMalloc(x, y, __FILE__, __LINE__, __func__)
static inline void * _SafeMalloc(size_t nmemb,
                                 size_t size,
                                 const char * fileName,
                                 int lineNumber,
                                 const char * functionName)
{
        if (nmemb > 0 && SIZE_MAX / nmemb < size) 
                _ExitWithError("Requested allocation exceeds SIZE_MAX",
                               fileName, 
                               lineNumber, 
                               functionName);
        void * pointer = malloc(nmemb * size);
        if (!pointer)
                _ExitWithError("Unable to allocate memory as requested",
                               fileName,
                               lineNumber,
                               functionName);
        return pointer;
}

/*
 * Allocate and zero memory with bounds checking. Exit with error if allocation
 * fails.
 */
#define SafeCalloc(x, y) _SafeCalloc(x, y, __FILE__, __LINE__, __func__)
static inline void * _SafeCalloc(size_t nmemb, 
                                 size_t size,
                                 const char * fileName,
                                 int lineNumber,
                                 const char * functionName)
{
        if (SIZE_MAX / nmemb < size) 
                _ExitWithError("Requested allocation exceeds SIZE_MAX",
                               fileName, 
                               lineNumber, 
                               functionName);
        void * pointer = calloc(nmemb, size);
        if (!pointer)
                _ExitWithError("Unable to allocate memory as requested",
                               fileName,
                               lineNumber,
                               functionName);
        return pointer;
}

/*
 * Frees a pointer and sets it to NULL.
 */
#define FreeAndClear(x) free(x); (x) = NULL;
