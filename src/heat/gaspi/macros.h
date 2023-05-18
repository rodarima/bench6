#ifndef MACROS_H
#define MACROS_H

#include <GASPI.h>
#include <stdlib.h>
#include <stdio.h>

#define CHECK(f...)                                                       \
    {                                                                     \
        const gaspi_return_t __r = f;                                     \
        if (__r != GASPI_SUCCESS) {                                       \
            printf("Error: '%s' [%s:%i]: %i\n",#f,__FILE__,__LINE__,__r); \
            exit (EXIT_FAILURE);                                          \
        }                                                                 \
    }

#define WAIT_FOR_ENTRIES(queue, requested_entries)                   \
    {                                                                \
        gaspi_queue_id_t __queue = queue;                            \
        gaspi_number_t __requested_entries = requested_entries;      \
        gaspi_number_t __queue_size_max;                             \
        gaspi_number_t __queue_size;                                 \
                                                                     \
        CHECK(gaspi_queue_size_max(&__queue_size_max));              \
        CHECK(gaspi_queue_size(__queue, &__queue_size));             \
                                                                     \
        if (__queue_size + __requested_entries > __queue_size_max) { \
            CHECK(gaspi_wait(__queue, GASPI_BLOCK));                 \
        }                                                            \
    }

#endif // MACROS_H

