#include <errno.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define PRECISION 70
#define PARALLEL_DEPTH 4

typedef struct PowerSeries { uint32_t coeffs[PRECISION]; } PowerSeries;

size_t CountBinDigits(size_t n) {
  size_t result = 0;
  while (n != 0) {
    result += (n % 2);
    n = n >> 1;
  }
  return result;
}

void Reverse(PowerSeries *data, size_t length) {
  PowerSeries tmp;
  for (size_t i = 0; i < length / 2; ++i) {
    tmp = data[i];
    data[i] = data[length - i - 1];
    data[length - i - 1] = tmp;
  }
}

uint32_t AddMod(uint32_t a, uint32_t b, uint32_t prime) {
  size_t sum = (size_t)(a) + (size_t)(b);
  return (uint32_t)( sum % prime );
}

void UpdateEntryAtPosition(PowerSeries *data, size_t ind, size_t bit, uint32_t prime) {
  size_t nbit = bit >> 1;

  if ( (bit & ind) != 0 ) {
    for (size_t j = 0; j < PRECISION; ++j) {
      data[bit ^ ind].coeffs[j] = AddMod( data[ind].coeffs[j], data[bit ^ ind].coeffs[j], prime );
    }
    if ( (nbit & ind) != 0 ) {
      for (size_t j = 1; j < PRECISION; ++j) {
        data[ind ^ bit ^ nbit].coeffs[j] = AddMod(data[ind].coeffs[j - 1], data[ind ^ bit ^ nbit].coeffs[j], prime);
      }
    }
  } 
}

void MultiplyByTransferMatrixAtPosition(PowerSeries *data, size_t height, size_t pos, uint32_t prime) {
  size_t bit = (1ull) << (height - 1 - pos);
  size_t nbit = bit >> 1;
  size_t length = (1ull) << height;

  if (height <  2 * PARALLEL_DEPTH + 3) {
    for (size_t ind = 0; ind < length; ++ind) {
      UpdateEntryAtPosition(data, ind, bit, prime);
    }
  } else {
    size_t shard = 0;
    size_t shards_number = (1 << PARALLEL_DEPTH);
    size_t local_length = (length >> PARALLEL_DEPTH); 
    if (pos < PARALLEL_DEPTH) {
      #pragma omp parallel for private(shard)
      for (shard = 0; shard < shards_number; ++shard){
        size_t tail = shard;
        for (size_t i = 0; i < local_length; ++i) {
          UpdateEntryAtPosition(data, tail + (i << PARALLEL_DEPTH), bit, prime);
        }
      }
    } else {
      #pragma omp parallel for private(shard)
      for (shard = 0; shard < shards_number; ++shard){
        size_t head = (size_t)(shard << (height - PARALLEL_DEPTH));
        for (size_t i = 0; i < local_length; ++i) {
          UpdateEntryAtPosition(data, i + head, bit, prime);
        }
      }
    }
  }
}

void MultiplyByTransferMatrix(PowerSeries *data, size_t height, uint32_t prime) {
  size_t length = (1ull) << height;
  Reverse(data, length);

  for (size_t pos = 0; pos < height; ++pos) {
    MultiplyByTransferMatrixAtPosition(data, height, pos, prime);
  }

  for (size_t i = 1; i < length; ++i) {
    size_t bd = CountBinDigits(i);
    for (size_t j = PRECISION - 1; j >= bd; --j) {
      data[i].coeffs[j] = data[i].coeffs[j - bd];
    }
    for (size_t j = 0; j < bd; ++j) {
      data[i].coeffs[j] = 0;
    }
  }

}

PowerSeries* CountTilings(size_t height, size_t min_width, size_t max_width, uint32_t prime) {
  PowerSeries *result = (PowerSeries *)calloc(max_width - min_width + 1, sizeof(PowerSeries));
  if (result == NULL) {
    printf("Cannot calloc for the result\n");
    return result;
  }

  size_t items_number = ((size_t)1) << height;
  PowerSeries *data = (PowerSeries *)calloc( items_number, sizeof(PowerSeries) );
  if (data == NULL) {
    printf("Cannot calloc, %s \n", strerror(errno));
    return result;
  }

  data[0].coeffs[0] = 1; 

  for (size_t i = 1; i <= max_width; ++i) {
    MultiplyByTransferMatrix(data, height, prime);
    if (i >= min_width) {
      result[i - min_width] = data[0];
    }
  }

  free(data);
  return result;
}

void PrintPowerSeries(PowerSeries a) {
  for (size_t i = 0; i < PRECISION; ++i) {
    printf("%ju ", a.coeffs[i]);
  }
}

int main(int argc, char *argv[]) {

  size_t height = atoi(argv[1]);
  size_t min_width = atoi(argv[2]);
  size_t max_width = atoi(argv[3]);
  size_t prime = atoi(argv[4]);

  omp_set_num_threads(1 << PARALLEL_DEPTH);
  omp_set_dynamic(0);

  PowerSeries *result = CountTilings(height, min_width, max_width, prime);
  for (size_t i = min_width; i <= max_width; ++i) {
    PrintPowerSeries(result[i - min_width]);
    printf(";");
  }
  free(result);

  return 0;
}
