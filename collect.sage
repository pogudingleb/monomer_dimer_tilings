def RemoveBigOh(s):
  return "+".join( s.split("+")[:-1] )

def Reconstruction(data):
  """ data is a set of pairs (mod, residue) """
  try:
    n = crt([mod(p[1], p[0]) for p in data]).rational_reconstruction()
  except ArithmeticError:
    n = -1
  return n

def ReconstructSeries(num_primes, mask, precision):
  data = []

  for i in range(num_primes):
    with open(mask + str(i + 1), 'r') as f:
      f.readline()
      prime = int(f.readline().split()[-1])
      GroundField = GF(prime)
      S.<p> = PowerSeriesRing(GroundField, default_prec = precision)
      T.<z> = PowerSeriesRing(GroundField, default_prec = precision)
      # Read grand potential
      f.readline()
      s = f.readline()
      grand_potential = T( RemoveBigOh(s) )
      # Read activity
      f.readline()
      s = f.readline()
      activity = S( RemoveBigOh(s) )
      # Read final
      s = f.readline()
      energy = S( RemoveBigOh(s) )
      data.append( (prime, grand_potential.padded_list()[:precision], activity.padded_list()[:precision], energy.padded_list()[:precision]) )
     
  result = []
  for i in range(1, 4):
    reconstructed_series = []
    for j in range(precision):
      residues = [(p[0], p[i][j]) for p in data]
      reconstructed_series.append( Reconstruction(residues) )
    result.append(reconstructed_series)

  return result

def CompareSeries(s_left, s_right, precision):
  return min(filter(lambda i: (i >= precision) or (s_left[i] != s_right[i]) or (s_right[i] == -1), range(precision + 1)))

def CompareResults(res_left, res_right, precision):
  return [ CompareSeries(res_left[i], res_right[i], precision) for i in range(3)]

import sys
from os import listdir
from sage.rings.finite_rings.integer_mod_ring import crt

mask = sys.argv[1]
precision = int(sys.argv[2])

files = os.listdir("./")
num_primes = 1
for filename in files:
  if filename[:len(mask)] == mask:
    if num_primes < int(filename[len(mask):]):
      num_primes = int(filename[len(mask):])

print "Num primes " + str(num_primes)

for i in range(3, num_primes):
  print CompareResults( ReconstructSeries(i, mask, precision), ReconstructSeries(i + 1, mask, precision), precision )

for r in ReconstructSeries(num_primes, mask, precision):
  print r

