#! /home/staff/pogudin/SageMath/sage

@parallel(1)
def _ComputeTilingsNumber(n, height, widths):
  min_width = min(widths)
  max_width = max(widths)
  command = './tilings ' + str(height) + ' ' + str(min_width) + ' ' + str(max_width) + ' ' + str(PRIME)
  raw_result = subprocess.Popen(command, stdout = subprocess.PIPE, shell = True).communicate()
  result = map(
    lambda s: ActivityPowerSeries(s.strip().split(' ')).log(), 
    raw_result[0].strip().split(';')[:-1]
  )

  return result

def _ComputeTriangularPolynomials(n):
  result = [ [0] * (i + 1) for i in range(n + 1) ]
  max_height = (n + 2) // 2

  inputs = []
  for h in reversed(range(1, max_height + 1)):
    widths = filter(lambda x: x > 0, range(n - 1 - h, n + 3 - h))
    inputs.append((n, h, widths))

  outputs = _ComputeTilingsNumber(inputs)
  for out in outputs:
    h = out[0][0][1]
    widths = out[0][0][2]
    polynomials = out[1]
    for i in range(len(widths)):
      result[-h][-(i + 1)] = polynomials[-(i + 1)]

  for i in range(n + 1 - max_height):
    for j in range(i + 1):
      result[i][j] = result[-(j + 1)][n - i]  

  return result


def _ComputeTriangularCoefficients(n):
  if n == 0:
    return [[1]]
  else:
    result = _ComputeTriangularCoefficients(n - 1)
    result.append([0] + result[-1])
    w = 0
    for i in range(n + 1):
      for j in range(i + 1):
        w = w + (n + 1 - i) * (j + 1) * result[i][j]
    result[-1][0] = 1 - w
    return result


def ComputeGrandPotential(n, d):
  """
  Computes \Gamma(z) in terms of Gaunt paper, eq. (2.7)
  """
  coefficients = _ComputeTriangularCoefficients(n)
  logs = _ComputeTriangularPolynomials(n)
  result = 0
  for i in range(n + 1):
    for j in range(i + 1):
      result = result + coefficients[i][j] * logs[i][j]
  return GroundField(1/d) * result


def ComputeGrandPotentialNew(n, d):
  """
  Computation of the grand potential using corrections
  """
  n0 = n - 5
  result = ComputeGrandPotential(n0, d)
  result = result \
    + z^(n0 + 1) * (-1)^n0 * GroundField(1/6) * (27 * 2^(n0 + 2) - 5 * n0^3 - 21 * n0^2 - 70 * n0 - 102) \
    + z^(n0 + 2) * (-1)^(n0 + 1) * GroundField(1/6) * (81 * n0 * 2^(n0 + 3) + 9 * 2^(n0 + 5) - 49 * n0^4 - 115 * n0^3 - 521 * n0^2 - 815 * n0 - 258) \
    + z^(n0 + 3) * (-1)^n0 * GroundField(1 / 60) * ( 1215 * n0^2 * 2^(n0 + 4) + 5805 * n0 * 2^(n0 + 3) + 1545 * 2^(n0 + 3) - 3 * n0^6 - 2203 * n0^5 - 5240 * n0^4 - 21445 * n0^3 - 62347 * n0^2 - 44482 * n0 - 11820 )\
    + z^(n0 + 4) * (-1)^(n0 + 1) * GroundField(1 / 10080) * (-23758560 + 82845 * 2^(6 + n0) - 18537792 * n0 + 894915 * 2^(6 + n0) * n0 - 88211380 * n0^2 + 297675 * 2^(7 + n0) * n0^2 - 55584984 * n0^3 + 25515 * 2^(8 + n0) * n0^3 - 9604861 * n0^4 - 4256028 * n0^5 - 1068550 * n0^6 - 3276 * n0^7 - 9 * n0^8) \
    + z^(n0 + 5) * (-1)^n0 * GroundField(1/100800) * (-5423241600 - 6792975 * 2^(6+n0) + 1644194640 * n0 + 128949975 * 2^(5+n0) * n0 - 4398077956 * n0^2 + 58188375 * 2^(6+n0) * n0^2 - 6793265520 * n0^3 + 4124925 * 2^(8+n0) * n0^3 - 1624200675 * n0^4 + 382725 * 2^(8+n0) * n0^4 - 317981615 * n0^5 - 158415943 * n0^6 - 22867730 * n0^7 - 109425 * n0^8 - 575 * n0^9 - n0^10)

  return result
  

def ComputeDimerActivity(n, d):
  """
  Computes dimer activity as a power series (truncated at n + 1) in density
  """
  grand = ComputeGrandPotentialNew(n, d)
  print "GrandPotential: "
  print grand
  density = z * grand.derivative(z)
  return density.reverse().subs(z = p *  GroundField(1/(2 * d)) )

def ComputeFreeEnergy(n):
  """
  Computation of f_{\infty, \infty} via equation (19) from Kong's paper without -0.5p ln(p) + (0.5 + ln 2)p term
  """
  activity = ComputeDimerActivity(n, 2)
  print "Activity: "
  print activity
  return GroundField(-1/2) * (activity / (p * GroundField(1 / 4))).log().integral()


####################################################################

import sys
import time
import subprocess

print "Started at " + time.strftime("%a, %d %b %Y %H:%M:%S +0000")
start = time.time()

PRIME = sys.argv[1]
def_prec = 70
GroundField = GF(PRIME)
ActivityPowerSeries.<z> = PowerSeriesRing(GroundField, default_prec = def_prec)
DensityPowerSeries.<p> = PowerSeriesRing(GroundField, default_prec = def_prec)

print "Computing modulo {}".format(PRIME)

n = sys.argv[2]
f = ComputeFreeEnergy(int(n))
print f

print "Finished at " + time.strftime("%a, %d %b %Y %H:%M:%S +0000")
end = time.time()
print "It took " + str(end - start) + " seconds"
