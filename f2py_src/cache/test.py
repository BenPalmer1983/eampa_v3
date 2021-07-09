
from f_cache import cache
import numpy
import time



def g(x):
  y = 0.0
  for i in range(100):
    y = (i/100000) * numpy.exp(0.0005 * x * (i/10000)) + numpy.log(0.007 + i/10000 + abs(x))
  return y


def f(f_key, x, use_cache=False):
  if(use_cache):
    in_cache, y = cache.get(f_key, x)
    if(not in_cache):
      y = g(x)
      cache.set(f_key, x, y)
  else:
    y = g(x)
  return y

def test_a():
  f_key = 0
  y = f(f_key, 0.8951, True)


def test_b():
  xs = 10000
  x = numpy.linspace(-9000.0,19050.0,xs)
  y = numpy.zeros((xs,),)
  numpy.random.shuffle(x)

  l = 4
  for n in range(l):
    s = time.time()
    f_key = n % 2 + 1
    for m in range(len(x)):
      y[m] = f(f_key, x[m], True)
    print(sum(y), time.time() - s)


test_b()







