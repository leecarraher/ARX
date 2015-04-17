'''Copyright 2010 Lee Carraher. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Lee Carraher ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Lee Carraher OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of Lee Carraher.
'''
'''
Created on Jul 15, 2009

@author: lee
'''
from numpy import linalg
from numpy.ma.core import zeros


#gaussian elimination based solver
def solver(m, sol, eps=1.0 / (10 ** 10)):
      (h, w) = (len(m), len(m[0]))
      for y in range(0, h):
            m[y].append(sol[y])
            
      w = w + 1
      for y in range(0, h):
            maxrow = y
            for y2 in range(y + 1, h):# find pivot
                  if abs(m[y2][y]) > abs(m[maxrow][y]):
                        maxrow = y2
                        
            temp = m[y]
            m[y] = m[maxrow]
            m[maxrow]=temp

            if abs(m[y][y]) <= eps:#check singularity
                  return []
            for y2 in range(y + 1, h):#eliminate column y
                  c = m[y2][y] / m[y][y]
                  for x in range(y, w):
                        m[y2][x] -= m[y][x] * c
      for y in range(h - 1, 0 - 1, -1):#back substitution
            c = m[y][y]
            for y2 in range(0, y):
                  for x in range(w - 1, y - 1, -1):
                        m[y2][x] -= m[y][x] * m[y2][y] / c
            m[y][y] /= c
            for x in range(h, w):# normalize row y
                  m[y][x] /= c
                  
      rtr = []
      for t in range(0,h):
            rtr.append(m[t][w-1])
            
      return rtr



def ARLeastSquares(inputseries, degree):

      k=0
      length = len(inputseries)
      mat = zeros((degree,degree))
      coefficients = zeros(degree)
      
      i=degree-1
      while i<length-1:
            hi = i+1
            j=0
            while j<degree:
                  hj = i-j
                  coefficients[j] += inputseries[hi] * inputseries[hj]
                  k=j
                  
                  while k <degree:
                        mat[j][k] +=inputseries[hj] * inputseries[i-k]
                        k+=1
                  j+=1
            i+=1      
      for i in range(0,degree):
            coefficients[i] /= length - degree
            for j in range(i, degree):
                  mat[i][j] /=(length - degree)
                  mat[j][i] = mat[i][j]
                  
      M= linalg.lstsq(mat,coefficients)[0]
      return M
      
      
def testAR(series, coefs):
      ests=[]
      err = 0.0;
      for i in range(0,len(series)):
            est = 0;
            if i>len(coefs):
                  for j in range(0,len(coefs)):
                        est+=coefs[j]*series[i-j-1]
                  err += (est-series[i])**2
            ests.append(est)
      return ests,err
      

def readFile(S):
      import csv
      reader = csv.reader(open(S), delimiter=',')
      X = []
      for row in reader:
            for col in row:
                  X.append(float(col))
            
      return X

if __name__ == '__main__':
      import sys
      if len(sys.argv)<3:
            print "usage: python autoRegressive.py filename degree [--no]"
            print "\t filename - a line delimitted data file\n\t degree - the degree or histor of the yule-walker matrix\n\t --no - means no output graph, or no pylab extension"
            sys.exit(1)
      filename = sys.argv[1]
      inputseries = readFile(filename)
      mean = 0.0
      for t in range(0,len(inputseries)):
            mean+=inputseries[t]
      mean/=len(inputseries)
      
      for t in range(0,len(inputseries)):
            inputseries[t]= inputseries[t] - mean
            
      C = ARLeastSquares(inputseries,int(sys.argv[2]))
      V,err = testAR(inputseries,C)
      print err
      if len(sys.argv)<4 or sys.argv[3]!="--no":
            import pylab
            pylab.plot(range(0,len(inputseries)),inputseries)
            pylab.plot(range(0,len(inputseries)),V,linestyle='--')
            pylab.show()
         
