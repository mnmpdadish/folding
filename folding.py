
from numpy import *

def main():
   test1 = Fold([3,2,0], array([[0,0,0],[1,0,0],[0,1,0],[1,1,0]]), array([[2,0,0],[0,2,0],[0,0,1]]))
   print(test1)
   

def Fold(position,sites,e):
   vol = (dot(e[0],cross(e[1],e[2])))   # equivalent to determinant, but keep the number integers
   if vol < 0:
      vol *= -1
      e[0] *= -1 
   e = transpose(e)
   found = False
   
   for siteIndex in range(0,len(sites)):
      site = sites[siteIndex]
      r = position - site
      R = CramersRuleInt(e,r)
      Q = R % vol
      if all(Q == array([0,0,0])):
         found= True
         break   
   if found:
      R2 = dot(e,R)//vol
      return R2,site
   else:
      print('Unable to fold the position',position, 'on the original cluster with this superlattice')
      assert()



def CramersRuleInt(e,r):
   #solve the equation for x in e*x = r using Cramer's rule, where r and x are 3-coordinates vectors, 
   #and e is an invertible 3x3 matrix, all input number should be integers. note that this function
   #return R and not x, where x = R / det(e)
   R = zeros((3),dtype=int)
   R[0] =  r[0]*(e[1][1]*e[2][2] - e[2][1]*e[1][2]) - r[1]*(e[0][1]*e[2][2] - e[2][1]*e[0][2]) + r[2]*(e[0][1]*e[1][2] - e[1][1]*e[0][2])
   R[1] = -r[0]*(e[1][0]*e[2][2] - e[2][0]*e[1][2]) + r[1]*(e[0][0]*e[2][2] - e[2][0]*e[0][2]) - r[2]*(e[0][0]*e[1][2] - e[1][0]*e[0][2])
   R[2] =  r[0]*(e[1][0]*e[2][1] - e[1][1]*e[2][0]) - r[1]*(e[0][0]*e[2][1] - e[0][1]*e[2][0]) + r[2]*(e[0][0]*e[1][1] - e[0][1]*e[1][0])
   return R

main()
