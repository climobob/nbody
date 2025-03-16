import sys
from math import *
import pstats
import copy

import numpy as np

from astronomy import *
 
#-----------------------------------------------------------------
class point_mass: 

  def __init__(self):
     self.x = np.zeros((3), dtype=np.float64)
     self.u = np.zeros((3), dtype=np.float64)
     self.a = np.zeros((3), dtype=np.float64)
     self.dx = np.zeros((3), dtype=np.float64) # auxiliary variable
     self.m = 0.
     self.k = 0.      # keplerian gravitational constant

  def init_loc(self, x1, x2, x3):
    self.x[0] = x1
    self.x[1] = x2
    self.x[2] = x3

  def init_vel(self, v1, v2, v3):
    self.u[0] = v1
    self.u[1] = v2
    self.u[2] = v3

  def ke(self):
    return 0.5*self.m*(self.u[0]*self.u[0] + self.u[1]*self.u[1] + self.u[2]*self.u[2])

  def show(self, l = 1., fout = sys.stdout):
    print(self.x/l, self.u, self.m, file = fout)
     
  #################################################################################
  # Most of time is spent in update_loc, update_vel, gravity
  # Gravity ~= 3*(_loc + _vel) time
  # explicit listing of elements is about 2x faster than using vector notation
  def update_loc(self, dt = 1.):
    #slower: self.x  += dt*self.u
    self.x[0]  += dt*self.u[0]
    self.x[1]  += dt*self.u[1]
    self.x[2]  += dt*self.u[2]

  def update_vel(self, dt = 1.):
    self.u[0]  += dt*self.a[0]
    self.u[1]  += dt*self.a[1]
    self.u[2]  += dt*self.a[2]

    self.a[0]  = 0.
    self.a[1]  = 0.
    self.a[2]  = 0.

  def update(self, dt = 1.):
    #Euler symplectic -- https://www.mgaillard.fr/2021/07/11/euler-integration.html
    self.u[0]  += dt*self.a[0]
    self.u[1]  += dt*self.a[1]
    self.u[2]  += dt*self.a[2]

    self.a[0]  = 0.
    self.a[1]  = 0.
    self.a[2]  = 0.

    self.x[0]  += dt*self.u[0] 
    self.x[1]  += dt*self.u[1]
    self.x[2]  += dt*self.u[2]

  def gravity(self, y):
    #maybe 7% slower: self.dx = y.x - self.x
    self.dx[0] = y.x[0] - self.x[0]
    self.dx[1] = y.x[1] - self.x[1]
    self.dx[2] = y.x[2] - self.x[2]

    # These two are equivalent in time
    r = sqrt(self.dx[0]*self.dx[0] + self.dx[1]*self.dx[1] + self.dx[2]*self.dx[2])
    a0 = y.k / r / r / r
    #e2 r = pow(self.dx[0]*self.dx[0] + self.dx[1]*self.dx[1] + 
    #           self.dx[2]*self.dx[2], 1.5)
    #e2 a0 = y.k / r

    #adds ~60% to run time: self.a += a0*self.dx
    self.a[0] += a0*self.dx[0]
    self.a[1] += a0*self.dx[1]
    self.a[2] += a0*self.dx[2]
    
    # This demands that caller realize that action = -reaction and not call both i,j and j,i
    a0 *= (self.k / y.k)
    y.a[0] -= a0*self.dx[0]
    y.a[1] -= a0*self.dx[1]
    y.a[2] -= a0*self.dx[2]
    

#--------------------------------------------------------------------
# class real_mass: -- bodies with size + shape --> moments of inertia
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# something to work on collections of point masses
class point_system: 
  com = point_mass()
  tke = 0.0
  tpe = 0.0

  def __init__(self, nbody = 3):
     self.nbody = nbody
     self.body  = []
     for i in range(0, nbody):
       tmp = point_mass()
       self.body.append(tmp)
       del tmp

  # PE defined by sum over all other bodies, 
  #    don't do both i,j and j,i -- double counts
  def PE(self):
    sum = 0.0
    tmp = 0.0
    for i in range(0, self.nbody):
      #debug: self.body[i].show()
      for j in range(i+1, self.nbody):
        #debug: self.body[j].show()
        #debug: print( self.body[i].m, self.body[j].k, dist(self.body[i], self.body[j]), flush=True)
        tmp = self.body[i].m * self.body[j].k / dist(self.body[i], self.body[j])
        sum -= tmp
    return sum

  def ke(self):
    sum = 0.0
    for i in range(0, self.nbody):
      sum += 0.5*self.body[i].m*(self.body[i].u[0]*self.body[i].u[0] + 
                                 self.body[i].u[1]*self.body[i].u[1] + 
                                 self.body[i].u[2]*self.body[i].u[2]   )
    return sum

  def center_of_mass(self):
    com = point_mass()
    self.com = copy.deepcopy(com)
    for i in range (0, self.nbody):
      self.com.m += self.body[i].m
      self.com.x += self.body[i].x*self.body[i].m
      self.com.u += self.body[i].u*self.body[i].m
    self.com.x /= self.com.m  
    self.com.u /= self.com.m
    del com
    return self.com

  def comoving(self):
    self.tke = point_system.ke(self)
    self.tpe = point_system.PE(self)
    # adjust so as to reference comoving with center of mass:
    for i in range(0, self.nbody):
      self.body[i].x -= self.com.x
      self.body[i].u -= self.com.u
    self.com.x -= self.com.x
    self.com.u -= self.com.u
    self.tke = point_system.ke(self)
    self.tpe = point_system.PE(self)

  def conservation(self):
    self.center_of_mass()
    self.tke = point_system.ke(self)
    self.tpe = point_system.PE(self)
    
  def com2(first, second):
    com = point_mass()
    com.m = first.m + second.m
    if (com.m > 0):
      com.x = first.m*first.x + second.m*second.x 
      com.u = first.m*first.u + second.m*second.u 
      com.x /= com.m
      com.u /= com.m
    com.k = com.m*astronomy.G
    return com

  def collision(self, tolerance = 0.01*3.84399e8, fout = sys.stdout ):
     #RG: Note that current version is preserving momentum, but not PE, KE
     newcount = self.nbody
     for i in range (0, self.nbody-1):
       for j in range (i+1, self.nbody):
         if (dist(self.body[i], self.body[j]) < tolerance):
           # merge masses -- center of mass for the two
           #debug:
           print("collision between bodies ",i,j, "masses ",self.body[i].m, self.body[j].m, flush = True, file = fout)
           self.body[i].show(fout = fout)
           self.body[j].show(fout = fout)
           tmp = point_system.com2(self.body[i], self.body[j])
           self.body[i] = copy.deepcopy(tmp)
           self.body[j].m = 0.
           newcount -= 1
           #debug: 
           print("after collision: ", self.body[i].m, flush = True, file = fout)

     # clean up step decrementing nbody and removing massless bodies
     for i in range (0, self.nbody-1):
       if (self.body[i].m == 0.):
         #debug: 
         print("body ",i," is massless, swap in a massy body", file = fout, flush=True)
         for j in range (self.nbody-1, i, -1):
           if (self.body[j].m != 0.):
             #debug:
             print(" swapping body ",j,"in to massless body ",i, flush=True, file = fout)
             self.body[i]   = copy.deepcopy(self.body[j])
             self.body[j].m = 0.
             break
     #debug: 
     print("change in body count: ",self.nbody, newcount, flush=True, file = fout)
     self.nbody = newcount

     # update system tpe, tke:
     #debug: print("tke tpe e before collision:",self.tke, self.tpe, self.tpe+self.tke, flush=True, file = fout)
     tmpke = self.tke
     tmppe = self.tpe
     self.tke = point_system.ke(self)
     self.tpe = point_system.PE(self)
     print("delta tke tpe e after collision: ",tmpke-self.tke, tmppe-self.tpe, 
                                  tmpke+tmppe-self.tpe-self.tke, flush=True, file = fout)
     # center of mass before and after:
     #debug: exit(1)
  
  def closest(self):
    nearest = 10.*astronomy.rmoon
    for i in range(0, self.nbody):
      for j in range(i+1, self.nbody):
        #about 50% faster to use tmp
        tmp = dist(self.body[i], self.body[j])
        if (tmp < nearest):
          nearest = tmp
        #if (dist(self.body[i], self.body[j]) < nearest):
        #  nearest = dist(self.body[i], self.body[j])
    return nearest

  def step(self, dt):
    #RG: is ~2x more efficient given action == -reaction
    for k in range(0, self.nbody-1):
      for l in range(k+1, self.nbody):
        self.body[k].gravity(self.body[l])
  
    for j in range (0, self.nbody):
      #self.body[j].update_loc(dt)
      #self.body[j].update_vel(dt)
      self.body[j].update(dt)

  def show(self, fout = sys.stdout ):
    print("\n system, nbodies = ",self.nbody, file = fout)
    for i in range (0, self.nbody):
      self.body[i].show(fout = fout) 

#--------------------------------------------------------------------
