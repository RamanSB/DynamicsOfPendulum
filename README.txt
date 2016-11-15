Project A - README File 

To replicate the results of the Single Pendulum system, one must first instantiate the
single pendulum with parameter values that indicate; The scaled damping coefficient (D),
the initial angle (theta_0), the minimum and max time (tMin, tMax) and finally the step size (h).

This is done as follows:
instanceName = SinglePendulum(D, theta_0, tMin, tMax, h)

the following example demonstrates a single pendulum with:
**D=0, theta0 = 0.01, tMin = 1, tMax = 20, h=0.1**

pend = SinglePendulum(0.2, 0.01, 0, 100, 0.05)

***

To solve the Single pendulum System with a FDM we simply call the method of our FDM by using the 
. operator, i.e. pend.RK4() This will solve the single pendulum system using the Runge-Kutta method of global error O(h^4).

pend.RK4() ---> Runge Kutta (4th order)
pend.explicitEuler() ---> Explicit Euler
pend.leapfrog() ---> Leapfrog
pend.implicitEuler() ---> Implicit Euler 
_________________________________________________________________________


 To replicate the results of the Double Pendulum system, one must first instantiate the
single pendulum with parameter values that indicate; The step size (h), the Mass ratio (R), the damping coefficient (G), and the finalTime (maxTime).

This is done as follows:
instanceName = DoublePendulum(h, R, G, maxTime)

the following example demonstrates a double pendulum with:
**h =0.1, R = 100 ,G =1, tMax = 100*

pendInstance = DoublePendulum(0.1, 100, 1, 100)

To solve the double pendulum system I have used the explicit RK4 method, this is accessed via a method using the dot operator i.e.
pendInstance.DoubleRK4()
