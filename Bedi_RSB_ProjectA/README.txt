{\rtf1\ansi\ansicpg1252\cocoartf1404\cocoasubrtf470
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 Project A - README File \
\
To replicate the results of the Single Pendulum system, one must first instantiate the\
single pendulum with parameter values that indicate; The scaled damping coefficient (D),\
the initial angle (theta_0), the minimum and max time (tMin, tMax) and finally the step size (h).\
\
This is done as follows:\
instanceName = SinglePendulum(D, theta_0, tMin, tMax, h)\
\
the following example demonstrates a single pendulum with:\
**D=0, theta0 = 0.01, tMin = 1, tMax = 20, h=0.1**\
\
pend = SinglePendulum(0.2, 0.01, 0, 100, 0.05)\
\
***\
\
To solve the Single pendulum System with a FDM we simply call the method of our FDM by using the \
. operator, i.e. \'91pend.RK4()\'92. This will solve the single pendulum system using the Runge-Kutta method of global error O(h^4).\
\
pend.RK4() \'97> Runge Kutta (4th order)\
pend.explicitEuler() \'97> Explicit Euler\
pend.leapfrog() \'97> Leapfrog\
pend.implicitEuler() \'97> Implicit Euler \
_________________________________________________________________________\
\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 To replicate the results of the Double Pendulum system, one must first instantiate the\
single pendulum with parameter values that indicate; The step size (h), the Mass ratio (R), the damping coefficient (G), and the finalTime (maxTime).\
\
This is done as follows:\
instanceName = DoublePendulum(h, R, G, maxTime)\
\
the following example demonstrates a double pendulum with:\
**h =0.1, R = 100 ,G =1, tMax = 100*\
\
pendInstance = DoublePendulum(0.1, 100, 1, 100)\
\
To solve the double pendulum system I have used the explicit RK4 method, this is accessed via a method using the dot operator\'85 i.e.\
\
pendInstance.DoubleRK4()\
}