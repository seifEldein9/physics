
from math import*
#mechanics physics

#Speed calculation
def velocity(distance,time):
    return distance/time
#time calculation 
def time(distance,velocity):
    return distance/velocity
#Distance calculation
def distance(velocity,time):
    return velocity*time
#Calculating acceleratio
def acceleration(velocity,time):   
    return velocity/time
#volume
def volume(length,Height,width):   
    return length*Height*width
#density
def density(mass,volume):
    return mass/volume
#force
def force(mass,acceleration):
    return mass*acceleration
#work
def work(force,distance):
    return force*distance
#work with cos theta
def workTheta(force,distance,theta):
    return force*distance*cos(theta)
#The resultant of two forces between them is theta
def resultantOftwoForces(force1,force2,theta):
    return sqrt(pow(force1,2)+pow(force2,2)+2*force1*force2*cos(theta))
#The direction of the resultant force between two forces is theta angle
def directionOfTheResultant(force1,force2,theta):
    return atan((force2*sin(theta))/(force1+force2*cos(theta)))
#Force analysis in two orthogonal directions
def forceAnalysisInTwoOrthogonalDirections(force,theta):
    fx = force*cos(theta)
    Rx = Rx + fx
    fy = force*sin(theta)
    Ry = Ry + fy
    return sqrt(pow(Rx,2)+pow(Ry,2))
#The direction of the resultant two-force analysis
def directionOfTheResultantTwoForce(Rx,Ry):
    return atan((Ry)/(Rx))
#gravitational constant
def g():
    return 9.8
#linear distance
def linearDistance(velocity,time,acceleration):
    return (velocity*time+(1/2)*acceleration*pow(time,2))
#linear time
def linearTime(velocity,initialvelocity,acceleration):
    return (initialvelocity-velocity)/(acceleration)
#circular time
def circularTime(theta,Omega):
    return (theta)/(Omega)
#periodic displacement
def periodicDisplacement(Omega,time):
    return (Omega)*(time)
#periodic speed
def omega(theta,time):
    return (theta)/(time)
#frequency
def omega(Omega):
    return (Omega)/(2*3.141592654)
#Angular velocity
def angularVelocity(theta,time,frequency):
    return (theta)/(time)*((2*3.141592654)*frequency)
#linear wheel
def linearAcceleration(Omega,time):
    return (Omega)/(time)
#tangential velocity
def tangentialVelocity(Omega,length):
    return (Omega)*(length)
#torque 
def torque(force,length):
    return (force)*(length)
#torque with theta
def torqueWithTheta(force,length,theta):
    return (force)*(length)*sin(theta)
#Friction force
def frictionforce(u,fn):
    return (u)*(fn)
#friction coefficient
def frictionCoefficient(Ff,fn):
    return (Ff)/(fn)
#gravitational force
def gravitationalForce(Mass):
    return (Mass)*(9.80665)
#The general law of attraction
def generalLawOfAttraction(Mass1,Mass2,length):
    return (6.67428*pow(10,-11))*((Mass1*Mass2)/(pow(length,2)))
#Kinetic energy
def kineticEnergy(Mass,velocity):
    return (1/2)*(Mass)*pow(velocity,2)
#Energy situation
def energySituation(work,length):
    return (work)*(length)
#HookeIsLaw
def hookeIsLaw(k,x):
    return (-k)*(x)
#periodic time
def periodicTime(Omega):
    return (2*(3.141592654))/(Omega)
#Power
def power(work,time):
    return (work)/(time)
#Conversion from horsepower to joules
def conversionFromHorsepowerToJoules(hp):
    return (hp)*(735)
#Graphic ability
def graphicAbility(fp,Bp):
    return (fp)+(Bp)