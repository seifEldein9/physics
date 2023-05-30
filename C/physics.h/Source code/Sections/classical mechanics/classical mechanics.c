#include<math.h>

double velocity(double distance,double time);
double time(double distance,double velocity);
double distance(double velocity,double time);
double acceleration(double velocity,double time);
double volume(double length,double Height,double width);
double density(double mass,double volume);
double force(double mass,double acceleration);
double work(double force,double distance);
double workTheta(double force,double distance,double theta);
double resultantOftwoForces(double force1,double force2,double theta);
double directionOfTheResultant(double force1,double force2,double theta);
double forceAnalysisInTwoOrthogonalDirections(double force,double theta);
double directionOfTheResultantTwoForce(double Rx,double Ry);
double g(void);
double linearDistance(double velocity,double time,double acceleration);
double linearTime(double velocity,double initialvelocity,double acceleration);
double circularTime(double theta,double Omega);
double periodicDisplacement(double Omega,double time);
double omega(double theta,double time);
double frequency(double Omega);
double angularVelocity(double theta,double time,double frequency);
double linearAcceleration(double Omega,double time);
double tangentialVelocity(double Omega,double length);
double torque(double force,double length);
double torqueWithTheta(double force,double length,double theta);
double frictionforce(double u,double fn);
double frictionCoefficient(double Ff,double fn);
double gravitationalForce(double Mass);
double generalLawOfAttraction(double Mass1,double Mass2,double length);
double kineticEnergy(double Mass,double velocity);
double energySituation(double work,double length);
double hookeIsLaw(double k,double x);
double periodicTime(double Omega);
double power(double work,double time);
double conversionFromHorsepowerToJoules(double hp);
double graphicAbility(double fp,double Bp);



int main(){

velocity(0/*distance*/,0/*time*/);
time(0/*distance*/,0/*velocity*/);
distance(0/*velocity*/,0/*time*/);
acceleration(0/*velocity*/,0/*time*/);
volume(0/*length*/,0/*height*/,0/*Width*/);
density(0/*mass*/,0/*volume*/);
force(0/*mass*/,0/*acceleration*/);
work(0/*force*/,0/*distance*/);
workTheta(0/*force*/,0/*distance*/,0/*theta*/);
resultantOftwoForces(0/*force1*/,0/*force2*/,0/*theta*/);
directionOfTheResultant(0/*force1*/,0/*force2*/,0/*theta*/);
forceAnalysisInTwoOrthogonalDirections(0/*force*/,0/*theta*/);
directionOfTheResultantTwoForce(0/*Rx*/,0/*Ry*/);
g(); //the wheel of gravity
linearTime(0/*velocity*/,0/*initialvelocity*/,0/*acceleration*/);
circularTime(0/*theta*/,0/*Omega*/);
periodicDisplacement(0/*Omega*/,0/*time*/);
omega(0/*theta*/,0/*time*/);
frequency(0/*Omega*/);
angularVelocity(0/*theta*/,0/*time*/,0/*frequency*/);
linearAcceleration(0/*Omega*/,0/*time*/);
tangentialVelocity(0/*Omega*/,0/*height*/);
torque(0/*force*/,0/*height*/);
torqueWithTheta(0/*force*/,0/*height*/,0/*theta*/);
frictionforce(0/*u*/,0/*fn*/);
frictionCoefficient(0/*Frictionforce*/,0/*fn*/);
gravitationalForce(0/*Mass*/);
generalLawOfAttraction(0/*Mass1*/,0/*Mass2*/,0/*height*/);
kineticEnergy(0/*Mass*/,0/*velocity*/);
energySituation(0/*work*/,0/*height*/);
hookeIsLaw(0/*k*/,0/*x*/);
periodicTime(0/*Omega*/);
power(0/*work*/,0/*time*/);
conversionFromHorsepowerToJoules(0/*hp*/);
graphicAbility(0/*fp*/,0/*Bp*/);


return 0;

}
//Speed calculation
double velocity(double distance,double time)
{
double v = (distance/time);

return v;
}
//time calculation
double time(double distance,double velocity)
{
double t = (distance/velocity);

return t;
}
//Distance calculation
double distance(double velocity,double time)
{
double d = (velocity*time);

return d;
}
//Calculating acceleration
double acceleration(double velocity,double time)
{
double a = (velocity/time);

return a;
}
//volume
double volume(double length,double Height,double width)
{
double v = (length*Height*width);

return v;
}
//density
double density(double mass,double volume)
{
double p = (mass/volume);

return p;
}
//force
double force(double mass,double acceleration)
{
double f = (mass*acceleration);

return f;
}
//work
double work(double force,double distance)
{
double w = (force*distance);

return w;
}
//work with cos theta
double workTheta(double force,double distance,double theta)
{
double w = (force*distance*cos(theta));

return w;
}
//The resultant of two forces between them is theta
double resultantOftwoForces(double force1,double force2,double theta)
{
double R = sqrt(pow(force1,2)+pow(force2,2)+2*force1*force2*cos(theta));

return R;
}
//The direction of the resultant force between two forces is theta angle
double directionOfTheResultant(double force1,double force2,double theta)
{
double a = atan((force2*sin(theta))/(force1+force2*cos(theta)));

return a;
}
//Force analysis in two orthogonal directions
double forceAnalysisInTwoOrthogonalDirections(double force,double theta)
{

    double Rx,Ry;
    double  fx = force*cos(theta);
    Rx = Rx + fx;
    double   fy = force*sin(theta);
    Ry = Ry + fy;


double R = sqrt(pow(Rx,2)+pow(Ry,2));

return R;
}
//The direction of the resultant two-force analysis
double directionOfTheResultantTwoForce(double Rx,double Ry)
{

double a = atan((Ry)/(Rx));

return a;
}
//gravitational constant
double g(void)
{
  double  g = 9.8;

 return g;
}
//linear distance
double linearDistance(double velocity,double time,double acceleration)
{

double s = (velocity*time+(1/2)*acceleration*pow(time,2));

return s;
}
//linear time
double linearTime(double velocity,double initialvelocity,double acceleration)
{

double t = (initialvelocity-velocity)/(acceleration);

return t;
}
//circular time
double circularTime(double theta,double Omega)
{

double t = (theta)/(Omega);

return t;
}
//periodic displacement
double periodicDisplacement(double Omega,double time)
{

double th = (Omega)*(time);

return th;
}
//periodic speed
double omega(double theta,double time)
{


double w = (theta)/(time);

return w;
}
//frequency
double frequency(double Omega)
{


double f = (Omega)/(2*3.141592654);


return f;
}
//Angular velocity
double angularVelocity(double theta,double time,double frequency)
{


double w = (theta)/(time)*((2*3.141592654)*frequency);

return w;
}
//linear wheel
double linearAcceleration(double Omega,double time)
{


double at = (Omega)/(time);

return at;
}
//tangential velocity
double tangentialVelocity(double Omega,double length)
{


double vt = (Omega)*(length);

return vt;
}
//torque
double torque(double force,double length)
{


double Mo = (force)*(length);

return Mo;
}
//torque with theta
double torqueWithTheta(double force,double length,double theta)
{


double Mo = (force)*(length)*sin(theta);

return Mo;
}
//Friction force
double frictionforce(double u,double fn)
{


double Ff = (u)*(fn);

return Ff;
}
//friction coefficient
double frictionCoefficient(double Ff,double fn)
{


double u = (Ff)/(fn);

return u;
}
//gravitational force
double gravitationalForce(double Mass)
{


double w = (Mass)*(9.80665);

return w;
}
//The general law of attraction
double generalLawOfAttraction(double Mass1,double Mass2,double length)
{


double f = (6.67428*pow(10,-11))*((Mass1*Mass2)/(pow(length,2)));

return f;
}
//Kinetic energy
double kineticEnergy(double Mass,double velocity)
{


double ke = (1/2)*(Mass)*pow(velocity,2);

return ke;
}
//Energy situation
double energySituation(double work,double length)
{


double ke = (work)*(length);

return ke;
}
//HookeIsLaw
double hookeIsLaw(double k,double x)
{


double f = (-k)*(x);

return f;
}
//periodic time
double periodicTime(double Omega)
{


double T = (2*(3.141592654))/(Omega);

return T;
}
//Power
double power(double work,double time)
{


double P = (work)/(time);

return P;
}
//Conversion from horsepower to joules
double conversionFromHorsepowerToJoules(double hp)
{


double P = (hp)*(735);

return P;
}
//Graphic ability
double graphicAbility(double fp,double Bp)
{


double ip = (fp)+(Bp);

return ip;
}


