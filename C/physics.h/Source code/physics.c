#include<math.h>

//mechanics physics
////////////////////////////////////
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
////////////////////////////////////
//Electrical and magnetic physics
//////////////////////////////////
double coulombIsEqualToHowManyElectrons(double coulomb);
double electronsEqualToTheNumberOfCoulombs(double electrons);
double electricWork(double coulomb,double volt);
double coulomb(double work,double volt);
double current(double volt,double resistance);
double electricalPower(double volt,double current);
double volt(double current,double resistance);
double resistance(double volt,double current);
double conductorResistance(double conductiveMaterial,double height,double area);
double conductorMaterialType(double resistance,double height,double area);
double electricalConductivity(double conductorMaterialType);
double magneticField(double magneticFlux,double area);
double magneticFlux(double magneticField,double area);
double magneticFluxArea(double magneticFlux,double magneticField);
double magneticFluxDensity(double current,double numberOfTurns,double magneticPermeability,double RadiusLength);
double magneticFieldForce(double magneticField,double current,double height);
double magneticFieldForceWithTheta(double magneticField,double current,double height,double theta);
double fluxDensity(double magneticFieldForce,double current,double height,double theta);
double magneticMoment(double current,double area,double numberOfTurns);
double averageElectromotiveForceGeneratedByAChargedCoil(double magneticFlux,double time,double numberOfTurns);
double mutualInduction(double factorAffectingTheInductanceCoefficient,double current,double time);
double factorAffectingTheInductanceCoefficient(double mutualInduction,double current,double time);
double capacitorCapacitance(double coulomb,double volt);
double capacitorCharge(double absoluteInsulationConstant,double plateSpace,double TheDistanceBetweenThePanels);
double waveSpeed(double frequencie,double waveLength);
double ConvertTemperatureFromCelsiusToKelvin(double Celsius);
double FinnIslaw(double KelvinTemperature);
double bulbFilamentTemperature(double FinnIslaw);
double photonEnergy(double frequencie);
double PlanckIsConstant(void);
double conversionFromElectronVoltsToJoules(double ElectronVolts);
double conversionFromJoulesToElectronVolts(double joules);
double electronVolt(double mass,double velocity);
double futonSpeed(double mass,double energy);
double LawOfConservationOfMassAndEnergy(double mass,double TheSpeedOfLight);
double photonMomentum(double mass,double TheSpeedOfLight);
double PhotonEnergyInMotion(double Energy,double TheSpeedOfLight);
double waveLength(double PhotonEnergyInMotion);
double CoulombIsLaw(double Coulomb1,double Coulomb2,double length);
///////////////////////////////////////////////////////


int main()
{
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
linearDistance(0/*velocity*/,0/*time*/,0/*acceleration*/);
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
coulombIsEqualToHowManyElectrons(0/*coulomb*/);
electronsEqualToTheNumberOfCoulombs(0/*electrons*/);
electricWork(0/*coulomb*/,0/*volt*/);
coulomb(0/*work*/,0/*volt*/);
current(0/*volt*/,0/*resistance*/);
electricalPower(0/*work*/,0/*time*/);
volt(0/*current*/,0/*resistance*/);
resistance(0/*volt*/,0/*current*/);
conductorResistance(0/*conductiveMaterial*/,0/*height*/,0/*space*/);
conductorMaterialType(0/*resistance*/,0/*height*/,0/*space*/);
electricalConductivity(0/*conductorMaterialType*/);
magneticField(0/*magneticFlux*/,0/*space*/);
magneticFlux(0/*magneticField*/,0/*space*/);
magneticFluxArea(0/*magneticFlux*/,0/*magneticField*/);
magneticFluxDensity(0/*current*/,0/*numberOfTurns*/,0/*magneticPermeability*/,0/*RadiusLength*/);
magneticFieldForce(0/*magneticField*/,0/*current*/,0/*height*/);
magneticFieldForceWithTheta(0/*magneticField*/,0/*current*/,0/*height*/,0/*theta*/);
fluxDensity(0/*magneticFieldForce*/,0/*current*/,0/*height*/,0/*theta*/);
magneticMoment(0/*current*/,0/*area*/,0/*numberOfTurns*/);
averageElectromotiveForceGeneratedByAChargedCoil(0/*magneticFlux*/,0/*time*/,0/*numberOfTurns*/);
mutualInduction(0/*factorAffectingTheInductanceCoefficient*/,0/*current*/,0/*time*/);
factorAffectingTheInductanceCoefficient(0/*mutualInduction*/,0/*current*/,0/*time*/);
capacitorCapacitance(0/*coulomb*/,0/*volt*/);
capacitorCharge(0/*absoluteInsulationConstant*/,0/*plateSpace*/,0/*TheDistanceBetweenThePanels*/);
waveSpeed(0/*frequencie*/,0/*waveLength*/);
ConvertTemperatureFromCelsiusToKelvin(0/*Celsius temperature*/);
FinnIslaw(0/*Kelvin Temperature*/);
bulbFilamentTemperature(0/*FinnIslaw*/);
photonEnergy(0/*frequencie*/);
PlanckIsConstant();
conversionFromElectronVoltsToJoules(0/*ElectronVolts*/);
conversionFromJoulesToElectronVolts(0/*joules*/);
electronVolt(0/*mass*/,0/*velocity*/);
futonSpeed(0/*mass*/,0/*energy*/);
LawOfConservationOfMassAndEnergy(0/*mass*/,0/*TheSpeedOfLight*/);
photonMomentum(0/*mass*/,0/*TheSpeedOfLight*/);
PhotonEnergyInMotion(0/*Energy*/,0/*TheSpeedOfLight*/);
waveLength(0/*PhotonEnergyInMotion*/);
CoulombIsLaw(0/*Coulomb1*/,0/*Coulomb2*/,0/*length*/);



    return 0;
}


//mechanics physics
///////////////////////////////////////////
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
////////////////////////////////////
//Electrical and magnetic physics
///////////////////////////////////
//A coulomb is equal to how many electrons
double coulombIsEqualToHowManyElectrons(double coulomb)
{
double c = (coulomb)*(6.25*pow(10,18));

return c;
}
//Electrons equal to the number of coulombs
double electronsEqualToTheNumberOfCoulombs(double electrons)
{
double e = (electrons)*(1.6*pow(10,-19));

return e;
}
//The work done to move the electric charge
double electricWork(double coulomb,double volt)
{
double w = (coulomb)*(volt);

return w;
}
//coulomb
double coulomb(double work,double volt)
{
double c = (work)/(volt);

return c;
}
//current
double current(double volt,double resistance)
{
double I = (volt)/(resistance);

return I;
}
//Electrical Power
double electricalPower(double volt,double current)
{
double p = (volt)*(current);

return p;
}
//volt
double volt(double current,double resistance)
{
double v = (current)*(resistance);

return v;
}
//resistance
double resistance(double volt,double current)
{
double R = (volt)/(current);

return R;
}
//electrical resistance of the conductor
double conductorResistance(double conductiveMaterial,double height,double area)
{
double R = (conductiveMaterial)*((height)/(area));

return R;
}
//Conductor material type
double conductorMaterialType(double resistance,double height,double area)
{
double pe = (resistance)*((area)/(height));

return pe;
}
//electrical conductivity
double electricalConductivity(double conductorMaterialType)
{
double a = 1/(conductorMaterialType);

return a;
}
//the magnetic field
double magneticField(double magneticFlux,double area)
{
double B = (magneticFlux)/(area);

return B;
}
//magneticFlux
double magneticFlux(double magneticField,double area)
{
double fm = (magneticField)*(area);

return fm;
}
//magnetic flux area
double magneticFluxArea(double magneticFlux,double magneticField)
{
double A = (magneticFlux)/(magneticField);

return A;
}
//magnetic flux density
double magneticFluxDensity(double current,double numberOfTurns,double magneticPermeability,double RadiusLength)
{
double B = (magneticPermeability*numberOfTurns*current)/(2*RadiusLength);

return B;
}
//magnetic field force
double magneticFieldForce(double magneticField,double current,double height)
{
double Fb = (magneticField)*(current)*(height);

return Fb;
}
//Magnetic field strength with theta
double magneticFieldForceWithTheta(double magneticField,double current,double height,double theta)
{
double Fb = (magneticField)*(current)*(height)*(sin(theta));

return Fb;
}
//flux density
double fluxDensity(double magneticFieldForce,double current,double height,double theta)
{
double B = (magneticFieldForce)/((current)*(height)*(sin(theta)));

return B;
}
//magnetic moment
double magneticMoment(double current,double area,double numberOfTurns)
{
double Md = (current)*(area)*(numberOfTurns);

return Md;
}
//The average electromotive force generated by a charged coil
double averageElectromotiveForceGeneratedByAChargedCoil(double magneticFlux,double time,double numberOfTurns)
{
double emf = (-numberOfTurns)*((magneticFlux)/(time));

return emf;
}
//mutual induction
double mutualInduction(double factorAffectingTheInductanceCoefficient,double current,double time)
{
double emf2 = (-factorAffectingTheInductanceCoefficient)*((current)/(time));

return emf2;
}
//factorAffectingTheInductanceCoefficient
double factorAffectingTheInductanceCoefficient(double mutualInduction,double current,double time)
{
double M = (mutualInduction)/((current)/(time));

return M;
}
//capacitor capacitance
double capacitorCapacitance(double coulomb,double volt)
{
double C = (coulomb)/(volt);

return C;
}
//capacitor charge
double capacitorCharge(double absoluteInsulationConstant,double plateSpace,double TheDistanceBetweenThePanels)
{
double C = (absoluteInsulationConstant)*(plateSpace)/(TheDistanceBetweenThePanels);

return C;
}
//wave speed
double waveSpeed(double frequencie,double waveLength)
{
double V = (frequencie)*(waveLength);

return V;
}
//Convert temperature from Celsius to Kelvin
double ConvertTemperatureFromCelsiusToKelvin(double Celsius)
{
double K = (Celsius)+(273);

return K;
}
//Finn Is law
double FinnIslaw(double KelvinTemperature)
{
double lm = (2.898*pow(10,-3))/(KelvinTemperature);

return lm;
}
//bulb filament temperature
double bulbFilamentTemperature(double FinnIslaw)
{
double Tk = (2.898*pow(10,-3))/(FinnIslaw);

return Tk;
}
//photon energy
double photonEnergy(double frequencie)
{
double E = (6.625*pow(10,-34))*(frequencie);

return E;
}
//Planck's constant
double PlanckIsConstant(void)
{
double h = (6.625*pow(10,-34));

return h;
}
//Conversion from electron volts to joules
double conversionFromElectronVoltsToJoules(double ElectronVolts)
{
double J = (ElectronVolts)*(1.6*pow(10,-19));

return J;
}
//Conversion from joules to electron volts
double conversionFromJoulesToElectronVolts(double joules)
{
double Ve = (joules)*(1.25*pow(10,18));

return Ve;
}
//electron volt
double electronVolt(double mass,double velocity)
{
double Ve = (1/2)*(mass*pow(velocity,2));

return Ve;
}
//futon speed
double futonSpeed(double mass,double energy)
{
double c = (sqrt(energy/mass));

return c;
}
//Law of conservation of mass and energy
double LawOfConservationOfMassAndEnergy(double mass,double TheSpeedOfLight)
{
double E = (mass*pow(TheSpeedOfLight,2));

return E;
}
//photon momentum
double photonMomentum(double mass,double TheSpeedOfLight)
{
double Pl = (mass)*(TheSpeedOfLight);

return Pl;
}
//Photon energy in motion
double PhotonEnergyInMotion(double Energy,double TheSpeedOfLight)
{
double E = (Energy)/(pow(TheSpeedOfLight,2));

return E;
}
//wave length
double waveLength(double PhotonEnergyInMotion)
{
double lm = (6.625*pow(10,-34))/(PhotonEnergyInMotion);

return lm;
}
//Coulomb's law
double CoulombIsLaw(double Coulomb1,double Coulomb2,double length)
{


double f = (8.9875517873681764*pow(10,9))*((Coulomb1*Coulomb2)/(pow(length,2)));

return f;
}
