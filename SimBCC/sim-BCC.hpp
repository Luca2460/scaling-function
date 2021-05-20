#ifndef SIM_HPP
#define SIM_HPP

#include <valarray>
#include <random>
#include <fstream>
#include <string>

#include "utils.hpp"
#include "config.hpp"

class Sim
{
public: //edit below added ds
    Sim(size_t N, size_t Nmeasure, double Nthermal, double stride
	, std::string const& filename, double temperature, double kb, double J, double ds
	, Vec H, bool ferroStart, std::string const& outstate, std::string const& intstate);
    Sim(Config const& params);

    void run();

    double energy() const;
    Spin magnetization() const;

private:
    // Return the right index using periodic boundary conditions
    size_t index(int i, int j, int k) const;
    size_t index(std::valarray<size_t> v) const;

    // Run it without storing any informations
    void quietRun(size_t Nsteps);

    // Compute only the difference in energy
    double deltaE(std::valarray<size_t> site, Spin const& newSpin, int l) const; // ADDED l IN BCC VERSION

    // Store and load the state of the simulation
    void storeState(std::string const& filename) const;
    void loadState(std::string const& filename);
    
    // Simulation params
    size_t _N; // Size of our system along one dimension
    size_t _Nmeasure;
    size_t _Nthermal;
    size_t _stride;
    std::ofstream _out;
    std::string _outstate;

    // Physical params
    double _kbT; // Temperature (in energy units)
    double _J;
    double _ds; // differential scaling for spins in the same plane/different planes
    Vec _H;
    // double _Vcoupling // Coupling between different planes for tetragonal lattice
    
    Grid3D<Spin> _spins;

    double _energy;
    Vec _magnetization;
};

#endif
