#include <cmath>
#include <iostream>
#include <iomanip>

# include <cassert>

#include "sim.hpp"
#include "utils.hpp"

using namespace std;


Sim::Sim(size_t N, size_t Nmeasure, double Nthermal
	 , double stride, string const& filename, double temperature
	 , double kb, double J, double Jz, double Jxy, Vec H, bool ferroStart, std::string const& outstate
	 , std::string const& instate)
    : _N{N}, _Nmeasure{Nmeasure}, _Nthermal{static_cast<size_t>(Nthermal * N * N * N)}
    , _stride{static_cast<size_t>(stride * N * N * N)}, _out{filename}, _outstate{outstate}, _kbT{temperature * kb}
    , _J{J}, _Jz{Jz}, _Jxy{Jxy}, _H{H}, _energy{0.0}, _magnetization{0.0, 0.0, 0.0}
{      
    // Make N^3 spins
    _spins.resize(N * N * N);
    if(instate != ""){
	loadstate(instate);
    }else if(ferroStart) // Put all spins in the Z direction
	fill(begin(_spins), end(_spins), Vec{0, 0, 1});
    else // Generate random starting configuration
	generate(begin(_spins), end(_spins), [](){ return uniform_on_sphere(); });

    // Compute the energy and magnetization
    _energy = energy();
    _magnetization = magnetization();
}

Sim::Sim(Config const& params) :
    Sim(params.get<size_t>("N"), params.get<size_t>("Nmeasure"), params.get<double>("Nthermal")
	, params.get<double>("stride"), params.get<string>("filename"), params.get<double>("temperature")
	, params.get<double>("kb"), params.get<double>("J"), params.get<double>("Jz"), params.get<double>("Jxy"), params.get<Vec>("H"), params.get<bool>("ferroStart")
	, params.get<string>("outstate"), params.get<string>("instate")) {}

void Sim::run()
{
    // We may have to put this somewhere else later
    double V = _N * _N * _N;
    // Thermalization step
    
    cout << "Thermalization ..." << endl << endl;
    quietRun(_Nthermal);

    size_t progressStride = _Nmeasure / 20;

    // Run, skipping _stride steps each time to reduce correlation
    for(size_t i = 0 ; i < _Nmeasure ; ++i)
    {
	if(i % progressStride == 0)
	    cout << "iteration " << (_Nthermal + i * _stride) << " / " << _Nthermal + _Nmeasure * _stride << " (" << i * 100.0 / _Nmeasure<< " %)" << endl;
	quietRun(_stride);
	_out << i << " " << _energy << " " << printRaw(_magnetization);
	// Energy and magnetization per site
        _out << " " << _energy / V << " " << printRaw<double>((1 / V) * _magnetization);
	// Sum of each spin component square per site
	double Sx = 0, Sy = 0, Sz = 0;
	for (auto const& spin : _spins)
	{
            Sx += spin[0] * spin[0];
	    Sy += spin[1] * spin[1];
	    Sz += spin[2] * spin[2];
	}
	_out << " " << printRaw<double>((1 / V) * Vec{Sx, Sy, Sz}) << endl;
    }
    if(_outstate != "")
	storeState(_outstate);
}

void Sim::storeState(string const& filename) const
{
    ofstream o{filename};
    for(auto const& spin : _spins)
	o << spin;
}

void Sim::loadstate(string const& filename)
{
    ifstream i{filename};
    for(auto & spin : _spins)
	i >> spin;
}

// // SPOSTA I DATASET PRIMA DI CAMBIARLO 
// // MAKE DI NUOVO
// double vs = 2; // vertical scaling factor for the coupling in the z axis

double Sim::energy() const
{
    Vec magnetization{0, 0, 0};
    double interaction = 0;
    
    for(size_t i = 0 ; i < _N ; ++i)
    {
	for(size_t j = 0 ; j < _N ; ++j)
	{
	    for(size_t k = 0 ; k < _N ; ++k)
	    {
		auto const& s = _spins[index(i, j, k)];
		
		// Get half of all the nearest neighbors
		auto const& n1 = _spins[index(i + 1, j, k)];
		auto const& n2 = _spins[index(i, j + 1, k)];
		auto const& n3 = _spins[index(i, j, k + 1)]; // These are the neighbors on the z axis which would need differen coupling

		magnetization += s;
		interaction += ((_Jxy*(n1 + n2) + _Jz*n3) * s).sum(); // I might introduce the coupling here as a scaling factor multiplying n3
	    // interaction += ((n1 + n2 + n3) * s).sum();
        }
	}
    }
    
    return -_J * interaction - (_H * magnetization).sum();
}

Spin Sim::magnetization() const
{
    /* We cannot use _spins.sum() as the cluster has an older version
     * of GCC where this bug
     * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87641 is not
     * fixed */
    return accumulate(begin(_spins), end(_spins), Vec{0.0, 0.0, 0.0});
}

size_t Sim::index(int i, int j, int k) const
{
    int N = (int)_N;
    i = (i % N + N) % N;
    j = (j % N + N) % N;
    k = (k % N + N) % N;
    
    return _N * _N * i + _N * j + k;
}

size_t Sim::index(std::valarray<size_t> v) const
{
    return index(v[0], v[1], v[2]);
}

void Sim::quietRun(size_t Nsteps)
{
    // Assuming _energy is up to date
    // As well as _magnetization
    
    for(size_t i = 0 ; i < Nsteps ; ++i)
    {
	// Take a random site
        auto rs = random_site(_N);
	// Change its spin
	auto newSpin = uniform_on_sphere();
	// Compute the difference in energy
	auto delta = deltaE(rs, newSpin);
	// Accept or not ; if T=0 then only accept if we decrease the energy
	if(delta < 0 || (_kbT != 0 && (uniform_unit() <= exp(- (1.0 / _kbT) * delta))))
	{
	    // cout << "accept" << endl;
	    // Accept
	    _energy += delta;
	    _magnetization += (newSpin - _spins[index(rs)]);
	    Spin oldspin = _spins[index(rs)];
	    _spins[index(rs)] = newSpin;

	    // if(abs(_energy - energy()) > 1e-6)
	    // {
	    // 	cout << "Mistmatch : " << endl;
	    // 	cout << "site" << rs << endl;
	    // 	cout << "spin" << oldspin << endl;
	    // 	cout << "newSpin " << newSpin << endl;
	    // 	cout << "delta " << delta << endl;
	    // }
	}
	// else
	    // cout << "reject" << endl;
    }
}

double Sim::deltaE(valarray<size_t> site, Spin const& newSpin) const
{
    int i = site[0];
    int j = site[1];
    int k = site[2];

    // Get all nearest neighbors
    auto const& n1 = _spins[index(i + 1, j, k)];
    auto const& n2 = _spins[index(i, j + 1, k)];
    auto const& n3 = _spins[index(i, j, k + 1)]; // Upper plane
    auto const& n4 = _spins[index(i - 1, j, k)];
    auto const& n5 = _spins[index(i, j - 1, k)];
    auto const& n6 = _spins[index(i, j, k - 1)]; // Lower plane
    
    auto const& spin = _spins[index(i, j, k)];

    Vec spinSum = _Jxy*(n1 + n2 + n4 + n5) + _Jz*(n3 + n6); // n3 and n6 be multiplied by the same scaling factor as in energy funct
    // Vec spinSum = n1 + n2 + n3 + n4 + n5 + n6; before adding Jz
    Vec coeff = _J * spinSum + _H;
    return -(coeff * (newSpin - spin)).sum();
}
