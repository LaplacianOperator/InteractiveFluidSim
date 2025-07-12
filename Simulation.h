#pragma once
#include <vector>
#include <glm.hpp>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <omp.h>

struct Particle {

	glm::vec3 position;

	glm::vec3 normal;

	float rho;
	float mass;
	float pressure;
	float colorfield;
	float curvature;
	int hash;
	float temperature;

	glm::vec3 PressureForce;
	glm::vec3 ViscousForce;
	glm::vec3 ExternalForce;
	glm::vec3 SurfaceTension;
	glm::vec3 BouyantForce;
	glm::vec3 Acceleration;
	glm::vec3 Velocity;
	
	bool isBottom;
	bool isSurface;

	Particle(glm::vec3  Position = glm::vec3(0.0f, 0.0f, 0.0f), float Mass = 1.0f)
		: position(Position), rho(1000.0f), mass(Mass), pressure(0.0f), PressureForce(0.0f), ViscousForce(0.0f), Acceleration(0.0f), Velocity(glm::vec3(0.0f, 0.0f, 0.0f)), ExternalForce(0.0f), hash(0), SurfaceTension(glm::vec3(0.0f, 0.0f, 0.0f)), normal(glm::vec3(0.0f)), colorfield(0.0f), curvature(0.0f), temperature(0.5f), isBottom(false), isSurface(true), BouyantForce(glm::vec3(0.0f))
	{

	};
};

enum class kerneltype {
	poly6,
	viscous,
	spiky,
	c2
};

struct Kernel {
	
	kerneltype type;
	float h;

	
	Kernel(kerneltype Type, float H)
		: type(Type), h(H){}

	float EvaluateKernel(glm::vec3& ri, glm::vec3& rj) const;
	glm::vec3 EvaluateGradKernel(glm::vec3& ri, glm::vec3& rj) const;
	float EvaluateLaplacianKernel(glm::vec3& ri, glm::vec3& rj) const;
};

struct FluidProperties{
	
	int NumberOfParticles;
	float spacing;
	int rows;
	int columns;
	float viscosity;
	float rhoi;
	float tempi;
	float c;
	float gamma;
	float sigma;
	float threshold;
	float kappa;
	float tempOfBottom;
	float tempOfTop;
	float beta;

	FluidProperties(int n, float Spacing, int Rows, int Columns, float Viscosity, float Rhoi, float C, float Gamma, float Sigma, float Kappa)
		: NumberOfParticles(n), spacing(Spacing), rows(Rows), columns(Columns), viscosity(Viscosity), rhoi(Rhoi), c(C), gamma(Gamma), sigma(Sigma), kappa(Kappa), tempi(0.50f), tempOfBottom(0.50f), beta(0.10f), tempOfTop(0.10f), threshold(0.010f)
	{}

};

struct SimulationParameters {

	float dt;
	float t;
	float width;
	float height;
	
	SimulationParameters(float Dt, float T, float Width, float Height)
	: dt(Dt), t(T), width(Width), height(Height)
	{}

};

class Simulation
{
private:
	
	SimulationParameters m_parameters;
	FluidProperties m_properties;

	std::vector <Particle> m_Particles;
	Kernel m_kernel;


	bool m_isRunning;
	bool m_interactive;
 	
	// Spatial Hashing 
	int nx = ceil(2 * m_parameters.width / m_kernel.h);
	int ny = ceil(2 * m_parameters.height / m_kernel.h);

	std::vector <int> m_CellStart;
	std::vector <int> m_CellEnd;
	
	void ComputeHashParticle();
	void Sort();
	void GetEndAndStart();
	std::vector <int> NeighborSearch(int index);
	
	

	// Steps

	void CalculateDensity();
	void CalculateTemperature();
	void CalculatePressureForce();
	void CalculateViscousForce();
	void CalculateSurfaceForce();
	void CalculateExternalForce();
	void CalculateBouyantForce();
	void ComputeAcceleration();
	void BoundaryBox();


		



public:
	
	
	Simulation(SimulationParameters Parameters, FluidProperties Properties, std::vector <Particle> Particles, Kernel kernel);
	~Simulation();

	std::vector <Particle> get_particles() { return m_Particles; }
	
	void set_h(float h) { m_kernel.h = h; }

	void set_v(float v) { m_properties.viscosity = v; }
	void set_rho(float rho) { m_properties.rhoi = rho; }
	void set_sigma(float sigma) { m_properties.sigma = sigma; }
	void set_c(float c) { m_properties.c = c; }
	void set_t(float t) { m_properties.tempi = t; }
	void set_kappa(float k) { m_properties.kappa = k; }
	void set_tempOfBottom(float t) { m_properties.tempOfBottom = t; }
	void set_beta(float b) { m_properties.beta = b; }
	void set_threshold(float tr) { m_properties.threshold = tr; }
	void EnableBottomHeat();



	void step();
	void pause() { m_isRunning = false; }
	void unpause() { m_isRunning = true; }
	void set_interactive(bool interactive) { m_interactive = interactive; }

	void AddForceAtCursor(float CursorPosx, float CursorPosy, bool mouseleft, bool mouseright);
	void reset();


};