#include "Simulation.h"


static float pi = 3.14159;
static float g = 9.8;

inline static float power(float val, int power) {

	float value = 1;

	for (int i = 0; i < power; i++) {
		value *= val;
	}

	return value;
}

float Kernel::EvaluateKernel(glm::vec3& ri, glm::vec3& rj) const
{	
	float w;

	glm::vec3 rij = ri - rj;
	
	float r = glm::length(rij);

	if (type == kerneltype::poly6) {
		if (r >= 0 && r <= h) {
			
			float temp = h * h - r * r;

			
			w = 315.0f * power(temp, 3) / (64.0f * pi * power(h, 9));
			return w;

		}
		else {
			w = 0.0f;
			return w;
		}
	
	}

	else if (type == kerneltype::spiky) {
		if (r >= 0 && r <= h) {

			float temp = h - r;


			w = 15.0f * power(temp, 3) / (pi * power(h, 6));
			return w;

		}
		else {
			w = 0.0f;
			return w;
		}
	}

	else if (type == kerneltype::viscous) {
		if (r > 0 && r <= h) {
			if (r == 0) {
				return 0.0f;
			}

			float temp = r * r / (h * h) + h / (2.0f * r) - power(r,3) / (2.0f * power(h,3)) - 1.0f;


			w = 15.0f * temp / (2.0f * pi * power(h, 3));

			return w;

		}
		else {
			w = 0.0f;
			return w;
		}
	}
	else if (type == kerneltype::c2) {
		float q = r / h;

		if (q >= 0 && q <= 1) {

			float w = 7.0f / (4.0f * pi * h * h) * (1.0f - q) * (1.0f - q) * (1.0f - q) * (1.0f - q) * (1.0f + 4.0f * q);


			return w;

		}
		else {

			w = 0.0f;

			return w;
		}





	}

	return 404.0f;
}

glm::vec3 Kernel::EvaluateGradKernel(glm::vec3& ri, glm::vec3& rj) const
{
	float w;

	glm::vec3 rij = ri - rj;
	float r = glm::length(rij);

	if (type == kerneltype::spiky) {
		if (r >= 0 && r <= h) {

			w = -45.0f * (h - r) * (h - r) / (pi * power(h, 6));
			glm::vec3 W = w * glm::normalize(rij);
			return W;
		}
		else {
			w = 0.0f;
			glm::vec3 W = w * glm::normalize(rij);
			return W;
		}

	}
	if (type == kerneltype::c2) {
		float q = r / h;
		if (q >= 0 && q <= 1) {

			glm::vec3 w = -140.0f / (4.0f * pi * h * h * h) *  q * (1.0f - q) * (1.0f - q) * (1.0f - q) * rij / r;

			return w;
		}
		else {

			return glm::vec3(0.0f);
		}


	}

	return glm::vec3(4, 0, 4);
}

float Kernel::EvaluateLaplacianKernel(glm::vec3& ri, glm::vec3& rj) const 
{
	float w;
	glm::vec3 rij = ri - rj;

	float r = glm::length(rij);

	if (type == kerneltype::viscous) {
		if (r >= 0 && r <= h) {

			w = 45.0f * (h - r) / (pi * power(h, 6));
			return w;
		}

		else {
			w = 0.0f;
			return w;
		}

	}

	if (type == kerneltype::c2) {

		float q = r / h;

		if (q >= 0 && q <= 1) {

			float w = 7.0f / (4.0f * pi * h * h * h * h) * (-20.0f * (1.0f - q) * (1.0f - q) * (1.0f - 4.0f * q) + (-20.0f * (1.0f - q) * (1.0f - q) * (1.0f - q) / q));


			return w;

		}
		else {

			w = 0.0f;

			return w;
		}



	}

	return 404.0f;
	
}

void Simulation::ComputeHashParticle()
{
	for (int i = 0; i < m_Particles.size(); i++) {

		int ix = floor((m_Particles[i].position.x + m_parameters.width) / m_kernel.h);
		
		int iy = floor((m_Particles[i].position.y + m_parameters.height) / m_kernel.h);	

		ix = glm::clamp(ix, 0, nx - 1);
		iy = glm::clamp(iy, 0, ny - 1);

		m_Particles[i].hash = iy * nx  + ix;

	}

}

void Simulation::Sort()
{
	std::sort(m_Particles.begin(), m_Particles.end(), [](const Particle& a, const Particle& b) {
		return a.hash < b.hash;
		});



}

void Simulation::GetEndAndStart()
{

	int maxhash = m_Particles.back().hash;

	std::vector <int> CellEnd(nx * ny, -1);
	std::vector <int> CellStart(nx * ny, -1);
			
	int currentHash = m_Particles[0].hash;

	CellStart[currentHash] = 0;

	for (int i= 1; i < m_Particles.size(); i++) {
		
		int h = m_Particles[i].hash;


		if (h != currentHash) {

			CellEnd[currentHash] = i - 1;
			currentHash = m_Particles[i].hash;
			CellStart[currentHash] = i;

		}

	}

	CellEnd[currentHash] = m_Particles.size() - 1;
	 
	m_CellStart = CellStart;
	m_CellEnd = CellEnd;


}

std::vector<int> Simulation::NeighborSearch(int index)
{
	std::vector<int> neighbors;

	Particle& p = m_Particles[index];

	int ix = p.hash % nx;
	int iy = p.hash / nx;

	for (int dx = -1; dx <= 1; dx++) {
		for (int dy = -1; dy <= 1; dy++) {

			int nix = ix + dx;
			int niy = iy + dy;

			if (nix < 0 || nix >= nx || niy < 0 || niy >= ny)
				continue;

			int neighborHash = niy * nx + nix;

			if (neighborHash < 0 || neighborHash >= m_CellStart.size())
				continue;

			int start = m_CellStart[neighborHash];
			int end = m_CellEnd[neighborHash];
			
			if (start == -1 || end == -1) continue;

			for (int j = start; j <= end; j++) {
				if (j == index) continue;

				float dist = glm::distance(p.position, m_Particles[j].position);
				if (dist < m_kernel.h) {
					neighbors.push_back(j);
				}
			}
		}
	}

	return neighbors;
}


void Simulation::CalculateDensity()
{
	
	m_kernel.type = kerneltype::poly6;
	
	#pragma omp parallel for
	for (int i = 0; i < m_Particles.size(); i++) {
		
		float rho = 0;

		std::vector<int> neighbors = NeighborSearch(i);

		for (int j : neighbors) {
			
			
			float p = m_Particles[j].mass * m_kernel.EvaluateKernel(m_Particles[i].position, m_Particles[j].position);
			rho += p;

				

		}

		m_Particles[i].rho = rho * (1.0f - m_properties.beta * (m_Particles[i].temperature - m_properties.tempi));
		
		float B = m_properties.rhoi * m_properties.c * m_properties.c / m_properties.gamma;

		m_Particles[i].pressure = B * (power((rho / m_properties.rhoi), 7) - 1.0f);

		//m_Particles[i].pressure = m_k * (rho - m_properties.rhoi);

	}

	
}

void Simulation::CalculateTemperature()
{
	m_kernel.type = kerneltype::spiky;


	#pragma omp parallel for
	for (int i = 0; i < m_Particles.size(); i++) {

		
		float Temperature = 0.0f;
		std::vector<int> neighbors = NeighborSearch(i);

		for (int j : neighbors) {
			
			glm::vec3 rij = m_Particles[i].position - m_Particles[j].position;
			float r = glm::length(rij);
			if (r < 1e-6f) continue;

			if (i != j) {

				Temperature += m_Particles[j].mass * m_properties.kappa * (m_Particles[j].temperature - m_Particles[i].temperature) * (glm::dot(rij, m_kernel.EvaluateGradKernel(m_Particles[i].position, m_Particles[j].position))) / ((m_Particles[i].rho + m_Particles[j].rho) * (glm::length(rij) * glm::length(rij) + 1e-4f * m_kernel.h * m_kernel.h));

				


			}

			else {

			}

		}

		m_Particles[i].temperature += Temperature * m_parameters.dt;


	}
}

void Simulation::CalculatePressureForce()
{

	m_kernel.type = kerneltype::spiky;
	
		#pragma omp parallel for
		for (int i = 0; i < m_Particles.size(); i++) {

			glm::vec3 PressureForce(0.0f);
			std::vector<int> neighbors = NeighborSearch(i);

			for (int j: neighbors) {

				if (i != j) {

					//float scalar = m_Particles[j].mass * (m_Particles[i].pressure / (m_Particles[i].rho * m_Particles[i].rho)) + (m_Particles[j].pressure/ (m_Particles[j].rho * m_Particles[j].rho));

					float scalar = m_Particles[j].rho * m_Particles[j].mass * (m_Particles[i].pressure / m_Particles[i].rho + m_Particles[j].pressure / m_Particles[j].rho) ;

					glm::vec3 w = m_kernel.EvaluateGradKernel(m_Particles[i].position, m_Particles[j].position);

			
					glm::vec3 pf = scalar * w;

					PressureForce -= pf;
					

				}

				else {

				}

			}

			m_Particles[i].PressureForce = PressureForce;


		}


	
}

void Simulation::CalculateViscousForce()
{
	
	m_kernel.type = kerneltype::viscous;

	#pragma omp parallel for
	for (int i = 0; i < m_Particles.size(); i++) {

		glm::vec3 ViscousForce(0.0f);
		std::vector<int> neighbors = NeighborSearch(i);

		for (int j: neighbors) {

			if (i != j) {

				float scalar = m_Particles[j].mass * m_kernel.EvaluateLaplacianKernel(m_Particles[i].position, m_Particles[j].position) / m_Particles[j].rho;

				glm::vec3 v = m_Particles[j].Velocity - m_Particles[i].Velocity;


				glm::vec3 vf = scalar * v;

				ViscousForce += vf;


			}

			else {

			}

		}

		m_Particles[i].ViscousForce = m_properties.viscosity * ViscousForce;


	}

}

void Simulation::CalculateSurfaceForce()
{

	m_kernel.type = kerneltype::viscous;


	#pragma omp parallel for
	for (int i = 0; i < m_Particles.size(); i++) {
		float cs = 0;
		glm::vec3 n(0.0f);
		float kappa = 0;

		std::vector<int> neighbors = NeighborSearch(i);

		for (int j : neighbors) {


			cs = m_Particles[j].mass * m_kernel.EvaluateKernel(m_Particles[i].position, m_Particles[j].position) / m_Particles[j].rho;
			m_Particles[i].colorfield += cs;
			
			n = m_Particles[j].mass * m_kernel.EvaluateGradKernel(m_Particles[i].position, m_Particles[j].position) / m_Particles[j].rho;
			m_Particles[i].normal += n;

			kappa = - m_Particles[j].mass * m_kernel.EvaluateLaplacianKernel(m_Particles[i].position, m_Particles[j].position) / (m_Particles[j].rho *glm::length(n));
			m_Particles[i].curvature += kappa;



			m_Particles[i].SurfaceTension += m_properties.sigma * kappa * n;
		}



	}


}

void Simulation::CalculateExternalForce()
{
	#pragma omp parallel for
	for (int i = 0; i < m_Particles.size(); i++) {
		
		float fg = g * m_Particles[i].mass;
		m_Particles[i].ExternalForce = glm::vec3(0.0f, -fg, 0.0f);

	}
}

void Simulation::CalculateBouyantForce()
{
	for (int i = 0; i < m_Particles.size(); i++) {

		m_Particles[i].BouyantForce.y = m_Particles[i].rho * m_properties.beta * (m_Particles[i].temperature - m_properties.tempi) * g;

	}
}

void Simulation::ComputeAcceleration()
{
	#pragma omp parallel for
	for (int i = 0; i < m_Particles.size(); i++) {

		glm::vec3 pf = m_Particles[i].PressureForce;
		glm::vec3 vf = m_Particles[i].ViscousForce;
		glm::vec3 ef = m_Particles[i].ExternalForce;
		glm::vec3 bf = m_Particles[i].BouyantForce;
		

		m_Particles[i].Acceleration = (pf + vf + ef + bf) / m_Particles[i].mass;



	}


}

void Simulation::BoundaryBox()
{
	#pragma omp parallel for
	for (int i = 0; i < m_Particles.size(); i++) {

		if (m_Particles[i].position.x > m_parameters.width) {
			m_Particles[i].Velocity.x *= -0.2f;
			m_Particles[i].position.x = m_parameters.width - 0.1f;
		}
		if (m_Particles[i].position.x < -m_parameters.width) {
			m_Particles[i].Velocity.x *= -0.2f;
			m_Particles[i].position.x = -m_parameters.width + 0.1f;
		}

		if (m_Particles[i].position.y > m_parameters.height) {

			m_Particles[i].Velocity.y *= -0.2f;
			m_Particles[i].position.y = m_parameters.height - 0.1f;
		}
		if (m_Particles[i].position.y < -m_parameters.height) {

			m_Particles[i].Velocity.y *= -0.2f;
			m_Particles[i].position.y = -m_parameters.height + 0.1f;
		}
	}
}


static void SetupGrid(std::vector<Particle>& p, float spacing, int cols, int rows) {
	int index = 0;

	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < cols; x++) {
			if (index >= p.size()) return;

			p[index].position.x = -10.0f + x * spacing;
			p[index].position.y = 7.0f - y * spacing;
			index++;
		}
	}
}
static void SetupMass(std::vector<Particle>& p, float rhoi, float spacing) {

	for (int i = 0; i < p.size(); i++) {

		p[i].mass = rhoi * spacing * spacing;

	}


}

Simulation::Simulation(SimulationParameters Parameters, FluidProperties Properties, std::vector<Particle> Particles, Kernel kernel)
	: m_parameters(Parameters), m_properties(Properties), m_Particles(Particles), m_kernel(kernel), m_isRunning(true), m_interactive(false)
{
	SetupGrid(m_Particles, m_properties.spacing, m_properties.columns, m_properties.rows);
	SetupMass(m_Particles, m_properties.rhoi, m_properties.spacing);
}

Simulation::~Simulation()
{
}


void Simulation::step()
{
	if (m_isRunning) {

		ComputeHashParticle();
		Sort();
		GetEndAndStart();

		CalculateDensity();
		CalculateTemperature();
		CalculatePressureForce();
		CalculateViscousForce();
		//CalculateSurfaceForce();
		CalculateExternalForce();
		CalculateBouyantForce();
		ComputeAcceleration();
	

		for (int i = 0; i < m_Particles.size(); i++) {

			m_Particles[i].Velocity += m_parameters.dt  * m_Particles[i].Acceleration;

		
			m_Particles[i].position += m_parameters.dt  * m_Particles[i].Velocity;


		}



		BoundaryBox();
	}
	else
	{ 

	}
}

void Simulation::EnableBottomHeat()
{

	for (int i = 0; i < m_Particles.size(); i++) {

		if (m_Particles[i].position.y < -m_parameters.height * 0.98) {

			m_Particles[i].isBottom = true;
		}
		if (m_Particles[i].isBottom) {


			m_Particles[i].temperature = m_properties.tempOfBottom;
		}


		if (m_Particles[i].position.x < -m_parameters.width * 0.98 || m_Particles[i].position.x > m_parameters.width * 0.98) {


			m_Particles[i].temperature = m_properties.tempi;

		}
		

	}

}


void Simulation::AddForceAtCursor(float CursorPosx, float CursorPosy, bool mouseleft, bool mouseright)
{
	if (m_isRunning && m_interactive) {
		
		glm::vec3 cursor(CursorPosx, CursorPosy, 0.0f);
		
		float h = 6.0f;
	

		for (int i = 0; i < m_Particles.size(); i++) {

			float distance = glm::distance(cursor, m_Particles[i].position);

			if (distance < h) {
				if (mouseright) {

					m_Particles[i].Velocity += glm::normalize((cursor - m_Particles[i].position)) * 0.2f;
				}
				else if (mouseleft){

					m_Particles[i].Velocity -= glm::normalize((cursor - m_Particles[i].position)) * 0.2f;
				}
			}
			

		}

	}
}

void Simulation::reset()
{
	std::vector <Particle> temp(m_properties.NumberOfParticles);

	SetupGrid(temp, m_properties.spacing, m_properties.columns, m_properties.rows);
	SetupMass(temp, m_properties.rhoi, m_properties.spacing);

	m_Particles = temp;



}


