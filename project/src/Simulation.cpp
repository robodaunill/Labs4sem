#include "Simulation.hpp"
#include <random>
#include <glm/glm.hpp>

Simulation::Simulation(size_t particleCount) {
    particles.resize(particleCount);
    initializeParticles();
}

void Simulation::initializeParticles() {
    std::mt19937 gen(42);
    std::uniform_real_distribution<float> dist(-0.5f, 0.5f);
    
    for(auto& p : particles) {
        p.position = glm::vec3(dist(gen), 5.0f + dist(gen), dist(gen));
        p.velocity = glm::vec3(0.0f);
    }
}

void Simulation::update(float dt) {
    const glm::vec3 gravity(0.0f, -9.81f, 0.0f);
    
    for(auto& p : particles) {
        p.velocity += gravity * dt;
        p.position += p.velocity * dt;
        
        if(p.position.y < 0.0f) {
            p.position.y = 0.0f;
            p.velocity.y *= -0.8f;
        }
    }
}