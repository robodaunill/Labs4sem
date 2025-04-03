#pragma once
#include <vector>
#include <glm/glm.hpp>

struct Particle {
    glm::vec3 position;
    glm::vec3 velocity;
};

class Simulation {
public:
    Simulation(size_t particleCount = 5000, float sphereRadius = 1.0f);
    void update(float dt);
    const std::vector<Particle>& getParticles() const { return particles; }

private:
    std::vector<Particle> particles;
    float sphereRadius;
    glm::vec3 sphereCenter = glm::vec3(0.0f, 5.0f, 0.0f); // Центр начальной сферы
    
    void initializeParticles();
    void checkCollisions();
};