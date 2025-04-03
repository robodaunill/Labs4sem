#pragma once
#include <vector>
#include <glm/glm.hpp>

struct Particle {
    glm::vec3 position;
    glm::vec3 velocity;
    float pressure = 0.0f;
    float density = 1.0f;
};
class Simulation {
public:
    Simulation(size_t particleCount = 1000);
    void update(float dt);
    const std::vector<Particle>& getParticles() const { return particles; }

private:
    std::vector<Particle> particles;
    void initializeParticles();
};