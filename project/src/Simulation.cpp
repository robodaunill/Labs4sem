#include "Simulation.hpp"
#include <random>
#include <glm/gtx/norm.hpp> // Для glm::length2

Simulation::Simulation(size_t particleCount, float radius) 
    : sphereRadius(radius) {
    particles.resize(particleCount);
    initializeParticles();
}

void Simulation::initializeParticles() {
    std::mt19937 gen(42);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    
    for(auto& p : particles) {
        glm::vec3 point;
        float lengthSq;
        
        // Генерация точек внутри единичной сферы
        do {
            point = glm::vec3(dist(gen), dist(gen), dist(gen));
            lengthSq = glm::length2(point);
        } while(lengthSq > 1.0f || lengthSq < 0.01f); // Исключаем точки близко к центру
        
        // Масштабирование и смещение
        p.position = sphereCenter + point * sphereRadius;
        p.velocity = glm::vec3(0.0f);
    }
}

void Simulation::update(float dt) {
    const glm::vec3 gravity(0.0f, -9.81f, 0.0f);
    
    for(auto& p : particles) {
        p.velocity += gravity * dt;
        p.position += p.velocity * dt;
    }
    
    checkCollisions();
}

void Simulation::checkCollisions() {
    const float groundLevel = 0.0f;
    const float damping = 0.8f;
    
    for(auto& p : particles) {
        // Столкновение с полом
        if(p.position.y < groundLevel) {
            p.position.y = groundLevel;
            p.velocity.y *= -damping;
        }
        
        // Столкновение со сферой (сохранение формы)
        glm::vec3 toCenter = sphereCenter - p.position;
        float distance = glm::length(toCenter);
        
        if(distance > sphereRadius) {
            glm::vec3 normal = toCenter / distance;
            p.position = sphereCenter - normal * sphereRadius * 0.95f;
            p.velocity = glm::reflect(p.velocity, normal) * damping;
        }
    }
}