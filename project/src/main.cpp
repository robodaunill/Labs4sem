#include "Simulation.hpp"
#include "VTKWriter.hpp"
#include <vulkan/vulkan.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <cstring>

#define VK_MESA_HEADLESS_EXTENSION_NAME "VK_MESA_headless"
#define VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME "VK_KHR_portability_enumeration"
#define VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR 0x00000001

class VulkanContext {
public:
    VulkanContext() { initialize(); }
    ~VulkanContext() { cleanup(); }

private:
    VkInstance instance;
    VkPhysicalDevice physicalDevice = VK_NULL_HANDLE;

    void initialize() {
        createInstance();
        selectPhysicalDevice();
    }

    void createInstance() {
        // Проверка доступных расширений
        uint32_t extCount = 0;
        vkEnumerateInstanceExtensionProperties(nullptr, &extCount, nullptr);
        std::vector<VkExtensionProperties> availableExts(extCount);
        vkEnumerateInstanceExtensionProperties(nullptr, &extCount, availableExts.data());
    
        // Собираем список доступных расширений
        std::vector<const char*> extensions;
        bool hasHeadless = false;
        bool hasPortability = false;
    
        for(const auto& ext : availableExts) {
            if(strcmp(ext.extensionName, VK_MESA_HEADLESS_EXTENSION_NAME) == 0) hasHeadless = true;
            if(strcmp(ext.extensionName, VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME) == 0) hasPortability = true;
        }
    
        if(hasPortability) extensions.push_back(VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME);
        if(hasHeadless) extensions.push_back(VK_MESA_HEADLESS_EXTENSION_NAME);
    
        VkInstanceCreateInfo createInfo{};
        createInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
        createInfo.enabledExtensionCount = static_cast<uint32_t>(extensions.size());
        createInfo.ppEnabledExtensionNames = extensions.data();
    
        if(hasPortability) {
            createInfo.flags |= VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR;
        }
    
        if(vkCreateInstance(&createInfo, nullptr, &instance) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create Vulkan instance!");
        }
    }

    void selectPhysicalDevice() {
        uint32_t deviceCount = 0;
        vkEnumeratePhysicalDevices(instance, &deviceCount, nullptr);
        
        std::vector<VkPhysicalDevice> devices(deviceCount);
        vkEnumeratePhysicalDevices(instance, &deviceCount, devices.data());

        for(const auto& device : devices) {
            VkPhysicalDeviceProperties props;
            vkGetPhysicalDeviceProperties(device, &props);
            
            if(props.deviceType == VK_PHYSICAL_DEVICE_TYPE_CPU) {
                physicalDevice = device;
                break;
            }
        }

        if(physicalDevice == VK_NULL_HANDLE) {
            throw std::runtime_error("Failed to find CPU device!");
        }
    }

    void cleanup() {
        vkDestroyInstance(instance, nullptr);
    }
};
int main() {
    try {
        VulkanContext vulkanContext;
        Simulation simulation(5000, 1.5f); // 5000 частиц, радиус 1.5
        
        const float dt = 0.016f;
        int frame = 0;
        
        auto startTime = std::chrono::high_resolution_clock::now();
        
        for(int i = 0; i < 300; ++i) {
            simulation.update(dt);
            
            if(i % 5 == 0) {
                VTKWriter::writeFrame(
                    simulation.getParticles(), 
                    "output/frame", 
                    frame++
                );
            }
        }
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        
        std::cout << "Simulation completed in " << duration.count() << " ms\n";
        std::cout << "Saved " << frame << " VTK frames\n";
        
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}