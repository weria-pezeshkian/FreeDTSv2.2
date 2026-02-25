#pragma once
#include <unordered_map>
#include <string>
#include <functional>
#include <istream>

class AbstractVertexPositionIntegrator;

/**
 * ===============================================================
 * VertexPositionIntegratorFactory
 * ===============================================================
 *
 * Singleton factory responsible for creating concrete
 * AbstractVertexPositionIntegrator objects at runtime based on
 * a string identifier.
 *
 * This enables:
 *  - Decoupling integrator selection from simulation code
 *  - Runtime configuration via input files
 *  - Easy extensibility without modifying existing logic
 *
 * Each derived integrator registers a creator function with
 * this factory, typically at static initialization time.
 */
class VertexPositionIntegratorFactory {
public:
    /**
     * Creator function type.
     *
     * Each creator:
     *  - Reads constructor parameters from an input stream
     *  - Returns a newly allocated derived integrator
     */
    using Creator =
        std::function<AbstractVertexPositionIntegrator*(std::istream&)>;

    /**
     * ---------------------------------------------------------------
     * Instance
     * ---------------------------------------------------------------
     * Returns the singleton instance of the factory.
     *
     * The static local instance ensures:
     *  - Exactly one factory exists
     *  - Thread-safe initialization (C++11 and later)
     */
    static VertexPositionIntegratorFactory& Instance() {
        static VertexPositionIntegratorFactory instance;
        return instance;
    }

    /**
     * ---------------------------------------------------------------
     * Register
     * ---------------------------------------------------------------
     * Registers a new integrator type with the factory.
     *
     * @param name    String identifier used in input files
     * @param creator Function that constructs the integrator
     */
    void Register(const std::string& name, Creator creator) {
        m_Creators[name] = creator;
    }

    /**
     * ---------------------------------------------------------------
     * Create
     * ---------------------------------------------------------------
     * Creates a concrete integrator based on its registered name.
     *
     * @param name   String identifier of the integrator type
     * @param input  Stream providing constructor parameters
     *
     * @return Pointer to a newly created integrator,
     *         or nullptr if the name is not registered.
     *
     * @note The caller is responsible for managing the lifetime
     *       of the returned object.
     */
    AbstractVertexPositionIntegrator* Create(
        const std::string& name,
        std::istream& input)
    {
        auto it = m_Creators.find(name);
        if (it == m_Creators.end())
            return nullptr;

        return it->second(input);
    }

private:
    /// Private constructor enforces singleton pattern
    VertexPositionIntegratorFactory() = default;

    /// Disable copy and assignment
    VertexPositionIntegratorFactory(const VertexPositionIntegratorFactory&) = delete;
    VertexPositionIntegratorFactory& operator=(const VertexPositionIntegratorFactory&) = delete;

    /// Registry mapping names to creator functions
    std::unordered_map<std::string, Creator> m_Creators;
};
