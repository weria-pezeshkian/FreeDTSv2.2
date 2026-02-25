/*
===============================================================
FactoryInclusionConversionMethod
===============================================================

This class implements a **singleton factory** for creating instances of
`AbstractInclusionConversion` derived classes (e.g., active exchnage, equliburim exchange etc.) dynamically at runtime based on a string
identifier.

The factory allows decoupling of object creation from the rest of the code,
enabling polymorphic instantiation without hardcoding derived types.
This is particularly useful when reading simulation configurations from files
and creating the appropriate volume coupling objects on-the-fly.

---------------------------------------------------------------
How it works:
---------------------------------------------------------------
1. Each derived class defines a static "creator" function that reads
   the required parameters from an input stream and constructs an instance.
2. Each derived class registers itself with the factory using its default
   name and creator function. Typically, this is done using a static object
   at global scope.
3. When the simulation code apply Conversion object, it calls
   FactoryInclusionConversionMethod::Instance().Create(name, input),
   which:
   - Looks up the name in the registry.
   - If found, invokes the creator function to produce a new object.
   - Returns a pointer to the newly created AbstractInclusionConversion.
---------------------------------------------------------------
*/

class FactoryInclusionConversionMethod {
public:
    // Define the type of the creator function
    // It takes:
    //   - std::istream& input: input stream to read parameters
    // Returns a pointer to a newly constructed AbstractInclusionConversion
    using Creator = std::function<AbstractInclusionConversion*(std::istream&)>;

    /*
    ---------------------------------------------------------------
    Instance
    ---------------------------------------------------------------
    Returns a reference to the singleton instance of the factory.
    The static local variable ensures:
      - Only one instance exists during program execution.
      - Thread-safe initialization in C++11 and later.
    */
    static FactoryInclusionConversionMethod& Instance() {
        static FactoryInclusionConversionMethod instance;
        return instance;
    }

    /*
    ---------------------------------------------------------------
    Register
    ---------------------------------------------------------------
    Registers a new creator function with a given string identifier.

    Parameters:
      - name: The string name used to identify this derived class.
      - c: The creator function responsible for constructing objects of
           this type.

    Example:
      FactoryInclusionConversionMethod::Instance().Register(
          "ActiveTwoStateInclusion",
            ActiveTwoStateInclusion::Create
      );
    */
    void Register(const std::string& name, Creator c) {
        m_Creators[name] = c;
    }

    /*
    ---------------------------------------------------------------
    Create
    ---------------------------------------------------------------
    Creates a new AbstractInclusionConversion object based on a string name.

    Parameters:
      - name: The string identifier for the desired derived class.
      - input: Input stream from which to read constructor parameters.

    Returns:
      - A pointer to a newly created AbstractInclusionConversion, or nullptr
        if the name is not registered.

    Notes:
      - The caller is responsible for managing the lifetime of the returned
        object (typically using smart pointers is recommended).
      - This function decouples the code that creates objects from the
        knowledge of derived classes, enabling dynamic and extensible
        design.
    */
    AbstractInclusionConversion* Create(
        const std::string& name,
        std::istream& input)
    {
        auto it = m_Creators.find(name);
        if(it == m_Creators.end())
            return nullptr; // name not found in registry

        // Call the registered creator function to construct the object
        return it->second(input);
    }

private:
    // Private constructor to enforce singleton pattern
    FactoryInclusionConversionMethod() = default;

    // Disable copy constructor and assignment operator
    FactoryInclusionConversionMethod(const FactoryInclusionConversionMethod&) = delete;
    FactoryInclusionConversionMethod& operator=(const FactoryInclusionConversionMethod&) = delete;

    // Registry mapping string names to creator functions
    std::unordered_map<std::string, Creator> m_Creators;
};
