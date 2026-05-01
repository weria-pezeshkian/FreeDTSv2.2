/*
 ===============================================================
 FactoryForceonVerticesfromInclusionsMethod
 ===============================================================

 This class implements a singleton factory for creating instances of
 `AbstractForceonVerticesfromInclusions` derived classes dynamically at runtime
 based on a string identifier.

 Each derived class constructor must take:
     (parameters...)


 ---------------------------------------------------------------
 How it works:
 ---------------------------------------------------------------
 1. Each derived class defines a static creator function:

        static AbstractForceonVerticesfromInclusions* Create(
            std::istream& input);

 2. The creator reads parameters from the input stream and constructs
    the object, passing the State pointer as the final argument.

 3. The derived class registers the creator with the factory.


 FactoryForceonVerticesfromInclusionsMethod::Instance()
            .Create(name, input);

 ---------------------------------------------------------------
*/

class State;   // forward declaration

class FactoryForceonVerticesfromInclusionsMethod {
public:

    /*
    ---------------------------------------------------------------
    Creator Function Type
    ---------------------------------------------------------------
    Every creator must accept:

        std::istream& input
    */
    using Creator =
        std::function<AbstractForceonVerticesfromInclusions*(std::istream&)>;

    /*
    ---------------------------------------------------------------
    Singleton Instance
    ---------------------------------------------------------------
    */
    static FactoryForceonVerticesfromInclusionsMethod& Instance()
    {
        static FactoryForceonVerticesfromInclusionsMethod instance;
        return instance;
    }

    /*
    ---------------------------------------------------------------
    Register
    ---------------------------------------------------------------
    Registers a creator function with a string identifier.
    */
    void Register(const std::string& name, Creator c)
    {
        m_Creators[name] = c;
    }

    /*
    ---------------------------------------------------------------
    Create
    ---------------------------------------------------------------
    Creates a derived evolution method.

    Parameters:
        name  : identifier of derived class
        input : stream containing parameters

    Returns:
        pointer to AbstractForceonVerticesfromInclusions or nullptr
    */
    AbstractForceonVerticesfromInclusions* Create(
        const std::string& name,
        std::istream& input)
    {
        auto it = m_Creators.find(name);

        if (it == m_Creators.end()){
            return nullptr;
        }

        return it->second(input);
    }

private:

    FactoryForceonVerticesfromInclusionsMethod() = default;

    FactoryForceonVerticesfromInclusionsMethod(
        const FactoryForceonVerticesfromInclusionsMethod&) = delete;

    FactoryForceonVerticesfromInclusionsMethod& operator=(
        const FactoryForceonVerticesfromInclusionsMethod&) = delete;

    std::unordered_map<std::string, Creator> m_Creators;
};
