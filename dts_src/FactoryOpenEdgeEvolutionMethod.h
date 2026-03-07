/*
 ===============================================================
 FactoryOpenEdgeEvolutionMethod
 ===============================================================

 This class implements a singleton factory for creating instances of
 `AbstractOpenEdgeEvolution` derived classes dynamically at runtime
 based on a string identifier.

 Each derived class constructor must take:
     (parameters..., State*)

 The State pointer is always passed as the LAST constructor argument.

 ---------------------------------------------------------------
 How it works:
 ---------------------------------------------------------------
 1. Each derived class defines a static creator function:

        static AbstractOpenEdgeEvolution* Create(
            std::istream& input,
            State* state);

 2. The creator reads parameters from the input stream and constructs
    the object, passing the State pointer as the final argument.

 3. The derived class registers the creator with the factory.

 4. The State class calls:

        FactoryOpenEdgeEvolutionMethod::Instance()
            .Create(name, input, this);

 ---------------------------------------------------------------
*/

class State;   // forward declaration

class FactoryOpenEdgeEvolutionMethod {
public:

    /*
    ---------------------------------------------------------------
    Creator Function Type
    ---------------------------------------------------------------
    Every creator must accept:

        std::istream& input
        State* state

    The state pointer is always passed LAST so derived constructors
    follow the rule:

        Derived(...parameters..., State* state)
    */
    using Creator =
        std::function<AbstractOpenEdgeEvolution*(std::istream&, State*)>;

    /*
    ---------------------------------------------------------------
    Singleton Instance
    ---------------------------------------------------------------
    */
    static FactoryOpenEdgeEvolutionMethod& Instance()
    {
        static FactoryOpenEdgeEvolutionMethod instance;
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
        state : pointer to the State object creating the evolution

    Returns:
        pointer to AbstractOpenEdgeEvolution or nullptr
    */
    AbstractOpenEdgeEvolution* Create(
        const std::string& name,
        std::istream& input,
        State* state)
    {
        auto it = m_Creators.find(name);

        if (it == m_Creators.end())
            return nullptr;

        return it->second(input, state);
    }

private:

    FactoryOpenEdgeEvolutionMethod() = default;

    FactoryOpenEdgeEvolutionMethod(
        const FactoryOpenEdgeEvolutionMethod&) = delete;

    FactoryOpenEdgeEvolutionMethod& operator=(
        const FactoryOpenEdgeEvolutionMethod&) = delete;

    std::unordered_map<std::string, Creator> m_Creators;
};
